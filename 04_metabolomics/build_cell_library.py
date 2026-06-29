#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_cell_library.py - 基于显著特征从 KEGG API 构建细胞/哺乳动物代谢物本地库

说明：
- HMDB 官方并未提供公开稳定的 REST API（直接访问会被 Cloudflare 拦截）。
- 本脚本改用 KEGG REST API（学术使用）的 exact_mass 范围查询，获取与中性质量匹配的化合物。
- 支持多加合物并行查询（POS: [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+；NEG: [M-H]-, [M+Cl]-, [M+CH3COO]-）。
- 对每个特征，若多种加合物均命中，选择 |ppm_diff| 最小的匹配写入库；所有命中写入 expanded 文件。
- 生成的 CSV 可直接作为 run_analysis.py 的 --lib 参数使用（列：Name, Mz, RT, KEGG_ID, NeutralMass, Adduct, PPM_Diff）。
"""
import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

# 加合物转换：neutral_mass = observed_mz + shift
ADDUCT_SHIFTS = {
    "POS": {
        "[M+H]+": -1.00727646688,
        "[M+Na]+": -22.989218,
        "[M+K]+": -38.963158,
        "[M+NH4]+": -18.033823,
    },
    "NEG": {
        "[M-H]-": 1.00727646688,
        "[M+Cl]-": -34.968853,
        "[M+CH3COO]-": -59.013851,
    },
}

KEGG_GET_URL = "http://rest.kegg.jp/get/{entry}"


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    reraise=True,
)
def kegg_get_name(kegg_id: str, timeout: int = 10):
    """通过 KEGG get 获取化合物名称。"""
    url = KEGG_GET_URL.format(entry=f"cpd:{kegg_id}")
    try:
        resp = requests.get(url, timeout=timeout)
    except (requests.Timeout, requests.ConnectionError):
        raise
    if resp.status_code != 200:
        return kegg_id
    for line in resp.text.splitlines():
        if line.startswith("NAME"):
            try:
                return line.split(maxsplit=1)[1].split(";")[0].strip()
            except IndexError:
                return kegg_id
    return kegg_id


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    reraise=True,
)
def query_kegg_by_mass(neutral_mass: float, ppm: float = 5.0, timeout: int = 10):
    """
    通过 KEGG exact_mass 范围查询化合物。
    返回第一个匹配的 (KEGG_ID, Name, ExactMass)。
    """
    if neutral_mass <= 0:
        return None
    tol = neutral_mass * ppm / 1e6
    low, high = neutral_mass - tol, neutral_mass + tol
    url = f"http://rest.kegg.jp/find/compound/{low:.6f}-{high:.6f}/exact_mass"
    try:
        resp = requests.get(url, timeout=timeout)
    except (requests.Timeout, requests.ConnectionError):
        raise

    if resp.status_code != 200:
        logging.warning(f"KEGG HTTP {resp.status_code} for mass {neutral_mass:.6f}")
        return None

    text = resp.text.strip()
    if not text:
        return None

    # KEGG find 返回格式：cpd:ID\texact_mass
    first_line = text.splitlines()[0]
    parts = first_line.split("\t")
    if not parts:
        return None
    kegg_id = parts[0].strip().replace("cpd:", "")
    exact_mass = float(parts[1].strip()) if len(parts) > 1 else None
    name = kegg_get_name(kegg_id)
    return kegg_id, name, exact_mass


def process_hits(hits_path: Path, mode: str, ppm: float):
    """
    读取显著特征表，按模式下的所有加合物查询 KEGG。
    返回 (library_records, expanded_records)。
    """
    if not hits_path.exists():
        logging.warning(f"File not found: {hits_path}")
        return [], []

    df = pd.read_excel(hits_path)
    if "Mz" not in df.columns:
        raise ValueError(f"{hits_path} must contain 'Mz' column")

    adducts = ADDUCT_SHIFTS[mode]
    mz_values = df["Mz"].dropna().unique()
    total = len(mz_values)
    logging.info(f"Processing {mode}: {total} unique m/z values, adducts={list(adducts.keys())}")

    library_records = []
    expanded_records = []

    for i, mz in enumerate(mz_values, 1):
        mz = float(mz)
        candidates = []
        for adduct, shift in adducts.items():
            neutral = mz + shift
            if neutral <= 0:
                continue
            logging.info(f"[{mode}] ({i}/{total}) m/z={mz:.4f} {adduct} -> neutral={neutral:.6f}")
            try:
                result = query_kegg_by_mass(neutral, ppm=ppm)
            except Exception as e:
                logging.error(f"KEGG query failed for {adduct} m/z {mz}: {e}")
                result = None

            if result:
                kegg_id, name, exact_mass = result
                ppm_diff = (exact_mass - neutral) / neutral * 1e6 if exact_mass is not None else None
                # 取该 m/z 对应的所有 RT（可能有多个特征共享 m/z）
                for rt in df.loc[df["Mz"] == mz, "RT"].dropna().unique():
                    rec = {
                        "Name": name,
                        "Mz": mz,
                        "RT": float(rt),
                        "KEGG_ID": kegg_id,
                        "NeutralMass": neutral,
                        "Mode": mode,
                        "Adduct": adduct,
                        "PPM_Diff": ppm_diff,
                        "KEGG_ExactMass": exact_mass,
                    }
                    candidates.append(rec)
                    expanded_records.append(rec)
            time.sleep(0.1)

        if candidates:
            # 对每个特征，选择 |ppm_diff| 最小的匹配写入库
            best = min(candidates, key=lambda r: abs(r["PPM_Diff"]) if r["PPM_Diff"] is not None else 1e9)
            library_records.append(best)
            logging.info(f"[{mode}] Best match for m/z={mz:.4f}: {best['KEGG_ID']} ({best['Name']}) {best['Adduct']} ppm_diff={best['PPM_Diff']:.3f}")
        else:
            logging.info(f"[{mode}] No KEGG match for m/z={mz:.4f}")

    return library_records, expanded_records


def main():
    parser = argparse.ArgumentParser(
        description="Build a cell/mammalian metabolite library from significant hits via KEGG API (multi-adduct)"
    )
    parser.add_argument("--pos-hits", default="output/pos_analysis/significant_hits.xlsx",
                        help="POS significant hits Excel")
    parser.add_argument("--neg-hits", default="output/neg_analysis_qcctrl/significant_hits.xlsx",
                        help="NEG significant hits Excel")
    parser.add_argument("--output", default="library/cell_library.csv", help="Output library CSV")
    parser.add_argument("--expanded-output", default=None,
                        help="Optional CSV with all adduct-level KEGG matches for inspection")
    parser.add_argument("--ppm", type=float, default=5.0, help="Mass tolerance in ppm")
    args = parser.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_library = []
    all_expanded = []

    if Path(args.pos_hits).exists():
        lib, exp = process_hits(Path(args.pos_hits), "POS", args.ppm)
        all_library.extend(lib)
        all_expanded.extend(exp)
    else:
        logging.warning(f"POS hits not found: {args.pos_hits}")

    if Path(args.neg_hits).exists():
        lib, exp = process_hits(Path(args.neg_hits), "NEG", args.ppm)
        all_library.extend(lib)
        all_expanded.extend(exp)
    else:
        logging.warning(f"NEG hits not found: {args.neg_hits}")

    if not all_library:
        logging.warning("No KEGG matches found. Library not written.")
    else:
        lib_df = pd.DataFrame(all_library)
        lib_df.to_csv(out_path, index=False, encoding="utf-8-sig")
        logging.info(f"[OK] Library saved: {out_path} ({len(lib_df)} entries, {lib_df['KEGG_ID'].nunique()} unique KEGG IDs)")

    if args.expanded_output and all_expanded:
        exp_path = Path(args.expanded_output)
        exp_path.parent.mkdir(parents=True, exist_ok=True)
        exp_df = pd.DataFrame(all_expanded)
        exp_df.to_csv(exp_path, index=False, encoding="utf-8-sig")
        logging.info(f"[OK] Expanded matches saved: {exp_path} ({len(exp_df)} records)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
