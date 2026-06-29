#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用 PubChem PUG REST API 稳健注释显著特征
（含重试、指数退避、速率限制、404 快速跳过）
"""
import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

REQUEST_INTERVAL = 0.5  # 秒，PubChem 建议不超过 5 req/s

# 常见加合物中性质量转换（ observed m/z -> neutral monoisotopic mass ）
ADDUCT_SHIFTS = {
    # 正离子模式
    "[M+H]+": -1.00727646688,
    "[M+Na]+": -22.989218,
    "[M+K]+": -38.963158,
    "[M-H2O+H]+": -19.017841,  # 脱水 + 质子
    "[M+NH4]+": -18.033823,
    # 负离子模式
    "[M-H]-": 1.00727646688,
    "[M+Cl]-": -34.968853,
    "[M+HCOO]-": -44.998214,
    "[M+CH3COO]-": -59.013851,
}


def parse_adducts(adduct_str: str):
    """把逗号分隔的加合物字符串解析为 (label, shift) 列表。"""
    if not adduct_str:
        return []
    items = [a.strip() for a in adduct_str.split(",") if a.strip()]
    out = []
    for a in items:
        if a in ADDUCT_SHIFTS:
            out.append((a, ADDUCT_SHIFTS[a]))
        else:
            raise ValueError(f"Unknown adduct '{a}'. Supported: {list(ADDUCT_SHIFTS.keys())}")
    return out


def default_adducts(mode: str):
    """根据模式返回默认加合物。"""
    mode = mode.lower()
    if mode == "pos":
        return [("[M+H]+", ADDUCT_SHIFTS["[M+H]+"])]
    if mode == "neg":
        return [("[M-H]-", ADDUCT_SHIFTS["[M-H]-"])]
    return []


class PubChemNoMatch(Exception):
    """Raised when PubChem returns 404 (no compounds in mass range)."""
    pass


@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=30),
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    reraise=True,
)
def query_pubchem_by_mass(neutral_mass: float, ppm_tol: float = 5.0, timeout: int = 15):
    """
    通过中性精确质量查询 PubChem 化合物，返回前 3 个匹配名称（IUPAC/Title）。
    若 404 则视为无匹配；其他 HTTP 错误会重试。
    """
    tol = neutral_mass * ppm_tol / 1e6
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"monoisotopic_mass/range/{neutral_mass - tol:.6f}/{neutral_mass + tol:.6f}/cids/JSON"
    )
    try:
        resp = requests.get(url, timeout=timeout)
    except (requests.Timeout, requests.ConnectionError):
        raise

    if resp.status_code == 404:
        # No compounds in this mass range
        return None
    if resp.status_code != 200:
        logging.warning(f"Unexpected HTTP {resp.status_code} for mass {neutral_mass:.6f}")
        raise requests.HTTPError(f"HTTP {resp.status_code}")

    data = resp.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None

    cid_str = ",".join(map(str, cids[:3]))
    prop_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"cid/{cid_str}/property/IUPACName,Title/JSON"
    )
    prop_resp = requests.get(prop_url, timeout=timeout)
    prop_resp.raise_for_status()
    prop_data = prop_resp.json()
    props = prop_data.get("PropertyTable", {}).get("Properties", [])
    names = [p.get("IUPACName") or p.get("Title") for p in props]
    names = [n for n in names if n]
    return "; ".join(names[:3]) if names else None


def annotate_hits(
    input_excel: Path,
    output_excel: Path,
    mode: str = "pos",
    adducts=None,
    ppm: float = 5.0,
    max_run: int = None,
):
    hits = pd.read_excel(input_excel)
    if "Mz" not in hits.columns:
        raise ValueError("Input Excel must contain 'Mz' column")

    if adducts is None:
        adducts = default_adducts(mode)
    if not adducts:
        raise ValueError("No adducts specified. Use --mode pos|neg or --adducts.")

    total = len(hits)
    n = min(total, max_run) if max_run else total
    annotations = []
    success = 0

    for idx in range(n):
        mz = hits.at[idx, "Mz"]
        if pd.isna(mz):
            annotations.append("")
            continue
        mz = float(mz)
        per_adduct = []
        for label, shift in adducts:
            neutral_mass = mz + shift
            if neutral_mass <= 0:
                continue
            logging.info(f"Querying ({idx + 1}/{n}): m/z={mz:.4f} {label} -> neutral={neutral_mass:.6f}")
            try:
                name = query_pubchem_by_mass(neutral_mass, ppm_tol=ppm)
            except Exception as e:
                logging.error(f"Failed after retries for {label} m/z {mz}: {e}")
                name = None
            if name:
                per_adduct.append(f"{label}:{name}")
            time.sleep(REQUEST_INTERVAL)

        if per_adduct:
            annotations.append(" | ".join(per_adduct))
            success += 1
        else:
            annotations.append("Unmatched")

    # Fill remaining rows if max_run < total
    if n < total:
        annotations += [""] * (total - n)

    hits["PubChem_Annotation"] = annotations
    hits.to_excel(output_excel, index=False, engine="openpyxl")
    logging.info(f"Saved annotated result: {output_excel}")
    logging.info(f"Successfully annotated {success}/{n} queried features")


def main():
    parser = argparse.ArgumentParser(
        description="PubChem API annotation for significant hits (adduct-aware)"
    )
    parser.add_argument("--hits", required=True, help="Significant hits Excel file")
    parser.add_argument("--output", default="significant_annotated_api.xlsx", help="Output Excel file")
    parser.add_argument("--mode", choices=["pos", "neg"], default="pos",
                        help="Ionization mode (default: pos -> [M+H]+)")
    parser.add_argument("--adducts", default=None,
                        help="Comma-separated adducts, e.g. '[M+H]+,[M+Na]+'")
    parser.add_argument("--ppm", type=float, default=5.0, help="Mass tolerance in ppm")
    parser.add_argument("--max", type=int, default=None, help="Max number of features to annotate")
    args = parser.parse_args()

    adducts = parse_adducts(args.adducts) if args.adducts else None
    annotate_hits(Path(args.hits), Path(args.output), mode=args.mode,
                  adducts=adducts, ppm=args.ppm, max_run=args.max)
    return 0


if __name__ == "__main__":
    sys.exit(main())
