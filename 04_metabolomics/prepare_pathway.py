#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_pathway.py - 基于 KEGG ID 为显著特征准备通路富集输入文件

输入：
  - 本地库（build_cell_library.py 生成的 cell_library.csv，含 KEGG_ID）
  - POS/NEG 显著特征表（significant_hits.xlsx）
输出：
  - pathway_input.csv：MetaboAnalyst / KEGG Mapper 可直接使用的化合物- fold change 表
  - pathway_mapping.csv：每个特征的 KEGG_ID、KEGG Name、通路列表
  - pathway_summary.csv：每个通路覆盖的显著特征数
"""
import argparse
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

KEGG_LINK_URL = "https://rest.kegg.jp/link/pathway/{target}"
KEGG_GET_URL = "https://rest.kegg.jp/get/{entry}"


@retry(retry=retry_if_exception_type((requests.exceptions.RequestException,)),
       stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=15), reraise=True)
def kegg_link_pathways(kegg_id: str, timeout: int = 20):
    """查询某个 KEGG compound 所在的所有 pathway。"""
    if not kegg_id or not kegg_id.startswith("C"):
        return []
    url = KEGG_LINK_URL.format(target=f"cpd:{kegg_id}")
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    pathways = []
    for line in resp.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].startswith("path:"):
            pathways.append(parts[1].replace("path:", ""))
    return pathways


@retry(retry=retry_if_exception_type((requests.exceptions.RequestException,)),
       stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=15), reraise=True)
def kegg_get_name(kegg_id: str, timeout: int = 20):
    """通过 KEGG get 获取化合物名称。"""
    if not kegg_id or not kegg_id.startswith("C"):
        return ""
    url = KEGG_GET_URL.format(entry=f"cpd:{kegg_id}")
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    for line in resp.text.splitlines():
        if line.startswith("NAME"):
            return line.split(maxsplit=1)[1].split(";")[0].strip()
    return ""


def main():
    parser = argparse.ArgumentParser(
        description="Prepare pathway enrichment input from significant hits and KEGG library"
    )
    parser.add_argument("--lib", default="library/cell_library.csv", help="Local library with KEGG_ID")
    parser.add_argument("--pos-hits", default=None, help="POS significant hits Excel (optional)")
    parser.add_argument("--neg-hits", default=None, help="NEG significant hits Excel (optional)")
    parser.add_argument("--hits", default=None, help="Single hits Excel to process instead of pos/neg")
    parser.add_argument("--output-dir", default="output/pathway", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 加载库
    lib_path = Path(args.lib)
    if not lib_path.exists():
        logging.error(f"Library not found: {lib_path}. Please run build_cell_library.py first.")
        return 1
    lib = pd.read_csv(lib_path)
    if "KEGG_ID" not in lib.columns or "Mz" not in lib.columns or "RT" not in lib.columns:
        logging.error("Library must contain KEGG_ID, Mz, RT columns")
        return 1
    # 兼容不同库列名
    if "Compound_Name" in lib.columns and "Name" not in lib.columns:
        lib = lib.rename(columns={"Compound_Name": "Name"})
    if "Pathways" not in lib.columns:
        lib["Pathways"] = ""

    # 合并显著特征
    hits_list = []
    if args.hits:
        if Path(args.hits).exists():
            hits_list.append(pd.read_excel(args.hits))
        else:
            logging.error(f"Hits file not found: {args.hits}")
            return 1
    else:
        for p, label in [(args.pos_hits, "POS"), (args.neg_hits, "NEG")]:
            if p and Path(p).exists():
                hits_list.append(pd.read_excel(p))
            elif p:
                logging.warning(f"{label} hits not found: {p}")
    if not hits_list:
        logging.error("No significant hits files found")
        return 1
    hits = pd.concat(hits_list, ignore_index=True)

    # 按 Mz/RT 与库匹配，获取 KEGG_ID 与 Name
    hits["KEGG_ID"] = None
    hits["KEGG_Name"] = ""
    hits["Pathways"] = ""

    pathway_counter = defaultdict(list)  # pathway_id -> list of feature names
    pathway_cache = {}  # kegg_id -> pathways

    for idx, row in hits.iterrows():
        mz, rt = row["Mz"], row["RT"]
        # 按 Mz/RT 容差匹配
        match = lib[
            (abs(lib["Mz"] - mz) < 1e-4) &
            (abs(lib["RT"] - rt) < 0.01)
        ]
        if not match.empty:
            kegg_id = str(match.iloc[0]["KEGG_ID"])
            if not kegg_id or kegg_id.lower() in ("nan", "none", ""):
                continue
            kegg_name = str(match.iloc[0]["Name"]) if pd.notna(match.iloc[0]["Name"]) else ""
            lib_pw = str(match.iloc[0]["Pathways"]) if "Pathways" in match.columns and pd.notna(match.iloc[0]["Pathways"]) else ""
            hits.at[idx, "KEGG_ID"] = kegg_id
            hits.at[idx, "KEGG_Name"] = kegg_name
            logging.info(f"Mapping feature m/z={mz:.4f} RT={rt:.2f} -> {kegg_id} ({kegg_name})")

            if lib_pw:
                pathways = [p.strip() for p in lib_pw.split(";") if p.strip()]
            else:
                if kegg_id not in pathway_cache:
                    try:
                        pathway_cache[kegg_id] = kegg_link_pathways(kegg_id)
                        time.sleep(0.1)
                    except Exception as e:
                        logging.warning(f"KEGG link failed for {kegg_id}: {e}")
                        pathway_cache[kegg_id] = []
                pathways = pathway_cache[kegg_id]
            if pathways:
                hits.at[idx, "Pathways"] = ";".join(pathways)
                for pw in pathways:
                    pathway_counter[pw].append(kegg_name or f"m/z_{mz:.4f}")
            else:
                logging.info(f"  {kegg_id} has no linked pathways")
        else:
            logging.info(f"No KEGG match for feature m/z={mz:.4f} RT={rt:.2f}")

    # 输出 1：MetaboAnalyst / KEGG Mapper 输入
    # 第一列化合物名，第二列 Log2FC，第三列 KEGG_ID（可选）
    input_df = hits.copy()
    input_df["Compound_Name"] = input_df["KEGG_Name"].replace("", "Unknown")
    input_df = input_df[["Compound_Name", "Log2FC", "KEGG_ID"]]
    input_path = out_dir / "pathway_input.csv"
    input_df.to_csv(input_path, index=False, encoding="utf-8-sig")
    logging.info(f"[OK] Pathway input saved: {input_path}")

    # 输出 2：完整映射表
    mapping_path = out_dir / "pathway_mapping.csv"
    hits.to_csv(mapping_path, index=False, encoding="utf-8-sig")
    logging.info(f"[OK] Pathway mapping saved: {mapping_path}")

    # 输出 3：通路汇总
    if pathway_counter:
        summary_rows = []
        for pw_id, compounds in pathway_counter.items():
            summary_rows.append({
                "Pathway_ID": pw_id,
                "Compound_Count": len(compounds),
                "Compounds": ";".join(compounds),
            })
        summary_df = pd.DataFrame(summary_rows).sort_values("Compound_Count", ascending=False)
        summary_path = out_dir / "pathway_summary.csv"
        summary_df.to_csv(summary_path, index=False, encoding="utf-8-sig")
        logging.info(f"[OK] Pathway summary saved: {summary_path}")
        logging.info(f"成功映射 KEGG ID: {hits['KEGG_ID'].notna().sum()} / {len(hits)}")
        logging.info(f"涉及通路数: {len(summary_df)}")
    else:
        logging.warning("未找到任何通路链接")

    return 0


if __name__ == "__main__":
    sys.exit(main())
