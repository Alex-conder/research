#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kegg_ora.py - 基于 KEGG pathway 的 Over-Representation Analysis（ORA）

输入：prepare_pathway.py 生成的 pathway_mapping.csv
输出：pathway_enrichment.csv（含 Fisher 精确检验 p 值 / BH-FDR）

使用规则：
- 将非空 KEGG_ID 的显著特征作为 query set；
- 通过 KEGG REST API 获取每个 pathway 所含化合物总数作为背景；
- 以参与分析的 pathway 化合物并集为 universe；
- 单尾 Fisher 精确检验，alternative='greater'。
"""
import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

KEGG_LINK_CPD_PW_URL = "https://rest.kegg.jp/link/cpd/{pathway_id}"
KEGG_GET_URL = "https://rest.kegg.jp/get/{entry}"


def fetch_pathway_compounds(pathway_id: str, timeout: int = 20):
    """返回某 pathway 中所有 cpd 的 KEGG ID 集合。"""
    url = KEGG_LINK_CPD_PW_URL.format(pathway_id=pathway_id)
    try:
        resp = requests.get(url, timeout=timeout)
    except Exception as e:
        logging.warning(f"Failed to fetch compounds for {pathway_id}: {e}")
        return set()
    if resp.status_code != 200:
        logging.warning(f"KEGG returned {resp.status_code} for {pathway_id}")
        return set()
    cpds = set()
    for line in resp.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].startswith("cpd:"):
            cpds.add(parts[1].replace("cpd:", ""))
    return cpds


def fetch_pathway_name(pathway_id: str, timeout: int = 10):
    """获取 pathway 名称。"""
    url = KEGG_GET_URL.format(entry=pathway_id)
    try:
        resp = requests.get(url, timeout=timeout)
    except Exception as e:
        logging.warning(f"Failed to fetch name for {pathway_id}: {e}")
        return ""
    if resp.status_code != 200:
        return ""
    for line in resp.text.splitlines():
        if line.startswith("NAME"):
            # NAME     pathway name - ...
            return line.split(maxsplit=1)[1].split(" - ")[0].strip()
    return ""


def main():
    parser = argparse.ArgumentParser(description="KEGG over-representation analysis")
    parser.add_argument("--mapping", required=True, help="Path to pathway_mapping.csv")
    parser.add_argument("--summary", default=None, help="Optional pathway_summary.csv (unused)")
    parser.add_argument("--output", default="output/pathway_enrichment.csv", help="Output CSV")
    args = parser.parse_args()

    mapping_path = Path(args.mapping)
    if not mapping_path.exists():
        logging.error(f"Mapping file not found: {mapping_path}")
        return 1

    df = pd.read_csv(mapping_path)
    if "KEGG_ID" not in df.columns or "Pathways" not in df.columns:
        logging.error("Mapping file must contain KEGG_ID and Pathways columns")
        return 1

    # 1. 构建 query set（非空 KEGG_ID，按唯一 ID 计算）
    query_df = df.dropna(subset=["KEGG_ID"]).copy()
    query_ids = {x for x in query_df["KEGG_ID"].astype(str) if x and x.lower() not in ("nan", "none")}
    logging.info(f"Query KEGG IDs: {len(query_ids)}")

    # 2. 收集涉及的 pathway
    pathway_to_hits = {}
    for _, row in query_df.iterrows():
        kegg_id = str(row["KEGG_ID"])
        pw_str = str(row["Pathways"]) if pd.notna(row["Pathways"]) else ""
        for pw in pw_str.split(";"):
            pw = pw.strip()
            if pw and not pw.startswith("path:"):
                pathway_to_hits.setdefault(pw, set()).add(kegg_id)

    pathway_ids = sorted(pathway_to_hits.keys())
    if not pathway_ids:
        logging.warning("No pathways found for query KEGG IDs")
        return 0

    logging.info(f"Pathways to evaluate: {len(pathway_ids)}")

    # 3. 获取每个 pathway 的化合物背景
    pathway_compounds = {}
    for pw in pathway_ids:
        cpds = fetch_pathway_compounds(pw)
        pathway_compounds[pw] = cpds
        logging.info(f"{pw}: {len(cpds)} compounds, {len(pathway_to_hits[pw])} hits")
        time.sleep(0.1)

    # 4. Universe：所有被评估 pathway 的化合物并集 + query set
    universe = set(query_ids)
    for cpds in pathway_compounds.values():
        universe |= cpds
    universe_size = len(universe)
    logging.info(f"Universe size: {universe_size}")

    # 5. Fisher exact test
    rows = []
    for pw in pathway_ids:
        hits = pathway_to_hits[pw]
        a = len(hits)
        b = len(query_ids) - a
        c = len(pathway_compounds[pw] - query_ids)
        d = universe_size - len(pathway_compounds[pw]) - b
        if a == 0:
            continue
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "Pathway_ID": pw,
            "Hits": a,
            "Query_Total": len(query_ids),
            "Pathway_Total": len(pathway_compounds[pw]),
            "Universe_Size": universe_size,
            "Odds_Ratio": odds,
            "P_value": pval,
        })

    if not rows:
        logging.warning("No pathways with hits")
        return 0

    result = pd.DataFrame(rows)
    # BH FDR
    _, fdr, _, _ = multipletests(result["P_value"].values, method="fdr_bh")
    result["FDR_BH"] = fdr

    # 获取 pathway 名称
    names = []
    for pw in result["Pathway_ID"]:
        names.append(fetch_pathway_name(pw))
        time.sleep(0.1)
    result["Pathway_Name"] = names

    # 重排列
    result = result[["Pathway_ID", "Pathway_Name", "Hits", "Query_Total",
                     "Pathway_Total", "Universe_Size", "Odds_Ratio", "P_value", "FDR_BH"]]
    result = result.sort_values("P_value")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, index=False, encoding="utf-8-sig")
    logging.info(f"[OK] Pathway enrichment saved: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
