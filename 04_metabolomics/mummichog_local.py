#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mummichog_local.py - 本地 m/z-to-pathway 兜底分析

不依赖化合物名称映射，直接用显著 m/z 的中性质量匹配 KEGG compound，
统计通路命中并使用超几何检验估计显著性。
"""
import argparse
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

from adaptive_annotation import ADDUCT_SHIFTS, APICache, query_kegg_by_mass_range

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def load_hits(pos_path: Path, neg_path: Path):
    dfs = []
    for path, mode in [(pos_path, "POS"), (neg_path, "NEG")]:
        if not path.exists():
            logging.warning("%s 不存在，跳过", path)
            continue
        df = pd.read_excel(path)
        df["Mode"] = mode
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_kegg_pathways(cpd_id: str, cache: APICache):
    info = cache.get_kegg_cpd(cpd_id)
    if info:
        return info.get("pathways", [])
    return []


def run_mummichog(hits: pd.DataFrame, ppm_tol: float, cache: APICache):
    grouped = hits.groupby(["Mz", "RT", "Mode"], as_index=False).first()
    feature_hits = []
    all_cpd_pathways = {}

    for mode, adducts in (("POS", list(ADDUCT_SHIFTS["POS"].keys())),
                          ("NEG", list(ADDUCT_SHIFTS["NEG"].keys()))):
        sub = grouped[grouped["Mode"] == mode]
        if sub.empty:
            continue
        for adduct in adducts:
            neutral_masses = sub["Mz"] + ADDUCT_SHIFTS[mode][adduct]
            nm_arr = neutral_masses[neutral_masses > 0]
            if nm_arr.empty:
                continue
            low = max(0.0, nm_arr.min() * (1 - ppm_tol / 1e6))
            high = nm_arr.max() * (1 + ppm_tol / 1e6)
            logging.info("%s %s: KEGG mass range %.4f - %.4f", mode, adduct, low, high)
            try:
                kegg_hits = query_kegg_by_mass_range(low, high, cache)
            except Exception as e:
                logging.warning("KEGG 查询失败 %s %s: %s", mode, adduct, e)
                continue
            for _, row in sub.iterrows():
                neutral = row["Mz"] + ADDUCT_SHIFTS[mode][adduct]
                if neutral <= 0:
                    continue
                best = None
                for h in kegg_hits:
                    ppm = (h["exact_mass"] - neutral) / neutral * 1e6
                    if abs(ppm) > ppm_tol:
                        continue
                    if best is None or abs(ppm) < abs(best["ppm"]):
                        best = {"cpd_id": h["cpd_id"], "ppm": ppm,
                                "exact_mass": h["exact_mass"], "adduct": adduct}
                if best:
                    pws = get_kegg_pathways(best["cpd_id"], cache)
                    all_cpd_pathways[best["cpd_id"]] = pws
                    feature_hits.append({
                        "Feature_Mz": row["Mz"],
                        "Feature_RT": row["RT"],
                        "Mode": mode,
                        "Adduct": adduct,
                        "KEGG_ID": best["cpd_id"],
                        "PPM_Error": best["ppm"],
                        "Pathways": ";".join(pws),
                    })

    hits_df = pd.DataFrame(feature_hits)
    if hits_df.empty:
        return hits_df, pd.DataFrame()

    # 通路统计
    pw_counter = Counter()
    pw_features = defaultdict(list)
    for rec in feature_hits:
        for pw in all_cpd_pathways.get(rec["KEGG_ID"], []):
            pw_counter[pw] += 1
            pw_features[pw].append(rec["KEGG_ID"])

    # 超几何检验背景：所有匹配到的 KEGG cpd
    bg_cpds = set(hits_df["KEGG_ID"].unique())
    M = len(bg_cpds)
    n = len(hits_df.drop_duplicates(["Feature_Mz", "Feature_RT", "Mode"]))
    rows = []
    for pw, k in pw_counter.items():
        N = len([c for c in bg_cpds if pw in all_cpd_pathways.get(c, [])])
        if N == 0:
            continue
        pval = hypergeom.sf(k - 1, M, N, n)
        rows.append({
            "Pathway_ID": pw,
            "Hits": k,
            "Pathway_Size_in_Background": N,
            "Background_Size": M,
            "Query_Size": n,
            "P_value": pval,
            "Hit_KEGG_IDs": ";".join(sorted(set(pw_features[pw]))),
        })
    pw_df = pd.DataFrame(rows).sort_values("P_value")
    if not pw_df.empty:
        from statsmodels.stats.multitest import multipletests
        pw_df["FDR_BH"] = multipletests(pw_df["P_value"], method="fdr_bh")[1]
    return hits_df, pw_df


def main():
    parser = argparse.ArgumentParser(description="Mummichog-style local m/z-to-pathway analysis")
    parser.add_argument("--pos-hits", default="output/pos_analysis/significant_hits.xlsx", type=Path)
    parser.add_argument("--neg-hits", default="output/neg_analysis_qcctrl/significant_hits.xlsx", type=Path)
    parser.add_argument("--output-dir", default="output/mummichog", type=Path)
    parser.add_argument("--cache-path", default="cache/adaptive_cache.sqlite", type=Path)
    parser.add_argument("--ppm", type=float, default=20.0)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = APICache(args.cache_path)
    hits = load_hits(args.pos_hits, args.neg_hits)
    if hits.empty:
        logging.error("未读取到显著特征")
        sys.exit(1)
    logging.info("显著特征 %d 行，唯一 m/z %d", len(hits), hits["Mz"].nunique())

    hits_df, pw_df = run_mummichog(hits, args.ppm, cache)
    hits_df.to_csv(args.output_dir / "mummichog_hits.csv", index=False)
    pw_df.to_csv(args.output_dir / "mummichog_pathways.csv", index=False)
    logging.info("Mummichog 完成：%d 个通路", len(pw_df))


if __name__ == "__main__":
    main()
