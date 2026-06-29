#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
consensus_annotation.py - 多数据库融合打分

读取 adaptive_annotation.py 生成的 all_candidates.csv，对每个显著特征在
KEGG/PubChem/ChEBI 候选间打分，输出共识注释表与通路输入表。
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

SOURCE_WEIGHT = {"KEGG": 1.2, "ChEBI": 1.1, "PubChem": 1.0}
ADDUCT_BONUS = {
    "[M+H]+": 0.10, "[M-H]-": 0.10,
    "[M+Na]+": 0.05, "[M+K]+": 0.05, "[M+NH4]+": 0.05,
    "[M+Cl]-": 0.05, "[M+HCOO]-": 0.05, "[M+CH3COO]-": 0.05,
    "[M-H2O+H]+": 0.03,
}


def score_candidate(row: pd.Series, max_ppm: float = 50.0) -> float:
    """共识打分：来源权重 * ppm 接近度 + 加合物先验。"""
    src = row["Source"]
    w = SOURCE_WEIGHT.get(src, 1.0)
    ppm = abs(float(row["PPM_Diff"]))
    ppm_score = max(0.0, 1.0 - ppm / max_ppm)
    adduct_bonus = ADDUCT_BONUS.get(row["Adduct"], 0.0)
    return w * ppm_score + adduct_bonus


def build_consensus(candidates_path: Path, max_ppm: float = 50.0):
    cand = pd.read_csv(candidates_path)
    if cand.empty:
        return cand, cand.copy()
    cand["AbsPPM"] = cand["PPM_Diff"].abs()
    cand["Score"] = cand.apply(lambda r: score_candidate(r, max_ppm), axis=1)

    # 每个特征按来源保留最佳，再跨来源排序
    best_per_source = (
        cand.sort_values(["Feature_Mz", "Feature_RT", "Mode", "Source", "AbsPPM"])
        .groupby(["Feature_Mz", "Feature_RT", "Mode", "Source"], as_index=False)
        .first()
    )
    ranked = best_per_source.sort_values(
        ["Feature_Mz", "Feature_RT", "Mode", "Score", "AbsPPM"],
        ascending=[True, True, True, False, True],
    )
    top = ranked.groupby(["Feature_Mz", "Feature_RT", "Mode"], as_index=False).first()

    support = (
        cand.groupby(["Feature_Mz", "Feature_RT", "Mode"])["Source"]
        .nunique()
        .reset_index(name="Supporting_DBs")
    )
    top = top.merge(support, on=["Feature_Mz", "Feature_RT", "Mode"], how="left")
    top["Consensus_Level"] = top["Supporting_DBs"].apply(
        lambda n: "High" if n >= 3 else ("Medium" if n == 2 else "Single")
    )
    return cand, top


def to_cell_library(top_df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "Mz": top_df["Feature_Mz"],
        "RT": top_df["Feature_RT"],
        "Mode": top_df["Mode"],
        "Adduct": top_df["Adduct"],
        "Neutral_Mass": top_df["Neutral_Mass"],
        "Library_Match": top_df["Source"],
        "KEGG_ID": top_df.apply(lambda r: r["ID"] if r["Source"] == "KEGG" else "", axis=1),
        "Compound_Name": top_df["Name"],
        "Exact_Mass": top_df["Exact_Mass"],
        "PPM_Error": top_df["PPM_Diff"],
        "Formula": top_df["Formula"],
        "SMILES": top_df["SMILES"],
        "PubChem_CID": top_df["PubChem_CID"],
        "ChEBI_ID": top_df["ChEBI_ID"],
        "Consensus_Score": top_df["Score"],
        "Consensus_Level": top_df["Consensus_Level"],
        "Pathways": "",
    })


def main():
    parser = argparse.ArgumentParser(description="多数据库共识注释")
    parser.add_argument("--candidates", default="output/adaptive/all_candidates.csv", type=Path)
    parser.add_argument("--output-dir", default="output/consensus", type=Path)
    parser.add_argument("--output-library", default="library/cell_library_consensus.csv", type=Path)
    parser.add_argument("--max-ppm", type=float, default=50.0)
    args = parser.parse_args()

    if not args.candidates.exists():
        print(f"候选文件不存在：{args.candidates}", file=sys.stderr)
        sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cand, top = build_consensus(args.candidates, args.max_ppm)

    cand.to_csv(args.output_dir / "consensus_all_candidates_scored.csv", index=False)
    top.to_csv(args.output_dir / "consensus_annotation.csv", index=False)

    cell_lib = to_cell_library(top)
    cell_lib.to_csv(args.output_library, index=False)

    kegg_n = cell_lib["KEGG_ID"].replace("", pd.NA).dropna().nunique()
    pc_n = cell_lib["PubChem_CID"].replace("", pd.NA).dropna().nunique()
    cb_n = cell_lib["ChEBI_ID"].replace("", pd.NA).dropna().nunique()
    print(
        f"共识注释完成：{len(top)} 个特征，KEGG={kegg_n}, PubChem={pc_n}, ChEBI={cb_n}"
    )


if __name__ == "__main__":
    main()
