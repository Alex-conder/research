#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neg_blankonly_fc.py - NEG Blank-only Fold Change 分析

由于 NEG 仅有一个 Blank Control，无法进行稳健的统计检验。
本脚本剔除 mixed std，以 Blank 为基准，计算 Treatment 样本的平均倍数变化：
    Log2FC = log2(mean(Treatment) + 1) - log2(Blank + 1)
筛选 |Log2FC| > 2 的特征作为“候选差异特征”，用于后续通路富集。
"""
import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def parse_mz_rt(name: str):
    m = re.search(r"M(\d+\.?\d*)_T(\d+\.?\d*)", name, re.I)
    if m:
        return float(m.group(1)), float(m.group(2))
    return np.nan, np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="output/feature_matrix_NEG_qcctrl_filtered.csv",
                        help="NEG filtered matrix CSV with Sample/Group columns")
    parser.add_argument("--output", default="output/neg_blankonly_fc_hits.xlsx",
                        help="Output Excel of FC-filtered hits")
    parser.add_argument("--fc", type=float, default=2.0, help="|Log2FC| threshold")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if "Sample" not in df.columns or "Group" not in df.columns:
        raise ValueError("Input CSV must contain 'Sample' and 'Group' columns")

    # 剔除 mixed std
    df = df[~df["Sample"].str.contains("std", case=False, na=False)].copy()

    blank_df = df[df["Sample"].str.contains("Blank00", case=False, na=False)]
    trt_df = df[~df["Sample"].str.contains("Blank00", case=False, na=False)]

    logging.info(f"Blank samples: {len(blank_df)}, Treatment samples: {len(trt_df)}")
    if len(blank_df) != 1:
        raise ValueError(f"Expected exactly 1 Blank, found {len(blank_df)}")
    if len(trt_df) < 1:
        raise ValueError("No Treatment samples found")

    meta_cols = ["Sample", "Group"]
    feature_cols = [c for c in df.columns if c not in meta_cols]

    blank_vals = blank_df[feature_cols].iloc[0].astype(float)
    trt_mean = trt_df[feature_cols].mean(axis=0).astype(float)

    log2fc = np.log2(trt_mean + 1.0) - np.log2(blank_vals + 1.0)

    records = []
    for feat in feature_cols:
        mz, rt = parse_mz_rt(feat)
        records.append({
            "Mz": mz,
            "RT": rt,
            "Log2FC": log2fc[feat],
            "Blank_Intensity": blank_vals[feat],
            "Treatment_Mean_Intensity": trt_mean[feat],
        })

    result = pd.DataFrame(records)
    sig = result[result["Log2FC"].abs() > args.fc].copy()
    sig = sig.sort_values("Log2FC", key=abs, ascending=False)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sig.to_excel(out_path, index=False, engine="openpyxl")

    print(f"[OK] NEG Blank-only FC analysis saved: {out_path}")
    print(f"    Total features: {len(result)}")
    print(f"    |Log2FC| > {args.fc}: {len(sig)}")
    print(f"    Up-regulated: {(sig['Log2FC'] > 0).sum()}")
    print(f"    Down-regulated: {(sig['Log2FC'] < 0).sum()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
