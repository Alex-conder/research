#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_neg_blankonly.py - 生成仅保留 Blank 作为 Control 的 NEG 矩阵

原始 NEG 矩阵：
  - 1 个 Blank (Control)
  - 2 个 Treatment Sample
  - 2 个 mixed std (QC)

本脚本只保留 Blank + Treatment，剔除 mixed std，用于敏感性分析。
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import pyreadstat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="output/feature_matrix_NEG_qcctrl_filtered.csv",
                        help="NEG feature matrix CSV (with Sample column)")
    parser.add_argument("--output", default="output/feature_matrix_NEG_blankonly_filtered.sav",
                        help="Output SPSS .sav file")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if "Sample" not in df.columns:
        raise ValueError(f"{args.input} must contain a 'Sample' column")

    # 仅保留 Blank（Control）和真正的 Treatment Sample，剔除 mixed std
    keep_mask = ~df["Sample"].str.contains("std", case=False, na=False)
    filtered = df[keep_mask].copy()

    # 去掉 Sample 列，避免被当作特征
    filtered = filtered.drop(columns=["Sample"])

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pyreadstat.write_sav(filtered, str(out_path))
    print(f"[OK] Blank-only NEG matrix saved: {out_path} (samples={len(filtered)}, features={filtered.shape[1]-1})")
    print(f"   Groups: {filtered['Group'].unique().tolist()}")
    print(f"   Group counts: {filtered['Group'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
