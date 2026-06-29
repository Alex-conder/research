#!/usr/bin/env python3
"""
Filter a built feature matrix before differential analysis:
  - Remove features with too few detections
  - Remove isotope peaks (m/z ~ +1.003 Da, similar RT)
  - Keep top-N features by mean intensity
  - Save as CSV and SPSS .sav
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat


def remove_isotopes(meta: pd.DataFrame, mean_int: np.ndarray,
                    mz_tol: float = 0.005, rt_tol: float = 0.05) -> pd.DataFrame:
    """Remove isotope features, keeping the one with highest mean intensity."""
    df = meta.copy()
    df["mean_intensity"] = mean_int
    df = df.sort_values("Mz").reset_index(drop=True)
    keep = np.ones(len(df), dtype=bool)
    for i in range(len(df)):
        if not keep[i]:
            continue
        mz_i = df.at[i, "Mz"]
        rt_i = df.at[i, "RT"]
        int_i = df.at[i, "mean_intensity"]
        j = i + 1
        while j < len(df) and df.at[j, "Mz"] - mz_i <= 1.003 + mz_tol:
            if abs(abs(df.at[j, "Mz"] - mz_i) - 1.003) <= mz_tol and abs(df.at[j, "RT"] - rt_i) <= rt_tol:
                if df.at[j, "mean_intensity"] > int_i:
                    keep[i] = False
                    break
                else:
                    keep[j] = False
            j += 1
    return df.loc[keep].drop(columns="mean_intensity")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, help="feature_matrix_{MODE}.csv")
    parser.add_argument("--meta", required=True, help="feature_metadata_{MODE}.csv")
    parser.add_argument("--output-prefix", required=True, help="Output file prefix (no extension)")
    parser.add_argument("--min-samples", type=int, default=3, help="Min non-zero samples")
    parser.add_argument("--top-n", type=int, default=10000, help="Keep top-N features by mean intensity")
    parser.add_argument("--isotope-mz-tol", type=float, default=0.005)
    parser.add_argument("--isotope-rt-tol", type=float, default=0.05)
    args = parser.parse_args()

    mat = pd.read_csv(args.matrix, index_col=0)
    meta = pd.read_csv(args.meta)

    X = mat.drop(columns=["Group"])
    group = mat["Group"]

    if len(meta) != X.shape[1]:
        raise ValueError(f"Metadata rows ({len(meta)}) != matrix columns ({X.shape[1]})")

    # Detection filter
    detections = (X > 0).sum(axis=0)
    keep = detections >= args.min_samples
    X = X.loc[:, keep]
    meta = meta.loc[keep.values].reset_index(drop=True)
    print(f"After detection filter (>= {args.min_samples}): {X.shape[1]} features")

    # Isotope removal
    mean_int = X.mean(axis=0).values
    meta = remove_isotopes(meta, mean_int, args.isotope_mz_tol, args.isotope_rt_tol)
    # meta rows are a subset; align X columns by index position
    kept_idx = meta.index
    X = X.iloc[:, kept_idx]
    X.columns = [f"M{r['Mz']:.4f}_T{r['RT']:.4f}" for _, r in meta.iterrows()]
    print(f"After isotope removal: {X.shape[1]} features")

    # Top-N by mean intensity
    mean_int = X.mean(axis=0)
    top = mean_int.sort_values(ascending=False).head(args.top_n).index
    X = X[top]
    print(f"After top-{args.top_n}: {X.shape[1]} features")

    out = X.copy()
    out.insert(0, "Group", group)

    out_csv = Path(args.output_prefix + ".csv")
    out_sav = Path(args.output_prefix + ".sav")
    out.to_csv(out_csv, encoding="utf-8-sig")
    pyreadstat.write_sav(out, str(out_sav))
    meta.to_csv(args.output_prefix + "_meta.csv", index=False, encoding="utf-8-sig")
    print(f"Saved {out_csv} and {out_sav}")


if __name__ == "__main__":
    sys.exit(main())
