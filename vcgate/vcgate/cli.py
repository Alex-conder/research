"""vcgate CLI：给一批留出扰动的预测谱打门控分数。

示例：
    python -m vcgate oof \
        --pred-expr outputs/pred_expr_seed0.npz \
        --predictions outputs/predictions_seed0.csv \
        --h5ad data/gears_demo/replogle_k562_essential/perturb_processed.h5ad \
        --out gate_scores_seed0.csv

pred_expr npz 需含 pred_expr / pred_perturbations / pred_genes 三个键
（GEARS step1_run.py 的固化格式）。ctrl 均值从 --h5ad 的 ctrl 细胞算，
也可 --ctrl-mean ctrl.npy 直接给（跳过 anndata 依赖）。
"""
import argparse
import sys

import numpy as np
import pandas as pd

from vcgate.features import SPECTRUM_FEATURES, compute_features
from vcgate.probe import cv_gate_auroc, median_split_labels, oof_gate_scores


def _ctrl_mean_from_h5ad(h5ad_path, pred_genes):
    import anndata as ad
    adata = ad.read_h5ad(h5ad_path, backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    adata.file.close()
    if var_genes == list(pred_genes):
        return ctrl_mean
    m = dict(zip(var_genes, ctrl_mean))
    return np.array([m[g] for g in pred_genes])


def main(argv=None):
    p = argparse.ArgumentParser(prog="vcgate",
                                description="reliability gate for "
                                            "perturbation-response predictions")
    sub = p.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("oof", help="cross-validated out-of-fold gate scores")
    o.add_argument("--pred-expr", required=True, help="npz: pred_expr/perturbations/genes")
    o.add_argument("--predictions", required=True,
                   help="CSV with columns perturbation, pearson_delta, n_cells")
    o.add_argument("--h5ad", default=None, help="processed h5ad (ctrl mean source)")
    o.add_argument("--ctrl-mean", default=None, help="npy with ctrl mean (skips h5ad)")
    o.add_argument("--out", required=True, help="output CSV path")
    args = p.parse_args(argv)

    z = np.load(args.pred_expr)
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()

    if args.ctrl_mean:
        ctrl_mean = np.load(args.ctrl_mean)
    elif args.h5ad:
        ctrl_mean = _ctrl_mean_from_h5ad(args.h5ad, pred_genes)
    else:
        o.error("--h5ad or --ctrl-mean is required")

    pred_df = pd.read_csv(args.predictions).set_index("perturbation").loc[perts]
    feat = compute_features(pred_expr, perts, pred_genes, ctrl_mean)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    names = SPECTRUM_FEATURES + ["log_n_cells"]
    X = feat[names].values.astype(np.float64)

    error = 1.0 - pred_df["pearson_delta"].values
    y = median_split_labels(error)
    scores = oof_gate_scores(X, y)
    auroc_mean, auroc_std = cv_gate_auroc(X, y)

    out = pd.DataFrame({
        "perturbation": perts,
        "gate_score": scores,
        "error": error,
        "y_worst_half": y,
    }).sort_values("gate_score", ascending=False)
    out.to_csv(args.out, index=False)
    print(f"[vcgate] n={len(y)}, CV AUROC = {auroc_mean:.4f} ± {auroc_std:.4f}"
          f" (criterion > 0.6)")
    print(f"[vcgate] -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
