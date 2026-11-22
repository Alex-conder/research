# step6b_precision_coverage.py — 门控探针的 precision-coverage 操作点分析
# 产出：outputs{,_rpe1}/precision_coverage.json + outputs/precision_coverage.png
# 口径与 step6_run.py 完全一致（同特征、同标签、同 CV 划分种子），OOF 分数经
# cross_val_predict 获得，无泄漏。标签为中位数分割，故随机基线 precision=0.5。
import json
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features

SEED = 0
REPO = os.path.dirname(os.path.abspath(__file__))

DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1"},
]

COVERAGES = np.linspace(0.02, 1.0, 50)


def precision_coverage(scores, y, cov):
    """top-cov 分数覆盖带内，y=1（最差一半）的占比。"""
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def run_dataset(data_name, out_dir):
    out_path = os.path.join(REPO, out_dir)
    z = np.load(f"{out_path}/pred_expr_seed{SEED}.npz")
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()

    adata = ad.read_h5ad(
        f"{REPO}/data/gears_demo/{data_name}/perturb_processed.h5ad",
        backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    adata.file.close()
    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))

    pred_df = pd.read_csv(
        f"{out_path}/predictions_seed{SEED}.csv").set_index("perturbation")
    pred_df = pred_df.loc[perts]

    feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    feature_names = SPECTRUM_FEATURES + \
        [c for c in COVARIATES if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    y = (1.0 - pred_df["pearson_delta"].values
         > np.median(1.0 - pred_df["pearson_delta"].values)).astype(int)

    # OOF 分数（与 step6 同 CV 结构）
    pipe = LogisticRegression(C=1.0, max_iter=2000)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    Xs = StandardScaler().fit_transform(X)
    oof = cross_val_predict(pipe, Xs, y, cv=cv, method="predict_proba")[:, 1]

    curve = {f"{c:.2f}": precision_coverage(oof, y, c) for c in COVERAGES}

    def at(cv, c):
        key = f"{min(COVERAGES, key=lambda v: abs(v - c)):.2f}"
        return cv[key]
    # 单特征对照：最强单特征 delta_max 的 OOF 曲线
    j = feature_names.index("delta_max")
    oof_max = cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                                Xs[:, [j]], y, cv=cv,
                                method="predict_proba")[:, 1]
    curve_max = {f"{c:.2f}": precision_coverage(oof_max, y, c)
                 for c in COVERAGES}

    res = {
        "task": "precision-coverage operating-point analysis (step6b)",
        "dataset": data_name, "n": int(len(y)), "base_rate": 0.5,
        "probe": {"precision_at_20": at(curve, 0.20), "precision_at_33":
                  at(curve, 0.33), "precision_at_10": at(curve, 0.10),
                  "curve": curve},
        "delta_max_single": {"precision_at_20": at(curve_max, 0.20),
                             "precision_at_33": at(curve_max, 0.33),
                             "precision_at_10": at(curve_max, 0.10),
                             "curve": curve_max},
    }
    with open(f"{out_path}/precision_coverage.json", "w",
              encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[{data_name}] probe P@10/20/33 = "
          f"{res['probe']['precision_at_10']:.3f}/"
          f"{res['probe']['precision_at_20']:.3f}/"
          f"{res['probe']['precision_at_33']:.3f} (baseline 0.5); "
          f"delta_max P@20 = {res['delta_max_single']['precision_at_20']:.3f}")
    return res


results = [run_dataset(d["data"], d["out"]) for d in DATASETS]

# ---- 双面板图（英文，300 dpi）----
fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
titles = ["Norman 2019 (n=217)", "Replogle 2022 RPE1 (n=306)"]
for ax, res, ti in zip(axes, results, titles):
    cov = sorted(float(c) for c in res["probe"]["curve"])
    ax.plot(cov, [res["probe"]["curve"][f"{c:.2f}"] for c in cov],
            color="#1f77b4", lw=1.8, label="full probe (9 features)")
    ax.plot(cov, [res["delta_max_single"]["curve"][f"{c:.2f}"] for c in cov],
            color="#9467bd", lw=1.4, ls="--",
            label="best single feature (delta_max)")
    ax.axhline(0.5, color="gray", ls=":", lw=1.2, label="random gate (0.5)")
    ax.axvline(0.20, color="red", ls="-", lw=0.8, alpha=0.5)
    ax.annotate(f"P@20% = {res['probe']['precision_at_20']:.2f}",
                xy=(0.20, res["probe"]["precision_at_20"]),
                xytext=(0.30, res["probe"]["precision_at_20"] - 0.02),
                fontsize=9, color="#1f77b4")
    ax.set_xlabel("coverage (fraction of queries flagged)")
    ax.set_title(ti)
    ax.set_xlim(0, 1)
    ax.set_ylim(0.45, 1.0)
    ax.grid(alpha=0.25)
axes[0].set_ylabel("precision (fraction flagged that are worst-half)")
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle("Gate operating point: precision vs coverage "
             "(median-split labels, random baseline = 0.5)")
fig.tight_layout()
fig.savefig(f"{REPO}/outputs/precision_coverage.png", dpi=300)
print("[SAVED]", f"{REPO}/outputs/precision_coverage.png")
