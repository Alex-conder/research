# step6c_strata_robustness.py — 门控探针的分层鲁棒性/失效场景分析
# 按真值效应强度三分位与细胞数三分位分层，报每层 OOF AUROC。
# 回答审稿问题："门控在什么查询子集上失效？"（轻量分布偏移检验，
# 无需重训 GEARS；重训练版 multi-split 按预注册纪律放补充材料。）
# 产出：outputs{,_rpe1}/strata_robustness.json + outputs/strata_robustness.png
import json
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features

SEED = 0
REPO = os.path.dirname(os.path.abspath(__file__))
DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs", "tag": "K562"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1", "tag": "RPE1"},
]


def oof_auc(X, y, seeds=(0, 1, 2)):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    aucs = []
    for s in seeds:
        oof = cross_val_predict(LogisticRegression(C=1.0, max_iter=2000,
                                                   random_state=s),
                                X, y, cv=cv, method="predict_proba")[:, 1]
        aucs.append(float(roc_auc_score(y, oof)))
    return aucs


def tertile_labels(v):
    q = pd.qcut(pd.Series(v), 3, labels=["low", "mid", "high"])
    return q.astype(str).values


def run_dataset(data_name, out_dir, tag):
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
    X = StandardScaler().fit_transform(
        feat[feature_names].values.astype(np.float64))
    err = 1.0 - pred_df["pearson_delta"].values
    y = (err > np.median(err)).astype(int)

    res = {"task": "strata robustness / failure-mode analysis (step6c)",
           "dataset": tag, "n": int(len(y)),
           "overall_auroc_3seeds": oof_auc(X, y)}
    for cov_name, cov in [("true_strength_tertile",
                           tertile_labels(pred_df["strength"].values)),
                          ("n_cells_tertile",
                           tertile_labels(pred_df["n_cells"].values))]:
        strata = {}
        for lvl in ["low", "mid", "high"]:
            m = cov == lvl
            if m.sum() < 30 or len(np.unique(y[m])) < 2:
                strata[lvl] = {"n": int(m.sum()), "auroc": None}
                continue
            strata[lvl] = {"n": int(m.sum()),
                           "auroc": float(np.mean(oof_auc(X[m], y[m]))),
                           "worst_half_rate": float(y[m].mean())}
        res[cov_name] = strata
        line = "  ".join(
            f"{k}: n={v['n']}" +
            (f" AUROC={v['auroc']:.3f}" if v["auroc"] else " AUROC=NA")
            for k, v in strata.items())
        print(f"[{tag}] {cov_name}: {line}")
    with open(f"{out_path}/strata_robustness.json", "w",
              encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    return res


results = [run_dataset(d["data"], d["out"], d["tag"]) for d in DATASETS]

# ---- 图：分层 AUROC 条形图（英文，300 dpi）----
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), sharey=True)
for ax, res in zip(axes, results):
    labels, vals, ns = [], [], []
    overall = float(np.mean(res["overall_auroc_3seeds"]))
    labels.append("overall")
    vals.append(overall)
    ns.append(res["n"])
    for cov in ["true_strength_tertile", "n_cells_tertile"]:
        for lvl in ["low", "mid", "high"]:
            s = res[cov][lvl]
            labels.append(f"{cov.split('_')[1]}_{lvl}")
            vals.append(s["auroc"] if s["auroc"] else np.nan)
            ns.append(s["n"])
    colors = ["#1f77b4"] + ["#9ecae1"] * 3 + ["#c7e9c0"] * 3
    bars = ax.bar(range(len(vals)), vals, color=colors, edgecolor="gray",
                  lw=0.5)
    ax.axhline(0.5, color="red", ls=":", lw=1.2)
    ax.axhline(overall, color="#1f77b4", ls="--", lw=1.0, alpha=0.6)
    for b, v, n in zip(bars, vals, ns):
        if not np.isnan(v):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.01, f"{v:.2f}",
                    ha="center", fontsize=7)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylim(0.4, 1.0)
    ax.set_title(f"{res['dataset']} (n={res['n']})")
    ax.grid(axis="y", alpha=0.25)
axes[0].set_ylabel("OOF AUROC (3-seed mean)")
fig.suptitle("Gate robustness across query strata "
             "(dashed = overall; red dotted = chance)")
fig.tight_layout()
fig.savefig(f"{REPO}/outputs/strata_robustness.png", dpi=300)
print("[SAVED]", f"{REPO}/outputs/strata_robustness.png")
