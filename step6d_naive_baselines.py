# step6d_naive_baselines.py — 门控探针 vs 朴素不确定性基线
# 回答审稿问题："探针是否只是比随机/朴素方法好？"
# 基线1：|delta| 预测幅度本身（delta_l2）作为误差排序分（无训练）；
# 基线2：delta_max 单特征 logistic（与 single_feature_auroc 同口径复算）；
# 基线3：探针集成方差——同一探针用 5 个不同 CV 划分种子训练，
#        取逐查询 OOF 概率的方差作为置信分（lbfgs 求解器对
#        random_state 不敏感，故以划分种子作为集成多样性来源）；
# 主方法：现有 9 特征 gating probe（与 step6b 完全同结构）。
# 统一口径：OOF AUROC（高分=高误差）+ P@20（随机基线 0.5）。
# 产出：outputs{,_rpe1}/naive_baselines.json + outputs/naive_baselines.png
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
ENSEMBLE_SEEDS = [0, 1, 2, 3, 4]
REPO = os.path.dirname(os.path.abspath(__file__))
DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs", "tag": "K562"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1", "tag": "RPE1"},
]


def precision_at(scores, y, cov):
    """top-cov 分数覆盖带内，y=1（最差一半）的占比（与 step6b 同）。"""
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def load_dataset(data_name, out_dir):
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
    err = 1.0 - pred_df["pearson_delta"].values
    y = (err > np.median(err)).astype(int)
    return feat, feature_names, y, pred_df, out_path


def oof_scores(X, y, split_seed, single_col=None):
    """同一探针在不同 CV 划分种子下的 OOF 概率（标准化 + L2 logistic）。"""
    Xm = X if single_col is None else X[:, [single_col]]
    Xs = StandardScaler().fit_transform(Xm)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=split_seed)
    return cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                             Xs, y, cv=cv, method="predict_proba")[:, 1]


def run_dataset(data_name, out_dir, tag):
    feat, feature_names, y, pred_df, out_path = load_dataset(data_name, out_dir)
    X = feat[feature_names].values.astype(np.float64)

    # 主方法：9 特征探针（SEED=0 划分，与 step6b 同）
    probe = oof_scores(X, y, SEED)

    # 基线2：delta_max 单特征（同口径）
    j_max = feature_names.index("delta_max")
    delta_max = oof_scores(X, y, SEED, single_col=j_max)

    # 基线1：|delta| 幅度（delta_l2），无训练的纯排序分
    j_l2 = feature_names.index("delta_l2")
    magnitude = X[:, j_l2]  # 特征本身即 L2 范数，与 y 无标准化依赖

    # 基线3：5 个划分种子的探针集成方差（高分=高分歧=存疑）
    ens = np.column_stack([oof_scores(X, y, s) for s in ENSEMBLE_SEEDS])
    ens_var = ens.var(axis=1)

    methods = {
        "probe_full": {"scores": probe,
                       "desc": "9-feature gating probe (main)"},
        "delta_max_single": {"scores": delta_max,
                             "desc": "delta_max single-feature logistic"},
        "magnitude_l2": {"scores": magnitude,
                         "desc": "|delta| L2 norm, no training"},
        "ensemble_variance": {"scores": ens_var,
                              "desc": "variance of 5 split-seed probe OOF "
                                      "probabilities"},
    }
    res = {"task": "naive uncertainty baselines vs gating probe (step6d)",
           "dataset": tag, "n": int(len(y)), "base_rate": 0.5,
           "label": "y = (1 - pearson_delta) > median, worst half",
           "metrics": "OOF AUROC (higher score = more likely worst-half) "
                      "and P@20 (top-20% flagged, worst-half fraction)",
           "methods": {}}
    for name, m in methods.items():
        auroc = float(roc_auc_score(y, m["scores"]))
        res["methods"][name] = {
            "desc": m["desc"], "auroc": auroc,
            "p_at_20": precision_at(m["scores"], y, 0.20)}
    with open(f"{out_path}/naive_baselines.json", "w",
              encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    line = "  ".join(f"{k}: AUROC={v['auroc']:.3f} P@20={v['p_at_20']:.3f}"
                     for k, v in res["methods"].items())
    print(f"[{tag}] {line}")
    return res


def main():
    results = [run_dataset(d["data"], d["out"], d["tag"]) for d in DATASETS]

    # ---- 图：每数据集一组，方法为 x，AUROC 与 P@20 双指标（英文，300 dpi）----
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.7), sharey=True)
    order = ["probe_full", "delta_max_single", "magnitude_l2",
             "ensemble_variance"]
    labels = ["probe (9 feat.)", "delta_max\n(single)", "|delta| L2\n(no train)",
              "ensemble\nvariance"]
    titles = ["Replogle K562 (n=217)", "Replogle 2022 RPE1 (n=306)"]
    for ax, res, ti in zip(axes, results, titles):
        aucs = [res["methods"][k]["auroc"] for k in order]
        p20s = [res["methods"][k]["p_at_20"] for k in order]
        x = np.arange(len(order))
        w = 0.36
        b1 = ax.bar(x - w / 2, aucs, w, color="#1f77b4", edgecolor="gray",
                    lw=0.5, label="OOF AUROC")
        b2 = ax.bar(x + w / 2, p20s, w, color="#ff7f0e", edgecolor="gray",
                    lw=0.5, label="P@20")
        for b in list(b1) + list(b2):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.008,
                    f"{b.get_height():.2f}", ha="center", fontsize=7)
        ax.axhline(0.5, color="red", ls=":", lw=1.2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(ti)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("score (chance = 0.5)")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Gating probe vs naive uncertainty baselines "
                 "(same OOF protocol, median-split labels)")
    fig.tight_layout()
    fig.savefig(f"{REPO}/outputs/naive_baselines.png", dpi=300)
    print("[SAVED]", f"{REPO}/outputs/naive_baselines.png")


if __name__ == "__main__":
    main()
