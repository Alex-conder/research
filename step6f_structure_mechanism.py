# step6f_structure_mechanism.py — 门控误差"结构机制"双检验（纯计算，不重训模型）
# 回答审稿问题："探针探测到的误差结构是任务内在难度还是指标 artifact？"
#
# 检验一（解释 2：任务内在难度）：
#   对每个 held-out 基因，用"训练集逐基因平均 delta"（mean predictor，
#   即所有训练单扰动真值 delta 的逐基因均值）作最朴素基线，按与 step1/6
#   完全相同的口径（top-50 |true_delta| 基因上 Pearson）算误差；
#   算 mean-baseline 误差 vs GEARS 误差的逐基因 Spearman/Pearson 相关；
#   顺带给出 step6 归档 single_feature_auroc 对应特征与误差的 Spearman
#   相关矩阵一角。
#   判定线（本脚本约定）：|rho| > 0.5 → 结构部分来自任务内在难度；否则难度非主因。
#
# 检验二（解释 3：指标 artifact）：
#   对同一批 held-out 基因，用备选指标重算误差并排序：
#     all-Pearson   全基因 Pearson（predictions CSV 的 pearson_delta_all 列）
#     all-Spearman  全基因 Spearman（pred_delta vs true_delta，现算）
#     top50-MSE     top-50 |true_delta| 基因上的均方误差（现算）
#   与主口径误差排序算 Spearman 相关；并按 step6c 完全相同的分层协议
#   （strength 三分位、5-fold OOF、3 划分种子均值、层内重拟合）重算每层
#   探针 AUROC，检查 "Norman low/mid 强度层 AUROC 低于随机" 是否只在
#   top-DEG Pearson 口径下出现。
#
# 数据来源（全部现有归档）：outputs{,_rpe1}/pred_expr_seed0.npz、
#   predictions_seed0.csv、data/gears_demo/<data>/perturb_processed.h5ad。
# 产出：outputs{,_rpe1}/structure_mechanism.json + outputs/structure_mechanism.png
#   （verdict 文：outputs/structure_mechanism_verdict.md，由实跑数字人工撰写）
import json
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6d_naive_baselines import SEED, load_dataset

REPO = os.path.dirname(os.path.abspath(__file__))
DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs", "tag": "Norman"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1", "tag": "RPE1"},
]
TOP_N_DEG = 50                # 主误差口径（step1_run.py）
DIFFICULTY_RHO_THRESH = 0.5   # 检验一判定线（|rho| 超此线即难度为结构来源之一）
STRATA_SEEDS = (0, 1, 2)      # 与 step6c 相同
ERROR_VARIANTS = ["top50_Pearson", "all_Pearson", "all_Spearman",
                  "top50_MSE", "mean_baseline_top50_Pearson"]


def build_truth_and_mean_baseline(data_name, holdout_genes, pred_genes):
    """返回 (true_delta, mean_delta, ctrl_aligned)（均为 pred_genes 顺序）。

    true_delta:  (n_heldout, n_genes) 留出扰动组真值均值 - ctrl 训练均值
    mean_delta:  (n_genes,) 训练集所有单扰动真值 delta 的逐基因均值
                 （mean predictor：对任何查询都预测同一个 delta）
    """
    adata = ad.read_h5ad(
        f"{REPO}/data/gears_demo/{data_name}/perturb_processed.h5ad",
        backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))

    conditions = adata.obs.condition.unique().tolist()
    singles = [c for c in conditions if ("ctrl" in c and c != "ctrl")]
    genes_all = sorted({c.replace("ctrl+", "").replace("+ctrl", "")
                        for c in singles})
    genes_all = [g for g in genes_all if "_" not in g]
    missing = set(holdout_genes) - set(genes_all)
    assert not missing, f"留出基因不在 split 基因名单: {sorted(missing)[:5]}"
    holdout_set = set(holdout_genes)
    train_conds = [c for c in singles
                   if c.replace("ctrl+", "").replace("+ctrl", "")
                   not in holdout_set]
    assert len(train_conds) > 0

    # 训练集逐基因 mean delta（增量累加，避免整表 densify）
    acc = np.zeros(len(pred_genes), dtype=np.float64)
    for c in train_conds:
        m = np.asarray(adata[adata.obs.condition == c].X.mean(axis=0)).ravel()
        acc += m - ctrl_aligned
    mean_delta = acc / len(train_conds)

    # 留出扰动真值 delta（npz perts 均在扰动图内、n_cells>0，与 step1 同假设）
    true_rows = []
    for g in holdout_genes:
        cand = [c for c in singles if c.startswith(g + "+")
                or c.endswith("+" + g)]
        assert cand, f"留出基因 {g} 找不到 condition"
        m = np.asarray(
            adata[adata.obs.condition == cand[0]].X.mean(axis=0)).ravel()
        true_rows.append(m - ctrl_aligned)
    adata.file.close()
    return np.vstack(true_rows), mean_delta, ctrl_aligned, len(train_conds)


def top_deg_pearson(pred_delta, true_delta, n_top=TOP_N_DEG):
    """step1_run.py 口径：top-N |true_delta| 基因上 Pearson(pred_delta, true_delta)。"""
    idx = np.argsort(-np.abs(true_delta))[:n_top]
    return float(np.corrcoef(pred_delta[idx], true_delta[idx])[0, 1])


def all_gene_spearman(pred_delta, true_delta):
    return float(sps.spearmanr(pred_delta, true_delta)[0])


def top_deg_mse(pred_delta, true_delta, n_top=TOP_N_DEG):
    idx = np.argsort(-np.abs(true_delta))[:n_top]
    d = pred_delta[idx] - true_delta[idx]
    return float((d * d).mean())


def oof_auc(X, y, seeds=STRATA_SEEDS):
    """与 step6c 完全相同：全局已标准化的 X，5-fold（random_state=0），
    logistic random_state 取 3 个划分种子，返回各 pooled OOF AUROC。"""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    aucs = []
    for s in seeds:
        oof = cross_val_predict(
            LogisticRegression(C=1.0, max_iter=2000, random_state=s),
            X, y, cv=cv, method="predict_proba")[:, 1]
        aucs.append(float(roc_auc_score(y, oof)))
    return aucs


def tertile_labels(v):
    q = pd.qcut(pd.Series(v), 3, labels=["low", "mid", "high"])
    return q.astype(str).values


def strata_by_strength(X, strength, err_by_metric):
    """逐指标的 strength 三分层 OOF AUROC（step6c 协议：层内重拟合、3 种子均值）。

    返回 {metric: {overall: 3seed均值, strata: {low/mid/high: (auroc, n)}}}。
    """
    out = {}
    for metric, err in err_by_metric.items():
        y = (err > np.median(err)).astype(int)
        res = {"overall_auroc_3seeds": float(np.mean(oof_auc(X, y))),
               "n": int(len(y)), "strata": {}}
        tert = tertile_labels(strength)
        for lvl in ["low", "mid", "high"]:
            m = tert == lvl
            if m.sum() < 30 or len(np.unique(y[m])) < 2:
                res["strata"][lvl] = {"n": int(m.sum()), "auroc": None}
                continue
            res["strata"][lvl] = {
                "n": int(m.sum()),
                "auroc": float(np.mean(oof_auc(X[m], y[m]))),
                "worst_half_rate": float(y[m].mean())}
        out[metric] = res
    return out


def corr_entry(x, y):
    rho, p_rho = (float(v) for v in sps.spearmanr(x, y))
    r, p_r = (float(v) for v in sps.pearsonr(x, y))
    return {"spearman_rho": rho, "spearman_p": p_rho,
            "pearson_r": r, "pearson_p": p_r, "n": int(len(x))}


def run_dataset(data_name, out_dir, tag):
    out_path = os.path.join(REPO, out_dir)
    z = np.load(f"{out_path}/pred_expr_seed{SEED}.npz")
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()
    feat, feature_names, y_main, pred_df, _ = load_dataset(
        data_name, out_dir)
    assert feat.index.tolist() == perts
    X = StandardScaler().fit_transform(
        feat[feature_names].values.astype(np.float64))

    # ---------- 真值向量 + mean predictor（检验一） ----------
    true_delta, mean_delta, ctrl_aligned, n_train_conds = \
        build_truth_and_mean_baseline(data_name, perts, pred_genes)
    pred_delta = pred_expr.astype(np.float64) - ctrl_aligned

    err = {}
    err["top50_Pearson"] = 1.0 - pred_df["pearson_delta"].values
    err["all_Pearson"] = 1.0 - pred_df["pearson_delta_all"].values
    err["all_Spearman"] = 1.0 - np.array(
        [all_gene_spearman(pred_delta[i], true_delta[i])
         for i in range(len(perts))])
    err["top50_MSE"] = np.array(
        [top_deg_mse(pred_delta[i], true_delta[i])
         for i in range(len(perts))])
    err["mean_baseline_top50_Pearson"] = 1.0 - np.array(
        [top_deg_pearson(mean_delta, true_delta[i])
         for i in range(len(perts))])

    # 主口径自检：从向量重算 vs CSV（step6 标签来源）一致性
    recheck = 1.0 - np.array([top_deg_pearson(pred_delta[i], true_delta[i])
                              for i in range(len(perts))])
    csv_check = {
        "max_abs_diff_top50": float(
            np.max(np.abs(recheck - err["top50_Pearson"]))),
        "max_abs_diff_all": float(np.max(np.abs(
            (1.0 - np.array([np.corrcoef(pred_delta[i], true_delta[i])[0, 1]
                             for i in range(len(perts))]))
            - err["all_Pearson"])))}

    # ---------- 检验一：任务内在难度 ----------
    diff = corr_entry(err["mean_baseline_top50_Pearson"], err["top50_Pearson"])
    diff["verdict_rule"] = f"|spearman rho| > {DIFFICULTY_RHO_THRESH} " \
        "=> 结构部分来自任务内在难度"
    diff["verdict_intrinsic_difficulty"] = bool(
        abs(diff["spearman_rho"]) > DIFFICULTY_RHO_THRESH)

    # 相关矩阵一角：误差 vs step6 single_feature 特征（+ 特征间）
    probe_cols = feature_names  # 9 特征（与 gating_probe.json 同名单）
    corner = {}
    for c in probe_cols:
        corner[c] = corr_entry(feat[c].values, err["top50_Pearson"])
    feat_feat = {}
    for i, a in enumerate(probe_cols):
        for b in probe_cols[i + 1:]:
            feat_feat[f"{a}__{b}"] = float(
                sps.spearmanr(feat[a].values, feat[b].values)[0])

    # ---------- 检验二：指标 artifact ----------
    metric_corr = {}
    for m in ["all_Pearson", "all_Spearman", "top50_MSE"]:
        metric_corr[m] = corr_entry(err["top50_Pearson"], err[m])
    err_rank_corr_matrix = pd.DataFrame(
        {a: [float(sps.spearmanr(err[a], err[b])[0]) for b in ERROR_VARIANTS]
         for a in ERROR_VARIANTS}, index=ERROR_VARIANTS)

    strata = strata_by_strength(X, pred_df["strength"].values, err)

    res = {"task": "structure mechanism: intrinsic difficulty vs metric "
                   "artifact (step6f_structure_mechanism)",
           "dataset": tag, "n_heldout": int(len(perts)),
           "n_train_conditions_for_mean_baseline": int(n_train_conds),
           "label": "y = (1 - metric) > median, worst half; main metric = "
                    "top-50 |true_delta| Pearson (step1/step6)",
           "csv_recheck_vs_vector_recompute": csv_check,
           "test1_intrinsic_difficulty": {
               "mean_baseline": "per-gene mean true delta over all training "
                                "single perturbations (no model)",
               "error_metric": "1 - top50_DEG_Pearson, same as main",
               "mean_baseline_vs_gears_error": diff,
               "feature_error_corner": corner,
               "feature_feature_spearman": feat_feat},
           "test2_metric_artifact": {
               "error_rank_spearman_vs_main": metric_corr,
               "error_rank_corr_matrix": err_rank_corr_matrix.to_dict(),
               "strength_strata_auroc_by_metric": strata},
           "norman_lowmid_below_chance_check": {
               m: {lvl: strata[m]["strata"][lvl]["auroc"]
                   for lvl in ["low", "mid"]}
               for m in ERROR_VARIANTS}}
    with open(f"{out_path}/structure_mechanism.json", "w",
              encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)

    print(f"[{tag}] n={len(perts)} train_conds={n_train_conds}")
    print(f"[{tag}] csv recheck: top50 max|diff|={csv_check['max_abs_diff_top50']:.2e}"
          f" all max|diff|={csv_check['max_abs_diff_all']:.2e}")
    print(f"[{tag}] T1 difficulty: spearman rho={diff['spearman_rho']:.3f}"
          f" (p={diff['spearman_p']:.2e})  pearson r={diff['pearson_r']:.3f}"
          f" (p={diff['pearson_p']:.2e})  -> intrinsic_difficulty="
          f"{diff['verdict_intrinsic_difficulty']}")
    for m, c in metric_corr.items():
        print(f"[{tag}] T2 err-rank spearman(main vs {m}) rho={c['spearman_rho']:.3f}"
              f" (p={c['spearman_p']:.2e})")
    for m in ERROR_VARIANTS:
        s = strata[m]
        line = "  ".join(
            f"{lvl}: n={s['strata'][lvl]['n']}"
            + (f" AUROC={s['strata'][lvl]['auroc']:.3f}"
               if s["strata"][lvl]["auroc"] is not None else " AUROC=NA")
            for lvl in ["low", "mid", "high"])
        print(f"[{tag}] T2 strata[{m}] overall={s['overall_auroc_3seeds']:.3f} {line}")
    return res, err


def main():
    all_res, all_err = [], []
    for d in DATASETS:
        res, err = run_dataset(d["data"], d["out"], d["tag"])
        all_res.append(res)
        all_err.append(err)

    # ---- 图：左面板难度散点，右面板指标间误差排序相关热图（英文，300 dpi）----
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    for row, (res, err) in enumerate(zip(all_res, all_err)):
        short = "Norman" if row == 0 else "RPE1"
        ax = axes[row, 0]
        x = err["mean_baseline_top50_Pearson"]
        yv = err["top50_Pearson"]
        ax.scatter(x, yv, s=14, alpha=0.55, edgecolor="none")
        d1 = res["test1_intrinsic_difficulty"]["mean_baseline_vs_gears_error"]
        ax.set_xlabel(f"mean-baseline error (1 - top-{TOP_N_DEG} Pearson)")
        ax.set_ylabel("GEARS error (1 - top-50 Pearson)")
        ax.set_title(f"{short}: rho={d1['spearman_rho']:.2f} "
                     f"(p={d1['spearman_p']:.1e}), r={d1['pearson_r']:.2f}",
                     fontsize=11)
        ax.grid(alpha=0.25)

        ax = axes[row, 1]
        M = pd.DataFrame(res["test2_metric_artifact"]
                         ["error_rank_corr_matrix"]).loc[
            ERROR_VARIANTS, ERROR_VARIANTS]
        im = ax.imshow(M.values, vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(M)))
        ax.set_xticklabels(["top50\nPearson", "all\nPearson", "all\nSpearman",
                            "top50\nMSE", "mean-b.\ntop50"], fontsize=8)
        ax.set_yticks(range(len(M)))
        ax.set_yticklabels(["top50 Pearson", "all Pearson", "all Spearman",
                            "top50 MSE", "mean-baseline top50"], fontsize=8)
        for i in range(len(M)):
            for j in range(len(M)):
                ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center",
                        va="center", fontsize=8,
                        color="white" if M.values[i, j] < 0.6 else "black")
        ax.set_title(f"{short}: error-ranking Spearman", fontsize=11)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Structure mechanism: intrinsic difficulty (left) vs "
                 "metric sensitivity (right)")
    fig.tight_layout()
    fig.savefig(f"{REPO}/outputs/structure_mechanism.png", dpi=300)
    print("[SAVED]", f"{REPO}/outputs/structure_mechanism.png")


if __name__ == "__main__":
    main()
