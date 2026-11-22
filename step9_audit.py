"""Step9 运行器：门控信号的机制审计（2026-10-14，审稿加固，不改判据只加检验）。

背景：step6 门控探针的全特征 AUROC 0.66-0.80，但最好特征是 delta_max
（预测谱幅度）。存在一种机械解释：预测 delta 越小，top-50 基因上的 Pearson
相关系数数值上越不稳定 → 误差天然变大。若信号全部来自这个"效应量几何"，
那探针检测的是"评测病态性"而非"模型自知"——结论的表述必须不同。

本脚本在**完全不触碰预注册判据**的前提下回答四个问题：
1. 幅度家族消融：剔除 delta_l2/l1/max/top50_mean 后，非幅度残余特征
   （n_de_025/frac_de_025/neg_frac/delta_entropy + log_n_cells）还剩多少信号？
   （Q1 的答案决定论文的真实身份）
2. 置换 null：打乱标签的 AUROC 分布（检验 0.66 是否显著异于随机）。
3. Bootstrap CI：按扰动 cluster 重采样的 CV AUROC 置信区间（替代 fold std）。
4. 覆盖-风险曲线：过滤掉探针认为最可疑的 k% 查询，能捕获多少比例的误差、
   保留预测的均误差降到多少（把 AUROC 翻译成部署语言）。
5. 强度分解（诊断用，非部署特征）：strength 是真值派生（step1_run.py:96-105），
   只用于回答"门控信号有多少是效应量代理"。

用法：
    source venv/Scripts/activate
    python step9_audit.py                          # 默认 outputs (norman)
    DATA=replogle_rpe1_essential OUT_DIR=outputs_rpe1 python step9_audit.py
输出（$OUT_DIR/）：audit.json, audit_coverage_risk.png, audit_verdict.md
"""
import json
import os
import time

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import compute_features

OUT = os.environ.get("OUT_DIR", "outputs")
DATA = os.environ.get("DATA", "replogle_k562_essential")
SEED = 0
N_PERM = 1000        # 置换 null 次数
N_BOOT = 1000        # bootstrap 重采样次数
os.makedirs(OUT, exist_ok=True)

# 幅度家族 vs 非幅度残余（消融的两臂）
MAGNITUDE_FAMILY = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean"]
RESIDUAL_FEATURES = ["n_de_025", "frac_de_025", "neg_frac", "delta_entropy"]

PROBE = lambda: make_pipeline(  # noqa: E731  与 step6 完全相同的探针设置
    StandardScaler(),
    LogisticRegression(C=1.0, max_iter=2000, random_state=SEED))


def cv_auroc(X, y):
    """与 step6 run_probe 相同的分层 5-fold CV（shuffle, seed=0）。
    返回 (fold 均值 AUROC, pooled out-of-fold 决策分)。"""
    from sklearn.model_selection import cross_validate
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = cross_val_predict(PROBE(), X, y, cv=cv, method="decision_function")
    sc = cross_validate(PROBE(), X, y, cv=cv, scoring="roc_auc")["test_score"]
    return float(sc.mean()), oof


def load_all():
    """镜像 step6_run.py 第 54-103 行的数据口径。"""
    z = np.load(f"{OUT}/pred_expr_seed{SEED}.npz")
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()
    adata = ad.read_h5ad(
        f"data/gears_demo/{DATA}/perturb_processed.h5ad", backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    adata.file.close()
    del adata
    ctrl_aligned = ctrl_mean if var_genes == pred_genes else np.array(
        dict(zip(var_genes, ctrl_mean))[g] for g in pred_genes)

    pred_df = pd.read_csv(
        f"{OUT}/predictions_seed{SEED}.csv").set_index("perturbation").loc[perts]
    feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    error = 1.0 - pred_df["pearson_delta"].values
    y = (error > np.median(error)).astype(int)
    return feat, pred_df, error, y


def main():
    t0 = time.time()
    feat, pred_df, error, y = load_all()
    n = len(y)
    print(f"[step9] n={n}, 正类={y.sum()}, OUT={OUT}, DATA={DATA}")

    X_full = feat[MAGNITUDE_FAMILY + RESIDUAL_FEATURES + ["log_n_cells"]].values
    X_mag = feat[MAGNITUDE_FAMILY + ["log_n_cells"]].values
    X_res = feat[RESIDUAL_FEATURES + ["log_n_cells"]].values

    # ---- 1. 幅度家族消融 ----
    auroc_full, oof_full = cv_auroc(X_full, y)
    auroc_mag, _ = cv_auroc(X_mag, y)
    auroc_res, _ = cv_auroc(X_res, y)
    pooled_full = roc_auc_score(y, oof_full)
    print(f"[step9] 消融: full={auroc_full:.4f}  幅度家族={auroc_mag:.4f}  "
          f"非幅度残余={auroc_res:.4f}")

    # ---- 2. 置换 null（同 CV 流程，打乱标签）----
    rng = np.random.default_rng(0)
    perm_aurocs = np.empty(N_PERM)
    for i in range(N_PERM):
        yp = rng.permutation(y)
        perm_aurocs[i] = roc_auc_score(yp, cross_val_predict(
            PROBE(), X_full, yp, cv=StratifiedKFold(
                n_splits=5, shuffle=True, random_state=SEED),
            method="decision_function"))
        if (i + 1) % 200 == 0:
            print(f"[step9] perm {i+1}/{N_PERM} "
                  f"({time.time()-t0:.0f}s)")
    p_value = float((perm_aurocs >= pooled_full).mean())
    print(f"[step9] 置换 null: mean={perm_aurocs.mean():.4f} "
          f"std={perm_aurocs.std():.4f}  p(>=observed)={p_value:.4f}")

    # ---- 3. Bootstrap CI（按扰动 cluster 重采样 pooled OOF 对）----
    boot_aurocs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            boot_aurocs[i] = np.nan
            continue
        boot_aurocs[i] = roc_auc_score(y[idx], oof_full[idx])
    boot_clean = boot_aurocs[~np.isnan(boot_aurocs)]
    ci_lo, ci_hi = np.percentile(boot_clean, [2.5, 97.5])
    print(f"[step9] bootstrap CI (pooled OOF): "
          f"[{ci_lo:.4f}, {ci_hi:.4f}] (n={len(boot_clean)})")

    # ---- 4. 覆盖-风险曲线（pooled OOF 分数）----
    order = np.argsort(-oof_full)          # 分数高 = 探针认为更可能错
    err_sorted = error[order]
    cum_err = np.cumsum(err_sorted)
    total_err = err_sorted.sum()
    frac = np.arange(1, n + 1) / n
    captured = cum_err / total_err
    cov = {"frac_filtered": frac.tolist(), "err_captured": captured.tolist()}
    grid = [0.10, 0.20, 0.30]
    cov_points = {}
    for g_ in grid:
        k = int(g_ * n)
        retained = n - k
        cov_points[f"filter_{int(g_*100)}pct"] = {
            "err_captured_frac": float(captured[k - 1]),
            "retained_mean_error": float(err_sorted[k:].mean()),
            "overall_mean_error": float(error.mean()),
        }

    # ---- 5. 强度分解（诊断，strength 为真值派生，非部署特征）----
    strength = pred_df["strength"].values
    delta_max = feat["delta_max"].values
    decomp = {
        "spearman_error_strength": float(spearmanr(error, strength).statistic),
        "spearman_error_delta_max": float(spearmanr(error, delta_max).statistic),
        "spearman_strength_delta_max":
            float(spearmanr(strength, delta_max).statistic),
        "auroc_strength_groundtruth":
            float(roc_auc_score(y, -strength)),
        "auroc_delta_max_single":
            float(roc_auc_score(y, -feat["delta_max"].values)),
    }
    print(f"[step9] 分解: rho(err,strength)={decomp['spearman_error_strength']:.3f} "
          f"rho(err,delta_max)={decomp['spearman_error_delta_max']:.3f} "
          f"AUROC(strength)={decomp['auroc_strength_groundtruth']:.4f}")

    # ---- 图：覆盖-风险曲线 ----
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(frac, captured, label="gating probe (OOF)")
    ax.plot(frac, frac, "k--", alpha=0.5, label="random filter")
    for g_ in grid:
        k = int(g_ * n)
        ax.plot(g_, captured[k - 1], "o")
        ax.annotate(f"{int(g_*100)}%: {captured[k-1]:.0%} errors",
                    (g_, captured[k - 1]), textcoords="offset points",
                    xytext=(6, -4), fontsize=8)
    ax.set_xlabel("fraction of queries filtered")
    ax.set_ylabel("fraction of total error captured")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{OUT}/audit_coverage_risk.png", dpi=300)

    # ---- 产物 ----
    result = {
        "task": "mechanism audit of the gating signal (step9)",
        "date": "2026-10-14",
        "out_dir": OUT, "data": DATA, "seed": SEED,
        "n_perts": int(n), "n_positive": int(y.sum()),
        "ablation": {
            "full_9feat_cv_auroc": float(auroc_full),
            "magnitude_family_cv_auroc": float(auroc_mag),
            "residual_nonmagnitude_cv_auroc": float(auroc_res),
            "magnitude_family": MAGNITUDE_FAMILY,
            "residual_features": RESIDUAL_FEATURES,
            "interpretation_rule": (
                "residual >= 0.6 -> 幅度之外仍有模型侧信号; "
                "residual < 0.6 -> 信号基本=效应量几何(评测病态性), "
                "论文结论改写为部署可检测的难度代理"),
        },
        "permutation_null": {
            "n_perm": N_PERM, "mean": float(perm_aurocs.mean()),
            "std": float(perm_aurocs.std()),
            "p_value_ge_observed": p_value,
        },
        "bootstrap_ci_pooled_oof": {
            "n_boot": N_BOOT, "ci95": [float(ci_lo), float(ci_hi)],
        },
        "coverage_risk": {
            "points": cov_points,
            "note": "filter = discard highest-OOF-score queries first",
        },
        "strength_decomposition_groundtruth_diagnostic": decomp,
        "runtime_sec": time.time() - t0,
    }
    with open(f"{OUT}/audit.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    verdict = f"""# audit_verdict.md — 门控信号机制审计（Step9，2026-10-14）

## 四个问题的答案（{DATA}，n={n}，seed 0）

1. **幅度家族消融**：全特征 CV AUROC = {auroc_full:.4f}；只留幅度家族 =
   {auroc_mag:.4f}；**非幅度残余 = {auroc_res:.4f}**。
2. **置换 null**（{N_PERM} 次）：null 均值 {perm_aurocs.mean():.4f} ±
   {perm_aurocs.std():.4f}，观察到 {pooled_full:.4f} 的 p = {p_value:.4f}。
3. **Bootstrap 95% CI**（pooled OOF，{N_BOOT} 次重采样）：
   [{ci_lo:.4f}, {ci_hi:.4f}]。
4. **覆盖-风险**：过滤最可疑 {int(grid[1]*100)}% 查询 → 捕获
   {captured[int(grid[1]*n)-1]:.1%} 的总误差，保留预测均误差
   {err_sorted[int(grid[1]*n):].mean():.4f}（全体 {error.mean():.4f}）。
5. **强度分解**（真值诊断）：rho(error, strength) =
   {decomp['spearman_error_strength']:.3f}；rho(error, delta_max) =
   {decomp['spearman_error_delta_max']:.3f}；strength 单独 AUROC =
   {decomp['auroc_strength_groundtruth']:.4f}（不可部署，仅诊断）。

## 判读（与预注册判据无关，不改写 step6 归档结论）

见 audit.json 的 interpretation_rule。产物：audit.json、
audit_coverage_risk.png、本文件。
"""
    with open(f"{OUT}/audit_verdict.md", "w", encoding="utf-8") as f:
        f.write(verdict)
    print(f"[step9] done in {time.time()-t0:.0f}s -> {OUT}/audit.json")


if __name__ == "__main__":
    main()
