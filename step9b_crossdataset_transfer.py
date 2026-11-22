# step9b_crossdataset_transfer.py — 同一 K562 就地探针跨数据集迁移到 RPE1
#
# 背景（论文核心发现的缺口实验）：
#   "误差签名是模型×数据集涌现性质"目前跨数据集证据只有签名侧
#   （scGenePT 同模型换数据集签名余弦 0.01 / GEARS 0.044，step6f、step9），
#   检测端只有各自就地 AUROC（K562 0.66 / RPE1 0.80），还缺
#   "同一探针跨数据集"的迁移数字。本脚本把归档的 K562 就地探针
#   原样（不重训、不 refit scaler）打分到 RPE1 全部 306 扰动上：
#     - 迁移 AUROC ≥ 0.6 且 CI 下界 > 0.5 → "检测可迁移、签名方向不可迁移"
#     - 0.5–0.6 → 弱迁移，如实报告
#     - < 0.5 → 跨数据集全面不可迁移
#
# 数据与口径（全部归档，无新下载）：
#   - K562 拟合探针：outputs/gating_probe_transfer.json（deployment 版，
#     in-sample 拟合于全部 217 留出扰动；feature_names / scaler_mean /
#     scaler_scale / coef / intercept），直接加载，不重新拟合。
#   - RPE1 特征与标签：step6d_naive_baselines.load_dataset(
#     "replogle_rpe1_essential", "outputs_rpe1")（step6d/6f 同款口径，
#     y = (1 - pearson_delta) > RPE1 全集中位数，与就地探针完全同标签口径）。
#   - 9 特征按探针的 feature_names 重排，保证列序一致。
#   - RPE1 就地参照：outputs_rpe1/gating_probe.json 的
#     cv_auroc_mean = 0.7965（引用归档值，不复算）；
#     K562 就地参照：outputs/gating_probe.json 的 cv_auroc_mean = 0.6556。
#
# 统计：AUROC 的扰动级 bootstrap 95% CI（B=1000，主 seed 20261022
#   派生独立 Generator，与 step8b 同体系）；迁移 vs 随机的单侧置换 p
#   （打乱标签 B=1000 次）。
#
# 对照组：K562 探针中 delta_max 系数符号直接用于 RPE1 原始 delta_max
#   排序的单特征迁移 AUROC（与 step6d delta_max 单特征基线同思路）。
#
# 产物：outputs_crossmodel/analysis/crossdataset_transfer.{json,png}
#       （300 dpi，英文，三柱 + CI：K562 就地 / K562→RPE1 迁移 / RPE1 就地）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step9b_crossdataset_transfer.py
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")

BOOT_MASTER_SEED = 20261022  # bootstrap / 置换检验专属主 seed（step8b 同体系）
B_BOOT = 1000
B_PERM = 1000

# 预注册判据（写死，不事后调整）
PREREG_CRITERIA = {
    "transfer_ok": "transfer AUROC >= 0.6 AND bootstrap 95% CI lower bound > 0.5",
    "weak": "0.5 <= transfer AUROC < 0.6",
    "no_transfer": "transfer AUROC < 0.5",
    "verdict_ok": ("检测跨数据集可迁移，签名方向不可迁移"
                   "（detection transfers across datasets, "
                   "signature direction does not）"),
    "verdict_weak": "弱迁移（weak transfer），如实报告",
    "verdict_no": "跨数据集全面不可迁移（nothing transfers across datasets）",
}


def load_rpe1_xy():
    """RPE1 (X, y, feature_names, perts)，口径与 step6d/6f 完全一致。"""
    from step6d_naive_baselines import load_dataset
    feat, feature_names, y, pred_df, _ = load_dataset(
        "replogle_rpe1_essential", "outputs_rpe1")
    X = feat[feature_names].values.astype(np.float64)
    return X, y, feature_names, list(pred_df.index)


def transfer_scores(rec, X):
    """归档探针打分：(X - mean)/scale @ coef + intercept（不 refit）。"""
    Xs = (X - np.array(rec["scaler_mean"])) / np.array(rec["scaler_scale"])
    return Xs @ np.array(rec["coef"]) + rec["intercept"]


def bootstrap_auroc_ci(scores, y, rng, b_boot=B_BOOT):
    """扰动级 bootstrap：有放回重采样 → AUROC 的 2.5/97.5 分位数 CI。"""
    n = len(y)
    aucs = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:  # 极端重采样单类，重抽
            idx = rng.integers(0, n, n)
        aucs[b] = roc_auc_score(y[idx], scores[idx])
    return (float(np.percentile(aucs, 2.5)),
            float(np.percentile(aucs, 97.5)))


def permutation_p_value(scores, y, observed, rng, b_perm=B_PERM):
    """单侧置换检验：打乱标签 B 次，p = (1 + #{perm >= obs})/(B+1)。"""
    n = len(y)
    cnt = 0
    for _ in range(b_perm):
        perm = rng.permutation(n)
        if roc_auc_score(y, scores[perm]) >= observed:
            cnt += 1
    return (cnt + 1) / (b_perm + 1)


def main():
    # --- 1. 加载 K562 就地探针（deployment 版，不重训） ---
    with open(os.path.join(REPO, "outputs", "gating_probe_transfer.json"),
              encoding="utf-8") as f:
        probe = json.load(f)
    assert probe["dataset"] == "replogle_k562_essential"

    # --- 2. RPE1 特征与标签（step6d 口径，与就地探针同标签定义） ---
    X, y, feature_names, perts = load_rpe1_xy()
    assert list(probe["feature_names"]) == list(feature_names), \
        f"特征列不一致: {probe['feature_names']} vs {feature_names}"
    n = len(y)
    assert n == 306, f"RPE1 扰动数应为 306，实际 {n}"
    print(f"[RPE1] n={n}，9 特征列序与 K562 探针一致 OK")

    # --- 3. 迁移打分（归档 scaler，不 refit）→ AUROC ---
    s_transfer = transfer_scores(probe, X)
    auroc_transfer = float(roc_auc_score(y, s_transfer))

    # --- 4. 对照组：delta_max 系数符号迁移（单特征，无训练） ---
    j_max = list(probe["feature_names"]).index("delta_max")
    sign = float(np.sign(probe["coef"][j_max]))
    s_dm = sign * X[:, j_max]  # 标准化保序，符号直接用原始 delta_max
    auroc_dm = float(roc_auc_score(y, s_dm))

    # --- 5. 就地参照（引用归档，不复算） ---
    with open(os.path.join(REPO, "outputs_rpe1", "gating_probe.json"),
              encoding="utf-8") as f:
        rpe1_insitu = json.load(f)
    with open(os.path.join(REPO, "outputs", "gating_probe.json"),
              encoding="utf-8") as f:
        k562_insitu = json.load(f)
    auroc_k562_insitu = float(k562_insitu["cv_auroc_mean"])
    auroc_rpe1_insitu = float(rpe1_insitu["cv_auroc_mean"])
    assert abs(auroc_rpe1_insitu - 0.7965) < 1e-3, \
        f"RPE1 就地归档 AUROC 与任务口径 0.7965 不符: {auroc_rpe1_insitu}"

    # --- 6. 统计：bootstrap CI + 单侧置换 p（seed 20261022 体系） ---
    rng_master = np.random.default_rng(BOOT_MASTER_SEED)
    rng_boot = np.random.default_rng(rng_master.integers(0, 2**31))
    lo, hi = bootstrap_auroc_ci(s_transfer, y, rng_boot)
    rng_perm = np.random.default_rng(rng_master.integers(0, 2**31))
    p_value = permutation_p_value(s_transfer, y, auroc_transfer, rng_perm)

    # --- 7. 预注册判据落档 ---
    if auroc_transfer >= 0.6 and lo > 0.5:
        verdict = PREREG_CRITERIA["verdict_ok"]
        bucket = "transfer_ok"
    elif auroc_transfer >= 0.5:
        verdict = PREREG_CRITERIA["verdict_weak"]
        bucket = "weak"
    else:
        verdict = PREREG_CRITERIA["verdict_no"]
        bucket = "no_transfer"

    res = {"task": "cross-dataset transfer: same K562 in-situ probe applied "
                   "verbatim to RPE1 (step9b)",
           "probe": {"source": "outputs/gating_probe_transfer.json",
                     "dataset": probe["dataset"], "n_fit": probe["n"],
                     "note": "deployment-style in-sample fit on all 217 "
                             "held-out K562 perturbations; NOT refitted "
                             "here; archived scaler applied without refit"},
           "target": {"dataset": "replogle_rpe1_essential", "n": n,
                      "label": "y = (1 - pearson_delta) > median over all "
                               "306 RPE1 perturbations (same definition as "
                               "the in-situ RPE1 probe, via "
                               "step6d_naive_baselines.load_dataset)"},
           "transfer": {"auroc": auroc_transfer,
                        "boot_ci95": [lo, hi], "boot_b": B_BOOT,
                        "permutation_p_one_sided_vs_random": p_value,
                        "perm_b": B_PERM,
                        "seed": BOOT_MASTER_SEED},
           "delta_max_sign_transfer": {"auroc": auroc_dm,
                                       "rule": "sign(K562 coef[delta_max]) "
                                               "* raw RPE1 delta_max"},
           "references_archived": {
               "k562_in_situ_cv_auroc_mean": auroc_k562_insitu,
               "rpe1_in_situ_cv_auroc_mean": auroc_rpe1_insitu,
               "note": "cited from outputs{,_rpe1}/gating_probe.json, "
                       "not recomputed"},
           "preregistered_criteria": PREREG_CRITERIA,
           "verdict_bucket": bucket, "verdict": verdict}

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, "crossdataset_transfer.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[RESULT] K562->RPE1 transfer AUROC = {auroc_transfer:.4f} "
          f"95% CI [{lo:.4f}, {hi:.4f}], p(one-sided) = {p_value:.4f}")
    print(f"[RESULT] delta_max sign transfer AUROC = {auroc_dm:.4f}")
    print(f"[VERDICT] {bucket}: {verdict}")
    print("[SAVED]", out_json)
    make_figure(auroc_transfer, lo, hi, auroc_k562_insitu, auroc_rpe1_insitu,
                auroc_dm, n)


def make_figure(auroc_transfer, lo, hi, auroc_k562, auroc_rpe1, auroc_dm, n):
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    labels = ["GEARS probe\nK562 in-situ",
              "K562 probe →\nRPE1 (transfer)",
              "GEARS probe\nRPE1 in-situ"]
    aucs = [auroc_k562, auroc_transfer, auroc_rpe1]
    colors = ["#8c8c8c", "#1f77b4", "#2ca02c"]
    bars = ax.bar([0, 1, 2], aucs, 0.55, color=colors, edgecolor="gray",
                  lw=0.5)
    # 迁移柱加 bootstrap CI 须线；就地柱引用归档 CV 标准误（不重复画图误差线，
    # 只标注 CI 属于迁移探针）
    ax.errorbar([1], [auroc_transfer],
                yerr=[[auroc_transfer - lo], [hi - auroc_transfer]],
                fmt="none", ecolor="black", elinewidth=1.3, capsize=4)
    for b, a in zip(bars, aucs):
        ax.text(b.get_x() + b.get_width() / 2, a + 0.012, f"{a:.3f}",
                ha="center", fontsize=9)
    ax.text(1, 0.32, f"95% CI\n[{lo:.3f},\n{hi:.3f}]", ha="center",
            fontsize=7.5, color="white")
    ax.axhline(0.5, color="red", ls=":", lw=1.2)
    ax.text(2.35, 0.505, "chance", color="red", fontsize=7, ha="right")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("AUROC (median-split error labels)", fontsize=9)
    ax.set_title("Cross-dataset transfer of one error probe "
                 f"(RPE1 n={n})", fontsize=10)
    ax.text(0.01, 0.97, f"delta_max sign-only transfer: {auroc_dm:.3f}",
            transform=ax.transAxes, fontsize=7.5, color="gray", va="top")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path = os.path.join(ANALYSIS_DIR, "crossdataset_transfer.png")
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
