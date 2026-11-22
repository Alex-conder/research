# step8c_sample_size_curve.py — 就地校准的最小样本量曲线
#
# 回答的问题："就地校准（对目标模型就地拟合门控探针）需要多少个扰动
# 才够？"
#
# 设计（先写死再跑）：
#   1. 数据：与 step8b 完全相同的 4 个（模型×数据集）签名数据集
#      （加载函数直接 import step8b.load_xy，口径不复制、不发明）。
#   2. 标签：y = (1 − pearson_delta) > median(error)，在**全量数据**上一次性
#      定义（step6 定义性标签），划分后不再重算，保证各 n 的测试标签固定。
#   3. 划分：与 step7/step7b 完全相同的扰动级 70/30 生成逻辑
#      （master seed 20261016，split_seeds=[0,1,2]，逐 seed 由主 seed 派生
#      Generator、permutation、前 30% 为测试集）——同一套划分，结果可比。
#   4. 扫描：训练子集大小 n ∈ {10,15,20,30,45,60,90,120,152}， capped 在
#      各数据集 70% 训练池大小（K562 217→152，RPE1 306→214，scF 45→31）。
#      每个 n 从训练池随机抽 20 个子集（不放回、各 n 独立抽样），子集内
#      in-sample 拟合同规格探针（StandardScaler + LogisticRegression
#      C=1.0, max_iter=2000），在固定 30% 测试集上算 AUROC。
#      聚合：3 个划分种子 × 20 子集 = 60 个 AUROC → 均值 + 95% percentile CI。
#   5. 判据线：AUROC = 0.6（step6 定义性阈值），图上标注。
#   6. 读数：各模型 mean AUROC 首次 ≥ 0.6 的最小 n（如实记录，含"全程
#      未达 0.6"的情形）。
#
# 随机源：SUBSET_MASTER_SEED = 20261021（子集抽样专属 seed；划分种子仍用
# step7 体系的 20261016，保证与 step7/7b 同一测试集）。
#
# 产物：outputs_crossmodel/analysis/sample_size_curve.json + .png
#       （300 dpi，英文，x=n log 刻度，y=AUROC，四模型均值±CI 曲线）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step8c_sample_size_curve.py
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from step8b_signature_ci import SIG_LABELS, SIG_ORDER, load_xy

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")

SPLIT_MASTER_SEED = 20261016   # 与 step7/step7b 相同 → 相同 70/30 划分
SPLIT_SEEDS = [0, 1, 2]
TEST_FRAC = 0.30
SUBSET_MASTER_SEED = 20261021  # 子集抽样专属 seed
N_SUBSETS = 20
AUROC_CRITERION = 0.6

N_GRID = [10, 15, 20, 30, 45, 60, 90, 120, 152]


def step7_splits(n):
    """复刻 step7b 的 70/30 划分生成（同一 master seed → 同一 permutation）。"""
    rng_master = np.random.default_rng(SPLIT_MASTER_SEED)
    out = []
    for _sp in SPLIT_SEEDS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        perm = rng.permutation(n)
        n_te = int(round(TEST_FRAC * n))
        out.append((perm[n_te:], perm[:n_te]))  # (train, test)
    return out


def fit_score(X_tr, y_tr, X_te):
    """子集内 in-sample 拟合同规格探针，返回测试集 decision score。"""
    sc = StandardScaler().fit(X_tr)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(
        sc.transform(X_tr), y_tr)
    return sc.transform(X_te) @ clf.coef_[0] + clf.intercept_[0]


def main():
    rng_master = np.random.default_rng(SUBSET_MASTER_SEED)
    res = {"task": "minimum sample size for in-situ calibration "
                   "(step8c)",
           "label": "y = (1 - pearson_delta) > median over FULL dataset "
                    "(step6 definition, fixed before splitting)",
           "split": "perturbation-level 70/30, split_seeds=[0,1,2], "
                    f"identical generation as step7/step7b (master seed "
                    f"{SPLIT_MASTER_SEED})",
           "sampling": f"{N_SUBSETS} random subsets per (n, split seed), "
                       "drawn without replacement from the 70% train pool; "
                       f"dedicated master seed {SUBSET_MASTER_SEED}",
           "probe": "StandardScaler + LogisticRegression(C=1.0, "
                    "max_iter=2000), in-sample fit on each subset "
                    "(deployment-style)",
           "auroc_criterion": AUROC_CRITERION,
           "curves": {}}

    fig_x, fig_y, fig_lo, fig_hi = {}, {}, {}, {}
    for sig_name in SIG_ORDER:
        X, y, _fn = load_xy(sig_name)
        n = len(y)
        splits = step7_splits(n)
        train_pool = min(len(tr) for tr, _te in splits)
        n_grid = [k for k in N_GRID if k <= train_pool]
        curve = []
        auroc_mean, auroc_lo, auroc_hi = [], [], []
        for k in n_grid:
            vals = []
            for tr, te in splits:
                # 子集抽样专用 Generator：按 (k, split) 从主 seed 派生，
                # 20 次重采样共用一条流（顺序确定，可复现）
                rng = np.random.default_rng(
                    rng_master.integers(0, 2**31))
                for _rep in range(N_SUBSETS):
                    sub = rng.choice(tr, size=k, replace=False)
                    while len(np.unique(y[sub])) < 2:
                        sub = rng.choice(tr, size=k, replace=False)
                    s = fit_score(X[sub], y[sub], X[te])
                    vals.append(float(roc_auc_score(y[te], s)))
            v = np.array(vals)
            m, lo, hi = float(v.mean()), float(np.percentile(v, 2.5)), \
                float(np.percentile(v, 97.5))
            curve.append({"n": k, "auroc_mean": m,
                          "ci_lo": lo, "ci_hi": hi,
                          "n_evals": int(len(v))})
            auroc_mean.append(m)
            auroc_lo.append(lo)
            auroc_hi.append(hi)
        reached = next((c["n"] for c in curve
                        if c["auroc_mean"] >= AUROC_CRITERION), None)
        res["curves"][sig_name] = {
            "label": SIG_LABELS[sig_name], "n_available": int(n),
            "train_pool": int(train_pool), "n_grid": n_grid,
            "curve": curve,
            "min_n_mean_auroc_ge_0.6": reached}
        fig_x[sig_name], fig_y[sig_name] = n_grid, auroc_mean
        fig_lo[sig_name], fig_hi[sig_name] = auroc_lo, auroc_hi
        print(f"[{sig_name}] n={n} train_pool={train_pool} "
              f"grid={n_grid} | 最小达标 n={reached} | "
              f"末端 AUROC={auroc_mean[-1]:.3f}")

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, "sample_size_curve.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, fig_x, fig_y, fig_lo, fig_hi,
                os.path.join(ANALYSIS_DIR, "sample_size_curve.png"))


def make_figure(res, xs, ys, los, his, path):
    colors = {"GEARS(K562)": "#1f77b4", "scgenept_k562": "#d62728",
              "scgenept_rpe1": "#ff7f0e", "scf_k562": "#2ca02c",
              "scgpt_k562": "#9467bd"}
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    for sig_name in SIG_ORDER:
        c = colors[sig_name]
        x, yv = np.array(xs[sig_name]), np.array(ys[sig_name])
        lo, hi = np.array(los[sig_name]), np.array(his[sig_name])
        ax.plot(x, yv, "-o", color=c, ms=4, lw=1.6,
                label=res["curves"][sig_name]["label"])
        ax.fill_between(x, lo, hi, color=c, alpha=0.15, lw=0)
        reached = res["curves"][sig_name]["min_n_mean_auroc_ge_0.6"]
        if reached is not None:
            j = xs[sig_name].index(reached)
            ax.plot(reached, yv[j], "s", color=c, ms=8, mfc="none",
                    mec=c, mew=1.5)
    ax.axhline(AUROC_CRITERION, color="red", ls="--", lw=1.2)
    ax.text(0.995, AUROC_CRITERION + 0.008, "AUROC = 0.6 criterion",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8,
            color="red")
    ax.axhline(0.5, color="gray", ls=":", lw=1.0)
    ax.text(0.995, 0.505, "chance = 0.5",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=7,
            color="gray")
    ax.set_xscale("log")
    ax.set_xticks([10, 15, 20, 30, 45, 60, 90, 120, 152])
    ax.set_xticklabels([10, 15, 20, 30, 45, 60, 90, 120, 152], fontsize=8)
    ax.set_xlim(9, 170)
    ax.set_ylim(0.35, 0.85)
    ax.set_xlabel("number of perturbations in calibration subset "
                  "(train side only, log scale)")
    ax.set_ylabel("test-set AUROC (fixed 30% holdout)")
    ax.set_title("In-situ calibration sample-size curve", fontsize=11)
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, loc="lower right")
    fig.text(0.5, 0.005, "square marker = smallest n with mean AUROC "
                          ">= 0.6", ha="center", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
