# step7_routing.py — 跨模型逐扰动路由（routing）：误差探针能否逐扰动预测
# "GEARS 和 scGenePT 谁的预测更可信"，并按此路由提升整体预测质量？
#
# 预注册式分析设计（先写死再跑）：
#   1. 扰动级 70/30 划分，3 个划分种子（SPLIT_SEEDS = [0, 1, 2]），
#      30% 作 held-out 路由测试集；
#   2. 在 70% 上：分别对 GEARS 特征、scGenePT 特征就地训同规格探针
#      （StandardScaler + L2 logistic C=1.0, max_iter=2000；
#      标签 y = (1 - pearson_delta) > 训练集中位数，与 step6 同口径）；
#   3. 在 30% 上：逐扰动比较两探针的预测误差概率，路由 = 选预测误差
#      较低的模型（平手取 GEARS）。对照组：
#        a) all-GEARS   b) all-scGenePT   c) oracle（逐扰动选真实 pearson 更高者，上界）
#   4. 指标：四组的 mean top-50 pearson_delta（30% 测试扰动），
#      3 seeds 报均值±std；路由组相对最优单模型的提升；
#   5. 显著性：配对 bootstrap（扰动为抽样单位，1000 次/种子），
#      路由 vs 最优单模型（逐种子取两单模型中测试均值更高者）的单侧 p。
#
# 数据（全部本地归档，不重跑服务器）：
#   GEARS 预测：outputs/pred_expr_seed0.npz（217 扰动 × 5000 基因，归档恢复；
#     与 step6 探针同一 artifact，SEED=0 与 step6d.load_dataset 一致）；
#   GEARS 标签：outputs/predictions_seed0.csv（pearson_delta = top-50 |delta| Pearson，
#     step1_run.py:160-162 口径）；
#   scGenePT 预测/标签：outputs_crossmodel/scgenept/scgenept_k562_predictions.npz
#     + scgenept_k562_labels.csv（经 step6f.load_crossmodel 对齐）；
#   特征：step6_gating_probe.compute_features（9 特征 = 8 谱统计 + log_n_cells，
#     纯 pred_expr/ctrl_mean/n_cells 函数，泄漏纪律同 step6）。
#
# 产物：outputs_crossmodel/analysis/routing.json + routing.png（300dpi 英文）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7_routing.py --dry-run
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7_routing.py
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features
from step6d_naive_baselines import load_dataset as load_gears_dataset
from step6f_crossmodel import load_crossmodel

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")

SPLIT_SEEDS = [0, 1, 2]
TEST_FRAC = 0.30
N_BOOT = 1000

GEARS_NPZ = os.path.join(REPO, "outputs", "pred_expr_seed0.npz")
GEARS_CSV = os.path.join(REPO, "outputs", "predictions_seed0.csv")
SG_NPZ = os.path.join(REPO, "outputs_crossmodel", "scgenept",
                      "scgenept_k562_predictions.npz")
SG_CSV = os.path.join(REPO, "outputs_crossmodel", "scgenept",
                      "scgenept_k562_labels.csv")
DATA_NAME = "replogle_k562_essential"


def check_products():
    missing = [p for p in (GEARS_NPZ, GEARS_CSV, SG_NPZ, SG_CSV)
               if not os.path.exists(p)]
    if missing:
        for p in missing:
            print(f"[MISSING] {p}")
        raise SystemExit(1)
    print(f"[OK] {os.path.relpath(GEARS_NPZ, REPO)} + predictions_seed0.csv（GEARS，归档恢复）")
    print(f"[OK] {os.path.relpath(SG_NPZ, REPO)} + scgenept_k562_labels.csv")


def fit_probe(X_tr, err_tr):
    """就地同规格探针：y = err > 训练中位数；StandardScaler + L2 logistic。"""
    y_tr = (err_tr > np.median(err_tr)).astype(int)
    if len(np.unique(y_tr)) < 2:
        return None
    scaler = StandardScaler().fit(X_tr)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(scaler.transform(X_tr), y_tr)
    return scaler, clf


def run_routing(dry_run=False):
    # GEARS：与 step6d.load_dataset 完全同路径（特征 + 标签）
    feat_g, feature_names, _y_g, pred_df_g, _ = load_gears_dataset(DATA_NAME, "outputs")
    X_g = feat_g[feature_names].values.astype(np.float64)
    perts_g = feat_g.index.tolist()
    err_g = 1.0 - pred_df_g["pearson_delta"].values  # 与 perts_g 对齐（step6d loc 过）

    # scGenePT：与 step6f.load_crossmodel 同路径（含 NaN 行剔除留痕）
    X_s, _y_s, fn_s, perts_s, extra = load_crossmodel(SG_NPZ, SG_CSV)
    assert fn_s == feature_names, f"特征列不一致: {fn_s} vs {feature_names}"
    lab_s = pd_read_labels()
    err_s = 1.0 - lab_s.loc[perts_s, "pearson_delta"].values
    pear_s = lab_s.loc[perts_s, "pearson_delta"].values

    # 对齐到共同扰动名单（两边都已剔 NaN；按 GEARS 名单取交集保序）
    common = [p for p in perts_g if p in set(perts_s)]
    idx_g = [perts_g.index(p) for p in common]
    idx_s = [perts_s.index(p) for p in common]
    X_g, err_g = X_g[idx_g], err_g[idx_g]
    X_s, err_s = X_s[idx_s], err_s[idx_s]
    pear_g = pred_df_g["pearson_delta"].values[[perts_g.index(p) for p in common]]
    n = len(common)
    print(f"[data] 共同扰动 n={n}（GEARS 特征空间 5000 基因 / "
          f"scGenePT {extra['n_genes_feature_space']} 基因）")

    if dry_run:
        return None

    rng_master = np.random.default_rng(20261016)
    seed_results = []
    boot_diffs_all = []  # 各 seed 的 bootstrap 差值（路由 − 最优单模型）
    for sp_seed in SPLIT_SEEDS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        perm = rng.permutation(n)
        n_te = int(round(TEST_FRAC * n))
        te, tr = perm[:n_te], perm[n_te:]

        probes = {}
        for name, X, err in (("gears", X_g, err_g), ("scgenept", X_s, err_s)):
            fit = fit_probe(X[tr], err[tr])
            probes[name] = (fit[0].transform(X[te]) @ fit[1].coef_[0]
                            + fit[1].intercept_[0]) if fit else None
        if probes["gears"] is None or probes["scgenept"] is None:
            raise RuntimeError("训练集标签无两类，无法拟合探针")

        # 预测误差概率（logistic 决策分转概率，单调等价于分数比较）
        def sigmoid(z):
            return 1.0 / (1.0 + np.exp(-z))
        p_err_g = sigmoid(probes["gears"])
        p_err_s = sigmoid(probes["scgenept"])

        pear_g_te, pear_s_te = pear_g[te], pear_s[te]
        route_to_s = p_err_s < p_err_g  # 路由：选预测误差较低者；平手取 GEARS
        routed = np.where(route_to_s, pear_s_te, pear_g_te)
        oracle = np.maximum(pear_g_te, pear_s_te)
        mean_g, mean_s = float(pear_g_te.mean()), float(pear_s_te.mean())
        best_single = max(mean_g, mean_s)
        best_name = "gears" if mean_g >= mean_s else "scgenept"

        # 配对 bootstrap：扰动为抽样单位
        boot = []
        for _ in range(N_BOOT):
            b = rng.integers(0, n_te, n_te)
            bs = max(pear_g_te[b].mean(), pear_s_te[b].mean())
            boot.append(routed[b].mean() - bs)
        boot = np.array(boot)
        p_val = float((boot <= 0).mean())  # 单侧：路由不优于最优单模型的概率
        boot_diffs_all.append(boot)

        seed_results.append({
            "split_seed": sp_seed, "n_train": int(len(tr)), "n_test": int(len(te)),
            "mean_pearson_delta": {
                "all_gears": mean_g, "all_scgenept": mean_s,
                "routed": float(routed.mean()), "oracle": float(oracle.mean())},
            "best_single": best_name,
            "routed_minus_best_single": float(routed.mean() - best_single),
            "routed_minus_oracle": float(routed.mean() - oracle.mean()),
            "routed_to_scgenept_frac": float(route_to_s.mean()),
            "oracle_to_scgenept_frac": float((pear_s_te > pear_g_te).mean()),
            "probe_test_auroc": {
                # 探针自身在 30% 上预测"本模型是否后一半"的 AUROC（诊断）
                "gears": float(auc((err_g[te] > np.median(err_g[tr])),
                                   p_err_g)),
                "scgenept": float(auc((err_s[te] > np.median(err_s[tr])),
                                      p_err_s))},
            "bootstrap": {"n": N_BOOT,
                          "diff_mean": float(boot.mean()),
                          "diff_std": float(boot.std(ddof=1)),
                          "p_one_sided": p_val},
        })
        print(f"[seed {sp_seed}] test n={n_te} | GEARS={mean_g:.4f} "
              f"scGenePT={mean_s:.4f} routed={routed.mean():.4f} "
              f"oracle={oracle.mean():.4f} | routed-best={routed.mean()-best_single:+.4f} "
              f"p={p_val:.3f} | route→sg {route_to_s.mean():.2%}")

    boot_pool = np.concatenate(boot_diffs_all)
    pooled = {"n_resamples": int(len(boot_pool)),
              "diff_mean": float(boot_pool.mean()),
              "diff_std": float(boot_pool.std(ddof=1)),
              "p_one_sided": float((boot_pool <= 0).mean())}

    def mstd(key):
        v = [r["mean_pearson_delta"][key] for r in seed_results]
        return {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)),
                "per_seed": [float(x) for x in v]}

    res = {
        "task": "cross-model per-perturbation routing via error probes (step7)",
        "dataset": DATA_NAME, "n_perturbations": n,
        "split": f"perturbation-level {int(round((1-TEST_FRAC)*100))}/"
                 f"{int(round(TEST_FRAC*100))}, split_seeds={SPLIT_SEEDS}",
        "probe_spec": "StandardScaler + LogisticRegression(C=1.0, max_iter=2000), "
                      "in-situ per model on the 70% train; y = (1 - pearson_delta) "
                      "> train median (step6 label)",
        "routing_rule": "per perturbation, route to the model whose probe predicts "
                        "lower error probability; ties -> GEARS",
        "metric": "mean top-50 |delta| Pearson (pearson_delta, step1 metric) on "
                  "the 30% held-out routing test set",
        "gears_source": "recovered from archive outputs/pred_expr_seed0.npz "
                        "(same artifact as step6 probe; no recompute)",
        "scgenept_note": "138/217 perturbations have target NOT in model input "
                         "space (pred ~ ctrl, data fact kept as-is)",
        "groups": {k: mstd(k) for k in
                   ["all_gears", "all_scgenept", "routed", "oracle"]},
        "routed_minus_best_single": {
            "mean": float(np.mean([r["routed_minus_best_single"]
                                   for r in seed_results])),
            "std": float(np.std([r["routed_minus_best_single"]
                                 for r in seed_results], ddof=1)),
            "per_seed": [r["routed_minus_best_single"] for r in seed_results]},
        "routed_minus_oracle": {
            "mean": float(np.mean([r["routed_minus_oracle"] for r in seed_results])),
            "std": float(np.std([r["routed_minus_oracle"] for r in seed_results],
                                ddof=1))},
        "routed_to_scgenept_frac": {
            "mean": float(np.mean([r["routed_to_scgenept_frac"]
                                   for r in seed_results])),
            "per_seed": [r["routed_to_scgenept_frac"] for r in seed_results]},
        "oracle_to_scgenept_frac": {
            "mean": float(np.mean([r["oracle_to_scgenept_frac"]
                                   for r in seed_results]))},
        "probe_test_auroc": {
            m: {"mean": float(np.mean([r["probe_test_auroc"][m]
                                        for r in seed_results])),
                "per_seed": [r["probe_test_auroc"][m] for r in seed_results]}
            for m in ["gears", "scgenept"]},
        "bootstrap_pooled": pooled,
        "bootstrap_per_seed": [r["bootstrap"] for r in seed_results],
        "seed_results": seed_results,
    }

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, "routing.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[SAVED] {out_json}")

    make_figure(res, seed_results, boot_pool,
                os.path.join(ANALYSIS_DIR, "routing.png"))
    return res


def pd_read_labels():
    import pandas as pd
    return pd.read_csv(SG_CSV).set_index("perturbation")


def auc(y, s):
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(y, s)


def make_figure(res, seed_results, boot_pool, path):
    groups = [("all_gears", "All GEARS", "#1f77b4"),
              ("all_scgenept", "All scGenePT", "#d62728"),
              ("routed", "Probe-routed", "#2ca02c"),
              ("oracle", "Oracle (upper bound)", "#9467bd")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    vals = [res["groups"][k]["mean"] for k, _c, _l in groups]
    errs = [res["groups"][k]["std"] for k, _c, _l in groups]
    bars = ax1.bar(range(len(groups)), vals,
                   color=[g[2] for g in groups], alpha=0.85,
                   edgecolor="gray", lw=0.5,
                   yerr=errs, capsize=4,
                   error_kw=dict(lw=1.0, ecolor="black"))
    for b, v in zip(bars, vals):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.008, f"{v:.3f}",
                 ha="center", fontsize=9)
    ax1.set_xticks(range(len(groups)))
    ax1.set_xticklabels([g[1] for g in groups], fontsize=9)
    ax1.set_ylabel("mean top-50 pearson_delta (30% held-out)")
    ax1.set_title("Per-perturbation cross-model routing")
    ax1.grid(axis="y", alpha=0.25)
    ax1.text(0.99, 0.02, "error bars: std over 3 split seeds",
             transform=ax1.transAxes, ha="right", fontsize=7, color="gray")

    ax2.hist(boot_pool, bins=40, color="#2ca02c", alpha=0.75,
             edgecolor="white")
    ax2.axvline(0.0, color="red", ls="--", lw=1.2)
    ax2.axvline(boot_pool.mean(), color="black", ls=":", lw=1.2)
    p = res["bootstrap_pooled"]["p_one_sided"]
    ax2.set_xlabel("routed − best single model (bootstrap resamples)")
    ax2.set_ylabel("count")
    ax2.set_title(f"Paired bootstrap (perturbation unit, 3 seeds × 1000)\n"
                  f"mean={boot_pool.mean():+.4f}, one-sided p={p:.3f}")
    ax2.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="只检查产物齐备性并打印分析计划，不跑分析")
    args = ap.parse_args()
    check_products()
    if args.dry_run:
        print("\n[DRY-RUN] 产物齐备。分析计划：")
        print(f"  - 扰动级 70/30 划分，seeds={SPLIT_SEEDS}，30% 为路由测试集")
        print("  - 70% 上就地对 GEARS / scGenePT 特征各训同规格探针"
              "（StandardScaler + L2 logistic C=1.0）")
        print("  - 30% 上逐扰动路由到预测误差较低模型；"
              "对照：全 GEARS / 全 scGenePT / oracle")
        print(f"  - 配对 bootstrap {N_BOOT} 次/种子，扰动为抽样单位")
        print("  - 产物：outputs_crossmodel/analysis/routing.json + routing.png")
        return
    run_routing()


if __name__ == "__main__":
    main()
