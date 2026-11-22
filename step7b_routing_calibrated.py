# step7b_routing_calibrated.py — 挽救实验：共同锚点（绝对 pearson_delta 回归）校准后的
# 跨模型逐扰动路由复测。
#
# 背景：step7_routing.py 发现失败机制——二分类探针分数是"模型内相对排序"
# （各自的中位数阈值对应完全不同的绝对误差水平），跨模型直接比较分数会
# 把 ~58% 扰动路由到整体差一个量级的 scGenePT（oracle 只需 9.7%）。
#
# 预注册式设计（先写死再跑）：
#   1. 数据与划分：沿用 step7 的 K562 217 扰动、3 个 70/30 划分种子
#      （split 生成逻辑与 step7 完全一致，保证基线可比）；GEARS 预测 =
#      outputs/pred_expr_seed0.npz，scGenePT =
#      outputs_crossmodel/scgenept/scgenept_k562_predictions.npz。
#   2. 统一误差口径：scGenePT 的 full 空间有 187/5000 基因整列 NaN，故在
#      4813 基因公共空间（两模型基因顺序一致的基因名交集子集）上用 h5ad
#      真值（step1 口径）对**两模型同口径重算** top-50 |delta| Pearson；
#      同时与两份归档标签 csv 比对留痕（一致性检查）。
#   3. 校准器（在 70% 训练侧，9 特征 → 绝对 pearson_delta 连续回归）：
#      a) Ridge（alpha=1.0，直接拟合训练集、预测测试集，无泄漏）；
#      b) isotonic 单调校准：Ridge 在训练侧做 5-fold OOF 预测，
#         isotonic 拟合 (OOF 预测 → 训练真值)，测试侧先 Ridge 再过
#         isotonic（out_of_bounds="clip"）——OOF 输入保证 isotonic 校准
#         本身不接触训练点的 in-sample 拟合值，防泄漏。
#   4. 路由规则：30% 上选校准后预测 pearson_delta 更高的模型。
#   5. 对照组：全 GEARS、全 scGenePT、未校准探针路由（step7 基线复算）、oracle。
#   6. 指标：mean top-50 pearson_delta，3 seeds 均值±std；校准路由 vs 全 GEARS
#      配对 bootstrap（1000 次/种子，pooled）；校准路由 vs oracle 差距。
#   7. 机制验证：校准路由实际路由给 scGenePT 的比例（oracle=9.7%，
#      step7 未校准=57.9%）是否回到合理区间。
#
# 产物：outputs_crossmodel/analysis/routing_calibrated.json + _calibrated.png（300dpi 英文）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7b_routing_calibrated.py --dry-run
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7b_routing_calibrated.py
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES
from step6d_naive_baselines import load_dataset as load_gears_dataset
from step6f_crossmodel import load_crossmodel

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")

SPLIT_SEEDS = [0, 1, 2]
TEST_FRAC = 0.30
N_BOOT = 1000
RIDGE_ALPHA = 1.0
ISO_FOLDS = 5
BOOT_MASTER_SEED = 20261016  # 与 step7 相同 → 相同 70/30 划分

GEARS_NPZ = os.path.join(REPO, "outputs", "pred_expr_seed0.npz")
GEARS_CSV = os.path.join(REPO, "outputs", "predictions_seed0.csv")
SG_NPZ = os.path.join(REPO, "outputs_crossmodel", "scgenept",
                      "scgenept_k562_predictions.npz")
SG_CSV = os.path.join(REPO, "outputs_crossmodel", "scgenept",
                      "scgenept_k562_labels.csv")
H5AD = os.path.join(REPO, "data", "gears_demo",
                    "replogle_k562_essential", "perturb_processed.h5ad")
TOP_N = 50


# ---------------------------------------------------------------- 数据加载
def load_all():
    """加载两模型特征与预测，统一口径重算 pearson_delta（h5ad 真值，
    4813 基因公共空间）。返回统一的加载字典。"""
    # 特征（与 step6/step7 同路径）
    feat_g, feature_names, _y, pred_df_g, _ = load_gears_dataset(
        "replogle_k562_essential", "outputs")
    X_g = feat_g[feature_names].values.astype(np.float64)
    perts_g = feat_g.index.tolist()
    X_s, _ys, fn_s, perts_s, extra = load_crossmodel(SG_NPZ, SG_CSV)
    assert fn_s == feature_names

    common = [p for p in perts_g if p in set(perts_s)]
    idx_g = [perts_g.index(p) for p in common]
    idx_s = [perts_s.index(p) for p in common]

    # 原始预测谱（5000 基因，两模型基因顺序一致，见 step7b dry-run 检查）
    zg = np.load(GEARS_NPZ, allow_pickle=True)
    zs = np.load(SG_NPZ, allow_pickle=True)
    g_genes = [str(g) for g in zg["pred_genes"]]
    s_genes = [str(g) for g in zs["full_gene_space"]]
    assert g_genes == s_genes, "两模型 full 空间基因顺序不一致，需按名对齐"
    pred_g_full = np.asarray(zg["pred_expr"], dtype=np.float64)[idx_g]
    pred_s_full = np.asarray(zs["pred_expr_full"], dtype=np.float64)[idx_s]
    nan_col = np.isnan(pred_s_full).any(axis=0)  # 187 个整列 NaN 基因
    keep = ~nan_col
    print(f"[data] 公共空间 {int(keep.sum())} 基因（剔除整列 NaN "
          f"{int(nan_col.sum())}）")

    # h5ad 真值（step1 口径）：ctrl 均值 + 各 condition 均值，对齐基因顺序
    import anndata as ad
    adata = ad.read_h5ad(H5AD, backed="r")
    var_genes = adata.var["gene_name"].tolist()
    obs_cond = adata.obs["condition"].astype(str).values
    ctrl_mean = np.asarray(
        adata[obs_cond == "ctrl"].X.mean(axis=0)).ravel()
    gene_pos = {g: i for i, g in enumerate(var_genes)}
    align = np.array([gene_pos[g] for g in g_genes])
    ctrl_mean = ctrl_mean[align]
    true_means = np.empty((len(common), len(g_genes)))
    cond_set = set(obs_cond)
    missing = []
    for i, g in enumerate(common):
        cond = next((c for c in (f"ctrl+{g}", f"{g}+ctrl") if c in cond_set),
                    None)
        if cond is None:
            missing.append(g)
            true_means[i] = np.nan
            continue
        true_means[i] = np.asarray(
            adata[obs_cond == cond].X.mean(axis=0)).ravel()[align]
    adata.file.close()
    if missing:
        print(f"[data] 警告：{len(missing)} 个扰动无真值 condition: {missing[:5]}")

    # 统一口径重算 top-50 |delta| Pearson（公共 4813 基因空间）
    pear_g = pearson_top50(pred_g_full[:, keep], true_means[:, keep],
                           ctrl_mean[keep])
    pear_s = pearson_top50(pred_s_full[:, keep], true_means[:, keep],
                           ctrl_mean[keep])

    # 口径一致性留痕：与归档标签 csv 比对
    lab_g = pred_df_g["pearson_delta"].values[[perts_g.index(p) for p in common]]
    lab_s = pd.read_csv(SG_CSV).set_index("perturbation").loc[
        common, "pearson_delta"].values
    consistency = {
        "gears_recomputed_vs_archive": label_diff(pear_g, lab_g),
        "scgenept_recomputed_vs_archive": label_diff(pear_s, lab_s),
        "note": "两模型统一在 4813-gene 公共空间、同一真值重算；"
                "归档标签各自在不同基因空间计算（GEARS 5000 / "
                "scGenePT matched 4813），允许小幅口径差",
    }
    ok = ~np.isnan(pear_g) & ~np.isnan(pear_s)
    print(f"[data] 共同扰动 n={len(common)}，有效 {int(ok.sum())}，"
          f"NaN pearson {int((~ok).sum())}")

    return {
        "X_g": X_g[idx_g][ok], "X_s": X_s[idx_s][ok],
        "pear_g": pear_g[ok], "pear_s": pear_s[ok],
        "perts": np.array(common)[ok], "feature_names": feature_names,
        "consistency": consistency,
        "n_feature_space": int(keep.sum()),
    }


def pearson_top50(pred, true, ctrl):
    """top-N |true delta| 基因上的 Pearson（step1_run.py:160-162 口径）。"""
    p_delta = pred - ctrl
    t_delta = true - ctrl
    out = np.full(len(pred), np.nan)
    for i in range(len(pred)):
        top = np.argsort(-np.abs(t_delta[i]))[:TOP_N]
        if np.std(t_delta[i][top]) == 0 or np.std(p_delta[i][top]) == 0:
            continue
        out[i] = np.corrcoef(p_delta[i][top], t_delta[i][top])[0, 1]
    return out


def label_diff(a, b):
    m = ~np.isnan(a) & ~np.isnan(b)
    return {"n": int(m.sum()),
            "pearson_corr": float(np.corrcoef(a[m], b[m])[0, 1]),
            "mean_abs_diff": float(np.mean(np.abs(a[m] - b[m]))),
            "max_abs_diff": float(np.max(np.abs(a[m] - b[m])))}


# ---------------------------------------------------------------- 校准器
def ridge_predict(X_tr, y_tr, X_te):
    """Ridge（alpha=1.0）在训练侧拟合，预测测试侧。训练侧内部再标准化
    （特征已由探针口径构造，此处直接缩放进 Ridge 管线最简单可靠）。"""
    sc = StandardScaler().fit(X_tr)
    m = Ridge(alpha=RIDGE_ALPHA).fit(sc.transform(X_tr), y_tr)
    return m.predict(sc.transform(X_te)), (sc, m)


def iso_calibrate(X_tr, y_tr, X_te, fold_seed):
    """isotonic 单调校准：Ridge 训练侧 5-fold OOF 预测 → isotonic 拟合
    (OOF 预测 → 训练真值)；测试侧 Ridge 预测后过 isotonic（clip 外推）。"""
    sc = StandardScaler().fit(X_tr)
    Xs = sc.transform(X_tr)
    kf = KFold(n_splits=ISO_FOLDS, shuffle=True, random_state=fold_seed)
    oof = np.full(len(y_tr), np.nan)
    for tr_i, va_i in kf.split(Xs):
        m = Ridge(alpha=RIDGE_ALPHA).fit(Xs[tr_i], y_tr[tr_i])
        oof[va_i] = m.predict(Xs[va_i])
    iso = IsotonicRegression(out_of_bounds="clip").fit(oof, y_tr)
    m_full = Ridge(alpha=RIDGE_ALPHA).fit(Xs, y_tr)
    return iso.predict(m_full.predict(sc.transform(X_te))), iso, oof


def uncalibrated_route(pairs):
    """step7 基线复算：同规格二分类探针（y=err>训练中位数）+ sigmoid
    概率比较，选预测误差较低者。pairs = [(X_tr, err_tr, X_te) per model]，
    返回 route_to_scgenept 布尔（平手取 GEARS）。"""
    probs = []
    for Xtr, etr, Xte in pairs:
        y_tr = (etr > np.median(etr)).astype(int)
        sc = StandardScaler().fit(Xtr)
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(
            sc.transform(Xtr), y_tr)
        z = sc.transform(Xte) @ clf.coef_[0] + clf.intercept_[0]
        probs.append(1.0 / (1.0 + np.exp(-z)))
    p_g, p_s = probs
    return p_s < p_g


# ---------------------------------------------------------------- 主分析
def run(dry_run=False):
    for p in (GEARS_NPZ, GEARS_CSV, SG_NPZ, SG_CSV, H5AD):
        if not os.path.exists(p):
            print(f"[MISSING] {p}")
            raise SystemExit(1)
    print(f"[OK] GEARS npz + csv；scGenePT npz + csv；h5ad 真值")

    d = load_all()
    if dry_run:
        return None
    X_g, X_s = d["X_g"], d["X_s"]
    pear_g, pear_s = d["pear_g"], d["pear_s"]
    n = len(pear_g)
    print(f"[data] n={n} 扰动 × {d['n_feature_space']} 基因公共空间；"
          f"全量 mean pearson: GEARS={pear_g.mean():.4f} "
          f"scGenePT={pear_s.mean():.4f}")

    rng_master = np.random.default_rng(BOOT_MASTER_SEED)  # 与 step7 同划分
    seed_results = []
    boot_store = {"ridge": [], "iso": []}
    for sp_seed in SPLIT_SEEDS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        perm = rng.permutation(n)
        n_te = int(round(TEST_FRAC * n))
        te, tr = perm[:n_te], perm[n_te:]

        pg_te, ps_te = pear_g[te], pear_s[te]

        # --- 校准器（两模型分别） ---
        def cal_pair(X):
            rid, _ = ridge_predict(X[tr], pear_g_or_s_tr, X[te])
            isoc, _, _ = iso_calibrate(X[tr], pear_g_or_s_tr, X[te], sp_seed)
            return rid, isoc

        preds = {}
        for name, X in (("gears", X_g), ("scgenept", X_s)):
            pear_g_or_s_tr = pear_g[tr] if name == "gears" else pear_s[tr]
            rid, isoc = cal_pair(X)
            preds[name] = {"ridge": rid, "iso": isoc}

        def routed_pearson(score_g, score_s):
            to_s = score_s > score_g  # 选校准后预测 pearson 更高者
            return np.where(to_s, ps_te, pg_te), to_s

        r_ridge, to_s_ridge = routed_pearson(preds["gears"]["ridge"],
                                             preds["scgenept"]["ridge"])
        r_iso, to_s_iso = routed_pearson(preds["gears"]["iso"],
                                         preds["scgenept"]["iso"])
        r_unc = uncalibrated_route(
            [(X_g[tr], 1 - pear_g[tr], X_g[te]),
             (X_s[tr], 1 - pear_s[tr], X_s[te])])
        r_unc_pear = np.where(r_unc, ps_te, pg_te)
        oracle = np.maximum(pg_te, ps_te)

        mean_g, mean_s = float(pg_te.mean()), float(ps_te.mean())
        groups = {"all_gears": mean_g, "all_scgenept": mean_s,
                  "probe_uncalibrated": float(r_unc_pear.mean()),
                  "ridge_routed": float(r_ridge.mean()),
                  "isotonic_routed": float(r_iso.mean()),
                  "oracle": float(oracle.mean())}

        boots = {}
        for key, arr in (("ridge", r_ridge), ("iso", r_iso)):
            b = []
            for _ in range(N_BOOT):
                bi = rng.integers(0, n_te, n_te)
                b.append(arr[bi].mean() - pg_te[bi].mean())
            b = np.array(b)
            boots[key] = {"diff_mean": float(b.mean()),
                          "diff_std": float(b.std(ddof=1)),
                          "p_one_sided": float((b <= 0).mean())}
            boot_store[key].append(b)

        seed_results.append({
            "split_seed": sp_seed, "n_train": int(len(tr)),
            "n_test": int(len(te)),
            "mean_pearson_delta": groups,
            "routed_to_scgenept_frac": {
                "ridge_routed": float(to_s_ridge.mean()),
                "isotonic_routed": float(to_s_iso.mean()),
                "probe_uncalibrated": float(r_unc.mean()),
                "oracle": float((ps_te > pg_te).mean())},
            "calibrated_minus_all_gears": {
                "ridge_routed": groups["ridge_routed"] - mean_g,
                "isotonic_routed": groups["isotonic_routed"] - mean_g},
            "calibrated_minus_oracle": {
                "ridge_routed": groups["ridge_routed"] - groups["oracle"],
                "isotonic_routed": groups["isotonic_routed"] - groups["oracle"]},
            "bootstrap_vs_all_gears": boots,
        })
        print(f"[seed {sp_seed}] test n={n_te} | GEARS={mean_g:.4f} "
              f"sg={mean_s:.4f} | uncal={groups['probe_uncalibrated']:.4f} "
              f"ridge={groups['ridge_routed']:.4f} iso={groups['isotonic_routed']:.4f} "
              f"oracle={groups['oracle']:.4f} | →sg% "
              f"uncal={r_unc.mean():.1%} "
              f"ridge={to_s_ridge.mean():.1%} iso={to_s_iso.mean():.1%} "
              f"oracle={(ps_te > pg_te).mean():.1%} | p(ridge)="
              f"{boots['ridge']['p_one_sided']:.3f} "
              f"p(iso)={boots['iso']['p_one_sided']:.3f}")

    def agg(key):
        v = [r["mean_pearson_delta"][key] for r in seed_results]
        return {"mean": float(np.mean(v)), "std": float(np.std(v, ddof=1)),
                "per_seed": [float(x) for x in v]}

    pooled = {}
    for key in ("ridge", "iso"):
        bp = np.concatenate(boot_store[key])
        pooled[key] = {"n_resamples": int(len(bp)),
                       "diff_mean": float(bp.mean()),
                       "diff_std": float(bp.std(ddof=1)),
                       "p_one_sided": float((bp <= 0).mean())}

    def frac_agg(name):
        v = [r["routed_to_scgenept_frac"][name] for r in seed_results]
        return {"mean": float(np.mean(v)), "per_seed": [float(x) for x in v]}

    res = {
        "task": "calibrated cross-model per-perturbation routing (step7b)",
        "dataset": "replogle_k562_essential", "n_perturbations": n,
        "metric_space": "unified 4813-gene public space, top-50 |delta| "
                        "Pearson recomputed from h5ad truth for BOTH models "
                        "(step1 metric)",
        "label_consistency_check": d["consistency"],
        "split": "perturbation-level 70/30, split_seeds=[0,1,2]; identical "
                 "split generation as step7_routing.py (master seed "
                 f"{BOOT_MASTER_SEED})",
        "calibrators": {
            "ridge": f"StandardScaler + Ridge(alpha={RIDGE_ALPHA}) regressing "
                     "absolute pearson_delta, fit on 70% train, predict 30% "
                     "test (no leakage)",
            "isotonic": "Ridge 5-fold OOF predictions on the 70% train "
                        "(shuffle, split-seed) -> IsotonicRegression(OOF_pred "
                        "-> train truth, out_of_bounds='clip'); test = "
                        "iso(Ridge_pred). OOF input prevents the calibrator "
                        "from seeing in-sample Ridge fits (leakage guard).",
        },
        "routing_rule": "per perturbation, route to the model with higher "
                        "calibrated predicted pearson_delta",
        "gears_source": "recovered from archive outputs/pred_expr_seed0.npz "
                        "(no recompute); scGenePT pred_expr_full from "
                        "scgenept_k562_predictions.npz",
        "groups": {k: agg(k) for k in
                   ["all_gears", "all_scgenept", "probe_uncalibrated",
                    "ridge_routed", "isotonic_routed", "oracle"]},
        "routed_to_scgenept_frac": {k: frac_agg(k) for k in
                                    ["probe_uncalibrated", "ridge_routed",
                                     "isotonic_routed", "oracle"]},
        "calibrated_minus_all_gears": {
            k: {"mean": float(np.mean(
                    [r["calibrated_minus_all_gears"][k] for r in seed_results])),
                "std": float(np.std(
                    [r["calibrated_minus_all_gears"][k] for r in seed_results],
                    ddof=1))}
            for k in ["ridge_routed", "isotonic_routed"]},
        "calibrated_minus_oracle": {
            k: {"mean": float(np.mean(
                    [r["calibrated_minus_oracle"][k] for r in seed_results])),
                "std": float(np.std(
                    [r["calibrated_minus_oracle"][k] for r in seed_results],
                    ddof=1))}
            for k in ["ridge_routed", "isotonic_routed"]},
        "bootstrap_pooled_vs_all_gears": pooled,
        "bootstrap_per_seed_vs_all_gears": [
            r["bootstrap_vs_all_gears"] for r in seed_results],
        "seed_results": seed_results,
    }

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, "routing_calibrated.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[SAVED] {out_json}")
    make_figure(res, boot_store, os.path.join(
        ANALYSIS_DIR, "routing_calibrated.png"))
    return res


def make_figure(res, boot_store, path):
    order = [("all_gears", "All GEARS", "#1f77b4"),
             ("all_scgenept", "All scGenePT", "#d62728"),
             ("probe_uncalibrated", "Probe-routed\n(uncalibrated, step7)", "#8c564b"),
             ("ridge_routed", "Ridge-calibrated\nrouting", "#2ca02c"),
             ("isotonic_routed", "Isotonic-calibrated\nrouting", "#17becf"),
             ("oracle", "Oracle (upper\nbound)", "#9467bd")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.4))

    vals = [res["groups"][k]["mean"] for k, _l, _c in order]
    errs = [res["groups"][k]["std"] for k, _l, _c in order]
    colors = [c for _k, _l, c in order]
    bars = ax1.bar(range(len(order)), vals, color=colors, alpha=0.88,
                   edgecolor="gray", lw=0.5, yerr=errs, capsize=4,
                   error_kw=dict(lw=1.0, ecolor="black"))
    for b, v in zip(bars, vals):
        ax1.text(b.get_x() + b.get_width() / 2, v + 0.008, f"{v:.3f}",
                 ha="center", fontsize=8.5)
    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels([l for _k, l, _c in order], fontsize=8)
    ax1.set_ylabel("mean top-50 pearson_delta (30% held-out)")
    ax1.set_title("Calibrated per-perturbation cross-model routing")
    ax1.grid(axis="y", alpha=0.25)
    ax1.text(0.99, 0.02, "error bars: std over 3 split seeds",
             transform=ax1.transAxes, ha="right", fontsize=7, color="gray")

    colors2 = {"ridge": "#2ca02c", "iso": "#17becf"}
    labels2 = {"ridge": "Ridge-calibrated", "iso": "Isotonic-calibrated"}
    for key in ("ridge", "iso"):
        bp = np.concatenate(boot_store[key])
        ax2.hist(bp, bins=40, color=colors2[key], alpha=0.55,
                 edgecolor="white",
                 label=f"{labels2[key]} (p={res['bootstrap_pooled_vs_all_gears'][key]['p_one_sided']:.3f})")
    ax2.axvline(0.0, color="red", ls="--", lw=1.2)
    ax2.set_xlabel("calibrated routed − all GEARS (bootstrap resamples)")
    ax2.set_ylabel("count")
    ax2.set_title("Paired bootstrap vs all GEARS\n(perturbation unit, "
                  "3 seeds × 1000)")
    ax2.grid(alpha=0.25)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="只检查产物齐备性并打印分析计划，不跑分析")
    args = ap.parse_args()
    if args.dry_run:
        for p in (GEARS_NPZ, GEARS_CSV, SG_NPZ, SG_CSV, H5AD):
            if not os.path.exists(p):
                print(f"[MISSING] {p}")
                raise SystemExit(1)
        print("[DRY-RUN] 产物齐备。分析计划：")
        print("  - 统一口径：4813 基因公共空间，两模型 pearson_delta 均从 "
              "h5ad 真值重算，与归档标签比对留痕")
        print("  - 划分：与 step7 完全相同的 3 个 70/30（master seed "
              f"{BOOT_MASTER_SEED}），保证基线可比")
        print(f"  - 校准器：Ridge(alpha={RIDGE_ALPHA})；isotonic（Ridge "
              f"{ISO_FOLDS}-fold OOF 输入防泄漏）")
        print("  - 对照：全 GEARS / 全 scGenePT / 未校准探针路由 / oracle")
        print("  - 产物：outputs_crossmodel/analysis/routing_calibrated."
              "json + routing_calibrated.png")
        return
    run(dry_run=False)


if __name__ == "__main__":
    main()
