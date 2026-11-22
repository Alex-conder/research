# step13_s4_s3_fixes.py — 红队挂账实算：S4 闭环配对检验 + S3 预测缩放对照
#
# 背景（redteam_20261024/gating.md）：
#   S4（建议级）：gate-vs-random 差 0.0065 无正式检验 -> 本脚本构造
#     random-rejection 基线（与门同 coverage 下随机拒答，3 个随机种子），
#     对 10 个可执行组合逐点配对 bootstrap（B=1000）报差值 CI；
#     判定：门的 risk 显著低于随机（CI 不含 0）的组合占比。
#     配对口径：同一 bootstrap 重采样内，随机拒答按 3 个种子各抽一次
#     （每 coverage 取 k=round(coverage x n) 个均匀随机分），三门分数在
#     同一重采样上重算，差值逐重采样取——公共基线噪声在对差中抵消。
#     稳健性复核：期望随机臂（= 重采样基率，随机拒答的期望精度）配对
#     给出同 verdict。
#   S3（建议级）：探针平凡性检验缺"预测缩放对照" -> 对 GEARS-K562/RPE1，
#     探针输入特征统一缩放（StandardScaler，全数据一次 fit）后重训重评
#     AUROC，与原口径（fold 内 StandardScaler 管线）对比；
#     补充：预测 delta 乘常数因子（0.5/2.0）重算特征后的 AUROC
#     （红队原文"恒定/缩放预测对照"，pearson 标签对正缩放不变）。
#
# 口径复用（不自造逻辑）：
#   - 组合加载/探针/符号层/两层门：step10_signature_resolution + step10_two_layer_gate
#   - GEARS 特征与标签：step6d_naive_baselines.load_dataset / step6_gating_probe
#   - risk 口径：risk = 1 - precision（top-coverage 带内最差一半占比，step6b 同）
#
# 随机源（专属，全部注明）：
#   S4_BOOT_MASTER_SEED = 20261103（逐组合由其派生 Generator）
#   S4_RANDOM_BASELINE_SEEDS = [910101, 910102, 910103]（随机拒答基线 3 种子；
#     点估计基线曲线每 种子 x coverage 抽签 S4_BASELINE_DRAWS 次；
#     配对 bootstrap 内每重采样每种子逐 coverage 各抽 1 次）
#   S3 复用 step6 原协议 seed=0（CV 划分与 LogisticRegression 同源）
#
# 产物：outputs_rebuild/s4_s3_fixes.json / s4_s3_fixes.md
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step13_s4_s3_fixes.py
import json
import os
import sys
import time

import numpy as np

for _s in ("stdout", "stderr"):  # GBK 控制台纪律：不崩、不丢字符
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")

from step10_signature_resolution import COMBOS, fit_coef, load_combo
from step10_two_layer_gate import (FULL_COV, N_FOLDS, OOF_SEED, SIBLING,
                                   oof_probe_scores, precision_at,
                                   symbol_scores, two_layer_scores)
from step6_gating_probe import run_probe
from step6d_naive_baselines import load_dataset

S4_BOOT_MASTER_SEED = 20261103
S4_RANDOM_BASELINE_SEEDS = [910101, 910102, 910103]
S4_BASELINE_DRAWS = 1000
S4_B_BOOT = 1000

S3_SEED = 0  # step6 原协议种子
S3_PRED_SCALE_FACTORS = [0.5, 2.0]
# 存档的 step6 原口径 AUROC（outputs{,_rpe1}/gating_probe.json），用于复现校验
S3_ARCHIVED = {
    "gears_k562": {"cv_auroc_mean": 0.6556473829201102,
                   "pooled_cv_auroc": 0.6546890927624872},
    "gears_rpe1": {"cv_auroc_mean": 0.7965175164758932,
                   "pooled_cv_auroc": 0.7894399589901321},
}

S4_COMBO_KEYS = [k for k, v in COMBOS.items()
                 if v["kind"] != "missing" and SIBLING.get(k) is not None]
S4_ARMS = ["in_situ", "symbol", "two_layer"]


# ------------------------------------------------------- S4 点估计基线曲线
def random_rejection_baseline(y, coverages, seeds, n_draws):
    """随机拒答基线（全数据点估计）：每 (seed, coverage) 做 n_draws 次等
    大概率拒答抽签，precision = 被标记查询中 y=1 占比；risk = 1-precision。
    返回逐 seed 与跨 seed 平均的 risk_by_coverage。"""
    n = len(y)
    per_seed = {}
    for s in seeds:
        g = np.random.default_rng(s)
        row = {}
        for c in coverages:
            k = max(1, int(round(c * n)))
            r = g.random((n_draws, n))
            top = np.argpartition(r, -k, axis=1)[:, -k:]
            prec = y[top].mean()
            row[f"{c:.2f}"] = {"k": int(k),
                               "precision_mean": float(prec),
                               "risk": float(1.0 - prec)}
        per_seed[str(s)] = row
    mean_row = {}
    for c in coverages:
        key = f"{c:.2f}"
        risks = [per_seed[str(s)][key]["risk"] for s in seeds]
        precs = [per_seed[str(s)][key]["precision_mean"] for s in seeds]
        mean_row[key] = {"risk_mean_over_seeds": float(np.mean(risks)),
                         "risk_min_over_seeds": float(np.min(risks)),
                         "risk_max_over_seeds": float(np.max(risks)),
                         "risk_spread_over_seeds": float(np.max(risks) -
                                                         np.min(risks)),
                         "precision_mean": float(np.mean(precs))}
    return {"per_seed": per_seed, "mean_over_seeds": mean_row,
            "n_draws_per_seed_coverage": n_draws}


# ------------------------------------------------------- S4 配对 bootstrap
def gate_arm_scores(Xb, yb, signs):
    """同一重采样上的三门分数（与 step10 arm_scores 同结构）。"""
    s_a = oof_probe_scores(Xb, yb)
    s_b = symbol_scores(Xb, signs)
    s_c = two_layer_scores(s_b, s_a, yb)
    return {"in_situ": s_a, "symbol": s_b, "two_layer": s_c}


def s4_paired_bootstrap(X, y, signs, coverages, b_boot, rng, stream_rngs):
    """扰动级配对 bootstrap：每次重采样重算三门分数；随机拒答臂在同一
    重采样上按 3 个种子各抽一次（每 coverage 取 k=round(c*n) 个均匀随机
    分），主口径 mc3 取种种子平均；复核口径 expected 取重采样基率
    （随机拒答的期望精度，蒙特卡洛噪声为零）。差 = risk_gate - risk_rand。
    """
    n = len(y)
    diffs = {mode: {a: np.zeros((b_boot, len(coverages)))
                    for a in S4_ARMS}
             for mode in ("mc3", "expected")}
    for rep in range(b_boot):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb = X[idx], y[idx]
        sc = gate_arm_scores(Xb, yb, signs)
        rand_scores = np.column_stack([g.random(n) for g in stream_rngs])
        base_rate = yb.mean()
        for j, c in enumerate(coverages):
            k = max(1, int(round(c * n)))
            top = np.argpartition(rand_scores, -k, axis=0)[-k:]
            rand_prec_mc3 = float(yb[top].mean())
            for a, s in sc.items():
                gate_risk = 1.0 - precision_at(s, yb, c)
                diffs["mc3"][a][rep, j] = gate_risk - (1.0 - rand_prec_mc3)
                diffs["expected"][a][rep, j] = gate_risk - (1.0 - base_rate)
    return diffs


def summarize_diffs(diffs, coverages):
    """逐 coverage 差值均值 + 95% percentile CI；全 coverage 均值差的 CI。"""
    out = {}
    for arm, d in diffs.items():
        per_cov = {}
        for j, c in enumerate(coverages):
            lo, hi = np.percentile(d[:, j], [2.5, 97.5])
            per_cov[f"{c:.2f}"] = {
                "diff_mean": float(d[:, j].mean()),
                "ci_lo": float(lo), "ci_hi": float(hi),
                "ci_excludes_zero": bool(lo > 0 or hi < 0)}
        overall = d.mean(axis=1)
        lo, hi = np.percentile(overall, [2.5, 97.5])
        out[arm] = {
            "diff_risk_gate_minus_random_by_coverage": per_cov,
            "mean_diff_over_coverage": float(overall.mean()),
            "ci_lo": float(lo), "ci_hi": float(hi),
            "ci_excludes_zero": bool(lo > 0 or hi < 0),
            "gate_significantly_lower_risk": bool(hi < 0),
        }
    return out


def run_s4_combo(key, coverages, b_boot, rng_master):
    c = COMBOS[key]
    t0 = time.perf_counter()
    X, y, feature_names = load_combo(key)
    sib_key = SIBLING[key]
    Xs, ys, _ = load_combo(sib_key)
    signs = np.sign(fit_coef(Xs, ys))
    t_load = time.perf_counter() - t0

    arms = gate_arm_scores(X, y, signs)
    point_risk = {a: {f"{c:.2f}": float(1.0 - precision_at(s, y, c))
                      for c in coverages}
                  for a, s in arms.items()}

    base = random_rejection_baseline(y, coverages,
                                     S4_RANDOM_BASELINE_SEEDS,
                                     S4_BASELINE_DRAWS)

    rng = np.random.default_rng(rng_master.integers(0, 2**31))
    stream_rngs = [np.random.default_rng(s + b_boot + 1)
                   for s in S4_RANDOM_BASELINE_SEEDS]
    t1 = time.perf_counter()
    diffs = s4_paired_bootstrap(X, y, signs, coverages, b_boot, rng,
                                stream_rngs)
    t_boot = time.perf_counter() - t1
    summ = {mode: summarize_diffs(d, coverages) for mode, d in diffs.items()}

    res = {"engine": c["engine"], "dataset": c["dataset"], "cls": c["cls"],
           "n": int(len(y)), "symbol_signs_from": sib_key,
           "point_risk_by_coverage": point_risk,
           "random_rejection_baseline": base,
           "paired_bootstrap_gate_minus_random": {
               "b_boot": b_boot, "primary_mode": "mc3",
               "modes": summ},
           "timing_s": {"load": round(t_load, 3), "bootstrap": round(t_boot, 3)}}
    m = summ["mc3"]["symbol"]
    e = summ["expected"]["symbol"]
    print(f"[S4][{key}] n={len(y)} symbol(mc3): {m['mean_diff_over_coverage']:+.4f} "
          f"CI=[{m['ci_lo']:+.4f},{m['ci_hi']:+.4f}] win={m['gate_significantly_lower_risk']} "
          f"| symbol(expected): {e['mean_diff_over_coverage']:+.4f} "
          f"CI=[{e['ci_lo']:+.4f},{e['ci_hi']:+.4f}] win={e['gate_significantly_lower_risk']} "
          f"| boot {t_boot:.1f}s")
    return res


def run_s4():
    rng_master = np.random.default_rng(S4_BOOT_MASTER_SEED)
    combos, times = {}, {}
    for key in S4_COMBO_KEYS:
        t0 = time.perf_counter()
        combos[key] = run_s4_combo(key, FULL_COV, S4_B_BOOT, rng_master)
        times[key] = round(time.perf_counter() - t0, 3)
    summary = {}
    for mode in ("mc3", "expected"):
        summary[mode] = {}
        for arm in S4_ARMS:
            wins = sum(1 for v in combos.values()
                       if v["paired_bootstrap_gate_minus_random"]["modes"]
                       [mode][arm]["gate_significantly_lower_risk"])
            summary[mode][arm] = {"n_win": int(wins), "n": len(combos),
                                  "fraction": wins / len(combos)}
    out = {"spec": {
               "coverages": FULL_COV, "b_boot": S4_B_BOOT,
               "n_folds": N_FOLDS, "oof_seed": OOF_SEED,
               "random_baseline_seeds": S4_RANDOM_BASELINE_SEEDS,
               "baseline_draws_per_seed_coverage": S4_BASELINE_DRAWS,
               "boot_master_seed": S4_BOOT_MASTER_SEED,
               "arms": S4_ARMS,
               "gate_product_arm": "symbol (R2 verdict tier=symbol_layer, "
                                   "outputs_rebuild/two_layer_gate.json)",
               "metric": "risk = 1 - precision in top-coverage flagged band "
                         "(step6b caliber); diff = risk_gate - risk_random",
               "excluded": "scf_k562 (in-situ arm only, no sibling product) "
                           "and scf_rpe1 (product missing) per R2 amendment A2",
               "pairing": "both arms recomputed on the same bootstrap "
                          "resample; random arm = 3 fixed seeds x 1 draw per "
                          "coverage per resample (primary, mc3), or the "
                          "resample base rate = expected precision of random "
                          "rejection (robustness, expected)",
               "random_baseline_point_curve": "on the full data, 3 seeds x "
                                              "1000 draws per coverage, "
                                              "averaged over seeds"},
           "combos": combos, "summary": summary, "timing_s": times}
    pm = summary["mc3"]
    print(f"[S4] mc3 wins: in_situ {pm['in_situ']['n_win']}/10, "
          f"symbol {pm['symbol']['n_win']}/10, "
          f"two_layer {pm['two_layer']['n_win']}/10")
    return out


# ---------------------------------------------------------------- S3 探针 AUROC
def probe_cv(X, y, in_fold_scale):
    """step6 run_probe 同协议（5-fold shuffle seed=0, L2 logistic C=1.0），
    in_fold_scale=True 时为原管线（fold 内 StandardScaler），
    False 时为纯 LR（输入特征已在外部统一缩放或不缩放）。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import (StratifiedKFold, cross_val_predict,
                                         cross_validate)
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    if in_fold_scale:
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=1.0, max_iter=2000, random_state=S3_SEED))
    else:
        clf = LogisticRegression(C=1.0, max_iter=2000, random_state=S3_SEED)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=S3_SEED)
    aurocs = cross_validate(clf, X, y, cv=cv, scoring="roc_auc")["test_score"]
    pooled = cross_val_predict(clf, X, y, cv=cv, method="decision_function")
    return {"auroc_mean": float(aurocs.mean()),
            "auroc_std": float(np.std(aurocs, ddof=1)),
            "auroc_per_fold": [float(a) for a in aurocs],
            "pooled_auroc": float(roc_auc_score(y, pooled))}


def load_gears_pred_bundle(data_name, out_dir):
    """复刻 step6d load_dataset 的取数，另返回 pred_expr/ctrl/perts/genes，
    供预测级缩放对照重算特征。"""
    import anndata as ad
    import pandas as pd
    from step6_gating_probe import SPECTRUM_FEATURES

    out_path = os.path.join(REPO, out_dir)
    z = np.load(f"{out_path}/pred_expr_seed{S3_SEED}.npz")
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
        f"{out_path}/predictions_seed{S3_SEED}.csv").set_index("perturbation")
    pred_df = pred_df.loc[perts]

    from step6_gating_probe import compute_features
    feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    feature_names = SPECTRUM_FEATURES + ["log_n_cells"]
    err = 1.0 - pred_df["pearson_delta"].values
    y = (err > np.median(err)).astype(int)
    return {"feat": feat, "feature_names": feature_names, "y": y,
            "pred_expr": pred_expr, "ctrl_aligned": ctrl_aligned,
            "perts": perts, "pred_genes": pred_genes,
            "pred_df": pred_df, "n": int(len(y))}


def run_s3_combo(key, data_name, out_dir):
    from step6_gating_probe import compute_features
    from sklearn.preprocessing import StandardScaler

    b = load_gears_pred_bundle(data_name, out_dir)
    X = b["feat"][b["feature_names"]].values.astype(np.float64)
    y = b["y"]
    names = b["feature_names"]

    original = run_probe(X, y, list(names), seed=S3_SEED)
    original = {k: original[k] for k in
                ("auroc_mean", "auroc_std", "auroc_per_fold", "pooled_auroc")}

    X_uni = StandardScaler().fit_transform(X)  # 统一缩放：全数据一次 fit
    unified = probe_cv(X_uni, y, in_fold_scale=False)
    raw = probe_cv(X, y, in_fold_scale=False)  # 诊断：完全不缩放

    archived = S3_ARCHIVED[key]
    validation = {
        "archived_cv_auroc_mean": archived["cv_auroc_mean"],
        "recomputed_cv_auroc_mean": original["auroc_mean"],
        "abs_diff": abs(original["auroc_mean"] -
                        archived["cv_auroc_mean"]),
        "archived_pooled_cv_auroc": archived["pooled_cv_auroc"],
        "recomputed_pooled_cv_auroc": original["pooled_auroc"],
        "pooled_abs_diff": abs(original["pooled_auroc"] -
                               archived["pooled_cv_auroc"]),
        "match": bool(abs(original["auroc_mean"] -
                          archived["cv_auroc_mean"]) < 1e-9)}

    pred_scale = {}
    for cscale in S3_PRED_SCALE_FACTORS:
        pred_c = (b["ctrl_aligned"] +
                  cscale * (b["pred_expr"] - b["ctrl_aligned"]))
        feat_c = compute_features(pred_c.astype(np.float32), b["perts"],
                                  b["pred_genes"], b["ctrl_aligned"])
        feat_c["log_n_cells"] = b["feat"]["log_n_cells"].values
        X_c = feat_c[names].values.astype(np.float64)
        # 缩放不变性自查：neg_frac / delta_entropy 对正缩放应严格不变
        inv = np.allclose(b["feat"]["neg_frac"].values,
                          feat_c["neg_frac"].values, atol=1e-6) and \
            np.allclose(b["feat"]["delta_entropy"].values,
                        feat_c["delta_entropy"].values, atol=1e-6)
        # pearson 标签对正缩放数学不变：标签直接用原 y（留痕断言口径）
        r = probe_cv(X_c, y, in_fold_scale=True)  # 原协议重评
        pred_scale[f"delta_x{cscale}"] = {
            **r, "scale_invariant_features_unchanged": bool(inv),
            "label": "y unchanged (pearson_delta invariant to positive "
                     "rescaling of predicted delta)"}

    res = {
        "n": b["n"], "n_features": len(names), "feature_names": list(names),
        "original_protocol": {
            **original,
            "desc": "step6 run_probe: fold 内 StandardScaler + L2 logistic "
                    "(C=1.0), stratified 5-fold CV shuffle seed=0; headline "
                    "= CV AUROC mean"},
        "unified_scaler_control": {
            **unified,
            "desc": "StandardScaler 全数据一次 fit（统一缩放）后纯 LR "
                    "重训，同一 CV 划分重评",
            "cv_auroc_diff_vs_original": unified["auroc_mean"] -
            original["auroc_mean"],
            "pooled_auroc_diff_vs_original": unified["pooled_auroc"] -
            original["pooled_auroc"]},
        "no_scaling_diagnostic": {
            **raw,
            "desc": "诊断（非主对照）：原始量纲特征直接进 LR，无任何缩放"},
        "prediction_scaling_supplementary": pred_scale,
        "validation_vs_archived": validation,
    }
    d = res["unified_scaler_control"]["cv_auroc_diff_vs_original"]
    dp = res["unified_scaler_control"]["pooled_auroc_diff_vs_original"]
    print(f"[S3][{key}] n={b['n']} orig={original['auroc_mean']:.6f} "
          f"unified={unified['auroc_mean']:.6f} diff={d:+.6f} "
          f"(pooled diff={dp:+.6f}) raw={raw['auroc_mean']:.6f} "
          f"| archived match={validation['match']}")
    return res


def run_s3():
    combos = {}
    for key, data_name, out_dir in [
            ("gears_k562", "replogle_k562_essential", "outputs"),
            ("gears_rpe1", "replogle_rpe1_essential", "outputs_rpe1")]:
        combos[key] = run_s3_combo(key, data_name, out_dir)
    return {"spec": {
                "seed": S3_SEED,
                "primary_control": "probe input features StandardScaler "
                                   "(fit once on full data) -> retrain -> "
                                   "re-evaluate AUROC, vs step6 original "
                                   "in-fold-scaler pipeline (same folds)",
                "supplementary": "predicted delta scaled by constant factor "
                                 f"{S3_PRED_SCALE_FACTORS}, features "
                                 "recomputed, labels unchanged (pearson "
                                 "invariant), original protocol re-eval",
                "caveat": "全数据 fit 的 scaler 在口径上微跨测试折；该对照的"
                          "目的正是把缩放射到 CV 口径外观察 AUROC 是否移动"},
            "combos": combos}


def main():
    t_all = time.perf_counter()
    res = {"task": "redteam S4 paired test + S3 scaling control "
                   "(step13, real computation)",
           "redteam_ref": "redteam_20261024/gating.md S3/S4",
           "r2_product": "outputs_rebuild/two_layer_gate.json verdict "
                         "tier=symbol_layer",
           "date": time.strftime("%Y-%m-%d"),
           "s4": run_s4(),
           "s3": run_s3(),
           "timing_total_s": round(time.perf_counter() - t_all, 2)}
    os.makedirs(OUT_DIR, exist_ok=True)
    out_json = os.path.join(OUT_DIR, "s4_s3_fixes.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    print(f"[DONE] total {res['timing_total_s']:.1f}s")


if __name__ == "__main__":
    main()
