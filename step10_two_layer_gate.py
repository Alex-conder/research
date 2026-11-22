# step10_two_layer_gate.py — R2 两层门设计对决（prereg_rebuild_v2.md）
#
# 锁定口径（prereg R2）：
#   - 三臂：(a) 纯就地探针（5-fold OOF，step6 规格）；
#           (b) 纯符号层 = 同引擎另一数据集就地探针的幅度特征系数符号，
#               零训练，跨数据集直接携带；目标侧只用无监督标准化（不碰标签）；
#               分数 = 原始幅度排序（修正案 1 A1 明确不动）
#           (c) 两层门 = 符号层 + 就地校准（修正案 1 A1 裁定）：
#               第二层 = 就地探针拟合概率在符号组（符号层方向分组）内的
#               isotonic 重校准（5-fold OOF 防泄漏）；
#               臂 (c) 分数 = 符号层方向 x 重校准探针概率（连续）
#   - 评测：risk-coverage 轨迹 coverage 1.0 -> 0.5 步长 0.05（11 档），
#     带内 precision = 被标记查询中确为最差一半的占比（step6b 同口径），
#     risk = 1 - precision；配对 bootstrap（B=1000）比臂间差异
#   - 判定（锁定，两档制；修正案 1 A2 分母修正）：两层臂在 >=8/11 有效组合
#     优于符号臂（CI 不含 0）= 产品定为两层门；
#     否则 = 产品定为符号层、探针降为辅助
#     （锁定 >=8/12 按等严格度上取整，11 有效组合，更保守）
#
# 口径风险留痕（smoke 必查，写进 smoke_report）：
#   旧操作化（对符号分本身做 OOF isotonic）被 smoke 证为结构性退化：
#   isotonic 单调，对同符号分数做 OOF 校准后排序变化只来自跨折非单调性
#   （折划分噪声，smoke |delta|~0.002，CI 含 0）。修正案 1 A1 已裁定改用
#   "符号组内 isotonic 重校准探针概率 x 方向"；符号分与两层分的秩相关 /
#   最大秩位移诊断列保留，供复核第二层是否带来真实重排。
#
# 随机源（专属）：BOOT_MASTER_SEED = 20261026；OOF 划分种子 = 0（step6 同）。
#
# 产物：outputs_rebuild/two_layer_gate.{json,png}（全量）
#       outputs_rebuild/smoke/two_layer_gate.smoke.{json,png}（--smoke）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step10_two_layer_gate.py --smoke
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

for _s in ("stdout", "stderr"):  # GBK 控制台纪律
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")

BOOT_MASTER_SEED = 20261026
OOF_SEED = 0
N_FOLDS = 5

FULL_B = 1000
FULL_COV = [round(1.0 - 0.05 * k, 2) for k in range(11)]  # 1.0 -> 0.5
SMOKE_B = 100
SMOKE_COV = [1.0, 0.75, 0.5]

VERDICT_RULE = ("two-layer beats symbol in >=8/11 evaluable combos "
                "(amendment A2: locked >=8/12 at equal strictness; scf_k562 "
                "lacks sibling product so 11 evaluable combos) "
                "(paired bootstrap CI of coverage-averaged precision diff "
                "excludes 0) -> product = two-layer gate; "
                "else product = symbol layer")

from step10_signature_resolution import (  # 口径单一来源：加载器与组合表
    COMBOS, SMOKE_COMBOS, fit_coef, load_combo)

# 符号层符号来源：同引擎另一数据集（就地全量拟合取符号）
SIBLING = {
    "gears_k562": "gears_rpe1", "gears_rpe1": "gears_k562",
    "scgenept_k562": "scgenept_rpe1", "scgenept_rpe1": "scgenept_k562",
    "scf_k562": None,  # 无 scFoundation RPE1 产物，符号臂不可得
    "scgpt_k562": "scgpt_rpe1", "scgpt_rpe1": "scgpt_k562",
    "geneformer_k562": "geneformer_rpe1", "geneformer_rpe1": "geneformer_k562",
    "uce_k562": "uce_rpe1", "uce_rpe1": "uce_k562",
}


def precision_at(scores, y, cov):
    """top-cov 分数带内 y=1（最差一半）占比（与 step6b/6d 同口径）。"""
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def oof_probe_scores(X, y, seed=OOF_SEED):
    """同规格就地探针 5-fold OOF 概率（step6f oof_scores 同结构）。"""
    Xs = StandardScaler().fit_transform(X)
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    return cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                             Xs, y, cv=cv, method="predict_proba")[:, 1]


def oof_isotonic_scores(s, y, seed=OOF_SEED):
    """5-fold OOF isotonic 校准：训练折拟合、测试折预测（防泄漏）。"""
    out = np.empty_like(s, dtype=np.float64)
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    for tr, te in cv.split(s.reshape(-1, 1), y):
        iso = IsotonicRegression(out_of_bounds="clip").fit(s[tr], y[tr])
        out[te] = iso.predict(s[te])
    return out


def symbol_scores(X, signs):
    """纯符号层打分：目标侧无监督 z-score 后与符号向量点积（零训练）。
    臂 (b) = 原始幅度排序（修正案 1 A1：本函数不动）。"""
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return ((X - mu) / sd) @ signs


def two_layer_scores(s_sym, s_probe, y, seed=OOF_SEED):
    """修正案 1 A1：第二层 = 就地探针拟合概率在符号组内的 isotonic
    重校准（符号组 = 符号层方向组；组内 5-fold OOF 防泄漏）。
    臂 (c) 分数 = 符号层方向 x 重校准探针概率（连续）。
    小组 / 单类组无法做 5-fold 校准时退回原始探针概率（该组方向
    仍为 0 时分数归 0，不影响排序）。"""
    sign = np.sign(s_sym)
    cal = np.empty_like(s_probe, dtype=np.float64)
    for g in np.unique(sign):
        idx = np.where(sign == g)[0]
        yy = y[idx]
        if (g == 0 or len(idx) < 2 * N_FOLDS
                or np.unique(yy).size < 2
                or np.bincount(yy).min() < N_FOLDS):
            cal[idx] = s_probe[idx]  # 无可校准结构：退回原始探针概率
            continue
        cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
        for tr, te in cv.split(s_probe[idx].reshape(-1, 1), yy):
            iso = IsotonicRegression(out_of_bounds="clip").fit(
                s_probe[idx[tr]], yy[tr])
            cal[idx[te]] = iso.predict(s_probe[idx[te]])
    return sign * cal


def arm_scores(X, y, signs):
    """三臂点估计分数。signs=None 时只算就地臂（符号臂不可得）。"""
    s_a = oof_probe_scores(X, y)
    if signs is None:
        return {"in_situ": s_a}
    s_b = symbol_scores(X, signs)
    s_c = two_layer_scores(s_b, s_a, y)
    return {"in_situ": s_a, "symbol": s_b, "two_layer": s_c}


def prec_table(scores_dict, y, coverages):
    return {arm: {f"{c:.2f}": precision_at(s, y, c) for c in coverages}
            for arm, s in scores_dict.items()}


def paired_bootstrap(X, y, signs, coverages, b_boot, rng):
    """扰动级配对 bootstrap：每次重采样重算三臂分数与带内 precision，
    记录逐 coverage 的臂间配对差（键 = "{a}_minus_{b}" = a - b）。"""
    n = len(y)
    if signs is not None:
        pairs = [("two_layer", "symbol"), ("two_layer", "in_situ"),
                 ("symbol", "in_situ")]
    else:
        pairs = []
    diffs = {f"{a}_minus_{b}": np.zeros((b_boot, len(coverages)))
             for a, b in pairs}
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb = X[idx], y[idx]
        s_a = oof_probe_scores(Xb, yb)
        sc = {"in_situ": s_a}
        if signs is not None:
            sc["symbol"] = symbol_scores(Xb, signs)
            sc["two_layer"] = two_layer_scores(sc["symbol"], s_a, yb)
        for i, c in enumerate(coverages):
            pa = {arm: precision_at(s, yb, c) for arm, s in sc.items()}
            for a, bb in pairs:
                diffs[f"{a}_minus_{bb}"][b, i] = pa[a] - pa[bb]
    return diffs


def summarize_diffs(diffs, coverages):
    """逐 coverage 均值差 + 全 coverage 均值差的 95% percentile CI。"""
    out = {}
    for name, d in diffs.items():
        per_cov = d.mean(axis=0)
        overall = d.mean(axis=1)  # 配对：同一重采样跨 coverage 平均
        lo, hi = np.percentile(overall, [2.5, 97.5])
        out[name] = {
            "mean_diff_by_coverage": {f"{c:.2f}": float(per_cov[i])
                                      for i, c in enumerate(coverages)},
            "mean_over_coverage": float(overall.mean()),
            "ci_lo": float(lo), "ci_hi": float(hi),
            "ci_excludes_zero": bool(lo > 0 or hi < 0),
            "favors": ("first_arm" if lo > 0 else
                       "second_arm" if hi < 0 else "undetermined"),
        }
    return out


def max_rank_displacement(s1, s2):
    r1 = np.argsort(np.argsort(s1))
    r2 = np.argsort(np.argsort(s2))
    return int(np.max(np.abs(r1 - r2)))


def run_combo(key, coverages, b_boot, rng_master):
    c = COMBOS[key]
    t0 = time.perf_counter()
    X, y, feature_names = load_combo(key)
    n = len(y)
    sib_key = SIBLING.get(key)
    signs = None
    sib_note = None
    if sib_key is None:
        sib_note = "symbol arm unavailable: no sibling product for this combo"
        print(f"[R2][{key}] no sibling -> symbol/two_layer arms skipped")
    else:
        Xs, ys, _ = load_combo(sib_key)
        signs = np.sign(fit_coef(Xs, ys))
        sib_note = f"signs from {sib_key} in-situ fit (zero training here)"
    t_load = time.perf_counter() - t0

    t1 = time.perf_counter()
    sc = arm_scores(X, y, signs)
    t_point = time.perf_counter() - t1

    table = prec_table(sc, y, coverages)

    diag = {}
    if signs is not None:
        rho, _ = spearmanr(sc["symbol"], sc["two_layer"])
        diag["spearman_symbol_vs_two_layer"] = float(rho)
        diag["max_rank_displacement_symbol_vs_two_layer"] = \
            max_rank_displacement(sc["symbol"], sc["two_layer"])
        # 诊断（非锁定臂）：就地探针 OOF isotonic 校准，观察校准对探针臂的影响
        t2 = time.perf_counter()
        s_a_cal = oof_isotonic_scores(sc["in_situ"], y)
        t_diag = time.perf_counter() - t2
        diag["probe_oof_isotonic_precision"] = {
            f"{cov:.2f}": precision_at(s_a_cal, y, cov) for cov in coverages}
        diag["probe_oof_isotonic_time_s"] = round(t_diag, 3)

    rng = np.random.default_rng(rng_master.integers(0, 2**31))
    t3 = time.perf_counter()
    diffs = paired_bootstrap(X, y, signs, coverages, b_boot, rng)
    t_boot = time.perf_counter() - t3
    summ = summarize_diffs(diffs, coverages)

    res = {
        "engine": c["engine"], "dataset": c["dataset"], "class": c["cls"],
        "n": int(n), "n_features": len(feature_names),
        "feature_names": feature_names,
        "symbol_layer": sib_note, "missing": False,
        "arms_available": list(sc.keys()),
        "precision_table": table,
        "paired_bootstrap": {"b_boot": b_boot, **summ},
        "diagnostics": diag,
        "timing_s": {"load": round(t_load, 3), "point_arms": round(t_point, 3),
                     "bootstrap": round(t_boot, 3)},
    }
    if "two_layer_minus_symbol" in summ:
        d = summ["two_layer_minus_symbol"]
        print(f"[R2][{key}] n={n} two_layer-symbol mean={d['mean_over_coverage']:+.4f} "
              f"CI=[{d['ci_lo']:+.4f},{d['ci_hi']:+.4f}] excl0={d['ci_excludes_zero']} "
              f"| boot {t_boot:.1f}s")
    else:
        print(f"[R2][{key}] n={n} in_situ only | boot {t_boot:.1f}s")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="迷你规格（B=100, coverage={1.0,0.75,0.5}, 2 组合）")
    ap.add_argument("--combos", nargs="*", default=None)
    args = ap.parse_args()

    smoke = args.smoke
    b_boot = SMOKE_B if smoke else FULL_B
    coverages = SMOKE_COV if smoke else FULL_COV
    combo_keys = args.combos or (SMOKE_COMBOS if smoke else list(COMBOS))

    boot_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "R2 two-layer gate duel (step10, prereg_rebuild_v2.md)",
           "prereg": "prereg_rebuild_v2.md#R2",
           "smoke": smoke,
           "spec": {"b_boot": b_boot, "coverage_grid": coverages,
                    "n_folds": N_FOLDS, "oof_seed": OOF_SEED,
                    "metric": "precision in top-coverage flagged band "
                              "(step6b caliber), risk = 1 - precision",
                    "full_run_spec": {"b_boot": FULL_B,
                                      "coverage_grid": FULL_COV}},
           "boot_master_seed": BOOT_MASTER_SEED,
           "verdict_rule": VERDICT_RULE,
           "combos": {}, "missing": []}
    t_all = time.perf_counter()
    for key in combo_keys:
        c = COMBOS[key]
        if c["kind"] == "missing":
            res["missing"].append(key)
            res["combos"][key] = {"engine": c["engine"], "dataset": c["dataset"],
                                  "class": c["cls"], "missing": True}
            continue
        res["combos"][key] = run_combo(key, coverages, b_boot, boot_master)
    total = time.perf_counter() - t_all

    evaluated = [k for k, v in res["combos"].items()
                 if not v.get("missing") and "two_layer" in v["arms_available"]]
    n_two_better = sum(
        1 for k in evaluated
        if res["combos"][k]["paired_bootstrap"]["two_layer_minus_symbol"]
        ["ci_excludes_zero"]
        and res["combos"][k]["paired_bootstrap"]["two_layer_minus_symbol"]
        ["favors"] == "first_arm")
    if smoke:
        verdict = {"tier": "smoke_only_not_adjudicated",
                   "note": "smoke numbers do not enter the locked criteria",
                   "n_two_better_smoke_numbers": int(n_two_better),
                   "n_combos_with_symbol_arm": len(evaluated)}
    else:
        verdict = {"rule": VERDICT_RULE,
                   "n_two_better": int(n_two_better),
                   "n_combos_with_symbol_arm": len(evaluated),
                   "threshold": ">=8/11 (locked; amendment A2: >=8/12 at "
                                "equal strictness, rounded up on 11 "
                                "evaluable combos = 72.7%, conservative)",
                   "denominator_note": "amendment A2 (locked before full "
                                       "run): scf_k562 lacks sibling "
                                       "product -> 11 evaluable combos; "
                                       "scf_k562 in-situ arm only, "
                                       "excluded from arm comparison"}
        verdict["tier"] = ("two_layer_gate" if n_two_better >= 8
                           else "symbol_layer")
    res["verdict"] = verdict
    res["timing_total_s"] = round(total, 2)

    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    stem = "two_layer_gate.smoke" if smoke else "two_layer_gate"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    print(f"[DONE] total {total:.1f}s | verdict {verdict.get('tier')}")


def make_figure(res, path):
    keys = [k for k, v in res["combos"].items()
            if not v.get("missing") and "two_layer" in v.get("arms_available", [])]
    if not keys:
        return
    n_p = len(keys)
    fig, axes = plt.subplots(1, n_p, figsize=(3.6 * n_p, 4.0), sharey=True)
    if n_p == 1:
        axes = [axes]
    colors = {"in_situ": "#1f77b4", "symbol": "#d62728", "two_layer": "#2ca02c"}
    for ax, key in zip(axes, keys):
        c = res["combos"][key]
        for arm in ("in_situ", "symbol", "two_layer"):
            if arm not in c["precision_table"]:
                continue
            xs = sorted(float(k) for k in c["precision_table"][arm])
            ys = [c["precision_table"][arm][f"{x:.2f}"] for x in xs]
            ax.plot(xs, ys, "-o", ms=3.5, lw=1.5, color=colors[arm], label=arm)
        ax.axhline(0.5, color="gray", ls=":", lw=1.0)
        ax.set_xlabel("coverage (fraction flagged)", fontsize=8)
        ax.set_title(f"{c['engine']} / {c['dataset']} (n={c['n']})", fontsize=9)
        ax.set_xlim(0.45, 1.02)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, loc="lower left")
    axes[0].set_ylabel("precision (flagged & worst-half)", fontsize=8)
    fig.suptitle(f"R2 gate arms precision-coverage "
                 f"(B={res['spec']['b_boot']}, smoke={res['smoke']})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
