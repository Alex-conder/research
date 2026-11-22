# step10_signature_resolution.py — R1 签名解析度-样本量曲线（prereg_rebuild_v2.md）
#
# 锁定口径（prereg R1）：
#   - 对象：6 引擎 x 2 数据集就地探针（GEARS / scGenePT / scFoundation /
#     scGPT / Geneformer / UCE x K562 / RPE1）
#   - 方法：逐特征 bootstrap（B=500）探针系数 95% percentile CI；
#     "稳定签名" = CI 不含零的特征集合；子采样曲线 n in {15,30,60,120,全量}
#     x 3 子采样种子；读数"解析首个稳定签名所需扰动数"
#   - embedding 类口径（2026-10-16 锁定）：误差 = err_emb_dist（模型空间距离），
#     特征 = 嵌入几何（pred_emb_norm / delta_emb_norm + 设计协变量），
#     不混入表达谱口径；true_emb_norm 为真值侧量（泄漏纪律）不进特征
#   - 判定（锁定，两档制）：>=4/12 组合在全量 n 有 >=2 个稳定特征 = A 档；
#     <4/12 = B 档。smoke 数字不进判据。
#
# 口径先例复用（不自造加载逻辑）：
#   - GEARS: step6d_naive_baselines.load_dataset（seed0 预测，9 特征，pearson 标签）
#   - scGenePT / scFoundation: step6f_crossmodel.load_crossmodel
#   - 就地探针规格：step8b fit_coef（StandardScaler + LogisticRegression C=1.0）
#   - wave 六份产物消费方式：outputs_crossmodel/wave/READY.md
#
# 随机源（专属）：BOOT_MASTER_SEED = 20261024；子采样 SUBSET_MASTER_SEED = 20261025。
# 每组合由主 seed 派生独立 Generator，点序确定可复现。
#
# 产物：outputs_rebuild/signature_resolution.{json,png}（全量）
#       outputs_rebuild/smoke/signature_resolution.smoke.{json,png}（--smoke）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step10_signature_resolution.py --smoke
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

for _s in ("stdout", "stderr"):  # GBK 控制台纪律：不崩、不丢字符
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")

BOOT_MASTER_SEED = 20261024
SUBSET_MASTER_SEED = 20261025

# 锁定规格（prereg R1）与 smoke 迷你规格
FULL_B = 500
FULL_N_GRID = [15, 30, 60, 120, None]   # None = 全量
FULL_SUB_SEEDS = [0, 1, 2]
SMOKE_B = 50
SMOKE_N_GRID = [15, 60, None]
SMOKE_SUB_SEEDS = [0]

STABLE_CRIT = 2     # 判定用：全量 n 稳定特征 >= 2
VERDICT_RULE = ("A: >=4/12 combos with >=2 stable features at full n; "
                "else B (two-tier, no middle)")

EMB_FEATURES = ["pred_emb_norm", "delta_emb_norm", "log_n_cells",
                "target_in_input_space"]

# ---------------------------------------------------------------- 组合注册表
WAVE_DIR = os.path.join(REPO, "outputs_crossmodel", "wave")
SG_DIR = os.path.join(REPO, "outputs_crossmodel", "scgenept")
SCF_DIR = os.path.join(REPO, "outputs_crossmodel", "scfoundation")

COMBOS = {
    # 表达谱类：误差 = 1 - pearson_delta，特征 = 9 谱特征（step6 口径）
    "gears_k562":   {"engine": "GEARS", "dataset": "K562", "cls": "expression",
                     "kind": "gears", "data": "replogle_k562_essential",
                     "out": "outputs"},
    "gears_rpe1":   {"engine": "GEARS", "dataset": "RPE1", "cls": "expression",
                     "kind": "gears", "data": "replogle_rpe1_essential",
                     "out": "outputs_rpe1"},
    "scgenept_k562": {"engine": "scGenePT", "dataset": "K562", "cls": "expression",
                      "kind": "crossmodel",
                      "npz": os.path.join(SG_DIR, "scgenept_k562_predictions.npz"),
                      "labels": os.path.join(SG_DIR, "scgenept_k562_labels.csv")},
    "scgenept_rpe1": {"engine": "scGenePT", "dataset": "RPE1", "cls": "expression",
                      "kind": "crossmodel",
                      "npz": os.path.join(SG_DIR, "scgenept_rpe1_predictions.npz"),
                      "labels": os.path.join(SG_DIR, "scgenept_rpe1_labels.csv")},
    "scf_k562":     {"engine": "scFoundation", "dataset": "K562", "cls": "expression",
                     "kind": "crossmodel",
                     "npz": os.path.join(SCF_DIR, "scf_k562_predictions.npz"),
                     "labels": os.path.join(SCF_DIR, "scf_k562_labels.csv")},
    "scf_rpe1":     {"engine": "scFoundation", "dataset": "RPE1", "cls": "expression",
                     "kind": "missing"},
    # embedding 类：误差 = err_emb_dist，特征 = 嵌入几何（READY.md 口径）
    "scgpt_k562":   {"engine": "scGPT", "dataset": "K562", "cls": "embedding",
                     "kind": "wave", "wave": "scgpt"},
    "scgpt_rpe1":   {"engine": "scGPT", "dataset": "RPE1", "cls": "embedding",
                     "kind": "wave", "wave": "scgpt"},
    "geneformer_k562": {"engine": "Geneformer", "dataset": "K562", "cls": "embedding",
                        "kind": "wave", "wave": "geneformer"},
    "geneformer_rpe1": {"engine": "Geneformer", "dataset": "RPE1", "cls": "embedding",
                        "kind": "wave", "wave": "geneformer"},
    "uce_k562":     {"engine": "UCE", "dataset": "K562", "cls": "embedding",
                     "kind": "wave", "wave": "uce"},
    "uce_rpe1":     {"engine": "UCE", "dataset": "RPE1", "cls": "embedding",
                     "kind": "wave", "wave": "uce"},
}

SMOKE_COMBOS = ["gears_k562", "geneformer_k562"]


# ---------------------------------------------------------------- 加载（口径复用）
def load_gears(data_name, out_dir):
    from step6d_naive_baselines import load_dataset
    feat, feature_names, y, _pred_df, _ = load_dataset(data_name, out_dir)
    X = feat[feature_names].values.astype(np.float64)
    return X, y, list(feature_names)


def load_crossmodel(npz_path, labels_path):
    from step6f_crossmodel import load_crossmodel
    X, y, feature_names, _perts, _extra = load_crossmodel(npz_path, labels_path)
    return X, y, list(feature_names)


def load_wave(wave_name, ds_key):
    """wave 产物消费（READY.md）：embedding 类误差 = err_emb_dist，
    特征 = 嵌入几何 pred_emb_norm / delta_emb_norm + log_n_cells +
    target_in_input_space。true_emb_norm 为真值侧量，不进特征。"""
    import pandas as pd
    lab = pd.read_csv(os.path.join(
        WAVE_DIR, wave_name, f"{wave_name}_wave_{ds_key}_labels.csv"))
    err = lab["err_emb_dist"].values.astype(np.float64)
    assert np.isfinite(err).all(), "err_emb_dist 含非有限值"
    y = (err > np.median(err)).astype(int)
    in_space = lab["target_in_input_space"].astype(bool).astype(np.float64)
    X = np.column_stack([
        lab["pred_emb_norm"].values.astype(np.float64),
        lab["delta_emb_norm"].values.astype(np.float64),
        np.log1p(lab["n_cells"].values.astype(np.float64)),
        in_space,
    ])
    assert np.isfinite(X).all(), "特征矩阵含 NaN"
    return X, y, list(EMB_FEATURES)


def load_combo(key):
    c = COMBOS[key]
    if c["kind"] == "gears":
        return load_gears(c["data"], c["out"])
    if c["kind"] == "crossmodel":
        return load_crossmodel(c["npz"], c["labels"])
    if c["kind"] == "wave":
        return load_wave(c["wave"], "k562" if c["dataset"] == "K562" else "rpe1")
    raise FileNotFoundError(f"{key}: 产物缺失（注册为 missing）")


# ---------------------------------------------------------------- 探针与 bootstrap
def fit_coef(X, y):
    """同规格就地探针（in-sample）：StandardScaler + L2 logistic, C=1.0。"""
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xs, y)
    return clf.coef_[0]


def bootstrap_ci(X, y, rng, b_boot):
    """扰动级 bootstrap：有放回重采样 -> 重拟合 -> 逐系数 95% percentile CI。"""
    n = len(y)
    coefs = np.empty((b_boot, X.shape[1]))
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:  # 极端重采样单类，重抽替代
            idx = rng.integers(0, n, n)
        coefs[b] = fit_coef(X[idx], y[idx])
    lo = np.percentile(coefs, 2.5, axis=0)
    hi = np.percentile(coefs, 97.5, axis=0)
    return lo, hi


def stable_set(lo, hi):
    """稳定签名 = CI 不含零的特征集合。"""
    return [j for j in range(len(lo)) if not (lo[j] <= 0 <= hi[j])]


def first_n_reaching(curve, thresh):
    """curve: [(n, mean_stable_count)]；返回首个 mean >= thresh 的 n（升序）。"""
    for n, m in sorted(curve, key=lambda t: (t[0] is None, t[0])):
        if m >= thresh:
            return n
    return None


# ---------------------------------------------------------------- 主流程
def run_combo(key, b_boot, n_grid, sub_seeds, rng_master):
    c = COMBOS[key]
    t0 = time.perf_counter()
    X, y, feature_names = load_combo(key)
    n = len(y)
    t_load = time.perf_counter() - t0

    rng = np.random.default_rng(rng_master.integers(0, 2**31))

    # --- 全量 n：bootstrap CI（判定输入） ---
    t1 = time.perf_counter()
    lo, hi = bootstrap_ci(X, y, rng, b_boot)
    t_full = time.perf_counter() - t1
    st = stable_set(lo, hi)
    point = fit_coef(X, y)
    feats = {}
    for j, f in enumerate(feature_names):
        feats[f] = {"point": float(point[j]),
                    "ci_lo": float(lo[j]), "ci_hi": float(hi[j]),
                    "stable": j in st}

    # --- 子采样曲线：每 (n, seed) 点做子采样内 bootstrap ---
    t2 = time.perf_counter()
    grid = [k for k in n_grid if k is None or k <= n]
    per_point = []   # (n_label, seed, stable_count)
    curve_rows = []
    for k in grid:
        for seed in sub_seeds:
            kk = n if k is None else k
            if kk >= n:  # 全量点不重复抽，直接用全量结果
                cnt = len(st)
            else:
                rng_sub = np.random.default_rng(rng_master.integers(0, 2**31))
                sub = rng_sub.choice(n, size=kk, replace=False)
                while len(np.unique(y[sub])) < 2:
                    sub = rng_sub.choice(n, size=kk, replace=False)
                lo_s, hi_s = bootstrap_ci(X[sub], y[sub], rng, b_boot)
                cnt = len(stable_set(lo_s, hi_s))
            per_point.append((kk, seed, cnt))
    t_curve = time.perf_counter() - t2

    for k in grid:
        kk = n if k is None else k
        cnts = [c_ for (nn, _s, c_) in per_point if nn == kk]
        curve_rows.append({"n": kk, "is_full": k is None or kk >= n,
                           "stable_count_mean": float(np.mean(cnts)),
                           "stable_count_min": int(min(cnts)),
                           "stable_count_max": int(max(cnts)),
                           "n_seeds": len(cnts)})

    sub_curve = [(r["n"], r["stable_count_mean"]) for r in curve_rows
                 if not r["is_full"]]
    res = {
        "engine": c["engine"], "dataset": c["dataset"], "class": c["cls"],
        "n": int(n), "n_features": len(feature_names),
        "feature_names": feature_names,
        "missing": False,
        "full": {"b_boot": b_boot,
                 "n_stable_features": len(st),
                 "stable_features": [feature_names[j] for j in st],
                 "features": feats},
        "curve": curve_rows,
        "n_first_stable_ge1": first_n_reaching(sub_curve, 1),
        "n_first_stable_ge2": first_n_reaching(sub_curve, STABLE_CRIT),
        "timing_s": {"load": round(t_load, 3),
                     "full_bootstrap": round(t_full, 3),
                     "curve": round(t_curve, 3)},
    }
    print(f"[R1][{key}] n={n} full-stable={len(st)}/{len(feature_names)} "
          f"curve=" + ",".join(f"{r['n']}:{r['stable_count_mean']:.1f}"
                              for r in curve_rows) +
          f" | full_boot {t_full:.1f}s curve {t_curve:.1f}s")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="迷你规格（B=50, n={15,60,full}, 单种子, 2 组合）")
    ap.add_argument("--combos", nargs="*", default=None, help="覆盖组合清单")
    args = ap.parse_args()

    smoke = args.smoke
    b_boot = SMOKE_B if smoke else FULL_B
    n_grid = SMOKE_N_GRID if smoke else FULL_N_GRID
    sub_seeds = SMOKE_SUB_SEEDS if smoke else FULL_SUB_SEEDS
    combo_keys = args.combos or (SMOKE_COMBOS if smoke else list(COMBOS))

    boot_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "R1 signature resolution vs sample size (step10, "
                   "prereg_rebuild_v2.md)",
           "prereg": "prereg_rebuild_v2.md#R1",
           "smoke": smoke,
           "spec": {"b_boot": b_boot, "n_grid": n_grid,
                    "sub_seeds": sub_seeds,
                    "stable_feature_def": "95% bootstrap CI excludes zero",
                    "full_run_spec": {"b_boot": FULL_B, "n_grid": FULL_N_GRID,
                                      "sub_seeds": FULL_SUB_SEEDS}},
           "boot_master_seed": BOOT_MASTER_SEED,
           "subset_master_seed": SUBSET_MASTER_SEED,
           "label": {"expression": "y = (1 - pearson_delta) > median "
                                   "(step6 definition)",
                     "embedding": "y = err_emb_dist > median (2026-10-16 "
                                  "locked embedding-error caliber)"},
           "features_note": "embedding class: pred_emb_norm / delta_emb_norm "
                            "+ log_n_cells + target_in_input_space; "
                            "true_emb_norm excluded (truth-derived, "
                            "leakage discipline)",
           "verdict_rule": VERDICT_RULE,
           "combos": {}, "missing": []}
    t_all = time.perf_counter()
    for key in combo_keys:
        c = COMBOS[key]
        if c["kind"] == "missing":
            res["missing"].append(key)
            res["combos"][key] = {"engine": c["engine"], "dataset": c["dataset"],
                                  "class": c["cls"], "missing": True}
            print(f"[R1][{key}] MISSING -> recorded")
            continue
        res["combos"][key] = run_combo(key, b_boot, n_grid, sub_seeds,
                                       boot_master)
    total = time.perf_counter() - t_all

    # --- 判定（仅全量模式且非 smoke 时生效） ---
    evaluated = [k for k, v in res["combos"].items() if not v.get("missing")]
    n_stable2 = sum(1 for k in evaluated
                    if res["combos"][k]["full"]["n_stable_features"] >= STABLE_CRIT)
    if smoke:
        verdict = {"tier": "smoke_only_not_adjudicated",
                   "note": "smoke numbers do not enter the locked criteria",
                   "n_combos_stable_ge2_smoke_numbers": int(n_stable2),
                   "n_combos_evaluated": len(evaluated)}
    else:
        tier = "A" if n_stable2 >= 4 else "B"
        verdict = {"tier": tier, "rule": VERDICT_RULE,
                   "n_combos_stable_ge2": int(n_stable2),
                   "n_combos_evaluated": len(evaluated),
                   "denominator_note": "locked rule says 4/12; evaluated "
                                       "combos here exclude missing products"}
    res["verdict"] = verdict
    res["timing_total_s"] = round(total, 2)

    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    stem = "signature_resolution.smoke" if smoke else "signature_resolution"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    print(f"[DONE] total {total:.1f}s | verdict {verdict.get('tier')}")


def make_figure(res, path):
    keys = [k for k, v in res["combos"].items() if not v.get("missing")]
    if not keys:
        return
    n_p = len(keys)
    fig, axes = plt.subplots(1, n_p, figsize=(3.2 * n_p, 3.8), sharey=True)
    if n_p == 1:
        axes = [axes]
    for ax, key in zip(axes, keys):
        c = res["combos"][key]
        rows = [r for r in c["curve"] if not r["is_full"]]
        xs = [r["n"] for r in rows]
        ys = [r["stable_count_mean"] for r in rows]
        ax.plot(xs, ys, "-o", color="#1f77b4", ms=4.5, lw=1.6)
        full_y = c["full"]["n_stable_features"]
        ax.plot([c["n"]], [full_y], "s", color="#d62728", ms=7,
                label=f"full n={c['n']}")
        ax.axhline(STABLE_CRIT, color="gray", ls="--", lw=1.0)
        ax.text(c["n"], STABLE_CRIT + 0.15, f"crit >= {STABLE_CRIT}",
                ha="right", fontsize=7, color="gray")
        ax.set_xscale("log")
        xt = sorted(xs + [c["n"]])
        ax.set_xticks(xt)
        ax.set_xticklabels(xt, fontsize=7)
        ax.set_xlabel("perturbations in subset (log)", fontsize=8)
        ax.set_title(f"{c['engine']} / {c['dataset']}\n(n={c['n']}, "
                     f"full-stable={full_y}/{c['n_features']})", fontsize=9)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=7, loc="upper left")
    axes[0].set_ylabel("mean # stable features (CI excludes 0)", fontsize=8)
    fig.suptitle(f"R1 signature resolution curve "
                 f"(B={res['spec']['b_boot']}, smoke={res['smoke']})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
