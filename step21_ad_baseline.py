# step21_ad_baseline.py -- AD (additive disagreement) baseline + dual-metric audit
#
# 预注册：prereg_ad_baseline.md（判据先于计算锁定，两档制）。
#
# 分析 A（次级证据，偏离口径 AD-mean）：5 数据集单基因留出 GEARS 产物上，
#   AD-mean = mean_g |pred_delta - mean_train_delta|（论文声明 AD 对单基因
#   OOD 无定义，见 prereg 偏离 D-AD1），与就地探针（step6 规格，三口径
#   y 重拟合，OOF 评估）与符号门（two_layer_gate 同构造，兄弟符号源锁定）
#   并排报三口径（full / top20 / top50）失败检测 AUROC。
# 分析 B（主判定场，论文忠实口径）：加载归档 gears_model_seed4（CPU、零训练），
#   对 norman 含 >=1 留出基因的组合条件做 on-file 推断；AD-combo =
#   mean_g |pred_delta - (delta_A + delta_B)|（观测单基因效应加和，prereg
#   偏离 D-AD2 披露）；探针 = 留出单基因上每口径全量拟合后外推打分；
#   符号门 = adamson(s3+s5) 符号零训练。
#
# 统计：AUROC（秩检验，平局平均秩）；配对差 AD-符号门 / AD-探针的扰动级
#   bootstrap（固定分数不重训），B=2000，95% percentile CI，seed=20261030。
#
# 产物：outputs_rebuild/ad_baseline.{json,png} + ad_baseline_verdict.md
#       + ad_baseline_norman_combo_pred.npz（组合预测留档）
#       smoke 模式 -> outputs_rebuild/smoke/ad_baseline.smoke.*
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step21_ad_baseline.py --smoke
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

for _s in ("stdout", "stderr"):  # GBK 控制台纪律
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

REPO = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(REPO, "data", "gears_demo")
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")

BOOT_MASTER_SEED = 20261030
B_BOOT_FULL = 2000
B_BOOT_SMOKE = 200
OOF_SEED = 0
METRICS = ["full", "top20", "top50"]
HOLDOUT_FRAC = 0.2
SPLIT_SEED = 0

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features  # noqa: E402

# ------------------------------------------------------------- 组合注册（分析 A）
A_COMBOS = [
    {"key": "gears_k562_s0", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 0},
    {"key": "gears_k562_s1", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 1},
    {"key": "gears_k562_s2", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 2},
    {"key": "gears_rpe1_s0", "ds": "replogle_rpe1_essential", "out": "outputs_rpe1",
     "seed": 0},
    {"key": "gears_norman_s4", "ds": "norman",
     "out": "outputs_multidataset/norman", "seed": 4},
    {"key": "gears_adamson_s3", "ds": "adamson",
     "out": "outputs_multidataset/adamson", "seed": 3},
    {"key": "gears_adamson_s5", "ds": "adamson",
     "out": "outputs_multidataset/adamson", "seed": 5},
    {"key": "gears_dixit_s4", "ds": "dixit",
     "out": "outputs_multidataset/dixit", "seed": 4},
    {"key": "gears_dixit_s5", "ds": "dixit",
     "out": "outputs_multidataset/dixit", "seed": 5},
]

# 符号源（锁定，prereg）：兄弟数据集就地全量拟合取符号（标准标签）
SIGN_SOURCES = {
    "replogle_k562_essential": [("replogle_rpe1_essential", "outputs_rpe1", 0)],
    "replogle_rpe1_essential": [("replogle_k562_essential", "outputs", 0)],
    "norman": [("adamson", "outputs_multidataset/adamson", 3),
               ("adamson", "outputs_multidataset/adamson", 5)],
    "adamson": [("norman", "outputs_multidataset/norman", 4)],
    "dixit": [("norman", "outputs_multidataset/norman", 4)],
}

_BUNDLES = {}


# ------------------------------------------------------------- 数据加载
def load_bundle(ds):
    """数据集级缓存：ctrl 均值、条件均值、留出基因重放、单基因效应。"""
    if ds in _BUNDLES:
        return _BUNDLES[ds]
    import anndata as ad
    t0 = time.time()
    h5 = os.path.join(DATA_DIR, ds, "perturb_processed.h5ad")
    adata = ad.read_h5ad(h5)
    var_genes = adata.var["gene_name"].astype(str).tolist()
    conds = adata.obs["condition"].astype(str).to_numpy()
    ctrl_mean = np.asarray(adata.X[conds == "ctrl"].mean(axis=0)).ravel()
    cond_means, cond_cells = {}, {}
    for c in np.unique(conds):
        m = conds == c
        cond_means[c] = np.asarray(adata.X[m].mean(axis=0)).ravel()
        cond_cells[c] = int(m.sum())
    adata = adata[:, :0]  # 释放表达矩阵
    del adata

    singles = [c for c in cond_means if ("ctrl" in c and c != "ctrl")]
    genes = sorted({c.replace("ctrl+", "").replace("+ctrl", "") for c in singles})
    genes = [g for g in genes if "_" not in g]
    rng = np.random.default_rng(SPLIT_SEED)
    holdout = set(rng.choice(genes, size=int(HOLDOUT_FRAC * len(genes)),
                             replace=False).tolist())
    single_cands = {}
    for g in genes:
        cand = [c for c in singles if c.startswith(g + "+") or
                c.endswith("+" + g)]
        if cand:
            single_cands[g] = sorted(cand)
    single_cond = {g: cands[0] for g, cands in single_cands.items()}
    single_delta = {g: (cond_means[c] - ctrl_mean)
                    for g, c in single_cond.items()}
    train_deltas = [d for g, d in single_delta.items() if g not in holdout]
    train_delta_mean = np.mean(train_deltas, axis=0)
    combos_all = [c for c in cond_means if ("ctrl" not in c and "+" in c)]
    b = {"var_genes": var_genes, "ctrl_mean": ctrl_mean,
         "cond_means": cond_means, "cond_cells": cond_cells,
         "holdout": holdout, "single_cond": single_cond,
         "single_cands": single_cands,
         "single_delta": single_delta, "train_delta_mean": train_delta_mean,
         "combos_all": sorted(combos_all)}
    _BUNDLES[ds] = b
    print(f"[load] {ds}: genes={len(var_genes)} conds={len(cond_means)} "
          f"holdout={len(holdout)} singles={len(single_cond)} "
          f"combos={len(combos_all)} ({time.time()-t0:.0f}s)", flush=True)
    return b


def align(ref_genes, values, pred_genes):
    if list(ref_genes) == list(pred_genes):
        return values
    pos = {g: i for i, g in enumerate(ref_genes)}
    return np.array([values[pos[g]] for g in pred_genes])


def load_product(ds, out_dir, seed, bundle):
    """产物级：预测谱、特征、真值 delta、标准标签。返回 dict 或 None。"""
    out_path = os.path.join(REPO, out_dir)
    z = np.load(os.path.join(out_path, f"pred_expr_seed{seed}.npz"))
    pred_expr = np.asarray(z["pred_expr"], dtype=np.float64)
    perts = [str(p) for p in z["pred_perturbations"].tolist()]
    pred_genes = [str(g) for g in z["pred_genes"].tolist()]
    pred_df = pd.read_csv(
        os.path.join(out_path, f"predictions_seed{seed}.csv")
    ).set_index("perturbation").loc[perts]

    ctrl = align(bundle["var_genes"], bundle["ctrl_mean"], pred_genes)
    # 每个留出基因可能有多条单基因条件名（如 'g+ctrl'/'ctrl+g'）；
    # 用产物 csv 的 n_cells 与 h5ad 条件计数对账选定真值条件（同 step1_run 语义）
    use_conds = []
    for p in perts:
        n_csv = int(pred_df.loc[p, "n_cells"])
        cands = bundle["single_cands"][p]
        match = [c for c in cands if bundle["cond_cells"][c] == n_csv]
        use_conds.append(match[0] if match else cands[0])
    true = np.stack([align(bundle["var_genes"], bundle["cond_means"][c],
                           pred_genes) for c in use_conds])
    n_cells = pred_df["n_cells"].values.astype(np.float64)
    log_n_cells = np.log1p(n_cells)

    feat = compute_features(pred_expr.astype(np.float32), perts, pred_genes,
                            ctrl.astype(np.float32))
    feat["log_n_cells"] = log_n_cells
    excluded = []
    if pred_df["in_pert_graph"].nunique() > 1:
        feat["in_pert_graph"] = pred_df["in_pert_graph"].astype(int).values
    else:
        excluded.append("in_pert_graph (constant)")
    feature_names = SPECTRUM_FEATURES + [c for c in COVARIATES
                                         if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    assert np.isfinite(X).all() and np.isfinite(pred_expr).all()

    err_std = 1.0 - pred_df["pearson_delta"].values
    y_std = (err_std > np.median(err_std)).astype(int)
    return {"perts": perts, "pred_genes": pred_genes, "pred_expr": pred_expr,
            "ctrl": ctrl, "true": true, "X": X, "feature_names": feature_names,
            "y_std": y_std, "excluded": excluded,
            "pred_delta": pred_expr - ctrl, "true_delta": true - ctrl}


# ------------------------------------------------------------- 误差与信号
def error_metrics(pred_delta, true_delta):
    absd = np.abs(pred_delta - true_delta)
    e_full = absd.mean(axis=1)
    order = np.argsort(-np.abs(true_delta), axis=1)
    e20 = np.take_along_axis(absd, order[:, :20], axis=1).mean(axis=1)
    e50 = np.take_along_axis(absd, order[:, :50], axis=1).mean(axis=1)
    return {"full": e_full, "top20": e20, "top50": e50}


def y_from_e(e):
    return (e > np.median(e)).astype(int)


def auroc_fast(y, s):
    """Mann-Whitney 秩 AUROC（平局平均秩）。调用方保证两类各 >=2。"""
    r = rankdata(s)
    n1 = int(y.sum())
    n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def oof_probe_scores(X, y):
    """step6 规格 OOF 决策分；少数类 <2 返回 None；5<n<10 用 LOO。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import (LeaveOneOut, StratifiedKFold,
                                         cross_val_predict)
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    if np.bincount(y).min() < 2:
        return None
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=OOF_SEED))
    if len(y) < 10:
        loo = LeaveOneOut()
        return cross_val_predict(clf, X, y, cv=loo, method="decision_function")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=OOF_SEED)
    return cross_val_predict(clf, X, y, cv=cv, method="decision_function")


def full_fit_scores(X_fit, y_fit, X_score):
    """全量拟合（部署外推用，分析 B 探针）。少数类 <2 返回 None。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    if np.bincount(np.asarray(y_fit)).min() < 2:
        return None
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=OOF_SEED))
    clf.fit(X_fit, y_fit)
    return clf.decision_function(X_score)


def fit_signs(X, y):
    """two_layer_gate 同构造：就地全量拟合取符号。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xs, y)
    return np.sign(clf.coef_[0])


def symbol_scores(X, signs):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return ((X - mu) / sd) @ signs


def ad_scores(pred_delta, ref_delta):
    return np.abs(pred_delta - ref_delta).mean(axis=1)


def paired_boot(y, s_ad, s_gate, rng, b_boot):
    """扰动级配对 bootstrap（固定分数）：AUROC_AD - AUROC_gate 的 95% CI。"""
    n = len(y)
    diffs = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        while np.bincount(y[idx]).min() < 2:
            idx = rng.integers(0, n, n)
        diffs[b] = auroc_fast(y[idx], s_ad[idx]) - \
            auroc_fast(y[idx], s_gate[idx])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return {"mean_diff": float(diffs.mean()), "ci_lo": float(lo),
            "ci_hi": float(hi),
            "ci_excludes_zero": bool(lo > 0 or hi < 0),
            "favors": ("AD" if lo > 0 else "gate" if hi < 0 else "undetermined")}


def eval_signals(y, scores, rng, b_boot):
    """scores: dict name->score（None 跳过）。返回 AUROC 表 + 配对 CI。"""
    out_auroc, out_boot = {}, {}
    ok = {k: s for k, s in scores.items() if s is not None}
    if np.bincount(y).min() < 2:
        return None, None
    for k, s in ok.items():
        out_auroc[k] = auroc_fast(y, s)
    if "ad" in ok:
        for k in ("symbol", "probe"):
            if k in ok:
                out_boot[f"ad_minus_{k}"] = paired_boot(y, ok["ad"], ok[k],
                                                       rng, b_boot)
    return out_auroc, out_boot


# ------------------------------------------------------------- 分析 A
def run_analysis_a(combos, rng, b_boot):
    signs_cache = {}
    res = {}
    for c in combos:
        key = c["key"]
        t0 = time.time()
        bundle = load_bundle(c["ds"])
        prod = load_product(c["ds"], c["out"], c["seed"], bundle)
        perts = prod["perts"]
        n = len(perts)

        # 符号源（按数据集缓存）
        if c["ds"] not in signs_cache:
            srcs = SIGN_SOURCES[c["ds"]]
            Xs, ys = [], []
            for s_ds, s_out, s_seed in srcs:
                sb = load_bundle(s_ds)
                sp = load_product(s_ds, s_out, s_seed, sb)
                Xs.append(sp["X"])
                ys.append(sp["y_std"])
            signs_cache[c["ds"]] = fit_signs(np.vstack(Xs), np.concatenate(ys))
        signs = signs_cache[c["ds"]]

        ad_ref = align(bundle["var_genes"], bundle["train_delta_mean"],
                       prod["pred_genes"])
        s_ad = ad_scores(prod["pred_delta"], ad_ref)
        s_sym = symbol_scores(prod["X"], signs)
        errs = error_metrics(prod["pred_delta"], prod["true_delta"])

        entry = {"dataset": c["ds"], "seed": c["seed"], "n": int(n),
                 "features": prod["feature_names"],
                 "excluded_features": prod["excluded"],
                 "symbol_source": [f"{d}/{o}/s{s}"
                                   for d, o, s in SIGN_SOURCES[c["ds"]]],
                 "signals_note": "ad = AD-mean (protocol deviation D-AD1)",
                 "per_metric": {}}
        for m in METRICS:
            y = y_from_e(errs[m])
            s_probe = oof_probe_scores(prod["X"], y)
            aurocs, boots = eval_signals(
                y, {"ad": s_ad, "symbol": s_sym, "probe": s_probe}, rng, b_boot)
            if aurocs is None:
                entry["per_metric"][m] = {
                    "evaluable": False,
                    "reason": "minority class < 2 (amendment-1 rule)"}
                continue
            entry["per_metric"][m] = {
                "evaluable": True, "n_positive": int(y.sum()),
                "auroc": aurocs, "paired_bootstrap": boots,
                "probe_cv": ("5-fold" if n >= 10 else "loo")}
        res[key] = entry
        a_str = {m: (entry["per_metric"][m].get("auroc")
                     if entry["per_metric"][m].get("evaluable") else None)
                 for m in METRICS}
        print(f"[A][{key}] n={n} " +
              " ".join(f"{m}={a_str[m]}" for m in METRICS) +
              f" ({time.time()-t0:.0f}s)", flush=True)
    return res


def summarize_a(res_a):
    """数据集级（种子均值）汇总表。"""
    groups = {}
    for key, e in res_a.items():
        ds = e["dataset"]
        groups.setdefault(ds, []).append(e)
    out = {}
    for ds, entries in groups.items():
        row = {"seeds": [e["seed"] for e in entries],
               "n_per_seed": [e["n"] for e in entries], "per_metric": {}}
        for m in METRICS:
            vals = {k: [] for k in ("ad", "symbol", "probe")}
            evaluable = all(e["per_metric"][m].get("evaluable")
                            for e in entries)
            if evaluable:
                for e in entries:
                    for k in vals:
                        vals[k].append(e["per_metric"][m]["auroc"][k])
                row["per_metric"][m] = {
                    "evaluable": True,
                    "auroc_seed_mean": {k: float(np.mean(v))
                                        for k, v in vals.items()},
                    "auroc_seed_std": {k: (float(np.std(v, ddof=1))
                                           if len(v) > 1 else 0.0)
                                       for k, v in vals.items()}}
            else:
                row["per_metric"][m] = {"evaluable": False}
        out[ds] = row
    return out


# ------------------------------------------------------------- 分析 B
def run_analysis_b(rng, b_boot, combo_cap=None):
    import torch
    from gears import GEARS, PertData

    ds = "norman"
    bundle = load_bundle(ds)
    holdout = bundle["holdout"]
    prod = load_product(ds, "outputs_multidataset/norman", 4, bundle)

    # 候选组合：>=1 留出基因
    cand = []
    for c in bundle["combos_all"]:
        parts = sorted(c.split("+"))
        if len(parts) == 2 and any(g in holdout for g in parts):
            cand.append((c, parts))
    seen, combos_sel = set(), []
    for c, parts in cand:
        t = tuple(parts)
        if t not in seen:
            seen.add(t)
            combos_sel.append((c, parts))
    if combo_cap:
        combos_sel = combos_sel[:combo_cap]

    # 资格：两成分均有观测单基因条件（论 文资格）
    eligible, excl = [], {"no_single_condition": 0}
    for c, parts in combos_sel:
        if not all(p in bundle["single_delta"] for p in parts):
            excl["no_single_condition"] += 1
            continue
        eligible.append((c, parts))

    t0 = time.time()
    pert_data = PertData(DATA_DIR)
    pert_data.load(data_name=ds)
    pert_data.prepare_split(split="simulation", seed=SPLIT_SEED,
                            test_perts="_".join(sorted(holdout)),
                            only_test_set_perts=True)
    pert_data.get_dataloader(batch_size=32, test_batch_size=128)
    model = GEARS(pert_data, device="cpu")
    model.model_initialize(hidden_size=64)
    sd = torch.load(os.path.join(REPO, "outputs_multidataset/norman",
                                 "gears_model_seed4/model.pt"),
                    map_location="cpu", weights_only=True)
    model.model.load_state_dict(sd)
    model.best_model.load_state_dict(sd)
    model.saved_pred = {}
    model.saved_logvar_sum = {}
    print(f"[B] model ready ({time.time()-t0:.0f}s)", flush=True)

    queries = [list(parts) for _, parts in eligible]
    pred_map, ok_pairs, fail = {}, [], 0

    def lookup(out, parts):
        want = set(parts)
        return next((kk for kk in out if set(kk.split("_")) == want), None)

    try:
        out = model.predict(queries)
        for (c, parts) in eligible:
            k = lookup(out, parts)
            if k is None:
                fail += 1
                continue
            pred_map[c] = np.asarray(out[k]).ravel().astype(np.float64)
            ok_pairs.append((c, parts))
    except Exception as exc:  # 逐个降级
        print(f"[B] batch predict failed ({exc!r}); falling back", flush=True)
        for (c, parts) in eligible:
            try:
                out = model.predict([list(parts)])
                k = lookup(out, parts)
                if k is None:
                    fail += 1
                    continue
                pred_map[c] = np.asarray(out[k]).ravel().astype(np.float64)
                ok_pairs.append((c, parts))
            except Exception:
                fail += 1
    excl["predict_failed"] = fail
    print(f"[B] combos: candidates={len(combos_sel)} eligible={len(eligible)} "
          f"predicted={len(ok_pairs)} excluded={excl} ({time.time()-t0:.0f}s)",
          flush=True)

    pred_genes = prod["pred_genes"]
    ctrl = prod["ctrl"]
    perts_c = [c for c, _ in ok_pairs]
    pred_expr_c = np.stack([align(bundle["var_genes"], pred_map[c], pred_genes)
                            for c in perts_c])
    true_c = np.stack([align(bundle["var_genes"], bundle["cond_means"][c],
                             pred_genes) for c in perts_c])
    np.savez_compressed(
        os.path.join(OUT_DIR, "ad_baseline_norman_combo_pred.npz"),
        pert=np.array(perts_c), pred_expr=pred_expr_c.astype(np.float32),
        pred_genes=np.array(pred_genes))

    # 组合特征（与单基因产物同特征列）
    feat = compute_features(pred_expr_c.astype(np.float32), perts_c,
                            pred_genes, ctrl.astype(np.float32))
    feat["log_n_cells"] = np.log1p(
        np.array([bundle["cond_cells"][c] for c in perts_c], dtype=np.float64))
    feature_names = list(prod["feature_names"])
    missing = [f for f in feature_names if f not in feat.columns]
    assert not missing, f"combo features missing {missing}"
    X_c = feat[feature_names].values.astype(np.float64)

    # 加和参考 delta_A + delta_B（观测单基因效应；偏离 D-AD2 披露）
    add_ref = np.stack([
        align(bundle["var_genes"],
              bundle["single_delta"][p1] + bundle["single_delta"][p2],
              pred_genes) for _, (p1, p2) in ok_pairs])
    s_ad = ad_scores(pred_expr_c - ctrl, add_ref)

    # 符号源：adamson(s3+s5)
    Xs, ys = [], []
    for s_ds, s_out, s_seed in SIGN_SOURCES["norman"]:
        sb = load_bundle(s_ds)
        sp = load_product(s_ds, s_out, s_seed, sb)
        Xs.append(sp["X"])
        ys.append(sp["y_std"])
    signs = fit_signs(np.vstack(Xs), np.concatenate(ys))
    s_sym = symbol_scores(X_c, signs)

    errs = error_metrics(pred_expr_c - ctrl, true_c - ctrl)
    per_metric = {}
    for m in METRICS:
        y = y_from_e(errs[m])
        s_probe = full_fit_scores(prod["X"], y_from_e(
            error_metrics(prod["pred_delta"], prod["true_delta"])[m]),
            X_c)
        aurocs, boots = eval_signals(
            y, {"ad": s_ad, "symbol": s_sym, "probe": s_probe}, rng, b_boot)
        if aurocs is None:
            per_metric[m] = {"evaluable": False,
                             "reason": "minority class < 2"}
            continue
        per_metric[m] = {"evaluable": True, "n_positive": int(y.sum()),
                         "auroc": aurocs, "paired_bootstrap": boots,
                         "probe_mode": "full-fit on held-out singles, "
                                       "deployed on combos (OOD)"}
    return {"model": "gears norman seed4 (archived weights, CPU, no training)",
            "n_candidates": len(combos_sel), "n_eligible": len(eligible),
            "n_predicted": len(ok_pairs), "excluded": excl,
            "n_perts": len(perts_c), "per_metric": per_metric,
            "signals_note": "ad = AD-combo (paper-faithful; deviation D-AD2)",
            "feature_names": feature_names}


# ------------------------------------------------------------- 判定（锁定）
def adjudicate(b):
    rule = ("CONFIRMED iff (a) full-metric paired bootstrap CI of "
            "(AD-combo - symbol) excludes 0 favoring AD, AND (b) top20-metric "
            "CI includes 0 or favors symbol. Everything else = NOT_CONFIRMED. "
            "Secondary: same for AD-combo - probe.")
    pm = b["per_metric"]
    if not all(pm[m].get("evaluable") for m in ("full", "top20")):
        return {"tier": "NOT_EVALUABLE", "rule": rule,
                "note": "analysis B not evaluable on full/top20"}
    d_full = pm["full"]["paired_bootstrap"]["ad_minus_symbol"]
    d_top = pm["top20"]["paired_bootstrap"]["ad_minus_symbol"]
    cond_a = d_full["ci_excludes_zero"] and d_full["favors"] == "AD"
    cond_b = (not d_top["ci_excludes_zero"]) or d_top["favors"] == "gate"
    tier = "CONFIRMED" if (cond_a and cond_b) else "NOT_CONFIRMED"
    sec = {}
    for pair in ("ad_minus_symbol", "ad_minus_probe"):
        if pair in pm["full"]["paired_bootstrap"]:
            sec[pair] = {
                "full": pm["full"]["paired_bootstrap"][pair],
                "top20": pm["top20"]["paired_bootstrap"].get(pair),
                "top50": pm["top50"]["paired_bootstrap"].get(pair)}
    return {"tier": tier, "rule": rule,
            "primary": {"ad_minus_symbol_full": d_full,
                        "ad_minus_symbol_top20": d_top,
                        "cond_a_full_favors_AD": bool(cond_a),
                        "cond_b_top20_reverses_or_nsig": bool(cond_b)},
            "all_pairs_by_metric": sec}


# ------------------------------------------------------------- 图与报告
def make_figure(res, path):
    panels = []
    ds_label = {"replogle_k562_essential": "K562 (GEARS, 3 seeds)",
                "replogle_rpe1_essential": "RPE1 (GEARS)",
                "norman": "Norman singles (GEARS s4)",
                "adamson": "Adamson (GEARS, 2 seeds)",
                "dixit": "Dixit (GEARS)"}
    seen_ds = []
    for key, e in res["analysis_A"]["combos"].items():
        ds = e["dataset"]
        if ds in seen_ds:
            continue
        seen_ds.append(ds)
        row = res["analysis_A"]["dataset_summary"].get(ds)
        panels.append((ds_label[ds], row))
    if res["analysis_B"].get("per_metric"):
        panels.append(("Norman combos (AD home setting)",
                       res["analysis_B"]))
    panels = [(t, r) for t, r in panels if r and any(
        r["per_metric"][m].get("evaluable") for m in METRICS
        if m in r["per_metric"])]
    if not panels:
        return
    n_p = len(panels)
    fig, axes = plt.subplots(1, n_p, figsize=(3.4 * n_p, 3.8), sharey=True)
    if n_p == 1:
        axes = [axes]
    colors = {"ad": "#d62728", "symbol": "#1f77b4", "probe": "#2ca02c"}
    names = {"ad": "AD", "symbol": "symbol gate", "probe": "in-situ probe"}
    for ax, (title, row) in zip(axes, panels):
        x = np.arange(len(METRICS))
        w = 0.26
        for i, sig in enumerate(("ad", "symbol", "probe")):
            vals, marks = [], []
            for m in METRICS:
                cell = row["per_metric"][m]
                if not cell.get("evaluable") or sig not in cell.get(
                        "auroc_seed_mean", cell.get("auroc", {})):
                    vals.append(np.nan)
                    marks.append("n/a")
                    continue
                if "auroc_seed_mean" in cell:
                    vals.append(cell["auroc_seed_mean"][sig])
                else:
                    vals.append(cell["auroc"][sig])
                pm_b = cell.get("paired_bootstrap", {}).get(
                    f"ad_minus_{sig}") if sig != "ad" else None
                marks.append(f"{vals[-1]:.3f}")
            bars = ax.bar(x + (i - 1) * w, vals, w, color=colors[sig],
                          edgecolor="gray", lw=0.5, label=names[sig])
            for b, mk in zip(bars, marks):
                if np.isnan(b.get_height()):
                    continue
                ax.text(b.get_x() + b.get_width() / 2,
                        b.get_height() + 0.008, mk, ha="center",
                        fontsize=6, rotation=90)
        ax.axhline(0.5, color="gray", ls=":", lw=1.0)
        ax.set_xticks(x)
        ax.set_xticklabels(METRICS, fontsize=8)
        ax.set_ylim(0.35, 1.02)
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.25, axis="y")
        ax.legend(fontsize=7, loc="lower left")
    axes[0].set_ylabel("failure-detection AUROC (Top50% error)", fontsize=8)
    fig.suptitle("AD baseline vs our gates, dual error metrics "
                 f"(paired bootstrap B={res['spec']['b_boot']})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}", flush=True)


def fmt_a(v):
    return "  n/a " if v is None else f"{v:.3f}"


def write_verdict_md(res, path):
    L = []
    L.append("# AD baseline dual-metric audit - verdict "
             "(step21, prereg_ad_baseline.md)\n")
    L.append(f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- smoke: {res['smoke']}")
    L.append("- signals: AD = additive disagreement (AD-mean on single-gene "
             "holdout products per deviation D-AD1; AD-combo on Norman "
             "combos, paper-faithful per deviation D-AD2); symbol gate = "
             "two_layer_gate construction; probe = step6 spec, y refit per "
             "metric")
    L.append("- error metrics: full = transcriptome-wide MAE of delta; "
             "top20/top50 = MAE over top-20/50 genes by |true delta| "
             "(per perturbation); failure = error > median (paper Top50%)")
    L.append("- CIs: perturbation-level paired bootstrap of AUROC "
             f"differences, B={res['spec']['b_boot']}, "
             f"seed={res['spec']['boot_seed']}\n")
    L.append("## Analysis B (primary): Norman combos, archived GEARS seed4, "
             "CPU, no training\n")
    b = res["analysis_B"]
    L.append(f"- candidates={b['n_candidates']} eligible={b['n_eligible']} "
             f"predicted={b['n_predicted']} excluded={b['excluded']}")
    L.append("")
    L.append("| metric | AD-combo | symbol gate | in-situ probe | "
             "AD - symbol [95% CI] | AD - probe [95% CI] |")
    L.append("|---|---|---|---|---|---|")
    for m in METRICS:
        cell = b["per_metric"][m]
        if not cell.get("evaluable"):
            L.append(f"| {m} | n/e | | | | |")
            continue
        a = cell["auroc"]
        bs = cell["paired_bootstrap"]
        dsym = bs.get("ad_minus_symbol", {})
        dpro = bs.get("ad_minus_probe", {})
        L.append(f"| {m} | {a['ad']:.3f} | {a['symbol']:.3f} | "
                 f"{a['probe']:.3f} | "
                 f"{dsym.get('mean_diff', float('nan')):+.3f} "
                 f"[{dsym.get('ci_lo', float('nan')):+.3f},"
                 f"{dsym.get('ci_hi', float('nan')):+.3f}] {dsym.get('favors','')} | "
                 f"{dpro.get('mean_diff', float('nan')):+.3f} "
                 f"[{dpro.get('ci_lo', float('nan')):+.3f},"
                 f"{dpro.get('ci_hi', float('nan')):+.3f}] {dpro.get('favors','')} |")
    v = res["verdict"]
    L.append("")
    L.append("## Verdict (locked two-tier)\n")
    L.append(f"- tier: **{v['tier']}**")
    L.append(f"- rule: {v['rule']}")
    if "primary" in v:
        p = v["primary"]
        L.append(f"- (a) full AD-symbol CI excludes 0 favoring AD: "
                 f"{p['cond_a_full_favors_AD']} "
                 f"(mean {p['ad_minus_symbol_full']['mean_diff']:+.3f}, "
                 f"CI [{p['ad_minus_symbol_full']['ci_lo']:+.3f}, "
                 f"{p['ad_minus_symbol_full']['ci_hi']:+.3f}])")
        L.append(f"- (b) top20 reverses or non-significant: "
                 f"{p['cond_b_top20_reverses_or_nsig']} "
                 f"(mean {p['ad_minus_symbol_top20']['mean_diff']:+.3f}, "
                 f"CI [{p['ad_minus_symbol_top20']['ci_lo']:+.3f}, "
                 f"{p['ad_minus_symbol_top20']['ci_hi']:+.3f}], "
                 f"favors {p['ad_minus_symbol_top20']['favors']})")
    L.append("")
    L.append("## Analysis A (secondary): 5-dataset single-gene holdout "
             "products, AD-mean (deviation D-AD1)\n")
    L.append("| dataset | seed | metric | AD-mean | symbol gate | "
             "in-situ probe | AD - symbol [95% CI] | AD - probe [95% CI] |")
    L.append("|---|---|---|---|---|---|---|---|")
    for key, e in res["analysis_A"]["combos"].items():
        for m in METRICS:
            cell = e["per_metric"][m]
            tag = f"{e['dataset']}" if m == METRICS[0] else ""
            if not cell.get("evaluable"):
                L.append(f"| {tag} | {e['seed']} | {m} | n/e | | | | |")
                continue
            a = cell["auroc"]
            bs = cell["paired_bootstrap"]
            dsym = bs.get("ad_minus_symbol", {})
            dpro = bs.get("ad_minus_probe", {})
            L.append(
                f"| {tag} | {e['seed']} | {m} | {a['ad']:.3f} | "
                f"{a['symbol']:.3f} | {a['probe']:.3f} | "
                f"{dsym.get('mean_diff', float('nan')):+.3f} "
                f"[{dsym.get('ci_lo', float('nan')):+.3f},"
                f"{dsym.get('ci_hi', float('nan')):+.3f}] {dsym.get('favors','')} | "
                f"{dpro.get('mean_diff', float('nan')):+.3f} "
                f"[{dpro.get('ci_lo', float('nan')):+.3f},"
                f"{dpro.get('ci_hi', float('nan')):+.3f}] {dpro.get('favors','')} |")
    L.append("")
    L.append("## Deviations (locked disclosures)\n")
    L.append("- D-AD1: paper states AD is undefined for single-gene OOD. On "
             "single-gene holdout products we use AD-mean = |pred_delta - "
             "mean train single-gene delta| (mean-effect additive "
             "reference). Not the paper's per-gene AD; secondary evidence "
             "only, not in the verdict rule.")
    L.append("- D-AD2: our split makes 'both single-gene conditions in model "
             "training' and 'combo unseen by model' mutually exclusive. "
             "AD-combo baseline uses all observed single-gene conditions "
             "(incl. held-out singles; model never saw them). Evaluated "
             "combo truth never enters the baseline (no circularity). "
             "Paper home setting has both components in training; our "
             "setting is harder for the model.")
    L.append("- Dixit n=3: not evaluable under the minority>=2 rule "
             "(amendment-1), recorded as n/e like the existing probes.")
    L.append("")
    L.append("## Where this goes in the paper (one paragraph)\n")
    L.append("Gating paper appendix 'Reliability signal comparison: additive "
             "disagreement (AD) baseline and dual-metric robustness'. "
             "Norman-combo panel is the AD paper's home setting run on our "
             "archived model (no retraining); the 5-dataset table shows our "
             "probe and symbol gate recomputed under full / top20 / top50 "
             "error definitions. Both metrics reported side by side; the "
             "Top20-MAE trap (paper 0.990 -> 0.619) is checked on our "
             "gates, and the dual-metric report is presented as the "
             "standard practice regardless of verdict tier.")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"[SAVED] {path}", flush=True)


# ------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    smoke = args.smoke
    b_boot = B_BOOT_SMOKE if smoke else B_BOOT_FULL
    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)

    combos = [c for c in A_COMBOS if c["key"] in (
        {"gears_k562_s0", "gears_norman_s4", "gears_adamson_s3"}
        if smoke else {c2["key"] for c2 in A_COMBOS})]
    combo_cap = 10 if smoke else None

    rng = np.random.default_rng(BOOT_MASTER_SEED)
    t_all = time.time()
    res_a = run_analysis_a(combos, rng, b_boot)
    summary_a = summarize_a(res_a)
    res_b = run_analysis_b(rng, b_boot, combo_cap=combo_cap)
    verdict = adjudicate(res_b)

    res = {"task": "AD baseline + dual-metric robustness audit (step21)",
           "prereg": "prereg_ad_baseline.md",
           "smoke": smoke,
           "spec": {"b_boot": b_boot, "boot_seed": BOOT_MASTER_SEED,
                    "oof_seed": OOF_SEED, "metrics": METRICS,
                    "failure_label": "error > median (paper Top50%)",
                    "topk": "per-perturbation by |true delta|",
                    "ad_direction": "larger AD = more likely failure "
                                    "(paper)"},
           "analysis_A": {"combos": res_a, "dataset_summary": summary_a},
           "analysis_B": res_b,
           "verdict": verdict,
           "timing_total_s": round(time.time() - t_all, 1)}
    stem = "ad_baseline.smoke" if smoke else "ad_baseline"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[SAVED] {out_json}", flush=True)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    if not smoke:
        write_verdict_md(res, os.path.join(out_dir, "ad_baseline_verdict.md"))
    print(f"[DONE] verdict={verdict['tier']} "
          f"total={res['timing_total_s']}s", flush=True)


if __name__ == "__main__":
    main()
