# -*- coding: utf-8 -*-
"""Unified external reliability protocol for MCF7 and Norman 2019.

Arms include source-probe transfer, sign layer, target-z continuous
coefficients, single delta-max sign, in-situ OOF, 30-label few-shot probe,
and ensemble disagreement when multiple seed predictions are available.
"""
import argparse
import json
import os
import warnings

import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)
BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)
FEATURES = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
            "n_de_025", "frac_de_025", "neg_frac", "delta_entropy", "log_n_cells"]
COVERAGES = np.round(np.arange(1.0, 0.49, -0.05), 2)
B_BOOT = 1000
BOOT_SEED = 20261029
OOF_SEED = 0


def feature_row(delta, n_cells):
    abs_d = np.abs(delta)
    p = abs_d / (abs_d.sum() + 1e-12)
    return {
        "delta_l2": float(np.linalg.norm(delta)),
        "delta_l1": float(abs_d.sum()),
        "delta_max": float(abs_d.max()),
        "delta_top50_mean": float(np.sort(abs_d)[-50:].mean()),
        "n_de_025": int((abs_d > 0.25).sum()),
        "frac_de_025": float((abs_d > 0.25).mean()),
        "neg_frac": float((delta < 0).mean()),
        "delta_entropy": float(-(p * np.log(p + 1e-12)).sum() / np.log(len(delta))),
        "log_n_cells": float(np.log1p(n_cells)),
    }


def mcf7_data():
    pred = np.load(os.path.join(BASE, "gse281048_testbed/outputs/gears_predictions.npz"),
                   allow_pickle=True)
    genes = pred["genes"].tolist()
    pred_arr = pred["predictions"]
    ctrl = pred["ctrl_mean"].astype(np.float64)
    error_mask = pred["error_mask"].astype(bool)
    errors = pd.read_csv(os.path.join(BASE, "gse281048_testbed/outputs/gears_gene_errors.csv"))
    error_map = errors.set_index("gene")["median_pearson_r"].to_dict()
    pert = pd.read_csv(os.path.join(BASE, "gse281048_testbed/perturbation_list.csv"))
    kd = pert[~pert["is_NT"]]
    n_cells = kd.groupby("perturb_gene")["n_cells"].sum().to_dict()
    rows = []
    for i, g in enumerate(genes):
        if not error_mask[i] or g not in error_map:
            continue
        row = feature_row(pred_arr[i].astype(np.float64) - ctrl, n_cells.get(g, 0))
        row.update({"gene": g, "median_pearson_r": float(error_map[g]),
                    "n_cells": int(n_cells.get(g, 0))})
        rows.append(row)
    df = pd.DataFrame(rows)
    df["error"] = 1 - df["median_pearson_r"]
    df["y"] = (df["error"] > df["error"].median()).astype(int)
    return df, [pred_arr], ctrl


def norman_data():
    files = [os.path.join(OUT, f"norman_predictions_seed{s}.npz") for s in [0, 1, 2]]
    files = [f for f in files if os.path.exists(f)]
    if not files:
        raise FileNotFoundError("Run step15_norman_predict.py first")
    zs = [np.load(f, allow_pickle=True) for f in files]
    genes = zs[0]["genes"].tolist()
    preds = [z["predictions"] for z in zs]
    ctrl = zs[0]["ctrl_mean"].astype(np.float64)
    gene_names = zs[0]["gene_names"].tolist()

    adata = anndata.read_h5ad("E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/norman/perturb_processed.h5ad")
    conditions = adata.obs["condition"].astype(str).to_numpy()
    var_map = {g: i for i, g in enumerate(adata.var["gene_name"].astype(str))}
    present = [i for i, g in enumerate(gene_names) if g in var_map]
    idx = [var_map[gene_names[i]] for i in present]
    X = adata.X[:, idx]
    ctrl_present = ctrl[present]
    ctrl_cells = conditions == "ctrl"
    single_by_gene = {}
    for g in genes:
        cs = []
        for c in np.unique(conditions):
            parts = c.split("+")
            if len(parts) == 2 and "ctrl" in parts and g in parts and c != "ctrl":
                cs.append(c)
        single_by_gene[g] = cs
    needed_conditions = sorted({c for cs in single_by_gene.values() for c in cs})
    profiles = {}
    for c in needed_conditions:
        profiles[c] = np.asarray(X[conditions == c].mean(axis=0)).ravel()
    rows = []
    for i, g in enumerate(genes):
        pdelta_full = preds[0][i].astype(np.float64) - ctrl
        pdelta = pdelta_full[present]
        rs, nc = [], 0
        for c in single_by_gene[g]:
            mask = conditions == c
            nc += int(mask.sum())
            tdelta = profiles[c] - ctrl_present
            top = np.argsort(np.abs(tdelta))[-50:]
            if np.std(tdelta[top]) == 0 or np.std(pdelta[top]) == 0:
                continue
            rs.append(float(pearsonr(tdelta[top], pdelta[top])[0]))
        if not rs:
            continue
        row = feature_row(pdelta, nc)
        row.update({"gene": g, "median_pearson_r": float(np.median(rs)),
                    "n_cells": int(nc), "n_conditions": len(rs)})
        rows.append(row)
    df = pd.DataFrame(rows)
    df["error"] = 1 - df["median_pearson_r"]
    df["y"] = (df["error"] > df["error"].median()).astype(int)
    return df, preds, ctrl


def source_scores(X, source):
    Xs = (X - np.asarray(source["scaler_mean"])) / np.asarray(source["scaler_scale"])
    return Xs @ np.asarray(source["coef"]) + float(source["intercept"])


def target_z(X):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def in_situ_scores(X, y, seed=OOF_SEED):
    pipe = Pipeline([("scale", StandardScaler()),
                     ("lr", LogisticRegression(C=1.0, max_iter=2000))])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    return cross_val_predict(pipe, X, y, cv=cv, method="predict_proba")[:, 1]


def fewshot_scores(X, y, n_calib=30, repeats=10, seed=20261030):
    rng = np.random.default_rng(seed)
    n = len(y)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    sums = np.zeros(n)
    counts = np.zeros(n)
    for _ in range(repeats):
        np_, nn_ = min(n_calib // 2, len(pos)), min(n_calib - n_calib // 2, len(neg))
        calib = np.concatenate([rng.choice(pos, np_, replace=False),
                                rng.choice(neg, nn_, replace=False)])
        test = np.setdiff1d(np.arange(n), calib)
        pipe = Pipeline([("scale", StandardScaler()),
                         ("lr", LogisticRegression(C=1.0, max_iter=2000))])
        pipe.fit(X[calib], y[calib])
        s = pipe.predict_proba(X[test])[:, 1]
        sums[test] += s
        counts[test] += 1
    out = np.full(n, np.nan)
    ok = counts > 0
    out[ok] = sums[ok] / counts[ok]
    if not np.all(ok):
        out[~ok] = in_situ_scores(X, y)[~ok]
    return out


def all_scores(X, y, source, ens_disagreement=None, fewshot_repeats=10):
    coef = np.asarray(source["coef"])
    j = FEATURES.index("delta_max")
    out = {
        "source_probe": source_scores(X, source),
        "sign_layer": target_z(X) @ np.sign(coef),
        "coef_target_z": target_z(X) @ coef,
        "delta_max_sign": np.sign(coef[j]) * X[:, j],
        "in_situ": in_situ_scores(X, y),
        "fewshot_30": fewshot_scores(X, y, repeats=fewshot_repeats),
    }
    if ens_disagreement is not None:
        out["ensemble_disagreement"] = ens_disagreement
    return out


def precision_at(scores, y, cov):
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def eval_point(scores, y, error):
    out = {}
    for arm, s in scores.items():
        out[arm] = {
            "auroc": float(roc_auc_score(y, s)),
            "spearman_error": float(spearmanr(s, error).statistic),
            "precision": {f"{c:.2f}": precision_at(s, y, c) for c in COVERAGES},
        }
    out["random"] = {"precision": {f"{c:.2f}": float(y.mean()) for c in COVERAGES}}
    return out


def paired_bootstrap(X, y, error, source, ens_disagreement=None):
    rng = np.random.default_rng(BOOT_SEED)
    n = len(y)
    arms = list(all_scores(X, y, source, ens_disagreement, fewshot_repeats=2).keys())
    diffs = {f"{a}_minus_random": np.zeros((B_BOOT, len(COVERAGES))) for a in arms}
    pairs = [("sign_layer", "source_probe"), ("sign_layer", "in_situ"),
             ("coef_target_z", "sign_layer"), ("coef_target_z", "source_probe"),
             ("coef_target_z", "in_situ"), ("fewshot_30", "sign_layer"),
             ("fewshot_30", "in_situ"), ("delta_max_sign", "coef_target_z")]
    for a, b in pairs:
        diffs[f"{a}_minus_{b}"] = np.zeros((B_BOOT, len(COVERAGES)))
    metrics = {a: {"auroc": [], "spearman": []} for a in arms}
    for b in range(B_BOOT):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb, eb = X[idx], y[idx], error[idx]
        ensb = ens_disagreement[idx] if ens_disagreement is not None else None
        sc = all_scores(Xb, yb, source, ensb, fewshot_repeats=5)
        base = float(yb.mean())
        for arm, s in sc.items():
            metrics[arm]["auroc"].append(roc_auc_score(yb, s))
            metrics[arm]["spearman"].append(spearmanr(s, eb).statistic)
        prec = {a: np.array([precision_at(s, yb, c) for c in COVERAGES]) for a, s in sc.items()}
        rand = np.full(len(COVERAGES), base)
        for a in arms:
            diffs[f"{a}_minus_random"][b] = rand - prec[a]
        for a, bb in pairs:
            diffs[f"{a}_minus_{bb}"][b] = prec[a] - prec[bb]
        if (b + 1) % 200 == 0:
            print(f"bootstrap {b+1}/{B_BOOT}", flush=True)
    return metrics, diffs


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def summarize(metrics, diffs):
    out = {a: {"auroc_mean": float(np.mean(v["auroc"])),
               "auroc_ci95": ci(v["auroc"]),
               "spearman_mean": float(np.mean(v["spearman"])),
               "spearman_ci95": ci(v["spearman"])} for a, v in metrics.items()}
    dout = {}
    for name, arr in diffs.items():
        mc = arr.mean(axis=1)
        dout[name] = {"mean_over_coverage": float(mc.mean()), "ci95": ci(mc)}
    return out, dout


def main(dataset):
    if dataset == "norman":
        df, preds, ctrl = norman_data()
    else:
        df, preds, ctrl = mcf7_data()
    source = json.load(open(os.path.join(BASE, "outputs/gating_probe_transfer.json"), encoding="utf-8"))
    X = df[FEATURES].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    error = df["error"].to_numpy(dtype=float)
    ens = None
    if len(preds) > 1:
        deltas = np.stack([p - ctrl for p in preds], axis=0)
        ens = np.linalg.norm(np.std(deltas, axis=0), axis=1)
    print(dataset, "n=", len(df), "seed predictions=", len(preds), flush=True)
    scores = all_scores(X, y, source, ens, fewshot_repeats=20)
    point = eval_point(scores, y, error)
    for a, s in scores.items():
        print(a, f"AUROC={point[a]['auroc']:.3f}", f"rho={point[a]['spearman_error']:.3f}", flush=True)
    metrics, diffs = paired_bootstrap(X, y, error, source, ens)
    summary, diff_summary = summarize(metrics, diffs)
    result = {
        "task": f"Unified external reliability protocol: {dataset}",
        "n_genes": int(len(df)),
        "n_positive": int(y.sum()),
        "features": FEATURES,
        "b_boot": B_BOOT,
        "boot_seed": BOOT_SEED,
        "point": point,
        "bootstrap": summary,
        "paired": diff_summary,
        "arms": list(scores.keys()),
    }
    prefix = os.path.join(OUT, f"{dataset}_protocol")
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    export = df.copy()
    for a, s in scores.items():
        export[f"score_{a}"] = s
    export.to_csv(prefix + "_scores.csv", index=False)

    arms = list(scores.keys())
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    xs = np.arange(len(arms))
    means = [summary[a]["auroc_mean"] for a in arms]
    lo = [summary[a]["auroc_ci95"][0] for a in arms]
    hi = [summary[a]["auroc_ci95"][1] for a in arms]
    ax.bar(xs, means, yerr=[np.array(means)-np.array(lo), np.array(hi)-np.array(means)], capsize=3)
    ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(xs, arms, rotation=25, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title(f"A. {dataset}: all {len(df)} eligible genes")
    ax = axes[1]
    for a in arms:
        ax.plot(COVERAGES, [point[a]["precision"][f"{c:.2f}"] for c in COVERAGES],
                marker="o", markersize=2.5, label=a)
    ax.plot(COVERAGES, [y.mean()] * len(COVERAGES), "k--", label="random")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Precision in flagged band")
    ax.set_title("B. Precision-coverage")
    ax.invert_xaxis()
    ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()

    lines = [f"# Unified external protocol: {dataset}", "",
             f"n={len(df)}; positive={int(y.sum())}; seed predictions={len(preds)}", "",
             "| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |", "|---|---|---|---|"]
    for a in arms:
        b, d = summary[a], diff_summary[f"{a}_minus_random"]
        lines.append(f"| {a} | {b['auroc_mean']:.3f} [{b['auroc_ci95'][0]:.3f}, {b['auroc_ci95'][1]:.3f}] | {b['spearman_mean']:.3f} | {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}] |")
    lines += ["", "Key pairwise precision differences:", ""]
    for name in diff_summary:
        if name.endswith("_minus_random"):
            continue
        d = diff_summary[name]
        lines.append(f"- {name}: {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}]")
    with open(prefix + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", prefix, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", choices=["norman", "mcf7"])
    args = parser.parse_args()
    main(args.dataset)
