# -*- coding: utf-8 -*-
"""Warm-start target calibration from archived K562 coefficients."""
import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

import step16_external_protocol as base

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = os.path.join(base.BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)
COVERAGES = base.COVERAGES
B_BOOT = 500
BOOT_SEED = 20261029
FEATURES = base.FEATURES


def warm_fit(X, y, init, C=1.0):
    init = np.asarray(init, dtype=float)
    n, d = X.shape
    X1 = np.c_[X, np.ones(n)]

    def obj(w):
        z = X1 @ w
        p = expit(z)
        eps = 1e-12
        loss = -np.mean(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))
        reg = 0.5 * np.sum(w[:-1] ** 2) / C
        return loss + reg

    def jac(w):
        z = X1 @ w
        p = expit(z)
        g = X1.T @ (p - y) / n
        g[:-1] += w[:-1] / C
        return g

    res = minimize(obj, init, jac=jac, method="L-BFGS-B", options={"maxiter": 2000})
    return res.x


def split_scores(X, y, source, n_calib=30, repeats=20, seed=20261030):
    rng = np.random.default_rng(seed)
    n = len(y)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    warm_sum = np.zeros(n)
    warm_cnt = np.zeros(n)
    scratch_sum = np.zeros(n)
    scratch_cnt = np.zeros(n)
    for _ in range(repeats):
        np_, nn_ = min(n_calib // 2, len(pos)), min(n_calib - n_calib // 2, len(neg))
        calib = np.concatenate([rng.choice(pos, np_, replace=False),
                                rng.choice(neg, nn_, replace=False)])
        test = np.setdiff1d(np.arange(n), calib)
        sc = StandardScaler().fit(X[calib])
        Xc = sc.transform(X[calib])
        Xt = sc.transform(X[test])
        init = np.r_[source["coef"], source["intercept"]]
        w = warm_fit(Xc, y[calib], init, C=1.0)
        warm = expit(np.c_[Xt, np.ones(len(test))] @ w)
        warm_sum[test] += warm
        warm_cnt[test] += 1
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xc, y[calib])
        scratch = clf.predict_proba(Xt)[:, 1]
        scratch_sum[test] += scratch
        scratch_cnt[test] += 1
    warm = np.full(n, np.nan)
    scratch = np.full(n, np.nan)
    okw, oks = warm_cnt > 0, scratch_cnt > 0
    warm[okw] = warm_sum[okw] / warm_cnt[okw]
    scratch[oks] = scratch_sum[oks] / scratch_cnt[oks]
    if not np.all(okw):
        warm[~okw] = base.in_situ_scores(X, y)[~okw]
    if not np.all(oks):
        scratch[~oks] = base.in_situ_scores(X, y)[~oks]
    return warm, scratch


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def paired_bootstrap(X, y, error, source):
    rng = np.random.default_rng(BOOT_SEED)
    n = len(y)
    diffs = {"warmstart_30_minus_fewshot_30": np.zeros((B_BOOT, len(COVERAGES))),
             "warmstart_30_minus_coef_target_z": np.zeros((B_BOOT, len(COVERAGES))),
             "warmstart_30_minus_random": np.zeros((B_BOOT, len(COVERAGES)))}
    metrics = {a: {"auroc": [], "spearman": []} for a in ["warmstart_30", "fewshot_30", "coef_target_z"]}
    for b in range(B_BOOT):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb, eb = X[idx], y[idx], error[idx]
        warm, scratch = split_scores(Xb, yb, source, repeats=5)
        z = base.target_z(Xb)
        sc = {"warmstart_30": warm, "fewshot_30": scratch,
              "coef_target_z": z @ np.asarray(source["coef"])}
        base_rate = float(yb.mean())
        for a, s in sc.items():
            metrics[a]["auroc"].append(roc_auc_score(yb, s))
            metrics[a]["spearman"].append(spearmanr(s, eb).statistic)
        prec = {a: np.array([base.precision_at(s, yb, c) for c in COVERAGES]) for a, s in sc.items()}
        diffs["warmstart_30_minus_fewshot_30"][b] = prec["warmstart_30"] - prec["fewshot_30"]
        diffs["warmstart_30_minus_coef_target_z"][b] = prec["warmstart_30"] - prec["coef_target_z"]
        diffs["warmstart_30_minus_random"][b] = base_rate - prec["warmstart_30"]
        if (b + 1) % 100 == 0:
            print(f"bootstrap {b+1}/{B_BOOT}", flush=True)
    return metrics, diffs


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
        df, preds, ctrl = base.norman_data()
    else:
        df, preds, ctrl = base.mcf7_data()
    source = json.load(open(os.path.join(base.BASE, "outputs/gating_probe_transfer.json"), encoding="utf-8"))
    X = df[FEATURES].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    error = df["error"].to_numpy(dtype=float)
    warm, scratch = split_scores(X, y, source, repeats=20)
    z = base.target_z(X)
    sc = {"warmstart_30": warm, "fewshot_30": scratch,
          "coef_target_z": z @ np.asarray(source["coef"]),
          "in_situ": base.in_situ_scores(X, y)}
    point = base.eval_point(sc, y, error)
    for a in sc:
        print(dataset, a, f"{point[a]['auroc']:.3f}", flush=True)
    metrics, diffs = paired_bootstrap(X, y, error, source)
    summary, diff_summary = summarize(metrics, diffs)

    # Label-efficiency curve
    curve = []
    for nc in [10, 15, 30, 60]:
        w, s = split_scores(X, y, source, n_calib=nc, repeats=10)
        curve.append({"n_calib": nc,
                      "warmstart_auroc": float(roc_auc_score(y, w)),
                      "fewshot_auroc": float(roc_auc_score(y, s))})
    result = {"task": f"Warm-start transfer calibration: {dataset}",
              "n_genes": int(len(df)), "point": point, "bootstrap": summary,
              "paired": diff_summary, "label_curve": curve,
              "b_boot": B_BOOT, "boot_seed": BOOT_SEED}
    prefix = os.path.join(OUT, f"{dataset}_warmstart")
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    export = df.copy()
    for a, s in sc.items():
        export[f"score_{a}"] = s
    export.to_csv(prefix + "_scores.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    arms = ["warmstart_30", "fewshot_30", "coef_target_z"]
    xs = np.arange(len(arms))
    means = [summary[a]["auroc_mean"] for a in arms]
    lo = [summary[a]["auroc_ci95"][0] for a in arms]
    hi = [summary[a]["auroc_ci95"][1] for a in arms]
    ax.bar(xs, means, yerr=[np.array(means)-np.array(lo), np.array(hi)-np.array(means)], capsize=4,
           color=["#2b8a3e", "#7c8aa0", "#74a892", "#b23a48"])
    ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(xs, arms, rotation=20, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title(f"A. {dataset}: warm-start calibration")
    ax = axes[1]
    nc = [r["n_calib"] for r in curve]
    ax.plot(nc, [r["warmstart_auroc"] for r in curve], "o-", label="warm start")
    ax.plot(nc, [r["fewshot_auroc"] for r in curve], "s--", label="from scratch")
    ax.set_xlabel("Number of target labels")
    ax.set_ylabel("AUROC")
    ax.set_title("B. Label efficiency")
    ax.legend()
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()

    lines = [f"# Warm-start transfer calibration: {dataset}", "",
             "| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |", "|---|---|---|---|"]
    for a in ["warmstart_30", "fewshot_30", "coef_target_z"]:
        b, d = summary[a], diff_summary["warmstart_30_minus_random"]
        lines.append(f"| {a} | {b['auroc_mean']:.3f} [{b['auroc_ci95'][0]:.3f}, {b['auroc_ci95'][1]:.3f}] | {b['spearman_mean']:.3f} | {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}] |")
    for name, d in diff_summary.items():
        lines.append(f"- {name}: {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}]")
    lines += ["", "Label efficiency:", ""]
    for r in curve:
        lines.append(f"- n={r['n_calib']}: warm={r['warmstart_auroc']:.3f}, scratch={r['fewshot_auroc']:.3f}")
    with open(prefix + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", prefix, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", choices=["norman", "mcf7"])
    args = parser.parse_args()
    main(args.dataset)
