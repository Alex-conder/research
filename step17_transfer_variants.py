# -*- coding: utf-8 -*-
"""Route B: robust target-normalized transfer variants on Norman and MCF7."""
import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

import step16_external_protocol as base

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = os.path.join(base.BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)
B_BOOT = 1000
BOOT_SEED = 20261029
COVERAGES = base.COVERAGES
FEATURES = base.FEATURES
STABLE2 = [FEATURES.index("delta_max"), FEATURES.index("delta_entropy")]
MAG4 = [FEATURES.index(x) for x in ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean"]]
JMAX = FEATURES.index("delta_max")
NO_MAX = [i for i in range(len(FEATURES)) if i != JMAX]


def median_iqr(X):
    med = np.median(X, axis=0)
    iqr = np.percentile(X, 75, axis=0) - np.percentile(X, 25, axis=0)
    iqr[iqr == 0] = 1.0
    return (X - med) / iqr


def rankdata(X):
    return np.apply_along_axis(lambda c: pd.Series(c).rank(pct=True).to_numpy(), 0, X)


def scores(X, y, source):
    coef = np.asarray(source["coef"])
    z = base.target_z(X)
    mz = median_iqr(X)
    rz = rankdata(X)
    return {
        "sign_layer": z @ np.sign(coef),
        "coef_target_z": z @ coef,
        "coef_median_iqr": mz @ coef,
        "coef_rank": rz @ coef,
        "stable2_z": z[:, STABLE2] @ coef[STABLE2],
        "stable2_rank": rz[:, STABLE2] @ coef[STABLE2],
        "delta_max_rank": np.sign(coef[JMAX]) * rz[:, JMAX],
        "magnitude4_z": z[:, MAG4] @ coef[MAG4],
        "no_delta_max_z": z[:, NO_MAX] @ coef[NO_MAX],
        "in_situ": base.in_situ_scores(X, y),
        "fewshot_30": base.fewshot_scores(X, y, repeats=10),
    }


def paired_bootstrap(X, y, error, source):
    rng = np.random.default_rng(BOOT_SEED)
    n = len(y)
    arms = list(scores(X, y, source).keys())
    diffs = {f"{a}_minus_random": np.zeros((B_BOOT, len(COVERAGES))) for a in arms}
    pairs = [("coef_rank", "coef_target_z"), ("coef_rank", "delta_max_rank"),
             ("stable2_z", "coef_target_z"), ("stable2_rank", "coef_target_z"),
             ("stable2_rank", "delta_max_rank"), ("coef_median_iqr", "coef_target_z"),
             ("delta_max_rank", "delta_max_rank"), ("no_delta_max_z", "coef_target_z")]
    # delta_max_sign is included as an alias for delta_max_rank when useful.
    for a, b in pairs:
        diffs[f"{a}_minus_{b}"] = np.zeros((B_BOOT, len(COVERAGES)))
    metrics = {a: {"auroc": [], "spearman": []} for a in arms}
    for b in range(B_BOOT):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb, eb = X[idx], y[idx], error[idx]
        sc = scores(Xb, yb, source)
        base_rate = float(yb.mean())
        for a, s in sc.items():
            metrics[a]["auroc"].append(roc_auc_score(yb, s))
            metrics[a]["spearman"].append(spearmanr(s, eb).statistic)
        prec = {a: np.array([base.precision_at(s, yb, c) for c in COVERAGES]) for a, s in sc.items()}
        rand = np.full(len(COVERAGES), base_rate)
        for a in arms:
            diffs[f"{a}_minus_random"][b] = rand - prec[a]
        for a, bb in pairs:
            diffs[f"{a}_minus_{bb}"][b] = prec[a] - prec[bb]
        if (b + 1) % 250 == 0:
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
        df, preds, ctrl = base.norman_data()
    else:
        df, preds, ctrl = base.mcf7_data()
    source = json.load(open(os.path.join(base.BASE, "outputs/gating_probe_transfer.json"), encoding="utf-8"))
    X = df[FEATURES].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    error = df["error"].to_numpy(dtype=float)
    sc = scores(X, y, source)
    point = base.eval_point(sc, y, error)
    for a in sc:
        print(dataset, a, f"{point[a]['auroc']:.3f}", flush=True)
    metrics, diffs = paired_bootstrap(X, y, error, source)
    summary, diff_summary = summarize(metrics, diffs)
    result = {"task": f"Robust transfer variants: {dataset}",
              "n_genes": int(len(df)), "n_positive": int(y.sum()),
              "point": point, "bootstrap": summary, "paired": diff_summary,
              "b_boot": B_BOOT, "boot_seed": BOOT_SEED}
    prefix = os.path.join(OUT, f"{dataset}_variants")
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    export = df.copy()
    for a, s in sc.items():
        export[f"score_{a}"] = s
    export.to_csv(prefix + "_scores.csv", index=False)

    arms = list(sc.keys())
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    ax = axes[0]
    xs = np.arange(len(arms))
    means = [summary[a]["auroc_mean"] for a in arms]
    lo = [summary[a]["auroc_ci95"][0] for a in arms]
    hi = [summary[a]["auroc_ci95"][1] for a in arms]
    ax.bar(xs, means, yerr=[np.array(means)-np.array(lo), np.array(hi)-np.array(means)], capsize=3)
    ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(xs, arms, rotation=25, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title(f"A. {dataset}: robust transfer variants")
    ax = axes[1]
    for a in ["sign_layer", "coef_target_z", "coef_rank", "stable2_rank", "delta_max_rank", "in_situ", "fewshot_30"]:
        ax.plot(COVERAGES, [point[a]["precision"][f"{c:.2f}"] for c in COVERAGES],
                marker="o", markersize=2.5, label=a)
    ax.plot(COVERAGES, [y.mean()] * len(COVERAGES), "k--", label="random")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Precision in flagged band")
    ax.set_title("B. Precision-coverage")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()

    lines = [f"# Robust transfer variants: {dataset}", "",
             "| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |", "|---|---|---|---|"]
    for a in arms:
        b, d = summary[a], diff_summary[f"{a}_minus_random"]
        lines.append(f"| {a} | {b['auroc_mean']:.3f} [{b['auroc_ci95'][0]:.3f}, {b['auroc_ci95'][1]:.3f}] | {b['spearman_mean']:.3f} | {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}] |")
    lines += ["", "Key paired precision differences:", ""]
    for name, d in diff_summary.items():
        if name.endswith("_minus_random"):
            continue
        lines.append(f"- {name}: {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}]")
    with open(prefix + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", prefix, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", choices=["norman", "mcf7"])
    args = parser.parse_args()
    main(args.dataset)
