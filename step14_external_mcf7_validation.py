# -*- coding: utf-8 -*-
"""Step 14: external GSE281048 MCF7 stress test, all eligible genes.

Primary analysis uses every unique gene with a computable top-50 |delta|
Pearson error (no confidence panel or calibration-set exclusion).
"""
import json
import os
import warnings

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
TESTBED = os.path.join(BASE, "gse281048_testbed")
OUT = os.path.join(BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)

FEATURES = [
    "delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
    "n_de_025", "frac_de_025", "neg_frac", "delta_entropy", "log_n_cells",
]
COVERAGES = np.round(np.arange(1.0, 0.49, -0.05), 2)
B_BOOT = 1000
BOOT_SEED = 20261029
OOF_SEED = 0


def load_gene_data():
    pred = np.load(os.path.join(TESTBED, "outputs/gears_predictions.npz"),
                   allow_pickle=True)
    genes = pred["genes"].tolist()
    pred_arr = pred["predictions"]
    ctrl_mean = pred["ctrl_mean"].astype(np.float64)
    error_mask = pred["error_mask"].astype(bool)

    gene_errors = pd.read_csv(
        os.path.join(TESTBED, "outputs/gears_gene_errors.csv"))
    error_map = gene_errors.set_index("gene")["median_pearson_r"].to_dict()

    pert = pd.read_csv(os.path.join(TESTBED, "perturbation_list.csv"))
    kd = pert[~pert["is_NT"]].copy()
    pathways = kd.groupby("perturb_gene")["pathway"].apply(list).to_dict()
    n_cells = kd.groupby("perturb_gene")["n_cells"].sum().to_dict()

    rows = []
    for i, gene in enumerate(genes):
        if not error_mask[i] or gene not in error_map:
            continue
        delta = pred_arr[i].astype(np.float64) - ctrl_mean
        abs_d = np.abs(delta)
        p = abs_d / (abs_d.sum() + 1e-12)
        rows.append({
            "gene": gene,
            "pathways": ";".join(pathways.get(gene, [])),
            "n_cells": int(n_cells.get(gene, 0)),
            "median_pearson_r": float(error_map[gene]),
            "delta_l2": float(np.linalg.norm(delta)),
            "delta_l1": float(abs_d.sum()),
            "delta_max": float(abs_d.max()),
            "delta_top50_mean": float(np.sort(abs_d)[-50:].mean()),
            "n_de_025": int((abs_d > 0.25).sum()),
            "frac_de_025": float((abs_d > 0.25).mean()),
            "neg_frac": float((delta < 0).mean()),
            "delta_entropy": float(-(p * np.log(p + 1e-12)).sum()
                                   / np.log(len(delta))),
            "log_n_cells": float(np.log1p(n_cells.get(gene, 0))),
        })
    df = pd.DataFrame(rows)
    if len(df) != 99:
        raise RuntimeError(f"Expected 99 eligible genes, got {len(df)}")
    df["error"] = 1.0 - df["median_pearson_r"]
    df["y"] = (df["error"] > df["error"].median()).astype(int)
    df["y_abs_r_lt_0.5"] = (df["median_pearson_r"] < 0.5).astype(int)
    return df


def source_scores(X, source):
    Xs = (X - np.asarray(source["scaler_mean"])) / np.asarray(source["scaler_scale"])
    return Xs @ np.asarray(source["coef"]) + float(source["intercept"])


def sign_scores(X, signs):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return ((X - mu) / sd) @ signs


def delta_max_sign_scores(X, source):
    j = FEATURES.index("delta_max")
    return np.sign(source["coef"][j]) * X[:, j]


def in_situ_scores(X, y, seed=OOF_SEED):
    pipe = Pipeline([
        ("scale", StandardScaler()),
        ("lr", LogisticRegression(C=1.0, max_iter=2000)),
    ])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    return cross_val_predict(pipe, X, y, cv=cv, method="predict_proba")[:, 1]


def precision_at(scores, y, cov):
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def all_scores(X, y, source):
    signs = np.sign(np.asarray(source["coef"]))
    return {
        "source_probe": source_scores(X, source),
        "sign_layer": sign_scores(X, signs),
        "delta_max_sign": delta_max_sign_scores(X, source),
        "in_situ": in_situ_scores(X, y),
    }


def eval_point(scores, y, error, abs_y):
    out = {}
    for arm, s in scores.items():
        out[arm] = {
            "auroc": float(roc_auc_score(y, s)),
            "auroc_abs_r_lt_0.5": float(roc_auc_score(abs_y, s)),
            "spearman_error": float(spearmanr(s, error).statistic),
            "precision": {f"{c:.2f}": precision_at(s, y, c) for c in COVERAGES},
        }
    out["random"] = {
        "precision": {f"{c:.2f}": float(y.mean()) for c in COVERAGES}
    }
    return out


def paired_bootstrap(X, y, error, abs_y, source):
    rng = np.random.default_rng(BOOT_SEED)
    n = len(y)
    arms = ["source_probe", "sign_layer", "delta_max_sign", "in_situ"]
    diffs = {}
    for a in ["sign_layer", "in_situ", "source_probe", "delta_max_sign"]:
        diffs[f"{a}_minus_random"] = np.zeros((B_BOOT, len(COVERAGES)))
    for a, b in [("sign_layer", "source_probe"),
                 ("sign_layer", "in_situ"),
                 ("source_probe", "in_situ")]:
        diffs[f"{a}_minus_{b}"] = np.zeros((B_BOOT, len(COVERAGES)))
    metrics = {a: {"auroc": [], "spearman": [], "auroc_abs": []} for a in arms}

    for b in range(B_BOOT):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2 or len(np.unique(abs_y[idx])) < 2:
            idx = rng.integers(0, n, n)
        Xb, yb, eb, ayb = X[idx], y[idx], error[idx], abs_y[idx]
        sc = all_scores(Xb, yb, source)
        base_rate = float(yb.mean())
        for arm, s in sc.items():
            metrics[arm]["auroc"].append(roc_auc_score(yb, s))
            metrics[arm]["spearman"].append(spearmanr(s, eb).statistic)
            metrics[arm]["auroc_abs"].append(roc_auc_score(ayb, s))
        prec = {arm: np.array([precision_at(s, yb, c) for c in COVERAGES])
                for arm, s in sc.items()}
        random_prec = np.full(len(COVERAGES), base_rate)
        for arm in ["sign_layer", "in_situ", "source_probe", "delta_max_sign"]:
            # risk = 1 - precision; arm - random = random_prec - arm_prec
            diffs[f"{arm}_minus_random"][b] = random_prec - prec[arm]
        for a, bb in [("sign_layer", "source_probe"),
                      ("sign_layer", "in_situ"),
                      ("source_probe", "in_situ")]:
            diffs[f"{a}_minus_{bb}"][b] = prec[a] - prec[bb]
        if (b + 1) % 100 == 0:
            print(f"bootstrap {b + 1}/{B_BOOT}", flush=True)
    return metrics, diffs


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def summarize(metrics, diffs):
    out = {}
    for arm, vals in metrics.items():
        out[arm] = {
            "auroc_mean": float(np.mean(vals["auroc"])),
            "auroc_ci95": ci(vals["auroc"]),
            "spearman_mean": float(np.mean(vals["spearman"])),
            "spearman_ci95": ci(vals["spearman"]),
            "auroc_abs_r_lt_0.5_mean": float(np.mean(vals["auroc_abs"])),
            "auroc_abs_r_lt_0.5_ci95": ci(vals["auroc_abs"]),
        }
    diff_out = {}
    for name, arr in diffs.items():
        mean_cov = arr.mean(axis=1)
        diff_out[name] = {
            "mean_over_coverage": float(mean_cov.mean()),
            "ci95": ci(mean_cov),
            "mean_by_coverage": {f"{c:.2f}": float(v) for c, v in
                                 zip(COVERAGES, arr.mean(axis=0))},
        }
    return out, diff_out


def pathway_sensitivity(df, scores):
    errors = pd.read_csv(os.path.join(TESTBED, "outputs/gears_errors.csv"))
    work = errors.merge(df[["gene", "y", "error"]], on="gene", how="inner")
    rows = []
    for pathway, sub in work.groupby("pathway"):
        if len(sub) < 10 or sub["y"].nunique() < 2:
            continue
        idx = sub.index.to_numpy()
        # Map condition rows back to per-gene score vectors.
        gene_to_pos = {g: i for i, g in enumerate(df["gene"])}
        pos = np.array([gene_to_pos[g] for g in sub["gene"]])
        row = {"pathway": pathway, "n_conditions": int(len(sub)),
               "n_genes": int(sub["gene"].nunique())}
        for arm, s in scores.items():
            row[arm] = float(roc_auc_score(sub["y"].to_numpy(), s[pos]))
        rows.append(row)
    return rows


def main():
    df = load_gene_data()
    source = json.load(open(os.path.join(BASE, "outputs/gating_probe_transfer.json"),
                            encoding="utf-8"))
    X = df[FEATURES].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    error = df["error"].to_numpy(dtype=float)
    abs_y = df["y_abs_r_lt_0.5"].to_numpy(dtype=int)

    scores = all_scores(X, y, source)
    point = eval_point(scores, y, error, abs_y)
    print("Point AUROC:")
    for arm in scores:
        print(arm, point[arm]["auroc"], "Spearman", point[arm]["spearman_error"])

    metrics, diffs = paired_bootstrap(X, y, error, abs_y, source)
    summary, diff_summary = summarize(metrics, diffs)
    pathways = pathway_sensitivity(df, scores)
    j_max = FEATURES.index("delta_max")
    diagnostics = {
        "source_vs_delta_max_pearson": float(np.corrcoef(
            scores["source_probe"], scores["delta_max_sign"])[0, 1]),
        "source_delta_max_scaler_mean": float(source["scaler_mean"][j_max]),
        "mcf7_delta_max_median": float(df["delta_max"].median()),
        "mcf7_delta_max_range": [float(df["delta_max"].min()),
                                 float(df["delta_max"].max())],
        "absolute_label_note": (
            f"The absolute r<0.5 label has {int(abs_y.sum())}/{len(abs_y)} positives "
            "and is therefore not used for interpretation."),
    }

    result = {
        "task": "External GSE281048 MCF7 validation, all eligible genes",
        "plan": "prereg_external_mcf7_v1.md",
        "n_genes": int(len(df)),
        "n_positive_median_split": int(y.sum()),
        "n_positive_abs_r_lt_0.5": int(abs_y.sum()),
        "features": FEATURES,
        "b_boot": B_BOOT,
        "boot_seed": BOOT_SEED,
        "oof_seed": OOF_SEED,
        "coverage_grid": [float(c) for c in COVERAGES],
        "point": point,
        "bootstrap": summary,
        "paired": diff_summary,
        "pathway_sensitivity_condition_level": pathways,
        "diagnostics": diagnostics,
        "caveats": [
            "GSE281048 MCF7 uses pathway-targeted Perturb-seq with 24 h cytokine stimulation; it is an external stress test, not a same-protocol Replogle validation.",
            "Primary analysis uses all 99 eligible unique genes; no confidence panel or calibration-set exclusion.",
            "The in_situ arm uses MCF7 labels and is a comparator, not a zero-label deployment method.",
            "Intervals are nominal perturbation-level bootstrap intervals; in_situ is refit per resample.",
        ],
    }
    with open(os.path.join(OUT, "mcf7_validation.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    export = df.copy()
    for arm, s in scores.items():
        export[f"score_{arm}"] = s
    export.to_csv(os.path.join(OUT, "mcf7_validation_scores.csv"), index=False)

    # Figure: AUROC, precision-coverage, pathway sensitivity.
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    arms = ["source_probe", "sign_layer", "delta_max_sign", "in_situ"]
    labels = ["K562 probe transfer", "sign layer", "delta-max sign", "MCF7 in-situ"]
    ax = axes[0]
    xs = np.arange(len(arms))
    means = [summary[a]["auroc_mean"] for a in arms]
    lows = [summary[a]["auroc_ci95"][0] for a in arms]
    highs = [summary[a]["auroc_ci95"][1] for a in arms]
    ax.bar(xs, means, yerr=[np.array(means) - np.array(lows),
                            np.array(highs) - np.array(means)], capsize=4,
           color=["#7c8aa0", "#2b8a3e", "#74a892", "#b23a48"])
    ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(xs, labels, rotation=20, ha="right")
    ax.set_ylabel("AUROC (median-split error)")
    ax.set_title("A. External MCF7 validation (all 99 genes)")

    ax = axes[1]
    for arm, label in zip(arms, labels):
        vals = [point[arm]["precision"][f"{c:.2f}"] for c in COVERAGES]
        ax.plot(COVERAGES, vals, marker="o", markersize=3, label=label)
    ax.plot(COVERAGES, [y.mean()] * len(COVERAGES), "k--", label="random")
    ax.set_xlabel("Coverage (fraction flagged)")
    ax.set_ylabel("Precision in flagged band")
    ax.set_title("B. Precision-coverage")
    ax.legend(fontsize=8)
    ax.invert_xaxis()

    ax = axes[2]
    if pathways:
        pw = [r["pathway"] for r in pathways]
        x = np.arange(len(pw))
        width = 0.27
        for k, (arm, label) in enumerate([("source_probe", "K562 probe"),
                                          ("sign_layer", "sign layer"),
                                          ("in_situ", "in-situ")]):
            vals = [r[arm] for r in pathways]
            ax.bar(x + (k - 1) * width, vals, width, label=label)
        ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
        ax.set_xticks(x, pw, rotation=0)
        ax.set_ylabel("Condition-level AUROC")
        ax.set_title("C. Pathway sensitivity (descriptive)")
        ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "mcf7_validation.png"), dpi=300)
    plt.close()

    # Compact Markdown report.
    lines = [
        "# External MCF7 validation (GSE281048)",
        "",
        f"Primary set: all {len(df)} eligible unique genes; median-split positives={int(y.sum())}; absolute r<0.5 positives={int(abs_y.sum())}.",
        "",
        "| Arm | AUROC (bootstrap mean [95% CI]) | Spearman(error) | sign-minus-random risk (coverage mean [95% CI]) |",
        "|---|---|---|---|",
    ]
    for arm, label in zip(arms, labels):
        b = summary[arm]
        d = diff_summary[f"{arm}_minus_random"]
        lines.append(f"| {label} | {b['auroc_mean']:.3f} [{b['auroc_ci95'][0]:.3f}, {b['auroc_ci95'][1]:.3f}] | {b['spearman_mean']:.3f} | {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}] |")
    lines += [
        "",
        "Pairwise coverage-averaged precision differences:",
        "",
    ]
    for name in ["sign_layer_minus_source_probe", "sign_layer_minus_in_situ",
                 "source_probe_minus_in_situ"]:
        d = diff_summary[name]
        lines.append(f"- {name}: {d['mean_over_coverage']:.4f} [{d['ci95'][0]:.4f}, {d['ci95'][1]:.4f}]")
    lines += [
        "",
        f"Diagnostic: source-probe and delta-max-sign score correlation = {diagnostics['source_vs_delta_max_pearson']:.6f}; source delta_max scaler mean = {diagnostics['source_delta_max_scaler_mean']:.3f}, MCF7 delta_max median = {diagnostics['mcf7_delta_max_median']:.2f}. The transferred full probe is therefore effectively dominated by delta_max on MCF7.",
        "",
        diagnostics["absolute_label_note"],
        "",
        "Caveat: GSE281048 MCF7 contains 24 h cytokine stimulation and is an external stress test, not a same-protocol Replogle validation. The in-situ arm uses MCF7 labels. Intervals are nominal and unadjusted for multiplicity.",
    ]
    with open(os.path.join(OUT, "mcf7_validation.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("Saved outputs to", OUT)


if __name__ == "__main__":
    main()
