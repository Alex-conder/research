# -*- coding: utf-8 -*-
"""Reliability-aware conformal gating on Norman 2019.

Combines split-half ground-truth reliability with split-conformal p-values.
"""
import json
import os
import warnings

import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

import step16_external_protocol as base

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = os.path.join(base.BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)
COVERAGES = base.COVERAGES
FEATURES = base.FEATURES
B_SPLIT = 20


def split_half_reliability(adata, pred_genes, gene_names, ctrl_mean, n_split=B_SPLIT):
    conditions = adata.obs["condition"].astype(str).to_numpy()
    var_map = {g: i for i, g in enumerate(adata.var["gene_name"].astype(str))}
    present = [i for i, g in enumerate(gene_names) if g in var_map]
    idx = [var_map[gene_names[i]] for i in present]
    X = adata.X[:, idx]
    ctrl_present = ctrl_mean[present]
    ctrl = conditions == "ctrl"
    ctrl_idx = np.where(ctrl)[0]
    rng = np.random.default_rng(20261030)
    rows = []
    for g in pred_genes:
        cs = []
        for c in np.unique(conditions):
            parts = c.split("+")
            if len(parts) == 2 and "ctrl" in parts and g in parts and c != "ctrl":
                cs.append(c)
        if not cs:
            continue
        # full effect for specificity
        full_effects = []
        for c in cs:
            full_effects.append(np.asarray(X[conditions == c].mean(axis=0)).ravel() - ctrl_present)
        full_effect = np.median(np.stack(full_effects), axis=0)
        rs = []
        for _ in range(n_split):
            perm = rng.permutation(ctrl_idx)
            ctrl_a, ctrl_b = perm[: len(perm) // 2], perm[len(perm) // 2:]
            for c in cs:
                idx_c = np.where(conditions == c)[0]
                perm_c = rng.permutation(idx_c)
                ca, cb = perm_c[: len(perm_c) // 2], perm_c[len(perm_c) // 2:]
                ea = np.asarray(X[ca].mean(axis=0)).ravel() - np.asarray(X[ctrl_a].mean(axis=0)).ravel()
                eb = np.asarray(X[cb].mean(axis=0)).ravel() - np.asarray(X[ctrl_b].mean(axis=0)).ravel()
                if np.std(ea) == 0 or np.std(eb) == 0:
                    continue
                rs.append(float(pearsonr(ea, eb)[0]))
        rows.append({"gene": g, "effect": full_effect, "reliability": float(np.median(rs) if rs else np.nan)})
    # shared axis: mean effect across reliable perturbations
    reliable_effects = np.stack([r["effect"] for r in rows if r["reliability"] >= 0.5])
    shared_axis = reliable_effects.mean(axis=0)
    shared_axis /= np.linalg.norm(shared_axis) + 1e-12
    for r in rows:
        e = r["effect"]
        cos = float(np.dot(e, shared_axis) / (np.linalg.norm(e) + 1e-12))
        r["phi"] = cos * cos
        if r["reliability"] < 0.5:
            r["quality"] = "unreliable"
        elif r["phi"] >= 0.5:
            r["quality"] = "shared"
        else:
            r["quality"] = "specific"
        del r["effect"]
    return pd.DataFrame(rows)


def conformal_scores(cal_scores, test_scores):
    cal = np.sort(np.asarray(cal_scores))
    n = len(cal)
    out = np.empty(len(test_scores))
    for i, s in enumerate(test_scores):
        out[i] = (1 + np.sum(cal >= s)) / (n + 1)
    return out


def main():
    df, preds, ctrl = base.norman_data()
    source = json.load(open(os.path.join(base.BASE, "outputs/gating_probe_transfer.json"), encoding="utf-8"))
    adata = anndata.read_h5ad(os.path.join(base.BASE, "data/gears_demo/norman/perturb_processed.h5ad"))
    pred = np.load(os.path.join(OUT, "norman_predictions_seed0.npz"), allow_pickle=True)
    genes = pred["genes"].tolist()
    gene_names = pred["gene_names"].tolist()
    rel = split_half_reliability(adata, genes, gene_names, ctrl)
    work = df.merge(rel, on="gene", how="left")
    work["quality"] = work["quality"].fillna("unreliable")

    # Conformal p-values for each arm, with n=30 calibration per split.
    rng = np.random.default_rng(20261030)
    y = work["y"].to_numpy(dtype=int)
    X = work[FEATURES].to_numpy(dtype=float)
    z = base.target_z(X)
    coef = np.asarray(source["coef"])
    j = FEATURES.index("delta_max")
    arms = {
        "source_probe": base.source_scores(X, source),
        "sign_layer": z @ np.sign(coef),
        "coef_target_z": z @ coef,
        "delta_max_sign": np.sign(coef[j]) * X[:, j],
    }
    result_scores = {}
    for arm, scores in arms.items():
        pvals = np.full(len(y), np.nan)
        for _ in range(20):
            # Conformal calibration must be unstratified to preserve exchangeability.
            calib = rng.choice(len(y), 60, replace=False)
            test = np.setdiff1d(np.arange(len(y)), calib)
            normal = scores[calib][y[calib] == 0]
            pvals[test] = conformal_scores(normal, scores[test])
        pvals[np.isnan(pvals)] = 1.0  # least anomalous if never held out
        result_scores[arm] = 1.0 - pvals  # higher = more anomalous

    point = base.eval_point(result_scores, y, work["error"].to_numpy(dtype=float))
    for a in result_scores:
        print(a, f"AUROC={point[a]['auroc']:.3f}", flush=True)

    # Reliability diagnostics
    alpha_diag = {}
    for a, s in result_scores.items():
        p = 1.0 - s
        row = {}
        for alpha in [0.05, 0.10, 0.20, 0.30, 0.50]:
            flag = p <= alpha
            row[str(alpha)] = {
                "coverage": float(flag.mean()),
                "flagged_fdp": float(y[flag].mean()) if flag.any() else None,
            }
        alpha_diag[a] = row
    diag = {
        "quality_counts": work["quality"].value_counts().to_dict(),
        "reliability_median": float(work["reliability"].median()),
        "reliability_quantiles": [float(x) for x in work["reliability"].quantile([0.1, 0.25, 0.5, 0.75, 0.9]).to_list()],
        "error_by_quality": work.groupby("quality")["error"].agg(["mean", "median", "count"]).to_dict(),
        "score_reliability_spearman": {a: float(spearmanr(result_scores[a], work["reliability"]).statistic) for a in result_scores},
        "error_reliability_spearman": float(spearmanr(work["error"], work["reliability"]).statistic),
        "conformal_alpha": alpha_diag,
    }
    result = {"task": "Reliability-aware conformal gating on Norman",
              "n_genes": int(len(work)), "point": point, "diagnostics": diag,
              "quality_counts": work["quality"].value_counts().to_dict()}
    prefix = os.path.join(OUT, "norman_reliability_conformal")
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    export = work.copy()
    for a, s in result_scores.items():
        export[f"score_{a}"] = s
    export.to_csv(prefix + "_scores.csv", index=False)

    # Figure: reliability histogram + conformal precision-coverage
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    ax = axes[0]
    ax.hist(work["reliability"], bins=24, color="#4c78a8")
    ax.axvline(0.5, color="k", linestyle="--", label="reliability threshold 0.5")
    ax.set_xlabel("split-half reliability")
    ax.set_ylabel("number of perturbations")
    ax.set_title("A. Ground-truth reliability on Norman")
    ax.legend(fontsize=8)
    ax = axes[1]
    for a in result_scores:
        ax.plot(COVERAGES, [point[a]["precision"][f"{c:.2f}"] for c in COVERAGES], marker="o", markersize=2.5, label=a)
    ax.plot(COVERAGES, [y.mean()] * len(COVERAGES), "k--", label="random")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Precision in flagged band")
    ax.set_title("B. Conformal risk-coverage")
    ax.invert_xaxis()
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()

    lines = ["# Reliability-aware conformal gating on Norman", "",
             f"Quality counts: {diag['quality_counts']}", "",
             f"Error-reliability Spearman: {diag['error_reliability_spearman']:.3f}", "",
             "| Arm | conformal AUROC | reliability Spearman |", "|---|---|---|"]
    for a in result_scores:
        lines.append(f"| {a} | {point[a]['auroc']:.3f} | {diag['score_reliability_spearman'][a]:.3f} |")
    lines += ["", "Conformal alpha=0.20 diagnostics:", ""]
    for a in result_scores:
        d = alpha_diag[a]["0.2"]
        lines.append(f"- {a}: coverage={d['coverage']:.3f}, flagged FDP={d['flagged_fdp']:.3f}")
    with open(prefix + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", prefix)


if __name__ == "__main__":
    main()
