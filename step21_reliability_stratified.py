# -*- coding: utf-8 -*-
"""Reliability-stratified gate evaluation on MCF7 and Replogle datasets."""
import argparse
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
from step10_signature_resolution import load_gears

warnings.filterwarnings("ignore", category=RuntimeWarning)
OUT = os.path.join(base.BASE, "outputs_external")
os.makedirs(OUT, exist_ok=True)
FEATURES = base.FEATURES
B_SPLIT = 20


def reliability(adata, pred_genes, gene_names, ctrl_mean, mode):
    conditions = adata.obs["condition"].astype(str).to_numpy()
    var_map = {g: i for i, g in enumerate(adata.var["gene_name"].astype(str))}
    present = [i for i, g in enumerate(gene_names) if g in var_map]
    idx = [var_map[gene_names[i]] for i in present]
    X = adata.X[:, idx]
    ctrl_present = ctrl_mean[present]
    rng = np.random.default_rng(20261030)
    rows = []
    for g in pred_genes:
        cs = []
        for c in np.unique(conditions):
            parts = c.split("+") if mode == "plus" else c.split("@")
            if len(parts) == 2 and "ctrl" in parts and g in parts and c != "ctrl":
                cs.append(c)
        if not cs:
            continue
        full_effects = []
        for c in cs:
            full_effects.append(np.asarray(X[conditions == c].mean(axis=0)).ravel() - ctrl_present)
        full_effect = np.median(np.stack(full_effects), axis=0)
        rs = []
        for _ in range(B_SPLIT):
            for c in cs:
                idx_c = np.where(conditions == c)[0]
                perm_c = rng.permutation(idx_c)
                ca, cb = perm_c[: len(perm_c) // 2], perm_c[len(perm_c) // 2:]
                if mode == "plus":
                    ctrl_idx = np.where(conditions == "ctrl")[0]
                else:
                    ctrl_idx = np.where(conditions == f"NT@{c.split('@')[1]}")[0]
                perm = rng.permutation(ctrl_idx)
                ctrl_a, ctrl_b = perm[: len(perm) // 2], perm[len(perm) // 2:]
                ea = np.asarray(X[ca].mean(axis=0)).ravel() - np.asarray(X[ctrl_a].mean(axis=0)).ravel()
                eb = np.asarray(X[cb].mean(axis=0)).ravel() - np.asarray(X[ctrl_b].mean(axis=0)).ravel()
                if np.std(ea) == 0 or np.std(eb) == 0:
                    continue
                rs.append(float(pearsonr(ea, eb)[0]))
        rows.append({"gene": g, "effect": full_effect, "reliability": float(np.median(rs) if rs else np.nan)})
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


def evaluate(df, source, dataset):
    X = df[FEATURES].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=int)
    z = base.target_z(X)
    coef = np.asarray(source["coef"])
    j = FEATURES.index("delta_max")
    arms = {
        "source_probe": base.source_scores(X, source),
        "sign_layer": z @ np.sign(coef),
        "coef_target_z": z @ coef,
        "delta_max_sign": np.sign(coef[j]) * X[:, j],
        "in_situ": base.in_situ_scores(X, y),
        "fewshot_30": base.fewshot_scores(X, y, repeats=10),
    }
    out = {}
    for q in ["specific", "shared", "unreliable"]:
        idx = df["quality"].to_numpy() == q
        if idx.sum() < 10 or len(np.unique(y[idx])) < 2:
            continue
        out[q] = {a: float(roc_auc_score(y[idx], s[idx])) for a, s in arms.items()}
    return out


def main(dataset):
    if dataset == "mcf7":
        df, preds, ctrl = base.mcf7_data()
        adata = anndata.read_h5ad(os.path.join(base.BASE, "gse281048_testbed/h5ad/mcf7_gse281048.h5ad"))
        pred = np.load(os.path.join(base.BASE, "gse281048_testbed/outputs/gears_predictions.npz"), allow_pickle=True)
        genes = pred["genes"].tolist()
        gene_names = pred["gene_names"].tolist()
        rel = reliability(adata, genes, gene_names, ctrl, mode="at")
    else:
        pred = np.load(os.path.join(base.BASE, f"outputs_{dataset}/pred_expr_seed0.npz"), allow_pickle=True)
        genes = pred["pred_perturbations"].tolist()
        gene_names = pred["pred_genes"].tolist()
        adata = anndata.read_h5ad(os.path.join(base.BASE, f"data/gears_demo/replogle_{dataset}_essential/perturb_processed.h5ad"))
        ctrl = np.asarray(adata.X[adata.obs["condition"].astype(str) == "ctrl"].mean(axis=0)).ravel().astype(np.float64)
        rel = reliability(adata, genes, gene_names, ctrl, mode="plus")
        X, y, feature_names = load_gears(f"replogle_{dataset}_essential", f"outputs_{dataset}" if dataset == "rpe1" else "outputs")
        df = pd.DataFrame(X, columns=feature_names)
        df["y"] = y
        df["error"] = y.astype(float)
    df = df.merge(rel, on="gene", how="left")
    df["quality"] = df["quality"].fillna("unreliable")
    source = json.load(open(os.path.join(base.BASE, f"outputs_{dataset}/gating_probe_transfer.json" if dataset == "rpe1" else "outputs/gating_probe_transfer.json"), encoding="utf-8"))
    perf = evaluate(df, source, dataset)
    result = {"task": f"Reliability-stratified gate evaluation: {dataset}",
              "n_genes": int(len(df)), "quality_counts": df["quality"].value_counts().to_dict(),
              "performance_by_quality": perf,
              "error_by_quality": df.groupby("quality")["error"].agg(["mean", "median", "count"]).to_dict()}
    prefix = os.path.join(OUT, f"{dataset}_reliability_stratified")
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    df.to_csv(prefix + "_scores.csv", index=False)

    # Plot per-quality AUROC by arm
    arms = ["source_probe", "sign_layer", "coef_target_z", "delta_max_sign", "in_situ", "fewshot_30"]
    fig, ax = plt.subplots(figsize=(11, 4.8))
    qs = [q for i, q in enumerate(["specific", "shared", "unreliable"]) if q in perf]
    x = np.arange(len(arms))
    width = 0.25
    for k, q in enumerate(qs):
        vals = [perf[q][a] for a in arms]
        ax.bar(x + (k - (len(qs)-1)/2) * width, vals, width, label=q)
    ax.axhline(0.5, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(x, arms, rotation=25, ha="right")
    ax.set_ylabel("AUROC")
    ax.set_title(f"{dataset}: gate performance by ground-truth quality")
    ax.legend()
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()

    lines = [f"# Reliability-stratified gate evaluation: {dataset}", "",
             f"Quality counts: {result['quality_counts']}", "",
             "| Quality | Arm | AUROC |", "|---|---|---|"]
    for q in qs:
        for a in arms:
            lines.append(f"| {q} | {a} | {perf[q][a]:.3f} |")
    with open(prefix + ".md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", choices=["mcf7", "k562", "rpe1"])
    args = parser.parse_args()
    main(args.dataset)
