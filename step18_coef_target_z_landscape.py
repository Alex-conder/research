# -*- coding: utf-8 -*-
"""Evaluate target-normalized continuous coefficients across the R2 landscape."""
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from step10_signature_resolution import COMBOS, fit_coef, load_combo
from step10_two_layer_gate import SIBLING, precision_at

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs_rebuild")
COVERAGES = np.round(np.arange(1.0, 0.49, -0.05), 2)
B_BOOT = 1000
BOOT_SEED = 20261026
OOF_SEED = 0


def target_z(X):
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def in_situ(X, y):
    pipe = Pipeline([("scale", StandardScaler()),
                     ("lr", LogisticRegression(C=1.0, max_iter=2000))])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=OOF_SEED)
    return cross_val_predict(pipe, X, y, cv=cv, method="predict_proba")[:, 1]


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def main():
    rows = []
    for key in COMBOS:
        sib = SIBLING.get(key)
        if sib is None:
            continue
        X, y, feature_names = load_combo(key)
        Xs, ys, _ = load_combo(sib)
        coef = fit_coef(Xs, ys)
        signs = np.sign(coef)
        z = target_z(X)
        scores = {"sign_layer": z @ signs, "coef_target_z": z @ coef,
                  "in_situ": in_situ(X, y)}
        point = {a: float(roc_auc_score(y, s)) for a, s in scores.items()}
        prec = {a: np.array([precision_at(s, y, c) for c in COVERAGES]) for a, s in scores.items()}
        rng = np.random.default_rng(BOOT_SEED)
        n = len(y)
        d_sign = np.zeros((B_BOOT, len(COVERAGES)))
        d_ins = np.zeros((B_BOOT, len(COVERAGES)))
        d_rand = np.zeros((B_BOOT, len(COVERAGES)))
        aurocs = {a: [] for a in scores}
        spears = {a: [] for a in scores}
        for b in range(B_BOOT):
            idx = rng.integers(0, n, n)
            while len(np.unique(y[idx])) < 2:
                idx = rng.integers(0, n, n)
            Xb, yb = X[idx], y[idx]
            zb = target_z(Xb)
            sb = {"sign_layer": zb @ signs, "coef_target_z": zb @ coef,
                  "in_situ": in_situ(Xb, yb)}
            for a, s in sb.items():
                aurocs[a].append(roc_auc_score(yb, s))
                spears[a].append(spearmanr(s, yb).statistic)
            pb = {a: np.array([precision_at(s, yb, c) for c in COVERAGES]) for a, s in sb.items()}
            d_sign[b] = pb["coef_target_z"] - pb["sign_layer"]
            d_ins[b] = pb["coef_target_z"] - pb["in_situ"]
            d_rand[b] = float(yb.mean()) - pb["coef_target_z"]
        row = {"combo": key, "n": int(n), "point": point,
               "coef_target_z_minus_sign_layer": {
                   "mean": float(d_sign.mean()), "ci95": ci(d_sign.mean(axis=1))},
               "coef_target_z_minus_in_situ": {
                   "mean": float(d_ins.mean()), "ci95": ci(d_ins.mean(axis=1))},
               "coef_target_z_minus_random": {
                   "mean": float(d_rand.mean()), "ci95": ci(d_rand.mean(axis=1))},
               "bootstrap": {a: {"auroc_mean": float(np.mean(aurocs[a])),
                                 "auroc_ci95": ci(aurocs[a]),
                                 "spearman_mean": float(np.mean(spears[a])),
                                 "spearman_ci95": ci(spears[a])} for a in scores},
               "precision": {a: {f"{c:.2f}": float(v) for c, v in zip(COVERAGES, prec[a])} for a in prec}}
        rows.append(row)
        print(key, point, "coefz-sign", row["coef_target_z_minus_sign_layer"], flush=True)

    result = {"task": "Target-normalized continuous coefficients across R2 landscape",
              "b_boot": B_BOOT, "boot_seed": BOOT_SEED, "combos": rows}
    with open(os.path.join(OUT, "coef_target_z_landscape.json"), "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    # Summary figure: coefz-sign paired difference by combo.
    fig, ax = plt.subplots(figsize=(10, 4.8))
    names = [r["combo"] for r in rows]
    vals = [r["coef_target_z_minus_sign_layer"]["mean"] for r in rows]
    lo = [r["coef_target_z_minus_sign_layer"]["ci95"][0] for r in rows]
    hi = [r["coef_target_z_minus_sign_layer"]["ci95"][1] for r in rows]
    x = np.arange(len(names))
    ax.bar(x, vals, yerr=[np.array(vals)-np.array(lo), np.array(hi)-np.array(vals)],
           color=["#2b8a3e" if lo[i] > 0 else "#7c8aa0" for i in range(len(vals))], capsize=4)
    ax.axhline(0, color="k", linewidth=1)
    ax.set_xticks(x, names, rotation=30, ha="right")
    ax.set_ylabel("coef_target_z - sign_layer\n(coverage-averaged precision difference)")
    ax.set_title("Target-normalized continuous coefficients vs sign layer")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "coef_target_z_landscape.png"), dpi=300)
    plt.close()

    lines = ["# Target-normalized continuous coefficients across R2 landscape", "",
             "| Combo | sign AUROC | coef_z AUROC | coef_z - sign [95% CI] | coef_z - in-situ [95% CI] | coef_z - random [95% CI] |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        b = r["bootstrap"]
        lines.append(f"| {r['combo']} | {b['sign_layer']['auroc_mean']:.3f} | {b['coef_target_z']['auroc_mean']:.3f} | {r['coef_target_z_minus_sign_layer']['mean']:.4f} [{r['coef_target_z_minus_sign_layer']['ci95'][0]:.4f}, {r['coef_target_z_minus_sign_layer']['ci95'][1]:.4f}] | {r['coef_target_z_minus_in_situ']['mean']:.4f} [{r['coef_target_z_minus_in_situ']['ci95'][0]:.4f}, {r['coef_target_z_minus_in_situ']['ci95'][1]:.4f}] | {r['coef_target_z_minus_random']['mean']:.4f} [{r['coef_target_z_minus_random']['ci95'][0]:.4f}, {r['coef_target_z_minus_random']['ci95'][1]:.4f}] |")
    with open(os.path.join(OUT, "coef_target_z_landscape.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("saved", os.path.join(OUT, "coef_target_z_landscape.json"))


if __name__ == "__main__":
    main()
