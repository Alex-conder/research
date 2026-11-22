"""vcgate — cheap reliability gating for perturbation-response predictions.

Turn a model's own predicted expression spectra into per-query gate scores:
"how likely is this prediction to be wrong?" Features are model-side only
(no ground-truth leakage); the probe is an L2 logistic regression with
stratified 5-fold CV, exactly as pre-registered in the accompanying paper
(see ../manuscript/main.tex and ../preregistration.md).

典型用法：
    from vcgate import compute_features, oof_gate_scores, gate_auroc

    feat = compute_features(pred_expr, perts, pred_genes, ctrl_mean)
    scores = oof_gate_scores(feat, error_above_median_labels)
    auroc = gate_auroc(labels, scores)

CLI：
    python -m vcgate oof --pred-expr pred_expr_seed0.npz \
        --predictions predictions_seed0.csv \
        --h5ad perturb_processed.h5ad --out gate_scores.csv
"""
from vcgate.features import COVARIATES, SPECTRUM_FEATURES, compute_features
from vcgate.probe import (GATE_CRITERION, cv_gate_auroc, fit_gate, gate_auroc,
                          median_split_labels, oof_gate_scores)

__version__ = "0.1.0"
__all__ = [
    "compute_features", "SPECTRUM_FEATURES", "COVARIATES",
    "oof_gate_scores", "gate_auroc", "fit_gate", "cv_gate_auroc",
    "median_split_labels", "GATE_CRITERION",
]
