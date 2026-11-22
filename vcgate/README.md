# vcgate

Cheap reliability gating for single-cell perturbation-response predictions.

Turns a model's **own** predicted expression spectra into per-query gate scores
("how likely is this prediction to be wrong?"), using model-side features only —
no ground-truth leakage. The probe is an L2 logistic regression under stratified
5-fold CV, exactly the pre-registered protocol of the accompanying case study
(GEARS; AUROC 0.66–0.75 on Norman 2019, 0.80 on Replogle 2022 RPE1; see
`../manuscript/main.tex`, `../preregistration.md`).

## Install

```bash
pip install .          # runtime deps only
pip install .[h5ad,dev]  # + anndata (ctrl mean from h5ad) + pytest
```

## Library use

```python
import numpy as np
from vcgate import compute_features, oof_gate_scores, gate_auroc, median_split_labels

feat = compute_features(pred_expr, perts, pred_genes, ctrl_mean)
feat["log_n_cells"] = np.log1p(n_cells)
y = median_split_labels(1.0 - pearson_delta)      # worst-half labels
scores = oof_gate_scores(feat.values, y)           # high = likely wrong
print(gate_auroc(y, scores))
```

## CLI

```bash
vcgate oof --pred-expr pred_expr_seed0.npz \
    --predictions predictions_seed0.csv \
    --h5ad perturb_processed.h5ad \
    --out gate_scores.csv
```

Input `npz` keys: `pred_expr` (n_pert × n_genes), `pred_perturbations`,
`pred_genes` — the format produced by GEARS runs in the parent repo
(`step1_run.py`).

## Leakage discipline

Accepted inputs: predicted expression, ctrl training mean, cell-count
covariates. Explicitly rejected: any ground-truth-derived column (e.g. a
"strength" computed from held-out true expression) — see
`../step6_gating_probe.py` docstring for the audit trail.

## Tests

```bash
python -m pytest tests/ -v
```
