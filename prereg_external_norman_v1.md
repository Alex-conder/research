# Norman 2019 external protocol and robust-transfer plan v1

- Date: 2026-10-29
- Status: recorded before running `step16_external_protocol.py norman` on the new predictions.
- Routes: A = reusable reliability evaluation protocol; B = target-normalized continuous-coefficient transfer.

## Data and inclusion

- Norman 2019 A549, single-gene perturbations only (conditions with exactly one target gene plus ctrl).
- Use every predicted unique gene with a computable top-50 |delta| Pearson error; no confidence panel or gate-score selection.
- Predictions come from the main K562 GEARS models (seeds 0/1/2). Seed 0 is the primary predictor; seeds 1/2 support ensemble disagreement.

## Arms

- `source_probe`: archived K562 scaler and coefficients, no refit.
- `sign_layer`: target-side unlabeled z-scored features dot the coefficient signs.
- `coef_target_z`: target-side unlabeled z-scored features dot the full continuous coefficients (route B).
- `delta_max_sign`: single-feature sign rule.
- `in_situ`: target-side 5-fold OOF logistic probe with fold-internal scaling.
- `fewshot_30`: repeated 30-label calibration probe, out-of-calibration averaged scores (uses target labels; comparator).
- `ensemble_disagreement`: across-seed predicted-delta variability (when available).

## Primary endpoints

- AUROC and continuous-error Spearman correlation, with perturbation-level bootstrap intervals.
- Coverage-averaged risk difference versus expected random rejection over coverage 1.0 to 0.5.
- Paired differences: `coef_target_z - sign_layer`, `coef_target_z - source_probe`, `coef_target_z - in_situ`, `fewshot_30 - sign_layer`, and `delta_max_sign - coef_target_z`.

## Decision criteria

- Route B is supported if `coef_target_z - sign_layer` has a positive coverage-averaged precision difference with nominal 95% CI excluding zero, and `coef_target_z` remains above chance. It is NOT required to beat in-situ or few-shot calibration.
- If sign layer fails but continuous coefficients recover, report this as evidence that coefficient magnitudes carry information beyond signs; do not claim universal deployment.
- All eligible genes are reported; negative and pathway/seed-sensitive results are retained.
- Bootstrap B=1000, master seed 20261029; OOF split seed 0; few-shot calibration size 30.
