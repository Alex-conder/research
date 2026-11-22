# External MCF7 validation plan v1 (internally dated before the new run)

- Date: 2026-10-29
- Status: analysis plan recorded before running `step14_external_mcf7_validation.py`.
- Purpose: replace the withdrawn threshold-only MCF7/TTA claim with a valid external stress test.

## Data and inclusion

- Dataset: GSE281048 MCF7, pathway-targeted Perturb-seq with 24 h cytokine stimulation.
- Primary analysis: all unique genes for which GEARS predictions and the same top-50 |delta| Pearson error are available (expected n=99). No confidence panel, no calibration-set exclusion, no selection by gate score or error.
- Labels: `error = 1 - median_pearson_r` across pathways; primary label is the dataset-relative median split. Continuous error and an absolute `r < 0.5` label are sensitivity analyses.
- Features: the same nine model-side features on the identical 5000-gene space as K562.

## Arms

1. `source_probe`: archived K562 scaler + coefficients + intercept, transferred without refitting.
2. `sign_layer`: signs of the archived K562 coefficients, dot product with target-side unlabeled z-scored features (zero target-label fitting).
3. `delta_max_sign`: sign of the K562 `delta_max` coefficient times raw MCF7 `delta_max` (the original single-feature sign rule).
4. `in_situ`: target-side 5-fold OOF L2 logistic probe with fold-internal standardization (uses MCF7 labels; comparator, not zero-label deployment).
5. `random`: expected random-rejection baseline at each coverage.

## Primary endpoints and interpretation

- AUROC for every arm, with perturbation-level bootstrap 95% intervals.
- Coverage-averaged risk difference `sign_layer - random` and `in_situ - random` over coverage 1.0 to 0.5 (11 points), paired bootstrap B=1000; a gate is favored only if the nominal 95% CI excludes zero in the harmful-for-random direction.
- Pairwise `sign_layer - source_probe` and `sign_layer - in_situ` differences. Sign superiority is NOT required and will not be claimed unless the CI supports it.

## Secondary/sensitivity analyses

- Spearman correlation of each score with continuous error.
- Absolute-label sensitivity: `r < 0.5` rather than median split.
- Condition-level (gene x pathway) descriptive AUROC by pathway to expose stimulation/pathway dependence; primary inference remains at unique-gene level.

## Reproducibility

- Bootstrap master seed: 20261029. OOF split seed: 0. B=1000.
- No model retraining on K562/RPE1, no new wet-lab data, no use of the old 60-gene panel or 15-gene calibration set.
- The MCF7 dataset is stimulus-confounded and is therefore an external stress test, not proof of general deployment across all cell lines.
