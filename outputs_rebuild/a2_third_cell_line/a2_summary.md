# A2 Third Cell-Line Validation (MCF7 GSE281048)

- MCF7 conditions: 189
- In-situ calibrated gate AUROC on MCF7: 0.581 (95% CI: [0.4651123037169549, 0.702322287350555])
- 3-fold CV AUROC on MCF7: 0.600
- K562 sign-layer threshold applied to MCF7: 0.3768729610671435
- Fraction of MCF7 conditions flagged by K562 threshold: 1.000

## Interpretation
The in-situ gate on MCF7 achieves AUROC ~0.58-0.60, comparable to K562/RPE1. This weakly supports cross-cell-line generalization, though the effect size is modest and the 95% CI includes chance. A true zero-parameter transfer AUROC requires per-condition ground-truth error labels, which are not stored in the MCF7 prediction NPZ.
