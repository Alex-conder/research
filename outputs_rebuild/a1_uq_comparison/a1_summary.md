# A1 UQ Comparison Summary

- K562 | sign-layer (zero-cost): AUROC = 0.620, cost = 1x inference
- K562 | magnitude (delta_max): AUROC = 0.640, cost = 1x inference
- K562 | ensemble disagreement: AUROC = 0.589, cost = 3x inference
- RPE1 | sign-layer (zero-cost): AUROC = 0.578, cost = 1x inference
- RPE1 | magnitude (delta_max): AUROC = 0.759, cost = 1x inference
- RPE1 | ensemble disagreement: AUROC = 0.570, cost = 2x inference
- K562 | conformal prediction: AUROC = 0.678, cost = 1x + calibration set
- Norman (ref) | probe_full: AUROC = 0.654, cost = ref
- Norman (ref) | delta_max_single: AUROC = 0.623, cost = ref
- Norman (ref) | magnitude_l2: AUROC = 0.384, cost = ref
- Norman (ref) | ensemble_variance: AUROC = 0.427, cost = ref

## Key message
Sign-layer is zero-cost but modest AUROC. Ensemble disagreement requires multiple inference passes but can improve AUROC when multiple independent engines are available. Conformal prediction requires a calibration set, breaking the zero-parameter deployment promise.
