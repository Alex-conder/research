# External MCF7 validation (GSE281048)

Primary set: all 99 eligible unique genes; median-split positives=49; absolute r<0.5 positives=98.

| Arm | AUROC (bootstrap mean [95% CI]) | Spearman(error) | sign-minus-random risk (coverage mean [95% CI]) |
|---|---|---|---|
| K562 probe transfer | 0.632 [0.521, 0.738] | 0.145 | -0.0607 [-0.1036, -0.0176] |
| sign layer | 0.482 [0.369, 0.595] | -0.078 | -0.0069 [-0.0538, 0.0394] |
| delta-max sign | 0.633 [0.522, 0.739] | 0.146 | -0.0612 [-0.1042, -0.0188] |
| MCF7 in-situ | 0.635 [0.465, 0.760] | 0.199 | -0.0553 [-0.1039, 0.0069] |

Pairwise coverage-averaged precision differences:

- sign_layer_minus_source_probe: -0.0538 [-0.1068, -0.0033]
- sign_layer_minus_in_situ: -0.0484 [-0.1238, 0.0259]
- source_probe_minus_in_situ: 0.0054 [-0.0415, 0.0495]

Diagnostic: source-probe and delta-max-sign score correlation = 0.999995; source delta_max scaler mean = 0.495, MCF7 delta_max median = 195.24. The transferred full probe is therefore effectively dominated by delta_max on MCF7.

The absolute r<0.5 label has 98/99 positives and is therefore not used for interpretation.

Caveat: GSE281048 MCF7 contains 24 h cytokine stimulation and is an external stress test, not a same-protocol Replogle validation. The in-situ arm uses MCF7 labels. Intervals are nominal and unadjusted for multiplicity.
