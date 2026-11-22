# Unified external protocol: mcf7

n=99; positive=49; seed predictions=1

| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |
|---|---|---|---|
| source_probe | 0.523 [0.402, 0.631] | -0.073 | 0.0016 [-0.0461, 0.0489] |
| sign_layer | 0.465 [0.353, 0.582] | -0.196 | 0.0092 [-0.0375, 0.0537] |
| coef_target_z | 0.515 [0.393, 0.626] | -0.100 | 0.0017 [-0.0466, 0.0518] |
| delta_max_sign | 0.523 [0.402, 0.629] | -0.073 | 0.0019 [-0.0459, 0.0489] |
| in_situ | 0.599 [0.438, 0.740] | 0.188 | -0.0350 [-0.0906, 0.0272] |
| fewshot_30 | 0.576 [0.415, 0.709] | 0.152 | -0.0266 [-0.0822, 0.0371] |

Key pairwise precision differences:

- sign_layer_minus_source_probe: -0.0076 [-0.0458, 0.0305]
- sign_layer_minus_in_situ: -0.0442 [-0.1267, 0.0319]
- coef_target_z_minus_sign_layer: 0.0075 [-0.0264, 0.0412]
- coef_target_z_minus_source_probe: -0.0002 [-0.0121, 0.0126]
- coef_target_z_minus_in_situ: -0.0367 [-0.1215, 0.0383]
- fewshot_30_minus_sign_layer: 0.0358 [-0.0407, 0.1179]
- fewshot_30_minus_in_situ: -0.0083 [-0.0554, 0.0301]
- delta_max_sign_minus_coef_target_z: -0.0001 [-0.0136, 0.0126]
