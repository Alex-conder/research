# Unified external protocol: mcf7

n=99; positive=49; seed predictions=1

| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |
|---|---|---|---|
| source_probe | 0.628 [0.517, 0.742] | 0.137 | -0.0596 [-0.1037, -0.0158] |
| sign_layer | 0.486 [0.384, 0.597] | -0.067 | -0.0079 [-0.0548, 0.0359] |
| coef_target_z | 0.626 [0.515, 0.738] | 0.131 | -0.0588 [-0.1046, -0.0149] |
| delta_max_sign | 0.629 [0.518, 0.743] | 0.138 | -0.0600 [-0.1042, -0.0158] |
| in_situ | 0.630 [0.464, 0.761] | 0.186 | -0.0533 [-0.1070, 0.0119] |
| fewshot_30 | 0.624 [0.451, 0.758] | 0.175 | -0.0505 [-0.1041, 0.0180] |

Key pairwise precision differences:

- sign_layer_minus_source_probe: -0.0517 [-0.1031, 0.0010]
- sign_layer_minus_in_situ: -0.0454 [-0.1167, 0.0305]
- coef_target_z_minus_sign_layer: 0.0509 [0.0025, 0.0988]
- coef_target_z_minus_source_probe: -0.0008 [-0.0095, 0.0081]
- coef_target_z_minus_in_situ: 0.0056 [-0.0450, 0.0530]
- fewshot_30_minus_sign_layer: 0.0426 [-0.0348, 0.1172]
- fewshot_30_minus_in_situ: -0.0028 [-0.0454, 0.0328]
- delta_max_sign_minus_coef_target_z: 0.0012 [-0.0076, 0.0099]
