# Unified external protocol: norman

n=102; positive=51; seed predictions=3

| Arm | AUROC [95% CI] | Spearman | risk minus random [95% CI] |
|---|---|---|---|
| source_probe | 0.547 [0.432, 0.651] | 0.072 | -0.0218 [-0.0674, 0.0237] |
| sign_layer | 0.479 [0.372, 0.581] | -0.081 | -0.0058 [-0.0503, 0.0377] |
| coef_target_z | 0.581 [0.472, 0.679] | 0.117 | -0.0362 [-0.0798, 0.0076] |
| delta_max_sign | 0.629 [0.519, 0.737] | 0.173 | -0.0627 [-0.1066, -0.0178] |
| in_situ | 0.705 [0.576, 0.810] | 0.388 | -0.0803 [-0.1228, -0.0285] |
| fewshot_30 | 0.695 [0.556, 0.800] | 0.368 | -0.0759 [-0.1194, -0.0213] |
| ensemble_disagreement | 0.430 [0.315, 0.544] | -0.096 | 0.0212 [-0.0248, 0.0660] |

Key pairwise precision differences:

- sign_layer_minus_source_probe: -0.0160 [-0.0491, 0.0164]
- sign_layer_minus_in_situ: -0.0745 [-0.1408, -0.0040]
- coef_target_z_minus_sign_layer: 0.0305 [0.0022, 0.0594]
- coef_target_z_minus_source_probe: 0.0145 [-0.0052, 0.0376]
- coef_target_z_minus_in_situ: -0.0441 [-0.1033, 0.0170]
- fewshot_30_minus_sign_layer: 0.0701 [-0.0065, 0.1367]
- fewshot_30_minus_in_situ: -0.0044 [-0.0311, 0.0191]
- delta_max_sign_minus_coef_target_z: 0.0265 [-0.0059, 0.0604]
