# Reliability-stratified gate evaluation: mcf7

Quality counts: {'shared': 46, 'unreliable': 39, 'specific': 14}

Error-reliability Spearman: 0.169

| Quality | Arm | AUROC [95% CI] | n |
|---|---|---|---|
| specific | source_probe | 0.644 [0.311, 0.917] | 14 |
| specific | sign_layer | 0.667 [0.364, 0.958] | 14 |
| specific | coef_target_z | 0.644 [0.356, 0.917] | 14 |
| specific | delta_max_sign | 0.622 [0.292, 0.919] | 14 |
| specific | in_situ | 0.600 [0.200, 0.950] | 14 |
| specific | fewshot_30 | 0.444 [0.060, 0.848] | 14 |
| shared | source_probe | 0.513 [0.341, 0.671] | 46 |
| shared | sign_layer | 0.515 [0.359, 0.691] | 46 |
| shared | coef_target_z | 0.498 [0.326, 0.665] | 46 |
| shared | delta_max_sign | 0.513 [0.352, 0.685] | 46 |
| shared | in_situ | 0.481 [0.306, 0.650] | 46 |
| shared | fewshot_30 | 0.284 [0.138, 0.441] | 46 |
| unreliable | source_probe | 0.514 [0.325, 0.704] | 39 |
| unreliable | sign_layer | 0.364 [0.188, 0.556] | 39 |
| unreliable | coef_target_z | 0.492 [0.300, 0.689] | 39 |
| unreliable | delta_max_sign | 0.511 [0.329, 0.696] | 39 |
| unreliable | in_situ | 0.726 [0.550, 0.886] | 39 |
| unreliable | fewshot_30 | 0.503 [0.321, 0.683] | 39 |

| Quality | target-z minus sign [95% CI] | delta-max minus sign [95% CI] |
|---|---|---|
| specific | -0.0222 [-0.2292, +0.1500] | -0.0444 [-0.3333, +0.2224] |
| shared | -0.0170 [-0.1418, +0.1078] | -0.0019 [-0.1424, +0.1437] |
| unreliable | +0.1277 [-0.0030, +0.2790] | +0.1467 [-0.0237, +0.3159] |
