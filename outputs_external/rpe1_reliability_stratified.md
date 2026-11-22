# Reliability-stratified gate evaluation: rpe1

Quality counts: {'shared': 170, 'specific': 86, 'unreliable': 50}

Error-reliability Spearman: -0.505

| Quality | Arm | AUROC [95% CI] | n |
|---|---|---|---|
| specific | source_probe | 0.966 [0.923, 0.999] | 86 |
| specific | sign_layer | 0.899 [0.822, 0.965] | 86 |
| specific | coef_target_z | 0.966 [0.914, 0.998] | 86 |
| specific | delta_max_sign | 0.924 [0.822, 0.997] | 86 |
| specific | in_situ | 0.953 [0.888, 0.998] | 86 |
| specific | fewshot_30 | 0.952 [0.891, 0.995] | 86 |
| shared | source_probe | 0.795 [0.712, 0.866] | 170 |
| shared | sign_layer | 0.800 [0.717, 0.873] | 170 |
| shared | coef_target_z | 0.795 [0.718, 0.868] | 170 |
| shared | delta_max_sign | 0.773 [0.678, 0.861] | 170 |
| shared | in_situ | 0.775 [0.695, 0.854] | 170 |
| shared | fewshot_30 | 0.772 [0.673, 0.857] | 170 |
| unreliable | source_probe | 0.773 [0.446, 1.000] | 50 |
| unreliable | sign_layer | 0.782 [0.483, 1.000] | 50 |
| unreliable | coef_target_z | 0.773 [0.448, 1.000] | 50 |
| unreliable | delta_max_sign | 0.733 [0.398, 0.959] | 50 |
| unreliable | in_situ | 0.760 [0.445, 1.000] | 50 |
| unreliable | fewshot_30 | 0.742 [0.408, 0.962] | 50 |

| Quality | target-z minus sign [95% CI] | delta-max minus sign [95% CI] |
|---|---|---|
| specific | +0.0667 [+0.0178, +0.1364] | +0.0242 [-0.0780, +0.1165] |
| shared | -0.0046 [-0.0629, +0.0591] | -0.0265 [-0.1114, +0.0646] |
| unreliable | -0.0089 [-0.1088, +0.1064] | -0.0489 [-0.2068, +0.1460] |
