# Cross-dataset reliability-stratified gate evaluation

AUROC uses dataset-relative worst-half error labels. Intervals are perturbation-level bootstrap 95% CIs.

| Dataset | Quality | n | Arm | AUROC [95% CI] |
|---|---|---:|---|---|
| Replogle K562 | specific | 109 | source_probe | 0.686 [0.578, 0.790] |
| Replogle K562 | specific | 109 | sign_layer | 0.532 [0.423, 0.644] |
| Replogle K562 | specific | 109 | coef_target_z | 0.686 [0.573, 0.785] |
| Replogle K562 | specific | 109 | delta_max_sign | 0.662 [0.553, 0.761] |
| Replogle K562 | specific | 109 | in_situ | 0.654 [0.547, 0.754] |
| Replogle K562 | specific | 109 | fewshot_30 | 0.654 [0.542, 0.759] |
| Replogle K562 | unreliable | 100 | source_probe | 0.588 [0.455, 0.714] |
| Replogle K562 | unreliable | 100 | sign_layer | 0.409 [0.284, 0.548] |
| Replogle K562 | unreliable | 100 | coef_target_z | 0.588 [0.457, 0.707] |
| Replogle K562 | unreliable | 100 | delta_max_sign | 0.563 [0.430, 0.686] |
| Replogle K562 | unreliable | 100 | in_situ | 0.543 [0.415, 0.667] |
| Replogle K562 | unreliable | 100 | fewshot_30 | 0.558 [0.425, 0.676] |
| Replogle RPE1 | specific | 86 | source_probe | 0.966 [0.923, 0.999] |
| Replogle RPE1 | specific | 86 | sign_layer | 0.899 [0.822, 0.965] |
| Replogle RPE1 | specific | 86 | coef_target_z | 0.966 [0.914, 0.998] |
| Replogle RPE1 | specific | 86 | delta_max_sign | 0.924 [0.822, 0.997] |
| Replogle RPE1 | specific | 86 | in_situ | 0.953 [0.888, 0.998] |
| Replogle RPE1 | specific | 86 | fewshot_30 | 0.952 [0.891, 0.995] |
| Replogle RPE1 | shared | 170 | source_probe | 0.795 [0.712, 0.866] |
| Replogle RPE1 | shared | 170 | sign_layer | 0.800 [0.717, 0.873] |
| Replogle RPE1 | shared | 170 | coef_target_z | 0.795 [0.718, 0.868] |
| Replogle RPE1 | shared | 170 | delta_max_sign | 0.773 [0.678, 0.861] |
| Replogle RPE1 | shared | 170 | in_situ | 0.775 [0.695, 0.854] |
| Replogle RPE1 | shared | 170 | fewshot_30 | 0.772 [0.673, 0.857] |
| Replogle RPE1 | unreliable | 50 | source_probe | 0.773 [0.446, 1.000] |
| Replogle RPE1 | unreliable | 50 | sign_layer | 0.782 [0.483, 1.000] |
| Replogle RPE1 | unreliable | 50 | coef_target_z | 0.773 [0.448, 1.000] |
| Replogle RPE1 | unreliable | 50 | delta_max_sign | 0.733 [0.398, 0.959] |
| Replogle RPE1 | unreliable | 50 | in_situ | 0.760 [0.445, 1.000] |
| Replogle RPE1 | unreliable | 50 | fewshot_30 | 0.742 [0.408, 0.962] |
| GSE281048 MCF7 | specific | 14 | source_probe | 0.644 [0.311, 0.917] |
| GSE281048 MCF7 | specific | 14 | sign_layer | 0.667 [0.364, 0.958] |
| GSE281048 MCF7 | specific | 14 | coef_target_z | 0.644 [0.356, 0.917] |
| GSE281048 MCF7 | specific | 14 | delta_max_sign | 0.622 [0.292, 0.919] |
| GSE281048 MCF7 | specific | 14 | in_situ | 0.600 [0.200, 0.950] |
| GSE281048 MCF7 | specific | 14 | fewshot_30 | 0.444 [0.060, 0.848] |
| GSE281048 MCF7 | shared | 46 | source_probe | 0.513 [0.341, 0.671] |
| GSE281048 MCF7 | shared | 46 | sign_layer | 0.515 [0.359, 0.691] |
| GSE281048 MCF7 | shared | 46 | coef_target_z | 0.498 [0.326, 0.665] |
| GSE281048 MCF7 | shared | 46 | delta_max_sign | 0.513 [0.352, 0.685] |
| GSE281048 MCF7 | shared | 46 | in_situ | 0.481 [0.306, 0.650] |
| GSE281048 MCF7 | shared | 46 | fewshot_30 | 0.284 [0.138, 0.441] |
| GSE281048 MCF7 | unreliable | 39 | source_probe | 0.514 [0.325, 0.704] |
| GSE281048 MCF7 | unreliable | 39 | sign_layer | 0.364 [0.188, 0.556] |
| GSE281048 MCF7 | unreliable | 39 | coef_target_z | 0.492 [0.300, 0.689] |
| GSE281048 MCF7 | unreliable | 39 | delta_max_sign | 0.511 [0.329, 0.696] |
| GSE281048 MCF7 | unreliable | 39 | in_situ | 0.726 [0.550, 0.886] |
| GSE281048 MCF7 | unreliable | 39 | fewshot_30 | 0.503 [0.321, 0.683] |

| Dataset | Quality | target-z minus sign [95% CI] | delta-max minus sign [95% CI] |
|---|---|---|---|
| Replogle K562 | specific | +0.1541 [+0.0272, +0.2939] | +0.1301 [-0.0390, +0.3067] |
| Replogle K562 | unreliable | +0.1796 [+0.0652, +0.3100] | +0.1543 [-0.0553, +0.3700] |
| Replogle RPE1 | specific | +0.0667 [+0.0178, +0.1364] | +0.0242 [-0.0780, +0.1165] |
| Replogle RPE1 | shared | -0.0046 [-0.0629, +0.0591] | -0.0265 [-0.1114, +0.0646] |
| Replogle RPE1 | unreliable | -0.0089 [-0.1088, +0.1064] | -0.0489 [-0.2068, +0.1460] |
| GSE281048 MCF7 | specific | -0.0222 [-0.2292, +0.1500] | -0.0444 [-0.3333, +0.2224] |
| GSE281048 MCF7 | shared | -0.0170 [-0.1418, +0.1078] | -0.0019 [-0.1424, +0.1437] |
| GSE281048 MCF7 | unreliable | +0.1277 [-0.0030, +0.2790] | +0.1467 [-0.0237, +0.3159] |
