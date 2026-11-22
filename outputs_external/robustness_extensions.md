# Robustness extensions

## Error-caliber robustness

| Dataset | Error caliber | AUROC [95% CI] | permutation p |
|---|---|---|---|
| Replogle K562 | top50_true_pearson | 0.655 [0.580, 0.727] | 0.0020 |
| Replogle K562 | top100_true_pearson | 0.655 [0.577, 0.724] | 0.0020 |
| Replogle K562 | all_gene_pearson | 0.743 [0.671, 0.809] | 0.0020 |
| Replogle K562 | all_gene_spearman | 0.773 [0.706, 0.832] | 0.0020 |
| Replogle RPE1 | top50_true_pearson | 0.789 [0.737, 0.840] | 0.0020 |
| Replogle RPE1 | top100_true_pearson | 0.798 [0.745, 0.845] | 0.0020 |
| Replogle RPE1 | all_gene_pearson | 0.856 [0.815, 0.894] | 0.0020 |
| Replogle RPE1 | all_gene_spearman | 0.858 [0.818, 0.895] | 0.0020 |

## Magnitude-free ablation

| Dataset | Arm | AUROC [95% CI] | permutation p |
|---|---|---|---|
| Replogle K562 | full_probe | 0.655 [0.578, 0.724] | 0.0020 |
| Replogle K562 | magnitude_free | 0.524 [0.447, 0.602] | 0.2954 |
| Replogle RPE1 | full_probe | 0.789 [0.735, 0.840] | 0.0020 |
| Replogle RPE1 | magnitude_free | 0.776 [0.727, 0.827] | 0.0020 |
