# Reliability-stratified gate evaluation: k562

Quality counts: {'specific': 109, 'unreliable': 100, 'shared': 8}

Error-reliability Spearman: -0.533

| Quality | Arm | AUROC [95% CI] | n |
|---|---|---|---|
| specific | source_probe | 0.686 [0.578, 0.790] | 109 |
| specific | sign_layer | 0.532 [0.423, 0.644] | 109 |
| specific | coef_target_z | 0.686 [0.573, 0.785] | 109 |
| specific | delta_max_sign | 0.662 [0.553, 0.761] | 109 |
| specific | in_situ | 0.654 [0.547, 0.754] | 109 |
| specific | fewshot_30 | 0.654 [0.542, 0.759] | 109 |
| unreliable | source_probe | 0.588 [0.455, 0.714] | 100 |
| unreliable | sign_layer | 0.409 [0.284, 0.548] | 100 |
| unreliable | coef_target_z | 0.588 [0.457, 0.707] | 100 |
| unreliable | delta_max_sign | 0.563 [0.430, 0.686] | 100 |
| unreliable | in_situ | 0.543 [0.415, 0.667] | 100 |
| unreliable | fewshot_30 | 0.558 [0.425, 0.676] | 100 |

| Quality | target-z minus sign [95% CI] | delta-max minus sign [95% CI] |
|---|---|---|
| specific | +0.1541 [+0.0272, +0.2939] | +0.1301 [-0.0390, +0.3067] |
| unreliable | +0.1796 [+0.0652, +0.3100] | +0.1543 [-0.0553, +0.3700] |
