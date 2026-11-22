# G-B2 定量互补性分析 verdict（step22，prereg_gb2_residual.md）

判定：主分析跨家族有序对 24 对中，AUROC ≥ 0.6 且 95% CI 完全在 0.5 之上的对数 = 0（判据 ≥1）-> **family_structure_is_the_whole_story**。

- CI 完全在 0.5 之上（正向显著）：0 对；完全在下（负向显著，'互反'定量形态）：6 对。
- 负向显著的对按预注册 §4 记为描述性'互反可预测'，不进主判定。

## 残差 AUROC 矩阵（K562，n=217）

行 = 源引擎探针，列 = 目标引擎残差；* = 95% CI 不含 0.5。

| src\\tgt | gears | scgenept | scgpt | geneformer | uce |
|---|---|---|---|---|---|
| gears | - | 0.596* | 0.497 | 0.357* | 0.387* |
| scgenept | 0.487 | - | 0.487 | 0.501 | 0.458 |
| scgpt | 0.448 | 0.536 | - | 0.557 | 0.493 |
| geneformer | 0.515 | 0.511 | 0.513 | - | 0.450 |
| uce | 0.512 | 0.551 | 0.534 | 0.445 | - |

| family pair | pairs | mean AUROC | mean |dev| | sig(+) | sig(-) | pass(locked) |
|---|---|---|---|---|---|---|
| cross_embedding_to_expression | 6 | 0.512 | 0.029 | 0 | 0 | 0 |
| cross_expression_to_embedding | 6 | 0.448 | 0.053 | 0 | 2 | 0 |
| within_embedding | 6 | 0.499 | 0.036 | 0 | 0 | 0 |
| within_expression | 2 | 0.542 | 0.055 | 1 | 0 | 0 |

共享难度占比（留自回归 R²）：gears_k562=0.024；scgenept_k562=0.082；scgpt_k562=0.006；geneformer_k562=0.106；uce_k562=0.090

## 残差 AUROC 矩阵（RPE1，n=306）

行 = 源引擎探针，列 = 目标引擎残差；* = 95% CI 不含 0.5。

| src\\tgt | gears | scgenept | scgpt | geneformer | uce |
|---|---|---|---|---|---|
| gears | - | 0.484 | 0.356* | 0.279* | 0.265* |
| scgenept | 0.497 | - | 0.509 | 0.493 | 0.513 |
| scgpt | 0.409* | 0.545 | - | 0.549 | 0.544 |
| geneformer | 0.440 | 0.482 | 0.486 | - | 0.480 |
| uce | 0.500 | 0.508 | 0.502 | 0.499 | - |

| family pair | pairs | mean AUROC | mean |dev| | sig(+) | sig(-) | pass(locked) |
|---|---|---|---|---|---|---|
| cross_embedding_to_expression | 6 | 0.481 | 0.037 | 0 | 1 | 0 |
| cross_expression_to_embedding | 6 | 0.403 | 0.105 | 0 | 3 | 0 |
| within_embedding | 6 | 0.510 | 0.022 | 0 | 0 | 0 |
| within_expression | 2 | 0.491 | 0.009 | 0 | 0 | 0 |

共享难度占比（留自回归 R²）：gears_rpe1=0.049；scgenept_rpe1=0.023；scgpt_rpe1=0.118；geneformer_rpe1=0.347；uce_rpe1=0.362

## 残差 AUROC 矩阵（scf_k562_diag（诊断，不进判定），n=45）

行 = 源引擎探针，列 = 目标引擎残差；* = 95% CI 不含 0.5。

| src\\tgt | gears | scgenept | scgpt | geneformer | uce | scf |
|---|---|---|---|---|---|---|
| gears | - | 0.516 | 0.530 | 0.425 | 0.367 | 0.536 |
| scgenept | 0.298* | - | 0.490 | 0.554 | 0.607 | 0.284* |
| scgpt | 0.411 | 0.598 | - | 0.582 | 0.537 | 0.494 |
| geneformer | 0.492 | 0.502 | 0.367 | - | 0.431 | 0.516 |
| uce | 0.519 | 0.578 | 0.567 | 0.522 | - | 0.580 |
| scf | 0.482 | 0.488 | 0.605 | 0.508 | 0.514 | - |

| family pair | pairs | mean AUROC | mean |dev| | sig(+) | sig(-) | pass(locked) |
|---|---|---|---|---|---|---|
| cross_embedding_to_expression | 9 | 0.521 | 0.044 | 0 | 0 | 0 |
| cross_expression_to_embedding | 9 | 0.511 | 0.060 | 0 | 0 | 0 |
| within_embedding | 6 | 0.501 | 0.068 | 0 | 0 | 0 |
| within_expression | 6 | 0.434 | 0.083 | 0 | 2 | 0 |

共享难度占比（留自回归 R²）：gears_k562=0.133；scgenept_k562=0.112；scgpt_k562=0.081；geneformer_k562=0.101；uce_k562=0.091；scf_k562=0.157

## 收敛表述（锁定阴性条款触发）

主分析跨家族对无一达 AUROC ≥ 0.6 且 CI 在 0.5 之上：家族结构即全部故事，如实收敛；R2 结论（产品 = 符号层）不变，与 error_signatures 深潜'签名噪声'互证。上表负向显著格（若有）为 kappa 结构'跨家族互反'的定量回声，按预注册为描述性记录。

注：主检验共 24 对（每数据集 12 = 表达谱 2 × embedding 3 × 双向 2），多重性未进判据（任务书锁定 '>=1 对'），显著对数与家族块均值已在上方如实列出。
