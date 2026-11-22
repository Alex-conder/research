# Family-Consequence verdict（step24，prereg_family_consequence.md）

三条判定互相独立；全部 on-file 回放，零基础模型训练。

## P1 省钱原则（家族内全引擎集成 vs 单代表门的效用差）

判定：主判定单元 12 个（K562/RPE1 × embedding 家族 × (T,R)），Δ 的 95% CI 含 0 的 12 个 -> **within_family_redundancy_has_cost_consequence**。

| dataset | family | target | rep | mean Δ | 95% CI | 含0 |
|---|---|---|---|---|---|---|
| K562 | expression | gears | scgenept | 0.000 | [0.000,0.000] | 结构性 |
| K562 | expression | scgenept | gears | 0.000 | [0.000,0.000] | 结构性 |
| K562 | embedding | scgpt | geneformer | -0.729 | [-2.247,0.745] | 是 |
| K562 | embedding | scgpt | uce | 0.111 | [-1.893,2.006] | 是 |
| K562 | embedding | geneformer | scgpt | -0.297 | [-2.253,1.112] | 是 |
| K562 | embedding | geneformer | uce | -0.302 | [-1.896,1.436] | 是 |
| K562 | embedding | uce | scgpt | -0.049 | [-0.261,0.116] | 是 |
| K562 | embedding | uce | geneformer | -0.024 | [-0.167,0.119] | 是 |
| RPE1 | expression | gears | scgenept | 0.000 | [0.000,0.000] | 结构性 |
| RPE1 | expression | scgenept | gears | 0.000 | [0.000,0.000] | 结构性 |
| RPE1 | embedding | scgpt | geneformer | -0.129 | [-1.475,1.146] | 是 |
| RPE1 | embedding | scgpt | uce | 0.295 | [-0.891,1.411] | 是 |
| RPE1 | embedding | geneformer | scgpt | -1.067 | [-3.974,1.651] | 是 |
| RPE1 | embedding | geneformer | uce | -1.120 | [-3.710,0.933] | 是 |
| RPE1 | embedding | uce | scgpt | -0.080 | [-0.518,0.376] | 是 |
| RPE1 | embedding | uce | geneformer | -0.117 | [-0.487,0.219] | 是 |
| scf_k562_diag（诊断） | expression | gears | scgenept | -0.442 | [-1.328,0.350] | 是 |
| scf_k562_diag（诊断） | expression | gears | scf | 0.273 | [-0.063,0.665] | 是 |
| scf_k562_diag（诊断） | expression | scgenept | gears | -0.287 | [-1.027,0.427] | 是 |
| scf_k562_diag（诊断） | expression | scgenept | scf | -0.188 | [-0.585,0.163] | 是 |
| scf_k562_diag（诊断） | expression | scf | gears | 0.181 | [-0.424,0.736] | 是 |
| scf_k562_diag（诊断） | expression | scf | scgenept | -0.564 | [-1.334,0.111] | 是 |
| scf_k562_diag（诊断） | embedding | scgpt | geneformer | -0.634 | [-2.000,0.810] | 是 |
| scf_k562_diag（诊断） | embedding | scgpt | uce | 0.817 | [-0.008,1.767] | 是 |
| scf_k562_diag（诊断） | embedding | geneformer | scgpt | 0.060 | [-0.341,0.499] | 是 |
| scf_k562_diag（诊断） | embedding | geneformer | uce | 0.156 | [-0.234,0.719] | 是 |
| scf_k562_diag（诊断） | embedding | uce | scgpt | -0.018 | [-0.109,0.065] | 是 |
| scf_k562_diag（诊断） | embedding | uce | geneformer | 0.042 | [-0.044,0.154] | 是 |

- 两臂绝对效用（对 random 的成本节省，4 档预算均值，非结构性格）：fam_full mean=-1.114，single mean=-0.922。P1 判据只依赖 Δ 的 CI；绝对效用水平用于解读'省钱'含义（两臂皆近 0 时，结论读作'加家族成员不买来门价值'）。

## P2 异质原则（家族分歧度作为第三信号）

判定：AUROC(d, 家族平均误差最差一半) 点估计 — K562=0.432，RPE1=0.537；判据双数据集 ≥ 0.6 -> **family_disagreement_not_predictive**。

| dataset | n | AUROC | 95% CI | sig>0.5 | sig<0.5 | Spearman(d, r̄) | 变体 expr/all | 变体 emb/all |
|---|---|---|---|---|---|---|---|---|
| K562 | 217 | 0.432 | [0.353,0.506] | 否 | 否 | -0.039 | 0.478 | 0.412 |
| RPE1 | 306 | 0.537 | [0.471,0.596] | 否 | 否 | 0.049 | 0.504 | 0.514 |
| scf_k562_diag（诊断） | 45 | 0.336 | [0.191,0.490] | 否 | 是 | -0.227 | 0.564 | 0.306 |

## P3 选择原则（家族代表性的元学习雏形，LODO 回放）

判定：R1 top-1 命中 2/2 折；R1 regret < R0 regret（双折）：否 -> **selection_principle_not_supported**。（明示：2 个全引擎折 + 11 行 on-file 档案，雏形级功效。）

| train -> test | 真实最优 | R1 预测 | R1 hit | R1 regret | R0 预测 | R0 regret | random regret |
|---|---|---|---|---|---|---|---|
| RPE1 -> K562 | GEARS (0.6543) | GEARS | Y | 0.0000 | GEARS | 0.0000 | 0.1242 |
| K562 -> RPE1 | GEARS (0.7897) | GEARS | Y | 0.0000 | GEARS | 0.0000 | 0.2057 |

- scf 描述性：scf_k562(diag) 就地 AUROC=0.6680；K562 训练下 R1 选 GEARS、R0 选 GEARS；RPE1 训练下 R1 选 GEARS、R0 选 GEARS。

## 档案表（P3 输入，全 on-file）

| key | engine | family | unit | n | balance | shared R² | own in-situ AUROC |
|---|---|---|---|---|---|---|---|
| gears_k562 | GEARS | expression | K562 | 217 | 0.498 | 0.133 | 0.6543 |
| gears_rpe1 | GEARS | expression | RPE1 | 306 | 0.500 | 0.049 | 0.7897 |
| scgenept_k562 | scGenePT | expression | K562 | 217 | 0.498 | 0.112 | 0.6016 |
| scgenept_rpe1 | scGenePT | expression | RPE1 | 306 | 0.500 | 0.023 | 0.4593 |
| scf_k562 | scFoundation | expression | scf_k562_diag | 45 | 0.489 | 0.157 | 0.6680 |
| scgpt_k562 | scGPT | embedding | K562 | 217 | 0.498 | 0.081 | 0.4783 |
| scgpt_rpe1 | scGPT | embedding | RPE1 | 306 | 0.500 | 0.118 | 0.5730 |
| geneformer_k562 | Geneformer | embedding | K562 | 217 | 0.498 | 0.101 | 0.4672 |
| geneformer_rpe1 | Geneformer | embedding | RPE1 | 306 | 0.500 | 0.347 | 0.5561 |
| uce_k562 | UCE | embedding | K562 | 217 | 0.498 | 0.091 | 0.4493 |
| uce_rpe1 | UCE | embedding | RPE1 | 306 | 0.500 | 0.362 | 0.5417 |

## 汇总（三条独立）

- P1: within_family_redundancy_has_cost_consequence
- P2: family_disagreement_not_predictive
- P3: selection_principle_not_supported

一句话原则（论文可用句按判定裁剪）：
- P1 过线：'同家族引擎的门错误冗余到可以只部署一台代表——对任一代表选择、任一预算档，单引擎门相对家族全集成门的效用损失 CI 均含 0。'
- P2 过线：'两家族门分数强烈分歧的扰动，家族平均误差确实更高——分歧度本身是可与两家族信号并列的第三信号。'
- P3 过线：'预算只容一台引擎时，家族归属 + 数据集属性可预测就地 AUROC 最高的引擎（LODO 回放）。'
- 任一未过线：对应句子不使用；阴性数字全部保留如上。

与既有档案关系：G-B 定性 kappa -> G-B2 跨家族残差互补阴性 -> 本步回答'家族结构能用来做什么'。三条判定独立，任一不过只禁对应原则句，不动 R2 结论（产品 = 符号层）。
