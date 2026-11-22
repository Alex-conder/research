# 红队挂账实算：S4 闭环配对检验 + S3 预测缩放对照

- 日期：2026-10-26；任务：redteam S4 paired test + S3 scaling control (step13, real computation)
- 红队来源：redteam_20261024/gating.md S3/S4；门产品口径：outputs_rebuild/two_layer_gate.json verdict tier=symbol_layer
- 实算脚本：`step13_s4_s3_fixes.py`；本报告由 `step13b_s4_s3_report.py` 生成
- S4 随机源：bootstrap 主种子 20261103；随机拒答基线种子 [910101, 910102, 910103]（点估计曲线每 种子 x coverage 抽签 1000 次；配对 bootstrap 内每重采样每种子逐 coverage 各抽 1 次）；OOF 划分种子 0（与 step6/step10 同源）
- S3 随机源：step6 原协议种子 0（CV 划分与 LR 初始化）
- 口径复用：组合加载/三门分数来自 step10_signature_resolution 与 step10_two_layer_gate；三门点估计已与 outputs_rebuild/two_layer_gate.json 对账一致（最大差 0.0）

## S4 闭环配对检验：门 vs 随机拒答基线

### 方法

对 R2 判定产品的门（符号层臂，另报就地探针臂与两层门臂作对照），在 coverage 网格 1.0 -> 0.5（步长 0.05，11 档）上构造随机拒答基线：每档 coverage 以等大概率拒答（标记）k=round(coverage x n) 个查询，risk = 1 - 被标记查询中最差一半的占比。点估计基线曲线取 3 个随机种子、每种子每档 1000 次抽签的平均。
配对 bootstrap（B=1000，扰动级有放回重采样）：同一重采样内重算三门分数，随机拒答臂按 3 个固定种子各抽一次取平均（主口径 mc3，差值的公共基线噪声在对差中抵消）；稳健性复核用期望随机臂（重采样基率，即随机拒答的期望精度，蒙特卡洛噪声为零）。每次重采样逐档记录 risk(门) - risk(随机)，报逐档差值均值与 95% percentile CI，以及 coverage 平均差值的 95% CI。判定：CI 不含 0 且方向为负（门 risk 显著低于随机）。scf_k562（仅就地臂）与 scf_rpe1（产物缺失）按 R2 修正案 A2 留痕排除，可执行组合 = 10。

### 表 S4-1 门产品臂（符号层，R2 判定 tier=symbol_layer），mc3 配对口径

| 组合 | 引擎 | 数据集 | n | 门 risk@cov0.50 | 随机 risk@cov0.50 | 配对差均值(coverage 1.0-0.5 平均) | 95% CI | 门显著更优 |
|---|---|---|---|---|---|---|---|---|
| gears k562 | GEARS | K562 | 217 | 0.4074 | 0.5019 | -0.0532 | [-0.0877, -0.0183] | 是 |
| gears rpe1 | GEARS | RPE1 | 306 | 0.3922 | 0.4999 | -0.0588 | [-0.0850, -0.0309] | 是 |
| scgenept k562 | scGenePT | K562 | 217 | 0.4907 | 0.5019 | -0.0235 | [-0.0576, +0.0160] | 否 |
| scgenept rpe1 | scGenePT | RPE1 | 306 | 0.4837 | 0.4992 | -0.0061 | [-0.0368, +0.0254] | 否 |
| scgpt k562 | scGPT | K562 | 217 | 0.5370 | 0.5019 | +0.0092 | [-0.0257, +0.0452] | 否 |
| scgpt rpe1 | scGPT | RPE1 | 306 | 0.4967 | 0.4996 | -0.0072 | [-0.0369, +0.0222] | 否 |
| geneformer k562 | Geneformer | K562 | 217 | 0.4815 | 0.5024 | +0.0067 | [-0.0294, +0.0441] | 否 |
| geneformer rpe1 | Geneformer | RPE1 | 306 | 0.5033 | 0.4998 | +0.0204 | [-0.0118, +0.0492] | 否 |
| uce k562 | UCE | K562 | 217 | 0.5000 | 0.5020 | +0.0092 | [-0.0267, +0.0456] | 否 |
| uce rpe1 | UCE | RPE1 | 306 | 0.5294 | 0.4999 | +0.0172 | [-0.0134, +0.0479] | 否 |

判定（该臂，mc3 口径）：2/10 组合门的 risk 显著低于随机（CI 不含 0 且方向为负）。

### 表 S4-2 对照臂：就地探针，mc3

| 组合 | 引擎 | 数据集 | n | 门 risk@cov0.50 | 随机 risk@cov0.50 | 配对差均值(coverage 1.0-0.5 平均) | 95% CI | 门显著更优 |
|---|---|---|---|---|---|---|---|---|
| gears k562 | GEARS | K562 | 217 | 0.3981 | 0.5019 | -0.0700 | [-0.1053, -0.0292] | 是 |
| gears rpe1 | GEARS | RPE1 | 306 | 0.2876 | 0.4999 | -0.1139 | [-0.1384, -0.0887] | 是 |
| scgenept k562 | scGenePT | K562 | 217 | 0.4537 | 0.5019 | -0.0322 | [-0.0732, +0.0154] | 否 |
| scgenept rpe1 | scGenePT | RPE1 | 306 | 0.5490 | 0.4992 | -0.0071 | [-0.0478, +0.0363] | 否 |
| scgpt k562 | scGPT | K562 | 217 | 0.5093 | 0.5019 | +0.0008 | [-0.0335, +0.0394] | 否 |
| scgpt rpe1 | scGPT | RPE1 | 306 | 0.4510 | 0.4996 | -0.0211 | [-0.0529, +0.0129] | 否 |
| geneformer k562 | Geneformer | K562 | 217 | 0.5185 | 0.5024 | +0.0045 | [-0.0337, +0.0429] | 否 |
| geneformer rpe1 | Geneformer | RPE1 | 306 | 0.4510 | 0.4998 | -0.0118 | [-0.0477, +0.0262] | 否 |
| uce k562 | UCE | K562 | 217 | 0.5463 | 0.5020 | +0.0047 | [-0.0287, +0.0427] | 否 |
| uce rpe1 | UCE | RPE1 | 306 | 0.4837 | 0.4999 | -0.0143 | [-0.0494, +0.0236] | 否 |

判定（该臂，mc3 口径）：2/10 组合门的 risk 显著低于随机（CI 不含 0 且方向为负）。

### 表 S4-3 对照臂：两层门，mc3

| 组合 | 引擎 | 数据集 | n | 门 risk@cov0.50 | 随机 risk@cov0.50 | 配对差均值(coverage 1.0-0.5 平均) | 95% CI | 门显著更优 |
|---|---|---|---|---|---|---|---|---|
| gears k562 | GEARS | K562 | 217 | 0.3889 | 0.5019 | -0.0361 | [-0.0756, +0.0034] | 否 |
| gears rpe1 | GEARS | RPE1 | 306 | 0.3725 | 0.4999 | -0.0125 | [-0.0428, +0.0175] | 否 |
| scgenept k562 | scGenePT | K562 | 217 | 0.5278 | 0.5019 | +0.0037 | [-0.0346, +0.0440] | 否 |
| scgenept rpe1 | scGenePT | RPE1 | 306 | 0.4967 | 0.4992 | -0.0077 | [-0.0406, +0.0252] | 否 |
| scgpt k562 | scGPT | K562 | 217 | 0.5185 | 0.5019 | +0.0049 | [-0.0308, +0.0404] | 否 |
| scgpt rpe1 | scGPT | RPE1 | 306 | 0.4967 | 0.4996 | -0.0027 | [-0.0388, +0.0317] | 否 |
| geneformer k562 | Geneformer | K562 | 217 | 0.5093 | 0.5024 | +0.0074 | [-0.0284, +0.0434] | 否 |
| geneformer rpe1 | Geneformer | RPE1 | 306 | 0.5294 | 0.4998 | +0.0130 | [-0.0184, +0.0427] | 否 |
| uce k562 | UCE | K562 | 217 | 0.5278 | 0.5020 | +0.0104 | [-0.0251, +0.0475] | 否 |
| uce rpe1 | UCE | RPE1 | 306 | 0.5163 | 0.4999 | +0.0152 | [-0.0123, +0.0418] | 否 |

判定（该臂，mc3 口径）：0/10 组合门的 risk 显著低于随机（CI 不含 0 且方向为负）。

### 稳健性复核（期望随机臂配对，产品臂 symbol）

| 组合 | 配对差均值 | 95% CI | 门显著更优 | 与 mc3 判定一致 |
|---|---|---|---|---|
| gears_k562 | -0.0533 | [-0.0843, -0.0225] | 是 | 是 |
| gears_rpe1 | -0.0592 | [-0.0824, -0.0363] | 是 | 是 |
| scgenept_k562 | -0.0233 | [-0.0547, +0.0097] | 否 | 是 |
| scgenept_rpe1 | -0.0063 | [-0.0325, +0.0200] | 否 | 是 |
| scgpt_k562 | +0.0092 | [-0.0198, +0.0397] | 否 | 是 |
| scgpt_rpe1 | -0.0073 | [-0.0332, +0.0183] | 否 | 是 |
| geneformer_k562 | +0.0067 | [-0.0235, +0.0383] | 否 | 是 |
| geneformer_rpe1 | +0.0204 | [-0.0064, +0.0463] | 否 | 是 |
| uce_k562 | +0.0096 | [-0.0234, +0.0418] | 否 | 是 |
| uce_rpe1 | +0.0172 | [-0.0094, +0.0434] | 否 | 是 |

期望随机臂口径判定：2/10 组合显著（与 mc3 主口径计数一致：是）。

### 附：产品臂逐档差值（risk_gate - risk_random，均值 [95% CI]，mc3）

| 组合 | 1.00 | 0.95 | 0.90 | 0.85 | 0.80 | 0.75 | 0.70 | 0.65 | 0.60 | 0.55 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gears_k562 | +0.000 [+0.000,+0.000]* | -0.016 [-0.031,+0.002]* | -0.017 [-0.043,+0.010]* | -0.023 [-0.056,+0.009]* | -0.037 [-0.073,+0.000]* | -0.056 [-0.098,-0.008] | -0.077 [-0.127,-0.028] | -0.091 [-0.144,-0.035] | -0.093 [-0.156,-0.026] | -0.091 [-0.163,-0.017] | -0.085 [-0.164,+0.000]* |
| gears_rpe1 | +0.000 [+0.000,+0.000]* | -0.020 [-0.031,-0.008] | -0.043 [-0.061,-0.025] | -0.062 [-0.086,-0.035] | -0.060 [-0.091,-0.029] | -0.060 [-0.094,-0.023] | -0.068 [-0.109,-0.030] | -0.072 [-0.119,-0.025] | -0.075 [-0.125,-0.024] | -0.082 [-0.141,-0.024] | -0.105 [-0.168,-0.039] |
| scgenept_k562 | +0.000 [+0.000,+0.000]* | -0.013 [-0.029,+0.006]* | -0.033 [-0.055,-0.009] | -0.046 [-0.074,-0.014] | -0.049 [-0.086,-0.010] | -0.041 [-0.090,+0.008]* | -0.027 [-0.079,+0.024]* | -0.018 [-0.073,+0.043]* | -0.015 [-0.079,+0.051]* | -0.011 [-0.081,+0.064]* | -0.004 [-0.083,+0.080]* |
| scgenept_rpe1 | +0.000 [+0.000,+0.000]* | -0.002 [-0.017,+0.013]* | +0.002 [-0.021,+0.023]* | +0.010 [-0.019,+0.037]* | +0.010 [-0.026,+0.044]* | -0.001 [-0.042,+0.042]* | -0.013 [-0.058,+0.031]* | -0.018 [-0.064,+0.032]* | -0.015 [-0.065,+0.038]* | -0.016 [-0.079,+0.044]* | -0.023 [-0.089,+0.046]* |
| scgpt_k562 | +0.000 [+0.000,+0.000]* | -0.002 [-0.019,+0.015]* | +0.001 [-0.026,+0.027]* | +0.006 [-0.029,+0.036]* | +0.011 [-0.027,+0.050]* | +0.016 [-0.031,+0.061]* | +0.014 [-0.037,+0.064]* | +0.012 [-0.043,+0.066]* | +0.013 [-0.049,+0.072]* | +0.014 [-0.053,+0.084]* | +0.015 [-0.065,+0.093]* |
| scgpt_rpe1 | +0.000 [+0.000,+0.000]* | +0.004 [-0.013,+0.019]* | +0.008 [-0.015,+0.029]* | +0.001 [-0.026,+0.027]* | -0.003 [-0.034,+0.029]* | -0.007 [-0.043,+0.029]* | -0.019 [-0.062,+0.026]* | -0.026 [-0.075,+0.022]* | -0.023 [-0.078,+0.031]* | -0.008 [-0.069,+0.054]* | -0.005 [-0.070,+0.065]* |
| geneformer_k562 | +0.000 [+0.000,+0.000]* | +0.008 [-0.011,+0.026]* | +0.006 [-0.024,+0.032]* | +0.001 [-0.034,+0.036]* | -0.001 [-0.040,+0.038]* | +0.009 [-0.037,+0.057]* | +0.014 [-0.037,+0.064]* | +0.011 [-0.045,+0.073]* | +0.009 [-0.054,+0.074]* | +0.009 [-0.059,+0.081]* | +0.008 [-0.071,+0.086]* |
| geneformer_rpe1 | +0.000 [+0.000,+0.000]* | +0.018 [+0.006,+0.030] | +0.029 [+0.007,+0.049] | +0.033 [+0.006,+0.059] | +0.030 [-0.003,+0.061]* | +0.031 [-0.010,+0.067]* | +0.025 [-0.017,+0.070]* | +0.022 [-0.023,+0.067]* | +0.020 [-0.036,+0.072]* | +0.013 [-0.050,+0.071]* | +0.003 [-0.063,+0.065]* |
| uce_k562 | +0.000 [+0.000,+0.000]* | +0.001 [-0.018,+0.019]* | +0.007 [-0.022,+0.032]* | +0.007 [-0.027,+0.040]* | +0.008 [-0.029,+0.046]* | +0.009 [-0.033,+0.053]* | +0.015 [-0.035,+0.064]* | +0.016 [-0.040,+0.071]* | +0.014 [-0.051,+0.077]* | +0.013 [-0.059,+0.081]* | +0.011 [-0.068,+0.086]* |
| uce_rpe1 | +0.000 [+0.000,+0.000]* | +0.018 [+0.002,+0.031] | +0.017 [-0.005,+0.036]* | +0.017 [-0.010,+0.042]* | +0.017 [-0.015,+0.048]* | +0.018 [-0.020,+0.055]* | +0.019 [-0.025,+0.064]* | +0.019 [-0.030,+0.069]* | +0.020 [-0.034,+0.072]* | +0.021 [-0.038,+0.077]* | +0.023 [-0.044,+0.085]* |

（* = 该档 CI 含 0；coverage=1.0 时门与随机标记全体，差恒为 0。）
随机拒答基线的种子稳定性：逐组合 3 种子基线 risk 最大极差见 JSON s4.combos[*].random_rejection_baseline.mean_over_seeds.risk_spread_over_seeds（量级 1e-3 或更小）。

## S3 预测缩放对照：探针 AUROC 的缩放稳健性

### 方法

主对照：对 GEARS-K562/RPE1 的 9 特征探针输入做统一缩放（StandardScaler 全数据一次 fit），以同一 CV 划分（stratified 5-fold, shuffle, seed=0）重训纯 L2 logistic 并重评 AUROC，与 step6 原口径（fold 内 StandardScaler 管线，CV AUROC 均值为主读数、pooled OOF AUROC 为辅读数）对比。重训前已对原口径做复现校验，与存档 outputs{,_rpe1}/gating_probe.json 一致（差 < 1e-9）。
补充对照（红队原文'恒定/缩放预测'）：把预测 delta 乘常数因子 0.5 与 2.0 重算特征（n_de_025/frac_de_025 经 0.25 阈值非线性变化，其余线性或不变），pearson 标签对正缩放数学不变（y 不动），按原协议重评。另报完全不缩放的诊断行。

### 结果

| 组合 | 口径 | CV AUROC 均值 | pooled OOF AUROC | 相对原口径差 |
|---|---|---|---|---|
| gears_k562 | 原口径（fold 内缩放） | 0.655647 | 0.654689 | - |
| gears_k562 | 统一缩放对照（全数据一次 fit） | 0.656080 | 0.654349 | +0.000433 |
| gears_k562 | 诊断：完全不缩放 | 0.628788 | 0.623853 | -0.026860 |
| gears_k562 | 补充：预测 delta x0.5 | 0.656553 | 0.653500 | +0.000905 |
| gears_k562 | 补充：预测 delta x2.0 | 0.656513 | 0.652990 | +0.000866 |
| gears_rpe1 | 原口径（fold 内缩放） | 0.796518 | 0.789440 | - |
| gears_rpe1 | 统一缩放对照（全数据一次 fit） | 0.796518 | 0.789654 | +0.000000 |
| gears_rpe1 | 诊断：完全不缩放 | 0.803184 | 0.795036 | +0.006667 |
| gears_rpe1 | 补充：预测 delta x0.5 | 0.796934 | 0.788543 | +0.000416 |
| gears_rpe1 | 补充：预测 delta x2.0 | 0.797364 | 0.788714 | +0.000846 |

### 判定

- gears_k562：统一缩放前后 CV AUROC 差 +0.000433（pooled 差 -0.000340）；预测级缩放 delta x0.5 / x2.0 的 AUROC 差分别为 +0.000905 / +0.000866；原口径复现与存档一致：True。
- gears_rpe1：统一缩放前后 CV AUROC 差 +0.000000（pooled 差 +0.000214）；预测级缩放 delta x0.5 / x2.0 的 AUROC 差分别为 +0.000416 / +0.000846；原口径复现与存档一致：True。

结论：探针 AUROC 不是缩放步骤的产物——把缩放射到 CV 口径外（全数据统一缩放）AUROC 不动；把预测本身缩放常数倍，AUROC 的移动在 1e-3 量级以内（见表），不改变任何判定。

## 附：论文附录英文草稿（主代理决定入稿位置；未改 main.tex）

### S4 段（Formal paired test against random rejection）

```
\textbf{Formal paired test against random rejection.} 
We replaced the descriptive gate-versus-random delta (Sec.~
\textit{closed-loop}) with a paired test on the risk--coverage 
grid. For each of the 10 evaluable gate combinations (R2 protocol; 
scFoundation K562 lacks the sibling product and the RPE1 product 
is missing), a random-rejection baseline flags a uniformly random 
k = round(coverage x n) queries at each of 11 coverage points 
(1.0 to 0.5, step 0.05); its risk is averaged over 3 seeds x 1000 
draws (seed-to-seed spread <= 0.0028). A paired perturbation-level 
bootstrap (B=1000; both arms recomputed on each resample, random 
arm redrawn per resample under the 3 fixed seeds; master seed 
20261103) tests the coverage-averaged risk difference per 
combination. The product gate (symbol layer) significantly beats 
random rejection in 2/10 combinations -- both GEARS 
expression-class combinations: K562 -0.0532 [-0.0877, -0.0183], 
RPE1 -0.0588 [-0.0850, -0.0309] (95% percentile CI) -- and is 
statistically indistinguishable from random rejection on the 
eight embedding-class wave combinations. The in-situ probe arm 
shows the same pattern (2/10; GEARS only), and the two-layer arm 
reaches significance in no combination (0/10). A zero-
Monte-Carlo-noise pairing (expected-random arm = resample base rate) reproduces the product-gate (symbol-arm) verdict in all ten combinations; for the two-layer arm the same pairing flips only the knife-edge GEARS-K562 call (CI upper +0.0034 vs -0.0002), leaving its count at 0--1/10 either way. The gate's 
closed-loop value is therefore confirmed on the expression class 
and bounded elsewhere; we report the non-significant combinations 
rather than averaging them into a cleaner story.
```

### S3 段（Prediction-scaling control）

```
\textbf{Prediction-scaling control.} 
To exclude the possibility that the probe AUROC is an artifact of 
the scaling stage, we retrained the 9-feature GEARS probe after a 
unified StandardScaler fit (full data, once) and re-evaluated 
under the identical 5-fold CV protocol (seed 0): the CV AUROC 
mean moves by +0.000433 on K562 (0.655647 -> 0.656080) and by 
0.000000 on RPE1 (0.796518 -> 0.796518); pooled OOF AUROC moves 
by at most 0.00034. Scaling the predicted perturbation deltas by 
constant factors 0.5 and 2.0, recomputing all probe features (the 
DE-count features change non-linearly through the 0.25 threshold; 
Pearson labels are invariant to positive rescaling), shifts AUROC 
by at most +0.0009 (K562) and +0.0008 (RPE1). The probe's signal 
is not a product of either feature scaling or prediction scale.
```

