# 结构机制双检验判定（step6f_structure_mechanism）

日期：2026-10-16。数据来源全部为现有归档（`pred_expr_seed0.npz`、`predictions_seed0.csv`、
`perturb_processed.h5ad`），未训练任何模型。一键复现：

```bash
cd E:/清华类脑计算研究中心/virtual_cell_constraint_mvp
PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step6f_structure_mechanism.py
```

口径锁定（与 step1/step6 相同）：误差 = 1 − pearson_delta（top-50 |true_delta| 基因上
Pearson(pred_delta, true_delta)）；标签 y = 误差 > 全体留出扰动误差中位数。
脚本内置自检：从向量重算与 CSV 的 max|diff| ≈ 1e-9（两数据集），口径复现无漂移。

## 检验一：任务内在难度（解释 2）—— 成立

mean predictor = 训练集全部单扰动真值 delta 的逐基因均值（Norman 875 个训练扰动，
RPE1 1237 个），对任何查询预测同一个 delta；误差仍按 top-50 Pearson 口径计算。

| 数据集 | n | Spearman ρ | p | Pearson r | p | 判定（\|ρ\|>0.5） |
|---|---|---|---|---|---|---|
| Norman (K562) | 217 | **0.762** | 2.2e-42 | 0.780 | 1.2e-45 | 难度是结构来源之一 |
| RPE1 | 306 | **0.910** | 5.0e-118 | 0.934 | 4.8e-138 | 难度是结构来源之一 |

- 一个不含任何模型信息、对全部查询输出常数的基线，其逐基因误差与 GEARS 误差
  强相关（尤其 RPE1 达 0.91）：**"哪些扰动基因难"在很大程度上是任务属性，任何模型
  都会在那里出错**。探针探测到的误差结构至少有相当部分来自任务内在难度，
  而非 GEARS 特有的失效模式。
- 相关矩阵一角（探针 9 特征 vs 误差，Spearman）：Norman 各谱特征 |ρ|≈0.18–0.31
  （delta_max −0.312 最强，log_n_cells −0.065 不显著）；RPE1 谱幅度类特征
  |ρ|≈0.50–0.51，log_n_cells +0.426、delta_entropy +0.481。即探针特征与误差是
  中等以下相关——探针并非简单地"重新发现难度"，难度也解释不了全部结构。

**限制**：mean predictor 是单一朴素基线，强相关不等于因果；它也可能部分反映
"真实效应小的扰动所有方法都差"（strength 与难度纠缠）。ρ>0.5 为脚本内约定阈值。

## 检验二：指标 artifact（解释 3）—— 部分成立，需如实标注

误差排序跨指标 Spearman 相关（主口径 = top-50 Pearson）：

| 对比 | Norman ρ (p) | RPE1 ρ (p) |
|---|---|---|
| all-Pearson | 0.821 (2.4e-54) | 0.913 (1.3e-120) |
| all-Spearman | 0.775 (1.1e-44) | 0.876 (3.0e-98) |
| top-50 MSE | **−0.254** (1.6e-04) | **−0.357** (1.3e-10) |

- Pearson 家族内部排序高度一致：**误差结构不是 top-DEG Pearson 独有的 artifact**。
- 但 top-50 MSE 与主口径仅弱负相关——MSE 是幅度耦合指标（未归一化），其误差排序
  确实携带不同信息；凡依赖"误差排序"的结论不应外推到 MSE 类指标。

关键结论复核（Norman strength 三分层探针 AUROC，step6c 协议：5-fold OOF、3 划分种子
均值、层内重拟合；每层 n≈72–73）：

| 误差口径 | low | mid | high |
|---|---|---|---|
| top-50 Pearson（原口径） | **0.391** | **0.301** | 0.703 |
| all-Pearson | **0.332** | 0.741 | 0.814 |
| all-Spearman | **0.470** | 0.727 | 0.862 |
| top-50 MSE | **0.396** | 0.644 | 0.648 |
| mean-baseline 误差（top-50 Pearson） | 0.503 | 0.510 | 0.470 |

- **"low 强度层 AUROC 低于随机"跨 Pearson 家族与 MSE 口径均保持**（0.33–0.47），
  该现象较稳；RPE1 无此现象（各层 0.75–0.90，与 step6c 一致）。
- **"mid 强度层低于随机"只在 top-50 Pearson 口径下出现**，换 all-Pearson /
  all-Spearman / top-50 MSE 即消失（0.64–0.74）。→ **该 half-sentence 必须标注为
  指标敏感，正文只宜引用 low 层现象。**
- 用 mean-baseline 误差做标签时分层 AUROC 全部回到随机附近（0.47–0.51），
  与检验一呼应：分层结构里"模型侧可探"的成分确实压在难度上。

**限制**：分层 AUROC 未给 CI（沿用 step6c 协议）；mid 层结论的敏感性检验是
事后（post-hoc）的，应视为降级而非推翻；单 seed（seed0）预测；MSE 未按基因
数/效应量归一化。

## 总判定

1. 探针误差结构**主要来自任务内在难度**（检验一 ρ=0.76/0.91，超 0.5 判定线）；
2. 误差排序在 Pearson 家族内**稳健**，不是 top-DEG Pearson artifact；但对 MSE 类
   指标不可外推；
3. Norman "low/mid 强度层低于随机"需拆句引用：**low 层跨指标稳健，mid 层仅
   top-50 Pearson 口径成立，属指标敏感结论**。

产物：`outputs{,_rpe1}/structure_mechanism.json`、`outputs/structure_mechanism.png`、
本文；脚本 `step6f_structure_mechanism.py`。
