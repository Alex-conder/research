# Step 7 判定：跨模型逐扰动路由（routing）——误差探针的"画像预测力"边界

日期：2026-10-16 ｜ 分析：`step7_routing.py` ｜ 产物：`outputs_crossmodel/analysis/routing.json` + `routing.png`

## 结论（一句话）

**在当前设计下，误差探针不能逐扰动可信地路由 GEARS 与 scGenePT 的预测：路由组显著劣于最优单模型（全 GEARS），p=1.000（bootstrap，扰动为抽样单位，3 seeds × 1000 次）。** 这不是探针失效——两个探针各自在模型内预测误差都有效（GEARS 测试 AUROC≈0.70，scGenePT≈0.58）——而是**跨模型概率不可比**：两个就地探针各自回答"在本模型的误差分布里这是不是后一半"，它们的 0.5 阈值对应完全不同的绝对误差水平，直接比较两探针输出的"预测误差概率"会把大量扰动错误路由到整体差一个量级的模型。本结果构成"可靠性画像预测力"的边界标定：画像信号是**模型内相对排序**信号，不是**跨模型绝对校准**信号。

## 四组指标（mean top-50 pearson_delta，30% held-out 路由测试集，3 划分种子 均值±std）

| 组 | mean ± std | per-seed |
|---|---|---|
| 全 GEARS（最优单模型） | **0.553 ± 0.035** | 0.591 / 0.548 / 0.521 |
| 全 scGenePT | 0.257 ± 0.020 | 0.265 / 0.271 / 0.234 |
| 探针路由 | 0.405 ± 0.030 | 0.439 / 0.385 / 0.391 |
| Oracle（逐扰动取真实更优，上界） | 0.568 ± 0.031 | 0.599 / 0.569 / 0.537 |

- 路由 − 最优单模型 = **−0.148 ± 0.017**（per-seed：−0.152 / −0.164 / −0.130）
- 配对 bootstrap（扰动为抽样单位，1000 次/种子，pooled 3000 次）：差值分布 100% 落在 0 以下，**单侧 p = 1.000**（即"路由不劣于最优单模型"这一零假设被完全拒绝的反方向：路由以确定性幅度劣于全 GEARS）
- 路由 − oracle = **−0.164 ± 0.019**：路由离上界的差距比 oracle 本身相对全 GEARS 的可提升空间（+0.015）大一个数量级
- oracle 仅需把 9.7% 的扰动路由给 scGenePT（oracle_to_scgenept_frac=0.097），而探针实际把 **57.9%**（per-seed 53.8% / 64.6% / 55.4%）路由给了 scGenePT —— 路由方向系统性错误

## 机制诊断

1. **探针本身有效（模型内）**：GEARS 就地探针在 30% 测试上预测"本模型后一半误差"AUROC=0.745/0.670/0.694（mean 0.703）；scGenePT 探针 AUROC=0.553/0.574/0.608（mean 0.578）。与 step6/step6f 的结论一致。
2. **跨模型不可比（核心失败点）**：标签 y = 误差 > 训练集中位数，中位数在两个模型间差近一倍（scGenePT 整体 pearson≈0.254 vs GEARS≈0.567，且 138/217 扰动的 target 不在 scGenePT 输入空间、预测≈ctrl，误差系统性高）。因此"scGenePT 探针说这不是后一半"远弱于"GEARS 探针说这不是后一半"，但路由规则把两者当成了同一量纲的分数。
3. **上界本身很窄**：即使 oracle 也只能比全 GEARS 提升 +0.015（scGenePT 在全部 217 扰动中仅 16 个真实更优）。在"一个模型全面占优"的数据集上，路由价值的先验上限就极低——本实验同时量化了这一上限。

## 边界标定（对"可靠性画像预测力"的含义）

- 已确认的能力：误差探针是**模型内**的相对可靠性排序器（AUROC 0.58–0.70，跨种子稳定），可用于"同一模型输出中先信哪些扰动"（precision-coverage，step6b 已验证）。
- 标定出的边界：**逐扰动跨模型路由需要跨模型可比的误差预测（绝对校准或共同参照）**，当前的二分类中位数探针不具备这一性质；直接把两探针分数相比较会因模型间基线误差水平不同而系统性偏向弱模型。
- 可行的后续方向（超出本次范围，仅记录）：(a) 用共同锚点校准（如对两模型预测同一扰动的谱特征差做回归而非各自中位数二分类）；(b) 路由规则改为"仅当强模型探针判'高风险'时才切换到弱模型"（单向门控），避免双向误路由；(c) 在模型水平相当（无全面占优）的模型对（如 GEARS 多种子/多架构）上重测，oracle 上限才可能张开。

## 口径与数据来源说明

- **GEARS 预测：从归档恢复，未重算。** 使用 `outputs/pred_expr_seed0.npz`（217 扰动 × 5000 基因，split_seed=0 留出名单），即 step6 探针的同一 artifact；标签用 `outputs/predictions_seed0.csv` 的 `pearson_delta`（top-50 |delta| Pearson，step1_run.py:160-162 口径）。特征经 `step6d_naive_baselines.load_dataset` 计算，与 step6 完全一致。
- **scGenePT**：`outputs_crossmodel/scgenept/scgenept_k562_predictions.npz` + `scgenept_k562_labels.csv`，经 `step6f_crossmodel.load_crossmodel` 对齐（4813 基因匹配空间；无 NaN 行剔除）。`target_in_input_space=False` 的 138 个扰动预测≈ctrl 为数据事实，保留原样。
- **探针规格**：StandardScaler + L2 logistic C=1.0（max_iter=2000），在 70% 训练扰动上就地对各模型特征分别拟合；标签 y = (1 − pearson_delta) > 训练集中位数（step6 口径）。
- **划分**：扰动级 70/30，划分种子 [0,1,2]（测试集各 65 扰动）；bootstrap 主种子 20261016。
- 显著性：配对 bootstrap，扰动为抽样单位，每种子 1000 次，pooled 3000 次，单侧（H1：路由 > 最优单模型）。

## 复现

```bash
cd "E:/清华类脑计算研究中心/virtual_cell_constraint_mvp"
PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7_routing.py --dry-run   # 产物检查
PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7_routing.py             # 全分析
```

产物：`outputs_crossmodel/analysis/routing.json`（全部数字与 per-seed 明细）、`outputs_crossmodel/analysis/routing.png`（300dpi，四组柱状图 + bootstrap 分布）。

---

# Step 7b 判定：共同锚点校准后的路由复测——"复活"为中性，未产生增益

日期：2026-10-16 ｜ 分析：`step7b_routing_calibrated.py` ｜ 产物：`outputs_crossmodel/analysis/routing_calibrated.json` + `routing_calibrated.png`

## 结论（一句话）

**校准成功消除了 step7 的有害误路由（路由给 scGenePT 的比例从 57.9% → 0.0%），但复活后的路由与全 GEARS 逐点相同（增益恒等于 0），未能回收 oracle 上界（−0.015）。** 判定：**半复活——校准修复了机制层面的跨模型不可比，但模型侧 9 特征不含识别 scGenePT 局部优胜扰动（16/217）的信号，路由价值仍为 0。** 结合 step7，对"画像预测力"的完整边界标定：特征能预测**模型内相对可靠性**（精度覆盖、门控有效），也能校准出**跨模型绝对水平差**（正确判断"scGenePT 全面更差"），但预测不了**跨模型的局部反差点**。

## 标签口径核验（先验证再分析）

- scGenePT 的 matched 4813 基因空间 = full 5000 空间中非 NaN 的 4813 列（逐元素一致），其归档标签在统一口径下重算 **pearson 相关 = 1.000**（mean |diff| ≈ 2e-9），口径完全一致。
- GEARS 归档标签在 5000 基因空间计算，统一口径（4813 公共空间）重算后 pearson 相关 = 0.986、mean |diff| = 0.025（top-50 基因集随空间略变，预期内的小幅口径差）。
- 因此本实验对**两模型统一在 4813 基因公共空间、同一 h5ad 真值重算** top-50 pearson_delta；全量均值 GEARS=0.5771、scGenePT=0.2542。217 扰动无一缺失真值。

## 六组指标（mean top-50 pearson_delta，30% held-out，3 seeds 均值±std；与 step7 完全相同的划分）

| 组 | mean ± std |
|---|---|
| 全 GEARS | 0.561 ± 0.038 |
| 全 scGenePT | 0.257 ± 0.020 |
| 未校准探针路由（step7 基线复算） | 0.410 ± 0.027 |
| Ridge 校准路由 | 0.561 ± 0.038 |
| Isotonic 校准路由 | 0.561 ± 0.038 |
| Oracle（上界） | 0.576 ± 0.031 |

- 校准路由 − 全 GEARS = **0.000 ± 0.000**（两种校准器、全部 3 种子逐点等于全 GEARS）；bootstrap 差值恒为 0，p 退化（无信息）——正确表述是"零增益"，而非"不显著"。
- 校准路由 − oracle = **−0.015 ± 0.007**（oracle 上界未被回收）。
- **路由比例（机制验证关键证据）**：未校准探针路由给 scGenePT = 56.9%（53.8/63.1/53.8），校准后 Ridge/isotonic = **0.0%**，oracle = 8.7%（6.2/10.8/9.2）。校准把路由方向从系统性错误修正为保守正确，但修正方式是"永不信任 scGenePT"，而非"识别出那 8.7% 该信任的场景"。

## 校准器方案与泄漏防护

- Ridge：StandardScaler + Ridge(alpha=1.0)，9 特征 → 绝对 pearson_delta 连续回归；70% 训练侧拟合、30% 测试侧预测（标准流程，无泄漏）。
- Isotonic：Ridge 在训练侧内部做 5-fold OOF 预测（shuffle，fold seed = 划分种子），`IsotonicRegression(out_of_bounds="clip")` 拟合 (OOF 预测 → 训练真值)；测试侧先 Ridge 再过 isotonic。OOF 输入保证单调校准器不接触训练点的 in-sample 拟合值（若用 in-sample 预测拟合 isotonic 会高估校准置信度）。
- 两种校准器结果逐点一致（isotonic 校准未改变任何路由决策），说明结果对校准器选择稳健，失败不在校准器而在特征信息量。

## 机制解释

校准器学到了"scGenePT 的绝对预测质量全面低于 GEARS"这一真实水平差（在全部 3 个测试集上 65/65 个扰动都预测 GEARS 更优），这正是 step7 缺失的共同锚点——因此有害误路由被根除。但 scGenePT 的 16 个真实优胜扰动（oracle 仅路由 8.7%）在 9 个模型侧特征（预测谱统计 + log_n_cells）上没有可辨识的签名：这些扰动上 scGenePT 的预测谱看起来"同样差"， Ridge/isotonic 都无法把它们排到 GEARS 前面。换言之，**跨模型局部反转点的预测需要模型侧特征所不包含的信息**（可能需要模型内部量、基因集先验或组合扰动结构知识）。

## 复现

```bash
cd "E:/清华类脑计算研究中心/virtual_cell_constraint_mvp"
PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7b_routing_calibrated.py --dry-run
PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7b_routing_calibrated.py
```

产物：`outputs_crossmodel/analysis/routing_calibrated.json`（全部数字、per-seed、口径核验）、`outputs_crossmodel/analysis/routing_calibrated.png`（300dpi，六组柱状图 + 双校准器 bootstrap 分布）。
