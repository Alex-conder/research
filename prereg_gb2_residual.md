# 预注册：G-B2 —— 错误家族结构的定量互补性分析（残差 AUROC）

- Date: 2026-10-31
- Status: recorded **before** running `step22_gb2_residual.py`。纯 on-file 分析：
  不训练任何新基础模型；唯一拟合对象是 (a) 家族既定规格的浅层 logistic
  探针（G-B 已算，直接复用其 OOF 分数），(b) 本步的误差秩 OLS 残差提取
  （p ≤ 5 + 截距，n = 217/306，in-sample 拟合——见 §3 记录）。
- 上位依据：任务书 G-B2；G-B 结论（`outputs_rebuild/ensemble_gate_verdict.md`，
  step17）：两两 kappa 呈口径家族结构（家族内错误冗余、跨家族错误互反），
  **定性**。本步把该结构发现升级为定量互补性分析。组合注册表、标签口径、
  探针规格复用 `prereg_expansion_gb.md` / `step10_signature_resolution.py`，
  不另造。

## 1. 数据与家族划分（锁定）

- 误差连续量（rank 变换的原料，同锁定二值标签的出处）：
  - 表达谱组合：`err = 1 − pearson_delta`（step6d / step6f 口径）。
  - embedding 组合：`err = err_emb_dist`（wave 标签 CSV，step17 口径）。
- **rank 口径（锁定）**：组合内秩百分位 `r_i = rank(err_i)/n`（升序，
  误差越大秩越高）。作用：跨引擎尺度无关（pearson 距离与 embedding 距离
  单位不可比），且与家族"最差一半 = 秩上半"的口径一致。
- 探针分数：复用 G-B 基础块——各引擎就地探针 5-fold OOF 概率
  （StratifiedKFold seed=0，fold 内 StandardScaler + L2 logistic C=1.0；
  由 `step10_two_layer_gate.oof_probe_scores` 重新生成，规格同 step17）。
  方向约定：分数越高 = 该引擎自身组合越可能落最差一半（与 G-B AUROC
  计算方向一致）。
- 家族（锁定）：表达谱 = {GEARS, scGenePT, scFoundation}；
  embedding = {scGPT, Geneformer, UCE}。
- 分析集：
  - **主分析**：K562 全覆盖 5 引擎（gears / scgenept / scgpt / geneformer
    / uce，n=217，同 G-B 基础块，scf 不入主分析——172 行无分数不插补）；
    RPE1 全覆盖 5 引擎（n=306）。scf_rpe1 产物缺失（R2 起已登记）。
  - **诊断**：scf_k562（n=45 子集，6 引擎按扰动名对齐），与 G-B 同处理
    （诊断不进锁定判定）。

## 2. 共享难度成分提取（锁定）

对每个有序对 (A→B)（A ≠ B，同数据集）：

1. 取数据集 D 的公共扰动集上，B 的误差秩 `r_B` 与全部引擎的 OOF 探针分数。
2. **OLS（含截距）**：`r_B ~ {probe_E : E ∈ engines(D), E ≠ A}`。
3. 残差 `resid_B^{(−A)}` = 拟合残差，即"扣掉其他所有引擎（含 B 自身）线性
   可解释成分后的 B 特异误差"。

**口径解释（锁定，任务书"所有引擎"的操作化）**：任务书原文为"对所有引擎的
探针分数做回归"。若把被检验的源引擎 A 放入回归元，则 OLS 正交性使
`corr(probe_A, resid_B) = 0` 按构造成立，互补性检验退化为只捡二值化非线性
残差，检验 vacuous。因此锁定操作化为**逐对留源（leave-source-out）**：
对有序对 (A→B)，回归元 = 除 A 外全部引擎（含 B 自身探针——最严格的
"增量信息"表述：A 相对包括 B 自己在内的所有其他引擎的增量）。此口径与
任务书判据的意图（"A 对 B 的残差误差有增量预测信息"）一致。

- 描述性分解（同步产出，不判定）：每组合"共享难度占比" = 留自回归
  （`r_B ~ 其他引擎探针`，不含 B 自身）的 R²；`1 − R²` = 引擎特异误差占比。

## 3. 统计方法与防泄漏（锁定）

- 残差提取为 in-sample OLS：p ≤ 5 + 截距，n ≥ 45，过拟合空间可忽略；
  逐组合记录 R² 备查。该选择在判定前公开记录。
- 检验：对 `resid_B^{(−A)}` 取二值 `z = 1[resid > 0]`（残差为正 = 比所有人
  预期更差），AUROC(probe_A, z)。
- 显著性：扰动级 bootstrap（B=1000，master seed **20261031**，逐对由主 seed
  派生 Generator），95% percentile CI；**CI 不含 0.5** 判显著。只重采样
  扰动索引，残差与分数向量固定（同 G-B 记录的理由）。
- 方向：AUROC 按"分数越高 → 残差为正概率越高"原向计算；`AUROC < 0.5` 的
  显著格 = 负向可预测（A 高分处 B 相对更好），即 G-B kappa 结构里"跨家族
  互反"的定量形态。矩阵全符号如实报告。

## 4. 判据（锁定，主判定）

- **主判定（scope：主分析两数据集的跨家族有序对，n=12×2=24 对）**：
  存在 ≥1 对跨家族组合（表达谱→embedding 或 embedding→表达谱，K562 或
  RPE1）残差 AUROC ≥ 0.6 且显著（95% CI 完全在 0.5 之上）
  → **"跨家族互补"定量成立**，G-B 节升级为"集成受家族结构限制但互补信息
  存在（残差证据）"。
  （更正记录，运行前：初稿写"n=18×2=36 对"，系把 scf 误并入主分析表达谱
  引擎数；主分析每数据集表达谱仅 GEARS/scGenePT 两台，scf_k562 按 §1 为
  诊断集。判据 ">=1 对"不变，此处仅更正计数。）
- **否则（全无）**：判**"家族结构即全部故事"**，如实收敛，不改写 R2 结论
  （产品 = 符号层）。
- 多重性记录（描述性，不改判定）：24 对主检验，报告显著对数与家族块
  平均偏差 |AUROC−0.5|；判据 ">=1 对"为任务书锁定阈值，照原文执行。
- 次解读（预注册，描述性）：负向显著（CI 完全在 0.5 之下）的跨家族格
  记为"互反可预测"（kappa 结构的定量回声）；家族内格（冗余方向）同矩阵
  报告。两者均不进主判定。

## 5. 产物与复跑

- `outputs_rebuild/gb2_residual.json`（全 spec、种子、逐对残差 AUROC/CI、
  家族块汇总、共享难度分解、判定）
- `outputs_rebuild/gb2_residual.md`（判定留痕 + 残差 AUROC 矩阵）
- `outputs_rebuild/gb2_residual.png`（逐数据集残差 AUROC 方阵热图，300dpi）
- `outputs_rebuild/smoke/gb2_residual.smoke.{json,png}`（--smoke：B=100、
  仅 K562；smoke 数字不进判据）
- 脚本：`step22_gb2_residual.py`；环境：本机 venv，CPU；UTF-8 纪律
  （不打印 emoji）。

## 6. 与既有档案的关系

- G-B（step17）kappa 矩阵给家族结构的定性描述；G-B2 用"扣掉共享难度后
  还能否互相预测"定量回答"跨家族互补信息是否存在"。
- 判定任一支均保留全矩阵与全部数字，不删对、不换阈值、不做事后重定义；
  阴性结果与 R2 结论（产品 = 符号层）兼容，与 error_signatures 深潜
  "签名噪声"互证关系在 verdict 中注明。
