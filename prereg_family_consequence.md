# 预注册：Family-Consequence —— "错误按架构家族组织"的三条有后果设计原则

- Date: 2026-10-30
- Status: recorded **before** running `step24_family_consequence.py`。纯 on-file
  回放：不训练任何新基础模型；唯一拟合对象与 G-B/G-B2/G-C 相同——家族既定
  规格的浅层 logistic 探针（5-fold OOF，StratifiedKFold seed=0，fold 内
  StandardScaler + L2 C=1.0，由 `step10_two_layer_gate.oof_probe_scores`
  确定性重算），与既有各步同一函数同一数据。
- 上位依据：任务书（gating 论文第二跳）；G-B（`outputs_rebuild/
  ensemble_gate_verdict.md`，step17）kappa 家族结构的**定性**观察
  （家族内错误冗余：embedding 三引擎 RPE1 最强 kappa +0.88；跨家族互反：
  GEARS×UCE −0.46/−0.63）；G-B2（`outputs_rebuild/gb2_residual.json`，
  step22）：跨家族残差互补为阴性（family_structure_is_the_whole_story），
  共享难度 R² 分层（表达谱误差高度引擎特异 0.02–0.08）。本步把家族结构从
  观察升级为**有实操后果的设计原则**，三条判定互相独立。
- 组合注册表、误差口径、家族划分、探针规格全部复用
  `prereg_gb2_residual.md` / `step10_signature_resolution.py`，不另造。

## 1. 数据与家族划分（锁定）

- 误差连续量：表达谱组合 `err = 1 − pearson_delta`（step6d/step6f 口径）；
  embedding 组合 `err = err_emb_dist`（wave 标签 CSV，step17 口径）。
- 误差秩：`r = rank(err)/n`（组合内升序百分位，与 G-B2 同口径）。
- 探针分数：复用 G-B/G-B2 基础块（OOF 概率；方向：越高 = 该引擎自身组合
  越可能落最差一半）。
- 家族（锁定）：表达谱 = {GEARS, scGenePT, scFoundation}；
  embedding = {scGPT, Geneformer, UCE}。
- 分析集（锁定）：
  - **主分析**：K562 全覆盖 5 引擎（n=217）、RPE1 全覆盖 5 引擎（n=306）。
  - **诊断**：scf_k562（n=45 子集，6 引擎按扰动名对齐），与 G-B2 同处理
    （诊断不进锁定判定）。

## 2. P1 —— 省钱原则（家族内全引擎集成 vs 家族代表单引擎的效用差）

**部署语义（锁定）**：门系统在家族 F 内可部署全部引擎的探针（全引擎集成门）
或只部署一台代表引擎的探针（单引擎门）；被门控的消费目标固定为家族内某
目标引擎 T 的预测误差 `err_T`（同一误差向量上比较两门臂，单位天然一致；
c·(n−B) 项两臂相同，效用差只取决于 top-B 采信的误差和）。

- 对数据集 D、家族 F、目标 T ∈ F：
  - 臂 `fam_full`：分数 = F∖{T} 全部成员探针的均值（家族对其余部分的
    集成门；**排除 T 自身探针**——检验的是"家族冗余信息"而非"T 自己的门
     trivially 最好"，避免判据被自门优势污染）。
  - 臂 `single_R`：分数 = R 的探针，R ∈ F∖{T}（任一单代表选择）。
- 效用（G-C 既定口径，任务书指定"回放效用差"）：预算 B ∈ {5, 10, 20, 50}
  （G-C 锁定网格），`U(s; B) = (k/n)·Σ err_T − topk_sum(err_T, s, B)`
  （相对 random 采信的成本节省，k = min(B, n)）。
- 效用损失：`Δ(T,R,B) = C(single_R) − C(fam_full)
  = topk_sum(err_T, s_R, B) − topk_sum(err_T, s_full, B)`；
  Δ > 0 = 单代表更差。
- 统计（锁定）：扰动级配对 bootstrap（两臂在同一重采样索引上重算），
  B=1000，master seed **20261101**（逐格由主 seed 派生 Generator），
  95% percentile CI。
- **主判定格（锁定）**：embedding 家族（3 成员，F∖{T} 为 2 台，非退化），
  K562 与 RPE1 两数据集：每数据集 (T,R) = 3×2 = 6 格 × 4 个预算档。
  **主判定单元 = (D, T, R)**：Δ 取 4 个预算档的均值（每个 bootstrap 重采样
  内同步重算 4 档再取均值，保持配对），共 **12 个主判定单元**。
- **判据 P1（锁定）**：12 个主判定单元的 Δ 95% CI **全部含 0** →
  "部署家族内任一代表单引擎相对家族内全引擎集成无显著门效用损失"，
  家族内冗余的实操推论（省钱原则）成立；任一单元 CI 完全在 0 之上 → P1
  不成立（如实报告该格）。48 个 (T,R,B) 细格全部列出；每 (D,F) 的大均值
  Δ 与 CI 作头条描述量。
- 结构性记录（锁定）：表达谱家族在主分析每数据集仅 2 成员（scf 不入主
  分析，同 G-B2 先例），F∖{T} 为单台，`fam_full ≡ single_R`，Δ≡0 由构造
  成立，记为**结构性平凡**、不作信息性判读；表达谱 3 成员的非退化检验由
  诊断集 scf_k562（{gears, scgenept, scf}）承担，结果描述性报告。
- 描述性（不判定）：每格两臂各自对 random 的效用 U（绝对效用语境——若
  两臂绝对效用均近 0，"省钱"结论的解读为"加成员不买来门价值"，verdict 中
  必须注明两臂绝对效用水平）。

## 3. P2 —— 异质原则（家族分歧作为第三信号）

**动机（锁定）**：G-B2 阴性 = 跨家族无正向增量互补信息；但 G-B 互反结构
（GEARS×UCE kappa −0.46/−0.63）提示跨家族**分歧**本身可能是信号：两个
口径家族对同一扰动的"难度估计"强烈不一致时，家族平均误差确实更高。

- 家族门分数（每数据集、每扰动，公共扰动集上对齐）：
  `s_expr = 表达谱家族成员探针均值`；`s_emb = embedding 家族成员探针均值`。
  （主分析：表达谱 2 台、embedding 3 台；诊断 scf_k562：表达谱 3 台。）
- **分歧度（锁定口径）**：`d = |rank(s_expr)/n − rank(s_emb)/n|`（家族均值
  分数各自做数据集内秩百分位后取差绝对值；秩口径免尺度假设，与 G-B2 纪律
  一致）。方向：d 越大 → 家族平均误差越高。
- **目标（锁定）**：家族平均误差秩 `r̄ = 数据集全部引擎的 r 的均值`；
  标签 `z = 1[r̄ ≥ median(r̄)]`（家族口径下"最差一半"）。
- 统计：AUROC(d, z)，扰动级 bootstrap B=1000（master seed 同上），
  95% CI；CI 与 0.5 的关系作描述性记录（不进门限）。
- **判据 P2（锁定，任务书原文）**：K562 与 RPE1 两数据集的 AUROC(d, z)
  **点估计均 ≥ 0.6** → 家族结构有预测性后果（异质原则：分歧可作第三
  信号）；任一数据集 < 0.6 → P2 不成立。诊断 scf_k562 描述性报告。
- 描述性（不判定）：两家族各自与全体引擎均值探针的 |秩差| 变体；
  分歧度与 r̄ 的 Spearman 相关。

## 4. P3 —— 选择原则（家族代表性选择的元学习雏形）

**动机（锁定）**：给定预算只训/部署一个引擎时，家族代表性如何选。回放：
按引擎类型 + 数据集属性预测哪个引擎就地 AUROC 最高。

- 训练档案（**全部 on-file，零重算**，来源 `outputs_rebuild/
  ensemble_gate.json` 的 `auroc.own_in_situ` 列 + `n`/`label_balance`，
  及 G-B2 `shared_difficulty.r2_leave_self_out`）：
  | 行 | 数据集单元 | 引擎数 |
  |---|---|---|
  | K562 | replogle_k562_essential | 5（gears/scgenept/scgpt/geneformer/uce） |
  | RPE1 | replogle_rpe1_essential | 5（同上） |
  | scf_k562(diag) | K562 45 扰动子集 | 仅 scf 一行（其余引擎该子集就地 AUROC 无 on-file 记录，不补算） |
- 特征（锁定）：引擎家族（表达谱/embedding 哑变量）、数据集单元
  （n_perts、label_balance）、引擎×单元共享难度 R²（G-B2 on-file）。
  就地 AUROC = 目标；每数据集单元的决策 = argmax 引擎。
- 锁定规则（零超参拟合的**原型规则**，二选一于判定前公开）：
  - **R1（家族感知）**：在训练单元上按家族聚合就地 AUROC 均值，预测家族
    = 均值较高者；家族内取训练单元均值最高者（ ties 取跨家族总均值高者）。
  - **R0（家族无感基线）**：同上但不在家族内取——直接取训练单元总均值
    最高的引擎（家族信息被抹除的对照）。
  - 参考：随机选择的期望 regret（对训练单元均值取 best−随机期望）。
- 评估（锁定）：留一数据集单元交叉验证（LODO）。主折：K562↔RPE1 双向
  （每折训练 5–6 行、测试 5 行）；scf_k562 单 scf 行不参与数值判定
  （稀疏，避免引擎池不对齐污染），其描述性记录：R1/R0 在 K562 训练下是否
  选中 scf。
- 指标：top-1 命中（预测引擎 == 该单元真实 argmax）；regret =
  best_AUROC − predicted_AUROC。
- **判据 P3（锁定）**：R1 在两主折上 top-1 全中 **且** 两折 regret 均 <
  R0 同折 regret → 选择原则成立（雏形级）；否则不成立。
  （明示：仅 2 个全引擎折 + 11 行档案，统计功效极低，判据为雏形级门槛，
  数字无论过不过均全量报告；不过线则论文不得使用该原则，只报告数字。）

## 5. 统计与防泄漏（锁定）

- 全部 replay：分数向量固定（OOF 探针，确定性重算同 G-B2），bootstrap
  只重采样扰动索引；不重拟合任何门。
- P3 完全不 bootstrap（on-file 点估计的确定性回放）；regret 为点值。
- 多重性记录（描述性，不改判定）：P1 主判定 12 单元（细格 48），P2 两个
  数据集点估计 + CI，P3 两折；判据阈值均为任务书/本节锁定值，照原文执行。
- GBK 控制台纪律：UTF-8 重配置，不打印 emoji。

## 6. 产物与复跑

- `outputs_rebuild/family_consequence.json`（全 spec、种子、P1 全部 Δ/CI、
  P2 两数据集 AUROC/CI、P3 档案表与 LODO 结果、三条独立判定）
- `outputs_rebuild/family_consequence.png`（三联图：P1 每 (D,F) Δ 均值+CI、
  P2 分歧度 AUROC、P3 LODO regret 对比，300dpi）
- `outputs_rebuild/family_consequence_verdict.md`（判定留痕）
- `outputs_rebuild/smoke/family_consequence.smoke.{json,png}`（--smoke：
  B=100、仅 K562；smoke 数字不进判据）
- 脚本：`step24_family_consequence.py`；环境：本机 venv，CPU。
- 用法：
  `PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step24_family_consequence.py --smoke`
  `PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step24_family_consequence.py`

## 7. 与既有档案的关系

- G-B 给家族结构的定性 kappa 图；G-B2 给"跨家族无残差互补"的定量阴性；
  本步回答"那家族结构**能用来做什么**"：省部署钱（P1）、当第三信号
  （P2）、指导选型（P3）。
- 三条判定独立：任一不过线只禁掉对应原则句，不动其余两条，不改 R2
  结论（产品 = 符号层），不删格、不换阈值、不做事后重定义。
- 阴性结果同样归档：原则句只写判定过的条款。
