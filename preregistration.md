# 预注册文档：守恒违例分数 OOD 信号验证（第一阶段）

> 锁定日期：2026-10-12　　锁定后任何判据修改须留痕并标注理由。

## 假设与判据

| ID | 假设 | 成立判据 | 不成立处置 |
|---|---|---|---|
| H2a | 守恒违例分数与 OOD 预测误差正相关 | 三个条件**同时**满足（见下方分解） | 任一不满足按对应条款归档（见下） |
| H2a-null | 信号非任意高维映射的伪相关 | 真实 r 高于置换 null 分布的 95 分位（≥100 次打乱基因-反应映射） | 未超过 null：信号为伪相关，归档并归因，方案终止 |
| H2a-partial | 信号非扰动强度代理 | 控制扰动强度后的偏相关仍显著（r_partial 的 p < 0.05；强度协变量 = 敲低效率或扰动的全基因组效应量） | 偏相关掉光（p ≥ 0.05）：信号是强度代理，归档归因，方案终止 |
| H2b | 违例分数提供独立于不确定性的信号 | corr(违例分数, ensemble 方差) < 0.3 | ≥ 0.3：信号冗余，转"仪表盘"路线（只检测不校正），本阶段结论照常归档 |

**H2b 阈值依据（防"数字无出处"追问）**：0.3 是我们的**判定阈值**，不是对标值。《陌生度》（Nat Mach Intell 2026）实测陌生度×不确定性相关性 0.10±0.05，那是它的结果；我们的期望是做到类似水平（≈0.1），超过 0.3 判定为无独立价值。判据与期望分离写明。

## H1（二阶段判据，与边界条件同条锁定，防"事后补后路"质疑）

| ID | 判据 | 边界条件（预期结果，非失败） |
|---|---|---|
| H1 | o.o.d. 场景下，硬约束校正的 top-DEG Pearson delta 增益 > 未校正，且 > 软约束（KG）增益 | i.i.d. 场景下硬约束净收益 ≤ 0 为**事前锁定的预期**（数据同化标准结论：强投影在观测强时损害拟合）。H1 只在 H2a、H2b 均通过后启动 |

## 评测设置（锁定）

- 数据：Replogle-Nadig K562 CRISPRi 子集（公开）。
- OOD 划分：扰动留出（perturbation holdout），与 scPerturBench perturbation-generalization 设置对齐；留出比例与 scPerturBench 一致。
- 引擎：GEARS，3 个随机种子（ensemble 不确定性 = 种子间方差）。
- 预测误差：top-DEG Pearson delta（对齐 scPerturBench / MAP 口径）。
- 代谢网络：Human-GEM 最新 release；只用高变代谢基因子集（表达变异 top 25% 的代谢基因）。
- GPR 映射：一阶段**只用 AND 关系反应**（多基因同时表达才激活），跳过 OR 关系以降低映射噪声（技术决策，留痕）。
- 违例分数：一阶段**不做归一化**，直接报告 ‖S·v̂‖₂ 的分布（量纲近似是已知简化，先看分布再定归一化方式）。

## 多重比较锁定与功效声明（2026-10-12 增补，出数据前封死）

- **主判据只有一个**：全网络违例分数 ‖S·v̂‖ 的 H2a/H2b 检验。任何子网络（糖酵解、TCA 等）分析均为**探索性**，不进入主判据推断，不报告其 p 值作为证据。
- **多重比较纪律**：除上表锁定的检验外，新增任何假设检验须先在此留痕，否则视为探索性。
- **功效声明**：扰动留出 n 若 < 30，p 值仅供辅助，主报告效应量 r 及其自助置信区间；判据不因样本量放宽，样本不足本身按协议偏离条款记录。
- **尺度控制**：‖S·v̂‖ 随表达整体缩放，分析时以表达总量为协变量或改用相对值（具体方式在 Step 4 实现时锁定并留痕）。

## 排除条款与协议偏离（防"数据不支持就做不了"质疑）

- 不按结果事后挑选代谢基因子集、不更换误差指标、不更换相关性类型；
- **协议偏离条款**：若因数据缺失、基因名映射丢失或 API 限制必须更换口径，须在 `outputs/verdict.md` 中明确记录"原口径不可用，换用新口径，原因是 XX"，标注为**协议偏离**而非协议内结果；偏离前后的结果并列报告；
- 判失败的结论与判成功同等价值，归档格式相同。

## 边界声明

H1 启动条件：H2a 与 H2b 均通过。校正数学形式届时二选一并在此留痕：QP 投影（真硬约束）或守恒正则（软形式；若选此，对比改为"投影 vs 正则"）。

## 探索性分析登记（E1，登记日期 2026-10-13，出数据前）

- **背景**：Day3 主判据 H2a 阴性（r=0.0226, p=0.74），违例分数 CV 仅 5.9%。
  归因假设：绝对表达作酶表达代理引入共同基线（所有扰动的 v̂ 都含同一 ctrl
  背景项），稀释了扰动特异信号。
- **E1 假设**：改用**扰动 vs 对照的预测表达 delta** 作为酶表达代理输入
  （v̂_r = AND 基因 delta 表达均值），违例分数的变异度恢复（CV 明显高于
  5.9%），且 corr(violation, error) 的效应量大于 Day3 绝对表达代理。
- **管线锁定**（与 Day3 同口径，保证可比）：同一 GPR AND-only 清洗
  （3,943 反应）、同一高变 67 基因子集、同一 167 保留反应、同一 S_sub、
  ‖S·v̂‖₂ 不归一化；delta = 预测表达 − 训练数据 ctrl 细胞均值（ctrl 均值
  口径与 step1_run.py 一致）；error = 1 − pearson_delta（同一
  outputs/predictions.csv）。
- **四项读出**（同时报告，缺一不可）：
  1. 违例分数 describe()/CV（对照 Day3 的 5.9%）；
  2. r（violation_delta vs error）的 Pearson 效应量；
  3. 200 次置换 null（同 Day3 方法）；
  4. 控制 strength 的偏相关（残差法，t df=n-3）。
- **多重比较纪律**：E1 按"多重比较锁定与功效声明"条款处理——**不复活 H2a
  主判据、不进主判据推断**；报告效应量与置信区间，p 值一律标注"探索性"。
- **决策边界**：E1 结果只作归因证据写清；若四项全活，是否升级二阶段由
  用户决策，脚本与报告不自行展开。

## 探索性登记（2026-10-16 入口A，出数据前）

- **背景**：误差签名（就地 9 特征 logistic 探针的标准化系数向量）已被项目
  确认为跨模型稳定现象（delta_max/delta_entropy 主导）。归档模型
  `outputs/gears_model_seed0/model.pt` 只有 3-epoch 最终 checkpoint，训练
  中签名何时成形无数据：形态（系数向量）与功能（探针 AUROC / pearson_delta）
  谁先稳定？
- **性质**：协议外探索（exploratory），不进主判据推断；single-model、
  single-seed、单一数据集（n=1），一切时序结论只作假设生成。
- **入口A 假设**：
  - T1（形态时点）：签名对终值签名余弦 trajectory crossing 0.8 的最早
    epoch k₁；
  - T2（功能时点）：探针 AUROC（3 split 种子均值）进入终值 ±0.02 带内最早
    epoch k₂；pearson_delta 均值进入终值 ±0.01 带内最早 epoch k₃；
  - T3（时序判定）：k₁ < k₂ 且 k₁ < k₃ → "签名早期定形、功能后至"
    （支持"误差结构是训练早期由数据/架构决定的功能性属性"叙事）；
    k₁ ≥ k₂ → "后期涌现"。
- **管线锁定**：训练配置与 step1_run 归档逐行对齐（SPLIT_SEED=0、
  HOLDOUT_FRAC=0.2、hidden=64、Adam(1e-3, 5e-4)+StepLR(0.5)、batch 32）；
  训练循环镜像 venv gears.train()（不改包文件）；每 epoch 末快照
  state_dict + optimizer/scheduler + RNG 状态（支持 --resume-snap）。
  阶段1（管线验证）与阶段2（12-epoch 轨迹）合并为一次训练：StepLR 衰减只
  依赖已训练 epoch 数，e3 快照 == 独立 3-epoch 终模型（数学等价，特此声明）。
  每快照分析口径全部复用现有脚本：pearson_delta（step1 top-50 口径）、
  9 特征（step6 compute_features + log_n_cells）、探针（step7 同规格
  StandardScaler + LogisticRegression(C=1.0)，70/30 扰动级划分，
  split_seeds=[0,1,2]，master seed 20261016）。
- **阶段1 判据（gate）**：e3 快照的全量就地拟合系数（= 归档
  error_signatures.json GEARS(K562) 的同口径拟合）与归档签名余弦
  > 0.8 → 管线有效；< 0.8 → 先查 split/数据口径，不解读轨迹。
- **判据纪律**：余弦一律用 L2 归一化系数向量（正则 + logistic + 共线特征
  使原始尺度不稳，已知 caveat）；任何 early/late 结论标注 n=1 限制；
  附带记录每 epoch Train/Validation Overall MSE。
- **降级条款**：12 epochs 若因显存/时间不可行，降级 6 epochs 并在 verdict
  如实说明。
- **产物**：step11_training_dynamics.py；outputs/training_dynamics/
  （快照 + per_epoch.json）；outputs/training_dynamics.{json,png}；
  outputs/training_dynamics_verdict.md。


## 探索性登记（2026-10-16 入口二）

- **背景**：误差签名是模型×数据集特有的涌现性质（GEARS 由 delta_max/delta_entropy 主导、scFoundation 由 log_n_cells 主导、scGenePT 跨数据集完全不迁移且 RPE1 就地探针 AUROC≈0.50）。机制问题：数据集属性能否解释签名差异？
- **D1**：签名主导特征 ⇔ 误差来源构成——log_n_cells 主导的组合其采样噪声占误差比例应更高（排序对照）。判定：**证伪**（n=4 ρ=−0.4；n=5 τ=0；采样噪声仅占误差 1–8%）。
- **D2**：就地 AUROC ⇔ err 对 9 特征回归的 CV R² 正相关。判定：**不足以判定**（n=5 ρ=0.7 为弱正相关但不单调；scFoundation n=45 下 R² 发散）。
- **D3**：跨数据集签名消失 ⇔ 数据集属性漂移（效应强度 KS 0.435、每条件细胞数比 1.63；scGenePT 余弦 0.006、GEARS 0.044）。判定：**方向一致**（纯描述，n=2）。
- **多重比较纪律**：三条均按探索性处理，不进主判据推断，不做显著性检验（n=4~5 只做排序/单调/描述）。
- **脚本**：`step9_dataset_property.py`（主 seed 20261022）；产物 `outputs/dataset_property.json` / `.png` / `_verdict.md`。
