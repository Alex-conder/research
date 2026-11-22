# 预注册：门控论文升层 G-B —— 跨引擎集成门（stacked meta-probe）

- Date: 2026-10-30
- Status: recorded **before** running `step17_ensemble_gate_gb.py`. 纯 on-file 分析：
  不训练任何新基础模型；唯一拟合对象是家族既定规格的浅层 logistic 探针
  （StandardScaler + L2 logistic, C=1.0，同 step6/step6f/step8b/step10 全族）。
- 上位依据：任务书 G-B（跨引擎集成门）；口径与组合注册表复用
  `step10_signature_resolution.py` / `prereg_rebuild_v2.md`（R2），不另造。

## 1. 数据盘点（锁定，运行前逐文件核过）

| 组合 | 引擎 | 数据集 | 口径 | 扰动数 | 状态 |
|---|---|---|---|---|---|
| gears_k562 / gears_rpe1 | GEARS | K562 / RPE1 | 表达谱 pearson | 217 / 306 | 在 |
| scgenept_k562 / scgenept_rpe1 | scGenePT | K562 / RPE1 | 表达谱 pearson | 217 / 306（与 GEARS 同名单同序） | 在 |
| scf_k562 | scFoundation | K562 | 表达谱 pearson | **45**（强度三层抽样 seed=0，为 217 的子集） | 在，仅子集 |
| scf_rpe1 | scFoundation | RPE1 | — | — | **缺失**（R2 起已登记） |
| scgpt / geneformer / uce × k562/rpe1 | 三 embedding 引擎 | K562 / RPE1 | embedding 距离 | 217 / 306（同名单同序） | 在（wave/ 已校验） |

- 标签（家族锁定口径，不改）：表达谱 `y = (1 − pearson_delta) > median`（step6 定义）；
  embedding `y = err_emb_dist > median`（2026-10-16 锁定口径）。scf_k562 标签在其 45 扰动
  内部取中位数。
- 可评测组合数（锁定）：**10**。= 12 个注册组合 − scf_rpe1（产物缺失）− scf_k562
  （无同引擎另一数据集产物 → 符号臂不可得，与 R2 修正案 A2 同处理，仅入诊断，
  不进锁定计数）。任务书"≥8/10"与本盘点一致，不改。

## 2. 集成门（stacked meta-probe）设计（锁定）

对每个评测组合 `E@D`（引擎 E 在数据集 D）：

**底层特征（三块）**
1. 基础块：数据集 D 上**每个可得引擎的就地探针 OOF 分数**（5-fold StratifiedKFold
   shuffle, seed=0，fold 内 StandardScaler + L2 logistic C=1.0；同 R2 `oof_probe_scores`
   规格）。该分数对该引擎自身标签是 OOF 的；对目标标签 `y_E` 而言，
   非 E 引擎的分数完全不接触 `y_E`，E 自身分数对 `y_E` 也是 OOF 的。
   - K562 目标（217 扰动）：基础块 = 5 个全覆盖引擎（gears / scgenept / scgpt /
     geneformer / uce）。**scf 不入基础块**（172 行无分数，不用插补）。
   - scf_k562 目标（45 扰动）：基础块 = 6 引擎分数按扰动名对齐到该 45 子集。
   - RPE1 目标（306 扰动）：基础块 = 5 引擎。
2. 符号块：R2 符号层分数 1 维（同引擎另一数据集就地拟合系数符号 · 本组合无监督
   z-score 特征；scf_k562 无符号块）。
3. 幅度块：本组合自身的幅度族特征原值（表达谱 9 谱统计 / embedding 4 几何特征，
   由注册表加载器给出；上层 fold 内标准化）。

**上层模型**：logistic（同规格），**leave-one-perturbation-out（LOPO）**——逐扰动
留出：在其余 n−1 行上 fit StandardScaler + L2 logistic C=1.0，预测留出扰动的
概率。n=217/306/45 次拟合。主读数为 LOPO 全样本 AUROC。

**防泄漏纪律（锁定）**
- 底层分数全部 OOF（5-fold, seed=0），且跨引擎分数互不经手目标标签；
- 上层 LOPO，留出扰动的标签从不进任何拟合；
- bootstrap 只重采样扰动索引、不重拟合模型（分数向量固定；LOPO 级重拟合
  ~2.2×10^6 次在 CPU 上不可行，此选择在判定前公开记录）。

## 3. 对照臂（锁定）

- **主对照 = GEARS 就地探针**（任务书指定"最优单引擎（GEARS 就地）"）：
  GEARS@D 的 OOF 探针分数作为 `y_E` 的预测器（跨引擎使用，分数本身不接触 `y_E`）。
- **次对照 = 符号层**（R2 符号臂，同上符号块）。
- 诊断（不进锁定计数）：本组合自身就地探针（own in-situ，AUROC 口径同 step6）。

## 4. 主终点与判定（锁定）

- 主终点：每组合 AUROC（LOPO OOF 概率 vs `y_E`）。
- 臂间差异：扰动级配对 bootstrap（B=1000，master seed 20261030，逐组合由主 seed
  派生 Generator），记录 `ensemble − GEARS_in_situ` 与 `ensemble − symbol` 的
  AUROC 差 95% percentile CI。
- **判定：集成 AUROC 在 ≥8/10 组合上优于 GEARS 就地臂（CI 不含 0）→ 集成有效，
  进主 claim。**
- **阴性条款（锁定，与任务书一致）：未达阈值则如实报告"集成无增益"**——
  这本身是层次结论：说明各引擎误差签名互不携带信息，与深潜
  `outputs_crossmodel/analysis/error_signatures.json` 的"签名噪声"判定互证。
  阴性不改写 R2 结论（产品仍为符号层）。

## 5. 跨引擎一致性矩阵（同步产出，定义锁定）

- 对象：每数据集内引擎两两组合（K562：6 引擎 15 对；RPE1：5 引擎 10 对）。
- 单元：两引擎误差二值（各自口径 y=1 = 最差一半）在**两引擎可得行的交集**
  （K562 全覆盖引擎对 n=217；涉 scf 的对 n=45；RPE1 n=306）上的 2×2 列联：
  `n11`（同错）/ `n10` / `n01` / `n00`（同对），+ 一致率、Cohen's kappa、phi。
- 解读（锁定）：kappa≈0 → 各引擎错误近似独立（集成增益有空间，但需实测；
  独立≠可预测）；kappa 显著 > 0 → 错误共享（集成上限低）。该矩阵回答
  "不同引擎的错误是否独立"——比单引擎门高一层的科学问题。

## 6. 产物与复跑

- `outputs_rebuild/ensemble_gate.json`（全 spec、种子、逐组合 AUROC/CI、一致性矩阵、verdict）
- `outputs_rebuild/ensemble_gate.png`（逐组合臂间 AUROC 对比 + kappa 一致性热图，300dpi）
- `outputs_rebuild/ensemble_gate_verdict.md`（判定留痕）
- 脚本：`step17_ensemble_gate_gb.py`（`--smoke` 迷你规格先行，smoke 数字不进判据）
- 环境：本机 venv，CPU；控制台 UTF-8 纪律沿用（不打印 emoji）。

## 7. 与既有档案的关系

- R2（prereg_rebuild_v2.md）判"产品 = 符号层"；G-B 是**升层**检验：若集成有效，
  主 claim 增加跨引擎集成门一层；若无效，R2 结论不变，G-B 阴性结果与
  error_signatures 深潜互证。
- 所有阴性、混合结果保留在产物中，不删组合、不换阈值、不做事后重定义。
