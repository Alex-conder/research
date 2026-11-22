# 预注册：门控论文升层 G-C —— 门引导的实验预算分配效用实验

- 日期：2026-10-30
- 状态：**已锁定**（全部判据、臂构成、成本模型锁定于任何效用计算运行之前）
- 起因：`meeting_prep/同行扫描_虚拟细胞可靠性_20261030.md` §3 第 6 条——
  任务级效用是 2026 年该领域的层次水位线（VCC"替代一次真实 Perturb-seq
  实验"叙事）；本实验把门从"评分工具"升级为"**虚拟细胞预测何时能替代一次
  真实 Perturb-seq 实验**"。
- 纪律：判定标准事先锁死，跑出来哪档写哪档；无效用也是结果，如实报；
  不删组合、不换阈值、不做事后重定义；GBK 控制台不打 emoji。

## 1. 回放设定（锁定）

模拟"实验预算 B 个扰动名额"：在留出扰动集合（GEARS 单基因留出重放产物，
与 step21 分析 A 同构造同数据）上，每个臂选 top-B 个扰动**采信**其预测
（不做真实实验，支付预测错误），其余 n−B 个扰动的真值被**测定**（每个支付
测定成本 c）。臂间只在"采信谁"上不同。

**总成本（锁定）**：`C_arm(B) = Σ_{i ∈ 采信集} err_i + c · (n − B)`

- `err_i` = 扰动 i 的 full 口径误差 = 全转录组 delta MAE
  `mean_g |pred_delta − true_delta|`（log-normalized 空间；step21 "full"）。
- `c` = 本组合 err 的中位数（一次真实 Perturb-seq 测定与该数据集典型预测
  误差同量级的锁定换算）。**固定 B 下 c·(n−B) 与臂无关，臂间比较与 c 无关；
  c 只进"相对全测定"的绝对节省次级读数。**
- 次级读数（不进判定）：相对全测定（c·n）的节省 = c·B − Σ_{采信} err；
  误差口径敏感性 top20/top50 MAE（同 step21 定义），只进 JSON 不进判定。

## 2. 臂构成（锁定）

**门引导采信臂（臂 a，三个门变体，各自独立判定）**
- 采信方向（锁定，2026-10-30 运行前补记）：探针与符号门的分数是**失败
  似然**分数（step6/step21 约定，高分 = 更可能失败），采信 = 取分数
  **最低**的 top-B（= 最可靠）；conformal_width 采信 σ **最窄**的 top-B；
  幅度基线采信 delta_max **最高**的 top-B（任务书锁定稻草人）。
1. `probe`：就地探针 OOF 决策分（step6 规格，5-fold StratifiedKFold seed=0；
   n<10 用 LOO；标签 y = err > 中位数按 full/top20/top50 三个口径各自重拟合）。
   主臂 = `probe_full`（与成本口径一致）；`probe_top20/top50` 照报。
   已知乐观性（OOF 与评测同集），与家族 step6/step21 惯例一致，如实声明；
   `symbol` 为零训练部署臂与之对照。
2. `symbol`：符号门（two_layer_gate 同构造：兄弟数据集就地全量拟合取系数
   符号 · 本组合无监督 z 分数；符号源锁定 = step21 SIGN_SOURCES）。
3. `conformal_width`：自适应 split-conformal 弃权分数，**宽度 = σ_i =
   ‖delta_pred_i‖₂/√G**（G-D prereg_conformal_arm.md §3 锁定构造，纯模型侧
   无真值泄漏）；采信 σ 最小（区间最窄）的 top-B。

**基线臂（臂 b/c，判定对手）**
- `random`：随机采信，成本取精确期望 = B·mean(err)（线性性；bootstrap 内
  同重采样集配对计算）。
- `magnitude`：幅度基线，采信 `delta_max`（= max_g |pred_delta_ig|，
  step6 特征，预测侧可部署）最高的 top-B。

诊断（不进判定）：各门变体对三口径 err 的失败检测 AUROC；AD-mean 信号
（step21 构造）作第四门变体仅入 JSON 诊断区。

## 3. 数据组合（锁定，运行前逐文件核过）

| 数据集 | 产物 | seed | n 留出扰动 | B∈{5,10,20,50} 可评测性 |
|---|---|---|---|---|
| replogle_k562_essential | `outputs/` | 0,1,2 | 217 ×3 | 全部 4 档 |
| replogle_rpe1_essential | `outputs_rpe1/` | 0 | 306 | 全部 4 档 |
| norman | `outputs_multidataset/norman/` | 4 | 20 | B=5,10；B=20,50 退化 |
| adamson | `outputs_multidataset/adamson/` | 3,5 | 16 ×2 | B=5,10；B=20,50 退化 |
| dixit | `outputs_multidataset/dixit/` | 4,5 | 3 ×2 | 全部退化（n=3 < min B） |

- 退化定义（锁定）：B_eff = min(B, n)；B_eff = n 时所有臂采信全集，
  成本恒等，差恰为 0。退化档不进判定，照实记录。
- 任务书判据为"≥3/5 数据集"：5 个注册数据集如上；dixit 因退化自动记
  **非正**（档档差为 0，CI 必含 0）。等效于 4 个可评测数据集至少 3 个正。
  该解读在本文件锁定，不改。

## 4. 统计与判定（锁定）

- 配对 bootstrap：扰动级有放回重采样（同重采样集内重排序重算两臂成本，
  配对差），B=1000，master seed=20261030，逐组合由
  SeedSequence(master).spawn 派生 Generator（同 G-B 纪律）；bootstrap 只
  重采样索引，不重拟合任何模型（分数向量固定，预先公开记录）。
- 比较量：`diff = cost_baseline − cost_gate`，95% percentile CI；
  CI 不含 0 且为正当且仅当门引导臂总成本显著低于基线臂。
- 组合级判定：对基线 X ∈ {random, magnitude}，
  `combo_positive_vs_X` 当且仅当**每个可评测 B 档** CI 均排除 0 且偏向门。
- 数据集级判定：`dataset_positive(g)` 当且仅当该数据集全部 seed 组合
  `combo_positive_vs_X` 对 X=random 与 magnitude 同时成立（多种子严格制）。
- **总判定（锁）**：门变体 g 在 ≥3/5 注册数据集上 dataset_positive(g)
  → "门 g 在实验预算约束下有正效用"，进主 claim（VCC 叙事）；
  否则如实报"无显著效用/部分效用"，阴性不改写 R2 结论（产品仍为符号层），
  主臂 = `probe_full`。
- 效率注记：n≤306 且成本为 top-B 求和，bootstrap 全组合重排序为 O(n)，
  CPU 分钟级；不占在跑训练资源。

## 5. 产物与复跑

- `outputs_rebuild/task_utility.json`（全 spec、逐组合逐 B 成本/差值/CI、
  AUROC 诊断、判定）
- `outputs_rebuild/task_utility.png`（逐数据集成本-B 曲线：三门 + 两基线 +
  全测定参考线；判定汇总面板；300dpi）
- `outputs_rebuild/task_utility_verdict.md`（判定留痕 + 论文定位段）
- 脚本：`step22_task_utility_gc.py`（`--smoke` 迷你规格先行，smoke 数字
  不进判据；smoke 组合 = gears_k562_s0 + gears_norman_s4）
- 环境：本机 venv，CPU，on-file 零训练（不重训 GEARS、不加载 torch 权重）；
  控制台 UTF-8 纪律沿用。

## 6. 与既有档案的关系

- 数据与信号构造复用 `step21_ad_baseline.py`（分析 A 同盘同协议）：
  load_bundle/load_product/error_metrics/oof_probe_scores/fit_signs/
  symbol_scores/SIGN_SOURCES，不另造。
- conformal_width 构造锚定 `prereg_conformal_arm.md`（G-D）§3 的 σ_i 定义，
  与 G-D 臂同族；本实验回答的是 G-D 不问的**排序-预算效用**问题。
- 阴性/混合结果全保留：这本身是层次结论（门只擅评分、不擅替实验）。

## 7. 偏离记录（锁定披露）

- D-GC1（2026-10-30，运行前）：smoke 首跑启动后发现采信方向实现错误
  （对失败似然分数取了高分采信，与 §2 采信方向相反），**在产出任何判据
  数字前**中断、修正脚本并补记 §2 采信方向句。无任何判据数字来自错误
  实现。
