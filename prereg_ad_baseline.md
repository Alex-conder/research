# 预注册：AD（additive disagreement）基线实算 + 双口径稳健性（gating 论文附录补强）

- 日期：2026-10-30
- 状态：**已锁定**（判据锁定于任何 AD 计算运行之前）
- 起因：`meeting_prep/同行扫描_虚拟细胞可靠性_20261030.md` §3 第 4 条——Gao et al.
  (arXiv 2610.18880) 提出 additive disagreement（AD）为失败检测信号：Norman 组合 OOD 上
  failure-detection AUROC 全转录组 MAE 口径 0.990，但 Top20-MAE 口径掉到 0.619——
  "口径敏感性"是该领域已知陷阱；我们的符号门/探针必须并排跑 AD 基线、且双口径全报。
- 纪律：判定两档制，事先锁死，跑出来哪档写哪档，无中间态；失败也归档。如实报全部数字。

## 引用论文的操作化口径（从原文锁定）

- 效应量：δ_c = μ_c − μ_ctrl（条件均值 − ctrl 均值，log-normalized 空间）
- 组合加和参考：δ_A + δ_B；资格 = 两成分基因的**观测**单基因条件都存在
- AD(A,B) = (1/G)·||f_θ(A,B) − (δ_A+δ_B)||_1，即预测 delta 与加和参考的逐基因 L1 距离取均值
- 失败检测 AUROC 在 Top50%/Top25%/Top10% 最高误差档报告；**方向 = AD 越大越可能失败**
  （原文："predictions with larger deviations from the additive baseline tend to have
  larger errors"）
- Top20 误差 = 每留出扰动按 |观测效应| 取前 20 基因上的 MAE
- 原文明确："AD is undefined for single-gene OOD on Adamson"（单基因留出设定下
  逐基因加和参考无定义）

## 误差三口径（本实验，锁定）

对每测试扰动 i，pred_delta = pred_expr − ctrl_mean，true_delta = 条件真值均值 − ctrl_mean
（ctrl_mean = h5ad 全部 ctrl 细胞均值，与既有产品同口径）：

- `full`：e = mean_g |pred_delta − true_delta|，全部基因
- `top20` / `top50`：按 |true_delta| 降序取每扰动前 20/50 基因上的同式 MAE
  （top20 = 论文口径；top50 桥接本项目产品历史口径 TOP_N_DEG=50）

失败标签（三口径各自独立）：y_m = 1[e_m > median(e_m)]——等价论文 Table 2 的
Top50% 档，且与本项目探针中位数约定一致。

## 分析 A：5 数据集单基因留出产物（AD-mean，偏离口径，次级证据）

组合（GEARS 产物，split 一律 SPLIT_SEED=0 / HOLDOUT_FRAC=0.2 重放核对）：
k562 seed0/1/2（outputs/，n=217）、rpe1 seed0（outputs_rpe1/，n=306）、
norman seed4（outputs_multidataset/norman/，n=20）、adamson seed3/5（n=16）、
dixit seed4/5（n=3，预期两少数类 <2 不可评估，照实记）。

- **协议偏离 D-AD1（锁定）**：论文声明 AD 对单基因留出无定义（成分基因自身的
  单基因效应随留出进了测试集，训练集中不存在）。操作化为 **AD-mean**：
  加和基线退化为训练集单基因效应均值 mean_train_delta（全体训练单基因条件
  delta 的等权均值；"训练条件" = 不含任何留出基因的条件，重放 prepare_split
  simulation 语义），AD-mean_i = mean_g |pred_delta_i − mean_train_delta|。
  AD-mean 不是论文逐基因 AD，只作次级证据，**不进判定档线**。
- 就地探针：step6 规格（9 谱特征 + log_n_cells，in_pert_graph 非常数才入），
  y_m 三口径重拟合；评估 = 同组合 OOF（n≥10 分层 5-fold seed=0；5<n<10 LOO；
  少数类 <2 不可评估——修正案 1 同规则）。
- 符号门：two_layer_gate 同构造（兄弟数据集就地全量拟合系数取符号；目标侧
  无监督 z-score 点积）。**符号源锁定为兄弟数据集标准标签
  （1−pearson_delta>median，step6 口径）拟合，三口径共用同一符号向量**——
  测的是"同一门"跨口径的稳定性，不做逐口径符号重拟合。
  兄弟映射（锁定）：k562←rpe1(s0)；rpe1←k562(s0)；norman←adamson(s3+s5 合并拟合)；
  adamson←norman(s4)；dixit←norman(s4)（dixit 预期不可评估，符号只留痕）。

## 分析 B：Norman 组合重放（论文忠实口径，主判定场）

- 零训练：加载已归档 `outputs_multidataset/norman/gears_model_seed4` 权重（CPU，
  CUDA_VISIBLE_DEVICES="" 强制 CPU，不碰在跑训练），对 norman h5ad 中
  **≥1 成分基因 ∈ 留出基因**的组合条件做 on-file 推断（模型对组合 OOD）。
- 入组资格（锁定，四条全满足）：(i) h5ad 真实组合条件；(ii) ≥1 成分 ∈ 21 个留出基因；
  (iii) 两成分基因均有观测单基因条件（论 文资格）；(iv) 两基因可 predict
  （∈ model.pert_list；预测失败的组合留痕剔除）。
- **AD-combo（论文忠实）**：AD_i = mean_g |pred_delta_i − (δ_A+δ_B)|，
  δ_g = 观测单基因条件均值 − ctrl_mean。**基线用全部观测单基因条件**（含留出
  单基因——论文语言 "observed single-gene effects"；模型从未见过留出基因；
  组合真值不进基线，无循环）。此为对论文"两成分均在训练中"设定的**偏离
  D-AD2**，照实披露（我们的留出划分使"模型见过两成分单基因"与"组合 OOD"不可兼得）。
- 就地探针：留出单基因（n=20）上每口径全量拟合 → 打分组合（部署外推，留痕）；
- 符号门：adamson(s3+s5) 符号，零训练。
- 预检：57 个组合含 ≥1 留出基因。

## 统计（锁定）

- 每（组合×口径×信号）报 AUROC（分析 A 逐种子，数据集汇总 = 种子均值，同本项目惯例）。
- 配对差 AD−符号门、AD−探针的扰动级 bootstrap：**固定分数**重采样（探针 OOF/
  部署分数拟合一次后不再重训；AD/符号门本就无训练），B=2000，95% percentile CI，
  BOOT_MASTER_SEED=20261030。

## 判定（锁定，两档制，主判定场 = 分析 B）

- **CONFIRMED（口径敏感性 = 领域通用现象）**：同时满足
  (a) full 口径 (AD-combo − 符号门) 配对 bootstrap CI 不含 0 且偏 AD；
  (b) top20 口径 (AD-combo − 符号门) CI 含 0 或偏符号门（反转或不显著）。
  附带同规则对照 (AD-combo − 探针)（次级）。
- **NOT_CONFIRMED**：其余一切情形（AD 在 full 不显著优于我们，或 top20 仍显著优）。
- 分析 A 五数据集表（AD-mean vs 探针 vs 符号门 × full/top20/top50）为次级证据，
  照实全报，不进档线。
- 结果无论哪档，论文动作相同：双口径（含 top50 桥接）全报 + 与 AD 并排。

## 产物与附录定位

- `outputs_rebuild/ad_baseline.json`（全部数字）、`outputs_rebuild/ad_baseline.png`
  （分组柱状：数据集 × 口径 × 三信号 AUROC + norman 组合面板）
- `outputs_rebuild/ad_baseline_verdict.md`（两档判定 + 偏离记录 D-AD1/D-AD2 +
  README 式一句话：进 gating 论文附录"可靠性信号对照：加性偏离（AD）基线与
  双口径稳健性"）
- `outputs_rebuild/ad_baseline_norman_combo_pred.npz`（组合预测留档，可复跑核对）
- 脚本 `step21_ad_baseline.py`（`--smoke` 迷你规格 = 每数据集首种子 + 组合上限 10，
  不进判定）；本预注册。

## 执行约束

- 本机 venv、CPU、单进程；**零训练**；唯一推断 = 分析 B 组合预测（归档权重，CPU）。
- 归档 outputs*/、data/ 只读；新产物只写 outputs_rebuild/。
- GBK 控制台纪律：stdout 只打 ASCII，文件全 utf-8，不打印 emoji。
- smoke 先行，全量随后；同规格可复跑（种子全锁定）。
