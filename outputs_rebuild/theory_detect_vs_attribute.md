# 检测–归因不对称（detect–attribute asymmetry）：虚拟细胞可靠性信号的半形式化理论框架

> theory section 素材（gating 论文升层用）。2026-10-30。
> 纪律：本文每一个数字均引用仓库内归档产物（见 §5 索引表），不依赖任何新计算；
> 凡属假设/推导而非归档事实处，显式标注（A1/A2、D1–D3）。
> 地位声明：这是**半形式化**框架——给出可操作定义、可检验推导和经验锚点，
> 不声称定理级严格性；高斯基准与 √n 律均为"校准用基准"，非证明。

---

## 0. 一句话框架

虚拟细胞预测的误差 E 可分解为**共享难度成分 S**（跨引擎可互相预测的"这道题难不难"）
与**引擎特异成分 U**（每家的失效方式）；**检测（detection）只需要对 S 的可排序性，
归因（attribution）则需要对 U 的结构识别**——前者便宜且可迁移，后者昂贵且按上下文逐格付费。
我们的全部归档实证（R² 分层、余弦零分布、迁移不对称、样本量-分辨率曲线、
双口径翻转、同难同易、家族互反）都是这一不对称的投影。

*English key sentence: Virtual-cell prediction errors decompose into a shared difficulty
component S and an engine-specific component U; detecting failure requires only
orderability with respect to S, whereas attributing failure requires structural
identification of U — the former is cheap and transferable, the latter is expensive
and must be paid per context.*

---

## 1. 形式化骨架

### 1.1 分解

设引擎 e（GEARS / scGenePT / scGPT / Geneformer / UCE / scFoundation）、数据集 d、
扰动 i ∈ {1..n}。取定一个误差口径（caliber）后，

> **E_{e,d,i} = S_{d,i} + U_{e,d,i}**　（A1：加性分解；在该口径下近似成立）

- **S（shared difficulty，共享难度）**：由上下文（数据集 × 扰动）决定、与引擎无关的
  "这题有多难"。操作化定义（归档口径，step22）：用留一引擎回归（leave-one-engine-out
  regression）从**其余引擎的误差**预测本引擎误差，其 R² 即 S 的可预测占比。
- **U（engine-specific，引擎特异残差）**：回归剩下的部分——该引擎在这道题上的
  "自家失效方式"。

### 1.2 检测 = 对 E 的可排序性（rankability）

> **定义（检测）**：无真值信号 g(x_i)（探针分数、置信度、AD、门分数……）对误差口径 E
> 检测成功 ⟺ AUROC(g, median-split(E)) > 1/2 且 bootstrap 95% CI 完全在 1/2 之上
> （= 本仓库全部探针实验的锁定判据口径）。

**充分条件**：任何与 S 相关的信号即可检测。因为 S 是 E 的成分，corr(g, S) ≠ 0
⟹ g 对 E 有非零秩相关（A2：不要求识别 S 本身，只要求排序相关）。
**代价声明**：检测不需要知道误差从哪来、是谁的错——只需要"这道题大概率难"。

高斯基准（校准用，非定理）：若 (g, S, U) 联合高斯、corr(g,S)=r、
ρ ≡ Var(S)/Var(E)，则

> AUROC(g, E) = 1/2 + arcsin(r·√ρ)/π

即**任何"难度型"信号的检测上限由共享方差占比 ρ 单调约束**：ρ→0 时一切难度型
信号退化为抛硬币；ρ 越大，越简单的 S 代理（甚至输入侧难度特征：效应强度、细胞数）
就能达到越高的 AUROC。

### 1.3 归因 = 对 U 的结构识别（structural identification）

> **定义（归因）**：给出误差**为什么**发生——即识别 U 对可观测输入特征的依赖结构
> （本仓库操作化为误差签名 β：logistic 探针在 9 个（表达谱）/4 个（embedding）
> 预测侧特征上的系数向量），或等价地识别 U 所在的低维坐标系。

**必要条件**：(i) U 的信号必须在统计上可分辨（特征系数 bootstrap 95% CI 不含零，
即 signature_resolution 的 stable-feature 口径）；(ii) 该结构必须在跨上下文
（跨数据集、跨引擎）中稳定，否则只是单上下文的描述而非"归因"。
我们的归档证据显示 (ii) 系统性失败（见 1.4）：**U 的坐标系是 engine×dataset
涌现性质**，归因因此只能在每个上下文内独立购买，且单价随特征空间维数上升。

### 1.4 数字填充（每条都指向归档文件，详见 §5 索引）

**(a) S 的方差占比分层（共享难度 R²）——框架的"物理常数"。**
留一引擎回归 R²（gb2_residual）：

| 引擎（家族） | K562 | RPE1 |
|---|---|---|
| GEARS（表达谱） | 0.024 | 0.049 |
| scGenePT（表达谱） | 0.082 | 0.023 |
| scGPT（embedding） | 0.006 | 0.118 |
| Geneformer（embedding） | 0.106 | 0.347 |
| UCE（embedding） | 0.090 | 0.362 |

关键读数：表达谱引擎误差中 S 只占 **~2%–8%**（scGPT K562 低至 0.6%），
embedding 引擎在 K562 为 6%–11%、在 RPE1 跳升到 **35%–36%**
（Geneformer 0.347 / UCE 0.362）——同一数据集上，家族的"误差可互相解释程度"
差一个量级。这直接预告了（且被 (b)(c) 证实）：谁的 ρ 高，谁的误差就越"公共"，
越容易被外来信号排序；谁的 ρ 低，谁的误差越"私有"。

**(b) U 不可跨家族预测——检测的"家族边界"。**
跨家族有序对 24 对（2 数据集 × 12 对）中，AUROC ≥ 0.6 且 CI 在 0.5 之上的对数 = **0**
（判据 ≥1）→ 锁定收敛 **family_structure_is_the_whole_story**（gb2_residual verdict）。
家族块均值全部贴 0.5：within_embedding 0.499/0.510，cross 块 0.512/0.448/0.481/0.403。
**一引擎的 U 对另一家族的误差没有正的可预测结构**——归因信息不外流。

**(c) 检测便宜且可迁移；归因不可迁移——签名余弦零分布。**
- 检测：就地 AUROC 0.656（K562 cv）/ 0.797（RPE1 cv）；**同一 K562 探针零重训直打
  RPE1：AUROC 0.758 [0.705, 0.809]，置换 p=0.001**；与就地差距 −0.044
  [−0.078, −0.010]（"small but significant degradation"，step9c 检验 B）；
  delta_max 系数符号迁移仍 0.768。GEARS 检测样本量操作曲线：**n=15 即达
  AUROC ≥ 0.6**（step8c）。
- 归因：跨数据集签名方向余弦 **≈ 0**（scGenePT 0.006，GEARS 0.044——零分布水平；
  dataset_property_verdict D3）；同模型跨数据集签名主导特征互换
  （K562 由 delta_max 主导、RPE1 由 log_n_cells 主导）。**检测迁移、归因不迁移**。

**(d) 归因按上下文付费——签名分辨率成本曲线（signature_resolution）。**
满 n 时 9/4 特征中稳定（bootstrap CI 不含零）的个数：GEARS 2/2（K562/RPE1）、
scGenePT 4/0、scFoundation 1、scGPT 0/3、Geneformer 1/1、UCE 0/1；
verdict A 仅 4/11 格 ≥2 个稳定特征。成本曲线：scGenePT-K562 到 **n=120 才出现
第 2 个稳定特征**（满 n=217 共 4 个）；scGPT-RPE1 **到满 n=306 才稳定 3/4**；
scf-K562 网格内永远 <2；且曲线非单调（GEARS-RPE1 从 n=15 的均值 3.0 起起伏伏，
满 n 只剩 2）。签名 CI 深潜（step8b）：GEARS 9 系数 7 个跨零、scGenePT-RPE1 9/9 全跨零。
**归因分辨率随 n 亚线性增长，且必须每个 engine×dataset 格各付一次。**

**(e) 符号迁移的条件性（two_layer_gate / step9b）。**
符号层（跨数据集零训练移植的系数符号门）vs 就地重排的两层门：
two-layer 在 0/10 格击败符号层 → 锁定产物 **symbol_layer**；但符号 vs 两层排序
Spearman 仅 0.58–0.85、最大秩位移 115–160。**粗粒度符号（谁大概更危险）可迁移；
细粒度幅度分层不可迁移**——这正是"S 级信息可携、U 级信息不可携"在门控结构上的投影。
且符号门在最结构化口径（top20）表现最好（0.583，ad_baseline Analysis B），
说明它吃到的是 U 型结构而非 S。

**(f) 同难同易与局部反转有界（step9c 检验 A）。**
GEARS vs scGenePT 逐扰动优胜 16/217，**远少于独立零分布均值 ~34**（B=10⁴，
双侧 p=2e-4）→ 两引擎误差正相关：S 真实存在且在引擎间共享；
oracle−全 GEARS headroom 仅 0.0152±0.0069（占全 GEARS 0.561 的 2.7%）→
局部反转不是可利用结构。**S 共享、U 私有、且反转空间有界**——与分解一致。

**(g) 效用与检测同源（task_utility）。**
门引导选 top-B 替代真实实验的成本节约：在探针 AUROC 高的格（RPE1 就地 0.797）
随预算 B 单调放大（B=50 对随机门 +0.941 [0.700,1.177]、对幅度门 +1.886）；
在探针弱的格（Adamson/Norman）节约 ≈0 且 CI 跨零；主口径 1/5 数据集阳性、
symbol 臂 0/5。**门能省多少实验，由它的排序质量决定；排序质量由 S 可排序性 +
同上下文 U 可预测性供给**——效用是检测的派生物，不是归因的派生物。

**(h) U 不是标签噪声（dataset_property）。**
采样噪声只占误差 1%–8%（noise floor），证伪"签名是噪声构成"假设；
log_n_cells 主导反映**模型侧数据支持结构**（少细胞条件下的预测回缩）。
即 U 是模型×数据集的结构性失效，不是测量噪声——归因问题真实存在且无法靠
"多测几个细胞"消解。

---

## 2. 三个可检验推导

> 每条给出：陈述 → 机制论证 → 已有归档证据 → 检验方式
> （优先接既有预注册判据；标注【现有】/【待做】）。

### D1. 检测 AUROC 的上限受共享方差占比 ρ 约束；误差口径切换会按 S/U 构成重排信号排名

**陈述**：对任何仅与 S 相关的无真值信号 g，AUROC(g) ≤ 1/2 + arcsin(√ρ)/π（高斯基准，
ρ = Var(S)/Var(E)）；改变误差口径 = 改变 E 的 S/U 构成 = 改变所有检测器的相对排名。
特别地，纯 S 型信号在 U 富集口径下必须退化，U 型结构信号在 S 富集口径下必须退化。

**机制**：AUROC 是 g 与 E 的秩相关的单调函数；E = S + U 中，难度型信号只通过 S
与 E 相关，其可榨取的相关性随 √ρ 封顶。误差口径（full-MAE vs top20-MAE）不是
中性选择：full-MAE 由全转录组基线难度主导（S 富集），top20-MAE 聚焦大效应基因、
各家在这些基因上的分配错误才是差异来源（U 富集）。

**已有证据**：(i) ρ 分层本身（1.4a）；(ii) AD 论文的口径漂移（同行扫描 §1d：
AUROC 0.990 full → 0.619 top20）；(iii) 我们自己的双口径审计复现同构翻转
（ad_baseline Analysis B，同一批 Norman 组合、同一归档模型：
full 口径 AD-combo 0.978 ≫ 符号门 0.335（Δ +0.642 [+0.471,+0.797]），
top20 口径翻成 AD 0.550 vs 符号门 0.583（Δ −0.034，CI 跨零）；
5 数据集表同向：full 口径 AD 在 5 格中 4 格高于符号门（replogle
0.719–0.863 vs 0.094–0.215、norman 0.880 vs 0.360、adamson 0.625–0.688
vs 0.375–0.438；rpe1 例外，0.469 vs 0.445 打平），
top20 口径符号门反超或打平；锁定 two-tier 规则 CONFIRMED 正是"双口径必报"的产物）。

**检验**：
- 【现有】ad_baseline 的 two-tier CONFIRMED 判据（(a) full CI 排除 0 偏向 AD 且
  (b) top20 CI 含 0 或偏向符号门）——D1 预测这应成为该领域**默认结果形态**而非意外；
- 【待做，接 gb2 判据】构造"纯难度探针"（仅用输入侧难度特征：目标基因效应强度、
  n_cells、target_in_input_space，零引擎输出），在 10 个 engine×dataset 格上测误差
  排序 AUROC；预测其 (a) AUROC 被对应格的 √ρ 基准约束（表达谱格 ≤~0.6 量级、
  RPE1 embedding 格允许更高），(b) AUROC 与 ρ 排序一致。若纯难度探针在 ρ≈0.02 的格
  上显著超过 √ρ 基准，D1 的高斯校准失效，框架需修正；
- 【待做】口径扫描（full / top50 / top20 / top10 / DEG-方向口径）× 信号面板
  （AD、难度探针、符号门、就地探针），预测排名按口径单调旋转。

### D2. 归因的样本成本随 U 坐标维数增长、且不可跨上下文合并

**陈述**：把归因操作化为在 d 维候选特征中稳定分辨 k 个 U 坐标（系数 CI 不含零）。
每个坐标的稳定化需要 n ≳ (z/β_j)²（β_j 为标准化效应）；由于 U 坐标系跨上下文
不相关（余弦≈0，1.4c），每个 engine×dataset 格的这笔样本**必须独立支付**，
不存在"合并 5 个数据集一次归因"的捷径。推论：扩大特征空间（d↑）在固定 n 下
**降低**可归因坐标数（多重性惩罚），而不是增加。

**机制**：签名系数的 bootstrap CI 半宽 ~ 1/√n；稳定 = 半宽 < |β_j|/2。U 的方向
是涌现的（跨上下文余弦 0.006–0.044），所以任何先验（别的格的签名）不提供信息。

**已有证据**：signature_resolution 曲线（1.4d：scGenePT-K562 到 n=120 才 2 稳定、
scGPT-RPE1 到 n=306 才 3 稳定、scf 网格内永远 <2、曲线非单调）；
signature_ci 跨零率（GEARS 7/9、scGenePT-RPE1 9/9）；同数据集跨模型签名
余弦 +0.58 vs 同模型跨数据集 +0.01（manuscript README 归档）——上下文间坐标系
近似随机，数据集边界比模型边界更根本。

**检验**：
- 【现有】R1 判据（满 n ≥2 稳定特征的格数 ≥4/12 → tier A；否则 B）——框架说
  tier A 应是该领域的常态 ceiling，tier B（分辨不出稳定结构）是默认；
- 【待做，直接复用 R1 口径】固定 n_grid，把特征空间从 9 扩到通路级（d: 9→50+），
  预测稳定坐标数不增或下降；若显著上升则 D2 被削弱；
- 【待做】跨格合并训练签名（池化 5 数据集）再逐格评估——预测合并签名在逐格
  上的可分辨性**不优于**就地签名（余弦≈0 的直接推论）；这同时是"大模型大数据
  能否归因"这一行业问题的可操作回答。

### D3. 跨家族互反 = U 的家族特定坐标系 + 互补失效模式

**陈述**：设 U_e = A_{f(e)}·θ_e（A_f = 家族装载矩阵，θ_e = 引擎在家族坐标系内的
失效载荷）。跨家族残差预测由 U 协方差的非对角块决定：若两家族 U 独立 → AUROC=0.5；
若家族间存在**互补失效模式**（一族在另一族失效处反而准确）→ 非对角块为负 →
AUROC 显著 <0.5，即归档中的"互反"。且负向显著应对称地出现在双向块中按失效
几何不对称分布。

**机制**：表达谱家族（显式图/加性结构）与 embedding 家族（隐式预训练表示）的
失效几何天然互补：embedding 家族在稀疏证据下 mode-collapse（预测回缩到均值，
误差小但结构错），恰在表达谱家族过度分散（大效应题上噪声大）处显得"准"，反之
亦然。S 在两者间正相关（同难同易，1.4f），U 互补——所以**互反是 U 的指纹，
不是 S 的**。

**已有证据**：gb2 矩阵的方向性负显著：RPE1 cross_expression_to_embedding 块
均值 0.403、3 对 CI 全在 0.5 之下（K562 同向 0.448、2 对）；而
cross_embedding_to_expression 0.481/0.512、within_embedding 0.499/0.510、
within_expression 0.542/0.491 均不显著——不对称的负向结构，纯噪声无法解释
（噪声应对称贴 0.5）。负向格按预注册记为描述性"互反可预测"。

**检验**：
- 【现有】gb2 负向格的块不对称性（上表）+ step9c 检验 A（16/217 < 零分布 34：
  S 正相关，与 D3 互补成立）；
- 【待做，接同行扫描 §4 机制建议】用 AD 论文的交互项分解把每格残差拆成
  "加性偏离 + 交互偏离"：预测加性偏离跨引擎同号（S 的本质）、交互偏离跨家族
  反号（U 坐标系不同的直接读出）；这是对"机制深度"审稿质疑的同行语言回答；
- 【待做】在家族块模型中显式加入 S 的部分回归：预测跨家族系数从"显著负"
  被吸收为 ≈0，即互反全部经 U 通道、不经 S 通道。

---

## 3. 与同行语言对接

### 3.1 AD 论文的口径漂移（同行扫描 §1d）：从事故升级为预测

AD 论文（Gao, arXiv 2610.18880）的核心警告：additive disagreement 在
full-MAE 口径 AUROC 0.990，换 Top20-MAE 掉到 0.619——"可靠性随误差定义漂移"。
**框架的回应**：这不是该论文的局部瑕疵，而是 D1 的必然结果——AD 是 S 型信号
（加性偏离 ≈ 难度代理），full 口径 E 中 S 富集故其排序有效；top20 口径把 E
换成 U 富集构成，任何 S 型信号必须退化，同时 U 型结构信号（我们的符号门）浮现。
我们在自己归档模型上复现了同构翻转（1.4/D1 证据 iii），并把"双口径必报"
固化为 two-tier 判定规则（ad_baseline CONFIRMED）。**对接输出**：给 AD 论文
（及整个 reliability 线）的理论语言——"detector ranking is a function of the
error caliber's S/U composition"；把"报双口径"从 checklist 条目升级为可推导的要求。
反向也成立：符号门在 top20 最优（0.583）提示"U 富集口径才是结构化信号的用武之地"。

### 3.2 SBB 三原则（Vollenweider & Bühlmann, bioRxiv 2026.04.20.719650）

- **Signal（诊断元指标）**：我们的探针/门正是诊断元指标。框架给它加上信息论
  说明书：元指标可榨取的信息 = 它能接触到的成分。难度型元指标（AD、置信度、
  幅度） ceiling = √ρ 律；同上下文拟合的结构化元指标（就地探针、符号门）可超出
  ceiling，但超出部分不可迁移（0.758 vs 0.797；符号层 0.58–0.85 Spearman +
  115–160 秩位移）。
- **Bounds（可识别性上界）**：SBB 的 Bounds 原则在我们的框架里有可测实现：
  归因的可识别性上界 = 每上下文稳定坐标数（signature_resolution 的 R1 判据
  A/B 两级 + step8b 的跨零率）。这把"可识别性"从口号变成可报告的数：
  本矩阵的 11 格中 4 格可分辨 ≥2 个坐标、7 格不可——**归因上界是逐格属性**。
- **Baselines（最强简单基线）**：框架说明"最难打败的基线"取决于口径的 S/U 构成：
  S 富集口径的最难基线是线性加性/AD（VCC 教训：纯 AI 打不过统计基线）；
  U 富集口径的最难基线变成符号门这类零训练结构信号。完整报告 = 双基线 × 双口径。

### 3.3 路线图五问（Frontiers Cell Dev Biol 2026 报告清单）与验证阶梯

| 路线图问题 | 框架语言下的回答（归档锚点） |
|---|---|
| held-out 是什么 | 预注册 split 族（5 数据集 × 多种子，preregistration 附录） |
| 每测试点剩余支持量 | **共享难度 R² 的逐格补数**：specific_fraction 0.64–0.99——"该引擎误差中多少是公共的"即支持量诊断（gb2_residual.json） |
| 最难打败的基线 | 口径依赖：full→AD（0.978）；top20→符号门（0.583）（ad_baseline） |
| 不确定性如何量化 | 两层回答：检测层 = 门排序（可排序性）；保证层 = conformal 臂（G-D，归档欠覆盖 0.8pp 如实报）。框架告诫：**校准声明是检测层声明**，归因声明需要 D2 的额外昂贵证据 |
| 是否有正交验证 | 跨数据集迁移 = 检测正交（AUROC 0.758 保持）+ 归因不正交（余弦≈0）的精确组合（step9b）——"正交验证"必须区分检测/归因两层 |

对路线图"可信虚拟细胞五条标准"中的**校准的不确定性**：框架给出澄清——
可校准的是 S 型不确定性（可排序、可覆盖保证），U 型不确定性的归因校准在
当前样本规模下不可购得（D2）；claim 阶梯（重建<响应<行为）每上移一层，
要求的成分从 S 扩展到 S+U，这正是"claim 越高对 UQ 要求越强"的机制解释。

---

## 4. 论文化定位

### 4.1 推荐：gating 论文 Discussion 的理论小节（主线）

**定位**：gating 论文当前是方法学短文（6 页 + 附录，arXiv q-bio.QM 主、cs.LG 副；
manuscript/README）。本框架最适合作为 Discussion 中一个命名小节
（"A detect–attribute asymmetry in virtual-cell reliability"），约 0.75–1 页：
- 形式化骨架压缩为 1 段（E = S + U；两个操作化定义）；
- 三个推导各 1 段，每段末句给 Predictions（接 §2 的【待做】，接口已备好）；
- 1 个汇总表（= 本文 1.4 的投影表）+ 1 个图建议（ρ 分层 vs 迁移/余弦的散点：
  每格标 ρ、转移 AUROC、签名余弦——一图讲完不对称）。
**理由**：(i) 全部数字已归档，零新计算，符合该文"预注册 + 诚实阴性结果"的卖点
排序；(ii) 它直接回应打回意见里的"机制问题未触及"——把七个分散的阴性/阳性
结果统一成一个机制性陈述，这是从"工具"升"优秀"的关键一跳；(iii) 推导以
Predictions 形式给出，把审稿人的"so what"转成"可证伪的下一步"。
**风险与对策**：短文篇幅装不下完整推导 → 骨架进正文、证明性论证与数字索引
（本文 §5 表）进附录。

### 4.2 备选：独立短文

**何时值得做**：若 gating 论文投稿后审稿人（或 AD 论文正式发表后的社区讨论）
对该框架有实质引用/争议，则扩展为独立短文；素材已够（纯写作工程，无需新计算）。
**体裁**：Perspective/Commentary，不是 research article。
**目标 venue（按预期收益排序）**：
1. **arXiv (q-bio.QM) + ML4H / AIDrugX（NeurIPS 或 ICML workshop）**——最快、
   与 gating 论文互相引流，workshop 读者正是可靠性评测社区；
2. **PLOS Computational Biology "Perspective"** 或 **Bioinformatics "Perspective"**——
   对口"评测方法论"读者，SBB/路线图对话自然；
3. **Cell Systems / Nature Machine Intelligence 的 Comment 类**——需与正式发表的
   AD 论文直接对话，命中率低但影响面最大；建议留作 2 之后的升级选项。
**独立短文的增量素材**（相对 Discussion 小节）：D1–D3 的完整推导、与
SBB/路线图的系统对表（§3.2/3.3 扩展）、以及"行业含义"段（D2 对"大数据能否
解决归因"的回答 = 不能跨上下文合并付费）。

### 4.3 建议时序

gating 论文带 Discussion 小节投稿（本轮）→ 审稿意见回收后决定是否在修订中
扩写 → 论文接收/预印本有引用后，再启动独立短文（venue 按 4.2 排序）。

---

## 5. 数字索引表（归档文件 → 引用值）

| 量 | 值 | 来源 |
|---|---|---|
| 共享难度 R²（K562/RPE1 × 5 引擎） | 0.006–0.118 / 0.023–0.362（表达谱 0.02–0.08 段 vs embedding 0.35–0.36 段） | `outputs_rebuild/gb2_residual.{json,md}`（step22） |
| 跨家族残差预测 | 24 对中 0 对 AUROC≥0.6 且 CI>0.5；块均值 0.403–0.542；负向显著合计 6 对（K562 2、RPE1 4）；verdict family_structure_is_the_whole_story | 同上 |
| 就地检测 AUROC | K562 cv 0.656；RPE1 cv 0.797 | `outputs{,_rpe1}/gating_probe.json`（经 step9b 引用） |
| 跨数据集迁移 | AUROC 0.758 [0.705,0.809]，置换 p=0.001；迁移−就地 −0.0441 [−0.0776,−0.0099]；delta_max 符号迁移 0.768 | `outputs_crossmodel/analysis/crossdataset_transfer{,_verdict}.{json,md}`（step9b）；`local_reversal_verdict.md` 检验 B（step9c） |
| 同难同易 / 局部反转 | 16/217 优胜 < 零分布均值 33.73（B=10⁴，双侧 p=2e-4）；oracle headroom 0.0152±0.0069（2.7%） | `outputs_crossmodel/analysis/local_reversal_verdict.md` 检验 A（step9c） |
| 签名余弦（跨上下文 ≈ 0） | scGenePT K562–RPE1 0.006；GEARS 0.044；同数据集跨模型 +0.58 | `outputs/dataset_property_verdict.md`（D3）；`manuscript/README.md` v2.1；`crossdataset_transfer_verdict.md` |
| 签名分辨率曲线 | 满 n 稳定特征：4/11 格 ≥2（verdict A）；scGenePT-K562 n=120 才 ≥2；scGPT-RPE1 满 n 306 才 3/4；曲线非单调 | `outputs_rebuild/signature_resolution.json`（step10，R1） |
| 签名 CI 跨零 | GEARS 7/9；scGenePT-RPE1 9/9；scFoundation 8/9 | `outputs_crossmodel/analysis/signature_ci.json`（step8b） |
| 符号 vs 两层门 | two-layer 0/10 格胜符号 → symbol_layer；Spearman 0.58–0.85，最大秩位移 115–160 | `outputs_rebuild/two_layer_gate.json`（step10，R2） |
| AD 双口径审计 | Norman 组合：full AD 0.978 vs 符号 0.335（Δ +0.642 [+0.471,+0.797]）；top20 AD 0.550 vs 符号 0.583（Δ −0.034）；探针 full 0.999 / top20 0.475；two-tier **CONFIRMED**；5 数据集表同向 | `outputs_rebuild/ad_baseline_verdict.md`（step21） |
| 任务效用 | 主口径 probe_full 1/5 阳性；RPE1 B=50：对随机 +0.941 [0.700,1.177]、对幅度 +1.886；K562 B=50 同向；Adamson/Norman ≈0；symbol 0/5 | `outputs_rebuild/task_utility_verdict.md`（step22） |
| 采样噪声占比 | 误差中 1–8%（noise floor），证伪标签噪声假设 | `outputs/dataset_property_verdict.md` |
| 检测样本量 | GEARS n=15 即 AUROC ≥ 0.6（终值 0.672） | `outputs_crossmodel/analysis/sample_size_curve.json`（step8c，经 manuscript README 引用） |
| AD 论文口径漂移（同行） | AUROC 0.990（full）→ 0.619（Top20） | `meeting_prep/同行扫描_虚拟细胞可靠性_20261030.md` §1(d)（arXiv 2610.18880） |

## 附：英文关键句（论文级，可直接引用）

1. **框架**：*"We propose a detect–attribute asymmetry: prediction errors of virtual-cell
   engines decompose into a shared difficulty component S and an engine-specific component U.
   Detecting failure requires only orderability with respect to S; attributing failure
   requires structural identification of U. The two goals differ in cost by orders of magnitude
   and in transferability by sign."*
2. **D1**：*"The AUROC of any difficulty-typed reliability signal is bounded by the shared
   fraction of error variance; an error-caliber switch re-weights the S/U composition and
   therefore can reverse the ranking of detectors — metric drift is a prediction, not an
   accident."*
3. **D2**：*"Attribution has a per-context sample cost: the number of stably resolved
   error-signature coordinates grows sub-linearly with n and cannot be pooled across
   engine×dataset contexts, because signature directions across contexts sit at the null
   (cosine 0.006–0.044)."*
4. **D3**：*"Cross-family reciprocal predictability of residuals is the fingerprint of
   family-specific coordinate systems of U under complementary failure modes — negative
   where coordinate systems anti-align, at chance where they are unrelated."*
5. **对接**：*"Detection transfers across datasets while attribution does not: a probe
   fitted on K562 retains AUROC 0.758 on RPE1, yet the direction of its error signature
   vanishes (cosine ≈ 0). Calibration claims are therefore detection-level claims;
   attribution claims require separate, per-context evidence."*
