# step8 分水岭实验 verdict：内部隐层表示 vs 输出几何特征预测 GEARS 自身误差

## 设计（预注册，脚本头部同文）

K562 Replogle essential，217 个留出扰动（与 step6/7 同一批）。对每个扰动复刻
GEARS `predict` 前向（`gears/gears.py:349-361`、`utils.py:457-477`：ctrl 细胞
有放回抽 300 个、eval + no_grad），用 hook 取 decoder 前的最终节点嵌入
（`gears/model.py:186-190`，`bn_pert_base`→ReLU→`recovery_w` 的输入），
每个扰动三种 64 维内部读出（hidden_size=64，300 细胞平均）：meanpool（全基因
节点均值）、pertnode（被扰动基因节点）、crossgene（`model.py:197`
`cross_gene_state` 的全局输出）。与输出几何 9 特征（step6 同函数、按模型种子
用 `pred_expr_seed{m}` 重算）在**完全相同**的 split（扰动级 70/30，seeds 0/1/2，
master 20261016，逐行复刻 step7_routing.py:112-119）上用**完全相同**的探针
（StandardScaler + L2 logistic C=1.0/max_iter=2000，70% 就地拟合，30% AUROC，
y = err > 训练侧中位数）对比。GEARS 模型种子 0（主）与 1（复现），不重训，
仅前向。判据预写死：mean diff ≥ +0.05 且 3 个 split 的配对 bootstrap 95% CI
下界全 > 0 → "藏在内部表示"；|diff| < 0.02 → "输出几何已足够"；其余如实报告。

## 数字（test AUROC，3 split 均值±std；diff = 该特征集 − 同样本 geom9）

| 特征集 | seed0（n=217） | seed1（n=217） |
|---|---|---|
| geom9（基线） | **0.703**±0.039（0.746/0.670/0.694，复现 step6 的 0.669/0.694/0.746） | **0.753**±0.058（0.742/0.815/0.703） |
| meanpool64 | 0.499±0.043（diff −0.204） | 0.700±0.033（diff −0.054） |
| crossgene64 | 0.724±0.008（diff +0.020） | 0.797±0.041（diff +0.044） |
| geom9+meanpool64 | 0.727±0.015（diff +0.024） | 0.760±0.038（diff +0.007，tie） |
| geom9+crossgene64 | 0.724±0.016（diff +0.021） | 0.762±0.047（diff +0.008，tie） |
| pertnode64（n=83*） | 0.622±0.040（diff +0.052） | 0.779±0.068（diff +0.109） |
| geom9+pertnode64（n=83*） | 0.662±0.043（diff +0.092） | 0.805±0.150（diff +0.134） |

\* 134/217 个扰动的靶基因不在 5000 维表达空间内（step7 已记录的"pred≈ctrl"
数据事实），cell graph 中无对应节点，pertnode 对它们无定义；仅在 83 个可评估
扰动上评测，geom9 在同一子集重算配对（geom9 子集本身不稳：seed0 0.465/0.603/
0.640，seed1 0.693/0.500/0.819）。

bootstrap：全部非基线特征集 × 2 模型种子 × 3 split 的配对 diff 95% CI **均跨 0**
（最接近的是 seed1 geom9+pertnode64 的 split2 CI [+0.000, +0.347]，split0/1 仍
跨 0）。无任何特征集满足"win"判据（mean diff ≥ +0.05 且全部 CI 下界 > 0）；
仅 seed1 的两个拼接落在 tie 档（|diff| < 0.02）。**全部判定：中间档。**

## 判据判定与结论

落在预写死的**中间档**，证据方向明确偏向"输出几何已是足够近似"：

1. 最强公平读出 crossgene64（扰动信息首次跨基因混合后的全局内部状态）与
   geom9 统计上不可区分（diff +0.020/+0.044，所有 CI 跨 0）——内部表示没有
   携带比输出几何更多的误差签名；
2. 拼接在两种子上无一致增益（+0.007~+0.024），seed0 上甚至无提升；
3. meanpool64 在 seed0 上恰为随机（0.499）——这是**架构事实**而非提取 bug：
   扰动 embedding 只在 `model.py:166-184` 被加到被扰动的那一个节点上，
   随后的 `bn_pert_base`（eval 模式逐通道仿射，无跨节点混合）+ ReLU 即进
   decoder，因此提取点处其余 4999 个节点嵌入**逐扰动不变**，mean-pool 把信号
   稀释 1/5000（前向正确性由 sanity check 背书：复刻预测 vs 归档 pred_expr
   逐扰动 mean pearson = 0.9979）；
4. pertnode 在 83 子集上的正 diff（+0.05~+0.13）主要由 geom9 在该小子集上的
   自身不稳定驱动，CI 全部跨 0，不构成"win"证据。

**一句话结论：在 GEARS 上，误差签名不藏在 decoder 前的内部隐层表示里——
即使取扰动信息首次全局混合的 cross-gene 读出，其误差预测力也不优于（约等于）
输出几何 9 特征；此前路由实验中"弱模型的赢点可能在内部表示"的猜想，在 GEARS
自身被证伪，当前"输出几何门控"的论文边界真实。**

## 关键代码依据与坑（备查）

- 提取点：`venv/Lib/site-packages/gears/model.py:186-190`（hook =
  `recovery_w` 的 forward_pre_hook）；crossgene：`model.py:197`；
  扰动单节点注入：`model.py:166-184`。
- 输入构造：`gears/gears.py:349-361`（predict 路径）+
  `gears/utils.py:457-477`；`predict` 内 ctrl 抽样（`np.random.randint`）
  未播种，本实验固定 `np.random.seed(20261016)`（分布相同，可复现）。
- 加载：不用 PertData，直接 `GEARS_Model(config.pkl)` + `model.pt`
  （config 已含 G_go/G_coexpress）；pert_list 复现自
  `pertdata.py:126-128`（essential_all_data_pert_genes.pkl ∩ gene2go_all.pkl，
  9853 节点，与 config num_perts 一致）。
- 坑：①pert_list（9853 扰动图）索引 vs 5000 维表达空间基因索引不可混淆；
  ②134/217 靶基因不在表达空间 → pertnode 子集协议；③mean-pool 结构性
  不变（见上）；④step6d.load_dataset 硬编码 SEED=0，几何特征按模型种子
  重写加载；⑤GEARS `predict` 有 `saved_pred` 缓存，提取自行构造输入绕过。

产物：`step8_internal_repr.py`、`outputs/internal_repr.json`（全部数字：
7 特征组 × 2 模型种子 × 3 split 的 AUROC、AUROC CI、配对 diff CI）、
`outputs/internal_repr.png`（300 dpi，英文）、
`outputs/internal_repr_cache_v2_seed{0,1}.npz`（内部表示缓存，可复现）。
