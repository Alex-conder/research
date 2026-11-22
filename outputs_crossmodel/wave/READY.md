# READY — 多模型波产物消费就绪说明

> 状态：传输与校验已完成（详见 `pull_validate.log.md`）。本文件说明六份产物是什么、
> 误差口径是什么、step6f 式迁移测试/签名管线怎么消费、建议的第一步分析。
> 生成日期：2026-10-24。

## 1. 文件清单

6 份产物 = 每引擎 × {K562, RPE1} 三件套（`_predictions.npz` + `_labels.csv` + `_meta.json`），
本机根目录：`virtual_cell_constraint_mvp/outputs_crossmodel/wave/`

| 引擎 | 数据集 | 扰动数 | 嵌入维度 | npz 大小 | npz 本机路径 |
|---|---|---|---|---|---|
| scGPT whole-human | K562 | 217 | 512 (CLS) | 27.2MB | `wave/scgpt/scgpt_wave_k562_predictions.npz` |
| scGPT whole-human | RPE1 | 306 | 512 (CLS) | 38.1MB | `wave/scgpt/scgpt_wave_rpe1_predictions.npz` |
| Geneformer V2-104M | K562 | 217 | 768 (mean-pool) | 40.7MB | `wave/geneformer/geneformer_wave_k562_predictions.npz` |
| Geneformer V2-104M | RPE1 | 306 | 768 (mean-pool) | 57.0MB | `wave/geneformer/geneformer_wave_rpe1_predictions.npz` |
| UCE 33L | K562 | 217 | 1280 (CLS) | 68.0MB | `wave/uce/uce_wave_k562_predictions.npz` |
| UCE 33L | RPE1 | 306 | 1280 (CLS) | 95.2MB | `wave/uce/uce_wave_rpe1_predictions.npz` |

配套 `_labels.csv`（探针列）与 `_meta.json`（权重 sha256、协议、命令行、耗时）同目录。
服务器副本：`ah-server:/opt/research/outputs_crossmodel/{scgpt,geneformer,uce}/`，
md5 逐一比对一致（见 `pull_validate.log.md` 哈希表）。

## 2. 口径（与 2026-10-16 口径决定一致）

- 三引擎输出均为 **per-cell cell embedding，不是表达谱**。各自的 `output_type` 键与
  meta 均已标注。
- **误差口径：embedding 类用模型空间距离，不与表达谱类混比。** 具体探针列：
  - `pred_emb_norm` / `true_emb_norm` / `delta_emb_norm`：嵌入与 delta 的 L2 范数；
  - `err_emb_dist`：pred 与 true 嵌入距离（embedding 空间的核心误差量）；
  - `pearson_delta` / `pearson_delta_all` 列**全空（=NaN），设计如此**。
- 池化约定（与 scGenePT/scFoundation 首轮对齐）：ctrl 池 32 细胞 seed 0；
  true 细胞 ≤32（采样 rng 777000+pi，labels `n_cells` = 实际细胞数）；
  `pred_delta_emb = pred_emb_mean - ctrl_emb`（基线已在 npz 内给出，无需重算）。
- **NaN 填充规则（已验证）**：条件真实细胞数 <32 时，`true_emb_cells` 的尾部整细胞
  为 NaN 填充槽；逐行满足 `NaN 槽数 == 32 - true_n_cells`，`true_emb_mean` 只用
  有效细胞平均、恒为有限值。消费方按 `np.isfinite` 掩码使用，勿对原数组整体求均值。
- 扰动标记列：scGPT 为 `target_masked`，Geneformer/UCE 为 `target_deleted`；
  `target_in_input_space` 口径各引擎略有差异（vocab/映射不同），逐引擎留痕。

## 3. step6f 式管线消费方式

- npz 键同构（六份一致，validate_wave.py 强制校验）：
  `pred_perturbations (n,) / pred_emb_cells (n,32,d) / pred_emb_mean (n,d) /
  pred_delta_emb (n,d) / true_emb_cells (n,32,d) / true_emb_mean (n,d) /
  true_n_cells (n,) / ctrl_emb (d,) / ctrl_emb_cells (32,d) /
  ctrl_pool_indices (32,) / target_in_input_space (n,) /
  target_masked|target_deleted (n,) / output_type`。
- 签名管线的最小消费单元是 `pred_delta_emb`（每扰动一条 d 维签名）与
  `err_emb_dist`（签名质量排序）；需要 per-cell 分布时用 `pred_emb_cells`。
- 跨引擎比较只在**同一误差口径内**做（embedding 距离三引擎互比；表达谱
  pearson_delta 只用 scGenePT/scFoundation 首轮产物）。跨口径比较仅限秩相关/
  一致性层面，并需在报告中标注口径差异。

## 4. 建议的第一步分析（分析前需新预注册）

**crossgene 签名在 embedding 类模型上的就地拟合**：以 `ctrl_emb` 为原点、
`pred_delta_emb` 为扰动签名，检验 (a) 同一扰动跨细胞子样本的签名稳定性、
(b) K562/RPE1 间签名秩一致性、(c) `err_emb_dist` 与扰动强度（target_in_input_space、
n_cells）的关系。三引擎同池同 seed，可直接配对。属新分析假设，**动手前须先预注册**
（假设、指标、阈值、排除准则），结果解读以预注册为准。

## 5. 复验与留痕

- 校验脚本：`wave/validate_wave.py`（只读；`python validate_wave.py wave/`）。
- 传输与校验日志：`wave/pull_validate.log.md`（含逐文件 md5、行数、全部通过记录）。
