# 化学/药物扰动小规模探针 · 数据源侦察报告

日期：2026-10-16 · 范围：仅侦察与数据落盘，未跑任何模型训练/推理

## 1. 结论速览

- **首选数据源：Tahoe-100M**（tahoebio/Tahoe-100M，HF 官方仓库，CC0-1.0）。已本机落盘 19 个 parquet 分片 + 4 个元数据表（约 1.6 GB，覆盖 plate 1/4/5/8/10/12/13，229 药物、50 细胞系、含 DMSO 对照与 0.05/0.5/5.0 µM 剂量维度），并已完成 **100 药物探针子集**（`tahoe_probe_cells.h5ad`，265,002 细胞，Morgan FP MaxMin 选样，MOA 24/24 类全覆盖）。
- **sci-Plex3（首选候选 #1）数据本身理想**（188 药物 × 4 剂量 × 3 细胞系，单细胞，GEARS 系论文常用化学 OOD 基准），但其公开镜像（figshare 22122701、GEO GSE139444）在本机网络均被阻断（figshare API/ndownloader 返回 403/202，NCBI 连接失败）。留作后续走代理/服务器侧的备选。
- **LINCS L1000 不适合**做本探针底库：bulk 而非单细胞、仅 978 landmark 基因，与 Replogle 单细胞基线不可比；仍是药物分子特征/scaffold 信息的辅助资源。
- **GEARS 不支持化学扰动类型**（pert 必须能在 GO 基因图中查到），本探针的基模型需用 chemCPA / biolord / scFoundation(嵌入+GEARS 改造) 等替代路径，详见 §6。

## 2. 数据源对比表

| 维度 | sci-Plex3 (Srivatsan 2020) | Tahoe-100M (Zhang/Vevo 2025) | LINCS L1000 |
|---|---|---|---|
| 类型 | 单细胞化学扰动 (sci-RNA-seq) | 单细胞化学扰动 (Mosaic/parse) | bulk 表达谱 (L1000 检测) |
| 规模 | ~650k 细胞，188 药物，3 细胞系 (A549/MCF7/K562) | ~100M 细胞 (95.6M 已发布)，~380–1,138 药物×剂量组合，50 细胞系 | ~1.3M 谱，~20k 化合物 |
| 剂量维度 | 4 剂量 (10nM/100nM/1µM/10µM) + vehicle | 3 剂量 (0.05/0.5/5.0 µM) + DMSO（按 plate 分组） | 多剂量但 bulk 汇总签名 |
| 时间维度 | 24h (固定) | 固定处理时长 | 6h/24h 等 |
| 批次 | 多次测序实验混合 | 14 plates，明确的 plate/batch 列 | 两期 (Phase I/II) |
| 格式 | 原始 count matrix (GEO)；curated h5ad/MuData (figshare 22122701, pertpy) | h5ad (GCS, 13–37 GB/plate) 或 HF parquet 分片 (3,388 shards × ~28k 细胞) | GCTX (level 3/4/5) |
| 公开可得性 | GEO GSE139444；figshare 22122701；pertpy `srivatsan_2020_sciplex3()` | CC0-1.0，免申请；HF / GCS / DNAnexus / lamin 均可 | 完全公开 (GEO GSE92742/GSE70138) |
| **本机可达性** | **✗ figshare 403、NCBI 不可达** | **✓ hf-mirror.com 与 GCS 直连均通** | 未测（底库定位，优先级低） |
| GEARS 兼容性 | 同左：GEARS 不原生支持药物 pert；需改 condition 语义 | 同左 | N/A（bulk） |
| 许可 | GEO 使用条款（公开科研用） | CC0-1.0 | 公开 |

## 3. 选定源（Tahoe-100M）数据概况

已落盘（`E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/chemical_probe/`）：

| 文件 | 大小 | 内容 |
|---|---|---|
| `tahoe100m_hf_shard0..12,242,968,1452,1936,2420,2904.parquet` | ~72–122 MB × 19 ≈ 1.6 GB | 19 个分片 × 28,225 细胞：稀疏表达 (genes idx + expressions) + drug/cell_line_id/plate/moa-fine/canonical_smiles；覆盖 plate 1/4/5/8/10/12/13（即 0.05/0.5/5.0 µM 三档剂量 + plate13 混合） |
| `tahoe_drug_metadata.parquet` | 40 KB | 379 药物：targets, moa-broad/fine, 是否上市/临床, canonical_smiles, pubchem_cid |
| `tahoe_cell_line_metadata.parquet` | 19 KB | 50 细胞系 × 驱动基因注释 (DepMap ID, organ, driver variant) |
| `tahoe_gene_metadata.parquet` | 1.3 MB | 62,710 基因 symbol ↔ ensembl ↔ token_id（分片 genes 列 = 该表行号） |
| `tahoe_sample_metadata.parquet` | 66 KB | 1,344 样本：plate, QC 均值, drug, **drugname_drugconc（剂量）** |
| `tahoe100m_plate3_tutorial_2k.h5ad` | 7.3 MB | 官方 tutorial 子集（2k 细胞），用于核对 h5ad schema |
| **`tahoe_probe_cells.h5ad`** | **810 MB** | **探针子集：265,002 细胞 × 62,710 基因；100 药物 + DMSO 对照 13,607 细胞；50 细胞系；剂量 0.05/0.5/5.0 µM（216,831 / 24,491 / 10,073 细胞）**；obs 含 condition/cell_type（GEARS 兼容列名） |
| `probe_drug_selection.csv` | 13 KB | 100 个选定药物（MaxMin scaffold 多样性）：moa、SMILES、pubchem_cid、细胞数 |
| `probe_summary.json` | <1 KB | 子集汇总指标 |
| `build_probe_subset.py` | — | 可复跑的子集生成脚本（从分片 glob 自动发现输入） |

关键数字（19 分片全量统计）：

- 536,271 细胞中 229 个独特药物（含 DMSO 标签）；候选过滤（≥300 细胞且有 SMILES）后 222 个
- **100 药物选样完成**：RDKit Morgan FP (radius 2, 1024b) + MaxMin 多样性贪心，**moa-fine 覆盖 24/24 类**（222 候选即覆盖全部 24 类）
- 50 个细胞系在每一分片均有覆盖；每 drug×cell_line×dose 组上限 200 细胞
- **对照**：`sample_metadata` 中 28 个 DMSO 样本（每 plate 2 个），探针子集含 13,607 个 DMSO 细胞
- **剂量分布**：0.05 µM × 424 样本、0.5 µM × 402、5.0 µM × 490、DMSO × 28；plate→剂量映射：plate 1/4/7/10=0.05，2/5/8/11=0.5，3/6/9/12/14=5.0，13=混合；**HF 分片按 plate 连续排布**（shard 0–12=plate4，242=plate5，968=plate8，1452=plate10，1936=plate12，2420=plate13，2904=plate1）
- 表达矩阵：每细胞 ~0.7–1.9k 个非零基因（Parse 试剂盒基因数偏少，注意与 Replogle ~5–8k 基因检出对比）

## 4. 与现有基因扰动数据（Replogle K562/RPE1）的关键差异

| 维度 | Replogle (基因 KO) | Tahoe-100M (化学) | 对门控框架的影响 |
|---|---|---|---|
| 扰动标识 | 基因 symbol（有 GO 图先验） | 药物名 + SMILES（无语义先验，除非接分子编码器） | GEARS 的 pert 图先验失效；门控信号需改从嵌入不确定性/剂量-反应一致性来 |
| 剂量 | 无（单剂量 KO） | 0.05/0.5/5.0 µM 三档 | **新增维度**：可做剂量-反应单调性作为可靠性信号 |
| 时间 | 固定 | 固定 | 无差异 |
| 批次 | 单一实验 | 14 plates + sample QC 列 | 门控需把 plate 作为 batch covariate，否则可靠性分数会吸收批次效应 |
| 细胞系 | K562 / RPE1（各自单独数据集） | 50 细胞系混合在统一格式中 | 与 Replogle 的单系设定不对齐；探针建议**固定 1–3 个细胞系** |
| 检出深度 | ~5–8k 基因/细胞 | ~1.9k 基因/细胞 | 表达向量维度/稀疏度不同，门控阈值需重新标定 |

## 5. 100 扰动选样方案（已执行，结果见 `probe_drug_selection.csv`）

原则：scaffold 多样性优先（化学探针的核心），MOA 覆盖兜底。

1. **候选过滤**：19 分片中细胞数 ≥300 且有 SMILES 的 **222 个药物**（536k 细胞池）
2. **scaffold 多样性**：RDKit Morgan FP (radius 2, 1024b) 两两 Tanimoto 距离 + MaxMin 贪心选 100 个；实现为 numpy 贪心（rdkit 2026.03/py3.14 的 `MaxMinPicker.Pick` 会段错误，已绕开）
3. **MOA 覆盖校验**：选中药物 moa-fine 覆盖 **24/24 类**（候选池本身即全类覆盖，MaxMin 未造成 MOA 偏斜）
4. **剂量**：选样跨剂量池执行，探针子集天然含 0.05/0.5/5.0 µM 三档（细胞数 216,831 / 24,491 / 10,073；5.0 µM 档偏少——plate12/14 各只有 1 个探测分片）。**剂量-反应分析建议以 0.05 µM 为主档**，0.5/5.0 用于单调性校验
5. **细胞系**：50 细胞系全保留在子集中；正式实验建议固定 A549（lung，与 LINCS/sci-Plex 可比性最好）为主设定，另保留全 50 系做"跨系泛化"子分析
6. **对照**：DMSO 细胞 13,607 个（同 plate 配比），`condition` 列已按 GEARS 习惯标 `ctrl`

## 6. 模型兼容性评估

### GEARS：✗ 不原生支持化学扰动
源码级证据（`gears/pertdata.py`）：
- `PertData.load()` 原生键只有 `norman/adamson/dixit/replogle_k562_essential/replogle_rpe1_essential`，**无 `sci-plex3` 或任何化学数据集**
- 自定义数据走 `new_data_process()`，但 `condition` 语义必须是基因名：`set_pert_genes()` 将 condition 与 GO 基因图 (`gene2go`) 求交，`filter_pert_in_go()` 把不在 GO 图中的扰动**直接丢弃**；`get_pert_idx()` 用基因 symbol 查 pert 节点索引，药物名查不到会返回 None 并破坏图构建
- README 自述 "multiple perturbation types" 指 KO/CRISPRi/过表达等基因型，均要求 condition 为基因名

结论：GEARS 用于化学扰动必须改造（如把 pert embedding 换成药物分子编码器输出），工作量超出 1 周探针范围。

### 替代路径（✓ 支持药物输入）
| 模型 | 药物扰动支持 | 依据 | 适配本探针的代价 |
|---|---|---|---|
| **chemCPA** (scvi-tools) | ✓ 原生：drug one-hot/embedding + dose covariate，训练于 LINCS 与 sci-Plex | theislab/chemCPA README | 低：与本探针 h5ad 格式直接兼容，scvi 生态成熟 |
| **biolord** | ✓ 在 sci-Plex3 上做过 unseen-drug counterfactual | Nat Methods 2024 (PMC11554562) | 低-中 |
| **PerturbNet** | ✓ "unseen chemical and genetic perturbations"，化学侧用药物结构编码 | Yu et al. 2025 (PMC12322087) | 中：需药物结构→latent 编码器 |
| **scFoundation** | △ 扰动预测走 GEARS 集成（基因）；药物仅为 DeepCDR(bulk IC50)/SCAD(敏感性)，非单细胞药物响应 | biomap-research/scFoundation README | 高：需自行把药物编码接进 GEARS 分支 |
| **scGenePT/scGPT** | △ 主要为基因扰动微调 (Replogle/Norman/Adamson)；scDrugMap 等用 LoRA 适配药物响应但非 perturbation-response 标准设定 | SC-Arena 2026 评测、scDrugMap 2025 | 高 |

**建议**：探针基模型用 **chemCPA**（剂量协变量原生支持，与门控的"误差-可靠性"耦合最自然：剂量作为可靠性条件变量）；GEARS 只保留为基因扰动主线的对照。

## 7. 1 周探针执行清单

- [x] Day 1：跑 `build_probe_subset.py` 生成 `tahoe_probe_cells.h5ad`（100 药物 + DMSO，≤200 细胞/drug×cell_line×dose，265k 细胞 × 62,710 基因）+ `probe_drug_selection.csv` + `probe_summary.json` ✅ 2026-10-17 完成
- [ ] Day 1-2：chemCPA 环境搭建（scvi-tools），探针 h5ad → chemCPA 数据格式（AnnData with `drug`/`dose` covariates）
- [ ] Day 2-4：chemCPA 在探针数据上训练（单细胞系 A549 固定 vs 50 系混合两设定）；**此处开始才涉及模型训练，本侦察阶段不涉及**
- [ ] Day 4-5：门控信号定义：per-drug 预测误差（held-out cells）vs 可靠性分数（剂量-反应一致性、per-plate 批次残差、嵌入不确定性）
- [ ] Day 5-6：与基因扰动主线对齐：同一门控阈值标定流程在 Tahoe 探针上复现
- [ ] Day 7：补 sci-Plex3（走可访问 figshare/GEO 的网络环境下载 figshare 22122701 curated h5ad），做 188 药物小规模交叉验证，写 NC 扩展数据段
- [ ] 风险：Tahoe 探针基因检出浅（~1.9k/细胞）→ 门控阈值不可与 Replogle 直接互换，需在报告中明确声明"阈值按数据集重标定"

## 附：网络可达性记录（本机，2026-10-16）

- 通：`hf-mirror.com`（HF 数据集镜像）、`storage.googleapis.com`（Arc VCA GCS 桶）、`github.com/raw.githubusercontent.com`、`ebi.ac.uk`
- 不通：`figshare.com`（API 403 / ndownloader 202 空响应）、`ncbi.nlm.nih.gov`、`zenodo.org`、`huggingface.co`（主站）、`backblazeb2.com`
- 含义：sci-Plex3 的两个标准入口均需换网络环境；Tahoe-100M 有两条独立可用通道（HF 镜像 + GCS）
