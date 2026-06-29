# 脑缺血后 Pirb 阳性细胞单细胞图谱 — 多数据集初步分析报告

## 一、已完成的数据处理

| 数据集 | 样本 | 状态 | 关键结果 |
|--------|------|------|----------|
| **GSE174574** | 6 个 10x 样本（3 Sham + 3 MCAO, 24 h） | ✅ 分析完成 | 58,287 个细胞；Pirb 在 MCAO 中上调约 3 倍；免疫细胞 Pirb+ 率最高 |
| **GSE171169** | 4 个 10x CD45high 样本（5d + 14d） | ✅ 分析完成 | 10,295 个免疫细胞；Pirb+ 率约 27–28% |
| **GSE227651** | 4 个 10x 样本（1d / 3d / 7d / sham） | ⏳ 重新下载 matrix | 1d 样本完整；3d/7d/sham matrix 损坏，后台下载中 |
| **GSE225948** | 29 个 CSV 样本（Brain + Blood） | ⏳ 重新下载 14 个 counts | 14 个 counts.csv.gz 损坏，后台下载中 |
| **GSE233815** | 57 个 bulk count files | ✅ 已下载 | 为 bulk RNA-seq count，非单细胞；空间子系列未下载 |

---

## 二、GSE174574 主要发现

### 2.1 数据质控
- 原始细胞：58,528
- 质控后：58,287（过滤标准：n_genes ≥ 500, total_counts ≥ 1000, pct_mt < 20%）

### 2.2 Pirb 在 MCAO 后显著上调
| 条件 | Pirb 阳性率 | 平均表达（log1p norm） |
|------|------------|----------------------|
| Sham | 2.7–2.8%   | 0.040–0.043          |
| MCAO | 8.2–9.8%   | 0.119–0.137          |

### 2.3 细胞类型注释与 Pirb 跨细胞类型表达
基于 marker 基因将细胞注释为 7 大类：Astrocyte、Microglia、Immune、Endothelial、Pericyte、Oligodendrocyte、Ependymal。

| 细胞类型 | MCAO Pirb+ % | Sham Pirb+ % | 主要结论 |
|----------|-------------|-------------|----------|
| **Immune** | **31.37%** | 24.24% | Pirb 阳性率最高，提示外周/脑浸润免疫细胞是 Pirb 主要表达群体 |
| **Microglia** | **6.16%** | 4.59% | 小胶质细胞 Pirb 表达显著，MCAO 后升高 |
| **Ependymal** | **4.41%** | 0.88% | MCAO 后明显上调 |
| **Astrocyte** | **3.13%** | 0.66% | 星形胶质细胞 Pirb+ 比例升高约 4.7 倍，支持原报告 |
| Oligodendrocyte | 2.52% | 0.92% | 少量表达 |
| Pericyte | 1.83% | 0.61% | 少量表达 |
| Endothelial | 1.77% | 0.58% | 少量表达 |

**关键洞察**：Pirb 上调不仅限于星形胶质细胞，在免疫细胞和小胶质细胞中表达更高。这提示 Pirb 可能参与缺血后的神经炎症调控。

### 2.4 Pirb+ vs Pirb− 星形胶质细胞差异表达
- 共鉴定 **164 个显著差异基因**（adj. p < 0.05, |log2FC| > 0.5）。
- **上调基因**（Pirb+ vs Pirb−）富含免疫/髓系/溶酶体标志：
  - **Cd14、Ly86**：脂多糖受体复合物，炎症激活标志
  - **Csf1r**：巨噬细胞/小胶质细胞谱系标志
  - **C1qa、C1qb**：补体 C1q，A1 神经毒性星形胶质细胞诱导因子
  - **Spp1（Osteopontin）**：炎症与组织修复相关
  - **Ctsb、Ctsz、Ctss、Ctsd**：溶酶体/组织蛋白酶
  - **Fcer1g、Fcgr3**：Fc 受体，免疫细胞特征

这与原报告提出的"Pirb+ 星形胶质细胞具有免疫/髓系特征"高度一致，可能反映：
1. A1-like 神经毒性星形胶质细胞状态；
2. 免疫细胞 doublet / 环境 RNA 污染（需进一步严格质控验证）。

---

## 三、GSE171169 主要发现

CD45high 浸润免疫细胞数据集（tMCAO 后 5d、14d）：

| 时间 | 细胞数 | Pirb 阳性率 | 平均表达 |
|------|--------|------------|----------|
| 5d   | 5,412  | **28.16%** | 0.378    |
| 14d  | 4,883  | **27.18%** | 0.386    |

- 免疫细胞中 Pirb 高表达且在不同时间点保持稳定，与 GSE174574 中 Immune 细胞结果一致。
- 提示 Pirb 在卒中后外周免疫细胞浸润/驻留中具有持续功能意义。

---

## 四、当前局限与待完成工作

1. **GSE227651** 和 **GSE225948** 的多个原始文件损坏，正在后台重新下载。
2. **细胞类型注释**基于简单 marker 得分，未使用 SingleR/CellTypist 等自动化注释，可作为后续优化。
3. **Pirb+ 星形胶质细胞身份**需要进一步排除 doublet（如 Scrublet/DoubletFinder）和环境 RNA 污染。
4. **机制深挖**：尚未完成配体-受体分析（Nogo-A/MAG/OMgp → Pirb）、通路富集（RhoA-ROCK、NF-κB）、细胞通讯（CellChat/CellPhoneDB）。
5. **GSE233815** 空间转录组数据尚未下载，无法验证 Pirb+ 细胞的空间定位。

---

## 五、下一步计划

1. 等待 GSE227651 / GSE225948 后台下载完成，补全时间序列与免疫/内皮分析。
2. 对 GSE174574 进行 doublet 检测，重新评估 Pirb+ 星形胶质细胞的纯度。
3. 通路富集分析（GO/KEGG/Reactome）Pirb+ vs Pirb− 差异基因。
4. 构建配体-受体互作假设：小胶质细胞 IL-1α/TNF/C1q → 诱导 A1 星形胶质细胞 → Pirb 上调。
5. 若空间数据可获取，验证 Pirb+ 细胞是否富集于梗死周围/半暗带。

---

## 六、输出文件索引

```
04_reports/figures/GSE174574/
├── GSE174574_qc.h5ad
├── GSE174574_annotated.h5ad
├── qc_summary.csv
├── pirb_expression_summary.csv
├── pirb_by_celltype_condition.csv
├── cell_metadata_annotated.csv
├── umap_overview.png
├── umap_celltype.png
├── marker_dotplot.png
├── pirb_violin.png
├── pirb_violin_celltype.png
├── pirb_positive_rate_celltype.png
└── de_pirb/
    ├── DE_PirbPos_vs_Neg_Astrocyte.csv
    ├── DE_PirbPos_vs_Neg_Microglia.csv
    ├── DE_PirbPos_vs_Neg_Immune.csv
    ├── volcano_Astrocyte.png
    ├── volcano_Microglia.png
    └── volcano_Immune.png

04_reports/figures/GSE171169/
├── GSE171169_qc.h5ad
├── qc_summary.csv
├── pirb_summary.csv
├── cell_metadata.csv
└── umap_overview.png
```
