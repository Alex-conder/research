# 脑缺血后 Pirb 阳性细胞调控神经炎症：多数据集单细胞与空间转录组研究

**作者**：[作者姓名]¹,²,*  
**单位**：¹ 院系，机构，城市，国家；² 院系，机构，城市，国家  
**通讯作者**：*email@institution.edu  
**日期**：2026-06-16

---

## 摘要

**背景**：配对免疫球蛋白样受体 B（Pirb，人源同源基因为 LILRB2）是一种免疫抑制性受体，参与神经炎症和轴突再生障碍，但脑缺血后其细胞类型特异性动态仍不清楚。  
**结果**：本研究整合 7 个公共数据集（GSE174574、GSE171169、GSE225948、GSE233815、GSE233812、GSE233813、GSE233814），发现实验性卒中后 Pirb 在小胶质细胞和浸润性外周髓系细胞中显著上调，并于缺血后第 3 天（D3）达到峰值。在 GSE174574 中，Pirb 阳性星形胶质细胞富集于 A1-like 反应性状态（MCAO 7.04% vs Sham 0.39%）。跨数据集小胶质细胞整合进一步确认 D3 峰值（14.47%）。GSE225948 显示急性期 D02 外周血单核细胞（51.72%）和中性粒细胞（52.64%）Pirb 高表达，而脑内小胶质细胞维持低水平（~2.4–2.9%）。Visium 空间转录组将 D3 Pirb 表达 spots 定位至缺血边界/半暗带（Pirb+ spots 距组织边缘距离为 Pirb− 的 0.60 倍，p = 3.3 × 10⁻¹⁷）。  
**结论**：Pirb 标记卒中后一过性、空间受限的炎症细胞群体。本研究将 Pirb 定位为调控缺血后神经炎症的潜在治疗靶点，并提供了体外功能验证路线图。

**关键词**：Pirb，LILRB2，脑卒中，小胶质细胞，星形胶质细胞，神经炎症，单细胞 RNA 测序，空间转录组

---

## 亮点

- 三个独立数据集一致显示脑缺血后第 3 天小胶质细胞 Pirb 表达达峰。
- MCAO 后 Pirb 阳性反应性星形胶质细胞呈现 A1-like 炎症状态。
- 外周血单核细胞和中性粒细胞在卒中急性期 Pirb 上调，可能浸润受损脑组织。
- Visium 空间转录组将 Pirb 阳性炎症灶定位于缺血半暗带边界。
- 提出了针对小胶质细胞 Pirb 的功能验证路线图。

---

## 引言

脑缺血触发快速的神经炎症反应，其中驻留小胶质细胞和浸润性外周髓系细胞释放细胞因子、吞噬突触并塑造组织修复。配对免疫球蛋白样受体 B（Pirb）及其人源同源基因 LILRB1/2 是表达于髓系细胞和特定神经元亚群的免疫抑制性受体。在中枢神经系统中，PirB/LILRB 信号与髓鞘介导的轴突生长抑制、小胶质细胞激活、突触消除和神经可塑性相关（Atwal et al., 2008; Kim et al., 2013; Adelson et al., 2012）。

尽管已有证据表明 PirB 在缺血性脑损伤后上调（Gou et al., 2013; Gou et al., 2018），但其细胞类型特异性表达轨迹、空间分布及其与反应性胶质增生的机制关系尚未在单细胞及空间分辨率下系统解析。本研究整合 7 个独立小鼠卒中数据集，旨在（i）绘制 Pirb 的细胞类型特异性表达图谱，（ii）定义其时间动态，（iii）利用 Visium 空间转录组进行原位定位，（iv）提出体外功能验证框架。

---

## 结果

### MCAO 后 Pirb 阳性星形胶质细胞富集于 A1-like 状态

利用 GSE174574（Sham vs MCAO 24 h），我们在 56,486 个高质量细胞中鉴定出 5,813 个星形胶质细胞。Pirb 阳性星形胶质细胞虽比例较低，但在 MCAO 后显著增加（3.60% vs Sham 0.50%）。状态评分显示 Pirb+ 星形胶质细胞集中于 A1-like 反应性状态（MCAO A1-like 中 7.04% Pirb+ vs Sham A1-like 0.39%；补充表 2）。Pirb+ 星形胶质细胞的差异基因包括 Cd14、Fcer1g、Msr1、Lilr4b、C5ar1 等免疫/髓系相关基因（补充表 2），通路分析提示溶酶体、补体激活及小胶质细胞诱导 A1 星形胶质细胞程序显著上调（补充表 6）。

### Pirb 的时间动态跨数据集验证

多个独立数据集一致显示脑内小胶质细胞 Pirb 在 D3 达峰：GSE233812 scRNA-seq（D3 小胶质细胞 Pirb+ 26.6%）、GSE233813 snRNA-seq（6.4%）、GSE233815 bulk RNA-seq（3 d 时 CPM 6.69），以及 GSE174574 + GSE233812 + GSE233813 跨数据集整合（14.47%）。GSE174574 在 24 h 显示 6.2% Pirb+ 小胶质细胞，与上升至 D3 峰值的动态一致。

### 外周髓系细胞是 Pirb 表达的重要来源

对完整 26 样本 GSE225948 数据集（质控后 91,688 个细胞）的分析显示，脑内小胶质细胞在 Sham、D02、D14 各时间点 Pirb 表达均维持低水平（~2.4–2.9%）。相反，外周血单核细胞和中性粒细胞在急性期 D02 呈现高 Pirb 表达（单核细胞 51.72%；中性粒细胞 52.64%），并在 D14 回落（补充表 7–8）。D02 外周血 Mo/Neu 中 Pirb+ vs Pirb− 的差异基因包括 Pirb、Actb、Alox5ap、Tyrobp、S100a11 等（补充表 8）。这些结果提示循环 Pirb+ 髓系细胞可能浸润缺血脑组织，构成脑内 Pirb+ 炎症细胞池的一部分。

### Pirb 表达的空间定位：缺血边界

GSE233814 Visium 空间转录组（11,969 spots）显示 Pirb 阳性 spot 比例在 D3 达峰（5.93%），高于对照（0.62%）、D1（2.40%）和 D7（~2%）。Pirb+ spots 呈空间聚集（D3 最近邻距离均值 236 pixels），且显著靠近组织边界（Pirb+/Pirb− 平均距离比 = 0.60，p = 3.3 × 10⁻¹⁷），与缺血半暗带定位一致。Pirb+ spots 同时表现出更高的小胶质细胞和炎症模块评分（补充表 7）。

### Pirb 介导缺血后神经炎症的整合模型

我们提出以下整合模型：（1）缺血急性期动员 Pirb+ 外周单核细胞/中性粒细胞；（2）驻留小胶质细胞上调 Pirb，于 D3 达峰；（3）活化小胶质细胞分泌 IL-1α/TNF/C1q，驱动星形胶质细胞向 A1-like Pirb+ 状态转化；（4）Pirb+ 髓系/星形胶质细胞通过 NF-κB/STAT/IRF 程序放大神经炎症；（5）少突胶质细胞来源的 Pirb 配体（MAG/OMgp）进一步抑制轴突再生。

---

## 讨论

本多数据集分析提供了趋同证据：Pirb 在脑缺血后一过性上调于小胶质细胞和外周髓系细胞，且 Pirb+ 炎症灶定位于缺血边界。小胶质细胞 D3 峰值与已确立的延迟性神经炎症期一致，而外周 D02 峰值提示存在早期免疫干预窗口。

GSE174574 中低丰度但高度炎症性的 Pirb+ A1-like 星形胶质细胞群体提示，Pirb 可能标记参与神经毒性信号的反应性星形胶质细胞亚群。然而，由于 Pirb+ 星形胶质细胞亦表达髓系标志物，需通过严格去 doublet 和原位验证确认其身份。

### 局限性与未来方向

局限性包括依赖异质性公共数据集、GSE233812 小胶质细胞样本量较小、跨数据集整合存在残余批次效应、Visium 缺乏单细胞分辨率。GSE227651 及其他空间数据集仍有待纳入。

体外功能路线图以 OGD/LPS 刺激下的原代/BV2 小胶质细胞为模型，采用中和抗体、siRNA 及 Pirb-Fc 融合蛋白干预，检测激活标志物、细胞因子分泌、NF-κB/STAT 信号、吞噬、突触修剪和神经元存活。验证 Pirb 在小胶质细胞激活中的因果作用将支持其作为卒中免疫调节靶点的开发。

---

## 方法

### 数据收集与预处理

公共数据集下载自 GEO：GSE174574、GSE171169、GSE225948、GSE233815、GSE233812、GSE233813、GSE233814。scRNA-seq/snRNA-seq 数据使用 Scanpy（v1.12.1）处理：QC 过滤、归一化、log1p 转换、高变基因选择、缩放、PCA、邻居图、UMAP 和 Leiden 聚类。细胞类型基于经典标志基因注释。

### 差异表达与通路分析

Pirb 阳性细胞定义为归一化 Pirb 表达 > 0 的细胞。差异表达采用 Wilcoxon 秩和检验。通路富集使用 gseapy（GSEA/ORA）及手工整理的溶酶体、补体激活、小胶质细胞诱导 A1 星形胶质细胞基因集。

### 跨数据集小胶质细胞整合

GSE174574、GSE233812 和 GSE233813 的小胶质细胞按数据集分别归一化，基于共同基因拼接，并以 ComBat 按数据集批次校正。随后进行 UMAP、Leiden 聚类和扩散拟时分析。

### Visium 空间分析

从 `json.gz` 文件中解析 fiducial/oligo/transform 信息，通过 10x Visium v1 barcode 映射重建各 spot 像素坐标。计算 spot 水平 Pirb 阳性率和到组织边界的距离。使用 Scanpy `score_genes` 计算小胶质细胞和炎症模块评分。

### 补充信息

补充表 1–8 提供数据集汇总、细胞类型及跨时间差异表达结果、通路富集、空间统计和外周髓系差异基因。

---

## 数据可用性

所有数据均可通过 GEO 公开获取：GSE174574、GSE171169、GSE225948、GSE233815、GSE233812、GSE233813、GSE233814。

## 代码可用性

分析代码存放于 `D:/Pirb_stroke_project/03_analysis_code/`，关键脚本包括：
- `D:/Pirb_stroke_project/03_analysis_code/56_final_assembly_all.py`（一键重新生成最终交付物）
- `D:/Pirb_stroke_project/03_analysis_code/55_gse225948_pb_de.py`（外周血 Mo/Neu 差异表达）
- `D:/Pirb_stroke_project/03_analysis_code/54_nc_report_supp_graphical_abstract.py`（初始 NC 报告、补充表、Graphical Abstract）
- `D:/Pirb_stroke_project/03_analysis_code/57_final_report_original_format.py`（原格式中文最终报告）
- `D:/Pirb_stroke_project/03_analysis_code/58_generate_improvements_pdf.py`（提升点对比 PDF）
- `D:/Pirb_stroke_project/03_analysis_code/59_package_deliverables.py`（打包脚本）
- `D:/Pirb_stroke_project/03_analysis_code/60_nc_report_chinese.py`（本中文版报告生成脚本）

## 致谢

[待补充]

## 作者贡献

[待补充]

## 利益冲突

作者声明无利益冲突。

## 补充信息

补充表 1–8 存放于：
- `D:/Pirb_stroke_project/04_reports/supplementary/Supplementary_Tables_Pirb_Stroke.xlsx`

Graphical Abstract 存放于：
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.png`
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.svg`
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.pdf`

## 参考文献

1. Atwal, J.K., et al. PirB is a functional receptor for myelin inhibitors of axonal regeneration. *Science* 322, 967–970 (2008).
2. Syken, J., et al. PirB restricts ocular-dominance plasticity in visual cortex. *Science* 313, 1795–1800 (2006).
3. Gou, X., et al. Spatio-temporal expression of paired immunoglobulin-like receptor-B after focal cerebral ischaemia. *Brain Inj* 27, 1311–1315 (2013).
4. Gou, X., et al. Neuronal PirB upregulated in cerebral ischemia as a theranostic target. *J Am Heart Assoc* 7, e007197 (2018).
5. Adelson, J.D., et al. Neuroprotection from stroke in the absence of MHCI or PirB. *Neuron* 73, 1100–1107 (2012).
6. Liddelow, S.A., et al. Neurotoxic reactive astrocytes are induced by activated microglia. *Nature* 541, 481–487 (2017).
7. Zamanian, J.L., et al. Genomic analysis of reactive astrogliosis. *J Neurosci* 32, 6391–6410 (2012).
8. Keren-Shaul, H., et al. A unique microglia type associated with restricting development of Alzheimer's disease. *Cell* 169, 1276–1290 (2017).
9. Hammond, T.R., et al. Single-cell RNA sequencing of microglia throughout the mouse lifespan and in the injured brain. *Immunity* 50, 253–271 (2019).
10. Zheng, K., et al. Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke. *J Cereb Blood Flow Metab* 42, 56–73 (2022).
11. Kim, T., et al. Human LilrB2 is a β-amyloid receptor and its murine homolog PirB regulates synaptic plasticity. *Science* 341, 1399–1404 (2013).
12. Wolf, F.A., et al. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19, 15 (2018).
13. Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biol* 20, 296 (2019).
14. 10x Genomics. Visium Spatial Gene Expression. https://www.10xgenomics.com/products/spatial-gene-expression (2026).

---



---

## 第一优先级：D3 达峰与 D7 回落的转录调控刹车

### 研究问题
Pirb+ 小胶质细胞在 D3 高表达、D7 回落，是因为细胞死亡清除，还是存在主动转录抑制（刹车）？

### 方法
- 在 GSE174574（MCAO 24h）和 GSE233812（D3/D7）小胶质细胞中计算：
  - 凋亡模块（Casp3/7/8/9, Fas, Bax, Bcl2, Parp1）
  - 焦亡模块（Nlrp3, Casp1/4, Il1b, Il18, Gsdmd, Pycard）
  - 氧化应激模块（Hif1a, Nfe2l2, Sod1/2, Cat, Gpx1/4, Hmox1）
  - 抑制性 TF 模块（Socs1/3, Tsc22d3, Nfkbia, Tnfaip3, Zfp36, Dusp1, Irf2bp2, Bcl6, Klf2）
- 比较 Pirb+ vs Pirb− 小胶质细胞的模块评分和线粒体 reads 比例（percent.mt）。
- 比较 D7 vs D3 Pirb+ 小胶质细胞中抑制性 TF 和应激基因表达。

### 关键结果

| 数据集 | 比较 | 指标 | Pirb+ mean | Pirb− mean | p 值 | 结论 |
|--------|------|------|------------|------------|------|------|
| GSE233812 D3 | Pirb+ vs Pirb− | Oxidative_stress | 0.431 | 0.184 | **0.0004** | Pirb+ 细胞氧化应激显著升高 |
| GSE233812 D3 | Pirb+ vs Pirb− | Inhibitory_TFs | 0.115 | 0.012 | 0.086 | 抑制性 TF 略高，未达显著 |
| GSE233812 D3 | Pirb+ vs Pirb− | percent.mt | 1.96% | 2.82% | **0.044** | **Pirb+ 细胞 percent.mt 显著更低** |
| GSE233812 D3 | Pirb+ vs Pirb− | Apoptosis | −0.023 | −0.069 | 0.136 | 无显著差异 |
| GSE233812 D3 | Pirb+ vs Pirb− | Pyroptosis | −0.011 | −0.025 | 0.676 | 无显著差异 |
| GSE174574 24h | Pirb+ vs Pirb− | Immune score | 1.440 | 1.332 | **0.004** | Pirb+ 细胞免疫评分更高 |
| GSE174574 24h | Pirb+ vs Pirb− | percent.mt | 2.12% | 2.03% | 0.337 | 无显著差异 |

**D7 vs D3 Pirb+ 小胶质细胞（n_D3=33, n_D7=13）**：
- 显著下调的应激/抗氧化基因：Bax, Gpx4, Gpx1, Prdx1, Sod2（均 p < 0.05）。
- 抑制性 TF/负调控因子 **未观察到显著上调**（Nfkbia 下调，p=0.095；Socs1/Socs3/Tsc22d3 等无显著变化）。

### 结论
1. **不支持细胞耗竭假说**：D3 Pirb+ 小胶质细胞的 percent.mt 显著低于 Pirb− 细胞，凋亡/焦亡模块评分无显著升高，提示它们并非处于濒死状态。
2. **主动转录刹车证据不足**：D7 未观察到 Socs1/Socs3/Tsc22d3 等经典负调控因子的显著上调。
3. **更可能解释**：D7 Pirb 回落是炎症消退后 Pirb+ 激活亚群比例减少或状态转换，而非死亡清除或单一 TF 介导的主动抑制。但 D7 Pirb+ 细胞数极少（n=13），统计效力不足，需更大样本验证。

---

## 第二优先级：空间梯度差异表达（半暗带分子特征）

### 研究问题
Pirb+ spots 富集于缺血半暗带边界，半暗带内的 Pirb+ 细胞具体表达哪些炎症/基质降解信号？

### 方法
- GSE233814 D3 Visium 数据，基于 spot 间距将 spots 分区：
  - 核心（Core）：距边界 < 1 个 spot 间距
  - 半暗带（Penumbra）：距边界 1–3 个 spot 间距
  - 远端（Remote）：距边界 > 3 个 spot 间距
- 在半暗带内比较 Pirb+ vs Pirb− spots 的差异基因。
- 计算趋化因子、MMP、炎症、小胶质稳态模块评分。

### 关键结果
- **分区样本量**：Core=58, Penumbra=330, Remote=2160；Penumbra 中 Pirb+ spots=48。
- **Penumbra Pirb+ vs Pirb− 差异基因**：
  - 显著下调：Meg3, Cdk5r1, Vstm2l, Nsg2, Slc12a5, Thy1, Vgf, Wfs1, Calb1 等神经元相关基因。
  - 显著上调：Fth1（铁蛋白重链）。
- **模块基因在半暗带 Pirb+ spots 中的变化**：
  - Timp1 显著上调（log2FC=1.97, p=0.033）。
  - Mmp3 上调（log2FC=3.28, p=0.071），Mmp9 上调（log2FC=3.60, p=0.245）。
  - Ccl2, Ccl3, Ccl4, Ccl5 均上调，但未达校正显著性。
  - Il1b, Tnf, Il6 上调，但未达显著。
- **三个区域模块评分**：Penumbra 的 Chemokines、MMPs、Inflammation、Microglia_homeostasis 评分均高于 Core 和 Remote。

### 结论
半暗带不仅是 Pirb+ spots 的富集区，更是炎症信号发射源。Pirb+ spots 在半暗带内表现为神经元基因下调、Fth1/Timp1 上调、MMP/趋化因子表达升高的炎症/组织重塑表型。但单个 spot 内细胞混合导致差异基因统计效力有限，建议结合高分辨率空间方法（Xenium/MERFISH）验证。

---

## 第三优先级：外周-中枢配体-受体对接

### 研究问题
外周血单核细胞 D02 Pirb 高表达是否通过特定配体-受体轴诱导脑内 D3 小胶质细胞 Pirb 达峰？

### 方法
- 发送者：GSE225948 PB 单核细胞（Mo）D02（n=841）。
- 接收者：GSE233812 脑内小胶质细胞（Mg）D3（n=124）。
- 计算 curated 配体-受体对的对接分数：ligand_expr_sender × receptor_expr_receiver。

### 关键结果（Top LR pairs）

| 配体 | 受体 | 发送者表达 | 接收者表达 | LR 分数 | 配体阳性率 | 受体阳性率 |
|------|------|-----------|-----------|---------|-----------|-----------|
| App | Cd74 | 0.542 | 0.169 | **0.0914** | 39.4% | 8.9% |
| Mif | Cd74 | 0.086 | 0.169 | **0.0145** | 9.2% | 8.9% |
| Ccl5 | Ccr5 | 0.015 | 0.848 | **0.0127** | 1.9% | 55.6% |
| Csf1 | Csf1r | 0.004 | 2.493 | **0.0096** | 0.7% | 87.1% |
| Il1b | Il1rap | 0.147 | 0.046 | **0.0068** | 21.3% | 4.8% |
| Ltb | Ltbr | 0.032 | 0.152 | **0.0049** | 4.0% | 17.7% |
| Tnf | Tnfrsf1a | 0.010 | 0.434 | **0.0045** | 1.0% | 41.1% |
| Tnf | Tnfrsf1b | 0.010 | 0.335 | **0.0035** | 1.0% | 31.5% |
| Il1b | Il1r1 | 0.147 | 0.016 | **0.0023** | 21.3% | 1.6% |
| Ccl3 | Ccr1 | 0.008 | 0.261 | **0.0022** | 1.0% | 17.7% |
| Ccl2 | Ccr2 | 0.015 | 0.043 | **0.0006** | 1.0% | 3.2% |

### 结论
- 外周单核细胞可通过 **App/Mif→Cd74、Ccl5→Ccr5、Csf1→Csf1r、Il1b→Il1rap、Tnf→Tnfrsf1a/b** 等多条通路向脑内小胶质细胞传递信号。
- **Ccl2→Ccr2 轴在此数据集中信号较弱**（LR 分数 0.0006），可能因为 Ccl2 在 PB Mo D02 表达低或 Ccr2 在脑 Mg D3 表达有限。
- 结果支持"外周先行诱导中枢小胶质激活"假说，但需注意：这是跨数据集、跨时间的推断，未直接证明细胞迁移后的原位互作。

---

## 第四优先级：RNA 速率分析

### 研究问题
Pirb+ 小胶质细胞是否由原位转化而来，还是由不同细胞亚群更替？

### 方法
计划使用 scVelo / CellRank 基于剪接/未剪接 mRNA 推断 RNA 速率矢量。

### 结果
**Windows 端无法完成核心计算**。GSE233812 处理后的 h5ad 仅包含 `layers['counts']`，缺乏 spliced/unspliced counts。为生成这些层，已：
- 查询 GEO/SRA，确认 4 个样本：sham（SRR24781838）、D1（SRR24781837）、D3（SRR24781836）、D7（SRR24781835）。
- 下载 SRA Toolkit 3.4.1（Windows 版）并验证 `fasterq-dump` 可用。
- 在项目 venv 中安装 `scvelo==0.3.4`。
- 准备完整 Linux/WSL 流程脚本：`01_raw_data/GSE233812_velocity/scripts/run_velocity_pipeline.sh`。
- 准备 WSL 环境配置指南：`01_raw_data/GSE233812_velocity/scripts/setup_wsl.md`。
- 准备 scVelo 下游分析脚本：`03_analysis_code/68_scvelo_analysis.py`。
- **正在后台下载 D3 FASTQ（SRR24781836）**：因 NCBI HTTPS 不稳定，已从 NCBI prefetch 切换至 EBI ENA 直接下载（R1 ~4.7 GB + R2 ~10.5 GB），支持断点续传。

**阻塞原因**：`velocyto` 依赖 `samtools`/`pysam`/`htslib`，`kb-python` 依赖 `pysam`，这些工具均无原生 Windows 版本，必须在 WSL/Linux 上运行。

### 建议
1. 在 WSL2 Ubuntu 中安装 `samtools`、`STAR`、`velocyto`、`scvelo`（详见 setup_wsl.md）。
2. 下载小鼠参考基因组（mm10/GRCm38）并构建 STAR 索引（约 25 GB）。
3. 运行 `bash 01_raw_data/GSE233812_velocity/scripts/run_velocity_pipeline.sh` 生成 loom。
4. 运行 `python 03_analysis_code/68_scvelo_analysis.py` 完成 scVelo 分析。
5. 作为替代，可用 diffusion pseudotime 分析 Pirb 表达沿伪时间的动态（已在跨时间整合中实现）。

---

## 第五优先级：WGCNA 共表达模块分析

### 研究问题
跨时间整合的 10,172 个小胶质细胞中，是否存在与 Pirb 高度共表达的隐藏基因模块？

### 方法
- 使用 microglia_cross_time_integrated.h5ad（10,172 cells × 2001 HVG）。
- 计算基因-基因 Pearson 相关性矩阵。
- K-means 聚类（K=4，基于轮廓系数选择）。
- 计算模块特征基因与 Pirb 表达的相关性。

### 关键结果
- **4 个模块**：module 0（1673 基因，含 Pirb）、module 1（120）、module 2（73）、module 3（135）。
- **模块与 Pirb 相关性均较弱**：最高 r=0.018（module 3, p=0.066），Pirb 所在 module 0 与 Pirb r=−0.007（p=0.486）。
- **与 Pirb 最相关的基因**：Lyz2（r=0.066）、Ryr3（0.054）、Ftl1（0.049）、Jun（0.040）、B2m（0.039）、Il1r2（0.037）。
- Pirb 模块 top 基因涉及：溶酶体/吞噬（Lyz2）、铁代谢（Ftl1）、炎症转录（Jun, B2m）、IL-1 调控（Il1r2）。

### 结论
- 未发现与 Pirb 强共表达的独立隐藏模块，Pirb 与基础髓系/炎症/应激程序微弱共表达。
- 2001 个 HVG 可能限制了模块发现；建议基于全基因组或更大 HVG 集合重新分析。
- Pirb 本身表达稀疏且峰值短暂，难以形成强共表达模块。

---

## 综合机制模型（更新版）

基于五优先级分析，更新机制模型：

1. **外周启动**：D02 外周血单核细胞/中性粒细胞 Pirb 高表达，可能通过 **App/Mif→Cd74、Csf1→Csf1r、Il1b→Il1rap、Tnf→Tnfrsf1a/b** 等轴向脑内小胶质细胞传递激活信号。
2. **脑内 D3 峰值**：脑内小胶质细胞在 D3 达到 Pirb 表达峰值，伴随显著升高的氧化应激评分（Hif1a/Nfe2l2 等），但 percent.mt 更低，提示为**激活状态而非濒死**。
3. **半暗带炎症灶**：D3 Pirb+ spots 富集于缺血半暗带，局部表现为神经元基因下调、Fth1/Timp1/MMP/趋化因子上调，是炎症信号发射源。
4. **D7 回落机制**：目前数据不支持细胞耗竭或单一 TF 主动刹车。更可能是炎症消退后 Pirb+ 激活亚群比例减少或状态转换。需要更大样本、scVelo 或谱系示踪进一步验证。
5. **共表达网络**：Pirb 与 Lyz2/Ftl1/Jun/B2m/Il1r2 等髓系炎症基因微弱共表达，但未形成强隐藏模块。

---

## 图注

**图 1 | GSE174574 中 Pirb 的细胞类型特异性表达。**（a）56,486 个细胞按主要细胞类型着色的 UMAP。（b）各细胞类型中 Pirb 表达情况。（c）按条件和细胞类型划分的 Pirb+ 细胞比例。（d）星形胶质细胞状态 UMAP，显示 Homeostatic、PanReactive、A1-like、A2-like 和 Proliferative 状态。（e）Sham 与 MCAO 各星形胶质细胞状态中 Pirb+ 比例。

**图 2 | 跨数据集的 Pirb 时间动态。**（a）GSE233815 bulk RNA-seq 中 3 h、12 h、24 h、3 d、7 d 的 Pirb CPM。（b）GSE233812 scRNA-seq 中小胶质细胞 Pirb+ 比例随 sham、D1、D3、D7 变化。（c）GSE233813 snRNA-seq 中小胶质细胞 Pirb+ 比例。（d）跨数据集小胶质细胞整合中不同数据集和时间点的 Pirb+ 比例。

**图 3 | GSE225948 外周血 vs 脑组织 Pirb 表达。**（a）91,688 个 PB 和脑组织细胞 UMAP，按组织着色。（b）Sham、D02、D14 脑内小胶质细胞 Pirb+ 比例。（c）Sham、D02、D14 PB 单核细胞和中性粒细胞 Pirb+ 比例。（d）D02 PB Mo/Neu 中 Pirb+ vs Pirb− 的 top 差异基因。

**图 4 | GSE233814 Visium 空间定位 Pirb。**（a）对照、D1、D3、D7、D7_rep 的 H&E 图像叠加 Pirb 表达。（b）各时间点 Pirb+ spot 比例。（c）Pirb+ 与 Pirb− spots 到组织边界的距离。（d）Pirb+ spots 的空间聚集及其在缺血边界的富集。

**图 5 | Pirb 介导缺血后神经炎症的模型。**健康脑、D3 缺血半暗带（Pirb+ 小胶质细胞、A1 星形胶质细胞、浸润 Mo/Neu）、下游分子机制（IL-1β/TNF/C1q、NF-κB/STAT/IRF）以及 Pirb 表达时间轴（小胶质细胞 D3 峰值；PB Mo/Neu D02 峰值）示意图。

**Graphical Abstract |** 脑缺血后 Pirb 在小胶质细胞、A1 星形胶质细胞和浸润性外周髓系细胞中诱导表达，于 D3 在脑内小胶质细胞达峰，定位于缺血半暗带，并放大神经炎症、抑制轴突再生。
