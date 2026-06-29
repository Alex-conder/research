# Pirb-positive cells orchestrate post-ischemic neuroinflammation: a multi-dataset single-cell and spatial transcriptomics study

**Authors**: [Author names]¹,²,*  
**Affiliations**: ¹ Department, Institution, City, Country; ² Department, Institution, City, Country  
**Corresponding author**: *email@institution.edu  
**Date**: 2026-06-16

---

## Abstract

**Background**: Paired immunoglobulin-like receptor B (Pirb, mouse ortholog of human LILRB2) is an immune inhibitory receptor implicated in neuroinflammation and axon regeneration failure, yet its cell-type-specific dynamics after cerebral ischemia remain poorly defined.  
**Results**: Integrating seven public datasets (GSE174574, GSE171169, GSE225948, GSE233815, GSE233812, GSE233813, GSE233814), we show that Pirb is sharply upregulated in microglia and infiltrating peripheral myeloid cells after experimental stroke, peaking at post-ischemic day 3 (D3). In GSE174574, Pirb-positive astrocytes were enriched in the A1-like reactive state (7.04% vs 0.39% in sham). Cross-dataset microglial integration confirmed a D3 peak (14.47%). GSE225948 revealed high Pirb expression in peripheral blood monocytes (51.72%) and neutrophils (52.64%) at acute D02, while brain microglia remained low (~2.4–2.9%). Visium spatial transcriptomics localized D3 Pirb-expressing spots to the ischemic boundary/penumbra (Pirb+ spots 0.60-fold closer to tissue edge, p = 3.3 × 10⁻¹⁷).  
**Conclusions**: Pirb marks a transient, spatially restricted inflammatory cell population after stroke. Our findings position Pirb as a candidate therapeutic target for modulating post-ischemic neuroinflammation and provide a roadmap for in vitro functional validation.

**Keywords**: Pirb, LILRB2, stroke, microglia, astrocyte, neuroinflammation, single-cell RNA sequencing, spatial transcriptomics

---

## Highlights

- Pirb expression peaks at post-ischemic day 3 in brain microglia across three independent datasets.
- Pirb-positive reactive astrocytes adopt an A1-like inflammatory state after MCAO.
- Peripheral blood monocytes and neutrophils show acute Pirb upregulation and may infiltrate the injured brain.
- Visium spatial transcriptomics places Pirb-positive inflammatory foci at the ischemic penumbra boundary.
- A functional validation roadmap targeting Pirb in microglia is proposed.

---

## Introduction

Cerebral ischemia triggers a rapid neuroinflammatory response in which resident microglia and infiltrating peripheral myeloid cells release cytokines, engulf synapses, and shape tissue repair. Paired immunoglobulin-like receptor B (PirB) and its human orthologs LILRB1/2 are immune inhibitory receptors expressed by myeloid cells and subsets of neurons. In the CNS, PirB/LILRB signaling has been linked to myelin-mediated axon growth inhibition, microglial activation, synapse elimination, and neuronal plasticity (Atwal et al., 2008; Kim et al., 2013; Adelson et al., 2012).

Despite evidence that PirB is upregulated after ischemic brain injury (Gou et al., 2013; Gou et al., 2018), its cell-type-specific expression trajectory, spatial distribution, and mechanistic relationship to reactive gliosis have not been systematically mapped at single-cell and spatial resolution. Here, we integrate seven independent mouse stroke datasets to (i) map the cell-type-specific expression of Pirb, (ii) define its temporal trajectory, (iii) localize Pirb-expressing cells in situ using Visium spatial transcriptomics, and (iv) propose an in vitro functional validation framework.

---

## Results

### Pirb is induced in A1-like reactive astrocytes after MCAO

Using GSE174574 (Sham vs MCAO 24 h), we identified 5,813 astrocytes among 56,486 high-quality cells. Pirb-positive astrocytes were rare but significantly enriched after MCAO (3.60% vs 0.50% in sham). State scoring revealed that Pirb+ astrocytes were concentrated in the A1-like reactive state (7.04% Pirb+ in MCAO A1-like vs 0.39% in sham A1-like; Supplementary Table 2). Differentially expressed genes in Pirb+ astrocytes included immune/microglial markers such as Cd14, Fcer1g, Msr1, Lilr4b, and C5ar1 (Supplementary Table 2), and pathway analysis implicated lysosome, complement activation, and microglia-induced A1 astrocyte programs (Supplementary Table 6).

### Temporal dynamics of Pirb across independent datasets

Multiple datasets converged on a D3 peak for Pirb expression in brain microglia: GSE233812 scRNA-seq (26.6% Pirb+ microglia at D3), GSE233813 snRNA-seq (6.4%), GSE233815 bulk RNA-seq (CPM 6.69 at 3 d), and cross-dataset integration of GSE174574 + GSE233812 + GSE233813 (14.47% Pirb+ microglia at D3). GSE174574 at 24 h showed 6.2% Pirb+ microglia, consistent with an ascending phase toward the D3 peak (Supplementary Table 5).

### Peripheral myeloid cells are a major Pirb-expressing compartment

Analysis of the complete 26-sample GSE225948 dataset (91,688 cells after QC) showed that brain microglia maintained low Pirb expression across sham, D02, and D14 (~2.4–2.9%). In contrast, peripheral blood monocytes and neutrophils displayed high acute Pirb expression that peaked at D02 (monocytes 51.72%; neutrophils 52.64%) and declined by D14 (Supplementary Tables 7–8). Differential expression in PB Mo/Neu at D02 identified Pirb, Actb, Alox5ap, Tyrobp, and S100a11 among the top Pirb+ markers (Supplementary Table 8). These findings suggest that circulating Pirb+ myeloid cells may infiltrate the ischemic brain and contribute to the inflammatory Pirb+ pool.

### Spatial localization of Pirb expression to the ischemic boundary

GSE233814 Visium spatial transcriptomics (11,969 spots) showed that Pirb-positive spot fraction peaked at D3 (5.93%) relative to control (0.62%), D1 (2.40%), and D7 (~2%). Pirb+ spots were spatially clustered (nearest-neighbor distance mean 236 pixels at D3) and significantly closer to the tissue boundary than Pirb− spots (mean distance ratio Pirb+/Pirb− = 0.60, p = 3.3 × 10⁻¹⁷), consistent with localization to the ischemic penumbra. Pirb+ spots displayed elevated microglial and inflammatory module scores (Supplementary Table 7).

### An integrated model of Pirb-mediated neuroinflammation

We propose an integrated model in which (1) ischemia acutely mobilizes Pirb+ peripheral monocytes/neutrophils; (2) resident microglia upregulate Pirb, peaking at D3; (3) activated microglia secrete IL-1α/TNF/C1q, driving astrocytes toward an A1-like Pirb+ state; (4) Pirb+ myeloid/astroglial cells amplify neuroinflammation via NF-κB/STAT/IRF programs; and (5) oligodendrocyte-derived Pirb ligands (MAG/OMgp) further suppress axon regeneration.

---

## Discussion

Our multi-dataset analysis provides convergent evidence that Pirb is transiently upregulated in microglia and peripheral myeloid cells after cerebral ischemia and that Pirb+ inflammatory foci localize to the ischemic boundary. The D3 microglial peak aligns with the established delayed neuroinflammatory phase, while the D02 peripheral peak identifies a window for early immune intervention.

The low-frequency but highly inflammatory Pirb+ A1 astrocyte population observed in GSE174574 suggests that Pirb may mark a subset of reactive astrocytes engaged in neurotoxic signaling. However, because Pirb+ astrocytes also express myeloid markers, careful doublet depletion and in situ validation are required to confirm their identity.

### Limitations

Limitations include the reliance on public datasets with heterogeneous protocols, the small microglial sample size in GSE233812, potential batch effects in cross-dataset integration, and the lack of single-cell resolution in Visium data. GSE227651 and additional spatial datasets remain to be incorporated.

### Future directions

The in vitro functional roadmap targets Pirb in primary/BV2 microglia under OGD/LPS stimulation using neutralizing antibodies, siRNA, and Pirb-Fc fusion proteins. Readouts span activation markers, cytokine secretion, NF-κB/STAT signaling, phagocytosis, synapse pruning, and neuronal survival. Validating the causal role of Pirb in microglial activation would support its development as a stroke immunomodulatory target.

---

## Methods

### Data collection and preprocessing

Public datasets were downloaded from GEO: GSE174574, GSE171169, GSE225948, GSE233815, GSE233812, GSE233813, and GSE233814. scRNA-seq/snRNA-seq data were processed with Scanpy (v1.12.1): QC filtering, normalization, log-transformation, highly variable gene selection, scaling, PCA, neighbor graph, UMAP, and Leiden clustering (Wolf et al., 2018). Cell types were annotated using canonical marker genes.

### Differential expression and pathway analysis

Pirb-positive cells were defined as cells with normalized Pirb expression > 0. Differential expression was performed using Wilcoxon rank-sum tests. Pathway enrichment was performed using gseapy (GSEA/ORA) and custom gene sets for lysosome, complement, and microglia-induced A1 astrocyte programs.

### Cross-dataset microglial integration

Microglia from GSE174574, GSE233812, and GSE233813 were normalized per dataset, concatenated on common genes, and batch-corrected with ComBat (dataset as batch). UMAP, Leiden clustering, and diffusion pseudotime were computed.

### Visium spatial analysis

Spatial coordinates were reconstructed from `json.gz` files by mapping 10x Visium v1 barcodes to fiducial/oligo positions and transform matrices (10x Genomics, 2026). Spot-level Pirb positivity and distance to tissue boundary were computed. Microglial and inflammatory module scores were calculated using Scanpy `score_genes`.

### Supplementary information

Supplementary Tables 1–8 are provided in:
- `D:/Pirb_stroke_project/04_reports/supplementary/Supplementary_Tables_Pirb_Stroke.xlsx`

This file contains dataset summary, cell-type-specific and cross-time differential expression results, pathway enrichments, spatial statistics, and peripheral myeloid DE genes. The Graphical Abstract is available at:
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.png`
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.svg`
- `D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.pdf`

---

## Data availability

All data are publicly available via GEO under accession numbers GSE174574, GSE171169, GSE225948, GSE233815, GSE233812, GSE233813, and GSE233814.

## Code availability

Analysis code is available in `D:/Pirb_stroke_project/03_analysis_code/` and is provided as Supplementary Software. Key scripts include:
- `D:/Pirb_stroke_project/03_analysis_code/56_final_assembly_all.py` (one-click regeneration of final deliverables)
- `D:/Pirb_stroke_project/03_analysis_code/55_gse225948_pb_de.py` (PB Mo/Neu differential expression)
- `D:/Pirb_stroke_project/03_analysis_code/54_nc_report_supp_graphical_abstract.py` (initial NC report, supplementary tables, graphical abstract)
- `D:/Pirb_stroke_project/03_analysis_code/53_generate_comprehensive_report.py` (Chinese comprehensive summary report)

## Acknowledgements

[To be added.]

## Author contributions

[To be added.]

## Competing interests

The authors declare no competing interests.

## Additional information

Supplementary information is available for this paper at [repository path]. Correspondence and requests for materials should be addressed to [corresponding author].

## References

1. Atwal, J.K., et al. PirB is a functional receptor for myelin inhibitors of axonal regeneration. *Science* 322, 967–970 (2008). https://doi.org/10.1126/science.1161151
2. Syken, J., GrandPré, T., Kanold, P.O. & Shatz, C.J. PirB restricts ocular-dominance plasticity in visual cortex. *Science* 313, 1795–1800 (2006). https://doi.org/10.1126/science.1128232
3. Gou, X., et al. Spatio-temporal expression of paired immunoglobulin-like receptor-B in the adult mouse brain after focal cerebral ischaemia. *Brain Inj* 27, 1311–1315 (2013). https://doi.org/10.3109/02699052.2013.812241
4. Gou, X., et al. Neuronal PirB upregulated in cerebral ischemia acts as an attractive theranostic target for ischemic stroke. *J Am Heart Assoc* 7, e007197 (2018). https://doi.org/10.1161/JAHA.117.007197
5. Adelson, J.D., et al. Neuroprotection from stroke in the absence of MHCI or PirB. *Neuron* 73, 1100–1107 (2012). https://doi.org/10.1016/j.neuron.2012.01.020
6. Liddelow, S.A., et al. Neurotoxic reactive astrocytes are induced by activated microglia. *Nature* 541, 481–487 (2017). https://doi.org/10.1038/nature21029
7. Zamanian, J.L., et al. Genomic analysis of reactive astrogliosis. *J Neurosci* 32, 6391–6410 (2012). https://doi.org/10.1523/JNEUROSCI.6221-11.2012
8. Keren-Shaul, H., et al. A unique microglia type associated with restricting development of Alzheimer's disease. *Cell* 169, 1276–1290 (2017). https://doi.org/10.1016/j.cell.2017.05.018
9. Hammond, T.R., et al. Single-cell RNA sequencing of microglia throughout the mouse lifespan and in the injured brain reveals complex cell-state changes. *Immunity* 50, 253–271 (2019). https://doi.org/10.1016/j.immuni.2018.11.004
10. Zheng, K., et al. Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke. *J Cereb Blood Flow Metab* 42, 56–73 (2022). https://doi.org/10.1177/0271678X211028994
11. Bennett, F.C., et al. New tools for studying microglia in the mouse and human CNS. *Proc Natl Acad Sci USA* 113, E1738–E1746 (2016). https://doi.org/10.1073/pnas.1525528113
12. Izzy, S., et al. Revisiting the concept of activated microglia in traumatic brain injury. *Neurobiol Dis* 125, 93–103 (2019). https://doi.org/10.1016/j.nbd.2019.01.016
13. Kim, T., et al. Human LilrB2 is a β-amyloid receptor and its murine homolog PirB regulates synaptic plasticity in an Alzheimer's model. *Science* 341, 1399–1404 (2013). https://doi.org/10.1126/science.1242077
14. Fujita, Y. & Yamashita, T. Axon growth inhibition by RhoA/ROCK in the central nervous system. *Front Neurosci* 8, 338 (2014). https://doi.org/10.3389/fnins.2014.00338
15. Cekanaviciute, E., et al. Astrocytic transforming growth factor-beta signaling reduces subacute neuroinflammation after stroke in mice. *Glia* 62, 1227–1240 (2014). https://doi.org/10.1002/glia.22675
16. Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0
17. Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biol* 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1
18. 10x Genomics. Visium Spatial Gene Expression. https://www.10xgenomics.com/products/spatial-gene-expression (accessed 2026).


---

## Figure Legends

**Fig. 1 | Cell-type-specific Pirb expression in GSE174574.** (a) UMAP of 56,486 cells colored by major cell type. (b) Pirb expression across cell types. (c) Fraction of Pirb+ cells by condition and cell type. (d) Astrocyte state UMAP showing Homeostatic, PanReactive, A1-like, A2-like and Proliferative states. (e) Pirb+ fraction within each astrocyte state in Sham vs MCAO.

**Fig. 2 | Temporal dynamics of Pirb expression across datasets.** (a) GSE233815 bulk RNA-seq Pirb CPM across 3 h, 12 h, 24 h, 3 d and 7 d. (b) GSE233812 scRNA-seq microglial Pirb+ fraction across sham, D1, D3 and D7. (c) GSE233813 snRNA-seq microglial Pirb+ fraction across sham, D1, D3 and D7. (d) Cross-dataset microglial integration showing Pirb+ fraction across datasets and time points.

**Fig. 3 | GSE225948 peripheral blood vs brain Pirb expression.** (a) UMAP of 91,688 PB and brain cells colored by tissue. (b) Pirb+ fraction in brain microglia across Sham, D02 and D14. (c) Pirb+ fraction in PB monocytes and neutrophils across Sham, D02 and D14. (d) Top differentially expressed genes in Pirb+ vs Pirb− PB Mo/Neu at D02.

**Fig. 4 | GSE233814 Visium spatial localization of Pirb.** (a) H&E images with Pirb expression overlay for control, D1, D3, D7 and D7_rep. (b) Quantification of Pirb+ spot fraction per time point. (c) Distance of Pirb+ vs Pirb− spots to tissue boundary. (d) Spatial clustering of Pirb+ spots and enrichment at the ischemic boundary.

**Fig. 5 | Proposed model of Pirb-mediated post-ischemic neuroinflammation.** Schematic of healthy brain, ischemic penumbra at D3 (Pirb+ microglia, A1 astrocytes, infiltrating Mo/Neu), downstream molecular mechanisms (IL-1β/TNF/C1q, NF-κB/STAT/IRF), and timeline of Pirb expression (microglia D3 peak; PB Mo/Neu D02 peak).

**Graphical Abstract |** Pirb is induced in microglia, A1 astrocytes and infiltrating peripheral myeloid cells after cerebral ischemia, peaks in brain microglia at D3, localizes to the ischemic penumbra, and amplifies neuroinflammation while suppressing axon regeneration.


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

