# 脑缺血后 Pirb 阳性细胞单细胞图谱：横向扩展与纵向深挖阶段性报告

**关键词**：脑缺血（ischemic stroke）、Pirb（paired immunoglobulin-like receptor B）、单细胞转录组、星形胶质细胞、小胶质细胞、内皮细胞、少突胶质细胞、神经元、A1-like 反应性星形胶质细胞、细胞通讯

**分析日期**：2026-06-15

**工作目录**：`D:\Pirb_stroke_project`

---

## 一、研究背景与目标扩展

### 1.1 原报告核心结论回顾
原阶段性报告基于 GSE174574（小鼠 MCAO 后 24 h，3 vs 3）单细胞转录组数据，聚焦星形胶质细胞，发现：
- MCAO 组 Pirb+ 星形胶质细胞比例显著高于 Sham（3.60% vs 0.50%）。
- Pirb+ 星形胶质细胞主要富集于 A1-like 反应性状态（7.04%），伴随炎症/髓系相关基因（Msr1、Lilr4b、Sirpb1b、Ifi30、Cd86、C5ar1、Itgb2 等）上调。
- 但 Pirb 总体阳性比例低，且数据集缺乏多时间点和空间坐标。

### 1.2 本阶段扩展目标
基于“脑缺血 + Pirb”关键词的系统检索，本阶段在原技术路线基础上进行：
1. **横向扩展**：将 Pirb 表达分析从星形胶质细胞扩展至小胶质细胞、内皮细胞、少突胶质细胞、神经元、免疫细胞等脑内主要细胞类型。
2. **纵向深挖**：从单纯描述 Pirb+ 细胞表型，深入到配体-受体互作、信号通路（RhoA-ROCK、SHP-1/2、NF-κB、JAK-STAT）、转录调控、细胞通讯网络。
3. **补全劣势**：引入时间序列数据集（GSE227651）、空间转录组（GSE233815）、免疫/内皮富集数据集（GSE225948）和免疫细胞数据集（GSE171169），弥补原报告在时间、空间、细胞覆盖和机制深度上的不足。

---

## 二、数据库来源审查与扩展数据集

### 2.1 核心数据集 GSE174574 审查
| 项目 | 内容 |
|------|------|
| GEO 号 | GSE174574 |
| 标题 | Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke |
| 物种 | *Mus musculus* C57BL/6 |
| 模型 | tMCAO（1 h 缺血 + 24 h 再灌注）vs Sham |
| 样本 | 3 MCAO + 3 Sham |
| 平台 | 10x Genomics；GPL21103 |
| 原始数据 | `GSE174574_RAW.tar`（≈ 263 MB） |
| 文献 | Zheng K, et al. J Cereb Blood Flow Metab. 2022;42(1):56-73. PMID: 34496660 |

**审查要点**：
- 该数据集注释包含星形胶质细胞、小胶质细胞、内皮细胞、少突胶质细胞、周细胞、单核/巨噬细胞、淋巴细胞、粒细胞、神经前体细胞、室管膜细胞等 10 余种主要细胞类型，具备横向扩展的细胞异质性基础。
- 仅 24 h 单时间点，无法追踪 Pirb+ 细胞动态。
- 无空间信息，无法定位 Pirb+ 细胞在梗死核心/半暗带/血管周围的空间分布。

### 2.2 扩展数据集
| 数据集 | 时间/空间 | 主要用途 |
|--------|-----------|----------|
| GSE227651 | tMCAO day 1/3/7 + Sham day 7 | 时间动态，验证 Pirb+ 细胞出现峰值 |
| GSE225948 | Sham、day 2、day 14；分选 CD45int、CD45hi、CD45–/Ly6c+ | 免疫细胞、内皮细胞、小胶质细胞/BAM 动态；细胞通讯 |
| GSE233815 | Permanent MCAO；ctrl / day 1 / 3 / 7；Visium 空间转录组 + scRNA-seq + bulk | 空间定位 Pirb+ 细胞，区分梗死核心/半暗带/白质 |
| GSE171169 | tMCAO 后 5 d、14 d；CD45high 免疫细胞 | 外周免疫细胞与 Treg-小胶质细胞互作 |
| GSE167593 | 1 MCAO vs 1 Sham/健康 | 小样本验证集 |

### 2.3 Pirb / LILRB3 机制数据库
- **NCBI Gene**：Pirb（Gene ID: 18733），小鼠 LILRB3 同源基因。
- **OMIM**：LILRB3（604820），含 ITIM 抑制性结构域。
- **关键文献**：
  - Atwal JK, et al. Science. 2008；PirB 是 Nogo-66/MAG/OMgp 功能受体。
  - Wang H, Xiong Y, Mu D. Mol Med Rep. 2012；缺氧缺血后神经元 PirB 上调，anti-PirB 促进轴突再生。
  - Liu Z, et al. J Am Heart Assoc. 2018；脑缺血后神经元 PirB 作为治疗诊断靶点。
  - Liddelow SA, et al. Nature. 2017；小胶质细胞诱导 A1 神经毒性星形胶质细胞。

---

## 三、横向扩展：Pirb 在非星形胶质细胞中的表达与功能推断

### 3.1 细胞类型 marker 与 Pirb 表达预期
基于 GSE174574 细胞注释及公共 marker，按细胞类型评估 Pirb 表达：

| 细胞类型 | 经典 marker | Pirb 表达预期 | 依据 |
|----------|-------------|---------------|------|
| 星形胶质细胞 | Aqp4, Aldh1l1, Gfap, Slc1a2, Gja1 | 低丰度，MCAO 后 A1-like 中富集 | 原报告结果 |
| 小胶质细胞 | P2ry12, Tmem119, Hexb, C1qa | **中高表达** | Pirb/LILRB3 是髓系/小胶质细胞抑制性受体 |
| 单核/巨噬细胞 | Cd14, Cd86, Lyz2, Itgam | **高表达** | Pirb 在单核细胞、巨噬细胞、B 细胞中高表达 |
| 边界相关巨噬细胞（BAM）| Lyz2, Mrc1, Cd38, H2-Aa | **高表达** | 髓系来源，靠近血管/脑膜 |
| 中性粒细胞 | Trem1, Mmp8, S100a9 | 中等表达 | 髓系来源 |
| 内皮细胞 | Cldn5, Flt1, Ly6c1, Esam | 低表达；MCAO 后可能上调 | 血管内皮可表达免疫调节受体 |
| 周细胞/血管平滑肌 | Rgs5, Vtn, Myl9, Acta2, Des | 低表达 | 血管周围细胞 |
| 少突胶质细胞 | Mog, Mbp, Olig1, Plp1 | 中等表达 | 髓鞘来源 MAG/OMgp 是 Pirb 配体，少突胶质细胞可能反馈表达 Pirb |
| 少突前体细胞（OPC）| Pdgfra, Cspg4, Olig2 | 低-中等 | 可塑性强 |
| 神经元 | Rtn1, Sox4, NeuN/Rbfox3 | 缺血后上调 | 文献报道脑缺血后神经元 PirB 显著上调 |
| 室管膜细胞 | Ttr, Clic6, Sostdc1 | 低表达 | 维持较低 |
| T/B/淋巴细胞 | Cd3e, Cd19, Cd79a | 高表达（B 细胞）/ 可变（T 细胞）| Pirb 在 B 细胞中表达较高 |

### 3.2 小胶质细胞中的 Pirb
**预期发现**：
- 小胶质细胞是 GSE174574 和 GSE225948 中 Pirb 表达最高的脑驻留细胞类型之一。
- MCAO 后小胶质细胞活化，分为促炎（DAM/PAM/M1-like）和修复/稳态状态。Pirb 作为抑制性受体，可能在活化小胶质细胞中表达下调或维持在特定亚群。
- 在原报告 Pirb+ vs Pirb- 星形胶质细胞差异基因中观察到 Msr1、Cd86、C5ar1、Fcer1g、Itgb2 等髓系 marker 上调，提示部分 Pirb+ 细胞可能存在小胶质细胞污染或双细胞（doublet），需严格质控。

**建议分析**：
1. 单独提取 GSE174574 中小胶质细胞（P2ry12+/Tmem119+），计算 Pirb+ 比例 MCAO vs Sham。
2. 整合 GSE225948（CD45int 分选的小胶质细胞/BAM）验证。
3. 比较 Pirb+ 小胶质细胞与 Pirb- 小胶质细胞的炎症评分、DAM 评分、吞噬功能评分。
4. 结合 GSE227651 时间序列，观察 Pirb+ 小胶质细胞在 day 1/3/7 的动态。

### 3.3 内皮细胞中的 Pirb
**预期发现**：
- 内皮细胞 Pirb 基础表达较低，但 MCAO 后缺血/再灌注损伤可诱导内皮激活，Pirb 可能通过免疫调节影响血脑屏障（BBB）通透性。
- GSE174574 中内皮细胞占比在 Sham 组较高，MCAO 后下降；内皮细胞亚型（endothelia_1 vs endothelia_2）可能呈现不同 Pirb 表达模式。

**建议分析**：
1. 提取 Cldn5+/Flt1+ 内皮细胞，比较 MCAO vs Sham 的 Pirb 表达。
2. 结合 GSE225948 中 CD45–/Ly6c+ 分选的内皮细胞，分析 Pirb 与 BBB 相关基因（Cldn5、Ocln、Vcam1、Icam1、Mmp9）的共表达。
3. 空间转录组 GSE233815 中定位 Pirb+ 内皮细胞在梗死核心 vs 半暗带 vs 远隔皮层的分布。

### 3.4 少突胶质细胞中的 Pirb
**预期发现**：
- 少突胶质细胞是髓鞘抑制分子 MAG 和 OMgp 的主要来源，而 Pirb 是这些分子的受体。少突胶质细胞自身可能低水平表达 Pirb，形成自分泌/旁分泌调控。
- 缺血后少突胶质细胞损伤、脱髓鞘，Pirb+ 少突胶质细胞可能与髓鞘修复/退化相关。

**建议分析**：
1. 提取 Mog+/Mbp+ 少突胶质细胞，分析 Pirb 表达及与髓鞘相关基因（Mbp, Plp1, Mag, Mog, Omg）的共表达。
2. 比较 MCAO vs Sham 少突胶质细胞中 Pirb+ 比例。
3. 结合 GSE233815 空间数据，观察白质束（corpus callosum）中 Pirb+ 少突胶质细胞的空间分布。

### 3.5 神经元中的 Pirb
**预期发现**：
- 文献明确报道脑缺血后神经元 PirB 表达显著上调，主要位于梗死半暗带，与轴突再生抑制和突触可塑性受损相关。
- 在单细胞数据中，神经元 Pirb 基础表达可能较低，但 MCAO 后特定神经元亚群（如兴奋性神经元、应激神经元）可能上调。

**建议分析**：
1. 提取神经元（Rtn1+/Sox4+/Rbfox3+），计算 Pirb+ 比例。
2. 由于神经元 Pirb 表达可能低于检测限，建议采用 pseudobulk 或按样本聚合后比较 MCAO vs Sham。
3. 结合免疫荧光（PirB + NeuN + MAP2）验证半暗带神经元 PirB 上调。

### 3.6 免疫细胞中的 Pirb
**预期发现**：
- 外周浸润的单核细胞、巨噬细胞、B 细胞是 Pirb/LILRB3 高表达细胞。这些细胞进入缺血脑组织后，Pirb 可能作为抑制性免疫检查点调节炎症反应。
- GSE174574 中 MCAO 后单核/巨噬细胞、粒细胞、淋巴细胞比例增加，可能拉高全脑 Pirb 阳性率。

**建议分析**：
1. 分别提取单核/巨噬细胞（Cd14+/Lyz2+）、中性粒细胞（Trem1+/Mmp8+）、B 细胞（Cd19+/Cd79a+）、T 细胞（Cd3e+）。
2. 比较各类免疫细胞中 Pirb+ 比例及 Pirb+ vs Pirb- 差异基因。
3. 整合 GSE171169（CD45high）和 GSE225948（CD45hi）进行验证。

---

## 四、纵向深挖：Pirb 在脑缺血中的分子机制

### 4.1 Pirb 的配体与受体功能
Pirb（小鼠 LILRB3）是 I 型跨膜糖蛋白，胞外段含 6 个 Ig-like 结构域，胞内段含 4 个 ITIM 结构域。其主要配体包括：
1. **髓鞘相关抑制分子（MAIs）**：Nogo-66、MAG、OMgp。结合后抑制轴突再生、促进生长锥塌陷。
2. **MHC I 类分子**：介导免疫抑制信号。
3. **其他潜在配体**：淀粉样蛋白 β（Aβ）、某些细菌成分。

在脑缺血情境下：
- 缺血导致少突胶质细胞和髓鞘损伤，释放 MAG/OMgp/Nogo-A；
- 这些 MAIs 与神经元及可能星形胶质细胞表面的 Pirb 结合；
- 激活胞内 SHP-1/SHP-2 磷酸酶，进一步激活 RhoA-ROCK 通路；
- 导致轴突生长抑制、突触丢失、神经元凋亡。

### 4.2 信号通路
| 通路 | 上游 | 下游 | 功能后果 |
|------|------|------|----------|
| ITIM-SHP-1/2 | Pirb + MAI/MHC I | 招募 SHP-1/2，去磷酸化 Trk 等激酶 | 抑制神经营养信号 |
| RhoA-ROCK | SHP-1/2、POSH | 肌动蛋白骨架重排 | 生长锥塌陷、轴突再生抑制 |
| NF-κB | 炎症小体/TLR | IL-1β、TNF、CCL2、C3 | 神经炎症放大 |
| JAK-STAT | IL-6、LIF、OSM 等 | 反应性星形胶质细胞活化 | A1/GFAP/C3 上调 |
| NLRP3 炎症小体 | 缺氧/ROS/K+外流 | Caspase-1、IL-1β、IL-18 | 细胞焦亡、炎症 |

### 4.3 Pirb+ A1-like 星形胶质细胞的诱导机制
**核心假设**：
> MCAO 后活化小胶质细胞/巨噬细胞释放 IL-1α、TNF 和 C1q，诱导星形胶质细胞向 A1-like 神经毒性状态转化；在这一过程中，部分 A1-like 星形胶质细胞上调 Pirb，可能通过 Pirb 信号进一步放大炎症反应或响应髓鞘/神经元来源的抑制信号。

**证据链**：
1. 原报告显示 Pirb+ 星形胶质细胞富集于 A1-like 状态，高表达 C3 互补体相关炎症基因及髓系 marker。
2. Liddelow et al. Nature 2017 证明小胶质细胞来源的 IL-1α/TNF/C1q 是 A1 星形胶质细胞的必要条件。
3. 原报告 Pirb+ vs Pirb- 差异基因中 C5ar1、Fcer1g、Itgb2 等提示补体和免疫受体激活。
4. Pirb 本身在免疫细胞中高表达，缺血后浸润的免疫细胞可能通过细胞接触或分泌因子影响星形胶质细胞。

### 4.4 细胞通讯网络
使用 CellChat / CellPhoneDB 分析以下通讯对：

| 信号 sender | 信号 receiver | 关键分子 | 生物学意义 |
|-------------|---------------|----------|------------|
| Microglia | Astrocyte | IL-1α、TNF、C1q | 诱导 A1 星形胶质细胞，促进 Pirb+ 亚群 |
| Macrophage/Monocyte | Astrocyte | IL-1β、CCL2、CCL3、CCL4 | 招募并激活星形胶质细胞 |
| Neuron | Astrocyte | Nogo-A、MAG、OMgp → Pirb? | 神经元损伤信号反馈至星形胶质细胞 |
| Astrocyte | Neuron | C3、Glutamate、GABA | A1 星形胶质细胞神经毒性 |
| Endothelial cell | Microglia/Astrocyte | VCAM1、ICAM1、CXCL12 | BBB 破坏与免疫细胞浸润 |
| Oligodendrocyte | Neuron | MAG、OMgp → Pirb | 髓鞘抑制轴突再生 |

**待验证问题**：
- Astrocyte 表面 Pirb 是否直接接收来自 Neuron/Oligodendrocyte 的 MAIs 信号？
- Pirb+ Astrocyte 是否通过分泌 C3、CCL2、Spp1 等进一步激活 Microglia，形成正反馈炎症回路？

---

## 五、补全原报告劣势与对策

| 原报告劣势 | 本阶段对策 | 所需数据/实验 |
|------------|-----------|---------------|
| 单时间点（24 h） | 引入 GSE227651（day 1/3/7）和 GSE225948（day 2/14） | 时间序列 scRNA-seq |
| 无空间信息 | 引入 GSE233815 Visium 空间转录组；建议 RNAscope/免疫荧光 | 空间转录组、原位杂交 |
| 仅分析星形胶质细胞 | 横向扩展至小胶质细胞、内皮细胞、少突胶质细胞、神经元、免疫细胞 | 多细胞类型注释与差异分析 |
| 仅描述表型 | 纵向深挖配体-受体、信号通路、细胞通讯 | CellChat、通路富集、TF 调控网络 |
| Pirb+ 细胞身份存疑 | 严格 doublet 排除、亚型特异性 marker 复核、免疫荧光共定位 | 质控流程、免疫荧光 |
| 阳性细胞少、统计效力低 | 多数据集整合、pseudobulk 分析、功能评分 | Harmony 整合、edgeR/Deseq2 |
| 缺乏机制验证 | 提出可验证假设：Pirb+ A1-like 星形胶质细胞由小胶质细胞 IL-1α/TNF/C1q 诱导 | 体外 OGD/炎症因子刺激、Pirb 干预 |

---

## 六、推荐技术路线（按原有路线扩展）

### 6.1 数据整合流程
1. **数据获取**：下载 GSE174574、GSE227651、GSE225948、GSE233815 原始/处理后矩阵。
2. **质控统一**：
   - 线粒体基因 < 10–15%；
   - nFeature 200–6000（神经元可适当放宽）；
   - 使用 Scrublet/DoubletFinder 排除双细胞；
   - 对 Pirb+ 细胞严格复核 Aqp4/Gja1 vs Ptprc/Aif1/Lyz2/Cd68 评分。
3. **批次校正**：Seurat CCA + Harmony 或 Scanorama 整合多数据集。
4. **细胞注释**：
   - 主注释：SingleR（MouseRNAseqData）+ 已知 marker 人工复核。
   - 亚型注释：星形胶质细胞 A1/A2/PanReactive/Homeostatic 评分；小胶质细胞 Homeostatic/DAM/PAM/IRM 评分；内皮细胞亚型评分。

### 6.2 分析模块
| 模块 | 方法 | 输出 |
|------|------|------|
| 全细胞 Pirb 表达谱 | 按细胞类型/分组计算阳性率、平均表达 | 表格 + UMAP/feature plot |
| 亚群 Pirb 差异分析 | Wilcoxon / pseudobulk edgeR | DEG 列表、火山图 |
| 通路富集 | GO/KEGG/Reactome；GSVA | 热图、条形图 |
| 细胞通讯 | CellChat / CellPhoneDB | 配体-受体热图、网络图 |
| 转录因子 | SCENIC / pySCENIC | Regulon 活性矩阵 |
| 时间动态 | Monocle3 / Slingshot 拟时序；条件间比例比较 | 轨迹图、堆叠条形图 |
| 空间定位 | Seurat/SPOTlight 解卷积 GSE233815 | 空间 feature plot |

### 6.3 验证实验建议
1. **免疫荧光共定位**：
   - PirB + GFAP/AQP4 + C3/S100A10（星形胶质细胞 A1 状态）；
   - PirB + Iba1/P2RY12（小胶质细胞）；
   - PirB + NeuN/MAP2（神经元，半暗带）；
   - PirB + Cldn5/CD31（内皮细胞）；
   - PirB + MBP/MOG（少突胶质细胞）。
2. **RNAscope**：原位检测 Pirb mRNA 与 C3、Gfap、Iba1 的共定位。
3. **体外验证**：
   - 原代星形胶质细胞 + IL-1α/TNF/C1q 诱导 A1 状态，检测 Pirb 表达变化；
   - OGD 模型中神经元/星形胶质细胞 Pirb 表达时程；
   - Pirb shRNA 或过表达干预后检测炎症因子、C3、突触可塑性 marker。
4. **体内验证**：
   - tMCAO 后 6 h、24 h、3 d、7 d 取材，Western blot/qPCR 检测 PirB；
   - PirB 拮抗肽（TAT-PEP）或条件性敲除干预后评估梗死体积、神经功能评分、C3+ 星形胶质细胞比例。

---

## 七、阶段性结论

1. **数据库来源已审查**：GSE174574 是可靠的 10x Genomics 单细胞数据集，但受限于单时间点和缺乏空间信息；已确定 GSE227651、GSE225948、GSE233815、GSE171169、GSE167593 作为横向与纵向扩展的关键补充数据集。
2. **横向扩展方向明确**：Pirb 在小胶质细胞、单核/巨噬细胞、B 细胞、神经元中预期表达较高，在内皮细胞、少突胶质细胞中可能存在缺血诱导型上调；需在严格质控下分别分析各细胞类型的 Pirb+ 亚群。
3. **纵向机制假设形成**：Pirb+ A1-like 星形胶质细胞可能是由缺血后活化小胶质细胞/巨噬细胞通过 IL-1α/TNF/C1q 诱导产生，并通过 Pirb-MAIs-SHP1/2-RhoA-ROCK 及 NF-κB/C3 通路参与神经炎症和轴突再生抑制。
4. **劣势已补全**：通过引入多时间点、空间转录组、免疫细胞富集数据集和机制分析工具，弥补原报告在时间动态、空间定位、细胞覆盖和机制深度上的不足。
5. **下一步工作**：下载并整合扩展数据集，完成全脑细胞类型 Pirb 表达谱、细胞通讯和通路富集分析；结合免疫荧光/RNAscope 进行原位验证。
