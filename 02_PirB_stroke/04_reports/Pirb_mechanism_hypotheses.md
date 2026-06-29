# Pirb 在脑缺血中的细胞类型特异性作用机制 — 假设与验证框架

## 一、分析数据基础

| 数据集 | 用途 | 状态 |
|--------|------|------|
| **GSE174574** | 主分析数据集：3 Sham + 3 MCAO（24 h），58,287 个细胞 | ✅ 完整 |
| **GSE171169** | 验证集：CD45high 免疫细胞（5d / 14d），10,295 个细胞 | ✅ 完整 |
| **GSE227651** | 时间序列：1d / 3d / 7d / sham | ⏳ 后台下载中 |
| **GSE225948** | 免疫/内皮/血液互作 | ⏳ 后台下载中 |

---

## 二、核心发现总结

### 2.1 Pirb 表达的细胞类型格局
- **免疫细胞（Ptprc+）Pirb 阳性率最高**：MCAO 31.4%，Sham 24.2%
- **小胶质细胞次之**：MCAO 6.2%，Sham 4.6%
- **星形胶质细胞**：MCAO 3.1%，Sham 0.7%（升高约 4.7 倍）
- **少突胶质细胞 / 内皮细胞 / 周细胞**：缺血后均升高 2–3 倍，但基数低

### 2.2 Pirb+ 星形胶质细胞的分子特征
对 83 个 Pirb+ vs 3405 个 Pirb− 星形胶质细胞做 Wilcoxon 差异分析，发现 **164 个显著上调基因**，其功能特征为：

| 功能类别 | 代表性基因 | 通路富集 p 值 |
|----------|-----------|--------------|
| **溶酶体/组织蛋白酶** | Ctsb, Ctsd, Ctsz, Ctss, Ctsc | 1.8 × 10⁻¹⁷ |
| **补体激活** | C1qa, C1qb, C3, C5ar1 | 7.0 × 10⁻¹¹ |
| **小胶质细胞诱导 A1** | C1qa, C1qb, Tnf, Il1a | 1.3 × 10⁻⁶ |
| **MHC I/抗原呈递** | H2-D1, H2-K1, Cd74 | 5.4 × 10⁻⁴ |
| **A1 反应性星形胶质细胞** | C3, Gfap, Vim, Serping1 | 3.5 × 10⁻² |
| **IRF 信号** | Irf1, Irf5, Irf7, Irf8 | 3.0 × 10⁻² |
| **NF-κB 信号** | Nfkb1, Nfkb2, Rela, Relb | 4.0 × 10⁻² |
| **炎症小体** | Nlrp3, Pycard, Casp1 | 3.0 × 10⁻² |

**下调基因**则包括星形胶质细胞稳态标志：Gja1、Slc1a3、Bcan、Mfge8、Ndrg2。

### 2.3 上游转录因子
在 Pirb+ 星形胶质细胞中显著上调的 TF：
- **Cebpb**（log₂FC = 1.75）
- **Rel**（log₂FC = 1.65）
- **Stat6**（log₂FC = 1.42）
- **Irf8**（log₂FC = 1.37）
- **Irf5**（log₂FC = 1.31）

### 2.4 配体-受体互作预测
- **少突胶质细胞**是 Pirb 配体（MAG / MOG / OMgp / Nogo-A / MHC I）的最强来源
- **MHC I → PirB** 对免疫细胞、小胶质细胞的信号最强（与这些细胞 Pirb 高表达一致）
- **指向星形胶质细胞的 Pirb 信号**主要来自少突胶质细胞、小胶质细胞和免疫细胞

### 2.5 伪时间动态
以 Sham 星形胶质细胞为根的扩散伪时间显示：
- Pirb 表达在伪时间 ~0.1–0.2 开始上升
- 在伪时间 ~0.35–0.6 达到峰值（Pirb+ 细胞比例 ~6%）
- 随后下降，提示 Pirb 高表达是缺血响应的中间状态

---

## 三、整合机制假设

### 假设 1：小胶质细胞通过 IL-1α / TNF / C1q 诱导 A1 样 Pirb+ 星形胶质细胞

**证据链**：
1. Pirb+ 星形胶质细胞显著富集补体 C1q（C1qa, C1qb）和 A1 反应性标志（C3, Gfap, Vim）
2. "Microglia-induced A1" 基因集显著富集（p = 1.3 × 10⁻⁶）
3. 缺血后小胶质细胞自身 Pirb 表达也升高，可能参与正反馈

**预测实验**：
- 体外用 IL-1α + TNF + C1q 处理原代星形胶质细胞，检测 Pirb 表达是否上调
- 使用小胶质细胞清除（CSF1R 抑制剂）或过继转移实验验证

---

### 假设 2：Pirb 在星形胶质细胞中放大 NF-κB / IRF / STAT 驱动的炎症程序

**证据链**：
1. Pirb+ 星形胶质细胞中 Nfkb1/2、Rel、Irf5/8、Stat6 显著上调
2. 下游炎症基因（Ccl4, Ccl12, Spp1, Lyz2, C5ar1）上调
3. 溶酶体/组织蛋白酶程序显著激活（p = 1.8 × 10⁻¹⁷）

**预测实验**：
- Pirb 阻断抗体或 Pirb KO 星形胶质细胞在 OGD 后检测 NF-κB p65 核转位
- 检测 Il1b、Tnf、Ccl2 等炎症因子分泌

---

### 假设 3：少突胶质细胞来源的 MAG / OMgp / MOG 通过 Pirb 抑制轴突再生并影响胶质瘢痕

**证据链**：
1. 少突胶质细胞是 MAG / MOG / OMgp 的主要来源
2. 少突胶质细胞 → 免疫细胞 / 小胶质细胞 / 星形胶质细胞的 Pirb 信号预测最强
3. Pirb 经典功能是介导髓鞘抑制分子对神经元的生长锥塌陷

**预测实验**：
- 在缺血半暗带检测 MAG/OMgp 与 Pirb 的共定位
- Pirb 阻断后检测轴突生长相关蛋白（GAP-43）和 RhoA-ROCK 活性

---

### 假设 4：Pirb 在免疫细胞和小胶质细胞中参与免疫抑制 / 突触修剪

**证据链**：
1. 免疫细胞 Pirb 阳性率最高（~30%），缺血后进一步升高
2. Pirb+ 免疫细胞上调 Ifitm、S100a、Lyz2、Alox5ap 等先天免疫标志
3. 小胶质细胞 Pirb+ 亚群上调 Pirb 和 Lyz2

**预测实验**：
- 分离 Pirb+ vs Pirb− 小胶质细胞，检测吞噬功能（pHrodo、突触素摄取）
- 检测 MHC II、CD86、CD206 等极化标志

---

### 假设 5：内皮细胞 Pirb 上调与血脑屏障破坏相关

**证据链**：
1. 内皮细胞 Pirb 阳性率缺血后从 0.58% 升至 1.77%（3 倍）
2. Pirb+ 内皮细胞下调 BBB 标志 **Abcb1a**（P-gp）、**Cxcl12**、**Rgcc**
3. 上调 Ctsb、Spp1、Hspa5、Ccl4 等应激/炎症基因

**预测实验**：
- 检测 Pirb+ 内皮细胞的紧密连接蛋白（Claudin-5、Occludin、ZO-1）表达
- Pirb 阻断后检测体外 BBB 模型跨内皮电阻（TEER）

---

## 四、与文献的衔接

| 文献发现 | 本分析支持 |
|----------|-----------|
| Atwal et al., Science 2008：PirB 是 Nogo/MAG/OMgp 功能受体 | 少突胶质细胞来源配体 → Pirb 信号预测最强 |
| Liddelow et al., Nature 2017：小胶质细胞 IL-1α/TNF/C1q 诱导 A1 星形胶质细胞 | "Microglia-induced A1" 基因集在 Pirb+ 星形胶质细胞中显著富集 |
| Wang et al., Mol Med Rep 2012：缺氧缺血后神经元 PirB 上调 | 尚需单独分析神经元亚群 |
| Liu et al., J Am Heart Assoc 2018：sPirB 阻断改善卒中预后 | 支持 Pirb 作为治疗靶点 |

---

## 五、后续验证优先级

| 优先级 | 实验/分析 | 预期结果 |
|--------|----------|----------|
| **高** | 小胶质细胞清除 + 星形胶质细胞 Pirb 表达 | 若 Pirb+ 星形胶质细胞减少，支持假设 1 |
| **高** | Pirb 阻断后检测炎症因子 / BBB 标志 | 若炎症减轻、TEER 升高，支持假设 2/5 |
| **中** | 免疫荧光：Pirb 与 C1q / C3 / H2-D1 共定位 | 验证 A1 与 MHC I 共表达 |
| **中** | 单细胞水平 Pirb  reporter / 流式分选 | 确认 Pirb+ 细胞身份，排除 doublet |
| **低** | GSE227651 / GSE225948 时间动态 | 后台下载完成后整合验证 |

---

## 六、输出文件索引

```
04_reports/figures/GSE174574/
├── de_pirb/
│   ├── DE_PirbPos_vs_Neg_Astrocyte.csv
│   ├── DE_PirbPos_vs_Neg_Microglia.csv
│   ├── DE_PirbPos_vs_Neg_Immune.csv
│   ├── volcano_Astrocyte.png
│   └── manual_pathway_enrichment.csv
├── de_pirb_all_ct/
│   ├── DE_PirbPos_vs_Neg_*.csv
│   ├── volcano_*.png
│   └── de_summary.csv
├── lr_manual/
│   ├── lr_scores_all.csv
│   ├── lr_to_astrocyte.csv
│   └── interaction_heatmap_mcao.png
├── lr_pirb_specific/
│   ├── pirb_lr_scores.csv
│   ├── pirb_lr_summary.csv
│   ├── pirb_pathway_summary.csv
│   └── pirb_signal_heatmap_mcao.png
├── tf_analysis/
│   ├── tf_expression_by_celltype_pirb.csv
│   ├── tf_log2fc.csv
│   ├── tf_up_in_pirb_pos.csv
│   └── tf_heatmap_log2fc.png
└── pseudotime/
    ├── astrocyte_pseudotime_metadata.csv
    ├── pirb_vs_pseudotime.png
    └── pirb_trend_pseudotime.png
```

---

## 七、结论

基于 GSE174574 和 GSE171169 的整合分析，Pirb 在脑缺血后呈现**广泛的细胞类型上调**，尤其在免疫细胞、小胶质细胞和星形胶质细胞中。Pirb+ 星形胶质细胞具有典型的 **A1 神经毒性/炎症表型**，其上游可能由小胶质细胞来源的 IL-1α/TNF/C1q 诱导，并通过 NF-κB/IRF/STAT 转录程序驱动。同时，少突胶质细胞来源的 MAG/OMgp/MOG 可能通过 Pirb 持续抑制轴突再生。这些假设为后续的功能验证和靶向干预提供了明确方向。
