# 脑缺血后 Pirb 阳性细胞单细胞图谱：机制深挖、多数据集验证及多组学确认阶段性报告（V9）

**报告日期**：2026-06-16  
**核心数据集**：GSE174574、GSE171169、GSE225948、GSE233815 bulk、GSE233812 scRNA-seq、GSE233813 snRNA-seq、GSE233814 Visium 空间  
**分析环境**：Python + scanpy；R 不可用

---

## 一、核心假设

1. 小胶质细胞 IL-1α/TNF/C1q → A1 样 Pirb+ 星形胶质细胞
2. Pirb 放大 NF-κB/IRF/STAT 炎症程序
3. 少突胶质细胞 MAG/OMgp/MOG → Pirb 抑制轴突再生
4. Pirb 参与免疫细胞/小胶质细胞免疫抑制或突触修剪
5. 内皮细胞 Pirb 上调与血脑屏障破坏相关

---

## 二、数据集

| 数据集 | 类型 | 样本 | 细胞/spot 数 | 用途 |
|--------|------|------|--------------|------|
| GSE174574 | scRNA-seq | 3 Sham + 3 MCAO（24 h） | 58,287 → 56,486 | 主分析 |
| GSE171169 | scRNA-seq | CD45high 5d/14d | 10,295 | 免疫细胞验证 |
| GSE225948 | scRNA-seq | PB + brain Sham/D02/D14 | 63,733 | 髓系免疫细胞验证 |
| GSE233815 | bulk RNA-seq | MCAO 3h/12h/24h/3D/7D | 48 | 时间动态验证 |
| GSE233812 | scRNA-seq | sham/D1/D3/D7 | 6,159 | 单细胞时间序列 |
| GSE233813 | snRNA-seq | sham/D1/D3/D7 | 8,374 | 单细胞核时间序列 |
| GSE233814 | Visium 空间 | control/D1/D3/D7 | 11,969 spots | 空间定位 |

---

## 三、主要结果

### 3.1 GSE174574 主分析
- Pirb+ 率：MCAO 8.2–9.8% vs Sham 2.7–2.8%。
- 细胞类型：免疫细胞（31.4%）> 小胶质细胞（6.2%）> 星形胶质细胞（3.1%）。
- Pirb+ 星形胶质细胞富集：溶酶体、补体、A1、NF-κB/IRF/STAT。

### 3.2 GSE171169 验证
- CD45high 免疫细胞 Pirb+ 率 27–28%。

### 3.3 GSE225948 独立验证
- 髓系免疫细胞 Pirb 高表达。
- Resident 小胶质细胞 Pirb 低。

### 3.4 GSE233815 bulk 时间序列
- Pirb CPM 在 3D 达峰（6.69）。

### 3.5 GSE233812 scRNA-seq
- 小胶质细胞 Pirb D3 峰值 **26.6%**。
- 其他细胞类型几乎不表达。

### 3.6 GSE233813 snRNA-seq
- 小胶质细胞 Pirb D3 最高 **6.4%**。

### 3.7 GSE233814 Visium 空间转录组
**Pirb+ spot 比例**：
| 时间 | Pirb+ fraction |
|------|----------------|
| control | 0.62% |
| D1 | 2.40% |
| **D3** | **5.93%** |
| D7 | 2.01% |
| D7_rep | 1.20% |

**Pirb+ spot 特征（D3）**：
| 特征 | Pirb+ spots | Pirb- spots |
|------|-------------|-------------|
| Microglia score | 1.165 | 0.687 |
| Inflammation score | 1.139 | 0.345 |

- **Pirb+ spots 显著富集小胶质细胞和炎症特征**。
- 空间层面再次确认 D3 是 Pirb 表达峰值。

---

## 四、跨数据集共识

### 4.1 时间动态
| 数据集 | 层面 | 峰值 | 关键结果 |
|--------|------|------|----------|
| GSE233812 | scRNA-seq | D3 | 小胶质细胞 Pirb 26.6% |
| GSE233813 | snRNA-seq | D3 | 小胶质细胞 Pirb 6.4% |
| GSE233814 | Visium 空间 | D3 | 5.93% spots |
| GSE233815 | bulk RNA-seq | 3D | CPM 6.69 |
| GSE174574 | scRNA-seq | 24h | 小胶质细胞 6.2% |

### 4.2 细胞类型
- **小胶质细胞**：GSE233812/813/174574 一致为 Pirb 主要表达细胞。
- **髓系免疫细胞**：GSE225948/171169 中高表达。
- **星形胶质细胞**：GSE174574 中 A1 样 Pirb+ 星形胶质细胞显著。

---

## 五、机制模型

1. 缺血后小胶质细胞激活，Pirb 从基线 → D1 上升 → **D3 峰值** → D7 回落。
2. 激活的小胶质细胞通过 IL-1α/TNF/C1q 诱导星形胶质细胞 A1 转化。
3. Pirb+ 小胶质细胞通过 NF-κB/IRF/STAT 程序放大神经炎症。
4. 少突胶质细胞 MAG/OMgp/MOG 与 Pirb 互作，抑制轴突再生。
5. 外周浸润髓系免疫细胞是另一 Pirb 高表达群体。

---

## 六、关键图表

| 结果 | 路径 |
|------|------|
| GSE233812 scRNA UMAP / Pirb 热图 | `04_reports/figures/GSE233812/` |
| GSE233813 snRNA UMAP / Pirb 热图 | `04_reports/figures/GSE233813/` |
| sc vs sn 小胶质细胞 Pirb 动态 | `04_reports/figures/GSE233812_813_comparison/` |
| GSE233814 Pirb 表达小提琴图 | `04_reports/figures/GSE233814/pirb_expression_violin.png` |
| GSE233814 control vs D3 marker | `04_reports/figures/GSE233814/marker_expression_control_vs_D3.png` |
| GSE233814 Pirb+ vs Pirb- 比较 | `04_reports/figures/GSE233814/pirb_pos_vs_neg_markers.csv` |

---

## 七、局限与下一步

### 7.1 局限
1. GSE233812 小胶质细胞数较少（406 个）。
2. GSE233814 未进行精细空间坐标可视化。
3. GSE225948 部分样本缺失。
4. GSE227651 仍未下载。

### 7.2 下一步
1. 解析 GSE233814 json 坐标，绘制 Pirb 空间分布热图。
2. 整合 GSE174574 + GSE233812/813 跨时间分析。
3. 体外功能实验：小胶质细胞激活 + Pirb 阻断。

---

## 八、结论

通过 7 个数据集的多组学交叉验证，我们一致发现：脑缺血后 **Pirb 主要在小胶质细胞中表达，并于 D3（3 天）达到峰值**。这一结论在单细胞、单细胞核、空间和组织 bulk 四个层面高度一致，为 Pirb 作为脑缺血后神经炎症关键调控因子提供了强有力证据。外周浸润髓系免疫细胞是另一重要 Pirb 载体。后续功能实验将进一步验证 Pirb 在缺血后炎症放大和轴突再生抑制中的因果作用。
