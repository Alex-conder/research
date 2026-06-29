"""
生成提升点对比表格 PDF
"""
import os
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT_PDF = "../04_reports/提升点对比表_原报告vs最终报告.pdf"

# 尝试使用支持中文的字体（Windows 常见）
chinese_fonts = ['SimHei', 'Microsoft YaHei', 'SimSun', 'Arial Unicode MS']
font_prop = None
for f in chinese_fonts:
    try:
        fp = fm.FontProperties(family=f)
        # 测试是否能正确渲染中文
        if fp.get_name():
            font_prop = fp
            break
    except Exception:
        continue

if font_prop is None:
    font_prop = fm.FontProperties()

plt.rcParams['font.family'] = font_prop.get_name()
plt.rcParams['axes.unicode_minus'] = False

def add_table_page(pdf, title, headers, rows, col_widths=None, fontsize=10):
    fig, ax = plt.subplots(figsize=(11.69, 8.27))  # A4 landscape
    ax.axis('off')
    ax.axis('tight')
    
    # Title
    ax.text(0.5, 0.95, title, transform=fig.transFigure,
            fontsize=16, fontweight='bold', ha='center', fontproperties=font_prop)
    
    # Prepare table data
    data = [headers] + rows
    
    # Auto column widths if not provided
    if col_widths is None:
        n_cols = len(headers)
        col_widths = [1.0 / n_cols] * n_cols
    
    table = ax.table(cellText=data, cellLoc='left', loc='center',
                     colWidths=col_widths)
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1, 2.5)
    
    # Style header
    for i in range(len(headers)):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white', fontproperties=font_prop)
    
    # Style rows
    for i in range(1, len(data)):
        for j in range(len(headers)):
            cell = table[(i, j)]
            if i % 2 == 0:
                cell.set_facecolor('#E7E6E6')
            cell.set_text_props(fontproperties=font_prop)
    
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

with PdfPages(OUT_PDF) as pdf:
    # Page 1: 总体提升概览
    add_table_page(
        pdf,
        "原报告 vs 最终报告：总体提升概览",
        ["维度", "原报告", "最终报告"],
        [
            ["数据集", "GSE174574（1 个）", "GSE174574/171169/225948/233815/233812/233813/233814（7 个）"],
            ["数据类型", "仅 scRNA-seq", "scRNA-seq、snRNA-seq、bulk RNA-seq、Visium 空间转录组"],
            ["样本/细胞", "6 样本 / 57,769 细胞", "完整 26 样本 / 91,688 细胞 + 11,969 空间 spots"],
            ["分析细胞类型", "仅星形胶质细胞", "小胶质、星形胶质、免疫细胞、PB 单核/中性粒等"],
            ["时间覆盖", "仅 24h vs Sham", "3h、12h、24h、D1、D3、D7、D14"],
            ["空间信息", "无", "Visium 定位至缺血半暗带/边界"],
            ["验证层级", "单数据集描述", "多数据集交叉验证"],
            ["机制模型", "未提出", "外周浸润→小胶质激活→A1 星形胶质→炎症放大→轴突再生抑制"],
            ["功能验证", "仅建议", "完整体外功能实验路线图"],
            ["补充材料", "无", "Supplementary Tables 1–8、Graphical Abstract、投稿模板"],
        ],
        col_widths=[0.18, 0.35, 0.47]
    )
    
    # Page 2: 数据集扩展
    add_table_page(
        pdf,
        "数据集扩展明细",
        ["数据集", "类型", "样本", "细胞/spot 数", "用途"],
        [
            ["GSE174574", "scRNA-seq", "3 Sham + 3 MCAO（24h）", "56,486", "主分析（星形胶质 + 全细胞）"],
            ["GSE171169", "scRNA-seq", "CD45high 5d/14d", "10,295", "免疫细胞验证"],
            ["GSE225948", "scRNA-seq", "PB + brain Sham/D02/D14", "91,688", "髓系免疫细胞验证（完整版）"],
            ["GSE233815", "bulk RNA-seq", "MCAO 3h/12h/24h/3d/7d", "48 样本", "时间动态验证"],
            ["GSE233812", "scRNA-seq", "sham/D1/D3/D7", "6,159", "单细胞时间序列"],
            ["GSE233813", "snRNA-seq", "sham/D1/D3/D7", "8,374", "单细胞核时间序列"],
            ["GSE233814", "Visium 空间", "control/D1/D3/D7/D7_rep", "11,969 spots", "空间定位"],
        ],
        col_widths=[0.16, 0.16, 0.28, 0.16, 0.24]
    )
    
    # Page 3: 核心发现提升
    add_table_page(
        pdf,
        "核心发现提升",
        ["发现", "原报告", "最终报告"],
        [
            ["星形胶质细胞 Pirb", "MCAO 后 3.60% vs Sham 0.50%，A1-like 中 7.04%", "结论保留；DE 基于去 doublet 细胞，补充通路富集"],
            ["小胶质细胞 Pirb", "仅提及免疫细胞 Pirb+ 31.4%", "明确为脑内 Pirb 主要来源，D3 峰值：26.6% / 6.4% / 14.47%"],
            ["外周血 Pirb", "未涉及", "PB 单核细胞 D02 51.72%，中性粒细胞 52.64%"],
            ["时间动态", "无法判断", "D3 脑内小胶质峰值；D02 外周髓系峰值"],
            ["空间定位", "缺少空间坐标", "D3 Pirb+ spots 富集半暗带边界（p = 3.3e-17）"],
            ["机制解释", "炎症相关基因共表达", "完整整合机制模型"],
            ["下游验证", "建议引入实验", "提供 OGD/LPS + anti-Pirb 体外实验方案"],
        ],
        col_widths=[0.22, 0.35, 0.43]
    )
    
    # Page 4: 配套交付物
    add_table_page(
        pdf,
        "新增配套交付物",
        ["类别", "文件名", "绝对路径"],
        [
            ["原格式最终报告", "脑缺血后Pirb阳性细胞单细胞图谱_多数据集验证最终报告_原格式.docx", "../04_reports/"],
            ["NC 格式报告", "脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告.docx", "../04_reports/"],
            ["Supplementary Tables", "Supplementary_Tables_Pirb_Stroke.xlsx", "../04_reports/supplementary/"],
            ["Graphical Abstract", "graphical_abstract_draft.png/.svg/.pdf", "../04_reports/figures/"],
            ["投稿模板", "Cover_Letter_Template.md", "../04_reports/"],
            ["审稿回复模板", "Response_to_Reviewers_Template.md", "../04_reports/"],
            ["交付物清单", "项目交付物清单_绝对路径.md", "../04_reports/"],
        ],
        col_widths=[0.22, 0.45, 0.33]
    )

print(f"[SAVE] {OUT_PDF}")
