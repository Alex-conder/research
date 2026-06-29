"""
将所有相关交付文件打包成 zip，放置于桌面。
"""
import os
import zipfile
from datetime import datetime

DESKTOP = "C:/Users/ASUS/Desktop"
ZIP_NAME = "Pirb_stroke_project_最终交付包.zip"
ZIP_PATH = os.path.join(DESKTOP, ZIP_NAME)

# 文件列表：(源绝对路径, zip 内目标路径)
files_to_pack = [
    # 最终报告
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_多数据集验证最终报告_原格式.docx", "01_最终报告_原格式/脑缺血后Pirb阳性细胞单细胞图谱_多数据集验证最终报告_原格式.docx"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_多数据集验证最终报告_原格式.md", "01_最终报告_原格式/脑缺血后Pirb阳性细胞单细胞图谱_多数据集验证最终报告_原格式.md"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告.docx", "02_NC格式报告/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告.docx"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告.md", "02_NC格式报告/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告.md"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告_中文版.docx", "02_NC格式报告/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告_中文版.docx"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告_中文版.md", "02_NC格式报告/脑缺血后Pirb阳性细胞单细胞图谱_NatureCommunications_完整报告_中文版.md"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_多数据集完整汇总报告.docx", "03_中文汇总报告/脑缺血后Pirb阳性细胞单细胞图谱_多数据集完整汇总报告.docx"),
    ("D:/Pirb_stroke_project/04_reports/脑缺血后Pirb阳性细胞单细胞图谱_多数据集完整汇总报告.md", "03_中文汇总报告/脑缺血后Pirb阳性细胞单细胞图谱_多数据集完整汇总报告.md"),

    # 深度机制分析专题报告
    ("D:/Pirb_stroke_project/04_reports/深度机制分析_5优先级结果汇总.docx", "04_深度机制分析/深度机制分析_5优先级结果汇总.docx"),
    ("D:/Pirb_stroke_project/04_reports/深度机制分析_5优先级结果汇总.md", "04_深度机制分析/深度机制分析_5优先级结果汇总.md"),

    # Supplementary Tables
    ("D:/Pirb_stroke_project/04_reports/supplementary/Supplementary_Tables_Pirb_Stroke.xlsx", "05_Supplementary_Tables/Supplementary_Tables_Pirb_Stroke.xlsx"),

    # Graphical Abstract
    ("D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.png", "06_Graphical_Abstract/graphical_abstract_draft.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.svg", "06_Graphical_Abstract/graphical_abstract_draft.svg"),
    ("D:/Pirb_stroke_project/04_reports/figures/graphical_abstract_draft.pdf", "06_Graphical_Abstract/graphical_abstract_draft.pdf"),

    # 投稿模板
    ("D:/Pirb_stroke_project/04_reports/Cover_Letter_Template.md", "07_投稿模板/Cover_Letter_Template.md"),
    ("D:/Pirb_stroke_project/04_reports/Response_to_Reviewers_Template.md", "07_投稿模板/Response_to_Reviewers_Template.md"),

    # 提升点对比与交付物清单
    ("D:/Pirb_stroke_project/04_reports/提升点对比表_原报告vs最终报告.pdf", "08_提升点对比/提升点对比表_原报告vs最终报告.pdf"),
    ("D:/Pirb_stroke_project/04_reports/项目交付物清单_绝对路径.md", "09_交付物清单/项目交付物清单_绝对路径.md"),

    # 关键差异基因与通路 CSV
    ("D:/Pirb_stroke_project/04_reports/figures/GSE174574/doublet_qc/DE_PirbPos_vs_Neg_Astrocyte_no_doublet.csv", "10_DE与通路/GSE174574_Astrocyte_PirbPos_vs_Neg_no_doublet.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb/DE_PirbPos_vs_Neg_Microglia.csv", "10_DE与通路/GSE174574_Microglia_PirbPos_vs_Neg.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb/DE_PirbPos_vs_Neg_Immune.csv", "10_DE与通路/GSE174574_Immune_PirbPos_vs_Neg.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/microglia_cross_time/pirb_pos_vs_neg_microglia_de.csv", "10_DE与通路/CrossTime_Microglia_PirbPos_vs_Neg.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE225948/DE_PirbPos_vs_Neg_PB_Mo_Neu_D02.csv", "10_DE与通路/GSE225948_PB_MoNeu_D02_PirbPos_vs_Neg.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb/manual_pathway_enrichment.csv", "10_DE与通路/GSE174574_Pathway_Enrichment.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/pirb_spatial_pattern_summary.csv", "10_DE与通路/GSE233814_Spatial_Summary.csv"),

    # 深度机制分析 CSV
    ("D:/Pirb_stroke_project/04_reports/figures/regulatory_brake/GSE174574_MCAO24h_pirb_pos_neg_comparison.csv", "10_DE与通路/RegulatoryBrake_GSE174574_MCAO24h.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/regulatory_brake/GSE233812_D3_pirb_pos_neg_comparison.csv", "10_DE与通路/RegulatoryBrake_GSE233812_D3.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/regulatory_brake/GSE233812_D7_pirb_pos_neg_comparison.csv", "10_DE与通路/RegulatoryBrake_GSE233812_D7.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/regulatory_brake/GSE233812_D7_vs_D3_inhibitory_stress_TFs.csv", "10_DE与通路/RegulatoryBrake_D7_vs_D3_inhibitory_TFs.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/gradient_de/D3_spatial_zones.csv", "10_DE与通路/GSE233814_D3_spatial_zones.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/gradient_de/penumbra_PirbPos_vs_Neg_DE.csv", "10_DE与通路/GSE233814_Penumbra_PirbPos_vs_Neg_DE.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/gradient_de/penumbra_PirbPos_vs_Neg_module_genes.csv", "10_DE与通路/GSE233814_Penumbra_module_genes.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/peripheral_central_lr/peripheral_PB_Mo_D02_to_brain_Mg_D3_lr_scores.csv", "10_DE与通路/Peripheral_Central_LR_scores.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/wgcna_microglia/gene_module_assignment.csv", "10_DE与通路/WGCNA_gene_module_assignment.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/wgcna_microglia/module_correlation_with_Pirb.csv", "10_DE与通路/WGCNA_module_correlation_with_Pirb.csv"),
    ("D:/Pirb_stroke_project/04_reports/figures/wgcna_microglia/gene_correlation_with_Pirb.csv", "10_DE与通路/WGCNA_gene_correlation_with_Pirb.csv"),

    # 关键图片
    ("D:/Pirb_stroke_project/04_reports/figures/GSE174574/marker_dotplot.png", "11_关键图片/GSE174574_marker_dotplot.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE225948/pirb_Mo_time.png", "11_关键图片/GSE225948_PB_Mo_Pirb_dynamics.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/microglia_cross_time/pirb_fraction_timeline.png", "11_关键图片/CrossTime_Microglia_Pirb_fraction_timeline.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/pirb_spatial_combined_panel.png", "11_关键图片/GSE233814_Pirb_spatial_combined_panel.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/regulatory_brake/boxplot_GSE233812_D3.png", "11_关键图片/RegulatoryBrake_GSE233812_D3.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/GSE233814/gradient_de/zone_module_scores_D3.png", "11_关键图片/GSE233814_zone_module_scores_D3.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/peripheral_central_lr/top_lr_pairs.png", "11_关键图片/Peripheral_Central_top_LR_pairs.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/peripheral_central_lr/lr_heatmap.png", "11_关键图片/Peripheral_Central_LR_heatmap.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/wgcna_microglia/module_corr_with_Pirb.png", "11_关键图片/WGCNA_module_corr_with_Pirb.png"),
    ("D:/Pirb_stroke_project/04_reports/figures/wgcna_microglia/top_gene_correlation_heatmap.png", "11_关键图片/WGCNA_top_gene_correlation_heatmap.png"),

    # 分析代码脚本
    ("D:/Pirb_stroke_project/03_analysis_code/53_generate_comprehensive_report.py", "12_分析代码/53_generate_comprehensive_report.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/54_nc_report_supp_graphical_abstract.py", "12_分析代码/54_nc_report_supp_graphical_abstract.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/55_gse225948_pb_de.py", "12_分析代码/55_gse225948_pb_de.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/56_final_assembly_all.py", "12_分析代码/56_final_assembly_all.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/57_final_report_original_format.py", "12_分析代码/57_final_report_original_format.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/58_generate_improvements_pdf.py", "12_分析代码/58_generate_improvements_pdf.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/59_package_deliverables.py", "12_分析代码/59_package_deliverables.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/60_nc_report_chinese.py", "12_分析代码/60_nc_report_chinese.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/61_regulatory_brake_analysis.py", "12_分析代码/61_regulatory_brake_analysis.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/62_spatial_gradient_de.py", "12_分析代码/62_spatial_gradient_de.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/63_peripheral_central_lr.py", "12_分析代码/63_peripheral_central_lr.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/64_rna_velocity.py", "12_分析代码/64_rna_velocity.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/65_wgcna_microglia.py", "12_分析代码/65_wgcna_microglia.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/66_deep_mechanism_summary.py", "12_分析代码/66_deep_mechanism_summary.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/67_append_deep_mechanism.py", "12_分析代码/67_append_deep_mechanism.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/68_scvelo_analysis.py", "12_分析代码/68_scvelo_analysis.py"),
    ("D:/Pirb_stroke_project/03_analysis_code/69_regenerate_docx.py", "12_分析代码/69_regenerate_docx.py"),

    # 原始参考报告
    ("C:/Users/ASUS/Desktop/脑缺血后Pirb阳性星形胶质细胞单细胞图谱_阶段性报告.docx", "13_原始参考报告/脑缺血后Pirb阳性星形胶质细胞单细胞图谱_阶段性报告.docx"),

    # RNA velocity 流程文档
    ("D:/Pirb_stroke_project/01_raw_data/GSE233812_velocity/scripts/run_velocity_pipeline.sh", "15_RNA_velocity流程/run_velocity_pipeline.sh"),
    ("D:/Pirb_stroke_project/01_raw_data/GSE233812_velocity/scripts/setup_wsl.md", "15_RNA_velocity流程/setup_wsl.md"),
    ("D:/Pirb_stroke_project/01_raw_data/GSE233812_velocity/scripts/download_D3_fastq.bat", "15_RNA_velocity流程/download_D3_fastq.bat"),
    ("D:/Pirb_stroke_project/01_raw_data/GSE233812_velocity/scripts/download_D3_fastq.sh", "15_RNA_velocity流程/download_D3_fastq.sh"),
    ("D:/Pirb_stroke_project/01_raw_data/GSE233812_velocity/scripts/download_with_retry.sh", "15_RNA_velocity流程/download_with_retry.sh"),
    ("D:/Pirb_stroke_project/04_reports/figures/rna_velocity/rna_velocity_status.txt", "15_RNA_velocity流程/rna_velocity_status.txt"),
]

# 检查所有文件存在
missing = []
for src, dst in files_to_pack:
    if not os.path.exists(src):
        missing.append(src)

if missing:
    print("[WARN] 以下文件缺失，将跳过：")
    for m in missing:
        print(f"  - {m}")

# 打包
with zipfile.ZipFile(ZIP_PATH, 'w', zipfile.ZIP_DEFLATED) as zf:
    for src, dst in files_to_pack:
        if os.path.exists(src):
            zf.write(src, dst)
            print(f"[ADD] {dst}")
        else:
            print(f"[SKIP] {dst}")

    # 递归添加 05_software_template 目录
    sw_dir = "D:/Pirb_stroke_project/05_software_template"
    if os.path.exists(sw_dir):
        for root, dirs, files in os.walk(sw_dir):
            for file in files:
                fpath = os.path.join(root, file)
                arcname = os.path.join("14_软件模板", os.path.relpath(fpath, sw_dir))
                zf.write(fpath, arcname)
                print(f"[ADD] {arcname}")

# 添加一个 README 说明
readme_content = f"""Pirb Stroke Project 最终交付包
生成日期：{datetime.now().strftime('%Y-%m-%d')}

目录结构：
- 01_最终报告_原格式：按原 Desktop 报告格式生成的中文最终报告（已含深度机制章节）
- 02_NC格式报告：Nature Communications 格式完整报告（英文 + 中文版）
- 03_中文汇总报告：多数据集完整汇总报告
- 04_深度机制分析：5 项优先级任务的结果汇总专题报告
- 05_Supplementary_Tables：Supplementary Tables 1-8（Excel）
- 06_Graphical_Abstract：Graphical Abstract 草图（PNG/SVG/PDF）
- 07_投稿模板：Cover Letter 和 Response to Reviewers 模板
- 08_提升点对比：原报告 vs 最终报告提升点 PDF 表格
- 09_交付物清单：所有交付物绝对路径清单
- 10_DE与通路：关键差异基因、通路富集、调控刹车、空间梯度、LR 对接、WGCNA CSV
- 11_关键图片：核心结果可视化图片
- 12_分析代码：生成上述结果的所有 Python 脚本
- 13_原始参考报告：Desktop 上的原始阶段性报告
- 14_软件模板：可复用的报告自动生成模板
- 15_RNA_velocity流程：scVelo 原始 FASTQ 下载、velocyto 流程脚本与 WSL 配置指南

项目根目录：D:/Pirb_stroke_project
"""
with zipfile.ZipFile(ZIP_PATH, 'a', zipfile.ZIP_DEFLATED) as zf:
    zf.writestr('README.txt', readme_content)

print(f"\n[SAVE] {ZIP_PATH}")
print(f"[SIZE] {os.path.getsize(ZIP_PATH) / 1024 / 1024:.2f} MB")
print("[DONE]")
