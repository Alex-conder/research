#!/usr/bin/env Rscript
# =====================================================
# MetaboAnalystR 本地通路富集分析（修正版）
# 用法: Rscript run_metaboanalyst_local.R [5ppm|10ppm]
# 输入: output/pathway_combined_<ppm>/pathway_input.csv
# 输出: output/metaboanalyst_results/
#
# 环境要求：
# 1. R 4.5.x + Rtools45（gcc 在 x86_64-w64-mingw32.static.posix/bin）
# 2. Bioconductor 3.22；若默认镜像慢，可在 R 中设置
#    options(BioC_mirror = "https://bioconductor.riken.jp")
# 3. MetaboAnalystR 通过 remotes::install_github("xia-lab/MetaboAnalystR") 安装
# 4. 首次运行会下载 compound_db.qs / syn_nms.qs 到工作目录；
#    如网络慢，可手动从 https://www.metaboanalyst.ca/resources/libs/ 下载。
# =====================================================

# 选择输入版本
args <- commandArgs(trailingOnly = TRUE)
ppm <- ifelse(length(args) > 0, args[1], "10ppm")
# 可传入第二个参数作为输入目录；否则按 ppm 推断
input_dir <- ifelse(length(args) > 1, args[2], paste0("./output/pathway_combined_", ppm))

# 0. 检查并安装 MetaboAnalystR（从 GitHub）
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
  message("正在从 GitHub 安装 MetaboAnalystR，可能需要数分钟...")
  remotes::install_github("xia-lab/MetaboAnalystR", build_vignettes = FALSE)
}
library(MetaboAnalystR)

# 1. 设置工作目录与输出目录
setwd("d:/vscode workplace/metabolomics_pipeline")
output_dir <- "./output/metaboanalyst_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. 读取通路输入文件
input_file <- file.path(input_dir, "pathway_input.csv")
if (!file.exists(input_file)) {
  stop("找不到输入文件: ", input_file)
}
data <- read.csv(input_file, stringsAsFactors = FALSE, encoding = "UTF-8")

# ORA 只需要化合物名称列表；去除 Unknown / 空值
cmpd_names <- data$Compound_Name
cmpd_names <- cmpd_names[!is.na(cmpd_names) & nzchar(trimws(cmpd_names)) & cmpd_names != "Unknown"]
cmpd_names <- unique(trimws(cmpd_names))

if (length(cmpd_names) < 3) {
  stop("有效化合物名称少于 3 个，无法进行 ORA。")
}
message("读取到 ", length(cmpd_names), " 个唯一化合物名称用于 ORA")

# 3. 初始化 MetaboAnalyst 对象（pathora = 通路 ORA）
mSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 72)

# 4. 设置化合物名称并做名称映射（KEGG/HMDB/PubChem/CHEBI/METLIN）
mSet <- Setup.MapData(mSet, cmpd_names)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)

# 查看映射结果（会被写入 name_map.csv）
map_table <- as.data.frame(mSet$dataSet$map.table, stringsAsFactors = FALSE)
message("成功映射 ", sum(map_table$Comment == "1"), " / ", nrow(map_table), " 个化合物")

# 5. 设置 KEGG 通路库：小鼠 "mmu"（Neuro-2a 来源）
mSet <- SetKEGG.PathLib(mSet, "mmu", "current")

# 6. 设置参考代谢组背景（可选）
# TRUE 表示以 HMDB 化合物为参考背景；FALSE 表示以 KEGG 通路库中所有化合物为背景
mSet <- SetMetabolomeFilter(mSet, FALSE)

# 7. 执行 ORA：rbc = relative-betweenness centrality；hyperg = 超几何检验
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")

# 8. 读取并保存结果
# MetaboAnalystR 默认会写入 pathway_results.csv；这里再复制到输出目录
ora_result <- mSet$analSet$ora.mat
if (is.null(ora_result) || nrow(ora_result) == 0) {
  warning("未得到 ORA 结果，请检查化合物名称映射是否为空。")
} else {
  out_csv <- file.path(output_dir, paste0("pathway_enrichment_metaboanalyst_", ppm, ".csv"))
  write.csv(ora_result, out_csv)
  message("[OK] ORA 结果已保存: ", out_csv)
}

# 9. 保存名称映射表
out_map <- file.path(output_dir, paste0("name_mapping_metaboanalyst_", ppm, ".csv"))
write.csv(map_table, out_map, row.names = FALSE)
message("[OK] 名称映射表已保存: ", out_map)

# 10. 绘制通路摘要图（本地未加载 current.kegglib，暂不启用）
# tryCatch({
#   PlotPathSummary(mSet, show.grid = TRUE, imgName = file.path(output_dir, "pathway_summary"), format = "png")
#   message("[OK] 通路摘要图已保存: ", file.path(output_dir, "pathway_summary.png"))
# }, error = function(e) {
#   message("通路摘要图绘制失败：", conditionMessage(e))
# })

message("\n[OK] MetaboAnalystR ORA 流程结束，结果目录: ", normalizePath(output_dir))

# =====================================================
# 附录：如果你想用 MSEA（基于连续 Log2FC，不截断显著性），
# 请将第 3-7 步替换为以下代码：
# =====================================================
# mSet2 <- InitDataObjects("conc", "pathqea", FALSE)
# mSet2 <- Read.TextData(mSet2, input_file, "rowu", "disc")
# mSet2 <- SanityCheckData(mSet2)
# mSet2 <- ReplaceMin(mSet2)
# mSet2 <- PreparePrenormData(mSet2)
# mSet2 <- Normalization(mSet2, "NULL", "NULL", "NULL", "NULL", ratio = FALSE, ratioNum = 20)
# mSet2 <- CrossReferencing(mSet2, "name")
# mSet2 <- CreateMappingResultTable(mSet2)
# mSet2 <- SetKEGGPathLib(mSet2, "mmu")
# mSet2 <- CalculateQeaScore(mSet2, "rbc", "gt")   # gt = global test
# write.csv(mSet2$analSet$qea.mat, file.path(output_dir, "pathway_msea_metaboanalyst.csv"))
