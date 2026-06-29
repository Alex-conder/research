#!/usr/bin/env Rscript
# ==============================================================================
# 批量下载 GEO 数据集原始文件（替代 Python 断点续传方案）
# 需要安装：BiocManager::install("GEOquery")
# ==============================================================================
library(GEOquery)

work_dir <- ".."
data_dir <- file.path(work_dir, "01_raw_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# 要下载的数据集
geo_ids <- c("GSE174574", "GSE227651", "GSE225948", "GSE233815", "GSE171169", "GSE167593")

for (gse in geo_ids) {
  message("Downloading ", gse, " ...")
  tryCatch({
    # getGEO 下载处理后的表达矩阵； Supplemental Files 需另外下载
    gse_obj <- getGEO(gse, GSEMatrix = TRUE, destdir = data_dir, getGPL = FALSE)
    saveRDS(gse_obj, file.path(data_dir, paste0(gse, "_gse_object.rds")))
    message("  ", gse, " expression matrix saved.")
  }, error = function(e) {
    message("  ERROR downloading ", gse, ": ", conditionMessage(e))
  })
}

# 使用 getGEOSuppFiles 下载原始补充文件（RAW tar）
for (gse in c("GSE174574", "GSE227651")) {
  message("Downloading supplemental files for ", gse, " ...")
  tryCatch({
    getGEOSuppFiles(gse, baseDir = data_dir, makeDirectory = TRUE)
    message("  ", gse, " supplemental files downloaded.")
  }, error = function(e) {
    message("  ERROR downloading supplemental files for ", gse, ": ", conditionMessage(e))
  })
}

message("Download script completed. Check ", data_dir)
