# ==============================================================================
# 脑缺血后 Pirb 阳性细胞单细胞图谱 — 横向扩展与纵向深挖分析流程
# 数据集：GSE174574（主） + GSE227651 / GSE225948 / GSE233815（扩展）
# ==============================================================================

# ---------------------------- 0. 环境准备 ------------------------------------
required_pkgs <- c("Seurat", "SeuratData", "SeuratDisk", "harmony", "sctransform",
                   "GEOquery", "dplyr", "ggplot2", "patchwork", "cowplot",
                   "SingleR", "celldex", "scran", "scater", "BiocParallel",
                   "clusterProfiler", "org.Mm.eg.db", "ReactomePA", "DOSE",
                   "CellChat", "nichenetr", "monocle3", "slingshot",
                   "SCENIC", "AUCell", "RcisTarget", "ComplexHeatmap", "circlize")

# install_if_missing <- function(pkgs) {
#   BiocManager::install(pkgs[!pkgs %in% installed.packages()[, "Package"]])
# }
# install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(SingleR)
  library(celldex)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(CellChat)
})

set.seed(20260615)

# 工作路径
work_dir <- "D:/Pirb_stroke_project"
data_dir <- file.path(work_dir, "01_raw_data")
out_dir <- file.path(work_dir, "04_reports")
fig_dir <- file.path(work_dir, "05_supplementary", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- 1. 下载数据（可选） ----------------------------
# 若网络通畅，可用 GEOquery 直接下载；否则先用 Python 脚本断点续传。
# library(GEOquery)
# gse <- getGEO("GSE174574", GSEMatrix = TRUE, getGPL = FALSE)
# pData(gse[[1]])  # 查看样本信息

# ---------------------------- 2. 读取 10x 数据 -------------------------------
read_10x_sample <- function(sample_dir, sample_id, condition) {
  counts <- Read10X(data.dir = sample_dir)
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
  obj$sample <- sample_id
  obj$condition <- condition
  obj$dataset <- "GSE174574"
  return(obj)
}

# 示例：假设已解压 GSE174574_RAW.tar 到 01_raw_data/GSE174574/
# 实际样本目录需根据解压后的文件名调整
gse174574_samples <- list(
  list(dir = file.path(data_dir, "GSE174574", "GSM5319987"), id = "GSM5319987", cond = "Sham"),
  list(dir = file.path(data_dir, "GSE174574", "GSM5319988"), id = "GSM5319988", cond = "Sham"),
  list(dir = file.path(data_dir, "GSE174574", "GSM5319989"), id = "GSM5319989", cond = "Sham"),
  list(dir = file.path(data_dir, "GSE174574", "GSM5319990"), id = "GSM5319990", cond = "MCAO"),
  list(dir = file.path(data_dir, "GSE174574", "GSM5319991"), id = "GSM5319991", cond = "MCAO"),
  list(dir = file.path(data_dir, "GSE174574", "GSM5319992"), id = "GSM5319992", cond = "MCAO")
)

# 仅当目录存在时读取
obj_list <- list()
for (s in gse174574_samples) {
  if (dir.exists(s$dir)) {
    obj_list[[s$id]] <- read_10x_sample(s$dir, s$id, s$cond)
  }
}

if (length(obj_list) == 0) {
  stop("未找到 GSE174574 样本目录。请先解压 RAW.tar 并检查路径。")
}

# ---------------------------- 3. 质控 ----------------------------------------
qc_and_filter <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                  nCount_RNA > 500 & percent.mt < 15)
  return(obj)
}
obj_list <- lapply(obj_list, qc_and_filter)

# ---------------------------- 4. 标准化与整合 --------------------------------
# 使用 SCTransform + Harmony 整合
obj_list <- lapply(obj_list, SCTransform, vars.to.regress = "percent.mt", verbose = FALSE)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                   anchor.features = features)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# PCA + Harmony 去批次
combined <- RunPCA(combined, features = features, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = c("sample", "condition"), assay.use = "integrated")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.8)

# ---------------------------- 5. 细胞注释 ------------------------------------
# 5.1 SingleR 自动注释
ref <- MouseRNAseqData()
pred <- SingleR(test = GetAssayData(combined, slot = "data"),
                ref = ref,
                labels = ref$label.main)
combined$singleR_label <- pred$labels[match(colnames(combined), rownames(pred))]

# 5.2 人工 marker 复核（主要脑内细胞类型）
markers <- list(
  Astrocyte = c("Aqp4", "Aldh1l1", "Gfap", "Slc1a2", "Gja1"),
  Microglia = c("P2ry12", "Tmem119", "Hexb", "C1qa", "Cx3cr1"),
  Macrophage_Monocyte = c("Cd14", "Cd86", "Lyz2", "Itgam", "Fcgr1"),
  Neuron = c("Rtn1", "Rbfox3", "Syt1", "Map2"),
  Oligodendrocyte = c("Mog", "Mbp", "Olig1", "Plp1", "Mag"),
  OPC = c("Pdgfra", "Cspg4", "Olig2"),
  Endothelial = c("Cldn5", "Flt1", "Ly6c1", "Esam", "Pecam1"),
  Pericyte_VSMC = c("Rgs5", "Vtn", "Acta2", "Des", "Myl9"),
  Ependymal = c("Ttr", "Clic6", "Sostdc1"),
  Lymphocyte = c("Cd3e", "Cd19", "Cd79a", "Ms4a1")
)

# 计算每个细胞类型的模块评分
for (ct in names(markers)) {
  genes <- markers[[ct]]
  genes <- intersect(genes, rownames(combined))
  if (length(genes) >= 3) {
    combined <- AddModuleScore(combined, features = list(genes), name = paste0(ct, "_score"))
  }
}

# 基于最高评分与 SingleR 结果综合注释（简化版）
# 实际项目中建议对每群人工复核 dotplot
combined$celltype <- combined$singleR_label
# 将 SingleR 的 Astrocytes, Microglia, Oligodendrocytes 等映射为统一命名
combined$celltype <- gsub("Astrocytes", "Astrocyte", combined$celltype)
combined$celltype <- gsub("Microglia", "Microglia", combined$celltype)
combined$celltype <- gsub("Oligodendrocytes", "Oligodendrocyte", combined$celltype)
combined$celltype <- gsub("Neurons", "Neuron", combined$celltype)
combined$celltype <- gsub("Endothelial cells", "Endothelial", combined$celltype)
combined$celltype <- gsub("Pericytes", "Pericyte", combined$celltype)

# ---------------------------- 6. Pirb 表达分析（横向扩展） -------------------
# 6.1 全细胞类型 Pirb 阳性率
pirb_expr <- FetchData(combined, vars = c("Pirb", "celltype", "condition"))
pirb_expr$Pirb_positive <- pirb_expr$Pirb > 0
pirb_summary <- pirb_expr %>%
  group_by(celltype, condition) %>%
  summarise(
    n_cells = n(),
    Pirb_positive_n = sum(Pirb_positive),
    Pirb_positive_pct = mean(Pirb_positive) * 100,
    mean_Pirb_expr = mean(Pirb),
    .groups = "drop"
  )
write.csv(pirb_summary, file.path(out_dir, "Pirb_positive_rate_by_celltype.csv"), row.names = FALSE)

# 6.2 可视化
p1 <- DimPlot(combined, reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
p2 <- FeaturePlot(combined, features = "Pirb", order = TRUE, pt.size = 0.1)
p3 <- FeaturePlot(combined, features = "Pirb", split.by = "condition", order = TRUE, pt.size = 0.1)
ggsave(file.path(fig_dir, "UMAP_celltype.pdf"), p1, width = 8, height = 6)
ggsave(file.path(fig_dir, "FeaturePlot_Pirb.pdf"), p2, width = 6, height = 5)
ggsave(file.path(fig_dir, "FeaturePlot_Pirb_split.pdf"), p3, width = 12, height = 5)

# ---------------------------- 7. 亚群分析 ------------------------------------
# 以星形胶质细胞为例，复现并扩展原报告
celltypes_to_analyze <- c("Astrocyte", "Microglia", "Endothelial", "Oligodendrocyte", "Neuron")

for (ct in celltypes_to_analyze) {
  ct_obj <- subset(combined, subset = celltype == ct)
  if (ncol(ct_obj) < 100) next
  ct_obj <- RunPCA(ct_obj, verbose = FALSE)
  ct_obj <- RunHarmony(ct_obj, group.by.vars = "sample")
  ct_obj <- RunUMAP(ct_obj, reduction = "harmony", dims = 1:20)
  ct_obj <- FindNeighbors(ct_obj, reduction = "harmony", dims = 1:20)
  ct_obj <- FindClusters(ct_obj, resolution = 0.6)

  # Pirb 阳性率
  ct_expr <- FetchData(ct_obj, vars = c("Pirb", "condition"))
  ct_expr$Pirb_positive <- ct_expr$Pirb > 0
  ct_summary <- ct_expr %>%
    group_by(condition) %>%
    summarise(n = n(), pct = mean(Pirb_positive) * 100, mean_expr = mean(Pirb), .groups = "drop")
  write.csv(ct_summary, file.path(out_dir, paste0("Pirb_positive_rate_", ct, ".csv")), row.names = FALSE)

  # Pirb+ vs Pirb- 差异分析
  if (sum(ct_expr$Pirb_positive) >= 10 && sum(!ct_expr$Pirb_positive) >= 10) {
    ct_obj$Pirb_binary <- ifelse(ct_expr$Pirb > 0, "Pirb_pos", "Pirb_neg")
    Idents(ct_obj) <- "Pirb_binary"
    degs <- FindAllMarkers(ct_obj, ident.1 = "Pirb_pos", ident.2 = "Pirb_neg",
                           only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(degs, file.path(out_dir, paste0("DEG_Pirbpos_vs_neg_", ct, ".csv")), row.names = FALSE)
  }

  # MCAO vs Sham 差异分析（pseudobulk 更稳健）
  # 此处用 Wilcoxon 作为示例
  Idents(ct_obj) <- "condition"
  degs_cond <- FindMarkers(ct_obj, ident.1 = "MCAO", ident.2 = "Sham",
                           min.pct = 0.25, logfc.threshold = 0.25)
  degs_cond$gene <- rownames(degs_cond)
  write.csv(degs_cond, file.path(out_dir, paste0("DEG_MCAO_vs_Sham_", ct, ".csv")), row.names = FALSE)
}

# ---------------------------- 8. 星形胶质细胞反应状态评分（复现+扩展） -------
astro <- subset(combined, subset = celltype == "Astrocyte")
astro <- RunPCA(astro, verbose = FALSE)
astro <- RunHarmony(astro, group.by.vars = "sample")
astro <- RunUMAP(astro, reduction = "harmony", dims = 1:20)
astro <- FindNeighbors(astro, reduction = "harmony", dims = 1:20)
astro <- FindClusters(astro, resolution = 0.5)

# A1 / A2 / PanReactive / Homeostatic 评分
astro_state_genes <- list(
  AstroA1_like = c("C3", "Gfap", "Vim", "H2.T23", "Ggta1", "Srgn"),
  AstroA2_like = c("S100a10", "S100a6", "Clcf1", "Tgm1", "Ptgs2"),
  AstroPanReactive = c("Gfap", "Vim", "Cst3", "Timp1", "Serpina3n"),
  AstroHomeostatic = c("Aqp4", "Aldh1l1", "Slc1a2", "Gja1", "Clu", "Ntsr2")
)
for (st in names(astro_state_genes)) {
  g <- intersect(astro_state_genes[[st]], rownames(astro))
  if (length(g) >= 3) {
    astro <- AddModuleScore(astro, features = list(g), name = st)
  }
}
astro$Pirb_positive <- FetchData(astro, vars = "Pirb")$Pirb > 0

# 按条件 + Pirb 阴阳性输出评分
astro_scores <- FetchData(astro, vars = c("condition", "Pirb_positive",
                                           "AstroA1_like1", "AstroA2_like1",
                                           "AstroPanReactive1", "AstroHomeostatic1"))
score_summary <- astro_scores %>%
  group_by(condition, Pirb_positive) %>%
  summarise(across(ends_with("1"), mean, na.rm = TRUE), .groups = "drop")
write.csv(score_summary, file.path(out_dir, "Astrocyte_state_scores_by_Pirb.csv"), row.names = FALSE)

# ---------------------------- 9. 细胞通讯（CellChat） -------------------------
# 构建 CellChat 对象，分析 MCAO 组细胞通讯
mcao <- subset(combined, subset = condition == "MCAO")
data.input <- GetAssayData(mcao, assay = "SCT", slot = "data")
meta <- data.frame(cell = colnames(data.input),
                   celltype = mcao$celltype,
                   condition = mcao$condition,
                   row.names = colnames(data.input))

cc <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cc@DB <- CellChatDB.mouse
cellchat_major <- subsetCellChat(cc, idents.use = c("Astrocyte", "Microglia", "Endothelial",
                                                    "Oligodendrocyte", "Neuron", "Macrophage_Monocyte"))
cellchat_major <- identifyOverExpressedGenes(cellchat_major)
cellchat_major <- identifyOverExpressedInteractions(cellchat_major)
cellchat_major <- computeCommunProb(cellchat_major, raw.use = TRUE)
cellchat_major <- filterCommunication(cellchat_major, min.cells = 10)
cellchat_major <- computeCommunProbPathway(cellchat_major)
cellchat_major <- aggregateNet(cellchat_major)

cellchat_major <- netAnalysis_computeCentrality(cellchat_major, slot.name = "netP")
netVisual_chord_gene(cellchat_major, sources.use = c("Microglia", "Macrophage_Monocyte"),
                     targets.use = c("Astrocyte"), lab.cex = 0.5,
                     title.name = "Microglia/Macrophage -> Astrocyte")
ggsave(file.path(fig_dir, "CellChat_microglia_to_astrocyte.pdf"), width = 8, height = 8)

# ---------------------------- 10. 通路富集（以 Astrocyte DEG 为例） ----------
astro_degs <- read.csv(file.path(out_dir, "DEG_MCAO_vs_Sham_Astrocyte.csv"))
sig <- astro_degs %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25)
genes <- sig$gene
entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

if (nrow(entrez) > 0) {
  ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
                  readable = TRUE, pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  write.csv(as.data.frame(ego), file.path(out_dir, "GO_BP_Astrocyte_MCAO_up.csv"), row.names = FALSE)

  kk <- enrichKEGG(gene = entrez$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
  write.csv(as.data.frame(kk), file.path(out_dir, "KEGG_Astrocyte_MCAO_up.csv"), row.names = FALSE)
}

# ---------------------------- 11. 保存结果 -----------------------------------
saveRDS(combined, file.path(data_dir, "GSE174574_integrated.rds"))
saveRDS(astro, file.path(data_dir, "GSE174574_astrocyte_subset.rds"))

message("分析完成。结果保存在：", out_dir)
message("图表保存在：", fig_dir)
