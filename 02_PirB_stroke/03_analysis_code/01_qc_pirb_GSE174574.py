"""
GSE174574 单细胞数据质控与 Pirb 表达初探
- 读取 6 个 10x 样本（文件平铺在 GSE174574/ 目录下，legacy genes.tsv.gz 命名）
- 基础质控
- 计算 Pirb 阳性率与表达水平
- 输出 CSV / 图片到 04_reports/figures/GSE174574
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 路径
DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE174574"
OUT_DIR  = "D:/Pirb_stroke_project/04_reports/figures/GSE174574"
os.makedirs(OUT_DIR, exist_ok=True)

# 从文件名识别样本前缀
barcode_files = sorted(glob.glob(os.path.join(DATA_DIR, "*_barcodes.tsv.gz")))
prefixes = sorted(set(os.path.basename(f).replace("_barcodes.tsv.gz", "") for f in barcode_files))
print(f"[INFO] Found {len(prefixes)} samples under {DATA_DIR}")

adatas = []
for prefix in prefixes:
    bc_file = os.path.join(DATA_DIR, f"{prefix}_barcodes.tsv.gz")
    gene_file = os.path.join(DATA_DIR, f"{prefix}_genes.tsv.gz")
    mtx_file = os.path.join(DATA_DIR, f"{prefix}_matrix.mtx.gz")
    
    # 手动读取 10x 三文件
    barcodes = pd.read_csv(bc_file, header=None, sep="\t")[0].values
    genes = pd.read_csv(gene_file, header=None, sep="\t")
    if genes.shape[1] == 2:
        gene_ids = genes[0].values
        gene_symbols = genes[1].values
    else:
        gene_symbols = genes[0].values
        gene_ids = gene_symbols
    
    X = mmread(mtx_file).T.tocsr()
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_symbols))
    ad.var_names_make_unique()
    ad.var["gene_ids"] = gene_ids
    ad.obs["sample"] = prefix
    
    if "sham" in prefix.lower():
        ad.obs["condition"] = "Sham"
    elif "mcao" in prefix.lower():
        ad.obs["condition"] = "MCAO"
    else:
        ad.obs["condition"] = "Unknown"
    adatas.append(ad)
    print(f"  {prefix}: {ad.n_obs} cells × {ad.n_vars} genes")

# 合并
adata = sc.concat(adatas, label="batch", keys=[a.obs["sample"][0] for a in adatas], index_unique="-")
adata.obs_names_make_unique()
print(f"[INFO] Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# 质控
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

# 过滤
n_before = adata.n_obs
adata = adata[adata.obs.n_genes_by_counts >= 500]
adata = adata[adata.obs.total_counts >= 1000]
adata = adata[adata.obs.pct_counts_mt < 20]
n_after = adata.n_obs
print(f"[INFO] QC: {n_before} -> {n_after} cells")

# 保存 QC 表
qc = adata.obs.groupby("sample").agg(
    condition=("condition", "first"),
    n_cells=("sample", "size"),
    median_n_genes=("n_genes_by_counts", "median"),
    median_total_counts=("total_counts", "median"),
    median_pct_mt=("pct_counts_mt", "median")
).reset_index()
qc.to_csv(os.path.join(OUT_DIR, "qc_summary.csv"), index=False)
print(qc)

# 标准化 + log1p
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Pirb 表达
pirb = "Pirb"
if pirb not in adata.var_names:
    print(f"[WARN] {pirb} not found in gene list")
else:
    adata.obs["Pirb_counts"] = adata[:, pirb].X.toarray().flatten()
    adata.obs["Pirb_detected"] = adata.obs["Pirb_counts"] > 0

    pirb_stats = adata.obs.groupby(["sample", "condition"]).agg(
        n_cells=("sample", "size"),
        pirb_positive=("Pirb_detected", "sum"),
        pirb_positive_pct=("Pirb_detected", "mean"),
        pirb_mean_counts=("Pirb_counts", "mean"),
        pirb_median_counts=("Pirb_counts", "median")
    ).reset_index()
    pirb_stats["pirb_positive_pct"] *= 100
    pirb_stats.to_csv(os.path.join(OUT_DIR, "pirb_expression_summary.csv"), index=False)
    print("\nPirb expression summary:")
    print(pirb_stats)

    # 小提琴图
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.violin(adata, keys=["Pirb_counts"], groupby="condition", ax=axes[0], show=False)
    axes[0].set_title("Pirb expression (log1p norm) by condition")
    sc.pl.violin(adata, keys=["Pirb_counts"], groupby="sample", ax=axes[1], show=False)
    axes[1].set_title("Pirb expression by sample")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_violin.png"), dpi=150)
    plt.close()

# 保存完整 log1p 数据（供下游使用）
adata.write_h5ad(os.path.join(OUT_DIR, "GSE174574_qc.h5ad"))

# 高变基因、降维、聚类（用于初步查看）
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack")
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5)

# 把 leiden 标签和降维坐标写回完整 adata
adata.obs["leiden"] = adata_hvg.obs["leiden"].values
for key in ["X_pca", "X_umap"]:
    if key in adata_hvg.obsm:
        adata.obsm[key] = adata_hvg.obsm[key]
for key in ["pca", "neighbors"]:
    if key in adata_hvg.uns:
        adata.uns[key] = adata_hvg.uns[key]
adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata.csv"))
adata.write_h5ad(os.path.join(OUT_DIR, "GSE174574_qc.h5ad"))

# UMAP 图
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(adata_hvg, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata_hvg, color="condition", ax=axes[1], show=False, title="Condition")
if pirb in adata_hvg.var_names:
    sc.pl.umap(adata_hvg, color="Pirb_counts", ax=axes[2], show=False, title="Pirb expression")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "umap_overview.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
