"""
处理 GSE233812 scRNA-seq 数据（sham / D1 / D3 / D7）
genes.tsv.gz 只有 2 列，需要自定义读取
"""
import os, glob, gc, tempfile
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import io, sparse
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

DATA_DIR = "../01_raw_data/GSE233812"
OUT_DIR = "../04_reports/figures/GSE233812"
os.makedirs(OUT_DIR, exist_ok=True)

samples = {
    "GSM7437213_sc-sham": "sham",
    "GSM7437214_sc-D1": "D1",
    "GSM7437215_sc-D3": "D3",
    "GSM7437216_sc-D7": "D7",
}

def read_10x_custom(prefix):
    """读取只有 2 列 genes 文件的 10x 数据"""
    barcodes = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None, compression='gzip')[0].values
    genes = pd.read_csv(f"{prefix}_genes.tsv.gz", header=None, compression='gzip', sep='\t')
    genes.columns = ['gene_id', 'gene_name']
    genes['feature_type'] = 'Gene Expression'
    
    # 读取 matrix
    mat = io.mmread(f"{prefix}_matrix.mtx.gz")
    mat = sparse.csr_matrix(mat.T)  # 转置为 cell x gene
    
    var = pd.DataFrame(index=genes['gene_name'].values)
    var['gene_id'] = genes['gene_id'].values
    var['feature_type'] = genes['feature_type'].values
    
    ad = anndata.AnnData(X=mat, obs=pd.DataFrame(index=barcodes), var=var)
    # 处理重复 gene symbol
    ad.var_names_make_unique()
    return ad

adatas = []
for prefix, time_point in samples.items():
    full_prefix = os.path.join(DATA_DIR, prefix)
    print(f"[READ] {prefix} -> {time_point}")
    ad = read_10x_custom(full_prefix)
    ad.obs['sample'] = prefix
    ad.obs['time_point'] = time_point
    ad.obs['time_point_num'] = {'sham': 0, 'D1': 1, 'D3': 3, 'D7': 7}[time_point]
    print(f"  Loaded {ad.n_obs} droplets x {ad.n_vars} genes")
    adatas.append(ad)
    gc.collect()

print("[INFO] Concatenating...")
adata = anndata.concat(adatas, axis=0, join='outer', fill_value=0)
adata.X = adata.X.tocsr()
print(f"[INFO] Merged: {adata.n_obs} droplets x {adata.n_vars} genes")

# 计算 MT%
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 初始过滤
print(f"[INFO] Before cell filter: {adata.n_obs}")
sc.pp.filter_cells(adata, min_genes=200)
print(f"[INFO] After min_genes=200: {adata.n_obs}")

# QC 分布图
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(adata.obs['n_genes_by_counts'], bins=100)
axes[0].set_title('n_genes distribution')
axes[1].hist(adata.obs['total_counts'], bins=100)
axes[1].set_title('total_counts distribution')
axes[2].hist(adata.obs['pct_counts_mt'], bins=100)
axes[2].set_title('pct_mt distribution')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "qc_distributions.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

# 宽松 QC
n_genes_low = 500
n_genes_high = 8000
total_counts_high = 40000
pct_mt_high = 20

print(f"[INFO] QC thresholds: n_genes {n_genes_low}-{n_genes_high}, total_counts < {total_counts_high}, pct_mt < {pct_mt_high}")
adata = adata[
    (adata.obs['n_genes_by_counts'] >= n_genes_low) &
    (adata.obs['n_genes_by_counts'] <= n_genes_high) &
    (adata.obs['total_counts'] <= total_counts_high) &
    (adata.obs['pct_counts_mt'] <= pct_mt_high)
]
print(f"[INFO] After QC: {adata.n_obs} cells")

# 过滤基因
sc.pp.filter_genes(adata, min_cells=10)
print(f"[INFO] After gene filter: {adata.n_vars} genes")

# 保存 counts
adata.layers['counts'] = adata.X.copy()
adata.raw = adata.copy()

# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 高变基因
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
adata_hvg = adata[:, adata.var.highly_variable].copy()

# 缩放、PCA
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')

# neighbors, UMAP, leiden
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.8, flavor='igraph', n_iterations=2)

# 传回
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obs['leiden'] = adata_hvg.obs['leiden']

# 可视化
sc.settings.set_figure_params(dpi=150, facecolor='white')
for color in ['time_point', 'leiden']:
    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, size=5)
    fig.savefig(os.path.join(OUT_DIR, f"umap_{color}.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

# Pirb 表达
if 'Pirb' in adata.var_names:
    adata.obs['Pirb_expr'] = adata[:, 'Pirb'].X.toarray().flatten()
    adata.obs['Pirb_positive'] = (adata.obs['Pirb_expr'] > 0).astype(int)
    
    summary = adata.obs.groupby('time_point').agg(
        n=('Pirb_positive', 'size'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean')
    ).reset_index()
    summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_overall.csv"), index=False)
    print(summary)
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    order = ['sham', 'D1', 'D3', 'D7']
    sns.violinplot(data=adata.obs, x='time_point', y='Pirb_expr', order=order, ax=ax)
    ax.set_title('Pirb expression across time course')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_expression_time.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    fig = sc.pl.umap(adata, color='Pirb_expr', show=False, return_fig=True, color_map='viridis', size=5)
    fig.savefig(os.path.join(OUT_DIR, "umap_Pirb_expr.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

# 保存
adata.write_h5ad(os.path.join(OUT_DIR, "../GSE233812_processed.h5ad"))
print("[SAVE] GSE233812_processed.h5ad")
print("[DONE]")
