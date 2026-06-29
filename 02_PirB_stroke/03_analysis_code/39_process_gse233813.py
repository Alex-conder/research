"""
处理 GSE233813 snRNA-seq 数据（sham / D1 / D3 / D7）
"""
import os, glob, gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import io, sparse
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE233813"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233813"
os.makedirs(OUT_DIR, exist_ok=True)

samples = {
    "GSM7437217_sn-sham": "sham",
    "GSM7437218_sn-D1": "D1",
    "GSM7437219_sn-D3": "D3",
    "GSM7437220_sn-D7": "D7",
}

def read_10x_custom(prefix):
    barcodes = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None, compression='gzip')[0].values
    genes = pd.read_csv(f"{prefix}_genes.tsv.gz", header=None, compression='gzip', sep='\t')
    genes.columns = ['gene_id', 'gene_name']
    mat = io.mmread(f"{prefix}_matrix.mtx.gz")
    mat = sparse.csr_matrix(mat.T)
    var = pd.DataFrame(index=genes['gene_name'].values)
    var['gene_id'] = genes['gene_id'].values
    ad = anndata.AnnData(X=mat, obs=pd.DataFrame(index=barcodes), var=var)
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

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"[INFO] Before cell filter: {adata.n_obs}")
sc.pp.filter_cells(adata, min_genes=200)
print(f"[INFO] After min_genes=200: {adata.n_obs}")

# snRNA-seq QC 阈值（核内 mt 通常较低）
n_genes_low = 500
n_genes_high = 8000
total_counts_high = 30000
pct_mt_high = 5

adata = adata[
    (adata.obs['n_genes_by_counts'] >= n_genes_low) &
    (adata.obs['n_genes_by_counts'] <= n_genes_high) &
    (adata.obs['total_counts'] <= total_counts_high) &
    (adata.obs['pct_counts_mt'] <= pct_mt_high)
]
print(f"[INFO] After QC: {adata.n_obs} cells")

sc.pp.filter_genes(adata, min_cells=10)
print(f"[INFO] After gene filter: {adata.n_vars} genes")

adata.layers['counts'] = adata.X.copy()
adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.8, flavor='igraph', n_iterations=2)

adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obs['leiden'] = adata_hvg.obs['leiden']

# 细胞类型注释
markers = {
    'Neuron': ['Snap25', 'Rbfox3', 'Syn1'],
    'Astrocyte': ['Gfap', 'Aqp4', 'Aldh1l1'],
    'Microglia': ['C1qa', 'C1qb', 'P2ry12', 'Tmem119'],
    'Oligodendrocyte': ['Mbp', 'Mog', 'Plp1'],
    'OPC': ['Pdgfra', 'Cspg4'],
    'Endothelial': ['Pecam1', 'Cldn5'],
    'Pericyte': ['Pdgfrb', 'Rgs5'],
}
for ct, genes in markers.items():
    available = [g for g in genes if g in adata.var_names]
    if len(available) == 0:
        continue
    sc.tl.score_genes(adata, gene_list=available, score_name=f'{ct}_score')

score_cols = [f'{ct}_score' for ct in markers.keys() if f'{ct}_score' in adata.obs.columns]
cluster_ct = adata.obs[score_cols].groupby(adata.obs['leiden']).mean().idxmax(axis=1).str.replace('_score', '').to_dict()
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_ct)

print("[INFO] Cell types:")
print(adata.obs.groupby('cell_type').size())

# 可视化
sc.settings.set_figure_params(dpi=150, facecolor='white')
for color in ['time_point', 'cell_type', 'leiden']:
    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, size=5)
    fig.savefig(os.path.join(OUT_DIR, f"umap_{color}.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

# Pirb
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
    
    ct_summary = adata.obs.groupby(['cell_type', 'time_point']).agg(
        n=('Pirb_positive', 'size'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean')
    ).reset_index()
    ct_summary = ct_summary[ct_summary['n'] >= 10]
    ct_summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_celltype_time.csv"), index=False)
    print(ct_summary)
    
    pivot = ct_summary.pivot_table(index='cell_type', columns='time_point', values='pirb_frac', fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot[['sham', 'D1', 'D3', 'D7']], annot=True, fmt='.3f', cmap='YlOrRd', ax=ax)
    ax.set_title('Pirb+ fraction by cell type and time (snRNA-seq)')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_fraction_celltype_time.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    fig = sc.pl.umap(adata, color='Pirb_expr', show=False, return_fig=True, color_map='viridis', size=5)
    fig.savefig(os.path.join(OUT_DIR, "umap_Pirb_expr.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

adata.write_h5ad(os.path.join(OUT_DIR, "../GSE233813_processed.h5ad"))
print("[SAVE] GSE233813_processed.h5ad")
print("[DONE]")
