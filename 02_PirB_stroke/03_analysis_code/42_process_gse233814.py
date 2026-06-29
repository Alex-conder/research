"""
处理 GSE233814 Visium 空间转录组（control / D1 / D3 / D7）
基础分析：Pirb 表达与空间分布
"""
import os, json, gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import io, sparse
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE233814"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233814"
os.makedirs(OUT_DIR, exist_ok=True)

samples = {
    "GSM7437221_C1-control": "control",
    "GSM7437222_B1-D1": "D1",
    "GSM7437223_D1-D3": "D3",
    "GSM7437224_C1-D7": "D7",
    "GSM7437225_D1-D7": "D7_rep",
}

def read_10x_custom(prefix):
    barcodes = pd.read_csv(f"{prefix}_barcodes.tsv.gz", header=None, compression='gzip')[0].values
    genes = pd.read_csv(f"{prefix}_features.tsv.gz", header=None, compression='gzip', sep='\t')
    genes.columns = ['gene_id', 'gene_name', 'feature_type']
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
    print(f"  Loaded {ad.n_obs} spots x {ad.n_vars} genes")
    adatas.append(ad)
    gc.collect()

print("[INFO] Concatenating...")
adata = anndata.concat(adatas, axis=0, join='outer', fill_value=0)
adata.X = adata.X.tocsr()
print(f"[INFO] Merged: {adata.n_obs} spots x {adata.n_vars} genes")

# 基础 QC
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"[INFO] Before QC: {adata.n_obs}")
adata = adata[(adata.obs['n_genes_by_counts'] >= 200) & (adata.obs['total_counts'] >= 500)]
print(f"[INFO] After QC: {adata.n_obs}")

sc.pp.filter_genes(adata, min_cells=5)
print(f"[INFO] After gene filter: {adata.n_vars} genes")

# 标准化
adata.layers['counts'] = adata.X.copy()
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Pirb 分析
if 'Pirb' in adata.var_names:
    adata.obs['Pirb_expr'] = adata[:, 'Pirb'].X.toarray().flatten()
    adata.obs['Pirb_positive'] = (adata.obs['Pirb_expr'] > 0).astype(int)
    
    summary = adata.obs.groupby('time_point').agg(
        n=('Pirb_positive', 'size'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean'),
    ).reset_index()
    summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_by_time.csv"), index=False)
    print(summary)
    
    # 小提琴图
    fig, ax = plt.subplots(figsize=(8, 5))
    order = ['control', 'D1', 'D3', 'D7', 'D7_rep']
    sns.violinplot(data=adata.obs, x='time_point', y='Pirb_expr', order=order, ax=ax)
    ax.set_title('Pirb expression in Visium spots across time')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_expression_violin.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    # 每个样本的 Pirb 阳性 fraction
    sample_summary = adata.obs.groupby('sample').agg(
        n=('Pirb_positive', 'size'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean'),
    ).reset_index()
    sample_summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_by_sample.csv"), index=False)
    print(sample_summary)

# 保存
adata.write_h5ad(os.path.join(OUT_DIR, "../GSE233814_processed.h5ad"))
print("[SAVE] GSE233814_processed.h5ad")
print("[DONE]")
