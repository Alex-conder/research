"""
处理 GSE225948 单细胞数据（缺血后 CD45hi + Microglia + EC）
数据格式：CSV 行=基因，列=细胞，值为 normalized expression
"""
import os, glob, gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

RAW_DIR = "../01_raw_data/GSE225948"
OUT_DIR = "../04_reports/figures/GSE225948"
os.makedirs(OUT_DIR, exist_ok=True)

# 找到所有成对的 counts + metadata
counts_files = sorted(glob.glob(os.path.join(RAW_DIR, "*_counts.csv.gz")))
print(f"[INFO] Found {len(counts_files)} counts files")

adatas = []

for cf in counts_files:
    base = os.path.basename(cf).replace("_counts.csv.gz", "")
    mf = os.path.join(RAW_DIR, base + "_metadata.csv.gz")
    if not os.path.exists(mf):
        print(f"[SKIP] {base}: metadata missing")
        continue
    
    print(f"[READ] {base}")
    try:
        df_counts = pd.read_csv(cf, index_col=0)
    except Exception as e:
        print(f"[SKIP] {base}: corrupt or unreadable counts ({e})")
        continue
    
    # 转置：行=细胞，列=基因
    df_counts = df_counts.T
    
    # 读取 metadata
    meta = pd.read_csv(mf, index_col=0)
    
    # metadata 的 index 需要与 counts 的 index 匹配
    common_cells = df_counts.index.intersection(meta.index)
    if len(common_cells) == 0:
        df_counts.index = [str(x).strip('"') for x in df_counts.index]
        meta.index = [str(x).strip('"') for x in meta.index]
        common_cells = df_counts.index.intersection(meta.index)
    
    if len(common_cells) == 0:
        print(f"[WARN] {base}: no common cells")
        continue
    
    df_counts = df_counts.loc[common_cells]
    meta = meta.loc[common_cells]
    
    # 构建稀疏矩阵
    X = sparse.csr_matrix(df_counts.values.astype(np.float32))
    var = pd.DataFrame(index=df_counts.columns)
    
    ad = anndata.AnnData(X=X, obs=meta, var=var)
    ad.obs['sample'] = base
    ad.obs['GSM'] = base.split('_')[0]
    
    col_map = {
        'tissue': 'tissue',
        'treatment': 'treatment',
        'age': 'age',
        'sex': 'sex',
        'sub.celltype': 'sub_celltype',
        'parent': 'parent_celltype',
        'nCount_RNA': 'nCount_RNA',
        'nFeature_RNA': 'nFeature_RNA',
        'percent.mt': 'percent_mt',
        'Sample_description': 'sample_description',
    }
    for old, new in col_map.items():
        if old in ad.obs.columns:
            ad.obs[new] = ad.obs[old].astype(str)
    
    adatas.append(ad)
    print(f"[DONE] {base}: {ad.n_obs} cells x {ad.n_vars} genes")
    del df_counts, meta, X
    gc.collect()

if len(adatas) == 0:
    raise ValueError("No valid samples loaded")

print("[INFO] Concatenating...")
adata = anndata.concat(adatas, axis=0, join='outer', fill_value=0)
adata.X = sparse.csr_matrix(adata.X)
print(f"[INFO] Merged: {adata.n_obs} cells x {adata.n_vars} genes")

# 时间推断
adata.obs['time_point'] = 'Unknown'

def parse_desc(desc):
    desc = str(desc).lower()
    if '2d' in desc or 'd02' in desc or 'day 2' in desc:
        return 'D02'
    elif '14d' in desc or 'd14' in desc or 'day 14' in desc:
        return 'D14'
    elif 'sham' in desc:
        return 'Sham'
    return 'Unknown'

tp_desc = adata.obs['sample_description'].apply(parse_desc)
adata.obs.loc[tp_desc != 'Unknown', 'time_point'] = tp_desc[tp_desc != 'Unknown']

tp_treat = adata.obs['treatment'].astype(str).str.upper()
adata.obs.loc[(adata.obs['time_point'] == 'Unknown') & (tp_treat.isin(['D02', 'D2'])), 'time_point'] = 'D02'
adata.obs.loc[(adata.obs['time_point'] == 'Unknown') & (tp_treat == 'D14'), 'time_point'] = 'D14'
adata.obs.loc[(adata.obs['time_point'] == 'Unknown') & (tp_treat == 'SHAM'), 'time_point'] = 'Sham'

adata.obs['tissue_time'] = adata.obs['tissue'].astype(str) + '_' + adata.obs['time_point'].astype(str)
adata.obs['age_group'] = adata.obs['age'].astype(str).apply(lambda x: 'aged' if x.startswith('M') else 'young')

print("[INFO] Time point distribution:")
print(adata.obs['time_point'].value_counts())
print("[INFO] Tissue distribution:")
print(adata.obs['tissue'].value_counts())
print("[INFO] Parent cell type distribution:")
print(adata.obs['parent_celltype'].value_counts())

# QC
print("[INFO] QC filtering...")
adata.obs['nFeature_RNA'] = pd.to_numeric(adata.obs['nFeature_RNA'], errors='coerce')
adata.obs['percent_mt'] = pd.to_numeric(adata.obs['percent_mt'], errors='coerce')
print(f"  Before QC: {adata.n_obs} cells")
adata = adata[(adata.obs['nFeature_RNA'] >= 200) & 
              (adata.obs['nFeature_RNA'] <= 8000) &
              (adata.obs['percent_mt'] <= 25)]
print(f"  After QC: {adata.n_obs} cells")

# 过滤低表达基因：至少在 50 个细胞中表达
sc.pp.filter_genes(adata, min_cells=50)
print(f"[INFO] After gene filter: {adata.n_vars} genes")

# 保存原始合并对象
adata.write_h5ad(os.path.join(OUT_DIR, "../GSE225948_merged_raw.h5ad"))
print("[SAVE] GSE225948_merged_raw.h5ad")

# 手动高变基因选择（避免 scanpy flavor 的 expm1 overflow）
print("[INFO] Selecting highly variable genes manually...")
means = np.array(adata.X.mean(axis=0)).flatten()
# 使用非零值的二阶矩近似
Xsq = adata.X.copy()
Xsq.data = Xsq.data ** 2
vars_ = np.array(Xsq.mean(axis=0)).flatten() - means ** 2
# 排除 mean 为 0 或 var 为 NaN/inf 的基因
valid = (means > 0) & np.isfinite(vars_) & (vars_ > 0)
# 按 mean 分层选高变基因
n_bins = 20
mean_bins = pd.qcut(means[valid], q=n_bins, duplicates='drop', labels=False)
dispersion = vars_ / (means + 1e-10)
# 在每个 bin 中选 dispersion 最高的
gene_rank = np.full(len(means), -1)
gene_rank[valid] = dispersion[valid]
# 选择 top 2000
threshold = np.sort(gene_rank[valid])[-2000]
adata.var['highly_variable'] = gene_rank >= threshold
print(f"[INFO] Selected {adata.var['highly_variable'].sum()} HVGs")

adata_hvg = adata[:, adata.var['highly_variable']].copy()

print("[INFO] Scaling and PCA...")
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')

print("[INFO] Computing neighbors and UMAP...")
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.8)

adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obs['leiden'] = adata_hvg.obs['leiden']

# 可视化
print("[INFO] Plotting UMAPs...")
sc.settings.set_figure_params(dpi=150, facecolor='white')
for color in ['parent_celltype', 'sub_celltype', 'tissue', 'time_point', 'tissue_time', 'leiden']:
    if color in adata.obs.columns:
        try:
            fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, size=5, legend_loc='on data' if color == 'leiden' else 'right margin')
            fig.savefig(os.path.join(OUT_DIR, f"umap_{color}.png"), dpi=200, bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            print(f"[WARN] UMAP {color} plot failed: {e}")

# Pirb 表达分析
print("[INFO] Pirb expression analysis...")
if 'Pirb' in adata.var_names:
    pirb_expr = np.array(adata[:, 'Pirb'].X.toarray()).flatten()
    adata.obs['Pirb_expr'] = pirb_expr
    adata.obs['Pirb_positive'] = (pirb_expr > 0).astype(int)
    
    # 总体 summary
    summary = adata.obs.groupby(['tissue', 'time_point', 'parent_celltype']).agg(
        n_cells=('Pirb_positive', 'size'),
        pirb_positive=('Pirb_positive', 'sum'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean'),
        pirb_median=('Pirb_expr', 'median'),
    ).reset_index()
    summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_by_group.csv"), index=False)
    print(summary.head(30))
    
    # 小提琴图：按 parent cell type
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    sns.violinplot(data=adata.obs, x='parent_celltype', y='Pirb_expr', ax=axes[0], order=sorted(adata.obs['parent_celltype'].unique()))
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90)
    axes[0].set_title('Pirb expression by parent cell type')
    sns.violinplot(data=adata.obs, x='tissue_time', y='Pirb_expr', ax=axes[1])
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90)
    axes[1].set_title('Pirb expression by tissue and time')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_expression_violin.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    # 关键细胞类型的时间动态
    key_cts = ['Mg', 'EC', 'MdC', 'Mo', 'Bc', 'Neu', 'Tc']
    for ct in key_cts:
        subset = adata[adata.obs['parent_celltype'] == ct]
        if subset.n_obs < 100:
            continue
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.violinplot(data=subset.obs, x='time_point', y='Pirb_expr', ax=ax, order=['Sham', 'D02', 'D14'])
        ax.set_title(f'Pirb expression in {ct} across time')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"pirb_{ct}_time.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)
    
    # Brain 年轻 vs 年老
    for ct in ['Mg', 'EC']:
        subset = adata[(adata.obs['parent_celltype'] == ct) & (adata.obs['tissue'] == 'brain')]
        if subset.n_obs < 100:
            continue
        subset.obs['age_time'] = subset.obs['age_group'].astype(str) + '_' + subset.obs['time_point'].astype(str)
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.violinplot(data=subset.obs, x='age_time', y='Pirb_expr', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title(f'Pirb in brain {ct}: young vs aged')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"pirb_{ct}_brain_age_time.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)
    
    # UMAP 上 Pirb 表达
    fig = sc.pl.umap(adata, color='Pirb_expr', show=False, return_fig=True, color_map='viridis', size=5)
    fig.savefig(os.path.join(OUT_DIR, "umap_Pirb_expr.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

# 保存最终对象
adata.write_h5ad(os.path.join(OUT_DIR, "../GSE225948_processed.h5ad"))
print("[SAVE] GSE225948_processed.h5ad")
print("[DONE]")
