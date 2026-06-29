"""
从 GSE225948_merged_raw.h5ad 加载，完成降维、聚类、Pirb 分析
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

IN_H5AD = "../04_reports/figures/GSE225948_merged_raw.h5ad"
OUT_DIR = "../04_reports/figures/GSE225948"
os.makedirs(OUT_DIR, exist_ok=True)

print("[INFO] Loading merged raw data...")
adata = sc.read_h5ad(IN_H5AD)
print(adata)

# 手动高变基因选择
print("[INFO] Selecting highly variable genes manually...")
means = np.array(adata.X.mean(axis=0)).flatten()
Xsq = adata.X.copy()
Xsq.data = Xsq.data ** 2
vars_ = np.array(Xsq.mean(axis=0)).flatten() - means ** 2
valid = (means > 0) & np.isfinite(vars_) & (vars_ > 0)
dispersion = vars_ / (means + 1e-10)
gene_rank = np.full(len(means), -1.0)
gene_rank[valid] = dispersion[valid]
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
sc.tl.leiden(adata_hvg, resolution=0.8, flavor='igraph', n_iterations=2)

adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']
adata.obs['leiden'] = adata_hvg.obs['leiden']

# 可视化
print("[INFO] Plotting UMAPs...")
sc.settings.set_figure_params(dpi=150, facecolor='white')
for color in ['parent_celltype', 'sub_celltype', 'tissue', 'time_point', 'tissue_time', 'leiden']:
    if color in adata.obs.columns:
        try:
            if color == 'leiden':
                fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, size=5, legend_loc='on data')
            else:
                fig = sc.pl.umap(adata, color=color, show=False, return_fig=True, size=5, legend_loc='right margin')
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
    
    summary = adata.obs.groupby(['tissue', 'time_point', 'parent_celltype'], observed=False).agg(
        n_cells=('Pirb_positive', 'size'),
        pirb_positive=('Pirb_positive', 'sum'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean'),
        pirb_median=('Pirb_expr', 'median'),
    ).reset_index()
    summary = summary[summary['n_cells'] > 0]
    summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_by_group.csv"), index=False)
    print(summary.head(40))
    
    # 小提琴图：按 parent cell type
    fig, axes = plt.subplots(1, 2, figsize=(20, 6))
    sns.violinplot(data=adata.obs, x='parent_celltype', y='Pirb_expr', ax=axes[0], order=sorted(adata.obs['parent_celltype'].unique()))
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=90)
    axes[0].set_title('Pirb expression by parent cell type')
    sns.violinplot(data=adata.obs, x='tissue_time', y='Pirb_expr', ax=axes[1], order=sorted(adata.obs['tissue_time'].unique()))
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=90)
    axes[1].set_title('Pirb expression by tissue and time')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_expression_violin.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    # 关键细胞类型的时间动态
    key_cts = ['Mg', 'EC', 'MdC', 'Mo', 'Bc', 'Neu', 'Tc', 'DC', 'NK']
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
        order = sorted(subset.obs['age_time'].unique())
        sns.violinplot(data=subset.obs, x='age_time', y='Pirb_expr', ax=ax, order=order)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        ax.set_title(f'Pirb in brain {ct}: young vs aged')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"pirb_{ct}_brain_age_time.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)
    
    # UMAP 上 Pirb 表达
    fig = sc.pl.umap(adata, color='Pirb_expr', show=False, return_fig=True, color_map='viridis', size=5)
    fig.savefig(os.path.join(OUT_DIR, "umap_Pirb_expr.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    # 热图：Pirb+ 比例
    pivot = summary.pivot_table(index='parent_celltype', columns=['tissue', 'time_point'], values='pirb_frac', fill_value=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax)
    ax.set_title('Pirb+ fraction by cell type, tissue and time')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_fraction_heatmap.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)

# 保存最终对象
adata.write_h5ad(os.path.join(OUT_DIR, "../GSE225948_processed.h5ad"))
print("[SAVE] GSE225948_processed.h5ad")
print("[DONE]")
