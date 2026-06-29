"""
对 GSE233812 进行细胞类型注释并分析 Pirb 时间动态
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE233812_processed.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233812"

adata = sc.read_h5ad(IN_H5AD)

# 小鼠脑细胞类型 marker
markers = {
    'Neuron': ['Snap25', 'Rbfox3', 'Syn1', 'Stmn2'],
    'Astrocyte': ['Gfap', 'Aqp4', 'Aldh1l1', 'Slc1a3'],
    'Microglia': ['C1qa', 'C1qb', 'P2ry12', 'Tmem119', 'Cx3cr1'],
    'Oligodendrocyte': ['Mbp', 'Mog', 'Plp1', 'Olig1'],
    'OPC': ['Pdgfra', 'Cspg4'],
    'Endothelial': ['Pecam1', 'Cldn5', 'Slc2a1'],
    'Pericyte': ['Pdgfrb', 'Rgs5', 'Cspg4'],
    'Immune': ['Ptprc', 'Lyz2', 'Cd3e', 'Cd19'],
}

# 计算每个 cluster 的 marker 得分
for ct, genes in markers.items():
    available = [g for g in genes if g in adata.var_names]
    if len(available) == 0:
        print(f"[WARN] No marker genes for {ct}")
        continue
    sc.tl.score_genes(adata, gene_list=available, score_name=f'{ct}_score')

# 为每个 leiden cluster 分配细胞类型
scores = adata.obs[[f'{ct}_score' for ct in markers.keys() if f'{ct}_score' in adata.obs.columns]]
cluster_ct = scores.groupby(adata.obs['leiden']).mean().idxmax(axis=1).str.replace('_score', '')
cluster_ct = cluster_ct.to_dict()
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_ct)

print("[INFO] Cell type assignment:")
print(adata.obs.groupby('cell_type').size())

# 可视化
sc.settings.set_figure_params(dpi=150, facecolor='white')
fig = sc.pl.umap(adata, color='cell_type', show=False, return_fig=True, size=10, legend_loc='right margin')
fig.savefig(os.path.join(OUT_DIR, "umap_cell_type.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

# 时间分布
fig, ax = plt.subplots(figsize=(10, 6))
ct_time = pd.crosstab(adata.obs['cell_type'], adata.obs['time_point'], normalize='columns')
sns.heatmap(ct_time, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax)
ax.set_title('Cell type composition by time point')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "celltype_composition_time.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

# Pirb 表达分析
if 'Pirb' in adata.var_names:
    # 按细胞类型和时间
    pirb_summary = adata.obs.groupby(['cell_type', 'time_point']).agg(
        n=('Pirb_positive', 'size'),
        pirb_frac=('Pirb_positive', 'mean'),
        pirb_mean=('Pirb_expr', 'mean'),
    ).reset_index()
    pirb_summary = pirb_summary[pirb_summary['n'] >= 10]
    pirb_summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_celltype_time.csv"), index=False)
    print(pirb_summary)
    
    # 热图
    pivot = pirb_summary.pivot_table(index='cell_type', columns='time_point', values='pirb_frac', fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot[['sham', 'D1', 'D3', 'D7']], annot=True, fmt='.3f', cmap='YlOrRd', ax=ax)
    ax.set_title('Pirb+ fraction by cell type and time')
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_fraction_celltype_time.png"), dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    # 小提琴图
    key_cts = pirb_summary[pirb_summary['time_point'] == 'D3'].sort_values('pirb_frac', ascending=False)['cell_type'].tolist()
    for ct in key_cts:
        subset = adata[adata.obs['cell_type'] == ct]
        if subset.n_obs < 20:
            continue
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.violinplot(data=subset.obs, x='time_point', y='Pirb_expr', order=['sham', 'D1', 'D3', 'D7'], ax=ax)
        ax.set_title(f'Pirb expression in {ct}')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"pirb_{ct}_time.png"), dpi=200, bbox_inches='tight')
        plt.close(fig)

# 保存
adata.write_h5ad(IN_H5AD)
print("[DONE]")
