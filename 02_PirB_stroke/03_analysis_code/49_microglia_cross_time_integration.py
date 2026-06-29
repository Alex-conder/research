"""
路线二：跨时间小胶质细胞整合分析
整合 GSE174574（Sham / MCAO 24h）+ GSE233812 scRNA-seq（sham/D1/D3/D7）
+ GSE233813 snRNA-seq（sham/D1/D3/D7）中的小胶质细胞，构建 Pirb+ 状态转移图谱。
"""
import os, gc
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')

OUT_DIR = "../04_reports/figures/microglia_cross_time"
os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. 加载各数据集并提取小胶质细胞
# -----------------------------------------------------------------------------
print("[INFO] Loading microglia from GSE174574...")
adata_174 = sc.read_h5ad("../04_reports/figures/GSE174574/GSE174574_annotated.h5ad")
mg_174 = adata_174[adata_174.obs['cell_type'] == 'Microglia'].copy()
mg_174.obs['dataset'] = 'GSE174574'
mg_174.obs['time_point'] = mg_174.obs['condition'].map({'Sham': 'sham', 'MCAO': 'MCAO_24h'})
mg_174.obs['time_ordered'] = mg_174.obs['time_point'].map({'sham': 0, 'MCAO_24h': 1})
mg_174.obs['Pirb_expr'] = mg_174.obs['Pirb_counts'].values
mg_174.obs['Pirb_positive'] = (mg_174.obs['Pirb_expr'] > 0).astype(int)
print(f"  GSE174574 microglia: {mg_174.n_obs}")
del adata_174; gc.collect()

print("[INFO] Loading microglia from GSE233812...")
adata_812 = sc.read_h5ad("../04_reports/figures/GSE233812_processed.h5ad")
mg_812 = adata_812[adata_812.obs['cell_type'] == 'Microglia'].copy()
mg_812.obs['dataset'] = 'GSE233812_sc'
mg_812.obs['time_ordered'] = mg_812.obs['time_point'].map({'sham': 0, 'D1': 1, 'D3': 2, 'D7': 3})
print(f"  GSE233812 microglia: {mg_812.n_obs}")
del adata_812; gc.collect()

print("[INFO] Loading microglia from GSE233813...")
adata_813 = sc.read_h5ad("../04_reports/figures/GSE233813_processed.h5ad")
mg_813 = adata_813[adata_813.obs['cell_type'] == 'Microglia'].copy()
mg_813.obs['dataset'] = 'GSE233813_sn'
mg_813.obs['time_ordered'] = mg_813.obs['time_point'].map({'sham': 0, 'D1': 1, 'D3': 2, 'D7': 3})
print(f"  GSE233813 microglia: {mg_813.n_obs}")
del adata_813; gc.collect()

# -----------------------------------------------------------------------------
# 2. 统一变量名并合并
# -----------------------------------------------------------------------------
print("[INFO] Merging microglia datasets...")
# 取共有基因
common_genes = list(set(mg_174.var_names) & set(mg_812.var_names) & set(mg_813.var_names))
print(f"  Common genes: {len(common_genes)}")

mg_174 = mg_174[:, common_genes].copy()
mg_812 = mg_812[:, common_genes].copy()
mg_813 = mg_813[:, common_genes].copy()

# 分别归一化：GSE233812/813 使用 counts layer；GSE174574 的 X 已经是 log1p normalized，直接使用
print("[INFO] Normalizing each dataset...")
for ad, name in [(mg_812, 'GSE233812_sc'), (mg_813, 'GSE233813_sn')]:
    if 'counts' in ad.layers:
        ad.X = ad.layers['counts'].copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
# GSE174574 已经是 log1p-normalized，保持不变

# 保留关键 obs 列
keep_obs = ['dataset', 'time_point', 'time_ordered', 'Pirb_expr', 'Pirb_positive',
            'n_genes_by_counts', 'total_counts', 'pct_counts_mt']
for ad in [mg_174, mg_812, mg_813]:
    for col in keep_obs:
        if col not in ad.obs.columns:
            ad.obs[col] = np.nan
    ad.obs = ad.obs[keep_obs].copy()

# 合并：所有数据集 X 均已是 log1p-normalized
mg_all = anndata.concat([mg_174, mg_812, mg_813], axis=0, join='outer', label='dataset',
                        keys=['GSE174574', 'GSE233812_sc', 'GSE233813_sn'], fill_value=0)
mg_all.X = mg_all.X.tocsr() if hasattr(mg_all.X, 'toarray') else mg_all.X
print(f"[INFO] Merged microglia: {mg_all.n_obs} cells x {mg_all.n_vars} genes")

# 统一 time_point
mg_all.obs['time_point'] = mg_all.obs['time_point'].astype(str)
mg_all.obs['time_point'] = mg_all.obs['time_point'].replace('nan', 'sham')

# 创建统一的时间标签：sham / D1 / D3 / D7 / MCAO_24h
mg_all.obs['time_label'] = mg_all.obs['time_point']

# -----------------------------------------------------------------------------
# 3. HVG + 降维聚类
# -----------------------------------------------------------------------------
print("[INFO] HVG, scaling...")

sc.pp.highly_variable_genes(mg_all, n_top_genes=2000, flavor='seurat_v3')
# 强制保留 Pirb 基因用于后续分析
if 'Pirb' in mg_all.var_names:
    mg_all.var.loc['Pirb', 'highly_variable'] = True
hvg_mask = mg_all.var['highly_variable'].values
mg_all = mg_all[:, hvg_mask].copy()
print(f"[INFO] HVG selected: {mg_all.n_vars} (including Pirb)")

# 统一时间类别顺序
order_time = ['sham', 'D1', 'MCAO_24h', 'D3', 'D7']
mg_all.obs['time_point'] = pd.Categorical(mg_all.obs['time_point'], categories=order_time, ordered=True)

sc.pp.scale(mg_all, max_value=10)
sc.tl.pca(mg_all, svd_solver='arpack')

# 批次校正：combat
print("[INFO] Running ComBat batch correction...")
sc.pp.combat(mg_all, key='dataset')

sc.pp.neighbors(mg_all, n_neighbors=15, n_pcs=20)
sc.tl.umap(mg_all)
sc.tl.leiden(mg_all, resolution=0.6)

# -----------------------------------------------------------------------------
# 4. Pirb 表达与细胞状态分析
# -----------------------------------------------------------------------------
print("[INFO] Pirb expression analysis...")
mg_all.obs['Pirb_expr'] = mg_all[:, 'Pirb'].X.toarray().flatten()
mg_all.obs['Pirb_positive'] = (mg_all.obs['Pirb_expr'] > 0).astype(int)

# 按数据集和时间汇总
summary = mg_all.obs.groupby(['dataset', 'time_point']).agg(
    n=('Pirb_positive', 'size'),
    pirb_pos=('Pirb_positive', 'sum'),
    pirb_frac=('Pirb_positive', 'mean'),
    pirb_mean=('Pirb_expr', 'mean'),
).reset_index()
summary.to_csv(os.path.join(OUT_DIR, "microglia_pirb_summary_by_dataset_time.csv"), index=False)
print(summary.to_string(index=False))

# 全局按时间汇总
global_summary = mg_all.obs.groupby('time_point').agg(
    n=('Pirb_positive', 'size'),
    pirb_pos=('Pirb_positive', 'sum'),
    pirb_frac=('Pirb_positive', 'mean'),
    pirb_mean=('Pirb_expr', 'mean'),
).reset_index()
# 按时间顺序
order_time = ['sham', 'D1', 'MCAO_24h', 'D3', 'D7']
global_summary['time_point'] = pd.Categorical(global_summary['time_point'].astype(str), categories=order_time, ordered=True)
global_summary = global_summary.sort_values('time_point')
global_summary.to_csv(os.path.join(OUT_DIR, "microglia_pirb_global_summary.csv"), index=False)
print("\nGlobal summary:")
print(global_summary.to_string(index=False))

# -----------------------------------------------------------------------------
# 5. 可视化
# -----------------------------------------------------------------------------
print("[INFO] Plotting...")
order_time = ['sham', 'D1', 'MCAO_24h', 'D3', 'D7']

# (a) UMAP 按数据集着色
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(mg_all, color='dataset', ax=ax, show=False, title='Microglia by dataset')
fig.savefig(os.path.join(OUT_DIR, "umap_dataset.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (b) UMAP 按时间着色
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(mg_all, color='time_point', ax=ax, show=False, title='Microglia by time point',
           palette='tab10')
fig.savefig(os.path.join(OUT_DIR, "umap_time_point.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (c) UMAP 按 Pirb 表达着色
cmap_pirb = LinearSegmentedColormap.from_list("pirb", ["#000000", "#1f4e79", "#ffd700", "#ff0000"])
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(mg_all, color='Pirb_expr', cmap=cmap_pirb, ax=ax, show=False,
           title='Pirb expression (log1p CPM)', vmin=0, vmax='p99')
fig.savefig(os.path.join(OUT_DIR, "umap_pirb_expr.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (d) UMAP 按 leiden 聚类着色
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(mg_all, color='leiden', ax=ax, show=False, title='Leiden clusters')
fig.savefig(os.path.join(OUT_DIR, "umap_leiden.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (e) 每个 cluster 的时间组成
cluster_time = pd.crosstab(mg_all.obs['leiden'], mg_all.obs['time_point'], normalize='index')
time_order_cols = ['sham', 'D1', 'MCAO_24h', 'D3', 'D7']
cluster_time = cluster_time[time_order_cols] if all(t in cluster_time.columns for t in time_order_cols) else cluster_time
fig, ax = plt.subplots(figsize=(10, 6))
sns.heatmap(cluster_time, cmap='YlOrRd', annot=True, fmt='.2f', ax=ax)
ax.set_title("Microglia cluster composition by time point")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "cluster_time_composition.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (f) 每个 cluster 的 Pirb 阳性率
cluster_pirb = mg_all.obs.groupby('leiden')['Pirb_positive'].agg(['sum', 'count', 'mean']).reset_index()
cluster_pirb.columns = ['leiden', 'pirb_pos_n', 'total_n', 'pirb_frac']
cluster_pirb.to_csv(os.path.join(OUT_DIR, "cluster_pirb_fraction.csv"), index=False)
fig, ax = plt.subplots(figsize=(8, 5))
sns.barplot(data=cluster_pirb, x='leiden', y='pirb_frac', order=cluster_pirb.sort_values('pirb_frac', ascending=False)['leiden'], ax=ax)
ax.set_title("Pirb+ fraction by microglia cluster")
ax.set_ylabel("Pirb+ fraction")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "cluster_pirb_fraction_barplot.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# (g) 时间动态折线图
fig, ax = plt.subplots(figsize=(8, 5))
sns.pointplot(data=global_summary, x='time_point', y='pirb_frac', order=order_time, color='darkred', ax=ax)
ax.set_title("Cross-time Pirb+ fraction in microglia")
ax.set_ylabel("Pirb+ fraction")
ax.set_xlabel("Time point")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_fraction_timeline.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# -----------------------------------------------------------------------------
# 6. 拟时序分析（diffmap + dpt）
# -----------------------------------------------------------------------------
print("[INFO] Pseudotime analysis...")
sc.tl.diffmap(mg_all, n_comps=15)

# 以 sham 中位数细胞为根
sham_cells = mg_all.obs_names[mg_all.obs['time_point'] == 'sham']
if len(sham_cells) > 0:
    # 使用 diffmap 坐标找 sham 细胞中最中心的作为根
    diff_coords = mg_all.obsm['X_diffmap']
    sham_idx = np.where(mg_all.obs['time_point'].values == 'sham')[0]
    sham_medoid = sham_idx[np.argmin(np.linalg.norm(diff_coords[sham_idx] - diff_coords[sham_idx].mean(axis=0), axis=1))]
    mg_all.uns['iroot'] = sham_medoid
    sc.tl.dpt(mg_all, n_dcs=10)

    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(mg_all, color='dpt_pseudotime', ax=ax, show=False, title='Diffusion pseudotime')
    fig.savefig(os.path.join(OUT_DIR, "umap_pseudotime.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Pirb 表达沿拟时序的变化
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.scatterplot(data=mg_all.obs, x='dpt_pseudotime', y='Pirb_expr', hue='time_point',
                    hue_order=order_time, alpha=0.3, s=10, ax=ax)
    sns.regplot(data=mg_all.obs, x='dpt_pseudotime', y='Pirb_expr', scatter=False,
                color='black', ax=ax)
    ax.set_title("Pirb expression along microglia activation pseudotime")
    plt.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "pirb_vs_pseudotime.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

# -----------------------------------------------------------------------------
# 7. 差异基因：Pirb+ vs Pirb- 小胶质细胞
# -----------------------------------------------------------------------------
print("[INFO] DE analysis: Pirb+ vs Pirb- microglia...")
mg_all.obs['Pirb_positive_str'] = mg_all.obs['Pirb_positive'].astype(str).astype('category')
sc.tl.rank_genes_groups(mg_all, groupby='Pirb_positive_str', method='wilcoxon', n_genes=100)

# 保存结果
result = mg_all.uns['rank_genes_groups']
groups = result['names'].dtype.names
de_list = []
for group in groups:
    n_genes = len(result['names'][group])
    for i in range(n_genes):
        de_list.append({
            'group': group,
            'gene': result['names'][group][i],
            'logfoldchange': result['logfoldchanges'][group][i],
            'pval': result['pvals'][group][i],
            'pval_adj': result['pvals_adj'][group][i],
            'score': result['scores'][group][i],
        })
de_df = pd.DataFrame(de_list)
de_df.to_csv(os.path.join(OUT_DIR, "pirb_pos_vs_neg_microglia_de.csv"), index=False)

# 热图：top DE genes（手动绘制，避免 ax 参数冲突）
n_top = 15
top_genes = []
for group in groups:
    genes = result['names'][group][:n_top]
    top_genes.extend(genes)
top_genes = list(dict.fromkeys(top_genes))  # 去重保持顺序

# 计算每个 Pirb 组的平均表达
expr_mat = []
for group in sorted(groups, key=lambda x: int(x)):
    mask = mg_all.obs['Pirb_positive_str'] == group
    expr = np.array(mg_all[mask, top_genes].X.mean(axis=0)).flatten()
    expr_mat.append(expr)
expr_mat = np.array(expr_mat)

fig, ax = plt.subplots(figsize=(12, 4))
im = ax.imshow(expr_mat, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
ax.set_xticks(range(len(top_genes)))
ax.set_xticklabels(top_genes, rotation=90, fontsize=8)
ax.set_yticks(range(len(groups)))
ax.set_yticklabels([f"Pirb+" if g == '1' else f"Pirb-" for g in sorted(groups, key=lambda x: int(x))])
ax.set_title("Top DE genes: Pirb+ vs Pirb- microglia")
plt.colorbar(im, ax=ax, label='Mean scaled expression')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_pos_vs_neg_heatmap.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# -----------------------------------------------------------------------------
# 8. 保存整合后的 h5ad
# -----------------------------------------------------------------------------
mg_all.write_h5ad(os.path.join(OUT_DIR, "microglia_cross_time_integrated.h5ad"))
print(f"[SAVE] {OUT_DIR}/microglia_cross_time_integrated.h5ad")

print("[DONE]")
