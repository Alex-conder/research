"""
GSE225948 关键结果总结图
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE225948_processed.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE225948"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)

# 确保 Pirb_expr 存在
if 'Pirb_expr' not in adata.obs.columns:
    adata.obs['Pirb_expr'] = np.array(adata[:, 'Pirb'].X.toarray()).flatten()
    adata.obs['Pirb_positive'] = (adata.obs['Pirb_expr'] > 0).astype(int)

# 关键细胞类型
key_cts = ['Neu', 'Mo', 'MdC', 'Gran', 'DC', 'Mg', 'EC', 'Bc', 'Tc', 'NK']
subset = adata[adata.obs['parent_celltype'].isin(key_cts)].copy()

# 1. 各组织/时间下关键细胞类型的 Pirb+ fraction 点图
summary = subset.obs.groupby(['tissue', 'time_point', 'parent_celltype'], observed=False).agg(
    n=('Pirb_positive', 'size'),
    pirb_frac=('Pirb_positive', 'mean'),
    pirb_mean=('Pirb_expr', 'mean'),
).reset_index()
summary = summary[summary['n'] >= 20]
summary['tissue_time'] = summary['tissue'].astype(str) + '_' + summary['time_point'].astype(str)
summary['parent_celltype'] = pd.Categorical(summary['parent_celltype'], categories=key_cts, ordered=True)
summary = summary.sort_values(['tissue', 'parent_celltype', 'time_point'])

# 按 tissue 分开绘制
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
for ax, tissue in zip(axes, ['PB', 'brain']):
    df = summary[summary['tissue'] == tissue]
    sns.pointplot(data=df, x='parent_celltype', y='pirb_frac', hue='time_point',
                  order=key_cts, hue_order=['Sham', 'D02', 'D14'],
                  palette=['gray', 'red', 'blue'], dodge=0.3, ax=ax, markers='o', capsize=0.1)
    ax.set_ylabel('Pirb+ fraction')
    ax.set_title(f'{tissue} - Pirb+ fraction by cell type and time')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(title='Time')
    ax.set_ylim(0, 0.7)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_fraction_pointplot.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 2. 小提琴图：关键细胞类型在 brain 中的时间动态
brain_subset = subset[subset.obs['tissue'] == 'brain']
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()
for idx, ct in enumerate(['Neu', 'Mo', 'MdC', 'Gran', 'Mg', 'DC']):
    ct_sub = brain_subset[brain_subset.obs['parent_celltype'] == ct]
    if ct_sub.n_obs < 20:
        axes[idx].text(0.5, 0.5, f'{ct}: too few cells', ha='center', va='center')
        axes[idx].set_title(ct)
        continue
    sns.violinplot(data=ct_sub.obs, x='time_point', y='Pirb_expr', ax=axes[idx],
                   order=['Sham', 'D02', 'D14'])
    axes[idx].set_title(f'{ct} (n={ct_sub.n_obs})')
    axes[idx].set_ylabel('Pirb expression')
plt.suptitle('Pirb expression in brain immune cells across stroke time course', y=1.02, fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_brain_key_celltypes_violin.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 3. PB 中的时间动态
pb_subset = subset[subset.obs['tissue'] == 'PB']
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()
for idx, ct in enumerate(['Neu', 'Mo', 'DC', 'Bc', 'Tc', 'NK']):
    ct_sub = pb_subset[pb_subset.obs['parent_celltype'] == ct]
    if ct_sub.n_obs < 20:
        axes[idx].text(0.5, 0.5, f'{ct}: too few cells', ha='center', va='center')
        axes[idx].set_title(ct)
        continue
    sns.violinplot(data=ct_sub.obs, x='time_point', y='Pirb_expr', ax=axes[idx],
                   order=['Sham', 'D02', 'D14'])
    axes[idx].set_title(f'{ct} (n={ct_sub.n_obs})')
    axes[idx].set_ylabel('Pirb expression')
plt.suptitle('Pirb expression in peripheral blood cells across stroke time course', y=1.02, fontsize=14)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_PB_key_celltypes_violin.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 4. 保存总结表格
summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_key_celltypes.csv"), index=False)
print("[SAVE] pirb_summary_key_celltypes.csv")
print(summary)
print("[DONE]")
