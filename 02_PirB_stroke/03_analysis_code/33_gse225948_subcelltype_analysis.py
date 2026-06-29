"""
GSE225948 sub.celltype 水平 Pirb 分析
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE225948_processed.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE225948"

adata = sc.read_h5ad(IN_H5AD)

summary = adata.obs.groupby(['tissue', 'time_point', 'parent_celltype', 'sub_celltype'], observed=False).agg(
    n=('Pirb_positive', 'size'),
    pirb_frac=('Pirb_positive', 'mean'),
    pirb_mean=('Pirb_expr', 'mean'),
).reset_index()
summary = summary[summary['n'] >= 30]

# 保存完整表格
summary.to_csv(os.path.join(OUT_DIR, "pirb_summary_subcelltype.csv"), index=False)

# 每个 parent cell type 中取 top subcelltype
for tissue in ['PB', 'brain']:
    for tp in ['Sham', 'D02', 'D14']:
        sub = summary[(summary['tissue']==tissue) & (summary['time_point']==tp)].copy()
        if len(sub) == 0:
            continue
        sub = sub.sort_values('pirb_frac', ascending=False).head(15)
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(data=sub, x='pirb_frac', y='sub_celltype', hue='parent_celltype', dodge=False, ax=ax)
        ax.set_title(f'{tissue} {tp} - Top Pirb+ sub.celltypes')
        ax.set_xlabel('Pirb+ fraction')
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"pirb_top_subcelltypes_{tissue}_{tp}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

# 选择关键 parent cell type 的 sub.celltype 热图
key_parents = ['Neu', 'Mo', 'MdC', 'Gran', 'DC', 'Mg']
sub_key = summary[summary['parent_celltype'].isin(key_parents)].copy()
pivot = sub_key.pivot_table(index=['parent_celltype', 'sub_celltype'], columns=['tissue', 'time_point'], values='pirb_frac', fill_value=0)
fig, ax = plt.subplots(figsize=(10, 12))
sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Pirb+ fraction'})
ax.set_title('Pirb+ fraction by sub.celltype, tissue and time')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_subcelltype_heatmap.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

print("[DONE] Sub.celltype analysis saved.")
