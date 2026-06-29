"""
GSE233814 空间转录组补充分析：Pirb 与炎症/小胶质 marker 的关系
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE233814_processed.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233814"

adata = sc.read_h5ad(IN_H5AD)

# 计算小胶质/炎症 marker 表达
markers = {
    'microglia': ['C1qa', 'C1qb', 'P2ry12', 'Tmem119', 'Cx3cr1'],
    'astrocyte_reactive': ['Gfap', 'C3', 'Vim'],
    'inflammation': ['Tnf', 'Il1a', 'Il1b', 'Spp1', 'Lyz2'],
    'macrophage': ['Adgre1', 'Cd68', 'Lyz2'],
}

for name, genes in markers.items():
    available = [g for g in genes if g in adata.var_names]
    if len(available) == 0:
        continue
    adata.obs[name] = np.array(adata[:, available].X.mean(axis=1)).flatten()

# Pirb+ vs Pirb- spots 的 marker 比较
adata.obs['Pirb_positive'] = (adata.obs['Pirb_expr'] > 0).astype(int)
compare = adata.obs.groupby('time_point').apply(
    lambda df: pd.Series({
        'pirb_pos_n': (df['Pirb_positive'] == 1).sum(),
        'pirb_neg_n': (df['Pirb_positive'] == 0).sum(),
        'pirb_pos_microglia': df.loc[df['Pirb_positive'] == 1, 'microglia'].mean(),
        'pirb_neg_microglia': df.loc[df['Pirb_positive'] == 0, 'microglia'].mean(),
        'pirb_pos_inflammation': df.loc[df['Pirb_positive'] == 1, 'inflammation'].mean(),
        'pirb_neg_inflammation': df.loc[df['Pirb_positive'] == 0, 'inflammation'].mean(),
        'pirb_pos_Gfap': df.loc[df['Pirb_positive'] == 1, 'Gfap'].mean() if 'Gfap' in df.columns else np.nan,
        'pirb_neg_Gfap': df.loc[df['Pirb_positive'] == 0, 'Gfap'].mean() if 'Gfap' in df.columns else np.nan,
    })
)
compare.to_csv(os.path.join(OUT_DIR, "pirb_pos_vs_neg_markers.csv"))
print(compare)

# 小提琴图：D3 vs control 的 marker 表达
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (name, label) in zip(axes, [
    ('microglia', 'Microglia score'),
    ('inflammation', 'Inflammation score'),
    ('Gfap', 'Gfap expression')
]):
    if name not in adata.obs.columns:
        ax.text(0.5, 0.5, f'{name} not available', ha='center')
        continue
    subset = adata[adata.obs['time_point'].isin(['control', 'D3'])]
    sns.violinplot(data=subset.obs, x='time_point', y=name, order=['control', 'D3'], ax=ax)
    ax.set_title(label)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "marker_expression_control_vs_D3.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

# Pirb 与 microglia score 的相关性
fig, ax = plt.subplots(figsize=(6, 5))
sns.scatterplot(data=adata.obs, x='microglia', y='Pirb_expr', hue='time_point', alpha=0.3, ax=ax)
ax.set_title('Pirb vs microglia score')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_vs_microglia_scatter.png"), dpi=200, bbox_inches='tight')
plt.close(fig)

print("[DONE]")
