"""
计算 GSE225948 外周血 (PB) 单核细胞 (Mo) 和中性粒细胞 (Neu)
在 D02 急性期的 Pirb+ vs Pirb- 差异基因，用于补充 Supplementary Table 8。
"""
import os
import scanpy as sc
import pandas as pd
import numpy as np

OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE225948"
OUT_CSV = os.path.join(OUT_DIR, "DE_PirbPos_vs_Neg_PB_Mo_Neu_D02.csv")

ad = sc.read_h5ad('04_reports/figures/GSE225948_processed.h5ad')

# 选取 PB, D02, parent_celltype 为 Mo 或 Neu 的细胞
mask = (
    (ad.obs['tissue'] == 'PB') &
    (ad.obs['time_point'] == 'D02') &
    (ad.obs['parent_celltype'].isin(['Mo', 'Neu']))
)
sub = ad[mask].copy()
print(f"Selected PB D02 Mo/Neu cells: {sub.shape[0]}")
print(sub.obs['parent_celltype'].value_counts())
print("Pirb+ fraction:", sub.obs['Pirb_positive'].mean())

# 确保 Pirb_positive 是字符串类别，用于 DE
sub.obs['Pirb_group'] = sub.obs['Pirb_positive'].astype(str).map({'1': 'Pirb_pos', '0': 'Pirb_neg'})

# rank_genes_groups
sc.tl.rank_genes_groups(sub, groupby='Pirb_group', groups=['Pirb_pos'], reference='Pirb_neg',
                        method='wilcoxon', n_genes=sub.shape[1], use_raw=False)

# 提取结果
res = sub.uns['rank_genes_groups']
# res['names'] 是 recarray
genes = res['names']['Pirb_pos']
df = pd.DataFrame({
    'gene': genes,
    'log2FC': res['logfoldchanges']['Pirb_pos'],
    'pval': res['pvals']['Pirb_pos'],
    'pval_adj': res['pvals_adj']['Pirb_pos'],
    'scores': res['scores']['Pirb_pos'],
})
# manually compute percentage expressed
pos_mask = sub.obs['Pirb_group'] == 'Pirb_pos'
neg_mask = sub.obs['Pirb_group'] == 'Pirb_neg'
X = sub.X.toarray() if hasattr(sub.X, 'toarray') else sub.X
expr = (X > 0).astype(float)
var_idx = {g: i for i, g in enumerate(sub.var_names)}
gene_idx = [var_idx[g] for g in genes]
df['pct_Pirb_pos'] = expr[pos_mask.values][:, gene_idx].mean(axis=0)
df['pct_Pirb_neg'] = expr[neg_mask.values][:, gene_idx].mean(axis=0)
df['cell_subset'] = 'PB_Mo_Neu_D02'
# 排序
df = df.sort_values(['pval_adj', 'log2FC'], ascending=[True, False]).reset_index(drop=True)

df.to_csv(OUT_CSV, index=False)
print(f"[SAVE] {OUT_CSV}")
print(df.head(10).to_string())
