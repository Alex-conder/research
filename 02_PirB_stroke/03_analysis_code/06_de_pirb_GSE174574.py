"""
GSE174574 Pirb+ vs Pirb- 差异表达分析（按细胞类型）
- 对 Astrocyte / Microglia / Immune 分别比较 Pirb+ vs Pirb-
- 使用 Wilcoxon rank-sum
- 输出差异基因 CSV 和火山图
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_H5AD = "../04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "../04_reports/figures/GSE174574/de_pirb"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# 确保 Pirb 表达列存在
pirb = "Pirb"
adata.obs["Pirb_detected"] = (adata[:, pirb].X.toarray().flatten() > 0).astype(int)

# 关注的细胞类型
cell_types = ["Astrocyte", "Microglia", "Immune"]

for ct in cell_types:
    ad_ct = adata[adata.obs["cell_type"] == ct].copy()
    n_pos = (ad_ct.obs["Pirb_detected"] == 1).sum()
    n_neg = (ad_ct.obs["Pirb_detected"] == 0).sum()
    print(f"\n[DE] {ct}: Pirb+ = {n_pos}, Pirb- = {n_neg}")
    if n_pos < 10 or n_neg < 10:
        print(f"[WARN] {ct}: too few cells, skipping")
        continue
    
    ad_ct.obs["Pirb_group"] = ad_ct.obs["Pirb_detected"].astype("category")
    sc.tl.rank_genes_groups(ad_ct, groupby="Pirb_group", method="wilcoxon", pts=True)
    result = ad_ct.uns["rank_genes_groups"]
    groups = result['names'].dtype.names
    
    # 提取 group 1 (Pirb+) vs 0 (Pirb-)
    group_key = '1'
    df = pd.DataFrame({
        "gene": result["names"][group_key],
        "logfoldchange": result["logfoldchanges"][group_key],
        "pvals": result["pvals"][group_key],
        "pvals_adj": result["pvals_adj"][group_key],
        "scores": result["scores"][group_key],
    })
    df.to_csv(os.path.join(OUT_DIR, f"DE_PirbPos_vs_Neg_{ct}.csv"), index=False)
    
    # 火山图
    df["-log10_padj"] = -np.log10(df["pvals_adj"].replace(0, 1e-300))
    df["significant"] = (df["pvals_adj"] < 0.05) & (abs(df["logfoldchange"]) > 0.5)
    
    plt.figure(figsize=(8, 6))
    plt.scatter(df["logfoldchange"], df["-log10_padj"], c="gray", s=5, alpha=0.5)
    sig = df[df["significant"]]
    plt.scatter(sig["logfoldchange"], sig["-log10_padj"], c="red", s=10, alpha=0.7)
    plt.xlabel("Log2 fold change (Pirb+ vs Pirb-)")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"{ct}: Pirb+ vs Pirb-")
    plt.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=0.8)
    plt.axvline(0, color="black", linestyle="-", linewidth=0.8)
    # 标注 top 10 上调/下调基因
    top_up = sig.nlargest(10, "logfoldchange")
    top_down = sig.nsmallest(10, "logfoldchange")
    for _, row in pd.concat([top_up, top_down]).iterrows():
        plt.text(row["logfoldchange"], row["-log10_padj"], row["gene"], fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"volcano_{ct}.png"), dpi=150)
    plt.close()
    
    print(f"  Significant DEGs: {df['significant'].sum()}")

print(f"\n[DONE] Results saved to {OUT_DIR}")
