"""
GSE174574 各细胞类型 Pirb+ vs Pirb- 差异表达分析
- Astrocyte, Microglia, Immune, Endothelial, Oligodendrocyte
- 输出 CSV 和火山图摘要
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb_all_ct"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
adata.obs["Pirb_detected"] = (adata[:, "Pirb"].X.toarray().flatten() > 0).astype(int)
adata.obs["Pirb_group"] = adata.obs["Pirb_detected"].astype("category")

# 关注的细胞类型
cell_types = ["Astrocyte", "Microglia", "Immune", "Endothelial", "Oligodendrocyte", "Pericyte"]

summary = []
for ct in cell_types:
    ad_ct = adata[adata.obs["cell_type"] == ct].copy()
    n_pos = (ad_ct.obs["Pirb_detected"] == 1).sum()
    n_neg = (ad_ct.obs["Pirb_detected"] == 0).sum()
    print(f"\n[DE] {ct}: Pirb+ = {n_pos}, Pirb- = {n_neg}")
    if n_pos < 10 or n_neg < 10:
        print(f"[WARN] {ct}: too few cells, skipping")
        continue
    
    sc.tl.rank_genes_groups(ad_ct, groupby="Pirb_group", method="wilcoxon", pts=True)
    result = ad_ct.uns["rank_genes_groups"]
    group_key = '1'
    df = pd.DataFrame({
        "gene": result["names"][group_key],
        "logfoldchange": result["logfoldchanges"][group_key],
        "pvals": result["pvals"][group_key],
        "pvals_adj": result["pvals_adj"][group_key],
        "scores": result["scores"][group_key],
    })
    df.to_csv(os.path.join(OUT_DIR, f"DE_PirbPos_vs_Neg_{ct}.csv"), index=False)
    
    sig_up = df[(df["pvals_adj"] < 0.05) & (df["logfoldchange"] > 0.5)]
    sig_down = df[(df["pvals_adj"] < 0.05) & (df["logfoldchange"] < -0.5)]
    summary.append({
        "cell_type": ct,
        "n_pirb_pos": int(n_pos),
        "n_pirb_neg": int(n_neg),
        "n_sig_up": len(sig_up),
        "n_sig_down": len(sig_down),
        "top_up_genes": ", ".join(sig_up.head(10)["gene"].tolist()),
        "top_down_genes": ", ".join(sig_down.head(10)["gene"].tolist()),
    })
    
    # 火山图
    df["-log10_padj"] = -np.log10(df["pvals_adj"].replace(0, 1e-300))
    df["significant"] = (df["pvals_adj"] < 0.05) & (abs(df["logfoldchange"]) > 0.5)
    plt.figure(figsize=(8, 6))
    plt.scatter(df["logfoldchange"], df["-log10_padj"], c="gray", s=5, alpha=0.5)
    sig = df[df["significant"]]
    plt.scatter(sig["logfoldchange"], sig["-log10_padj"], c="red", s=10, alpha=0.7)
    top = pd.concat([sig.nlargest(8, "logfoldchange"), sig.nsmallest(8, "logfoldchange")])
    for _, row in top.iterrows():
        plt.text(row["logfoldchange"], row["-log10_padj"], row["gene"], fontsize=7)
    plt.xlabel("Log2 fold change (Pirb+ vs Pirb-)")
    plt.ylabel("-log10 adjusted p-value")
    plt.title(f"{ct}: Pirb+ vs Pirb-")
    plt.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=0.8)
    plt.axvline(0, color="black", linestyle="-", linewidth=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"volcano_{ct}.png"), dpi=150)
    plt.close()
    
    print(f"  Sig up: {len(sig_up)}, Sig down: {len(sig_down)}")

pd.DataFrame(summary).to_csv(os.path.join(OUT_DIR, "de_summary.csv"), index=False)
print("\nDE summary:")
print(pd.DataFrame(summary))
print(f"[DONE] Results saved to {OUT_DIR}")
