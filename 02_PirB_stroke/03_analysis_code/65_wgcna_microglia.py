"""
第五优先级：小胶质细胞跨时间整合矩阵的共表达网络分析
使用 K-means 对基因相关性向量聚类，识别与 Pirb 相关的模块。
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import pearsonr

OUT_DIR = "../04_reports/figures/wgcna_microglia"
os.makedirs(OUT_DIR, exist_ok=True)


def module_eigengene(X, module_genes_idx):
    """计算模块特征基因（第一主成分）。"""
    sub = X[:, module_genes_idx]
    sub = (sub - sub.mean(axis=0)) / (sub.std(axis=0) + 1e-10)
    u, s, vt = np.linalg.svd(sub, full_matrices=False)
    me = u[:, 0] * s[0]
    return me


def main():
    ad = sc.read_h5ad('04_reports/figures/microglia_cross_time/microglia_cross_time_integrated.h5ad')
    ad.obs_names_make_unique()
    print(f"[INFO] Loaded integrated microglia: {ad.shape}")

    genes = ad.var_names.tolist()
    n_genes = len(genes)
    print(f"[INFO] Using {n_genes} highly variable genes")

    X = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.asarray(ad.X)
    if X.max() > 100:
        X = np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4)

    # 计算相关性矩阵
    print("[INFO] Computing gene-gene correlation matrix...")
    corr = np.corrcoef(X.T)
    np.fill_diagonal(corr, 1)
    corr_df = pd.DataFrame(corr, index=genes, columns=genes)
    corr_df.to_csv(os.path.join(OUT_DIR, "gene_gene_correlation_matrix.csv"))

    # 选择最佳 K
    print("[INFO] Selecting optimal K for gene correlation K-means...")
    best_k = 6
    best_silhouette = -1
    for k in range(4, 13):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(corr)
        sil = silhouette_score(corr, labels)
        print(f"  K={k}: silhouette={sil:.4f}")
        if sil > best_silhouette:
            best_silhouette = sil
            best_k = k
    print(f"[INFO] Selected K={best_k} (silhouette={best_silhouette:.4f})")

    kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(corr)

    gene_module_df = pd.DataFrame({"gene": genes, "module": labels})
    gene_module_df.to_csv(os.path.join(OUT_DIR, "gene_module_assignment.csv"), index=False)
    print("\nModule sizes:")
    print(gene_module_df["module"].value_counts().sort_index())

    # Pirb 所在模块
    pirb_module = gene_module_df.loc[gene_module_df["gene"] == "Pirb", "module"].values[0]
    print(f"\n[Pirb] located in module {pirb_module}")

    # 模块特征基因与 Pirb 相关性
    pirb_expr = ad[:, "Pirb"].X.toarray().flatten() if hasattr(ad[:, "Pirb"].X, 'toarray') else np.asarray(ad[:, "Pirb"].X).flatten()
    module_records = []
    for mod in sorted(gene_module_df["module"].unique()):
        mod_genes = gene_module_df[gene_module_df["module"] == mod]["gene"].tolist()
        mod_idx = [genes.index(g) for g in mod_genes]
        me = module_eigengene(X, mod_idx)
        r, p = pearsonr(me, pirb_expr)
        module_records.append({
            "module": mod,
            "n_genes": len(mod_genes),
            "corr_with_Pirb": r,
            "p_value": p,
            "is_pirb_module": (mod == pirb_module),
        })
    module_corr_df = pd.DataFrame(module_records).sort_values("corr_with_Pirb", ascending=False)
    module_corr_df.to_csv(os.path.join(OUT_DIR, "module_correlation_with_Pirb.csv"), index=False)
    print("\n=== Module correlation with Pirb ===")
    print(module_corr_df.to_string())

    # 保存 Pirb 模块基因
    pirb_module_genes = gene_module_df[gene_module_df["module"] == pirb_module]["gene"].tolist()
    pd.DataFrame({"gene": pirb_module_genes}).to_csv(
        os.path.join(OUT_DIR, f"Pirb_module_{pirb_module}_genes.csv"), index=False
    )

    # Pirb 与所有基因的相关性
    pirb_corr = corr_df["Pirb"].drop("Pirb").sort_values(ascending=False)
    pirb_corr.to_csv(os.path.join(OUT_DIR, "gene_correlation_with_Pirb.csv"))
    print(f"\nPirb module {pirb_module}: {len(pirb_module_genes)} genes")
    print("Top 20 genes by correlation with Pirb:")
    print(pirb_corr.head(20).to_string())

    # 保存模块内 top 相关基因
    with open(os.path.join(OUT_DIR, "module_top_genes.txt"), "w") as f:
        for mod in sorted(gene_module_df["module"].unique()):
            mod_genes = gene_module_df[gene_module_df["module"] == mod]["gene"].tolist()
            mod_corr = pirb_corr[pirb_corr.index.isin(mod_genes)].head(10)
            f.write(f"\nModule {mod} ({len(mod_genes)} genes) - top 10 by Pirb correlation:\n")
            for g, r in mod_corr.items():
                f.write(f"  {g}\t{r:.4f}\n")

    # 可视化
    # 1. 模块-Pirb 相关性
    plt.figure(figsize=(10, 6))
    colors = ["red" if x else "gray" for x in module_corr_df["is_pirb_module"]]
    sns.barplot(data=module_corr_df, x="module", y="corr_with_Pirb", hue="module", palette=colors, legend=False)
    plt.title(f"Module eigengene correlation with Pirb expression (K={best_k})")
    plt.ylabel("Pearson r")
    plt.xlabel("Module")
    plt.axhline(0, color='black', linestyle='--')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "module_corr_with_Pirb.png"), dpi=300)
    plt.close()

    # 2. Pirb 相关性分布
    plt.figure(figsize=(8, 5))
    sns.histplot(pirb_corr.values, bins=40, kde=True)
    plt.axvline(0, color='black', linestyle='--')
    plt.title(f"Distribution of gene-Pirb correlations (n={len(pirb_corr)})")
    plt.xlabel("Pearson r with Pirb")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_correlation_distribution.png"), dpi=300)
    plt.close()

    # 3. 模块大小
    plt.figure(figsize=(10, 5))
    sizes = gene_module_df["module"].value_counts().sort_index()
    sns.barplot(x=sizes.index, y=sizes.values)
    plt.title(f"Module sizes (K={best_k})")
    plt.xlabel("Module")
    plt.ylabel("Number of genes")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "module_sizes.png"), dpi=300)
    plt.close()

    # 4. 相关性矩阵按模块排序的热图（仅 top 500 基因）
    top_genes = pirb_corr.head(500).index.tolist() + ["Pirb"]
    top_genes = list(dict.fromkeys(top_genes))  # 去重
    if "Pirb" not in top_genes:
        top_genes.append("Pirb")
    top_genes = [g for g in top_genes if g in genes]
    sub_corr = corr_df.loc[top_genes, top_genes]
    plt.figure(figsize=(12, 10))
    sns.heatmap(sub_corr, cmap="RdBu_r", center=0, vmin=-0.3, vmax=0.3, xticklabels=False, yticklabels=False)
    plt.title("Gene-gene correlation heatmap (top 500 Pirb-correlated genes)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "top_gene_correlation_heatmap.png"), dpi=200)
    plt.close()

    print(f"\n[SAVED] all results to {OUT_DIR}")


if __name__ == "__main__":
    main()
