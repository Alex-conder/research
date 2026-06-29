"""
GSE174574 Pirb+ vs Pirb- 转录因子表达分析
- 手动检查关键免疫/炎症相关 TF 的表达差异
- 输出 CSV 和箱线图
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "../04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "../04_reports/figures/GSE174574/tf_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
adata.obs["Pirb_detected"] = (adata[:, "Pirb"].X.toarray().flatten() > 0).astype(int)

# 关键转录因子列表
tf_list = [
    # NF-kB 家族
    "Nfkb1", "Nfkb2", "Rela", "Relb", "Rel",
    # STAT 家族
    "Stat1", "Stat2", "Stat3", "Stat4", "Stat5a", "Stat5b", "Stat6",
    # IRF 家族
    "Irf1", "Irf2", "Irf3", "Irf4", "Irf5", "Irf7", "Irf8", "Irf9",
    # AP-1
    "Jun", "Fos", "Jund", "Fosb", "Junb",
    # 其他炎症/应激 TF
    "Nfkbia", "Ikbkb", "Mapk14", "Cebpb", "Cebpd", "Atf3", "Egr1",
    # A1 星形胶质细胞相关
    "Stat3", "Nfkb1", "Irf1",
]
tf_list = [g for g in tf_list if g in adata.var_names]
print(f"[INFO] Available TFs: {tf_list}")

# 按细胞类型和 Pirb 状态计算 TF 平均表达
records = []
for ct in ["Astrocyte", "Microglia", "Immune", "Endothelial", "Oligodendrocyte"]:
    ad_ct = adata[adata.obs["cell_type"] == ct]
    for tf in tf_list:
        expr = ad_ct[:, tf].X.toarray().flatten()
        df = pd.DataFrame({"expr": expr, "pirb": ad_ct.obs["Pirb_detected"].values})
        for pval in [0, 1]:
            sub = df[df["pirb"] == pval]
            records.append({
                "cell_type": ct,
                "tf": tf,
                "pirb_group": "Pirb-" if pval == 0 else "Pirb+",
                "n": len(sub),
                "mean_expr": sub["expr"].mean(),
                "median_expr": sub["expr"].median(),
                "pct_detected": (sub["expr"] > 0).mean() * 100
            })

tf_stats = pd.DataFrame(records)
tf_stats.to_csv(os.path.join(OUT_DIR, "tf_expression_by_celltype_pirb.csv"), index=False)

# 计算 Pirb+ vs Pirb- 的 fold change
fc = tf_stats.pivot_table(index=["cell_type", "tf"], columns="pirb_group", values="mean_expr").reset_index()
fc["log2FC"] = np.log2((fc["Pirb+"] + 1e-6) / (fc["Pirb-"] + 1e-6))
fc.to_csv(os.path.join(OUT_DIR, "tf_log2fc.csv"), index=False)

# 筛选显著上调的 TF（log2FC > 0.5）
sig_tf = fc[fc["log2FC"] > 0.5].sort_values("log2FC", ascending=False)
print("\nTop up-regulated TFs in Pirb+ cells:")
print(sig_tf.head(20))
sig_tf.to_csv(os.path.join(OUT_DIR, "tf_up_in_pirb_pos.csv"), index=False)

# 热图：TF log2FC across cell types
pivot = fc.pivot(index="tf", columns="cell_type", values="log2FC").fillna(0)
plt.figure(figsize=(10, 12))
sns.heatmap(pivot, cmap="RdBu_r", center=0, vmin=-2, vmax=2, linewidths=0.5)
plt.title("TF log2FC (Pirb+ vs Pirb-) across cell types")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "tf_heatmap_log2fc.png"), dpi=150)
plt.close()

# 箱线图：Astrocyte 中关键 TF 的 Pirb+ vs Pirb- 表达
astro = adata[adata.obs["cell_type"] == "Astrocyte"]
astro_df = pd.DataFrame({tf: astro[:, tf].X.toarray().flatten() for tf in tf_list})
astro_df["Pirb_group"] = astro.obs["Pirb_detected"].map({0: "Pirb-", 1: "Pirb+"}).values

# 选择变化最大的 top 12 TF
astro_fc = fc[fc["cell_type"] == "Astrocyte"].sort_values("log2FC", ascending=False).head(12)["tf"].tolist()
plot_df = astro_df.melt(id_vars="Pirb_group", value_vars=astro_fc, var_name="TF", value_name="Expression")
plt.figure(figsize=(14, 6))
sns.boxplot(data=plot_df, x="TF", y="Expression", hue="Pirb_group")
plt.xticks(rotation=45, ha="right")
plt.title("Top TF expression in Astrocyte Pirb+ vs Pirb-")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "tf_boxplot_astrocyte.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
