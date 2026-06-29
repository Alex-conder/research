"""
第一优先级：解析 D3 达峰与 D7 回落的转录调控刹车
1. 计算应激/凋亡/焦亡/氧化应激模块评分
2. 比较 Pirb+ vs Pirb- 小胶质细胞的 percent.mt 和模块评分
3. 追踪 D3 vs D7 Pirb+ 小胶质细胞中抑制性 TF/负调控因子的表达
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

OUT_DIR = "../04_reports/figures/regulatory_brake"
os.makedirs(OUT_DIR, exist_ok=True)

# 基因集
GENE_SETS = {
    "Apoptosis": ["Casp3", "Casp7", "Casp8", "Casp9", "Fas", "Faslg", "Bax", "Bak1", "Bcl2", "Bcl2l1", "Parp1", "Casp6", "Casp2"],
    "Pyroptosis": ["Nlrp3", "Casp1", "Casp4", "Casp5", "Il1b", "Il18", "Gsdmd", "Pycard", "Naip2", "Nlrc4"],
    "Oxidative_stress": ["Hif1a", "Nfe2l2", "Sod1", "Sod2", "Cat", "Gpx1", "Gpx4", "Hmox1", "Nqo1", "Prdx1"],
    "Inhibitory_TFs": ["Socs1", "Socs3", "Tsc22d3", "Nfkbia", "Tnfaip3", "Zfp36", "Dusp1", "Irf2bp2", "Bcl6", "Klf2"],
}

STRESS_GENES = [g for genes in [GENE_SETS["Apoptosis"], GENE_SETS["Pyroptosis"], GENE_SETS["Oxidative_stress"]] for g in genes]

def load_and_subset(path, celltype_col="cell_type", microglia_label="Microglia", condition_col=None, condition_val=None):
    ad = sc.read_h5ad(path)
    ad.obs_names_make_unique()
    mask = ad.obs[celltype_col] == microglia_label
    if condition_col and condition_val:
        mask = mask & (ad.obs[condition_col] == condition_val)
    sub = ad[mask].copy()
    print(f"Loaded {path}: microglia {sub.shape[0]}")
    return sub


def score_genes_safe(ad, gene_set_name, gene_list):
    """计算模块评分，过滤缺失基因。"""
    avail = [g for g in gene_list if g in ad.var_names]
    missing = [g for g in gene_list if g not in ad.var_names]
    print(f"  {gene_set_name}: {len(avail)}/{len(gene_list)} genes available, missing: {missing}")
    if len(avail) >= 3:
        sc.tl.score_genes(ad, gene_list=avail, score_name=f"score_{gene_set_name}", use_raw=False, random_state=42)
    else:
        ad.obs[f"score_{gene_set_name}"] = np.nan


def compare_pirb_pos_neg(ad, dataset_name, time_label):
    """比较 Pirb+ vs Pirb- 小胶质细胞的评分，返回 ad（含 Pirb_group）和 df。"""
    if "Pirb" not in ad.var_names:
        print(f"  Warning: Pirb not in {dataset_name} var_names")
        return ad, None
    ad.obs["Pirb_positive"] = (ad[:, "Pirb"].X.toarray().flatten() > 0).astype(int)
    ad.obs["Pirb_group"] = ad.obs["Pirb_positive"].map({1: "Pirb_pos", 0: "Pirb_neg"})

    records = []
    score_cols = [c for c in ad.obs.columns if c.startswith("score_")] + ["pct_counts_mt"]
    for col in score_cols:
        pos = ad.obs.loc[ad.obs["Pirb_group"] == "Pirb_pos", col].dropna()
        neg = ad.obs.loc[ad.obs["Pirb_group"] == "Pirb_neg", col].dropna()
        if len(pos) > 0 and len(neg) > 0:
            stat, pval = stats.mannwhitneyu(pos, neg, alternative='two-sided')
            records.append({
                "dataset": dataset_name,
                "time": time_label,
                "metric": col,
                "Pirb_pos_mean": pos.mean(),
                "Pirb_pos_median": pos.median(),
                "Pirb_pos_n": len(pos),
                "Pirb_neg_mean": neg.mean(),
                "Pirb_neg_median": neg.median(),
                "Pirb_neg_n": len(neg),
                "p_value": pval,
                "log2FC_mean": np.log2((pos.mean() + 1e-10) / (neg.mean() + 1e-10)),
            })
    df = pd.DataFrame(records)
    return ad, df


def d3_vs_d7_inhibitory(ad233812):
    """分析 GSE233812 D3 vs D7 Pirb+ 小胶质细胞中抑制性 TF 和应激基因表达。"""
    ad = ad233812.copy()
    ad.obs["Pirb_positive"] = (ad[:, "Pirb"].X.toarray().flatten() > 0).astype(int)
    mask = ad.obs["Pirb_positive"] == 1
    sub = ad[mask].copy()
    groups = sub.obs.groupby("time_point")
    d3_cells = sub[sub.obs["time_point"] == "D3"]
    d7_cells = sub[sub.obs["time_point"] == "D7"]

    genes = GENE_SETS["Inhibitory_TFs"] + STRESS_GENES
    genes = [g for g in genes if g in sub.var_names]

    results = []
    for g in genes:
        d3_expr = d3_cells[:, g].X.toarray().flatten()
        d7_expr = d7_cells[:, g].X.toarray().flatten()
        if d3_expr.std() == 0 and d7_expr.std() == 0:
            continue
        stat, pval = stats.mannwhitneyu(d7_expr, d3_expr, alternative='two-sided')
        results.append({
            "gene": g,
            "D3_mean": d3_expr.mean(),
            "D7_mean": d7_expr.mean(),
            "log2FC_D7_vs_D3": np.log2((d7_expr.mean() + 1e-10) / (d3_expr.mean() + 1e-10)),
            "p_value": pval,
            "D3_n": len(d3_expr),
            "D7_n": len(d7_expr),
        })
    df = pd.DataFrame(results).sort_values("p_value")
    return df


def plot_boxplots(ad, dataset_name, time_label, outdir):
    """绘制 Pirb+ vs Pirb- 的评分箱线图。"""
    score_cols = [c for c in ad.obs.columns if c.startswith("score_")] + ["pct_counts_mt"]
    plot_data = []
    for col in score_cols:
        for group in ["Pirb_pos", "Pirb_neg"]:
            vals = ad.obs.loc[ad.obs["Pirb_group"] == group, col].dropna()
            for v in vals:
                plot_data.append({"metric": col.replace("score_", "").replace("pct_counts_mt", "percent.mt"), "group": group, "value": v})
    df_plot = pd.DataFrame(plot_data)
    if df_plot.empty:
        return
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=df_plot, x="metric", y="value", hue="group", palette={"Pirb_pos": "#d62728", "Pirb_neg": "#1f77b4"})
    plt.xticks(rotation=45, ha='right')
    plt.title(f"{dataset_name} {time_label}: Pirb+ vs Pirb- microglia")
    plt.ylabel("Score / Percentage")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"boxplot_{dataset_name}_{time_label}.png"), dpi=300)
    plt.close()


def main():
    # 1. GSE174574 MCAO 24h microglia
    ad174 = load_and_subset(
        '04_reports/figures/GSE174574/GSE174574_annotated.h5ad',
        condition_col="condition", condition_val="MCAO"
    )
    for name, genes in GENE_SETS.items():
        score_genes_safe(ad174, name, genes)
    ad174, df174 = compare_pirb_pos_neg(ad174, "GSE174574", "MCAO_24h")
    plot_boxplots(ad174, "GSE174574", "MCAO_24h", OUT_DIR)

    # 2. GSE233812 D3/D7 microglia
    ad812 = load_and_subset('04_reports/figures/GSE233812_processed.h5ad')
    for name, genes in GENE_SETS.items():
        score_genes_safe(ad812, name, genes)
    ad812, df812_all = compare_pirb_pos_neg(ad812, "GSE233812", "all")
    ad812_d3 = ad812[ad812.obs["time_point"] == "D3"].copy()
    ad812_d7 = ad812[ad812.obs["time_point"] == "D7"].copy()
    _, df812_d3 = compare_pirb_pos_neg(ad812_d3, "GSE233812", "D3")
    _, df812_d7 = compare_pirb_pos_neg(ad812_d7, "GSE233812", "D7")
    plot_boxplots(ad812_d3, "GSE233812", "D3", OUT_DIR)
    plot_boxplots(ad812_d7, "GSE233812", "D7", OUT_DIR)

    # 3. D3 vs D7 抑制性 TF
    df_inhib = d3_vs_d7_inhibitory(ad812)

    # 4. 保存结果
    df174.to_csv(os.path.join(OUT_DIR, "GSE174574_MCAO24h_pirb_pos_neg_comparison.csv"), index=False)
    df812_d3.to_csv(os.path.join(OUT_DIR, "GSE233812_D3_pirb_pos_neg_comparison.csv"), index=False)
    df812_d7.to_csv(os.path.join(OUT_DIR, "GSE233812_D7_pirb_pos_neg_comparison.csv"), index=False)
    df_inhib.to_csv(os.path.join(OUT_DIR, "GSE233812_D7_vs_D3_inhibitory_stress_TFs.csv"), index=False)

    # 5. 打印关键结论
    print("\n=== GSE174574 MCAO 24h Pirb+ vs Pirb- ===")
    print(df174.to_string())
    print("\n=== GSE233812 D3 Pirb+ vs Pirb- ===")
    print(df812_d3.to_string())
    print("\n=== GSE233812 D7 Pirb+ vs Pirb- ===")
    print(df812_d7.to_string())
    print("\n=== GSE233812 D7 vs D3 Pirb+ cells: top inhibitory/stress genes ===")
    print(df_inhib.head(20).to_string())

    # 6. 抑制性 TF 显著上调的判断
    inhib_up = df_inhib[(df_inhib["gene"].isin(GENE_SETS["Inhibitory_TFs"])) & (df_inhib["log2FC_D7_vs_D3"] > 0) & (df_inhib["p_value"] < 0.05)]
    print(f"\n显著上调的抑制性 TF/负调控因子（D7 vs D3, p<0.05）: {len(inhib_up)} 个")
    print(inhib_up[["gene", "D3_mean", "D7_mean", "log2FC_D7_vs_D3", "p_value"]].to_string())

    # 7. percent.mt 判断
    mt_d3_pos = ad812[(ad812.obs["time_point"] == "D3") & (ad812.obs["Pirb_group"] == "Pirb_pos")].obs["pct_counts_mt"]
    mt_d7_pos = ad812[(ad812.obs["time_point"] == "D7") & (ad812.obs["Pirb_group"] == "Pirb_pos")].obs["pct_counts_mt"]
    mt_d3_neg = ad812[(ad812.obs["time_point"] == "D3") & (ad812.obs["Pirb_group"] == "Pirb_neg")].obs["pct_counts_mt"]
    print(f"\npercent.mt (Pirb+ microglia):")
    print(f"  D3 Pirb+: mean={mt_d3_pos.mean():.2f}, median={mt_d3_pos.median():.2f}")
    print(f"  D7 Pirb+: mean={mt_d7_pos.mean():.2f}, median={mt_d7_pos.median():.2f}")
    print(f"  D3 Pirb-: mean={mt_d3_neg.mean():.2f}, median={mt_d3_neg.median():.2f}")
    stat, pval = stats.mannwhitneyu(mt_d3_pos, mt_d3_neg, alternative='two-sided')
    print(f"  D3 Pirb+ vs Pirb- percent.mt p-value: {pval:.4g}")

    print(f"\n[SAVED] all results to {OUT_DIR}")


if __name__ == "__main__":
    main()
