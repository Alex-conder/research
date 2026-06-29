"""
第二优先级：GSE233814 Visium 空间梯度差异表达
1. 计算每个 spot 到组织边界的距离
2. 按距离分为梗死核心 / 半暗带 / 远端
3. 在半暗带内比较 Pirb+ vs Pirb- spots 的差异基因
4. 比较三个区域的整体炎症/趋化因子/基质金属蛋白酶模块评分
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy import stats

DATA_DIR = "../04_reports/figures/GSE233814"
OUT_DIR = os.path.join(DATA_DIR, "gradient_de")
os.makedirs(OUT_DIR, exist_ok=True)

# 模块基因集
MODULES = {
    "Chemokines": ["Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Cxcl10", "Cxcl2", "Cxcl9"],
    "MMPs": ["Mmp3", "Mmp9", "Mmp2", "Mmp14", "Mmp13", "Timp1"],
    "Inflammation": ["Il1b", "Tnf", "Il6", "Nfkbia", "Tnfaip3", "Socs3", "Cd68", "Itgam", "Aif1"],
    "Microglia_homeostasis": ["P2ry12", "Tmem119", "Cx3cr1", "Tgfbr1", "Selplg"],
}


def estimate_spot_spacing(xy):
    """通过计算最近邻距离的中位数估计 spot 间距（≈100μm）。"""
    if len(xy) < 2:
        return np.nan
    d = cdist(xy, xy)
    np.fill_diagonal(d, np.inf)
    nn = d.min(axis=1)
    return np.median(nn)


def compute_boundary_distance(spatial_df):
    """计算每个 spot 到组织凸包边界的距离。"""
    spatial_df = spatial_df.copy()
    spatial_df["dist_to_boundary"] = np.nan
    for sample in spatial_df["sample"].unique():
        sub = spatial_df[spatial_df["sample"] == sample]
        xy = sub[["imageX", "imageY"]].values
        if len(xy) < 4:
            continue
        hull = ConvexHull(xy)
        hull_points = xy[hull.vertices]
        dist = cdist(xy, hull_points).min(axis=1)
        spatial_df.loc[sub.index, "dist_to_boundary"] = dist
    return spatial_df


def classify_zones(spatial_df):
    """按估计 spot 间距将 spots 分为核心 / 半暗带 / 远端。"""
    spatial_df = spatial_df.copy()
    spatial_df["zone"] = "Unknown"
    for sample in spatial_df["sample"].unique():
        sub = spatial_df[spatial_df["sample"] == sample]
        xy = sub[["imageX", "imageY"]].values
        spacing = estimate_spot_spacing(xy)
        dist = sub["dist_to_boundary"].values
        # 核心：< 1 个间距；半暗带：1-3 个间距；远端：> 3 个间距
        core = dist < spacing
        penumbra = (dist >= spacing) & (dist < 3 * spacing)
        remote = dist >= 3 * spacing
        spatial_df.loc[sub.index[core], "zone"] = "Core"
        spatial_df.loc[sub.index[penumbra], "zone"] = "Penumbra"
        spatial_df.loc[sub.index[remote], "zone"] = "Remote"
        print(f"  {sample}: spacing={spacing:.1f} px; Core={core.sum()}, Penumbra={penumbra.sum()}, Remote={remote.sum()}")
    return spatial_df


def score_modules(ad, modules):
    """计算模块评分。"""
    for name, genes in modules.items():
        avail = [g for g in genes if g in ad.var_names]
        missing = [g for g in genes if g not in ad.var_names]
        if len(avail) >= 3:
            sc.tl.score_genes(ad, gene_list=avail, score_name=f"score_{name}", use_raw=False, random_state=42)
            print(f"  {name}: {len(avail)}/{len(genes)} available, missing: {missing}")
        else:
            ad.obs[f"score_{name}"] = np.nan


def spatial_de_in_penumbra(ad, spatial_df):
    """在半暗带内比较 Pirb+ vs Pirb- spots 的差异基因。"""
    # 确保 adata obs 与 spatial_df 对齐
    ad = ad[spatial_df["barcode"].values].copy()
    ad.obs["zone"] = spatial_df["zone"].values
    ad.obs["dist_to_boundary"] = spatial_df["dist_to_boundary"].values
    ad.obs["Pirb_positive"] = spatial_df["Pirb_positive"].values.astype(int)
    ad.obs["Pirb_group"] = ad.obs["Pirb_positive"].map({1: "Pirb_pos", 0: "Pirb_neg"})

    penumbra = ad[ad.obs["zone"] == "Penumbra"].copy()
    print(f"\nPenumbra spots: {penumbra.shape[0]} (Pirb+ {(penumbra.obs['Pirb_positive']==1).sum()})")

    if (penumbra.obs["Pirb_positive"] == 1).sum() < 10:
        print("  Too few Pirb+ spots in penumbra for DE.")
        return None, penumbra

    sc.tl.rank_genes_groups(
        penumbra,
        groupby="Pirb_group",
        groups=["Pirb_pos"],
        reference="Pirb_neg",
        method="wilcoxon",
        n_genes=penumbra.shape[1],
        use_raw=False,
    )
    res = penumbra.uns["rank_genes_groups"]
    df = pd.DataFrame({
        "gene": res["names"]["Pirb_pos"],
        "log2FC": res["logfoldchanges"]["Pirb_pos"],
        "pval": res["pvals"]["Pirb_pos"],
        "pval_adj": res["pvals_adj"]["Pirb_pos"],
        "score": res["scores"]["Pirb_pos"],
    })
    df = df.sort_values(["pval_adj", "log2FC"], ascending=[True, False])
    return df, penumbra


def compare_zones(ad, spatial_df):
    """比较三个区域的模块评分。"""
    ad = ad[spatial_df["barcode"].values].copy()
    ad.obs["zone"] = spatial_df["zone"].values
    score_cols = [c for c in ad.obs.columns if c.startswith("score_")]
    records = []
    for zone in ["Core", "Penumbra", "Remote"]:
        sub = ad.obs[ad.obs["zone"] == zone]
        if len(sub) == 0:
            continue
        for col in score_cols:
            records.append({
                "zone": zone,
                "metric": col,
                "mean": sub[col].mean(),
                "median": sub[col].median(),
                "n": len(sub),
            })
    return pd.DataFrame(records)


def plot_zone_scores(ad, spatial_df, outdir):
    """绘制三个区域模块评分箱线图。"""
    ad = ad[spatial_df["barcode"].values].copy()
    ad.obs["zone"] = spatial_df["zone"].values
    score_cols = [c for c in ad.obs.columns if c.startswith("score_")]
    plot_data = []
    for col in score_cols:
        for zone in ["Core", "Penumbra", "Remote"]:
            vals = ad.obs.loc[ad.obs["zone"] == zone, col].dropna()
            for v in vals:
                plot_data.append({"metric": col.replace("score_", ""), "zone": zone, "value": v})
    df_plot = pd.DataFrame(plot_data)
    if df_plot.empty:
        return
    plt.figure(figsize=(12, 6))
    order = ["Core", "Penumbra", "Remote"]
    sns.boxplot(data=df_plot, x="metric", y="value", hue="zone", hue_order=order, palette={"Core": "#d62728", "Penumbra": "#ff7f0e", "Remote": "#2ca02c"})
    plt.xticks(rotation=45, ha='right')
    plt.title("Module scores across spatial zones (GSE233814 D3)")
    plt.ylabel("Module score")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "zone_module_scores_D3.png"), dpi=300)
    plt.close()


def main():
    # 读取数据
    ad = sc.read_h5ad(os.path.join(DATA_DIR, "../GSE233814_processed.h5ad"))
    ad.obs_names_make_unique()
    spatial_df = pd.read_csv(os.path.join(DATA_DIR, "spot_pixel_coordinates_pirb.csv"))
    print(f"[INFO] Loaded {len(spatial_df)} spots, adata {ad.shape}")

    # 仅分析 D3
    spatial_df = spatial_df[spatial_df["time_point"] == "D3"].copy()
    print(f"[INFO] D3 spots: {len(spatial_df)}")

    # 计算边界距离和分区
    spatial_df = compute_boundary_distance(spatial_df)
    spatial_df = classify_zones(spatial_df)

    # 对齐 barcode
    spatial_df = spatial_df.set_index("barcode")
    common = list(set(spatial_df.index) & set(ad.obs_names))
    spatial_df = spatial_df.loc[common].reset_index()
    ad = ad[common].copy()
    print(f"[INFO] Aligned {len(common)} spots")

    # 计算模块评分
    score_modules(ad, MODULES)

    # 1. 空间梯度差异表达
    de_df, penumbra_ad = spatial_de_in_penumbra(ad, spatial_df)
    if de_df is not None:
        de_df.to_csv(os.path.join(OUT_DIR, "penumbra_PirbPos_vs_Neg_DE.csv"), index=False)
        print("\n=== Penumbra Pirb+ vs Pirb- top DE genes ===")
        print(de_df.head(20).to_string())

        # 筛选趋化因子/MMP/炎症基因
        markers = []
        for mod_genes in MODULES.values():
            markers.extend(mod_genes)
        marker_de = de_df[de_df["gene"].isin(markers)]
        marker_de.to_csv(os.path.join(OUT_DIR, "penumbra_PirbPos_vs_Neg_module_genes.csv"), index=False)
        print("\n=== Penumbra Pirb+ vs Pirb- module gene DE ===")
        print(marker_de.to_string())

    # 2. 三个区域模块评分比较
    zone_scores = compare_zones(ad, spatial_df)
    zone_scores.to_csv(os.path.join(OUT_DIR, "zone_module_score_summary.csv"), index=False)
    print("\n=== Zone module score summary ===")
    print(zone_scores.to_string())

    # 3. 可视化
    plot_zone_scores(ad, spatial_df, OUT_DIR)

    # 4. 保存空间分区信息
    spatial_df.to_csv(os.path.join(OUT_DIR, "D3_spatial_zones.csv"), index=False)

    print(f"\n[SAVED] all results to {OUT_DIR}")


if __name__ == "__main__":
    main()
