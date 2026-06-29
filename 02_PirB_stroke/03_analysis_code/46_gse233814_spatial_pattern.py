"""
GSE233814 Visium 空间模式分析：
1. Pirb+ spot 的最近邻距离与空间聚类程度
2. Pirb+ spot 到组织边界的距离分布
3. 按空间分区的 Pirb 阳性率
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull, distance
from scipy.spatial.distance import cdist
from scipy.stats import mannwhitneyu

DATA_DIR = "../04_reports/figures/GSE233814"
OUT_DIR = DATA_DIR

# 读取空间坐标 + Pirb 表达
spatial_df = pd.read_csv(os.path.join(DATA_DIR, "spot_pixel_coordinates_pirb.csv"))
print(f"[INFO] Loaded {len(spatial_df)} spots")

order = ['control', 'D1', 'D3', 'D7', 'D7_rep']

# -----------------------------------------------------------------------------
# 1. 每个样本内 Pirb+ 与 Pirb- spot 到最近 Pirb+ spot 的距离
# -----------------------------------------------------------------------------
print("[INFO] Computing distance to nearest Pirb+ spot...")

sample_results = []
for sample in spatial_df['sample'].unique():
    sub = spatial_df[spatial_df['sample'] == sample].copy()
    tp = sub['time_point'].iloc[0]
    xy = sub[['imageX', 'imageY']].values
    pos_mask = sub['Pirb_positive'].values.astype(bool)

    if pos_mask.sum() < 2:
        continue

    pos_xy = xy[pos_mask]

    # 所有 spot 到最近 Pirb+ spot 的距离
    dist_to_pos = cdist(xy, pos_xy).min(axis=1)
    sub['dist_to_pirb_pos'] = dist_to_pos

    # Pirb+ spot 之间的最近邻距离（排除自身）
    if pos_mask.sum() >= 2:
        d_pos = cdist(pos_xy, pos_xy)
        np.fill_diagonal(d_pos, np.inf)
        nn_dist_pos = d_pos.min(axis=1)
        mean_nn_pos = nn_dist_pos.mean()
        median_nn_pos = np.median(nn_dist_pos)
    else:
        mean_nn_pos = median_nn_pos = np.nan

    neg_nn_dist = dist_to_pos[~pos_mask]
    sample_results.append({
        'sample': sample,
        'time_point': tp,
        'n_total': len(sub),
        'n_pirb_pos': pos_mask.sum(),
        'frac_pirb_pos': pos_mask.mean(),
        'mean_nn_dist_pirb_pos': mean_nn_pos,
        'median_nn_dist_pirb_pos': median_nn_pos,
        'mean_dist_to_pirb_pos_neg': neg_nn_dist.mean() if len(neg_nn_dist) > 0 else np.nan,
        'median_dist_to_pirb_pos_neg': np.median(neg_nn_dist) if len(neg_nn_dist) > 0 else np.nan,
    })

    # 保存到 sub dataframe
    spatial_df.loc[sub.index, 'dist_to_pirb_pos'] = dist_to_pos

sample_results_df = pd.DataFrame(sample_results)
sample_results_df.to_csv(os.path.join(OUT_DIR, "pirb_spatial_pattern_summary.csv"), index=False)
print("[SAVE] pirb_spatial_pattern_summary.csv")
print(sample_results_df.to_string(index=False))

# -----------------------------------------------------------------------------
# 2. 组织边界距离：计算每个 spot 到组织凸包边界的距离
# -----------------------------------------------------------------------------
print("[INFO] Computing distance to tissue boundary (convex hull)...")

boundary_results = []
for sample in spatial_df['sample'].unique():
    sub = spatial_df[spatial_df['sample'] == sample].copy()
    tp = sub['time_point'].iloc[0]
    xy = sub[['imageX', 'imageY']].values

    if len(xy) < 4:
        continue

    hull = ConvexHull(xy)
    hull_points = xy[hull.vertices]

    # 计算每个 spot 到凸包顶点的最小距离（近似边界距离）
    dist_to_boundary = cdist(xy, hull_points).min(axis=1)
    sub['dist_to_boundary'] = dist_to_boundary
    spatial_df.loc[sub.index, 'dist_to_boundary'] = dist_to_boundary

    pos_mask = sub['Pirb_positive'].values.astype(bool)
    if pos_mask.sum() > 0 and (~pos_mask).sum() > 0:
        stat, pval = mannwhitneyu(
            sub.loc[pos_mask, 'dist_to_boundary'].values,
            sub.loc[~pos_mask, 'dist_to_boundary'].values,
            alternative='two-sided'
        )
    else:
        stat = pval = np.nan

    boundary_results.append({
        'sample': sample,
        'time_point': tp,
        'n_pirb_pos': pos_mask.sum(),
        'mean_dist_boundary_pos': sub.loc[pos_mask, 'dist_to_boundary'].mean() if pos_mask.sum() > 0 else np.nan,
        'mean_dist_boundary_neg': sub.loc[~pos_mask, 'dist_to_boundary'].mean() if (~pos_mask).sum() > 0 else np.nan,
        'mannwhitney_u': stat,
        'p_value': pval,
    })

boundary_df = pd.DataFrame(boundary_results)
boundary_df.to_csv(os.path.join(OUT_DIR, "pirb_boundary_distance_summary.csv"), index=False)
print("[SAVE] pirb_boundary_distance_summary.csv")
print(boundary_df.to_string(index=False))

# -----------------------------------------------------------------------------
# 3. 可视化
# -----------------------------------------------------------------------------
print("[INFO] Plotting spatial pattern results...")

# (a) Pirb+ spot 最近邻距离与阴性 spot 到 Pirb+ 距离对比
plot_df_list = []
for sample in spatial_df['sample'].unique():
    sub = spatial_df[spatial_df['sample'] == sample]
    tp = sub['time_point'].iloc[0]
    pos_mask = sub['Pirb_positive'].values.astype(bool)

    if pos_mask.sum() < 2:
        continue

    # 阳性 spot 之间的 NN 距离
    pos_xy = sub.loc[pos_mask, ['imageX', 'imageY']].values
    d_pos = cdist(pos_xy, pos_xy)
    np.fill_diagonal(d_pos, np.inf)
    nn_pos = d_pos.min(axis=1)

    for d in nn_pos:
        plot_df_list.append({'time_point': tp, 'distance': d, 'group': 'NN among Pirb+'})

    # 阴性 spot 到最近阳性 spot 的距离
    neg_dist = sub.loc[~pos_mask, 'dist_to_pirb_pos'].values
    for d in neg_dist:
        plot_df_list.append({'time_point': tp, 'distance': d, 'group': 'Pirb- to nearest Pirb+'})

plot_df = pd.DataFrame(plot_df_list)

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=plot_df, x='time_point', y='distance', hue='group',
            order=order, ax=ax)
ax.set_title("Spatial clustering of Pirb+ spots\n(smaller NN distance = more clustered)")
ax.set_ylabel("Distance (pixels)")
ax.set_xlabel("Time point")
plt.legend(title='')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_spatial_clustering_boxplot.png"), dpi=300, bbox_inches='tight')
print("  [SAVE] pirb_spatial_clustering_boxplot.png")
plt.close(fig)

# (b) Pirb+ vs Pirb- 到组织边界的距离
boundary_plot_list = []
for sample in spatial_df['sample'].unique():
    sub = spatial_df[spatial_df['sample'] == sample].dropna(subset=['dist_to_boundary'])
    tp = sub['time_point'].iloc[0]
    for _, row in sub.iterrows():
        boundary_plot_list.append({
            'time_point': tp,
            'distance_to_boundary': row['dist_to_boundary'],
            'group': 'Pirb+' if row['Pirb_positive'] else 'Pirb-'
        })

boundary_plot_df = pd.DataFrame(boundary_plot_list)

fig, ax = plt.subplots(figsize=(10, 6))
sns.violinplot(data=boundary_plot_df, x='time_point', y='distance_to_boundary',
               hue='group', order=order, split=True, ax=ax)
ax.set_title("Distance of Pirb+ vs Pirb- spots to tissue boundary")
ax.set_ylabel("Distance to convex hull boundary (pixels)")
ax.set_xlabel("Time point")
plt.legend(title='')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_boundary_distance_violin.png"), dpi=300, bbox_inches='tight')
print("  [SAVE] pirb_boundary_distance_violin.png")
plt.close(fig)

# (c) 边界富集散点：x=时间，y=平均边界距离比值（阳性/阴性）
ratio_df = boundary_df.copy()
ratio_df['boundary_distance_ratio'] = ratio_df['mean_dist_boundary_pos'] / ratio_df['mean_dist_boundary_neg']
ratio_df['time_point'] = pd.Categorical(ratio_df['time_point'], categories=order, ordered=True)

fig, ax = plt.subplots(figsize=(8, 5))
sns.stripplot(data=ratio_df, x='time_point', y='boundary_distance_ratio', size=10, color='darkred', ax=ax)
ax.axhline(1.0, color='black', linestyle='--', linewidth=1)
ax.set_title("Pirb+ spots boundary enrichment\nratio < 1: Pirb+ closer to boundary")
ax.set_ylabel("Mean distance ratio (Pirb+ / Pirb-)")
ax.set_xlabel("Time point")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_boundary_enrichment_ratio.png"), dpi=300, bbox_inches='tight')
print("  [SAVE] pirb_boundary_enrichment_ratio.png")
plt.close(fig)

# -----------------------------------------------------------------------------
# 4. 按边界距离分箱，看 Pirb 阳性率是否随边界距离变化
# -----------------------------------------------------------------------------
print("[INFO] Boundary distance binning analysis...")
bin_results = []
for sample in spatial_df['sample'].unique():
    sub = spatial_df[spatial_df['sample'] == sample].dropna(subset=['dist_to_boundary']).copy()
    tp = sub['time_point'].iloc[0]
    sub['boundary_quintile'] = pd.qcut(sub['dist_to_boundary'], q=5, labels=['Q1_edge','Q2','Q3','Q4','Q5_center'])
    qsum = sub.groupby('boundary_quintile', observed=False).agg(
        n=('Pirb_positive', 'size'),
        pirb_pos=('Pirb_positive', 'sum'),
        pirb_frac=('Pirb_positive', 'mean'),
    ).reset_index()
    qsum['sample'] = sample
    qsum['time_point'] = tp
    bin_results.append(qsum)

bin_df = pd.concat(bin_results, ignore_index=True)
bin_df.to_csv(os.path.join(OUT_DIR, "pirb_boundary_quintile_summary.csv"), index=False)
print("[SAVE] pirb_boundary_quintile_summary.csv")

# 可视化：每个时间点边界分箱的 Pirb 阳性率
fig, ax = plt.subplots(figsize=(10, 6))
sns.pointplot(data=bin_df, x='boundary_quintile', y='pirb_frac', hue='time_point',
              hue_order=order, dodge=0.2, ax=ax)
ax.set_title("Pirb positivity across boundary-to-center quintiles")
ax.set_ylabel("Pirb+ fraction")
ax.set_xlabel("Boundary distance quintile")
plt.legend(title='Time point')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "pirb_boundary_quintile_fraction.png"), dpi=300, bbox_inches='tight')
print("  [SAVE] pirb_boundary_quintile_fraction.png")
plt.close(fig)

print("[DONE]")
