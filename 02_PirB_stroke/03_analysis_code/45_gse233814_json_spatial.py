"""
GSE233814 Visium 空间精确解析：
从 json.gz 解析 fiducial/oligo 坐标与变换矩阵，将 barcode 映射到 (row,col)，
再对齐 oligo 的 imageX/imageY 像素坐标，最终把 Pirb 表达映射到组织空间位置，
生成高分辨率空间分布图。
"""
import os, gzip, json, gc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

sc.settings.verbosity = 3

DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE233814"
WHITELIST = "D:/Pirb_stroke_project/01_raw_data/visium-v1.txt"
COORD_FILE = "D:/Pirb_stroke_project/01_raw_data/visium-v1_coordinates.txt"
H5AD_PATH = "D:/Pirb_stroke_project/04_reports/figures/GSE233814_processed.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233814"
os.makedirs(OUT_DIR, exist_ok=True)

# 样本到 json 与显示名称的映射
sample_info = {
    "GSM7437221_C1-control": {
        "time_point": "control",
        "json": "GSM7437221_mouse_control_V11A20-380-C1.json.gz",
    },
    "GSM7437222_B1-D1": {
        "time_point": "D1",
        "json": "GSM7437222_V10M09-040-B1.json.gz",
    },
    "GSM7437223_D1-D3": {
        "time_point": "D3",
        "json": "GSM7437223_V10M09-040-D1.json.gz",
    },
    "GSM7437224_C1-D7": {
        "time_point": "D7",
        "json": "GSM7437224_V10M09-040-C1.json.gz",
    },
    "GSM7437225_D1-D7": {
        "time_point": "D7_rep",
        "json": "GSM7437225_mouse_D7_V11A20-380-D1.json.gz",
    },
}

order = ['control', 'D1', 'D3', 'D7', 'D7_rep']

# -----------------------------------------------------------------------------
# 1. 构建 barcode -> (row, col) 映射（coordinates 文件是 1-based，转为 0-based）
# -----------------------------------------------------------------------------
print("[INFO] Loading Visium v1 barcode coordinate whitelist...")
coord_df = pd.read_csv(COORD_FILE, sep='\t', header=None,
                       names=['barcode', 'col_1based', 'row_1based'])
coord_df['row'] = coord_df['row_1based'] - 1
coord_df['col'] = coord_df['col_1based'] - 1
barcode_to_rc = dict(zip(coord_df['barcode'].values,
                         zip(coord_df['row'].values, coord_df['col'].values)))
print(f"[INFO] Whitelist barcodes: {len(barcode_to_rc)}")

# -----------------------------------------------------------------------------
# 2. 读取已处理 h5ad（包含 QC 后的表达矩阵）
# -----------------------------------------------------------------------------
print(f"[INFO] Reading {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"[INFO] Loaded: {adata.n_obs} spots x {adata.n_vars} genes")

# 确保 Pirb 表达列存在
if 'Pirb_expr' not in adata.obs.columns:
    adata.obs['Pirb_expr'] = adata[:, 'Pirb'].X.toarray().flatten()
if 'Pirb_positive' not in adata.obs.columns:
    adata.obs['Pirb_positive'] = (adata.obs['Pirb_expr'] > 0).astype(int)

# -----------------------------------------------------------------------------
# 3. 为每个样本添加空间像素坐标并校验
# -----------------------------------------------------------------------------
print("[INFO] Mapping barcodes to image pixel coordinates...")
spatial_records = []  # 用于保存每个 spot 的元信息

for sample, info in sample_info.items():
    json_name = info['json']
    time_point = info['time_point']
    json_path = os.path.join(DATA_DIR, json_name)

    # 读取 json.gz
    with gzip.open(json_path, 'rt') as f:
        json_data = json.load(f)

    # oligo: (row, col) -> (imageX, imageY)
    oligo_map = {(o['row'], o['col']): (o['imageX'], o['imageY'])
                 for o in json_data['oligo']}

    # fiducial 信息（用于报告/可视化参考）
    fiducial_df = pd.DataFrame(json_data['fiducial'])
    fiducial_df['sample'] = sample

    # 取该样本的 spots
    mask = adata.obs['sample'] == sample
    barcodes = adata.obs.index[mask].str.replace("-1", "", regex=False).values

    coords = []
    unmatched = 0
    for bc in barcodes:
        rc = barcode_to_rc.get(bc)
        if rc is None or rc not in oligo_map:
            coords.append((np.nan, np.nan))
            unmatched += 1
        else:
            coords.append(oligo_map[rc])

    coords = np.array(coords, dtype=float)
    if 'spatial' not in adata.obsm:
        adata.obsm['spatial'] = np.full((adata.n_obs, 2), np.nan)
    adata.obsm['spatial'][mask] = coords

    print(f"  {sample} ({time_point}): {mask.sum()} spots, "
          f"{mask.sum() - unmatched} matched to image coords")

    # 保存记录
    sample_obs = adata.obs[mask].copy()
    sample_obs['imageX'] = coords[:, 0]
    sample_obs['imageY'] = coords[:, 1]
    sample_obs = sample_obs.reset_index().rename(columns={'index': 'barcode'})
    spatial_records.append(sample_obs[['barcode', 'sample', 'time_point',
                                       'imageX', 'imageY',
                                       'Pirb_expr', 'Pirb_positive']])

spatial_df = pd.concat(spatial_records, ignore_index=True)
spatial_df.to_csv(os.path.join(OUT_DIR, "spot_pixel_coordinates_pirb.csv"),
                  index=False)
print(f"[SAVE] spot_pixel_coordinates_pirb.csv ({len(spatial_df)} spots)")

# -----------------------------------------------------------------------------
# 4. 可视化：每个样本的 Pirb 表达空间热图
# -----------------------------------------------------------------------------
print("[INFO] Drawing high-resolution Pirb spatial maps...")

# 自定义热图颜色：黑 -> 蓝 -> 黄 -> 红
cmap_expr = LinearSegmentedColormap.from_list(
    "pirb", ["#000000", "#1f4e79", "#ffd700", "#ff0000"])

for sample, info in sample_info.items():
    time_point = info['time_point']
    mask = adata.obs['sample'] == sample
    sub = adata[mask].copy()

    xy = sub.obsm['spatial']
    valid = ~np.isnan(xy).any(axis=1)
    sub = sub[valid].copy()
    xy = xy[valid]

    expr = sub.obs['Pirb_expr'].values
    pos = sub.obs['Pirb_positive'].values.astype(bool)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # (a) 所有 spots 的 Pirb 表达热图
    ax = axes[0]
    scat = ax.scatter(xy[:, 0], xy[:, 1], c=expr, cmap=cmap_expr,
                      s=35, edgecolors='none', vmin=0, vmax=max(expr.max(), 1.0))
    ax.set_title(f"{sample} | {time_point}\nPirb expression (log1p CPM)")
    ax.set_xlabel("imageX (pixels)")
    ax.set_ylabel("imageY (pixels)")
    ax.invert_yaxis()  # 图像坐标原点在左上角
    plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04)

    # (b) 仅 Pirb+ spots（红色）
    ax = axes[1]
    ax.scatter(xy[:, 0], xy[:, 1], c='#cccccc', s=8, edgecolors='none', label='Pirb-')
    if pos.sum() > 0:
        ax.scatter(xy[pos, 0], xy[pos, 1], c='#e41a1c', s=45,
                   edgecolors='black', linewidths=0.3, label=f'Pirb+ (n={pos.sum()})')
    ax.set_title(f"Pirb+ spots\n{pos.sum()}/{len(pos)} ({pos.mean()*100:.2f}%)")
    ax.set_xlabel("imageX (pixels)")
    ax.set_ylabel("imageY (pixels)")
    ax.invert_yaxis()
    ax.legend(loc='upper right')

    # (c) Pirb 表达的分位数染色（便于观察弱信号空间分布）
    ax = axes[2]
    q90 = np.percentile(expr, 90)
    scat = ax.scatter(xy[:, 0], xy[:, 1], c=expr, cmap='magma',
                      s=35, edgecolors='none', vmin=0, vmax=max(q90, 0.5))
    ax.set_title(f"Pirb expression (capped at 90th percentile={q90:.2f})")
    ax.set_xlabel("imageX (pixels)")
    ax.set_ylabel("imageY (pixels)")
    ax.invert_yaxis()
    plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04)

    plt.tight_layout()
    outname = os.path.join(OUT_DIR, f"pirb_spatial_{sample}_{time_point}.png")
    fig.savefig(outname, dpi=300, bbox_inches='tight')
    print(f"  [SAVE] {outname}")
    plt.close(fig)

    gc.collect()

# -----------------------------------------------------------------------------
# 5. 组合图：5 个样本并排，统一颜色尺度
# -----------------------------------------------------------------------------
print("[INFO] Drawing combined spatial panel...")
vmax_global = max(1.0, spatial_df['Pirb_expr'].max())

fig, axes = plt.subplots(2, len(order), figsize=(22, 12))
for col_idx, tp in enumerate(order):
    sample = [s for s, v in sample_info.items() if v['time_point'] == tp][0]
    mask = adata.obs['sample'] == sample
    sub = adata[mask].copy()
    xy = sub.obsm['spatial']
    valid = ~np.isnan(xy).any(axis=1)
    sub = sub[valid]
    xy = xy[valid]
    expr = sub.obs['Pirb_expr'].values
    pos = sub.obs['Pirb_positive'].values.astype(bool)

    # 上行：表达热图
    ax = axes[0, col_idx]
    scat = ax.scatter(xy[:, 0], xy[:, 1], c=expr, cmap=cmap_expr,
                      s=30, edgecolors='none', vmin=0, vmax=vmax_global)
    ax.set_title(f"{tp}\nn={len(sub)} spots")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    if col_idx == len(order) - 1:
        plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04,
                     label='Pirb log1p CPM')

    # 下行：Pirb+ 分布
    ax = axes[1, col_idx]
    ax.scatter(xy[:, 0], xy[:, 1], c='#dddddd', s=6, edgecolors='none')
    if pos.sum() > 0:
        ax.scatter(xy[pos, 0], xy[pos, 1], c='#e41a1c', s=40,
                   edgecolors='black', linewidths=0.3)
    ax.set_title(f"Pirb+ {pos.sum()}/{len(pos)} ({pos.mean()*100:.2f}%)")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])

axes[0, 0].set_ylabel("Pirb expression", fontsize=12)
axes[1, 0].set_ylabel("Pirb+ spots", fontsize=12)
plt.suptitle("GSE233814 Visium: Pirb spatial distribution across MCAO time points",
             fontsize=14, y=1.02)
plt.tight_layout()
combined_path = os.path.join(OUT_DIR, "pirb_spatial_combined_panel.png")
fig.savefig(combined_path, dpi=300, bbox_inches='tight')
print(f"  [SAVE] {combined_path}")
plt.close(fig)

# -----------------------------------------------------------------------------
# 6. 提取并保存 json 元信息（变换矩阵、fiducial 摘要）
# -----------------------------------------------------------------------------
print("[INFO] Summarizing json spatial metadata...")
meta_records = []
for sample, info in sample_info.items():
    json_path = os.path.join(DATA_DIR, info['json'])
    with gzip.open(json_path, 'rt') as f:
        data = json.load(f)
    transform = np.array(data['transform'])
    meta_records.append({
        'sample': sample,
        'time_point': info['time_point'],
        'serialNumber': data.get('serialNumber'),
        'area': data.get('area'),
        'n_fiducial': len(data['fiducial']),
        'n_oligo': len(data['oligo']),
        'transform_a11': transform[0, 0],
        'transform_a12': transform[0, 1],
        'transform_a13': transform[0, 2],
        'transform_a21': transform[1, 0],
        'transform_a22': transform[1, 1],
        'transform_a23': transform[1, 2],
    })
meta_df = pd.DataFrame(meta_records)
meta_df.to_csv(os.path.join(OUT_DIR, "spatial_json_metadata.csv"), index=False)
print(f"[SAVE] spatial_json_metadata.csv")
print(meta_df[['sample', 'time_point', 'serialNumber', 'area',
               'n_fiducial', 'n_oligo']].to_string(index=False))

print("[DONE]")
