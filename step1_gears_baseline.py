"""Step 1: GEARS OOD baseline on Replogle-Nadig K562 (perturbation holdout).

产出 outputs/predictions.csv，列：
    perturbation, pearson_delta, uncertainty, strength, pred_seed0..2 (npy 另存)

其中 strength = 扰动强度协变量（敲低效率或扰动组 vs 对照组的全基因组效应量，
如 L2 距离），供 Step 4 的 H2a-partial 偏相关用（预注册 2026-10-12 增补，
在预注册锁定之前定义，属协议内）。

数据来源（公开，二选一）：
  A) GEARS 官方 repo 的数据下载脚本（github.com/snap-stanford/GEARS）
  B) scPerturBench 复现包（github.com/bm2-lab/scPerturBench）——推荐，
     其 o.o.d. 划分协议与本预注册锁定设置对齐
注意：本脚本为骨架，GEARS 具体 API 以 0.1.x 文档为准，未经端到端测试。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

OUT = "outputs"
os.makedirs(OUT, exist_ok=True)

# ---- 1. 加载 Replogle K562（AnnData，.obs 含 'gene' 扰动列）----
# adata = sc.read_h5ad("data/replogle_k562.h5ad")
raise SystemExit("请先将 Replogle K562 数据放入 data/ 并取消上面的注释")

# ---- 2. 扰动留出划分（对齐 scPerturBench perturbation generalization）----
# holdout_frac = 0.2  # 与 scPerturBench 一致
# test_genes = rng.choice(unique_genes, size=int(holdout_frac * n_genes), replace=False)
# train_adata / test_adata 按 test_genes 切分

# ---- 3. 三种子训练 GEARS ----
from gears import GEARS

seeds = [0, 1, 2]
predictions = {}
for seed in seeds:
    gears = GEARS(train_adata, {"perturbation": "gene"})
    gears.model_initialize()
    gears.train()                      # 默认超参，与预注册一致
    predictions[seed] = gears.predict(test_adata)   # dict: gene -> predicted expr vector

# ---- 4. top-DEG Pearson delta（对齐 scPerturBench/MAP 口径）----
def top_deg_pearson_delta(pred, true, n_top=50):
    """预测与真实表达在 top-n DEG 上的 Pearson 相关（delta 口径按 scPerturBench 实现）"""
    raise NotImplementedError("按 scPerturBench 的 PCC-delta 实现补齐")

rows = []
for gene in test_genes:
    deltas = [top_deg_pearson_delta(predictions[s][gene], true_expr[gene]) for s in seeds]
    rows.append({"perturbation": gene, "pearson_delta": np.mean(deltas),
                 "uncertainty": np.var(deltas)})
    for s in seeds:
        np.save(f"{OUT}/pred_{gene}_seed{s}.npy", predictions[s][gene])

pd.DataFrame(rows).to_csv(f"{OUT}/predictions.csv", index=False)
print("Step 1 done -> outputs/predictions.csv")
