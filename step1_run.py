"""Day 2: Replogle-Nadig K562 CRISPRi + 扰动留出（OOD）+ 首批预测。

对应预注册评测设置：
- 数据：GEARS 提供的 replogle_k562_essential 子集（dataverse 7458695，
  "not the complete dataset and has been filtered" —— Day 2 先用它跑通；
  完整 Replogle 数据 ~60GB，后续再评估）。
- OOD 划分：扰动留出。GEARS prepare_split(split='simulation',
  only_test_set_perts=True, test_perts=<留出基因>) —— 留出基因的全部细胞进测试集，
  训练集完全不含这些扰动。留出比例 0.2（scPerturBench perturbation
  generalization 的确切划分名单未能从 Figshare 复现包直接取得，记协议偏离，
  见 setup_log.md Day2 章节）。
- 引擎：GEARS，先单种子（seed 0）跑通；3 种子 ensemble 是预注册锁定设置，
  Day 2 时间盒内先出种子 0 的数字。
- 误差：每个留出扰动计算两种口径并列报告——
    pearson_delta_all   全基因 Pearson(pred_delta, true_delta)
    pearson_delta_top50 top-50 |true_delta| 基因上的 Pearson（近似 top-DEG 口径；
    GEARS 内部用 top-20 DE，scPerturBench/MAP 的 PCC-delta 为全基因，
    故两口径并列，孰进 verdict 待议）

用法：
    python step1_run.py            # 默认 10 epochs, seed 0
    EPOCHS=5 SEED=1 python step1_run.py
"""
import os
import time

import numpy as np
import pandas as pd

DATA_DIR = "./data/gears_demo"
# [D8 跨数据集验证] 数据集与输出目录参数化：默认 K562（归档不变），
# 跨数据集验证用 DATA=replogle_rpe1 OUT_DIR=outputs_rpe1，避免覆盖已归档产物。
DATA = os.environ.get("DATA", "replogle_k562_essential")
OUT = os.environ.get("OUT_DIR", "outputs")
os.makedirs(OUT, exist_ok=True)

HOLDOUT_FRAC = 0.2
SEED = int(os.environ.get("SEED", "0"))        # 模型初始化种子（ensemble 成员间变化来源）
# [修复 D5] 原实现 SEED 同时驱动留出划分与模型初始化，导致不同种子的
# 留出集不同（seed0∩seed1 仅 45/217），种子间方差失去意义。
# 现拆分：SPLIT_SEED 固定留出划分（ensemble 成员必须在同一考场），
# SEED 只驱动模型初始化。重训 seed1/2 时 SPLIT_SEED 必须 = 0（与 Day2 seed0 同考场）。
SPLIT_SEED = int(os.environ.get("SPLIT_SEED", "0"))
EPOCHS = int(os.environ.get("EPOCHS", "10"))
HIDDEN = 64
TOP_N_DEG = 50

t0 = time.time()
from gears import PertData, GEARS  # noqa: E402

# ---------- 1. 加载数据（zip 需提前 curl 下载并解压，见 setup_log.md E9/E11） ----------
pert_data = PertData(DATA_DIR)
pert_data.load(data_name=DATA)
adata = pert_data.adata
print(f"[data] cells={adata.n_obs}, genes={adata.n_vars}, "
      f"conditions={adata.obs.condition.nunique()} "
      f"({time.time()-t0:.0f}s)")

# ---------- 2. 扰动留出划分 ----------
conditions = adata.obs.condition.unique().tolist()
singles = [c for c in conditions if ("ctrl" in c and c != "ctrl")]
genes = sorted({c.replace("ctrl+", "").replace("+ctrl", "") for c in singles})
bad = [g for g in genes if "_" in g]
if bad:
    print(f"[split] 警告：{len(bad)} 个基因名含 '_'，prepare_split 的 test_perts "
          f"按 '_' 拼接会有歧义，已剔除: {bad[:5]}")
    genes = [g for g in genes if "_" not in g]

rng = np.random.default_rng(SPLIT_SEED)
holdout_genes = rng.choice(genes, size=int(HOLDOUT_FRAC * len(genes)),
                           replace=False).tolist()
print(f"[split] 总单扰动基因={len(genes)}, 留出={len(holdout_genes)} "
      f"(frac={HOLDOUT_FRAC}, split_seed={SPLIT_SEED}, model_seed={SEED})")

test_perts_arg = "_".join(holdout_genes)
pert_data.prepare_split(split="simulation", seed=SPLIT_SEED,
                        test_perts=test_perts_arg,
                        only_test_set_perts=True)
split_info = {k: len(v) for k, v in pert_data.set2conditions.items()}
print(f"[split] set2conditions: {split_info}")

# 留出扰动的 condition 名（'GENE+ctrl' 形式）
holdout_conds = [c for c in singles
                 if c.replace("ctrl+", "").replace("+ctrl", "") in holdout_genes]
# batch 调小：run3 OOM 教训（D2-3），单进程 + 小 batch 保内存安全
BATCH = int(os.environ.get("BATCH", "32"))
TEST_BATCH = int(os.environ.get("TEST_BATCH", "128"))
pert_data.get_dataloader(batch_size=BATCH, test_batch_size=TEST_BATCH)

# ---------- 3. 训练 ----------
import torch  # noqa: E402
torch.manual_seed(SEED)  # [修复 D5] 固定模型初始化种子，ensemble 成员间仅初始化不同
np.random.seed(SEED)
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"[train] device={device}, epochs={EPOCHS}, hidden={HIDDEN}")
model = GEARS(pert_data, device=device)
model.model_initialize(hidden_size=HIDDEN)
model.train(epochs=EPOCHS)
print(f"[train] done ({time.time()-t0:.0f}s)")

# ---------- 4. 预测留出扰动 + pearson delta + strength ----------
# 契约列：perturbation, pearson_delta, uncertainty, strength（见 setup_log.md D2-4）
ctrl = adata[adata.obs.condition == "ctrl"]
ctrl_mean = np.asarray(ctrl.X.mean(axis=0)).ravel()

# strength 协变量（不重训，直接从数据算）：
#   每个留出扰动组内，各细胞到对照均值的 L2 距离取平均 —— 全基因组效应量
strength_map = {}
for g in holdout_genes:
    cand = [c for c in singles if c.startswith(g + "+") or c.endswith("+" + g)]
    if not cand:
        continue
    sub = adata[adata.obs.condition == cand[0]]
    d = sub.X.toarray() - ctrl_mean          # (n_cells, n_genes) dense
    strength_map[g] = float(np.linalg.norm(d, axis=1).mean())

# predict 对不在扰动图（model.pert_list）中的基因会直接 raise，先过滤并留痕
in_graph = [g for g in holdout_genes if g in model.pert_list]
skipped = [g for g in holdout_genes if g not in model.pert_list]
if skipped:
    print(f"[eval] {len(skipped)}/{len(holdout_genes)} 个留出基因不在扰动图中，"
          f"无法 predict（记录为 NaN），例如: {skipped[:5]}")
if not in_graph:
    raise SystemExit("留出基因全部不在扰动图中，无法预测")

query = [[g] for g in in_graph]
pred = model.predict(query)  # dict: gene -> predicted expression vector (n_genes,)

# ---- Day3 增补：固化模型与预测表达向量（酶表达代理的输入）----
if os.environ.get("SAVE_MODEL", "1") == "1":
    model.save_model(f"{OUT}/gears_model_seed{SEED}")
    print(f"[save] model -> {OUT}/gears_model_seed{SEED}")
# 预测向量按 model.gene_list 顺序对齐（与 adata.var.gene_name 同序）
pred_mat = np.stack([np.asarray(pred[g]).ravel() for g in in_graph])
np.savez_compressed(
    f"{OUT}/pred_expr_seed{SEED}.npz",
    pred_expr=pred_mat.astype(np.float32),
    pred_perturbations=np.array(in_graph),
    pred_genes=np.array(list(model.gene_list)),
)
print(f"[save] pred_expr {pred_mat.shape} -> {OUT}/pred_expr_seed{SEED}.npz")


rows = []
for g in holdout_genes:
    if g not in pred:
        rows.append({"perturbation": g, "pearson_delta": np.nan,
                     "pearson_delta_all": np.nan, "uncertainty": np.nan,
                     "strength": strength_map.get(g, np.nan),
                     "n_cells": 0, "in_pert_graph": g in model.pert_list})
        continue
    p = np.asarray(pred[g]).ravel()
    cond = [c for c in singles if c.startswith(g + "+") or c.endswith("+" + g)][0]
    sub = adata[adata.obs.condition == cond]
    true_mean = np.asarray(sub.X.mean(axis=0)).ravel()
    n_cells = sub.n_obs
    p_delta = p - ctrl_mean
    t_delta = true_mean - ctrl_mean
    r_all = np.corrcoef(p_delta, t_delta)[0, 1]
    # pearson_delta（主列）：top-N |真实 delta| 差异基因上 Pearson（对齐 top-DEG 口径）
    top_idx = np.argsort(-np.abs(t_delta))[:TOP_N_DEG]
    r_top = np.corrcoef(p_delta[top_idx], t_delta[top_idx])[0, 1]
    rows.append({"perturbation": g, "pearson_delta": r_top,
                 "pearson_delta_all": r_all, "uncertainty": np.nan,  # 单种子，见 D2-4
                 "strength": strength_map.get(g, np.nan),
                 "n_cells": n_cells, "in_pert_graph": True})

df = pd.DataFrame(rows)
out_csv = f"{OUT}/predictions_seed{SEED}.csv"
df.to_csv(out_csv, index=False)
valid = df.dropna(subset=["pearson_delta"])
print(f"[eval] 留出扰动={len(df)}, 有效预测={len(valid)}")
print(f"[eval] pearson_delta(top{TOP_N_DEG}DEG): mean={valid.pearson_delta.mean():.4f} "
      f"median={valid.pearson_delta.median():.4f}")
print(f"[eval] pearson_delta_all: mean={valid.pearson_delta_all.mean():.4f} "
      f"median={valid.pearson_delta_all.median():.4f}")
print(f"[eval] strength: mean={valid.strength.mean():.2f}")
print(f"[eval] -> {out_csv}  (total {time.time()-t0:.0f}s)")
