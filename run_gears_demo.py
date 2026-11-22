"""GEARS 官方 demo 打通脚本（Day 1 链路验证）。

两级策略：
  1) 官方路径：PertData 下载 norman 小数据集 → prepare_split → model_initialize → train → predict
  2) 降级路径：若官方数据下载失败，构造合成 AnnData → new_data_process → 同流程

通过标准：GEARS.predict() 返回非空预测（官方 demo 或合成数据均可，如实记录用的是哪条）。
"""
import os
import sys
import traceback

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

DATA_DIR = "./data/gears_demo"
os.makedirs(DATA_DIR, exist_ok=True)


def try_official_demo():
    """官方 tutorial 路径：下载 norman 数据集。"""
    from gears import PertData, GEARS

    pert_data = PertData(DATA_DIR)
    pert_data.load(data_name="norman")  # 需从外网下载，可能失败
    pert_data.prepare_split(split="simulation", seed=1)
    pert_data.get_dataloader(batch_size=32, test_batch_size=128)

    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[demo] device = {device}")
    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=64)
    model.train(epochs=3)  # demo 只跑 3 个 epoch 验证链路
    pred = model.predict([["FEV"], ["CBL", "CNN1"]])
    return model, pred


def make_synthetic_adata(n_genes=300, n_cells_per_cond=40, seed=0):
    """构造满足 GEARS 格式要求的合成 AnnData。

    GEARS new_data_process 要求：adata.var 有 'gene_name'；adata.obs 有
    'condition' 与 'cell_type'。
    condition 编码约定（GEARS 内部约定，见 gears/utils.py parse_single_pert）：
      'ctrl'           = 对照
      'GENE+ctrl'      = 单基因扰动（必须含 '+'，否则会被当成 combo 解析而报 IndexError）
      'GENE1+GENE2'    = 组合扰动
    """
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    # 给一部分基因真实符号名，测试 graph 构建对符号名的要求
    for i, g in enumerate(["TP53", "MYC", "KRAS", "CBL", "FEV", "CNN1", "PTEN", "AKT1"]):
        genes[i] = g

    conds = (["ctrl"]
             + [f"{g}+ctrl" for g in genes[:20]]   # 单扰动：GEARS 约定 'GENE+ctrl'
             + ["TP53+MYC", "CBL+CNN1"])           # 组合扰动
    X, obs_cond, obs_ct = [], [], []
    base = rng.normal(0, 1, n_genes)
    for ci, cond in enumerate(conds):
        for _ in range(n_cells_per_cond):
            x = base + rng.normal(0, 0.3, n_genes)
            if cond != "ctrl":
                pert_genes = [g for g in cond.split("+") if g != "ctrl"]
                idx = [genes.index(g) for g in pert_genes if g in genes]
                x[idx] -= 2.0  # 敲低效应
            X.append(x)
            obs_cond.append(cond)
            obs_ct.append("synth")
    adata = sc.AnnData(sparse.csr_matrix(np.asarray(X)))
    # 注意：不要传 pd.Categorical。pandas 3 下行 apply(逐行) 时 categorical 列
    # 会产出 unhashable 的 Categorical 标量（GEARS get_DE_genes 内部会踩炸），
    # 用普通 object 列即可，GEARS 自己会 astype('category')。
    adata.obs["condition"] = obs_cond
    adata.obs["cell_type"] = obs_ct
    adata.var["gene_name"] = genes
    return adata


def try_synthetic_demo():
    """降级路径：合成 AnnData 走通 new_data_process + train + predict。"""
    from gears import PertData, GEARS

    adata = make_synthetic_adata()
    pert_data = PertData(DATA_DIR)
    pert_data.new_data_process(dataset_name="synth", adata=adata)
    pert_data.prepare_split(split="simulation", seed=1)
    pert_data.get_dataloader(batch_size=32, test_batch_size=128)

    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[demo-synth] device = {device}")
    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=64)
    model.train(epochs=3)
    pred = model.predict([["FEV"], ["CBL", "CNN1"]])
    return model, pred


if __name__ == "__main__":
    which = None
    if "--synthetic" in sys.argv:
        print("=== 直接走合成数据降级路径（--synthetic）===")
        model, pred = try_synthetic_demo()
        which = "synthetic"
    else:
        try:
            print("=== 尝试官方 demo（norman 数据集）===")
            model, pred = try_official_demo()
            which = "official(norman)"
        except Exception:
            print("!!! 官方 demo 失败，报错如下，转入合成数据降级路径 !!!")
            traceback.print_exc()
            print("=== 尝试合成数据降级路径 ===")
            model, pred = try_synthetic_demo()
            which = "synthetic"

    print(f"\n=== demo 跑通，路径 = {which} ===")
    print(f"predict 返回类型: {type(pred)}")
    if isinstance(pred, dict):
        for k, v in list(pred.items())[:5]:
            arr = np.asarray(v)
            print(f"  {k}: shape={arr.shape}, finite={np.isfinite(arr).all()}")
    assert pred is not None and len(pred) > 0, "predict 返回为空！"
    print("PASS: model_initialize + train + predict 全流程完成")
