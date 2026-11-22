"""Model-side feature construction for the reliability gate.

预测谱特征（每行一个留出扰动），与虚拟细胞门控研究的预注册口径一致
（step6_gating_probe.py 的发布版）。红线：只接收模型侧输入
（预测表达谱、ctrl 训练均值、协变量），绝不接收真值派生列。
"""
import numpy as np
import pandas as pd

# 预测 DE 阈值（log-normalized 表达空间）：|delta| > 0.25 ≈ log2 倍数变化 19%
DE_THRESH = 0.25

SPECTRUM_FEATURES = [
    "delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
    "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
]
COVARIATES = ["log_n_cells", "in_pert_graph"]


def compute_features(pred_expr, perts, pred_genes, ctrl_mean):
    """从预测表达谱构造模型侧特征。

    pred_expr: (n_pert, n_genes) float32，预测表达（log-normalized）
    perts:     (n_pert,) 扰动名，与 pred_expr 行对齐
    pred_genes:(n_genes,) 基因名，与 pred_expr 列对齐
    ctrl_mean: (n_genes,) ctrl 训练均值（调用方负责按 pred_genes 对齐）
    返回 DataFrame：index=perturbation，列为 SPECTRUM_FEATURES。
    """
    pred_expr = np.asarray(pred_expr, dtype=np.float32)
    ctrl_mean = np.asarray(ctrl_mean, dtype=np.float32)
    if pred_expr.ndim != 2:
        raise ValueError("pred_expr must be 2-D (n_pert, n_genes)")
    if pred_expr.shape[1] != len(pred_genes):
        raise ValueError("pred_expr columns != len(pred_genes)")
    if len(ctrl_mean) != len(pred_genes):
        raise ValueError("ctrl_mean length != len(pred_genes)")

    delta = pred_expr - ctrl_mean                     # (n_pert, n_genes)
    abs_d = np.abs(delta)

    delta_l2 = np.linalg.norm(delta, axis=1)
    delta_l1 = abs_d.sum(axis=1)
    delta_max = abs_d.max(axis=1)
    top50 = np.sort(abs_d, axis=1)[:, -50:]
    delta_top50_mean = top50.mean(axis=1)
    n_de = (abs_d > DE_THRESH).sum(axis=1).astype(np.float32)
    frac_de = n_de / delta.shape[1]
    neg_frac = (delta < 0).mean(axis=1)

    # 归一化 Shannon 熵（1 = 变化均匀散布全基因组，0 = 集中于单个基因）
    p = abs_d / abs_d.sum(axis=1, keepdims=True)
    H = -(p * np.log(p + 1e-12)).sum(axis=1) / np.log(delta.shape[1])

    return pd.DataFrame({
        "delta_l2": delta_l2,
        "delta_l1": delta_l1,
        "delta_max": delta_max,
        "delta_top50_mean": delta_top50_mean,
        "n_de_025": n_de,
        "frac_de_025": frac_de,
        "neg_frac": neg_frac,
        "delta_entropy": H,
    }, index=pd.Index(list(perts), name="perturbation"))
