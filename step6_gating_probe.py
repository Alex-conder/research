"""Day5 Step6（空位 3：可靠性门控）：门控探针的纯计算核心。

目的：用模型自身的廉价特征（预测表达谱派生 + 实验设计协变量）训练一个
logistic 回归探针，预测 "GEARS 对某个留出扰动会不会答错"
（error = 1 - pearson_delta 高于中位数）。

泄漏纪律（本文件的红线）：
- 只允许「模型侧」特征：预测表达谱 pred_expr、ctrl 训练均值、n_cells、
  in_pert_graph。
- **strength 已确认排除**：来源核验见 step1_run.py 第 96-105 行——
  strength_map 由留出扰动组的真实表达 sub.X 到 ctrl 均值的 L2 距离算出，
  是**真值派生**协变量，用作门控特征会造成泄漏（它本质上与 error 共享
  真实 delta 信息）。本模块不接收 strength 列。
- pearson_delta / pearson_delta_all / uncertainty 同理只作标签或弃用，
  绝不进特征矩阵。

标签定义（定义性选择，写死）：二分类 y = (error > median(error))，
即"答错"= 该留出扰动的预测质量落在全体留出扰动的后一半。

判据（定义性阈值）：5-fold CV AUROC 均值 > 0.6 视为"门控信号存在"
（0.5 = 随机；该线为项目内约定阈值，与 CoT benchmark ECE 探针达标线
同一量级，见 outputs/gating_verdict.md）。

用法：被 step6_run.py 调用；也可
    python step6_gating_probe.py   # 自测：合成数据冒烟测试
"""
import numpy as np
import pandas as pd

# 预测 DE 阈值（log-normalized 表达空间）：|delta| > 0.25 ≈ log2 倍数变化 19%
DE_THRESH = 0.25


def compute_features(pred_expr, perts, pred_genes, ctrl_mean):
    """从预测表达谱构造模型侧特征（每行一个留出扰动）。

    pred_expr: (n_pert, n_genes) float32，预测表达（log-normalized）
    perts:     (n_pert,) 扰动基因名，与 pred_expr 行对齐
    pred_genes:(n_genes,) 基因名，与 pred_expr 列对齐
    ctrl_mean: (n_genes,) ctrl 训练均值（adata.var 顺序，调用方负责按
               pred_genes 对齐）
    返回 DataFrame：index=perturbation，列为下方 FEATURES 中纯谱特征。
    """
    delta = pred_expr - ctrl_mean.astype(np.float32)        # (n_pert, n_genes)
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

    feat = pd.DataFrame({
        "perturbation": list(perts),
        "delta_l2": delta_l2,
        "delta_l1": delta_l1,
        "delta_max": delta_max,
        "delta_top50_mean": delta_top50_mean,
        "n_de_025": n_de,
        "frac_de_025": frac_de,
        "neg_frac": neg_frac,
        "delta_entropy": H,
    }).set_index("perturbation")
    return feat


SPECTRUM_FEATURES = [
    "delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
    "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
]
# 协变量（实验设计，非表达真值）：n_cells 取 log1p；in_pert_graph 在本数据
# 全为 True（常数无方差），调用方验证后自动剔除并留痕。
COVARIATES = ["log_n_cells", "in_pert_graph"]


def run_probe(X, y, feature_names, seed=0):
    """分层 5-fold CV 的 L2 logistic 门控探针。

    X: (n, d) 特征（本函数内做标准化，fold 内 fit 防泄漏）
    y: (n,) 二分类标签
    返回 dict：auroc_mean/std/per_fold、pooled 决策分与 ROC 坐标、
    单特征 AUROC 表。
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.model_selection import (StratifiedKFold, cross_val_predict,
                                         cross_validate)
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    # 注：sklearn>=1.8 中 penalty="l2" 已弃用（lbfgs 默认即 L2），不显式传参
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=seed))

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    scores = cross_validate(clf, X, y, cv=cv, scoring="roc_auc")
    aurocs = scores["test_score"]

    pooled_score = cross_val_predict(clf, X, y, cv=cv,
                                     method="decision_function")
    fpr, tpr, _ = roc_curve(y, pooled_score)
    pooled_auroc = roc_auc_score(y, pooled_score)

    # 单特征消融：每个特征单独跑同一 CV（标准化后）
    single = {}
    for j, name in enumerate(feature_names):
        s = cross_validate(make_pipeline(
            StandardScaler(),
            LogisticRegression(C=1.0, max_iter=2000, random_state=seed)),
            X[:, [j]], y, cv=cv, scoring="roc_auc")["test_score"]
        single[name] = {"auroc_mean": float(s.mean()),
                        "auroc_std": float(s.std())}

    return {
        "auroc_mean": float(aurocs.mean()),
        "auroc_std": float(np.std(aurocs, ddof=1)),
        "auroc_per_fold": [float(a) for a in aurocs],
        "pooled_auroc": float(pooled_auroc),
        "roc_fpr": fpr.tolist(), "roc_tpr": tpr.tolist(),
        "single_feature": single,
    }


if __name__ == "__main__":
    # 冒烟测试：合成数据，真信号特征应显著 > 0.5
    rng = np.random.default_rng(0)
    n = 200
    sig = rng.normal(size=n)
    X = np.column_stack([sig, rng.normal(size=(n, 3))])
    y = (sig + 0.8 * rng.normal(size=n) > 0).astype(int)
    r = run_probe(X, y, ["sig", "x1", "x2", "x3"])
    print(f"[smoke] auroc_mean={r['auroc_mean']:.3f} "
          f"(应 > 0.7), best single={max(r['single_feature'], key=lambda k: r['single_feature'][k]['auroc_mean'])}")
