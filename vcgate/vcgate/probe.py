"""门控探针：L2 logistic 回归 + 分层 5-fold CV（预注册口径的发布版）。

与论文/仓库 step6_gating_probe.py 设置完全一致：
StandardScaler + LogisticRegression(C=1.0)，分层 5-fold（shuffle, seed=0）。
"""
import numpy as np

# 预注册判据线（项目内约定阈值，非外部标准；见 preregistration.md）
GATE_CRITERION = 0.6


def _pipeline(seed=0):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=seed))


def oof_gate_scores(X, y, n_splits=5, seed=0):
    """pooled out-of-fold 门控分数（高 = 更可能答错），无泄漏（fold 内 fit）。"""
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return cross_val_predict(_pipeline(seed), np.asarray(X), np.asarray(y),
                             cv=cv, method="decision_function")


def cv_gate_auroc(X, y, n_splits=5, seed=0):
    """逐 fold AUROC 均值 ± std（预注册报告口径）。"""
    from sklearn.model_selection import StratifiedKFold, cross_validate
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    sc = cross_validate(_pipeline(seed), np.asarray(X), np.asarray(y),
                        cv=cv, scoring="roc_auc")["test_score"]
    return float(sc.mean()), float(sc.std(ddof=1))


def fit_gate(X, y, seed=0):
    """在全数据上拟合最终门控模型（部署形态：fit 一次，score 任意新查询）。"""
    return _pipeline(seed).fit(np.asarray(X), np.asarray(y))


def gate_auroc(y, scores):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(np.asarray(y), np.asarray(scores)))


def median_split_labels(error):
    """定义性标签：error 高于中位数 = 1（'答错' = 落在后一半）。"""
    error = np.asarray(error, dtype=float)
    return (error > np.median(error)).astype(int)
