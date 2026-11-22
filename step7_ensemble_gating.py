"""Day5 Step7（补全检验）：3 种子 ensemble 不确定性 → H2b 档案检验 + 门控探针对照。

两件事：
1. **H2b（预注册判据，守恒线遗留）**：corr(违例分数, ensemble 不确定性)。
   守恒约束路线已按预注册判阴性关闭，本检验**只为档案完整性**，
   结论标注"路线关闭后的补全检验，不改变已归档结论"。
2. **门控探针集成对照**：兑现 gating_verdict.md 的承诺——种子 1/2 训完后
   用 ensemble 方差补全 uncertainty，与 step6 纯廉价特征探针（AUROC 0.6556）对照。

ensemble 不确定性公式（docstring 留痕，口径与 step6 一致）：
    对每个留出扰动 i 和种子 s ∈ {0,1,2}：
        delta_{s,i} = pred_expr_{s,i} − ctrl_train_mean      （共同 ctrl 基线）
    ens_var_i     = mean_g var_s(delta_{s,i,g})              （ddof=1，逐基因
                  跨种子方差再对基因取均值）
    特征用 log_ens_var = log(ens_var_i)                      （方差右偏，取对数）

判据（照 preregistration.md，原文照录）：
    H2b: corr(违例分数, ensemble 方差) < 0.3 → 通过（违例提供独立信号）
         ≥ 0.3 → 信号冗余（转仪表盘路线，本阶段结论照常归档）
    0.3 为判定阈值非对标值（依据见 preregistration.md 阈值依据条款）。
    门控探针判据沿用 step6 定义性阈值：5-fold CV AUROC > 0.6。

用法：被 step7_run.py 调用；也可
    python step7_ensemble_gating.py   # 合成数据冒烟测试
"""
import numpy as np
import pandas as pd

SEEDS = (0, 1, 2)


def ensemble_uncertainty(pred_mats, ctrl_mean):
    """3 种子预测谱 → 逐扰动 ensemble 不确定性（见模块 docstring 公式）。

    pred_mats: list of (n_pert, n_genes) float32/float64，种子顺序对应 SEEDS，
               且行顺序（扰动）必须一致（调用方校验）。
    ctrl_mean: (n_genes,) ctrl 训练均值。
    返回 (n_pert,) ens_var（跨种子逐基因方差的基因均值，ddof=1）。
    """
    deltas = np.stack([p - ctrl_mean.astype(np.float32)
                       for p in pred_mats])          # (3, n_pert, n_genes)
    ens_var = deltas.var(axis=0, ddof=1).mean(axis=1)   # (n_pert,)
    return ens_var


def run_h2b(violation, uncertainty):
    """H2b：违例分数 × ensemble 不确定性的 Pearson r（定义性判据见 docstring）。"""
    from scipy import stats
    r, p = stats.pearsonr(violation, uncertainty)
    return {"r": float(r), "p": float(p), "n": int(len(violation))}


def corr_features(a, b):
    from scipy import stats
    r, p = stats.pearsonr(a, b)
    return {"r": float(r), "p": float(p)}


if __name__ == "__main__":
    # 冒烟测试：合成 3 种子预测谱 + 违例/标签
    rng = np.random.default_rng(0)
    n_pert, n_genes = 100, 200
    base = rng.normal(0, 1, (n_pert, n_genes)).astype(np.float32)
    ctrl = np.zeros(n_genes, dtype=np.float32)
    # 种子 1/2 = 种子 0 + 异方差噪声：扰动 0..49 噪声大 → 高 ens_var
    noise_scale = np.where(np.arange(n_pert) < 50, 2.0, 0.2)[:, None]
    m1 = base + rng.normal(0, noise_scale, base.shape).astype(np.float32)
    m2 = base + rng.normal(0, noise_scale, base.shape).astype(np.float32)
    ev = ensemble_uncertainty([base, m1, m2], ctrl)
    expect_high, expect_low = ev[:50].mean(), ev[50:].mean()
    assert expect_high > 10 * expect_low, \
        f"ens_var 未分离高/低噪声组: {expect_high:.3f} vs {expect_low:.3f}"
    print(f"[smoke] ens_var 高噪声组={expect_high:.3f} >> 低噪声组={expect_low:.3f} OK")

    # 违例与 ens_var 正相关 → H2b r 应 > 0
    viol = 100 + 50 * (ev - ev.min()) / np.ptp(ev) + rng.normal(0, 5, n_pert)
    h2b = run_h2b(viol, ev)
    assert h2b["r"] > 0.5, f"H2b 冒烟 r={h2b['r']:.3f} 应 > 0.5"
    print(f"[smoke] H2b r={h2b['r']:.3f} (p={h2b['p']:.3g}) OK")

    # 门控探针流程：log_ens_var 作特征应能恢复标签
    from step6_gating_probe import run_probe
    y = (np.arange(n_pert) < 50).astype(int)        # 高噪声 = "答错"
    log_ev = np.log(ev)
    r = run_probe(log_ev[:, None], y, ["log_ens_var"])
    assert r["auroc_mean"] > 0.9, f"探针冒烟 AUROC={r['auroc_mean']:.3f} 应 > 0.9"
    print(f"[smoke] 探针 AUROC={r['auroc_mean']:.3f} OK")
