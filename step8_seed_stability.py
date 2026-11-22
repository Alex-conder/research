"""Day5 Step8：门控探针的种子稳定性检查（为种子 2 落地收尾）。

背景：step6 的门控信号（seed0 上 5-fold CV AUROC = 0.6556）只在单种子预测上
验证过。本模块对 seed ∈ {0,1,2} 各自重跑 step6 完整流程
（特征构造 → L2 logistic 探针 → 分层 5-fold CV AUROC），检验信号对模型
初始化是否稳健。

稳定性判据（三条，层级：存在性 > 方向 > 量级；量级不过不否定存在性）：
- **信号存在**：至少两粒 AUROC > 0.6 且三粒均值 > 0.6；不过 → 归档阴性。
- **方向稳健**：三种子 AUROC 均 > 0.6 且 best 单特征一致；不过 → 归因。
- **量级稳健**：三种子 AUROC 极差 < 0.05（依据：0.05 ≈ step6 观测 fold 间
  std 0.048；n=217 小样本下 fold 级 std 约 0.05-0.06，极差落在该区间时
  优先归因小样本而非信号不稳）。不过 → 结论降级为"精确量级需更大样本确认"。

组合结论：三✅ → "信号稳健存在"；方向✅量级❌ → "信号存在，精确量级需更
大样本确认"；存在❌ → "信号不存在，归档阴性"；方向❌ → "信号不稳健，归因"。

不变量检查（[修复 D5] 的预期）：三个种子的留出扰动集合必须完全一致
（predictions_seed{s}.csv 第一列全等），不等即报错退出——这本身也是检查项。

复用纪律：特征构造与探针直接 import step6_gating_probe 的
compute_features / run_probe / SPECTRUM_FEATURES，不复制代码。

用法：被 step8_run.py 调用；也可
    python step8_seed_stability.py   # 合成数据冒烟测试
"""
import numpy as np
import pandas as pd

from step6_gating_probe import SPECTRUM_FEATURES, compute_features, run_probe

AUROC_EXIST = 0.6              # 信号存在/方向的阈值（沿用 step6 定义性阈值）
STABILITY_RANGE = 0.05          # 量级阈值：三种子 AUROC 极差 < 0.05 为稳健


def gating_for_seed(pred_expr, perts, pred_genes, ctrl_mean, pred_df,
                    probe_seed=0):
    """单种子完整门控流程（与 step6_run.py 第 2-4 节同口径）。

    pred_expr/parts/pred_genes: 该种子的 npz 内容
    ctrl_mean: (n_genes,) ctrl 训练均值（按 pred_genes 顺序）
    pred_df:   该种子的 predictions CSV（index=perturbation，与 perts 对齐）
    返回 dict：特征名、CV AUROC 均值/std/逐 fold/pooled、单特征表、最好单特征。
    """
    pred_df = pred_df.loc[list(perts)]
    feat = compute_features(pred_expr, perts, pred_genes, ctrl_mean)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    # in_pert_graph 全为 True（常数），step6 已留痕剔除，此处同口径不引入
    feature_names = SPECTRUM_FEATURES + ["log_n_cells"]
    X = feat[feature_names].values.astype(np.float64)
    assert not np.isnan(X).any(), "特征矩阵含 NaN"

    error = 1.0 - pred_df["pearson_delta"].values
    y = (error > np.median(error)).astype(int)

    res = run_probe(X, y, feature_names, seed=probe_seed)
    best = max(res["single_feature"],
               key=lambda k: res["single_feature"][k]["auroc_mean"])
    return {
        "n_perts": int(len(y)), "n_positive": int(y.sum()),
        "error_median": float(np.median(error)),
        "features": feature_names,
        "cv_auroc_mean": res["auroc_mean"],
        "cv_auroc_std": res["auroc_std"],
        "cv_auroc_per_fold": res["auroc_per_fold"],
        "pooled_cv_auroc": res["pooled_auroc"],
        "best_single_feature": best,
        "single_feature_auroc": res["single_feature"],
    }


def stability_summary(per_seed):
    """三种子结果 → 三判据逐项判定 + 组合结论。

    判据（层级：存在性 > 方向 > 量级；量级不过不否定存在性）：
    - 方向稳健 direction：全部种子 AUROC > 0.6 且 best 单特征一致。
      不过 → "信号不稳健，归因"。
    - 量级稳健 magnitude：AUROC 极差 < 0.05。阈值依据：0.05 ≈ step6 观测
      fold 间 std（0.048）。归因注释：n=217 小样本下 fold 级 std 约 0.05-0.06，
      极差落在该区间时优先归因小样本而非信号不稳。不过 → 不否定存在性，
      结论降级为"精确量级需更大样本确认"。
    - 信号存在 existence：至少两粒 AUROC > 0.6 且三粒均值 > 0.6。
      不过 → "信号不存在，归档阴性"。

    组合结论映射：
      方向✅量级✅存在✅ → "信号稳健存在"
      方向✅量级❌存在✅ → "信号存在，精确量级需更大样本确认"
      存在❌            → "信号不存在，归档阴性"
      方向❌            → "信号不稳健，归因"
    """
    aucs = {s: r["cv_auroc_mean"] for s, r in per_seed.items()}
    bests = {s: r["best_single_feature"] for s, r in per_seed.items()}
    rng = max(aucs.values()) - min(aucs.values())
    mean_auc = float(np.mean(list(aucs.values())))

    existence = sum(a > AUROC_EXIST for a in aucs.values()) >= 2 and \
        mean_auc > AUROC_EXIST
    direction = all(a > AUROC_EXIST for a in aucs.values()) and \
        len(set(bests.values())) == 1
    magnitude = rng < STABILITY_RANGE

    if not existence:
        verdict = "信号不存在，归档阴性"
    elif not direction:
        verdict = "信号不稳健，归因"
    elif not magnitude:
        verdict = "信号存在，精确量级需更大样本确认"
    else:
        verdict = "信号稳健存在"

    return {
        "auroc_by_seed": aucs,
        "best_single_feature_by_seed": bests,
        "mean_auroc": mean_auc,
        "range": float(rng),
        "criteria": {
            "existence_at_least_2_above_0.6_and_mean_above_0.6": bool(existence),
            "direction_all_above_0.6_and_same_best_feature": bool(direction),
            "magnitude_range_below_0.05": bool(magnitude),
            "existence_threshold": AUROC_EXIST,
            "magnitude_threshold": STABILITY_RANGE,
            "magnitude_attribution": "极差阈值 0.05 ≈ step6 fold std 0.048；"
                "n=217 小样本下 fold 级 std 约 0.05-0.06，极差落在该区间时"
                "优先归因小样本而非信号不稳；量级不过不否定存在性",
        },
        "verdict": verdict,
    }


if __name__ == "__main__":
    # ---- 端到端冒烟：3 个"种子"携带同一信号 → 各自 AUROC 都应高且极差应小 ----
    rng = np.random.default_rng(0)
    n_pert, n_genes = 200, 500
    sig = rng.normal(size=n_pert)                  # 共享信号：扰动幅度
    base = rng.normal(0, 1, (n_pert, n_genes)).astype(np.float32)
    ctrl = np.zeros(n_genes, dtype=np.float32)
    amp = np.abs(sig)[:, None]
    per_seed = {}
    for s in range(3):
        mat = base + amp * rng.normal(0, 1, base.shape).astype(np.float32) \
            + 0.1 * s * rng.normal(size=base.shape).astype(np.float32)
        # 预测质量随幅度衰减 → 大幅扰动"答错"
        quality = 1.0 / (1.0 + amp.ravel())
        pdf = pd.DataFrame({"pearson_delta": quality,
                            "n_cells": rng.integers(20, 400, n_pert)})
        pdf.index = [f"P{i}" for i in range(n_pert)]
        perts = pdf.index.tolist()
        per_seed[s] = gating_for_seed(mat, perts, [f"G{i}" for i in range(n_genes)],
                                      ctrl, pdf)
        print(f"[smoke] seed{s}: AUROC={per_seed[s]['cv_auroc_mean']:.3f}")
    summ = stability_summary(per_seed)
    assert all(a > 0.7 for a in summ["auroc_by_seed"].values()), \
        "合成信号下单种子 AUROC 应 > 0.7"
    assert summ["criteria"]["magnitude_range_below_0.05"], \
        f"极差 {summ['range']:.4f} 应 < {STABILITY_RANGE}"
    print(f"[smoke] 端到端: range={summ['range']:.4f}, verdict={summ['verdict']}")

    # ---- 组合结论映射单测：best 特征在合成随机数据上可能逐种子翻转，
    # 端到端冒烟不保证 direction=True，四档映射用构造结果逐项验证 ----
    def mk(aucs, best=("delta_max",) * 3):
        return {s: {"cv_auroc_mean": a, "best_single_feature": b}
                for s, (a, b) in enumerate(zip(aucs, best))}

    cases = [
        (mk([0.65, 0.66, 0.67]), "信号稳健存在"),                    # ✅✅✅
        (mk([0.65, 0.69, 0.72]), "信号存在，精确量级需更大样本确认"),  # ✅❌✅
        (mk([0.50, 0.55, 0.58]), "信号不存在，归档阴性"),            # 存在❌
        (mk([0.65, 0.58, 0.67]), "信号不稳健，归因"),                # 方向❌（一粒<0.6，
        #  另两粒>0.6 且均值 0.633>0.6 → 存在✅，按层级归方向档）
        (mk([0.65, 0.66, 0.67], ("delta_max", "delta_l2", "delta_max")),
         "信号不稳健，归因"),                                        # best 不一致
    ]
    for fake, expect in cases:
        got = stability_summary(fake)["verdict"]
        assert got == expect, f"映射错误: {got!r} != {expect!r}"
    print(f"[smoke] 组合结论映射 {len(cases)} 档全部通过")

