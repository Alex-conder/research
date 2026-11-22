"""Day3 Step4: H2a 判据检验（含 null 与偏相关），出图 + verdict.md。

判据（preregistration.md + 2026-10-12 增补版流程）：
- H2a:      corr(violation, error) 的 Pearson r 与 p；成立线 r > 0.3 且 p < 0.05
- H2a-null: 200 次置换 error 标签的 r 零分布；观测 r 超 95 分位为通过置换检验
- H2a-partial: 控制 strength 的偏相关 + t 检验（df = n-3）
- H2b:      本轮跳过（uncertainty 全 NaN，待种子 1/2），verdict.md 注明

输出：outputs/scatter_h2a.png, outputs/null_distribution.png, outputs/verdict.md
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

OUT = "outputs"
N_PERM = 200
rng = np.random.default_rng(0)

# ---- 1. 合并 violation 与 error/strength ----
pred = pd.read_csv(f"{OUT}/predictions.csv")          # perturbation, pearson_delta, strength, ...
viol = pd.read_csv(f"{OUT}/violation_score.csv")      # perturbation, violation
df = pred.merge(viol, on="perturbation", how="inner").dropna(
    subset=["pearson_delta", "violation", "strength"])
df["error"] = 1.0 - df["pearson_delta"]
n = len(df)
print(f"[step4] n={n}（合并后）")
assert n >= 20, f"样本太少 n={n}"

v, e, s = df.violation.values, df.error.values, df.strength.values

# ---- 2. H2a ----
r_obs, p_obs = stats.pearsonr(v, e)
print(f"[step4] H2a: r={r_obs:.4f}, p={p_obs:.4g}")

# ---- 3. H2a-null：置换 ----
r_null = np.array([stats.pearsonr(v, e[rng.permutation(n)])[0]
                   for _ in range(N_PERM)])
null_p95 = np.percentile(r_null, 95)
emp_p = (1 + np.sum(r_null >= r_obs)) / (1 + N_PERM)
print(f"[step4] H2a-null: r_obs={r_obs:.4f} vs null 95%={null_p95:.4f}, "
      f"empirical p={emp_p:.4f}")

# ---- 4. H2a-partial：控制 strength ----
def resid(y, x):
    b = np.polyfit(x, y, 1)
    return y - np.polyval(b, x)

r_partial = stats.pearsonr(resid(v, s), resid(e, s))[0]
t_stat = r_partial * np.sqrt((n - 3) / (1 - r_partial ** 2))
p_partial = 2 * stats.t.sf(abs(t_stat), n - 3)
print(f"[step4] H2a-partial|r(strength): r={r_partial:.4f}, t={t_stat:.3f}, "
      f"p={p_partial:.4g}")

# ---- 5. 图 ----
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.scatter(v, e, alpha=0.6)
b = np.polyfit(v, e, 1)
xs = np.linspace(v.min(), v.max(), 50)
ax.plot(xs, np.polyval(b, xs), "r--", lw=1)
ax.set_xlabel("conservation violation ||S*vhat||_2 (unnormalized)")
ax.set_ylabel("error = 1 - pearson_delta (top50 DEG)")
ax.set_title(f"H2a: r={r_obs:.3f}, p={p_obs:.3g}, n={n}")
fig.tight_layout()
fig.savefig(f"{OUT}/scatter_h2a.png", dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.hist(r_null, bins=30, alpha=0.7, label="null (200 perms)")
ax.axvline(r_obs, color="r", lw=2, label=f"observed r={r_obs:.3f}")
ax.axvline(null_p95, color="k", ls="--", label=f"null 95%={null_p95:.3f}")
ax.set_xlabel("Pearson r (violation, permuted error)")
ax.set_ylabel("count")
ax.legend()
fig.tight_layout()
fig.savefig(f"{OUT}/null_distribution.png", dpi=150)
plt.close(fig)

# ---- 6. verdict.md ----
h2a_pass = (r_obs > 0.3) and (p_obs < 0.05)
h2a_clear_fail = (r_obs < 0.1) or (p_obs >= 0.05)
setup = json.load(open(f"{OUT}/proxy_setup.json"))
lines = []
lines.append("# verdict.md — 第一阶段 H2a 检验结论（Day3，2026-10-13）\n")
lines.append("## 结果\n")
lines.append(f"- n = {n} 个留出测试扰动（seed 0，单种子；3 种子为预注册锁定设置，本轮时间盒偏离）")
lines.append(f"- **H2a：r = {r_obs:.4f}，p = {p_obs:.4g}**"
             f"（判据：r > 0.3 且 p < 0.05）→ {'**成立**' if h2a_pass else '**未成立**'}")
lines.append(f"- H2a-null：置换 200 次，观测 r 的 empirical p = {emp_p:.4f}"
             f"（null 95 分位 = {null_p95:.4f}）→ "
             f"{'通过' if r_obs > null_p95 else '未通过'}置换检验")
lines.append(f"- H2a-partial（控制 strength）：r = {r_partial:.4f}，"
             f"t({n-3}) = {t_stat:.3f}，p = {p_partial:.4g}")
lines.append("- **H2b：本轮跳过** —— uncertainty 全 NaN（种子 1/2 未训），"
             "待 ensemble 补齐后检验。\n")
lines.append("## 三条件分流\n")
if h2a_pass:
    lines.append("1. H2a 成立 → 进入 H2b（需先补种子 1/2 的 ensemble 不确定性）。")
    lines.append("2. H2b 通过（corr(violation, uncertainty) < 0.3）→ 启动 H1（二阶段）。")
    lines.append("3. H2b 不通过（≥ 0.3）→ 转仪表盘路线，结论照常归档。")
elif h2a_clear_fail:
    lines.append("1. H2a 不成立（r < 0.1 或 p ≥ 0.05）→ 按预注册**归档阴性**，"
                 "做归因分析（映射距离 / 守恒律不适用性），方案终止。")
    lines.append("2. 归档格式与阳性结论相同。\n")
    lines.append("## 归因线索（阴性时参考）\n")
    lines.append(f"- strength（扰动效应量）与 error 的相关："
                 f"r = {stats.pearsonr(s, e)[0]:.4f}")
    lines.append(f"- strength 与 violation 的相关："
                 f"r = {stats.pearsonr(s, v)[0]:.4f}")
    lines.append("- 代理通量为表达均值粗近似，一阶段噪声预期内（预注册边界声明）。")
else:
    lines.append("1. H2a 落在灰色区间（0.1 ≤ r ≤ 0.3）→ 不下结论；"
                 "补种子 1/2 + 完整 Replogle 数据后重测。")
    lines.append("2. H2b 待 uncertainty 补齐后并列检验。")
    lines.append("3. 若重测仍灰色 → 按预注册做归因分析并归档。\n")
lines.append("## 协议偏离记录（按偏离条款并列报告）\n")
lines.append("1. OOD 划分：scPerturBench 官方留出名单未取得，用 seed=0 自行留出 20%（D2-0）。")
lines.append("2. 数据：GEARS 过滤版 replogle_k562_essential 子集，非完整 Replogle K562。")
lines.append("3. pearson_delta 差异基因数：top-50 为自定（GEARS 内部 top-20，"
             "scPerturBench 未确证）；全基因口径 pearson_delta_all 已在 CSV 并列。")
lines.append("4. uncertainty 全 NaN：3 种子 ensemble 未跑（时间盒偏离）。")
lines.append(f"5. 高变子集解释：表达方差基于 217 个**预测**谱（非原始数据）；"
             f"保留反应 {setup['n_rxns_kept']}/{setup['n_rxns_total']}。")
lines.append("\n## 复现\n")
lines.append("- 脚本：step1_run.py / step2_run.py / step3_run.py / step4_run.py")
lines.append(f"- 中间产物：outputs/predictions.csv, violation_score.csv, "
             f"pred_expr_seed0.npz, proxy_setup.json")
lines.append("- 日志：setup_log.md（Day3 小节）")
open(f"{OUT}/verdict.md", "w", encoding="utf-8").write("\n".join(lines))
print(f"[step4] -> {OUT}/scatter_h2a.png, null_distribution.png, verdict.md")
print(f"[step4] 分流: h2a_pass={h2a_pass}, clear_fail={h2a_clear_fail}")
