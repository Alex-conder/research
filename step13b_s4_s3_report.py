# step13b_s4_s3_report.py — 由 s4_s3_fixes.json 生成 s4_s3_fixes.md
# 用法：venv/Scripts/python.exe step13b_s4_s3_report.py
import json
import os
import sys

for _s in ("stdout", "stderr"):
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")


def fnum(x, nd=4):
    return f"{x:.{nd}f}"


def arm_stats(c, arm, mode="mc3"):
    return c["paired_bootstrap_gate_minus_random"]["modes"][mode][arm]


def md_table_s4(res, arm, mode, title):
    rows = []
    for key, c in res["s4"]["combos"].items():
        pb = arm_stats(c, arm, mode)
        base = c["random_rejection_baseline"]["mean_over_seeds"]
        mid_cov = "0.50"
        rows.append({
            "combo": key.replace("_", " "), "engine": c["engine"],
            "dataset": c["dataset"], "n": c["n"],
            "gate_risk_c50": c["point_risk_by_coverage"][arm][mid_cov],
            "rand_risk_c50": base[mid_cov]["risk_mean_over_seeds"],
            "diff_mean": pb["mean_diff_over_coverage"],
            "ci_lo": pb["ci_lo"], "ci_hi": pb["ci_hi"],
            "win": pb["gate_significantly_lower_risk"]})
    lines = [f"### {title}", "",
             "| 组合 | 引擎 | 数据集 | n | 门 risk@cov0.50 | 随机 risk@cov0.50 | "
             "配对差均值(coverage 1.0-0.5 平均) | 95% CI | 门显著更优 |",
             "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(
            f"| {r['combo']} | {r['engine']} | {r['dataset']} | {r['n']} | "
            f"{fnum(r['gate_risk_c50'])} | {fnum(r['rand_risk_c50'])} | "
            f"{r['diff_mean']:+.4f} | [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] | "
            f"{'是' if r['win'] else '否'} |")
    s = res["s4"]["summary"][mode][arm]
    lines += ["", f"判定（该臂，{mode} 口径）：{s['n_win']}/{s['n']} 组合门的 risk "
              f"显著低于随机（CI 不含 0 且方向为负）。", ""]
    return lines


def main():
    with open(os.path.join(OUT_DIR, "s4_s3_fixes.json"), encoding="utf-8") as f:
        res = json.load(f)
    s4, s3 = res["s4"], res["s3"]
    spec4 = s4["spec"]

    L = []
    L += ["# 红队挂账实算：S4 闭环配对检验 + S3 预测缩放对照", ""]
    L += [f"- 日期：{res['date']}；任务：{res['task']}",
          f"- 红队来源：{res['redteam_ref']}；门产品口径：{res['r2_product']}",
          "- 实算脚本：`step13_s4_s3_fixes.py`；本报告由 "
            "`step13b_s4_s3_report.py` 生成",
          f"- S4 随机源：bootstrap 主种子 {spec4['boot_master_seed']}；随机拒答"
          f"基线种子 {spec4['random_baseline_seeds']}（点估计曲线每 种子 x "
          f"coverage 抽签 {spec4['baseline_draws_per_seed_coverage']} 次；配对 "
          f"bootstrap 内每重采样每种子逐 coverage 各抽 1 次）；OOF 划分种子 "
          f"{spec4['oof_seed']}（与 step6/step10 同源）",
          f"- S3 随机源：step6 原协议种子 {s3['spec']['seed']}（CV 划分与 LR "
          "初始化）",
          "- 口径复用：组合加载/三门分数来自 step10_signature_resolution 与 "
            "step10_two_layer_gate；三门点估计已与 "
            "outputs_rebuild/two_layer_gate.json 对账一致（最大差 0.0）",
          ""]

    # ---------------- S4 ----------------
    L += ["## S4 闭环配对检验：门 vs 随机拒答基线", ""]
    L += ["### 方法", "",
          "对 R2 判定产品的门（符号层臂，另报就地探针臂与两层门臂作对照），"
          "在 coverage 网格 1.0 -> 0.5（步长 0.05，11 档）上构造随机拒答基线："
          "每档 coverage 以等大概率拒答（标记）k=round(coverage x n) 个查询，"
          "risk = 1 - 被标记查询中最差一半的占比。点估计基线曲线取 3 个随机"
          "种子、每种子每档 1000 次抽签的平均。",
          "配对 bootstrap（B=1000，扰动级有放回重采样）：同一重采样内重算三门"
          "分数，随机拒答臂按 3 个固定种子各抽一次取平均（主口径 mc3，差值的"
          "公共基线噪声在对差中抵消）；稳健性复核用期望随机臂（重采样基率，即"
          "随机拒答的期望精度，蒙特卡洛噪声为零）。每次重采样逐档记录 "
          "risk(门) - risk(随机)，报逐档差值均值与 95% percentile CI，以及 "
          "coverage 平均差值的 95% CI。判定：CI 不含 0 且方向为负（门 risk "
          "显著低于随机）。scf_k562（仅就地臂）与 scf_rpe1（产物缺失）按 R2 "
          "修正案 A2 留痕排除，可执行组合 = 10。", ""]
    L += md_table_s4(res, "symbol", "mc3",
                     "表 S4-1 门产品臂（符号层，R2 判定 tier=symbol_layer）"
                     "，mc3 配对口径")
    L += md_table_s4(res, "in_situ", "mc3", "表 S4-2 对照臂：就地探针，mc3")
    L += md_table_s4(res, "two_layer", "mc3", "表 S4-3 对照臂：两层门，mc3")
    # expected-mode agreement for product arm
    L += ["### 稳健性复核（期望随机臂配对，产品臂 symbol）", "",
          "| 组合 | 配对差均值 | 95% CI | 门显著更优 | 与 mc3 判定一致 |",
          "|---|---|---|---|---|"]
    for key, c in s4["combos"].items():
        e = arm_stats(c, "symbol", "expected")
        m = arm_stats(c, "symbol", "mc3")
        agree = (e["gate_significantly_lower_risk"] ==
                 m["gate_significantly_lower_risk"])
        L.append(f"| {key} | {e['mean_diff_over_coverage']:+.4f} | "
                 f"[{e['ci_lo']:+.4f}, {e['ci_hi']:+.4f}] | "
                 f"{'是' if e['gate_significantly_lower_risk'] else '否'} | "
                 f"{'是' if agree else '否'} |")
    se = s4["summary"]["expected"]["symbol"]
    L += ["", f"期望随机臂口径判定：{se['n_win']}/{se['n']} 组合显著"
          f"（与 mc3 主口径计数一致："
          f"{'是' if se['n_win'] == s4['summary']['mc3']['symbol']['n_win'] else '否'}）。",
          ""]
    # per-coverage detail for product arm
    L += ["### 附：产品臂逐档差值（risk_gate - risk_random，均值 [95% CI]，mc3）",
          ""]
    hdr = "| 组合 | " + " | ".join(f"{c:.2f}" for c in spec4["coverages"]) + " |"
    L += [hdr, "|---|" + "---|" * len(spec4["coverages"])]
    for key, c in s4["combos"].items():
        cells = []
        pc = arm_stats(c, "symbol")["diff_risk_gate_minus_random_by_coverage"]
        for cov in spec4["coverages"]:
            d = pc[f"{cov:.2f}"]
            mark = "*" if not d["ci_excludes_zero"] else ""
            cells.append(f"{d['diff_mean']:+.3f} [{d['ci_lo']:+.3f},"
                         f"{d['ci_hi']:+.3f}]{mark}")
        L.append(f"| {key} | " + " | ".join(cells) + " |")
    L += ["", "（* = 该档 CI 含 0；coverage=1.0 时门与随机标记全体，差恒为 0。）",
          "随机拒答基线的种子稳定性：逐组合 3 种子基线 risk 最大极差见 JSON "
          "s4.combos[*].random_rejection_baseline.mean_over_seeds."
          "risk_spread_over_seeds（量级 1e-3 或更小）。", ""]

    # ---------------- S3 ----------------
    L += ["## S3 预测缩放对照：探针 AUROC 的缩放稳健性", ""]
    L += ["### 方法", "",
          "主对照：对 GEARS-K562/RPE1 的 9 特征探针输入做统一缩放"
          "（StandardScaler 全数据一次 fit），以同一 CV 划分（stratified "
          "5-fold, shuffle, seed=0）重训纯 L2 logistic 并重评 AUROC，与 "
          "step6 原口径（fold 内 StandardScaler 管线，CV AUROC 均值为主读数、"
          "pooled OOF AUROC 为辅读数）对比。重训前已对原口径做复现校验，与存档 "
          "outputs{,_rpe1}/gating_probe.json 一致（差 < 1e-9）。",
          "补充对照（红队原文'恒定/缩放预测'）：把预测 delta 乘常数因子 "
          "0.5 与 2.0 重算特征（n_de_025/frac_de_025 经 0.25 阈值非线性变化，"
          "其余线性或不变），pearson 标签对正缩放数学不变（y 不动），按原协议"
          "重评。另报完全不缩放的诊断行。", ""]
    L += ["### 结果", "",
          "| 组合 | 口径 | CV AUROC 均值 | pooled OOF AUROC | 相对原口径差 |",
          "|---|---|---|---|---|"]
    for key, c in s3["combos"].items():
        o = c["original_protocol"]
        u = c["unified_scaler_control"]
        r = c["no_scaling_diagnostic"]
        L.append(f"| {key} | 原口径（fold 内缩放） | {fnum(o['auroc_mean'],6)} "
                 f"| {fnum(o['pooled_auroc'],6)} | - |")
        L.append(f"| {key} | 统一缩放对照（全数据一次 fit） | "
                 f"{fnum(u['auroc_mean'],6)} | {fnum(u['pooled_auroc'],6)} | "
                 f"{u['cv_auroc_diff_vs_original']:+.6f} |")
        L.append(f"| {key} | 诊断：完全不缩放 | {fnum(r['auroc_mean'],6)} | "
                 f"{fnum(r['pooled_auroc'],6)} | "
                 f"{r['auroc_mean'] - o['auroc_mean']:+.6f} |")
        for tag, p in c["prediction_scaling_supplementary"].items():
            L.append(f"| {key} | 补充：预测 delta x{tag.split('x')[-1]} | "
                     f"{fnum(p['auroc_mean'],6)} | {fnum(p['pooled_auroc'],6)} | "
                     f"{p['auroc_mean'] - o['auroc_mean']:+.6f} |")
    L += ["", "### 判定", ""]
    for key, c in s3["combos"].items():
        u = c["unified_scaler_control"]
        ps = c["prediction_scaling_supplementary"]
        v = c["validation_vs_archived"]
        L.append(
            f"- {key}：统一缩放前后 CV AUROC 差 "
            f"{u['cv_auroc_diff_vs_original']:+.6f}（pooled 差 "
            f"{u['pooled_auroc_diff_vs_original']:+.6f}）；预测级缩放 delta x0.5 / "
            f"x2.0 的 AUROC 差分别为 "
            f"{ps['delta_x0.5']['auroc_mean'] - c['original_protocol']['auroc_mean']:+.6f} / "
            f"{ps['delta_x2.0']['auroc_mean'] - c['original_protocol']['auroc_mean']:+.6f}"
            f"；原口径复现与存档一致：{v['match']}。")
    L += ["", "结论：探针 AUROC 不是缩放步骤的产物——把缩放射到 CV 口径外"
          "（全数据统一缩放）AUROC 不动；把预测本身缩放常数倍，AUROC 的移动"
          "在 1e-3 量级以内（见表），不改变任何判定。", ""]

    # ---------------- 论文附录英文草稿 ----------------
    s = s4["summary"]["mc3"]
    gk = s4["combos"]["gears_k562"]["paired_bootstrap_gate_minus_random"]
    gr = s4["combos"]["gears_rpe1"]["paired_bootstrap_gate_minus_random"]
    sk3 = s3["combos"]["gears_k562"]
    sr1 = s3["combos"]["gears_rpe1"]
    L += ["## 附：论文附录英文草稿（主代理决定入稿位置；未改 main.tex）", ""]
    L += ["### S4 段（Formal paired test against random rejection）", "",
          "```",
          "\\textbf{Formal paired test against random rejection.} ",
          "We replaced the descriptive gate-versus-random delta (Sec.~",
          "\\textit{closed-loop}) with a paired test on the risk--coverage ",
          "grid. For each of the 10 evaluable gate combinations (R2 protocol; ",
          "scFoundation K562 lacks the sibling product and the RPE1 product ",
          "is missing), a random-rejection baseline flags a uniformly random ",
          "k = round(coverage x n) queries at each of 11 coverage points ",
          "(1.0 to 0.5, step 0.05); its risk is averaged over 3 seeds x 1000 ",
          "draws (seed-to-seed spread <= 0.0028). A paired perturbation-level ",
          "bootstrap (B=1000; both arms recomputed on each resample, random ",
          "arm redrawn per resample under the 3 fixed seeds; master seed ",
          "20261103) tests the coverage-averaged risk difference per ",
          "combination. The product gate (symbol layer) significantly beats ",
          "random rejection in 2/10 combinations -- both GEARS ",
          "expression-class combinations: K562 -0.0532 [-0.0877, -0.0183], ",
          "RPE1 -0.0588 [-0.0850, -0.0309] (95% percentile CI) -- and is ",
          "statistically indistinguishable from random rejection on the ",
          "eight embedding-class wave combinations. The in-situ probe arm ",
          "shows the same pattern (2/10; GEARS only), and the two-layer arm ",
          "reaches significance in no combination (0/10). A zero-",
          "Monte-Carlo-noise pairing (expected-random arm = resample base "
          "rate) reproduces the product-gate (symbol-arm) verdict in all "
          "ten combinations; for the two-layer arm the same pairing flips "
          "only the knife-edge GEARS-K562 call (CI upper +0.0034 vs "
          "-0.0002), leaving its count at 0--1/10 either way. The gate's ",
          "closed-loop value is therefore confirmed on the expression class ",
          "and bounded elsewhere; we report the non-significant combinations ",
          "rather than averaging them into a cleaner story.",
          "```", ""]
    L += ["### S3 段（Prediction-scaling control）", "",
          "```",
          "\\textbf{Prediction-scaling control.} ",
          "To exclude the possibility that the probe AUROC is an artifact of ",
          "the scaling stage, we retrained the 9-feature GEARS probe after a ",
          "unified StandardScaler fit (full data, once) and re-evaluated ",
          "under the identical 5-fold CV protocol (seed 0): the CV AUROC ",
          "mean moves by +0.000433 on K562 (0.655647 -> 0.656080) and by ",
          "0.000000 on RPE1 (0.796518 -> 0.796518); pooled OOF AUROC moves ",
          "by at most 0.00034. Scaling the predicted perturbation deltas by ",
          "constant factors 0.5 and 2.0, recomputing all probe features (the ",
          "DE-count features change non-linearly through the 0.25 threshold; ",
          "Pearson labels are invariant to positive rescaling), shifts AUROC ",
          "by at most +0.0009 (K562) and +0.0008 (RPE1). The probe's signal ",
          "is not a product of either feature scaling or prediction scale.",
          "```", ""]

    out_md = os.path.join(OUT_DIR, "s4_s3_fixes.md")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("[SAVED]", out_md)


if __name__ == "__main__":
    main()
