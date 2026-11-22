# -*- coding: utf-8 -*-
"""Step 12: 数据覆盖消融的签名对比分析（A/B/C/D）。

输入：outputs_ablation/{A_k562,B_rpe1,C_mixed,D_zeroshot}/<V>_predictions.npz
      + _labels.csv（服务器 train_gears_ablation.py 产物）
分析单元：
  A_k562, B_rpe1, D_zeroshot(rpe1), C_mixed 按 row_dataset 拆成 C_k562 / C_rpe1。
对每个单元复用 step6 口径：
  特征 = step6_gating_probe.compute_features（8 spectrum + log_n_cells 协变量）
  标签 y = (1 - pearson_delta) > median（最差一半）
  探针 = 分层 5-fold CV L2 logistic（step6_gating_probe.run_probe）
  签名 = 就地标准化 logistic 系数（step8b.fit_coef）
判定：
  C 的签名（相对 A / B）消失或显著改变 -> "能改"（签名由训练数据覆盖决定）
  C 的签名与 A/B 同构（主导特征相同、量级相近）-> "不能改"（由模型/任务结构决定）
产物：
  outputs_ablation/signature_comparison.json
  outputs_ablation/signature_comparison.png（300dpi，英文标注）
  outputs_ablation/verdict.md
"""
import json
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)
from step6_gating_probe import (COVARIATES, SPECTRUM_FEATURES,  # noqa: E402
                                compute_features, run_probe)
from step8b_signature_ci import fit_coef  # noqa: E402

OUT = os.path.join(REPO, "outputs_ablation")
ANALYSIS_DIR = OUT

BASELINE_SIGS_JSON = os.path.join(
    REPO, "outputs_crossmodel", "analysis", "error_signatures.json")


def load_unit(npz_path, labels_path, dataset_filter=None, restrict_genes=None):
    """加载一个分析单元，返回 (X, y, feature_names, meta)。

    restrict_genes: 若给出（如 2832 共有基因），把所有单元限定到该基因子空间
    再算特征——n_de/delta_l1 等计数/求和类特征受基因数影响，跨空间直接比
    AUROC/签名会混淆，统一空间后再比（A/B/D 存的是 5000 空间）。
    """
    z = np.load(npz_path, allow_pickle=True)
    pred_expr = z["pred_expr"]
    perts = [str(g) for g in z["pred_perturbations"]]
    pred_genes = [str(g) for g in z["pred_genes"]]
    row_ds = ([str(d) for d in z["row_dataset"]]
              if "row_dataset" in z.files else None)
    if "ctrl_mean_matched" in z.files:
        ctrl_mean = np.asarray(z["ctrl_mean_matched"], dtype=np.float64)
    else:
        ctrl_mean = np.asarray(z["ctrl_mean_full"], dtype=np.float64)
    if ctrl_mean.ndim > 1:
        # 逐行匹配的 ctrl（所有行通常相同）；折叠成 1D 均值向量，
        # 否则行数与 ok 过滤后的 pred_expr 不一致会广播失败（D 踩过）
        with np.errstate(all="ignore"):
            ctrl_mean = np.nanmean(ctrl_mean, axis=0)
        ctrl_mean = np.where(np.isnan(ctrl_mean), 0.0, ctrl_mean)
    if restrict_genes is not None:
        pos = {g: i for i, g in enumerate(pred_genes)}
        idx = [pos[g] for g in restrict_genes if g in pos]
        pred_expr = pred_expr[:, idx]
        ctrl_mean = ctrl_mean[:, idx] if ctrl_mean.ndim > 1 else ctrl_mean[idx]
        pred_genes = [pred_genes[i] for i in idx]
    # 全 NaN 列剔除
    nan_col = np.isnan(pred_expr).all(axis=0)
    if nan_col.any():
        pred_expr = pred_expr[:, ~nan_col]
        pred_genes = [g for g, k in zip(pred_genes, ~nan_col) if k]
        ctrl_mean = ctrl_mean[~nan_col]

    lab = pd.read_csv(labels_path)
    err = 1.0 - lab["pearson_delta"].values
    n_cells = lab["n_cells"].values.astype(float)
    in_graph = lab["in_pert_graph"].values.astype(bool) \
        if "in_pert_graph" in lab.columns else np.ones(len(lab), bool)

    if dataset_filter is not None and row_ds is not None:
        sel = np.array([d == dataset_filter for d in row_ds])
    else:
        sel = np.ones(len(err), bool)

    ok = sel & ~np.isnan(err) & in_graph & ~np.isnan(pred_expr).any(axis=1)
    dropped = int((sel & ~ok).sum())
    pred_expr, perts_arr = pred_expr[ok], np.array(perts)[ok]
    err, n_cells = err[ok], n_cells[ok]
    if len(err) < 20:
        raise SystemExit(f"有效行太少: {npz_path} filter={dataset_filter} n={len(err)}")
    y = (err > np.median(err)).astype(int)

    feat = compute_features(pred_expr, list(perts_arr), pred_genes, ctrl_mean)
    feat["log_n_cells"] = np.log1p(n_cells)
    feature_names = SPECTRUM_FEATURES + [c for c in COVARIATES
                                         if c in feat.columns]
    # 常数列（如 in_pert_graph）剔除
    keep = [f for f in feature_names if feat[f].std() > 0]
    dropped_const = sorted(set(feature_names) - set(keep))
    feature_names = keep
    X = feat[feature_names].values.astype(np.float64)
    assert not np.isnan(X).any()
    meta = {"n_total": int(sel.sum()), "n_used": int(len(y)),
            "n_dropped": dropped, "dropped_const": dropped_const,
            "n_genes": int(len(pred_genes)),
            "pearson_delta_mean": float(np.mean(1.0 - err))}
    return X, y, feature_names, meta


def cos_sim(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else 0.0


def dominant(feature_names, coef, k=3):
    order = np.argsort(-np.abs(coef))
    return [{"feature": feature_names[i],
             "coef": float(coef[i])} for i in order[:k]]


def main():
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    common_path = os.path.join(OUT, "common_genes.json")
    restrict = None
    if os.path.exists(common_path):
        with open(common_path, encoding="utf-8") as f:
            restrict = json.load(f)["common_genes_k562_order"]
        print(f"[cfg] 统一特征空间 -> {len(restrict)} 共有基因")
    units = {}
    specs = [
        ("A_k562", "A_k562", None),
        ("B_rpe1", "B_rpe1", None),
        ("C_k562", "C_mixed", "k562"),
        ("C_rpe1", "C_mixed", "rpe1"),
        ("D_zeroshot_rpe1", "D_zeroshot", None),
    ]
    for unit_name, variant, ds_filter in specs:
        d = os.path.join(OUT, variant)
        npz = os.path.join(d, f"{variant}_predictions.npz")
        lab = os.path.join(d, f"{variant}_labels.csv")
        if not (os.path.exists(npz) and os.path.exists(lab)):
            print(f"[skip] {unit_name}: 产物缺失 ({d})")
            continue
        X, y, fn, meta = load_unit(npz, lab, ds_filter, restrict)
        probe = run_probe(X, y, fn, seed=0)
        coef = fit_coef(X, y)
        units[unit_name] = {
            "variant": variant, "dataset_filter": ds_filter,
            "feature_names": fn, "meta": meta,
            "probe": {k: probe[k] for k in
                      ("auroc_mean", "auroc_std", "pooled_auroc",
                       "single_feature")},
            "signature": {f: float(c) for f, c in zip(fn, coef)},
            "dominant_top3": dominant(fn, coef),
        }
        print(f"[{unit_name}] n={meta['n_used']} "
              f"AUROC={probe['auroc_mean']:.3f} "
              f"top3={[(d['feature'], round(d['coef'],2)) for d in units[unit_name]['dominant_top3']]}")

    if len(units) < 4:
        print(f"[warn] 只有 {len(units)} 个单元可用，对比可能不完整")

    # ---- 签名相似度矩阵（就地系数，同特征序） ----
    names = list(units)
    sigs = np.array([[units[n]["signature"][f] for f in units[n]["feature_names"]]
                     for n in names])
    sim = {}
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            if j <= i:
                continue
            sim[f"{ni}||{nj}"] = cos_sim(sigs[i], sigs[j])

    # ---- 与基线归档签名对比（A vs GEARS(K562), B vs GEARS(RPE1)） ----
    baseline_cmp = {}
    if os.path.exists(BASELINE_SIGS_JSON):
        with open(BASELINE_SIGS_JSON, encoding="utf-8") as f:
            arch = json.load(f)
        want = {"A_k562": ["GEARS(K562)", "gears_k562", "GEARS_K562"],
                "B_rpe1": ["GEARS(RPE1)", "gears_rpe1", "GEARS_RPE1",
                           "gears_rpe1_seed0"]}
        for uname, cands in want.items():
            if uname not in units:
                continue
            hit = next((c for c in cands if c in arch), None)
            if hit is None:
                # 模糊匹配：名字含 gears 且含 k562/rpe1
                tag = "k562" if "k562" in uname else "rpe1"
                hit = next((k for k in arch if "gears" in k.lower()
                            and tag in k.lower()), None)
            if hit is None:
                continue
            bsig = arch[hit]
            fn = units[uname]["feature_names"]
            bv = np.array([bsig.get(f, np.nan) for f in fn])
            lv = np.array([units[uname]["signature"][f] for f in fn])
            okf = ~np.isnan(bv)
            if okf.sum() >= 3:
                baseline_cmp[uname] = {
                    "baseline_key": hit,
                    "cos_sim": cos_sim(lv[okf], bv[okf]),
                    "n_features": int(okf.sum()),
                }

    # ---- C 判定 ----
    verdict = {"question": "误差探针签名由训练数据覆盖决定，还是由模型架构决定？",
               "rule": "C 的签名消失/显著改变 -> 能改；C 与 A/B 同构 -> 不能改"}
    c_keys = [n for n in names if n.startswith("C_")]
    if c_keys:
        obs = {}
        for ck in c_keys:
            ds = "A_k562" if ck == "C_k562" else "B_rpe1"
            pair = f"{ck}||{ds}" if f"{ck}||{ds}" in sim else f"{ds}||{ck}"
            obs[ck] = {
                "vs_same_data_solo": sim.get(pair),
                "dominant_top3": units[ck]["dominant_top3"],
                "single_feature": {
                    k: round(v["auroc_mean"], 3)
                    for k, v in units[ck]["probe"]["single_feature"].items()},
                "auroc": units[ck]["probe"]["auroc_mean"],
            }
        verdict["C_observation"] = obs
        # 判据：主导特征是否保持（A: delta_max 主导；B: log_n_cells 主导）
        def dom_feat(u):
            return units[u]["dominant_top3"][0]["feature"] if u in units else None
        keep = all(dom_feat(ck) == dom_feat("A_k562" if ck == "C_k562" else "B_rpe1")
                   for ck in c_keys if (ck == "C_k562" and "A_k562" in units)
                   or (ck == "C_rpe1" and "B_rpe1" in units))
        verdict["dominant_feature_preserved_in_C"] = bool(keep)
        verdict["conclusion"] = (
            "不能改：混合训练后签名同构（主导特征与单数据集版本一致），"
            "签名由模型×任务结构决定" if keep else
            "能改：混合训练后签名主导特征改变，签名受训练数据覆盖影响")

    out = {
        "task": "data-coverage ablation: signature comparison (step12)",
        "label": "y = (1 - pearson_delta) > median (worst half), same as step6",
        "units": units, "signature_cos_sim": sim,
        "baseline_comparison": baseline_cmp, "verdict": verdict,
    }
    with open(os.path.join(ANALYSIS_DIR, "signature_comparison.json"), "w",
              encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"[out] -> {ANALYSIS_DIR}/signature_comparison.json")

    # ---- 图（300dpi，英文） ----
    make_figure(units, names, sigs)


def make_figure(units, names, sigs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    # 左：AUROC 条形（全探针 + 关键单特征）
    ax = axes[0]
    x = np.arange(len(names))
    w = 0.25
    ax.bar(x - w, [units[n]["probe"]["auroc_mean"] for n in names], w,
           yerr=[units[n]["probe"]["auroc_std"] for n in names],
           label="full probe (9-feat)", color="#4c72b0")
    sf = {n: units[n]["probe"]["single_feature"] for n in names}
    ax.bar(x, [sf[n].get("delta_max", {}).get("auroc_mean", np.nan)
               for n in names], w, label="delta_max alone", color="#dd8452")
    ax.bar(x + w, [sf[n].get("log_n_cells", {}).get("auroc_mean", np.nan)
                   for n in names], w, label="log_n_cells alone",
           color="#55a868")
    ax.axhline(0.5, ls="--", c="gray", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("5-fold CV AUROC")
    ax.set_ylim(0.4, 1.0)
    ax.set_title("Gating probe AUROC by training-data regime")
    ax.legend(fontsize=8, loc="lower right")

    # 右：签名热图（就地标准化系数）
    ax = axes[1]
    fn = units[names[0]]["feature_names"]
    im = ax.imshow(sigs, aspect="auto", cmap="RdBu_r",
                   vmin=-np.abs(sigs).max(), vmax=np.abs(sigs).max())
    ax.set_xticks(range(len(fn)))
    ax.set_xticklabels(fn, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_title("Error-probe signature (in-situ logistic coef)")
    fig.colorbar(im, ax=ax, shrink=0.8, label="standardized coefficient")
    for i in range(len(names)):
        j = int(np.argmax(np.abs(sigs[i])))
        ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1,
                                   fill=False, edgecolor="k", lw=1.5))
    fig.tight_layout()
    p = os.path.join(ANALYSIS_DIR, "signature_comparison.png")
    fig.savefig(p, dpi=300)
    fig.savefig(p.replace(".png", ".pdf"))
    print(f"[out] -> {p}")


if __name__ == "__main__":
    main()
