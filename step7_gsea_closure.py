# step7_gsea_closure.py — 应用价值闭环：门控过滤提升下游 GSEA 富集质量
# 回答应用章 claim："用门控分数过滤低置信度扰动预测后，下游通路富集的
# 假阳性下降、真通路保持"。
#
# 设计（先定 claim 再定口径）：
#   每个留出扰动两条基因排序——预测 delta 排序（pred_expr - ctrl 训练均值）
#   与真值 delta 排序（真实表达均值 - ctrl），分别做 prerank GSEA。
#   过滤只改变"保留哪些扰动"，GSEA 结果按扰动缓存，各臂复用同一份，
#   保证覆盖率混淆被完全控制（对比均在匹配覆盖率下进行）。
#
# 基因集：Enrichr GO_Biological_Process_2023（gseapy 在线获取，注明日期），
#   仅保留与 5000 基因宇宙交集 >= MIN_OVERLAP 的 term。
# GSEA 口径：gseapy.prerank，permutation_num=100，seed=0，min_size=10，
#   max_size=500；显著性 FDR q-val < 0.25（GSEA 惯例宽松线）。
#
# 真通路两种定义：
#   A（主，非循环）：该扰动基因自身所属 GO term 的并集（纯成员关系，
#     不依赖真值 DE）；recall_A = pred 显著集中命中 >=1 个自身 term 的扰动比例。
#   B（副，可能循环）：真值 delta 显著富集的 term 为真集；
#     recall_B = |pred_sig ∩ true_sig| / |true_sig|。
# 假阳性：FPR_pert = |pred_sig - true_sig| / |pred_sig|（pred_sig 为空则剔除留痕）。
#
# 臂（对比均在匹配覆盖率下）：
#   none      全部扰动（coverage=1）
#   gate 70/50 按门控分数升序保留最自信的前 70%/50%（高分=存疑，丢弃高分端）
#   random 70/50 随机丢同数量扰动（20 seed，阴性对照："只是丢数据"的解释）
# 统计：gate vs none 用交集扰动上的配对 Wilcoxon 符号秩检验；
#   gate vs random 用 gate 均值在 20 个随机种子均值分布中的分位数。
#
# 门控分数：与 step6b 完全同口径重算（9 特征谱特征 + log_n_cells 协变量，
#   5-fold CV OOF logistic 概率，SEED=0 划分），无泄漏。
#
# 产出：outputs{,_rpe1}/gsea_closure.json + outputs/gsea_closure.png（300dpi 英文）
# 用法：PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step7_gsea_closure.py
import json
import os
import pickle
import time

import anndata as ad
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features

SEED = 0
REPO = os.path.dirname(os.path.abspath(__file__))
DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs", "tag": "Norman"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1", "tag": "RPE1"},
]
LIBRARY = "GO_Biological_Process_2023"   # Enrichr，2026-10-16 获取
MIN_OVERLAP = 10        # term 与本数据 5000 基因宇宙的最小交集
FDR_CUT = 0.25          # GSEA 显著性阈值（惯例宽松线）
PERM_NUM = 100
COVERAGES = [1.0, 0.70, 0.50]
RANDOM_SEEDS = list(range(20))


def load_dataset(data_name, out_dir):
    """读 OOF 预测 + 真值表达，返回门控分数与逐扰动 delta 排序所需量。"""
    out_path = os.path.join(REPO, out_dir)
    z = np.load(f"{out_path}/pred_expr_seed{SEED}.npz")
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()

    adata = ad.read_h5ad(
        f"{REPO}/data/gears_demo/{data_name}/perturb_processed.h5ad",
        backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()

    # 真值：留出扰动组的真实表达均值（condition 形如 'GENE+ctrl' / 'ctrl+GENE'）
    true_delta = {}
    singles = [c for c in adata.obs.condition.unique()
               if "ctrl" in c and c != "ctrl"]
    cond_of = {}
    for c in singles:
        g = c.replace("ctrl+", "").replace("+ctrl", "")
        cond_of[g] = c
    for g in perts:
        sub = adata[adata.obs.condition == cond_of[g]]
        true_delta[g] = (np.asarray(sub.X.mean(axis=0)).ravel() - ctrl_mean)
    adata.file.close()

    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))
    pred_df = pd.read_csv(
        f"{out_path}/predictions_seed{SEED}.csv").set_index("perturbation")
    pred_df = pred_df.loc[perts]

    # 门控分数：step6b 同口径（高分 = 更可能是后一半 = 低置信）
    feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    feature_names = SPECTRUM_FEATURES + \
        [c for c in COVARIATES if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    err = 1.0 - pred_df["pearson_delta"].values
    y = (err > np.median(err)).astype(int)
    Xs = StandardScaler().fit_transform(X)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    gate = cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                             Xs, y, cv=cv, method="predict_proba")[:, 1]
    gate_scores = pd.Series(gate, index=perts)
    return perts, pred_genes, pred_expr, ctrl_aligned, true_delta, gate_scores


def fetch_gene_sets(pred_genes):
    """Enrichr 库 -> 过滤到与本宇宙有交集的 term dict。"""
    lib = gp.get_library(name=LIBRARY, organism="Human")
    gset = set(pred_genes)
    return {k: v for k, v in lib.items()
            if len(set(v) & gset) >= MIN_OVERLAP}


def gsea_sig_terms(rnk, gene_sets):
    """一次 prerank GSEA，返回 DataFrame(Term, NES, FDR)。"""
    r = gp.prerank(rnk=rnk, gene_sets=gene_sets, threads=-1,
                   permutation_num=PERM_NUM, min_size=MIN_OVERLAP,
                   max_size=500, seed=SEED, verbose=False)
    df = r.res2d[["Term", "NES", "FDR q-val"]].copy()
    df.columns = ["term", "nes", "fdr"]
    return df.reset_index(drop=True)


def run_gsea_all(perts, pred_genes, pred_expr, ctrl_aligned, true_delta,
                 gene_sets, cache_path):
    """逐扰动跑 pred/true 两份 GSEA；有缓存直接复用（同库同口径）。"""
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            cache = pickle.load(f)
        if cache.get("library") != LIBRARY:
            cache = {}
    res = cache.get("perts", {})
    t0 = time.time()
    for i, g in enumerate(perts):
        if g in res:
            continue
        pred_delta = pred_expr[i] - ctrl_aligned
        pred_rnk = pd.Series(pred_delta, index=pred_genes)
        true_rnk = pd.Series(true_delta[g], index=pred_genes)
        res[g] = {
            "pred": gsea_sig_terms(pred_rnk, gene_sets),
            "true": gsea_sig_terms(true_rnk, gene_sets),
        }
        if (i + 1) % 25 == 0:
            print(f"  [gsea] {i+1}/{len(perts)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            with open(cache_path, "wb") as f:
                pickle.dump({"library": LIBRARY, "perts": res}, f)
    with open(cache_path, "wb") as f:
        pickle.dump({"library": LIBRARY, "perts": res}, f)
    return res


def sig_terms(df):
    return set(df.loc[df["fdr"] < FDR_CUT, "term"])


def per_pert_metrics(g, gsea, pert_gene_sets):
    """单扰动的 (fpr, recall_A, recall_B)；不可定义处返回 None。"""
    pred_sig = sig_terms(gsea[g]["pred"])
    true_sig = sig_terms(gsea[g]["true"])
    fpr = (len(pred_sig - true_sig) / len(pred_sig)
           if pred_sig else None)
    own = pert_gene_sets.get(g, set())
    rec_a = float(len(pred_sig & own) > 0) if own else None
    rec_b = (len(pred_sig & true_sig) / len(true_sig)
             if true_sig else None)
    return fpr, rec_a, rec_b


def aggregate(sel, gsea, pert_gene_sets):
    """选中扰动子集上的三项指标均值（各自剔除不可定义者，留痕 n）。"""
    fprs, rec_as, rec_bs = [], [], []
    n_pred_empty = 0
    for g in sel:
        fpr, ra, rb = per_pert_metrics(g, gsea, pert_gene_sets)
        if fpr is None:
            n_pred_empty += 1
        else:
            fprs.append(fpr)
        if ra is not None:
            rec_as.append(ra)
        if rb is not None:
            rec_bs.append(rb)
    return {
        "coverage": len(sel) / len(gsea),
        "n": len(sel),
        "n_pred_sig_empty": n_pred_empty,
        "fpr": {"mean": float(np.mean(fprs)) if fprs else None,
                "median": float(np.median(fprs)) if fprs else None,
                "n": len(fprs),
                "per_pert": fprs},
        "recall_A": {"mean": float(np.mean(rec_as)) if rec_as else None,
                     "n": len(rec_as)},
        "recall_B": {"mean": float(np.mean(rec_bs)) if rec_bs else None,
                     "n": len(rec_bs)},
    }


def run_dataset(data_name, out_dir, tag):
    out_path = os.path.join(REPO, out_dir)
    print(f"=== [{tag}] load + gate scores ===", flush=True)
    perts, pred_genes, pred_expr, ctrl_aligned, true_delta, gate = \
        load_dataset(data_name, out_dir)

    print(f"=== [{tag}] fetch {LIBRARY} ===", flush=True)
    gene_sets = fetch_gene_sets(pred_genes)
    print(f"  terms after overlap>={MIN_OVERLAP}: {len(gene_sets)}", flush=True)

    # 扰动基因自身所属 term（真通路定义 A）
    pert_gene_sets = {g: {t for t, gs in gene_sets.items() if g in set(gs)}
                      for g in perts}

    cache_path = f"{out_path}/gsea_closure_cache.pkl"
    print(f"=== [{tag}] GSEA x2 x {len(perts)} perts ===", flush=True)
    gsea = run_gsea_all(perts, pred_genes, pred_expr, ctrl_aligned,
                        true_delta, gene_sets, cache_path)

    # ---- 各臂 ----
    order_conf = gate.sort_values().index.tolist()   # 低分=自信在前
    arms = {"none": {"desc": "no filtering (coverage=1)"}}
    for cov in COVERAGES[1:]:
        k = int(round(cov * len(perts)))
        keep = order_conf[:k]
        arms[f"gate_{cov:.2f}"] = {
            "desc": f"keep most-confident {cov:.0%} by gate score",
            "sel": keep}
        arms[f"random_{cov:.2f}"] = {"desc": f"random {cov:.0%} "
                                     f"({len(RANDOM_SEEDS)} seeds)",
                                     "seeds": {}}
        for s in RANDOM_SEEDS:
            r = np.random.default_rng(1000 + s)
            arms[f"random_{cov:.2f}"]["seeds"][s] = \
                list(r.permutation(perts)[:k])

    res = {"task": "GSEA downstream closure: gating filter improves "
                   "enrichment quality at matched coverage (step7)",
           "dataset": tag, "n_perturbations": len(perts),
           "library": LIBRARY, "library_date": "2026-10-16 (Enrichr via gseapy)",
           "gsea": f"prerank, permutation_num={PERM_NUM}, seed={SEED}, "
                   f"min_size={MIN_OVERLAP}, max_size=500",
           "significance": f"FDR q-val < {FDR_CUT}",
           "truth_pathways": {
               "A_primary_noncircular": "GO terms containing the perturbed "
                                        "gene itself (membership union)",
               "B_secondary": "terms significant (FDR<cut) on true-delta "
                              "ranking (circular with FPR reference, "
                              "secondary only)"},
           "fpr_definition": "|pred_sig - true_sig| / |pred_sig| per "
                             "perturbation, averaged",
           "gate_score": "same OOF protocol as step6b (9-feature probe, "
                         "high score = likely worst half)",
           "arms": {}}

    per_pert_all = {}
    for g in perts:
        fpr, ra, rb = per_pert_metrics(g, gsea, pert_gene_sets)
        per_pert_all[g] = {"fpr": fpr, "recall_A": ra, "recall_B": rb}

    for name, arm in arms.items():
        if name == "none":
            agg = aggregate(perts, gsea, pert_gene_sets)
            res["arms"][name] = {"desc": arm["desc"], **agg}
        elif name.startswith("gate"):
            agg = aggregate(arm["sel"], gsea, pert_gene_sets)
            res["arms"][name] = {"desc": arm["desc"], **agg}
        else:  # random: 逐 seed 聚合
            seed_aggs = [aggregate(sel, gsea, pert_gene_sets)
                         for sel in arm["seeds"].values()]
            fprs = [a["fpr"]["mean"] for a in seed_aggs]
            ra = [a["recall_A"]["mean"] for a in seed_aggs]
            rb = [a["recall_B"]["mean"] for a in seed_aggs]
            ex = seed_aggs[0]
            res["arms"][name] = {
                "desc": arm["desc"], "coverage": ex["coverage"],
                "n": ex["n"],
                "fpr": {"mean": float(np.mean(fprs)),
                        "std": float(np.std(fprs, ddof=1)),
                        "min": float(np.min(fprs)),
                        "max": float(np.max(fprs)),
                        "per_seed": fprs},
                "recall_A": {"mean": float(np.mean(ra)),
                             "std": float(np.std(ra, ddof=1))},
                "recall_B": {"mean": float(np.mean(rb)),
                             "std": float(np.std(rb, ddof=1))}}

    # ---- 统计检验 ----
    # 注意：单扰动的 FPR 由 (pred GSEA, true GSEA) 决定，与是否被过滤无关，
    # 因此"配对 none vs gate 的逐扰动 FPR 差"恒为 0（错误口径，已废弃）。
    # 正确口径是比较各臂保留/丢弃队列的组成：
    #  (a) 保留队列：gate 保留集 FPR 均值 vs 20 个随机保留集均值分布
    #      （分位数 + 逐扰点配对 Wilcoxon：d_g = fpr_g*(I(gate 丢 g) - I(随机丢 g))，
    #       检验门控是否系统性丢高 FPR 扰动；对 20 seed 取 p 中位数与最大值）；
    #  (b) 丢弃队列：gate 丢弃集 FPR 均值 vs 随机丢弃集分布
    #      （门控丢的是否正是假阳性多的扰动）。
    stats = {}
    fpr_vec = np.array([per_pert_all[g]["fpr"]
                        if per_pert_all[g]["fpr"] is not None else np.nan
                        for g in perts])
    for cov in COVERAGES[1:]:
        gname, rname = f"gate_{cov:.2f}", f"random_{cov:.2f}"
        gate_drop = np.array([g not in set(arms[gname]["sel"]) for g in perts])
        rand_kept_means = res["arms"][rname]["fpr"]["per_seed"]
        gate_mean = res["arms"][gname]["fpr"]["mean"]
        pct_kept = float(np.mean([m <= gate_mean for m in rand_kept_means]))
        x = np.where(gate_drop, fpr_vec, 0.0)
        pvals = []
        for sel in arms[rname]["seeds"].values():
            rdrop = np.array([g not in set(sel) for g in perts])
            y = np.where(rdrop, fpr_vec, 0.0)
            w = wilcoxon(x - y)          # 零差自动剔除
            pvals.append(float(w.pvalue))
        gate_dropped_mean = float(np.nanmean(fpr_vec[gate_drop]))
        rand_dropped = [
            float(np.nanmean(fpr_vec[np.array(
                [g not in set(sel) for g in perts])]))
            for sel in arms[rname]["seeds"].values()]
        pct_drop = float(np.mean([m >= gate_dropped_mean
                                  for m in rand_dropped]))
        stats[f"coverage_{cov:.2f}"] = {
            "kept_cohort": {
                "gate_mean_fpr": gate_mean,
                "random_mean_fpr": float(np.mean(rand_kept_means)),
                "random_std": float(np.std(rand_kept_means, ddof=1)),
                # P(random <= gate)，接近 0 支持门控优于随机丢数据
                "percentile_of_gate_among_random": pct_kept,
                "paired_wilcoxon_p_median_over_seeds":
                    float(np.median(pvals)),
                "paired_wilcoxon_p_max_over_seeds": float(np.max(pvals))},
            "dropped_cohort": {
                "gate_dropped_mean_fpr": gate_dropped_mean,
                "random_dropped_mean_fpr": float(np.mean(rand_dropped)),
                "random_dropped_std": float(np.std(rand_dropped, ddof=1)),
                # P(random >= gate)，接近 1 支持门控丢的正是高 FPR 扰动
                "percentile_of_gate_among_random": pct_drop}}
    res["stats"] = stats

    with open(f"{out_path}/gsea_closure.json", "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[{tag}] saved {out_path}/gsea_closure.json")
    for k, v in res["arms"].items():
        print(f"  {k}: FPR={v['fpr']['mean']:.3f} "
              f"recallA={v['recall_A']['mean']:.3f} "
              f"recallB={v['recall_B']['mean']:.3f} (n={v['n']})")
    return res


def main():
    results = [run_dataset(d["data"], d["out"], d["tag"]) for d in DATASETS]

    # ---- 图：FPR vs 覆盖率，三组对比（英文，300 dpi）----
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), sharey=True)
    titles = ["Replogle K562 (n=217)", "Replogle RPE1 (n=306)"]
    for ax, res, ti in zip(axes, results, titles):
        # 门控臂（两个覆盖点连线）
        gx, gy = [], []
        for cov in COVERAGES[1:]:
            a = res["arms"][f"gate_{cov:.2f}"]
            gx.append(a["coverage"])
            gy.append(a["fpr"]["mean"])
        a1 = res["arms"]["none"]
        ax.plot([1.0], [a1["fpr"]["mean"]], "o", color="#7f7f7f", ms=8,
                label="no filter")
        ax.plot(gx, gy, "-o", color="#1f77b4", ms=7, lw=1.8,
                label="gate filter (ours)")
        # 随机臂（均值 ± 1 std，带 jitter 的散点）
        for cov in COVERAGES[1:]:
            a = res["arms"][f"random_{cov:.2f}"]
            xs = [a["coverage"]] * len(a["fpr"]["per_seed"])
            ax.scatter(xs, a["fpr"]["per_seed"], s=14, alpha=0.45,
                       color="#d62728", zorder=2,
                       label="random drop (20 seeds)" if cov == 0.70 else None)
            ax.errorbar([a["coverage"]], [a["fpr"]["mean"]],
                        yerr=[a["fpr"]["std"]], fmt="s", color="#d62728",
                        ms=8, capsize=4, zorder=3)
        ax.set_xlabel("coverage (fraction of perturbations kept)")
        ax.set_title(ti)
        ax.set_xlim(0.4, 1.06)
        ax.set_ylim(0.38, 0.85)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("pathway FPR\n(sig. on prediction, not on truth)")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Gating filter reduces downstream GSEA false positives at "
                 "matched coverage", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.tight_layout()
    fig.savefig(f"{REPO}/outputs/gsea_closure.png", dpi=300)
    print("[SAVED]", f"{REPO}/outputs/gsea_closure.png")


if __name__ == "__main__":
    main()
