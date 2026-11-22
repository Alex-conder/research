# step9_dataset_property.py — 机制入口二（探索性登记，2026-10-16，不入主协议）
#
# 机制问题：为什么不同（模型×数据集）的误差签名长这样？数据集属性能否解释？
#
# ============================================================================
# 预注册假设（全部为 EXPLORATORY；n 小，不做 Pearson 显著性，只做排序/单调/描述）
# ============================================================================
# D1（签名主导特征 ⇔ 误差来源构成）：
#   逐扰动误差 err = 1 − pearson_delta 可分解为「采样噪声成分」与「模型成分」。
#   采样噪声成分操作化：对该扰动条件把 h5ad 真实细胞随机对半拆分，两半各自
#   估计 delta，top-50 |delta| 基因上两半 Pearson 均值为 ρ_half（B=50 次）。
#   全样本 delta 相对总体 delta 的估计噪声方差是半样本的一半，故"总体完美
#   预测器"（pred = 总体 delta）的期望误差下限：
#       noise_floor = 1 − 1/√(1 + (1−ρ_half)/(2ρ_half))
#   （split-half + Spearman–Brown 校正；n_cells 越大 ρ_half 越近 1，下限越低）。
#   sampling_share = noise_floor / err（逐扰动，再取组合均值）。
#   判据：若某签名由 log_n_cells 主导，则该组合的 sampling_share 排序应
#   高于 delta_max / delta_entropy 主导的组合（排序/单调对照，≥4 组合）。
#   若所有组合的采样占比排序与签名主导特征排序不一致 → D1 证伪。
# D2（就地可学性 ⇔ 误差的几何可预测性）：
#   就地探针 AUROC（归档值）应与「err 对 9 特征的 5-fold CV 回归 R²」正相关。
#   判据：Spearman ρ > 0 且至少单调；n=4~5 只做描述，不下显著性结论。
# D3（跨数据集签名消失 ⇔ 数据集属性漂移）：
#   量化 K562 vs RPE1 的属性漂移（真值效应强度分布 KS 距离/分位数比、
#   ctrl 逐基因表达离散度比、每条件细胞数分布比、逐条件细胞间异质性比），
#   与签名余弦对照（scGenePT K562 vs RPE1 ≈ 0.01 已知；GEARS K562 vs RPE1
#   在本脚本内按 step8b 同规格就地拟合后计算——GEARS RPE1 签名未归档，
#   此为探索性补算并留痕）。n=2 数据集对，纯描述。
#
# 数据集属性清单（每组合计算；逐扰动层面 + 数据集汇总层面）：
#   - 真值效应强度：top-50 |true_delta| 均值（true_delta = 条件全样本均值 −
#     ctrl 全样本均值，h5ad var 空间 5000 基因）
#   - ctrl 表达分布：逐基因均值/方差/零比例的汇总统计（中位数/IQR/均值，log 空间）
#   - 逐条件细胞数分布（h5ad 实测与标签值并列）
#   - 逐条件细胞间异质性（h5ad 含逐细胞表达，无需降级）：条件内逐基因
#     无偏方差的基因间均值
#
# 口径纪律（全部复用，不自造）：
#   - 特征/标签：step8b.load_xy（内部 = step6d.load_dataset / step6f.load_crossmodel，
#     特征 = step6_gating_probe.compute_features，y = (err > median)）
#   - 签名点估计核对：step8b.fit_coef 全量重拟合必须与
#     outputs_crossmodel/analysis/error_signatures.json 一致（<1e-6，assert）
#   - err 逐扰动值：GEARS 取 predictions_seed0.csv 的 pearson_delta；
#     跨模型取 _labels.csv 的 pearson_delta（与 step6f 同源文件）
#   - 就地 AUROC：GEARS 取 outputs{,_rpe1}/gating_probe.json 的 cv_auroc_mean；
#     跨模型取 outputs_crossmodel/analysis/{model}_{ds}.json 的 in_situ.auroc
#     （均为主 seed 0 / IN_SITU_SEEDS 均值归档值，本脚本不重跑探针 CV）
#
# 随机源：MASTER_SEED = 20261022（20261022 主 seed 体系专属 seed；每个组合的
# 蒙特卡洛与 CV 划分由主 seed 派生独立 Generator，互不复用同一随机流）。
#
# 产物：outputs/dataset_property.json + outputs/dataset_property.png（300dpi 英文）
#       outputs/dataset_property_verdict.md（人工撰写，本脚本只出数字）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step9_dataset_property.py
import json
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sstats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler

from step8b_signature_ci import (SOURCES, FEATURE_ORDER, fit_coef, load_xy)

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")
OUT_JSON = os.path.join(REPO, "outputs", "dataset_property.json")
OUT_PNG = os.path.join(REPO, "outputs", "dataset_property.png")

MASTER_SEED = 20261022          # 20261022 主 seed 体系（本脚本专属）
B_MC = 50                       # 噪声下限蒙特卡洛次数/扰动
TOP_N_DEG = 50                  # pearson_delta 口径（step1_run.py）

DATASETS = {
    "k562": {"data": "replogle_k562_essential", "h5ad": None},
    "rpe1": {"data": "replogle_rpe1_essential", "h5ad": None},
}
# 组合清单：4 个归档签名组合 + GEARS RPE1（签名未归档，探索性补算，D1/D2 第 5 点）
COMBOS = {
    "gears_k562":  {"sig": "GEARS(K562)",    "ds": "k562", "archived": True},
    "gears_rpe1":  {"sig": None,             "ds": "rpe1", "archived": False},
    "scgenept_k562": {"sig": "scgenept_k562", "ds": "k562", "archived": True},
    "scgenept_rpe1": {"sig": "scgenept_rpe1", "ds": "rpe1", "archived": True},
    "scf_k562":    {"sig": "scf_k562",       "ds": "k562", "archived": True},
}
ARCHIVED_AUROC = {   # 归档就地 AUROC（本脚本不重跑探针 CV，只做自洽核对）
    "gears_k562": ("outputs/gating_probe.json", "cv_auroc_mean"),
    "gears_rpe1": ("outputs_rpe1/gating_probe.json", "cv_auroc_mean"),
    "scgenept_k562": ("outputs_crossmodel/analysis/scgenept_k562.json", None),
    "scgenept_rpe1": ("outputs_crossmodel/analysis/scgenept_rpe1.json", None),
    "scf_k562": ("outputs_crossmodel/analysis/scfoundation_k562.json", None),
}


# ---------------------------------------------------------------- 组合加载
def load_combo(name):
    """加载组合的 (X, y, 签名, 逐扰动 err/n_cells, pred_expr/基因/ctrl_mean)。

    特征与标签口径 = step8b.load_xy（GEARS RPE1 走 step6d.load_dataset 同式）。
    返回 dict；perts 与 X 行严格对齐。
    """
    from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features
    info = COMBOS[name]
    if name == "gears_rpe1":
        from step6d_naive_baselines import load_dataset
        feat, feature_names, y, pred_df, out_path = load_dataset(
            "replogle_rpe1_essential", "outputs_rpe1")
        X = feat[feature_names].values.astype(np.float64)
        perts = feat.index.tolist()
        err = 1.0 - pred_df["pearson_delta"].values
        n_cells = pred_df["n_cells"].values.astype(int)
        z = np.load(os.path.join(out_path, "pred_expr_seed0.npz"))
        pred_expr = z["pred_expr"]
        pred_genes = [str(g) for g in z["pred_genes"].tolist()]
        adata = ad.read_h5ad(
            f"{REPO}/data/gears_demo/replogle_rpe1_essential/"
            "perturb_processed.h5ad", backed="r")
        var_genes = adata.var["gene_name"].tolist()
        cm = np.asarray(
            adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
        adata.file.close()
        ctrl_mean = (cm if var_genes == pred_genes else
                     np.array([dict(zip(var_genes, cm))[g] for g in pred_genes]))
        coef = fit_coef(X, y)
    elif info["sig"] is not None:
        X, y, feature_names = load_xy(info["sig"])
        assert list(feature_names) == FEATURE_ORDER
        coef = fit_coef(X, y)
        # 并行加载 pred_expr / 标签 err / n_cells（与 step6d/6f 同源文件）
        src = SOURCES[info["sig"]]
        kind = src[0]
        if kind == "gears":
            from step6d_naive_baselines import load_dataset
            _feat, _fn, _y, pred_df, out_path = load_dataset(src[1], src[2])
            perts = _feat.index.tolist()
            err = 1.0 - pred_df["pearson_delta"].values
            n_cells = pred_df["n_cells"].values.astype(int)
            z = np.load(os.path.join(out_path, "pred_expr_seed0.npz"))
            pred_expr, pred_genes = z["pred_expr"], [str(g) for g in z["pred_genes"].tolist()]
            adata = ad.read_h5ad(
                f"{REPO}/data/gears_demo/{src[1]}/perturb_processed.h5ad",
                backed="r")
            var_genes = adata.var["gene_name"].tolist()
            cm = np.asarray(
                adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
            adata.file.close()
            ctrl_mean = (cm if var_genes == pred_genes else
                         np.array([dict(zip(var_genes, cm))[g]
                                   for g in pred_genes]))
        else:
            from step6f_crossmodel import load_crossmodel
            _X, _y, _fn, perts, _extra = load_crossmodel(src[1], src[2])
            assert np.array_equal(_X, X) and np.array_equal(_y, y), \
                f"{name}: load_xy 与 load_crossmodel 结果不一致"
            lab = pd.read_csv(src[2]).set_index("perturbation").loc[perts]
            err = 1.0 - lab["pearson_delta"].values
            n_cells = lab["n_cells"].values.astype(int)
            z = np.load(src[1], allow_pickle=True)
            pred_expr = z["pred_expr"]
            pred_genes = [str(g) for g in z["pred_genes"].tolist()]
            key = ("ctrl_mean_matched" if "ctrl_mean_matched" in z.files
                   else "ctrl_mean_full")
            ctrl_mean = np.asarray(z[key], dtype=np.float64)
            nan_col = np.isnan(pred_expr).all(axis=0)
            if nan_col.any():  # 与 load_crossmodel 同式剔除全 NaN 列
                pred_expr = pred_expr[:, ~nan_col]
                pred_genes = [g for g, k in zip(pred_genes, ~nan_col) if k]
                ctrl_mean = ctrl_mean[~nan_col]
        # 剔除无真值行（与 load_crossmodel 同口径；GEARS load_dataset 已对齐）
        ok = ~np.isnan(err)
        if not ok.all():
            X, y = X[ok], y[ok]
            pred_expr, perts = pred_expr[ok], [p for p, k in zip(perts, ok) if k]
            err, n_cells = err[ok], n_cells[ok]
            coef = fit_coef(X, y)  # 重拟合以保持一致
    assert not np.isnan(X).any() and not np.isnan(err).any()
    return {"X": X, "y": y, "coef": coef, "perts": perts, "err": err,
            "n_cells": n_cells, "pred_expr": pred_expr,
            "pred_genes": pred_genes, "ctrl_mean": ctrl_mean}


# ------------------------------------------------------- 真值数据集属性（h5ad）
class TruthStats:
    """每个真值 h5ad 只打开一次：ctrl 逐基因汇总 + 各条件全样本 delta/异质性。"""

    def __init__(self, data_name):
        path = f"{REPO}/data/gears_demo/{data_name}/perturb_processed.h5ad"
        self.adata = ad.read_h5ad(path, backed="r")
        self.var_genes = self.adata.var["gene_name"].tolist()
        self.g2i = {g: i for i, g in enumerate(self.var_genes)}
        self.conds = self.adata.obs["condition"].astype(str).values
        self.cond_set = set(self.conds)
        ctrl = self.adata[self.conds == "ctrl"].X.tocsr()
        self.n_ctrl = ctrl.shape[0]
        self.ctrl_mean = np.asarray(ctrl.mean(axis=0)).ravel()
        m2 = np.asarray(ctrl.power(2).mean(axis=0)).ravel()
        self.ctrl_var = np.maximum(m2 - self.ctrl_mean ** 2, 0.0)
        nnz = np.diff(ctrl.indptr)  # 逐基因非零细胞数
        self.ctrl_zero_frac = 1.0 - nnz / self.n_ctrl
        self._cache = {}

    def find_cond(self, gene):
        for c in (f"ctrl+{gene}", f"{gene}+ctrl"):
            if c in self.cond_set:
                return c
        return None

    def condition_stats(self, gene):
        """(全样本 delta(5000,), 细胞数, 异质性=条件内逐基因方差的基因间均值)"""
        if gene in self._cache:
            return self._cache[gene]
        cond = self.find_cond(gene)
        if cond is None:
            self._cache[gene] = None
            return None
        X = self.adata[self.conds == cond].X.tocsr()
        n = X.shape[0]
        mean = np.asarray(X.mean(axis=0)).ravel()
        m2 = np.asarray(X.power(2).mean(axis=0)).ravel()
        var = np.maximum(m2 - mean ** 2, 0.0)
        delta = mean - self.ctrl_mean
        het = float(var.mean())
        self._cache[gene] = (delta, n, het, X)
        return self._cache[gene]

    def split_half_rho(self, gene, rng, b_mc, idx=None):
        """把条件细胞随机对半拆 B 次，两半 delta 在 top-50 |full delta| 基因上
        Pearson 的均值（该条件 delta 估计的 split-half reliability）。
        idx: 匹配基因空间的列索引（None = 全 5000 基因空间）。"""
        out = self.condition_stats(gene)
        if out is None:
            return None
        delta_full_all, n, het, X = out
        if n < 4:
            return None
        if idx is not None:
            delta_full = delta_full_all[idx]
        else:
            delta_full = delta_full_all
        top = np.argsort(-np.abs(delta_full))[:TOP_N_DEG]
        if np.allclose(delta_full[top], 0):
            return None
        n_half = n // 2
        rs = np.empty(b_mc)
        for b in range(b_mc):
            perm = rng.permutation(n)
            dA = np.asarray(X[perm[:n_half]].mean(axis=0)).ravel() - self.ctrl_mean
            dB = np.asarray(X[perm[n_half:2 * n_half]].mean(axis=0)).ravel() - self.ctrl_mean
            if idx is not None:
                dA, dB = dA[idx], dB[idx]
            rs[b] = np.corrcoef(dA[top], dB[top])[0, 1]
        return float(np.mean(rs))

    def close(self):
        self.adata.file.close()


# ---------------------------------------------------------------- 分析主体
def noise_floor_for_combo(combo, ts, rng, b_mc=B_MC):
    """逐扰动采样噪声下限：split-half ρ → Spearman–Brown 校正的期望误差下限。

    noise_floor = 1 − 1/√(1 + (1−ρ_half)/(2ρ_half))（ρ_half>0；否则记不可用）。
    top-50 |true_delta| 基因在组合的匹配基因空间上取（step1_run 口径）。
    返回 (noise_floor 数组, rho_half 数组, 不可用扰动数)。
    """
    g2i = ts.g2i
    idx = np.array([g2i.get(g, -1) for g in combo["pred_genes"]])
    n_gene_missing = 0
    if not (idx >= 0).all():
        # 匹配基因空间有基因不在 h5ad 5000 基因内：退回全空间（留痕）
        n_gene_missing = int((idx < 0).sum())
        idx = None
    nf = np.full(len(combo["perts"]), np.nan)
    rho = np.full(len(combo["perts"]), np.nan)
    n_bad = 0
    for i, g in enumerate(combo["perts"]):
        r = ts.split_half_rho(g, rng, b_mc, idx=idx)
        if r is None or r <= 0:
            n_bad += 1
            continue
        rho[i] = r
        nf[i] = 1.0 - 1.0 / np.sqrt(1.0 + (1.0 - r) / (2.0 * r))
    return nf, rho, n_bad, n_gene_missing


def cv_r2(X, err, rng):
    """err ~ 9 特征的 5-fold CV R²（标准化 + OLS，划分种子由主 seed 派生）。"""
    Xs = StandardScaler().fit_transform(X)
    cv = KFold(n_splits=5, shuffle=True, random_state=int(rng.integers(0, 2**31)))
    scores = cross_val_score(LinearRegression(), Xs, err, cv=cv, scoring="r2")
    return float(scores.mean())


def main():
    rng_master = np.random.default_rng(MASTER_SEED)
    with open(os.path.join(ANALYSIS_DIR, "error_signatures.json"),
              encoding="utf-8") as f:
        archived_sigs = json.load(f)

    res = {"task": "mechanism entry 2: dataset property vs error signature "
                   "(step9, EXPLORATORY, registered 2026-10-16)",
           "status": "exploratory, outside main protocol; n small — "
                     "ranking/monotonic/descriptive only, no significance tests",
           "master_seed": MASTER_SEED, "b_mc": B_MC,
           "label": "y = (1 - pearson_delta) > median; err = 1 - pearson_delta "
                    "(step1_run top-50 |true_delta| Pearson)",
           "noise_floor_def": "per-pert split-half: mean corr between two "
                              "random cell halves (B=50) on top-50 |delta| "
                              "genes; noise_floor = 1 - 1/sqrt(1+(1-rho)/"
                              "(2*rho)) (Spearman-Brown); "
                              "sampling_share = noise_floor / err",
           "combos": {}, "dataset_summary": {}, "d3_dataset_drift": {}}

    # ---- 组合循环 ----
    combos = {}
    for name in COMBOS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        combo = load_combo(name)
        # 签名自洽核对（归档组合）
        if COMBOS[name]["sig"] is not None:
            arch = np.array([archived_sigs[COMBOS[name]["sig"]][f]
                             for f in FEATURE_ORDER])
            dev = float(np.max(np.abs(combo["coef"] - arch)))
            assert dev < 1e-6, f"{name} 签名重拟合与归档不一致 (dev={dev:.2e})"
        # 归档 AUROC 自洽核对
        path, key = ARCHIVED_AUROC[name]
        with open(os.path.join(REPO, path), encoding="utf-8") as f:
            rec = json.load(f)
        auroc = float(rec[key]) if key else float(rec["in_situ"]["auroc"])
        ts = TruthStats(DATASETS[COMBOS[name]["ds"]]["data"])
        nf, rho_h, n_bad, n_gmiss = noise_floor_for_combo(combo, ts, rng)
        ts.close()
        r2 = cv_r2(combo["X"], combo["err"], rng)
        share = nf / combo["err"]
        w = np.abs(combo["coef"]) / np.abs(combo["coef"]).sum()
        dom = FEATURE_ORDER[int(np.argmax(np.abs(combo["coef"])))]
        combos[name] = {"combo": combo, "noise_floor": nf,
                        "auroc": auroc, "r2": r2}
        res["combos"][name] = {
            "dataset": COMBOS[name]["ds"], "n": int(len(combo["y"])),
            "signature_archived": COMBOS[name]["sig"] is not None,
            "signature_coef": {f: float(c) for f, c in
                               zip(FEATURE_ORDER, combo["coef"])},
            "signature_dominant_feature": dom,
            "signature_abs_weight": {f: float(x) for f, x in
                                     zip(FEATURE_ORDER, w)},
            "in_situ_auroc_archived": auroc,
            "err_regression_cv_r2": r2,
            "noise_floor": {
                "mean": float(np.nanmean(nf)), "median": float(np.nanmedian(nf)),
                "split_half_rho_mean": float(np.nanmean(rho_h)),
                "per_pert": [None if np.isnan(v) else float(v) for v in nf],
                "n_unusable": int(n_bad),
                "n_genes_not_in_h5ad": int(n_gmiss)},
            "sampling_share_mean": float(np.nanmean(share)),
            "sampling_share_median": float(np.nanmedian(share)),
            "err_mean": float(combo["err"].mean()),
            "n_cells_median": float(np.median(combo["n_cells"])),
            "per_pert": {
                "err": combo["err"].tolist(),
                "n_cells": combo["n_cells"].tolist(),
                "sampling_share": [None if np.isnan(v) else float(v)
                                   for v in share]},
        }
        print(f"[{name}] n={len(combo['y'])} AUROC={auroc:.3f} "
              f"R2={r2:.3f} noise_floor={np.nanmean(nf):.4f} "
              f"share={np.nanmean(share):.3f} dom={dom}")

    # ---- 数据集汇总属性（每 h5ad 一次；逐扰动属性用该数据集各组合 perts 并集）----
    for ds_key, ds in DATASETS.items():
        ts = TruthStats(ds["data"])
        logv = np.log1p(ts.ctrl_var)
        res["dataset_summary"][ds_key] = {
            "h5ad_shape": [int(ts.adata.n_obs), int(ts.adata.shape[1])],
            "n_conditions": int(ts.adata.obs["condition"].nunique()),
            "n_ctrl_cells": int(ts.n_ctrl),
            "ctrl_gene_var_log1p": {"median": float(np.median(logv)),
                                    "mean": float(logv.mean()),
                                    "iqr": float(np.percentile(logv, 75)
                                                 - np.percentile(logv, 25))},
            "ctrl_gene_mean": {"median": float(np.median(ts.ctrl_mean)),
                               "mean": float(ts.ctrl_mean.mean())},
            "ctrl_zero_frac": {"median": float(np.median(ts.ctrl_zero_frac)),
                               "mean": float(ts.ctrl_zero_frac.mean())},
        }
        # 逐扰动属性：用该数据集所有组合的 perts 并集
        genes = sorted({g for name, t in combos.items()
                        if COMBOS[name]["ds"] == ds_key
                        for g in t["combo"]["perts"]})
        strength, het, ncell = {}, {}, {}
        for g in genes:
            out = ts.condition_stats(g)
            if out is None:
                continue
            delta, n, h, _X = out
            top = np.sort(np.abs(delta))[-TOP_N_DEG:].mean()
            strength[g] = float(top)
            het[g] = float(h)
            ncell[g] = int(n)
        res["dataset_summary"][ds_key]["per_pert"] = {
            "n_genes": len(genes),
            "true_strength_top50": strength,
            "within_condition_heterogeneity": het,
            "h5ad_n_cells": ncell,
        }
        ts.close()

    # ---- D1：采样占比排序 vs 签名主导特征排序 ----
    names4 = ["gears_k562", "scgenept_k562", "scgenept_rpe1", "scf_k562"]
    names5 = names4[:1] + ["gears_rpe1"] + names4[1:]
    order_share, order_wnc = {}, {}
    d1 = {}
    for ns in (names4, names5):
        tag = f"n{len(ns)}"
        order_share[tag] = sorted(ns, key=lambda n: -res["combos"][n]
                                  ["sampling_share_mean"])
        order_wnc[tag] = sorted(ns, key=lambda n: -res["combos"][n]
                                ["signature_abs_weight"]["log_n_cells"])
        rho, _ = sstats.spearmanr(
            [res["combos"][n]["sampling_share_mean"] for n in ns],
            [res["combos"][n]["signature_abs_weight"]["log_n_cells"] for n in ns])
        # 排序一致度：Kendall tau between the two orderings
        tau, _ = sstats.kendalltau(
            [order_share[tag].index(n) for n in ns],
            [order_wnc[tag].index(n) for n in ns])
        d1[tag] = {"order_by_sampling_share": order_share[tag],
                   "order_by_log_n_cells_weight": order_wnc[tag],
                   "spearman_rho_share_vs_wnc": float(rho),
                   "kendall_tau_orderings": float(tau)}
    res["D1"] = {"criterion": "log_n_cells-dominated combos should rank higher "
                              "in sampling share; inconsistent orderings "
                              "everywhere => falsified",
                 "descriptive_only": True, **d1}

    # ---- D2：AUROC vs CV R² ----
    d2 = {}
    for ns in (names4, names5):
        tag = f"n{len(ns)}"
        aucs = [res["combos"][n]["in_situ_auroc_archived"] for n in ns]
        r2s = [res["combos"][n]["err_regression_cv_r2"] for n in ns]
        rho, _ = sstats.spearmanr(aucs, r2s)
        d2[tag] = {"auroc": aucs, "cv_r2": r2s, "combos": ns,
                   "spearman_rho": float(rho),
                   "monotone": bool(np.all(np.diff([
                       r2s[aucs.index(a)] for a in sorted(aucs)]) >= 0))}
    res["D2"] = {"criterion": "rho > 0 and at least monotone; n=4~5 "
                              "descriptive only", **d2}

    # ---- D3：K562 vs RPE1 属性漂移 + 签名余弦 ----
    def sigvec(n):
        return np.array([res["combos"][n]["signature_coef"][f]
                         for f in FEATURE_ORDER])
    cos = lambda a, b: float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    cos_sg = cos(sigvec("scgenept_k562"), sigvec("scgenept_rpe1"))
    cos_g = cos(sigvec("gears_k562"), sigvec("gears_rpe1"))
    assert 0.0 <= abs(cos_sg) < 0.05, f"scGenePT 跨数据集签名余弦偏离已知 0.01: {cos_sg}"
    drift = {}
    for prop in ["true_strength_top50", "within_condition_heterogeneity",
                 "h5ad_n_cells"]:
        a = res["dataset_summary"]["k562"]["per_pert"][prop]
        b = res["dataset_summary"]["rpe1"]["per_pert"][prop]
        common = sorted(set(a) & set(b))
        va, vb = np.array([a[g] for g in common]), np.array([b[g] for g in common])
        drift[prop] = {
            "ks_stat": float(sstats.ks_2samp(
                list(a.values()), list(b.values())).statistic),
            "median_ratio_k562_over_rpe1": float(np.median(va) / np.median(vb)),
            "q25_ratio": float(np.percentile(va, 25) / np.percentile(vb, 25)),
            "q75_ratio": float(np.percentile(va, 75) / np.percentile(vb, 75)),
            "n_common_genes": int(len(common))}
    va = res["dataset_summary"]["k562"]["ctrl_gene_var_log1p"]["median"]
    vb = res["dataset_summary"]["rpe1"]["ctrl_gene_var_log1p"]["median"]
    drift["ctrl_gene_var_log1p_median_ratio"] = float(va / vb)
    res["D3"] = {"criterion": "larger property drift co-occurs with lower "
                              "cross-dataset signature cosine; n=2 dataset "
                              "pairs, purely descriptive",
                 "signature_cosine": {"scgenept_k562_vs_rpe1": cos_sg,
                                      "gears_k562_vs_rpe1": cos_g},
                 "property_drift": drift}

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", OUT_JSON)
    make_figure(res, OUT_PNG)
    print(json.dumps({"D1_n4": res["D1"]["n4"],
                      "D1_n5": res["D1"]["n5"],
                      "D2_n4": res["D2"]["n4"],
                      "D2_n5": res["D2"]["n5"],
                      "D3_cosine": res["D3"]["signature_cosine"],
                      "D3_drift": res["D3"]["property_drift"]},
                     ensure_ascii=False, indent=2))


def make_figure(res, path):
    order = ["gears_k562", "gears_rpe1", "scgenept_k562", "scgenept_rpe1",
             "scf_k562"]
    labels = ["GEARS\nK562", "GEARS\nRPE1", "scGenePT\nK562", "scGenePT\nRPE1",
              "scFoundation\nK562"]
    dom_color = {"log_n_cells": "#2ca02c", "delta_max": "#1f77b4",
                 "delta_entropy": "#ff7f0e"}
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))

    # (a) sampling share vs dominant signature feature
    ax = axes[0]
    shares = [res["combos"][n]["sampling_share_mean"] for n in order]
    doms = [res["combos"][n]["signature_dominant_feature"] for n in order]
    cols = [dom_color.get(d, "gray") for d in doms]
    bars = ax.bar(range(len(order)), shares, color=cols, edgecolor="gray",
                  lw=0.5)
    for b, s, d in zip(bars, shares, doms):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.003,
                f"{s:.2f}\n({d.replace('delta_', '')})", ha="center",
                va="bottom", fontsize=6.5)
    ax.set_ylim(0, max(shares) * 1.35)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("mean sampling-noise share of error")
    ax.set_title("(a) sampling share (color = dominant\nsignature feature)",
                 fontsize=9)
    ax.grid(axis="y", alpha=0.25)

    # (b) signature |weight| of the three mechanism features
    ax = axes[1]
    feats = ["log_n_cells", "delta_max", "delta_entropy"]
    x = np.arange(len(order))
    w = 0.26
    for j, (f, c) in enumerate(zip(feats, [dom_color[f] for f in feats])):
        vals = [res["combos"][n]["signature_abs_weight"][f] for n in order]
        b = ax.bar(x + (j - 1) * w, vals, w, color=c, edgecolor="gray",
                   lw=0.4, label=f)
        for bb, v in zip(b, vals):
            ax.text(bb.get_x() + bb.get_width() / 2, bb.get_height() + 0.006,
                    f"{v:.2f}", ha="center", fontsize=6, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("|standardized coef| share")
    ax.set_title("(b) signature weight on\ndata-support vs shape features",
                 fontsize=9)
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(axis="y", alpha=0.25)

    # (c) in-situ AUROC vs err~feature CV R^2
    ax = axes[2]
    aucs = [res["combos"][n]["in_situ_auroc_archived"] for n in order]
    r2s = [res["combos"][n]["err_regression_cv_r2"] for n in order]
    ax.scatter(aucs, r2s, c=[dom_color.get(d, "gray") for d in doms], s=55,
               edgecolor="black", lw=0.6, zorder=3)
    # 手工错开的标注（顶部三个点很挤）
    offsets = {"gears_k562": (6, -14), "gears_rpe1": (-10, 8),
               "scgenept_k562": (6, 4), "scgenept_rpe1": (-88, 6),
               "scf_k562": (8, 2)}
    for a_, r_, l, n in zip(aucs, r2s, labels, order):
        ax.annotate(l.replace("\n", " "), (a_, r_), textcoords="offset points",
                    xytext=offsets[n], fontsize=7)
    ax.axhline(0.0, color="gray", ls=":", lw=0.8)
    ax.axvline(0.5, color="red", ls=":", lw=0.8)
    ax.set_xlabel("in-situ probe AUROC (archived)")
    ax.set_ylabel("CV R² of err ~ 9 features")
    ax.set_title("(c) learnability vs geometric\npredictability of error "
                 "(n=5, descriptive)", fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("Dataset properties vs error signatures (EXPLORATORY, "
                 f"seed {res['master_seed']})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
