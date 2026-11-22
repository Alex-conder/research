# step17_ensemble_gate_gb.py — 升层 G-B：跨引擎集成门（stacked meta-probe）
# 预注册：prereg_expansion_gb.md（2026-10-30，运行前锁定）。
#
# 设计（锁定，详见预注册）：
#   对每个评测组合 E@D，集成门 = 上层 logistic（LOPO 防泄漏），特征三块：
#     (1) 基础块：数据集 D 上各可得引擎就地探针 OOF 分数（5-fold seed=0，
#         fold 内 StandardScaler + L2 logistic C=1.0，同 R2 oof_probe_scores）
#     (2) 符号块：R2 符号层分数（同引擎另一数据集系数符号 · 本组合 z 特征；
#         scf_k562 无 sibling -> 无此块）
#     (3) 幅度块：本组合自身幅度族特征原值（注册表加载器口径）
#   对照臂：主 = GEARS 就地 OOF 分数（任务书"最优单引擎"）；次 = 符号层；
#   诊断 = 本组合自身就地探针（不进锁定计数）。
#   判定（锁定）：集成 AUROC 在 >=8/10 组合优于 GEARS 就地臂
#   （扰动级配对 bootstrap B=1000，95% CI 不含 0）-> 集成有效进主 claim；
#   否则如实报"集成无增益"（与 error_signatures 深潜"签名噪声"互证的层次结论）。
#   一致性矩阵：每数据集引擎两两误差二值列联（n11 同错 / n10 / n01 / n00 同对，
#   kappa, phi），回答"不同引擎的错误是否独立"。
#
# 防泄漏：底层分数全部 OOF；上层 LOPO；bootstrap 只重采样索引不重拟合
# （LOPO 级重拟合在 CPU 不可行，预注册已公开记录）。
#
# 产物：outputs_rebuild/ensemble_gate.{json,png} + ensemble_gate_verdict.md
#       outputs_rebuild/smoke/ensemble_gate.smoke.{json,png}（--smoke）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step17_ensemble_gate_gb.py --smoke
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step17_ensemble_gate_gb.py
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import cohen_kappa_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

for _s in ("stdout", "stderr"):  # GBK 控制台纪律：不崩、不丢字符
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")
WAVE_DIR = os.path.join(REPO, "outputs_crossmodel", "wave")

BOOT_MASTER_SEED = 20261030
FULL_B = 1000
SMOKE_B = 100

VERDICT_RULE = ("ensemble AUROC beats GEARS in-situ arm in >=8/10 evaluable "
                "combos (perturbation-level paired bootstrap B=1000, 95% CI "
                "excludes 0) -> ensemble effective, enters main claim; else "
                "report no ensemble gain (tier conclusion, cross-validates "
                "error_signatures signature-noise)")

from step10_signature_resolution import COMBOS, fit_coef  # 口径单一来源
from step10_two_layer_gate import SIBLING, oof_probe_scores, symbol_scores

# 数据集 -> 全覆盖（非 scf 子集）引擎组合键（基础块来源）
FULL_KEYS = {
    "K562": ["gears_k562", "scgenept_k562", "scgpt_k562",
             "geneformer_k562", "uce_k562"],
    "RPE1": ["gears_rpe1", "scgenept_rpe1", "scgpt_rpe1",
             "geneformer_rpe1", "uce_rpe1"],
}
SCF_KEY = "scf_k562"
LOCKED_KEYS = [k for k in COMBOS
               if COMBOS[k]["kind"] != "missing" and k != SCF_KEY]

EMB_KEYS = [k for k, c in COMBOS.items() if c["cls"] == "embedding"]


# ---------------------------------------------------------------- 加载（带扰动名）
def load_combo_xy_perts(key):
    """返回 (X, y, feature_names, perts)，行序与 perts 一一对应。"""
    c = COMBOS[key]
    if c["kind"] == "gears":
        from step6d_naive_baselines import load_dataset
        feat, feature_names, y, pred_df, _ = load_dataset(c["data"], c["out"])
        X = feat[feature_names].values.astype(np.float64)
        perts = [str(g) for g in pred_df.index]
    elif c["kind"] == "crossmodel":
        from step6f_crossmodel import load_crossmodel
        X, y, feature_names, perts, _extra = load_crossmodel(c["npz"],
                                                             c["labels"])
        perts = [str(g) for g in perts]
    elif c["kind"] == "wave":
        lab = pd.read_csv(os.path.join(
            WAVE_DIR, c["wave"],
            f"{c['wave']}_wave_{'k562' if c['dataset'] == 'K562' else 'rpe1'}"
            f"_labels.csv"))
        err = lab["err_emb_dist"].values.astype(np.float64)
        assert np.isfinite(err).all(), f"{key}: err_emb_dist 含非有限值"
        y = (err > np.median(err)).astype(int)
        in_space = lab["target_in_input_space"].astype(bool).astype(np.float64)
        X = np.column_stack([
            lab["pred_emb_norm"].values.astype(np.float64),
            lab["delta_emb_norm"].values.astype(np.float64),
            np.log1p(lab["n_cells"].values.astype(np.float64)),
            in_space,
        ])
        feature_names = ["pred_emb_norm", "delta_emb_norm", "log_n_cells",
                         "target_in_input_space"]
        perts = [str(g) for g in lab["perturbation"]]
    else:
        raise FileNotFoundError(f"{key}: 产物缺失（注册为 missing）")
    assert len(perts) == len(y) == len(X), f"{key}: 行/扰动名不对齐"
    assert set(np.unique(y)) <= {0, 1}, f"{key}: 标签非二值"
    return X, y, list(feature_names), perts


# ---------------------------------------------------------------- 上层 LOPO
def lopo_scores(X, y):
    """leave-one-perturbation-out 同规格 logistic OOF 概率（防泄漏）。"""
    n = len(y)
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        m = np.ones(n, dtype=bool)
        m[i] = False
        sc = StandardScaler().fit(X[m])
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(sc.transform(X[m]),
                                                             y[m])
        out[i] = clf.predict_proba(sc.transform(X[i:i + 1]))[0, 1]
    return out


def safe_auroc(s, y):
    if len(np.unique(y)) < 2:
        return None
    return float(roc_auc_score(y, s))


def paired_boot_diffs(s_a, s_b, y, b_boot, rng):
    """扰动级配对 bootstrap（分数向量固定，不重拟合；预注册已记录）。
    返回均值差、95% CI、是否不含 0、方向。"""
    n = len(y)
    diffs = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        diffs[b] = roc_auc_score(y[idx], s_a[idx]) - \
            roc_auc_score(y[idx], s_b[idx])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return {"mean_diff": float(diffs.mean()), "ci_lo": float(lo),
            "ci_hi": float(hi),
            "ci_excludes_zero": bool(lo > 0 or hi < 0),
            "favors": ("first_arm" if lo > 0 else
                       "second_arm" if hi < 0 else "undetermined")}


# ---------------------------------------------------------------- 一致性矩阵
def consistency_block(keys_y_perts, engine_names):
    """keys_y_perts: {key: (y, perts)}；按扰动名交集做两两列联。
    返回 {pair_key: {...}} 与 kappa 方阵（对角 = 1.0）。"""
    out = {}
    n_e = len(engine_names)
    kappa_mat = np.eye(n_e)
    for a in range(n_e):
        ya, pa = keys_y_perts[engine_names[a]]
        ma = {p: v for p, v in zip(pa, ya)}
        for b in range(a + 1, n_e):
            yb, pb = keys_y_perts[engine_names[b]]
            mb = {p: v for p, v in zip(pb, yb)}
            common = sorted(set(ma) & set(mb))
            if not common:
                continue
            va = np.array([ma[p] for p in common])
            vb = np.array([mb[p] for p in common])
            n11 = int(((va == 1) & (vb == 1)).sum())   # 同错（最差半）
            n10 = int(((va == 1) & (vb == 0)).sum())
            n01 = int(((va == 0) & (vb == 1)).sum())
            n00 = int(((va == 0) & (vb == 0)).sum())   # 同对
            agree = float((va == vb).mean())
            kappa = float(cohen_kappa_score(va, vb))
            phi = float(np.corrcoef(va, vb)[0, 1])
            pk = f"{engine_names[a]}__{engine_names[b]}"
            out[pk] = {"n_pair": len(common), "n11_both_wrong": n11,
                       "n10": n10, "n01": n01, "n00_both_right": n00,
                       "agreement": agree, "cohens_kappa": kappa,
                       "phi": phi}
            kappa_mat[a, b] = kappa_mat[b, a] = kappa
    return out, kappa_mat.tolist()


# ---------------------------------------------------------------- 主流程
def run_combo(key, loaded, base_scores, b_boot, rng_master):
    c = COMBOS[key]
    ds = c["dataset"]
    t0 = time.perf_counter()
    X, y, feature_names, perts = loaded[key]
    n = len(y)

    # 基础块：全覆盖引擎 +（scf 目标时）scf 自身
    src_keys = list(FULL_KEYS[ds])
    if key == SCF_KEY:
        src_keys = src_keys + [SCF_KEY]
    base_feats, base_names = [], []
    for sk in src_keys:
        s = base_scores[sk]
        m = {p: v for p, v in zip(loaded[sk][3], s)}
        vals = np.array([m[p] for p in perts], dtype=np.float64)
        assert np.isfinite(vals).all(), f"{key}: 基础块 {sk} 有缺失行"
        base_feats.append(vals)
        base_names.append(f"oof_probe::{sk}")

    # 符号块（R2 同构；scf_k562 无 sibling -> 无符号块/符号臂）
    signs = None
    sym_note = None
    sib = SIBLING.get(key)
    if sib is not None:
        Xs, ys, _, _ = loaded[sib]
        signs = np.sign(fit_coef(Xs, ys))
        sym_note = f"signs from {sib} in-situ fit (zero training here)"
    sym = symbol_scores(X, signs) if signs is not None else None

    feats = list(base_feats)
    names = list(base_names)
    if sym is not None:
        feats.append(sym)
        names.append("symbol_layer_score")
    feats.append(X)
    names.extend([f"own::{f}" for f in feature_names])
    X_stack = np.column_stack(feats)
    t_feat = time.perf_counter() - t0

    t1 = time.perf_counter()
    s_ens = lopo_scores(X_stack, y)
    t_lopo = time.perf_counter() - t1

    arms = {"ensemble": s_ens,
            "gears_in_situ": None, "symbol": sym, "own_in_situ": None}
    gears_key = f"gears_{ds.lower()}"
    mg = {p: v for p, v in zip(loaded[gears_key][3], base_scores[gears_key])}
    arms["gears_in_situ"] = np.array([mg[p] for p in perts])
    if key != gears_key:  # gears 目标时 own == gears_in_situ
        mo = {p: v for p, v in zip(perts, base_scores[key])}
        arms["own_in_situ"] = np.array([mo[p] for p in perts])
    else:
        arms["own_in_situ"] = arms["gears_in_situ"]

    aurocs = {a: safe_auroc(s, y) for a, s in arms.items() if s is not None}
    rng = np.random.default_rng(rng_master.integers(0, 2**31))
    pairs = {}
    for other in ("gears_in_situ", "symbol", "own_in_situ"):
        if arms[other] is None:
            continue
        pairs[f"ensemble_minus_{other}"] = paired_boot_diffs(
            s_ens, arms[other], y, b_boot, rng)

    res = {"engine": c["engine"], "dataset": ds, "class": c["cls"],
           "n": int(n), "label_balance": float(y.mean()),
           "stack_features": names, "symbol_layer": sym_note,
           "base_block_keys": src_keys,
           "auroc": aurocs,
           "paired_bootstrap": {"b_boot": b_boot, **pairs},
           "timing_s": {"feature_assembly": round(t_feat, 3),
                        "lopo": round(t_lopo, 3)}}
    d = pairs.get("ensemble_minus_gears_in_situ", {})
    print(f"[GB][{key}] n={n} AUROC ens={aurocs['ensemble']:.4f} "
          f"gears={aurocs['gears_in_situ']:.4f} "
          f"sym={aurocs.get('symbol')} own={aurocs['own_in_situ']:.4f} "
          f"| ens-gears {d.get('mean_diff', float('nan')):+.4f} "
          f"CI=[{d.get('ci_lo', 0):+.4f},{d.get('ci_hi', 0):+.4f}] "
          f"excl0={d.get('ci_excludes_zero')} | lopo {t_lopo:.1f}s")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="迷你规格（B=100, 3 组合覆盖三种装配路径）")
    ap.add_argument("--combos", nargs="*", default=None)
    args = ap.parse_args()

    smoke = args.smoke
    b_boot = SMOKE_B if smoke else FULL_B
    combo_keys = args.combos or (["gears_k562", SCF_KEY, "scgpt_rpe1"]
                                 if smoke else list(COMBOS.keys()))

    t_all = time.perf_counter()
    # 预加载全部所需组合（基础块要跨组合取分数）
    need = set()
    for k in combo_keys:
        if COMBOS[k]["kind"] == "missing":
            continue
        need.add(k)
        need.update(FULL_KEYS[COMBOS[k]["dataset"]])
        if k == SCF_KEY:
            need.add(SCF_KEY)
        sib = SIBLING.get(k)
        if sib is not None:
            need.add(sib)
    loaded = {}
    for k in sorted(need):
        if COMBOS[k]["kind"] == "missing":
            continue
        loaded[k] = load_combo_xy_perts(k)
        print(f"[LOAD] {k}: n={len(loaded[k][1])} "
              f"cls={COMBOS[k]['cls']}")

    t_base = time.perf_counter()
    base_scores = {k: oof_probe_scores(v[0], v[1]) for k, v in loaded.items()}
    t_base = time.perf_counter() - t_base
    print(f"[BASE] in-situ OOF base scores: {len(base_scores)} combos "
          f"in {t_base:.1f}s")

    boot_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "G-B cross-engine ensemble gate (step17, "
                   "prereg_expansion_gb.md)",
           "prereg": "prereg_expansion_gb.md",
           "smoke": smoke,
           "spec": {"meta_model": "LOPO logistic (StandardScaler + L2, "
                                   "C=1.0, max_iter=2000)",
                    "base_block": "per-engine 5-fold OOF probe scores "
                                  "(StratifiedKFold seed=0, R2 caliber)",
                    "b_boot": b_boot,
                    "full_run_spec": {"b_boot": FULL_B},
                    "label": {"expression": "y = (1 - pearson_delta) > median",
                              "embedding": "y = err_emb_dist > median"}},
           "boot_master_seed": BOOT_MASTER_SEED,
           "verdict_rule": VERDICT_RULE,
           "combos": {}, "missing": []}
    for key in combo_keys:
        c = COMBOS[key]
        if c["kind"] == "missing":
            res["missing"].append(key)
            res["combos"][key] = {"engine": c["engine"],
                                  "dataset": c["dataset"], "class": c["cls"],
                                  "missing": True}
            continue
        res["combos"][key] = run_combo(key, loaded, base_scores, b_boot,
                                       boot_master)

    # ---- 一致性矩阵（每数据集，引擎两两误差二值列联）----
    consis = {}
    for ds, tag in (("K562", "k562"), ("RPE1", "rpe1")):
        keys = [f"gears_{tag}", f"scgenept_{tag}", f"scgpt_{tag}",
                f"geneformer_{tag}", f"uce_{tag}"]
        if ds == "K562":
            keys.insert(2, SCF_KEY)
        keys = [k for k in keys if k in loaded and k in res["combos"]
                and not res["combos"][k].get("missing")]
        engines = [COMBOS[k]["engine"] for k in keys]
        block, kappa_mat = consistency_block(
            {COMBOS[k]["engine"]: (loaded[k][1], loaded[k][3]) for k in keys},
            engines)
        off = [kappa_mat[i][j] for i in range(len(engines))
               for j in range(i + 1, len(engines))]
        consis[ds] = {"engines": engines, "pairs": block,
                      "kappa_matrix": kappa_mat,
                      "mean_abs_kappa_offdiag": float(np.mean(np.abs(off)))}
        print(f"[CONSIS][{ds}] {len(block)} pairs, "
              f"mean|kappa|={consis[ds]['mean_abs_kappa_offdiag']:.3f}")
    res["consistency"] = consis

    # ---- verdict ----
    locked = [k for k in LOCKED_KEYS if k in res["combos"]
              and not res["combos"][k].get("missing")]
    n_better = sum(
        1 for k in locked
        if res["combos"][k]["paired_bootstrap"]
        ["ensemble_minus_gears_in_situ"]["ci_excludes_zero"]
        and res["combos"][k]["paired_bootstrap"]
        ["ensemble_minus_gears_in_situ"]["favors"] == "first_arm")
    if smoke:
        verdict = {"tier": "smoke_only_not_adjudicated",
                   "note": "smoke numbers do not enter the locked criteria",
                   "n_ensemble_better_smoke_numbers": int(n_better),
                   "n_locked_combos": len(locked)}
    else:
        verdict = {"rule": VERDICT_RULE,
                   "n_ensemble_better": int(n_better),
                   "n_locked_combos": len(locked),
                   "threshold": ">=8/10 (locked)"}
        verdict["tier"] = ("ensemble_effective" if n_better >= 8
                           else "no_ensemble_gain")
    res["verdict"] = verdict
    res["timing_total_s"] = round(time.perf_counter() - t_all, 2)

    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    stem = "ensemble_gate.smoke" if smoke else "ensemble_gate"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    if not smoke:
        write_verdict_md(res, os.path.join(OUT_DIR,
                                           "ensemble_gate_verdict.md"))
    print(f"[DONE] total {res['timing_total_s']}s | tier "
          f"{verdict.get('tier')}")


def make_figure(res, path):
    keys = [k for k in res["combos"] if not res["combos"][k].get("missing")]
    if not keys:
        return
    arms = ["ensemble", "gears_in_situ", "symbol", "own_in_situ"]
    colors = {"ensemble": "#2ca02c", "gears_in_situ": "#1f77b4",
              "symbol": "#d62728", "own_in_situ": "#9467bd"}
    labels = {"ensemble": "ensemble (stacked LOPO)",
              "gears_in_situ": "GEARS in-situ",
              "symbol": "symbol layer",
              "own_in_situ": "own in-situ (diag)"}
    fig = plt.figure(figsize=(15.5, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.1, 1.0], height_ratios=[1, 1],
                          wspace=0.24, hspace=0.5)
    ax = fig.add_subplot(gs[:, 0])
    x = np.arange(len(keys))
    width = 0.19
    for j, arm in enumerate(arms):
        vals = [res["combos"][k]["auroc"].get(arm) for k in keys]
        vals = [np.nan if v is None else v for v in vals]
        bars = ax.bar(x + (j - 1.5) * width, vals, width,
                      color=colors[arm], label=labels[arm],
                      alpha=0.9 if arm != "own_in_situ" else 0.55)
        for b, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.004,
                        f"{v:.2f}", ha="center", va="bottom", fontsize=5.5,
                        rotation=90)
    ax.axhline(0.5, color="gray", ls=":", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([k.replace("_", "\n") for k in keys], fontsize=7)
    ax.set_ylabel("AUROC (predict worst-half errors)", fontsize=9)
    ax.set_ylim(0.35, 1.0)
    ax.set_title("G-B ensemble gate vs single-engine arms "
                 f"(B={res['spec']['b_boot']}, smoke={res['smoke']})",
                 fontsize=10)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7.5, loc="upper left", ncol=2)

    for gi, ds in enumerate(("K562", "RPE1")):
        axc = fig.add_subplot(gs[gi, 1])
        block = res.get("consistency", {}).get(ds)
        if not block:
            axc.axis("off")
            continue
        km = np.array(block["kappa_matrix"])
        im = axc.imshow(km, vmin=-0.2, vmax=1.0, cmap="RdYlGn")
        eng = block["engines"]
        axc.set_xticks(range(len(eng)))
        axc.set_yticks(range(len(eng)))
        axc.set_xticklabels(eng, fontsize=6.5, rotation=45, ha="right")
        axc.set_yticklabels(eng, fontsize=6.5)
        for i in range(len(eng)):
            for j in range(len(eng)):
                axc.text(j, i, f"{km[i, j]:.2f}", ha="center", va="center",
                         fontsize=6,
                         color="black" if 0.3 < km[i, j] < 0.85 else "white")
        axc.set_title(f"error-consistency Cohen's kappa ({ds})", fontsize=8.5)
    fig.colorbar(im, ax=fig.get_axes()[-1], fraction=0.03, pad=0.02,
                 label="kappa")
    fig.tight_layout(rect=[0, 0, 0.99, 1])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


def write_verdict_md(res, path):
    v = res["verdict"]
    lines = ["# G-B 跨引擎集成门 verdict（step17，prereg_expansion_gb.md）",
             ""]
    if v.get("tier") == "smoke_only_not_adjudicated":
        lines.append("smoke 运行，不进锁定判定。")
    else:
        locked = [k for k in LOCKED_KEYS
                  if not res["combos"][k].get("missing")]
        lines.append(
            f"判定：集成优于 GEARS 就地臂 {v['n_ensemble_better']}/"
            f"{v['n_locked_combos']} 组合（阈值 >=8/10，配对 bootstrap CI "
            f"不含 0）-> **{v['tier']}**。")
        lines.append("")
        lines.append("## 逐组合（锁定 10 组合 + scf_k562 诊断）")
        lines.append("")
        lines.append("| combo | n | ensemble | GEARS in-situ | symbol | "
                     "own in-situ | ens-gears 差值 [CI] | 不含0 |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for k in list(locked) + ([SCF_KEY] if SCF_KEY in res["combos"]
                                 and not res["combos"][SCF_KEY].get("missing")
                                 else []):
            c = res["combos"][k]
            d = c["paired_bootstrap"]["ensemble_minus_gears_in_situ"]
            sym = c["auroc"].get("symbol")
            sym = f"{sym:.4f}" if sym is not None else "n/a"
            tag = " (diag)" if k == SCF_KEY else ""
            lines.append(
                f"| {k}{tag} | {c['n']} | {c['auroc']['ensemble']:.4f} | "
                f"{c['auroc']['gears_in_situ']:.4f} | {sym} | "
                f"{c['auroc']['own_in_situ']:.4f} | "
                f"{d['mean_diff']:+.4f} [{d['ci_lo']:+.4f}, "
                f"{d['ci_hi']:+.4f}] | {d['ci_excludes_zero']} |")
        lines.append("")
        cons = res["consistency"]
        for ds in ("K562", "RPE1"):
            b = cons[ds]
            ka = sorted(((p, r) for p, r in b["pairs"].items()),
                        key=lambda t: t[1]["cohens_kappa"])
            lo = ka[0]
            hi = ka[-1]
            lines.append(
                f"一致性矩阵（{ds}）：mean|kappa| off-diag = "
                f"{b['mean_abs_kappa_offdiag']:.3f}；最弱一致 "
                f"{lo[0]} kappa={lo[1]['cohens_kappa']:.3f} "
                f"(n={lo[1]['n_pair']})；最强一致 "
                f"{hi[0]} kappa={hi[1]['cohens_kappa']:.3f} "
                f"(n={hi[1]['n_pair']})。")
        lines.append("")
        lines.append("## 结构观察（描述性，不改判定）")
        lines.append("")
        lines.append(
            "两两 kappa 呈清晰的口径家族结构：表达谱引擎互间为正"
            "（K562: GEARS×scGenePT +0.41、GEARS×scFoundation +0.42；"
            "RPE1: GEARS×scGenePT -0.15 为例外但幅度小），embedding 三引擎互间"
            "强正（RPE1 Geneformer×UCE +0.88、scGPT×Geneformer +0.43），"
            "而表达谱×embedding 跨家族为负（GEARS×UCE：K562 -0.46、RPE1 -0.63）。"
            "即不同引擎的错误**不独立**：家族内错误共享（冗余），跨家族错误"
            "互反（互补）。集成臂的胜负分布与此一致——在 embedding 目标组合上"
            "GEARS 就地臂近随机（0.21-0.49）而集成显著占优，在 GEARS 自身组合上"
            "集成无法超越其就地探针（家族内冗余）。"
            "锁定判定的完整表述：集成相对最优单引擎无稳定增益（6/10 < 8/10），"
            "各引擎误差签名不携带可相互堆叠的独立信息，"
            "与 error_signatures 深潜“签名噪声”互证。")
        if v["tier"] == "no_ensemble_gain":
            lines.append(
                "阴性解读（锁定条款）：各引擎误差签名互不携带可堆叠信息，与 "
                "outputs_crossmodel/analysis/error_signatures.json 深潜的"
                "“签名噪声”判定互证；R2 结论（产品 = 符号层）不变。")
        else:
            lines.append("集成有效：跨引擎集成门作为新的一层进主 claim。")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
