# step6e_disease_strata.py — 疾病相关基因分层：探针误差探测的生物学解读
# 用 gseapy 从 Enrichr 取 DisGeNET 库，取并集得到疾病相关基因集；
# 把测试集扰动基因按"是否落在疾病基因集"分成两层，在各层内分别算
# 探针 AUROC 与 delta_max 基线 AUROC（与 step6d 同 OOF 口径，3 划分种子均值）。
# 层内 n<30 时如实标注小样本、不给强结论。
# 产出：outputs{,_rpe1}/disease_strata.json + outputs/disease_strata.png
import datetime
import json
import os

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features
from step6d_naive_baselines import SEED, load_dataset

REPO = os.path.dirname(os.path.abspath(__file__))
DATASETS = [
    {"data": "replogle_k562_essential", "out": "outputs", "tag": "K562"},
    {"data": "replogle_rpe1_essential", "out": "outputs_rpe1", "tag": "RPE1"},
]
SPLIT_SEEDS = [0, 1, 2]
SMALL_N = 30


def get_disease_genes():
    """Enrichr DisGeNET 库并集；记录库名与抓取日期。
    gseapy 失败时回退到 Enrichr GMT 直链（网络抖动时更稳），仍失败则报错。"""
    lib, via = None, "gseapy.get_library"
    for attempt in range(3):
        try:
            lib = gp.get_library(name="DisGeNET", organism="Human")
            break
        except Exception as e:
            print(f"[Enrichr] gseapy attempt {attempt + 1} failed: {e}")
            lib = None
    if lib is None:
        import urllib.request
        url = ("https://maayanlab.cloud/Enrichr/geneSetLibrary"
               "?mode=text&libraryName=DisGeNET")
        with urllib.request.urlopen(url, timeout=120) as r:
            text = r.read().decode("utf-8")
        lib = {}
        for line in text.strip().split("\n"):
            name, _, rest = line.partition("\t")
            lib[name] = [g for g in rest.split("\t") if g]
        via = "Enrichr GMT direct download"
    genes = set()
    for gs in lib.values():
        genes.update(gs)
    meta = {"source": "Enrichr", "library": "DisGeNET",
            "organism": "Human", "n_gene_sets": len(lib),
            "n_union_genes": len(genes),
            "fetched": datetime.date.today().isoformat(),
            "fetch_method": via}
    return genes, meta


def oof_scores(X, y, split_seed, single_col=None):
    Xm = X if single_col is None else X[:, [single_col]]
    Xs = StandardScaler().fit_transform(Xm)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=split_seed)
    return cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                             Xs, y, cv=cv, method="predict_proba")[:, 1]


def strata_auroc(scores, y, mask):
    """单层内 OOF AUROC（3 划分种子均值）；样本不足或单类返回 None。"""
    if mask.sum() < SMALL_N or len(np.unique(y[mask])) < 2:
        return None
    return float(np.mean([roc_auc_score(y[mask], scores[s][mask])
                          for s in range(len(SPLIT_SEEDS))]))


def run_dataset(data_name, out_dir, tag, disease_genes, meta):
    feat, feature_names, y, pred_df, out_path = load_dataset(
        data_name, out_dir)
    X = feat[feature_names].values.astype(np.float64)
    perts = feat.index.tolist()
    in_disease = np.array([p in disease_genes for p in perts])

    j_max = feature_names.index("delta_max")
    probe = [oof_scores(X, y, s) for s in SPLIT_SEEDS]
    dmax = [oof_scores(X, y, s, single_col=j_max) for s in SPLIT_SEEDS]

    layers = {}
    for name, mask in [("disease_genes", in_disease),
                       ("non_disease_genes", ~in_disease)]:
        layers[name] = {
            "n": int(mask.sum()),
            "worst_half_rate": float(y[mask].mean()),
            "probe_auroc": strata_auroc(probe, y, mask),
            "delta_max_auroc": strata_auroc(dmax, y, mask),
            "small_sample": bool(mask.sum() < SMALL_N),
        }
    res = {"task": "disease-gene stratification of gate AUROC (step6e)",
           "dataset": tag, "n": int(len(y)),
           "disease_gene_set": meta, "layers": layers,
           "caveat": "strata with n<30 flagged small_sample; "
                     "no strong claims on small strata"}
    with open(f"{out_path}/disease_strata.json", "w",
              encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    line = "  ".join(
        f"{k}: n={v['n']}"
        + (f" probe={v['probe_auroc']:.3f} dmax={v['delta_max_auroc']:.3f}"
           if v["probe_auroc"] is not None else " AUROC=NA")
        + (" [small n]" if v["small_sample"] else "")
        for k, v in layers.items())
    print(f"[{tag}] {line}")
    return res


disease_genes, meta = get_disease_genes()
print(f"[Enrichr] DisGeNET: {meta['n_gene_sets']} sets, "
      f"{meta['n_union_genes']} union genes, fetched {meta['fetched']}")
results = [run_dataset(d["data"], d["out"], d["tag"], disease_genes, meta)
           for d in DATASETS]

# ---- 图：分层 AUROC 条形图（英文，300 dpi）----
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), sharey=True)
titles = ["Replogle K562 (n=217)", "Replogle 2022 RPE1 (n=306)"]
for ax, res, ti in zip(axes, results, titles):
    layers = res["layers"]
    names = ["disease_genes", "non_disease_genes"]
    x = np.arange(len(names))
    w = 0.36
    pv = [layers[k]["probe_auroc"] for k in names]
    dv = [layers[k]["delta_max_auroc"] for k in names]
    ns = [layers[k]["n"] for k in names]
    b1 = ax.bar(x - w / 2, [v if v is not None else 0 for v in pv], w,
                color="#1f77b4", edgecolor="gray", lw=0.5, label="full probe")
    b2 = ax.bar(x + w / 2, [v if v is not None else 0 for v in dv], w,
                color="#9467bd", edgecolor="gray", lw=0.5,
                label="delta_max (single)")
    for b, v, n in zip(list(b1) + list(b2), pv + dv, ns * 2):
        if v is not None:
            ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                    f"{v:.2f}", ha="center", fontsize=8)
    ax.axhline(0.5, color="red", ls=":", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([f"disease genes\n(n={ns[0]})",
                        f"other genes\n(n={ns[1]})"], fontsize=8)
    ax.set_ylim(0.4, 1.0)
    ax.set_title(ti)
    ax.grid(axis="y", alpha=0.25)
axes[0].set_ylabel("OOF AUROC (3 split-seed mean)")
axes[0].legend(fontsize=8, loc="upper right")
fig.suptitle("Gate AUROC within disease-gene vs other-gene queries "
             "(Enrichr DisGeNET union)")
fig.tight_layout()
fig.savefig(f"{REPO}/outputs/disease_strata.png", dpi=300)
print("[SAVED]", f"{REPO}/outputs/disease_strata.png")
