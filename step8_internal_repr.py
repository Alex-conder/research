# step8_internal_repr.py — 分水岭实验：GEARS 内部隐层表示 vs 输出几何特征，
# 谁更能预测模型自身的扰动级误差？
#
# 背景与问题：step6 的输出几何探针（9 个预测谱统计量特征）在 K562 217 扰动上
# AUROC ≈ 0.70；step7 跨模型路由显示"弱模型的赢点不在输出几何、可能在内部表示"。
# 本实验直接检验：GEARS 的**最终节点嵌入**（GNN 融合扰动 embedding 后、decoder 前）
# 是否比输出几何特征更能预测"该扰动会被答错"。
#
# ─────────────────────────────────────────────────────────────────────────────
# 预注册式设计（先写死再跑）：
#
# 1. 数据：K562 Replogle essential，217 个留出扰动（与 step6/7 同一批）。
#    真值 h5ad：data/gears_demo/replogle_k562_essential/perturb_processed.h5ad；
#    标签：outputs/predictions_seed{m}.csv 的 pearson_delta
#    （top-50 |真实 delta| 基因上 Pearson，step1_run.py:160-162 口径）；
#    err = 1 - pearson_delta；y = err > **训练侧中位数**（step7 fit_probe 口径）。
#
# 2. 内部表示提取（不重训，只用已有 checkpoint 前向）：
#    - 加载 outputs/gears_model_seed{m}/config.pkl + model.pt，
#      直接构造 gears.model.GEARS_Model（config 已含 G_go/G_coexpress/num_genes/
#      num_perts，无需 PertData；pert_list 复现自 essential_all_data_pert_genes.pkl
#      ∩ gene2go_all.pkl，与 gears/pertdata.py:126-128 完全一致，9853 个）。
#    - 提取点：GEARS_Model.forward 中 decoder 输入，即 model.py:186-190
#      `base_emb = bn_pert_base(...); base_emb = transform(base_emb)` 之后、
#      `recovery_w(base_emb)` 之前的最终节点嵌入（hidden_size=64 维/节点），
#      用 forward_pre_hook 挂在 recovery_w 上捕获（shape: n_cells*5000 × 64）。
#    - 扰动输入构造与推理完全一致：复刻 gears.gears.predict 路径
#      （gears.py:349-361）——create_cell_graph_dataset_for_prediction 从 ctrl 细胞
#      有放回抽 300 个细胞、DataLoader 组 batch、eval + no_grad。
#      唯一区别：predict 内 np.random.randint 未播种，这里先
#      np.random.seed(20261016) 固定（可复现性纪律，分布相同）。
#    - 每个扰动三种定长向量（300 细胞平均）：
#        meanpool  : 全部 5000 基因节点嵌入做 mean-pooling → 64 维；
#        pertnode  : 仅取被扰动基因对应节点的嵌入 → 64 维；
#        crossgene : cross_gene_state（model.py:197，第一层 decoder 后对
#                    全基因状态的 5000→64 全局 MLP）输出 → 64 维。
#      架构事实（决定三种读出的信息量，写入 verdict）：扰动 embedding 只在
#      model.py:166-184 被加到**被扰动基因的那一两个节点**上，随后
#      bn_pert_base（eval 模式逐通道仿射，无跨节点混合）+ ReLU 即进入 decoder。
#      因此提取点处除被扰动节点外的所有节点嵌入**逐扰动不变**——meanpool
#      结构性接近于扰动不变（预期≈随机），pertnode 是提取点处唯一有信号的
#      读出；crossgene 是扰动信息首次经跨基因混合后的全局内部状态，
#      是"内部表示是否知道误差"的最强公平检验，一并记录。
#      注意：217 个留出扰动中 134 个的靶基因不在 5000 维表达空间内
#      （step7 已记录此数据事实：pred ≈ ctrl），这些扰动在 cell graph 中
#      没有对应节点，pertnode 对它们无定义（记 NaN）。因此 pertnode64 与
#      geom9+pertnode64 仅在 83 个"靶基因在表达空间内"的扰动子集上评估；
#      配对基线 geom9 在同一子集上重算，split 用与 step7 完全相同的生成
#      逻辑（step7_splits(n_sub)），此为数据驱动的必要协议偏差，如实记录。
#    - 提取同时记录该 forward 的预测谱，与归档 pred_expr_seed{m}.npz 逐扰动
#      Pearson 作 sanity check（应 ≈1，验证前向复刻正确）。
#    - 模型种子：m ∈ {0, 1}（seed0 主实验、seed1 复现）。
#
# 3. 特征集（每种都与输出几何用**完全相同**的 split 与探针）：
#    a. geom9    : 8 个谱统计量 + log_n_cells（step6_gating_probe.compute_features
#                  同函数，逐模型种子用 pred_expr_seed{m} 重算；seed0 应复现 ≈0.70）
#    b. meanpool64 / pertnode64 / crossgene64 : 内部表示（单独）
#    c. geom9+meanpool64 / geom9+pertnode64 / geom9+crossgene64 : 拼接
#    探针规格与 step6/7 完全相同：StandardScaler + LogisticRegression(C=1.0,
#    max_iter=2000)，70% 训练侧就地拟合（含中位数阈值），30% 测试侧 AUROC。
#
# 4. split（必须与 step7 完全一致以便对比）：扰动级 70/30，
#    SPLIT_SEEDS=[0,1,2]，master seed 20261016：
#    rng_master = default_rng(20261016)；每个 split_seed 用
#    rng = default_rng(rng_master.integers(0, 2**31))；perm = rng.permutation(n)；
#    n_te = int(round(0.30*n))；te, tr = perm[:n_te], perm[n_te:]
#    （逐行复刻 step7_routing.py:112-119 / step7b_routing_calibrated.py:237-244）。
#
# 5. bootstrap：每个（模型种子 × split 种子）在 30% 测试集上以扰动为抽样单位
#    重采样 2000 次：(i) 每个特征集 AUROC 的 95% CI；(ii) 特征集 − geom9 配对
#    AUROC 差的 95% CI（内部表示判据用）。
#
# 6. 判据（预先写死，对每个非基线特征集在其各自样本上判定，
#    对每个模型种子）：
#    - d = 3 个 split 种子上 (AUROC_set − AUROC_geom9) 的均值；
#    - "误差签名藏在内部表示"：d ≥ 0.05 且 3 个 split 种子的配对 bootstrap 95% CI
#      下界全部 > 0；
#    - "输出几何已是足够近似"：|d| < 0.02；
#    - 其余 → 中间档，如实报告。
#
# 产物：outputs/internal_repr.json（全部数字）、outputs/internal_repr.png
# （300 dpi，英文标签）、outputs/internal_repr_verdict.md（≤1 页，中文）。
#
# 纪律：不改任何现有文件；不重训 GEARS；所有随机源固定（split 20261016 体系、
# 提取 np.random.seed(20261016)、bootstrap 每组合独立 default_rng 派生）。
#
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step8_internal_repr.py
import json
import os
import pickle
import time

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import DataLoader

from step6_gating_probe import SPECTRUM_FEATURES, compute_features

REPO = os.path.dirname(os.path.abspath(__file__))
DATA_NAME = "replogle_k562_essential"
H5AD = os.path.join(REPO, "data", "gears_demo", DATA_NAME, "perturb_processed.h5ad")
PERT_GENES_PKL = os.path.join(REPO, "data", "gears_demo", "essential_all_data_pert_genes.pkl")
GENE2GO_PKL = os.path.join(REPO, "data", "gears_demo", "gene2go_all.pkl")

MODEL_SEEDS = [0, 1]          # GEARS checkpoint 种子（0 主实验，1 复现）
SPLIT_SEEDS = [0, 1, 2]       # 与 step7 相同的划分种子
TEST_FRAC = 0.30
MASTER_SEED = 20261016        # 与 step7 相同的 master seed
EXTRACT_SEED = 20261016       # ctrl 细胞抽样的固定种子
N_BOOT = 2000
N_CELLS_PER_PERT = 300        # 与 gears predict 默认 num_samples 一致

INTERNAL_VARIANTS = ["meanpool64", "pertnode64", "crossgene64"]
# 七种特征集：geom9 为基线；后三种为拼接
FEATURE_SETS = ["geom9", "meanpool64", "pertnode64", "crossgene64",
                "geom9+meanpool64", "geom9+pertnode64", "geom9+crossgene64"]
# 判据阈值（预注册）
THRESH_WIN, THRESH_TIE = 0.05, 0.02


def get_pert_list():
    """复现 gears/pertdata.py:126-128 的 pert_names（perturbation graph 节点名单）。"""
    essential = pickle.load(open(PERT_GENES_PKL, "rb"))
    gene2go = pickle.load(open(GENE2GO_PKL, "rb"))
    g2g = {i: gene2go[i] for i in essential if i in gene2go}
    return np.unique(list(g2g.keys()))


def build_model(model_seed):
    """从 config.pkl + model.pt 直接构造 GEARS_Model（不依赖 PertData）。"""
    from gears.model import GEARS_Model
    mdir = os.path.join(REPO, "outputs", f"gears_model_seed{model_seed}")
    cfg = pickle.load(open(os.path.join(mdir, "config.pkl"), "rb"))
    model = GEARS_Model(cfg).to(cfg["device"])
    sd = torch.load(os.path.join(mdir, "model.pt"),
                    map_location=cfg["device"])
    model.load_state_dict(sd)
    model.eval()
    return model, cfg


def extract_internal_representations(model_seed, perts, ctrl_adata, pert_list,
                                     gene_names_5000, cfg):
    """对每个扰动复刻 gears predict 前向，hook 取最终节点嵌入。

    返回：meanpool (n_pert, 64)、pertnode (n_pert, 64)、crossgene (n_pert, 64)、
          sanity 预测谱 (n_pert, 5000)。
    关键代码依据（venv/Lib/site-packages/gears/）：
      - 提取点 model.py:186-190：bn_pert_base → transform(ReLU) → recovery_w；
        hook 挂在 recovery_w（forward_pre_hook，捕获其输入）→ meanpool/pertnode；
      - crossgene：model.py:197 cross_gene_embed = cross_gene_state(out)，
        forward hook 捕获输出（每 graph 一个 64 维全局向量）。
      - 输入构造 gears.py:349-361 + utils.py:457-477：
        ctrl 细胞有放回抽 300 → create_cell_graph_for_prediction → DataLoader 组 batch。
    注意两套索引不要混淆：
      - pi = 扰动在 pert_list（9853 节点扰动图，pert_emb 行号）中的位置，
        决定"哪个扰动 embedding 被加到图里"（forward 的 pert_idx）；
      - gene_idx = 扰动基因在 5000 维表达空间（adata.var/pred_genes 顺序，
        即 cell graph 节点顺序）中的位置，决定 pertnode 取哪个节点。
    """
    from gears.utils import create_cell_graph_for_prediction
    device = cfg["device"]
    num_genes = cfg["num_genes"]
    model, _ = build_model(model_seed)
    gene_idx = {g: i for i, g in enumerate(gene_names_5000)}
    pert_arr = np.array(pert_list)

    np.random.seed(EXTRACT_SEED)  # predict 内 ctrl 抽样未播种，这里固定（分布相同）
    n = len(perts)
    meanpool = np.zeros((n, cfg["hidden_size"]), dtype=np.float32)
    pertnode = np.zeros((n, cfg["hidden_size"]), dtype=np.float32)
    crossgene = np.zeros((n, cfg["hidden_size"]), dtype=np.float32)
    pred_spectra = np.zeros((n, num_genes), dtype=np.float32)
    captured = {}

    def hook_base(module, inputs):
        captured["base"] = inputs[0].detach()

    def hook_cross(module, inputs, output):
        captured["cross"] = output.detach()

    h1 = model.recovery_w.register_forward_pre_hook(hook_base)
    h2 = model.cross_gene_state.register_forward_hook(hook_cross)
    t0 = time.time()
    with torch.no_grad():
        for i, g in enumerate(perts):
            pi = [int(np.where(g == pert_arr)[0][0])]  # pert_emb 行号（pert_list 索引）
            gj = gene_idx.get(g)                       # 5000 维表达空间基因索引（不在→None）
            Xs = ctrl_adata[np.random.randint(0, len(ctrl_adata),
                                              N_CELLS_PER_PERT), :].X.toarray()
            cg = [create_cell_graph_for_prediction(X, pi, [g]).to(device)
                  for X in Xs]
            batch = next(iter(DataLoader(cg, N_CELLS_PER_PERT, shuffle=False)))
            batch.to(device)
            out = model(batch)
            base = captured["base"].cpu().numpy()          # (n_cells*5000, 64)
            assert base.shape[0] == N_CELLS_PER_PERT * num_genes
            node = base.reshape(N_CELLS_PER_PERT, num_genes, -1)
            meanpool[i] = node.mean(axis=(0, 1))           # 全基因节点 mean-pool
            # 靶基因不在 5000 维表达空间的扰动无对应节点 → NaN（协议见头部注释）
            pertnode[i] = (node[:, gj, :].mean(axis=0)
                           if gj is not None else np.nan)
            cross = captured["cross"].cpu().numpy()        # (n_cells, 64)
            assert cross.shape[0] == N_CELLS_PER_PERT
            crossgene[i] = cross.mean(axis=0)
            pred_spectra[i] = out.mean(dim=0).cpu().numpy()
            if (i + 1) % 50 == 0:
                print(f"  [extract seed{model_seed}] {i+1}/{n} "
                      f"({time.time()-t0:.0f}s)", flush=True)
    h1.remove()
    h2.remove()
    return meanpool, pertnode, crossgene, pred_spectra


def load_geometry(model_seed, perts):
    """与 step6d.load_dataset 同口径的输出几何 9 特征 + err（按 npz 扰动名单对齐）。

    step6d 固定 SEED=0；这里按模型种子参数化（pred_expr_seed{m}.npz /
    predictions_seed{m}.csv），以便"同模型种子内部表示 vs 同模型种子输出几何"的
    公平对比；seed0 即复现 step6 基线（应 ≈0.70）。
    """
    out_path = os.path.join(REPO, "outputs")
    z = np.load(os.path.join(out_path, f"pred_expr_seed{model_seed}.npz"))
    pred_expr, npz_perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()
    adata = ad.read_h5ad(H5AD, backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    adata.file.close()
    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))
    pred_df = pd.read_csv(
        os.path.join(out_path, f"predictions_seed{model_seed}.csv")
    ).set_index("perturbation").loc[npz_perts]

    feat = compute_features(pred_expr, npz_perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    feature_names = SPECTRUM_FEATURES + ["log_n_cells"]

    # 按目标扰动名单对齐（与 step7 的扰动顺序无关——split 在本函数外按同一
    # perts 顺序生成，这里保证 geom / internal / err 三者严格对齐）
    X = feat.loc[perts, feature_names].values.astype(np.float64)
    err = 1.0 - pred_df.loc[perts, "pearson_delta"].values
    return X, err, feature_names


def step7_splits(n):
    """逐行复刻 step7_routing.py:112-119 的 split 生成（master 20261016）。"""
    rng_master = np.random.default_rng(MASTER_SEED)
    out = []
    for sp_seed in SPLIT_SEEDS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        perm = rng.permutation(n)
        n_te = int(round(TEST_FRAC * n))
        out.append((sp_seed, perm[n_te:], perm[:n_te]))  # (seed, tr, te)
    return out


def fit_probe(X_tr, err_tr):
    """step7 fit_probe 同规格：y = err > 训练中位数；StandardScaler + L2 logistic。"""
    y_tr = (err_tr > np.median(err_tr)).astype(int)
    if len(np.unique(y_tr)) < 2:
        return None
    scaler = StandardScaler().fit(X_tr)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(
        scaler.transform(X_tr), y_tr)
    return scaler, clf


def probe_scores(fit, X):
    scaler, clf = fit
    return scaler.transform(X) @ clf.coef_[0] + clf.intercept_[0]


def auroc_ci(y, s, rng, n_boot=N_BOOT):
    """以扰动为抽样单位的 AUROC bootstrap 95% CI。"""
    n = len(y)
    stat = np.empty(n_boot)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        yb, sb = y[idx], s[idx]
        # 重采样可能缺类：回退到全量统计量
        if len(np.unique(yb)) < 2:
            stat[b] = roc_auc_score(y, s)
            continue
        # 向量化 AUROC（避免 sklearn 在退化情形抛错）：P(s_pos > s_neg) 平局半计
        d = sb[pos][:, None] - sb[neg][None, :]
        stat[b] = (np.mean(d > 0) + 0.5 * np.mean(d == 0))
    return float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5))


def vec_auroc(y, s):
    """向量化 AUROC（P(s_pos>s_neg)，平局半计），与 roc_auc_score 等价。"""
    sp, sn = s[y == 1], s[y == 0]
    d = sp[:, None] - sn[None, :]
    return float(np.mean(d > 0) + 0.5 * np.mean(d == 0))


def diff_auroc_ci(y, s_new, s_ref, rng, n_boot=N_BOOT):
    """配对 bootstrap：AUROC(新特征集) − AUROC(geom9) 的 95% CI（测试集内）。"""
    n = len(y)
    stat = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            stat[b] = vec_auroc(y, s_new) - vec_auroc(y, s_ref)
            continue
        stat[b] = vec_auroc(yb, s_new[idx]) - vec_auroc(yb, s_ref[idx])
    return float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5))


def probe_pass(X_sets, err, names, splits, m):
    """在一组样本上对若干特征集跑 step7 同规格探针协议。

    返回 (per_split entries, summary dict)。names 中须含 "geom9" 作配对基线
    （pairs 中每个非 geom9 特征集都与它在同一测试集上比 AUROC 差）。
    bootstrap rng 按 (MASTER_SEED, model_seed, split_seed) 固定派生。
    """
    entries = []
    for sp_seed, tr, te in splits:
        fit_geom = fit_probe(X_sets["geom9"][tr], err[tr])
        assert fit_geom is not None
        s_geom_te = probe_scores(fit_geom, X_sets["geom9"][te])
        y_te = (err[te] > np.median(err[tr])).astype(int)
        rng = np.random.default_rng((MASTER_SEED, m, sp_seed))
        entry = {"split_seed": sp_seed,
                 "n_train": int(len(tr)), "n_test": int(len(te)),
                 "auroc": {}, "auroc_ci95": {}, "vs_geom9_diff_ci95": {}}
        for name in names:
            fit = fit_probe(X_sets[name][tr], err[tr])
            assert fit is not None
            s_te = probe_scores(fit, X_sets[name][te])
            a = vec_auroc(y_te, s_te)
            assert abs(a - roc_auc_score(y_te, s_te)) < 1e-9
            entry["auroc"][name] = a
            entry["auroc_ci95"][name] = auroc_ci(y_te, s_te, rng)
            if name != "geom9":
                entry["vs_geom9_diff_ci95"][name] = diff_auroc_ci(
                    y_te, s_te, s_geom_te, rng)
        entries.append(entry)
        line = "  ".join(f"{k}={v:.3f}" for k, v in entry["auroc"].items())
        print(f"[seed{m}/split{sp_seed}] n={len(tr)}+{len(te)} {line}")

    summary = {}
    for name in names:
        aucs = [e["auroc"][name] for e in entries]
        entry_sum = {"sample_n": int(len(err)),
                     "auroc_mean": float(np.mean(aucs)),
                     "auroc_std": float(np.std(aucs, ddof=1)),
                     "auroc_per_split": aucs,
                     "auroc_ci95_per_split": [e["auroc_ci95"][name]
                                              for e in entries]}
        if name != "geom9":
            diffs = [e["auroc"][name] - e["auroc"]["geom9"] for e in entries]
            los = [e["vs_geom9_diff_ci95"][name][0] for e in entries]
            entry_sum.update({
                "diff_vs_geom9_mean": float(np.mean(diffs)),
                "diff_vs_geom9_per_split": [float(d) for d in diffs],
                "vs_geom9_diff_ci95_per_split": [
                    e["vs_geom9_diff_ci95"][name] for e in entries],
                "all_ci_lower_bounds_above_0": bool(all(l > 0 for l in los))})
            d = entry_sum["diff_vs_geom9_mean"]
            if d >= THRESH_WIN and entry_sum["all_ci_lower_bounds_above_0"]:
                verdict = "win"
            elif abs(d) < THRESH_TIE:
                verdict = "tie"
            else:
                verdict = "intermediate"
            entry_sum["criterion_verdict"] = verdict
        summary[name] = entry_sum
    return entries, summary



def run_all():
    # ---- 数据与对齐 ----
    z0 = np.load(os.path.join(REPO, "outputs", "pred_expr_seed0.npz"))
    perts = z0["pred_perturbations"].tolist()          # 217 个扰动（与 step6 同批）
    n = len(perts)
    adata = ad.read_h5ad(H5AD, backed="r")
    ctrl_adata = adata[adata.obs.condition == "ctrl"]
    pert_list = get_pert_list()
    splits = step7_splits(n)
    print(f"[data] n={n} 扰动, ctrl cells={len(ctrl_adata)}, "
          f"pert graph nodes={len(pert_list)}, test sizes="
          f"{[len(te) for _s, _tr, te in splits]}")

    results = {"task": "watershed: internal hidden representation vs output "
                       "geometry for predicting GEARS self-error (step8)",
               "dataset": DATA_NAME, "n_perturbations": n,
               "model_seeds": MODEL_SEEDS, "split_seeds": SPLIT_SEEDS,
               "master_seed": MASTER_SEED, "extract_seed": EXTRACT_SEED,
               "split_protocol": "replicated step7_routing.py:112-119 "
                                 "(perturbation-level 70/30, master 20261016)",
               "label": "y = (1 - pearson_delta) > train-median (step6/step7)",
               "probe_spec": "StandardScaler + LogisticRegression(C=1.0, "
                             "max_iter=2000), in-situ on 70% train",
               "extraction_point": "GEARS_Model.forward decoder input "
                                   "(model.py:186-190): node embeddings after "
                                   "bn_pert_base+ReLU, before recovery_w; "
                                   "hook = forward_pre_hook on recovery_w",
               "input_construction": "replicates gears.gears.predict "
                                     "(gears.py:349-361, utils.py:457-477): "
                                     "300 bootstrapped ctrl cells per pert, "
                                     "eval+no_grad; np.random.seed(20261016) "
                                     "for reproducibility",
               "bootstrap": {"n_resamples": N_BOOT, "unit": "perturbation",
                             "ci": "95% percentile"},
               "criterion": {"win": f"mean diff >= {THRESH_WIN} and all 3 "
                                    "per-split paired bootstrap CIs lower "
                                    "bound > 0",
                             "tie": f"|mean diff| < {THRESH_TIE}",
                             "applied_to": ["meanpool64", "pertnode64",
                                            "crossgene64",
                                            "geom9+meanpool64",
                                            "geom9+pertnode64",
                                            "geom9+crossgene64"]},
               "model_seed_results": {}}

    X_sets_per_seed = {}
    for m in MODEL_SEEDS:
        print(f"[extract] model seed {m} ...", flush=True)
        model, cfg = build_model(m)
        del model  # 提取函数内部自行加载；此处只为尽早报错
        cache = os.path.join(REPO, "outputs",
                             f"internal_repr_cache_v2_seed{m}.npz")
        if os.path.exists(cache):
            cz = np.load(cache)
            meanpool, pertnode = cz["meanpool"], cz["pertnode"]
            crossgene, pred_spectra = cz["crossgene"], cz["pred_spectra"]
            print(f"  [cache] loaded {cache}")
        else:
            meanpool, pertnode, crossgene, pred_spectra = \
                extract_internal_representations(m, perts, ctrl_adata,
                                                 pert_list,
                                                 z0["pred_genes"].tolist(), cfg)
            np.savez_compressed(cache, meanpool=meanpool, pertnode=pertnode,
                                crossgene=crossgene, pred_spectra=pred_spectra,
                                perts=np.array(perts))
            print(f"  [cache] saved {cache}")

        # sanity：提取时的预测谱 vs 归档 pred_expr
        z = np.load(os.path.join(REPO, "outputs", f"pred_expr_seed{m}.npz"))
        ref = dict(zip(z["pred_perturbations"].tolist(), z["pred_expr"]))
        pears = [np.corrcoef(pred_spectra[i], ref[p])[0, 1]
                 for i, p in enumerate(perts)]
        sanity = float(np.mean(pears))
        print(f"  [sanity] forward replication vs archived pred_expr "
              f"mean pearson = {sanity:.4f} (应≈1)")
        assert sanity > 0.99, "前向复刻与归档预测不一致，停止"

        X_geom, err, geom_names = load_geometry(m, perts)
        X_sets = {
            "geom9": X_geom,
            "meanpool64": meanpool.astype(np.float64),
            "pertnode64": pertnode.astype(np.float64),
            "crossgene64": crossgene.astype(np.float64),
            "geom9+meanpool64": np.column_stack(
                [X_geom, meanpool.astype(np.float64)]),
            "geom9+pertnode64": np.column_stack(
                [X_geom, pertnode.astype(np.float64)]),
            "geom9+crossgene64": np.column_stack(
                [X_geom, crossgene.astype(np.float64)]),
        }
        X_sets_per_seed[m] = (X_sets, err)
        results["model_seed_results"][str(m)] = {
            "feature_dims": {k: int(v.shape[1]) for k, v in X_sets.items()},
            "geom_feature_names": geom_names,
            "sanity_pred_pearson_vs_archived": sanity,
            "internal_repr_dim": int(meanpool.shape[1])}

    adata.file.close()

    # ---- 探针：七特征组 × 2 模型种子 × 3 split 种子，两个样本 pass ----
    # pass_full（n=217，与 step7 完全相同的 split）：geom9 / meanpool64 /
    #   crossgene64 / 两个拼接；pass_pertnode（n=83，靶基因在表达空间内的子集，
    #   同一生成逻辑）：pertnode64 / geom9+pertnode64（geom9 子集重算，names 中
    #   必须显式含 geom9 作配对基线）
    gene_set = set(z0["pred_genes"].tolist())
    pert_sub = [i for i, g in enumerate(perts) if g in gene_set]
    splits_sub = step7_splits(len(pert_sub))
    for m in MODEL_SEEDS:
        X_sets, err = X_sets_per_seed[m]
        ms_res = results["model_seed_results"][str(m)]
        ms_res["sample"] = {
            "full_n": n, "pertnode_subset_n": len(pert_sub),
            "pertnode_subset_reason": "靶基因须在 5000 维表达空间内才有 "
                                      "cell-graph 节点（134/217 不在，"
                                      "pred≈ctrl，step7 已记录）",
            "pertnode_subset_split_protocol": "step7_splits(n_sub)，与 step7 "
                                              "同一生成逻辑（必要协议偏差，"
                                              "数据驱动）"}

        names_full = ["geom9", "meanpool64", "crossgene64",
                      "geom9+meanpool64", "geom9+crossgene64"]
        X_full = {k: X_sets[k] for k in names_full}
        entries_f, summary_f = probe_pass(X_full, err, names_full, splits, m)
        ms_res["split_results_full"] = entries_f

        names_sub = ["geom9", "pertnode64", "geom9+pertnode64"]
        X_sub = {k: X_sets[k][pert_sub] for k in names_sub}
        entries_s, summary_s = probe_pass(X_sub, err[pert_sub], names_sub,
                                          splits_sub, m)
        ms_res["split_results_pertnode"] = entries_s

        ms_res["summary"] = {}
        ms_res["summary"].update(summary_f)
        ms_res["summary"].update(summary_s)
        for name, s in ms_res["summary"].items():
            if name == "geom9":
                continue
            print(f"[summary seed{m}] {name} (n={s['sample_n']}): "
                  f"AUROC={s['auroc_mean']:.3f}±{s['auroc_std']:.3f}, "
                  f"diff vs geom9={s['diff_vs_geom9_mean']:+.3f} "
                  f"-> {s['criterion_verdict']}")

    # ---- 产物 ----
    out_json = os.path.join(REPO, "outputs", "internal_repr.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"[SAVED] {out_json}")
    make_figure(results, os.path.join(REPO, "outputs", "internal_repr.png"))
    return results


def make_figure(results, path):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6), sharey=True)
    labels = {"geom9": "output\ngeometry (9)",
              "meanpool64": "internal\nmean-pool (64)",
              "pertnode64": "internal\npert-node (64)*",
              "crossgene64": "internal\ncross-gene (64)",
              "geom9+meanpool64": "geom+\nmean-pool (73)",
              "geom9+pertnode64": "geom+\npert-node (73)*",
              "geom9+crossgene64": "geom+\ncross-gene (73)"}
    colors = {"geom9": "#7f7f7f", "meanpool64": "#1f77b4",
              "pertnode64": "#d62728", "crossgene64": "#2ca02c",
              "geom9+meanpool64": "#9467bd", "geom9+pertnode64": "#8c564b",
              "geom9+crossgene64": "#bcbd22"}
    # 每个特征集的 per-split 散点取自对应样本 pass 的结果
    entries_key = {k: "split_results_full" for k in FEATURE_SETS}
    for k in ("pertnode64", "geom9+pertnode64"):
        entries_key[k] = "split_results_pertnode"
    rng_j = np.random.default_rng(0)
    for ax, m in zip(axes, MODEL_SEEDS):
        ms = results["model_seed_results"][str(m)]
        x = np.arange(len(FEATURE_SETS))
        means = [ms["summary"][k]["auroc_mean"] for k in FEATURE_SETS]
        stds = [ms["summary"][k]["auroc_std"] for k in FEATURE_SETS]
        bars = ax.bar(x, means, 0.62, color=[colors[k] for k in FEATURE_SETS],
                      alpha=0.85, edgecolor="gray", lw=0.5)
        ax.errorbar(x, means, yerr=stds, fmt="none", ecolor="black",
                    elinewidth=1.0, capsize=3)
        for k, xi in zip(FEATURE_SETS, x):
            ens = ms[entries_key[k]]
            jit = rng_j.normal(0, 0.03, len(ens))
            ax.plot(xi + jit, [e["auroc"][k] for e in ens],
                    "o", ms=4, mfc="none", mec="black", alpha=0.7, ls="none")
        for xi, v in zip(x, means):
            ax.text(xi, v + 0.012, f"{v:.3f}", ha="center", fontsize=8)
        ax.axhline(0.5, color="red", ls=":", lw=1.1)
        ax.set_xticks(x)
        ax.set_xticklabels([labels[k] for k in FEATURE_SETS], fontsize=8,
                           rotation=28, ha="right")
        ax.set_ylim(0.35, 0.9)
        ax.set_title(f"GEARS model seed {m}", fontsize=10)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("test AUROC (30% held-out, y = worst-half error)")
    fig.suptitle("Internal hidden representation vs output geometry for "
                 "predicting GEARS self-error\n(K562 Replogle essential; "
                 "bars = mean over 3 split seeds, points = per split; "
                 "*perturbed-node readout evaluable on n=83 subset, see text)")
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")


if __name__ == "__main__":
    run_all()
