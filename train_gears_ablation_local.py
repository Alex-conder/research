# -*- coding: utf-8 -*-
"""GEARS 数据覆盖消融（A/B/C/D）：训练 + held-out 评估 + 兼容产物。
Windows/本机 GPU 版（RTX 5060 Laptop 8GB），由 train_gears_ablation.py
（服务器 / Linux 版）适配而来，协议完全一致：

- cell-gears 0.1.2（打补丁版，与基线同码）、hidden=64、batch=32、split_seed=0、
  holdout_frac=0.2、GEARS "simulation" split + only_test_set_perts。
- A: train k562 -> eval k562；B: train rpe1 -> eval rpe1
- C: train k562+rpe1 混合 -> 两者都测（同一 held-out 名单 = 两数据集名单并集）
- D: train k562 -> eval rpe1（零样本；复用 A 权重，仅评估）

与 Linux 版的差异仅在本机适配层：
- 路径全部指向本机项目（data/gears_demo/...，pyg 缓存与共表达 CSV 已就位）；
- dataset 目录直接用官方 demo 目录（replogle_k562_essential / replogle_rpe1_essential），
  不做 symlink；mixed 用 data/gears_demo/mixed/；
- resource/uname 等 POSIX 调用替换为 psutil/platform（loadavg 护栏在本机恒为 0，
  保留代码不动）；
- 每 epoch 日志附加 GPU 峰值显存。

用法：
  python train_gears_ablation_local.py --variant A_k562 --epochs 10
  python train_gears_ablation_local.py --variant D_zeroshot \
      --weights outputs_ablation/A_k562/gears_model/model.pt
  # 测速（不出产物）：
  python train_gears_ablation_local.py --variant A_k562 --max-steps 50
"""
import argparse
import hashlib
import json
import os
import pickle
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import torch

SPLIT_SEED = 0
HOLDOUT_FRAC = 0.2
HIDDEN = 64
TOP_N_DEG = 50
BATCH = 32
TEST_BATCH = 128
LOADAVG_GUARD = 200.0

PROJECT = os.path.dirname(os.path.abspath(__file__))
GEARS_DEMO = os.path.join(PROJECT, "data", "gears_demo")
H5AD = {"k562": os.path.join(GEARS_DEMO, "replogle_k562_essential",
                             "perturb_processed.h5ad"),
        "rpe1": os.path.join(GEARS_DEMO, "replogle_rpe1_essential",
                             "perturb_processed.h5ad")}
META_NPZ = {"k562": os.path.join(PROJECT, "outputs", "pred_expr_seed0.npz"),
            "rpe1": os.path.join(PROJECT, "outputs_rpe1",
                                 "pred_expr_seed0.npz")}
# PertData(AB_ROOT) 需要 gene2go_all.pkl / essential_all_data_pert_genes.pkl，
# 两者已在 data/gears_demo/ 下；共表达 CSV 写到
# <AB_ROOT>/<dataset_name>/simulation_0_...csv（目录已存在）。
AB_DATA = GEARS_DEMO
AB_ROOT = GEARS_DEMO
DS_DIR = {"k562": os.path.join(GEARS_DEMO, "replogle_k562_essential"),
          "rpe1": os.path.join(GEARS_DEMO, "replogle_rpe1_essential"),
          "mixed": os.path.join(GEARS_DEMO, "mixed")}
DEFAULT_OUT = os.path.join(PROJECT, "outputs_ablation")
VARIANTS = {
    "A_k562":     {"train": ["k562"],       "eval": ["k562"]},
    "B_rpe1":     {"train": ["rpe1"],       "eval": ["rpe1"]},
    "C_mixed":    {"train": ["k562", "rpe1"], "eval": ["k562", "rpe1"]},
    "D_zeroshot": {"train": ["k562"],       "eval": ["rpe1"]},
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def loadavg():
    try:
        with open("/proc/loadavg") as f:
            return float(f.read().split()[0])
    except Exception:
        return 0.0


def wait_for_load():
    while loadavg() > LOADAVG_GUARD:
        log(f"loadavg {loadavg():.0f} > {LOADAVG_GUARD:.0f}, wait 120s")
        time.sleep(120)


def vm_hwm_gb():
    """Linux 用 rusage 峰值 RSS；Windows 退回 psutil 当前 RSS。"""
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    except ImportError:
        try:
            import psutil
            return psutil.Process().memory_info().rss / 1e9
        except Exception:
            return 0.0


def gpu_peak_gb(reset=True):
    if not torch.cuda.is_available():
        return 0.0
    v = torch.cuda.max_memory_allocated() / 1e9
    if reset:
        torch.cuda.reset_peak_memory_stats()
    return v


def compute_holdout(adata):
    """与 step1_run.py 完全一致的留出名单协议（在 PertData.load 的 GO 过滤后
    的 adata 上调用——load 会剔除扰动基因不在 GO 图的条件）。返回 holdout。"""
    conditions = adata.obs.condition.unique().tolist()
    singles = [c for c in conditions if ("ctrl" in c and c != "ctrl")]
    genes = sorted({c.replace("ctrl+", "").replace("+ctrl", "") for c in singles})
    genes = [g for g in genes if "_" not in g]
    rng = np.random.default_rng(SPLIT_SEED)
    holdout = rng.choice(genes, size=int(HOLDOUT_FRAC * len(genes)),
                         replace=False).tolist()
    return genes, [str(g) for g in holdout]


def apply_go_filter(adata):
    """复刻 PertData.load 内的 GO 图条件过滤（需要 pert_names）。"""
    from gears.utils import filter_pert_in_go
    with open(os.path.join(AB_ROOT, "essential_all_data_pert_genes.pkl"),
              "rb") as f:
        essential = pickle.load(f)
    with open(os.path.join(AB_ROOT, "gene2go_all.pkl"), "rb") as f:
        gene2go = pickle.load(f)
    g2g = {i: gene2go[i] for i in essential if i in gene2go}
    pert_names = np.unique(list(g2g.keys()))
    keep = adata.obs.condition.apply(lambda x: filter_pert_in_go(x, pert_names))
    return adata[keep.index.values[keep.values]]


def ensure_dataset_dir(d):
    """本机：数据集目录 = 官方 demo 目录（已含 perturb_processed.h5ad、
    data_pyg 缓存、simulation_0 共表达 CSV、splits/），直接校验存在即可。"""
    ddir = DS_DIR[d]
    link = os.path.join(ddir, "perturb_processed.h5ad")
    if not os.path.exists(link):
        raise SystemExit(f"dataset dir missing h5ad: {link}")
    return ddir


def ensure_mixed_h5ad():
    """C 版本：两数据集拼接（唯一开发点）。

    batch effect 处理口径（如实记录，不过度工程化）：
    - 两数据集各自已完成 GEARS 标准预处理（各自 log-normalize）；
    - 拼接保留 obs['dataset'] 标签供分析留痕；同名列条件（如 ctrl+GENE 在
      K562 与 RPE1 都存在）按 GEARS 原生行为合并为一个 condition（ctrl 池
      亦合并），模型条件空间 = 并集；
    - 不引入额外 covariate 机制（GEARS 0.1.2 无此接口），batch effect 吸收
      能力 = 模型原生上限。此局限写入 meta 与 verdict。
    """
    import scanpy as sc
    import anndata as ad
    ddir = DS_DIR["mixed"]
    out = os.path.join(ddir, "perturb_processed.h5ad")
    os.makedirs(ddir, exist_ok=True)
    if os.path.exists(out):
        log(f"[data] mixed h5ad exists: {out}")
        return ddir
    log("[data] building mixed k562+rpe1 h5ad ...")
    a_k = sc.read_h5ad(H5AD["k562"])
    a_r = sc.read_h5ad(H5AD["rpe1"])
    g_k = a_k.var["gene_name"].astype(str).tolist()
    g_r = a_r.var["gene_name"].astype(str).tolist()
    if g_k != g_r:
        # 两数据集 HVG 选择不同（实测交集 2832/5000）：按共有基因交集对齐到
        # k562 基因序。位置式重索引（保持 var_names=Ensembl 对齐）。
        common = [g for g in g_k if g in set(g_r)]
        if len(common) < 0.5 * len(g_k):
            raise SystemExit(f"gene space overlap too small: {len(common)}")
        a_k = a_k[:, [g_k.index(g) for g in common]]
        a_r = a_r[:, [g_r.index(g) for g in common]]
        log(f"[data] gene spaces aligned to common {len(common)} genes "
            f"(k562 {len(g_k)}, rpe1 {len(g_r)})")
    # 注：uns['non_zeros_gene_idx'] 的取值是"保留全部基因"的列排列，对 loss
    # 数值是 no-op（对全列求和与顺序无关），两数据集各自的排列直接合并即可；
    # uns['rank_genes_groups_cov_all'] 以 condition_name（K562_/rpe1_ 前缀，
    # 两数据集不相交）为键、Ensembl ID 为值，与基因序无关，可直接合并。
    a_k.obs["dataset"] = "k562"
    a_r.obs["dataset"] = "rpe1"
    # 两数据集 obs_names 是 10x 细胞条形码，跨数据集必然重复；PertData.load
    # 的 GO 过滤用 adata[obs_names] 取值（anndata 要求索引唯一），重复会直接
    # 炸 InvalidIndexError。拼数据集前缀保证唯一（条件逻辑只用 obs.condition）。
    a_k.obs_names = [f"{n}-k562" for n in a_k.obs_names]
    a_r.obs_names = [f"{n}-rpe1" for n in a_r.obs_names]
    # uns 合并：rank_genes_groups_cov_all / non_zeros_gene_idx 取并集
    uns = {}
    for key in ("rank_genes_groups_cov_all", "non_zeros_gene_idx"):
        merged = {}
        for a in (a_k, a_r):
            v = a.uns.get(key)
            if isinstance(v, dict):
                merged.update(v)
        uns[key] = merged
    mixed = ad.concat([a_k, a_r], join="outer", label=None, index_unique=None)
    for k, v in uns.items():
        mixed.uns[k] = v
    log(f"[data] mixed: cells={mixed.n_obs}, genes={mixed.n_vars}, "
        f"conditions={mixed.obs.condition.nunique()}")
    mixed.write(out, compression=None)
    log(f"[data] mixed h5ad -> {out} ({os.path.getsize(out)/1e9:.1f}GB)")
    return ddir


def train_dir_name(variant, train_ds):
    return "mixed" if len(train_ds) > 1 else train_ds[0]


def build_pertdata(variant, train_ds):
    """加载（或建缓存）训练用 PertData，并按协议切分。"""
    from gears import PertData
    ds_key = train_dir_name(variant, train_ds)
    if ds_key == "mixed":
        ddir = ensure_mixed_h5ad()
    else:
        ddir = ensure_dataset_dir(ds_key)
    t0 = time.time()
    pert_data = PertData(AB_ROOT)  # gene2go_all.pkl 已在 data/gears_demo
    cache = os.path.join(ddir, "data_pyg", "cell_graphs.pkl")
    if os.path.exists(cache):
        log(f"[data] pyg cache hit: {cache} ({os.path.getsize(cache)/1e9:.1f}GB)")
    else:
        log("[data] pyg cache missing, will build (one-off, slow) ...")
    pert_data.load(data_path=ddir)
    log(f"[data] loaded {ds_key}: cells={pert_data.adata.n_obs}, "
        f"conds={pert_data.adata.obs.condition.nunique()} ({time.time()-t0:.0f}s)")
    return pert_data, ds_key


def get_holdouts(eval_ds):
    """按 step1 协议从原始 h5ad 计算各数据集留出名单，并与 meta npz 对齐校验。"""
    import scanpy as sc
    holdout = {}
    for d in eval_ds:
        adata = sc.read_h5ad(H5AD[d], backed="r")
        _, ho = compute_holdout(apply_go_filter(adata.to_memory()))
        adata.file.close()
        z = np.load(META_NPZ[d], allow_pickle=True)
        ref = [str(g) for g in z["pred_perturbations"]]
        in_ref = [g for g in ho]           # meta npz = in_graph 顺序
        assert set(ref) == set(in_ref), (
            f"{d} holdout mismatch vs {META_NPZ[d]}: "
            f"|symdiff|={len(set(ref) ^ set(in_ref))}")
        holdout[d] = ho
        log(f"[split] {d}: genes pool ok, holdout={len(ho)} (aligned with meta npz)")
    return holdout


def do_train(model, out_dir, epochs, resume=True, max_steps=0):
    """带 checkpoint 的 GEARS 训练循环（复刻 gears.train，补 ckpt/计时/显存日志）。"""
    from copy import deepcopy
    from gears.inference import evaluate, compute_metrics
    from gears.utils import loss_fct, uncertainty_loss_fct
    from torch import optim
    from torch.optim.lr_scheduler import StepLR

    train_loader = model.dataloader["train_loader"]
    val_loader = model.dataloader["val_loader"]
    model.model = model.model.to(model.device)
    optimizer = optim.Adam(model.model.parameters(), lr=1e-3, weight_decay=5e-4)
    scheduler = StepLR(optimizer, step_size=1, gamma=0.5)

    ckpt_last = os.path.join(out_dir, "ckpt_last.pt")
    ckpt_best = os.path.join(out_dir, "ckpt_best.pt")
    hist_path = os.path.join(out_dir, "history.jsonl")
    start_epoch = 0
    min_val = np.inf
    best_model = None
    if resume and os.path.exists(ckpt_last):
        ck = torch.load(ckpt_last, map_location=model.device, weights_only=False)
        model.model.load_state_dict(ck["model_state"])
        optimizer.load_state_dict(ck["opt_state"])
        scheduler.load_state_dict(ck["sched_state"])
        start_epoch = ck["epoch"] + 1
        min_val = ck["best_val"]
        if os.path.exists(ckpt_best):
            best_model = deepcopy(model.model)
        log(f"[ckpt] resume from epoch {start_epoch} (best_val_de={min_val:.4f})")

    history = []
    for epoch in range(start_epoch, epochs):
        wait_for_load()
        t_ep = time.time()
        model.model.train()
        n_steps = 0
        for step, batch in enumerate(train_loader):
            batch.to(model.device)
            optimizer.zero_grad()
            y = batch.y
            if model.config["uncertainty"]:
                pred, logvar = model.model(batch)
                loss = uncertainty_loss_fct(
                    pred, logvar, y, batch.pert, reg=model.config["uncertainty_reg"],
                    ctrl=model.ctrl_expression, dict_filter=model.dict_filter,
                    direction_lambda=model.config["direction_lambda"])
            else:
                pred = model.model(batch)
                loss = loss_fct(pred, y, batch.pert, ctrl=model.ctrl_expression,
                                dict_filter=model.dict_filter,
                                direction_lambda=model.config["direction_lambda"])
            loss.backward()
            torch.nn.utils.clip_grad_value_(model.model.parameters(), 1.0)
            optimizer.step()
            if step % 50 == 0:
                log(f"Epoch {epoch+1} Step {step+1} Train Loss: {loss.item():.4f}")
            n_steps = step + 1
            if max_steps and n_steps >= max_steps:
                break
        scheduler.step()
        train_res = evaluate(train_loader, model.model, model.config["uncertainty"],
                             model.device, max_batches=500)
        val_res = evaluate(val_loader, model.model, model.config["uncertainty"],
                           model.device)
        train_metrics, _ = compute_metrics(train_res)
        val_metrics, _ = compute_metrics(val_res)
        ep_sec = time.time() - t_ep
        log(f"Epoch {epoch+1}: train MSE {train_metrics['mse']:.4f} | "
            f"val MSE {val_metrics['mse']:.4f} | val DE20 MSE {val_metrics['mse_de']:.4f} | "
            f"{ep_sec:.0f}s ({ep_sec/max(n_steps,1):.2f}s/step, {n_steps} steps) | "
            f"HWM {vm_hwm_gb():.1f}GB | GPU peak {gpu_peak_gb():.2f}GB")
        history.append({"epoch": epoch + 1, "train_mse": train_metrics["mse"],
                        "val_mse": val_metrics["mse"],
                        "val_mse_de": val_metrics["mse_de"],
                        "sec": ep_sec, "steps": n_steps,
                        "vm_hwm_gb": round(vm_hwm_gb(), 2),
                        "gpu_peak_gb": round(gpu_peak_gb(reset=False), 2)})
        with open(hist_path, "a") as f:
            f.write(json.dumps(history[-1]) + "\n")
        is_best = val_metrics["mse_de"] < min_val
        if is_best:
            min_val = val_metrics["mse_de"]
            best_model = deepcopy(model.model)
        torch.save({"epoch": epoch, "model_state": model.model.state_dict(),
                    "opt_state": optimizer.state_dict(),
                    "sched_state": scheduler.state_dict(),
                    "best_val": min_val,
                    "val_mse_de": val_metrics["mse_de"]}, ckpt_last)
        if is_best:
            torch.save({"epoch": epoch, "model_state": model.model.state_dict(),
                        "opt_state": optimizer.state_dict(),
                        "sched_state": scheduler.state_dict(),
                        "best_val": min_val}, ckpt_best)
        if max_steps and n_steps >= max_steps:
            log(f"[smoke] max-steps {max_steps} reached, stop training early")
            break

    if best_model is not None:
        model.best_model = best_model
    else:
        model.best_model = deepcopy(model.model)
    return history


def predict_holdout(model, holdout_genes):
    in_graph = [g for g in holdout_genes if g in model.pert_list]
    skipped = [g for g in holdout_genes if g not in model.pert_list]
    if skipped:
        log(f"[eval] {len(skipped)}/{len(holdout_genes)} not in pert graph -> NaN")
    pred = model.predict([[g] for g in in_graph]) if in_graph else {}
    pred_map = {g: np.asarray(pred[g]).ravel() for g in in_graph}
    return pred_map, in_graph, skipped


def evaluate_and_save(model, variant, eval_ds, holdout, out_prefix, meta_extra):
    """对每个 eval 数据集计算 held-out 预测 + pearson_delta，写兼容产物。

    基因空间处理：模型预测空间 = model.gene_list（A/B: 5000；C: 共有 2832；
    D: k562 的 5000）。eval 数据集 h5ad 的 var 序与其不一致时，按基因名映射到
    模型空间；D（k562->rpe1）中 rpe1 缺失的基因整列剔除（存 2832 共有空间），
    所有指标与特征均在存储空间上计算，行内无 NaN。
    """
    import scanpy as sc
    gene_list = list(model.gene_list)

    # ---- 第一遍：确定每个 eval 数据集的基因映射 + ctrl 均值 ----
    ds_info = {}
    for d in eval_ds:
        adata = sc.read_h5ad(H5AD[d], backed="r")
        var_genes = adata.var["gene_name"].astype(str).tolist()
        obs_cond = adata.obs["condition"].astype(str).values
        is_ctrl = obs_cond == "ctrl"
        ctrl_full = np.asarray(adata[is_ctrl].X.mean(axis=0)).ravel()
        if var_genes == gene_list:
            remap = None                       # 同序，直接用
        else:
            pos = {g: i for i, g in enumerate(var_genes)}
            keep = [g for g in gene_list if g in pos]       # 模型序
            m_idx = [i for i, g in enumerate(gene_list)
                     if g in pos]                         # 模型空间位置
            e_idx = [pos[g] for g in keep]                 # eval 空间位置
            if variant == "C_mixed":
                assert len(keep) == len(gene_list), (
                    f"C: {d} 未覆盖模型全部基因空间（{len(keep)}/{len(gene_list)}）")
            # remap = (eval空间索引, 存储基因名, 模型空间索引)
            # 注意两边的基因序不同（mixed 合并序 ≠ k562/rpe1 发布序），
            # 模型侧预测 p 必须用 m_idx，eval 侧真值/ctrl 用 e_idx，
            # 不能共用一个索引（D 之前没炸纯属两数据集 var 序碰巧一致）。
            remap = (e_idx, keep, m_idx)
        ds_info[d] = {"var_genes": var_genes, "obs_cond": obs_cond,
                      "cond_set": set(obs_cond), "ctrl_full": ctrl_full,
                      "remap": remap, "n_cells_total": int(adata.n_obs)}
        adata.file.close()
    stored_genes = gene_list
    for d in eval_ds:
        r = ds_info[d]["remap"]
        if r is not None:
            stored_genes = r[1]                # D: 共有子空间
            break
    n_genes = len(stored_genes)

    # ---- 第二遍：预测（每数据集一次）+ 逐行真值/指标 ----
    rows, pred_rows = [], []
    pred_cache = {}
    n_rows = sum(len(holdout[d]) for d in eval_ds)
    true_mat = np.full((n_rows, n_genes), np.nan, np.float32)
    ctrl_mat = np.full((n_rows, n_genes), np.nan, np.float32)
    ri = 0
    for d in eval_ds:
        info = ds_info[d]
        remap = info["remap"]
        ctrl_d = (info["ctrl_full"] if remap is None
                  else info["ctrl_full"][remap[0]])
        adata = sc.read_h5ad(H5AD[d], backed="r")
        obs_cond = info["obs_cond"]
        cond_set = info["cond_set"]
        if d not in pred_cache:
            pred_cache[d] = predict_holdout(model, holdout[d])
        pred_map, in_graph, skipped = pred_cache[d]
        for j, g in enumerate(holdout[d]):
            cond = next((c for c in (f"ctrl+{g}", f"{g}+ctrl")
                         if c in cond_set), None)
            ctrl_mat[ri + j] = ctrl_d.astype(np.float32)
            if g in pred_map:
                p = pred_map[g]
                if remap is not None:
                    p = p[remap[2]]
                if cond is None:
                    rows.append({"perturbation": g, "dataset": d,
                                 "pearson_delta": np.nan,
                                 "pearson_delta_all": np.nan,
                                 "n_cells": 0, "in_pert_graph": True})
                    pred_rows.append(np.full(n_genes, np.nan, np.float32))
                    continue
                sub = adata[obs_cond == cond]
                tm = np.asarray(sub.X.mean(axis=0)).ravel()
                if remap is not None:
                    tm = tm[remap[0]]
                true_mat[ri + j] = tm.astype(np.float32)
                p_delta, t_delta = p - ctrl_d, tm - ctrl_d
                r_all = float(np.corrcoef(p_delta, t_delta)[0, 1])
                top = np.argsort(-np.abs(t_delta))[:TOP_N_DEG]
                r_top = float(np.corrcoef(p_delta[top], t_delta[top])[0, 1])
                rows.append({"perturbation": g, "dataset": d,
                             "pearson_delta": r_top, "pearson_delta_all": r_all,
                             "n_cells": int(sub.n_obs), "in_pert_graph": True})
                pred_rows.append(p.astype(np.float32))
            else:
                n = int((obs_cond == cond).sum()) if cond is not None else 0
                rows.append({"perturbation": g, "dataset": d,
                             "pearson_delta": np.nan, "pearson_delta_all": np.nan,
                             "n_cells": n, "in_pert_graph": False})
                pred_rows.append(np.full(n_genes, np.nan, np.float32))
        adata.file.close()
        ri += len(holdout[d])

    lab = pd.DataFrame(rows)
    pred_mat = np.stack(pred_rows)

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    np.savez_compressed(
        out_prefix + "_predictions.npz",
        pred_expr=pred_mat,
        pred_expr_full=pred_mat,
        pred_perturbations=np.array([r["perturbation"] for r in rows]),
        row_dataset=np.array([r["dataset"] for r in rows]),
        pred_genes=np.array(stored_genes),
        model_gene_space=np.array(gene_list),
        true_expr=true_mat,
        true_expr_full=true_mat,
        ctrl_mean_matched=ctrl_mat,
        ctrl_mean_full=ctrl_mat,
        full_gene_space=np.array(stored_genes),
        target_in_input_space=np.array([r["in_pert_graph"] for r in rows]),
    )
    lab.to_csv(out_prefix + "_labels.csv", index=False)

    valid = lab.dropna(subset=["pearson_delta"])
    summary = {}
    for d in eval_ds:
        vd = valid[valid.dataset == d]
        summary[d] = {"n": int(len(vd)),
                      "pearson_delta_mean": float(vd.pearson_delta.mean()),
                      "pearson_delta_median": float(vd.pearson_delta.median()),
                      "pearson_delta_all_mean": float(vd.pearson_delta_all.mean())}
    meta = {
        "task": f"GEARS data-coverage ablation variant {variant}",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "variant": variant,
        "protocol": {"split_seed": SPLIT_SEED, "holdout_frac": HOLDOUT_FRAC,
                     "hidden": HIDDEN, "batch": BATCH, "top_n_deg": TOP_N_DEG,
                     "gears_split": "simulation/only_test_set_perts",
                     "baseline_ref": "local step1_run.py cell-gears 0.1.2 (patched)"},
        "data": {d: {"h5ad": H5AD[d], "h5ad_sha256": sha256(H5AD[d]),
                     "heldout_source": META_NPZ[d]} for d in eval_ds},
        "mixed_batch_effect_note": (
            "C 版本：k562+rpe1 直接拼接；同名列条件按 GEARS 原生行为合并；"
            "obs['dataset'] 留痕；无额外 batch covariate（GEARS 0.1.2 无接口）"),
        "gene_space_note": (
            f"模型预测空间 {len(gene_list)} 基因；存储/评估空间 {n_genes} 基因"
            + ("（k562->rpe1 共有子空间，D 版本）" if n_genes != len(gene_list)
               else "")),
        "eval_protocol": {
            "label_definition": "pearson_delta = Pearson(pred_delta, true_delta) "
                                f"on top-{TOP_N_DEG} |true_delta| genes",
            "true_mean": "per-dataset condition cell mean",
            "ctrl_mean": "per-dataset ctrl cell mean"},
        "runtime": meta_extra,
        "host": platform.platform(),
        "summary": summary,
    }
    with open(out_prefix + "_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    log(f"[eval] summary={json.dumps(summary)}")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True, choices=list(VARIANTS))
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--out-root", default=DEFAULT_OUT)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-steps", type=int, default=0,
                    help="烟雾测速：训练 N step 即停，不出产物")
    ap.add_argument("--resume", type=int, default=1)
    ap.add_argument("--weights", default="",
                    help="D_zeroshot 用：A 版本 model.pt 路径")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    t_start = time.time()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.set_num_threads(args.threads)
    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    os.environ.setdefault("MKL_NUM_THREADS", str(args.threads))
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"[cfg] variant={args.variant} epochs={args.epochs} threads={args.threads} "
        f"device={device} torch={torch.__version__}")
    if device == "cuda":
        log(f"[cfg] GPU: {torch.cuda.get_device_name(0)} | "
            f"total {torch.cuda.get_device_properties(0).total_memory/1e9:.1f}GB")

    cfg = VARIANTS[args.variant]
    train_ds, eval_ds = cfg["train"], cfg["eval"]
    out_dir = os.path.join(args.out_root, args.variant)
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = os.path.join(out_dir, args.variant)

    from gears import GEARS

    holdout = get_holdouts(sorted(set(train_ds) | set(eval_ds)))

    if args.variant == "D_zeroshot":
        # 用训练数据集（k562）的 PertData 构建模型对象，加载 A 权重，
        # 在 rpe1 held-out 上预测（零样本；pert 图覆盖限制如实记录）。
        pert_data, ds_key = build_pertdata(args.variant, train_ds)
        test_perts = "_".join(holdout[train_ds[0]])
        pert_data.prepare_split(split="simulation", seed=SPLIT_SEED,
                                test_perts=test_perts,
                                only_test_set_perts=True)
        log(f"[split] {ds_key}: " +
            str({k: len(v) for k, v in pert_data.set2conditions.items()}))
        pert_data.get_dataloader(batch_size=BATCH, test_batch_size=TEST_BATCH)
        model = GEARS(pert_data, device=device)
        model.model_initialize(hidden_size=HIDDEN)
        wpath = args.weights or os.path.join(
            DEFAULT_OUT, "A_k562", "gears_model", "model.pt")
        sd = torch.load(wpath, map_location=device, weights_only=True)
        model.model.load_state_dict(sd)
        model.best_model = model.model
        log(f"[D] loaded A weights from {wpath}")
        summary = evaluate_and_save(
            model, args.variant, eval_ds, holdout, out_prefix,
            {"command": " ".join(sys.argv), "weights": wpath,
             "total_seconds": time.time() - t_start})
        log(f"[out] -> {out_prefix}_* (total {time.time()-t_start:.0f}s)")
        return

    pert_data, ds_key = build_pertdata(args.variant, train_ds)
    if args.variant == "C_mixed":
        test_perts = "_".join(sorted(set(holdout["k562"]) | set(holdout["rpe1"])))
    else:
        test_perts = "_".join(holdout[train_ds[0]])
    pert_data.prepare_split(split="simulation", seed=SPLIT_SEED,
                            test_perts=test_perts,
                            only_test_set_perts=True)
    log(f"[split] {ds_key}: " +
        str({k: len(v) for k, v in pert_data.set2conditions.items()}))
    pert_data.get_dataloader(batch_size=BATCH, test_batch_size=TEST_BATCH)

    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=HIDDEN)
    log(f"[model] initialized hidden={HIDDEN} uncertainty={model.config['uncertainty']}")

    history = do_train(model, out_dir, args.epochs,
                       resume=bool(args.resume), max_steps=args.max_steps)
    if args.max_steps:
        log(f"[smoke] step-rate done: {history[-1]['steps']} steps in "
            f"{history[-1]['sec']:.0f}s -> {history[-1]['sec']/history[-1]['steps']:.2f}s/step "
            f"| GPU peak {history[-1]['gpu_peak_gb']:.2f}GB")
        return

    model.save_model(os.path.join(out_dir, "gears_model"))
    log(f"[save] model -> {out_dir}/gears_model")
    summary = evaluate_and_save(
        model, args.variant, eval_ds, holdout, out_prefix,
        {"command": " ".join(sys.argv),
         "epochs": args.epochs,
         "epoch_seconds": [h["sec"] for h in history],
         "vm_hwm_gb": round(vm_hwm_gb(), 2),
         "total_seconds": time.time() - t_start})
    log(f"[out] -> {out_prefix}_* (total {time.time()-t_start:.0f}s)")


if __name__ == "__main__":
    main()
