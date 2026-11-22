# -*- coding: utf-8 -*-
"""Geneformer V2-104M wave 推理: 每细胞 embedding(pred 侧 + true 侧)。

口径(2026-10-18 多模型波, 两类分治):
- embedding 类误差 = 模型自身空间 d(pred_emb, true_emb);
- pred_emb = ctrl 池(pool=32, seed=0, 与 scgenept 同池)删除扰动基因 token 后
  的每细胞嵌入; true_emb = 真实扰动细胞(<=32/pert, seed=777000+pi 子采样)
  的每细胞嵌入;
- 签名特征(pred_emb 范数/delta_emb 范数/弦距 ctrl 等)由分析端从本产物复算。

产物 <out>_predictions.npz(每细胞, 非仅均值):
  pred_perturbations / ctrl_emb_cells (pool,d) / ctrl_emb / ctrl_pool_indices
  / pred_emb_cells (n,pool,d) / pred_emb_mean (n,d)
  / true_emb_cells (n,cap,d, NaN 填充) / true_emb_mean (n,d) / true_n_cells
  / target_in_input_space / target_deleted / ctrl_mean_full / full_gene_space
  / output_type
断点: 每扰动写 <out>_checkpoint.npz(指纹=基因空间+扰动名单+pool+true_cells)。

用法(容器内):
  python -u /research/scripts/run_geneformer_wave.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/outputs_crossmodel/geneformer/geneformer_wave_k562
"""
import argparse
import json
import os
import pickle
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import torch

MODEL_DIR = "/research/models/geneformer/"
TOK_DICT = MODEL_DIR + "token_dictionary_gc104M.pkl"
NAME2ID = MODEL_DIR + "gene_name_id_dict_gc104M.pkl"
WEIGHTS = MODEL_DIR + "Geneformer-V2-104M/"
DEFAULT_H5AD = "/research/data/perturb_processed_k562.h5ad"
DEFAULT_LOCAL_NPZ = "/research/meta/pred_expr_seed0.npz"
SEED = 0
MAX_LEN = 2048


def row_mean(a):
    """(n, k, d) 每行均值, 全 NaN 行保持 NaN 不告警。"""
    s = np.nansum(a, axis=1)
    c = np.isfinite(a).sum(axis=1)
    out = np.full(a.shape[0:1] + a.shape[2:], np.nan, dtype=np.float32)
    np.divide(s, np.maximum(c, 1), out=out, where=(c > 0))
    return out


def sha256(path, chunk=1 << 22):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_dicts():
    with open(TOK_DICT, "rb") as f:
        tok = pickle.load(f)
    with open(NAME2ID, "rb") as f:
        n2i = pickle.load(f)
    return tok, n2i


def cell_sentences(X, var_genes, tok, n2i):
    """(n_cells, n_genes) -> token id 序列 list + gene->pos 映射 list。"""
    sentences, pos_maps = [], []
    for i in range(X.shape[0]):
        x = X[i]
        nz = np.where(x > 0)[0]
        order = nz[np.argsort(-x[nz])][:MAX_LEN - 2]
        toks, pmap = [], {}
        for gi in order:
            eid = n2i.get(var_genes[gi])
            t = tok.get(eid) if eid else None
            if t is None:
                continue
            pmap[var_genes[gi]] = len(toks)
            toks.append(t)
        sentences.append(toks)
        pos_maps.append(pmap)
    return sentences, pos_maps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=DEFAULT_H5AD)
    ap.add_argument("--local-npz", default=DEFAULT_LOCAL_NPZ)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", type=int, default=32)
    ap.add_argument("--true-cells", type=int, default=32)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SCT_THREADS", "48")))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    t_start = time.time()
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.set_num_threads(args.threads)
    device = torch.device("cpu")
    print(f"[cfg] h5ad={args.h5ad} pool={args.pool} true_cells={args.true_cells} "
          f"batch={args.batch} threads={args.threads}", flush=True)

    z = np.load(args.local_npz)
    holdout_genes = [str(g) for g in z["pred_perturbations"]]
    local_genes = [str(g) for g in z["pred_genes"]]

    adata = ad.read_h5ad(args.h5ad, backed="r")
    var_genes = adata.var["gene_name"].astype(str).tolist()
    assert var_genes == local_genes
    obs_cond = adata.obs["condition"].astype(str).values
    is_ctrl = obs_cond == "ctrl"
    n_ctrl = int(is_ctrl.sum())
    ctrl_mean_full = np.asarray(
        adata[is_ctrl].X.mean(axis=0)).ravel().astype(np.float32)
    print(f"[data] cells={adata.n_obs}, ctrl_cells={n_ctrl}", flush=True)

    # ---- 断点 ----
    CKPT = args.out + "_checkpoint.npz"
    resume_rows, pool_idx_ckpt = [], None
    if args.limit == 0 and os.path.exists(CKPT):
        cz = np.load(CKPT, allow_pickle=True)
        if (list(cz["ckpt_gene_space"]) == var_genes
                and list(cz["ckpt_perts"]) == holdout_genes
                and int(cz["ckpt_pool"]) == args.pool
                and int(cz["ckpt_true"]) == args.true_cells):
            resume_rows = [int(i) for i in cz["ckpt_rows"]]
            pool_idx_ckpt = cz["ckpt_pool_idx"]
            print(f"[ckpt] 恢复 {len(resume_rows)} 个扰动", flush=True)
        else:
            print("[ckpt] 指纹不匹配, 从头开始", flush=True)

    rng = np.random.default_rng(SEED)
    pool_idx = pool_idx_ckpt if pool_idx_ckpt is not None else \
        rng.choice(np.where(is_ctrl)[0], size=min(args.pool, n_ctrl),
                   replace=False)
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)

    tok, n2i = load_dicts()
    pool_sents, pool_pos = cell_sentences(ctrl_pool, var_genes, tok, n2i)
    target_in_space = np.array([g in var_genes for g in holdout_genes])
    target_has_token = np.array(
        [tok.get(n2i.get(g)) is not None if n2i.get(g) else False
         for g in holdout_genes])

    from transformers import BertModel
    model = BertModel.from_pretrained(WEIGHTS, output_attentions=False)
    model.to(device)
    model.eval()
    print(f"[model] Geneformer loaded", flush=True)

    def embed_batch(sents):
        maxlen = max(len(s) for s in sents) + 2
        ids = np.zeros((len(sents), maxlen), dtype=np.int64)
        att = np.zeros((len(sents), maxlen), dtype=np.int64)
        for r, s in enumerate(sents):
            ids[r, 0] = tok["<cls>"]
            ids[r, 1:1 + len(s)] = s
            ids[r, 1 + len(s)] = tok["<eos>"]
            att[r, :2 + len(s)] = 1
        with torch.no_grad():
            out = model(input_ids=torch.from_numpy(ids).to(device),
                        attention_mask=torch.from_numpy(att).to(device))
        h = out.last_hidden_state
        mask = torch.from_numpy(att).to(device).unsqueeze(-1).float()
        gene_mask = mask.clone()
        idx = att.sum(1) - 1
        gene_mask[torch.arange(len(sents)), torch.from_numpy(idx)] = 0
        gene_mask[:, 0] = 0
        emb = (h * gene_mask).sum(1) / gene_mask.sum(1).clamp(min=1)
        return emb.cpu().numpy().astype(np.float32)

    def embed_sents(sents):
        outs = []
        for b in range(0, len(sents), args.batch):
            outs.append(embed_batch(sents[b:b + args.batch]))
        return np.concatenate(outs, axis=0)

    def embed_pool_dropped(drop_gene=None):
        if drop_gene is None:
            return embed_sents(pool_sents)
        sents = []
        for i, s in enumerate(pool_sents):
            p = pool_pos[i].get(drop_gene, None)
            sents.append(s if p is None else s[:p] + s[p + 1:])
        return embed_sents(sents)

    # ---- 基线 ----
    t0 = time.time()
    if os.path.exists(CKPT + ".ctrl_emb.npy") and args.limit == 0:
        ctrl_emb_cells = np.load(CKPT + ".ctrl_emb.npy")
        print("[emb] ctrl baseline loaded from ckpt", flush=True)
    else:
        ctrl_emb_cells = embed_pool_dropped()
        if args.limit == 0:
            np.save(CKPT + ".ctrl_emb.npy", ctrl_emb_cells)
    ctrl_emb = ctrl_emb_cells.mean(axis=0)
    print(f"[emb] ctrl baseline {ctrl_emb_cells.shape} in {time.time()-t0:.0f}s",
          flush=True)

    # ---- true 细胞选择(每 pert 确定性, 跨模型同规则) ----
    cond_set = set(obs_cond)
    true_idx = {}
    for pi, g in enumerate(holdout_genes):
        cond = next((c for c in (f"ctrl+{g}", f"{g}+ctrl") if c in cond_set),
                    None)
        if cond is None:
            true_idx[pi] = np.array([], dtype=np.int64)
            continue
        pos = np.where(obs_cond == cond)[0]
        r = np.random.default_rng(777000 + pi)
        k = min(args.true_cells, len(pos))
        true_idx[pi] = np.sort(r.choice(pos, size=k, replace=False))
    print(f"[data] true cells/pert: min={min(len(v) for v in true_idx.values())} "
          f"max={max(len(v) for v in true_idx.values())}", flush=True)

    n_all = len(holdout_genes)
    D = ctrl_emb.shape[0]
    cap = args.true_cells
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    pred_cells = np.full((n_all, args.pool, D), np.nan, dtype=np.float32)
    true_cells = np.full((n_all, cap, D), np.nan, dtype=np.float32)
    true_n = np.zeros(n_all, dtype=np.int32)
    target_deleted = np.zeros(n_all, dtype=bool)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        rr = np.array(resume_rows)
        pred_cells[rr] = cz["ckpt_pred_cells"]
        true_cells[rr] = cz["ckpt_true_cells"]
        true_n[rr] = cz["ckpt_true_n"]
        target_deleted[rr] = cz["ckpt_deleted"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for pi in run_rows:
        g = holdout_genes[pi]
        ti = time.time()
        deleted = bool(target_has_token[pi]) and any(
            g in pool_pos[i] for i in range(len(pool_pos)))
        pred_cells[pi] = embed_pool_dropped(g) if deleted \
            else ctrl_emb_cells.copy()
        target_deleted[pi] = deleted
        idx = true_idx[pi]
        if len(idx):
            Xs = np.asarray(adata[idx].X.todense()).astype(np.float32)
            sents, _ = cell_sentences(Xs, var_genes, tok, n2i)
            true_cells[pi, :len(idx)] = embed_sents(sents)
            true_n[pi] = len(idx)
        done_rows.append(pi)
        if not args.limit:
            rr = np.array(done_rows)
            np.savez(CKPT, ckpt_rows=rr,
                     ckpt_pred_cells=pred_cells[rr],
                     ckpt_true_cells=true_cells[rr],
                     ckpt_true_n=true_n[rr],
                     ckpt_deleted=target_deleted[rr],
                     ckpt_perts=np.array(holdout_genes),
                     ckpt_gene_space=np.array(var_genes),
                     ckpt_pool=args.pool, ckpt_true=args.true_cells,
                     ckpt_pool_idx=np.asarray(pool_idx))
        n_done = len(done_rows)
        if n_done <= 5 or n_done % 20 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} deleted={deleted} "
                  f"{time.time()-ti:.1f}s (eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done {n_all} perts in {t_infer/60:.1f} min", flush=True)
    adata.file.close()

    # ---- 产物 ----
    pred_mean = row_mean(pred_cells)
    true_mean = row_mean(true_cells)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez_compressed(
        args.out + "_predictions.npz",
        pred_perturbations=np.array(holdout_genes),
        pred_emb_cells=pred_cells, pred_emb_mean=pred_mean,
        pred_delta_emb=pred_mean - ctrl_emb[None, :],
        true_emb_cells=true_cells, true_emb_mean=true_mean,
        true_n_cells=true_n,
        ctrl_emb=ctrl_emb, ctrl_emb_cells=ctrl_emb_cells,
        ctrl_pool_indices=pool_idx.astype(np.int64),
        target_in_input_space=target_in_space,
        target_deleted=target_deleted,
        ctrl_mean_full=ctrl_mean_full,
        full_gene_space=np.array(var_genes),
        output_type=np.array(
            ["per_cell_embeddings (pred=ctrl pool w/ target token deleted; "
             "true=real pert cells); embedding-class model, pearson N/A"]),
    )
    labels = []
    for pi, g in enumerate(holdout_genes):
        pm, tm = pred_mean[pi], true_mean[pi]
        labels.append({
            "perturbation": g, "pearson_delta": np.nan,
            "pearson_delta_all": np.nan, "n_cells": int(true_n[pi]),
            "target_in_input_space": bool(target_in_space[pi]),
            "target_deleted": bool(target_deleted[pi]),
            "pred_emb_norm": float(np.linalg.norm(pm)),
            "true_emb_norm": float(np.linalg.norm(tm)),
            "delta_emb_norm": float(np.linalg.norm(pm - ctrl_emb)),
            "err_emb_dist": float(np.linalg.norm(pm - tm)),
        })
    pd.DataFrame(labels).to_csv(args.out + "_labels.csv", index=False)

    meta = {
        "task": "Geneformer V2-104M wave: per-cell embeddings K562/RPE1",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {"name": "Geneformer V2-104M",
                  "weights_sha256": sha256(WEIGHTS + "model.safetensors"),
                  "inference_protocol": "in silico token deletion; cell emb = "
                      "mean gene-token last_hidden_state; per-cell kept",
                  "output_type": "cell embedding (768-d)"},
        "data": {"h5ad": args.h5ad, "heldout_source": args.local_npz,
                 "ctrl_pool_size": int(len(ctrl_pool)),
                 "true_cells_cap": cap,
                 "true_sampling": "default_rng(777000+pi) w/o replacement"},
        "runtime": {"command": " ".join(sys.argv),
                    "threads": args.threads,
                    "infer_seconds_total": t_infer,
                    "total_seconds": time.time() - t_start},
    }
    with open(args.out + "_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    if not args.limit and os.path.exists(CKPT):
        os.remove(CKPT)
        p = CKPT + ".ctrl_emb.npy"
        if os.path.exists(p):
            os.remove(p)
    print(f"[out] -> {args.out}_predictions.npz (total "
          f"{time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
