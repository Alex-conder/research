# -*- coding: utf-8 -*-
"""Geneformer (ctheodoris/Geneformer, V2-104M) 对 Replogle essential held-out
扰动的零样本 in silico deletion 推理（K562）。输出类型：细胞嵌入（非表达谱，如实标注）。

协议：
- 基因名（h5ad var gene_name）→ Ensembl ID（gene_name_id_dict_gc104M.pkl）→
  token（token_dictionary_gc104M.pkl）。
- 每个 ctrl 池细胞（pool=32, seed=0）按表达降序取非零基因构成 cell sentence；
- 基线嵌入 ctrl_emb = 池内细胞嵌入均值（BERT last_hidden_state 对基因 token 取均值）；
- 扰动：从 cell sentence 删除扰动基因 token（Geneformer ISPy 风格的 in silico
  deletion），pred_emb[pi] = 池内均值；pred_delta_emb = pred_emb - ctrl_emb。
- 目标基因不在输入空间/不在字典：无法删除 → 该扰动 pred 退化为 ctrl（留痕
  target_in_input_space / target_deleted）。
- 原生输出为 embedding，不产生表达谱 → labels 的 pearson 列为 NaN（不硬造）。

断点续跑：每扰动写 <out>_checkpoint.npz。

用法（容器内）：
  python /research/scripts/run_geneformer_k562.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/results/geneformer_k562 [--pool 32] [--batch 16] [--threads 64]
产物：<out>_predictions.npz（同构键 + pred_emb/ctrl_emb/delta）+ _labels.csv + _meta.json
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
DEFAULT_OUT = "/research/results/geneformer_k562"
SEED = 0
MAX_LEN = 2048          # cell sentence 上限（非零基因数通常 < 2048）


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


def cell_sentences(ctrl_pool, var_genes, tok, n2i):
    """把 ctrl 池 (n_cells, n_genes) 转成 token id 序列（按表达降序，非零）。

    返回 list[list[int]]（去掉了无 token 的基因），以及
    gene_pos[i][gene] = 该基因在 cell i sentence 中的位置（-1=未表达/无 token）。
    """
    n_cells = ctrl_pool.shape[0]
    sentences = []
    pos_maps = []
    for i in range(n_cells):
        x = ctrl_pool[i]
        nz = np.where(x > 0)[0]
        order = nz[np.argsort(-x[nz])][:MAX_LEN - 2]
        toks, pmap = [], {}
        for gi in order:
            g = var_genes[gi]
            eid = n2i.get(g)
            t = tok.get(eid) if eid else None
            if t is None:
                continue
            pmap[g] = len(toks)
            toks.append(t)
        sentences.append(toks)
        pos_maps.append(pmap)
    return sentences, pos_maps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=DEFAULT_H5AD)
    ap.add_argument("--local-npz", default=DEFAULT_LOCAL_NPZ)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--pool", type=int, default=32)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SCT_THREADS", "64")))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    t_start = time.time()
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.set_num_threads(args.threads)
    device = torch.device("cpu")
    print(f"[cfg] pool={args.pool} batch={args.batch} threads={args.threads}",
          flush=True)

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

    CKPT = args.out + "_checkpoint.npz"
    resume_rows = []
    pool_idx_ckpt = None
    if args.limit == 0 and os.path.exists(CKPT):
        cz = np.load(CKPT, allow_pickle=True)
        if (list(cz["ckpt_gene_space"]) == var_genes
                and list(cz["ckpt_perts"]) == holdout_genes
                and int(cz["ckpt_pool"]) == args.pool):
            resume_rows = [int(i) for i in cz["ckpt_rows"]]
            pool_idx_ckpt = cz["ckpt_pool_idx"]
            print(f"[ckpt] 恢复 {len(resume_rows)} 个已完成扰动", flush=True)
        else:
            print("[ckpt] 指纹不匹配，从头开始", flush=True)

    rng = np.random.default_rng(SEED)
    pool_idx = pool_idx_ckpt if pool_idx_ckpt is not None else \
        rng.choice(np.where(is_ctrl)[0], size=min(args.pool, n_ctrl),
                   replace=False)
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)

    tok, n2i = load_dicts()
    print(f"[vocab] tokens={len(tok)}, name2id={len(n2i)}", flush=True)
    sentences, pos_maps = cell_sentences(ctrl_pool, var_genes, tok, n2i)
    lens = [len(s) for s in sentences]
    print(f"[tok] sentence len min/med/max = {min(lens)}/{int(np.median(lens))}/"
          f"{max(lens)}", flush=True)
    target_in_space = np.array([g in var_genes for g in holdout_genes])
    target_has_token = np.array(
        [tok.get(n2i.get(g)) is not None if n2i.get(g) else False
         for g in holdout_genes])
    print(f"[vocab] 目标不在输入空间: {int((~target_in_space).sum())}/217, "
          f"无 token: {int((~target_has_token).sum())}/217", flush=True)

    from transformers import BertModel
    model = BertModel.from_pretrained(WEIGHTS, output_attentions=False)
    model.to(device)
    model.eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"[model] Geneformer loaded, params={n_params/1e6:.0f}M", flush=True)

    def embed_batch(sents):
        """sents: list[list[int]] -> 每细胞基因 token 均值嵌入 (n, 768)。"""
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
        h = out.last_hidden_state              # (n, L, 768)
        mask = torch.from_numpy(att).to(device).unsqueeze(-1).float()
        # 基因 token 掩码：去掉 CLS(0) 与 EOS（每序列最后一个 att=1 位置）
        gene_mask = mask.clone()
        idx = att.sum(1) - 1
        gene_mask[torch.arange(len(sents)), torch.from_numpy(idx)] = 0
        gene_mask[:, 0] = 0
        emb = (h * gene_mask).sum(1) / gene_mask.sum(1).clamp(min=1)
        return emb.cpu().numpy().astype(np.float32)

    def embed_cells(drop_gene=None):
        sents = []
        for i, s in enumerate(sentences):
            if drop_gene is not None:
                p = pos_maps[i].get(drop_gene, None)
                if p is not None:
                    s = s[:p] + s[p + 1:]
            sents.append(s)
        outs = []
        for b in range(0, len(sents), args.batch):
            outs.append(embed_batch(sents[b:b + args.batch]))
        return np.concatenate(outs, axis=0)

    # ---- 基线 ctrl 嵌入 ----
    t0 = time.time()
    if os.path.exists(CKPT + ".ctrl_emb.npy") and args.limit == 0:
        ctrl_emb_cells = np.load(CKPT + ".ctrl_emb.npy")
        print(f"[emb] ctrl baseline loaded from ckpt", flush=True)
    else:
        ctrl_emb_cells = embed_cells()
        if args.limit == 0:
            np.save(CKPT + ".ctrl_emb.npy", ctrl_emb_cells)
    ctrl_emb = ctrl_emb_cells.mean(axis=0)
    print(f"[emb] ctrl baseline {ctrl_emb_cells.shape} in {time.time()-t0:.0f}s",
          flush=True)

    # ---- 逐扰动 deletion ----
    n_all = len(holdout_genes)
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    D = ctrl_emb.shape[0]
    pred_emb = np.full((n_all, D), np.nan, dtype=np.float32)
    target_deleted = np.zeros(n_all, dtype=bool)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        pred_emb[np.array(resume_rows)] = cz["ckpt_preds"]
        target_deleted[np.array(resume_rows)] = cz["ckpt_deleted"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for pi in run_rows:
        g = holdout_genes[pi]
        ti = time.time()
        deleted = bool(target_has_token[pi]) and any(
            g in pos_maps[i] for i in range(len(pos_maps)))
        if deleted:
            emb = embed_cells(drop_gene=g).mean(axis=0)
        else:
            emb = ctrl_emb.copy()
        pred_emb[pi] = emb
        target_deleted[pi] = deleted
        done_rows.append(pi)
        if not args.limit:
            np.savez(CKPT, ckpt_rows=np.array(done_rows),
                     ckpt_preds=pred_emb[np.array(done_rows)],
                     ckpt_deleted=target_deleted[np.array(done_rows)],
                     ckpt_perts=np.array(holdout_genes),
                     ckpt_gene_space=np.array(var_genes),
                     ckpt_pool=args.pool,
                     ckpt_pool_idx=np.asarray(pool_idx))
        n_done = len(done_rows)
        if n_done <= 5 or n_done % 20 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} deleted={deleted} "
                  f"{time.time()-ti:.1f}s (eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done in {t_infer/60:.1f} min", flush=True)

    # ---- 标签（pearson 留 NaN：embedding 级产物，不硬造表达谱）----
    labels = []
    cond_set = set(obs_cond)
    for pi, g in enumerate(holdout_genes):
        cands = [f"ctrl+{g}", f"{g}+ctrl"]
        cond = next((c for c in cands if c in cond_set), None)
        n_cells = int((obs_cond == cond).sum()) if cond else 0
        labels.append({"perturbation": g, "pearson_delta": np.nan,
                       "pearson_delta_all": np.nan, "n_cells": n_cells,
                       "target_in_input_space": bool(target_in_space[pi]),
                       "target_deleted": bool(target_deleted[pi]),
                       "delta_emb_norm": float(
                           np.linalg.norm(pred_emb[pi] - ctrl_emb))})
    adata.file.close()
    lab = pd.DataFrame(labels)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pred_delta = pred_emb - ctrl_emb[None, :]
    np.savez_compressed(
        args.out + "_predictions.npz",
        pred_perturbations=np.array(holdout_genes),
        pred_emb=pred_emb,                       # (217, 768) 扰动后细胞嵌入
        pred_delta_emb=pred_delta,               # (217, 768)
        ctrl_emb=ctrl_emb,                       # (768,)
        ctrl_emb_cells=ctrl_emb_cells,           # (32, 768) 每细胞基线嵌入
        ctrl_pool_indices=pool_idx.astype(np.int64),
        target_in_input_space=target_in_space,
        target_deleted=target_deleted,
        ctrl_mean_full=ctrl_mean_full,
        full_gene_space=np.array(var_genes),
        output_type=np.array(
            ["cell_embedding_via_in_silico_deletion (no expression prediction; "
             "pearson columns in _labels.csv are NaN by design)"]),
    )
    lab.to_csv(args.out + "_labels.csv", index=False)
    pd.DataFrame({"full_gene_space": var_genes}).to_csv(
        args.out + "_genes_full.csv", index=False)

    meta = {
        "task": "Geneformer V2-104M zero-shot in silico deletion: Replogle "
                "essential held-out perturbation (K562)",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {
            "name": "Geneformer V2-104M (ctheodoris/Geneformer)",
            "source": "hf-mirror ctheodoris/Geneformer Geneformer-V2-104M",
            "weights_dir": WEIGHTS,
            "weights_sha256": sha256(WEIGHTS + "model.safetensors"),
            "inference_protocol": "delete perturbed gene token from ranked cell "
                "sentence, cell embedding = mean of gene-token last_hidden_state, "
                "pred_emb = pool mean; delta vs ctrl baseline",
            "output_type": "cell embedding (768-d); NOT an expression profile",
        },
        "data": {
            "h5ad": args.h5ad,
            "heldout_source": args.local_npz,
            "split": "GEARS simulation, split_seed=0",
            "ctrl_pool_size": int(len(ctrl_pool)),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "torch_version": torch.__version__,
            "threads": args.threads,
            "infer_seconds_total": t_infer,
            "total_seconds": time.time() - t_start,
        },
    }
    with open(args.out + "_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    if not args.limit and os.path.exists(CKPT):
        os.remove(CKPT)
        for suf in (".ctrl_emb.npy",):
            if os.path.exists(CKPT + suf):
                os.remove(CKPT + suf)
    print(f"[out] -> {args.out}_predictions.npz (total "
          f"{time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
