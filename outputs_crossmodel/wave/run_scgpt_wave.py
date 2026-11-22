# -*- coding: utf-8 -*-
"""scGPT whole-human wave 推理: 每细胞 embedding(pred 侧 + true 侧)。

口径(2026-10-18 多模型波):
- embedding 类误差 = 模型自身空间 d(pred_emb, true_emb);
- pred_emb = ctrl 池(pool=32, seed=0)在 mask 协议(扰动基因值置 -1)下的
  每细胞 cell_emb(model.forward 自带 _get_cell_emb_from_layer);
- true_emb = 真实扰动细胞(<=32/pert, seed=777000+pi 子采样)的每细胞 cell_emb;
- 与 run_scgpt_k562.py 的差别: 该脚本输出 mlm_output 表达谱(mask 协议),
  本脚本输出 cell_emb(embedding 类口径), 两者互补不冲突。

产物 <out>_predictions.npz: 同 geneformer wave 同构键 + scgpt 特有
  cell_emb_style / vocab_size / pred_gene_indices_in_h5ad。
断点: 每扰动写 <out>_checkpoint.npz。

用法(容器内):
  python -u /research/scripts/run_scgpt_wave.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/outputs_crossmodel/scgpt/scgpt_wave_k562
"""
import argparse
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import torch

MODEL_DIR = "/research/models/scgpt_whole_human/"
MODEL_WEIGHTS = MODEL_DIR + "best_model.pt"
VOCAB_FILE = MODEL_DIR + "vocab.json"
DEFAULT_H5AD = "/research/data/perturb_processed_k562.h5ad"
DEFAULT_LOCAL_NPZ = "/research/meta/pred_expr_seed0.npz"

SEED = 0
MASK_VALUE = -1.0
PAD_VALUE = -2.0
CFG = dict(embsize=512, nhead=8, d_hid=512, nlayers=12,
           nlayers_cls=3, n_cls=1, mvc_decoder_style="inner product")


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


def load_model(device):
    from scgpt.model.model import TransformerModel
    from scgpt.tokenizer import GeneVocab

    vocab = GeneVocab.from_file(VOCAB_FILE)
    if "<pad>" not in vocab:
        vocab.append_token("<pad>")
    model = TransformerModel(
        ntoken=len(vocab), d_model=CFG["embsize"], nhead=CFG["nhead"],
        d_hid=CFG["d_hid"], nlayers=CFG["nlayers"],
        nlayers_cls=CFG["nlayers_cls"], n_cls=CFG["n_cls"], vocab=vocab,
        dropout=0.0, pad_token="<pad>", pad_value=PAD_VALUE,
        do_mvc=True, mvc_decoder_style=CFG["mvc_decoder_style"],
        use_fast_transformer=False)
    ckpt = torch.load(MODEL_WEIGHTS, weights_only=True, map_location=device)
    sd = ckpt.get("model_state_dict", ckpt) if isinstance(ckpt, dict) else ckpt
    sd = {k.replace(".self_attn.Wqkv.", ".self_attn.in_proj_"): v
          for k, v in sd.items()}
    missing, unexpected = model.load_state_dict(sd, strict=False)
    crit_miss = [k for k in missing
                 if not k.startswith(("cls_decoder", "batch_encoder",
                                      "mvc_decoder", "dab_decoder",
                                      "ecs_decoder"))]
    if unexpected or crit_miss:
        print(f"[model] WARN unexpected={len(unexpected)} "
              f"crit_missing={crit_miss}", flush=True)
    model.to(device)
    model.eval()
    return model, vocab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=DEFAULT_H5AD)
    ap.add_argument("--local-npz", default=DEFAULT_LOCAL_NPZ)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", type=int, default=32)
    ap.add_argument("--true-cells", type=int, default=32)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SCT_THREADS", "32")))
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
    assert var_genes == local_genes, "h5ad 基因序与 local npz 不一致"
    obs_cond = adata.obs["condition"].astype(str).values
    is_ctrl = obs_cond == "ctrl"
    n_ctrl = int(is_ctrl.sum())
    print(f"[data] cells={adata.n_obs}, ctrl_cells={n_ctrl}", flush=True)
    ctrl_mean_full = np.asarray(
        adata[is_ctrl].X.mean(axis=0)).ravel().astype(np.float32)

    CKPT = args.out + "_checkpoint.npz"
    resume_rows, pool_idx_ckpt = [], None
    if args.limit == 0 and os.path.exists(CKPT):
        cz = np.load(CKPT, allow_pickle=True)
        if (list(cz["ckpt_gene_space"]) == var_genes
                and list(cz["ckpt_perts"]) == holdout_genes
                and int(cz["ckpt_pool"]) == args.pool
                and int(cz["ckpt_true"]) == args.true_cells
                and int(cz["ckpt_batch"]) == args.batch):
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

    model, vocab = load_model(device)
    in_vocab = np.array([g in vocab for g in var_genes])
    keep = np.where(in_vocab)[0]
    matched_genes = [var_genes[i] for i in keep]
    gene_ids = np.array([vocab[g] for g in matched_genes], dtype=np.int64)
    ctrl_pool = ctrl_pool[:, keep]
    pos_in_matched = {g: i for i, g in enumerate(matched_genes)}
    target_in_space = np.array([g in pos_in_matched for g in holdout_genes])
    print(f"[vocab] matched {len(keep)}/{len(var_genes)} "
          f"(cell_emb_style={model.cell_emb_style})", flush=True)

    gene_ids_t = torch.tensor(gene_ids).long().unsqueeze(0).to(device)
    pad_mask = torch.zeros_like(gene_ids_t, dtype=torch.bool, device=device)

    def embed_cells(x, mask_target=None):
        """x (n, n_matched) -> 每细胞 cell_emb (n, d)。mask_target=列号或 None。"""
        xx = x.copy()
        if mask_target is not None:
            xx[:, mask_target] = MASK_VALUE
        xt = torch.from_numpy(xx).float()
        outs = []
        with torch.no_grad():
            for s in range(0, len(xt), args.batch):
                xb = xt[s:s + args.batch].to(device)
                gi = gene_ids_t.expand(len(xb), -1)
                # 快速路径: 只要 cell_emb, 跳过 mlm_decoder(512->60697 投影,
                # 每细胞 ~300 GFLOP)。与完整 forward 的 cell_emb 逐元素等价
                # (已验证 max abs diff = 0.0)。
                # src_key_padding_mask 传 None: 本协议无 padding(全 4813 基因
                # 在场), 全 False 掩码只会触发 CPU 上更慢的 nested-tensor
                # 路径(基准 118s vs 74.7s/batch8, 等价性 0.0)。
                enc = model._encode(gi.to(device), xb, None)
                outs.append(model._get_cell_emb_from_layer(enc, xb)
                            .float().cpu().numpy())
        return np.concatenate(outs, axis=0)

    # ---- 基线 ctrl 嵌入 ----
    t0 = time.time()
    if os.path.exists(CKPT + ".ctrl_emb.npy") and args.limit == 0:
        ctrl_emb_cells = np.load(CKPT + ".ctrl_emb.npy")
        print("[emb] ctrl baseline loaded from ckpt", flush=True)
    else:
        ctrl_emb_cells = embed_cells(ctrl_pool)
        if args.limit == 0:
            np.save(CKPT + ".ctrl_emb.npy", ctrl_emb_cells)
    ctrl_emb = ctrl_emb_cells.mean(axis=0)
    print(f"[emb] ctrl baseline {ctrl_emb_cells.shape} in {time.time()-t0:.0f}s",
          flush=True)

    # ---- true 细胞选择 ----
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
    target_masked = np.zeros(n_all, dtype=bool)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        rr = np.array(resume_rows)
        pred_cells[rr] = cz["ckpt_pred_cells"]
        true_cells[rr] = cz["ckpt_true_cells"]
        true_n[rr] = cz["ckpt_true_n"]
        target_masked[rr] = cz["ckpt_masked"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for pi in run_rows:
        g = holdout_genes[pi]
        ti = time.time()
        ppos = pos_in_matched.get(g, None)
        target_masked[pi] = ppos is not None
        pred_cells[pi] = embed_cells(ctrl_pool, mask_target=ppos)
        idx = true_idx[pi]
        if len(idx):
            Xs = np.asarray(adata[idx].X.todense()).astype(np.float32)[:, keep]
            true_cells[pi, :len(idx)] = embed_cells(Xs)
            true_n[pi] = len(idx)
        done_rows.append(pi)
        if not args.limit:
            rr = np.array(done_rows)
            np.savez(CKPT, ckpt_rows=rr,
                     ckpt_pred_cells=pred_cells[rr],
                     ckpt_true_cells=true_cells[rr],
                     ckpt_true_n=true_n[rr],
                     ckpt_masked=target_masked[rr],
                     ckpt_perts=np.array(holdout_genes),
                     ckpt_gene_space=np.array(var_genes),
                     ckpt_pool=args.pool, ckpt_true=args.true_cells,
                     ckpt_batch=args.batch,
                     ckpt_pool_idx=np.asarray(pool_idx))
        n_done = len(done_rows)
        if n_done <= 5 or n_done % 10 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} masked={ppos is not None} "
                  f"{time.time()-ti:.1f}s (eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done {n_all} perts in {t_infer/60:.1f} min", flush=True)
    adata.file.close()

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
        target_masked=target_masked,
        pred_gene_indices_in_h5ad=keep.astype(np.int32),
        ctrl_mean_full=ctrl_mean_full,
        full_gene_space=np.array(var_genes),
        cell_emb_style=np.array([str(model.cell_emb_style)]),
        output_type=np.array(
            ["per_cell_embeddings (pred=ctrl pool masked-target cell_emb; "
             "true=real pert cells cell_emb); embedding-class model"]),
    )
    labels = []
    for pi, g in enumerate(holdout_genes):
        pm, tm = pred_mean[pi], true_mean[pi]
        labels.append({
            "perturbation": g, "pearson_delta": np.nan,
            "pearson_delta_all": np.nan, "n_cells": int(true_n[pi]),
            "target_in_input_space": bool(target_in_space[pi]),
            "target_masked": bool(target_masked[pi]),
            "pred_emb_norm": float(np.linalg.norm(pm)),
            "true_emb_norm": float(np.linalg.norm(tm)),
            "delta_emb_norm": float(np.linalg.norm(pm - ctrl_emb)),
            "err_emb_dist": float(np.linalg.norm(pm - tm)),
        })
    pd.DataFrame(labels).to_csv(args.out + "_labels.csv", index=False)

    meta = {
        "task": "scGPT whole-human wave: per-cell cell_emb K562/RPE1",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {"name": "scGPT whole-human",
                  "weights_sha256": sha256(MODEL_WEIGHTS),
                  "cell_emb_style": str(model.cell_emb_style),
                  "inference_protocol": "mask target to -1 on ctrl pool; "
                      "true cells unmasked; cell_emb via _encode + "
                      "_get_cell_emb_from_layer (skips 60697-vocab MLM head; "
                      "verified elementwise equal to full-forward cell_emb); "
                      "per-cell kept",
                  "output_type": "cell embedding (512-d cell_emb)"},
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
