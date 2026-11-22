# -*- coding: utf-8 -*-
"""UCE 33L wave 推理: 每细胞 embedding(pred 侧 + true 侧)。

口径(2026-10-18 多模型波):
- embedding 类误差 = 模型自身空间 d(pred_emb, true_emb);
- pred_emb = ctrl 池(pool=32, seed=0)删除扰动基因(采样权重置 0)后的
  每细胞 CLS 嵌入; true_emb = 真实扰动细胞(<=32/pert, seed=777000+pi
  子采样)的每细胞 CLS 嵌入;
- 与 run_uce_k562.py 的区别: 基线脚本只留 pred/ctrl 均值, 本脚本留每细胞
  双端嵌入, 供 d(pred_emb, true_emb) 与签名特征复算。

产物 <out>_predictions.npz: 与 geneformer wave 同构 + uce_gene_space /
pred_gene_indices_in_h5ad。断点: 每扰动写 <out>_checkpoint.npz。

用法(容器内):
  python -u /research/scripts/run_uce_wave.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/outputs_crossmodel/uce/uce_wave_k562
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
import torch.nn as nn

UCE_REPO = "/research/uce_repo/UCE-main"
sys.path.insert(0, UCE_REPO)
sys.path.insert(0, UCE_REPO + "/data_proc")

MODEL_FILES = "/research/models/uce/model_files/"
MODEL_33L = MODEL_FILES + "33l_8ep_1024t_1280.torch"
TOKEN_FILE = MODEL_FILES + "all_tokens.torch"
OFFSET_PKL = MODEL_FILES + "species_offsets.pkl"
SPEC_CHROM_CSV = MODEL_FILES + "species_chrom.csv"
HUMAN_PE = (MODEL_FILES + "protein_embeddings/" +
            "Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt")

DEFAULT_H5AD = "/research/data/perturb_processed_k562.h5ad"
DEFAULT_LOCAL_NPZ = "/research/meta/pred_expr_seed0.npz"

SEED = 0
PAD_LENGTH = 1536
SAMPLE_SIZE = 1024
PAD_TOKEN_IDX = 0
CHROM_TOKEN_LEFT_IDX = 1
CHROM_TOKEN_RIGHT_IDX = 2
CLS_TOKEN_IDX = 3
CHROM_TOKEN_OFFSET = 143574
OUT_DIM = 1280


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


def build_mapping(var_genes):
    from data_utils import get_spec_chrom_csv
    human_pe = torch.load(HUMAN_PE)
    human_pe = {k.upper(): v for k, v in human_pe.items()}
    spec_pe_genes = list(human_pe.keys())
    del human_pe
    with open(OFFSET_PKL, "rb") as f:
        offset = pickle.load(f)["human"]
    gene_to_chrom_pos = get_spec_chrom_csv(SPEC_CHROM_CSV)
    spec_chrom = gene_to_chrom_pos[
        gene_to_chrom_pos["species"] == "human"].set_index("gene_symbol")

    keep, pe_idx, chroms, starts = [], [], [], []
    for i, g in enumerate(var_genes):
        gu = g.upper()
        if gu not in spec_pe_genes or gu not in spec_chrom.index:
            continue
        keep.append(i)
        pe_idx.append(spec_pe_genes.index(gu) + offset)
        chroms.append(spec_chrom.loc[gu, "spec_chrom"])
        starts.append(spec_chrom.loc[gu, "start"])
    chroms = pd.Categorical(
        chroms, categories=gene_to_chrom_pos["spec_chrom"].cat.categories
    ).codes
    print(f"[uce] mapped {len(keep)}/{len(var_genes)} genes to UCE human space",
          flush=True)
    return (np.array(keep, dtype=np.int64),
            torch.tensor(pe_idx, dtype=torch.long),
            np.array(chroms, dtype=np.int64),
            np.array(starts))


def sample_sentence(weights, pe_row_idxs, chroms, starts, rng):
    w = weights / weights.sum()
    choice_idx = rng.choice(np.arange(len(w)), size=SAMPLE_SIZE, p=w,
                            replace=True)
    choosen_chrom = chroms[choice_idx]
    chrom_sort = np.argsort(choosen_chrom, kind="stable")
    choice_idx = choice_idx[chrom_sort]
    new_chrom = chroms[choice_idx]
    choosen_starts = starts[choice_idx]

    ordered = np.full(PAD_LENGTH, CLS_TOKEN_IDX, dtype=np.int64)
    i = 1
    uq_chroms = np.unique(new_chrom)
    rng.shuffle(uq_chroms)
    for chrom in uq_chroms:
        ordered[i] = int(chrom) + CHROM_TOKEN_OFFSET
        i += 1
        loc = np.where(new_chrom == chrom)[0]
        sort_by_start = np.argsort(choosen_starts[loc], kind="stable")
        to_add = choice_idx[loc[sort_by_start]]
        ordered[i:i + len(to_add)] = pe_row_idxs[to_add]
        i += len(to_add)
        ordered[i] = CHROM_TOKEN_RIGHT_IDX
        i += 1
    mask = np.zeros(PAD_LENGTH, dtype=np.float32)
    mask[:i] = 1.0
    ordered[i:] = PAD_TOKEN_IDX
    return ordered, mask, i


def load_model(device):
    from model import TransformerModel
    model = TransformerModel(token_dim=5120, d_model=1280, nhead=20,
                             d_hid=5120, nlayers=33, output_dim=OUT_DIM,
                             dropout=0.05)
    empty_pe = torch.zeros(145469, 5120)
    empty_pe.requires_grad = False
    model.pe_embedding = nn.Embedding.from_pretrained(empty_pe)
    sd = torch.load(MODEL_33L, map_location="cpu")
    model.load_state_dict(sd, strict=True)
    all_pe = torch.load(TOKEN_FILE)
    if all_pe.shape[0] == 143574:
        torch.manual_seed(23)
        chrom_tensors = torch.normal(mean=0, std=1, size=(1895, 5120))
        all_pe = torch.vstack((all_pe, chrom_tensors))
        all_pe.requires_grad = False
    if all_pe.shape[0] != 145469:
        all_pe.requires_grad = False
        model.pe_embedding = nn.Embedding.from_pretrained(all_pe)
    print(f"[model] UCE-33L loaded (pe={tuple(all_pe.shape)})", flush=True)
    model.to(device)
    model.eval()
    return model


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=DEFAULT_H5AD)
    ap.add_argument("--local-npz", default=DEFAULT_LOCAL_NPZ)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", type=int, default=32)
    ap.add_argument("--true-cells", type=int, default=32)
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SCT_THREADS", "48")))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    t_start = time.time()
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.set_num_threads(args.threads)
    device = torch.device("cpu")
    print(f"[cfg] pool={args.pool} true_cells={args.true_cells} "
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

    rng_pool = np.random.default_rng(SEED)
    pool_idx = pool_idx_ckpt if pool_idx_ckpt is not None else \
        rng_pool.choice(np.where(is_ctrl)[0], size=min(args.pool, n_ctrl),
                        replace=False)
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)

    keep, pe_row, chroms, starts = build_mapping(var_genes)
    mapped_genes = [var_genes[i] for i in keep]
    pos_in_mapped = {g: i for i, g in enumerate(mapped_genes)}
    target_in_space = np.array([g in pos_in_mapped for g in holdout_genes])
    print(f"[uce] 目标不在 UCE 空间的扰动: {int((~target_in_space).sum())}/"
          f"{len(holdout_genes)}", flush=True)

    model = load_model(device)

    def embed_cells(weights, drop_pos=None, seed_base=0):
        """weights (n, n_mapped) -> 每细胞 CLS 嵌入 (n, 1280)。"""
        sents, masks = [], []
        for c in range(len(weights)):
            w = weights[c].copy()
            if drop_pos is not None:
                w[drop_pos] = 0.0
            if w.sum() <= 0:
                w = weights[c].copy()
            rng = np.random.default_rng(seed_base + c)
            s, m, _ = sample_sentence(w, pe_row, chroms, starts, rng)
            sents.append(s)
            masks.append(m)
        outs = []
        for b in range(0, len(sents), args.batch):
            bs = np.stack(sents[b:b + args.batch])
            bm = np.stack(masks[b:b + args.batch])
            tok = model.pe_embedding(torch.from_numpy(bs).long().to(device))
            src = tok.permute(1, 0, 2)
            with torch.no_grad():
                _, emb = model(src, torch.from_numpy(bm).to(device))
            outs.append(emb.cpu().numpy().astype(np.float32))
        return np.concatenate(outs, axis=0)

    # ctrl 池权重(固定 seed_base 与基线脚本同池同采样)
    pool_w = np.clip(ctrl_pool[:, keep], 0, None)
    expressed = pool_w > 0

    t0 = time.time()
    if os.path.exists(CKPT + ".ctrl_emb.npy") and args.limit == 0:
        ctrl_emb_cells = np.load(CKPT + ".ctrl_emb.npy")
        print("[emb] ctrl baseline loaded from ckpt", flush=True)
    else:
        ctrl_emb_cells = embed_cells(pool_w, seed_base=SEED * 1000003)
        if args.limit == 0:
            np.save(CKPT + ".ctrl_emb.npy", ctrl_emb_cells)
    ctrl_emb = ctrl_emb_cells.mean(axis=0)
    ctrl_emb = ctrl_emb / max(np.linalg.norm(ctrl_emb), 1e-12)
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
    cap = args.true_cells
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    pred_cells = np.full((n_all, args.pool, OUT_DIM), np.nan, dtype=np.float32)
    true_cells = np.full((n_all, cap, OUT_DIM), np.nan, dtype=np.float32)
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
        dpos = pos_in_mapped.get(g, None)
        deleted = bool(target_in_space[pi]) and dpos is not None and \
            bool(expressed[:, dpos].any())
        if deleted:
            emb = embed_cells(pool_w, drop_pos=dpos,
                              seed_base=SEED * 1000003)
        else:
            emb = ctrl_emb_cells.copy()
        pred_cells[pi] = emb
        target_deleted[pi] = deleted
        idx = true_idx[pi]
        if len(idx):
            Xs = np.clip(np.asarray(
                adata[idx].X.todense()).astype(np.float32)[:, keep], 0, None)
            true_cells[pi, :len(idx)] = embed_cells(
                Xs, seed_base=SEED * 2000003 + pi * 10007)
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
        if n_done <= 5 or n_done % 10 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} deleted={deleted} "
                  f"{time.time()-ti:.1f}s (eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done {n_all} perts in {t_infer/60:.1f} min", flush=True)
    adata.file.close()

    pred_mean = row_mean(pred_cells)
    pred_mean = pred_mean / np.maximum(
        np.linalg.norm(pred_mean, axis=1, keepdims=True), 1e-12)
    true_mean = row_mean(true_cells)
    true_mean = true_mean / np.maximum(
        np.linalg.norm(true_mean, axis=1, keepdims=True), 1e-12)
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
        ctrl_mean_uce_space=ctrl_mean_full[keep],
        uce_gene_space=np.array(mapped_genes),
        pred_gene_indices_in_h5ad=keep.astype(np.int32),
        full_gene_space=np.array(var_genes),
        output_type=np.array(
            ["per_cell_embeddings (UCE-33L CLS, 1280-d; pred=ctrl pool w/ "
             "target deletion; true=real pert cells); embedding-class model"]),
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
        "task": "UCE-33L wave: per-cell embeddings K562/RPE1",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {"name": "UCE 33L",
                  "weights_sha256": sha256(MODEL_33L),
                  "inference_protocol": "official eval pipeline + fixed seeds; "
                      "deletion = zero sampling weight; per-cell CLS kept",
                  "output_type": "cell embedding (1280-d CLS)"},
        "data": {"h5ad": args.h5ad, "heldout_source": args.local_npz,
                 "ctrl_pool_size": int(len(pool_w)),
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
