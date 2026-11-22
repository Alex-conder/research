-rw-r--r--. 1 root root 16986 Oct 16 12:13 /opt/research/scripts/run_uce_k562.py
# -*- coding: utf-8 -*-
"""UCE (snap-stanford/UCE, 33-layer) 对 Replogle essential held-out 扰动的
零样本 in silico deletion 推理（K562）。输出类型：细胞嵌入（CLS, 1280-d，
如实标注，不硬造表达谱）。

协议（复刻 UCE eval_single_anndata 官方流程，只加确定性 seed 与扰动删除）：
- 基因子集到 UCE human protein-embedding 空间（spec_pe_genes ∩ spec_chrom csv）；
- 细胞句：按 log1p(权重) 有放回采 sample_size=1024 个基因 token，按染色体
  （乱序）→ 基因组起点排序，包裹染色体左右 token，CLS 开头，PAD 到 1536；
  权重 ∝ 表达值（h5ad X 为 log1p-normalized，log1p(expm1(x))=x，等价于
  log1p(counts)，协议偏离已记录）；
- 扰动：采样前将该基因权重置 0（=从细胞句删除该基因）；
- 嵌入 = decoder CLS token 输出，L2 归一化；pred_emb[pi] = 池(32)均值；
  ctrl_emb = 未删除基线；pred_delta_emb = 差。
- 模型：33L TransformerModel(token_dim=5120, d_model=1280, nhead=20,
  d_hid=5120, output_dim=1280)，权重/ESM2 token 嵌入/protein embeddings
  均为 figshare 官方文件；pe_embedding 加载同官方（all_tokens.torch +
  torch.manual_seed(23) 生成 1895 个染色体 token，若 checkpoint 已含
  145469 行则保留 checkpoint 内训练好的染色体嵌入）。

断点续跑：每扰动写 <out>_checkpoint.npz。

用法（容器内）：
  python /research/scripts/run_uce_k562.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/results/uce_k562 [--pool 32] [--batch 8] [--threads 64]
产物：<out>_predictions.npz + _labels.csv + _meta.json（pearson 列为 NaN，by design）
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
DEFAULT_OUT = "/research/results/uce_k562"

SEED = 0
PAD_LENGTH = 1536
SAMPLE_SIZE = 1024
PAD_TOKEN_IDX = 0
CHROM_TOKEN_LEFT_IDX = 1
CHROM_TOKEN_RIGHT_IDX = 2
CLS_TOKEN_IDX = 3
CHROM_TOKEN_OFFSET = 143574
OUT_DIM = 1280


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
    """UCE human 空间映射：返回 keep_idx, pe_row_idxs, chroms, starts。

    pe_row_idxs[i] = species_offset + human_pe_genes.index(gene.upper())，
    与官方 adata_path_to_prot_chrom_starts 完全一致。
    """
    from data_utils import get_spec_chrom_csv
    human_pe = torch.load(HUMAN_PE)
    human_pe = {k.upper(): v for k, v in human_pe.items()}
    spec_pe_genes = list(human_pe.keys())
    del human_pe  # 只需要基因名单（token 嵌入在 all_tokens.torch）
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
    print(f"[uce] mapped {len(keep)}/{len(var_genes)} genes to UCE human space, "
          f"offset={offset}", flush=True)
    return (np.array(keep, dtype=np.int64),
            torch.tensor(pe_idx, dtype=torch.long),
            np.array(chroms, dtype=np.int64),
            np.array(starts))


def sample_sentence(weights, pe_row_idxs, chroms, starts, rng):
    """单细胞句采样（复刻 eval_data.sample_cell_sentences，确定性 rng）。"""
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
    if all_pe.shape[0] == 143574:          # 官方逻辑：补染色体 token
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
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--pool", type=int, default=32)
    ap.add_argument("--batch", type=int, default=8)
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

    rng_pool = np.random.default_rng(SEED)
    pool_idx = pool_idx_ckpt if pool_idx_ckpt is not None else \
        rng_pool.choice(np.where(is_ctrl)[0], size=min(args.pool, n_ctrl),
                        replace=False)
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)

    keep, pe_row, chroms, starts = build_mapping(var_genes)
    mapped_genes = [var_genes[i] for i in keep]
    # 权重 ∝ 表达（log-normalized；见 docstring 等价性说明），负值截 0
    pool_w = np.clip(ctrl_pool[:, keep], 0, None)
    expressed = pool_w > 0
    print(f"[uce] pool weights ready {pool_w.shape}", flush=True)

    pos_in_mapped = {g: i for i, g in enumerate(mapped_genes)}
    target_in_space = np.array([g in pos_in_mapped for g in holdout_genes])
    print(f"[uce] 目标不在 UCE 空间的扰动: {int((~target_in_space).sum())}/217",
          flush=True)

    model = load_model(device)

    def embed_pool(drop_gene=None):
        """32 个 ctrl 细胞的 UCE 嵌入 (32, 1280)。drop_gene=从句中删除。"""
        drop_pos = pos_in_mapped.get(drop_gene, None)
        sents, masks = [], []
        for c in range(len(pool_w)):
            w = pool_w[c].copy()
            if drop_pos is not None:
                w[drop_pos] = 0.0
            if w.sum() <= 0:
                w = pool_w[c].copy()
            rng = np.random.default_rng(SEED * 1000003 + c)   # 每细胞确定性
            s, m, _ = sample_sentence(w, pe_row, chroms, starts, rng)
            sents.append(s)
            masks.append(m)
        outs = []
        for b in range(0, len(sents), args.batch):
            bs = np.stack(sents[b:b + args.batch])
            bm = np.stack(masks[b:b + args.batch])
            tok = model.pe_embedding(
                torch.from_numpy(bs).long().to(device))   # (B, L, 5120)
            src = tok.permute(1, 0, 2)                    # (L, B, 5120)
            with torch.no_grad():
                _, emb = model(src, torch.from_numpy(bm).to(device))
            outs.append(emb.cpu().numpy().astype(np.float32))
        return np.concatenate(outs, axis=0)

    t0 = time.time()
    ctrl_emb_cells = embed_pool()
    ctrl_emb = ctrl_emb_cells.mean(axis=0)
    ctrl_emb = ctrl_emb / max(np.linalg.norm(ctrl_emb), 1e-12)
    print(f"[emb] ctrl baseline {ctrl_emb_cells.shape} in {time.time()-t0:.0f}s",
          flush=True)

    n_all = len(holdout_genes)
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    pred_emb = np.full((n_all, OUT_DIM), np.nan, dtype=np.float32)
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
        deleted = bool(target_in_space[pi]) and \
            bool(expressed[:, pos_in_mapped[g]].any())
        emb = embed_pool(drop_gene=g).mean(axis=0) if deleted \
            else ctrl_emb_cells.mean(axis=0)
        emb = emb / max(np.linalg.norm(emb), 1e-12)
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
        if n_done <= 5 or n_done % 10 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} deleted={deleted} "
                  f"{time.time()-ti:.1f}s (eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done in {t_infer/60:.1f} min", flush=True)

    labels = []
    cond_set = set(obs_cond)
    ctrl_mapped_m = ctrl_mean_full[keep]
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
    np.savez_compressed(
        args.out + "_predictions.npz",
        pred_perturbations=np.array(holdout_genes),
        pred_emb=pred_emb,                       # (217, 1280) L2-normalized
        pred_delta_emb=pred_emb - ctrl_emb[None, :],
        ctrl_emb=ctrl_emb,
        ctrl_emb_cells=ctrl_emb_cells,
        ctrl_pool_indices=pool_idx.astype(np.int64),
        target_in_input_space=target_in_space,
        target_deleted=target_deleted,
        ctrl_mean_full=ctrl_mean_full,
        ctrl_mean_uce_space=ctrl_mapped_m,
        uce_gene_space=np.array(mapped_genes),
        pred_gene_indices_in_h5ad=keep.astype(np.int32),
        full_gene_space=np.array(var_genes),
        output_type=np.array(
            ["cell_embedding (UCE-33L CLS, 1280-d) via in_silico gene "
             "deletion; no expression prediction; pearson columns in "
             "_labels.csv are NaN by design"]),
    )
    lab.to_csv(args.out + "_labels.csv", index=False)
    pd.DataFrame({"uce_gene_space": mapped_genes}).to_csv(
        args.out + "_genes_uce.csv", index=False)
    pd.DataFrame({"full_gene_space": var_genes}).to_csv(
        args.out + "_genes_full.csv", index=False)

    meta = {
        "task": "UCE 33-layer zero-shot in silico deletion: Replogle "
                "essential held-out perturbation (K562)",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {
            "name": "UCE 33L (snap-stanford/UCE)",
            "source": "figshare model files (see RUNBOOK)",
            "weights_file": MODEL_33L,
            "weights_sha256": sha256(MODEL_33L),
            "token_file_sha256": sha256(TOKEN_FILE),
            "inference_protocol": "official eval_single_anndata pipeline with "
                "fixed per-cell seeds; deletion = zero target gene sampling "
                "weight; embedding = normalized CLS of decoder",
            "input_note": "h5ad X is log-normalized; weights proportional to "
                "X values == log1p(recovered counts) up to library-size "
                "normalization, sampling probabilities identical",
            "output_type": "cell embedding (1280-d); NOT an expression profile",
        },
        "data": {
            "h5ad": args.h5ad,
            "heldout_source": args.local_npz,
            "split": "GEARS simulation, split_seed=0",
            "ctrl_pool_size": int(len(pool_w)),
            "n_genes_uce_space": int(len(mapped_genes)),
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
    print(f"[out] -> {args.out}_predictions.npz (total "
          f"{time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

