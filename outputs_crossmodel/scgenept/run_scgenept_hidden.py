# -*- coding: utf-8 -*-
"""scGenePT 隐层表示导出（内部表示 vs 输出几何 探针实验用）。

由 run_scgenept_k562.py 复制改写：数据/基因匹配/ctrl 池/模型加载与其完全一致
（同 seed、同 pool 采样、同指纹），保证与 scgenept_*_predictions.npz 逐扰动对齐。

改动点：
- 对 model.transformer_encoder 最后一个 block 挂 forward hook，捕获
  (batch, seq_len, embsize) 隐层输出（该输出是 post-norm 的 block 输出、
  未经 encoder 末端 final LayerNorm——torch 2.1 nn.TransformerEncoder 的
  norm 在 layers 循环之后施加）。
- 序列内无显式 [CLS] token（cell_emb_style="cls" 实际取位置 0 的基因 token），
  因此每个扰动同时存三种 pooled 表示：
    hidden_meanpool  全序列 token 均值（每细胞 mean-pool 后再对 32 细胞池取均值）
    hidden_firsttok  位置 0 token（scGenePT cls-style cell embedding 的定义）
    hidden_targettok 扰动目标基因所在 token（目标不在输入空间的行为 NaN 行）
- 断点落盘 <out>_hidden_checkpoint.npz（指纹=基因空间+扰动名单+pool+batch），
  重启自动续跑；完成后写 <out>_hidden.npz + <out>_hidden_meta.json 并清理断点。

用法（容器内）：
  python /research/scripts/run_scgenept_hidden.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/results/scgenept_k562 [--pool 32] [--batch 8] [--threads 64]
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

REPO = "/research/scgenept_repo/repo"
sys.path.insert(0, REPO)

MODELS_DIR = "/research/models/"
MODEL_TYPE = "scgenept_ncbi+uniprot_gpt"
MODEL_WEIGHTS = (MODELS_DIR + "finetuned/scgenept_ncbi+uniprot/norman/"
                 "best_model_gpt3.5_ada_rnd_seed_42.pt")
VOCAB_FILE = REPO + "/models/pretrained/scgpt/vocab.json"
DEFAULT_H5AD = "/research/data/perturb_processed_k562.h5ad"
DEFAULT_LOCAL_NPZ = "/research/meta/pred_expr_seed0.npz"
DEFAULT_OUT = "/research/results/scgenept_k562"

SEED = 0


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=DEFAULT_H5AD)
    ap.add_argument("--local-npz", default=DEFAULT_LOCAL_NPZ,
                    help="held-out 名单 + 本地基因序（step6 口径）")
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--pool", type=int, default=32,
                    help="ctrl 细胞池大小（与 predictions.npz 同一池）")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--threads", type=int,
                    default=int(os.environ.get("SCT_THREADS", "64")))
    ap.add_argument("--limit", type=int, default=0, help="只跑前 N 个扰动（测速用，禁用断点）")
    args = ap.parse_args()

    t_start = time.time()
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    torch.set_num_threads(args.threads)
    device = torch.device("cpu")
    print(f"[cfg] h5ad={args.h5ad} pool={args.pool} batch={args.batch} "
          f"threads={args.threads} device=cpu (hidden export)", flush=True)

    # ---- 1. held-out 名单（与 run_scgenept_k562.py / 本地 step6 完全一致）----
    z = np.load(args.local_npz)
    holdout_genes = [str(g) for g in z["pred_perturbations"]]
    local_genes = [str(g) for g in z["pred_genes"]]
    print(f"[data] held-out perts={len(holdout_genes)}, "
          f"local pred gene space={len(local_genes)}", flush=True)

    # ---- 2. 读 h5ad（backed，只物化需要的行）----
    adata = ad.read_h5ad(args.h5ad, backed="r")
    var_genes = adata.var["gene_name"].astype(str).tolist()
    assert var_genes == local_genes, (
        "服务器 h5ad 基因序与本地 npz pred_genes 不一致！")
    obs_cond = adata.obs["condition"].astype(str).values
    is_ctrl = obs_cond == "ctrl"
    n_ctrl = int(is_ctrl.sum())
    print(f"[data] cells={adata.n_obs}, genes={adata.n_vars}, ctrl_cells={n_ctrl}",
          flush=True)

    # ---- 断点检测（与 run_scgenept_k562.py 同指纹约定）----
    CKPT = args.out + "_hidden_checkpoint.npz"
    resume_rows = []
    pool_idx_ckpt = None
    if args.limit == 0 and os.path.exists(CKPT):
        cz = np.load(CKPT, allow_pickle=True)
        fingerprint_ok = (
            list(cz["ckpt_gene_space"]) == var_genes
            and list(cz["ckpt_perts"]) == holdout_genes
            and int(cz["ckpt_pool"]) == args.pool
            and int(cz["ckpt_batch"]) == args.batch)
        if fingerprint_ok:
            resume_rows = [int(i) for i in cz["ckpt_rows"]]
            pool_idx_ckpt = cz["ckpt_pool_idx"]
            print(f"[ckpt] 恢复 {len(resume_rows)} 个已完成扰动，"
                  f"剩余 {len(holdout_genes) - len(resume_rows)}", flush=True)
        else:
            print("[ckpt] 指纹不匹配，忽略断点从头开始", flush=True)

    # ctrl 池（确定性采样；与 predictions 同一 seed/池）
    rng = np.random.default_rng(SEED)
    if pool_idx_ckpt is not None:
        pool_idx = pool_idx_ckpt
    else:
        pool_idx = rng.choice(np.where(is_ctrl)[0],
                              size=min(args.pool, n_ctrl), replace=False)
    t0 = time.time()
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)
    print(f"[data] ctrl pool {ctrl_pool.shape} densified in {time.time()-t0:.0f}s",
          flush=True)

    # ---- 3. 基因匹配 scGPT vocab（与 run_scgenept_k562.py 相同）----
    from utils.data_loading import (get_embs_to_include,
                                    initialize_genept_embeddings,
                                    initialize_go_embeddings,
                                    match_genes_to_scgpt_vocab_from_adata)
    from utils.scgpt_config import (D_HID, EMBSIZE, NHEAD, N_LAYERS_CLS, N_CLS,
                                    NLAYERS, PAD_TOKEN, PAD_VALUE, PERT_PAD_ID,
                                    SPECIAL_TOKENS)
    from models.scGenePT import scGenePT

    stub = ad.AnnData(X=np.zeros((1, len(var_genes)), dtype=np.float32),
                      var=pd.DataFrame({"gene_name": var_genes}))
    vocab, all_gene_ids, dataset_genes, gene2idx = \
        match_genes_to_scgpt_vocab_from_adata(VOCAB_FILE, stub, SPECIAL_TOKENS)
    in_vocab = np.array([g in vocab for g in dataset_genes])
    keep = np.where(in_vocab)[0]
    n_matched = len(keep)
    print(f"[vocab] matched {n_matched}/{len(var_genes)} genes "
          f"(vocab size={len(vocab)})", flush=True)
    matched_genes = [var_genes[i] for i in keep]
    gene_ids = all_gene_ids[keep]
    ctrl_pool = ctrl_pool[:, keep]          # 与 gene_ids 对齐到匹配子空间
    print(f"[data] ctrl pool sliced to matched space {ctrl_pool.shape}",
          flush=True)

    target_in_space = np.array(
        [g in gene2idx and bool(in_vocab[gene2idx[g]]) for g in holdout_genes])
    print(f"[vocab] held-out 目标基因不在输入空间的扰动: "
          f"{int((~target_in_space).sum())}/{len(holdout_genes)}", flush=True)

    # ---- 4. 构建模型并加载权重（与 run_scgenept_k562.py 相同）----
    embs_to_include = get_embs_to_include(MODEL_TYPE)
    genept_embs, genept_emb_type, genept_emb_dim, _ = \
        initialize_genept_embeddings(embs_to_include, dataset_genes, vocab,
                                     MODEL_TYPE, MODELS_DIR)
    go_embs_to_include, go_emb_type, go_emb_dim, _ = \
        initialize_go_embeddings(embs_to_include, dataset_genes, vocab,
                                 MODEL_TYPE, MODELS_DIR)
    t0 = time.time()
    model = scGenePT(
        ntoken=len(vocab), d_model=EMBSIZE, nhead=NHEAD, d_hid=D_HID,
        nlayers=NLAYERS, nlayers_cls=N_LAYERS_CLS, n_cls=N_CLS, vocab=vocab,
        n_perturbagens=2, dropout=0.0, pad_token=PAD_TOKEN, pad_value=PAD_VALUE,
        pert_pad_id=PERT_PAD_ID, use_fast_transformer=False,
        embs_to_include=embs_to_include, genept_embs=genept_embs,
        genept_emb_type=genept_emb_type, genept_emb_size=genept_emb_dim,
        go_embs_to_include=go_embs_to_include, go_emb_type=go_emb_type,
        go_emb_size=go_emb_dim)
    params = torch.load(MODEL_WEIGHTS, weights_only=True, map_location=device)
    params = {k.replace("Wqkv.", "in_proj_"): v for k, v in params.items()}
    model.load_state_dict(params)
    model.to(device)
    model.eval()
    print(f"[model] loaded in {time.time()-t0:.0f}s "
          f"(weights={os.path.getsize(MODEL_WEIGHTS)/1e6:.0f}MB)", flush=True)

    # ---- 4b. 挂隐层 hook（最后一个 transformer block 的输出）----
    # torch 2.1 的 nn.TransformerEncoder 在 eval + 无 padding 时走 sparsity
    # fast path：逐层调用仍在，但层输出是 NestedTensor——hook 里转成 padded
    # 稠密张量（本任务 pad_mask 全 False，padded 结果即完整输出）。
    if hasattr(model.transformer_encoder, "enable_nested_tensor"):
        model.transformer_encoder.enable_nested_tensor = False
    last_layer = model.transformer_encoder.layers[-1]
    captured = []

    def _hook(module, inputs, output):
        out = output[0] if isinstance(output, tuple) else output
        if out.is_nested:
            out = torch.nested.to_padded_tensor(out, 0.0)
        captured.append(out.detach().float().cpu())

    last_layer.register_forward_hook(_hook)
    print(f"[hook] registered on transformer_encoder.layers[-1] "
          f"(nlayers={len(model.transformer_encoder.layers)}, embsize={EMBSIZE})",
          flush=True)

    # ---- 5. 逐扰动批式推理（隐层导出 + 断点落盘）----
    gene_ids_t = torch.tensor(gene_ids).long().unsqueeze(0).to(device)
    pad_mask = torch.zeros_like(gene_ids_t, dtype=torch.bool, device=device)
    pos_in_matched = {g: i for i, g in enumerate(matched_genes)}

    n_all = len(holdout_genes)
    d_hid = int(EMBSIZE)
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []          # 测速模式不用断点
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]

    hid_mean = np.full((n_all, d_hid), np.nan, dtype=np.float32)
    hid_first = np.full((n_all, d_hid), np.nan, dtype=np.float32)
    hid_target = np.full((n_all, d_hid), np.nan, dtype=np.float32)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        hid_mean[np.array(resume_rows)] = cz["ckpt_hidden_meanpool"]
        hid_first[np.array(resume_rows)] = cz["ckpt_hidden_firsttok"]
        hid_target[np.array(resume_rows)] = cz["ckpt_hidden_targettok"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for k, pi in enumerate(run_rows):
        g = holdout_genes[pi]
        ti = time.time()
        flags = np.zeros((len(ctrl_pool), n_matched), dtype=np.float32)
        tgt_pos = pos_in_matched.get(g, None)
        if tgt_pos is not None:
            flags[:, tgt_pos] = 1.0
        flags_t = torch.from_numpy(flags).long()
        sum_mean = np.zeros(d_hid, dtype=np.float64)
        sum_first = np.zeros(d_hid, dtype=np.float64)
        sum_target = np.zeros(d_hid, dtype=np.float64)
        n_cells = 0
        with torch.no_grad():
            for s in range(0, len(ctrl_pool), args.batch):
                xb = torch.from_numpy(ctrl_pool[s:s + args.batch]).float()
                gi = gene_ids_t.expand(len(xb), -1)
                pm = pad_mask.expand(len(xb), -1)
                del captured[:]
                model(gi.to(device), xb.to(device), flags_t[s:s+args.batch].to(device),
                      src_key_padding_mask=pm.to(device),
                      CLS=False, CCE=False, MVC=False, ECS=False,
                      do_sample=False)
                if not captured:
                    raise RuntimeError("隐层 hook 未触发——检查 "
                                       "enable_nested_tensor fast path")
                h = captured[0]                       # (batch, seq, d_hid)
                if h.shape[:2] != (len(xb), n_matched):
                    raise RuntimeError(f"hook 捕获 shape {tuple(h.shape)} "
                                       f"!= ({len(xb)}, {n_matched}, {d_hid})")
                b = len(xb)
                sum_mean += h.mean(dim=1).numpy().sum(axis=0)
                sum_first += h[:, 0, :].numpy().sum(axis=0)
                if tgt_pos is not None:
                    sum_target += h[:, tgt_pos, :].numpy().sum(axis=0)
                n_cells += b
        hid_mean[pi] = (sum_mean / n_cells).astype(np.float32)
        hid_first[pi] = (sum_first / n_cells).astype(np.float32)
        if tgt_pos is not None:
            hid_target[pi] = (sum_target / n_cells).astype(np.float32)
        done_rows.append(pi)
        if not args.limit:                      # 断点落盘（~1.7MB，秒级）
            idx = np.array(done_rows)
            np.savez(CKPT, ckpt_rows=idx,
                     ckpt_hidden_meanpool=hid_mean[idx],
                     ckpt_hidden_firsttok=hid_first[idx],
                     ckpt_hidden_targettok=hid_target[idx],
                     ckpt_perts=np.array(holdout_genes),
                     ckpt_gene_space=np.array(var_genes),
                     ckpt_pool=args.pool, ckpt_batch=args.batch,
                     ckpt_pool_idx=np.asarray(pool_idx))
        n_done = len(done_rows)
        if n_done <= 5 or n_done % 10 == 0:
            el = time.time() - t_inf0
            eta = el / max(1, n_done - len(resume_rows)) * (n_all - n_done)
            print(f"[inf] {n_done}/{n_all} {g} {time.time()-ti:.1f}s "
                  f"(eta {eta/60:.1f}min)", flush=True)
    t_infer = time.time() - t_inf0
    print(f"[inf] done {n_all} perts x {len(ctrl_pool)} cells "
          f"in {t_infer/60:.1f} min", flush=True)
    adata.file.close()

    # ---- 6. 产物 ----
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez(
        args.out + "_hidden.npz",
        hidden_meanpool=hid_mean,
        hidden_firsttok=hid_first,
        hidden_targettok=hid_target,
        perturbations=np.array(holdout_genes),
        target_in_input_space=target_in_space,
        hidden_layer=np.array("transformer_encoder.layers[-1] (post-norm, "
                              "pre-final-encoder-LN)"),
        embsize=np.int64(d_hid),
        n_layers=np.int64(NLAYERS),
        pool=np.int64(len(ctrl_pool)),
        matched_gene_space=np.array(matched_genes),
    )
    print(f"[out] {args.out}_hidden.npz keys: hidden_meanpool/hidden_firsttok/"
          f"hidden_targettok {hid_mean.shape}, perturbations ({n_all},)", flush=True)

    # ---- 7. 元信息 ----
    with open(os.path.join(REPO, "GIT_COMMIT")) as f:
        git_commit = f.read().strip()
    import subprocess
    meta = {
        "task": "scGenePT hidden-state export (last transformer block output) "
                "for internal-representation vs output-geometry probe",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {
            "name": "scGenePT (scgenept_ncbi+uniprot, Norman finetuned)",
            "model_type": MODEL_TYPE,
            "repo": "https://github.com/czi-ai/scGenePT",
            "git_commit": git_commit,
            "weights_file": MODEL_WEIGHTS,
            "weights_sha256": sha256(MODEL_WEIGHTS),
        },
        "data": {
            "h5ad": args.h5ad,
            "heldout_source": args.local_npz,
            "n_genes_matched_vocab": int(n_matched),
            "ctrl_pool_size": int(len(ctrl_pool)),
        },
        "hidden_spec": {
            "source": "forward hook on model.transformer_encoder.layers[-1]; "
                      "captured tensor (batch, seq_len, embsize) is the last "
                      "block's post-norm output, BEFORE the encoder's final "
                      "LayerNorm (torch 2.1 applies norm after the layer loop)",
            "note_no_cls_token": "序列内无显式 [CLS] token；cell_emb_style="
                                 "'cls' 实际取位置 0 的基因 token，故同时导出 "
                                 "hidden_meanpool / hidden_firsttok / "
                                 "hidden_targettok（目标基因 token，不在输入"
                                 "空间的行为全 NaN）",
            "pooling": "per-cell pooling over seq dim, then mean over the "
                       f"{len(ctrl_pool)} ctrl cells (same pool as "
                       "scgenept_*_predictions.npz)",
            "embsize": d_hid,
            "n_layers": int(NLAYERS),
        },
        "runtime": {
            "command": " ".join(sys.argv),
            "torch_version": torch.__version__,
            "threads": args.threads,
            "n_perts": int(n_all),
            "resumed_from_checkpoint": int(len(resume_rows)),
            "seconds_per_pert_mean": t_infer / max(1, len(run_rows)),
            "infer_seconds_total": t_infer,
            "total_seconds": time.time() - t_start,
        },
        "host": subprocess.check_output(["uname", "-a"]).decode().strip(),
    }
    with open(args.out + "_hidden_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    if not args.limit and os.path.exists(CKPT):
        os.remove(CKPT)          # 完成后清理断点
    print(f"[out] -> {args.out}_hidden.npz / _hidden_meta.json "
          f"(total {time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
