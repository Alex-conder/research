# -*- coding: utf-8 -*-
"""scGenePT (scgenept_ncbi+uniprot, Norman-finetuned) 对 Replogle essential
held-out 扰动的零样本跨数据集推理（K562 与 RPE1 通用）。

口径与本地 step1_run.py / step6_run.py 对齐：
- held-out 扰动名单 = 本地 outputs*/pred_expr_seed0.npz 的 pred_perturbations
  （split_seed=0 的 simulation test 集，与 predictions_seed0.csv 同序）。
- ctrl 均值 / 真值均值 = 同一 perturb_processed.h5ad 按 obs.condition 分组求均值。
- pearson_delta（主列）= top-50 |true_delta| 基因上的 Pearson(pred_delta, true_delta)；
  pearson_delta_all = 匹配基因全集上的 Pearson。

跨数据集适配说明（协议偏离，记录于 meta）：
- 模型在 Norman(5045 基因)上微调，此处输入为目标数据集的 ~5000 基因表达向量；
  基因经 scGPT whole-human vocab 匹配后输入。
- 推理为确定性均值（explicit_zero_prob=False，mlm_output=pred，与官方训练配置一致）。

断点续跑：每完成 1 个扰动即写 <out>_checkpoint.npz（行号+预测+ctrl 池+基因空间指纹），
重启时校验指纹（基因空间/池大小/batch/扰动名单）一致则跳过已完成行并复用同一 ctrl 池。

用法（容器内）：
  python /research/scripts/run_scgenept_k562.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/results/scgenept_k562 [--pool 32] [--batch 8] [--threads 64]
产物：<out>_predictions.npz + _labels.csv + _meta.json（+ 运行中 _checkpoint.npz）
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

TOP_N_DEG = 50
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
                    help="ctrl 细胞池大小（每扰动从池内均值）")
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
          f"threads={args.threads} device=cpu", flush=True)

    # ---- 1. held-out 名单（与本地 step6 完全一致）----
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

    # ctrl 均值（全基因，供下游 step6 特征用）
    ctrl_mean_full = np.asarray(
        adata[is_ctrl].X.mean(axis=0)).ravel().astype(np.float32)

    # ---- 断点检测（在采样 ctrl 池之前，保证续跑用同一池）----
    CKPT = args.out + "_checkpoint.npz"
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

    # ctrl 池（确定性采样；续跑时复用原池）
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

    # ---- 3. 基因匹配 scGPT vocab ----
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

    # 目标基因不在输入基因空间（h5ad var）的扰动：无法打 pert flag（数据事实：
    # GEARS 图含扰动节点，token 模型不含），模型实际看到未扰动输入；仍照常评估并留痕。
    target_in_space = np.array(
        [g in gene2idx and bool(in_vocab[gene2idx[g]]) for g in holdout_genes])
    n_no_target = int((~target_in_space).sum())
    print(f"[vocab] held-out 目标基因不在输入空间的扰动: {n_no_target}/"
          f"{len(holdout_genes)}（模型收不到扰动信号，照常评估并留痕）",
          flush=True)

    # ---- 4. 构建模型并加载权重 ----
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

    # ---- 5. 逐扰动批式推理（带断点落盘）----
    gene_ids_t = torch.tensor(gene_ids).long().unsqueeze(0).to(device)
    pad_mask = torch.zeros_like(gene_ids_t, dtype=torch.bool, device=device)
    pos_in_matched = {g: i for i, g in enumerate(matched_genes)}

    n_all = len(holdout_genes)
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []          # 测速模式不用断点
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    preds = np.full((n_all, n_matched), np.nan, dtype=np.float32)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        preds[np.array(resume_rows)] = cz["ckpt_preds"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for k, pi in enumerate(run_rows):
        g = holdout_genes[pi]
        ti = time.time()
        flags = np.zeros((len(ctrl_pool), n_matched), dtype=np.float32)
        if g in pos_in_matched:
            flags[:, pos_in_matched[g]] = 1.0
        flags_t = torch.from_numpy(flags).long()
        outs = []
        with torch.no_grad():
            for s in range(0, len(ctrl_pool), args.batch):
                xb = torch.from_numpy(ctrl_pool[s:s + args.batch]).float()
                gi = gene_ids_t.expand(len(xb), -1)
                pm = pad_mask.expand(len(xb), -1)
                od = model(gi.to(device), xb.to(device), flags_t[s:s+args.batch].to(device),
                           src_key_padding_mask=pm.to(device),
                           CLS=False, CCE=False, MVC=False, ECS=False,
                           do_sample=False)
                outs.append(od["mlm_output"].float().cpu().numpy())
        preds[pi] = np.concatenate(outs, axis=0).mean(axis=0)
        done_rows.append(pi)
        if not args.limit:                      # 断点落盘（uncompressed，~4MB，秒级）
            np.savez(CKPT, ckpt_rows=np.array(done_rows),
                     ckpt_preds=preds[np.array(done_rows)],
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

    # ---- 6. 真值与标签（口径同 step1_run.py）----
    ctrl_m = ctrl_mean_full[keep]
    labels = []
    true_expr = np.full((n_all, n_matched), np.nan, dtype=np.float32)
    cond_set = set(obs_cond)
    for pi, g in enumerate(holdout_genes):
        cands = [f"ctrl+{g}", f"{g}+ctrl"]
        cond = next((c for c in cands if c in cond_set), None)
        if cond is None:
            labels.append({"perturbation": g, "pearson_delta": np.nan,
                           "pearson_delta_all": np.nan, "n_cells": 0,
                           "target_in_input_space": bool(target_in_space[pi])})
            continue
        sub = adata[obs_cond == cond]
        tm = np.asarray(sub.X.mean(axis=0)).ravel().astype(np.float32)[keep]
        true_expr[pi] = tm
        p_delta = preds[pi] - ctrl_m
        t_delta = tm - ctrl_m
        r_all = float(np.corrcoef(p_delta, t_delta)[0, 1])
        top = np.argsort(-np.abs(t_delta))[:TOP_N_DEG]
        r_top = float(np.corrcoef(p_delta[top], t_delta[top])[0, 1])
        labels.append({"perturbation": g, "pearson_delta": r_top,
                       "pearson_delta_all": r_all, "n_cells": int(sub.n_obs),
                       "target_in_input_space": bool(target_in_space[pi])})
    adata.file.close()

    lab = pd.DataFrame(labels)
    valid = lab.dropna(subset=["pearson_delta"])
    v_flag = valid[valid.target_in_input_space]
    print(f"[eval] 有效={len(valid)}/{len(lab)} "
          f"pearson_delta(top{TOP_N_DEG}): mean={valid.pearson_delta.mean():.4f} "
          f"median={valid.pearson_delta.median():.4f} | all: "
          f"mean={valid.pearson_delta_all.mean():.4f}", flush=True)
    if len(v_flag):
        print(f"[eval] 目标在输入空间的子集 n={len(v_flag)}: "
              f"pearson_delta mean={v_flag.pearson_delta.mean():.4f} "
              f"median={v_flag.pearson_delta.median():.4f}", flush=True)

    # ---- 7. 产物（含全基因空间 NaN 填充版，便于与本地 step6 对齐）----
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    preds_full = np.full((n_all, len(var_genes)), np.nan, dtype=np.float32)
    preds_full[:, keep] = preds
    true_full = np.full((n_all, len(var_genes)), np.nan, dtype=np.float32)
    true_full[:, keep] = true_expr
    np.savez_compressed(
        args.out + "_predictions.npz",
        pred_expr=preds,
        pred_expr_full=preds_full,
        pred_perturbations=np.array(holdout_genes),
        pred_genes=np.array(matched_genes),
        pred_gene_indices_in_h5ad=keep.astype(np.int32),
        target_in_input_space=target_in_space,
        true_expr=true_expr,
        true_expr_full=true_full,
        ctrl_mean_matched=ctrl_m,
        ctrl_mean_full=ctrl_mean_full,
        full_gene_space=np.array(var_genes),
    )
    lab.to_csv(args.out + "_labels.csv", index=False)

    # ---- 8. 元信息 ----
    with open(os.path.join(REPO, "GIT_COMMIT")) as f:
        git_commit = f.read().strip()
    import subprocess
    meta = {
        "task": "scGenePT cross-model generalization: Replogle essential "
                "held-out perturbation prediction (zero-shot from Norman finetune)",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {
            "name": "scGenePT (scgenept_ncbi+uniprot, Norman finetuned)",
            "model_type": MODEL_TYPE,
            "repo": "https://github.com/czi-ai/scGenePT",
            "git_commit": git_commit,
            "weights_file": MODEL_WEIGHTS,
            "weights_sha256": sha256(MODEL_WEIGHTS),
            "weights_bytes": os.path.getsize(MODEL_WEIGHTS),
            "gene_embeddings_file": MODELS_DIR +
                "gene_embeddings/NCBI+UniProt_embeddings-gpt3.5-ada.pkl",
            "inference_mode": "deterministic mean (explicit_zero_prob=False)",
        },
        "data": {
            "h5ad": args.h5ad,
            "h5ad_sha256": sha256(args.h5ad),
            "heldout_source": args.local_npz,
            "split": "GEARS simulation, split_seed=0, holdout frac=0.2",
            "n_genes_h5ad": len(var_genes),
            "n_genes_matched_vocab": int(n_matched),
            "n_ctrl_cells_total": n_ctrl,
            "ctrl_pool_size": int(len(ctrl_pool)),
        },
        "eval_protocol": {
            "label_definition": "pearson_delta = Pearson(pred_delta, true_delta) "
                                "on top-50 |true_delta| genes; pearson_delta_all on "
                                "all matched genes (same as local step1_run.py)",
            "true_mean": "mean of condition cells in perturb_processed.h5ad",
            "ctrl_mean": "mean of ctrl cells in perturb_processed.h5ad",
            "top_n_deg": TOP_N_DEG,
            "cross_dataset_note": "model finetuned on Norman 5045-gene space; "
                                  "applied zero-shot via scGPT whole-human vocab "
                                  "gene matching",
            "target_not_in_input_space": (
                "部分 held-out 目标基因不在数据集输入基因内（GEARS 图含扰动节点而 "
                "token 模型无法标记缺失 token），这些扰动的预测≈ctrl，逐扰动留痕 "
                "target_in_input_space"),
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
        "summary": {
            "pearson_delta_mean": float(valid.pearson_delta.mean()),
            "pearson_delta_median": float(valid.pearson_delta.median()),
            "pearson_delta_all_mean": float(valid.pearson_delta_all.mean()),
            "n_valid": int(len(valid)),
        },
    }
    with open(args.out + "_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    if not args.limit and os.path.exists(CKPT):
        os.remove(CKPT)          # 完成后清理断点
    print(f"[out] -> {args.out}_predictions.npz / _labels.csv / _meta.json "
          f"(total {time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
