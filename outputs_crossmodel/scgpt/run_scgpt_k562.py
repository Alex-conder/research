# -*- coding: utf-8 -*-
"""scGPT whole-human 预训练权重 (bowang-lab/scGPT, perturblab/scgpt-human 镜像)
对 Replogle essential held-out 扰动的零样本推理（K562）。

协议（canonical scGPT 零样本扰动预测，对齐 scGPT 论文 / Tutorial_Perturbation）：
- 输入 = ctrl 细胞池（pool=32, seed=0）的 log-normalized 表达，连续值；
- 扰动基因在输入中 mask（值置 -1, mask_value），模型经 MLM 头 (mlm_output)
  预测全部 token 的值；pred_expr = 池内 mlm_output 均值；
- 未扰动位置输出≈输入（恒等重建），扰动位置为真实预测；同时单独留痕
  pred_pert_gene_value（扰动基因处的预测值，NaN=目标不在输入空间）。
- 基因经 scGPT whole-human vocab (60697) 匹配；与 scgenept 同序同空间。

口径与 run_scgenept_k562.py 对齐：held-out=meta/pred_expr_seed0.npz；
pearson_delta = top-50 |true_delta| 基因上 Pearson(pred_delta, true_delta)。

断点续跑：每扰动写 <out>_checkpoint.npz（指纹=基因空间+扰动名单+pool+batch）。

用法（容器内）：
  python /research/scripts/run_scgpt_k562.py \
      --h5ad /research/data/perturb_processed_k562.h5ad \
      --local-npz /research/meta/pred_expr_seed0.npz \
      --out /research/results/scgpt_k562 [--pool 32] [--batch 8] [--threads 64]
产物：<out>_predictions.npz + _labels.csv + _meta.json + _genes.csv（+ 运行中 _checkpoint.npz）
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
DEFAULT_OUT = "/research/results/scgpt_k562"

TOP_N_DEG = 50
SEED = 0
MASK_VALUE = -1.0          # scGPT mask_value（pad_value=-2，互不相同）
PAD_VALUE = -2.0

# whole-human args.json（perturblab/scgpt-human）
CFG = dict(embsize=512, nhead=8, d_hid=512, nlayers=12,
           nlayers_cls=3, n_cls=1, mvc_decoder_style="inner product")


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
    """构建 scGPT TransformerModel 并加载 whole-human 权重（自适应键名）。"""
    import json as _json
    from scgpt.model.model import TransformerModel
    from scgpt.tokenizer import GeneVocab

    with open(VOCAB_FILE) as f:
        vocab = GeneVocab.from_file(VOCAB_FILE)
    if "<pad>" not in vocab:
        vocab.append_token("<pad>")

    model = TransformerModel(
        ntoken=len(vocab),
        d_model=CFG["embsize"], nhead=CFG["nhead"], d_hid=CFG["d_hid"],
        nlayers=CFG["nlayers"], nlayers_cls=CFG["nlayers_cls"], n_cls=CFG["n_cls"],
        vocab=vocab, dropout=0.0, pad_token="<pad>", pad_value=PAD_VALUE,
        do_mvc=True, mvc_decoder_style=CFG["mvc_decoder_style"],
        use_fast_transformer=False)
    ckpt = torch.load(MODEL_WEIGHTS, weights_only=True, map_location=device)
    sd = ckpt.get("model_state_dict", ckpt) if isinstance(ckpt, dict) else ckpt
    # fast_transformer 训练保存的键名 Wqkv.* -> 标准 MultiheadAttention in_proj_*
    # （scgpt 0.2.1 用 nn.TransformerEncoder，键为 transformer_encoder.layers.N.self_attn.Wqkv.*）
    sd = {k.replace(".self_attn.Wqkv.", ".self_attn.in_proj_"): v
          for k, v in sd.items()}
    missing, unexpected = model.load_state_dict(sd, strict=False)
    crit_miss = [k for k in missing
                 if not k.startswith(("cls_decoder", "batch_encoder",
                                      "mvc_decoder", "dab_decoder",
                                      "ecs_decoder"))]
    if unexpected or crit_miss:
        print(f"[model] WARN unexpected={len(unexpected)} crit_missing={crit_miss}",
              flush=True)
        if unexpected[:5]:
            print(f"[model] unexpected sample: {unexpected[:5]}", flush=True)
    model.to(device)
    model.eval()
    return model, vocab


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
    print(f"[cfg] h5ad={args.h5ad} pool={args.pool} batch={args.batch} "
          f"threads={args.threads} device=cpu", flush=True)

    # ---- 1. held-out 名单 ----
    z = np.load(args.local_npz)
    holdout_genes = [str(g) for g in z["pred_perturbations"]]
    local_genes = [str(g) for g in z["pred_genes"]]
    print(f"[data] held-out perts={len(holdout_genes)}, "
          f"local pred gene space={len(local_genes)}", flush=True)

    # ---- 2. 读 h5ad ----
    adata = ad.read_h5ad(args.h5ad, backed="r")
    var_genes = adata.var["gene_name"].astype(str).tolist()
    assert var_genes == local_genes, "服务器 h5ad 基因序与本地 npz pred_genes 不一致！"
    obs_cond = adata.obs["condition"].astype(str).values
    is_ctrl = obs_cond == "ctrl"
    n_ctrl = int(is_ctrl.sum())
    print(f"[data] cells={adata.n_obs}, genes={adata.n_vars}, ctrl_cells={n_ctrl}",
          flush=True)
    ctrl_mean_full = np.asarray(
        adata[is_ctrl].X.mean(axis=0)).ravel().astype(np.float32)

    # ---- 断点 ----
    CKPT = args.out + "_checkpoint.npz"
    resume_rows = []
    pool_idx_ckpt = None
    if args.limit == 0 and os.path.exists(CKPT):
        cz = np.load(CKPT, allow_pickle=True)
        if (list(cz["ckpt_gene_space"]) == var_genes
                and list(cz["ckpt_perts"]) == holdout_genes
                and int(cz["ckpt_pool"]) == args.pool
                and int(cz["ckpt_batch"]) == args.batch):
            resume_rows = [int(i) for i in cz["ckpt_rows"]]
            pool_idx_ckpt = cz["ckpt_pool_idx"]
            print(f"[ckpt] 恢复 {len(resume_rows)} 个已完成扰动，"
                  f"剩余 {len(holdout_genes) - len(resume_rows)}", flush=True)
        else:
            print("[ckpt] 指纹不匹配，忽略断点从头开始", flush=True)

    rng = np.random.default_rng(SEED)
    pool_idx = pool_idx_ckpt if pool_idx_ckpt is not None else \
        rng.choice(np.where(is_ctrl)[0], size=min(args.pool, n_ctrl),
                   replace=False)
    t0 = time.time()
    ctrl_pool = np.asarray(adata[pool_idx].X.todense()).astype(np.float32)
    print(f"[data] ctrl pool {ctrl_pool.shape} densified in {time.time()-t0:.0f}s",
          flush=True)

    # ---- 3. 模型 + vocab 匹配 ----
    model, vocab = load_model(device)
    n_params = sum(p.numel() for p in model.parameters())
    in_vocab = np.array([g in vocab for g in var_genes])
    keep = np.where(in_vocab)[0]
    n_matched = len(keep)
    print(f"[vocab] matched {n_matched}/{len(var_genes)} genes "
          f"(vocab size={len(vocab)}, params={n_params/1e6:.0f}M)", flush=True)
    matched_genes = [var_genes[i] for i in keep]
    gene_ids = np.array([vocab[g] for g in matched_genes], dtype=np.int64)
    ctrl_pool = ctrl_pool[:, keep]
    pos_in_matched = {g: i for i, g in enumerate(matched_genes)}
    target_in_space = np.array([g in pos_in_matched for g in holdout_genes])
    print(f"[vocab] 目标基因不在输入空间的扰动: {int((~target_in_space).sum())}/"
          f"{len(holdout_genes)}", flush=True)

    # ---- 4. 逐扰动推理 ----
    gene_ids_t = torch.tensor(gene_ids).long().unsqueeze(0).to(device)
    pad_mask = torch.zeros_like(gene_ids_t, dtype=torch.bool, device=device)

    n_all = len(holdout_genes)
    if args.limit:
        run_rows = list(range(min(args.limit, n_all)))
        resume_rows = []
    else:
        done_set = set(resume_rows)
        run_rows = [i for i in range(n_all) if i not in done_set]
    preds = np.full((n_all, n_matched), np.nan, dtype=np.float32)
    pred_pert_val = np.full(n_all, np.nan, dtype=np.float32)
    if resume_rows:
        cz = np.load(CKPT, allow_pickle=True)
        preds[np.array(resume_rows)] = cz["ckpt_preds"]
        pred_pert_val[np.array(resume_rows)] = cz["ckpt_pert_val"]

    done_rows = list(resume_rows)
    t_inf0 = time.time()
    for k, pi in enumerate(run_rows):
        g = holdout_genes[pi]
        ti = time.time()
        x = ctrl_pool.copy()
        ppos = pos_in_matched.get(g, None)
        if ppos is not None:
            x[:, ppos] = MASK_VALUE       # 扰动基因 mask
        xt = torch.from_numpy(x).float()
        outs = []
        with torch.no_grad():
            for s in range(0, len(xt), args.batch):
                xb = xt[s:s + args.batch]
                gi = gene_ids_t.expand(len(xb), -1)
                pm = pad_mask.expand(len(xb), -1)
                od = model(gi.to(device), xb.to(device),
                           src_key_padding_mask=pm.to(device),
                           CLS=False, CCE=False, MVC=False, ECS=False,
                           do_sample=False)
                outs.append(od["mlm_output"].float().cpu().numpy())
        p = np.concatenate(outs, axis=0).mean(axis=0)
        preds[pi] = p
        if ppos is not None:
            pred_pert_val[pi] = p[ppos]
        done_rows.append(pi)
        if not args.limit:
            np.savez(CKPT, ckpt_rows=np.array(done_rows),
                     ckpt_preds=preds[np.array(done_rows)],
                     ckpt_pert_val=pred_pert_val[np.array(done_rows)],
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

    # ---- 5. 真值与标签 ----
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

    # ---- 6. 产物（与 scgenept_k562_predictions.npz 同构）----
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
        # scGPT 特有：扰动基因处预测值（NaN=不在输入空间）；输出类型留痕
        pred_pert_gene_value=pred_pert_val,
        output_type=np.array(["gene_expression_prediction_via_masked_perturbation"]),
    )
    lab.to_csv(args.out + "_labels.csv", index=False)
    pd.DataFrame({"matched_genes": matched_genes}).to_csv(
        args.out + "_genes_matched.csv", index=False)
    pd.DataFrame({"full_gene_space": var_genes}).to_csv(
        args.out + "_genes_full.csv", index=False)

    # ---- 7. 元信息 ----
    meta = {
        "task": "scGPT whole-human zero-shot: Replogle essential held-out "
                "perturbation prediction (K562)",
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model": {
            "name": "scGPT whole-human pretrained (bowang-lab/scGPT)",
            "source": "hf-mirror perturblab/scgpt-human (re-upload of "
                      "bowang-lab/scGPT whole_human best_model.pt)",
            "weights_file": MODEL_WEIGHTS,
            "weights_sha256": sha256(MODEL_WEIGHTS),
            "weights_bytes": os.path.getsize(MODEL_WEIGHTS),
            "args_json": CFG,
            "inference_protocol": "canonical zero-shot masking: perturbed gene "
                "value masked to -1 in ctrl cells, mlm_output read for all "
                "tokens; pred_expr = pool mean; unmasked positions are "
                "identity reconstructions",
            "output_type": "gene expression prediction (mlm_output)",
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
                "on top-50 |true_delta| genes; pearson_delta_all on all matched "
                "genes (same as local step1_run.py / scgenept run)",
            "true_mean": "mean of condition cells in perturb_processed.h5ad",
            "ctrl_mean": "mean of ctrl cells in perturb_processed.h5ad",
            "top_n_deg": TOP_N_DEG,
            "caveat": "零样本 mask 协议下未扰动位置预测≈输入，pearson_delta 主要反映"
                "扰动基因单点及其邻域；pred_pert_gene_value 单独留痕",
            "target_not_in_input_space": "目标基因不在输入空间的扰动无法 mask，"
                "预测≈ctrl（138/217，同 scgenept 数据事实）",
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
        os.remove(CKPT)
    print(f"[out] -> {args.out}_predictions.npz / _labels.csv / _meta.json "
          f"(total {time.time()-t_start:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
