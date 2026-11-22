# -*- coding: utf-8 -*-
"""scFoundation (xTrimoGene MAE, 100M) 单扰动 CPU 推理测速。

权重: perturblab/scfoundation-gene (model.pt, 477MB, Apache-2.0,
      源自 biomap-research/scFoundation, Nature Methods 2024)
代码: biomap-research/scFoundation @ main 的 model/ 目录（官方加载约定）。

测速口径：模拟真实推理 1 个 ctrl 细胞 → 重建全部 19264 基因表达
（与官方 get_embedding.py 路径一致：encoder 只过非零基因，decoder 过全基因）。
输入的非零模式取 K562 ctrl 均值映射到 19264 基因空间后的真实支撑集，
值为 1..99 的模拟分箱值（ timing 与数值语义无关）。

用法（容器内）：python /research/scripts/scf_speedtest.py
产物：/research/results/scfoundation_speedtest.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch

CODE_DIR = "/research/scfoundation_repo/model"
MODEL_PT = "/research/models/scfoundation_gene_model.pt"
H5AD = "/research/data/perturb_processed_k562.h5ad"
OUT = "/research/results/scfoundation_speedtest.json"
THREADS = int(os.environ.get("SCF_THREADS", "8"))

sys.path.insert(0, CODE_DIR)
from load import getEncoerDecoderData  # noqa: E402
from pretrainmodels import select_model  # noqa: E402

N_GENES = 19264


def main():
    torch.set_num_threads(THREADS)
    os.environ["OMP_NUM_THREADS"] = str(THREADS)
    device = torch.device("cpu")

    # ---- 1. 构建模型（配置来自 perturblab/scfoundation-gene config.json）----
    config = {
        "model": "mae_autobin",
        "n_class": N_GENES,
        "seq_len": 19266,   # = 权重 pos_emb 行数 19267 - 1（perturblab 打包的实际值）
        "bin_num": 100,
        "bin_alpha": 1.0,
        "pad_token_id": 103,
        "mask_token_id": 102,
        "encoder": {"module_type": "transformer", "hidden_dim": 768,
                    "depth": 12, "heads": 12},
        "decoder": {"module_type": "performer", "hidden_dim": 512,
                    "depth": 6, "heads": 8, "dim_head": 64},
    }
    model = select_model(config)
    sd = torch.load(MODEL_PT, map_location="cpu", weights_only=False)
    pos_rows = sd["pos_emb.weight"].shape[0]
    print(f"[model] pos_emb rows={pos_rows} -> max_seq_len={pos_rows-1}",
          flush=True)
    missing, unexpected = model.load_state_dict(sd, strict=False)
    print(f"[model] missing={len(missing)} unexpected={len(unexpected)}",
          flush=True)
    if missing[:5]:
        print("  missing sample:", missing[:5], flush=True)
    if unexpected[:5]:
        print("  unexpected sample:", unexpected[:5], flush=True)
    model.to(device)
    model.eval()

    # ---- 2. 构造 1 个细胞的输入（真实非零支撑集）----
    gene_df = pd.read_csv(os.path.join(CODE_DIR, "OS_scRNA_gene_index.19264.tsv"),
                          sep="\t")
    scf_genes = gene_df["gene_name"].astype(str).tolist()
    assert len(scf_genes) == N_GENES, len(scf_genes)

    import anndata as ad
    adata = ad.read_h5ad(H5AD, backed="r")
    var_genes = adata.var["gene_name"].astype(str).tolist()
    ctrl_mean = np.asarray(adata[adata.obs.condition == "ctrl"].X.mean(axis=0))
    adata.file.close()
    ctrl_mean = np.asarray(ctrl_mean).ravel()
    g2i = {g: i for i, g in enumerate(var_genes)}
    mapped = [g2i[g] for g in scf_genes if g in g2i]
    support = np.zeros(N_GENES, dtype=bool)
    support[mapped] = ctrl_mean[[g2i[g] for g in scf_genes if g in g2i]] > 0
    n_nz = int(support.sum())
    print(f"[data] scf gene space={N_GENES}, K562 映射基因={len(mapped)}, "
          f"非零支撑={n_nz}", flush=True)

    rng = np.random.default_rng(0)
    values = np.zeros((1, N_GENES), dtype=np.float32)
    values[0, support] = rng.integers(1, 100, size=n_nz).astype(np.float32)
    data = torch.from_numpy(values).to(device)
    data_raw = data.clone()

    (encoder_data, encoder_pos, encoder_pad, encoder_labels,
     decoder_data, decoder_pad, _, _, decoder_pos) = \
        getEncoerDecoderData(data, data_raw, config)
    enc_len = int(encoder_data.shape[1])
    print(f"[data] encoder seq len={enc_len}, decoder seq len={N_GENES}",
          flush=True)

    def forward_once():
        with torch.no_grad():
            return model(x=encoder_data, padding_label=encoder_pad,
                         encoder_position_gene_ids=encoder_pos,
                         encoder_labels=encoder_labels,
                         decoder_data=decoder_data, mask_gene_name=False,
                         mask_labels=None,
                         decoder_position_gene_ids=decoder_pos,
                         decoder_data_padding_labels=decoder_pad)

    # ---- 3. 计时 ----
    t0 = time.time()
    out = forward_once()
    t_warm = time.time() - t0
    n_out = int(out.numel())
    print(f"[run] warmup forward {t_warm:.1f}s, output shape={list(out.shape)}",
          flush=True)

    # encoder 单独计时
    with torch.no_grad():
        x = model.token_emb(torch.unsqueeze(encoder_data, 2), output_weight=0)
        x = x + model.pos_emb(encoder_pos)
        t0 = time.time()
        enc_out = model.encoder(x, padding_mask=encoder_pad)
        t_enc = time.time() - t0
    print(f"[run] encoder-only {t_enc:.1f}s", flush=True)

    times = []
    for _ in range(2):
        t0 = time.time()
        forward_once()
        times.append(time.time() - t0)
    t_full = float(np.mean(times))
    print(f"[run] full forward mean={t_full:.1f}s "
          f"(runs={[round(t, 1) for t in times]})", flush=True)

    meta = {
        "model": "scFoundation xTrimoGene MAE (gene model, 100M params)",
        "weights": "hf-mirror: perturblab/scfoundation-gene model.pt (477MB)",
        "source_repo": "biomap-research/scFoundation (Nature Methods 2024)",
        "hardware": "AMD EPYC 7702, torch CPU",
        "threads": THREADS,
        "note": ("测速时 scGenePT 全量任务同机运行（容器 32 线程配额）；"
                 "THREADS=8 避免互相干扰，绝对值偏保守"),
        "input": {"n_genes_total": N_GENES, "n_nonzero_encoder_input": n_nz,
                  "encoder_seq_len": enc_len,
                  "support": "K562 ctrl 均值>0 的真实基因支撑集(映射到19264)"},
        "timings_seconds": {"warmup_forward": t_warm, "encoder_only": t_enc,
                            "full_forward_mean": t_full,
                            "full_forward_runs": times},
        "extrapolation_per_perturbation": {
            "assumption": "每扰动需对 P 个 ctrl 细胞各 1 次 full forward 再取均值",
            "for_pool_32": t_full * 32,
            "for_pool_8": t_full * 8,
        },
        "output_shape": list(out.shape),
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"[out] -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
