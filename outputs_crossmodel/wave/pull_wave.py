# -*- coding: utf-8 -*-
"""Batch SFTP pull of wave artifacts from ah-server (172.16.92.2), resuming by size."""
import sys
import time
from pathlib import Path

import paramiko

HOST, PORT = "172.16.2.7", 22022
KEY = r"E:\animal_house_build\.ssh\animal_house_deploy2"


def connect():
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, port=PORT, username="root", key_filename=KEY, timeout=30)
    return c

REMOTE = "/opt/research/outputs_crossmodel"
LOCAL = Path(__file__).resolve().parent

FILES = {  # relpath -> expected size (bytes), from `stat -c %s` on server 2026-10-24
    "scgpt/scgpt_wave_k562_predictions.npz": 27218092,
    "scgpt/scgpt_wave_k562_labels.csv": 19351,
    "scgpt/scgpt_wave_k562_meta.json": 1062,
    "scgpt/scgpt_wave_rpe1_predictions.npz": 38084233,
    "scgpt/scgpt_wave_rpe1_labels.csv": 27518,
    "scgpt/scgpt_wave_rpe1_meta.json": 1169,
    "geneformer/geneformer_wave_k562_predictions.npz": 40746945,
    "geneformer/geneformer_wave_k562_labels.csv": 19366,
    "geneformer/geneformer_wave_k562_meta.json": 932,
    "geneformer/geneformer_wave_rpe1_predictions.npz": 57020338,
    "geneformer/geneformer_wave_rpe1_labels.csv": 27524,
    "geneformer/geneformer_wave_rpe1_meta.json": 1040,
    "uce/uce_wave_k562_predictions.npz": 68028636,
    "uce/uce_wave_k562_labels.csv": 18880,
    "uce/uce_wave_k562_meta.json": 896,
    "uce/uce_wave_rpe1_predictions.npz": 95191097,
    "uce/uce_wave_rpe1_labels.csv": 20789,
    "uce/uce_wave_rpe1_meta.json": 1003,
}

c = connect()
sftp = c.open_sftp()
t0 = time.time()
done_bytes = 0
for rel, size in FILES.items():
    dst = LOCAL / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size == size:
        print(f"SKIP {rel} (already complete)")
        done_bytes += size
        continue
    ta = time.time()
    sftp.get(f"{REMOTE}/{rel}", str(dst))
    dt = time.time() - ta
    got = dst.stat().st_size
    status = "OK" if got == size else f"SIZE_MISMATCH got={got} want={size}"
    print(f"GET  {rel} {got/1e6:.1f}MB {dt:.1f}s {got/max(dt,1e-9)/1e6:.2f}MB/s {status}")
    done_bytes += got
sftp.close()
c.close()
total = time.time() - t0
print(f"TOTAL {done_bytes/1e6:.1f}MB in {total:.1f}s ({done_bytes/max(total,1e-9)/1e6:.2f}MB/s incl. skips)")
