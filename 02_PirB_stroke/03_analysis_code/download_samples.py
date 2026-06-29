#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量下载 GEO 样本级别补充文件
"""
import sys
import os
import time
import urllib.request


def download(url: str, out_path: str, max_retries: int = 50):
    out_path = os.path.abspath(out_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    retries = 0
    while retries < max_retries:
        start_byte = os.path.getsize(out_path) if os.path.exists(out_path) else 0
        try:
            req = urllib.request.Request(url)
            if start_byte > 0:
                req.add_header("Range", f"bytes={start_byte}-")
            with urllib.request.urlopen(req, timeout=30) as resp:
                mode = "ab" if start_byte > 0 else "wb"
                with open(out_path, mode) as f:
                    while True:
                        chunk = resp.read(1024 * 1024)
                        if not chunk:
                            break
                        f.write(chunk)
            print(f"[SUCCESS] {out_path}")
            return
        except Exception as e:
            retries += 1
            wait = min(2 ** retries, 30)
            print(f"[ERROR] {url}: {e}. Retry {retries}/{max_retries} in {wait}s...")
            time.sleep(wait)
    print(f"[FAIL] {url}")


if __name__ == "__main__":
    # GSE227651 remaining samples
    base_dir = "../01_raw_data/GSE227651"
    samples = {
        "GSM7104633": "3d",
        "GSM7104634": "7d",
        "GSM7104635": "sham",
    }
    for gsm, label in samples.items():
        out_dir = os.path.join(base_dir, f"{gsm}_{label}")
        for ftype in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
            url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{gsm}_{label}-{ftype}"
            out = os.path.join(out_dir, ftype)
            download(url, out)
