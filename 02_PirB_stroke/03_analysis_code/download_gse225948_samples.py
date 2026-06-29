#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量下载 GSE225948 所有样本级别 CSV 文件
"""
import os
import re
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
    meta_path = "../01_raw_data/GSE225948_metadata_full.txt"
    out_dir = "../01_raw_data/GSE225948"
    with open(meta_path, 'r', encoding='utf-8') as f:
        text = f.read()
    supp_urls = re.findall(r'ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/samples/GSM\d+nnn/GSM\d+/suppl/[^\s]+\.gz', text)
    print(f"Found {len(supp_urls)} sample supplementary URLs")
    for url in supp_urls:
        filename = url.split('/')[-1]
        out = os.path.join(out_dir, filename)
        if os.path.exists(out):
            print(f"[SKIP] {out}")
            continue
        download(url, out)
