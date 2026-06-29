#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
断点续传下载 GEO 原始数据，下载完成后校验文件大小
支持 FTP/HTTP，可设置重试次数、块大小、超时
用法：
    python download_geo_resumable.py <url> <output_path> [chunk_size] [max_retries]
示例：
    python download_geo_resumable.py \
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227651/suppl/GSE227651_RAW.tar" \
        "../01_raw_data/GSE227651_RAW.tar" \
        1048576 100
"""
import sys
import os
import time
import urllib.request


def get_total_size(url: str) -> int:
    """通过 HEAD 请求获取总大小；FTP 返回 -1"""
    try:
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=30) as resp:
            return int(resp.headers.get('Content-Length', -1))
    except Exception:
        return -1


def download(url: str, out_path: str, chunk_size: int = 1024 * 1024, max_retries: int = 100):
    out_path = os.path.abspath(out_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    total_size = get_total_size(url)
    print(f"[INFO] Reported total size: {total_size / 1024 / 1024:.2f} MB" if total_size > 0 else "[INFO] Total size unknown")

    retries = 0
    while retries < max_retries:
        start_byte = os.path.getsize(out_path) if os.path.exists(out_path) else 0
        print(f"[INFO] Output: {out_path}")
        print(f"[INFO] Resume from byte: {start_byte}")

        # 如果已完成且大小匹配，直接返回
        if total_size > 0 and start_byte >= total_size:
            print("[SUCCESS] File already complete.")
            return

        try:
            req = urllib.request.Request(url)
            if start_byte > 0:
                req.add_header("Range", f"bytes={start_byte}-")
            with urllib.request.urlopen(req, timeout=30) as resp:
                # 获取剩余大小
                remaining = int(resp.headers.get("Content-Length", 0))
                expected_total = start_byte + remaining if start_byte > 0 else remaining
                print(f"[INFO] Expected total (this session): {expected_total / 1024 / 1024:.2f} MB")
                mode = "ab" if start_byte > 0 else "wb"
                with open(out_path, mode) as f:
                    downloaded = start_byte
                    last_time = time.time()
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        now = time.time()
                        if now - last_time >= 10:
                            pct = (downloaded / expected_total * 100) if expected_total > 0 else 0
                            print(f"[PROGRESS] {downloaded / 1024 / 1024:.2f} MB / {expected_total / 1024 / 1024:.2f} MB ({pct:.1f}%)")
                            last_time = now

            # 下载完成后校验
            final_size = os.path.getsize(out_path)
            print(f"[INFO] Final size: {final_size / 1024 / 1024:.2f} MB")
            if total_size > 0 and final_size < total_size:
                print(f"[WARN] Incomplete: {final_size}/{total_size}. Will retry...")
                retries += 1
                time.sleep(min(2 ** retries, 30))
                continue
            print("[SUCCESS] Download completed and verified.")
            return
        except Exception as e:
            retries += 1
            wait = min(2 ** retries, 30)
            print(f"[ERROR] {e}. Retrying {retries}/{max_retries} in {wait}s...")
            time.sleep(wait)

    print("[FAIL] Max retries exceeded.")
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    url = sys.argv[1]
    out = sys.argv[2]
    chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 1024 * 1024
    retries = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    download(url, out, chunk, retries)
