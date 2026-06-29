"""
GSE225948 损坏 counts 文件 HTTPS/FTP 断点续传重下载。
自动检测 n_genes=0 或读取失败的 counts 文件，从 NCBI 重新下载。
"""
import os, subprocess, time, sys
import pandas as pd
from pathlib import Path

DATA_DIR = Path("D:/Pirb_stroke_project/01_raw_data/GSE225948")
LOG_PATH = DATA_DIR / "redownload_log.txt"

def log(msg):
    print(msg)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}\n")

def download_file(gsm, filename, max_time=600):
    local_path = DATA_DIR / filename
    # 备份旧损坏文件
    if local_path.exists():
        backup_path = local_path.with_suffix(local_path.suffix + ".old")
        try:
            local_path.rename(backup_path)
            log(f"[BACKUP] {filename} -> {backup_path.name}")
        except Exception as e:
            log(f"[WARN] backup failed: {e}")

    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/{filename}"
    cmd = [
        "curl", "-C", "-", "-L", "-o", str(local_path), url,
        "--max-time", str(max_time),
        "-w", "\nHTTP: %{http_code}  Size: %{size_download}\n"
    ]
    log(f"[DOWNLOAD] {url} -> {local_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    log(result.stdout)
    if result.stderr:
        log(f"[STDERR] {result.stderr[:500]}")
    return result.returncode

def validate_file(filename):
    fpath = DATA_DIR / filename
    if not fpath.exists():
        return False, "missing"
    try:
        n_genes = len(pd.read_csv(fpath, compression='gzip', usecols=[0]))
        if n_genes < 10000:
            return False, f"too_few_genes_{n_genes}"
        return True, f"OK_{n_genes}"
    except Exception as e:
        return False, f"error_{type(e).__name__}"

# 检测需要重下的文件
files_to_redownload = []
for f in sorted(DATA_DIR.glob("*_counts.csv.gz")):
    ok, reason = validate_file(f.name)
    if not ok:
        gsm = f.name.split('_')[0]
        files_to_redownload.append((gsm, f.name, reason))

log(f"[INFO] Files to redownload: {len(files_to_redownload)}")
for gsm, fn, reason in files_to_redownload:
    log(f"  {fn} ({reason})")

# 逐个下载并校验
for gsm, filename, reason in files_to_redownload:
    log(f"\n[START] {filename}")
    for attempt in range(3):
        rc = download_file(gsm, filename, max_time=900)
        ok, reason2 = validate_file(filename)
        if ok:
            log(f"[SUCCESS] {filename} -> {reason2}")
            break
        else:
            log(f"[RETRY {attempt+1}] {filename} still invalid: {reason2}")
            time.sleep(10)
    else:
        log(f"[FAILED] {filename} after 3 attempts")

log("[DONE]")
