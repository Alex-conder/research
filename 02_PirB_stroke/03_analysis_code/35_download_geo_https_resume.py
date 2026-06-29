"""
使用 HTTPS + curl 断点续传下载 GEO 大文件
"""
import os, sys, subprocess, time

raw_dir = "D:/Pirb_stroke_project/01_raw_data"
os.makedirs(raw_dir, exist_ok=True)

def download_https(url, local_path, max_retries=100):
    local_path = os.path.abspath(local_path)
    print(f"[DOWNLOAD] {url} -> {local_path}")
    retries = 0
    while retries < max_retries:
        cmd = ["curl", "-C", "-", "-L", "-o", local_path, url,
               "--max-time", "600", "-w", "\\nHTTP: %{http_code}  Size: %{size_download}\\n"]
        print(f"[RUN] {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.returncode == 0:
            print(f"[SUCCESS] {local_path}")
            return True
        else:
            retries += 1
            wait = min(2 ** retries, 60)
            print(f"[ERROR] curl failed (code {result.returncode}). Retry {retries}/{max_retries} in {wait}s...")
            print(result.stderr[:500])
            time.sleep(wait)
    print(f"[FAILED] {url}")
    return False

# GSE233812 scRNA-seq RAW tar
f812 = os.path.join(raw_dir, "GSE233812_RAW.tar")
download_https("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233812/suppl/GSE233812_RAW.tar", f812)

# GSE233813 snRNA-seq RAW tar
f813 = os.path.join(raw_dir, "GSE233813_RAW.tar")
download_https("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233813/suppl/GSE233813_RAW.tar", f813)

# GSE233814 空间转录组：列出并下载非 tif 文件
print("[INFO] Listing GSE233814 spatial files...")
spatial_samples = {
    "GSM7437221": "C1-control",
    "GSM7437222": "B1-D1",
    "GSM7437223": "D1-D3",
    "GSM7437224": "C1-D7",
    "GSM7437225": "D1-D7",
}
out814 = os.path.join(raw_dir, "GSE233814")
os.makedirs(out814, exist_ok=True)

from ftplib import FTP
for gsm, suffix in spatial_samples.items():
    ftp_dir = f"/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/"
    try:
        ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=30)
        ftp.login()
        ftp.cwd(ftp_dir)
        files = ftp.nlst()
        ftp.quit()
    except Exception as e:
        print(f"[WARN] Cannot list {ftp_dir}: {e}")
        continue
    for fname in files:
        if fname.endswith(".tif.gz"):
            continue
        if any(fname.endswith(ext) for ext in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz", "json.gz"]):
            url = f"https://ftp.ncbi.nlm.nih.gov{ftp_dir}{fname}"
            local = os.path.join(out814, fname)
            download_https(url, local)

print("[DONE]")
