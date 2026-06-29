"""
下载 GSE233815 子系列 scRNA-seq (GSE233812) 和 snRNA-seq (GSE233813)
以及 GSE233814 空间转录组 matrix/json（跳过 tif 图像）
"""
import os, sys, subprocess
from ftplib import FTP

raw_dir = "D:/Pirb_stroke_project/01_raw_data"

def run_ftp(remote_path, local_path):
    print(f"[DOWNLOAD] {remote_path} -> {local_path}")
    subprocess.run([sys.executable, "03_analysis_code/download_geo_ftp_resumable.py", remote_path, local_path])

# GSE233812 scRNA-seq RAW tar (166 MB)
run_ftp("/geo/series/GSE233nnn/GSE233812/suppl/GSE233812_RAW.tar", f"{raw_dir}/GSE233812_RAW.tar")

# GSE233813 snRNA-seq RAW tar (308 MB)
run_ftp("/geo/series/GSE233nnn/GSE233813/suppl/GSE233813_RAW.tar", f"{raw_dir}/GSE233813_RAW.tar")

# GSE233814 空间转录组：只下载 matrix/barcodes/features/json，跳过 tif
spatial_samples = {
    "GSM7437221": "C1-control",
    "GSM7437222": "B1-D1",
    "GSM7437223": "D1-D3",
    "GSM7437224": "C1-D7",
    "GSM7437225": "D1-D7",
}
out814 = f"{raw_dir}/GSE233814"
os.makedirs(out814, exist_ok=True)

for gsm, sample_suffix in spatial_samples.items():
    ftp_dir = f"/geo/samples/{gsm[:7]}nnn/{gsm}/suppl/"
    # 列出文件
    ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=30)
    ftp.login()
    ftp.cwd(ftp_dir)
    files = ftp.nlst()
    ftp.quit()
    
    for fname in files:
        if fname.endswith(".tif.gz"):
            continue
        if any(fname.endswith(ext) for ext in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz", "json.gz"]):
            remote = f"{ftp_dir}{fname}"
            local = f"{out814}/{fname}"
            run_ftp(remote, local)

print("[DONE] GSE233812/813/814 download attempt finished.")
