"""
使用 ftplib 进行可续传的 GEO FTP 下载
- 支持断点续传（REST）
- 自动重试
"""
import os, sys, time
from ftplib import FTP, error_temp, error_perm

def download_ftp(host, remote_path, local_path, max_retries=50):
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    start_byte = os.path.getsize(local_path) if os.path.exists(local_path) else 0
    print(f"[INFO] Download {remote_path}")
    print(f"[INFO] Local: {local_path}, resume from {start_byte}")
    
    retries = 0
    while retries < max_retries:
        try:
            ftp = FTP(host, timeout=60)
            ftp.login()
            ftp.cwd(os.path.dirname(remote_path))
            filename = os.path.basename(remote_path)
            total = ftp.size(filename)
            print(f"[INFO] Total size: {total / 1024 / 1024:.2f} MB")
            
            if start_byte >= total:
                print("[SUCCESS] Already complete.")
                ftp.quit()
                return
            
            mode = "ab" if start_byte > 0 else "wb"
            with open(local_path, mode) as f:
                if start_byte > 0:
                    ftp.voidcmd(f"REST {start_byte}")
                
                downloaded = start_byte
                last_time = time.time()
                last_size = downloaded
                
                def callback(chunk):
                    nonlocal downloaded, last_time, last_size
                    f.write(chunk)
                    downloaded += len(chunk)
                    now = time.time()
                    if now - last_time >= 10:
                        pct = downloaded / total * 100
                        speed = (downloaded - last_size) / (now - last_time) / 1024
                        print(f"[PROGRESS] {downloaded/1024/1024:.2f}/{total/1024/1024:.2f} MB ({pct:.1f}%) {speed:.1f} KB/s")
                        last_time = now
                        last_size = downloaded
                
                ftp.retrbinary(f"RETR {filename}", callback, rest=start_byte)
            
            # 验证大小
            final_size = os.path.getsize(local_path)
            if final_size == total:
                print("[SUCCESS] Download completed and verified.")
                ftp.quit()
                return
            else:
                print(f"[WARN] Size mismatch: {final_size}/{total}")
                start_byte = final_size
                retries += 1
                ftp.quit()
                time.sleep(5)
        except Exception as e:
            retries += 1
            print(f"[ERROR] {e}. Retrying {retries}/{max_retries}...")
            start_byte = os.path.getsize(local_path) if os.path.exists(local_path) else 0
            time.sleep(min(2 ** retries, 60))
    
    print("[FAILED] Max retries reached.")

if __name__ == "__main__":
    remote = sys.argv[1]
    local = sys.argv[2]
    # remote should be like /geo/samples/GSM7104nnn/GSM7104633/suppl/GSM7104633_3d-matrix.mtx.gz
    download_ftp("ftp.ncbi.nlm.nih.gov", remote, local)
