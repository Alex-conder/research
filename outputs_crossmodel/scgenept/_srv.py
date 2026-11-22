# -*- coding: utf-8 -*-
"""服务器 SSH/SFTP 辅助（paramiko）。用法：
    python _srv.py exec '<shell 命令>' [timeout]
    python _srv.py put <本地路径> <远端路径>
    python _srv.py get <远端路径> <本地路径>
    python _srv.py probe   # 连接 + 基础信息
"""
import sys
import time

import paramiko

HOST, USER, PWD = "172.16.92.2", "root", "lab123lab"


def connect():
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, username=USER, password=PWD, timeout=20)
    return c


def exec_cmd(c, cmd, timeout=120, get_pty=False):
    i, o, e = c.exec_command(cmd, timeout=timeout, get_pty=get_pty)
    out = o.read().decode("utf-8", "replace")
    err = e.read().decode("utf-8", "replace")
    rc = o.channel.recv_exit_status()
    return rc, out, err


def progress(transferred, total):
    if total > 0:
        pct = transferred / total * 100
        print(f"\r  {transferred/1e6:.0f}/{total/1e6:.0f} MB ({pct:.1f}%)",
              end="", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    c = connect()
    if mode == "exec":
        cmd = sys.argv[2]
        timeout = int(sys.argv[3]) if len(sys.argv) > 3 else 120
        rc, out, err = exec_cmd(c, cmd, timeout)
        sys.stdout.write(out)
        if err.strip():
            sys.stderr.write(err)
        print(f"[rc={rc}]")
    elif mode == "put":
        t0 = time.time()
        sftp = c.open_sftp()
        sftp.put(sys.argv[2], sys.argv[3], callback=progress)
        dt = time.time() - t0
        import os
        sz = os.path.getsize(sys.argv[2])
        print(f"\nput {sz/1e6:.0f} MB in {dt:.0f}s = {sz/1e6/dt:.1f} MB/s")
    elif mode == "get":
        t0 = time.time()
        sftp = c.open_sftp()
        sftp.get(sys.argv[2], sys.argv[3], callback=progress)
        dt = time.time() - t0
        print(f"\nget done in {dt:.0f}s")
    elif mode == "probe":
        for cmd in [
            "cat /proc/loadavg",
            "df -h /opt / | tail -3",
            'docker ps --format "{{.Names}}\t{{.Status}}"',
            "python3 --version; docker --version",
            "curl -sI --max-time 10 https://czi-scgenept-public.s3.amazonaws.com/ "
            "-o /dev/null -w 's3_http=%{http_code} time=%{time_total}s\\n'",
            "curl -sI --max-time 10 https://hf-mirror.com -o /dev/null "
            "-w 'hf_http=%{http_code}\\n'",
            "nproc; free -g | head -2",
        ]:
            rc, out, err = exec_cmd(c, cmd)
            print(f"$ {cmd}\n{out.strip()}")
            if rc:
                print(f"  [rc={rc}] {err.strip()[:300]}")
    c.close()
