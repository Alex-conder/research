# -*- coding: utf-8 -*-
"""Tiny paramiko helper: run remote commands / SFTP on 172.16.92.2.
Usage: C:/Python314/python.exe _ssh.py "remote command" [timeout_sec]
       C:/Python314/python.exe _ssh.py put local remote
       C:/Python314/python.exe _ssh.py get remote local
"""
import sys

import paramiko

HOST, USER, PWD = "172.16.92.2", "root", "lab123lab"


def connect():
    c = paramiko.SSHClient()
    c.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    c.connect(HOST, username=USER, password=PWD, timeout=30)
    return c


def run(cmd, timeout=120):
    c = connect()
    i, o, e = c.exec_command(cmd, timeout=timeout)
    out = o.read().decode(errors="replace")
    err = e.read().decode(errors="replace")
    rc = o.channel.recv_exit_status()
    c.close()
    if err.strip():
        out += "<STDERR>" + err
    return rc, out


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "put":
        c = connect()
        sftp = c.open_sftp()
        sftp.put(sys.argv[2], sys.argv[3])
        sftp.close()
        c.close()
        print("PUT_OK", sys.argv[3])
    elif mode == "get":
        c = connect()
        sftp = c.open_sftp()
        sftp.get(sys.argv[2], sys.argv[3])
        sftp.close()
        c.close()
        print("GET_OK", sys.argv[3])
    else:
        rc, out = run(mode, int(sys.argv[2]) if len(sys.argv) > 2 else 120)
        print(out)
        sys.exit(rc)
