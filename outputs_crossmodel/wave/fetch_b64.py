# -*- coding: utf-8 -*-
"""Fetch server files via base64-over-SSH-stdout (SFTP write is broken).

Usage: C:/Python314/python.exe fetch_b64.py <remote_path> <local_path>
Writes remote `base64 -w0` stdout decoded to local_path, then compares
with remote md5sum. Exits nonzero on mismatch.
"""
import base64
import hashlib
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))  # _ssh.py lives in parent
from _ssh import connect  # noqa: E402


def main():
    remote, local = sys.argv[1], sys.argv[2]
    c = connect()
    i, o, e = c.exec_command(f"base64 -w0 {remote}; echo; md5sum {remote}",
                             timeout=600)
    raw = o.read().decode()
    err = e.read().decode()
    rc = o.channel.recv_exit_status()
    if rc != 0:
        print("REMOTE_ERR", err[:500])
        sys.exit(1)
    lines = raw.strip().split("\n")
    b64 = lines[0].strip()
    remote_md5 = lines[1].split()[0]
    data = base64.b64decode(b64)
    with open(local, "wb") as f:
        f.write(data)
    local_md5 = hashlib.md5(data).hexdigest()
    c.close()
    print(f"OK {local} {len(data)}B md5={'MATCH' if local_md5 == remote_md5 else 'MISMATCH ' + local_md5 + ' vs ' + remote_md5}")
    sys.exit(0 if local_md5 == remote_md5 else 2)


if __name__ == "__main__":
    main()
