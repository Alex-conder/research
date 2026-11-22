# -*- coding: utf-8 -*-
import sys

def load(path):
    d = {}
    for line in open(path, encoding="utf-8"):
        parts = line.split()
        if len(parts) == 2:
            d[parts[1].replace(chr(92), "/").lstrip("./").lstrip("*")] = parts[0]
    return d

loc = load(sys.argv[1])
srv = load(sys.argv[2])
allok = True
for k in sorted(srv):
    l = loc.get(k)
    ok = l == srv[k]
    allok = allok and ok
    print(("OK   " if ok else "FAIL ") + k + ("" if ok else " local=%s server=%s" % (l, srv[k])))
missing = [k for k in loc if k not in srv]
if missing:
    print("EXTRA_LOCAL:", missing)
print("MD5_OVERALL:", "PASS" if allok and not missing else "FAIL")
