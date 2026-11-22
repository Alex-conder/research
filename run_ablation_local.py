# -*- coding: utf-8 -*-
"""GEARS 数据覆盖消融本机一键驱动（Windows）：A -> B -> D -> C 顺序执行。

单版本失败不中断后续（错误记录进日志并汇总）；所有输出追加到
outputs_ablation/train_all.log。用 pythonw 无窗口 detached 启动：

  pythonw run_ablation_local.py [EPOCHS]     （默认 10）

监控：tail -f outputs_ablation/train_all.log
"""
import os
import subprocess
import sys
import time

PROJECT = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(PROJECT, "venv", "Scripts", "python.exe")
SCRIPT = os.path.join(PROJECT, "train_gears_ablation_local.py")
LOG = os.path.join(PROJECT, "outputs_ablation", "train_all.log")
EPOCHS = int(sys.argv[1]) if len(sys.argv) > 1 else 10

os.makedirs(os.path.dirname(LOG), exist_ok=True)

D_WEIGHTS = os.path.join(PROJECT, "outputs_ablation", "A_k562",
                         "gears_model", "model.pt")
STEPS = [
    ("A_k562", EPOCHS, []),
    ("B_rpe1", EPOCHS, []),
    ("D_zeroshot", 0, ["--weights", D_WEIGHTS]),
    ("C_mixed", EPOCHS, []),
]


def logline(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
    with open(LOG, "ab", buffering=0) as f:
        f.write(line.encode("utf-8", errors="replace"))


def run(variant, epochs, extra):
    cmd = [PY, "-u", SCRIPT, "--variant", variant, "--epochs", str(epochs),
           "--threads", "8", *extra]
    logline(f"===== START {variant} (epochs={epochs}) =====")
    logline("CMD " + " ".join(cmd))
    t0 = time.time()
    # CREATE_NO_WINDOW：子进程不挂控制台。上次 A/B 因 forrtl "window-CLOSE
    # event"（控制台被关闭/发 CTRL_CLOSE）被杀，无控制台即免疫控制台事件。
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    with open(LOG, "ab", buffering=0) as logf:
        rc = subprocess.call(cmd, cwd=PROJECT, stdout=logf,
                             stderr=subprocess.STDOUT,
                             creationflags=flags)
    logline(f"===== END {variant} rc={rc} "
            f"({(time.time()-t0)/60:.1f} min) =====")
    return rc


def main():
    logline(f"########## run_ablation_local START epochs={EPOCHS} "
            f"host={os.environ.get('COMPUTERNAME', '?')} ##########")
    failures = []
    for v, e, x in STEPS:
        try:
            rc = run(v, e, x)
        except Exception as exc:  # 驱动自身异常也要继续
            logline(f"!!!!! {v} driver exception: {exc!r}")
            rc = -1
        if rc != 0:
            failures.append((v, rc))
            logline(f"!!!!! {v} FAILED rc={rc}, continuing")
    logline(f"########## run_ablation_local DONE failures={failures} ##########")


if __name__ == "__main__":
    main()
