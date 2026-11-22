"""轮询服务器直到 K562/RPE1 推理双双结束，最长 20h。退出码 0 = 双双完成。"""
import sys, time
sys.path.insert(0, ".")
from _srv import connect, exec_cmd

c = connect()
deadline = time.time() + 20 * 3600
while time.time() < deadline:
    rc, out, _ = exec_cmd(c,
        "echo files=$(grep -l '\\[out\\]' /opt/research/results/full.log /opt/research/results/rpe1.log 2>/dev/null | wc -l);"
        "echo kproc=$(docker top research_scgenept 2>/dev/null | grep -c run_scgenept_k562);"
        "echo rproc=$(docker top research_scgenept_rpe1 2>/dev/null | grep -c run_scgenept_k562);"
        "grep '\\[inf\\]' /opt/research/results/full.log | tail -1;"
        "grep '\\[inf\\]' /opt/research/results/rpe1.log | tail -1;"
        "cat /proc/loadavg", 60)
    print(time.strftime("%H:%M:%S"), out.strip().replace("\n", " | "), flush=True)
    if "files=2" in out or ("kproc=0" in out and "rproc=0" in out):
        rc, out2, _ = exec_cmd(c, "ls -la /opt/research/results/ | grep -E 'scgenept_(k562|rpe1)'", 60)
        print(out2, flush=True)
        sys.exit(0 if "files=2" in out else 3)
    time.sleep(600)
sys.exit(2)
