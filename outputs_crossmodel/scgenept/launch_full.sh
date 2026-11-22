#!/bin/bash
# K562 全量推理（64 线程；BLAS 线程同步封顶防 pids-limit 击穿）
setsid nohup docker exec \
  -e SCT_THREADS=64 -e OMP_NUM_THREADS=64 -e MKL_NUM_THREADS=64 \
  -e OPENBLAS_NUM_THREADS=64 \
  research_scgenept \
  python /research/scripts/run_scgenept_k562.py --pool 32 --batch 8 \
  --out /research/results/scgenept_k562 \
  < /dev/null > /opt/research/results/full.log 2>&1 &
echo "launched pid $!"
