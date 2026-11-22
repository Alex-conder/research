#!/bin/bash
# RPE1 全量推理（独立容器 64 核配额；线程/BLAS 封顶 64）
setsid nohup docker exec \
  -e SCT_THREADS=64 -e OMP_NUM_THREADS=64 -e MKL_NUM_THREADS=64 \
  -e OPENBLAS_NUM_THREADS=64 \
  research_scgenept_rpe1 \
  python /research/scripts/run_scgenept_k562.py \
  --h5ad /research/data/perturb_processed_rpe1.h5ad \
  --local-npz /research/meta_rpe1/pred_expr_seed0.npz \
  --pool 32 --batch 8 \
  --out /research/results/scgenept_rpe1 \
  < /dev/null > /opt/research/results/rpe1.log 2>&1 &
echo "launched pid $!"
