# Redteam null tests (step23, 2026-10-30)

结果级红黑队：四个头条新数字的置换零检验（B=1000）。先与存档对账
（tol 1e-6）再跑置换；判读：p < 0.05 -> 数字非零假设产物，可回填；
p >= 0.05 -> 加 null-flag 如实报。

## norman 探针均值 0.7833

- 统计量：3 种子 cv_auroc_mean 的跨种子均值（锁定判据口径）
- seed3: obs cv_mean=0.6000, pooled OOF=0.5500, 置换 p(pooled)=0.4026
- seed4: obs cv_mean=0.8500, pooled OOF=0.7100, 置换 p(pooled)=0.1189
- seed5: obs cv_mean=0.9000, pooled OOF=0.9000, 置换 p(pooled)=0.0040
- obs mean=0.7833; null mean=0.5099 sd=0.1049
- **置换 p = 0.0040** (B=1000)

## adamson 探针均值 0.7000

- 统计量：3 种子 cv_auroc_mean 的跨种子均值（锁定判据口径）
- seed3: obs cv_mean=0.8000, pooled OOF=0.6406, 置换 p(pooled)=0.2178
- seed4: obs cv_mean=0.6000, pooled OOF=0.6250, 置换 p(pooled)=0.2298
- seed5: obs cv_mean=0.7000, pooled OOF=0.6406, 置换 p(pooled)=0.1708
- obs mean=0.7000; null mean=0.5083 sd=0.1224
- **置换 p = 0.0619** (B=1000)

## G-C RPE1 s0 B=50 probe_full 效用

- obs 点差 cost_random - cost_gate = 0.9387（存档 bootstrap 均值 0.9413 [0.6995,1.1773]）
- null mean=-0.0032 sd=0.1437
- **置换 p = 0.0010** (B=1000)
- 次级 B=10: obs=0.2268, 存档 boot=0.2439, p=0.0010
- 次级 B=20: obs=0.5028, 存档 boot=0.4926, p=0.0010
- 次级 B=5: obs=0.1428, 存档 boot=0.1387, p=0.0010

## G-D A_k562 conformal 单元覆盖 0.8921

- N=560000 单元，覆盖 499563，名义 0.9
- MC 标签打乱（名义率）B=1000: 双侧 p=0.0010，单侧欠覆盖 p=0.0010
- 精确二项：双侧 p=5.093e-85，单侧欠覆盖 p=2.609e-85
- 判读：覆盖可分辨于名义率零假设（p<0.05，方向为系统性欠覆盖 0.8pp）；锁定接收带 [0.85, 0.95] 内（Wilson95 存档 [0.8913,0.8929]），PASS 判据为带成员资格 -> 无 null-flag，欠覆盖幅度如实写。

