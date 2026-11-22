# RUNBOOK — 隔离研究环境 & scGenePT 跨模型泛化验证（Replogle K562）

> 服务器：172.16.92.2（Rocky Linux 9.8，256 线程 EPYC 7702 / 881GB RAM / 无 GPU）
> 工作区：`/opt/research`（与动物房系统完全隔离；自建容器 restart=no，显式资源配额）
> 硬约束：**绝不** stop/restart/remove `animal_house_*` 容器；动手前后看 `cat /proc/loadavg`。

## 0. 目录结构

```
/opt/research/
  data/perturb_processed_k562.h5ad   # Replogle K562 essential（GEARS 处理版, 2.9GB）
  data/perturb_processed_rpe1.h5ad   # Replogle RPE1（2.85GB）
  meta/          # K562 held-out 名单与本地 step6 口径文件（pred_expr_seed0.npz / predictions_seed0.csv / split_seed0.pkl）
  meta_rpe1/     # RPE1 同三件套（306 扰动）
  scgenept_repo/ # czi-ai/scGenePT @ 5847c57e4cccd9949a7f0b5f3775d5dcae5a5eac
  scfoundation_repo/  # biomap-research/scFoundation model/ 代码（测速用）
  models/        # scGenePT 微调权重(954MB) + NCBI+UniProt 嵌入(416MB) + scfoundation-gene model.pt(477MB)
  results/       # 推理产物 + 日志
  scripts/       # setup_env.sh / run_scgenept_k562.py / launch_full.sh / launch_rpe1.sh / scf_speedtest.py
```

## 1. 环境搭建（一次性）

```bash
bash /opt/research/scripts/setup_env.sh > /opt/research/scripts/setup.log 2>&1
```

要点（踩过的坑）：
- 基础镜像 `python:3.10-slim-bookworm`：**必须 3.10**——scgpt 0.2.1 要求 <3.11（3.11 只能用 0.2.2+，未验证）。
- torch 先装 CPU 版：`pip install torch==2.1.2+cpu torchtext==0.16.2+cpu --index-url https://download.pytorch.org/whl/cpu`，
  否则 pip 会拉 2.5GB CUDA 版。
- **numpy 必须钉死 1.26.4**（setup 脚本最后一条 pip 已做）：scgpt 的依赖链会把 numpy 顶到 2.x，
  torch 2.1.2 是 numpy1 ABI，混用直接炸。
- **torch 必须钉死 2.1.2+cpu**：一次依赖解析事故把 torch 顶到 2.14.1+cu130（全套 nvidia CUDA 包），
  torchtext 0.16.2 的 libtorchtext.so 立刻链接失败（libc10 not found）。修复：
  `pip install --no-cache-dir torch==2.1.2+cpu torchtext==0.16.2+cpu --index-url https://download.pytorch.org/whl/cpu`
  再 `pip uninstall` 13 个 nvidia-* 包。修复后已 commit 为镜像 `research:scgenept-py310-v2`（10.7GB）。
- 容器：`docker run -d --name research_scgenept --restart no --cpus 32 --memory 128g --pids-limit 512 -v /opt/research:/research -w /research python:3.10-slim-bookworm sleep infinity`
- 服务器访问 github.com 不稳定（HTTP/2 断流/超时）：repo 用本地下载 tarball 后 SFTP 上传，
  并把 commit 写进 `scgenept_repo/repo/GIT_COMMIT`。
- 服务器直连 S3 仅 ~200KB/s：scGenePT 权重(954MB)/嵌入(416MB)在本机下载后 SFTP 上传更快。
  本机与服务器之间 SFTP ~4-8MB/s。
- 宿主机跑 docker exec 时路径是**容器内**路径（/research/...），不是 /opt/research/...。

## 2. 数据上传（本机 → 服务器，paramiko SFTP）

| 文件 | 远端 | 说明 |
|---|---|---|
| `data/gears_demo/replogle_k562_essential/perturb_processed.h5ad` (2.9GB) | `data/perturb_processed_k562.h5ad` | 唯一需要的组学数据；ctrl/真值均值都从它算 |
| `outputs/pred_expr_seed0.npz` (3.9MB) | `meta/` | held-out 名单（217 扰动）+ 本地基因序 |
| `outputs/predictions_seed0.csv` | `meta/` | 本地 GEARS 标签 + n_cells（对齐校验用） |
| splits/...simulation_0_0.75_528ce902c6.pkl | `meta/split_seed0.pkl` | seed0 simulation split（校验用） |

单文件均 <5GB，未压缩。h5ad 上传后与本地字节数一致（2,898,717,719）。

## 3. scGenePT 推理

```bash
# K562 全量（217 扰动 × 32 ctrl 细胞池）：
bash /opt/research/scripts/launch_full.sh
# RPE1 全量（306 扰动，第二容器，NUMA socket1）：
bash /opt/research/scripts/launch_rpe1.sh
# 即：docker exec -e SCT_THREADS=64 -e OMP/MKL/OPENBLAS_NUM_THREADS=64 \
#       research_scgenept[_rpe1] python /research/scripts/run_scgenept_k562.py \
#       --h5ad <rpe1 时指定> --local-npz <rpe1 时指定> --pool 32 --batch 8 \
#       --out /research/results/scgenept_k562|scgenept_rpe1
```

- **断点续跑**：每扰动写 `<out>_checkpoint.npz`（指纹=基因空间+扰动名单+pool+batch），
  重启自动恢复已完成部分（ctrl 池同指纹可复用）。
- **NUMA 钉扎**（8 NUMA node / 2 socket，不钉会跨 socket 抢带宽）：
  `docker update --cpuset-cpus "0-63,128-191" research_scgenept`（socket0，配额 96 核）
  `docker update --cpuset-cpus "64-127,192-255" research_scgenept_rpe1`（socket1，配额 64 核）
  实测 RPE1 从 ~168s/扰动降到 ~112s/扰动；K562 因每扰动细胞量更大，稳态 ~171s/扰动。
- **杀任务**：`docker top <容器> | awk '/run_scgenept_k562/ {print $2}' | xargs -r kill`。
  千万别用 `pkill -f`——会匹配到自身命令行把 SSH 会话杀死。
- **docker commit 会 PAUSE 容器**（冻结推理进程），绝不在任务运行时 commit。
- 长跑一律 `setsid nohup docker exec ... < /dev/null > 日志 2>&1 &`，SSH 断连不杀进程。

评估口径（与本地 step1_run.py/step6_run.py 完全对齐）：
- held-out 名单 = `meta/pred_expr_seed0.npz` 的 pred_perturbations（217，split_seed=0 simulation test）。
- `pearson_delta` = top-50 |true_delta| 基因上 Pearson(pred_delta, true_delta)；`pearson_delta_all` = 匹配基因全集。
- true_mean = 该扰动 condition 细胞均值；ctrl_mean = ctrl 细胞均值（同一 h5ad，log-normalized 空间）。
- **跨数据集适配**（协议偏离，已记录）：模型在 Norman 5045 基因空间微调，此处经 scGPT whole-human vocab
  匹配后输入 K562 的 4813/5000 基因（187 个不在 vocab）。
- **数据事实**：217 个 held-out 中 138 个的目标基因不在 h5ad 的 5000 输入基因里（GEARS 图含扰动节点，
  token 模型无法标记不存在的 token）→ 这些扰动模型实际收到全零 pert flag，预测≈ctrl；
  逐扰动记录 `target_in_input_space`，汇总同时报全集与子集。

产物：`results/scgenept_k562_predictions.npz`（matched 4813 空间 + full 5000 NaN 填充版 +
true_expr + ctrl_mean + 基因序）、`_labels.csv`、`_meta.json`（模型/数据 sha256、命令行、耗时）。

## 4. scFoundation CPU 测速

```bash
docker exec -e SCF_THREADS=8 -e OMP_NUM_THREADS=8 -e OPENBLAS_NUM_THREADS=8 \
  -e MKL_NUM_THREADS=8 research_scgenept python /research/scripts/scf_speedtest.py
```
权重：`hf-mirror.com/perturblab/scfoundation-gene/resolve/main/model.pt`（477MB，Apache-2.0，
源自 biomap-research/scFoundation Nature Methods 2024；另有 genbio-ai/scFoundation models.ckpt 1.8GB）。
注意：与主任务同容器并行时**必须**用环境变量把 BLAS 线程压到 8 以下，否则 numpy C 扩展
在 pids-limit 下导入即崩（`PyCapsule_Import ... datetime` 假象报错）。

## 5. 收尾

```bash
docker stop research_scgenept research_scgenept_rpe1   # 不用即停；restart=no，宿主机重启不会自启
docker start research_scgenept                          # 复用环境
cat /proc/loadavg                      # 确认负载回落
```
动物房容器全程未触碰（`docker ps` 快照见汇报）。

## 6. 关键数字（2026-10-15 实测）

- scGenePT K562（64 线程, socket0 钉扎, pool32×batch8, seq=4813）：**~171s / 扰动**，
  217 扰动 ≈ 10.3h；RPE1（64 线程, socket1 钉扎）：**~112s / 扰动**，306 扰动 ≈ 9.5h；
  两容器并行总吞吐 ~2 扰动/143s。旧基准 32 线程不钉扎单任务 ~86-94s/扰动。
- scFoundation（xTrimoGene MAE 100M）：8 线程（与主任务并行，保守值）单细胞 full forward
  **≈17.2s**（encoder 4627 token 10.7s + performer decoder 全 19264 基因）；权重加载 0 missing/0 unexpected。
  换算：32 线程满负荷下单扰动（池 32 细胞）估 ~2.5-4 min；217 扰动全量 CPU 约 10-15h → 批量可行但接近上限。
- 冒烟（PRPF3, pool4）scGenePT pearson_delta(top50)=0.398；基准 3 扰动均值 0.02
  （跨数据集零样本信号弱，正式结果以全量为准）。
- 负载：双任务并行推理期间 load ~65-78（128 线程跑在 256 线程机器，内存充足）。
