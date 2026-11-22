# RUNBOOK — scGPT whole-human 跨模型泛化验证（Replogle K562）

> 服务器：172.16.92.2（Rocky Linux 9.8，256 线程 EPYC 7702 / 881GB RAM / 无 GPU）
> 工作区：`/opt/research`（与动物房系统完全隔离；自建容器 restart=no）
> 硬约束：**绝不** stop/restart/remove `animal_house_*` 容器。

## 0. 产物约定（与 scgenept 同构）

`/research/results/scgpt_k562_predictions.npz`：
`pred_expr (217, n_matched)` / `pred_expr_full (217, 5000)` / `pred_perturbations (217,)` /
`pred_genes (n_matched,)` / `pred_gene_indices_in_h5ad` / `target_in_input_space (217,)` /
`true_expr` / `true_expr_full` / `ctrl_mean_matched` / `ctrl_mean_full` / `full_gene_space (5000,)` /
外加 scGPT 特有：`pred_pert_gene_value (217,)`（扰动基因处预测值，NaN=不在输入空间）、
`output_type`（标注输出类型）。
另有 `_labels.csv` / `_meta.json` / `_genes_matched.csv` / `_genes_full.csv`。

## 1. 容器与环境

```bash
docker run -d --name research_scgpt --restart no --cpus 64 --memory 128g \
  --pids-limit 512 --cpuset-cpus "0-63,128-191" \
  -v /opt/research:/research -w /research \
  research:scgenept-py310-v2 sleep infinity
```

- 复用 `research:scgenept-py310-v2` 镜像（scGenePT 修复后的 env：torch 2.1.2+cpu /
  scgpt 0.2.1 / anndata 0.10.9 / numpy 1.26.4）。**不要**在这个容器里 pip 升级
  torch/numpy（ABI 事故史见 scgenept/RUNBOOK.md）。
- cpuset 与 research_scgenept（K562 任务）共用 socket0：等 scgenept 任务结束后
  独占带宽；如需 NUMA 调整：`docker update --cpuset-cpus "..." research_scgpt`。

## 2. 权重（hf-mirror 直下载，无需 S3/GDrive）

```bash
mkdir -p /opt/research/models/scgpt_whole_human && cd /opt/research/models/scgpt_whole_human
curl -sL --retry 5 -o best_model.pt https://hf-mirror.com/perturblab/scgpt-human/resolve/main/best_model.pt
curl -sL -o vocab.json https://hf-mirror.com/perturblab/scgpt-human/resolve/main/vocab.json
curl -sL -o args.json  https://hf-mirror.com/perturblab/scgpt-human/resolve/main/args.json
```

- `perturblab/scgpt-human` = bowang-lab/scGPT 官方 whole-human `best_model.pt`
  的 HF 镜像（ntoken=60697, 12 层, embsize=512, ~51.3M params, 205MB）。
- 下载实测 ~8.4MB/s，40s 内完成。
- **坑**：checkpoint 是 fast_transformer 训练保存的，attention 键名为
  `*.attn.Wqkv.*`；CPU 上 `use_fast_transformer=False` 的 MultiheadAttention
  需要 `in_proj_*`。脚本里已做 `".attn.Wqkv." → ".attn.in_proj_"` 重命名 + strict=False。
- checkpoint 里还有 `flag_encoder.weight` 等多余键（perturblab 内部框架用），
  load 时报 unexpected 属正常；`cls_decoder` 等头缺失不影响（CLS=False）。
- vocab.json 与 scgenept_repo 内的 scGPT vocab **逐字节一致**（60697 tokens）。

## 3. 推理协议（canonical scGPT 零样本 mask 协议）

输入 = ctrl 池（pool=32, seed=0）log-normalized 表达（连续值，与 scGenePT 同输入处理）；
扰动基因值 mask 为 -1（mask_value，pad_value=-2 不冲突）；MLM 头 `mlm_output`
读出全 token 预测；`pred_expr = 池均值`。未扰动位置为恒等重建，扰动位置为真实预测。
局限（写入 meta）：mask 协议下预测谱整体≈ctrl，pearson_delta 主要反映扰动基因单点，
`pred_pert_gene_value` 单独留痕供下游用。

## 4. 运行

```bash
# 冒烟（3 扰动）：
docker exec -e SCT_THREADS=32 research_scgpt python /research/scripts/run_scgpt_k562.py \
  --limit 3 --threads 32 --out /research/results/scgpt_smoke
# 全量（217 扰动，setsid nohup 后台，断 SSH 不杀）：
setsid nohup docker exec -e SCT_THREADS=64 -e OMP_NUM_THREADS=64 \
  -e MKL_NUM_THREADS=64 -e OPENBLAS_NUM_THREADS=64 \
  research_scgpt python /research/scripts/run_scgpt_k562.py \
  --threads 64 --out /research/results/scgpt_k562 \
  < /dev/null > /opt/research/results/scgpt_k562.log 2>&1 &
```

- 断点续跑：每扰动写 `scgpt_k562_checkpoint.npz`（指纹=基因空间+扰动名单+pool+batch），
  重启自动跳过已完成。
- 杀任务：`docker top research_scgpt | awk '/run_scgpt_k562/ {print $2}' | xargs -r kill`
  （**别用 pkill -f**，会匹配自身命令行把 SSH 会话杀死）。
- 数据事实（与 scgenept 相同）：217 个 held-out 中 138 个目标基因不在 h5ad 输入空间，
  无法 mask，预测≈ctrl，逐扰动留痕 `target_in_input_space`。

## 5. 日志/监控

```bash
tail -f /opt/research/results/scgpt_k562.log
cat /proc/loadavg   # 机器 256 线程；scgenept 双任务并行期间 load 80 左右正常
docker stats research_scgpt --no-stream
```

（耗时实测待补——见 _meta.json 与下方"关键数字"）

## 6. 关键数字（2026-10-16 实测）

- 冒烟（3 扰动, 32 线程, 与 scgenept 双任务争抢 socket0）：~222s/扰动；
  pearson_delta(top50) 3 扰动均值 **0.317**（目标在输入空间子集 n=2 均值 0.325，
  对照 scGenePT 冒烟 PRPF3 0.398，同量级零样本信号）。
- 全量已启动：`setsid nohup docker exec -e SCT_THREADS=32 ... research_scgpt
  python -u /research/scripts/run_scgpt_k562.py --threads 32
  --out /research/results/scgpt_k562`（日志 `results/scgpt_k562.log`）。
  高负载下 ETA ~13h；scgenept 任务结束后会显著加速。32 线程比 64 线程单扰动
  更快（scGenePT 同款现象：64 线程 171s/扰动 vs 32 线程 86-94s/扰动）。
- **坑（已修）**：scgpt 0.2.1 用 `nn.TransformerEncoder`，权重键为
  `transformer_encoder.layers.N.self_attn.Wqkv.*`（不是 scGenePT repo 的
  `encoder.attn.Wqkv.*`），重命名必须匹配 `.self_attn.Wqkv.` 前缀。

## 7. 收尾

```bash
docker stop research_scgpt   # 不用即停；restart=no
docker start research_scgpt  # 复用环境
```
