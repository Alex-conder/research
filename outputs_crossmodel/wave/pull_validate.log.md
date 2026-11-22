# 多模型波产物 传输 + 校验日志

> 日期：2026-10-24。源：`ah-server:/opt/research/outputs_crossmodel/`（read-only 拉取）。
> 目标：本机 `virtual_cell_constraint_mvp/outputs_crossmodel/wave/`。
> 共用机纪律遵守：纯网络传输，零计算占用；动物房四容器未触碰。

## 1. 传输过程

| 阶段 | 方式 | 结果 |
|---|---|---|
| 1 | `scp ah-server:...` 顺序拉 18 文件 | ~425s 完成 scgpt 6 文件（66.3MB），实测 ~0.16MB/s；判定路径瓶颈后主动停止 |
| 2 | `pull_wave.py`（paramiko SFTP，走 ah-server 别名同密钥 `E:\animal_house_build\.ssh\animal_house_deploy2`，172.16.2.7:22022） | 前 ~310s 完成 geneformer k562 三件套 + rpe1 部分；因后台任务时限切换 |
| 3 | `pull_wave.py` 断点续跑（按文件大小判定，`-u` 无缓冲） | 全部完成 |

最终一轮输出（摘）：

```
GET  geneformer/geneformer_wave_rpe1_predictions.npz 57.0MB 468.9s 0.12MB/s OK
GET  uce/uce_wave_k562_predictions.npz                 68.0MB 346.1s 0.20MB/s OK
GET  uce/uce_wave_rpe1_predictions.npz                 95.2MB 492.8s 0.19MB/s OK
TOTAL 326.4MB in 1309.2s (0.25MB/s incl. skips)
```

- 总量：18 文件 / 326.4MB（6 npz 311.2MB + 12 小文件）。
- 总墙钟：约 34 分钟（含两次方式切换的损耗）；瓶颈为 ssh 转发路径吞吐
  （实测 0.12-0.20MB/s，与 RUNBOOK 记录的直连 4-8MB/s 不同——本次经 172.16.2.7:22022 转发）。
- 逐个文件大小与服务器 `stat` 一致（全部 OK，无 SIZE_MISMATCH）。
- 传输脚本：`pull_wave.py`（可重跑，幂等续传）。

## 2. md5 校验（18/18 PASS）

本机 `md5sum` 与服务器 `md5sum` 逐一比对（`_local_md5.txt` / `_server_md5.txt`，
比对器 `_md5_check.py`）。全部一致：

```
OK  geneformer/geneformer_wave_k562_labels.csv
OK  geneformer/geneformer_wave_k562_meta.json
OK  geneformer/geneformer_wave_k562_predictions.npz
OK  geneformer/geneformer_wave_rpe1_labels.csv
OK  geneformer/geneformer_wave_rpe1_meta.json
OK  geneformer/geneformer_wave_rpe1_predictions.npz
OK  scgpt/scgpt_wave_k562_labels.csv
OK  scgpt/scgpt_wave_k562_meta.json
OK  scgpt/scgpt_wave_k562_predictions.npz
OK  scgpt/scgpt_wave_rpe1_labels.csv
OK  scgpt/scgpt_wave_rpe1_meta.json
OK  scgpt/scgpt_wave_rpe1_predictions.npz
OK  uce/uce_wave_k562_labels.csv
OK  uce/uce_wave_k562_meta.json
OK  uce/uce_wave_k562_predictions.npz
OK  uce/uce_wave_rpe1_labels.csv
OK  uce/uce_wave_rpe1_meta.json
OK  uce/uce_wave_rpe1_predictions.npz
MD5_OVERALL: PASS
```

过程注记：首轮手工誊录服务器 md5 时有一处笔误（geneformer k562 meta.json），
以服务器实时 md5sum 重取为准，重取后本地=服务器=e1bccfabe2acf6954cdee653c796a42c，
确认是记录笔误而非传输错误。

## 3. 内容校验（6/6 产物 PASS，各 27/27 项）

校验器：`validate_wave.py`（只读；`python validate_wave.py wave/`）。
每份产物检查 27 项，六份全部通过。检查类别：

- labels：行数（K562=217 / RPE1=306）、列齐全（base + 引擎 flag + 4 探针列）、
  pearson 两列全空（=NaN，设计如此）、探针列全为有限数。
- npz：必需键齐全；各嵌入数组维度 = 引擎维度（scGPT 512 / Geneformer 768 / UCE 1280）；
  pred 侧与 true_emb_mean / ctrl_emb 全有限；ctrl 池=32；true_n_cells 全正；
  output_type 与 meta 均标注 embedding。
- meta.json：output_type 为 embedding、h5ad 路径与数据集一致。

| 产物 | 行数 | 嵌入维度 | 校验 | 备注 |
|---|---|---|---|---|
| scgpt/k562 | 217 | 512 | 27/27 PASS | 4 行条件细胞数<32（NaN 填充槽） |
| scgpt/rpe1 | 306 | 512 | 27/27 PASS | 35 行有填充槽 |
| geneformer/k562 | 217 | 768 | 27/27 PASS | 4 行有填充槽 |
| geneformer/rpe1 | 306 | 768 | 27/27 PASS | 35 行有填充槽 |
| uce/k562 | 217 | 1280 | 27/27 PASS | 4 行有填充槽 |
| uce/rpe1 | 306 | 1280 | 27/27 PASS | 35 行有填充槽 |

## 4. 校验中发现并确认的两个设计口径（非缺陷）

1. **pearson 列为空串而非 "nan" 字面量**：labels.csv 的 pearson_delta / pearson_delta_all
   对全部 217/306 行是空单元格。按 2026-10-16 口径决定（embedding 类误差用模型空间
   距离，不与表达谱类混比），空 = 无值，属设计如此。消费脚本按空/NaN 两种写法都要能读。
2. **true_emb_cells 的整细胞 NaN 是短条件的填充槽**：6 份产物中受影响的行完全一致
   （K562 4 行：SF1/RPL39/WDR70/PRPF31；RPE1 35 行），逐行满足
   `NaN 槽数 == 32 - true_n_cells`，无半细胞（部分维度）NaN，labels 的 n_cells 列与
   true_n_cells 完全一致。已排除数据缺陷假设（服务器端复核：这 4 个条件的细胞在 h5ad
   内表达总量正常，NaN 与零表达细胞无关，就是细胞数不足 32 的 padding）。
   消费方务必按 `np.isfinite` 掩码取用 true_emb_cells。

## 5. 附加统计（供消费方参考）

| 产物 | target_in_input_space=True | 扰动生效 flag=True |
|---|---|---|
| scgpt/k562 | 79/217 | 79 (target_masked) |
| scgpt/rpe1 | 142/306 | 142 |
| geneformer/k562 | 83/217 | 79 (target_deleted；另 4 个有基因无 Ensembl/token 映射) |
| geneformer/rpe1 | 143/306 | 142 |
| uce/k562 | 79/217 | 79 (target_deleted) |
| uce/rpe1 | 142/306 | 142 |

## 6. 产物文件（本机最终路径）

- `wave/scgpt/scgpt_wave_k562_{predictions.npz,labels.csv,meta.json}`
- `wave/scgpt/scgpt_wave_rpe1_{predictions.npz,labels.csv,meta.json}`
- `wave/geneformer/geneformer_wave_k562_{predictions.npz,labels.csv,meta.json}`
- `wave/geneformer/geneformer_wave_rpe1_{predictions.npz,labels.csv,meta.json}`
- `wave/uce/uce_wave_k562_{predictions.npz,labels.csv,meta.json}`
- `wave/uce/uce_wave_rpe1_{predictions.npz,labels.csv,meta.json}`

辅助文件（可复跑/复验）：`pull_wave.py`、`validate_wave.py`、`_md5_check.py`、
`_local_md5.txt`、`_server_md5.txt`。消费说明见同目录 `READY.md`。
