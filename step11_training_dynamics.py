# step11_training_dynamics.py — 入口 A（训练动力学）：误差签名是"早期定形"还是"后期涌现"？
#
# ============================ 探索性预注册（2026-10-16，出数据前写死） ============================
#
# 【性质】协议外探索（exploratory），不进主判据推断；single-model single-seed（n=1）
# 的限制在所有结论中如实标注。
#
# 【科学问题】GEARS 的"误差签名"（就地 9 特征 logistic 探针的标准化系数向量）在
# 训练过程中何时成形——形态（系数向量）与功能（探针 AUROC / 预测质量 pearson_delta）
# 谁先稳定？归档模型 outputs/gears_model_seed0/model.pt 只有最终 checkpoint
# （3 epochs，4346s，见 day2_run5.log），本实验补出训练轨迹。
#
# 【阶段 1（管线验证）+ 阶段 2（主实验）合并为一次 12-epoch 训练的说明】
#   原设计为两次独立训练（3 epochs 复现 + 12 epochs 轨迹）。二者数学上等价：
#   StepLR(step_size=1, gamma=0.5) 的 lr 衰减只依赖"已训练的 epoch 数"，不依赖
#   总轮数；数据顺序由同一随机流决定；optimizer/scheduler 状态逐 epoch 演化相同。
#   因此 12-epoch 训练的 e3 快照 == 独立 3-epoch 训练的最终模型（模 CUDA 原子操作
#   非确定性）。为节省 ~75min GPU 时间，合并为一次训练；复现判据（gate）在 e3
#   快照分析后立即判定并留痕，不达标仍继续跑完（不销毁已耗算力），verdict 中
#    prominently 标注 gate 结果。
#   —— 此为执行层面的合并，不改变任何统计口径；特此声明。
#
# 【训练配置（必须与归档对齐，逐项核对）】
#   数据源/划分：step1_run.py 原逻辑逐行复用——data/gears_demo/replogle_k562_essential，
#     SPLIT_SEED=0, HOLDOUT_FRAC=0.2, np.random.default_rng(0) 选 217 个留出单扰动基因；
#     prepare_split(split="simulation", seed=0, test_perts=..., only_test_set_perts=True)
#     （split 文件已本地缓存，"Local copy of split is detected"，加载即归档同款划分）；
#     batch=32 / test_batch=128；hidden_size=64；模型初始化 torch.manual_seed(0) +
#     np.random.seed(0)（注意：gears 包 import 时有模块级 torch.manual_seed(0)，
#     MODEL_SEED=0 时与归档完全一致）。
#   训练循环：镜像 venv gears/gears.py GEARS.train()（不 fork 包源码、不改包文件；
#     调包内 loss_fct/evaluate/compute_metrics 同款入口）：
#     Adam(lr=1e-3, weight_decay=5e-4) + StepLR(step_size=1, gamma=0.5) +
#     clip_grad_value_(1.0)；train 侧评估用 venv 兼容补丁的 evaluate(..., max_batches=500)。
#     自写驱动的唯一动机：train() 每次调用重建 optimizer/scheduler，无法满足
#     --resume-snap 的"optimizer/scheduler 状态连续"要求。
#
# 【快照】每 epoch 结束保存 outputs/training_dynamics/snap_e{k}.pt：
#   model/optimizer/scheduler state_dict + torch CPU/CUDA RNG state + np RNG state +
#   epoch + min_val，保证 --resume-snap 续跑与不间断训练尽量一致（残留差异 =
#   CUDA 非确定性，如实标注）。epoch 0（model_initialize 后、训练前）为基线分析点，
#   不落权重快照。
#
# 【每快照分析（特征/标签口径全部复用现有脚本）】
#   预测：对 217 个留出扰动中在扰动图内的基因（in_graph，与 step1 相同过滤）
#     用该 epoch 末权重前向（ gears.predict 同款路径；分析前重置 saved_pred 缓存，
#     并将 best_model 指到当前 epoch 权重——本实验关心 epoch-k 权重本身，不是
#     best-checkpoint 语义）。
#   pearson_delta：step1_run.py:152-166 口径——top-50 |真实 delta| 基因上
#     Pearson(pred_delta, true_delta)；对每个快照的预测现算（标签随快照变）。
#   9 特征：step6_gating_probe.compute_features（8 谱统计）+ log_n_cells
#     （n_cells 取归档 outputs/predictions_seed0.csv，实验设计量，与快照无关）。
#   探针：step7_routing.fit_probe 同规格——StandardScaler + LogisticRegression(
#     C=1.0, max_iter=2000)；y = (1 - pearson_delta) > 训练中位数；
#     扰动级 70/30 划分、split_seeds=[0,1,2]、主 seed 20261016（与 step7 完全同
#     的 Generator 派生代码，n=217 复现 step7 的 perm）。
#   记录：3 个 split 种子各自的系数向量（标准化系数，即 scaler+LR pipeline 的
#     coef_）与 65 测试侧 AUROC（标签中位数取训练侧，与 step7 probe_test_auroc
#     同口径）；全 217 的 pearson_delta 均值；全量 217 就地拟合系数（= 归档
#     error_signatures.json GEARS(K562) 的口径，step8b 已验证 max|dev|<1e-6）。
#
# 【阶段 1 判据（gate，e3 时判定）】
#   e3 的全量就地拟合系数 与 outputs_crossmodel/analysis/error_signatures.json
#   的 "GEARS(K562)" 算余弦（未归一的原始系数直接算 cos；另报标准化后 cos）。
#   主判据：全量拟合 cos > 0.8（与归档完全同口径）；辅助：152 训练侧拟合
#   （3 split 种子均值）对归档签名的 cos。不达标 → verdict 中标注"先查口径"，
#   不解读后续轨迹形态。
#
# 【假设（先写死）】
#   T1（形态时点）：签名对最终签名（epoch 12，3 split 种子均值向量）的余弦
#     trajectory crossing 0.8 的最早 epoch k₁。
#   T2（功能时点）：探针 AUROC（3 split 种子均值）进入其终值 ±0.02 带内的最早
#     epoch k₂。
#   T3：pearson_delta 均值进入其终值 ±0.01 带内的最早 epoch k₃。
#   判定：k₁ < k₂ 且 k₁ < k₃ → "签名早期定形、功能后至"（支持"误差结构是训练
#     早期由数据/架构决定的功能性属性"）；否则 → "后期涌现"。
#   epoch 0 基线点参与全部 trajectory。
#
# 【判据纪律 / caveat（已知风险，写死）】
#   - 系数尺度不稳（L2 正则 + logistic + 共线特征：n_de_025 与 frac_de_025
#     近乎线性相关）：余弦一律用"系数除以自身 L2 范数"的标准化向量计算，
#     verdict 中写明此 caveat。
#   - 任何 early/late 结论标注 single-model、single-seed、单一数据集限制。
#   - 附带记录每 epoch Train/Validation Overall MSE（gears evaluate 口径）。
#
# 【降级条款】若 12 epochs 因显存/时间不可行，降级为 6 epochs 并在 verdict
#   如实说明（--epochs 参数即此开关）。
#
# 【续跑纪律】每 epoch 快照与分析记录立即落盘（per_epoch.json 全量重写，原子写）；
#   时间耗尽时报告进度与续跑命令：
#     PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step11_training_dynamics.py \
#         --resume-snap --epochs 12
#   绝不中途杀训练进程。
#
# 【产物】outputs/training_dynamics/snap_e{1..12}.pt + per_epoch.json；
#   outputs/training_dynamics.json + training_dynamics.png（300dpi 英文双面板）；
#   outputs/training_dynamics_verdict.md（跑完后另行撰写）。
#
# 【随机源】config seed：SPLIT_SEED=0（划分）、MODEL_SEED=0（初始化）；
#   快照分析专用 seed：MASTER_SEED=20261016（探针 70/30 划分，与 step7 同源）。
#
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step11_training_dynamics.py --epochs 12
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step11_training_dynamics.py --resume-snap --epochs 12
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step11_training_dynamics.py --plot-only   # 从 per_epoch.json 重出最终 json/png
# =================================================================================================
import argparse
import json
import os
import time

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(REPO, "data", "gears_demo")
DATA = "replogle_k562_essential"
OUT = os.path.join(REPO, "outputs", "training_dynamics")
PRED_CSV = os.path.join(REPO, "outputs", "predictions_seed0.csv")
ARCHIVED_SIG_JSON = os.path.join(REPO, "outputs_crossmodel", "analysis",
                                 "error_signatures.json")

HOLDOUT_FRAC = 0.2
SPLIT_SEED = 0          # 留出划分种子（与 step1_run 归档一致）
MODEL_SEED = 0          # 模型初始化种子（与 step1_run 归档一致）
HIDDEN = 64
BATCH = 32
TEST_BATCH = 128
TOP_N_DEG = 50          # pearson_delta 的 top-N 口径（step1_run）
EPOCHS_DEFAULT = 12
GATE_EPOCH = 3          # 阶段 1 复现判据 applied 的 epoch

# 探针 70/30 划分（与 step7_routing 完全同的派生代码）
SPLIT_SEEDS = [0, 1, 2]
TEST_FRAC = 0.30
MASTER_SEED = 20261016

FEATURE_ORDER = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
                 "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
                 "log_n_cells"]


# ---------------------------------------------------------------- 数据与模型初始化
def setup_data_and_model():
    """逐行镜像 step1_run.py 的数据/划分/初始化顺序（归档同款）。"""
    import torch
    from gears import PertData, GEARS

    t0 = time.time()
    pert_data = PertData(DATA_DIR)
    pert_data.load(data_name=DATA)
    adata = pert_data.adata
    print(f"[data] cells={adata.n_obs}, genes={adata.n_vars}, "
          f"conditions={adata.obs.condition.nunique()} ({time.time()-t0:.0f}s)",
          flush=True)

    # ---- 扰动留出划分（step1_run.py:60-80 原逻辑）----
    conditions = adata.obs.condition.unique().tolist()
    singles = [c for c in conditions if ("ctrl" in c and c != "ctrl")]
    genes = sorted({c.replace("ctrl+", "").replace("+ctrl", "") for c in singles})
    bad = [g for g in genes if "_" in g]
    if bad:
        genes = [g for g in genes if "_" not in g]

    rng = np.random.default_rng(SPLIT_SEED)
    holdout_genes = rng.choice(genes, size=int(HOLDOUT_FRAC * len(genes)),
                               replace=False).tolist()
    print(f"[split] 总单扰动基因={len(genes)}, 留出={len(holdout_genes)} "
          f"(frac={HOLDOUT_FRAC}, split_seed={SPLIT_SEED}, model_seed={MODEL_SEED})",
          flush=True)

    test_perts_arg = "_".join(holdout_genes)
    pert_data.prepare_split(split="simulation", seed=SPLIT_SEED,
                            test_perts=test_perts_arg,
                            only_test_set_perts=True)
    print(f"[split] set2conditions: "
          f"{ {k: len(v) for k, v in pert_data.set2conditions.items()} }",
          flush=True)

    pert_data.get_dataloader(batch_size=BATCH, test_batch_size=TEST_BATCH)

    import torch as _torch  # noqa: F401
    torch.manual_seed(MODEL_SEED)
    np.random.seed(MODEL_SEED)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = GEARS(pert_data, device=device)
    model.model_initialize(hidden_size=HIDDEN)
    print(f"[init] device={device}, hidden={HIDDEN} ({time.time()-t0:.0f}s)",
          flush=True)

    # ---- 真值侧静态量（快照无关，只算一次）----
    ctrl = adata[adata.obs.condition == "ctrl"]
    ctrl_mean = np.asarray(ctrl.X.mean(axis=0)).ravel()
    var_genes = adata.var["gene_name"].tolist()
    pred_genes = list(model.gene_list)
    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))

    holdout_conds = [c for c in singles
                     if c.replace("ctrl+", "").replace("+ctrl", "")
                     in holdout_genes]
    in_graph = [g for g in holdout_genes if g in model.pert_list]
    skipped = [g for g in holdout_genes if g not in model.pert_list]
    if skipped:
        print(f"[eval] {len(skipped)}/{len(holdout_genes)} 个留出基因不在扰动图中："
              f"{skipped[:5]}", flush=True)

    # 每个留出基因的真实 delta（对齐 pred_genes 顺序）与 n_cells
    true_delta, n_cells_map = {}, {}
    for g in holdout_genes:
        cand = [c for c in singles
                if c.startswith(g + "+") or c.endswith("+" + g)]
        if not cand:
            continue
        sub = adata[adata.obs.condition == cand[0]]
        t = np.asarray(sub.X.mean(axis=0)).ravel()
        if var_genes != pred_genes:
            t = np.array([dict(zip(var_genes, t))[x] for x in pred_genes])
        true_delta[g] = (t - ctrl_aligned).astype(np.float64)
        n_cells_map[g] = int(sub.n_obs)

    # n_cells 协变量以归档 predictions_seed0.csv 为准（实验设计量，两处应一致）
    if os.path.exists(PRED_CSV):
        pdf = pd.read_csv(PRED_CSV).set_index("perturbation")
        for g in n_cells_map:
            if g in pdf.index and not np.isnan(pdf.loc[g, "n_cells"]):
                n_cells_map[g] = int(pdf.loc[g, "n_cells"])

    ctx = {
        "model": model, "device": device, "adata": adata,
        "holdout_genes": holdout_genes, "in_graph": in_graph,
        "pred_genes": pred_genes, "ctrl_aligned": ctrl_aligned,
        "true_delta": true_delta, "n_cells_map": n_cells_map,
    }
    return ctx


# ---------------------------------------------------------------- 训练循环（镜像 gears.train）
def train_one_epoch(model, optimizer, scheduler, device, epoch_idx):
    """镜像 venv gears/gears.py GEARS.train() 的单 epoch（含兼容补丁 max_batches=500）。"""
    import torch
    import torch.nn as nn
    from gears.inference import evaluate, compute_metrics
    from gears.utils import loss_fct, print_sys

    train_loader = model.dataloader["train_loader"]
    val_loader = model.dataloader["val_loader"]

    model.model.train()
    for step, batch in enumerate(train_loader):
        batch.to(device)
        optimizer.zero_grad()
        y = batch.y
        pred = model.model(batch)
        loss = loss_fct(pred, y, batch.pert,
                        ctrl=model.ctrl_expression,
                        dict_filter=model.dict_filter,
                        direction_lambda=model.config["direction_lambda"])
        loss.backward()
        nn.utils.clip_grad_value_(model.model.parameters(), clip_value=1.0)
        optimizer.step()
        if step % 50 == 0:
            print_sys(f"Epoch {epoch_idx} Step {step + 1} "
                      f"Train Loss: {loss.item():.4f}")

    scheduler.step()
    train_res = evaluate(train_loader, model.model, model.config["uncertainty"],
                         device, max_batches=500)
    val_res = evaluate(val_loader, model.model, model.config["uncertainty"],
                       device)
    train_metrics, _ = compute_metrics(train_res)
    val_metrics, _ = compute_metrics(val_res)
    print_sys(f"Epoch {epoch_idx}: Train Overall MSE: "
              f"{train_metrics['mse']:.4f} Validation Overall MSE: "
              f"{val_metrics['mse']:.4f}. ")
    print_sys(f"Train Top 20 DE MSE: {train_metrics['mse_de']:.4f} "
              f"Validation Top 20 DE MSE: {val_metrics['mse_de']:.4f}. ")
    return {"train_mse": float(train_metrics["mse"]),
            "val_mse": float(val_metrics["mse"]),
            "train_mse_de": float(train_metrics["mse_de"]),
            "val_mse_de": float(val_metrics["mse_de"])}


# ---------------------------------------------------------------- 快照保存/恢复
def save_snapshot(model, optimizer, scheduler, epoch, min_val, path):
    import torch
    payload = {
        "model": {k: v.detach().cpu() for k, v in
                  model.model.state_dict().items()},
        "optimizer": optimizer.state_dict(),
        "scheduler": scheduler.state_dict(),
        "epoch": epoch,
        "min_val": float(min_val),
        "rng_torch": torch.get_rng_state(),
        "rng_cuda": (torch.cuda.get_rng_state() if torch.cuda.is_available()
                     else None),
        "rng_np": np.random.get_state(),
    }
    torch.save(payload, path)


def load_snapshot(model, optimizer, scheduler, path):
    import torch
    payload = torch.load(path, map_location=model.device, weights_only=False)
    model.model.load_state_dict(payload["model"])
    optimizer.load_state_dict(payload["optimizer"])
    scheduler.load_state_dict(payload["scheduler"])
    torch.set_rng_state(payload["rng_torch"])
    if payload.get("rng_cuda") is not None and torch.cuda.is_available():
        torch.cuda.set_rng_state(payload["rng_cuda"])
    np.random.set_state(payload["rng_np"])
    return int(payload["epoch"]), float(payload["min_val"])


def latest_snapshot(out_dir):
    snaps = [f for f in os.listdir(out_dir)
             if f.startswith("snap_e") and f.endswith(".pt")]
    if not snaps:
        return None
    epochs = [int(f[len("snap_e"):-len(".pt")]) for f in snaps]
    return max(epochs)


# ---------------------------------------------------------------- 探针分析（口径复用）
def analyze_snapshot(ctx):
    """对当前 epoch 权重量：217 预测 → pearson_delta → 9 特征 → 3-split 就地探针。

    返回 record dict（全部数字可 JSON 序列化）。
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    from step6_gating_probe import compute_features

    model = ctx["model"]
    in_graph = ctx["in_graph"]
    pred_genes = ctx["pred_genes"]
    ctrl = ctx["ctrl_aligned"].astype(np.float32)
    true_delta = ctx["true_delta"]

    # ---- 预测（用当前 epoch 权重，清缓存防串快照）----
    model.best_model = model.model
    model.saved_pred = {}
    pred = model.predict([[g] for g in in_graph])
    pred_mat = np.stack([np.asarray(pred[g]).ravel() for g in in_graph])

    # ---- pearson_delta（step1_run 口径：top-50 |true_delta| 上 Pearson）----
    pearson = {}
    for i, g in enumerate(in_graph):
        p_delta = pred_mat[i].astype(np.float64) - ctrl
        t_delta = true_delta[g]
        top_idx = np.argsort(-np.abs(t_delta))[:TOP_N_DEG]
        pearson[g] = float(np.corrcoef(p_delta[top_idx],
                                       t_delta[top_idx])[0, 1])

    # ---- 9 特征（step6 口径）----
    feat = compute_features(pred_mat, in_graph, pred_genes, ctrl)
    feat["log_n_cells"] = np.log1p(
        [ctx["n_cells_map"][g] for g in in_graph])
    X = feat[FEATURE_ORDER].values.astype(np.float64)
    perts = feat.index.tolist()
    err = np.array([1.0 - pearson[g] for g in perts])

    # ---- 3-split 就地探针（step7 同 Generator 派生）----
    n = len(perts)
    rng_master = np.random.default_rng(MASTER_SEED)
    per_seed = []
    for _sp in SPLIT_SEEDS:
        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        perm = rng.permutation(n)
        n_te = int(round(TEST_FRAC * n))
        te, tr = perm[:n_te], perm[n_te:]

        err_tr, err_te = err[tr], err[te]
        med_tr = np.median(err_tr)
        y_tr = (err_tr > med_tr).astype(int)
        y_te = (err_te > med_tr).astype(int)

        scaler = StandardScaler().fit(X[tr])
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(
            scaler.transform(X[tr]), y_tr)
        coef = clf.coef_[0].astype(np.float64)
        score_te = scaler.transform(X[te]) @ coef + clf.intercept_[0]
        auroc = float(roc_auc_score(y_te, score_te))
        per_seed.append({"coef": coef.tolist(), "auroc": auroc,
                         "n_train": int(len(tr)), "n_test": int(len(te))})

    # ---- 全量 217 就地拟合（= 归档 error_signatures.json 口径）----
    y_full = (err > np.median(err)).astype(int)
    scaler_full = StandardScaler().fit(X)
    clf_full = LogisticRegression(C=1.0, max_iter=2000).fit(
        scaler_full.transform(X), y_full)
    coef_full = clf_full.coef_[0].astype(np.float64)

    coef_mean = np.mean([s["coef"] for s in per_seed], axis=0)

    return {
        "n_perturbations": n,
        "pearson_delta": {g: pearson[g] for g in perts},
        "pearson_delta_mean": float(np.mean(list(pearson.values()))),
        "probe": {
            "per_split_seed": [
                {"split_seed": sp, "coef": s["coef"], "auroc": s["auroc"]}
                for sp, s in zip(SPLIT_SEEDS, per_seed)],
            "coef_mean": coef_mean.tolist(),
            "auroc_mean": float(np.mean([s["auroc"] for s in per_seed])),
        },
        "full_fit_coef": coef_full.tolist(),
    }


def cos_std(a, b):
    """标准化向量的余弦（系数尺度不稳的已知风险 → 一律用单位向量）。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return float("nan")
    return float(a @ b / (na * nb))


def load_archived_gears_sig():
    with open(ARCHIVED_SIG_JSON, encoding="utf-8") as f:
        sigs = json.load(f)
    return np.array([sigs["GEARS(K562)"][feat] for feat in FEATURE_ORDER],
                    dtype=np.float64)


def append_record(out_dir, record):
    path = os.path.join(out_dir, "per_epoch.json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    else:
        data = {"records": []}
    data["records"] = [r for r in data["records"]
                       if r["epoch"] != record["epoch"]]
    data["records"].append(record)
    data["records"].sort(key=lambda r: r["epoch"])
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)
    os.replace(tmp, path)


# ---------------------------------------------------------------- 主流程
def run_training(epochs, resume):
    import torch
    import torch.optim as optim
    from torch.optim.lr_scheduler import StepLR

    os.makedirs(OUT, exist_ok=True)
    ctx = setup_data_and_model()
    model = ctx["model"]
    device = ctx["device"]
    # 必须在创建 optimizer 之前上设备：.to(device) 会替换 Parameter 张量，
    # 先建 optimizer 会持有陈旧的 CPU 参数引用（epoch0 分析会把模型挪到 GPU）。
    model.model = model.model.to(device)

    optimizer = optim.Adam(model.model.parameters(), lr=1e-3,
                           weight_decay=5e-4)
    scheduler = StepLR(optimizer, step_size=1, gamma=0.5)
    min_val = np.inf
    start_epoch = 0

    if resume:
        last = latest_snapshot(OUT)
        if last is None:
            print("[resume] 未找到快照，从头开始", flush=True)
        else:
            path = os.path.join(OUT, f"snap_e{last}.pt")
            start_epoch, min_val = load_snapshot(
                model, optimizer, scheduler, path)
            print(f"[resume] 从 {path} 恢复：epoch={start_epoch}, "
                  f"min_val={min_val:.4f}", flush=True)

    archived_sig = load_archived_gears_sig()

    # ---- epoch 0 基线（仅初始化权重，未训练）----
    if start_epoch == 0:
        t_a = time.time()
        rec0 = analyze_snapshot(ctx)
        rec0.update({"epoch": 0, "phase": "baseline",
                     "elapsed_s": 0.0,
                     "gate_cos_full_vs_archived": None,
                     "gate_cos_152mean_vs_archived": None})
        append_record(OUT, rec0)
        print(f"[epoch 0/基线] pearson_delta_mean="
              f"{rec0['pearson_delta_mean']:.4f}, "
              f"probe AUROC={rec0['probe']['auroc_mean']:.3f} "
              f"(分析 {time.time()-t_a:.0f}s)", flush=True)

    t_train0 = time.time()
    for epoch in range(start_epoch + 1, epochs + 1):
        t_e = time.time()
        metrics = train_one_epoch(model, optimizer, scheduler, device, epoch)
        if metrics["val_mse_de"] < min_val:
            min_val = metrics["val_mse_de"]

        snap_path = os.path.join(OUT, f"snap_e{epoch}.pt")
        save_snapshot(model, optimizer, scheduler, epoch, min_val, snap_path)
        print(f"[snap] {snap_path}", flush=True)

        rec = analyze_snapshot(ctx)
        rec.update({"epoch": epoch, **metrics,
                    "elapsed_s": time.time() - t_train0})
        if epoch == GATE_EPOCH:
            cos_full = cos_std(rec["full_fit_coef"], archived_sig)
            cos_152 = cos_std(rec["probe"]["coef_mean"], archived_sig)
            rec["gate_cos_full_vs_archived"] = cos_full
            rec["gate_cos_152mean_vs_archived"] = cos_152
            verdict = "PASS (管线有效)" if cos_full > 0.8 else \
                "FAIL (<0.8，先查 split/数据口径)"
            print(f"[GATE e{GATE_EPOCH}] 复现判据：全量拟合 cos={cos_full:.4f}，"
                  f"152 均值拟合 cos={cos_152:.4f} → {verdict}", flush=True)
        append_record(OUT, rec)
        print(f"[epoch {epoch}/{epochs}] pearson_delta_mean="
              f"{rec['pearson_delta_mean']:.4f}, "
              f"probe AUROC={rec['probe']['auroc_mean']:.3f}, "
              f"val_mse={metrics['val_mse']:.4f} "
              f"(本 epoch 共 {time.time()-t_e:.0f}s)", flush=True)

    finalize(json_path=os.path.join(REPO, "outputs", "training_dynamics.json"),
             png_path=os.path.join(REPO, "outputs", "training_dynamics.png"))
    print("[done] 训练完成，最终产物已生成", flush=True)


# ---------------------------------------------------------------- 终态产物（json + png）
def finalize(json_path, png_path):
    per_epoch_path = os.path.join(OUT, "per_epoch.json")
    with open(per_epoch_path, encoding="utf-8") as f:
        data = json.load(f)
    records = sorted(data["records"], key=lambda r: r["epoch"])
    archived_sig = load_archived_gears_sig()

    final = next(r for r in records if r["epoch"] == max(x["epoch"]
                                                         for x in records))
    sig_final = np.array(final["probe"]["coef_mean"], dtype=np.float64)
    final_auroc = final["probe"]["auroc_mean"]
    final_pearson = final["pearson_delta_mean"]

    for r in records:
        cos_f = cos_std(r["probe"]["coef_mean"], sig_final)
        cos_seeds = [cos_std(s["coef"], sig_final)
                     for s in r["probe"]["per_split_seed"]]
        r["cos_to_final"] = cos_f
        r["cos_to_final_per_seed"] = cos_seeds
        if r["epoch"] != final["epoch"]:
            r.setdefault("gate_cos_full_vs_archived", None)
            r.setdefault("gate_cos_152mean_vs_archived", None)
        r["cos_full_vs_archived"] = cos_std(r["full_fit_coef"], archived_sig)

    def first_epoch(pred):
        for r in records:
            if pred(r):
                return int(r["epoch"])
        return None

    k1 = first_epoch(lambda r: r["cos_to_final"] >= 0.8)
    k2 = first_epoch(lambda r: abs(r["probe"]["auroc_mean"] - final_auroc)
                     <= 0.02)
    k3 = first_epoch(lambda r: abs(r["pearson_delta_mean"] - final_pearson)
                     <= 0.01)

    res = {
        "task": "training dynamics of the GEARS error signature "
                "(exploratory, entry A, 2026-10-16)",
        "dataset": DATA, "model_seed": MODEL_SEED, "split_seed": SPLIT_SEED,
        "n_perturbations": final["n_perturbations"],
        "epochs_total": max(r["epoch"] for r in records),
        "gate": {
            "epoch": GATE_EPOCH,
            "cos_full_fit_vs_archived":
                final_gate(records, "gate_cos_full_vs_archived"),
            "cos_152mean_vs_archived":
                final_gate(records, "gate_cos_152mean_vs_archived"),
            "criterion": "full-fit cos > 0.8 vs error_signatures.json "
                         "GEARS(K562) (same fitting protocol as archive)"},
        "timing": {
            "k1_form": k1, "k2_auroc": k2, "k3_pearson": k3,
            "T3_verdict": (
                "form-early (k1<k2 and k1<k3): signature shape stabilizes "
                "before probe function" if (k1 is not None and k2 is not None
                                          and k3 is not None
                                          and k1 < k2 and k1 < k3)
                else "emergent/late or indeterminate (k1>=k2 or k1>=k3)")},
        "caveats": [
            "cosine computed on L2-normalized coefficient vectors "
            "(regularized logistic + collinear features make raw scale "
            "unstable)",
            "single model, single seed, single dataset (n=1): timing "
            "conclusions are illustrative, not inferential",
            "epoch-0 baseline is the random-init model",
            "resume runs are bit-comparable to uninterrupted runs except "
            "for CUDA atomic nondeterminism"],
        "records": records,
    }
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    print(f"[SAVED] {json_path}", flush=True)

    make_figure(res, png_path)


def final_gate(records, key):
    for r in records:
        if r.get(key) is not None:
            return r[key]
    return None


def make_figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    recs = sorted(res["records"], key=lambda r: r["epoch"])
    ep = [r["epoch"] for r in recs]
    cos = [r["cos_to_final"] for r in recs]
    cos_min = [min(r["cos_to_final_per_seed"]) for r in recs]
    auroc = [r["probe"]["auroc_mean"] for r in recs]
    pear = [r["pearson_delta_mean"] for r in recs]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax1.plot(ep, cos, "o-", color="#1f77b4", lw=1.6, label="mean over 3 split seeds")
    ax1.fill_between(ep, cos_min, cos, color="#1f77b4", alpha=0.15,
                     label="min over 3 split seeds")
    ax1.axhline(0.8, color="red", ls="--", lw=1.2)
    ax1.text(ep[-1], 0.805, "0.8", color="red", fontsize=8, ha="right")
    ax1.set_xlabel("epoch")
    ax1.set_ylabel("cosine to final-epoch signature")
    ax1.set_title("Signature shape vs epoch (L2-normalized coefs)")
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8, loc="lower right")

    final_auroc = auroc[-1]
    final_pear = pear[-1]
    ax2.plot(ep, auroc, "s-", color="#2ca02c", lw=1.6,
             label="probe AUROC (mean, 3 seeds)")
    ax2.axhspan(final_auroc - 0.02, final_auroc + 0.02, color="#2ca02c",
                alpha=0.12)
    ax2.set_xlabel("epoch")
    ax2.set_ylabel("probe AUROC (65 held-out perts)", color="#2ca02c")
    ax2.tick_params(axis="y", labelcolor="#2ca02c")
    ax2.set_ylim(0.35, 1.0)

    ax2b = ax2.twinx()
    ax2b.plot(ep, pear, "^-", color="#d62728", lw=1.6,
              label="pearson_delta mean (217 perts)")
    ax2b.axhspan(final_pear - 0.01, final_pear + 0.01, color="#d62728",
                 alpha=0.10)
    ax2b.set_ylabel("mean top-50 pearson_delta", color="#d62728")
    ax2b.tick_params(axis="y", labelcolor="#d62728")

    lines = ax2.get_lines() + ax2b.get_lines()
    ax2.legend(lines, [l.get_label() for l in lines], fontsize=8,
               loc="center right")
    ax2.set_title("Probe function & prediction quality vs epoch\n"
                  "(bands = final $\\pm$0.02 AUROC / $\\pm$0.01 pearson)")
    ax2.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=EPOCHS_DEFAULT,
                    help="总训练轮数（降级档：6）")
    ap.add_argument("--resume-snap", action="store_true",
                    help="从 outputs/training_dynamics/ 中最后一个快照续跑")
    ap.add_argument("--plot-only", action="store_true",
                    help="只从 per_epoch.json 重新生成最终 json/png")
    args = ap.parse_args()

    if args.plot_only:
        finalize(json_path=os.path.join(REPO, "outputs",
                                        "training_dynamics.json"),
                 png_path=os.path.join(REPO, "outputs",
                                       "training_dynamics.png"))
        return
    run_training(args.epochs, args.resume_snap)


if __name__ == "__main__":
    main()
