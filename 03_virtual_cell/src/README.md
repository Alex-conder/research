# src 目录

本目录用于存放从 notebook 中抽离的可复用函数。

目前核心逻辑均内联在 `notebooks/` 中，以保证每个 notebook 可独立阅读。后续若项目扩大，建议将重复函数迁移至此。

示例结构（未来）：

```
src/
├── hh_utils.py       # HH 模型辅助函数
├── fba_utils.py      # FBA 分析与可视化
├── scvi_utils.py     # scVI 训练与扰动
└── interface.py      # 模块耦合接口
```
