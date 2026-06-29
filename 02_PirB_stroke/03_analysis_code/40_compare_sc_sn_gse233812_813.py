"""
比较 GSE233812 (scRNA-seq) 和 GSE233813 (snRNA-seq) 中小胶质细胞 Pirb 动态
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE233812_813_comparison"
os.makedirs(OUT_DIR, exist_ok=True)

# 读取 summary
sc = pd.read_csv("D:/Pirb_stroke_project/04_reports/figures/GSE233812/pirb_summary_celltype_time.csv")
sn = pd.read_csv("D:/Pirb_stroke_project/04_reports/figures/GSE233813/pirb_summary_celltype_time.csv")

sc['dataset'] = 'scRNA-seq (GSE233812)'
sn['dataset'] = 'snRNA-seq (GSE233813)'

combined = pd.concat([sc, sn], ignore_index=True)
micro = combined[combined['cell_type'] == 'Microglia'].copy()
micro['time_point_num'] = micro['time_point'].map({'sham': 0, 'D1': 1, 'D3': 3, 'D7': 7})
micro = micro.sort_values('time_point_num')

# 点图 + 连线
fig, ax = plt.subplots(figsize=(8, 5))
sns.pointplot(data=micro, x='time_point', y='pirb_frac', hue='dataset', 
              order=['sham', 'D1', 'D3', 'D7'], palette=['red', 'blue'], ax=ax, markers='o', dodge=0.1)
ax.set_ylabel('Pirb+ fraction in Microglia')
ax.set_title('Microglial Pirb expression across stroke time course')
ax.legend(title='Dataset')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "microglia_pirb_sc_vs_sn.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 保存对比表
micro[['dataset', 'time_point', 'n', 'pirb_frac', 'pirb_mean']].to_csv(
    os.path.join(OUT_DIR, "microglia_pirb_sc_vs_sn.csv"), index=False)
print(micro[['dataset', 'time_point', 'n', 'pirb_frac', 'pirb_mean']])
print("[DONE]")
