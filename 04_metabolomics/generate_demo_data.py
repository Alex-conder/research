#!/usr/bin/env python3
"""Generate a small synthetic SPSS .sav dataset for pipeline testing."""
import numpy as np
import pandas as pd
from pathlib import Path

np.random.seed(42)

n_control = 8
n_treatment = 8
n_features = 60

# Feature m/z and RT
mz = np.linspace(100, 800, n_features) + np.random.normal(0, 1, n_features)
rt = np.linspace(1, 15, n_features) + np.random.normal(0, 0.2, n_features)

# Base intensities
base = np.random.lognormal(mean=8, sigma=1.5, size=n_features)

# Make ~10 features differentially abundant
diff_idx = np.random.choice(n_features, size=10, replace=False)

def make_samples(n, group_factor):
    samples = []
    for _ in range(n):
        x = base * np.random.lognormal(0, 0.2, n_features)
        x[diff_idx] *= group_factor
        # Introduce some missing values
        mask = np.random.rand(n_features) < 0.03
        x[mask] = np.nan
        samples.append(x)
    return np.array(samples)

ctrl = make_samples(n_control, 1.0)
trt = make_samples(n_treatment, 3.0)  # treatment up-regulated

X = np.vstack([ctrl, trt])
groups = ["Control"] * n_control + ["Treatment"] * n_treatment

# Column names as Mz_RT format
col_names = [f"M{mz[i]:.4f}_T{rt[i]:.4f}" for i in range(n_features)]
df = pd.DataFrame(X, columns=col_names)
df.insert(0, "Group", groups)

Path("./data").mkdir(exist_ok=True)
import pyreadstat
pyreadstat.write_sav(df, "./data/raw_data.sav")
print("Created ./data/raw_data.sav")

# Demo local library
lib = pd.DataFrame({
    "Name": [f"Metabolite_{i}" for i in diff_idx],
    "Mz": mz[diff_idx],
    "RT": rt[diff_idx],
})
Path("./library").mkdir(exist_ok=True)
lib.to_csv("./library/plant_library.csv", index=False)
print("Created ./library/plant_library.csv")
