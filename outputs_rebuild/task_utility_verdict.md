# G-C task-level utility (step22, prereg_expansion_gc.md) - verdict

- generated: 2026-10-30 16:47:55
- smoke: False
- replay: GEARS single-gene holdout products (step21 analysis A identical data/protocol); zero training, CPU, on-file
- cost: C(B) = sum of trusted prediction errors (full transcriptome delta MAE) + c*(n-B); c = per-combo median err (arm comparisons at fixed B are c-free)
- statistics: perturbation-level paired bootstrap, B=1000, master seed=20261030; no model refits inside bootstrap (score vectors fixed, preregistered)
- locked rule: combo positive vs a baseline iff the paired 95% CI of (cost_baseline - cost_gate) excludes 0 favoring the gate at EVERY evaluable B; dataset positive iff ALL its seed combos are positive vs BOTH baselines; confirmed iff >=3/5 registered datasets (dixit degenerate by construction -> auto-negative)

## Gate arm: probe_full PRIMARY -> NOT_CONFIRMED (1/5)

| dataset | B | d(random-gate) [95% CI] | d(magnitude-gate) [95% CI] |
|---|---|---|---|
| gears_k562_s0 | 5 | +0.075 [+0.019,+0.127] gate | +0.309 [+0.213,+0.408] gate |
|  | 10 | +0.141 [+0.069,+0.218] gate | +0.563 [+0.428,+0.696] gate |
|  | 20 | +0.285 [+0.168,+0.392] gate | +0.977 [+0.747,+1.214] gate |
|  | 50 | +0.665 [+0.492,+0.827] gate | +1.861 [+1.552,+2.180] gate |
| gears_k562_s1 | 5 | +0.098 [+0.051,+0.147] gate | +0.336 [+0.231,+0.420] gate |
|  | 10 | +0.186 [+0.111,+0.256] gate | +0.555 [+0.342,+0.714] gate |
|  | 20 | +0.369 [+0.251,+0.472] gate | +0.894 [+0.620,+1.180] gate |
|  | 50 | +0.857 [+0.692,+1.014] gate | +1.924 [+1.587,+2.260] gate |
| gears_k562_s2 | 5 | -0.023 [-0.147,+0.065] undetermined | +0.090 [-0.020,+0.227] undetermined |
|  | 10 | -0.034 [-0.189,+0.102] undetermined | +0.110 [-0.052,+0.273] undetermined |
|  | 20 | +0.011 [-0.190,+0.212] undetermined | +0.327 [+0.069,+0.579] gate |
|  | 50 | +0.366 [+0.154,+0.537] gate | +0.847 [+0.534,+1.176] gate |
| gears_rpe1_s0 | 5 | +0.139 [+0.077,+0.165] gate | +0.309 [+0.197,+0.417] gate |
|  | 10 | +0.244 [+0.122,+0.323] gate | +0.561 [+0.400,+0.722] gate |
|  | 20 | +0.493 [+0.346,+0.610] gate | +1.034 [+0.789,+1.248] gate |
|  | 50 | +0.941 [+0.700,+1.177] gate | +1.886 [+1.515,+2.246] gate |
| gears_norman_s4 | 5 | +0.014 [-0.003,+0.030] undetermined | +0.043 [+0.005,+0.091] gate |
|  | 10 | +0.028 [+0.004,+0.055] gate | +0.053 [+0.011,+0.106] gate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s3 | 5 | +0.002 [-0.039,+0.040] undetermined | +0.009 [-0.046,+0.063] undetermined |
|  | 10 | +0.010 [-0.015,+0.037] undetermined | +0.029 [-0.032,+0.086] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s5 | 5 | +0.012 [-0.018,+0.041] undetermined | +0.030 [-0.014,+0.068] undetermined |
|  | 10 | +0.020 [-0.019,+0.051] undetermined | +0.038 [-0.022,+0.087] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s4 | 5 | n/e | n/e |
| gears_dixit_s4 | 10 | n/e | n/e |
| gears_dixit_s4 | 20 | n/e | n/e |
| gears_dixit_s4 | 50 | n/e | n/e |
| gears_dixit_s5 | 5 | n/e | n/e |
| gears_dixit_s5 | 10 | n/e | n/e |
| gears_dixit_s5 | 20 | n/e | n/e |
| gears_dixit_s5 | 50 | n/e | n/e |

dataset_positive: replogle_k562_essential=N, replogle_rpe1_essential=Y, norman=N, adamson=N, dixit=N

## Gate arm: probe_top20 -> NOT_CONFIRMED (1/5)

| dataset | B | d(random-gate) [95% CI] | d(magnitude-gate) [95% CI] |
|---|---|---|---|
| gears_k562_s0 | 5 | -0.011 [-0.077,+0.055] undetermined | +0.225 [+0.118,+0.320] gate |
|  | 10 | -0.016 [-0.108,+0.069] undetermined | +0.402 [+0.233,+0.551] gate |
|  | 20 | -0.041 [-0.178,+0.094] undetermined | +0.652 [+0.417,+0.881] gate |
|  | 50 | -0.189 [-0.402,+0.019] undetermined | +1.005 [+0.720,+1.336] gate |
| gears_k562_s1 | 5 | +0.062 [+0.005,+0.114] gate | +0.299 [+0.200,+0.382] gate |
|  | 10 | +0.156 [+0.085,+0.219] gate | +0.523 [+0.312,+0.683] gate |
|  | 20 | +0.265 [+0.154,+0.371] gate | +0.795 [+0.521,+1.071] gate |
|  | 50 | +0.478 [+0.277,+0.672] gate | +1.545 [+1.202,+1.909] gate |
| gears_k562_s2 | 5 | -0.040 [-0.164,+0.060] undetermined | +0.074 [-0.039,+0.192] undetermined |
|  | 10 | -0.037 [-0.194,+0.083] undetermined | +0.100 [-0.036,+0.242] undetermined |
|  | 20 | -0.036 [-0.207,+0.123] undetermined | +0.275 [+0.078,+0.496] gate |
|  | 50 | +0.250 [+0.003,+0.467] gate | +0.734 [+0.401,+1.082] gate |
| gears_rpe1_s0 | 5 | +0.132 [+0.064,+0.160] gate | +0.303 [+0.197,+0.405] gate |
|  | 10 | +0.224 [+0.111,+0.300] gate | +0.539 [+0.380,+0.694] gate |
|  | 20 | +0.363 [+0.182,+0.518] gate | +0.908 [+0.641,+1.178] gate |
|  | 50 | +0.699 [+0.462,+0.940] gate | +1.637 [+1.256,+2.038] gate |
| gears_norman_s4 | 5 | +0.002 [-0.024,+0.024] undetermined | +0.031 [-0.001,+0.079] undetermined |
|  | 10 | +0.015 [-0.012,+0.043] undetermined | +0.039 [+0.006,+0.087] gate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s3 | 5 | +0.010 [-0.030,+0.042] undetermined | +0.017 [-0.036,+0.058] undetermined |
|  | 10 | +0.011 [-0.016,+0.040] undetermined | +0.028 [-0.031,+0.088] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s5 | 5 | +0.012 [-0.021,+0.040] undetermined | +0.030 [-0.014,+0.070] undetermined |
|  | 10 | +0.020 [-0.017,+0.050] undetermined | +0.039 [-0.017,+0.088] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s4 | 5 | n/e | n/e |
| gears_dixit_s4 | 10 | n/e | n/e |
| gears_dixit_s4 | 20 | n/e | n/e |
| gears_dixit_s4 | 50 | n/e | n/e |
| gears_dixit_s5 | 5 | n/e | n/e |
| gears_dixit_s5 | 10 | n/e | n/e |
| gears_dixit_s5 | 20 | n/e | n/e |
| gears_dixit_s5 | 50 | n/e | n/e |

dataset_positive: replogle_k562_essential=N, replogle_rpe1_essential=Y, norman=N, adamson=N, dixit=N

## Gate arm: probe_top50 -> NOT_CONFIRMED (1/5)

| dataset | B | d(random-gate) [95% CI] | d(magnitude-gate) [95% CI] |
|---|---|---|---|
| gears_k562_s0 | 5 | +0.011 [-0.044,+0.072] undetermined | +0.247 [+0.152,+0.340] gate |
|  | 10 | +0.056 [-0.026,+0.138] undetermined | +0.474 [+0.324,+0.612] gate |
|  | 20 | +0.099 [-0.043,+0.219] undetermined | +0.792 [+0.546,+1.011] gate |
|  | 50 | +0.161 [-0.052,+0.373] undetermined | +1.362 [+1.041,+1.703] gate |
| gears_k562_s1 | 5 | +0.084 [+0.031,+0.132] gate | +0.322 [+0.227,+0.404] gate |
|  | 10 | +0.153 [+0.085,+0.212] gate | +0.513 [+0.288,+0.681] gate |
|  | 20 | +0.245 [+0.139,+0.348] gate | +0.775 [+0.516,+1.040] gate |
|  | 50 | +0.489 [+0.271,+0.695] gate | +1.552 [+1.196,+1.890] gate |
| gears_k562_s2 | 5 | -0.052 [-0.172,+0.045] undetermined | +0.061 [-0.055,+0.204] undetermined |
|  | 10 | -0.080 [-0.217,+0.040] undetermined | +0.065 [-0.066,+0.216] undetermined |
|  | 20 | -0.066 [-0.292,+0.138] undetermined | +0.248 [+0.048,+0.486] gate |
|  | 50 | +0.262 [+0.027,+0.459] gate | +0.747 [+0.464,+1.035] gate |
| gears_rpe1_s0 | 5 | +0.126 [+0.025,+0.165] gate | +0.299 [+0.166,+0.413] gate |
|  | 10 | +0.214 [+0.096,+0.299] gate | +0.527 [+0.364,+0.698] gate |
|  | 20 | +0.413 [+0.251,+0.549] gate | +0.955 [+0.727,+1.191] gate |
|  | 50 | +0.817 [+0.572,+1.045] gate | +1.751 [+1.350,+2.134] gate |
| gears_norman_s4 | 5 | +0.014 [-0.004,+0.031] undetermined | +0.042 [+0.003,+0.092] gate |
|  | 10 | +0.028 [+0.006,+0.054] gate | +0.053 [+0.012,+0.102] gate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s3 | 5 | +0.010 [-0.027,+0.041] undetermined | +0.016 [-0.035,+0.059] undetermined |
|  | 10 | +0.011 [-0.017,+0.039] undetermined | +0.030 [-0.033,+0.088] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s5 | 5 | +0.012 [-0.021,+0.041] undetermined | +0.029 [-0.015,+0.067] undetermined |
|  | 10 | +0.019 [-0.018,+0.051] undetermined | +0.037 [-0.016,+0.088] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s4 | 5 | n/e | n/e |
| gears_dixit_s4 | 10 | n/e | n/e |
| gears_dixit_s4 | 20 | n/e | n/e |
| gears_dixit_s4 | 50 | n/e | n/e |
| gears_dixit_s5 | 5 | n/e | n/e |
| gears_dixit_s5 | 10 | n/e | n/e |
| gears_dixit_s5 | 20 | n/e | n/e |
| gears_dixit_s5 | 50 | n/e | n/e |

dataset_positive: replogle_k562_essential=N, replogle_rpe1_essential=Y, norman=N, adamson=N, dixit=N

## Gate arm: symbol -> NOT_CONFIRMED (0/5)

| dataset | B | d(random-gate) [95% CI] | d(magnitude-gate) [95% CI] |
|---|---|---|---|
| gears_k562_s0 | 5 | -0.229 [-0.290,-0.175] baseline | +0.007 [-0.083,+0.090] undetermined |
|  | 10 | -0.443 [-0.544,-0.338] baseline | -0.020 [-0.136,+0.082] undetermined |
|  | 20 | -0.717 [-0.885,-0.567] baseline | -0.023 [-0.109,+0.052] undetermined |
|  | 50 | -1.236 [-1.451,-1.036] baseline | -0.036 [-0.149,+0.049] undetermined |
| gears_k562_s1 | 5 | -0.227 [-0.278,-0.181] baseline | +0.009 [-0.079,+0.077] undetermined |
|  | 10 | -0.442 [-0.539,-0.358] baseline | -0.077 [-0.275,+0.058] undetermined |
|  | 20 | -0.758 [-0.909,-0.609] baseline | -0.229 [-0.442,-0.058] baseline |
|  | 50 | -1.317 [-1.532,-1.072] baseline | -0.249 [-0.436,-0.090] baseline |
| gears_k562_s2 | 5 | -0.166 [-0.243,-0.103] baseline | -0.047 [-0.159,+0.067] undetermined |
|  | 10 | -0.314 [-0.420,-0.212] baseline | -0.172 [-0.307,-0.024] baseline |
|  | 20 | -0.506 [-0.660,-0.357] baseline | -0.200 [-0.375,-0.058] baseline |
|  | 50 | -0.868 [-1.065,-0.691] baseline | -0.388 [-0.573,-0.204] baseline |
| gears_rpe1_s0 | 5 | -0.102 [-0.239,+0.025] undetermined | +0.068 [-0.056,+0.225] undetermined |
|  | 10 | -0.236 [-0.382,-0.094] baseline | +0.077 [-0.066,+0.239] undetermined |
|  | 20 | -0.421 [-0.618,-0.228] baseline | +0.122 [+0.007,+0.248] gate |
|  | 50 | -0.697 [-1.017,-0.356] baseline | +0.242 [+0.064,+0.468] gate |
| gears_norman_s4 | 5 | -0.014 [-0.038,+0.007] undetermined | +0.015 [-0.009,+0.062] undetermined |
|  | 10 | -0.020 [-0.051,+0.005] undetermined | +0.005 [-0.008,+0.029] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s3 | 5 | -0.015 [-0.042,+0.010] undetermined | -0.008 [-0.037,+0.014] undetermined |
|  | 10 | -0.031 [-0.058,+0.004] undetermined | -0.013 [-0.049,+0.015] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s5 | 5 | -0.015 [-0.043,+0.011] undetermined | +0.002 [-0.029,+0.033] undetermined |
|  | 10 | -0.015 [-0.047,+0.026] undetermined | +0.003 [-0.045,+0.058] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s4 | 5 | degenerate | degenerate |
|  | 10 | degenerate | degenerate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s5 | 5 | degenerate | degenerate |
|  | 10 | degenerate | degenerate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |

dataset_positive: replogle_k562_essential=N, replogle_rpe1_essential=N, norman=N, adamson=N, dixit=N

## Gate arm: conformal_width -> NOT_CONFIRMED (1/5)

| dataset | B | d(random-gate) [95% CI] | d(magnitude-gate) [95% CI] |
|---|---|---|---|
| gears_k562_s0 | 5 | +0.047 [-0.048,+0.143] undetermined | +0.282 [+0.155,+0.399] gate |
|  | 10 | +0.041 [-0.079,+0.165] undetermined | +0.460 [+0.280,+0.632] gate |
|  | 20 | +0.087 [-0.078,+0.253] undetermined | +0.779 [+0.539,+1.041] gate |
|  | 50 | +0.513 [+0.273,+0.741] gate | +1.707 [+1.361,+2.079] gate |
| gears_k562_s1 | 5 | +0.112 [+0.055,+0.152] gate | +0.349 [+0.246,+0.434] gate |
|  | 10 | +0.192 [+0.116,+0.269] gate | +0.563 [+0.347,+0.738] gate |
|  | 20 | +0.369 [+0.258,+0.478] gate | +0.887 [+0.619,+1.162] gate |
|  | 50 | +0.839 [+0.670,+0.999] gate | +1.908 [+1.577,+2.253] gate |
| gears_k562_s2 | 5 | +0.029 [-0.042,+0.085] undetermined | +0.147 [+0.000,+0.304] gate |
|  | 10 | +0.084 [+0.002,+0.154] gate | +0.222 [+0.052,+0.398] gate |
|  | 20 | +0.198 [+0.102,+0.291] gate | +0.514 [+0.319,+0.727] gate |
|  | 50 | +0.517 [+0.381,+0.657] gate | +1.001 [+0.680,+1.334] gate |
| gears_rpe1_s0 | 5 | +0.075 [-0.060,+0.177] undetermined | +0.246 [+0.080,+0.396] gate |
|  | 10 | +0.093 [-0.057,+0.240] undetermined | +0.405 [+0.225,+0.597] gate |
|  | 20 | +0.206 [+0.045,+0.370] gate | +0.745 [+0.496,+0.999] gate |
|  | 50 | +0.348 [+0.089,+0.609] gate | +1.284 [+0.862,+1.675] gate |
| gears_norman_s4 | 5 | +0.025 [+0.014,+0.038] gate | +0.053 [+0.017,+0.102] gate |
|  | 10 | +0.035 [+0.017,+0.059] gate | +0.061 [+0.023,+0.111] gate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s3 | 5 | +0.018 [-0.023,+0.056] undetermined | +0.024 [-0.032,+0.076] undetermined |
|  | 10 | +0.007 [-0.019,+0.034] undetermined | +0.025 [-0.035,+0.085] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_adamson_s5 | 5 | +0.007 [-0.020,+0.037] undetermined | +0.025 [-0.014,+0.061] undetermined |
|  | 10 | +0.021 [-0.013,+0.051] undetermined | +0.038 [-0.013,+0.088] undetermined |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s4 | 5 | degenerate | degenerate |
|  | 10 | degenerate | degenerate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |
| gears_dixit_s5 | 5 | degenerate | degenerate |
|  | 10 | degenerate | degenerate |
|  | 20 | degenerate | degenerate |
|  | 50 | degenerate | degenerate |

dataset_positive: replogle_k562_essential=N, replogle_rpe1_essential=N, norman=Y, adamson=N, dixit=N

## Verdict (locked)
- tier: **NOT_CONFIRMED** (primary gate = probe_full)
- rule: primary = probe_full; per-gate verdicts also reported; negative result reported as-is (no utility is a result); R2 product unchanged (symbol layer)

## Where this goes in the paper (one paragraph)
Gating paper main text, 'task-level utility' panel: the gate is no longer only an error-scoring tool but a decision tool -- it decides which virtual-cell predictions may stand in for a real Perturb-seq measurement under an experimental budget (VCC 'replace one real experiment' narrative). Per-dataset utility tables above; positive verdict licenses the claim 'the gate has positive utility under experimental budget constraints'; a negative verdict is reported as-is (a gate that scores errors but cannot save experiments), which would itself be a tier-level conclusion about reliability signals in this field.
