# AD baseline dual-metric audit - verdict (step21, prereg_ad_baseline.md)

- generated: 2026-10-30 13:27:59
- smoke: False
- signals: AD = additive disagreement (AD-mean on single-gene holdout products per deviation D-AD1; AD-combo on Norman combos, paper-faithful per deviation D-AD2); symbol gate = two_layer_gate construction; probe = step6 spec, y refit per metric
- error metrics: full = transcriptome-wide MAE of delta; top20/top50 = MAE over top-20/50 genes by |true delta| (per perturbation); failure = error > median (paper Top50%)
- CIs: perturbation-level paired bootstrap of AUROC differences, B=2000, seed=20261030

## Analysis B (primary): Norman combos, archived GEARS seed4, CPU, no training

- candidates=57 eligible=57 predicted=57 excluded={'no_single_condition': 0, 'predict_failed': 0}

| metric | AD-combo | symbol gate | in-situ probe | AD - symbol [95% CI] | AD - probe [95% CI] |
|---|---|---|---|---|---|
| full | 0.978 | 0.335 | 0.999 | +0.642 [+0.471,+0.797] AD | -0.021 [-0.057,+0.001] undetermined |
| top20 | 0.550 | 0.583 | 0.475 | -0.034 [-0.279,+0.228] undetermined | +0.075 [+0.014,+0.140] AD |
| top50 | 0.616 | 0.462 | 0.541 | +0.157 [-0.090,+0.409] undetermined | +0.074 [+0.017,+0.135] AD |

## Verdict (locked two-tier)

- tier: **CONFIRMED**
- rule: CONFIRMED iff (a) full-metric paired bootstrap CI of (AD-combo - symbol) excludes 0 favoring AD, AND (b) top20-metric CI includes 0 or favors symbol. Everything else = NOT_CONFIRMED. Secondary: same for AD-combo - probe.
- (a) full AD-symbol CI excludes 0 favoring AD: True (mean +0.642, CI [+0.471, +0.797])
- (b) top20 reverses or non-significant: True (mean -0.034, CI [-0.279, +0.228], favors undetermined)

## Analysis A (secondary): 5-dataset single-gene holdout products, AD-mean (deviation D-AD1)

| dataset | seed | metric | AD-mean | symbol gate | in-situ probe | AD - symbol [95% CI] | AD - probe [95% CI] |
|---|---|---|---|---|---|---|---|
| replogle_k562_essential | 0 | full | 0.844 | 0.144 | 0.909 | +0.699 [+0.583,+0.806] AD | -0.066 [-0.122,-0.015] gate |
|  | 0 | top20 | 0.404 | 0.615 | 0.635 | -0.210 [-0.356,-0.062] gate | -0.230 [-0.361,-0.091] gate |
|  | 0 | top50 | 0.443 | 0.569 | 0.615 | -0.125 [-0.274,+0.026] undetermined | -0.171 [-0.293,-0.042] gate |
| replogle_k562_essential | 1 | full | 0.863 | 0.094 | 0.900 | +0.767 [+0.683,+0.843] AD | -0.037 [-0.078,+0.002] undetermined |
|  | 1 | top20 | 0.504 | 0.456 | 0.501 | +0.047 [-0.103,+0.195] undetermined | +0.003 [-0.093,+0.102] undetermined |
|  | 1 | top50 | 0.524 | 0.433 | 0.544 | +0.091 [-0.060,+0.236] undetermined | -0.019 [-0.112,+0.078] undetermined |
| replogle_k562_essential | 2 | full | 0.719 | 0.215 | 0.769 | +0.505 [+0.383,+0.625] AD | -0.050 [-0.107,+0.008] undetermined |
|  | 2 | top20 | 0.536 | 0.421 | 0.524 | +0.112 [-0.034,+0.253] undetermined | +0.012 [-0.070,+0.094] undetermined |
|  | 2 | top50 | 0.523 | 0.431 | 0.575 | +0.091 [-0.051,+0.232] undetermined | -0.053 [-0.141,+0.033] undetermined |
| replogle_rpe1_essential | 0 | full | 0.469 | 0.445 | 0.821 | +0.024 [-0.044,+0.095] undetermined | -0.349 [-0.432,-0.262] gate |
|  | 0 | top20 | 0.542 | 0.551 | 0.640 | -0.008 [-0.076,+0.062] undetermined | -0.100 [-0.173,-0.026] gate |
|  | 0 | top50 | 0.545 | 0.502 | 0.691 | +0.044 [-0.025,+0.109] undetermined | -0.146 [-0.229,-0.064] gate |
| norman | 4 | full | 0.880 | 0.360 | 0.680 | +0.518 [+0.150,+0.828] AD | +0.199 [+0.044,+0.404] AD |
|  | 4 | top20 | 0.720 | 0.480 | 0.440 | +0.241 [-0.202,+0.640] undetermined | +0.279 [+0.066,+0.517] AD |
|  | 4 | top50 | 0.880 | 0.360 | 0.680 | +0.515 [+0.154,+0.840] AD | +0.199 [+0.050,+0.394] AD |
| adamson | 3 | full | 0.688 | 0.438 | 0.531 | +0.252 [-0.267,+0.719] undetermined | +0.151 [-0.083,+0.383] undetermined |
|  | 3 | top20 | 0.609 | 0.172 | 0.672 | +0.435 [+0.000,+0.818] undetermined | -0.051 [-0.413,+0.270] undetermined |
|  | 3 | top50 | 0.609 | 0.172 | 0.672 | +0.437 [-0.017,+0.810] undetermined | -0.064 [-0.400,+0.270] undetermined |
| adamson | 5 | full | 0.625 | 0.375 | 0.766 | +0.255 [-0.209,+0.714] undetermined | -0.139 [-0.436,+0.156] undetermined |
|  | 5 | top20 | 0.625 | 0.375 | 0.766 | +0.256 [-0.219,+0.688] undetermined | -0.140 [-0.469,+0.156] undetermined |
|  | 5 | top50 | 0.625 | 0.375 | 0.766 | +0.247 [-0.218,+0.683] undetermined | -0.139 [-0.450,+0.145] undetermined |
| dixit | 4 | full | n/e | | | | |
|  | 4 | top20 | n/e | | | | |
|  | 4 | top50 | n/e | | | | |
| dixit | 5 | full | n/e | | | | |
|  | 5 | top20 | n/e | | | | |
|  | 5 | top50 | n/e | | | | |

## Deviations (locked disclosures)

- D-AD1: paper states AD is undefined for single-gene OOD. On single-gene holdout products we use AD-mean = |pred_delta - mean train single-gene delta| (mean-effect additive reference). Not the paper's per-gene AD; secondary evidence only, not in the verdict rule.
- D-AD2: our split makes 'both single-gene conditions in model training' and 'combo unseen by model' mutually exclusive. AD-combo baseline uses all observed single-gene conditions (incl. held-out singles; model never saw them). Evaluated combo truth never enters the baseline (no circularity). Paper home setting has both components in training; our setting is harder for the model.
- Dixit n=3: not evaluable under the minority>=2 rule (amendment-1), recorded as n/e like the existing probes.

## Where this goes in the paper (one paragraph)

Gating paper appendix 'Reliability signal comparison: additive disagreement (AD) baseline and dual-metric robustness'. Norman-combo panel is the AD paper's home setting run on our archived model (no retraining); the 5-dataset table shows our probe and symbol gate recomputed under full / top20 / top50 error definitions. Both metrics reported side by side; the Top20-MAE trap (paper 0.990 -> 0.619) is checked on our gates, and the dual-metric report is presented as the standard practice regardless of verdict tier.
