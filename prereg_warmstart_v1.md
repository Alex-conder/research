# Warm-start transfer calibration plan v1

- Date: 2026-10-30
- Status: recorded before running `step19_warmstart_transfer.py`.
- Motivation: target-normalized continuous coefficients help in some domains, but few-shot calibration can be stronger. We test whether source coefficients also provide a useful initialization for target calibration.

## Method

For each stratified calibration split of n=10, 15, 30, 60 labels:

- Fit feature scaler on calibration only.
- Optimize an L2 logistic probe (C=1.0) from the archived K562 coefficient vector and intercept, using L-BFGS-B.
- Compare against the same split's from-scratch L2 logistic probe, target-z transfer without labels, and full in-situ OOF.

Primary evaluation: AUROC, precision-coverage, and paired bootstrap for `warmstart_30 - fewshot_30` and `warmstart_30 - coef_target_z`. Also report label-efficiency curves (n=10, 15, 30, 60).

## Decision criteria

- Warm start is useful only if it improves over from-scratch few-shot at n<=30 with nominal 95% CI excluding zero, and remains above chance. It is not required to beat full in-situ.
- If warm start helps on one dataset but not the other, report it as domain-dependent. No universal claim.
- Bootstrap B=500, master seed 20261029; 20 repeated splits for point estimates, 5 for bootstrap.
