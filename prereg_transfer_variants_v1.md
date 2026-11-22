# Transfer variants plan v1 (route B)

- Date: 2026-10-30
- Status: recorded before running `step17_transfer_variants.py`.
- Motivation: sign-only transfer discards coefficient magnitudes, while z-scored continuous coefficients remain sensitive to heavy-tailed feature distributions. We test two principled robustifications, not a large subset search.

## Variants

1. `sign_layer`: target z-scored features dot coefficient signs.
2. `coef_target_z`: target z-scored features dot full coefficients.
3. `coef_median_iqr`: median/IQR-scaled features dot full coefficients (robust location/scale).
4. `coef_rank`: within-target percentile-ranked features dot full coefficients (rank-normalized transfer).
5. `stable2_z`: only `delta_max` and `delta_entropy`, the two K562-stable features from R1, with z-scoring.
6. `stable2_rank`: the same two features with rank normalization.
7. `delta_max_rank`: single `delta_max` feature with sign and rank normalization.
8. `magnitude4_z`: four continuous magnitude statistics with z-scoring.
9. `no_delta_max_z`: all features except `delta_max` with z-scoring (residual information diagnostic).

In-situ and 30-label few-shot probes remain comparators, not zero-label methods.

## Decision criteria

- A variant is considered an improvement only if its paired coverage-averaged precision difference versus `coef_target_z` is positive with nominal 95% CI excluding zero, and it remains above chance. We also report the contrast versus `delta_max_sign`; superiority over that stronger baseline is not required.
- If rank normalization helps on one dataset but not the other, report it as domain-dependent. No universal deployment claim.
- Bootstrap B=1000, master seed 20261029. All eligible genes are used, with no panel selection.
