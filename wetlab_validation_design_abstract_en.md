# Wet-lab Validation Design: Prospective Validation of the Gating Probe (One-page Summary)

**Status**: design document; funding and platform access approved (schedule to be confirmed with Prof. Wu). All effectiveness claims below are **hypotheses to be tested**, not results.

## What we are validating

We build gating probes that predict a virtual-cell model's own error from 9 geometric features of its outputs. In retrospective K562 benchmarks the probes reach AUROC 0.66–0.80, transfer across datasets, and can be calibrated on-site with as few as ~15 perturbations. **No wet-lab evidence exists yet.** This experiment asks: before running an experiment, can the gate flag which perturbation predictions will be wrong?

## Design in one paragraph

We build an OOD testbed from GSE281048 (MCF7 breast cancer cells, CRISPRi signaling-pathway perturbations, >1500 conditions across 5 pathways x 6 cell lines; pathway-targeted with 24h cytokine stimulation and no unstimulated arm — verified against the GEO record 2026-10-19; not present in mainstream pretraining corpora), run the models (GEARS primary; scFoundation secondary), calibrate the gate on ~15 held-out public perturbations, and select a final panel of 30–60 genes stratified into high-confidence / low-confidence / borderline tiers (10–20 each) plus non-targeting controls. The panel is **frozen** before any wet-lab work. The experiment runs the same panel in MCF7 (primary; MDA-MB-231 as a cross-cell-line exploration) with scRNA-seq readout (bulk as fallback), ≥3 biological replicates. Measured expression is then compared against model predictions, and errors are tested against the gate's pre-registered confidence tiers.

## Pre-registered hypotheses

- **H1**: gate-flagged low-confidence predictions are genuinely less accurate. Criterion: AUROC ≥ 0.65 (lower edge of the K562 range) with bootstrap 95% CI above 0.5.
- **H2**: confidence-tiered top-k selection improves hit rate. Criterion: median error of the high-confidence tier below the low-confidence tier, bootstrap CI excluding 0.
- **H1-null**: the gate is not merely a perturbation-strength proxy (partial correlation control). Failure outcomes are archived and reported with the same weight as successes; no post-hoc criterion changes.

## Practical notes

- Timeline: ~10 weeks from the freeze point (computational prep 3 weeks → pilot gating → perturb-seq window → analysis), plus platform queue time.
- Main risks: platform scheduling (mitigated by decoupling the freeze point from the lab window), CRISPRi efficiency, batch effects (tier-balanced batch assignment, frozen with the panel).
- Power honesty: with 10–20 perturbations per tier, the design resolves only large effects (Cohen's d ≳ 0.9); we report effect sizes with CIs rather than relying on small-sample p-values.

## Outputs

Results feed the prospective-validation section of the gating paper (supplement or main text depending on effect size; negative results reported regardless). A standalone short paper is considered only if H2 is strongly positive **and** replicated in the optional second wave (glioblastoma, GL261/U87 vs. Liu 2024 testbed).

**Pending decisions**: platform schedule & division of labor; final cell line; ethics/biosafety approval status; scRNA-seq vs bulk; go/no-go on the glioma wave.
