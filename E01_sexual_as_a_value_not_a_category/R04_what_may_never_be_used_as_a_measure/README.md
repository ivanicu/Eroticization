# R04 — what may never be used as a measure

## DECISION: which variables in this release are unusable as outcomes or predictors, and which comparisons need what matching?

**Closed.** SAFE: self-reported induction is unusable as an outcome (no timing signature at >=1yr with 8.8x power; no structural signature; it tracks breadth). Group comparisons must be coverage-matched (corr with coverage gap +0.815) -- but NEVER matched when the question is about breadth itself, which is what R08 got wrong.

**11 sub-rounds**, each one belief update.

| r | directory |
|---|---|
| `r01` | [`r01_modality_invariance`](r01_modality_invariance) |
| `r02` | [`r02_deficit_reference_class`](r02_deficit_reference_class) |
| `r03` | [`r03_coverage_confound`](r03_coverage_confound) |
| `r04` | [`r04_matched_modality`](r04_matched_modality) |
| `r05` | [`r05_induction_timing`](r05_induction_timing) |
| `r06` | [`r06_timing_mde`](r06_timing_mde) |
| `r07` | [`r07_single_interest_power`](r07_single_interest_power) |
| `r08` | [`r08_off_manifold`](r08_off_manifold) |
| `r09` | [`r09_manifold_mde`](r09_manifold_mde) |
| `r10` | [`r10_matching_removed_it`](r10_matching_removed_it) |
| `r11` | [`r11_breadth_or_attribution`](r11_breadth_or_attribution) |

Numbers and intervals live in the top-level [`README.md`](../../README.md).
