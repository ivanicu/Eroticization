# R04 — what may never be used as a measure

## DECISION: which variables in this release are unusable as outcomes or predictors, and which comparisons need what matching?

**Closed.** SAFE: self-reported induction is unusable as an outcome (no timing signature at >=1yr with 8.8x power; no structural signature; it tracks breadth). Group comparisons must be coverage-matched (corr with coverage gap +0.815) -- but NEVER matched when the question is about breadth itself, which is what R08 got wrong.

**11 sub-rounds**, each one belief update.

| R | directory |
|---|---|
| `r01` | [`R01_modality_invariance`](R01_modality_invariance) |
| `r02` | [`R02_deficit_reference_class`](R02_deficit_reference_class) |
| `r03` | [`R03_coverage_confound`](R03_coverage_confound) |
| `r04` | [`R04_matched_modality`](R04_matched_modality) |
| `r05` | [`R05_induction_timing`](R05_induction_timing) |
| `r06` | [`R06_timing_mde`](R06_timing_mde) |
| `r07` | [`R07_single_interest_power`](R07_single_interest_power) |
| `r08` | [`R08_off_manifold`](R08_off_manifold) |
| `r09` | [`R09_manifold_mde`](R09_manifold_mde) |
| `r10` | [`R10_matching_removed_it`](R10_matching_removed_it) |
| `r11` | [`R11_breadth_or_attribution`](R11_breadth_or_attribution) |

Numbers and intervals live in the top-level [`README.md`](../../README.md).
