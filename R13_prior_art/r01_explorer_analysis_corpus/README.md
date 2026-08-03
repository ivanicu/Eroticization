# R13·r01 — the explorer's own analysis corpus

Pre-registered in [`PREREGISTRATION.md`](../../PREREGISTRATION.md) as the round that runs before
any other, with a consequence fixed in advance: **any headline coded VERIFICATION is relabelled,
no exceptions, including the three I most want to keep.**

## Why it was pre-registered first

The CoVal programme spent 140 rounds discovering an algorithm documented in its own dataset card.
This release has a public explorer whose source repository ships an `analysis/` directory. **Fifty-two
sub-rounds ran without opening it.**

## Corpus

`github.com/austeane/aella-survey-site`, `analysis/swarm/` — 15 analyses dated 2026-02-13, five
months before this project started: `01-age-onset` · `02-mental-health` · `03-personality-kinks` ·
`04-politics-deep` · `05-gender-deep` · `06-taboo-clusters` · `07-relationships` ·
`08-orientation-identity` · `09-porn-media` · `10-surprises` · `11-bootstrap-cis` ·
`12-multivariate` · `13-interactions` · `14-missingness` · `15-clustering`. Plus
`META-FINDINGS.md`, `findings.json`, and `docs/schema/interesting-findings.md`. Also
`github.com/darknet-doll/BKS_survey`, which ships OCEAN scoring with acquiescence-bias notes.

## Coding of the eleven headlines

| headline | code | evidence |
|---|---|---|
| group comparisons must be coverage-matched | **VERIFICATION** | `14-missingness`: *"if missingness differs systematically by gender, orientation, or politics, then group comparisons on gated columns are comparing different subpopulations"* and *"any 'gender difference' in these kinks is partially an artifact of who answered"* |
| earlier onset ↔ higher intensity (my "leakage confound") | **VERIFICATION** | `01-age-onset` Finding 1, r = −0.12 to −0.17 per kink, n up to 8,741. Mine: −0.126 |
| breadth is part acquiescence | **VERIFICATION** | `BKS_Data_Review.md` ships acquiescence-bias notes for the OCEAN columns |
| the gated tree, P(enter\|parent>0)=0.99 | **VERIFICATION** | `14-missingness` §1, per-column missingness table |
| four coordinates / latent structure | **contested** | `15-clustering` runs k-means latent profiles, k=4, silhouette 0.17–0.18 — clusters of *people*, not a shared coordinate subspace tested for cross-domain transfer. Adjacent, different object |
| remaining seven | **unchecked** | `12-multivariate`, `13-interactions`, `09-porn-media`, `06-taboo-clusters` not yet read line by line |

## Verdict

Two headlines and two supporting facts are prior art. The `+0.815` correlation and the matching
correction are a novel *quantification* of a hazard that was already documented — the hazard is
theirs, the measurement of its effect on subspace congruence is mine. That distinction does not
rescue the sentence I actually wrote.

**Unchecked is not clean.** Four of the fifteen analyses bear on live claims and have not been read.
