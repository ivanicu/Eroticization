# Frozen lines

A line is frozen when further computation on it cannot identify what it is measuring. Freezing is
not a verdict that the line was wrong — several produced correct numbers. It is a verdict that
**the next round on this line would not change any decision**, because the ambiguity is in the
object rather than in the estimate. Each entry records what would unfreeze it, because a freeze
without an unfreeze condition is abandonment with better manners.

## 1. Direction of the consumption link — `r32`, `r48`, `r49`

**What it was trying to identify.** Whether consuming erotic content shapes which things get
eroticised, or an unusual set of coordinates drives consumption.

**Why it cannot separate.** The surviving effect is real and well-bounded (0.0439 after matching
coverage, breadth and sex; 3.6× the neuroticism reference). Both causal directions predict exactly
this number. No cross-sectional statistic distinguishes them, and three rounds spent on
increasingly careful matching improved the *estimate* without touching the *identification*.

**Unfreeze:** longitudinal measurement, or Ivan's phase 3 conditioning design. Nothing in this
release will do it.

## 2. Whether induction ever happens — `r34`–`r40`

**What it was trying to identify.** Whether the 82.7% who say porn induced new fetishes are
reporting an event.

**Why it cannot separate.** The timing signature is absent at ≥1 year with 8.8× power, and the
structural signature is bounded near zero from two directions. But absence of a *retrospective
self-report* signature is not absence of the phenomenon — it is a statement about the measurement.
Every remaining route through this release is another self-report.

**Unfreeze:** prospective onset measurement, or a conditioning experiment. Also: any release
carrying finer than 2-year onset bins would move the timing MDE from 1 year to months.

## 3. The A-vs-B separation by generalisation — `r07`–`r13`

**What it was trying to identify.** Whether a dedicated sexual-content system exists.

**Why it cannot separate.** Retraction #13: a dedicated module can be compositional, so success at
predicting from ordinary features is uninformative about A. These rounds measure compositionality
versus item-specificity, which is a different axis. The line is frozen *as an A/B test* and remains
live as a compositionality test.

**Unfreeze:** dissociation designs — selective interference, adaptation, cross-modal transfer with
a stronger modality contrast than this release carries.

## 4. Coordinate 3's identity — `r19`, `r20`

**What it was trying to identify.** The abjection ↔ blood/burning/weapons axis.

**Why it cannot separate.** It survives every attack (0.53 after deleting 16 of 32 contributing
blocks, 0.92 after deleting the fluid family) but its two poles are named from loadings, and the
option vocabulary is not a controlled feature space. Whether it is *disgust vs danger*, *interior
vs exterior*, or *contamination vs injury* cannot be decided from option text.

**Unfreeze:** a stimulus set built to vary those three contrasts orthogonally — which is a
collection problem, not an analysis problem.

## 5. Discriminant validity of the three role axes — `A02·R16`, `R17`, `R18`

**What it was trying to identify.** Whether POWER, GAZE and SUBSTANCE are three constructs or one
measured three ways — *without* the disattenuation that `ADVERSARY_FORECAST` #5 (p=0.40) correctly
flagged as fragile (GAZE's alpha is 0.163, and dividing by √0.163 is unstable exactly where it
matters).

**Why it cannot separate — and it took three rounds to find out.**

| round | what happened |
|---|---|
| `R16` | profile correlations against a 14-variable battery. Pre-registered kill fired ("101% of ceiling → one construct"). **Negative control returned 0.55–0.62 for a pure-noise axis**, so the verdict was inadmissible |
| `R17` | diagnosed it as battery intercorrelation pulling every profile toward the first PC; whitened the battery. **Gate failed again at 0.41–0.44.** So that diagnosis was wrong |
| `R18` | sampled the noise control 300 times instead of once. **Null is unbiased — mean −0.019 — with sd 0.356**, against a theoretical 1/√(14−3) = 0.302. The single draw R16/R17 relied on was the 90th percentile of that null |

**The instrument was never biased. It is unbiased and enormous.** A profile is 14 numbers, so the
correlation between two profiles has a 95% null band of **[−0.65, +0.66]** — and I twice read a
single draw from that band as a property of the design.

**What the design can still say, stated with its bound.** It detects profile similarity above
≈0.66 and nothing below:

| pair | profile r | p vs sampled null | verdict |
|---|---:|---:|---|
| POWER–SUBSTANCE | 0.934 | 0.000 | shares an external profile |
| GAZE–SUBSTANCE | 0.825 | 0.003 | shares an external profile |
| POWER–GAZE | 0.639 | 0.067 | indistinguishable from noise |

**And that is compatible with distinctness rather than against it.** Profile similarity is
necessary but not sufficient for identity: two orthogonal directions can share an external profile
by both loading on sex and age. POWER–SUBSTANCE sit at profile r 0.93 while their *direct*
correlation is 0.112 — which reads as **distinct axes sharing a demographic embedding**, not as one
construct. Deciding between those two readings is what this design cannot do.

**Unfreeze condition.** A battery of roughly **102 external variables** would bring the null sd
below 0.10; this release has **14**. Alternatively, abandon profile similarity and test
discriminance by incremental prediction — does each axis add held-out variance to a criterion the
other two already predict — which needs a criterion this release does not carry. Either way the
requirement is a different site, not more compute.
