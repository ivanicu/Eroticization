# Retractions

Every claim this repository made and then killed, scoped or corrected — in the order it
happened, with what did the killing.

The rounds are numbered by when they ran, not by what survived. Read in that order the
repository looks like a sequence of findings. It is not. **Twelve of these fourteen entries
are a later round destroying an earlier round's conclusion, and in all twelve both rounds are
mine.** No outside challenger has run yet, so every row here is self-inflicted, and §
`ADVERSARY_FORECAST.md` records what I expect an outsider to add — written before one arrives,
so it scores my calibration rather than my luck.

This file exists because the git log has all of it and nobody reads a git log.

**Reading order matters.** Entries 1–4 are one failure mode: an instrument too weak for the
question, read as an answer. Entries 5–8 are a second: a name assigned before a stability check.
Entries 9–12 are a third and the most expensive — a group comparison that was partly measuring
how much of the survey each group had walked through. Entries 13–14 are mine about my own
framework, not about the data. Where a later entry supersedes an earlier one the earlier text
is **annotated, never rewritten** — a ledger that edits its own history is the thing it exists
to prevent.

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 1 | **"The existing kink databases cannot separate Ivan's models A/B/C; phase 1 is hypothesis generation only"** — my opening position, stated before touching the data | `r13` leave-one-block-out. Person factors fitted on 31 domains predict which options a person endorses in a domain the factors never saw: **32/32 blocks positive, median gain +0.0340, permuted floor −0.0029, 31/32 above 3× floor** | Nothing of the claim. The correct half was narrower: the database cannot separate A from B *by generalisation tests*, because a dedicated module can be compositional (see #13) |
| 2 | **"Cross-domain transfer of individual structure is weak"** (`r08`, median \|r\| = 0.064 over 418 block pairs) | `r09`. PC1-vs-PC1 is blind to a shared subspace that is rotated differently in each block. Held-out CCA on the same residuals: **0.272 against a 0.055 permutation floor, 285/321 pairs ≥2× floor** | Nothing. Same data, same centering, opposite conclusion — the instrument was the finding |
| 3 | **Six named "factor coordinates"** printed by the joint rank-8 factorisation (`r11`) | `r11` itself, on the line below them. Held-out cell R² **+0.0497** against a gating-preserving permutation null of **+0.0504**. It tied its own null | The factors are inadmissible. `r12` established this was *not* instrument blindness — the person side recovers sex at r=0.509 — so the metric, not the instrument, was unfit |
| 4 | **"Multilingual duplicate options contaminate block 1, the largest block (n=15,250), and its zero overlap with the English twins shows two disjoint language subsamples"** (`r13`) | Counting. **One respondent. Four singleton options.** Zero overlap is what n=1 always looks like | A real hygiene item: singleton options exist and enter every factorisation as noise columns. Floor-filtered at n≥20 from `r17` onward (5 dropped) |
| 5 | **Four named coordinates from SVD on the concatenated blocks** (`r18`) | `r18`'s own stability control. All four named block 1 as their largest contributor; drop-top-block stability **0.39 / 0.70 / 0.25 / 0.38**. Variance-weighted factorisation returns the biggest block's coordinates | The method correction: MAXVAR generalized CCA, where block size enters only through a projector's rank. Stability went to 0.78 / 0.78 / 0.93 / 0.96 / 0.88 |
| 6 | **GCCA coordinate 1** (`r19`), agreement 0.299 | Its own per-block person-shuffle null: **0.302**. Observed sits *below* its null | Nothing. Reported dead rather than as a factor |
| 7 | **GCCA coordinate 5** (`r20`) | Deleting the 12 blocks contributing its top-25 loadings: stability **0.25** | Dropped. Coordinates 2–4 survive the same attack at 0.63 / 0.53 / 0.58 and survive deleting the whole 6-block fluid family at 0.94 / 0.92 / 0.73 |
| 8 | **"Schedule adherence is an individual trait"**, with correlates (biomale +0.093, own mean onset +0.100) — `r29` | Split-half over categories: r=+0.120, **Spearman–Brown 0.214**. Mostly measurement noise | The shared schedule itself: per-person rank agreement with the population ordering mean +0.232, **74.7% of 9,691 people positive against a 48.9% null, Cohen d 0.69**. The *adherence* correlates are correlates of nothing and were withdrawn before publication |
| 9 | **"Unit acquisition of a coordinate"** — my own hypothesis from `r26`, that categories sharing a GCCA coordinate were acquired together | `r27`. With preference similarity controlled: **COORD ΔR² +0.0008, t = −0.46**; TEMPO ΔR² **+0.0726, t = +4.59**. The two predictors correlate −0.028, so COORD is absent rather than hidden | Developmental tempo, which then survived its own binning null (observed t +4.19 vs null mean +0.50 ± 0.75, max +2.49) |
| 10 | **"The erotic grammar is largely modality-invariant"** (`r30`, deficit 0.0360 vs a 0.0965 sex deficit) | `r31`'s reference class. Nine person-variable splits measured identically: modality ranks **6th of 9**, median deficit 0.0313. It is typical, not exceptional | The animation split (drawn vs live-action, no real bodies) at 0.0123 raw / 0.0204 matched, 0.8 sd — bounded below sex-sized, which is the constraint on A that survives |
| 11 | **"How much porn you consume reorganises erotic grammar 2.4× more than sex does"** (`r31`, deficit **0.2285**, z=60) | `r32`: **corr(deficit, between-group gap in blocks entered) = +0.815** across the nine splits. Block-count-matched, pornhabit falls to **0.0871** — 38% survives. Then `r49` matched theta and sex as well: **0.0439** | The ordering **reversed**. Sex measured comparably is 0.0778. Anything published off the first pass would have said the opposite of what the data says |
| 12 | **"Self-reported porn-induced fetish acquisition has no structural signature"** (`r37`–`r38`, misfit rho +0.038, extremity +0.020, both null with stated MDEs) | `r39`, my own self-attack. Both measures live in within-block profile *shape*, and I had block-count matched — deleting the only place the effect could live. Unmatched: **n_categories 17.02 → 27.80 of 68, rho +0.2410, p 1.7e-177** | The null is correct **about shape**, with bounds (added misfit < 0.1 sd; extremity < 10%). It was a null about the wrong quantity, and the scope now says so |
| 13 | **"Model A predicts ordinary semantic features cannot predict what is sexual"** — my framing of Ivan's model A, used to rank the whole phase-1 design | Faces. A dedicated system (FFA, double dissociation) computes from ordinary image features and generalises perfectly to novel faces. Dedicated ⇏ non-compositional. The test is sound in one direction only: *failure to predict* supports A; *success* does not refute it | The reframe: rounds r07–r13 measure **compositionality vs item-specificity**, not content vs value. A/B separate on dissociation, acquisition and cross-modal transfer instead |
| 14 | **A norm/transgression parameter `n_i`** added to Ivan's operator, predicting that rarity-preference is a domain-general trait | `r09`. Cross-domain transfer of count-controlled rarity taste: **rho +0.065, 26 of 321 pairs above 0.2**. Not a domain-general trait as operationalised | The prediction that a normalisation manipulation should split subgroups by *sign* is untested and remains the cheapest available separator — it is in `PREREGISTRATION.md`, not here |

---

## What the pattern is

Three failure modes produced twelve of fourteen, and they are not independent:

**An instrument too weak, read as an answer** (#1, #2, #3, #6). In every case a null or a small
number was produced by a measurement that could not have shown the effect, and in three of four
the positive control that would have caught it existed and was run *afterwards*. `r12` is the
clearest: the person side of a factorisation recovers sex at r=0.509 while its held-out cell
metric ties its own null, because item covariance saturates that metric. **The instrument was
live and the metric was blind, and those are different verdicts.**

**A name assigned before a stability check** (#5, #7, #8). Naming is where the temptation is,
because a named coordinate reads as understanding. Every naming attempt in this repository that
was not immediately followed by a deletion test turned out to be reporting the largest block's
structure or an unreliable composite.

**A group comparison measuring survey coverage** (#10, #11, and it reached back into #12).
The release is a gated tree, so groups differ in how much of it they walked through, and
subspace congruence falls when coverage differs for reasons that have nothing to do with the
construct. This is undocumented anywhere I have read about this dataset and the correlation is
+0.815. It produced the single largest overstatement here — a factor of 5.2, and a reversed
ordering.

**And the fourth, which has only one entry and is the most dangerous** (#12): matching on the
confound deleted the effect. Match for a *shape* question, never for a *breadth* question, and
say which is being asked before matching. The same procedure is a control in one round and a
destroyer of evidence in the next.

---

## Entries 15–16, added by `R13·r01` — prior art, which the forecast predicted at p=0.75

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 15 | **"Any group comparison on this release must be block-count matched or it partly measures survey coverage. Nothing I read about this dataset says so."** — `R06`, and the sentence was in `README.md` as the transferable methodological result | The explorer's own `analysis/swarm/14-missingness.md`, dated **2026-02-13**, five months earlier: *"if missingness differs systematically by gender, orientation, or politics, then group comparisons on gated columns are comparing different subpopulations, not the full sample"*, with per-column missingness by gender, orientation and politics, responder-vs-non-responder Cohen's d > 0.7 on 17 of 18 columns, and *"any 'gender difference' in these kinks is partially an artifact of who answered"* | The **hazard is theirs**. What is mine is narrower and still stands: `corr(subspace-congruence deficit, between-group coverage gap) = +0.815` over nine splits, and the matched correction that took pornhabit from 0.2285 to 0.0871 to 0.0439 and reversed its ordering against sex. A quantification of a documented hazard is not a discovery of it, and the sentence claiming nothing documented it is simply false |
| 16 | **"Onset is a proxy for intensity"** presented as a confound I found and controlled (`R04`, corr(residual onset, residual preference) = −0.126) | `analysis/swarm/01-age-onset.md` **Finding 1**, same date: *"Earlier onset predicts higher current intensity"*, r = −0.12 to −0.17 across nine kinks, n up to 8,741 | The control, not the number. Stripping it and finding RSA 0.621 → 0.600 is still the right move and still holds. But I reported discovering a relationship that is the first published finding in the release's own analysis directory |

**What this pair actually costs.** Not much numerically — no interval moves. It costs the claim to
novelty on the one methodological result I had been treating as the project's most transferable
output, and it reproduces the CoVal failure exactly: **fifty-two sub-rounds against an object
whose own analysis directory was one `git ls-tree` away.** The forecast in
`ADVERSARY_FORECAST.md` put this at p=0.75 and named the coverage artifact as a likely casualty,
which is the one piece of calibration this file can report in its own favour.

**And the honest tail:** four of the fifteen published analyses bear on live claims and have not
been read line by line. Those claims are `[unchecked]`, which is not `[clean]`.

---

## Entry 17, added by `E01·R06·r02` — the one where the other analysis was right and I was wrong

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 17 | **"θ is nearly domain-encapsulated. All 15 non-sexual variables jointly R² = 0.012. Childhood adversity, adult sexual assault, corporal punishment and sexual repressiveness of upbringing are all bounded under 0.09 disattenuated — the family of theories rooting broad or unusual sexuality in adversity gets essentially no support here."** (`R05·r06`, and it was stated at a sample size where 0.023 would have been visible, which made it sound decisive) | A full `realstat` round built to adjudicate it against `12-multivariate` Test 7, which reported the opposite (**Any 10.665 vs None 9.750, d=0.151, survives demographic controls, verdict REAL**). Two errors on my side, and they compound: **(a)** my breadth was a count over the 68 gated rating columns, so survey progression sits between adversity and the outcome — but `totalfetishcategory` is asked of **everyone, 0% missing**, and on an ungated outcome progression is not a mediator to control, it is over-control. **(b)** The 24-cell specification grid shows exactly this: every cell with progression controlled is near zero (0.021 / 0.018 / −0.008), every cell without it is not. Headline cell — ungated outcome, demographics + acquiescence controlled, Pearson: **r = +0.0590**, effect/floor **7.5**, sham 0.006, placebo 0.028, positive control recovers biomale→(receivepain−givepain) at r = −0.407. Pre-registered kill fired: **world B refuted** | The acquiescence correction, which is mine and still stands (0.070 → 0.059 unadjusted → adjusted). And the narrower true statement: adversity's association with erotic breadth is **real and small**, r ≈ 0.06, about 0.35% of variance — not "essentially no support", and not large either. My error was reporting a specification-dependent near-zero as a property of the world |

**Why this one is different from 1–16.** Every earlier entry is a later round of mine killing an
earlier round of mine. This is the first where **someone else's analysis was right and mine was
wrong**, and I found it only because `R06·r01` went looking for prior art after 52 sub-rounds. The
disagreement was sitting in a public file the whole time, and my own specification grid already
contained the answer — 0.016 / 0.050 / 0.088 across three specs — which I reported and then
summarised by picking the one that agreed with the conclusion I had already written.

**That is the failure mode with no entry above it: publishing the specification curve and then
narrating the cell you liked.** The grid is not a robustness appendix. It is the result.

---

## Entry 18, added by `E01·A02·R14`–`R15` — the count was a property of how many I asked for

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 18 | **"Four coordinates survive a block split-half (7 above 3× floor); practically the grammar is ~4-dimensional"** — reported with held-out canonical r .357 / .290 / .193 / .130 and carried as a headline through nine subsequent rounds | A specification curve I had pre-registered against myself. **`R14`, 60 cells (option floor 5/10/20/40/80 × K 4/6/8/10 × 3 seeds):** the option floor is *completely inert* — every row identical, 32 blocks entering at every floor — while the count tracks **K** almost 1:1 (K=4→4, 6→5, 8→7, 10→8). **`R15`, sweeping K to 24:** count = 10 at K=12 and **16 at K=24**, +6, non-saturating. The pre-registered kill fired. Positive control held throughout — the permutation floor stayed flat at 0.008–0.011 across all cells, so this is a moving count against a fixed floor, not a moving floor | **The transfer, not its dimensionality.** Every retained coordinate sits 30–50× its own permutation floor, and the profile falls steeply (0.461 · 0.450 · 0.321 · 0.311 · 0.256 · 0.228 · 0.201 · 0.159 · 0.125 · 0.104 …). Cross-domain structure is real and large. **How much of it there is is not identified by this procedure**, and "~4-dimensional" is withdrawn as a number rather than replaced by a bigger one |

**And the replacement failed in the same round, which is the part worth keeping.** `R15`
pre-registered the magnitude profile as the K-free object that would survive if the count did not.
It is not K-free: the top canonical correlation itself rises from **0.304 at K=4 to 0.461 at K=24**,
because CCA with more components finds better-aligned combinations of the same halves. So the
procedure is monotone in K *everywhere* — in the count and in the magnitudes — and neither is
identified without fixing K on grounds outside the data. **I proposed the fallback in the same
docstring as the test, and the test killed both.** That is what a pre-registered replacement is
for; had I chosen it afterwards it would have looked like a finding.

**What this does NOT touch**, stated because a retraction that quietly widens is its own failure:
the three role axes (`A02`) come from *named indicators* — POWER, GAZE, SUBSTANCE — not from the
K-dependent CCA, and their mutual disattenuated r ≤ 0.362 at fixed n stands. The leave-one-block-out
result (32/32 blocks, median gain +0.0340) is a prediction into held-out domains and carries no K.

**Forecast scoring**, since `ADVERSARY_FORECAST.md` exists to be graded rather than admired.
Prediction #2, written before any of this: *"The four-coordinate count is unstable to the option
floor. `n≥20` was chosen once and never swept. A challenger sweeping it finds 3 or 6."* — **p=0.65,
direction CORRECT, mechanism WRONG.** The count is unstable, and the option floor has nothing to do
with it. I predicted the right casualty for the wrong reason, which is worth less than it looks:
had I swept only what I predicted, I would have found the floor inert and concluded the count was
robust.

---

## Entry 19, added by `E01·A07·R01` — the one K I picked was the only K where the ordering came out the way I reported it

`R14`–`R15` proved the CCA machinery is monotone in K. Three surviving headlines used the same
machinery at a fixed, never-swept **K=5**. Swept over K ∈ {2,3,5,8,12,16} × 3 seeds, with a
**known-null pool** added for the first time — a coin-flip sham and a row-parity placebo, 36 cells,
p95 = 0.0216 — against a measured ceiling sd of 0.0133.

| split | reported | median over K | range | vs null p95 | effect/floor |
|---|---:|---:|---|---|---:|
| sex | 0.0927 | **+0.1437** | .081–.206 | ABOVE | 10.8 |
| modality written/visual | 0.0334 | **+0.0546** | .017–.157 | ABOVE | 4.1 |
| consumption (pornhabit) | 0.0439 | **+0.0398** | −.012–.103 | ABOVE | 3.0 |
| neuroticism *(my reference)* | 0.0121 | **+0.0126** | −.020–.054 | **INSIDE THE NULL** | 1.0 |

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 19 | **"Consumption reorganises the grammar at 0.0439, 3.6× the neuroticism reference"** and **"modality is typical, ranking 6th of 9, matched deficit 0.0334"** — reported two rounds apart and used together to argue consumption matters and modality does not | The K sweep. **K=5 is the only value at which consumption outranks modality**; at K ∈ {2,3,8,12,16} the order is sex > modality > consumption. Decomposed: modality at K=5 is **0.0326 against 0.0670 at every other K — ratio 0.49**, so my K halved it; consumption at K=5 is **0.0463 against 0.0372 — ratio 1.24**, so my K inflated it. The single arbitrary choice moved the two claims in opposite directions, which is exactly what it takes to flip an ordering | Both effects, restated and reversed. **Modality is the larger of the two** (0.0546 vs 0.0398) and both clear the known-null p95 with effect/floor 4.1 and 3.0. Sex is robustly largest at every K, effect/floor 10.8 — the positive control held throughout |

**The reference was noise.** I validated consumption by saying it was "3.6× the neuroticism
reference measured the same way". Neuroticism's median deficit is **+0.0126, inside the known-null
pool, effect/floor 1.0.** So the argument was *3.6× nothing*. Consumption clears the null on its
own, so the conclusion survives — but the sentence that made it sound established was comparing
against a quantity I had never checked was non-null.

**And the reason this was findable only now: I had never run a placebo.** Fifty-plus rounds of
congruence comparisons, a nine-variable "reference class" in `A04` — and every single one of those
nine was a *real* variable. Not one was known-null. A reference class of real variables tells you
where a new variable ranks; it cannot tell you whether the bottom of the class is above zero.
**The first time I split on row-index parity, the bottom of my reference class turned out to be
sitting in it.**

**Correction made in flight, recorded because the near-miss is the point:** on first reading the
rank table I concluded modality was *inside* the placebo range, because the placebo's max (0.0374)
exceeds modality's K=5 value (0.0326). That comparison is wrong — a max against a single cell.
Against medians and the null p95, modality is comfortably above. I nearly retracted a real finding
by attacking it with a cheaper statistic than the one that produced it, which `realstat` §3 names
as the most expensive kind of error.
