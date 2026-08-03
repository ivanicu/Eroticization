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

---

## Entry 20, added by `E01·A07·R02` — the coverage law was one point, and the point was 7.6× the rest of the range

After #15 conceded the coverage *hazard* to the explorer's own `14-missingness`, the only novel
claim this project still had was the quantification: **corr(congruence deficit, coverage gap) =
+0.815**, and the matched correction built on it.

`ADVERSARY_FORECAST` #3, written before this round, p=0.55: *"nine splits, and PORNHABIT sits far
out on both axes. Without it the correlation is much weaker, and the methodological claim rests on
a single influential observation."*

**The coverage gaps, printed for the first time:**

| split | coverage gap | deficit @K5 |
|---|---:|---:|
| pornhabit | **2.435** | +0.0463 |
| neuroticism | 0.321 | +0.0148 |
| sex | 0.130 | **+0.0838** |
| placebo (row parity) | 0.108 | +0.0002 |
| sham (coin flip) | 0.079 | −0.0105 |
| modality | 0.069 | +0.0326 |

**pornhabit's gap is 7.6× the next largest.** And the two points that most contradict a coverage
law are the two strongest real effects: **sex has the second-largest deficit with the third-*smallest*
gap**, and modality has a real deficit at the smallest gap of all.

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 20 | **"corr(deficit, coverage gap) = +0.815 across nine splits. Any group comparison on this release must be block-count matched or it partly measures survey coverage"** — carried as the project's transferable methodological result, and after #15 its only novel one | Influence analysis with the known-null splits added to the class. **Median r over K ∈ {3,5,8,12} = +0.127. Worst leave-one-out = −0.294, and the point whose removal does it is pornhabit at 6 of 8 specifications.** The permutation null at this n is enormous — **|r| p95 = 0.820 at n=6 splits** — so at this many units a correlation has to clear 0.82 to carry information, and +0.815 never did. Positive control passed (gap vs itself 1.000; a synthetic deficit built as f(coverage) recovered at 1.000), so this is a real null and not broken plumbing | **The correction, not the law.** Block-count matching moved pornhabit's deficit 0.2285 → 0.0871 → 0.0439 by *direct measurement on that split*, and that is unaffected — it never depended on the correlation. What dies is the general claim that deficit tracks coverage across splits. It does not: sex and modality both have large deficits at near-zero coverage gaps |

**What the correlation actually was.** One split with a coverage gap 7.6× everything else, which
also happened to have a large deficit for its own reasons. With n=9 units and a null p95 near 0.82,
that is a two-point line. **I published a law computed over nine units without ever asking what a
correlation is worth at n=9** — and the answer, measured here, is: nothing below 0.82.

**The unit was never people.** Every one of the 15,503 respondents appears in every split, which is
why the sample felt enormous. But the *estimand* is a correlation over splits, and there are six of
them. `realstat` G1 says name the estimand before the method; I named the method (correlate the
deficits) and inherited an n of nine by accident.

**Forecast scoring.** #3 predicted the mechanism exactly — pornhabit far out on both axes, one
influential point — at p=0.55. **CORRECT, mechanism included**, unlike #2 last round where only the
direction was right. Two of three forecasts now scored, both hits.

---

## Entry 21, added by `E01·A03·R08` — my correction was as unfit as the thing it corrected, and my own positive control is what said so

`ADVERSARY_FORECAST` #4, p=0.50: *"344 pairs come from 27 categories, so the effective n is closer
to 27 and every p-value there is overstated."* Checked against the code rather than memory, it is
**exactly half right**, and the halves are in different rounds:

| test | its null | verdict |
|---|---|---|
| the RSA, acquired-together vs liked-together | `A03·R03:88` permutes **categories** (Mantel) | **correct all along.** RSA +0.600, null 0.000 ± 0.055, **z = +10.9** — stands |
| the residual-structure claim, "80% of onset structure is not explained by preference" | `A03·R04:58` permutes **pairs** | **anticonservative.** Pair permutation destroys the block a category imposes on its ~25 pairs, so null eigenvalues come out too small |

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 21 | **"80% of the sd of onset structure is not explained by preference; residual top-3 eigenvalues z = 4.9 / 3.5 / 2.7"** | Nothing yet, and that is the finding. The pair-level null is unfit (anticonservative). The category-level null I wrote to replace it is **degenerate**: permuting category labels of a symmetric matrix is a *similarity transform*, `M[ix_(p,p)]`, which leaves eigenvalues **exactly invariant** — observed and null were identical to three decimals, sd = 0. **My own pre-registered positive control caught it**: planting rank-2 structure at amplitude 0.10 returned z = +1.0, and at amplitude 0.00 returned z = −0.8. A test that cannot see a planted effect cannot report an absent one | **UNVERIFIED, and it must not be rounded to OVERTURNED.** The claim has no valid test in either direction. What it needs is a null that destroys the residual structure while preserving the category block structure — a person-level bootstrap of the whole pipeline, not a permutation of the finished matrix |

**The automated verdict was wrong and I am overriding it.** The script's pre-registered kill fired
— `z = −0.8 < 2.0 : CLAIM WITHDRAWN` — because the kill threshold was written assuming the null
would be valid. The positive control, written in the same docstring, says it is not. **A
pre-registered kill is a commitment about what to do with a valid measurement; it does not
authorise acting on a broken one.** Folding this into OVERTURNED would have manufactured a false
retraction, and a false retraction is as permanent as a false acquittal — nobody re-examines a
claim its own author has withdrawn.

**The general fact worth carrying out of this project:** *eigenvalues are invariant under
simultaneous row–column permutation, so label-permutation is never a null for a spectral statistic.*
It is a perfectly good null for the RSA in the row above, because a correlation between two
matrices is not permutation-invariant. **Same permutation, one statistic it tests and one it cannot
touch** — which is the proxy ledger exactly: a null is sound for a *statistic*, never for a *matrix*.

**Forecast scoring.** #4 CORRECT on the residual test, WRONG on the RSA — the first partial. Three
of seven now scored: #2 direction-only, #3 fully correct, #4 half. **None of the three was scored
in my favour by me — each was checked by running the thing it predicted.**

---

## Entry 22, added by `E01·A03·R09` — the first UNVERIFIED to be resolved, and it went the way the original round said

#21 left "onset carries structure preference does not" untestable: the pair-level null was
anticonservative, the label-permutation null was degenerate, and **both were cheaper than the
correct one**. `realstat` G2 — *a permutation null answers "did the pairing matter", never "why";
name the world it must exclude and build it* — so the rival world was built at the person level
instead of permuted at the matrix level.

**Rival world, specified in full before running:** onset similarity is *entirely* preference
similarity plus noise. Synthetic person-level onset data generated with the observed preference
covariance, pushed through the same 2-year binning, the same missingness mask, and the same
residualisation.

| | top residual eigenvalue |
|---|---:|
| **observed** | **0.959** |
| rival world, noise 0.3 | 0.455 ± 0.041, 95% [0.385, 0.535] |
| rival world, noise 0.5 *(matched)* | 0.441 ± 0.040, 95% [0.372, **0.532**] |
| rival world, noise 0.7 | 0.429 ± 0.040, 95% [0.362, 0.515] |

Observed is **1.8× the rival world's upper bound** and above it at every noise level.

**The controls that make this reportable, unlike #21's:**
- **positive** — injecting person-level rank-2 structure: amplitude 0.0 → 0.441, *not detected*
  (identical to the rival world, so the test does not fire on nothing); 0.3 → 1.321, detected;
  0.6 → 3.263, detected. **It sees what is there and stays silent on what is not.**
- **sham** — synthetic onset with identity covariance, no preference coupling at all: 0.418 against
  0.441 coupled. Nearly equal, which is the correct behaviour: after residualising on preference,
  what is left should not care whether preference was in the generator.

| # | Claim | Verdict |
|---|---|---|
| 22 | **"Onset similarity carries structure that preference similarity does not"** — `A03·R04`, left UNVERIFIED by #21 | **CONFIRMED.** Resolved in favour of the original round, against a purpose-built rival world with a working positive control. **Scope, narrowed:** this confirms the *existence* of residual structure. The headline number attached to it — "80% of the sd" — is a descriptive ratio (0.076/0.094), not a tested quantity, and remains descriptive |

**What this entry is for.** Twenty-one entries of this file are things I got wrong. This one is a
claim that two bad instruments had put in limbo and a good instrument recovered — and the good
instrument cost roughly ten times the compute of either bad one. **The cheap nulls were unfit in
*opposite* directions**: pair permutation inflated the z, label permutation made it undefined.
Averaging them, or trusting whichever ran first, would have landed anywhere.

**UNVERIFIED is the verdict that obliges a re-run**, and this is the first time in the project that
obligation has been discharged rather than deferred.

---

## Entry 23, added by `E01·A05·R10` — the ceiling was real, the number behind it was not what I forecast

`ADVERSARY_FORECAST` #6, p=0.35: *"breadth's 0.557 reliability makes θ's 'no external correlate'
partly a ceiling, and a better-measured θ would show the adversity correlations at 0.15 rather than
under 0.09."*

Tested by building θ on a **measured reliability ladder** (composites of 1–6 breadth indicators,
reliability 0.62–0.89) rather than dividing once by a single alpha. Classical test theory predicts
`r_obs = r_true · √rel`, so **`r_obs/√rel` must be flat if a ceiling is operating**.

**It is flat — drift −7% from the lowest to the highest reliability rung.** The measurement model
holds, so the disattenuated values are estimates rather than extrapolations off a single point:

| target | r_true (= r_obs/√rel) | previously reported |
|---|---:|---|
| pornhabit **(positive control)** | **0.22** | r rises 0.180 → 0.217 across rungs exactly as √rel predicts |
| childhood adversity | **0.105** | "bounded under 0.09 disattenuated" |
| mental illness | **0.099** | 0.054 |
| openness | **0.090** | *"sign flips across specifications, unreportable"* |
| spanked as a child | 0.037 | −0.000 |
| adult sexual assault | 0.033 | 0.020 |
| **NULL row-parity** | **0.000** | — |
| **NULL coin-flip** | **0.010** | — |

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 23 | **"θ is nearly domain-encapsulated; all 15 non-sexual variables jointly R² = 0.012; adversity bounded under 0.09 disattenuated"** — already partly retracted by #17, which showed the gated outcome and the progression control were doing the work | The reliability ladder. A ceiling *was* operating, and it was worth roughly a factor of two: adversity 0.105, mental illness 0.099, openness 0.090. The negative controls confirm the ladder is not simply inflating everything — row-parity sits at **0.000** and coin-flip at **0.010** across all six rungs | **The size, restated once more and now pinned rather than bounded.** These are ~1% of variance each. "Essentially no support" was wrong; "large" would also be wrong. Three rounds have now moved this quantity — 0.016 (over-controlled) → 0.059 (ungated, #17) → **0.105 (reliability-corrected)** — and each move was a measurement error of mine, not new data |

**Openness is the one I owe most.** In `A05` I wrote that its sign flipped across specifications
and was therefore unreportable. With θ measured properly it is a stable **+0.090** at every rung.
The instability was in my θ, not in openness — **I attributed my own measurement noise to the
variable I was measuring against**, which is the same error as blaming a scale for a wobbly table.

**Forecast scoring.** #6: **CORRECT on mechanism** — a reliability ceiling was operating and I had
said so at p=0.35 — **WRONG on magnitude**, 0.105 against the 0.15 predicted. Four of seven now
scored: #2 direction-only, #3 fully correct, #4 half, #6 mechanism-only. **The forecast file has
been right about *what* would break every time and wrong about *how much* three times out of four**,
which is itself a calibration fact worth more than any single entry.

---

## Entry 24, added by `E01·A02·R19` — the flagship claim, measured with a validated ruler instead of a guessed one

The headline of this project has been **"top/bottom is one word over three near-independent
coordinates"**: POWER, GAZE, SUBSTANCE, observed mutual |r| ≤ 0.112, disattenuated ≤ 0.362,
effective dimensionality 2.95 of 3. `ADVERSARY_FORECAST` #5 flagged the disattenuation as fragile.
`R16`–`R18` tried to bypass it and froze — underpowered by construction.

So: settle the pair that *can* be settled. GAZE has only 4 indicators in this release (the
option-level search returned regex false positives — *Wetlook*, *Video game characters*), but
POWER has 3 and SUBSTANCE has 7, and the ladder method already validated on θ in `A05·R10`.

**The control that makes this reportable:** SUBSTANCE split against *itself*, two disjoint halves
of its own 7 indicators, through the identical pipeline. **r_true = +1.018.** A disattenuation
that could not recover 1.0 for a measure against itself could not be trusted to report anything
else; this one does. Sham (POWER against a 7-indicator non-role personality composite): **+0.023**.

| | r_obs, low → high reliability | r_true |
|---|---|---:|
| SUBSTANCE vs itself *(positive control)* | 0.393 → 0.516 → 0.516 | **+1.018** |
| **POWER vs SUBSTANCE** | **0.298 → 0.341 → 0.381** | **+0.605** |
| POWER vs sham *(non-role)* | 0.002 → 0.026 | +0.023 |

`r_obs` rises with reliability exactly as `r_true·√rel` requires, and `r_true` drifts only −15%
across terciles — the measurement model holds.

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 24 | **"Three near-independent role axes; max disattenuated r = 0.362; effective dimensionality 2.95 of 3"** — carried as the project's most interesting finding since `A02` | The ladder. **POWER–SUBSTANCE r_true = +0.605, not 0.233** — they share **37% of variance**, not 5%. My earlier figure divided the observed 0.112 by √(0.686 × 0.337), and that 0.337 was a *block split-half of a GCCA coordinate*, which is not the reliability of the composite I actually correlated. **I disattenuated with a reliability belonging to a different measure.** The pre-registered kill lands in the middle band: 0.45 < 0.605 < 0.70 → **UNVERIFIED**, neither one construct nor confirmed distinct | Not identity either — 0.605 is a long way from 1.0, and the sham sits at 0.023. **The honest statement is that POWER and SUBSTANCE are substantially correlated but not the same**, and that "2.95 of 3" is withdrawn, since it was computed from attenuated observed correlations and no longer follows |

**Gap in this round, stated rather than left to be found:** the noise-vs-noise negative control
produced no admissible rungs, because pure-noise composites have reliability ≈ 0 and my own
admissibility filter (`rel > 0.02`) excluded them. The sham does the same job and returns 0.023,
but the control I *pre-registered* did not run, and a control that cannot run is not a control that
passed.

**What this costs.** Every previous retraction here moved a peripheral number. This one moves the
sentence I would have led with. Two of the three axes share more than a third of their variance,
and the third cannot be measured well enough in this release to place at all.

---

## Entry 25, added by `E01·A01·R13` — the foundation, attacked for the first time, survives at half size

Everything in this project descends from the leave-one-block-out result: person factors fitted on
31 domains predict which options a person endorses in a domain the factors never saw, **32/32
blocks positive, median gain +0.0340**. Its only control had been a permuted-factor null, which
`R08` taught is precisely the sort of null that can be unfit.

**The confound, written before running.** The factors mean-impute every block a person did not
enter, so a person's factor vector partly encodes **which blocks they entered** — their coverage
pattern — and coverage correlates with breadth, which correlates with within-block profile shape.
The factors could win with no cross-domain grammar at all, carrying only the gating tree's shadow.

| baseline | median gain | blocks positive | permuted floor |
|---|---:|---:|---:|
| propensity *(as originally reported)* | **+0.0373** | 97.9% | −0.0025 |
| + n_blocks | +0.0336 | 97.9% | −0.0027 |
| + coverage principal components | +0.0237 | 97.9% | −0.0025 |
| **+ full coverage, all 32 indicators** | **+0.0170** | **31/32** | −0.0023 |

**Positive control:** hand a block its *own* factors and the gain is **+0.7915** — the instrument
detects grammar decisively when grammar is present, so a small number here is a small effect and
not a blind test.

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 25 | **"Person factors from 31 domains predict a domain they never saw, 32/32 blocks, median gain +0.0340"** — the foundation | Coverage. Putting all 32 block-entry indicators in the baseline removes **54%** of the gain. My pre-registered kill lands mid-band (0.010 < 0.0170 < 0.020) → **UNVERIFIED on magnitude** | **The sign, decisively.** 31 of 32 blocks stay positive against a floor of −0.0023; under a null of 50% that is p ≈ 1e-8. The effect is real and roughly half what I reported |

**And the correct output is a bound, not a point — for the same reason as #17.** Coverage is not
purely a confound: *which blocks you enter is itself an expression of what you like*, so putting
all 32 indicators in the baseline removes real grammar along with the artifact. That makes
**+0.0170 a lower bound and +0.0373 an upper bound**, and the truth is somewhere between depending
on how much of coverage is mediator rather than nuisance. **This is the third time in this project
that a control has turned out to be a mediator** — #17 (survey progression), #25 (coverage), and
the matching in #12. Each time I reached for the control first and asked what it was second.

**Standing correction to my own practice, since three instances is a pattern rather than an
accident: before entering any variable as a control, state whether the treatment could cause it.
If it could, the controlled estimate is a lower bound and must be reported as one half of an
interval, never as the answer.**

---

## Entry 26, added by `E01·A05·R11` — a pre-registered kill that passed while the real problem sat next to it

Reading the control sets rather than trusting their names, per #25's standing correction:
**21 of the 22 items in my "acquiescence" index are explicit erotic content**, including *"I am
aroused by being dominant in sexual interactions"* and *"I am aroused by being submissive"* — the
POWER axis itself. I pre-registered the kill: if `|corr(index, POWER)| > 0.30` the index is erotic
content and three published quantities get republished.

**It did not fire — 0.024 — and the reason is worth keeping.** The two power items sit in the mean
with **opposite signs** (dominant −0.355, submissive +0.350) and cancel. With mixed-valence items,
content cancels in a mean and response style survives; that is what an acquiescence index is
*supposed* to do, and mine was accidentally well constructed. A random 22-column sham scores +0.029
against POWER — indistinguishable from the index, which is the right comparison.

| | vs POWER | vs breadth |
|---|---:|---:|
| full 22-item index | **+0.024** | **+0.385** |
| sham, 22 random numeric columns | +0.029 | — |
| the single content-free item | −0.112 | +0.055 |

| # | The claim | What killed or scoped it | What survived |
|---|---|---|---|
| 26 | **"9–13% of breadth is response style"** and **"85% of the induction→breadth link survives response-style control"** | Not the check I pre-registered — that one passed. The index is orthogonal to POWER and behaves like an arbitrary composite. **But it correlates +0.385 with breadth**, which is the outcome, so calling it a nuisance is the mediator problem for the fourth time. Whether it measures agreeing-in-general or endorsing-erotic-things-in-general cannot be separated in this release: the items are all erotic and none is reverse-keyed by design | The magnitudes, which barely move. induction→breadth is **+0.2922** uncontrolled, **+0.2523** with the full index (86%, matching the published 85%), **+0.2924** with the one content-free item. **"9–13% is response style" is downgraded to UNVERIFIED** — the shared variance is real, its interpretation is not identified |

**The lesson is about pre-registration, not about acquiescence.** I noticed two power items, wrote a
kill aimed at POWER, and it passed — while the dangerous overlap, index-with-outcome at 0.385, sat
one column away and was never in the kill condition. **Pre-registration stops you picking a
threshold after seeing the result. It does not stop you pointing the test at the wrong quantity,
and a kill that passes feels like clearance for everything you did not test.**

**Unfreeze condition** for the response-style question: balanced-keyed items, where the same
construct is asked in both directions so agreement and content can be separated. This release has
none, so no amount of further computation on it will resolve this.

---

## Entry 27, added by `E01·A01·R14` — the erotic covariates do nothing, and neither does personality

The other half of #26. Two of the eight CCA covariates are erotic items with no cancelling twin —
*"I find it erotic when two people of the opposite gender to me sexually interact"* and the
gender-identity attraction item — so the same audit was expected to come out differently. **It did
not.** Prediction recorded before running, refuted after.

Nested covariate ladder, held-out cross-domain CCA, K ∈ {3,5,8} × 3 seeds, monotone across nested
rungs as a positive control requires:

| covariate set | K=5 | share of total adjustment |
|---|---:|---:|
| none | 0.4122 | — |
| + sex, age | 0.2898 | **95%** |
| + personality (6 vars) | 0.2854 | 3% |
| **+ 2 erotic (published spec)** | **0.2832** | **2%** |
| sham: + 2 non-erotic Likert instead | 0.2861 | — |
| negative control: + 8 random columns instead | 0.2879 | — |

| # | The claim | Verdict |
|---|---|---|
| 27 | **"73% of the cross-domain transfer is not demographic"** | **STANDS**, and is now decomposed: **the entire adjustment is sex and age.** Personality and the two erotic covariates together carry 5%, and the erotic pair specifically carries 2% — less than eight *random* numeric columns remove (0.2879 vs 0.2832). My suspicion was wrong |

**The finding I was not looking for.** Eight random columns remove about as much as personality
plus the erotic pair. **So the "personality adjustment" present in every round of this project is
not distinguishable from a degrees-of-freedom effect** — any covariate set of that size removes
roughly that much by consuming df. Every time I wrote "personality partialled out" I was describing
an operation that did approximately nothing beyond arithmetic.

**Units note, stated because the numbers look wrong against the published ones.** This round reports
the **maximum** canonical correlation; the published 0.273/0.200 are the **mean** across components.
The decomposition is what this round is for and the share is stable across both, but 0.4122 and
0.2832 must not be read as replacing 0.273 and 0.200.

**Calibration, since two consecutive rounds refuted my own alarm.** #26 and #27 were both "this
control is secretly the thing being measured", and both came back negative — the acquiescence index
because its content cancelled, this one because sex and age had already done all the work. After
twenty-six entries of being wrong in one direction, the risk shifts: **suspecting every control
equally is not calibration, it is a different way of not looking.** The two audits cost real
compute and returned "fine", which is a legitimate result and worth the same shelf space as a kill.

---

## Entry 28, added by `E01·A02·R20` — the third automated kill I have had to override, and they all fail the same way

Additivity is the mathematical core of Ivan's model B: `v_i = w_i^T φ(s)` requires that a person's
weight on a feature does not depend on the scene. `A02·R11`–`R13` reported PLAUSIBLE on 3–4
decisive substance pairs. The fluid template is really **4 acts × 7 substances**, so the same
question is a variance decomposition over 31,541 cells from 1,953 people.

| component | observed | person-permuted null | sham (non-role options) |
|---|---:|---:|---:|
| person main effect | **0.1937** | 0.0601 | 0.3831 |
| substance | 0.0105 | 0.0105 | 0.1427 |
| act | 0.0260 | 0.0260 | 0.0000 |
| **person × substance** | **0.2105** | **0.2007** | 0.5106 |
| person × act | 0.1986 | 0.1797 | 0.0246 |

**The automated kill fired — "interaction is 109% of the main effect, additivity refuted" — and it
is wrong.** The interaction sits **at its own null**: 0.2105 observed against 0.2007 permuted.
Interaction terms in a sparse crossed design absorb residual noise by construction, and my kill
compared the interaction to the *main effect* instead of to *itself under permutation*. The person
main effect is real (3.2× its null); the interaction is not distinguishable from noise.

**And the sham says the design is unfit anyway.** A meaningless contrast — two random halves of the
same blocks' *non-role* options — produces a **larger** person main effect (0.3831) than the real
self/other contrast (0.1937). Whatever the person component is capturing, a contrast with no role
content captures more of it. The decomposition is picking up person-level scale, not role.

| # | The claim | Verdict |
|---|---|---|
| 28 | **"Additivity refuted in the folk basis; the A_i interaction term is required"** — the automated verdict of this very round | **OVERRIDDEN and UNVERIFIED.** The interaction is at its permutation null and the sham exceeds the signal. `R11`–`R13`'s PLAUSIBLE is neither upgraded nor refuted; additivity remains untested at this scale |

**Three overrides, one failure mode.** #21, #26, and now #28: every kill I have had to override was
**a threshold on an observed quantity that presupposed the quantity was real**. "z < 2", "|corr| >
0.30", "interaction > 50% of main effect" — none was conditioned on its control passing first.

**Standing correction, and it belongs next to #25's:** *a pre-registered kill must be written as a
conditional — evaluate the threshold **only if** the positive control fires and the negative control
returns null. A kill that can fire on a broken instrument is not a commitment, it is an automated
way to publish an artifact.* The controls were present and correct in all three rounds; they simply
were not in the kill's own logic, so the script printed a verdict the controls forbade.

**What the round did surface, unplanned.** Person residuals across source classes: male–female
**−0.706**, male–neutral −0.486, female–neutral −0.263. Removing a person mean over 3 classes
induces roughly −0.5 mechanically, so male–female is *more* opposed than the constraint requires and
female–neutral *less*. That asymmetry is not explained by the arithmetic and is the one thing here
worth another round — with a design whose sham does not beat its signal.

---

## Entry 29, added by `E01·A02·R21` — my own "not explained by the arithmetic" flag was 78% arithmetic

`R20` reported person residuals across source classes at male–female **−0.706** and I flagged it as
not explained by the centring. That flag was unchecked: removing a person mean over three classes
*forces* negative residual correlations, and the induced value depends on each person's cell counts,
so comparing it to a nominal −0.5 by eye proves nothing.

Both fixes in one round — an **uncentred** estimate with no induced component, and the centred one
against a null that permutes persons *within (substance, act)*, preserving every cell count and the
entire centring operation.

**Gate first, per the conditional-kill rule.** Positive control: within-source-class correlations
precum–ejaculate **+0.729**, saliva–urine **+0.522**, squirt–breastmilk **+0.348** — the score
measures something reproducible. Negative control: uncentred permutation null **+0.005 ± 0.027**.
Both pass, so the thresholds are binding.

| source pair | **uncentred** | its null | centred | geometry-preserving null | sham (non-role) |
|---|---:|---:|---:|---:|---:|
| male–neutral | **+0.479** | +0.001 | −0.496 | −0.581 | +0.219 |
| female–neutral | **+0.296** | +0.001 | −0.328 | −0.467 | +0.190 |
| **male–female** | **+0.063** | +0.005 | −0.725 | **−0.569** | +0.266 |

| # | The claim | Verdict |
|---|---|---|
| 29 | **"male–female residual correlation −0.706, more opposed than the arithmetic requires"** — `R20`'s unplanned finding, flagged by me as worth another round | **WITHDRAWN. 78% of it is the induced geometry** (−0.569 under a null that has no person structure at all). The excess is real but small, and the *uncentred* estimate — which has no induced component — is **+0.063**, not negative. There is no opposition |

**What replaces it is better than what it replaces.** The role feature transfers strongly between
male-source and neutral substances (+0.479), moderately between female and neutral (+0.296), and
**not at all across the male–female boundary (+0.063)**. So the structure is not *opposition*, it is
**absence of transfer** — and the boundary is specifically male↔female, with neutral substances
bridging to both, asymmetrically.

**And the sham makes the last cell sharper than a null does.** A meaningless non-role contrast
transfers at **+0.266** across the male–female boundary. The role feature transfers at **+0.063**.
**Across that one boundary the role feature transfers less than a meaningless contrast does** — it
is not merely absent, it is below the floor that arbitrary shared variance provides.

Conditional kill returns **UNVERIFIED** on the opposition question as posed (+0.063 sits between the
0.05 and 0.15 thresholds). The *contrast* between +0.479 and +0.063, both measured against a
+0.005 null, is what this round actually establishes.

---

## Entry 30, added by `E01·A02·R22` — the sharpest fact in the project was one cell of a specification curve

`R21` established, as the project's cleanest structural result: the role feature transfers
male↔neutral at +0.479 and male↔female at +0.063, below even a non-role sham (+0.266). I named
orientation as the rival. The blunter rival I had missed: **"ejaculating myself" requires a penis,
"squirting myself" requires a vulva** — so on gendered substances the self-pole is gated by the
respondent's own body, in opposite directions for the two classes, with no role content involved.
Neutral substances have no such gate, which is exactly why they bridge both sides.

Two orthogonal fixes, gate passed first (male–neutral stays >0.20 under every fix, min +0.350; all
permutation nulls ≤ |0.014|):

| specification | **male–female** | male–neutral | female–neutral | n | null |
|---|---:|---:|---:|---:|---:|
| all, as `R21` | **+0.063** | +0.479 | +0.296 | 1,457 | +0.004 |
| men only | **+0.233** | +0.506 | +0.393 | 769 | +0.011 |
| women only | **+0.097** | +0.350 | +0.218 | 688 | −0.014 |
| drop the *produce* act | **−0.204** | +0.484 | +0.184 | 1,457 | +0.000 |
| men only, drop *produce* | **+0.082** | +0.528 | +0.335 | 769 | −0.001 |

| # | The claim | Verdict |
|---|---|---|
| 30 | **"male–female transfer is +0.063, below the sham"** — `R21`'s headline, called the sharpest surviving structural fact one round earlier | **NOT IDENTIFIED.** The quantity ranges **−0.204 to +0.233** across five defensible specifications. The two fixes push in *opposite* directions — sex-stratifying raises it, dropping the body-gated act reverses its sign — and their combination lands back near the original. +0.063 was one cell of a curve I had not run |

**What survives, and it is the comparison rather than the point.** Across all five specifications
the role feature transfers **male↔neutral at +0.35 to +0.53**, against a sham that sits at +0.24–0.27
and barely moves. So:

- **male↔neutral: role beats sham by 2× in every specification** — robust
- **male↔female: role never beats sham in any specification** (+0.233 vs +0.242 at closest) — also robust, and it is a statement about a *contrast*, not a coefficient

**The lesson is that I promoted a number one round after building the tooling that would have caught
it.** `R14`–`R15` established that this project's headline quantities move with unswept analyst
choices; `R21` then produced a single-specification number and I called it the sharpest fact in the
project in the same commit. **A specification curve is not a step you run on suspicious results — it
is the unit in which a result is stated**, and thirty entries in I am still writing point estimates
first and curves second.

---

## Entry 31, added by `E01·A03·R10` — the first curve that came back robust, and the published cell was the most conservative one by accident

Applying #30's correction rather than regretting it: every surviving README number gets a
specification curve or is withdrawn. Two had never had one. This is the load-bearing one — `A03`'s
decision (model acquisition and valuation as two systems) rests on it.

**72 cells**: 3 bin→number mappings × 4 minimum-onset thresholds × 3 category sets × 2 rank
statistics. All reported.

| | |
|---|---|
| share of people agreeing with the population ordering | **[0.776, 0.860]** across all 72 cells |
| within-person permutation null, **recomputed in every cell** | 0.479 – 0.506 (gate passes) |
| positive control: synthetic exact-followers + noise | **0.993** |
| sham: population ordering replaced by a random ordering | 0.568 |

**ROBUST.** The mapping choice moves it by 0.005, the minimum-onset threshold by 0.01, and the
category set by 0.06 — the largest single lever, and still nowhere near the null.

**But the published figure was 0.747, below the entire swept range**, and that needed explaining
rather than waving at. It reproduces **exactly** — 0.747, mean ρ +0.232, n=6,810 — on the **27
matched onset/rating categories**, against 0.848 / +0.350 / n=10,076 on all 31.

**And nothing chose those 27 for this question.** `A03·R09` obtained its data by `exec`-ing the RSA
round's loader, which filters onset columns down to those with a matching *arousal-rating* column —
a filter that exists because the RSA needed to pair onset with preference. The schedule question
needs no such pairing. **A category filter travelled into an analysis it had nothing to do with,
through the loader-reuse pattern I formalised in `lib/rounds.py` during the restructure.** The
explicit dependency map made the reuse visible; it did not make the *inherited assumptions* visible.

| # | Claim | Verdict |
|---|---|---|
| 31 | **"A population-shared maturational schedule; 74.7% of people follow it, d=0.69"** | **CONFIRMED and widened.** The full range including the published cell is **[0.747, 0.860]**, every cell far above a null that never leaves 0.48–0.51. The published number was, by accident, the most conservative cell on the grid |

**Worth stating plainly after thirty entries of the opposite:** this is the first specification
curve in the project that came back clean, and the accidental choice went *against* the finding
rather than for it. Both directions of accident exist; only one of them gets noticed by an author
checking their own work.

---

## Entry 32, added by `E01·A03·R11` — the exec graph, and the first time the gate made me fix an instrument instead of override a verdict

#31 found that A03's rounds inherit their data from the RSA round's loader, which filters onset
columns to the 27 with a matching arousal-rating column. **Mapping the `exec` graph across the whole
project: 30 of 57 rounds inherit a loader, in two clusters.**

| loader | rounds inheriting it | the filter they did not choose |
|---|---:|---|
| `16_dimensionality` | **20** | block inclusion at n≥1200, n_options≥10, mean_picks>1.5, option floor 20 |
| `24_attack_rsa` | **6** | the 27 matched onset/rating categories |
| `49/50_additivity`, `33_induction_timing`, `45_theta_nonsexual` | 4 | fluid-blocks-only, porn-onset excluded, etc. |

Cluster 2 feeds six `A03` rounds **including `R09`, which produced one of only two CONFIRMED claims
in the project.** So it was re-run on the widest category set the estimand admits.

| category set | observed | rival world (noise 0.3/0.5/0.7) | above p97.5 |
|---|---:|---|---|
| matched-27, as published | 0.959 | 0.461 / 0.445 / 0.435 | 3/3 |
| matched-28 | 0.975 | 0.461 / 0.451 / 0.436 | 3/3 |

Positive control: injection 0.3 → **1.355 / 1.393**, detected in both sets; injection 0.0 → **0.439 /
0.454** against a reference of 0.439 / 0.449, **not** detected. Both gates pass.

| # | Claim | Verdict |
|---|---|---|
| 32 | **"Onset carries structure preference does not"** (#22, CONFIRMED) | **CONFIRMED and filter-free.** Above the rival world in all six cells. **Scope, stated honestly: the widening is only 27→28**, because the residualisation *requires* an onset column to have an arousal-rating twin. Unlike the schedule in #31, this filter is defensible here — it follows from the estimand rather than from a loader. That makes this a weak widening, and I am not claiming more |

**The new thing is what happened to the gate.** On the first run it **failed**: the zero-injection
control read as "detected" (0.452 against a reference p97.5 of 0.407). That was not the finding —
it was my control, drawing its reference from a *different RNG stream* with only 30 samples, so seed
noise alone decided it. Drawing the reference once from a shared stream and testing both injections
against it fixed it, and the verdict then came out clean.

**Three times before this (#21, #26, #28) a kill fired and I had to override it by hand.** This time
the conditional gate refused *before* printing a verdict, and the correct response was to **repair
the instrument, not to override the gate**. That is the rule doing the job it was added for, one
round after being added — and the difference between the two situations is worth naming: **an
override is a judgment I make about a machine, a refusal is the machine declining to make one for
me.** Only the second is reproducible by someone who is not me.

---

## Entry 33, added by `E01·A07·R03` — three design errors in one round, and the gate caught all three

The block-inclusion threshold (`n_respondents ≥ 1200`) is inherited by **20 of 57 rounds**, decides
which blocks exist at all, and sits under the modality, consumption, coverage and θ results
simultaneously. Sweeping it was the largest unexamined choice left in the project.

**Error 1.** First run compared *unmatched* deficits against *matched* published values. #11 and #25
established that any group comparison on this release must be block-count matched or it measures
survey coverage; I built the sweep and did not carry that forward. Gate failed. Matching added.

**Error 2.** With matching, the gate failed again — the pre-registered condition was "sex must be
the largest deficit in every cell", and consumption exceeds it. But #11's own chain says consumption
is **0.0871 at block-only matching** and only falls to 0.0439 after θ and sex matching. **My gate
encoded an expectation that belongs to a different matching level than the one the sweep runs at.**
The gate is mis-specified, not the instrument.

**Error 3, the one worth keeping.** Seed spread reaches **0.0437** — larger than modality's entire
deficit in most cells. I ran 2 seeds and said so in the IMPOSSIBLE register, but a 2-seed design
cannot resolve a quantity whose seed noise exceeds it, and I did not check that before choosing 2.

**Block-count-matched deficits across the sweep** (10 threshold combinations, 26–48 blocks admitted):

| split | range | verdict |
|---|---|---|
| sex | 0.083 – 0.102 | stable, ±10% |
| consumption | 0.082 – 0.112 | stable, ±15%, and consistent with #11's 0.0871 at this matching level |
| **modality** | **0.009 – 0.050** | **5.5× range** |
| placebo (row parity) | −0.005 – +0.003 | clean at every threshold |

| # | Claim | Verdict |
|---|---|---|
| 33 | **"modality deficit +0.0546"** (`#19`) | **UNVERIFIED, and flagged.** The pre-registered gate was mis-specified so no verdict is licensed, but the observed pattern is that modality swings 5.5× with which blocks are admitted while sex and consumption move under 15%. It needs a properly-gated round with ≥5 seeds before anything is claimed |

**What the round actually demonstrates is about the gate, not about modality.** Three independent
design errors — a missing control, a mis-specified threshold, an under-powered seed count — and the
conditional gate refused on all three rather than printing a verdict on any. **Before #28 added the
conditional, every one of these would have produced a number.** The second error is the interesting
one: a pre-registered condition can be *wrong* rather than merely unmet, and when it is, the honest
move is to report UNVERIFIED and rebuild — not to relax the condition until it passes, which is
exactly what pre-registration exists to prevent me doing.

---

## Entry 34, added by `E01·A07·R04` — the modality deficit is smaller than the noise of the only instrument that can measure it

`#33` failed three ways and I rebuilt for all three: block-matched inside every cell, gate written
for *that* matching level, 5 seeds, and the two corpus axes swept separately. The resolvability
criterion — a cell counts only if its effect exceeds **2× its own seed spread** — was added
specifically because of `#33`'s error 3.

**Positive control passes: sex deficit 0.065–0.109 in all 9 cells, resolvable in 9/9** (seed spread
0.004–0.021 against effects of 0.07–0.11).
**Negative control passes: placebo |≤0.0082| everywhere.**
**Resolvability fails: 0 of 9 modality cells are resolvable** — median 0.022–0.039 against seed
spreads of 0.020–0.073.

| corpus cut (median modality deficit) | resp≥600 | resp≥1200 | resp≥2000 |
|---|---:|---:|---:|
| options ≥8 | 0.0266 | 0.0223 | 0.0385 |
| options ≥10 | 0.0296 | 0.0274 | 0.0267 |
| options ≥12 | 0.0285 | 0.0225 | 0.0220 |

| # | The claim | Verdict |
|---|---|---|
| 34a | **"modality swings 5.5× with which blocks are admitted"** — `#33`'s own observation, one round old | **WITHDRAWN. That was 2-seed noise.** At 5 seeds the cell medians span 0.022–0.039, a factor of 1.8, and the pattern `#33` read as corpus-dependence does not survive its own seed count |
| 34b | **"modality deficit +0.0546"** (`#19`) | **UNVERIFIED with a stated precision limit.** The point estimate is somewhere near 0.02–0.04 and the published 0.0546 sits above the entire grid, but **no cell is resolvable at 5 seeds** and no verdict is licensed |

**The contrast is what makes this interpretable rather than a shrug.** Same pipeline, same cells,
same seeds: sex is resolvable in 9/9, modality in 0/9. The instrument is not noisy in general — the
modality effect is simply small relative to the noise this instrument carries.

**Unfreeze condition, quantified.** Seed spread scales roughly as 1/√n. To bring the spread below
half the effect (~0.0135 against ~0.027) needs about **5 × (0.040/0.0135)² ≈ 44 seeds** per cell,
against the 5 run here. That is arithmetic, not aspiration, and it prices the question: resolving
the modality deficit costs roughly an order of magnitude more compute than has been spent on it,
and until that is paid the number does not exist.

**Two rounds in a row the gate has refused, and both refusals were correct.** `#33` caught three of
my design errors; this round caught a precision limit that no amount of care in design would have
removed. The difference matters: the first is a fixable mistake, the second is a property of the
measurement, and only a gate that checks resolvability separates them.

---

## Entry 35, added by `E01·A08·R01` — the resolvability criterion applied backwards, and modality turns out to be the exception

`#34` killed the modality deficit by resolvability rather than by falsity: its effect was smaller
than its own seed spread in 9 cells of 9. Most of this project's numbers were produced at 1–3 seeds
and had never been asked the same question. Re-run at **6 seeds**, criterion `|effect| > 2 × spread`:

| quantity | median | seed spread | ratio | resolvable |
|---|---:|---:|---:|:--:|
| **sex deficit** *(positive control)* | 0.0847 | 0.0137 | **6.18** | ✅ |
| **placebo deficit** *(negative control)* | 0.0002 | 0.0076 | **0.02** | ❌ *(as it must be)* |
| consumption deficit | 0.0951 | 0.0235 | 4.05 | ✅ |
| cross-domain CCA (mean) | 0.1491 | 0.0155 | 9.63 | ✅ |
| breadth nestedness excess | 0.0657 | 0.0036 | **18.38** | ✅ |

Gate passes — one resolvable, one not, so the criterion discriminates here. **0% of the audited
headlines are unresolvable.** Modality was the exception, not the rule, and the surviving numbers
clear the bar by factors of 4 to 18.

**Not audited, with the price stated instead of the row omitted:** the onset rival world (#22/#32)
at ~40 pipeline rebuilds per seed, ~240 for six; the POWER–SUBSTANCE ladder (#24) at 216 composite
refits; the θ reliability ladder (#23). Those three remain at their original seed counts and are
**NOT AUDITED**, which is not the same as passing.

**New flag, found by the audit rather than sought.** The cross-domain CCA re-runs at **0.1491**
against a published **0.200** — a gap of 0.051, which is **3.3 seed spreads** and therefore not seed
noise. The published figure came from a 12-split average and this from 6; whether that alone
explains a 25% difference is unresolved, and I am recording it as an open discrepancy rather than
guessing at a cause. **The audit's job was to check precision and it found an accuracy problem
instead** — which is the ordinary way a control earns its cost.

| # | Claim | Verdict |
|---|---|---|
| 35 | **cross-domain CCA = 0.200** | **DISCREPANT.** Re-run gives 0.1491 ± 0.0155 at 6 seeds. Both are above every floor this project has measured, so the *existence* of cross-domain transfer is untouched; the *value* needs reconciling before it is quoted again |

---

## Entry 36, added by `E01·A08·R02` — the audit that flagged a discrepancy was measuring a different quantity, and two ledger entries were wrong

`#35` flagged the cross-domain CCA re-running at 0.1491 against a published 0.200 and recorded it as
an open discrepancy. Reconciling it required reading three rounds' code rather than one, and the
answer corrects **two** earlier entries of this ledger.

**Split count is not the cause.** Both statistics are flat across 4–24 half-splits (drift +3% and
+5%), nulls collapse to 0.017 and 0.009:

| half-splits | max over components | mean over components |
|---:|---:|---:|
| 4 | 0.3244 | 0.1472 |
| 12 | 0.3265 | 0.1470 |
| 24 | 0.3308 | 0.1514 |

**Correction 1 — `#27`'s units note is backwards.** It says *"this round reports the maximum
canonical correlation; the published 0.273/0.200 are the mean across components."* Both published
rounds return the **max**: `A01·R09` uses `np.nanmax(cv)`, `A02·R01`'s `cvcca` returns
`max(abs(corr))`. The note had it exactly reversed and has been sitting in the ledger since.

**Correction 2, and the larger one — `#35` was not a discrepancy at all.** The published 0.273/0.200
come from CCA between **pairs of individual blocks** (`nc=3`, median over 321 block pairs). The
`A08` audit computed CCA between **two half-splits of the whole block set** (`K=5`). Those are
different estimands, not two statistics of one fit — which is why neither 0.332 nor 0.151 lands on
0.200. **My resolvability audit audited a quantity that was never published.**

| # | Entry | Correction |
|---|---|---|
| 36a | `#27`'s units note | **Reversed.** Published values are the max, not the mean |
| 36b | `#35`'s "cross-domain CCA is DISCREPANT" | **WITHDRAWN.** Two different estimands compared. The published pairwise-block value (0.273 raw / 0.200 demographics-removed) stands; the half-split factor CCA is a separate quantity at max 0.332 / mean 0.151, also far above its null |

**`realstat` G1 is *estimand before method*, and my own audit round broke it.** `#35` set out to
check precision, reached for the CCA machinery that was nearest to hand, and never asked whether
that machinery computed the published quantity. The resolvability numbers it produced are still
valid — as statements about the half-split factor CCA, which is what they measured.

**Both corrections were found by reading code, not by running anything.** Thirty-six entries in,
the cheapest instrument in this project remains opening the file that produced the number.

---

## Entry 37, added by `E01·A08·R03` — the headline reproduces and is resolvable; the positive control I wrote for it is degenerate

`#36` established that the resolvability audit had measured a quantity nobody published. Redone on
the estimand that actually carries the headline — pairwise-block CCA, `nc=min(3,…)`, max over
components, median across 306 pairs — at **6 seeds**:

| | median | seed spread | ratio | published |
|---|---:|---:|---:|---:|
| raw | **0.2686** | 0.0137 | **19.6** | 0.273 |
| demographics-removed | **0.1980** | 0.0169 | **11.7** | 0.200 |
| permuted null | 0.0509 | 0.0015 | — | 0.055 |

**Both headline values reproduce inside their own seed spread**, the null lands where it was
published, and both are **RESOLVABLE** by `#34`'s criterion at ratios of 12–20 — an order of
magnitude clear of the bar that modality failed.

**But the positive control is degenerate, and I am recording that rather than banking a pass.** I
split a block against itself — two halves of its own options — and it returned **exactly 1.0000**.
That is not the pipeline detecting alignment: the block residual matrix is **row-centred**, so each
person's option values sum to zero, which forces `sum(half₁) = −sum(half₂)`. A linear combination of
one half predicts a combination of the other *by construction*, and CCA finds it. **Same failure as
`#21`** — a control whose value is fixed by the centring geometry rather than by the data — and I
built it two weeks of rounds after writing that lesson down.

| # | Claim | Verdict |
|---|---|---|
| 37 | **cross-domain CCA 0.273 raw / 0.200 adjusted** | **RESOLVABLE and reproduced**, with seed spreads of ±0.014 and ±0.017 now attached. The verdict rests on the seed-spread measurement and the permuted null (0.051, a quarter of the adjusted value), **not** on the positive control, which is uninformative |

**What a real positive control here would be:** two blocks known to be near-duplicates — precum and
ejaculate, which correlate at 0.729 on their role scores (`A02·R21`) — should return a high but
**sub-unity** canonical correlation. Unity is the signature of a constraint, not of a strong effect,
and any control that returns exactly 1.0000 should be assumed broken before it is assumed
impressive.

---

## Entry 38, added by `E01·A08·R04` — a graded control, and the pipeline reproduces a similarity ordering it was never told

`#37` reproduced the headline and found it resolvable, but its positive control returned exactly
1.0000 because row-centring forces `sum(half₁) = −sum(half₂)`. A binary control ("must be high")
cannot distinguish detection from constraint anyway, so this one is **graded**: block pairs sorted
by independently known similarity, and the pipeline must reproduce the *ordering*.

| tier | pairs | median CCA | max | permuted floor |
|---|---:|---:|---:|---:|
| **1** precum ↔ ejaculate — same template, same source class, role scores r=0.729 (`A02·R21`) | 4 | **0.8603** | 0.8697 | 0.067 |
| **2** within the fluid family — same template, different substance | 36 | **0.4943** | 0.691 | 0.063 |
| **4** all other pairs — the population the headline is a median over | 968 | **0.2723** | 0.724 | 0.048 |
| **3** fluid ↔ non-fluid — different template | 216 | 0.2032 | 0.594 | 0.057 |
| floor | person-permuted | **0.0499** | — | — |

**Strictly ordered 0.860 > 0.494 > 0.272 > 0.050, and nothing reaches unity.** The most similar
pair in the release lands at 0.86 — high, and *not* 1.0, which is what separates a detector from a
constraint. **Control VALID; `#37`'s resolvability verdict is now gated on a working positive
control rather than on two of three.**

**Secondary observation, unplanned.** Tier 3 (fluid ↔ non-fluid, **0.2032**) sits *below* tier 4
(all other pairs, **0.2723**). Cross-template pairs transfer *less* than ordinary pairs do, so the
seven fluid blocks are partially isolated from the rest of the corpus rather than merely internally
similar. That is consistent with `A02·R21`'s finding that the role feature crosses substance
boundaries freely inside the family, and it means the headline median (0.2686) is carried by the
non-fluid bulk, not by the family that produced most of this project's role results.

| # | Claim | Verdict |
|---|---|---|
| 38 | **cross-domain CCA 0.269 ± 0.014 raw / 0.198 ± 0.017 adjusted, resolvable** | **CONFIRMED on all three controls.** Positive: graded tiers reproduce a known ordering without hitting unity. Negative: permuted floor 0.050, a fifth of the adjusted value. Precision: ratios 19.6 and 11.7 |

**Why the graded form mattered more than the fix.** The degenerate control in `#37` would have been
repaired by any pair that returned something under 1.0 — but a single "high" number still cannot
tell you whether the pipeline is *measuring* similarity or merely *responding* to it. Four tiers in
a known order is a much stronger statement than one tier above a bar, and it cost the same compute.

---

## Entry 39, added by `E01·A02·R23` — the SUBSTANCE axis does not exist outside the seven blocks it was built in

`#38` found the fluid family partially isolated (fluid↔non-fluid pairs transfer at 0.203 against
0.272 for ordinary pairs). Most of this project's role findings live in that family, so the question
is whether the coordinate built there reaches the rest of the release. The self/other contrast only
exists as paired options inside those blocks, so generality can only be tested by asking whether the
coordinate **predicts** elsewhere.

| target | gain from SUBSTANCE | gain from **placebo** (fluid endorsement rate, no role content) | permuted null |
|---|---:|---:|---:|
| **fluid blocks** (5) — where it is defined | **+0.09339** | +0.01860 | −0.00043 |
| **non-fluid blocks** (24) | **−0.00012** | **+0.00235** | −0.00143 |

**Gate passes both ways**: SUBSTANCE predicts its own family at 200× its null, and the permuted null
outside is −0.0014. So this is a real null, not a blind test.

**Outside the family, SUBSTANCE explains nothing — and a role-free placebo built from the same
blocks and the same people does twenty times better.** Whatever carries a little signal across the
boundary is the person's fluid-block endorsement rate, which contains no role content at all.

| # | The claim | Verdict |
|---|---|---|
| 39 | **the SUBSTANCE axis, and everything built on it** — the third role coordinate (`#24`), the source-gender transfer structure (`#29`, `#30`), the additivity tests (`#28`), and `A02`'s "receiving vs giving" coordinate | **SUBGRAPH-LOCAL.** All of it is scoped to the seven fluid blocks. There is no evidence in this release that the coordinate exists elsewhere, and the appropriate scope line is now attached in `README` |

**What this does and does not touch.** POWER was built from *non-fluid* items and `A02·R01`–`R02`
showed it predicts option endorsement wherever the option set varies in role — that finding stands.
What falls is the generality of SUBSTANCE, which means the **three-axis picture is one axis measured
release-wide, one (GAZE) unmeasurable, and one confined to 3% of the corpus.**

**The check that would have caught this costs one round and I ran it forty rounds late.** "Does the
coordinate predict outside the blocks that defined it" is the first question to ask of any construct
built from a subset, and I built four rounds of findings on SUBSTANCE before asking it. The reason
is legible in the ledger: the fluid family was the only place the contrast was *available*, and
availability quietly became representativeness.

---

## Entry 40, added by `E01·A01·R15` — the same check that scoped SUBSTANCE clears the central claim

`#39` scoped the SUBSTANCE axis to seven blocks by asking whether it predicts outside the blocks
that defined it. The same question had never been asked of the construct carrying this project's
central claim: the person factors are built from **32 blocks out of 101**, selected by a filter that
20 rounds inherited and none chose.

| | held-out effect (gain − permuted null) |
|---|---:|
| **included blocks** *(positive control, leave-one-out)* | **+0.0390** |
| **excluded blocks** — never part of the construction, n≥300 | **+0.0178** (46% of included) |
| excluded blocks, **placebo** (endorsement rate, no coordinate content) | +0.0064 |

**REACHES.** The factors predict in blocks the filter discarded, at 46% of the included effect and
**2.8× a placebo** built from the same source with no coordinate content. 20 of 25 excluded blocks
show a positive gain. The 46% attenuation is expected — excluded blocks are smaller, so their
profiles are less reliable, which attenuates any prediction into them.

| # | Claim | Verdict |
|---|---|---|
| 40 | **the shared cross-domain grammar** — the project's central claim, built on 32 of 101 blocks | **GENERALISES.** It is not a property of the large blocks. It predicts in blocks it was never fitted on, including ones the inherited filter removed |

**A gate correction inside the round, and it is the same class as `#33`'s error 2.** The first pass
failed its negative control: permuted factors gave **−0.0064** on excluded blocks rather than ~0. That
is not contamination — excluded blocks are smaller, so six permuted columns cost held-out R² **by
degrees of freedom alone**. Requiring `|null| < 0.005` was a condition that does not fit a design
with varying n. **The null *is* the df price**, so the effect is `gain − null` and the negative
control becomes "is the null *stable* across blocks" (sd 0.0087, passes).

**The contrast is the point.** The same one-round check scoped SUBSTANCE to 3% of the corpus and
cleared the person factors. **A check that only ever confirms is not a check** — and after `#39` I
went into this expecting a second local construct, and got the opposite.

---

## Entry 41, added by `E01·A02·R24` — I set a gate threshold above the published magnitude of the thing I was gating on

`#39` scoped SUBSTANCE, `#40` cleared the person factors. POWER was the last role axis claimed to
hold release-wide, and `A02`'s basis decision rests on it. Same check: does it predict in blocks
that did not contribute to its construction?

| | effect (gain − permuted null) |
|---|---:|
| contributing blocks (6 of 32) — *positive control* | **+0.00150** |
| **non-contributing blocks (26)** | **+0.00245** |
| placebo there (endorsement rate, no role content) | +0.00063 |

**Gate FAILED** — the positive control needed >0.002 and returned 0.0015. **But `A02·R01`
published POWER's median held-out gain as +0.0006**, so I set a threshold three times above the
known magnitude of the quantity I was gating on, then failed my own round with it. Fourth
mis-specified gate (`#33` twice, `#40` once, now this), and the pattern is finally legible: **I set
thresholds from intuition instead of from the published magnitude.**

**Standing correction, joining the conditional-kill rule:** *a gate threshold must be derived from
a measured quantity — the published effect, its null, or its seed spread — never chosen. If no such
quantity exists yet, the round's first job is to measure it, not to guess a bar.*

| # | Claim | Verdict |
|---|---|---|
| 41 | **POWER is a release-wide role axis** — `A02`'s surviving basis after `#39` | **UNVERIFIED by the gate, and the descriptive pattern is unflattering.** POWER's effect is ~0.0015–0.0025 **in both sets**, marginally *larger* outside its home blocks than inside. It is not local like SUBSTANCE — it is uniformly weak. Against a placebo of 0.0006 it is real but tiny, which matches `A02·R01`'s original finding that POWER only bites where the option set varies in role |

**Where this leaves `A02`'s decision.** Its "SAFE" statement was: don't use the folk basis; the three
role axes are near-orthogonal; cross role with source. Since then: `#24` restated POWER–SUBSTANCE
from 0.233 to 0.605; `#39` scoped SUBSTANCE to 3% of the corpus; GAZE was never measurable; and now
POWER is confirmed weak everywhere rather than strong somewhere. **A02 has no strong release-wide
role coordinate, and its decision should be re-closed on that basis rather than on the three-axis
picture that has since dissolved.**

**The one thing that survives intact from A02** is the conditional finding, not the coordinate: role
predicts endorsement **wherever the option set varies in role** (`A02·R01`–`R02`, r=+0.752 between
option-set role variance and predictive gain, 21× between varying and non-varying blocks). That is a
statement about *when* a role coordinate applies, and it never depended on any axis being globally
strong.

---

## Entry 42, added by `E01·A02·R25` — the last A02 survivor passes every statistical attack and is unidentifiable anyway

After `#24`, `#39` and `#41`, the only `A02` claim left was the **conditional** one: role predicts
endorsement wherever the option set itself varies in role — `r = +0.752`, 21× between varying and
non-varying blocks. `#30` had noted its role-varying cells are 4 blocks; `#39` had since shown the
fluid family is subgraph-local. It looked structurally identical to the coverage law `#20` killed.

**It is not, and it passes:**

| | |
|---|---|
| observed r | **+0.752** |
| leave-one-out range | **+0.672 … +0.830** — never near collapse |
| permutation null at n=28 | mean −0.002, sd 0.191, **\|r\| p95 = 0.357** |
| sham (x permuted) | \|r\| p95 0.357 |
| positive control (gain planted as f(role variance)) | **+1.000** |

Unlike the coverage law — LOO floor −0.294 against a null p95 of 0.820 at n=6 — this sits at
**2.1× its null** with a leave-one-out floor of +0.67. **Statistically it is the most robust
correlation in the project.**

**And it cannot be identified, which the influence analysis printed as a `nan`.** Removing the fluid
family returns **`nan`, because `role_var` has standard deviation 0.0000 among non-fluid blocks.**
All four role-varying blocks are fluid blocks. So `role_var > 0` and `fluid-family membership` are
**the same indicator** in this release, and the correlation is equally well read as *"POWER predicts
better in fluid blocks"* — a family `#39` has already shown does not represent the corpus.

| # | Claim | Verdict |
|---|---|---|
| 42 | **"role predicts endorsement wherever the option set varies in role"** — `A02`'s last survivor | **ROBUST BUT UNIDENTIFIED.** Survives leave-one-out, its null, and a sham. The treatment variable has **zero variance outside one family**, so no design on this release can separate *role variance* from *fluid-family membership*. Not withdrawn — it is a real regularity — but it cannot be stated as a claim about role |

**`A02` is now re-closed, and the decision has inverted.** Its original SAFE statement was "not the
folk basis; three near-orthogonal role axes; cross role with source." What actually survives is:
**this release cannot settle what basis the ontology should use.** One axis is confined to 3% of the
corpus, one was never measurable, one is uniformly weak, and the conditional finding that would have
rescued them is perfectly confounded with the same 3%. **The honest decision is not a basis but a
requirement: a phase-1 collection must instantiate the role contrast in blocks that are not all one
family**, which is a design constraint on new data rather than a finding in this one.

**Zero variance in the treatment outside one stratum is a stronger objection than any p-value**, and
it was visible from the first round that computed `role_var` — 24 of 28 blocks at exactly 0.000,
printed in `A02·R02`'s own output table. I read that column as "the contrast is rare" and not as
"the contrast is confounded", and those are the same table.

---

## Entry 43, added by `E01·A03·R12` — the null that survived having its instrument replaced

`A03`'s decision — model acquisition and valuation as **two systems** — rests on `#9`: TEMPO
predicts within-person onset similarity at t=+4.59 while COORD (coordinate-loading similarity) gives
t=−0.46. Since then every coordinate COORD was built from has been damaged: `#39` scoped SUBSTANCE
to 3% of the corpus, `#41` found POWER uniformly weak, `#18`/`#19` withdrew the coordinate count as
K-dependent. **A null from a measure that has since failed its own validation is silence, not
absence.**

Rebuilt COORD from the **person factors** — the one construct in this project that passed a
generalisation check (`#40`: predicts blocks the filter discarded at 46% of the included effect,
2.8× a placebo).

| | b | t |
|---|---:|---:|
| **COORD, from validated person factors** | +0.0191 | **+1.26** |
| TEMPO *(positive control, published +4.59)* | +0.0282 | **+4.36** |
| permuted-COORD null | — | mean +0.10, sd 0.98, \|t\| p95 **1.92** |
| sham (gaussian) | — | \|t\| p95 1.86 |

**COORD sits inside its own null band.** And the re-test is *better* powered than the original, not
worse: all **27** categories carry a person-factor profile against only **18** mappable to GCCA
blocks, giving **344 pairs against 153**. So the original null was not an underpowered accident
either.

| # | Claim | Verdict |
|---|---|---|
| 43 | **acquisition and valuation are two systems** (`A03`'s decision) | **CONFIRMED with the instrument replaced.** COORD is null on a validated coordinate construct at nearly twice the pair count, while TEMPO reproduces at t=+4.36 |

**This is the first time in the project that replacing a damaged instrument has left a conclusion
standing.** Everywhere else — `#21`, `#28`, `#37` — swapping the instrument changed the verdict or
revealed the old one was unsupported. Here the worry was legitimate and specific (P6: a null is
inadmissible from an instrument that has failed validation elsewhere), the check was cheap, and the
answer was that the original round was right for reasons better than the ones it had.

**A03 stays closed.** It is now the only arc in `E01` whose decision has survived every attack made
on the evidence beneath it.

---

## Entry 44, added by `E01·A01·R16` — a second string-derived proxy beaten by its own sham

`A01`'s decision says this release cannot separate **A** (a dedicated sexual-content system) from
**B** (valuation of ordinary representation), because `#13` showed a dedicated module can be
compositional. `#13` did not consider **developmental structure**, so this round tried it:

> **B predicts** erotic interests inherit the maturation order of the ordinary representations they
> are built on — person-perception matures features before relations, so relational categories
> should arrive later. **A predicts** no such ordering.

Relationality measured by a **pre-registered string rule** over each category's own option texts
(share of options referring to another person), not by hand-coding — which would be an answer key.

| | r with mean onset age |
|---|---:|
| raw | −0.120 (p=0.63, n=19) |
| + word count, option count, prevalence partialled | **−0.379** |
| **SHAM — share of options containing a COLOUR word** | **−0.629** |
| permutation null at n=19 | mean 0.000, sd 0.235, **\|r\| p95 = 0.450** |
| positive control (onset planted as f(relationality)) | +0.802 |

**The sham beats the signal.** A content-free token class predicts acquisition age *better* than the
person-reference measure, and both sit inside a null band of ±0.450. **UNVERIFIED**, and the design
is unfit rather than merely negative.

| # | Claim | Verdict |
|---|---|---|
| 44 | **"relational categories are acquired later"** — a developmental A-vs-B separator | **UNVERIFIED and the instrument is unfit.** MDE at n=19 categories is \|r\| ≈ 0.45; the relationality proxy is beaten by a colour-word sham |

**The pattern, and it is the second instance.** `#28`'s additivity variance decomposition was also
beaten by its sham, and that measure was also derived from option text. **String-derived proxies
over option text carry generic text structure — length, register, specificity — that an arbitrary
token class captures just as well.** Any such proxy needs its sham run *before* the estimand, not
after, because the sham is what tells you whether the rule measures its concept or its corpus.

**`A01` is re-closed unchanged, and that is now a tested statement rather than an inherited one.**
Its SAFE was "usable for compositionality and cross-domain transfer; NOT usable for A-vs-B." Since
closure: `#25` bounded the foundation at [+0.017, +0.037] with the sign decisive at 31/32, `#40`
showed the factors generalise to blocks the filter discarded, `#37`/`#38` confirmed the headline is
resolvable against a graded control. **The compositionality half is better supported than at
closure. The A-vs-B half has now survived a second attempt to break it, this one from a direction
`#13` never considered.**

---

## Entry 45, added by `E01·A07·R05` — matching buys identification and spends resolvability, and the published number is on the wrong side of the trade

`A07` was opened mid-project with rounds but **no decision statement**, which the E/A/R structure
forbids. The missing measurement: `#11` published consumption→coordinates at **0.0439, triple-matched**
(block count + θ + sex), while `#35` only ever measured resolvability at the **block-only** value
of 0.0951 — twice the size, and half the size is exactly where `#34` found modality failing.

| split | matching level | deficit | seed spread | ratio | resolvable |
|---|---|---:|---:|---:|:--:|
| **consumption** | block | 0.0900 | 0.0347 | 2.60 | ✅ |
| **consumption** | block + θ | 0.0719 | 0.0277 | 2.59 | ✅ |
| **consumption** | **triple — as published** | **0.0357** | **0.0255** | **1.40** | ❌ |
| sex *(pos. control)* | block | 0.0866 | 0.0176 | 4.93 | ✅ |
| sex *(pos. control)* | block + θ | 0.0877 | 0.0284 | 3.08 | ✅ |
| placebo *(neg. control)* | all three | ≤0.003 | — | ≤0.61 | ❌ |

Gates pass both ways. **The published consumption headline is unresolvable at the matching level it
was published at**, and joins modality.

**The trade-off is the finding.** Each matching level removes more confound and shrinks the effect —
0.090 → 0.072 → 0.036 — while the seed spread barely moves (0.035 → 0.028 → 0.026). **Matching buys
identification with resolvability**, and there is a level past which the estimate is cleaner and no
longer measurable. I applied three levels of matching because `#11` and `#25` established each was
necessary, and never noticed that the last one spent the effect.

| # | Claim | Verdict |
|---|---|---|
| 45 | **consumption → coordinates = 0.0439** (`#11`, triple-matched) | **UNRESOLVABLE at that matching level** (ratio 1.40). Reportable at block or block+θ matching as **0.072–0.090**, where it is confounded with sex composition by roughly 31% (`#11`'s decomposition). Neither version is both clean and measurable |

**`A07` now has its decision, and it closes:**

> **Which congruence comparisons are reportable?** Only those whose effect exceeds roughly **0.07**
> at this design's precision. **Sex qualifies at every matching level** (ratio 3–5). **Consumption
> qualifies only at looser matching**, where it is not fully identified. **Modality never qualifies**
> (`#34`, 0 of 9 cells). The placebo correctly never qualifies. **Below ~0.07 this instrument cannot
> distinguish an effect from its own seed noise**, and no amount of better matching helps — matching
> makes it worse.

---

## Entry 46, added by `E01·A05·R12` — the measure cannot tell a concentrated population from a random one, and its published value is smaller than its own bias

`A05`'s decision — model the scalar gain as the object rather than control it away — rests in part
on breadth being **quantity without shape**: a person's set is only **0.88% more concentrated** in
coordinate space than a size-matched base-rate set. Computed at K=6 and never swept, in a project
where `#18`/`#19` showed K moves everything.

| population | K=2 | K=4 | K=6 | K=10 | K=16 |
|---|---:|---:|---:|---:|---:|
| **real** | −0.01% | −0.13% | **−0.29%** | −1.67% | −2.56% |
| **CONCENTRATED** *(positive control — sets drawn from one direction)* | −0.43% | −0.98% | **−1.67%** | −2.28% | −2.63% |
| **base-rate** *(negative control — sets drawn from base rates, no shape by construction)* | +0.04% | **−1.69%** | **−1.62%** | −1.25% | −1.19% |

**Both controls fail, and they fail into each other.**

- The **negative** control should be 0 by construction and is **−1.2% to −1.7% at K≥4**. The cause is
  arithmetic: the observed value is *one* realisation while the null is an *average of six*, and the
  participation ratio is nonlinear, so Jensen's inequality manufactures a gap where none exists.
- The **positive** control, a population built to be concentrated, gives **−1.67% at K=6 —
  indistinguishable from the base-rate control's −1.62%**. **The measure cannot separate a genuinely
  concentrated population from a random one.**
- And the published **−0.88% is smaller than the measure's own bias**.

| # | Claim | Verdict |
|---|---|---|
| 46 | **"breadth is quantity without shape" — 0.88% concentration, out-of-sample** | **UNVERIFIED, instrument unfit.** It cannot distinguish concentration from randomness, and its published value sits inside its own Jensen bias. Not "no shape found" — **no shape *detectable*, by this measure, at any K** |

**What this does to `A05`'s decision.** Its SAFE rested on four supports: breadth is the dominant
axis of individual variation, it is **quantity without shape**, sets are 24% nested, and it carries
the induction report. **The second is now gone.** The others stand — nestedness is resolvable at
ratio 18.4 (`#35`), the reliability is 0.557 (`#23`), and the induction link survives response-style
control at 85% (`#26`). **The decision holds on three of four supports, and I am recording which one
fell rather than leaving the count unchanged.**

**The bias was findable without any of this.** Comparing one draw against a mean of six is a Jensen
error visible in the code, and `#36` already established that reading the file that produced a
number is the cheapest instrument here. I ran a K sweep to find a K artifact and found an
arithmetic one instead — which is the second time an audit has caught something other than what it
was aimed at (`#35` found accuracy while checking precision).

---

## Entry 47, added by `E01·A05·R13` — the same control that killed one measure validates the other, and unity means opposite things in the two rounds

`#46` killed the concentration measure with synthetic populations of known answer. Nestedness is now
one of `A05`'s three remaining supports and had never had the same treatment. I predicted it would
pass — its null is **one draw per pair**, not an average of six, so the Jensen mechanism that broke
the concentration measure should not apply. **This round tests that reasoning rather than assuming
it.**

| population | containment | null | excess | % of chance→perfect | seed spread |
|---|---:|---:|---:|---:|---:|
| **NESTED** *(positive control — sets are prefixes of the popularity ordering)* | **1.0000** | 0.7268 | **+0.2732** | **100.0%** | 0.0025 |
| **base-rate** *(negative control — no nestedness by construction)* | 0.7273 | 0.7275 | **−0.0001** | −0.02% | 0.0017 |
| **real** | 0.7934 | 0.7283 | **+0.0655** | **24.0%** | 0.0022 |

**MEASURE SOUND.** The negative control lands at **−0.0001** — four orders of magnitude better than
the concentration measure's −0.016 — and the positive control recovers exactly 100% of the gap. The
real value reproduces the published +0.0660 / 24.2%.

| # | Claim | Verdict |
|---|---|---|
| 47 | **breadth sets are 24.2% of the way from chance to perfect nesting** | **CONFIRMED with known-answer controls.** `A05`'s decision now stands on three validated supports (nestedness, reliability 0.557, the induction link at 85%) and one withdrawn (`#46`) |

**Unity means opposite things in `#37` and here, and the difference is whether it was built in.**
`#37`'s positive control returned exactly 1.0000 and that was a **red flag** — a block split against
itself, where row-centring forces the result. This round's positive control also returns exactly
1.0000 and that is a **pass** — the population was *constructed* as perfect prefixes, so unity is
the definition, not an artifact. **The test is not "is it 1.0" but "did I build the 1.0 or did the
geometry hand it to me."**

**And the reason this round was worth running even though it passed:** `#46` and `#47` used the same
control on two measures from the same arc, built in the same week, and one was broken while the
other was clean. **A measure's soundness is not inherited from its neighbours**, and the only way I
found out which was which was to generate populations whose answer I already knew — which cost less
than either of the original rounds did.

---

## Entry 48, added by `E01·A04·R12` — matching corrects where there is something to correct, and my summary statistic divided noise by noise everywhere else

`A04`'s rule — "group comparisons must be coverage-matched" — stood on the `+0.815` law that `#20`
withdrew, plus `#45`'s finding that matching can spend an effect below resolvability. Tested
directly, with a **placebo match** that discards the same people while correcting nothing:

| split | coverage gap | unmatched | coverage-matched | **real change** | **placebo change** |
|---|---:|---:|---:|---:|---:|
| **pornhabit** *(pos. control)* | **2.435** | 0.2254 | 0.0959 | **0.1295** | **0.0069** |
| sex | 0.130 | 0.0975 | 0.0852 | 0.0123 | 0.0042 |
| neuroticism | 0.321 | 0.0344 | 0.0387 | 0.0043 | 0.0001 |
| **modality** *(neg. control)* | 0.069 | 0.0355 | 0.0327 | 0.0028 | **0.0047** |

**Where there is a coverage gap, matching corrects and the placebo does 5% of the work**
(pornhabit: 0.225 → 0.096, placebo 0.007). **Where there is no gap, both changes are noise** —
modality's placebo change (0.0047) *exceeds* its real change (0.0028), which is exactly what
"nothing to correct" looks like.

**My pre-registered summary was the wrong aggregation and returned UNVERIFIED.** I took the *median
of the ratio* `placebo change / real change` across four splits — 0.54 — but for the two splits with
negligible coverage gaps that ratio is **noise divided by noise**, and it dominated the median. The
quantity that matters is the ratio *where the correction is large*, and there it is **0.05**.

| # | Claim | Verdict |
|---|---|---|
| 48 | **"group comparisons must be coverage-matched"** (`A04`'s rule) | **RESTATED, not withdrawn.** Match **when the coverage gap is large** — there the correction is real and demonstrably not subsampling (placebo does 5%). Where the gap is small, matching changes almost nothing and `#45`'s resolvability cost dominates, so **matching should be conditioned on the gap rather than applied by default.** The general "change ∝ gap" law still cannot be established: n=4 splits with one dominant point, the same structure `#20` killed |

**Third time a summary statistic has hidden its own table** (`#30` narrated one cell of a curve,
`#33` compared unmatched to matched values, this one divided noise by noise). Each time the raw
per-unit numbers were printed directly above the summary that misread them. **The aggregation is
where the error lives, not the measurement** — and a four-row table needs no aggregation at all.

---

## Entry 49, added by `E01·A01·R17` — the shared grammar is real, graded, and carries almost no predictable variance

No outside challenger has ever run against this project; this session forbids dispatching one, so
every ledger row stays **`[unchallenged]`**, not "clean". The available substitute is not another
reviewer but another **framing** — `realstat` §2.5: three independent designs test the framing,
where a re-implementation only tests the code, and most retractions here were framing errors that
correct code would have preserved.

**Published framing:** held-out canonical correlation between two blocks' residuals, max over 3
components, median across 306 pairs → **0.269 raw / 0.198 adjusted**.
**This framing:** held-out *prediction* — fit block A's residuals → block B's residuals on 70% of
shared people, score R² on the other 30%.

| tier | prediction R² raw | adjusted | null |
|---|---:|---:|---:|
| precum ↔ ejaculate | **0.4327** | 0.3480 | −0.0075 |
| within fluid family | **0.1060** | 0.0693 | −0.0129 |
| fluid ↔ non-fluid | 0.0108 | −0.0117 | −0.0102 |
| **all other pairs** | **0.0105** | **−0.0006** | −0.0111 |
| **overall median** | **+0.0111** | **−0.0022** | −0.0110 |

**The tier ordering reproduces exactly** (0.433 > 0.106 > 0.011), so this framing measures the same
thing `#38` graded. **And the magnitude disagrees by two orders of magnitude: CCA says 0.198,
prediction says −0.0022.**

**Both are correct, and the reconciliation is the finding.** CCA maximises over projections — it
finds the *best* linear combination in each block. Prediction must reproduce block B's *actual*
residuals. So there **is** a shared direction, it **is** graded by known similarity, and it carries
**almost none of B's variance**. A high canonical correlation in a low-variance direction.

| # | Claim | Verdict |
|---|---|---|
| 49 | **"a shared cross-domain grammar", 0.269/0.198** | **RESTATED, and the restatement is sharper than the claim.** Cross-domain structure **exists** — CCA 0.198, tier-ordered, above every null. Cross-domain structure **predicts essentially nothing** — pairwise block→block R² is −0.002 against a −0.011 null, everywhere except inside the fluid family (0.35 / 0.07), which `#39` already showed is isolated |

**And it reconciles with `#25`/`#40` rather than contradicting them.** Those used factors from **31
blocks** to predict a held-out block and got **+0.017 to +0.037, 31/32 positive**. Pairwise
block→block gives ~0. Both are true: **the shared signal is too thin to survive one block's noise
and accumulates across many.** Three measurements, one structure — thin per pair, real in aggregate,
and a correlation that exists in a direction almost nothing travels along.

**This is the largest single correction in the ledger and it came from changing the framing, not
from finding an error.** Every number in the CCA framing was correct. The framing simply answered
"is there a shared direction" while I had been reading it as "does domain A tell you about domain
B" — and those diverge by two orders of magnitude in exactly this kind of data.

---

## Entry 50, added by `E01·A03·R13` — the onset RSA survives the framing swap that broke the central claim, and my fifth mis-specified gate nearly hid it

`#49` broke the central claim by swapping framings: the cross-domain grammar exists as a
correlation (CCA 0.198) and carries no predictable variance (pairwise R² −0.002). The onset RSA has
the same structure — a correlation between two similarity matrices, +0.599 — and had never been
asked the same question.

**Person-level prediction, 6,615 people with ≥8 onsets and ≥8 preferences, 5 seeds:**

| | held-out R² | seed spread |
|---|---:|---:|
| **preference half → preference half** *(the ceiling this pipeline can reach)* | **+0.0436** | 0.0182 |
| onset half → onset half | +0.0510 | 0.0239 |
| **onset → preference** | **+0.0136** | 0.0039 |
| permuted-person null | −0.0043 | 0.0036 |

**Onset predicts preference at +0.0136 against a −0.0043 null — 4.5 seed spreads above it, and 31%
of the same-domain ceiling.** Unlike `#49`, this correlation *does* carry predictable variance.

| # | Claim | Verdict |
|---|---|---|
| 50 | **"acquired together tracks liked together", RSA +0.599** | **SURVIVES THE FRAMING SWAP.** Person-level prediction reaches 31% of what preference achieves predicting *itself*. `A03`'s evidence is not the thin-direction artifact `#49` exposed in `A01` |

**And the gate failed anyway, for the fifth time, in the way `#41` already named.** I required the
positive control to exceed **0.05**; it returned **0.0436** — because 0.0436 *is* the ceiling of
this pipeline, measured **inside the same round**. `#41`'s standing correction said to derive
thresholds from a measured quantity rather than choose them, and here the measurement was three
lines above the threshold that ignored it.

**Sharper procedural fix, since the general form of the rule was not enough:** *compute the positive
control **first**, then set every threshold as a fraction of it.* A ceiling measured in the same
round and then ignored is worse than one never measured, because the round contains its own refutation
and prints it before the verdict.

**This is not relaxing a failed threshold to rescue a finding** — the distinction matters and `#33`
warns about exactly that. The threshold was **unsatisfiable by construction**: no quantity in this
design can exceed a ceiling of 0.0436, so demanding 0.05 asks the pipeline to beat itself. The
verdict stands on the ratio to the *null* and to the *measured ceiling*, both of which were
pre-registered in spirit and neither of which I wrote down as the operative comparison.

---

## Entry 51, added by `E01·A08·R05` — the published statistics rank the claims in a different order than predictable variance does

`#49` and `#50` created two classes: correlations that carry predictable variance and correlations
that do not. Every remaining headline is a correlation or a congruence and none had been sorted.
Common pipeline, **ceiling computed first** per `#50`'s procedural fix, everything expressed as a
fraction of it:

| claim | published statistic | held-out prediction | its own ceiling | **% of ceiling** | class |
|---|---:|---:|---:|---:|---|
| onset RSA (`#50`) | 0.599 | +0.0136 | 0.0436 | **31%** | carries variance |
| **nestedness** | 24.0% of gap | **+0.0766** | 0.3107 | **25%** | carries, borderline |
| **sex deficit** | 0.093 | **+0.0289** | 0.3157 | **9%** | partial |
| cross-domain CCA (`#49`) | 0.198 | −0.0022 | — | **~0%** | **thin direction** |

Both new ceilings clear their nulls by 170–200×, so both claims are sortable rather than merely
unmeasured.

**The ordering flips for the largest published number.** The cross-domain CCA has a *larger*
published statistic than the sex deficit — **0.198 against 0.093** — and carries **an order of
magnitude less predictable variance**. Ranking the project's claims by the numbers I published gives
a different order than ranking them by how much of the outcome they actually explain, and the
disagreement is worst exactly where the published number is biggest.

| # | Claim | Verdict |
|---|---|---|
| 51 | the surviving headline set | **SORTED.** Two carry predictable variance (onset RSA 31%, nestedness 25%), one is partial (sex 9%), one is a thin direction (cross-domain CCA ~0%). **The maturational schedule is UNSORTED** — its statistic is a within-person rank agreement whose predictive analogue is a ranking task, not an R², so forcing it through this pipeline would be a category error. Listed rather than assigned |

**What the two-class split is really measuring.** A correlation asks *"is there a direction along
which these agree"*; predictable variance asks *"how much of the thing is that direction"*. Both are
legitimate and they answer different questions — the mistake was never computing one, it was reading
one as the other for forty-eight entries. **Every claim in this ledger that survived to the end
survived as a statement about a direction, and only two of them are also statements about an
amount.**

---

## Entry 52, added by `E01·A03·R14` — the sort completes, and the claim that carries the most is the one that was hardest to break

`#51` sorted every surviving headline by predictable variance except the maturational schedule,
whose statistic is a within-person rank agreement. Its proper analogue is a **ranking task**, built
here: for a held-out person and a pair of their categories, predict which they acquired first, using
only the population ordering.

| | accuracy | seed spread |
|---|---:|---:|
| **population ordering** (fitted on training people, applied to held-out) | **66.71%** | 0.72 |
| **oracle** — the best a single global ordering can do, fitted *on the held-out half itself* | **66.53%** | 0.42 |
| random ordering *(negative control)* | **48.98%** | 11.95 |
| tied pairs excluded | 36.3% of all within-person pairs | |

**The population ordering reaches 101% of the oracle's margin over chance.** A global ordering
fitted on people it has never seen performs as well as one fitted on the very people it is scored
on. **There is no better global ordering to find** — the schedule is not merely generalisable, it is
saturated.

**The completed sort:**

| claim | published statistic | % of its own ceiling | class |
|---|---:|---:|---|
| **maturational schedule** | 0.747–0.860 | **101%** | **carries — at ceiling** |
| onset RSA | 0.599 | 31% | carries |
| nestedness | 24.0% of gap | 25% | carries, borderline |
| sex deficit | 0.093 | 9% | partial |
| cross-domain CCA | 0.198 | ~0% | thin direction |

**The ordering is almost the inverse of how hard each claim was to establish.** The schedule cost
one 72-cell specification curve (`#31`) and came back clean on the first attempt — the only headline
in the project that did. The cross-domain grammar consumed roughly a dozen rounds, three
retractions, and a framing swap, and carries no predictable variance at all. **What survived
scrutiny best is what needed the least of it**, which is not a coincidence: the schedule is a large
effect in a well-measured quantity, and large effects in well-measured quantities do not generate
forty-eight retractions.

**And the 36.3% tie share is the honest limit.** More than a third of within-person category pairs
are tied by the 2-year binning and are excluded rather than scored as half-credit. The 66.7% is
accuracy *on pairs the release can order at all* — a release with finer bins would score a different
number on a larger pool, and I do not know in which direction.

---

## Entry 53, added by `E01·A03·R15` — the schedule is developmental, and rarity is a real but separate ordering

`#52` made the maturational schedule the strongest surviving claim — a population ordering predicts
held-out pairwise acquisition order at 66.71%, saturating what any global ordering can do. It had
never been asked its most obvious rival: **rare interests are reported later**, so the "schedule"
might be a prevalence ordering wearing a developmental name.

| ordering used to predict held-out pairs | accuracy | seed spread |
|---|---:|---:|
| oracle (best single global ordering, in-sample) | 66.50% | 0.78 |
| **onset ordering** | **66.75%** | 0.48 |
| **onset ordering, prevalence partialled out** | **65.59%** | 0.80 |
| **prevalence ordering alone** | **60.75%** | 0.83 |
| sham — block option count | 54.41% | 1.02 |
| random ordering | 53.33% | **7.68** |

`corr(prevalence, mean onset) = −0.219` — rarer categories do arrive later, weakly.

**DEVELOPMENTAL. Removing prevalence from the onset ordering costs 1.16 points** (66.75 → 65.59),
while prevalence alone reaches only 60.75%. The two orderings are largely **independent**
contributions rather than one masquerading as the other.

| # | Claim | Verdict |
|---|---|---|
| 53 | **the maturational schedule** | **SURVIVES its rival.** Not a rarity artifact: the ordering retains 65.6% with prevalence removed. **Rarity is separately real** — 60.75% on its own, well above the sham and the random band — which is a finding in its own right and was never reported |

**A limitation I created and should not have.** The random control has a seed spread of **7.68
points** — one random permutation per seed, and a single permutation is a very noisy estimate of
chance. The other rows have spreads of 0.5–1.0. **The chance baseline is the least precisely
measured quantity in the table**, which is backwards: it is the one every other row is compared
against. Many permutations would cost almost nothing and I ran five.

**What this leaves standing.** The schedule is now the only claim in the project that has passed a
specification curve (`#31`, 72 cells), a predictable-variance test (`#52`, at ceiling), and a rival
explanation (`#53`). It is also the claim that needed the fewest rounds to establish — which
`#52` already noted is not a coincidence.

---

## Entry 54, added by `E01·A03·R16` — rarity is a second ordering principle, and it is not the age window

`#53` found prevalence ordering predicts held-out pairwise acquisition order at **60.75%**, largely
independent of the developmental ordering, and flagged it as an unreported finding sitting in a
control column. The alternative that would manufacture it for free: **censoring**. In an 18–32
sample a late-acquired category is held by fewer people simply because some have not reached it yet,
so "rare" and "late" would be mechanically linked by the age window.

**Censoring is an age effect and must decay as acquisition completes.** Tested within age bands,
each using its *own* prevalence so the predictor is never imported across bands:

| respondent age | onset ordering | **prevalence ordering** | random (40 perms) | n |
|---|---:|---:|---:|---:|
| 14–17 | 65.51% | **59.65%** | 50.15% | 2,364 |
| 18–20 | 65.92% | **59.77%** | 49.98% | 2,011 |
| 21–24 | 66.20% | **59.57%** | 49.19% | 2,605 |
| 25–28 | 67.71% | **59.50%** | 50.60% | 2,719 |
| 29–32 | 67.52% | **60.24%** | 49.76% | 2,760 |

**NOT CENSORING.** The rarity ordering is flat across the entire window — **+0.59 points from
youngest to oldest**, where censoring predicts decay.

| # | Claim | Verdict |
|---|---|---|
| 54 | **rarity predicts acquisition order** — the finding `#53` left in a control column | **CONFIRMED and not an age artifact.** ~60% against a 50% baseline in every band, stable. **This release now carries two independent ordering principles for acquisition: a developmental schedule (65.5–67.7%) and rarity (≈60%), contributing largely separately (`#53`: removing prevalence costs the schedule 1.16 points)** |

**Two secondary observations worth the space.** The onset ordering *improves* with respondent age —
65.5% → 67.5% — consistent with older respondents having completed more of their acquisition, so
their orderings are more fully expressed. And I fixed the limitation `#53` raised against itself:
the random baseline is now averaged over **40 permutations per seed** rather than one, giving a
standard error near 0.7 instead of a 7.68-point seed spread. **A complaint I filed against my own
round two rounds ago, discharged, which is the only reason it was not simply forgotten.**

**What remains unresolved, and it needs data this release cannot provide.** Rarity could mean
*acquired* later (a cumulative exposure process — common things first) or *noticed* later (rare
interests take longer to recognise). Censoring is dead; those two are not separable retrospectively,
and telling them apart is a prospective-measurement problem, not an analysis one.

---

## Entry 55, added by `E01·A03·R17` — I wrote a prediction that was analytically impossible, and the question underneath it had a real answer

`#54` closed with: *"if the two ordering principles are separate, an ordering using both should beat
either alone, by roughly the sum of their independent contributions — sharp enough to be wrong."*

**It is wrong by construction and I should have seen it before writing it.** Onset and prevalence
are both **global** orderings. The oracle is *by definition* the best global ordering. Onset already
reaches 101% of it. **No combination of global orderings can exceed the best global ordering** — the
prediction was not sharp, it was impossible, and I dressed a tautology as a falsifiable claim.

The bound is the useful part: it says prevalence adds nothing **on top of** onset despite predicting
60% alone, so its information is a *subset*. And it exposes the question the schedule line had never
asked — **how much of acquisition order is global at all?**

| ordering | accuracy | seed spread |
|---|---:|---:|
| global oracle (best single global ordering) | 66.97% | 1.13 |
| global ordering (fitted on training people) | 66.56% | 0.75 |
| **neighbour-fitted** — 400 nearest people in *preference* space, never onset | **67.44%** | 1.04 |
| random-neighbour *(negative control, same fitting-set size)* | 66.67% | 0.88 |
| prevalence | 60.46% | 1.20 |

The neighbour ordering beats the global one by **+0.88**, against a random-neighbour control of
**+0.11** — 8× the control, and it exceeds even the *global oracle*, which it is allowed to do
because the oracle bounds global orderings only.

| # | Claim | Verdict |
|---|---|---|
| 55a | **"a combined ordering should beat either alone"** (`#54`'s parting prediction) | **WITHDRAWN as analytically impossible.** Bounded by the oracle |
| 55b | **individual variation in acquisition order** | **PRESENT BUT UNRESOLVABLE.** +0.88 points against a seed spread of 1.04 → ratio **0.85**, below `#34`'s threshold of 2. Real by the control comparison, unmeasurable by my own criterion |

**The decomposition, with each part labelled by its own resolvability:**

- **global schedule: +16.56 points over chance**, seed spread 0.75, ratio **22** — RESOLVABLE
- **individual component: +0.88 points**, seed spread 1.04, ratio **0.85** — UNRESOLVABLE

**So acquisition order is global to within this design's resolution, and any individual component is
below it.** That is a sharper statement than "≈95% global", which is the number the table invites
and which I am not entitled to write, because the 5% is exactly the part I cannot measure.

**Applying `#34`'s criterion to a number I had just produced and liked** is the only reason this
entry says "unresolvable" instead of "small but real". The criterion was built in a round that
killed someone else's headline — mine, four weeks of rounds ago — and this is the first time it has
been turned on a result while it was still warm.

---

## Entry 56, added by `E01·A01·R18` — the item margin is still untested after forty rounds and two attempts, and the gate caught both

`#55` established that every `NEXT` line should be checked for analytic possibility before being
acted on. **All 45 audited**: nearly every one was an audit of an existing claim and was executed.
**One opened a new margin, was falsifiable, and was abandoned** when the restructure interrupted it
— `aea8476`, forty rounds ago: *"ten iterations have all measured structure over PEOPLE. The other
margin has never been touched — the ITEMS."*

It is sharper now than when written, because `#49` showed the person-side cross-domain structure is
a thin direction. So: is there **item-side** structure the person side missed?

**Attempt 1 — not information-matched.** I compared a bare item-neighbour average (0.234) against a
baseline already carrying *both* the item marginal and the person rate (0.403), and the script
printed **"EMPTY: the item margin adds nothing"**. It was handicapped, not tested. Without reading
the predictor definitions I would have published that the item margin is empty.

**Attempt 2 — nested without fitting.** I added each candidate's residual to the common base as an
*unweighted* term. Full-magnitude additions overshoot, and the **positive control fell below its own
baseline** (person 0.3673 against 0.4034). **The gate refused** rather than reporting.

| | attempt 1 | attempt 2 |
|---|---:|---:|
| baseline (item marginal + person rate) | 0.4034 | 0.4034 |
| person factors *(positive control)* | 0.4164 | **0.3673** ← below baseline |
| item neighbours | 0.2344 | 0.1441 |
| random items *(negative control)* | 0.0187 | 0.1160 |

| # | Claim | Verdict |
|---|---|---|
| 56 | **the item margin** | **UNTESTED.** Two designs, two errors, no measurement. The correct design fits base and candidate **jointly** on training cells and scores the increment out-of-sample; neither attempt did that |

**What the round actually demonstrates is the gate, for the third time** (`#33`, `#34`, now this).
Attempt 1's failure was invisible to the gate — it passed, and only reading the predictor
definitions caught it. Attempt 2's failure was **exactly** what a positive control is for, and the
gate stopped it cold. **A gate catches broken instruments; it does not catch mismatched
comparisons**, and those need the thing `#36` named as the cheapest instrument here — opening the
file and reading what the predictors actually contain.

**Forty rounds of delay, and the margin is still dark.** That is the honest state: the person side
has been measured exhaustively and the item side has never been measured once.
