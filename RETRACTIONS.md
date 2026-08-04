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

## Entry 18, added by `E01·A02·R038`–`R15` — the count was a property of how many I asked for

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

## Entry 19, added by `E01·A07·R102` — the one K I picked was the only K where the ordering came out the way I reported it

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

## Entry 20, added by `E01·A07·R103` — the coverage law was one point, and the point was 7.6× the rest of the range

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

## Entry 21, added by `E01·A03·R057` — my correction was as unfit as the thing it corrected, and my own positive control is what said so

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

## Entry 22, added by `E01·A03·R058` — the first UNVERIFIED to be resolved, and it went the way the original round said

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

## Entry 23, added by `E01·A05·R093` — the ceiling was real, the number behind it was not what I forecast

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

## Entry 24, added by `E01·A02·R043` — the flagship claim, measured with a validated ruler instead of a guessed one

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

## Entry 25, added by `E01·A01·R013` — the foundation, attacked for the first time, survives at half size

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

## Entry 26, added by `E01·A05·R094` — a pre-registered kill that passed while the real problem sat next to it

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

## Entry 27, added by `E01·A01·R014` — the erotic covariates do nothing, and neither does personality

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

## Entry 28, added by `E01·A02·R044` — the third automated kill I have had to override, and they all fail the same way

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

## Entry 29, added by `E01·A02·R045` — my own "not explained by the arithmetic" flag was 78% arithmetic

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

## Entry 30, added by `E01·A02·R046` — the sharpest fact in the project was one cell of a specification curve

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

## Entry 31, added by `E01·A03·R059` — the first curve that came back robust, and the published cell was the most conservative one by accident

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

## Entry 32, added by `E01·A03·R060` — the exec graph, and the first time the gate made me fix an instrument instead of override a verdict

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

## Entry 33, added by `E01·A07·R104` — three design errors in one round, and the gate caught all three

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

## Entry 34, added by `E01·A07·R105` — the modality deficit is smaller than the noise of the only instrument that can measure it

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

## Entry 35, added by `E01·A08·R107` — the resolvability criterion applied backwards, and modality turns out to be the exception

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

## Entry 36, added by `E01·A08·R108` — the audit that flagged a discrepancy was measuring a different quantity, and two ledger entries were wrong

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

## Entry 37, added by `E01·A08·R109` — the headline reproduces and is resolvable; the positive control I wrote for it is degenerate

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

## Entry 38, added by `E01·A08·R110` — a graded control, and the pipeline reproduces a similarity ordering it was never told

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

## Entry 39, added by `E01·A02·R047` — the SUBSTANCE axis does not exist outside the seven blocks it was built in

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

## Entry 40, added by `E01·A01·R015` — the same check that scoped SUBSTANCE clears the central claim

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

## Entry 41, added by `E01·A02·R048` — I set a gate threshold above the published magnitude of the thing I was gating on

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

## Entry 42, added by `E01·A02·R049` — the last A02 survivor passes every statistical attack and is unidentifiable anyway

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

## Entry 43, added by `E01·A03·R061` — the null that survived having its instrument replaced

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

## Entry 44, added by `E01·A01·R016` — a second string-derived proxy beaten by its own sham

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

## Entry 45, added by `E01·A07·R106` — matching buys identification and spends resolvability, and the published number is on the wrong side of the trade

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

## Entry 46, added by `E01·A05·R095` — the measure cannot tell a concentrated population from a random one, and its published value is smaller than its own bias

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

## Entry 47, added by `E01·A05·R096` — the same control that killed one measure validates the other, and unity means opposite things in the two rounds

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

## Entry 48, added by `E01·A04·R083` — matching corrects where there is something to correct, and my summary statistic divided noise by noise everywhere else

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

## Entry 49, added by `E01·A01·R017` — the shared grammar is real, graded, and carries almost no predictable variance

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

## Entry 50, added by `E01·A03·R062` — the onset RSA survives the framing swap that broke the central claim, and my fifth mis-specified gate nearly hid it

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

## Entry 51, added by `E01·A08·R111` — the published statistics rank the claims in a different order than predictable variance does

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

## Entry 52, added by `E01·A03·R063` — the sort completes, and the claim that carries the most is the one that was hardest to break

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

## Entry 53, added by `E01·A03·R064` — the schedule is developmental, and rarity is a real but separate ordering

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

## Entry 54, added by `E01·A03·R065` — rarity is a second ordering principle, and it is not the age window

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

## Entry 55, added by `E01·A03·R066` — I wrote a prediction that was analytically impossible, and the question underneath it had a real answer

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

## Entry 56, added by `E01·A01·R018` — the item margin is still untested after forty rounds and two attempts, and the gate caught both

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

---

## Entry 57, added by `E01·A01·R019` — the item margin is real, carries 71% of the person margin, and was never measured in fifty-seven rounds

Third attempt, built the way the first two should have been: **cell-level** masking so every feature
is computed from unmasked cells only, and all four models **nested on the same base and jointly
fitted** on training cells, scored on identical held-out cells.

| model | out-of-sample R² | seed spread | **increment over base** | resolvability |
|---|---:|---:|---:|---:|
| base — item marginal + person rate | +0.3131 | 0.0024 | — | — |
| **person factors** *(positive control)* | +0.3420 | 0.0024 | **+0.0289** | ratio 12 |
| **item neighbours** | +0.3336 | 0.0048 | **+0.0206** | ratio 4.3 |
| random neighbours *(negative control)* | +0.3137 | 0.0030 | +0.0006 | — |

**Both gates pass**: the established person margin reproduces (+0.0289), and random neighbours sit
at +0.0006. **578,989 held-out cells.**

| # | Claim | Verdict |
|---|---|---|
| 57 | **the item margin** | **REAL, at 71% of the person margin.** Both increments are resolvable by `#34`'s criterion. `#56`'s attempt-1 verdict of **"EMPTY"** was badly wrong — it came from a handicapped comparison, and the margin it dismissed is nearly three-quarters the size of the one this project spent fifty-seven rounds measuring |

**What is measured and what is not.** Each margin was tested **separately against the base**, so
71% means "item alone recovers 71% of what person alone recovers". It does **not** mean the item
margin adds 71% *on top of* the person margin — two margins can each be large and be the same
structure seen from two sides. **Fitting them jointly is the obvious next measurement and this round
did not do it**, so the redundancy question is open and I am not implying an answer.

**The cost of the delay is the entry.** `aea8476` flagged this margin forty rounds ago; the
restructure interrupted it; two later attempts failed on design rather than on data. In the interim
the person margin was swept for K, for block inclusion, for seeds, for coverage matching, for
framings — **five kinds of audit on the measured margin while the unmeasured one sat untouched.**
`#39` named this failure mode as *availability becoming representativeness*; this is the same error
one level up, where **the margin I had tooling for became the margin that existed.**

---

## Entry 58, added by `E01·A01·R020` — two-thirds the same structure, and the third that is not is below resolution

`#57` measured the item margin for the first time and left the decisive question open: does it add
anything **on top of** the person margin, or is it the same structure seen from the other side? Both
outcomes were possible here — unlike `#54`'s parting prediction, which an oracle bound made
impossible.

| model | out-of-sample R² | increment over base | seed spread |
|---|---:|---:|---:|
| base + **person** | +0.3422 | **+0.0293** | 0.0034 |
| base + **item** | +0.3331 | **+0.0202** | 0.0042 |
| **base + both** | **+0.3488** | **+0.0359** | 0.0040 |
| base + random neighbours *(neg. control)* | +0.3130 | +0.0001 | 0.0033 |
| base + person + random *(second neg. control)* | +0.3423 | +0.0294 | 0.0034 |

All three gates pass, including the second negative control: adding random neighbours **to the
person margin** changes nothing (+0.0294 vs +0.0293), so the joint fit is not inflating.

**PARTIAL OVERLAP — 67% shared.** Combined **+0.0359** = the larger margin **+0.0293** plus
**+0.0066**, where the item margin alone is +0.0202. **The item side contributes a third of itself
on top of the person side; two-thirds of it is the same structure.**

| # | Claim | Verdict |
|---|---|---|
| 58 | **person and item margins** | **SUBSTANTIALLY ONE STRUCTURE.** 67% of the item margin's explained variance is already in the person margin. **The unique item contribution is +0.0066 against a combined seed spread of 0.0040 — ratio 1.65, below `#34`'s threshold of 2, therefore UNRESOLVABLE.** So: mostly the same thing, and whatever is genuinely new on the item side is too small for this design to see |

**Applying `#34` to the number the whole round turns on.** It would have been easy to report
"+0.0066 unique item contribution" as the finding — it is positive, it survived three gates, and it
answers a question I had chased for forty rounds. It is also **1.65 seed spreads**, and the
criterion that killed the modality deficit does not become optional when the result is mine and
recent. **Third time the criterion has been turned on a warm result** (`#52`, `#55`, now this), and
the first time it has landed on the number a round was designed to produce.

**What the forty-round detour bought.** The item margin was worth opening: it is real (+0.0202,
ratio 4.8, well resolvable) and it was dismissed as "empty" by a broken comparison two rounds ago.
What it was not, is *new*. **The structure this project has been describing from the person side is
the structure — the other margin is a second view of it, not a second half of it.**

---

## Entry 59, added by `E01·A01·R021` — the aggregate cross-domain signal is larger than block-internal structure, and my verdict label confused largeness with sameness

`#58` closed on an apparent contradiction: within-block structure predicts at ratio 12 while `#49`
found cross-domain transfer carries no predictable variance. **Reading the code resolved half of it
immediately** — `R19`/`R20` fit the SVD *inside* each block, so their "person margin" was always
block-internal, never cross-block. I had quoted two different objects in one breath.

The decomposition, with cross-block factors fitted on **other blocks only**:

| model | out-of-sample R² | **increment** | seed spread |
|---|---:|---:|---:|
| base + within-block reconstruction | +0.3421 | **+0.0290** | 0.0034 |
| **base + cross-block factors** | +0.3540 | **+0.0409** | 0.0034 |
| **base + both** | **+0.3736** | **+0.0606** | 0.0032 |
| base + cross-block, person labels permuted *(neg. control)* | +0.3131 | **+0.0000** | 0.0042 |
| base + random neighbours *(neg. control)* | +0.3131 | **+0.0000** | 0.0042 |

Both nulls land at **exactly 0.0000**, which is the cleanest control pair in the project.

**The aggregate cross-domain signal is 141% of block-internal structure.** Factors built from 31
*other* blocks predict a person's endorsements in a held-out block **better than that block's own
low-rank structure does**. This is not in tension with `#49` — it is the reconciliation `#49` itself
proposed: pairwise block→block is ~0, and the signal **accumulates across blocks**. Now quantified:
one block gives nothing, thirty-one give +0.0409.

**And my verdict label was wrong.** The gate fired "ONE STRUCTURE" because cross-block exceeded 50%
of within-block — but **largeness is not sameness**. The overlap test is combined against the sum:
0.0290 + 0.0409 = 0.0699, combined = **0.0606**, so they are **87% additive** — two *largely
independent* contributions with ~13% shared. Seventh time a threshold has measured one thing and
labelled another (`#28` was the same error: comparing a component to another component instead of to
its own null).

| # | Claim | Verdict |
|---|---|---|
| 59 | **within-block vs cross-domain structure** | **TWO LARGELY INDEPENDENT CONTRIBUTIONS**, 87% additive, and the cross-domain one is the **larger** (+0.0409 vs +0.0290). The apparent contradiction with `#49` was a pairwise-vs-aggregate confusion in my own quoting, not in the data |

**This is the first round in the project where the cross-domain claim comes out stronger than
expected.** Forty-nine entries have shrunk, scoped or withdrawn it — thin direction, K-dependent
count, unresolvable modality, one-point coverage law. Measured as an *aggregate* against a
*marginals* base with two exactly-null controls, it is the largest single structural increment
recorded here.

---

## Entry 60, added by `E01·A01·R022` — the accumulation is √n, and that is the mechanism as well as the curve

`#59` showed pairwise block→block prediction is ~0 while 31 blocks give +0.0409, so the signal
accumulates. The curve had never been measured, and it is the one quantity a phase-1 collection
would need to size itself.

| source domains | increment | permuted-label null |
|---:|---:|---:|
| 1 | +0.0066 | 0.0000 |
| 2 | +0.0102 | 0.0000 |
| 4 | +0.0135 | 0.0000 |
| 8 | +0.0225 | 0.0000 |
| 16 | +0.0297 | 0.0000 |
| 31 | **+0.0415** | 0.0000 |

Gate passes both ways: n=31 reproduces `#59`'s +0.0409, and the permuted-label null is **exactly
0.0000 at every n** — a control that holds at one n is not a control, and this one holds at six.

**The functional form is √n, and it is not close:**

| fit | coefficient of variation across the six points |
|---|---:|
| **increment ∝ √n** | **6.4%** |
| increment ∝ log₂(n+1) | 11.2% |
| increment ∝ n | 51.9% |

**increment = 0.00723 × √(source domains).**

**The form is the mechanism, not just a curve.** √n is the signature of **one shared latent measured
with independent per-block noise**: averaging n noisy estimates cuts the noise by √n, so recoverable
signal grows as √n. That reconciles `#49` and `#59` exactly — a single block's estimate is buried in
its own noise (pairwise ≈ 0), and the latent emerges only as blocks accumulate. **The project spent
fifty rounds arguing about whether the shared structure was real; its accumulation exponent says
what kind of thing it is.**

| # | Claim | Verdict |
|---|---|---|
| 60 | **the cross-block signal accumulates as √n**, coefficient 0.00723 | **MEASURED.** Not saturated at 31 — still rising 40% from n=16 to n=31 — so this release cannot show a ceiling and I am not claiming one |

**What it costs to go further, stated as arithmetic rather than ambition:** reaching an increment of
0.06 needs **69** source domains, 0.08 needs **122**, 0.10 needs **191**. That is the price list for a
phase-1 collection, and it comes from a law fitted over n=1–31. **Extrapolating it to 122 assumes
the law holds four times beyond its measured range, which is an assumption, not a finding.**

---

## Entry 61, added by `E01·A01·R023` — blocks are not interchangeable, so √n does not mean what I said it meant one round ago

`#60` fitted `increment = 0.00723 × √n` at CV 6.4% and read √n as **the signature of one shared
latent measured with independent per-block noise** — "the form is the mechanism, not just a curve".
That reading makes a prediction inside this release: if every block is an independent noisy estimate
of the same thing, **blocks are interchangeable** and only *how many* matters, not *which*.

Fourteen random subsets of n=8, three seeds each. Both gates pass — the mean reproduces `#60`'s
n=8 value (+0.0222 vs +0.0225) and the permuted-label null is −0.00000 everywhere.

| | |
|---|---:|
| SD **across** subsets (which blocks) | **0.00376** |
| SD **within** a subset, across seeds | **0.00049** |
| **variance ratio** | **57.9×** |

**NOT INTERCHANGEABLE**, by a factor of 58. And the property that predicts which subsets do better
is not subtle:

| subset property | correlation with its increment |
|---|---:|
| **mean respondent count of the source blocks** | **+0.816** |
| mean option count | +0.103 |
| fluid-family share | +0.000 |

| # | Claim | Verdict |
|---|---|---|
| 61 | **"√n is the signature of one shared latent with independent per-block noise"** — `#60`'s mechanism reading, one round old | **UNSUPPORTED.** Blocks differ 58× more than seeds do, and **block size predicts subset quality at r = +0.816**. `#60` sampled subsets at random, so total respondents grew with n — **√(block count) and √(total sample) were perfectly confounded in that design**, and the same curve is produced by "more data helps". The √n **fit** stands as a description; the **mechanism** does not |

**The prediction I wrote was the right one and it cost one round to run.** `#60` closed by naming
exactly this test — *"any block subset that beats its own n would falsify the independence
assumption"* — and it did, immediately. **The difference between that and `#54`'s parting prediction
(bounded by an oracle, impossible to fail) is the difference `#55` was about**, and this is the first
NEXT line written after that lesson that turned out to be both falsifiable and false.

**What would settle the mechanism:** vary block count and total respondents *independently* — e.g.
subsample large blocks down to small-block size, so n rises while total N is held fixed. If the
increment still grows as √n, the per-block-latent reading survives. If it tracks √(total N), the
structure is one pooled estimate and block boundaries are incidental to it.

---

## Entry 62, added by `E01·A03·R067` — three levels of sub-global structure tested, all three below resolution

`#52` showed one global ordering saturates what any global ordering can do. `#55` found the
**individual** component is +0.88 against a 1.04 spread — unresolvable. Between "one for everyone"
and "one per person" sits the level never tested: **groups**, which should be far more detectable
because a group ordering is estimated from thousands of people.

| grouping | own-group ordering | other-group ordering | difference | seed spread | ratio |
|---|---:|---:|---:|---:|---:|
| **sex** | 67.47% | 65.96% | **+1.51** | 3.77 | **0.40** |
| **breadth** (median split) | 67.67% | 67.53% | **+0.14** | 5.56 | **0.02** |
| random ordering *(neg. control)* | 49.77–49.81% | — | — | ~3 | — |

Gates pass: own-group beats random by **>17 points** in every group, and the random ordering sits at
chance.

| # | Claim | Verdict |
|---|---|---|
| 62 | **group-specific acquisition schedules** | **NOT DETECTED, and this is a LOW-POWER null.** Within-group test sets are smaller, so seed spreads are 3.8–5.6 points and the MDE is roughly **7.5 points**. The observed sex difference is +1.5. **A difference of five points would have been invisible here**, and I am not reporting "one schedule" as though the design could have said otherwise |

**Three levels, three sub-resolution results:** individual (+0.88, ratio 0.85), group-by-sex (+1.51,
ratio 0.40), group-by-breadth (+0.14, ratio 0.02). **Everything below the global schedule lands
under this design's resolution**, and the global schedule itself sits 16.6 points above chance with
a ratio of 22.

**The honest shape of that.** It is not "acquisition order is global" — it is **"acquisition order is
global at the only resolution this release supports, and every finer structure tested has been
smaller than the noise of testing it"**. Those differ in what they license: the first invites a
mechanism, the second invites a better instrument. Given `#34`, `#52`, `#55` and now this, the
project has repeatedly found the same shape — one large well-measured global effect, and a fringe of
sub-resolution structure underneath it that a finer-binned release would be needed to see.

---

## Entry 63, added by `E01·A03·R068` — 66.7% was accuracy on the pairs the release can order; over all pairs it is bounded [60.5, 66.5]

`#52` reported the schedule at 66.71% and flagged, honestly but without quantifying it, that **36.3%
of within-person pairs are tied** by the 2-year binning and were excluded. Excluded pairs are not a
random third — a tie means both onsets fell in the same bin, so **ties are the pairs with the
smallest true gaps**, the hard ones.

| true onset gap | accuracy | random | share of all pairs |
|---|---:|---:|---:|
| **0 (tied — excluded)** | — | — | **36.3%** |
| 0.1 – 2.5 yr | **61.85%** | 50.39% | 25.0% |
| 2.5 – 4.5 yr | 66.68% | 50.82% | 16.9% |
| 4.5 – 8.5 yr | 70.85% | 51.81% | 15.5% |
| 8.5+ yr | **73.81%** | 49.74% | 6.4% |

Gates pass: every non-tied bin is above chance (min 61.9%) and the random ordering sits at chance in
all of them. **Accuracy is gap-driven — 12 points from the narrowest to the widest bin.**

| # | Claim | Verdict |
|---|---|---|
| 63 | **the schedule predicts held-out acquisition order at 66.7%** | **RESTATED AS BOUNDS: [60.5%, 66.5%].** The published figure is accuracy on the **63.7%** of pairs the binning can order. Scoring ties at half credit — what a predictor earns when it cannot break them — gives **60.52%** over all pairs |

**Why bounds and not a point.** A tie means both onsets landed in the same 2-year bin, but the true
order may still exist beneath the binning. If the ordering does better than chance on those sub-bin
gaps, the truth is above 60.5%; if the bins hide genuine simultaneity, it is near 60.5%. **66.5% is
the upper bound (ties don't count), 60.5% the lower (ties count as coin flips), and this release
cannot narrow the interval** — the information needed is exactly what the binning destroyed.

**The finding survives the restatement, and is smaller than advertised.** Every gap bin beats
chance, the gradient is orderly, and the schedule remains the strongest claim in the project. But
the number that has been quoted for eleven entries was measured on the subset where the measurement
is easiest, and I flagged the tie share in `#52` without ever asking what it cost. **A caveat noted
and not quantified is a caveat that has been filed, not paid.**

---

## Entry 64, added by `E01·A01·R024` — at equal data, fragmenting it across more blocks makes the signal worse, so the accumulation was never about blocks

`#60` fitted `increment = 0.00723 × √n_sources` and read it as one shared latent measured with
independent per-block noise. `#61` falsified the interchangeability that reading requires (subset
variance 58× seed variance; block size predicts subset quality at r=+0.816) and identified the
confound: random subsets grow in **total respondents** as n grows, so √(block count) and √(total
sample) were never separated.

The confound is removable by subsampling. **Hold total respondent-rows fixed; vary how many blocks
they are spread over.**

| budget (total rows) | n=8 | n=16 | n=31 |
|---:|---:|---:|---:|
| 2,000 | 0.0015 | — | — |
| 4,000 | 0.0025 | — | — |
| **≈6,400** | **0.0042** | **0.0028** | **0.0025** |

**Positive control passes** — at fixed n=8 the increment rises monotonically with budget
(0.0015 → 0.0025 → 0.0042). **Negative control passes** — permuted-label null ≤ |0.00001| in every
cell.

**At equal budget, MORE BLOCKS IS WORSE.** n=31 delivers **60%** of what n=8 delivers from the same
number of respondent-rows.

| # | Claim | Verdict |
|---|---|---|
| 64 | **"the cross-block signal accumulates with the number of source domains"** — the reading behind `#60`'s √n law | **REVERSED.** Accumulation tracks **total sample**, and block count is not merely incidental (`#61`'s hypothesis) but **actively costly**: fragmenting fixed data across more blocks degrades the estimate by 40%. The structure is **one pooled estimate**, and block boundaries are a tax on measuring it |

**What `#60`'s curve actually was.** Random subsets of n blocks contain ≈n × (mean block size)
respondents, so its x-axis was total sample in disguise. The √ shape is the ordinary √N of an
estimate from N observations — **the most common scaling law there is, and I read it as a statement
about latent structure.**

**Eighth verdict-label failure, and this one is mechanical:** the script printed `partial: nanx`
because budget 6,386 (n=31 × 206 rows) did not align with 6,400 in the pivot, so the ratio was
computed across an empty intersection. **A one-cell rounding mismatch, and the automated verdict
came out as `nan` while the table it was computed from is unambiguous.** Read the table, not the
verdict — for the eighth time.

**Scope, and it is real.** These increments (0.0042 at best) are an order of magnitude below the
full-data +0.0415, because every source is capped at ≤800 rows. The comparison at fixed budget is
internally valid; whether the same ordering holds at full scale is **not tested**, and the release
cannot test it — at full block sizes, n and total sample cannot be varied independently.

---

## Entry 65, added by `E01·A09·R112` — the epoch is named after a quantity the loader deletes on line 1

The epoch is `E01_sexual_as_a_value_not_a_category`. That title is a claim about the **relative
size** of the ITEM main effect (content is content) and the PERSON×ITEM interaction (the same
content carries different value for different people). **105 rounds were run inside this epoch and
none of them measured it**, because the shared loader's first two lines are

```python
R = M - M.mean(0, keepdims=True)     # <- the ITEM main effect, deleted
R = R - R.mean(1, keepdims=True)
```

**Every claim in this project is a claim about the interaction, made after its rival had been
removed from the data.** This round is the first to touch the raw uncentred matrix.

Held-out cell masking (15%, 3 seeds), Shapley decomposition over all 6 orderings, 32 blocks:

| K | ITEM | PERSON | INTERACTION | full R² |
|---|---:|---:|---:|---:|
| 1 | +0.1990 | +0.0893 | **+0.0169** | 0.314 |
| 2 | +0.1964 | +0.0924 | −0.0022 | 0.289 |
| 4 | +0.1934 | +0.0969 | −0.0584 | 0.230 |
| 8 | +0.1986 | +0.0984 | −0.1402 | 0.145 |

**The item main effect is real and it is the largest single component.** The within-person shuffle
kills it (share −0.008), so it is not row sums in disguise. It is also flat in K, as a main effect
must be.

**Two of my own gates were mis-specified, and one of them inverts the reading.**

| # | Claim | Verdict |
|---|---|---|
| 65a | **"the interaction contributes nothing (share −0.22 at K=4)"** | **WRONG, and the graded control is what says so.** The estimator carries a large negative bias: in a *purely additive* synthetic world with **no interaction at all**, the pipeline returns X_share = **−0.629**. Real data returns −0.22. Read against the dose curve (g=0 → −0.629, 0.25 → −0.552, 0.5 → −0.364, 1.0 → +0.058) the real data sits near **g ≈ 0.75**. A raw negative number here means a *substantial* interaction |
| 65b | **gate (b), "positive control fails at g=0"** — coded as `FLOOR < 0.02` | **A CHECK THAT CANNOT FAIL.** A floor of −0.629 satisfies `< 0.02` trivially. The intent was "returns ~zero at g=0"; what was written passes for any negative number, however large. Built 5×, caught 5× |
| 65c | **gate (d), the placebo** — coded as `abs(contribution) < 0.002` | **WRONG DIRECTION.** A permuted person component is not inert, it is *actively wrong*, so adding it must **hurt**. It contributed −0.136. A placebo *component* passes when its contribution is **≤ 0**, not when it is ≈ 0. Ninth verdict-label failure. The −0.136 is itself evidence: the person effect is genuinely person-specific, not a generic offset |

**Why the bias exists, measured not assumed.** X falls monotonically with K (+0.017 → −0.140) and
correlates with cells-per-parameter at **r = +0.517** across blocks. Rank-K soft-impute on masked
binary cells overfits, and the overfit scales with how few cells each parameter sees. Block 17
(1,367 × 10) returns X = −1.46; block 51 returns +0.38.

**Verdict on the epoch title: `UNVERIFIED`, and that is not an acquittal.** The A-vs-B comparison
requires the interaction on the same scale as the item effect, and the only calibration available
here is **one reference block's** dose curve applied to 32 blocks of wildly different shape — when
the bias is demonstrably shape-dependent. What is licensed today:

- **ITEM main effect: +0.199 held-out R², ~2/3 of everything explained, negative-control clean.**
- **PERSON main effect: +0.089, and person-specific (permuting it costs 0.136).**
- **INTERACTION: present** (real data is far above the additive floor) **but its magnitude is not
  yet on a comparable scale.**

The fix is cheap and specific: **a per-block additive synthetic control**, so every block is read
against a floor built from its own marginals and its own shape. That is `R02`.

---

## Entry 66, added by `E01·A09·R113` — the floor was built from a world nothing like the data, and the identification gate caught it

`R01` could not deliver a verdict because the interaction estimator's negative bias is
shape-dependent. `R02` gave every block its own floor, built from its own marginals — and then
**refused to report a verdict for any of the 32 blocks**, because the nuisance-matching check
failed everywhere.

| check | result |
|---|---|
| dose recovery (planted interaction recovered) | **32/32 PASS** |
| nuisance match (synthetic reproduces the block's own main effects) | **0/32 PASS** |

**Why.** The synthetic drew Bernoulli cells from an additive probability model. That world is far
**weaker** than the data on both main effects:

- item effect **0.157 lower** than real (median; range 0.031–0.396)
- person effect **0.0258 vs 0.0893 — 0.3× the real one**

A floor built from a world unlike the data is not a floor. Bernoulli-from-marginals cannot
reproduce real breadth clustering, and clipping probabilities to [0.02, 0.98] compresses the item
effect further.

| # | Claim | Verdict |
|---|---|---|
| 66 | **"a per-block additive synthetic gives each block its own floor"** | **UNIDENTIFIED — the synthetic is not the same world.** The correction `X_c = X_real − X_synth` is only meaningful if the two worlds differ *in the interaction and nothing else*. Here they differed in both main effects by more than the entire quantity being estimated |

**What this round is worth anyway.** The gate fired for a real reason, not as a formality, and it
fired *before* a verdict was published. `R02`'s uncorrected table shows the gap `X_c − I` negative
in **29 of 32 blocks** at K=1 — which, had the gate been a formality, would have been reported as
world A. **It may still be world A. It is not yet allowed to be.**

Also fixed here: `R01`'s placebo, re-run with the direction corrected. A permuted person component
contributes **−0.142**, which is `≤ 0` and therefore **PASS** — and its size is itself evidence that
the person effect is person-specific rather than a generic offset.

**The fix, and it is the right one for binary data:** fixed-margin randomisation (curveball)
preserves every row sum and every column sum **exactly**, so both main effects are matched *by
construction* rather than in expectation, and it destroys the interaction and nothing else. It also
buys a graded control `R02` could not have — running only a *fraction* of the mixing trades gives a
dose axis with exact margins at every point. That is `R03`.

---

## Entry 67, added by `E01·A09·R114` — the epoch is named after the smaller of the two components, and the loader deleted the larger one first

Fixed-margin randomisation (curveball) preserves **every row sum and every column sum exactly**
(asserted per draw, not assumed), so both main effects are matched *by construction* and only the
person×item interaction is destroyed. Running a *fraction* of the mixing trades gives a graded
control with exact margins at every dose.

**All gates pass, for the first time in this arc.**

| gate | K=1 | K=2 |
|---|---|---|
| graded control monotone in dose | **PASS** `[+0.021, +0.004, −0.010, −0.029, −0.069]` | **PASS** `[−0.003, −0.029, −0.053, −0.080, −0.132]` |
| not already at the floor at f=0 | **PASS** | **PASS** |
| margins matched, per block (\|dI\| ≤ 0.01) | **23/32** | **23/32** |
| median \|I_real − I_null\| | **0.0049** | **0.0061** |
| median \|P_real − P_null\| | **0.0022** | **0.0034** |

| | K=1 | K=2 |
|---|---:|---:|
| ITEM main effect `I` | **+0.2223** | **+0.2217** |
| interaction, bias-corrected `X_c` | **+0.0671** | **+0.1137** |
| ratio `X_c / I` | 0.302 | 0.513 |
| gap `X_c − I` | **−0.1366** | **−0.1029** |
| 2× seed spread | 0.0184 | 0.0242 |
| blocks with `X_c > 0` | **23/23** | **23/23** |
| blocks with `X_c > I` | 2/23 | 6/23 |

| # | Claim | Verdict |
|---|---|---|
| 67a | **the interaction is real** | **CONFIRMED and resolvable.** `X_c > 0` in **23 of 23** identified blocks, and the gap from the fixed-margin floor is 3.6–5.6× its own seed spread. The naive negative numbers in `R01` were estimator bias, exactly as `#65a` said |
| 67b | **`E01_sexual_as_a_value_not_a_category` — the epoch's own title** | **FALSE AS STATED.** The item main effect exceeds the bias-corrected interaction by −0.137 (K=1, **3.3×**) and −0.103 (K=2, **1.9×**), both **5–7× the seed spread**. It is not "a value, *not* a category". It is **both, with the category component the larger one** |
| 67c | **the 105 rounds preceding this arc** | **SCOPED, not retracted.** Every one of them ran on `R = M − M.mean(0) − M.mean(1)`, i.e. on the smaller component, after the larger had been deleted on line 1. Their findings stand *as findings about the interaction*. What none of them may claim is that the interaction is what this survey is mostly made of |

**The honest limit, and it is the live threat.** `X_c` is **K-dependent and rising** (0.067 → 0.114
from K=1 to K=2) because the estimator's overfit penalty grows faster on the null than on the real
matrix. The ordering holds at both K tested, by 4–7× the spread — but **an extrapolation to higher
K is not licensed by two points**, and if `X_c` kept rising at that rate it would cross `I` near
K≈4. That is not a caveat to file; it is the next round.

**Why 9 blocks fail the margin match when curveball preserves margins exactly.** They do not — the
per-draw assertion on row and column sums passed for all 32. `dI` is **sampling noise in the
masked-cell estimate** of the column means, and it is larger in small blocks. Dropping those 9 is
conservative: it removes noisy blocks, never adds favourable ones.

**The sentence I can no longer write:** *"in this data, sexual interest is a value assigned to
ordinary content rather than a content category"* — the content category is the bigger half.

---

## Entry 68, added by `E01·A09·R115` — the K-trend turns over, and the epoch title is false under every reading except a tie

`#67` left the ordering resting on two ranks with the corrected interaction still rising. The full
sweep, with **K=0 in the grid** so that "no interaction term at all" could win:

| K | full R² | ITEM | PERSON | INTERACTION | `X_null` | `X_c` |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.2911 | 0.2118 | 0.0793 | 0.0000 | 0.0000 | 0.0000 |
| **1** | **0.3124** | 0.2077 | 0.0871 | **+0.0177** | −0.0695 | 0.0872 |
| 2 | 0.2901 | 0.2050 | 0.0919 | −0.0069 | −0.1374 | 0.1305 |
| 3 | 0.2661 | 0.2031 | 0.0939 | −0.0309 | −0.1943 | 0.1634 |
| 4 | 0.2315 | 0.2014 | 0.0952 | −0.0651 | −0.2489 | 0.1838 |
| **6** | 0.1753 | 0.1997 | 0.0976 | −0.1220 | −0.3153 | **0.1932 ← peak** |
| 8 | 0.1415 | 0.2024 | 0.0953 | −0.1563 | −0.3041 | 0.1478 |
| 12 | 0.1164 | 0.2055 | 0.0862 | −0.1754 | −0.2591 | 0.0837 |

**`X_c` peaks at K=6 and turns over.** The extrapolation that threatened `#67` does not happen: the
corrected interaction never runs away, and its maximum over the entire sweep is **+0.1932**.

**Rank selected by held-out prediction** (rule fixed before the sweep was read): K\*=0 for **6 of 23
blocks** — for those, the best-predicting model contains **no interaction term at all**. K\*=1 for
13, and 4 blocks pick K ≥ 2.

| estimand | ITEM | INTERACTION | ratio | gap | 2× spread |
|---|---:|---:|---:|---:|---:|
| **PREDICTION** — what it delivers out of sample | +0.2223 | **+0.0189** | **11.7×** | −0.1969 | 0.0092 |
| **DETECTION** — structure present, bias-corrected | +0.2223 | **+0.0629** | **3.5×** | −0.1496 | 0.0092 |

Both gaps are **16–21× the seed spread**, and they agree. `#67`'s asymmetry (a raw item effect
against a corrected interaction) is now resolved by reporting both, and the verdict is the same
either way.

| # | Claim | Verdict |
|---|---|---|
| 68a | **`X_c` rises with K and may cross the item effect** — the live threat `#67` flagged | **DEAD.** `X_c` peaks at K=6 (+0.1932) and falls thereafter. It never exceeds the item main effect at any rank tested |
| 68b | **`E01_sexual_as_a_value_not_a_category`** | **FALSE AS STATED, under both estimands, at the pre-registered rank.** Item is 3.5–11.7× the interaction |
| 68c | **the most generous possible reading, chosen adversarially in the rival's favour** | at **K=6**, the rank that *maximises* the corrected interaction, the two are **+0.1932 vs +0.1997 — a tie.** Even here the title does not survive: the best case for "a value, **not** a category" is **"both, equally"** |

**The specification curve of the verdict itself.** Under the pre-registered rule: world A by
16–21× the spread. Under the rule most favourable to the epoch title: a tie. **There is no rank,
estimand, or selection rule in this grid under which the interaction is the larger component.**

**The sentence I can no longer write, in its last surviving form:** *"the interaction is at least as
large as the item effect somewhere in the rank sweep."* At its own maximum it reaches parity, and
that maximum costs 44% of the model's held-out R² (0.175 at K=6 vs 0.312 at K=1) to buy.

---

## Entry 69, added by `E01·A10·R116` — the item main effect is 1+1=2, and the algebra was checkable in one line

`A09` established the item main effect at 3.5–11.7× the interaction. Before that becomes a claim
about a **content category** it has to survive the arithmetic trap — *could this number have come
out otherwise?*

**ATTACK 1 (arithmetic).** For a binary matrix the held-out R² of column means has a closed form:
`I_hat = Var(p_j) / (E[p_j(1−p_j)] + Var(p_j))`. Regressing the **measured** `I` on it across the
23 identified blocks:

| | |
|---|---:|
| slope | **0.9878** |
| intercept | +0.0025 |
| **R²** | **0.9939** |
| median \|residual\| | 0.0031 |

**ATTACK 2 (gauge).** Replace each block's prevalences with a surrogate of the **same dispersion,
different shape** (dispersion matched to 0.0011). `I` is unchanged in **23 of 23 blocks** — median
|Δ| = 0.0074 against a 2× spread of 0.0260. **`I` cannot see which options are popular. It sees
only how spread out the prevalences are.**

| # | Claim | Verdict |
|---|---|---|
| 69a | **"the ITEM main effect is the larger component" (`#67b`, `#68b`)** | **NUMERICALLY INTACT, SEMANTICALLY REWORDED.** It is a **DERIVATION**, not a measurement — forced by prevalence dispersion at R² = 0.994, and blind to prevalence shape in 23/23 blocks. It must be stated as *"option base rates vary more than people vary around them"* |
| 69b | **the reading "there is a sexual-content category the survey detects"** | **UNSUPPORTED, and it always was.** In binary endorsement data, *base rate* and *content* are **the same number**. Nothing here distinguishes a content detector from the fact that some options are ticked more often. Model A was never tested by this arc; it was defined into it |
| 69c | **does this rescue the epoch title?** | **NO.** The comparison is unchanged — both quantities are still held-out R² on the same cells, and the item side still wins by 16–21× the spread. What changed is what the winning side *means* |

**What is actually established, stated at the size the evidence supports:** knowing **which option**
it is tells you far more about whether it is endorsed than knowing **who the person** is, once their
overall breadth is accounted for. That is a real and non-trivial answer to the A-vs-B question — and
it is a statement about *magnitude*, not about a mechanism.

**ATTACK 3 (presentation-order primacy) — and my tenth mis-specified gate.** The release
**alphabetises** multi-select answers (order-consistency 1.0000 across 119 pairs), so display order
is destroyed and this test only works under the unverifiable assumption *display == alphabetical*.

The round printed **`positive control FAIL — test is blind`**, and that was **my error, not the
test's**: I compared the planted effect against the **between-block SD (0.2828)** where the
**standard error of the mean (0.0590)** belongs. Corrected: the control at g=0.15 moves ρ to −0.314
vs a threshold of 0.229 → **FIRES**. The instrument is not blind.

Corrected reading of the observation: **ρ = −0.111, SE 0.059, |t| = 1.88** — *below* the 2-SE bar,
so **not resolvable**, but the point estimate has **exactly the sign primacy would produce** and sits
just inside the floor. This is **not a clean null**, and it is recorded as a residual threat with
what it would require: **the survey instrument's option order**, which this release does not carry.

**The sentence I can no longer write:** *"the survey shows sexual interest is organised by content
category."* It shows options differ in prevalence — and in this data that sentence and the previous
one are the same sentence.

---

## Entry 70, added by `E01·A10·R117` — the epoch title was two questions wearing one sentence, and the answer differs between them

`#69` showed the within-block contrast cannot adjudicate model A vs model B, because within a block
*base rate* and *content* are the same number. The level where they differ is **across** blocks:
base rates are block-local — an option in block 1 and an option in block 2 are different options —
so **no item-level quantity can transfer**, while a person-side readout weight can.

Four components, one scale, the same held-out cells, Shapley over all 24 orderings:

| component | K_within=0 | K_within=1 |
|---|---:|---:|
| `I` item prevalence — block-local, a **derivation** (`#69`) | **+0.2248** | +0.2222 |
| `P` person breadth — block-local | +0.0850 | +0.0898 |
| `C` **cross-block, person-side** (factors fit only on the other 31 blocks) | **+0.0222** | +0.0172 |
| `W` within-block interaction | 0 by design | +0.0078 |
| `C` corrected against the person-permutation null | **+0.0244** | +0.0176 |

**Gates.** Person-permutation null at **−0.0009**; `C` exceeds it by more than 2× seed spread in
**23 of 23 blocks** (median 2× spread 0.0048). Placebo passes: `I` and `P` move by |0.0006| and
|0.0013| under the permutation, as they must, since neither uses `C`.

| # | Claim | Verdict |
|---|---|---|
| 70a | **a person's deviation from prevalence transfers across content domains** | **CONFIRMED, 23/23 blocks**, at +0.0244 held-out R² with a clean permutation null. This is the *only* quantity in the release that cannot be a base rate |
| 70b | **`E01_sexual_as_a_value_not_a_category`** | **THE TITLE IS TWO QUESTIONS.** As a *variance* claim it is **FALSE** — prevalence dispersion is 9× the cross-block person signal (`#67`, `#68`, `#69`). As a claim about *what is domain-general* it is **SUPPORTED** — at the cross-block level the item side contributes nothing and the entire transferable signal is person-side |
| 70c | **"the item contributes 0 across blocks"** | **A DERIVATION, and labelled as one.** It is 0 because no option is shared between blocks, not because it was tested. A content-property route would give the item side a cross-block channel; that needs per-option content annotations, which this release does not carry and which `#28`/`#44` showed string proxies cannot substitute for |

**Not compute-matched, and it matters.** `C` uses 4 cross-block factors, `W` uses 1 within-block
factor. The apparent ordering `C (+0.0172) > W (+0.0078)` at K_within=1 is therefore **not
licensed** — it compares 4 parameters against 1. What *is* licensed is that `C` survives the
addition of `W` (0.0244 → 0.0176, a 28% overlap, not an absorption).

**The ontology shift.** The A-vs-B question was never *which is bigger*. Content explains the most
variance and explains it **block by block with no transfer**; the individualised readout explains
far less and is **the only thing that crosses domains**. Both facts are true simultaneously and the
epoch title asserted only one of them while sounding like it asserted both.

**What this means for the 105 rounds.** They ran on the double-centred residual — i.e. on `C` and
`W`, the domain-general and domain-specific person-side structure. That was the **right object for
the A-vs-B question** and the wrong object for a variance claim. `#67c` said their findings stand as
findings about the interaction; `#70` upgrades that: **the interaction is the part of this data that
distinguishes the two models at all.**

---

## Entry 71, added by `E01·A10·R118` — the domain-general structure keeps buying signal; the domain-specific one is one dimension deep and then noise

`#70` refused the ordering `C > W` because `C` used 4 factors and `W` used 1. The fix is not a rank
match but a **parameter count**, and here they are not the same thing:

- `C` at rank k fits **k×m** loadings on the target block. **The person scores are free** — they
  come entirely from the other 31 blocks and no target cell estimates them.
- `W` at rank k fits **k×(n+m)** parameters on the target, both scores and loadings.

Median across blocks: `df_C` = 18–144, `df_W` = 3,249–12,996. **`W` costs 22–800× more.**

| | Kw=0 | Kw=1 | Kw=2 | Kw=4 |
|---|---:|---:|---:|---:|
| **C** at Kc=1 | 0.0132 | 0.0086 | 0.0081 | 0.0086 |
| **C** at Kc=4 | 0.0288 | 0.0194 | 0.0183 | 0.0190 |
| **C** at Kc=8 | **0.0359** | 0.0253 | 0.0237 | 0.0242 |
| **W** (any Kc) | 0 | **+0.0078…+0.0138** | **−0.011…−0.018** | **−0.063…−0.070** |

**`C` corrected against the person-permutation null, and it does not saturate:** 0.0133 (Kc=1) →
0.0199 → 0.0296 → **0.0380 (Kc=8)**, with the null at −0.0001…−0.0021.

| # | Claim | Verdict |
|---|---|---|
| 71a | **domain-general exceeds domain-specific** | **CONFIRMED across the grid: 186 of 276 comparisons (67.4%) general, 32 (11.6%) specific, 58 (21.0%) not distinguishable** — every cell declared only above 2× its own seed spread |
| 71b | **the cell that disagrees, published rather than buried** | at **Kc=1, Kw=1 — the RANK-matched cell — `W` wins**: +0.0152 vs +0.0075, 8 blocks specific against 6 general. **A rank match is the one comparison that favours `W`, and it is the comparison `#70` was about to make.** It buys that win with 200× the parameters |
| 71c | **the shape of the difference, which is the actual finding** | **`W` is exhausted at rank 1 and goes negative at rank 2** (−0.011 → −0.070). **`C` grows monotonically to rank 8 and has not saturated.** The domain-general structure is multi-dimensional; the domain-specific structure is roughly one dimension deep and then noise |

**What is now sayable:** the person-side readout is **not assembled per domain**. A single
cross-domain structure, estimated without ever seeing the target block, predicts that block's cells
better than a structure fit on the block itself — while spending two orders of magnitude fewer
parameters there.

**What is still not sayable:** how many cross-block dimensions there are. `C` is still rising at
Kc=8 and the sweep stops there. That is the same non-saturation `#18` and `#49` hit from other
directions, and it is now the arc's residual gap.

---

## Entry 72, added by `E01·A10·R119` — there is no coordinate count, and that is why two rounds failed to find one

`#18` and `#49` both tried to settle how many domain-general coordinates exist, from *within* the
shared space, and both failed. This round asks by prediction instead: a dimension counts only if it
predicts held-out cells in a block whose data never touched its estimation. And the curve is read
**only against synthetic worlds of known rank**.

**The estimator counts, and it counts exactly.** Gain per dimension:

| world | knee | per-dimension increments |
|---|---:|---|
| control, **true rank 2** | **2** | **0.1277**, 0.0009, 0.0017, 0.0006, 0.0005, 0.0003, 0.0002, 0.00003 |
| control, **true rank 5** | **5** | 0.0753, 0.0721, **0.0774**, 0.0006, 0.0007, 0.0005, −0.0000, 0.00003 |
| **real data** | — | 0.0070, 0.0056, 0.0031, 0.0022, 0.0017, 0.0014, 0.0009, 0.0006 |

**The cliff ratio is the whole result.** At its true rank the r=2 control drops **147×** in one step
and the r=5 control drops **141×**. The sharpest drop anywhere in the real curve is **1.8×**.

| # | Claim | Verdict |
|---|---|---|
| 72a | **the estimator can count dimensions** | **CONFIRMED by dose-response.** Knee at 2 for a rank-2 world, at 5 for a rank-5 world, with a 141–147× cliff in both. This is the positive control `#18` and `#49` never had |
| 72b | **"how many domain-general coordinates are there?"** | **THE QUESTION PRESUPPOSES A CLIFF THAT DOES NOT EXIST.** The real spectrum decays smoothly — a factor of 12 across 32 dimensions with no step exceeding 1.8×. It is not *r factors plus noise*. `#18` and `#49` did not fail from weak instruments; **they were asking for a number the object does not have** |
| 72c | **`#70`'s cross-block magnitude, +0.0244** | **UNDERESTIMATED 2.6×.** It was measured at rank 4. Corrected `C` reaches **+0.0635 at Kc=32** and is still gaining at 1.45× the seed spread. That is **75% of the person-breadth effect (0.085)** and **29% of the item effect (0.222)** — the gap `#70` reported as 9× is nearer **3.5×** |
| 72d | **my saturation criterion** | **MIS-SPECIFIED — the eleventh.** `sat()` used *remaining gain as a fraction of total*, which returned Kc=5 for **both** controls and printed `FAIL — the estimator cannot count`. The correct statistic is gain **per dimension**, since the rank grid is uneven (…3, 5, 8, 12…). The table was unambiguous while the automated verdict was wrong, for the eleventh time |

**The ontology shift.** The project has spent rounds asking *which* coordinates and *how many*. The
object is a **slowly-decaying spectrum**, not a small basis — so "name the axes" was never going to
converge, and the several failures to name them (`#18` K-dependence, `#39` subgraph locality, `#49`
thin direction) are one fact seen three times, not three separate defeats.

**The sentence I can no longer write:** *"the domain-general readout has K coordinates for some
modest K."* Every K up to 32 buys real held-out prediction, and nothing in the sweep says where it
stops.

---

## Entry 73, added by `E01·A03·R069` — the effect estimate never moved; the noise did, and it was a cap I wrote myself

`#55` reported individual variation in acquisition order at **+0.88 points, seed spread 1.03, ratio
0.85** and logged it as a low-power null. `#72` then showed the 8-dimensional person embedding it
used discards most of the available signal. Re-running with the embedding widened and — the larger
fault — **the pair cap removed**:

`#55`'s accuracy function carried `cap=20000` with `if tot>=cap: break`, and that `break` exits the
**person** loop. At ~10 pairs per person it stopped after roughly 2,000 of the available people.

| | `#55` | `#73` |
|---|---:|---:|
| held-out pairs | 20,000 | **77,329** |
| held-out people | ~2,000 | **6,230** |
| neighbour − global | +0.88 | **+0.90** (median over 12 cells) |
| seed spread | 1.03 | **0.29** |
| ratio | 0.85 — unresolvable | **3.1 — resolvable** |

**The point estimate is the same to two decimals.** `#55` measured the effect correctly and could
not see it.

| dim | k=400 | k=1000 | k=2500 |
|---:|---:|---:|---:|
| 8 | +0.730 | +0.849 | +0.598 |
| 32 | **+1.015** | +1.007 | +0.698 |
| 64 | +0.977 | +1.007 | +0.667 |
| 68 (all) | +0.964 | +0.954 | +0.641 |

**Controls.** Random-neighbour gap is clean in **12/12 cells** (−0.13 to −0.46, all inside their own
spread). The graded positive control is monotone — a planted person-specific shift of 0.5 / 1 / 2
years produces gaps of +2.81 / +5.03 / +9.14 — **does not fire at g=0**, and gives an **MDE of 0.5
years**.

| # | Claim | Verdict |
|---|---|---|
| 73a | **`#55` "individual variation in acquisition order is not detectable"** | **OVERTURNED — it was a power failure, not a fact**, and the power was destroyed by a cap in my own code that silently truncated the person loop |
| 73b | **`#62`/`#55`'s reading, "one global acquisition schedule"** | **WITHDRAWN.** A neighbour-fitted ordering beats the global one by **+0.90 points** in 8 of 12 specifications with clean controls throughout. Acquisition order is *mostly* global but **not purely** |
| 73c | **the size, stated so it cannot be over-read** | Global ordering beats chance by **16.4 points** (50 → 66.4). The individual component adds **1.0 more** — about **6% of the orderable signal**, and equivalent to **under 0.3 years** of person-specific shift on the planted scale. Real, resolvable, and small |
| 73d | **`A03`'s "acquisition and valuation are two systems"** | **DOWNGRADED, not overturned.** The individual component is predicted from **preference space** — the valuation side — so the two are not independent. What survives is that they are *mostly* separable; what does not is *strictly* separable |

**Fourth instance of the same failure mode.** `#21`, `#26`, `#40`, `#50` were thresholds mis-set.
This one is different and worse: **a sampling cap that made the design blind, written by me, never
priced, and reported as a property of the world.** A null that comes from an unexamined `break` is
the exact shape of `L11` — untested is not the same as null-survives.

---

## Entry 74, added by `E01·A03·R070` — the cap bound in four rounds and not the fifth, and I had already generalised it to all five

`#73` found `cap=20000` with a `break` that exits the person loop, traced it to five rounds, and
wrote *"every one of them therefore ran on roughly 2,000 of 12,459 eligible people."* That sentence
was a generalisation from one round to four others **without measuring any of them**.

**Measured.** Realized pairs are fewer than attempted pairs, because `acc()` skips tied pairs — and
`#63` measured the tie rate at 36.3%:

| design | held-out people | attempted pairs | realized (×0.64) | cap binds? |
|---|---:|---:|---:|---|
| R14 / R15 / R16 / R17 — full pool | 6,229 | 55,537 | ~35,500 | **YES, ~43% of pairs dropped** |
| R18 — sex group A | 3,187 | 28,275 | ~18,100 | **no** |
| R18 — sex group B | 3,042 | 26,892 | ~17,200 | **no** |
| R18 — breadth high | 3,734 | 36,462 | ~23,300 | marginal |
| R18 — breadth low | 2,495 | 19,154 | ~12,300 | **no** |

The run confirms it directly: **both arms scored 16,903 vs 17,793 pairs over the same 3,115 people**,
and for the sex grouping the capped and uncapped numbers are **identical to three decimals** — shift
0.000, sd ratio 1.000.

| # | Claim | Verdict |
|---|---|---|
| 74a | **`#73`'s "every one of them ran on ~2,000 people"** | **OVERSTATED.** It bound in the four full-pool rounds (R14/R15/R16/R17) and **not** in R18, whose group halves are too small to reach the ceiling. I generalised a measured fact about one design to four I never checked — in an entry whose whole subject was an unpriced assumption about my own code |
| 74b | **`#62`, the group-level null** | **STANDS, and is not a cap artifact.** Re-run uncapped with 6 seeds: sex own−other **+1.183, spread 1.232, ratio 0.96**; breadth **+0.310, spread 0.825, ratio 0.38**. Still below 2× spread. Its low power comes from **small group sizes**, which no amount of uncapping fixes |
| 74c | **my own pair-count estimator, written in this round** | **WRONG by the tie rate.** It counted *attempts*, not *realized* pairs, and over-predicted by 1.57× — precisely `#63`'s 36.3% tie exclusion, a number this project had already measured and I did not apply |
| 74d | **the gate returning `UNVERIFIED`** | **CORRECT.** It required uncapping to shrink the spread, and on R18's design there was nothing to shrink. A gate that refuses when its premise is absent is working |

**What still needs re-pricing:** R14 (the 66.5% schedule headline), R15 (schedule or rarity), R16
(rarity or censoring). Their point estimates are unaffected — the truncation is a random subset, so
location is unbiased — but their **seed spreads are inflated**, which makes every resolvability
verdict in them **conservative**, not permissive. That is the safe direction, and it is the reason
this is a re-pricing rather than a retraction.

**The lesson, and it is the one `#73` was about:** an entry written to record an unpriced assumption
contained a fresh unpriced assumption in its own second sentence. Auditing a class of bug is not the
same as measuring each member of the class.

---

## Entry 75, added by `E01·A03·R071` — CLOSURE: the cap cost precision and cost no verdict

Labelled **Closure**, not Frontier: it protects an existing conclusion rather than separating worlds.
`#74` measured that the cap bound in R14/R15/R16. This re-prices them.

**The cap was binding here** — 20,000 pairs against 35,438 uncapped, over the same 6,230 people —
and the spread shrank **1.42×**, confirming `#74`'s measurement (the positive control this round
needed, and the one R18 could not provide).

| quantity | capped | uncapped | shift | sd capped → uncapped |
|---|---:|---:|---:|---|
| population ordering (**the 66.5% headline**) | 66.885 | **66.852** | −0.033 | 0.289 → **0.191** |
| oracle (best global ordering) | 66.779 | 66.694 | −0.085 | 0.488 → 0.344 |
| prevalence ordering | 60.551 | 60.457 | −0.094 | 0.428 → 0.325 |
| onset with prevalence projected out | 65.717 | 65.880 | +0.163 | 0.218 → 0.363 |

R16's censoring test, all four onset bands, onset vs prevalence:

| band | n | onset | prevalence | gap | ratio (uncapped) |
|---|---:|---:|---:|---:|---:|
| 0–11 | 261 | 55.87 | 51.56 | +4.31 | 2.63 |
| 11–14 | 1,525 | 63.90 | 57.98 | +5.92 | 8.26 |
| 14–18 | 2,904 | 67.74 | 61.06 | +6.68 | 12.56 |
| 18–99 | 1,541 | 70.46 | 64.05 | +6.42 | 8.12 |

| # | Claim | Verdict |
|---|---|---|
| 75a | **verdicts flipped by the cap** | **ZERO of 12.** Every comparison in R14/R15/R16 was resolvable in *both* arms. The unresolvable calls elsewhere in this project were unresolvable for reasons other than the cap |
| 75b | **the schedule headline** | **CONFIRMED at 66.852 ± 0.191** (8 seeds, 35,438 pairs, 6,230 people). `#63`'s bounds [60.5%, 66.5%] stand; the upper bound is now measured with a spread 1.4× tighter |
| 75c | **onset's information is not prevalence** | **CONFIRMED.** Projecting prevalence out of the onset ordering costs **0.97 points** (66.85 → 65.88), and onset beats prevalence in **all four censoring bands** by +4.3 to +6.7, resolvable in every band |

**Why this was worth running even though nothing moved.** `#73` and `#74` established that a cap I
wrote had silently truncated four designs. Leaving that unpriced would have meant every downstream
number carried an unquantified doubt — and an unquantified doubt is indistinguishable from a real
one. It cost 0.03 points of location and bought back 1.4× of precision.

---

## Entry 76, added by `E01·A10·R120` — my own p=0.70 self-overturn fails on survival and lands on magnitude, and no permutation could have caught it

`ADVERSARY_FORECAST.md` block 2, prediction #1 at **p = 0.70**: *the cross-block transfer `C` is
partly the gate, not the person.* Two independent attacks, because a permutation and a projection
fail differently.

| score type | real | perm (free) | perm (count-matched) | perm (pattern-matched) |
|---|---:|---:|---:|---:|
| **raw** | **+0.0290** | −0.0022 | −0.0007 | −0.0017 |
| **gate-free** (entry pattern + demographics projected out) | **+0.0115** | −0.0020 | −0.0018 | −0.0018 |

Projection verified non-vacuous: **17.3%** of person-score variance removed, residual orthogonal to
every entry indicator to **3.6 × 10⁻¹⁶**.

| # | Claim | Verdict |
|---|---|---|
| 76a | **my forecast, "`C` is the gate"** | **WRONG on survival.** `C` beats the strictest null in **22 of 23 blocks** at +0.0115. `#70`, `#71`, `#72` stand |
| 76b | **…and RIGHT on magnitude** | **`C` shrinks 60.5%** — +0.0290 → +0.0115 — once entry pattern and demographics are projected out. `#70`'s +0.0244 and `#72`'s +0.0635 are **inflated by ~60%** and are restated as **≈+0.010 and ≈+0.025** |
| 76c | **the permutation ladder — the control I actually forecast** | **TOOTHLESS.** Free, count-matched and pattern-matched nulls all return ≈ −0.002. Stratifying the permutation changed **nothing**. The entire correction came from the projection. This is `realstat` §G2 exactly: **a permutation null answers *did the pairing matter*, never *why*** — and I had forecast a permutation as the fix |
| 76d | **the pattern-matched arm's strictness** | **OVERSTATED BY ME, in the lenient direction.** **85.3% of people (13,154)** have no entry pattern shared with ≥20 others, so they are permuted freely. That arm is *mostly* the free permutation wearing a stricter name |

**The removed subspace was disproportionately the predictive part** — 17.3% of score variance
carried 60.5% of the transfer. That is not a nuisance being trimmed; it is a large fraction of what
`C` was measuring.

**The ambiguity this round does not resolve, and it changes what `76b` means.** The projection
bundles **entry indicators with demographics** (sex, age, five personality scales, powerlessness).
The A02-era loader always projected demographics out; the `A09`/`A10` scores did **not**. So the 60%
may be survey structure, may be ordinary demographic variance, and this design cannot say which.
Those are very different claims — *"the transfer is confounded with who answers which block"* versus
*"the transfer is partly sex and age"* — and separating them is one more projection.

**What it cost to run.** This was my own highest-probability prediction of self-overturn, named in a
file, with the control specified. Leaving it written down and unrun would have been `#57`'s pattern:
the reason a claim must fail, filed where nobody re-reads it.

---

## Entry 77, added by `E01·A10·R121` — the confound is the survey's own shape, not who answers it

`#76` shrank `C` by projecting out entry pattern **and** demographics together, and could not say
which mattered. Four arms, identical cells and masks, one variable changed at a time:

| projection | variance removed | `C` (mean) | shrinkage | beats free permutation |
|---|---:|---:|---:|---:|
| **none** — the `#70`/`#72` specification | 0% | +0.03785 | — | **23/23** |
| **demographics only** — sex, age, 5 personality scales, powerlessness | 5.9% | +0.02848 | 24.8% | **23/23** |
| **gate only** — entry indicators for the other 31 blocks | 14.1% | +0.02141 | **43.4%** | **22/23** |
| **both** — the `#76` specification | 17.3% | +0.01708 | 54.9% | **22/23** |

All residuals orthogonal to their design matrix to machine precision.

| # | Claim | Verdict |
|---|---|---|
| 77a | **`#76b`'s ambiguity, resolved: gate or demographics?** | **STRUCTURAL.** The gate alone removes **43.4%** of the transfer against demographics' **24.8%** — 1.75× more transfer from 2.4× more variance. Bounds, since the two are collinear and no split is identified: **gate [30.1%, 43.4%], demographics [11.5%, 24.8%]**, overlap 13.3% |
| 77b | **what that means for `#70`/`#71`/`#72`** | **A caveat they must now carry.** Part of what read as a shared cross-domain readout is **shared exposure to the same entry conditions** — people appear together in blocks because they cleared the same parent ratings, and that alone makes their residuals covary. `C` survives it (22/23 blocks) but at **45% of the published magnitude** |
| 77c | **the demographic share** | **Not a threat, but an inconsistency of mine.** The A02-era loader always projected demographics out; the `A09`/`A10` scores did not, so those rounds silently readmitted a nuisance the project had already decided to remove. 24.8% of `C` was that |

**Restated magnitudes, gate-and-demographics-free:**

| was | is |
|---|---|
| `#70` cross-block transfer **+0.0244** | **≈ +0.011** |
| `#72` at rank 32, **+0.0635** | **≈ +0.029** |
| `#72` "75% of the person-breadth effect" | **≈ 34%** |

**The ordering claims are unaffected.** `#71`'s domain-general-beats-domain-specific comparison used
the same scores on both sides, and `#72`'s no-cliff result is about the *shape* of the spectrum, not
its height. What shrinks is the size of the transfer, and the size was never what those entries
turned on.

**Third mis-specified accessor in three rounds** — `D.mode` resolves to `DataFrame.mode()`, so the
per-block table crashed after the main result had printed. Same family as `T.shift` in `#74`: a
column name that collides with a pandas method fails *loudly* here, which is the lucky case. The
unlucky case is a column named `count` or `size` silently returning a method and being truthy.

---

## Entry 78, added by `E01·A10·R122` — a random basis of identical rank and identical cost gets nothing, which is what an accounting argument cannot see

`ADVERSARY_FORECAST` block 2, prediction #4 (p=0.50): *"`#71`'s parameter-count argument is a bad
accounting. `C`'s person scores are called free because they are estimated elsewhere — but they are
estimated on **the same people**."*

The forecast asks for an accounting. An accounting is the wrong instrument, and `#76` had just cost
a round for exactly that — reaching for a proxy when the direct measurement exists. The direct
question: **does `C` survive when the people it is evaluated on are disjoint from the people whose
data built the basis?**

| arm | `C` (per-block median) | `C` > 0 | `W` | eval n |
|---|---:|---:|---:|---:|
| **raw** — the `#71` specification | +0.00870 | 20/23 | −0.038 | 3,908 |
| **person-holdout** — basis from fit-half only, eval-half scores by projection | **+0.00701** | **19/23** | −0.040 | 1,956 |
| **shuffled basis** — same rank, same parameter count, no structure | **+0.00027** | **3/23** | −0.052 | 1,956 |

| # | Claim | Verdict |
|---|---|---|
| 78a | **forecast #4, "the shared basis is unbilled"** | **WRONG.** `C` retains **81%** under person-holdout (+0.00870 → +0.00701) with every evaluation person disjoint from the basis. "Free" was the right accounting for the only reason that matters: **those parameters cannot overfit cells they were never near** |
| 78b | **the control that actually answers it** | **A random orthonormal basis of identical rank and identical parameter count gets +0.00027, positive in 3/23 blocks.** An accounting argument cannot distinguish it from the real basis — they cost the same. A prediction test separates them by **26×**. This is why the forecast's own framing was unanswerable in its own terms |
| 78c | **my pre-registered positive control** | **MIS-SPECIFIED — the twelfth, and this one could not PASS.** It demanded `W > 0` at K=4, and `#71`/`#72` had already measured `W` negative at every K ≥ 2. I wrote a gate whose criterion the instrument was *known* to violate. The mirror of a check that cannot fail, and equally useless |
| 78d | **the corrected control, run rather than asserted** | `W` at **rank 1** on the same eval halves: mean **+0.0214**, positive in **15/23** blocks, resolvable in **12/23**. The eval half **can** fit within-block structure, so the holdout arm is not blind — but at **65% of blocks**, not the 70% I had set, so this is `CONFIRMED at reduced power`, not clean |

**Scope, stated because `78d` requires it.** The person-holdout arm halves the evaluation sample
(3,908 → 1,956). `C` survives that halving; `W` at rank 1 survives it in two-thirds of blocks. A
challenger who wants to break `78a` should attack the power of the holdout arm, not its logic —
and that is a cheaper attack than the one I ran.

---

## Entry 79, added by `E01·A10·R123` — the loss was never swept, and log-loss punishes the interaction rather than rescuing it

`A09`/`A10` swept rank, estimand, null, block, ordering, projection and score type. Every number in
both arcs is **squared error on a binary cell** — an axis `realstat` §G4 requires and this project
had none of.

`ADVERSARY_FORECAST` block 2, prediction #3 (p=0.55): *"under log-loss the interaction's share rises,
because squared error under-weights confident-and-wrong predictions, which is exactly where a
person-specific readout would differ from a base rate."*

**Least-squares estimator, three losses, 23 blocks × 3 seeds:**

| component | Brier | log-loss | L1 |
|---|---:|---:|---:|
| `I` item | **+0.2013** | **+0.1972** | **+0.1863** |
| `P` person | +0.0958 | +0.0775 | +0.1084 |
| `C` cross-block | +0.0065 | +0.0025 | +0.0075 |
| `W` within-block | −0.0662 | **−0.3545** | **+0.0822** |
| full model | +0.2374 | −0.0774 | +0.3844 |

| estimator × loss | item wins | interaction wins | tied |
|---|---:|---:|---:|
| LS × Brier | **23/23** | 0 | 0 |
| LS × log-loss | **23/23** | 0 | 0 |
| LS × L1 | **19/23** | **3** | 1 |

| # | Claim | Verdict |
|---|---|---|
| 79a | **forecast #3, "log-loss raises the interaction's share"** | **WRONG, and wrong in its DIRECTION.** Under log-loss the interaction collapses to **−0.355** while the item effect barely moves (0.201 → 0.197). The *mechanism* the forecast named is real — log-loss punishes confident-and-wrong hardest — but **the overconfident predictor is the low-rank interaction estimator, not the base rate.** A base rate is never confident |
| 79b | **`A09`/`A10`'s ordering** | **SCALE-FREE across the losses tested.** Item wins 23/23 under Brier and log-loss. `#67`, `#68`, `#70` are not artefacts of squared error |
| 79c | **the loss most favourable to the interaction, published because it disagrees** | **L1.** Under absolute error `W` turns **positive (+0.082)**, the item:interaction ratio falls to **2.1×**, and **3 of 23 blocks flip to the interaction**. L1 down-weights exactly the confident-wrong cases log-loss magnifies. The ordering survives, but it is **least secure under L1**, and that is the cell a challenger should attack |
| 79d | **my logistic-estimator arm** | **INVALID, killed by its own negative control.** Person-permuted `C` returns **+0.042** under the logit estimator against −0.0006 under least squares — i.e. **80% of its apparent `C` (+0.053) is null**. The one-step IRLS lets the cross-block term absorb variance through the working weights regardless of what `U` contains. The whole arm is discarded; nothing in `79a`–`79c` uses it |
| 79e | **gate (b), the negative control** | **MIS-SPECIFIED — the thirteenth.** Written as `N[N.est=='ls']`, it checked **one of the two estimator arms** and printed `PASS` while the other arm's null was seven times its own effect. A control that examines part of the design is `#04`'s empty-population failure in miniature |
| 79f | **gate (a), the positive control** | **MIS-SPECIFIED — the fourteenth.** It required `C + W > 0.01` under every loss, bundling `W` — a component `#71`/`#72` had already measured as overfitting at K ≥ 2. Under log-loss `W` = −0.193 dragged the sum negative and the gate printed `FAIL`. **`C` alone is +0.154 to +0.203 in the planted world under all three losses**, so every metric can in fact see an interaction |

**Two mis-specified gates in one round, both from bundling.** `79e` bundled two estimators into one
check; `79f` bundled two components into one threshold. Neither threshold was wrong about its own
quantity — both were applied to a **sum that hides a sign**.

---

## Entry 80, added by `E01·A10·R124` — the metric that flatters the interaction before bias correction is the one that flatters it least after

`#79c` named its own weakest specification and sent this round at it: **L1 is the only loss under
which the within-block interaction `W` is positive**, so if the item-vs-interaction ordering flips
anywhere, it flips there. Rank swept under L1 with a **per-block fixed-margin (curveball) floor**,
margins asserted exact per draw.

| Kw | `W` L1 | L1 floor | **L1 corrected** | `W` Brier | Brier floor | **Brier corrected** |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | +0.0551 | +0.0005 | +0.0546 | +0.0149 | −0.0696 | +0.0845 |
| 2 | +0.0710 | −0.0017 | +0.0727 | −0.0101 | −0.1322 | +0.1221 |
| **4** | **+0.0833** | −0.0126 | **+0.0959** | −0.0648 | −0.2412 | **+0.1765** |
| 8 | +0.0519 | −0.0267 | +0.0786 | −0.1501 | −0.2866 | +0.1365 |

| loss | Kw | item | interaction | gap | blocks item / interaction / tied |
|---|---:|---:|---:|---:|---|
| Brier | 1 | +0.2192 | +0.0206 | +0.2155 | 21 / 1 / 1 |
| Brier | 8 | +0.2146 | −0.1484 | +0.3474 | 23 / 0 / 0 |
| L1 | 1 | +0.1968 | +0.0557 | +0.1525 | 21 / 2 / 0 |
| **L1** | **4** | **+0.1917** | **+0.0931** | **+0.1080** | **20 / 3 / 0** |
| L1 | 8 | +0.1873 | +0.0547 | +0.1244 | 22 / 0 / 1 |

| # | Claim | Verdict |
|---|---|---|
| 80a | **`#79b`, "the ordering is scale-free"** | **STANDS, now with the rank axis crossed against the loss axis.** The item effect wins at **every rank under both losses**. Its narrowest margin anywhere in this project is **L1 at Kw=4: 2.06×, and that gap is still 21.3× its own seed spread**, with 20 of 23 blocks on the item side |
| 80b | **`#79c`, "L1 is the loss most favourable to the interaction"** | **TRUE OF RAW NUMBERS, FALSE AFTER CORRECTION — and this reverses one round after I wrote it.** L1's fixed-margin floor sits near zero (+0.0005 to −0.027) while Brier's is strongly negative (−0.070 to −0.287). Corrected, **Brier reports 1.8× MORE within-block interaction than L1** (+0.1765 vs +0.0959). L1 looked favourable only because it barely penalises the estimator's overfit, so its raw number needed almost no correction |
| 80c | **`W`, the within-block interaction, itself** | **REAL at every rank under both losses.** It clears its own fixed-margin floor in all 8 cells. `#71`'s "exhausted at rank 1 and negative at rank 2" was a statement about the **uncorrected** Brier number; corrected, `W` peaks at **Kw=4** under both losses. The *ordering* in `#71` is unaffected — `C` was corrected there and `W` was not, which made `#71` conservative in the direction it concluded |
| 80d | **`worst.item` → `Series.item()`** | **THIRD pandas accessor collision** after `T.shift` (`#74`) and `D.mode` (`#77`). By `P7`, three of the same bug is infrastructure, not a third patch: **no column in this project may be named after a DataFrame or Series method.** `item`, `mode`, `shift`, `count`, `size`, `min`, `max`, `sum`, `mean`, `std`, `rank`, `pop`, `all`, `any`, `abs`, `where`, `mask`, `first`, `last`, `div`, `pow`, `T` |

**The general lesson in `80b`, which is worth more than the L1 result.** A loss that reports a
*larger raw* effect for a component is not a loss that is *more sensitive* to it — it may simply be
a loss that penalises the estimator less. **The comparison between metrics is only meaningful after
each is referred to its own null**, and this project spent two rounds reading raw cross-metric
numbers before doing that. The same error shape as `#65a`: a raw number compared across
specifications whose floors differ.

---

## Entry 81 — correcting `#80c`, which described `#71`'s internals from memory and got the direction backwards

`#80c` said: *"`C` was corrected there and `W` was not, which made `#71` conservative in the
direction it concluded."* **Both halves are false**, and I wrote them from memory one entry after a
round whose entire subject was reading tables instead of remembering them.

`A10/R118/run.py:152` is `gap = pb.C - pb.W`. **Neither side was corrected.** And the asymmetry that
does exist runs the **opposite way**:

| | its own null | so a raw number is… |
|---|---|---|
| `C` cross-block | person-permutation, **≈ −0.002** (`#80`: L1 +0.00005, Brier −0.0007) | almost unbiased |
| `W` within-block | fixed-margin, **−0.070 to −0.287** (`#80`) | **severely under-reported** |

| # | Claim | Verdict |
|---|---|---|
| 81a | **`#80c`'s account of `#71`** | **WRONG in both particulars.** Neither side was corrected, and the uncorrected comparison **handicaps the domain-SPECIFIC side**, not the general one |
| 81b | **`#71`'s conclusion, "domain-general beats domain-specific, 186/276"** | **NOW AT RISK, in the direction `#80c` claimed was safe.** At Kw=4 Brier, `W` raw is −0.065 and `W` corrected is **+0.177**; `C` at the same specification is ≈ +0.007 raw and ≈ +0.007 corrected. If those hold in one symmetric run, `#71` **inverts** |
| 81c | **why `#80c` happened** | I described a round I wrote nine entries earlier without opening it. The correction cost one `grep`. **Door ① in miniature: my evidence was a story about the object, and the object was 30 lines away** |

`#71` is the load-bearing claim of arc `A10` — the README's "domain-general" row cites it — so it is
**flagged `AT RISK` and not yet withdrawn**. The symmetric run is the next round, and it is designed
to be able to overturn `#71` rather than to defend it.

---

## Entry 82, added by `E01·A10·R125` — #71 inverts: the person-side readout is mostly assembled PER DOMAIN

`#81` flagged `#71` `AT RISK` because it compared a nearly-unbiased `C` against a severely
under-reported `W`. Symmetric run, both sides referred to their own nulls in the same run, same
`Kc × Kw` grid:

| Kc | Kw | `C` raw | `C` null | **`C` corrected** | `W` raw | `W` null | **`W` corrected** | gap raw | **gap corrected** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 0.0028 | −0.0001 | **0.0030** | 0.0168 | −0.0681 | **0.0849** | −0.014 | **−0.082** |
| 4 | 2 | 0.0061 | −0.0007 | **0.0068** | −0.0101 | −0.1347 | **0.1246** | +0.016 | **−0.118** |
| 8 | 4 | 0.0103 | −0.0014 | **0.0117** | −0.0666 | −0.2429 | **0.1763** | +0.077 | **−0.165** |

**Both gates pass.** Each null moves its own quantity (|ΔC| = 0.0073, |ΔW| = 0.1287); margins
asserted exact on every curveball draw. And the **doubly-destroyed world** — fixed-margin randomised
*and* person-permuted — returns **corrected C = −0.0003 and corrected W = +0.0001**, so neither
correction manufactures signal by subtracting a too-negative floor.

| | `#71` as published | `#82` symmetric |
|---|---|---|
| domain-**general** larger | **186** (67.4%) | **0** (0.0%) |
| domain-**specific** larger | 32 (11.6%) | **201** (97.1%) |
| not distinguishable | 58 | 6 |
| median gap | +0.016 | **−0.118** |

| # | Claim | Verdict |
|---|---|---|
| 82a | **`#71`, "domain-general beats domain-specific, 186/276"** | **WITHDRAWN, AND INVERTED.** Corrected, the within-block structure is **0.064–0.163** against cross-block **0.002–0.012** — **7× to 26× larger**, in **201 of 207** comparisons, none the other way. The 186/276 margin was entirely the **−0.135 handicap** `#71` handed the specific side |
| 82b | **`#71`'s parameter argument, which is not what was wrong with it** | **STILL TRUE, and now clearly a different claim.** `C` spends 18–144 target parameters to `W`'s 3,249–12,996. So the cross-block structure is far more **parameter-efficient** while being far **smaller**. `#71` asserted the second on evidence for the first — *efficient* and *larger* are not the same word |
| 82c | **the A-vs-B reading this changes** | **The person-side readout is mostly assembled PER DOMAIN.** `#70` said only the person side crosses blocks — still true, the item side is 0 there by construction. But the transferable part is **a small fraction of the person-side structure**, not its bulk. "A general eroticization operator" is substantially weaker than arc `A10` had been reporting |
| 82d | **what is NOT affected** | `#72`'s no-cliff result (about `C` alone, and `C` is corrected there) · `#69`'s derivation finding · `#79b`/`#80a`'s scale-freeness of the item-vs-interaction ordering · `#70`'s two-question split. Only the **general-vs-specific ranking** falls |

**This retraction exists because `#81` corrected `#80c`.** `#80c` asserted `#71` was conservative in
the direction it concluded; one `grep` showed the opposite; the round that followed overturned `#71`
outright. **A memory-sourced sentence about my own prior work was the only thing standing between
this project and a false headline it would have kept.**

---

## Entry 83, added by `E01·A10·R126` — the domain-specific spectrum has a knee at 5; the domain-general one has none

`#72` swept the cross-block rank and found no cliff. `#82` then showed the cross-block part is the
**thin** one, so `#72` characterised a residual. Same question, asked of the structure that carries
the mass, with the same known-rank discipline — **synthetic worlds whose within-block rank is 2 and 5**.

| Kw | **real, corrected** | ctrl r=2 | ctrl r=5 | seed spread |
|---:|---:|---:|---:|---:|
| 1 | 0.0917 | 0.2193 | 0.1833 | 0.0057 |
| 2 | 0.1327 | **0.3331** | 0.3192 | 0.0071 |
| 3 | 0.1658 | 0.3371 | 0.4314 | 0.0082 |
| **5** | **0.1995 ← peak** | 0.3308 | **0.5402 ← peak** | 0.0100 |
| 8 | 0.1444 | 0.2178 | 0.4347 | 0.0110 |
| 16 | 0.0294 | 0.0730 | 0.1297 | 0.0041 |

**Gain per dimension**, which is where the knee lives:

| world | knee | increments |
|---|---:|---|
| ctrl **r=2** | **2** | 0.1139, 0.0040, −0.0032, −0.0377, … |
| ctrl **r=5** | **5** | 0.1359, 0.1122, 0.0544, −0.0352, … |
| **real** | **5** | 0.0410, 0.0332, 0.0168, −0.0184, … |

| # | Claim | Verdict |
|---|---|---|
| 83a | **the estimator counts within-block rank** | **CONFIRMED by dose-response.** Knee at 2 for a rank-2 world, 5 for a rank-5 world, both with the gain going **negative** immediately after |
| 83b | **the domain-specific spectrum** | **KNEE AT 5, and the shape is indistinguishable from a true-rank-5 world.** Corrected skill peaks at +0.1995 and the per-dimension gain **changes sign** between rank 3 and 5 — exactly where the r=5 control's does |
| 83c | **the structural statement, which is the sharpest this project has made** | **The person-side readout is ~5 domain-specific factors per block, plus a general part with no rank limit** (`#72`: cross-block gain still positive and declining smoothly at rank 32, never changing sign). Two levels, two different kinds of object |
| 83d | **my cliff-ratio statistic** | **MIS-SPECIFIED — the fifteenth.** It divides consecutive per-dimension gains, and the denominator **changes sign** at the knee, so it returned 1.7 × 10¹⁰ and printed a "cliff at rank 3" that does not exist. Same family as `#79f` — an operator applied across a sign boundary. The **knee** is the correct statistic and it is calibrated; the ratio is discarded |
| 83e | **the top of the sweep** | **DEGENERATE, and I should have excluded it by design.** Blocks have m = 10–24 options, so rank 24 *is* full rank: the approximation reconstructs the observed residual exactly and real and null both return **0.0000**. Ranks above ~16 carry no information here |

**The confound `83c` must carry, and it is not small.** The within-block sweep lives in a space of at
most **m ≈ 10–24** dimensions; the cross-block sweep of `#72` lived in ~**500**. So *"specific is
low-rank, general is not"* is partly a statement about **how many dimensions were available to
each**. A knee at 5 out of ≤24 and no knee at 32 out of ~500 are not like-for-like. Settling it
needs the cross-block sweep re-run at the *same* dimensional budget — that is the residual gap, and
it is named here rather than left for a reader to notice.

---

## Entry 84, added by `E01·A10·R127` — the confound `#83c` named cannot be settled here, and the control is what says so

`#83c` claimed the within-block spectrum knees at 5 while the cross-block one does not, and named its
own confound: the two estimators search spaces of very different size (m = 10–24 columns vs ~500).
This round matched the budget — cross-block scores built from a **random subset of exactly m**
other-block columns, three independent subsets per block.

**The gate fired before the answer.** A synthetic world with a **known shared rank of 5**, run at the
same restricted budget:

| | per-dimension gains |
|---|---|
| known rank-5, **wide** basis (`#83`, within-block analogue) | 0.1359, 0.1122, 0.0544, **−0.0352** ← sign change at 5 |
| known rank-5, **budget-matched** basis (here) | 0.0288, 0.0234, 0.0126, 0.0063, 0.0038, 0.0001 — **no sign change at all** |

| # | Claim | Verdict |
|---|---|---|
| 84a | **the budget-matched comparison** | **UNIDENTIFIED, and that is not an acquittal.** At m columns the estimator cannot count rank even when the rank is known and equal to 5. An n × 18 SVD cannot separate 5 true factors from noise, so each added dimension recovers a fraction of the truth instead of exhausting it. **The instrument needs the wide basis to count, and with the wide basis the budgets are unmatched** |
| 84b | **`#83c`, "the two levels are different kinds of object"** | **DOWNGRADED from "the sharpest structural statement this project has made" to a contrast whose confound is now known to be UNRESOLVABLE ON THIS RELEASE.** It is not refuted; it is unsupported in a way I cannot fix here |
| 84c | **what the matched run does show, stated at the size it supports** | The real budget-matched curve rises monotonically (+0.0006 → +0.0097 over ranks 1→16) with **no sign change**, and so does the known-rank-5 control (+0.0445 → +0.1565). **At equal budget the real cross-block spectrum is shape-indistinguishable from a true-rank-5 world** — 16× smaller in magnitude, identical in the only feature the knee reads. That is weak evidence *against* `#83c`, and it is weak precisely because the instrument is blind here |
| 84d | **subset stability, published rather than averaged** | Knees per subset: **16, 16, 8**. Subset 2 disagrees, and its disagreement comes from a single noisy increment (−0.00006 at rank 12). One draw in three flips a knee, which is itself a reason not to read a knee off this arm |

**What would settle it** — the register entry, since `#83c` cannot be closed here:
- a release with **more options per block**, so the within-block estimator has a comparable budget; or
- **many more blocks**, so a matched subset is still large enough to resolve rank; or
- a within-block estimator that borrows strength across blocks without borrowing structure — which is
  the same identification problem in a different place.

**The round's real product is the control.** Without the known-rank-5 arm this would have read as
*"budget-matched general spectrum has no knee, so `#83c` stands"* — a confirmation, published, from
an instrument that had just been shown blind. That is `P5`'s star rule: **a `not found` is
inadmissible until the same instrument has passed a positive control**, and here it did not.

---

## Entry 85, added by `E01·A10·R128` — the three components are the same size, and "11.7×" was three raw numbers against one corrected one

`#82` inverted `#71` by referring both sides to their own nulls. The audit that followed reached the
README's **first** row, which still read *"item +0.222 vs person +0.085 vs interaction +0.019"* —
raw against raw, the exact comparison `#82` had just killed one row below.

Every component given the null that destroys **it** and preserves everything else, in one run:

| component | its null | invariant asserted | real → null |
|---|---|---|---|
| `I` item | within-person shuffle | row sums exact | **+0.2041 → −0.0004** |
| `P` person | within-column shuffle | column sums exact | **+0.0928 → −0.0525** |
| `W` within-block | fixed-margin curveball | **both** margins exact | −0.0536 → −0.1938 |
| `C` cross-block | person-permutation | — | +0.0063 → −0.0006 |

| K | `I` raw | **`I_c`** | `P` raw | **`P_c`** | `W` raw | **`W_c`** | **`C_c`** | interaction raw | **interaction corrected** | ratio raw | **ratio corrected** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.2081 | **0.2084** | 0.0864 | 0.1482 | 0.0168 | 0.0853 | 0.0030 | 0.0196 | **0.0883** | 10.6× | **2.36×** |
| 2 | 0.2054 | 0.2058 | 0.0914 | 0.1479 | −0.0087 | 0.1253 | 0.0040 | −0.0049 | 0.1293 | −41.7× | **1.59×** |
| 3 | 0.2034 | 0.2038 | 0.0936 | 0.1454 | −0.0325 | 0.1600 | 0.0056 | −0.0274 | 0.1656 | −7.4× | **1.23×** |
| **5** | 0.2003 | 0.2006 | 0.0969 | 0.1432 | −0.0924 | **0.1906** | 0.0084 | −0.0848 | **0.1989** | −2.4× | **1.01×** |
| 8 | 0.2034 | 0.2037 | 0.0958 | 0.1419 | −0.1513 | 0.1396 | 0.0135 | −0.1392 | 0.1531 | −1.5× | **1.33×** |

| # | Claim | Verdict |
|---|---|---|
| 85a | **"the item effect is 11.7× / 3.5× the interaction" (`#68`, README row ①)** | **WITHDRAWN.** Corrected, the ratio is **1.0×–2.4× across the whole rank sweep** and **1.05× at each side's own best rank**. The 11.7× was a corrected item effect against an uncorrected interaction carrying a −0.19 estimator bias |
| 85b | **the three components, symmetrically** | **THE SAME SIZE.** `I_c` **+0.208** · interaction `+0.199` · `P_c` **+0.148**. The person effect was also under-reported raw (0.093 → 0.148); only the item effect needed no correction, which is exactly what a main effect a margin-preserving null cannot touch should do |
| 85c | **`#68c`, "the best case for the epoch title is *both, equally*"** | **PROMOTED from adversarial edge case to CENTRAL ESTIMATE.** It was reached there by choosing the rank most generous to the rival; it is reached here by giving every component its own null and its own best rank — a symmetric procedure, not a generous one |
| 85d | **the per-block split, published because the aggregate hides it** | **item larger in 12/23, interaction larger in 7/23, tied in 4/23.** Median gap **+0.049** against a 2× spread of 0.020 — so the *typical* block does favour the item, while the *mean* is a tie, and **7 blocks go the other way**. "Tied" is the aggregate; the distribution is a genuine split |
| 85e | **what the epoch question now answers to** | Neither *"a value, not a category"* nor its converse. **Content and individualised valuation contribute equally**, with person breadth close behind. `A09`'s reversal of the epoch title stands as *"not larger"*; it never licensed *"much smaller"*, and for eighteen entries this project reported that it did |

**The N/A that makes this conservative toward the item.** No null destroys `I` while preserving the
interaction — a within-person shuffle kills both. So `I_c` is an **upper** estimate of the item
effect's uniquely-attributable part, and the true ratio is **at most** what is reported here.

---

## Entry 86, added by `E01·A05·R097` — a null correction credits a component with the DAMAGE its estimator does, and 87% of the person effect's correction was that

`#85b` reported the person main effect at **+0.148** corrected against **+0.093** raw, and explained
the gap as estimator noise: item effects are estimated from 1,200–15,000 observations per column,
person effects from 10–24 per row. **That explanation is testable — if it is noise, a shrunken
estimator should deliver the gap directly, with no correction at all.**

| world | raw row mean | James-Stein | empirical Bayes | shrink factor |
|---|---:|---:|---:|---:|
| **planted** person effect | **−0.0072** | **+0.0234** | +0.0210 | 0.536 |
| **no** person effect (within-column shuffle) | **−0.0642** | −0.0008 | **−0.0018** | 0.091 |
| **real** | +0.0866 | **+0.0959** | +0.0945 | 0.743 |

**Both gates pass, and the planted arm is stark: the raw row mean scores −0.0072 on a world
containing a genuine person effect** — worse than not using it — while shrinkage finds it. And on a
world with none, shrinkage returns −0.002 where the raw estimator returns −0.064.

| # | Claim | Verdict |
|---|---|---|
| 86a | **`#85b`'s "person +0.148"** | **OVERSTATED.** Shrinkage recovers only **12.7%** of the +0.0616 gap. The best point estimate is the **shrunken** one, **+0.0959**, and the honest report is the interval **[0.087, 0.148]** weighted toward its low end |
| 86b | **what a null correction actually measures** | **Skill PLUS the damage a misfit estimator does.** In the no-person world the raw estimator scores −0.064 — that is not a floor the component clears, it is a hole the *estimator* digs. Subtracting it credits the component with 0.064 it never earned. **87% of `#85b`'s person correction was that** |
| 86c | **the alarm, and it is not small** | **The same critique applies to `W` and `C`.** `W`'s fixed-margin floor is **−0.19** — an estimator digging the same hole — and `#82`, `#83` and `#85` all rest on `W_c = W_real − W_null`. If a *shrunken* low-rank estimator recovers only ~13% of that gap too, then `#82`'s inversion, `#83`'s knee and `#85`'s "three components are equal" all move, and the interaction's magnitude is far below **+0.199** |
| 86d | **`#85`'s comparison in the meantime** | **FLAGGED, not withdrawn.** Item +0.208 needs no correction (its null is −0.0004, so its estimator digs no hole). The interaction's +0.199 is the number now in doubt. **The direction of the doubt favours the item**, which is the opposite of the direction the last four entries moved |

**Why this was findable.** `#85` explained its own correction with a mechanism — *"person effects are
estimated from 18 observations"* — and a mechanism is a prediction. The prediction was that
shrinkage would close the gap. It closed 13% of it. **An entry that explains itself can be attacked
through its explanation; one that only reports a number cannot.**

**N/A:** separating breadth from acquiescence needs a reverse-keyed or forced-choice item. This
release has none, so nothing here distinguishes *"wants more"* from *"ticks more"*.

---

## Entry 87, added by `E01·A10·R129` — the −0.19 floor is an estimator artifact, and my own comparison let one arm peek at the test set

`#86c` predicted that if the interaction's null correction is crediting **damage** rather than
**noise**, a properly regularised estimator would recover only ~13% of the gap directly — and that
`#82`, `#83` and `#85` would then all move.

| world | best-rank **hard** truncation | tuned **soft** thresholding | λ |
|---|---:|---:|---:|
| fixed-margin (no interaction) | **−0.0615** | **+0.0032** | 3.45 |
| planted rank-5 | +0.2596 | +0.1884 | 3.97 |
| **real** | +0.0258 | **+0.0502** | 3.57 |

| # | Claim | Verdict |
|---|---|---|
| 87a | **my own gate (a)** | **FAILED, and the flaw is mine — the sixteenth mis-specified design element.** The hard arm takes `max over K of held-out skill`, i.e. **selects its rank on the test cells**, while the soft arm tunes λ on a validation split carved out of training cells. One arm peeked. The planted-world comparison is therefore void |
| 87b | **the fixed-margin floor, which survives the flaw** | **THE −0.06 HOLE IS AN ESTIMATOR ARTIFACT.** Hard truncation scores −0.0615 on a world with **no interaction at all**; soft thresholding scores **+0.0032**. Oracle rank-selection can only *help* hard, so its hole is real and not a selection artifact. **A regularised estimator has a floor of ≈ 0** |
| 87c | **the real-world comparison, which also survives** | **Soft beats oracle-selected hard: +0.0502 vs +0.0258.** A fortiori, since the loser had the unfair advantage |
| 87d | **`#86c`'s alarm** | **CONFIRMED in direction and magnitude.** Recovery = (0.0502 − 0.0258)/(0.1906 − 0.0258) = **14.8%**, against the person effect's 12.7%. The correction in `#82`/`#83`/`#85` was **crediting the interaction with the hole hard truncation digs** |
| 87e | **the interaction's honest magnitude** | **≈ +0.050**, not +0.199 — a **4× reduction** — measured with an estimator whose own null is ≈ 0 and therefore needs no correction at all. Against the item effect's **+0.208**, that is **≈ 4×**, not the 1.05× tie `#85` reported |

**`#82`, `#83`, `#85` are all FLAGGED `AT RISK`, not withdrawn.** `87a` means this round cannot
carry a retraction on its own: the honest comparison needs the hard arm to select K on training
cells like the soft arm does. That is one line, and it is the next round. But the direction is
already visible in `87b`, which does not depend on the flaw.

**The shape of the error, now seen twice in two entries.** `#86` and `#87` are the same finding at
two components: **subtracting a null does not correct an estimator, it credits the estimator's
failure to the thing being estimated.** The fix is not a better null — it is a better estimator,
whose null is zero because it does not overfit in the first place.

---

## Entry 88, added by `E01·A10·R130` — with nothing subtracted anywhere, the item effect is 3.5–5.6× the interaction, and `#85`'s tie was the correction machinery

`#86` and `#87` established that subtracting a null credits the estimator's failure to the component.
This round removes the machinery entirely: **shrunk column means · empirical-Bayes shrunk row means ·
per-column ridge on external scores · soft singular-value thresholding**, every hyperparameter chosen
on a validation split carved out of the **training** cells.

**Negative controls — each regularised component in the world that destroys it:**

| | value |
|---|---:|
| `I` under within-person shuffle | **−0.0000** |
| `P` under within-column shuffle | **−0.0002** |
| `C` under person-permutation | **−0.0021** |
| `W` under fixed-margin curveball | **−0.0219** ← the one that misses |

**`#87a` SETTLED — soft vs hard, both tuned honestly on training cells:**

| world | soft | hard |
|---|---:|---:|
| real | **+0.0273** | −0.0091 |
| no-interaction | **−0.0219** | **−0.0945** |

**The decomposition, no correction applied anywhere:**

| `I` | `P` | `C` | `W` | interaction (C+W) | **item : interaction** |
|---:|---:|---:|---:|---:|---:|
| **+0.2173** | +0.1034 | +0.0092 | +0.0306 | +0.0388 | **5.59×** |

Per block: **item larger in 21/23**, interaction larger in 1, tied 1. Median gap +0.178 against a 2×
spread of 0.014. Chosen hyperparameters: λ = 4.87, ridge α = 90.4.

| # | Claim | Verdict |
|---|---|---|
| 88a | **gate (b), "planted world recovers I, P and W"** | **FAILED ON MY PLANT, NOT THE ESTIMATOR.** I planted a person effect of sd 0.08. The graded ladder shows recovered `P` = −0.0002 / +0.0106 / +0.0571 / +0.2229 at sd 0/0.08/0.15/0.30 — monotone, silent at zero, and **my plant sat below its own MDE**. The real `P` of +0.1034 corresponds to a planted sd of ≈ 0.19. **Gate (b) PASSES when run above the magnitude it can detect** |
| 88b | **gate (a), `W` in the no-interaction world** | **FAILS by 0.0019** against a threshold I chose. A residual hole of −0.022 remains — 4.3× smaller than hard truncation's −0.0945, but not zero. **Its direction UNDER-reports the interaction**, so correcting it gives `W` ≈ +0.053, interaction ≈ +0.062, ratio **3.5×**. The conclusion survives the flaw in the direction that matters |
| 88c | **`#85`'s "the three components are the same size, 1.05×"** | **WITHDRAWN.** It was the correction machinery. On estimators that need no correction the ratio is **3.5×–5.6×** and the item is larger in **21 of 23 blocks**, not 12 |
| 88d | **`#82`'s inversion (`W` > `C`)** | **DIRECTION SURVIVES, MAGNITUDE DOES NOT.** Regularised: `W` +0.031 vs `C` +0.009 — **3.3×**, not the 7–26× `#82` reported. The domain-specific part is still the larger of the two person-side components, and both are small next to the item effect |
| 88e | **`A09`'s original direction** | **REINSTATED, on evidence it did not have at the time.** `#67`/`#68` said the item effect is the larger component and got there through an uncorrected comparison that happened to point the right way; `#85` overturned it with a correction that credited estimator damage; this round reaches it with **nothing subtracted at all** |

**Three reversals on one question in eleven entries** — `#68` item-dominant → `#85` tied → `#88`
item-dominant. Every step was driven by a methodological fix, not new data, and the **first and last
agree while the reasoning in between was wrong in both directions.** That is worth more than the
number: the project's error was never the direction, it was **comparing quantities whose estimators
fail differently**, and it took `#80`, `#82`, `#84`, `#86` and `#87` to name it.

---

## Entry 89, added by `E01·A10·R131` — the estimator with the honest floor cannot count dimensions, and the one that can count has a dishonest floor

`#88` showed hard rank truncation digs a −0.09 hole on a structureless world while soft thresholding
digs −0.02. The natural follow-up was to re-ask `#72`/`#83`'s dimensionality question on the better
estimator, where rank is not swept but **emerges** as the number of singular values surviving the
training-tuned threshold.

| world | tuned λ | **effective rank** | rank sd | skill |
|---|---:|---:|---:|---:|
| **no structure** (fixed-margin) | 3.32 | **16.10** | 4.40 | −0.011 |
| **known rank 2** | 4.04 | **15.80** | 4.59 | +0.127 |
| **known rank 5** | 4.30 | **15.78** | 4.57 | +0.166 |
| **real** | 3.33 | **16.00** | 4.51 | +0.035 |

| # | Claim | Verdict |
|---|---|---|
| 89a | **effective rank as a dimensionality estimator here** | **DEAD.** It reads ~16 on a rank-2 world, ~16 on a rank-5 world, and **~16 on a world with no structure at all**. Blocks have m = 10–24 columns, so this is near-full-rank everywhere: the count above λ is set by the **noise spectrum of the completed matrix**, not by the signal. Both gates fail, and they fail in the way that says the statistic is uninformative rather than the data is |
| 89b | **`#72` and `#83`, measured with hard truncation** | **NOT overturned by this.** Their knee was **calibrated** — a rank-2 world kneed at 2, a rank-5 world at 5, real at 5 — and `#88` showed hard truncation's **magnitude** is biased, not its **rank ordering**. A biased scale can still locate where added dimensions stop paying |
| 89c | **the resulting standoff, which is the finding** | **The two estimators are complementary and neither is sufficient.** Hard truncation can locate a knee and cannot be trusted on magnitude; soft thresholding is trustworthy on magnitude and carries no rank information. **Every dimensionality claim in this project rests on the estimator that `#88` disqualified for size**, and no estimator available here does both |
| 89d | **the skill column, which is a separate reading** | Real within-block skill at the tuned λ is **+0.035**, against **+0.127** for a planted rank-2 world and **+0.166** for rank-5 (both at loading scale 0.30). Whatever within-block structure exists is **far weaker than either control** — consistent with `#88`'s +0.031 and inconsistent with the +0.19 the correction machinery reported |

**What would settle dimensionality here** — the register entry, since `89c` closes no door quietly:
an estimator with a **selection-consistent** penalty (a rank-aware information criterion tuned by
cross-validation over *both* rank and shrinkage), or blocks with enough options that the noise
spectrum is not near-full-rank. This release has 10–24 options per block, and that is the binding
constraint on every rank question the project has asked.

---

## Entry 90, added by `E01·A10·R132` — on the probability scale the three components are the same size, and the item effect wins every skill measure because it is estimated 179× better

Every number in this project has been held-out skill — a quantity whose meaning depends on the
estimator, which is what `#86`–`#89` spent four entries discovering. This round inverts the plant:
**what per-cell probability perturbation, passed through the same estimator, reproduces the real
data's skill?**

| | percentage points |
|---|---:|
| option base rate — observed 22.6, binomial noise 0.7 | **± 22.6 pp** |
| person overall rate — observed 19.8, binomial noise **10.4** | **± 16.3 pp** |
| person × option interaction — inverted plant, family CV **16%** | **± 23.7 pp** |

**And the reconciliation, which is the entry.** Observations behind each estimate: **item n = 3,228 ·
person m = 18 — a 179× difference.** So:

| | magnitude (pp) | held-out skill (`#88`) |
|---|---:|---:|
| item | ±22.6 | **+0.217** |
| person | ±16.3 | +0.103 |
| interaction | **±23.7** | **+0.039** |

| # | Claim | Verdict |
|---|---|---|
| 90a | **the epoch question, answered on a scale that does not depend on an estimator** | **Content and individualised valuation move endorsement probability by the SAME AMOUNT** — ±22.6 pp vs ±23.7 pp. `#88`'s 5.6× skill advantage for the item effect is **entirely a learnability difference**, not a size difference: a base rate is estimated from 3,228 observations, a person's deviation from 18 |
| 90b | **what this makes of eleven entries of reversals** | **They were measuring predictability and calling it composition.** `#68` (item-dominant), `#85` (tied), `#88` (item-dominant) were all correct *about skill* and all silent about magnitude. **The interaction is the largest source of variation in this data and the smallest source of predictable variation**, and both halves of that sentence are measured |
| 90c | **my magnitude accounting, first pass** | **WRONG — the seventeenth mis-specified statistic.** I reported the planted perturbation's sd **before** clipping to [0.02, 0.98], overstating it by 5–56% depending on scale (at rank 10, scale 0.30: 94.6 pp pre-clip vs 41.7 pp realized). Corrected, the implied interaction falls from ±30.8 to **±23.7 pp** — and the family's CV **improves** from 22.4% to 15.9%, so the corrected quantity is also the better-identified one |
| 90d | **the person effect's noise share, which nobody had priced** | **18% of the observed person spread is binomial noise** (10.4 pp of 19.8). That is the same fact `#86` found as a 12.7% shrinkage recovery, now in interpretable units, and it is why the person component behaves differently from the item one in every round of this arc |

**Scope, and it is a shape assumption.** The inversion assumes the real interaction has the same
**Gaussian low-rank** shape as the plant. A sparser structure — a few people with strong specific
tastes rather than everyone deviating a little — would produce the same skill at a different sd.
Item and person spreads are **measured directly**; the interaction's is **inferred through a model**,
and that asymmetry is not removable here.

---

## Entry 91, added by `E01·A10·R133` — a structure carried by 5% of people is invisible at ANY strength, and that is what a fetish looks like

`#90` named its own shape assumption: the ±23.7 pp inversion assumed a **dense** structure, and a
sparse one — a few people with strong specific tastes — might imply a different magnitude. This
round planted sparse structures to find out. **It could not, and the reason is the finding.**

| carriers | scale 0.08 | 0.12 | 0.20 | 0.30 | **0.50** | per-carrier sd at 0.50 |
|---|---:|---:|---:|---:|---:|---:|
| **5%** | −0.013 | −0.013 | −0.012 | −0.012 | **−0.012** | **±50.4 pp** |
| **15%** | −0.010 | −0.007 | −0.003 | −0.001 | +0.002 | ±48.8 pp |
| **30%** | −0.007 | +0.001 | +0.014 | +0.024 | +0.033 | ±47.1 pp |
| **100%** | −0.002 | +0.024 | +0.089 | +0.158 | **+0.239** | ±42.5 pp |

The estimator's floor is **−0.012**. **At 5% carriers the skill never leaves the floor**, even when
those carriers' endorsement probabilities are moved by **±50 percentage points**.

| # | Claim | Verdict |
|---|---|---|
| 91a | **the shape test as designed** | **UNRUNNABLE — gate (b) fails.** Only the dense arm brackets the real skill of +0.036; no sparse arm reaches it at any magnitude in the ladder. The inversion cannot be compared across shapes because three of the four shapes cannot produce the observation at all |
| 91b | **the capability boundary, which is what the failure measures** | **This estimator — and every estimator used in this project — is BLIND to interaction structure carried by a minority.** 5% carriers: undetectable at ±50 pp. 15%: barely (+0.002 against a −0.012 floor). Detection needs **≈30% of people or more** |
| 91c | **what that does to `#90a`** | **STRENGTHENS the dense reading and bounds it at the same time.** Since sparse shapes cannot reproduce +0.036, the observed interaction **must** have a dense component — the shape assumption was not arbitrary. But a sparse structure could coexist **entirely invisibly**, and nothing in this project would show it |
| 91d | **the consequence for the project's own subject matter, and it is not a technicality** | **A minority of people with intense, specific attachments is precisely what a fetish is, phenomenologically.** Every method in these eleven arcs — factor analysis, CCA, low-rank completion, held-out skill — is a **variance-explained** method, and variance-explained methods weight by prevalence. **The eroticization operator, if it works the way kink actually presents, is exactly the object this entire toolkit cannot see** |
| 91e | **the dense arm as a replication check** | **PASSES.** It returns **±23.5 pp** against `#90`'s ±23.7 pp, on a separately-written pipeline. `#90`'s number replicates |

**This is a capability-boundary update, not a claim update.** Nothing measured so far becomes wrong.
What changes is the **scope of every negative result in the project**: 105 rounds of "no structure
found" now carry a quantified blind spot — *no structure found among the ≥30% of people who share
it*. `#62`, `#55`, `#49`, `#39` and every other null in this ledger inherit it.

**What would see a sparse structure**: a method that scores **per person** rather than per cell —
individual-level anomaly detection, or a mixture model with a small high-deviation component. Both
are buildable on this release, and neither has been tried.

---

## Entry 92, added by `E01·A11·R135` — the only quantile where the control fires is the only one where the data doesn't

`#91` measured that every method in this project is blind to a structure carried by a minority. This
round builds the missing one: a **per-person goodness-of-fit** statistic
`T_i = mean of (M−p̂)²/p̂(1−p̂)` over that person's held-out cells, pooled across all their blocks. No
per-person parameter is fitted, so nothing can overfit, and every person counts once regardless of
prevalence.

**Positive control — the exact world `#91` showed was invisible (5% of people, ±50 pp):**

| quantile | sparse5 | its null | diff | 2× spread | |
|---|---:|---:|---:|---:|---|
| share T>2 | 0.0080 | 0.0013 | +0.0067 | 0.0084 | **blind** |
| p90 | 1.3543 | 1.3182 | +0.0361 | 0.1372 | **blind** |
| p95 | 1.4710 | 1.4518 | +0.0192 | 0.2147 | **blind** |
| **p99** | 1.9462 | 1.7161 | **+0.2301** | 0.1781 | **SEES IT** |

**Real data vs the parametric null:**

| quantile | real | null | diff | 2× spread | |
|---|---:|---:|---:|---:|---|
| **p50** | 0.8353 | 0.8916 | **−0.0563** | 0.0538 | **RESOLVABLE — real is UNDER-dispersed** |
| p75 | 1.0856 | 1.0785 | +0.0071 | 0.1577 | no |
| p90 | 1.7196 | 1.4398 | +0.2798 | 0.3282 | no |
| **p95** | 2.0463 | 1.8120 | **+0.2343** | 0.1561 | RESOLVABLE |
| **p99** | 2.8271 | 2.2313 | +0.5958 | 0.9558 | **no** |
| **share T>2** | **0.0592** | **0.0287** | +0.0305 | 0.0217 | RESOLVABLE — **2.06×** |

| # | Claim | Verdict |
|---|---|---|
| 92a | **the verdict this round printed** | **WRONG — the eighteenth mis-specified verdict.** It keyed on **p99 alone**, the noisiest quantile, and printed *"no minority structure"* while p95 and the share-above-2 both resolved in the other direction |
| 92b | **…and the correct verdict is still `UNVERIFIED`, for a better reason** | **The only quantile where the positive control FIRES (p99) is the only one where the real data does NOT resolve.** At p95 and share-T>2 the data resolves but the instrument is **blind there** by its own control. `P5`'s star rule, in its sharpest form yet: a measurement is inadmissible where the instrument has not been shown to work, whichever way it comes out |
| 92c | **the null is mis-specified, independently** | **Real is UNDER-dispersed at the median (−0.056, resolvable).** Respondents pick a roughly stable *number* of options per block, which makes within-row picks negatively correlated — the independent-Bernoulli null cannot reproduce that. **A null that misses the median cannot be trusted in the tail.** The fix is a row-sum-preserving (curveball) null, per person |
| 92d | **the method itself** | **A real improvement, and a marginal one.** The old methods were blind to 5% carriers at every magnitude (`#91`); this one sees them at p99 only. That is the difference between *impossible* and *barely* |
| 92e | **`world='null'` parsed as NaN** | **The nineteenth gotcha, and a new family.** `null` is in pandas' default NA list, so the entire null arm vanished from the first analysis and every comparison returned `nan`. `#80d` banned column NAMES colliding with methods; this is a column **VALUE** colliding with a parser default. Same rule extended: **`null`, `NA`, `nan`, `None`, `N/A`, `inf` are never arm labels** |

**What this round is worth despite resolving nothing.** It converted `#91`'s blind spot from *"we
cannot see minorities"* to *"we can see them at p99, we have three seeds, and the null is wrong at
the median."* All three are fixable — more seeds, a row-sum-preserving null — and none of them was
knowable before the method existed.

---

## Entry 93, added by `E01·A11·R136` — with the right null the apparent tail vanishes, and so does the instrument's ability to see anything

`#92` left one admissible quantile and two fixable faults. Both are fixed here: the null becomes
**fixed-margin (curveball)**, reproducing the response format instead of mistaking it for signal, and
the seeds go from 3 to 8.

**Real vs fixed-margin null, 8 seeds — every quantile:**

| | real | null | diff | 2× spread | |
|---|---:|---:|---:|---:|---|
| share T>2 | 0.0590 | **0.0495** | +0.0095 | 0.0406 | no |
| p50 | 0.8373 | 0.8471 | −0.0099 | 0.0667 | no — **`#92c`'s fault is fixed** |
| p75 | 1.1338 | 1.0775 | +0.0562 | 0.1579 | no |
| p90 | 1.7796 | 1.7381 | +0.0415 | 0.2876 | no |
| p95 | 2.0558 | 2.0130 | +0.0428 | 0.3093 | no |
| **p99** | 2.6429 | 2.6350 | **+0.0079** | 0.8566 | no |

**Positive control — 5% carriers against its own fixed-margin null: blind at every quantile**, and
at p99 the null is *higher* than the plant (2.1673 vs 1.9208).

| # | Claim | Verdict |
|---|---|---|
| 93a | **`#92`'s apparent tail excess (p95 +0.234, share T>2 at 2.06×)** | **WAS THE NULL.** Against a null that reproduces the response format, the share above T=2 goes from 2.06× to **1.19×** and every quantile difference falls inside its own spread. `#92c` diagnosed the parametric null as mis-specified at the median; it was mis-specified in the tail too, and that is where the apparent finding lived |
| 93b | **the per-person misfit statistic itself** | **UNIDENTIFIED, and this is the real result.** `T_i` and the **row sum are confounded**: a person with strong specific structure picks a different *number* of options, the fitted person effect absorbs it, and a null that preserves row sums removes the signal along with the artifact. **Against the correct null the instrument is blind to the very world it was built to see** |
| 93c | **the verdict** | **UNVERIFIED, and emphatically not an acquittal.** Gate (a) passes, gate (b) fails, and the two together say: *the only null that is admissible is one this statistic cannot work against.* `#92`'s "we can see them at p99" is withdrawn one entry after it was written |
| 93d | **`D[...].T`** | **FOURTH accessor collision — and `T` is on the banned list I wrote in `#80d` two entries ago.** Writing the rule did not prevent the next violation of it. The rule needs to be enforced by something that is not me: a linter, or a naming convention with a prefix (`v_T`, `v_mode`) that cannot collide |

**What would actually work**, and it follows directly from `93b`: a statistic that **conditions on the
pick count**. Given that a person picked *k* options, are those *k* unusually concentrated on **rare**
ones? A surprisal score `S_i = −Σ log(base rate of each picked option)` is exactly that — it is
invariant to *how many* were picked and sensitive to *which*, and fixed-margin randomisation is its
natural and correct null. That is the design `#91` was reaching for and neither `#92` nor `#93` built.

---

## Entry 94, added by `E01·A11·R137`+`R04` — the most promising signal in the project, and a positive control that has failed to license it three times

`#93b` showed the per-person misfit statistic is confounded with the pick count, and named the
statistic that isn't: **mean surprisal** `S_i = mean over picked options of −log(base rate)` —
invariant to *how many* were picked, sensitive to *which*, and with fixed-margin randomisation as its
exactly-matched null.

**Real vs fixed-margin null, 6 seeds, base rates as a fixed external reference:**

| | real | null | diff | 2× spread | ratio |
|---|---:|---:|---:|---:|---:|
| p50 | 0.6988 | 0.7039 | **−0.0051** | 0.0014 | **3.6× — real picks MORE COMMON options** |
| p75 | 0.8076 | 0.7931 | +0.0146 | 0.0017 | 8.6× |
| p90 | 0.9105 | 0.8770 | +0.0336 | 0.0034 | 9.9× |
| **p95** | 0.9792 | 0.9289 | **+0.0503** | 0.0017 | **30×** |
| **p99** | 1.1018 | 1.0313 | **+0.0704** | 0.0049 | **14×** |

**The shape is the fetish signature**: the median person picks options **more common** than a
margin-preserving reassignment would give them, while the upper tail picks **rarer** ones. Bimodal,
in exactly the direction `#91` said the old toolkit could not see.

| # | Claim | Verdict |
|---|---|---|
| 94a | **the signal** | **PROVISIONAL — resolvable at 3.6–30× its own seed spread against an exactly-matched null**, with the pipeline self-consistent (the no-op arm reproduces the real comparison to +0.0499 vs +0.0503). It is the most promising positive result this project has produced |
| 94b | **…and it is NOT licensed, because the positive control has failed three times** | `R03` strength-0 arm **was the real data** (twentieth mis-specified design element) · `R03`'s plant was **self-cancelling** — adding picks to rare options raises their base rates and lowers their surprisal, so a *stronger* plant gave a *smaller* effect (+0.0496 → +0.0407) · `R04`'s swap plant is **below its own MDE**: 5% of people × 2 swaps, against ~35 pooled picks each, moves p95 by −0.0015 |
| 94c | **the pattern across three entries** | **Every one of these failures is a plant below the magnitude the design can detect** — the same fault as `#88`'s gate (b), which I diagnosed there and reproduced twice more here. **I keep choosing plant strengths by intuition instead of deriving them from the statistic's own resolution.** A plant magnitude should be computed, not guessed: it is the effect size that shifts the target quantile by more than 2× its seed spread |
| 94d | **what the failures do NOT do** | They do not weaken `94a`. A control that is too weak to fire says nothing about whether the real effect is real — it says the control was mis-sized. `UNVERIFIED`, and not an acquittal in either direction |

**The fix, derived rather than guessed.** p95 has a 2× spread of 0.0017 and the real effect is
+0.0503, so a control must plant an effect ≥ 0.0017 at p95 to be informative and ideally ≈ 0.05 to be
comparable. With ~35 pooled picks per person, a carrier needs their surprisal to move by ≈ 0.05 × 35
≈ **1.75 nats**, i.e. roughly **two swaps from a median-prevalence option to one ~6× rarer, in every
block they enter** — not two swaps in one block. That is the next round, and its magnitude is now a
calculation instead of a guess.

---

## Entry 95, added by `E01·A11·R139` — LICENSED: a minority whose picks concentrate on rare options, seen by the first statistic in this project whose null preserves exactly what it conditions on

`#94b` withheld the surprisal signal because three positive controls had failed, all by planting
below the detectable magnitude. `#94`'s closing paragraph derived the magnitude instead of guessing
it. Run at that magnitude:

| plant | achieved per-carrier shift | p95 real−null | 2× spread |
|---|---:|---:|---:|
| **0 swaps/block** (a genuine no-op on a separate matrix) | — | **+0.0494** | 0.0020 |
| 1 swap/block | +0.310 (12.9 nats) | **+0.0815** | 0.0065 |
| 2 swaps/block | +0.590 (24.4 nats) | **+0.1808** | 0.0119 |

**Both gates pass.** The 0-swap arm reproduces the real comparison (+0.0494 vs +0.0503). The plant
achieves a shift far above the threshold and **fires**: p95 effect **+0.1314** against a 2× spread of
0.0238.

| # | Claim | Verdict |
|---|---|---|
| 95a | **`#94a`, the surprisal signal** | **LICENSED.** The real data's upper tail picks **rarer options than a margin-preserving reassignment allows**: p95 **+0.0503 at 37× its own seed spread**, p99 +0.0704, p90 +0.0336 — and the median goes the **other way** (−0.0051), so the typical person picks *more common* options than chance. **Bimodal, and in the direction `#91` said was invisible** |
| 95b | **the magnitude, calibrated against the plant** | The real excess is equivalent to **≈1.5 median→rare swaps per block in 5% of people** — or proportionally more people making fewer swaps. Same one-to-many family as `#90`'s (rank, scale); the **product** is what is identified, not the split |
| 95c | **what this adds to eleven arcs of variance-explained methods** | **A component every one of them was structurally blind to.** `#91` measured that blindness (5% carriers undetectable at ±50 pp); this is the same release, the same people, and a statistic that conditions on the pick count instead of absorbing it |
| 95d | **my derived MDE, which was wrong in the safe direction** | I computed that a carrier needed ≈1.75 nats of shift. One swap per block delivers **12.9 nats — 7× more than needed.** The instrument is far more sensitive than my derivation, so the true MDE is well below one swap per block. Deriving the magnitude worked; **my arithmetic for it was off by an order of magnitude**, and only the achieved-shift column revealed that |

**What this is NOT, and the limits are real.** Surprisal cannot be attributed to a cause — an
idiosyncratic erotic attachment, a careless responder, and someone answering a block they barely
relate to all raise `S`. And this is a **population-tail** statement: it says a minority picks rare
options together, not that any particular person does.

**The one-sentence version of the whole project, now that this exists:** content and individualised
valuation move endorsement probability by about the same amount (`#90`), the item effect dominates
every *predictability* measure only because it is estimated 179× better (`#90b`), and **a minority
concentrating on uncommon options exists and was invisible to every method used before this arc**.

---

## Entry 96, added by `E01·A11·R140` — the separator is right, the execution was not mean-matched, and the two normalisations disagree

`#95` named its own limit: surprisal cannot distinguish an idiosyncratic attachment from a careless
responder. But the two make different predictions about **shape across blocks** — careless is **flat**
(high everywhere), attachment is **spiky** (high in a few blocks, ordinary elsewhere) — and that is
testable here. Both plants built, both run against their own fixed-margin nulls.

| world | mean S | between-block sd | **sd / mean** |
|---|---:|---:|---:|
| fixed-margin null, top-5% | 0.9567 | 0.2586 | 0.270 |
| **real, top-5%** | **1.0277** | **0.3004** | **0.292** |
| **flat** plant (careless) | 1.0815 | 0.3084 | **0.285** |
| **spiky** plant (attachment) | 0.7934 | 0.3339 | **0.421** |

| # | Claim | Verdict |
|---|---|---|
| 96a | **gate (a), "the two plants separate at comparable mean"** | **THE LABEL IS FALSE — the twenty-first mis-specified design element.** The plants differ in mean S by **36%** (1.08 vs 0.79). I wrote "at comparable mean" into the gate's text and never checked it. The gate passed on the sd gap alone, which is the quantity the mean confounds |
| 96b | **the verdict the round printed, `SPIKY`** | **NOT LICENSED.** On **raw sd** the real top-5% sits +0.0419 above chance and reads spiky. On **sd/mean** — the natural normalisation, since sd scales with mean — the real value **0.292 sits beside the FLAT plant's 0.285** and nowhere near the spiky plant's 0.421. **Two defensible normalisations, opposite answers** |
| 96c | **what is actually established** | **The design separates**: flat and spiky plants do produce different sd/mean (0.285 vs 0.421, a 48% gap). The instrument works. **The execution did not match means**, so the real data cannot be placed between them |
| 96d | **`#95` itself** | **UNAFFECTED.** `#95` licensed the *existence* of a rare-option minority against an exactly-matched null with a plant that fired monotonically. This round was about *what it is*, and it answers nothing |

**The fix is specific and cheap:** tune the two plants' swap counts until their **mean S matches**
(fewer swaps in the flat arm, more carriers or more blocks in the spiky arm), then compare sd at
equal mean. The 48% sd/mean gap between the plants says there is plenty of signal to read once the
confound is removed.

**Twenty-one mis-specified design elements now**, and this one is a new species: not a wrong
threshold, not a leaky arm, but **a gate whose prose asserted a condition the code never tested.**
The sentence "at comparable mean" was doing the work of a check.

---

## Entry 97, added by `E01·A11·R141`+`R08`+`R09` — the real minority is elevated more EVENLY than any discrete-carrier plant, and that is a different hypothesis than either of the two I was testing

`#96` left the attachment-vs-carelessness question open because the two plants were not mean-matched.
Three rounds to fix it, each caught by its own gate:

| round | what its gate caught |
|---|---|
| `R07` | plants were built by **adding to the real matrix**, so every arm has mean ≥ real by construction and the real point can never lie on a curve. **22nd mis-specified design element** |
| `R08` | rebased on the null — and the **spiky family saturates in mean at 0.976** (S2 0.9659 → S10 0.9759), because there are only 3 rare options per block and a carrier holding all of them cannot get more surprising. **I had parameterised the spiky arm on the axis that saturates. 23rd** |
| `R09` | reparameterised on **breadth**, which does not saturate, turning two ladders that could never meet into one family |

**One family, breadth `b` = blocks elevated per carrier, fixed intensity — all three gates pass:**

| b | mean | sd |
|---:|---:|---:|
| 0 (= the null, exactly) | 0.9556 | 0.2570 |
| 1 | 0.9613 | 0.2764 |
| 2 | 0.9714 | **0.3046** |
| 3 | 0.9839 | 0.3309 |
| 5 | 1.0199 | 0.3904 |
| 8 | 1.0954 | 0.4653 |
| 23 (all) | 1.6490 | 0.4593 |
| **real** | **1.0277** | **0.3004** |

| # | Claim | Verdict |
|---|---|---|
| 97a | **the implied breadth by mean-matching** | **b̂ = 5.31 of 23 blocks (23%)** |
| 97b | **…but the real point is OFF the curve** | At that mean the family predicts sd **0.3981**; the real sd is **0.3004** — short by 0.098, about **10× the family arms' own seed spread**. **The real minority has the MEAN of a b≈5 structure and the SD of a b≈2 one**: its elevation is spread far more evenly across a carrier's blocks than any uniform-intensity discrete-carrier plant produces |
| 97c | **my "lies on the family" tolerance** | **MIS-SPECIFIED — the twenty-fourth.** It compared against the *real arm's* seed spread, which is **structurally zero**: the real matrix does not vary with the seed, so nothing in that arm can have spread. A tolerance of 0 makes every point "off the curve" and every point "resolvable". The legitimate uncertainty for a fixed dataset is a **bootstrap over people**, not a spread over seeds |
| 97d | **the hypothesis this actually points at, which is neither of the two I was testing** | **A CONTINUUM, not a mixture.** In the plants the top-5% *is* the carrier group — a discrete subpopulation with block-specific elevation, hence high between-block sd. In the real data the top-5% is the **tail of a continuous distribution**, and selecting on the mean averages block-specific spikes away. Low sd at that mean is what a **gradient** looks like, not what a **distinct minority** looks like |

**The fork this opens is ontological and it is the right one.** `#95` licensed "a minority concentrating
on rare options". `#97d` says the word **minority** may be wrong — the same statistic is produced by a
*continuum of rare-option affinity* with no distinct subgroup at all. Those are different objects,
they imply different mechanisms, and the discriminating design is a **mixture plant vs a gradient
plant**, matched on (mean, sd) and compared on the **shape of the whole distribution** rather than one
tail quantile.

---

## Entry 98, added by `E01·A11·R144` — neither family fits, and the two failures bracket the answer from opposite sides

`#97d` opened the fork: is the rare-option signal a **distinct minority** or a **gradient across
everyone**? Both families tuned to reproduce the real p50 and p95, then tested on the held-out shape.
Uncertainty from a **bootstrap over people** (`#97c`), which is tiny at n ≈ 12,000: sd(p50) = 0.0018,
sd(p95) = 0.0035.

**Elevation above the fixed-margin null, by quantile:**

| | p50 | p90 | p95 | p99 |
|---|---:|---:|---:|---:|
| null | 0.7011 | 0.8549 | 0.8990 | 0.9859 |
| **real** | **0.6999** | 0.8893 | 0.9520 | 1.0706 |
| **real − null** | **−0.0012 (AT the null)** | **+0.034** | **+0.053** | **+0.085** |

| # | Claim | Verdict |
|---|---|---|
| 98a | **the GRADIENT family** | **EXCLUDED, decisively.** Its weakest setting misses the real (p50, p95) by **44 bootstrap sd**. A trait spread across everyone necessarily **lifts the median** — and the real median sits **exactly at the null** (−0.0012, well inside 1 sd). Whatever this is, **the typical person does not have it** |
| 98b | **the sharp MIXTURE family** | **ALSO EXCLUDED.** Its best member (3% of people, 4 swaps) matches p50 and p95 to 3.7 sd — genuinely close — and then **overshoots the held-out tail shape by 4×**: s99 = **6.21** against the real **1.87**. A 3% carrier group produces nothing until p97 and then a cliff; the real data is **already elevated at p90** |
| 98c | **what the two failures bracket** | **A BROAD, MILD minority.** The elevation is **zero at the median**, appears by **p75–p90**, and grows smoothly to +0.085 at p99. Too broad to be 3% of people; too concentrated to be everyone. `#95`'s word *minority* survives `98a`; its implied narrowness does not survive `98b` |
| 98d | **the shape, stated so it is not over-read** | The real profile `(s90, s95, s99) = (0.958, 1.275, 1.875)` sits very close to the **null's** `(0.933, 1.200, 1.727)`. The distribution is only **mildly** heavier-tailed than chance in shape; the signal is in the **location of the upper quantiles**, not in a dramatic change of form |
| 98e | **gate (a)'s failure, which is correct and not a nuisance** | With n ≈ 12,000 the bootstrap sd is 0.002–0.004, so a 3-sd tolerance is **0.6% relative**. Neither family reaches it. **That is the design working**: a tolerance set by the data's own precision refuses models that are merely close, and both models here are merely close |

**The bracketing is the finding, and it names the next sweep exactly.** Carrier fraction between 3%
and 100%, intensity at the low end (1 swap), tuned so that **p50 stays at the null while p90 lifts by
+0.034** — those two constraints together pin the fraction, because a larger fraction at lower
intensity trades off along exactly that axis.

---

## Entry 99, added by `E01·A11·R146` — the check eleven rounds never ran: the distribution is symmetrically WIDER, not one-sidedly elevated

Every round in this arc measured **p50 and upward**. `#94`/`#95` licensed *"a minority picks rarer
options"*; `#98` bracketed its size; `#R11` found no carrier model fits, because every random-carrier
plant lifts the median while the real p50 sits at zero. **A much simpler world produces all of that,
and measuring it costs one line: what does the LOWER tail do?**

**Real elevation above the fixed-margin null, across the whole range:**

| p1 | p5 | p10 | p25 | **p50** | p75 | p90 | p95 | p99 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **−0.0425** | **−0.0356** | **−0.0327** | **−0.0166** | **+0.0003** | +0.0160 | +0.0339 | +0.0527 | +0.0851 |

Bootstrap sd (300 resamples over people): 0.0048 · 0.0027 · 0.0031 · 0.0022 · 0.0019 · 0.0019 ·
0.0032 · 0.0037 · 0.0097. **The lower-tail depression is 11–13× its own sd** — as resolvable as the
upper elevation, and pointing the other way.

**Symmetry ratio (upper elevation ÷ lower depression) = +1.27.**

| # | Claim | Verdict |
|---|---|---|
| 99a | **`#95`'s "a minority concentrating on rare options"** | **WITHDRAWN AS WORDED.** The distribution is **symmetrically wider** than its null: the top picks rarer options *and the bottom picks commoner ones*, by comparable amounts, with the median exactly unmoved. That is **a person-level trait with variance in both directions**, not a subgroup |
| 99b | **what survives, and it is not smaller** | **A continuous person-level parameter of rare-option affinity, with real spread.** Some people systematically pick uncommon options; others systematically pick common ones; the fixed-margin null has neither. `#95`'s *existence* claim stands — the null cannot produce this — only the word *minority* fails |
| 99c | **why eleven rounds missed it** | **The statistic was designed to find a tail, so every round looked at a tail.** `#91` framed the target as a minority invisible to variance-explained methods, `#92`–`#98` inherited that frame, and none of them printed a quantile below the median. **The frame chose the measurement, and the measurement could only confirm the frame** |
| 99d | **my two controls in this round** | **BOTH BROKEN — the twenty-fifth mis-specified design element.** They selected a **random** quarter of people rather than a *ranked* one, so both lifted every quantile including p1 (+0.0105); and the width control's down-swap needed carriers to *hold* rare picks to give up, which people in a curveball null mostly do not. Gates (a) and (b) both fail, correctly |
| 99e | **…and why `99a` stands anyway** | **The symmetry is a DIRECT MEASUREMENT of the real data against its own null — no model is fitted and no control is needed to read it.** What the broken controls forbid is the *formal* comparison of 1.27 against a calibrated pure-width prediction of 1.0. What they do not touch is that both sides moved, by 11–13× their bootstrap sd, in opposite directions |

**And this lands on the project's original question.** A continuous person-level parameter governing
*which* options are endorsed, independent of *how many*, is exactly the shape of Ivan's model B —
`v_i(s,c,t) = w_i(c,t)ᵀ h(s)`, an individualised readout weight. Not a fetish-carrying subgroup; a
**dimension everyone has a value on**.

---

## Entry 100, added by `E01·A11·R147` — a reliable person-level trait in WHICH options are endorsed, independent of HOW MANY

`#99b` proposed *"a continuous person-level parameter of rare-option affinity"*. That is a claim about
a **trait**, and a trait must be **reliable** and must not be something already measured wearing a new
name. Split-half over disjoint sets of **blocks**, Spearman-Brown corrected, with a planted-trait
ceiling and a fixed-margin floor:

| world | raw split-half `r` | **raw reliability** | after residualising on pick count | **residualised reliability** | corr(S, picks) |
|---|---:|---:|---:|---:|---:|
| fixed-margin null | 0.366 | **0.536** | **−0.011** | **−0.022** | **+0.719** |
| **real** | 0.481 | 0.649 | **0.275** | **+0.432** | +0.608 |
| planted trait | 0.726 | 0.842 | 0.712 | **+0.832** | +0.242 |

n = 7,316 people with ≥6 blocks; 6 half-splits × 3 seeds.

| # | Claim | Verdict |
|---|---|---|
| 100a | **gate (b), "the null is unreliable"** | **FAILED ON THE WRONG COLUMN — the twenty-sixth mis-specified design element.** It read the *raw* reliability (0.536) when the design had **pre-specified** the residualised one precisely for this confound. **Curveball preserves each person's pick count exactly**, and surprisal correlates with count at **+0.719**, so the null's raw reliability is pick-count reliability. On the pre-specified column the null returns **−0.022** and the gate **PASSES** |
| 100b | **the trait** | **CONFIRMED.** Residualised split-half reliability **+0.432**, against a floor of **−0.022** and a planted-trait ceiling of **+0.832**. It is **52% of the way to a cleanly planted trait** and 20× its own floor |
| 100c | **is it just breadth?** | **NO.** 67% of the raw reliability survives removing pick count (0.649 → 0.432), and the correlation with pick count is **lower in the real data (+0.608) than in the null (+0.719)** — the association with count is mostly *structural*, not the signal |
| 100d | **what this establishes about the original question** | **A reliable individual-difference dimension in WHICH options a person endorses, independent of HOW MANY.** That is `v_i(s,c,t) = w_i(c,t)ᵀ h(s)` measured — an **individualised readout weight with real, reliable between-person variance**. Not a category detector, not a subgroup, not response style alone |

**What reliability cannot do, stated because `100d` is the strongest claim in the arc.** A stable
**response style** is also a reliable person property. Reliability establishes that *something stable
about which options a person picks* exists and is not pick count; it does not establish that the
stable thing is erotic preference rather than a way of answering questionnaires. **The external
correlate is the missing test** — and unlike everything else in this arc, that one needs variables
outside the endorsement matrix, which this release does carry.

---

## Entry 101, added by `E01·A11·R148`+`R15` — the trait tracks WHEN interests were acquired, more than it tracks personality

`#100` confirmed a reliable person-level trait in *which* options are endorsed and named what
reliability cannot settle: a stable **response style** is also a reliable person property. The two
readings make opposite predictions about variables **outside** the endorsement matrix.

`R14` ran it and both gates failed **on my errors**, both worth recording:

| # | error | class |
|---|---|---|
| 101a | the null's residualised affinity still correlated with agreeableness (−0.052) and onset (−0.048). Residualising on `picks + log(picks)` is **linear**; the leftover is non-linear and still tracks count | **27th mis-specified design element** |
| 101b | I compared **raw** breadth correlations (0.021 mean) against `#17`/`#23`'s **disattenuated** ~0.10. Different scales. On the same scale breadth–openness is **+0.0712 raw**, which matches — the gate failed on my arithmetic | **28th** |

`R15`'s fix is not a better regression: **curveball preserves each person's pick count exactly**
(asserted), so `affinity = S_real − mean(S_null)` is count-matched **per person, by construction**.

**Both gates now pass** — null-minus-null max |r| = **0.0187**; breadth reproduces `#17`/`#23` on the
same scale (**0.0712** vs ~0.075 implied).

| variable | affinity `r` | z | breadth `r` | null−null |
|---|---:|---:|---:|---:|
| **biomale** | **+0.0930** | **8.4** | +0.1141 | −0.008 |
| **mean onset age** | **−0.0838** | **7.1** | −0.1596 | +0.001 |
| agreeable | −0.0555 | 4.9 | +0.0001 | −0.013 |
| conscientious | −0.0470 | 4.2 | −0.0124 | −0.007 |
| extroversion | −0.0357 | 2.9 | +0.0176 | +0.011 |
| neuroticism | −0.0318 | 2.7 | −0.0055 | +0.019 |
| openness | +0.0234 | 2.0 | +0.0712 | −0.014 |
| age | −0.0157 | 1.4 | +0.0416 | +0.008 |
| powerlessness | −0.0058 | 0.5 | +0.0205 | −0.001 |

| # | Claim | Verdict |
|---|---|---|
| 101c | **the fork, `#100`'s open question** | **LEANS EROTIC, and only leans.** Onset age |r| = 0.0838 against a personality-block **mean** of 0.0332 (**2.5×**) — but against the **largest** personality correlate (agreeable, 0.0555) it is **1.5×**. A mean-versus-single comparison flatters the conclusion; the honest number is 1.5× |
| 101d | **the direction, which is interpretable** | **Higher rare-option affinity ↔ EARLIER acquisition** (r = −0.084; disattenuated −0.127). People whose endorsements concentrate on uncommon options report acquiring their interests younger. This is the first link between arc `A11` and the fifteen rounds `A03` spent on the acquisition schedule |
| 101e | **the largest correlate overall** | **Sex: +0.0930** (disattenuated +0.141) — men's endorsements concentrate on rarer options. Larger than onset and larger than anything in the personality block |
| 101f | **the size, so it is not over-read** | Every correlation here is |r| ≤ 0.093, i.e. **under 1% of variance**, and the *breadth* nuisance this project projects out everywhere has a **larger** correlation with onset (−0.160) than the trait does (−0.084). **The trait has an external anchor; it is not a strong one** |

**What is settled and what is not.** Settled: the trait is not pick count (`#100c`), it has external
correlates that a matched null does not (`101` gates), and its largest non-demographic correlate is
acquisition age rather than personality. **Not settled**: 1.5× is a lean, not a separation, and a
response style that happens to track onset age would produce the same table.

---

## Entry 102, added by `E01·A11·R150` — the onset link does not survive its own matched null, and my gate checked the null against a number instead of against the effect

`#101d` reported that rare-option affinity tracks earlier acquisition at **r = −0.0838, z = 7.1**, and
called it the first link between arc `A11` and `A03`'s fifteen rounds on the schedule. This round
tested it against the nuisance `#101f` flagged, on a population restricted to people with **≥6 onset
entries as well as ≥6 blocks**.

| measure | raw | bootstrap sd | partial (picks, log picks, answer count) |
|---|---:|---:|---:|
| **affinity** | **−0.0302** | 0.0132 | −0.0292 |
| breadth | +0.0111 | 0.0128 | +0.0118 |
| **null − null** | **−0.0275** | 0.0139 | −0.0272 |

Split-half (affinity from half the blocks vs onset from the other half of the categories):
**−0.0205**.

| # | Claim | Verdict |
|---|---|---|
| 102a | **the verdict this round printed, `SURVIVES`** | **WRONG — the twenty-ninth mis-specified design element, and the worst kind.** Gate (a) asked whether the null-minus-null correlation was **below an absolute 0.03**. It was, by 0.0025. **A negative control must be compared to the EFFECT it is supposed to be smaller than, never to a number I chose.** Against the effect it is **91% of it** |
| 102b | **`#101d`'s onset link** | **DOWNGRADED to `UNVERIFIED`.** At −0.0302 against a matched null of −0.0275 with overlapping bootstrap intervals, the effect and its own control are **not distinguishable** |
| 102c | **why `#101` got −0.0838 and this got −0.0302** | `#101` included everyone with ≥6 **blocks**, computing mean onset over however many onset entries they had — **including people with one or two**, whose mean onset is both noisy and confounded with how much they answered. Restricting to ≥6 onset entries **cuts the correlation to 36% of its size**. `#101`'s null−null was +0.0007 on that looser population and −0.0275 here, so the restricted estimate is the better-controlled one and the looser one was inflated |
| 102d | **what survives from `#101`** | `#100`'s **reliability** result is untouched (0.432 residualised, floor −0.022) — it uses no external variable. `#101e`'s **sex** correlate (+0.0930, z 8.4) is 3× the onset one and was not tested here. What falls is specifically **the onset link and the "erotic-parameter side" reading built on it** |
| 102e | **the pattern across `#101a`, `#101b`, `102a`** | **Three gate failures in two rounds, all in the comparison rather than the computation.** A control residualised on the wrong functional form; a reference compared across scales; a null compared against a constant. **The measurements were right every time and the comparisons were wrong every time** |

**The arc's honest state.** A reliable person-level trait in *which* options are endorsed exists
(`#100`, floor-controlled, 20× its floor). Its **only** external correlate that has survived a matched
null is **sex**. The onset link — the one finding that connected this arc to the rest of the project
— **is withdrawn pending a design where the null is compared to the effect.**

---

## Entry 103 — `lib/gates.py`: the comparison rules as code, validated by replaying every failure it was written for

`#102e` named the class behind twenty-nine mis-specified design elements: **every gate in this project
compares two numbers, and every failure was in HOW, not in WHAT.** Nine of the ten catalogued cases
are a comparison taking the **wrong second argument**. That is infrastructure, not a thirtieth patch
(`P7`).

`lib/gates.py` makes the second argument **required**:

| method | the failure it forecloses |
|---|---|
| `negative_control(name, null, effect=…)` | **`#102a`** — no absolute threshold exists; a null is judged against the effect it must be smaller than |
| `positive_control(name, planted, floor=…, spread=…)` | **`#88a`** (planted below its own MDE, 3×), **`#78c`** (a criterion the instrument was known to violate) |
| `same_scale(name, mine, theirs, scale=…)` | **`#101b`** — `scale` is a required string with no default; a comparison whose scale nobody wrote down is the one that compared raw to disattenuated |
| `resolvable(name, effect, spread)` | **`#97c`** — refuses a spread that is zero or non-finite, and says to bootstrap over units rather than seeds |
| `no_sign_crossing(name, series)` | **`#83d`** (a ratio across a sign change → 1.7 × 10¹⁰), **`#79f`** (a threshold on a sum that hides a sign) |
| `covers_every_arm(name, checked, arms)` | **`#79e`** — a control that examined one arm of a two-arm design and passed |
| `asserted(name, condition, detail)` | **`#96a`** — a condition stated in prose must be a boolean here, or it was never tested |

**Validated the only way infrastructure can be**: every historical failure replayed through it with
its real numbers.

| case | originally printed | library returns |
|---|---|---|
| `#102a` null 91% of the effect | **PASS** (|null| < 0.03) | **FAIL** — "null is 91% of the effect" |
| `#79e` one arm of two checked | **PASS** | **FAIL** — "MISSING ['logit']" |
| `#83d` cliff ratio across a sign change | **PASS** (1.7e10) | **FAIL** — "signs CROSS ZERO" |
| `#96a` "at comparable mean" in prose | **PASS** | **FAIL** — means differ by 0.2881 |
| `#97c` tolerance from a zero spread | verdict from a 0 tolerance | **FAIL** — "STRUCTURALLY ZERO" |
| `#78c` criterion the instrument violates | FAIL, undiagnosed | **FAIL** — "headroom −0.0497" |
| `#79f` threshold on a sign-crossing sum | FAIL, wrong reason | **FAIL** — "signs CROSS ZERO" |
| **`#101b`** raw vs disattenuated | **FAIL** (my arithmetic) | **PASS** — correct once the scale is declared |

**7 of 8 flipped to FAIL; the eighth flipped to PASS, and that is also correct** — `#101b`'s original
failure was mine, not the data's, and declaring the scale resolves it. A library that only ever said
FAIL would be a library that cannot be wrong, which is the failure mode this project has catalogued
twenty-nine times.

**What it does not fix.** It cannot know whether a control is the *right* control — `#101a`'s
residualisation was linear where the confound was not, and no signature check catches that. It
forecloses the comparison errors, which are nine of ten; **the tenth is still mine**.

---

## Entry 104, added by `E01·A11·R151` — CLOSURE: the arc's three surviving claims re-read through the instrument built from its own failures

Labelled **Closure** (`P0`): it protects existing conclusions rather than separating worlds. Its value
is that the protection is performed by `lib/gates.py`, which `#103` validated by catching **7 of 8**
historical failures I had passed by hand. Read from the committed artifacts — no recomputation, so
nothing can drift.

**`#95` — is the rare-option signal licensed?**

| gate | result |
|---|---|
| no-op arm reproduces the real comparison | **PASS** — \|+0.0494 − 0.0503\| = 0.0010 < 0.0027 |
| plant at 2 swaps/block clears the floor | **PASS** — +0.1808 > +0.0514, headroom +0.1293 |
| plant ladder monotone and one-signed | **PASS** — [0.0494, 0.0815, 0.1808] |
| real effect resolvable at p95 | **PASS** — **37.5× its own spread** |

**`#99` — is the widening symmetric?**

| gate | result |
|---|---|
| upper elevation (p90+p95) | **PASS** — +0.0865, **17.8×** its spread |
| lower depression (p5+p10) | **PASS** — +0.0683, **16.8×** its spread |
| median unmoved | **PASS** — \|+0.0003\| < 0.0038 |
| ratio within 2× of symmetry | **PASS** — up/dn = **1.27** |

**`#100` — is the trait reliable, and is it pick count?**

| gate | result |
|---|---|
| null reliability vs the effect | **PASS** — null is **5%** of the effect (not "below 0.15", which is what `#102a` taught) |
| planted trait recovered | **PASS** — +0.8317, headroom +0.8462 |
| real residualised reliability | **PASS** — +0.4317, **90.5×** its own spread |
| every arm read on the same column | **PASS** — the `#100a` failure, now structural |
| >50% of raw reliability survives removing picks | **PASS** — 0.4317 > 0.3247 |

**All three survive.** And they survive under the *corrected* comparison rules, not the ones I wrote
by hand: `#100`'s null is now judged against its effect (`#102a`'s lesson) rather than against a
threshold, and its arms are checked for completeness (`#79e`'s).

**Arc `A11`, final state:**

| claim | status |
|---|---|
| a rare-option signal exists that a margin-preserving null cannot produce | **LICENSED** (`#95`, `#104`) |
| it is a **symmetric widening**, not a minority — both tails move, median unmoved | **CONFIRMED** (`#99`, `#104`) |
| it is a **reliable person-level trait**, and it is not pick count | **CONFIRMED** (`#100`, `#104`) |
| its only external correlate surviving a matched null | **sex, +0.093** |
| it tracks acquisition age | **WITHDRAWN** (`#102`) |
| a *minority* structure carried by ≤30% of people | **still invisible to every variance-explained method** (`#91`) |

---

## Entry 105, added by `E01·A10·R134` — the direct measurement is algebraically impossible, which is why `#90` had to be indirect

`#90`'s **±23.7 pp** is the top row of the README and the only number in the standing table that is
**inferred through a model** rather than measured. `#90` named that asymmetry in its own scope
paragraph. This step attacked it with what looked like a direct route: binomial noise is independent
across cells, and curveball preserves **both** margins exactly, so `mean r²(real) − mean r²(null)`
should cancel the noise and the person-estimation error and leave the interaction variance.

**It returned exactly 0.00 pp — including on a planted world carrying 30.9 pp of interaction.**

The reason is an identity, not a weak instrument. For a binary matrix with `p_j` the column mean and
`b_i` an out-of-block person effect:

| term of `mean_ij (M_ij − p_j − b_i)²` | fixed by |
|---|---|
| `mean(M²) = mean(M)` (M is 0/1) | margins |
| `−(2/nm) Σ_j p_j · (column sum)_j` | **column sums** |
| `−(2/nm) Σ_i b_i · (row sum)_i` | **row sums** |
| `mean((p_j + b_i)²)` — contains no `M` | margins |

**Every term is a function of the margins alone.** Verified numerically: three trials with
31–35 pp of planted interaction, real MSR and null MSR agree to **machine precision, difference
exactly 0.000e+00**.

| # | Claim | Verdict |
|---|---|---|
| 105a | **this round's design** | **A DERIVATION, not a measurement — the `realstat` arithmetic trap, fired on my own step.** *Could this have come out otherwise?* No: the algebra forces it. Same class as `#69`, where the item main effect turned out to be forced by prevalence dispersion at R² = 0.994 |
| 105b | **`#90`'s indirectness** | **VINDICATED, and upgraded from a caveat to a theorem.** `#90` did not use plant inversion out of convenience — **a margin-cancelling direct route provably does not exist.** The scope paragraph should now read *"inferred because direct estimation is algebraically impossible"*, not *"inferred, and that asymmetry is not removable here"* |
| 105c | **the general law, which is worth more than the round** | **Any statistic that is a function of the margins alone is invisible to a margin-preserving null, at any interaction magnitude.** First moments of residuals are such statistics. Only statistics that are **not** margin-determined can see the interaction: a low-rank **fit** (`#88`), or the **distribution** of a per-person score across people (`#95`) — *its quantiles, never its mean* |
| 105d | **what `105c` retroactively explains** | Why `#93`'s per-person misfit `T` died against the correct null (its mean is margin-determined); why `#95` had to read **p90/p95/p99** rather than the mean surprisal; why every early round that compared averages against a permutation null found nothing. **Three separate dead ends in this project were one theorem** |
| 105e | **the gate library's report** | `lib/gates.py` printed 3 PASS and 1 FAIL on a run where every number was zero. `negative_control(0, effect=0)` passes vacuously, and `no_sign_crossing([0,0,0])` passes. **A degenerate-input guard is missing** — the library needs to refuse a comparison whose effect is exactly zero rather than evaluate it |

**The sentence I can no longer write:** *"the interaction magnitude is inferred rather than measured,
and a direct measurement would settle it."* **There is no direct measurement to run.**

---

## Entry 106, added by `E01·A12·R157` — 获取顺序留下了痕迹:两个最后都喜欢 A 和 B 的人,会因为「先得到哪个」而在其他方面不同

Ivan 最初的先验有三步:`普通表征 → 个体化性价值读出 → 递归性表征重构`。项目测了前两步,第三步被归档成"需要纵向数据"就再没碰过。**那个归档是错的** —— 递归有横断面签名,而且它是心理学问句不是统计问句:

> 如果一个已获得的兴趣**重塑**表征,那么两个最后都喜欢 A 和 B 的人,应该因为**先得到哪个**而在其他所有方面不同。纯读出下顺序不可能有影响 —— 你两个都有了,权重就是权重。

120 个类别对,中位 n = 1,184。用**其他所有偏好**预测"A 先还是 B 先",留出 AUC,基线里已含性别·年龄·广度·早熟度:

| | 值 |
|---|---:|
| 仅协变量 AUC | 0.5617 |
| + 其他所有偏好 | 0.5671 |
| 原始增量 | +0.0054 |
| 分层置换零 | **−0.0200** |
| **偏移校正后** | **+0.0254** |
| SE | 0.0023 → **11.1×** |
| 逐对为正 | **104/120 (86.7%)** |
| 自助 95% | [+0.0209, +0.0299] |

| # | 结论 | 判定 |
|---|---|---|
| 106a | **顺序携带信息** | **确认。** 偏移校正后 +0.0254,11.1× 自身 SE,120 对里 104 对为正,自助区间不含零 |
| 106b | **幅度,用种植标定** | 种植阶梯相对自己 g=0 地板是 **0 / +0.0064 / +0.0255**,真实效应 **+0.0254** —— **正好落在 g=0.25 上**。约 **四分之一**的顺序标签可以从"其他所有偏好"重构 |
| 106c | **我的 gate 用错了形状 —— 第三十个** | 我用 `negative_control` 问"零是否够小"。**这个零不该是零** —— 两臂模型容量不同(~70 个偏好预测器 vs 4 个协变量),零携带的是**过拟合代价**,系统性为负。用"零应≈0"的比较去判一个偏移零,会让一个真实效应 FAIL。库里加了 `offset_control`,选择规则写在 docstring:**问"零应该是零吗?"是 → `negative_control`;不是,它有已知的系统方向 → `offset_control`** |
| 106d | **最强替代解释,写在跑之前,现在必须付账** | **稳定的"类型"同时驱动顺序和内容。** 一个偏关系型的人,既会更早获得关系性兴趣,也会有关系性偏好 —— 顺序通过**类型**预测轮廓,而不是通过重塑。协变量里有早熟度,但没有"你是哪一类"。**这不能被这个设计排除** |
| 106e | **能排除它的判别设计,便宜且下一轮就做** | 把类别对分成**同类**(都具体 / 都关系)和**跨类**。类型驱动 → 效应集中在跨类对;重塑驱动 → 同类对里也在。**这是 `106d` 的分离器** |

**心理学上现在能说的**:*你按什么顺序获得性兴趣,和你最终喜欢什么,不是两件独立的事。* **不能说的**:是顺序**造成**了差异 —— 横断面数据不排除第三共因,而 `106d` 恰好指出了一个非常可信的第三共因。

---

## Entry 107, added by `E01·A12·R158` — 顺序的痕迹不是「你是哪一类人」:跨类对的效应被类型吸收殆尽,同类对的存活

`#106d` 在跑之前就写下了唯一严肃的替代解释:**一个稳定的类型同时驱动顺序和内容**。偏关系型的人既更晚获得关系性兴趣,也有关系性偏好。`#106e` 给了分离器,这一轮跑了它 —— 而且是两个分离器同时跑。

**种类划分写死在代码里,在看任何结果之前**:具体 = 物体/身体属性/物质,不需要建模他人意图(12 类);关系 = 需要建模他人意图、地位、同意或关系本身(14 类)。

| 对的种类 | 中位 n | 基线协变量下 | **直接吸收类型后** | 保留 | 逐对为正 |
|---|---:|---:|---:|---:|---:|
| **跨类**(78 对) | 1,321 | +0.0221 | **+0.0005** | **2%** | 83.3% |
| **同类**(62 对) | 1,337 | +0.0289 | **+0.0236** (6.8× SE) | **82%** | 88.7% |

"直接吸收类型"= 把每个人自己的**具体类平均起始年龄**、**关系类平均起始年龄**及其差放进基线协变量。

| # | 结论 | 判定 |
|---|---|---|
| 107a | **跨类对的顺序效应 = 类型** | **确认。** 加入类型协变量后从 +0.0221 掉到 **+0.0005**(0.2× SE)。`#106d` 的怀疑对跨类对完全正确 |
| 107b | **同类对的顺序效应不是类型** | **确认。** +0.0289 → **+0.0236**,保留 **82%**,**6.8× 自身 SE**。在两个关系性兴趣之间,先得到哪一个,仍然预测你其余的轮廓 —— 在已经吸收了"你有多关系型"之后 |
| 107c | **我的正对照第一次是错的 —— 第三十一个** | `g=0` 那一臂**就是真实数据本身**,种植是**替换**标签而不是叠加,所以在一个已有真实信号上种植不会增加。这是 `#94b` 那个自消种植的重演,我诊断过一次又犯了一次。重跑:先用分层置换毁掉真实顺序信号当地板,再在上面种植 |
| 107d | **修正后的正对照,并且方向重要** | 两类都单调、都能检测,而且**同类对更灵敏**(g=0.30 时相对地板 **+0.0656** vs 跨类 **+0.0564**)。**所以同类里如果是零,那会是一个可读的零** —— 它不是 |
| 107e | **这个设计仍排除不掉什么** | 一个作用在**比"具体/关系"更细的层次**上的第三共因。同类里仍可能有子类型。本设计只排除这一个层次,这是它的全部范围 |

**心理学上现在能说的,而且这是全项目最接近 Ivan 第三步的一句**:

> 在两个**同类**的性兴趣之间,你先获得哪一个,仍然预示着你其余的整个偏好轮廓 —— 而且这不是"你本来就是那种人"的重复表述,因为"你是哪种人"已经被直接减掉了。

**仍然不能说**:顺序**造成**了它。横断面数据不排除更细粒度的共因。

---

## Entry 108, added by `E01·A12·R159` — 剂量测不出来,而且原因是分层本身:随机砍掉 2/3 的人就能杀死这个效应

`#107e` 留下最后一个对手:一个作用在比"具体/关系"更细层次上的共因。这一轮不用更细的分类学,用**剂量** —— 共因不做的那个预测:

- 静态共因:"我本来就是 A 型" → 先得到 A,也长得像 A 型。**隔多久、几岁,都不该有影响**
- 表征重塑:A 在 B 到来之前的窗口里重塑表征 → **间隔越大重塑越多**,**发生越早可塑性越高**

**两个正对照都过。** 种一个真的随间隔递增的效应 → 测得出递增(−0.0144 → +0.0669 → +0.1008);种一个跟间隔无关的平效应 → **不产生假剂量**(0.0364 / 0.0669 / 0.0547)。所以分层本身不制造剂量。

**两个剂量都 FAIL:**

| | 低 | 中 | 高 | 高−低 | 2×spread |
|---|---:|---:|---:|---:|---:|
| 同类·间隔 | −0.0062 | −0.0029 | +0.0024 | **+0.0087** | 0.0232 |
| 同类·年龄(低龄−高龄) | +0.0045 | +0.0025 | +0.0049 | **−0.0004** | 0.0280 |

**在相信这个"无剂量"之前查了功率,而功率是灾难性的。** 同一批对、同一个 `eff()`,只改样本量,**随机取三分之一(与剂量完全无关)**:

| | 全样本 (~1,300) | 随机 1/3 (~430) | 保留 |
|---|---:|---:|---:|
| **同类** | **+0.0121 ± 0.0062** | **−0.0007 ± 0.0092** | **−6%** |
| 跨类 | −0.0167 ± 0.0040 | −0.0212 ± 0.0061 | — |

| # | 结论 | 判定 |
|---|---|---|
| 108a | **剂量问题** | **UNVERIFIED,不是"没有剂量"。** 分层格子在真实效应的量级上完全失明 |
| 108b | **第三十二个设计缺陷,而且是一个新种** | **一个分层设计,它的每格功率从未与合并效应量对照过。** 正对照在 **+0.053**(g=0.25)开火,而真实合并效应是 **+0.0236** —— 控制在 **2.2 倍**于真实量级的地方开火,所以它**没有为真实量级的零背书**。`#88a` 的镜像:那次是种植低于 MDE,这次是**控制高于待测效应** |
| 108c | **诊断是直接测出来的,不是推断的** | 一个与剂量无关的随机三分之一就能把效应杀成 −6%。**问题是每格 n,不是剂量不存在** |
| 108d | **顺带的一个不一致,记下来** | 本轮合并同类效应 **+0.0121**,`#107` 是 **+0.0236**。差别:本轮要求 n≥900(只剩 22 个同类对)且只用 1 次置换抽样,`#107` 用 n≥400(62 对)和 2 次。**估计量本身跨轮有方差**,量级约 2 倍 —— 这个数字比我此前引用它时更不稳 |
| 108e | **保住功率的修法,下一轮做** | 别在对内分层。**在对间做**:每个对贡献一个用全样本算出的效应,和它自己的典型间隔/典型首获年龄,然后把效应回归到剂量上。n 变成对数,而每个效应保留全样本功率。同类对内部本来就有间隔变异 |

**这一轮没有产生任何关于世界的信念更新,只有关于仪器的。** 按 §0.2,这是**成本回收,不是产出** —— 它花了预算移除一个错误的"无剂量"结论,并且把那个结论标成 UNVERIFIED 而不是让它进账本当发现。

---

## Entry 109, added by `E01·A12·R161` — 效应"缩小"是因为我在两轮之间悄悄换掉了零的种类,而换掉的那个是错的

同一个量,我自己重估三次,一路缩小:`#107` **+0.0236** → `#108` **+0.0121** → `#109` **+0.0058**。

**旋钮全扫,40 个同类对,只改估计量的三个旋钮:**

| ndraw | 1 | 2 | 4 | 8 |
|---|---:|---:|---:|---:|
| **分层置换** reps=5 | 0.0229 | 0.0215 | 0.0218 | **0.0237** |
| **分层置换** reps=8 | 0.0223 | 0.0205 | 0.0205 | **0.0224** |
| **平置换** reps=5 | 0.0072 | 0.0063 | 0.0048 | **0.0041** |
| **平置换** reps=8 | 0.0059 | 0.0052 | 0.0036 | **0.0031** |

**跟抽样数和 reps 都无关 —— 两条线各自完全平。漂移 100% 是置换的种类。** 已知强度种植在 16 种配置下全部测出且稳定(极差 0.0113,均值的 1/4),纯置换数据处处 ≤ 0.0033。所以估计量本身是好的。

**哪个零是对的 —— 造一个 y 只依赖协变量、P 没有任何额外信息的世界,直接量 (full − base):**

| | 值 | 与真值的误差 |
|---|---:|---:|
| **真值**(P 无额外信息时的 full − base) | **−0.0239 ± 0.0021** | — |
| **分层置换**给出的偏移 | −0.0199 ± 0.0031 | **+0.0040** |
| **平置换**给出的偏移 | −0.0032 ± 0.0022 | **+0.0207** |

| # | 结论 | 判定 |
|---|---|---|
| 109a | **`#107` 的 +0.0236 站住** | **确认。** 用真值偏移校正后同类效应 **+0.019 ~ +0.027**。`#107` 用的是对的零 |
| 109b | **`#108` 的 +0.0121 和 `#109` 的 +0.0058** | **撤回,是错零的伪影。** 平置换把偏移低估了 **+0.0207** —— 几乎正好是两条路线的全部差距 |
| 109c | **第三十三个设计缺陷,而且是新种** | **我在两轮之间换掉了一个承重的设计元素,而且没注意到。** `#107` 用分层置换,`#108` 写 `eff()` 时我为了简洁用了 `y[rp.permutation(len(y))]` —— 一行,看起来只是重构,实际上换掉了零假设本身。**"重构"和"改变估计量"在代码里长得一模一样** |
| 109d | **为什么分层是对的** | 被检验的命题是"P 在**协变量之上**是否携带信息"。零必须**毁掉 P→y 而保留 COV→y**。分层置换在协变量分层内打乱,正是这个;平置换把两者一起毁掉,于是它测的过拟合代价是在一个基线塌到 0.50 的**另一个 regime** 里的 —— 而真实数据的基线在 0.56 |
| 109e | **分层置换也不精确,记下来** | 它误差 +0.0040(分层是粗的:性别 × 广度三分位 × 早熟三分位 = 18 格,连续协变量的层内变异被毁掉)。**最好的做法是直接用合成无信号世界的偏移**,而不是任何置换 |
| 109f | **`#108`/`#109` 的剂量结论全部作废** | 它们都是在错误的零上算的。剂量问题回到完全未测状态 |

**按 §0.2,这一轮的产出**:一个立住的声明(同类顺序效应 ≈ **+0.02**,零由合成无信号世界直接量出),一个工具(**合成无信号世界作为偏移的标定物**,比任何置换都准),和一个关闭的决定(**平置换在这类增量设计里是错的零,不要再用**)。

---

## Entry 110, added by `E01·A12·R162` — 顺序效应在对的零上是 +0.0285(8.6× SE);剂量再次被自己的对照挡下,而这次机制查清了

`#109f` 把剂量问题打回完全未测(`#108`/`#109` 都用了错零)。这一轮用 `#109e` 认定的真值偏移 ——
**合成无信号世界**:把 `y` 回归到协变量上生成 `y_synth`(保留 COV→y,毁掉 P→y),在它上面直接量 `full − base`。

| | 中位 n | 平均间隔 | 平均首获年龄 | **效应** | SE | 纯置换臂 |
|---|---:|---:|---:|---:|---:|---:|
| **同类** (60 对) | 1,337 | 4.68 | 13.12 | **+0.0285** | 0.0033 | +0.0079 |
| 跨类 (70 对) | 1,362 | 4.82 | 13.04 | +0.0068 | 0.0014 | +0.0081 |

| # | 结论 | 判定 |
|---|---|---|
| 110a | **同类顺序效应,第三次独立估计** | **+0.0285 ± 0.0033(8.6× SE)。** 与 `#107` 的 +0.0236、`#109a` 的 +0.019~+0.027 一致。保守取纯置换臂 +0.0079 作残余偏移,则 **≈ +0.021**。这是这个量最好的估计 |
| 110b | **跨类效应几乎消失** | **+0.0068**,而同类是 **4.2 倍**。`#107` 的分离在对的零上重现 —— 跨类的顺序痕迹本来就是"你是哪一类",而类型已经在协变量里 |
| 110c | **剂量:再次 UNVERIFIED,但这次机制是查清的** | **平种植产生了 +0.0079 ± 0.0034 的假斜率**(2.3×),而观测到的首获年龄斜率是 **+0.0080** —— **一模一样**。间隔斜率 −0.0062(1.9×)不可分辨。所以剂量读不出来,不是因为没剂量,是因为设计本身有一个与间隔挂钩的伪影 |
| 110d | **伪影的机制,直接测出来了** | **间隔大的对,顺序标签更一边倒**:`corr(平均间隔, 类别平衡度) = −0.597`,小间隔对平衡度 0.379 vs 大间隔对 0.291。AUC 在类别不平衡时的估计 regime 不同 —— 这就是平种植也能产生斜率的来源,与剂量无关 |
| 110e | **能修,而且共线性可接受** | 把**类别平衡度**放进对间回归。`corr(gap, bal) = −0.597`,不到 0.8,可以同时估。`corr(首获年龄, bal) = +0.410` 也在可修范围 |
| 110f | **合成偏移本身有小残差** | 纯置换 `y` 上效应应为 0,实测 **+0.0079**(效应的 28%,过了 `negative_control`)。原因:置换后的 `y` 没有协变量依赖,合成世界退化成抛硬币 —— 偏移又回到了"平置换 regime"。**合成偏移只在真实 y 上是对的,在置换臂上不是** |

**心理学上这一轮加了什么**:上一轮说"同类之间先获得哪一个,预示其余轮廓"。这一轮把那个数在**对的零**上钉到 **+0.021~+0.029**,并且确认**跨类的那一半确实只是"你是哪种人"**(4.2 倍差距)。**"隔多久""几岁发生"是否有剂量 —— 仍然不知道,而且现在知道为什么问不出来。**

---

## Entry 111, added by `E01·A12·R163` — 剂量在这套数据上问不出来,这是能力边界不是零;而我为了追它偏离了自己写下的 NEXT 三轮

`#110e` 说把类别平衡度放进对间回归就能修伪影。做了,机制确认了,但不够:

| | 上一轮 | 控制平衡度后 |
|---|---:|---:|
| 平种植的**假斜率** | +0.0079 | **+0.0050**(降 37%) |
| 观测到的首获年龄斜率 | +0.0080 | **+0.0072** |
| **假斜率 / 观测斜率** | **0.99** | **0.69** |
| 同类·间隔斜率 | −0.0062 (1.9×) | −0.0056 (1.4×) |

| # | 结论 | 判定 |
|---|---|---|
| 111a | **`#110d` 的机制确认** | 控制类别平衡度把假斜率降了 37%,方向和量级都对。伪影的来源诊断正确 |
| 111b | **剂量问题:能力边界,不是零** | 用尽最好的控制之后,伪影仍是观测斜率的 **69%**。**"间隔多大""几岁发生"在这套数据的这族设计上问不出来** —— 这是关于仪器的陈述,不是关于世界的 |
| 111c | **我偏离了自己写下的 NEXT,三轮** | `#107` 的 NEXT 写的是**时间不对称**:「先得到的那个应该比后得到的那个更能预测其余轮廓」。然后我连跑三轮**剂量**,一个不同的设计,而且三轮全是 UNVERIFIED。**frontier §3 的盆地信号早就亮了 —— N 步同方向、损失不再下降 —— 我没读它** |
| 111d | **按 §0.2 结账** | `#108` `#110` `#111` 三轮:关于世界的信念更新 **0**,关于仪器的 **3**(每格功率、错零、平衡度伪影)。这是**成本回收,不是产出**。其中 `#109` 是真正的产出(错零被抓出来,`#107` 被救回),但那是意外收获,不是我追剂量的目标 |

**换方向,回到 `#107` 自己写下的那一个,而且它是一个完全不同的设计:**

> 对一个类别 A,把人按"A 在他自己的序列里来得早还是晚"分组,然后问:**当 A 来得早时,A 的评分是否更能被其余偏好轮廓预测?**

如果 A 到来时重塑了表征,早-A 的人其余轮廓应该更"A 形"。**这个设计里没有间隔,没有类别平衡度,所以 `#110d`/`#111a` 的伪影结构上不存在。** 混淆是"早-A 的人可能就是更喜欢 A"—— 用 A 的评分做匹配即可,而评分是直接观测的。

---

## Entry 112, added by `E01·A12·R164` — 一个兴趣来得越早,它与这个人其余偏好结合得越紧;而且不是"喜欢得更强"

`#111c` 承认我为了追剂量偏离了 `#107` 自己写下的 NEXT 三轮。回到那一个,而且它是一个完全不同的设计 —— **没有间隔,没有类别平衡度**,所以 `#110d`/`#111a` 的伪影结构上不存在。

问题一句话:**一个兴趣来得早的时候,它是否与这个人其余的偏好结合得更紧?**

对每个类别 A,按 A 在**这个人自己的序列**里来得早还是晚分组,然后用**其余所有偏好**预测 A 的评分,比较两组的留出 R²。24 个能对上评分列的类别,中位 n = 2,884。

| 臂 | 早-A 组 R² | 晚-A 组 R² | **差** | SE | Δ评分 |
|---|---:|---:|---:|---:|---:|
| **真实** | **0.0675** | 0.0437 | **+0.0237** | 0.0088 | **0.0000** |
| 早/晚标签打乱 | 0.0656 | 0.0586 | +0.0070 | 0.0042 | 0.0000 |
| 种植 0.10 | 0.0852 | 0.0439 | +0.0413 | 0.0086 | — |
| 种植 0.25 | 0.1266 | 0.0464 | +0.0802 | 0.0096 | — |

**全部 gate 通过:**

| gate | 结果 |
|---|---|
| A 的评分在两组间已匹配(断言在代码里,不是散文里) | **PASS** — Δ = **0.0000** |
| 早/晚标签在匹配格内打乱 | **PASS** — 零是效应的 **29%** |
| 种植的重塑被测出 | **PASS** — +0.0802,余量 +0.0541 |
| 种植阶梯单调不跨符号 | **PASS** — [0.007, 0.041, 0.080] |
| 早减晚的差值可分辨 | **PASS** — **2.7× 自身展布** |

| # | 结论 | 判定 |
|---|---|---|
| 112a | **早到来的兴趣与其余偏好结合更紧** | **确认。** 留出 R² 从 **0.0437 升到 0.0675**,提高 **54%**。零(标签打乱)是 +0.0070 |
| 112b | **不是"喜欢得更强"** | **两组对 A 的评分按格匹配,差值 0.0000。** 跑之前写下的最强混淆("兴趣强 → 来得早")被直接排除,而不是被论证掉 |
| 112c | **幅度用种植标定** | 真实 +0.0237 落在种植 0.10(+0.0413)之下,约相当于 **0.06 的植入重塑强度** |
| 112d | **与 `#110a` 是两个独立的签名,都指向同一件事** | `#110a`:同类之间**先获得哪一个**预示其余轮廓(+0.021~+0.029)。`#112a`:一个兴趣**来得越早**,与其余轮廓结合越紧(+0.0237)。两个设计共享的只有"起始年龄"这个变量,伪影结构完全不同 —— 前者被间隔/平衡度伪影困住,后者根本没有这两个量 |
| 112e | **仍然不能说** | 因果方向。以及**"来得早"与"记得来得早"分不开** —— 一个人若把某个兴趣当作自己的核心,可能既报告更早的起始年龄,也让其余轮廓围绕它组织。这是本设计无法排除的,而且它本身就是一个心理学假设 |

**心理学上这是全项目最接近 Ivan 第三步的一句**:

> **越早进入一个人性欲版图的东西,越处在那张版图的中心** —— 其余的一切更能被它预测。而这不是因为他更喜欢它:两组对它的喜欢程度是逐格匹配的。

---

## Entry 113, added by `E01·A12·R165` — "早"不是"相对更爱";而绑定只在真正在乎的兴趣上发生

`#112e` 留下唯一对手:**回忆偏差** —— 一个人把 X 当作核心,既报告更早的起始年龄,也让轮廓围绕 X 组织。两个分离器同跑。

### 分离器一:直接用核心性替代(这个真的判别)

若"早"只是"核心",那么用**直接测到的核心性**匹配之后效应应该消失。核心性 = A 的评分在**这个人自己所有评分**里的百分位(`#112` 匹配的是**绝对**评分,这一个固定的是"A 对这个人有多突出")。

| 匹配方式 | 早组 R² | 晚组 R² | **差** | SE | Δ匹配量 |
|---|---:|---:|---:|---:|---:|
| 绝对评分(`#112`) | 0.0675 | 0.0437 | **+0.0237** | 0.0088 | Δ评分 **0.0000** |
| **相对核心性** | 0.0710 | 0.0475 | **+0.0235** | 0.0070 | Δ核心性 **0.0007** |

**保留 99%。** 六个 gate 全过,两个零(打乱标签)分别是效应的 29% 和 38%,种植在核心性匹配下仍被测出(+0.0944,余量 +0.0698)。

### 分离器二:评分带交互(这个不判别 —— 是我自己写错了)

| 带 | 真实差 | SE | 打乱地板 | **种植** | 仪器 |
|---|---:|---:|---:|---:|---|
| **低评分 (0–2)** | **+0.0045** | 0.0080 | −0.0079 | **+0.0614** | **能用** |
| **高评分 (3–5)** | **+0.0231** | 0.0041 | −0.0041 | **+0.0851** | **能用** |

高−低 = **+0.0186 ± 0.0090(2.1×)**,而且**两带的仪器都通过了自己的正对照** —— 交互是真的。

| # | 结论 | 判定 |
|---|---|---|
| 113a | **"早"不是"相对更爱"** | **确认。** 按人内相对核心性逐格匹配后保留 **99%**。`#112a` 在这个版本的对手面前完全站住 |
| 113b | **交互是真的,而且两带仪器都验证过** | **确认。** 低带 R² 是负的(早 −0.0365,晚 −0.0410),但那不是"没效应"——种植 +0.0614 说明仪器在那里能用 |
| 113c | **但交互不判别 —— 第三十四个设计缺陷** | `#112e` 我写的是"回忆偏差预测差随评分增大,重塑不预测这个交互"。**错。重塑同样预测它** —— 一个人不在乎的兴趣本来就没什么可用来重塑的。**两个世界做同一个预测,frontier §2 的"平行"** —— 这个动作对那个世界什么都没说,我却把它当判别式写进了 NEXT |
| 113d | **`#112a` 现在的确切范围** | 挡住了**"早 = 相对评分更高"**这个版本的回忆偏差(核心性匹配,99% 保留)。**没挡住**"早 = 身份认同上更核心,而这一点超出评分排名能测到的范围"。核心性代理是**评分秩**,它是部分测量 |
| 113e | **交互本身是一个关于世界的发现,和判别无关** | **绑定只在真正在乎的兴趣上发生。** 低评分带效应 +0.0045(0.6×,不可分辨)而仪器能用;高评分带 +0.0231(5.6×)。**一个来得早但你并不很喜欢的兴趣,不绑定任何东西** |

**心理学上这一轮加的一句话**:

> 早到来只有在**你真的在乎它**的时候才把其余的东西组织到它周围。一个早早出现却无关痛痒的兴趣,不改变任何布局。

**而且这不是"你更爱它所以显得更早"** —— 两组按人内相对核心性逐格匹配,差 0.0007。

---

## Entry 114, added by `E01·A12·R166` — 记忆把心爱的兴趣往前拉,约 0.74 年;这是这个数据集里第一个被直接量出来的回忆偏差

`#113c` 承认我上一轮的"判别式"不判别(两个世界做同一个预测)。这一轮的分离器**不碰绑定统计量**,直接量**起始年龄报告本身**,而且两个世界做**相反**的预测:

- **回忆偏差**:记忆把心爱的往前拉 → 评分高的项应系统性**早于**人群时间表对它的预测
- **报告准确**:偏离与评分无关

用 A03 花了 22 轮建立的人群时间表,**留一人**计算基准(一个人不参与自己的基准),再减掉这个人自己的整体早熟度。

| 臂 | 斜率(年 / 评分标准差) | SE |
|---|---:|---:|
| **真实** | **−0.2000** | 0.0101 |
| 评分在**同类别内**打乱 | −0.0004 | 0.0096 |
| 人为往前拉 0.3 年 | −0.2804 | 0.0100 |
| 人为往前拉 0.8 年 | −0.4147 | 0.0101 |

**全部 gate 通过**(阶梯写法修正后):偏离与早熟度不相关(+0.0227);打乱零是效应的 **0%**;种植阶梯单调 [+0.0804, +0.2148];种植 0.8 年被测出,余量 +0.1945;真实斜率 **19.8× 自身展布**。

| # | 结论 | 判定 |
|---|---|---|
| 114a | **心爱的兴趣被报告得更早** | **确认,19.8× SE。** 用种植标定:相当于把高评分项**往前拉了约 0.74 年**。零(同类别内打乱评分)是效应的 0% |
| 114b | **这是本数据集里第一个被直接量出来的回忆偏差** | 此前 `#61` 只证明"色情诱导"自述没有时序签名。这一条不同:它是一个**可量化的、在评分上单调的报告畸变**,而且 A03 的时间表本身就是在这个畸变之上建立的 |
| 114c | **它对 `#112a` 的威胁,而且比我原以为的更微妙** | 直觉说"`#112` 已按评分逐格匹配,格内拉前量恒定,所以无害"。**这个直觉是错的。** `#112` 的"早"定义是 `rel = 本项起始 − 该人其余项均值`。回忆偏差把**本项**拉前 f(评分本项),也把**其余项均值**拉前 mean(f(评分其余项))。格内前者恒定,**后者不恒定** —— 一个"心爱之物很多"的人,其余项被整体拉前,于是他这一项显得**更晚**。**"早"因此与"你有多少个心爱之物"挂钩,而那与轮廓相关** |
| 114d | **`#112a`/`#113a` 现在的状态** | **FLAGGED,未撤回。** `#112` 的协变量里有**广度**(评分 >0 的项数),但没有**平均评分**,而 `114c` 的通路走的是后者 |
| 114e | **可直接修,而且是这个畸变第一次能被用来修东西** | `#114a` 给出了畸变的**大小和形状**(0.2 年 / 评分标准差)。用它反向校正起始年龄:`own_corrected = own + 0.2 × z(评分)`,再重跑 `#112`。**如果效应存活,回忆偏差被排除;如果崩塌,`#112a` 必须撤回** |

**心理学上这一轮的一句话**:

> **人会把自己最爱的性兴趣记得比实际更早 —— 平均往前拉约九个月。** 你性欲史的起点,部分是被你现在的偏好重写过的。

**IMPOSSIBLE(写在前面的):** "记忆拉前"与"早获得的东西后来更被珍视"分不开 —— 两者都产生负斜率。本轮只能测这个斜率在不在,不能定方向。

---

## Entry 115, added by `E01·A12·R167` — `#112a` 撤回:那不是"早到来的东西在中心",是"什么都给高分的人"

`#114c` 预测了一条通路,`#114e` 说可以直接测。测了,而且**结论与我上一轮提交的相反**。

| 臂 | 早组 R² | 晚组 R² | **差** | SE | 保留 |
|---|---:|---:|---:|---:|---:|
| **raw**(与 `#112` 可比的基线) | 0.0619 | 0.0438 | **+0.0181** | 0.0087 | — |
| **仅校正起始年龄**(+0.2×z(评分)) | 0.0647 | 0.0493 | **+0.0154** | 0.0076 | **85%** |
| **仅把人均评分放进协变量** | 0.2192 | 0.2173 | **+0.0019** | 0.0106 | **10%** |
| **both** | 0.2301 | 0.2222 | **+0.0079** | 0.0092 | 44%(0.9×,不可分辨) |

**三个 gate 通过,第四个 FAIL:** 评分匹配成立(Δ=0.0000)· 打乱零是效应的 32% · **种植在 both 臂仍被测出(+0.0513,余量 +0.0348)** · **but 真实效应 0.9× 自身展布**。

| # | 结论 | 判定 |
|---|---|---|
| 115a | **`#112a` 撤回** | "越早进入版图的东西越在中心" —— **不成立**。加入**人均评分**后效应降到原来的 **10%**,both 臂 0.9× 不可分辨,而**正对照证明 both 臂能测**(种植 +0.0513,余量 +0.0348)。这不是功率问题 |
| 115b | **`#113a` 同样撤回** | 它匹配的是**人内相对核心性**,也没有控制**人均评分**。`#113a` 说"保留 99%"是真的 —— 但那只说明相对核心性不是那条通路,`#114c` 指出的才是 |
| 115c | **杀死它的不是回忆偏差校正** | 仅校正起始年龄保留 **85%**。**是"人均评分"这一个协变量。** 我在 `#114c` 里把通路描述成"回忆偏差经由其余项均值",但实测下来,起始年龄的畸变只值 15%,**剩下的 90% 是一个更平凡的东西:什么都给高分的人** |
| 115d | **它的机制,从数字上直接看得出来** | 加入人均评分把两组的 R² 从 ~0.06 抬到 **~0.22** —— 人均评分是单个项目评分的极强预测器(给什么都打高分的人,这一项也高)。而 `rel = 本项 − 其余项均值` 让"早"与人均评分挂钩。**"早"一直在部分地测"这个人整体给分有多高"** |
| 115e | **`#107`/`#110` 现在也有嫌疑,而且是同一条通路** | 那个设计的结局是"A 先还是 B 先",协变量里有广度、早熟度、具体/关系倾向,**没有人均评分**。而 `y = (起始A < 起始B)` 经由回忆偏差依赖 `评分A − 评分B`,轮廓又预测评分。**同一条通路,没测过** |
| 115f | **A12 现在还立着的只有一条** | **`#114a`:记忆把心爱的兴趣往前拉约 0.74 年(19.8× SE)。** 它独立测量,不经过任何绑定统计量,不受这条通路影响 |

**心理学上必须收回的那句话**:上一轮我说"越早进入一个人性欲版图的东西,越处在那张版图的中心"。**收回。** 在把"这个人整体打分有多高"减掉之后,那个效应只剩下一成。

**A12 这个 arc 目前的净产出**:一个立住的发现(**人会把最爱的性兴趣记得早约九个月**),一个撤回(**"早=中心"**),和一个待查(`#107`/`#110` 的顺序效应是否走同一条通路)。

---

## Entry 116, added by `E01·A12·R168` — 顺序效应挺过了杀死它兄弟的那个混淆,而且原因是设计结构上的

`#115e` 把 `#107`/`#110` 的顺序对设计放到同一条通路下受审:`y = (起始A < 起始B)` 经由回忆偏差依赖 `评分A − 评分B`,而轮廓预测评分,协变量里从没有过人均评分。

| 协变量集 | 效应 | SE | 保留 |
|---|---:|---:|---:|
| **base**(`#110` 原协变量) | **+0.0246** | 0.0033 | — |
| **+人均评分** | **+0.0256** | 0.0033 | **104%** |
| **+对内两项各自评分及其差** | +0.0162 | 0.0027 | 66% |
| **+both** | **+0.0171** | 0.0027 | **70%** |
| 纯置换 y | +0.0083 | 0.0028 | — |
| 种植 0.20 | +0.0648 | 0.0036 | — |

**全部 gate 通过**:零是效应的 48% · 种植被测出(余量 +0.0492)· **+both 臂 6.4× 自身展布**,零的种类已命名(合成无信号世界)。

| # | 结论 | 判定 |
|---|---|---|
| 116a | **`#107`/`#110` 站住** | 加入人均评分**和**对内两项评分之后保留 **70%**,**6.4× 自身展布**。`#115e` 的怀疑对这个设计**不成立** |
| 116b | **人均评分对它完全没有作用** | +0.0246 → **+0.0256**,一点没动。杀死 `#112a` 的那个协变量在这里是**惰性的** |
| 116c | **原因是设计结构上的,这是本轮最有用的一条** | `#112` 的"早" = `本项 − 该人其余项均值` —— 直接与**人层面的整体给分水平**纠缠。`#107`/`#110` 的 `y = (A 先于 B)` 是**对内比较**,这个人的整体水平在两项之间**相消**。**成对设计对人层面作答水平的混淆是构造上免疫的** |
| 116d | **被拿走的 30% 是可解释的,而且正是 `#114a`** | 对内两项各自的评分拿走 34%,人均评分拿走 0。这正是回忆偏差经由 `评分A − 评分B` 的那一份 —— **`#114a` 量出来的畸变,在这里被量成了这个效应的三分之一** |
| 116e | **零偏高,记下来** | 纯置换 y 给 +0.0083,是效应的 **48%**,刚过 `negative_control` 的 0.5 线。`#110f` 已经诊断过原因(置换后的 y 没有协变量依赖,合成世界退化成抛硬币)。**它让 `+both` 的真实量级可能更接近 +0.009 而非 +0.017** |

**心理学上现在能说的**:

> 在两个**同类**的性兴趣之间,**你先获得哪一个,仍然预示着你其余的整个偏好轮廓** —— 在减掉你对这两样各自有多喜欢、以及你整体给分有多高之后。

**而它的兄弟命题("越早进入版图的东西越在中心",`#112a`)在同一个混淆下倒了。** 两个命题看起来是同一件事的两种说法,实际上一个是**对内**比较、一个是**对人内均值**比较,而后者结构上就在测作答水平。

---

## Entry 117, added by `E01·A12·R169` — 零修对后效应站得更稳;而"先来的把其余偏好拉向自己"刚过线

### 第一件:`#116e` 的零修好了

`#116e` 记下松头:纯置换零是效应的 **48%**,`#110f` 已诊断原因(置换后的 y 没有协变量依赖,`eff()` 内部的合成世界退化成抛硬币,偏移掉回平置换 regime)。正确的零是 **`eff(y_synth)`** —— `y_synth` 由**真实 y 的协变量拟合**生成,保留 COV→y、毁掉 P→y,而 `eff()` 内部对它再算一次同样的拟合,**regime 自洽**。

| 臂 | 值 | 占效应 |
|---|---:|---:|
| **真实** | **+0.0159**(6.1× SE) | — |
| **修正零 `eff(y_synth)`** | **+0.0023** | **14%** |
| 旧零 `eff(y_置换)` | +0.0096 | 60% |

| # | 结论 | 判定 |
|---|---|---|
| 117a | **偏移无偏,`#116e` 的顾虑解除** | 修正零是效应的 **14%**,不是 48%。`#116e` 担心"真实量级可能更接近 +0.009 而非 +0.017" —— **不成立,+0.0159 是真量级** |
| 117b | **旧零一直在吃掉真实信号** | 置换零 +0.0096 里,只有 +0.0023 是偏移,**其余 +0.0073 是被 regime 错位制造出来的**。这解释了为什么 `#110`/`#116` 的数字彼此有 2 倍的抖动 |

### 第二件:方向(标为 **DESCRIPTION**,不做因果主张)

跑之前做了 `#113c` 的判别检查:**方向本身不判别** —— 重塑和"你本来就是 A 型"都预测"A 先的人更 A 向"。所以这一段只描述已确立效应的方向。它的非平凡版本是:**在减掉这两项各自的评分之后,A 先的人是否仍在其余偏好上更 A 向**。

| | 值 |
|---|---:|
| A 先组减 B 先组,在「A 向 − B 向」方向上的标准化位移 | **+0.0266 ± 0.0127** |
| 分层置换零 | +0.0067(效应的 25%) |
| 位移为正的对 | **29/45** |
| 可分辨性 | **2.1×** —— 刚过线 |

| # | 结论 | 判定 |
|---|---|---|
| 117c | **方向是"拉向先来的那个",但只是刚过线** | +0.0266,**2.1× 自身展布**,29/45 对为正。**这是本 arc 最弱的一条,应当按这个强度引用** |
| 117d | **最大的几对是可读的** | genderplay 先于 incest · nonconsent 先于 pregnancy · specific roles 先于 sensory · clothing 先于 incest —— 位移 +0.23 ~ +0.29。反向的:preferred object 之于 sensory、sadomasochism 之于 appearance |
| 117e | **第五次 pandas 访问器撞名,而这次撞的是我自己禁用名单上的第一个词** | `G.shift` → `DataFrame.shift`。`#80d` 我写下"任何列名不得与 DataFrame/Series 方法同名",并把 **`shift` 列在第一位**;`#93d` 撞 `T` 时我说"规则需要由不是我的东西来强制执行";**然后我又撞了一次同名单上的词。写下规则第三次没有阻止我违反它** |

**心理学上这一轮加的**:上一轮说"先获得哪一个预示其余轮廓"。这一轮说方向 —— **其余的偏好被拉向先来的那一个** —— 但只有 2.1×,是本 arc 最弱的一条,而且它在设计上不能区分"被拉过去"和"你本来就在那一侧"。

---

## Entry 118, added by `E01·A12·R170` — 方向命题从 2.1× 升到 3.1×,但不是因为我预测的理由;而它原来的正对照是一个恒等式

`#117` 的 NEXT 预先承诺:双倍 n 若不过 3× 就按 2.1× 引用并停止。跑了,过了 —— 但**两个诊断都跟我的预测无关**。

| 门槛 | 对数 | 位移 | SE | **强度** | 位移为正 |
|---|---:|---:|---:|---:|---:|
| **n≥400** | **68** | +0.0339 | 0.0110 | **3.1×** | 46/68 |
| n≥250 | 73 | +0.0318 | 0.0106 | 3.0× | 48/73 |
| 分层置换零 (n≥400) | 68 | −0.0022 | 0.0076 | — | 效应的 **6%** |

| # | 结论 | 判定 |
|---|---|---|
| 118a | **加强成功,`#117c` 按 3.1× 引用** | 位移 **+0.0339,3.1× 自身展布**,零是效应的 6%,46/68 对为正 |
| 118b | **但降低门槛几乎没作用 —— 我预测错了原因** | 400→250 只把对数从 **68 加到 73**,强度 3.1× → 3.0×(**略降**)。我在 `#117` 的 NEXT 里说"对数大约翻倍",实际增加 **7%** |
| 118c | **真正的原因是我自己在 `#117` 里的一个 cap** | `R13:116` 有 `if len(dirs)>=45: break` —— **同一个门槛(n≥400)下本来就有 68 对,我截断到 45**。加上 `R13:104` 的 `pref` 阈值 80(本轮 60),`#117c` 的 2.1× 是在一个**被我任意截断的子集**上算的。**第六次同类:一个我自己写的静默 cap 改变了结论**(`#73` `#74` 之后) |
| 118d | **`#117c` 的正对照是一个恒等式 —— 它不可能失败** | 我的种植是 `proj = proj + g*lab`,直接加在**投影**上,所以 `plant15 − plant0 = 0.1500` **精确等于种进去的量**。这是 `realstat` 的算术陷阱:*这个数可能是别的吗?不可能,代数强制它*。所以 `#117c` 实际上**从未有过正对照** —— 我上一轮说"补上",补的是一个不能失败的检查 |
| 118e | **生成式正对照,这次是真的** | 在**偏好数据层面**推 A 先组,然后**从头重算** `w` 和投影:g = 0 / 0.10 / 0.30 给 +0.0339 / **+0.1140** / **+0.2723**,相对地板增量 **+0.0801 / +0.2383**,单调且远超 2×SE。**仪器确实能看见一个生成出来的方向位移** |

**心理学上的那句话现在可以按更高的强度说**:

> 在两个同类的性兴趣之间,**其余的偏好被拉向先来的那一个** —— 位移 3.1× 自身展布,68 对里 46 对为正,在减掉这两样各自的评分、它们的差、以及全部人层协变量之后。

**仍然标 DESCRIPTION**:`#113c` 查过,方向本身不判别 —— 重塑和"你本来就在那一侧"都预测它。强度提高了,性质没变。

---

## Entry 119, added by `E01·A13·R171` — 回忆畸变随年龄增强,所以它不是纯粹的"报告时刻建构";而"伪影是否为零"一直是错的问题

新 arc。`#114a` 量出的回忆偏差(心爱的被报告早约 0.74 年)我一直只当干扰用。它本身是个心理学对象,而且两个机制做**相反**的预测:

- **记忆重构**:畸变来自回忆过程 → 回忆间隔越长重构越多 → **年龄越大斜率越负**
- **当下叙事**:畸变来自报告时刻把现在最爱的整合进自传 → 与流逝时间无关 → **与年龄无关**

**第一次跑,gate 正确地拦下了:** 平种植(与年龄无关的畸变)产生假趋势 **−0.0691**,真实是 **−0.0660** —— 几乎一模一样,读不出来。

**机制查清了,而且是可测的**:`sd(dev)` 随年龄单调扩大(跨 5 层相关 **+0.988**,最老/最年轻 **1.54 倍**)。斜率的单位是"年 / 评分SD",**尺度本身在变**,所以任何非零畸变在老龄层都显得更大。

**修法:层内标准化,让斜率无量纲。**

| 受访者年龄 | 真实 | 打乱零 | 种植(年龄依赖) | **平种植** |
|---:|---:|---:|---:|---:|
| 15.5 | −0.0505 | −0.0086 | +0.0940 | −0.2066 |
| 19.0 | −0.0737 | +0.0020 | +0.0004 | −0.2219 |
| 22.5 | −0.0751 | +0.0075 | −0.0667 | −0.2053 |
| 26.5 | −0.0766 | +0.0036 | −0.1271 | −0.1959 |
| 30.5 | **−0.0917** | +0.0006 | −0.1838 | −0.1942 |
| **年龄趋势** | **−0.0122 ± 0.0035(3.5×)** | +0.0026(0.8×) | −0.0954(27.4×) | **+0.0070** |

| # | 结论 | 判定 |
|---|---|---|
| 119a | **畸变随年龄增强** | **确认,3.5×。** 从 15 岁组的 −0.0505 到 30 岁组的 −0.0917,**十五年里几乎翻倍** |
| 119b | **伪影现在指向相反方向** | 层内标准化把平种植的假趋势从 **−0.0691** 压到 **+0.0070**,而且**符号与真实效应相反** —— 校正只会让效应**更大**(−0.0192),不会更小 |
| 119c | **但判别范围比我框的窄,必须收口** | 我把它框成"记忆重构 vs 当下叙事"。实际判别的是**时间依赖 vs 时间无关**:纯粹的报告时刻建构被排除了,但"叙事在多年里被反复复述而固化"同样随时间增长,**它没被排除**。结论应写成:**这个畸变是随时间累积的,不是在回答问卷的那一刻凭空生成的** |
| 119d | **"伪影是否为零"一直是错的问题 —— 第三十五个** | 我用 `asserted(|伪影| < 2·SE)` 去判,它以 **0.0001** 之差 FAIL(0.0070 vs 阈值 0.0069)。而那个伪影**符号与效应相反**,校正只会放大效应。正确的问题不是"它是不是零",是"**它能不能造出这个效应**"。库里加了 `artifact_cannot_explain(artifact, effect, spread)`:符号相反 → 过;同号但小于效应一半 → 过;同号且同量级 → FAIL。**回放 `#110c`(伪影 +0.0079 vs 效应 +0.0080)→ 正确地 FAIL** |
| 119e | **`check_coverage` 第一次真实使用就开火** | 24/26 类别(2 个 n<400)。这是合法纳入标准不是成本控制,但守卫要求**声明**,于是轮次输出里现在有"纳入 24/26 类别"这一行。**上一轮刚做的工具,这一轮就改变了输出** |

**心理学上的一句话**:

> 人把最爱的性兴趣记得更早这件事,**随时间加深** —— 十五年里几乎翻倍。你的性欲起点不是一次性被重写的,它在被持续地往前挪。

---

## Entry 120, added by `E01·A13·R172` — 机制分不开,而且原因是真实曲率根本不可分辨;我的两个 gate 又各错了一次

`#119c` 留下残余:"随时间累积"对**记忆衰退**和**叙事固化**都成立。分离器是**形状** —— 固化只拉最爱的那些(顶端,凸),衰退沿评分整体变陡(线性)。

设计时先撞上一个硬约束:`dev` 的基准是样本内的,所以 `mean(dev)` 在任一分层里被钉住,**高低两组之差是被构造强制的**。能读的只有沿评分 0–5 的**曲线形状**。

| 臂 | 曲率(二次) | 线性 | 曲率/线性 |
|---|---:|---:|---:|
| **真实** | **−0.00038 ± 0.00131** | −0.00212 ± 0.00105 | +0.18 |
| 打乱零 | +0.00033 | −0.00052 | — |
| 顶端种植(固化) | −0.00647 (4.9×) | −0.01185 | **+0.55 ± 0.12** |
| 线性种植(衰退) | −0.00457 (3.5×) | −0.03102 | **+0.15 ± 0.04** |

| # | 结论 | 判定 |
|---|---|---|
| 120a | **机制分不开 —— 能力边界,不是零** | **真实曲率是 0.3× 自身展布**,根本没通过可分辨性。要把它推到 2×,需要**约 48 倍样本量** |
| 120b | **不是因为两个形状不可分 —— 我的阈值写错了** | 我写 `abs(差) > 0.5` 判"两个种植是否可分"。用它们**自己的 SE** 判:+0.546 ± 0.122 vs +0.147 ± 0.043,**差 0.399 ± 0.130 = 3.1×,其实可分**。**又一个随手写的常数阈值**,`#102a` 的同一类,第三十六个 |
| 120c | **种植又一次远高于真实量级** | 线性项:顶端种植是真实的 **5.6×**,线性种植 **14.6×**。`#88a` `#108b` 的老毛病第三次出现 —— 正对照在远高于待测量级处开火,所以它们的可分性**不为真实数据的读数背书** |
| 120d | **gate 之间有顺序,而我一直平铺着写 —— 第三十七个** | 我对一个 **0.3×** 的量跑了 `negative_control`("零能不能解释它")和形状比较("它像哪个形状")。**一个还没证明自己非零的量,没有形状可言,也不需要零来解释。** 库里加了 `require_resolvable_first()`:它不通过时,后续所有比较标成 **MOOT** 而不是 PASS/FAIL —— 因为对一个未分辨的效应说"零很小"是一个**假的通过** |
| 120e | **这一轮唯一带出来的实质信息** | 线性项 **−0.00212,2.0×**,方向为负 —— 沿评分越高越早,与 `#114a` 同向。**形状里只有线性成分是可分辨的,而线性成分对两个机制都相容** |

**心理学上必须停在哪里**:`#119a` 说畸变随时间加深(3.5×),这一条站住。**它是记忆变模糊还是故事被反复讲实,这套数据分不出来** —— 需要约 48 倍样本,或者一个能独立测到"回忆准确度"的外部锚。

**按 §0.2 结账**:这一轮关于世界的信念更新 **0**,关于仪器的 **2**(第三十六、三十七个)。**成本回收,不是产出。** 而 `#111c` 的规则说,同一个问题不追第三轮 —— 机制分离到此为止。

---

## Entry 121, added by `E01·A11·R152` — 条件于平均意外度,会把恋物子群的签名一起条件掉;而我的弥散对照不保边际,所以对比不成立

`#120` 的 NEXT 指向 `#104` 的"勾选数分层内的意外度秩"。重读 `#99` 后改问了一个更锋利的问题 —— **`#99` 测的是 S 的水平,没测它的形状,而恋物在现象学上不是"整体偏好罕见",是"对某一样异常强烈"**:

- **弥散(连续维度)**:稀有度均匀铺在一个人所有勾选上 → 集中度与零一致
- **集中(恋物子群)**:少数人的稀有度堆在极少数选项上 → 集中度出现右尾

**第一次跑:全部 gate 通过,判"更接近集中" —— 而我在相信它之前查了正对照有没有产生我设计的签名。它没有。**

| 条件集 | 零 | **集中种植** | 真实 |
|---|---:|---:|---:|
| **条件于 k 和 S**(第一次用的) | +2.006 | **+1.935(−0.071)** | +1.912(−0.094) |
| **只条件于 k** | +1.899 | **+1.981(+0.082)** | **+2.189(+0.290)** |

| # | 结论 | 判定 |
|---|---|---|
| 121a | **条件于 S 把签名一起条件掉了 —— 第三十八个,而且是一整类** | **一个集中的子群按定义就有高 S。** 在 S 分层内标准化,他们被拿去和同样高 S 的人比,**那个"异常"被定义成了正常**。集中种植的尾部超出因此从 **+0.082 变成 −0.071**,符号翻转。**条件于一个与被测结构共变的量,等于把结构条件掉** |
| 121b | **第一次跑的"更接近集中"作废** | 那个判断建立在一个签名已被破坏的正对照上。**gate 全过而结论无效**,因为我没有检查"正对照是否产生了预测的签名"—— `#113c` 的教训在一个新位置重演 |
| 121c | **只条件于 k 时真实确有尾部超出,7.0×** | 真实 **+0.290**,零的抖动 0.042。**但集中种植只给 +0.082,恰好等于 2× 零抖动 —— 它没有清晰开火**,所以这个超出不能被读成"子群" |
| 121d | **我的弥散对照不保边际,对比因此不成立 —— 第三十九个** | 弥散种植把常见勾选换成稀有的,**这改变列和**。而真实-vs-零的对比是在**精确保边际**下做的。两者不可比:所有弥散档把平均 S 抬高 +0.058~+0.184,而真实相对零是 **−0.0016** —— 因为 `#105c` 已经证明**平均 S 是边际决定量,curveball 下不可能移动**。我拿一个能移动它的种植去匹配一个不能移动它的对比 |
| 121e | **可修,而且修法是唯一正确的那个** | 保边际的弥散种植:**成对交换** —— A 得到一个稀有项、失去一个常见项,B 反向。列和与行和都精确不变,但人层面的异质性被造出来。这是唯一能与 curveball 对比放在同一尺度上的弥散世界 |

**按 §0.2 结账**:关于世界的信念更新 **0**;关于仪器 **2**(第三十八、三十九个)。**成本回收。** 而 `#111c` 的规则:`121e` 是这个问题的**第二次**尝试,如果它也失败,停止,不追第三轮。

**心理学上仍未答的那一句**:*那些人是"整体偏好罕见",还是"对某几样异常强烈"?* —— 真实数据在只条件于勾选数时确有 7.0× 的尾部超出,但我还不能说它是子群而不是重尾的连续维度。

---

## Entry 122, added by `E01·A11·R153` — 保边际的两个种植都太弱,问题停在能力边界;而尾部超出本身是真的

`#121e` 指定的修法做了:**成对交换** —— A 得稀有失常见、B 反向,行和列和都精确不变(每块每次抽样都断言)。两个世界**交换总数相同**(3000/块),只改分配给多少人:集中 = 5% 的人承担,弥散 = 60% 的人承担。

| 臂 | 尾−中 | 相对零 | 平均稀有项 |
|---|---:|---:|---:|
| 零 | +1.899 | 0.000 | 3.100 |
| 第二次独立零 | +1.857 | −0.042 | 3.102 |
| **退化种植(0 次交换)** | +1.899 | **0.000** | 3.100 |
| **集中种植(5% 的人)** | +1.885 | **−0.013** | 3.100 |
| **弥散种植(60% 的人)** | +1.909 | **+0.010** | 3.100 |
| **真实** | **+2.189** | **+0.290(7.0×)** | 3.090 |

| # | 结论 | 判定 |
|---|---|---|
| 122a | **两个种植都在零抖动之内,分不开** | 集中 −0.013、弥散 +0.010,零抖动 **0.042**。**3000 次交换/块远低于 MDE**。这是 `#88a` `#108c` `#120c` 的第四次,但方向相反:前三次种植**太强**,这次**太弱** |
| 122b | **我兑现了自己的预先承诺** | `#121e` 和本轮 docstring 都写明"这是第二次尝试,失败则停,不追第三轮"(`#111c`)。**它失败了,我停。** 预先承诺如果只在方便时兑现,就等于没有 |
| 122c | **但尾部超出本身是真的,而且这是一句关于人的话** | 真实 **+0.290,7.0× 零抖动**,而退化种植精确等于零(证明管线本身不造效应)。零精确保留每人的勾选数和每项的基率 —— **所以在"勾了多少"完全相同的前提下,现实里人们的稀有勾选比随机分配更不平均**:第 95 百分位的人比中位的人多勾罕见项的幅度,超出零允许的范围 |
| 122d | **它加了什么,没加什么** | 与 `#99` 的对称加宽**一致**,并把它推广到"稀有项计数、条件于勾选数"这个新坐标。**没有**加子群 —— "集中的少数"与"重尾的连续维度"在这套数据、这个种植尺度上分不开 |
| 122e | **能力边界的确切形式** | 需要一个能在 3000 次交换/块**之上**产生可分辨签名的种植尺度,而那意味着交换量已接近改变数据本身的程度。**更可能的出路不是更强的种植,是一个外部锚** —— 一个独立测到"这个人对这一样有多强烈"的变量,而这个 release 没有 |

**心理学上停在这里的那句话**:

> 在勾选数完全相同的前提下,**人们的罕见偏好比随机分配更不平均** —— 有些人的勾选明显更堆在罕见项上(7.0×)。**但这套数据分不出他们是"一小群对某几样异常强烈的人",还是"一条重尾的连续光谱"。**

`#91` 的盲区在 A11 里没有被打开。**它被精确地定位了,并且量出了打开它需要什么。**

---

## Entry 123 — 前页重写、一条被自己账本推翻却仍在前页的说法、以及一族用过 2 轮的仪器

三件事,一轮里。

### 一 · 前页按 §0.2 重写(Production)

`#122` 的 NEXT:**前页就是产品**,而它描述的是三个 arc 之前的世界。头部 **8,853 → 4,356 字符**。堆叠的警告块清掉(我在 `d070cb4` 修过一次,又长回来了),改成三段:**这说明了关于人的什么**(每条一句话,按实测强度,附账本收据)· **被撤回的**(留着,因为我曾经报告过)· **这套数据做不到的**(每条都量过边界)。

### 二 · 逐个数字回账本核,抓到一条过期的说法

按 `#80c` 的教训(我曾从记忆描述自己九条前的工作,方向搞反),每个头条数字都 grep 回 `RETRACTIONS.md`。**两处不符,一处实质**:

| 我写的 | 账本说 |
|---|---|
| 66.85 ± 0.19 | **66.852 ± 0.191**(`#75b`)。已改 |
| 色情诱导 rho **+0.2515** | **账本里根本没有这个数。** 它来自 A09 之前的 README 正文。`#26` 实测 **+0.2922** 未控制、**+0.2523** 加全指标,**并且把"其中 85% 是作答风格"降级为 UNVERIFIED** |

**前页一直在携带一个它自己的账本已经复杂化了的说法。** 已改,并把 UNVERIFIED 标记带上前页。

### 三 · 找外部锚时撞见的东西,以及一次当场自我纠正

`#122e` 说打开子群盲区需要**外部锚**。查 inventory:**逐兴趣的测量只有评分和起始年龄,没有强度/频率/实际行为** → `#122e` 的边界**由检查确认,不是由假设确认**。

但撞见两族没怎么用过的列:

| 族 | 列数 | 166 轮里用过 |
|---|---:|---:|
| `RATING_NEG_FIB` | 5 | 2 轮 |
| **`FORCED_CHOICE_MOST`** | **10** | **2 轮** |
| `RATING_BINNED_FIB` | 7 | 3 轮 |

| # | 结论 | 判定 |
|---|---|---|
| 123a | **`RATING_NEG_FIB` 看起来像反向计分 —— 我错了,当场纠正** | 值域 [−8, 0],与 0–5 的平均给分相关 **−0.194 ~ −0.318**(合并 −0.390)。我据此推断"这正是 `#26` 说没有的反向计分仪器"。**错。** 题面是**正向**措辞("I find dirtytalking erotic"),只有**数值刻度**是负的 —— 同一个构念翻了符号,不是反向计分的题。**`#26` 的说法站得住** |
| 123b | **但 `FORCED_CHOICE_MOST` 是真的,而且正是缺的那个仪器** | 强制单选**按构造消除作答水平**:必须且只能选一个,"什么都说是"无处施力。实测:各选项对应人群的**平均给分极差 0.324–0.506**,而给分本身的 sd 是 **0.694** —— 极差不到 sd 的 0.73 倍,**选择与作答水平基本正交** |
| 123c | **它 166 轮里只用过 2 轮** | 10 列,n 从 1,990 到 14,974,选项 9–23 个。**这是账本反复说"本 release 分不开"的那些界线上,唯一一个从没被认真用过的仪器** |

**按 §0.2 结账**:一个立住的产出(重写的前页 + 一条过期说法被抓掉),一个边界由检查确认(`#122e`),一次当场自我纠正(`123a`),和一个被发现的工具(`123b`)。

---

## Entry 124, added by `E01·A06·R101` — 自称"色情诱导了我"的人,在强制单选里选的东西和别人一模一样

`#123b` 找到了 `#26` 说本 release 没有的那个仪器:**强制单选按构造消除作答水平**(必须且只能选一个,"什么都说是"无处施力),实测各选项人群的平均给分极差 0.324–0.506 < 0.73 sd,**而 166 轮里只用过 2 轮**。

问题:**自称被色情诱导的人,他们选的东西和别人不一样吗?**(Yes 11,196 / No 2,334)

| 强制单选题 | n | 选项 | 效应 | 零 |
|---|---:|---:|---:|---:|
| youfeelmost | 4,344 | 23 | +0.0074 | −0.0008 |
| otherfeel1most | 4,200 | 23 | +0.0011 | +0.0077 |
| bondagemost | 2,148 | 19 | −0.0040 | −0.0096 |
| gentleness2most | 2,400 | 9 | −0.0052 | +0.0101 |
| nonconsent1most | 1,078 | 16 | −0.0056 | +0.0042 |
| powerdynamic2most | 1,922 | 18 | +0.0125 | +0.0006 |
| toys2most | 1,778 | 21 | +0.0318 | +0.0144 |
| humiliation2most | 566 | 21 | +0.0114 | +0.0061 |
| mentalalteration2most | 460 | 10 | −0.0101 | +0.0062 |
| **合并(9 题)** | | | **+0.0044 ± 0.0041(1.1×)** | **+0.0043(效应的 99%)** |

种植阶梯 **+0.0061 / +0.0231 / +0.0876**,单调,g=0.25 被测出,余量 +0.0733 → **MDE ≈ g=0.05**。

| # | 结论 | 判定 |
|---|---|---|
| 124a | **在广度固定的前提下,那句自述不预测你选什么** | 效应 **1.1× 自身展布**,零是效应的 **99%**。`require_resolvable_first` 正确地把后续比较全部标成 MOOT(`#120d` 的机制第一次真实生效) |
| 124b | **这是 `#61` 第一次拿到免于作答风格的确认** | `#61` 说那句自述"只追踪这个人整体勾了多少",`#26` 反对说"分不开泛泛同意与泛泛认可情欲"。**强制单选让这个反对失效** —— 在这个仪器上,"什么都说是"没有着力点,而结果是零 |
| 124c | **零的强度可以量化** | 种植 g=0.10 给 +0.0170 = 2× MDE。所以这个零许可到 **约 5% 的诱导标签由选择决定**的水平 —— 不是"没测到",是"如果有 5% 我会看到" |
| 124d | **覆盖率守卫第一次拦下一个真结论** | 第一次跑只有 **1/10** 个题通过纳入标准,`check_coverage` 抛错。诊断:五分层匹配后广度残差 **+0.58~+0.72**(诱导与广度 rho +0.29 而分层太粗)。**如果没有那个守卫,我会拿 1 个题的结果去回答一个 10 题的问题** |
| 124e | **修法同时修掉一个常数阈值** | 改成**卡尺 1:1 匹配**(每个 No 配广度最近的 Yes,卡尺 = 0.25 sd),容差从"|Δ| ≤ 0.5"改成"≤ 0.05 × 广度自身 sd"。覆盖率 1/10 → **9/10**。残差 −0.457(系统性地 Yes < No)—— **方向对零有利,是保守的** |
| 124f | **第四十个:我断言两个不同种子算出的量精确相等** | gate 里写 `g=0 必须精确等于真实臂`,实测 +0.0061 vs +0.0044。原因:两臂用了**不同的随机种子**(seed=5 vs seed=1),掩码不同。**退化种植的正确写法是复用同一个种子**,否则那条断言测的是种子不是设计 |

**心理学上的一句话**:

> **说"色情给了我这个癖好"的人,和不这么说的人,在"这些东西里哪一样最让你兴奋"上做出的选择完全一样。** 那句自述报告的是"我口味广",不是一段关于内容的因果史 —— 而这一次,是用一个"什么都说是"无处施力的仪器测出来的。

---

## Entry 125, added by `E01·A05·R098` — 广度携带内容,但只有一点点;而我给库开的口子当场被自己的回归测试堵回去

`#124` 的 NEXT 两件事:把 `#124f` 落地成代码,再用强制单选打 `#26` 自己那条 UNVERIFIED。

### 一 · 广度是内容还是作答风格

`#26` 把"广度里 9–13% 是作答风格"降级为 UNVERIFIED,理由不是测错而是**无法解释**:作答风格指标与广度相关 +0.385(而广度就是结局)。**强制单选让这个反对整个失效,而且这个设计根本不用那个指标 —— 中介问题因此不存在。**

- 广度 = 纯作答风格 → "什么都说是"在强制单选上无处施力 → **不该**预测你选哪一个
- 广度 = 真实内容 → 口味真广的人被迫只选一个时会选得**不一样**

| | 值 |
|---|---:|
| **效应(10 题合并)** | **+0.00173 ± 0.00031(5.5×)** |
| 逐题为正 | **9/10** |
| 打乱零 | −0.00147 ± 0.00140(**1.1×,与零无法区分**,且**与效应异号**) |
| 种植阶梯 | +0.0017 / +0.0032 / +0.0102,单调 |
| 退化臂 `g=0` | **精确复现真实臂**(+0.001726 = +0.001726) |
| 覆盖率 | **10/10** |

| # | 结论 | 判定 |
|---|---|---|
| 125a | **广度携带内容,不是纯作答风格** | 5.5× 自身展布,**9/10 题同向**,而零与效应**异号**(校正只会放大)。**在一个"什么都说是"无处施力的仪器上,广度仍然预测你选什么** |
| 125b | **但很小,必须按这个量级引用** | 用种植标定:g=0.15 给 +0.0085,真实 +0.0017 → **约 g=0.03**。**广度里约 3% 是可从选择中读出的内容**。这不是"广度是内容",是"广度里有一点内容,而且不是零" |
| 125c | **`#26` 的那条 UNVERIFIED 有了方向** | `#26` 说"是泛泛同意还是泛泛认可情欲事物,分不开"。**现在分开了一半**:至少有一部分**不是**泛泛同意,因为泛泛同意在强制单选上无处施力。剩下的那部分有多大,这个设计不测 |
| 125d | **`#124f` 的守卫第一次生效** | `degenerate_matches_reference` 要求 `g=0` 的退化臂**精确**复现真实臂。复用种子后:**+0.001726 vs +0.001726**,精确相等。上一轮那个 +0.0061 vs +0.0044 的差,确认是种子不是设计 |

### 二 · 我给 `negative_control` 开的口子,当场被回归测试堵回去

零是 **−0.00147 ± 0.00140** —— 距零仅 **1.1×**,**与零无法区分**。但 `negative_control` 用的是 `|零| < 0.5×|效应|`,而效应本身只有 +0.0017,所以门槛是 0.0009,它 FAIL 了。

**`#102a` 教的是"别拿零和常数比,要和效应比"。但当效应很小时,还必须问"零自己是否已经与零无法区分"。两者都要。**

于是我给它加了第二条通路(传 `null_spread`,若 |零| < 2×自身展布也算过)。**然后回放五个历史案例,发现这个口子让 `#102a` 那个 91% 的同号零通过了** —— 而那正是这个库存在的起因。

**修法**:自身展布的豁免**只在零帮不上忙时**成立 —— 与效应**异号**,或已小于效应的一半。回归五例全部正确:

| 案例 | 零/效应 | 同号? | 判定 |
|---|---|---|---|
| A05R15 | 85% | **异号** | **PASS** |
| `#102a` | 91% | 同号 | **FAIL** |
| `#110c` | 99% | 同号 | **FAIL** |
| `#124` | 99% | 同号 | **FAIL** |
| `#100` | 5% | 异号 | **PASS** |

**第四十一个,而且是一个新种:我为了修一个真实的漏判,给守卫开了一个口子,那个口子放走了它最初要抓的那个案例。** 抓住它的不是我,是我先写下来的那五个历史回放。

**心理学上的一句话**:

> **口味广的人,不只是"更爱说是"。** 当被迫在一堆东西里只选一个 —— 一个"什么都说是"完全用不上力的场合 —— 他们选的东西仍然和口味窄的人不一样。**只是这个差别很小:广度里大约 3% 能从选择里读出来。**

---

## Entry 126, added by `E01·A11·R154`+`R21` — 稀有偏好特质搬到另一个仪器上仍然在;而口味广的人恰恰选得更常见

`#125` 的 NEXT 指向跨题变异,设计时看到一个更锋利的版本 —— 它同时回答本项目最发达那条声明的最大弱点:

> `#95`/`#99`/`#100` 立起「稀有选项亲和」这条人格维度(信度 +0.432),**但它全部建立在 0-5 多选的勾选上,而那正是"什么都说是"有施力点的地方**。`#100c` 论证过它不是勾选数,但那是**同一族仪器内部**的论证。

强制单选按构造消除作答水平。搬过去:

| | 效应 | 零 | 逐题同向 |
|---|---:|---:|---:|
| **稀有亲和 S → 选中选项的冷门程度** | **+0.1449 ± 0.0191(7.6×)** | +0.0053(**4%**) | **10/10** |
| **广度 → 选中选项的冷门程度** | **−0.0200 ± 0.0061(3.3×)** | −0.0037(18%) | 9/10 |

两者**互相控制**。种植阶梯 +0.145 / +0.394 / +0.566 单调,退化臂 **精确复现**(+0.144930 = +0.144930)。

| # | 结论 | 判定 |
|---|---|---|
| 126a | **稀有亲和特质拿到跨仪器验证** | 在一个"什么都说是"完全无处施力的仪器上,它**仍然 7.6× 地预测**这个人被迫只选一个时选多冷门的。**`#100` 那条特质第一次走出它自己的仪器族** |
| 126b | **而广度指向相反** | 控制 S 之后,广度预测选**更常见**的(−0.0200,3.3×,9/10 同向)。**"喜欢的东西多"和"喜欢冷门的东西"是两个方向相反的维度** —— `#100c` 说它们不是同一件事,这里给出了方向 |
| 126c | **我跑之前没写下最强的那个混淆 —— 第四十二个,而且是新种** | 之前四十一个都是**写下了但写错了**;这一个是**根本没想到**:强制单选的选项与多选块的选项**重叠 89–100%**(10 题里 7 题),那样 S 与结局共享同一批 item,**是恒等式不是测量**。跑完才查出来 |
| 126d | **修法是留出块,而不是只留那 3 个干净题** | 对每个强制单选题,把与它重叠过半的多选块**从 S 里剔掉**再算。效应 **+0.1822 → +0.1449**(掉 20%),**10 题全部干净**而不是只剩 3 题 |
| 126e | **而重叠本来就不是驱动因素,这是可测的** | 重叠 0% 的 3 题给 **+0.1547**,重叠 ≥50% 的 7 题给 **+0.1408** —— **几乎相同**。所以那 20% 的下降是剔块带来的信息损失,不是伪影被去掉 |

**心理学上的两句话**:

> **有些人系统性地被不常见的东西吸引,而这不是"更爱说是"** —— 把他们放进一个必须二选一、"什么都说是"完全没用的场合,他们仍然挑更冷门的那个。
>
> **而"喜欢的东西多"是另一回事,方向还相反**:控制掉稀有亲和之后,口味广的人被迫只选一个时,选的是**更主流**的那个。

---

## Entry 127, added by `E01·A11·R156` — 两个现存声明在共享-item 守卫下都站住;而审计自己的正对照抓到了审计自己的 bug

`#126c` 是本项目第一个**在设计时被漏掉**的混淆(前四十一个都是写下了但写错了)。守卫 `check_disjoint_items` 已落地并回放验证。`#126` 的 NEXT:回头查同样组合 item 派生量的两个声明。

### 审计一 · `#100` 的信度(+0.432)

**通路**:分半按**不相交的块**分,这一层干净。但残差化用的 `picks` 是**全部块**的勾选总数 —— 两半的 item 都在里面。修法:各半用**自己**的勾选数残差化。

| | 原样(共享勾选数) | 修正(各半自己的) |
|---|---:|---:|
| 真实 | +0.4376 | **+0.4611(23.1×)** |
| 固定边际零 | +0.0044 | +0.0587(效应的 **13%**) |
| 种植特质 | — | **+0.8374**(余量 +0.7387) |

### 审计二 · `#116` 的顺序效应(+0.0159)

**通路**:预测器是完整的 68 列偏好矩阵,而**该对自己的两个评分列就在里面** —— 而 `#114` 已证明评分经由回忆偏差扭曲起始年龄。`#116` 把它们当协变量加了,但它们**仍在预测器里**。

| | 值 |
|---|---:|
| 含该对评分列 | +0.0171 ± 0.0026 |
| **剔除后** | **+0.0170 ± 0.0026(6.5×)** |
| **保留** | **99%** |

| # | 结论 | 判定 |
|---|---|---|
| 127a | **`#116` 干净** | 剔除该对自己的评分列后保留 **99%**,剔除带来的变化是效应的 **1%**。共享 item 通路存在但**不做工** |
| 127b | **`#100` 站住,而且修正后更大** | 各半用自己的勾选数残差化后信度 **+0.4611**(原 +0.4376),**23.1× 自身展布**,种植 +0.8374。**共享协变量的影响符号相反(−0.0236),校正只会放大** |
| 127c | **但修正后的零也变大了,记下来** | +0.0044 → **+0.0587**(效应的 13%)。半特异的勾选数更噪声,残差化更不彻底,留下更多由勾选数驱动的共享方差。**两条路线都让真实效应远在零之上,但修正路线的零不是零** |
| 127d | **审计自己的正对照抓到了审计自己的 bug —— 第四十三个** | 第一次跑时种植给 **−0.0569**(负的)。查代码:`aff[:M.shape[0]]` —— `aff` 按**全局人索引**,而 `M` 的行对应 `RAW[t]['ppl']`,不是 `ALLP` 的前若干个。**种植的亲和被逐块分给了错的人**,人层面信号被打散,信度自然为负。**改成 `aff[idx]` 后种植给 +0.8374** |
| 127e | **这就是为什么审计必须自带正对照** | 如果我只跑了真实臂和零臂,会看到 +0.4611 vs +0.0587 并宣布通过 —— 而那个结论恰好是对的,**但支撑它的仪器当时是坏的**。`P5` 的星号规则在审计层面同样成立 |

**按 §0.2 结账**:两个现存声明经受住了一个新守卫的回溯审查(**Closure**,标明),一个守卫落地并回放验证(**产出**),一个审计脚本的 bug 被它自己的正对照抓到(**成本回收**)。

---

## Entry 128, added by `E01·A14·R173` — 不寻常的东西一开始就在那里;而一遍双向去均值把这个结论的符号弄反了

`#127` 的 NEXT 是给六个守卫写一个 linter。写了,跑了,**然后必须给它降级**:171 个轮次里有 146
个早于 `lib/gates.py` 存在本身,所以 94% 的"缺失率"是剧场。有判别力的只有现存声明那 19 轮,
其中 18 轮被标记 —— **但那些标记大半是假阳性**(`degenerate` 只是在源码里看到 "plant" 字样)。

| # | 结论 | 判定 |
|---|---|---|
| 128a | **linter 是覆盖率仪器,不是缺陷探测器** | 它分不出"守卫缺席"和"等价的检查是手写的"。SAFE SIDE 只在**缺失**方向可读;绿色报成"未标记",永不报成"已保护"。`tools/guard_lint.py` 落地,连同它自己的代理账 |

然后本轮的真正内容:把本项目两个最大的现存发现第一次接上 —— `#75` 的**人群共享发育时间表**
和 `#100` 的**稀有亲和特质**,分别回答"什么时候"和"什么口味",从没人问过它们是不是同一条
轨迹的两端。

### 三个关于"一个人怎么变成他自己"的世界

- **A 共享辐射** — 每个人都从常见走向罕见,速率一样。
- **B 起点分歧** — 高稀有亲和的人一开始就在外围。
- **C 行程分歧** — 高稀有亲和的人是走得更久的人。

`ρ_i` = 扣掉类别固定效应(=人群时间表)与这个人整体早熟之后,他的起始年龄残差与类别
**稀有度**的人内相关。`S_i` 来自**多选题选项**,与起始年龄题目零重叠。

| # | 结论 | 判定 |
|---|---|---|
| 128b | **罕见的兴趣来得更早,不是更晚** | `mean(ρ) = −0.0328`,题内跨人置换零 +0.0018,**8.8× 自身自助展布**。人群层的题目回归几乎为零(斜率 +0.061 年/稀有度单位,r = +0.064,n=31)—— **这是纯人层效应:一个人自己那些别人很少有的兴趣,排在他自己那些大家都有的兴趣前面** |
| 128c | **越是最终口味罕见的人,这个提前量越大** | `corr(ρ,S) = −0.0502`(5.5×);**去掉年龄与类别数后 −0.0244(2.7×,保留 49%)**。可报的数是后者 |
| 128d | **世界 C 死,世界 B 活** | 版图不是从常见的中心慢慢长向罕见的边缘。**不寻常的东西一开始就在那里,常见的东西是后来才补上的** |
| 128e | **`#114` 的回忆偏差不能解释它,而且对主量符号相反** | 实际贡献(移除该通道前后的差):`mean(ρ)` **4%**,`corr(ρ,S)` **24%**。λ 规格曲线 0 → −0.2947 全程平坦(−0.0327 ~ −0.0329)。剥离的正对照:纯伪影臂剥完残余 **0.0%** |
| 128f | **反向工作的混淆,写出来因为它让结论更强** | 一个罕见类别之所以罕见是因为少有人获得它,而**近期才获得**的更可能被报成晚 —— 这条通路把罕见推向**更晚**,与观测方向相反 |

### 第四十四个错 —— 而且它把结论的**符号**弄反了

> **一遍「题目去均值 → 人内去均值」在缺失不平衡的面板上不是幂等的。**

人内去均值会把题目均值重新带回来,而**残余的题目主效应恰好沿稀有度排列**,于是它精确地
伪造出一个"题目属性 × 人"的交互:

| 交替投影次数 | `mean(ρ)` | `corr(ρ,S)` |
|---:|---:|---:|
| 1(我原来的做法) | **+0.0699** | −0.0900 |
| 2 | −0.0168 | −0.0560 |
| 3 | −0.0302 | −0.0511 |
| **收敛(13)** | **−0.0328** | −0.0502 |

| # | 结论 | 判定 |
|---|---|---|
| 128g | **这是 `#105` 边际决定量陷阱的镜像** | `#105` 说残差一阶矩是边际决定的,在保边际零下恒为零。这里是**反面**:去均值没做干净时,残余边际**不是**零,而它伪装成交互。**且置换零看不见这件事** —— 置换保留每题的值分布,所以它也带着同样的残余边际,两臂的残余不抵消,差值 18.2× 看上去无懈可击 |
| 128h | **是一次尺度扫描顺手抓到的,不是我看出来的** | λ 规格曲线在 λ=0 处给 −0.0168 而基线给 +0.0790 —— 同一批人同一个量不剥任何东西,只因为循环里 `demean` 被调用了两次。**若不发表整条规格曲线,这个 bug 不会浮现** |
| 128i | **`artifact_cannot_explain` 我喂错了尺度** | 单位幅度的伪影臂问的是"如果年龄 **100%** 由该通道构成会怎样"(给 −0.0410,是效应的 195%,门 FAIL);要的是"该通道**实际**贡献多少"= 移除它前后的差(−0.0051,24%,门 PASS)。**一个上界不是一个贡献。** 这是 `#119d`("是不是零"是错的问题)高一层的同一类错 |

**按 §0.2 结账**:一句关于人的新事实带着它的标度立住(**产出**),一个守卫库的用法错误被
定位并写清正确喂法(**产出**),一个会把结论符号弄反的分析错误被抓住(**成本回收**),
一个 linter 落地但必须自带降级说明(**产出,标明它不是什么**)。

**NEXT**:`128b` 说不寻常的东西一开始就在那里,但"一开始"在这份数据里最早只到 2 年分箱的
第一格。真正能分离的下一件事是 **`corr(ρ,S)` 的 49% 去了哪里** —— 扣掉年龄与类别数就掉一半,
而这两个都不是心理学量,是**仪器量**(#5 的覆盖度定律)。做法:把 `ρ_i` 按**块数匹配**重算
(A05 的标准做法),看剩下的 −0.0244 是不是全部由覆盖度携带。若是,那么"越罕见越早"是真的,
但"越罕见的人越早"只是"答题多的人报得细"。**这个分叉决定 128c 能不能进 README。**

---

## Entry 129, added by `E01·A14·R174` — 特质链接不是覆盖度;而一次尺度不匹配的残差化,把两轮的 `corr(·,S)` 变成了负的勾选数

`#128` 的 NEXT:`corr(ρ,S)` 扣掉年龄与类别数就掉一半,而这两个都是**仪器量**。按 `#5` 的
覆盖度定律,`128c` 进 README 之前必须做**块数匹配**——在设计上关掉这条通路,不靠回归假设线性。

### 第四十五个错,而且它是我这一轮所有"发现"的共同源头

`S`(稀有亲和特质)对勾选数的残差化写成了

```python
S[ok] = z(S[ok] - np.polyval(np.polyfit(z(PK[ok]), z(S[ok]), 1), z(PK[ok])))   # 错
S[ok] = z(S[ok]) - np.polyval(np.polyfit(z(PK[ok]), z(S[ok]), 1), z(PK[ok]))   # 对
```

**把一个 z 尺度的预测值,从原始尺度的变量里减掉。** 原始 `S` 的 sd 远小于预测值的 sd,
所以返回的几乎就是 **−z(勾选数)** 本身。

| # | 结论 | 判定 |
|---|---|---|
| 129a | **两轮的 `corr(·,S)` 测的不是特质,是负的勾选数** | `corr(所谓残差, PK) = **−0.9654**`。`#128c` 的 −0.0502、`#128` 网格里全部 `corr_S` 列、以及 R02 第一次跑出的"匹配后只剩 13–16%",**全部作废** |
| 129b | **我差点把这个 bug 当成一条关于整个 release 的重大发现** | 它顺带产出 `corr(类别数, S) = −0.8297`,读起来像"稀有亲和特质其实是覆盖度",一条会波及 `#100`/`#104` 的警告。**修正后是 −0.0218** —— 特质与覆盖度基本正交。**一个 bug 制造的假警报,比它制造的假发现更危险,因为假警报会去撤别人的东西** |
| 129c | **破绽只有一条,而且是免费的** | **残差与它所回归掉的协变量的相关,在构造上恒等于 0。** 一行断言。已落地为 `lib/gates.py::check_residualized`(第七个守卫) |
| 129d | **单位不匹配的代数不会抛异常** | 变量还叫 `S`,还被 z 标准化过,数值范围也对,产生的相关还是**显著的**。这是 `#117e`(pandas 访问器撞名)同一族:**代码合法、命名正确、结果有意义 —— 而对象已经被换掉了** |

### 修正后的判定

| 协变量集 | 匹配后 `corr(ρ,S)` | 保留 | 类别数残差 | 配对数 |
|---|---:|---:|---:|---:|
| 未匹配 | −0.0459(4.2×) | | | |
| 类别数 | −0.0389 | 85% | 0.003 sd | 4,463 |
| **类别数+年龄** | **−0.0417(3.1×)** | **91%** | 0.006 sd | 4,401 |
| 类别数+年龄+勾选数 | −0.0461 | 100% | 0.002 sd | 4,153 |

置换零 −0.0072(效应的 17%)· 种植正对照 +0.5878 · 三个规格一致。

| # | 结论 | 判定 |
|---|---|---|
| 129e | **`128c` 站住,而且比第一次报的更强** | 用真的特质,`corr(ρ,S) = −0.0459`;**块数匹配后保留 91%**,三个匹配规格 85–100% 一致。`#128` 里"扣掉年龄与类别数掉一半"是坏 S 的性质 —— 因为那个坏的 S **就是**覆盖度变量。修正后保留 **99%** |
| 129f | **`128b` 完全未受影响** | `mean(ρ) = −0.0328`(8.8×)不涉及 `S`。数字逐位不变 |
| 129g | **`#114` 对特质链接的实际贡献 39%** | 上界(单位幅度伪影臂)−0.1668,实际贡献 −0.0168 vs 效应 −0.0434。λ 规格曲线 0 → −0.2947 全程同号(−0.0459 → −0.0370) |

**按 §0.2 结账**:一条关于人的话立住并拿到了它的设计级控制(**产出**),第七个守卫落地
(**产出**),一个会污染两轮全部特质相关的 bug 被定位(**成本回收**),一条差点被我发出去的
假警报被撤回(**成本回收**)。

**NEXT**:`128b` 与 `129e` 现在说的是同一件事的两面 —— 罕见的兴趣来得更早(人层,8.8×),
而口味越罕见的人提前得越多(3.1×,匹配后)。**但"更早"最早只到 2 年分箱的第一格**,而
`#128f` 指出的反向通路(罕见=近期获得)只是被论证过,没有被测过。做法:把每个人的**最早
一格**单独拿出来,问「一个人报告的**第一个**兴趣,它的稀有度是否预测这个人最终的 S」——
这是同一个假说的**离散、不依赖残差尺度**的版本,而且它对双向去均值的做法完全免疫。
若第一个兴趣的稀有度就已经分化,世界 B 从"轨迹形状"升级为"起点本身"。

---

## Entry 130, added by `E01·A14·R175` — 人群从常见处开始,而每个人都在共享曲线上提前了自己那一份

`#129` 的 NEXT:把「罕见的来得更早」做成**离散版本** —— 一个人最早报告的那批兴趣,是不是
比从他**自己的曲目库**里随机抽的更罕见?零是他自己的库,所以这个检验对去均值的做法、
他喜欢多少东西、他整体早熟不早熟、以及覆盖度**全部免疫**。

### 两个数,反号,而且都是真的

| | Δ | 库内分位 | 对人内置换零 |
|---|---:|---:|---:|
| 原始年龄 | **−0.2345** | **0.3252** | 49.4× / 54.5× |
| 扣掉人群时间表 | **+0.0767** | **0.5314** | 14.5× |

| # | 结论 | 判定 |
|---|---|---|
| 130a | **人从大家都有的东西开始** | 最早那一格在他自己库里按罕见度排落在**第 33 百分位**。种植正对照单调,`#114` 实际贡献 17% |
| 130b | **但扣掉时间表后反过来,与 `#128b` 同向** | Δ = **+0.0767**(第 53 百分位),**14.5×**。这是第三个独立统计量,而且它对本弧踩过的每一个坑都免疫 |
| 130c | **口味越罕见的人提前得越多 —— 第二次,换了统计量** | `corr(Δ, S)` = +0.0510,按类别数匹配后 **+0.0532(5.1×,保留 104%)**,置换零 26%。与 `#129e` 的 −0.0417 同向(符号相反是因为 Δ 与 ρ 的定义方向相反) |

### 我的机制解释被自己的检验在同一轮里否定

写下的是「argmin 由**左尾**决定,而题目均值对左尾是盲的」。测了:

| 题目层 Spearman(稀有度, ·) | ρ | p |
|---|---:|---:|
| 起始年龄 **10 分位**(我的解释) | −0.091 | 0.628 |
| 起始年龄 **均值**(`#128` 用的) | +0.011 | 0.954 |
| 起始年龄 **中位数** | **+0.437** | **0.014** |

| # | 结论 | 判定 |
|---|---|---|
| 130d | **左尾解释死** | −0.091,p=0.63。写下来测,而不是叙述,是它当场死掉的唯一原因 |
| 130e | **真机制是中位数,而均值对它完全盲** | +0.437(p=0.014)vs +0.011(p=0.954),同一批 31 个题目。**同一组数据、同一个关系,换一个汇总统计量就从"完全没有"变成"显著"** |
| 130f | ⚠ **这修正 `#128` 的一句话** | 「人群层的题目回归几乎为零(r=+0.064)—— 这是纯人层效应」是在**均值**上算的。**`#128b` 的数字不受影响**(双向去均值精确移除题目均值),但它的**范围陈述**错了:人群时间表**确实**按常见→罕见排列。README 已改 |

### 第四十六个错 —— 在守卫库里

| # | 结论 | 判定 |
|---|---|---|
| 130g | **`require_resolvable_first` 的 MOOT 是 Gate 级的,会跨族传染** | 本轮一个 Gate 里放了两个**独立**的量:原始年龄的 `corr(Δ,S)`(1.1×,未分辨)与扣掉时间表的同名量(5.1×,分辨得很好)。前者把后者也标成了 MOOT。**一个未分辨的量只应让依赖它的行 MOOT。** 修法:`family=` 参数,MOOT 按族计。`#120d` 建这个门是为了防止"对未分辨的量问形状",不是为了让一个 Gate 里的量互相拖累 |

**按 §0.2 结账**:一句关于人的话在**第三个互不依赖的统计量**上立住(**产出**),一个关于
人群时间表的新事实(中位数 +0.437)立住并修正了一句已发布的范围陈述(**产出**),守卫库的
一个跨族传染缺陷被修(**产出**),我自己的机制解释被当场证伪(**成本回收**)。

**NEXT**:`130e` 是本轮最便宜也最可迁移的东西 —— 同一批 31 个题目,均值说 ρ=+0.011(p=0.95),
中位数说 ρ=+0.437(p=0.014)。**这条链上还有多少个"人群层几乎为零"是在均值上算出来的?**
做法:把本项目所有**题目层**的结论(`#69` 的题目主效应、`#88`/`#90` 的三成分、`#63`/`#75` 的
时间表上界)各自用中位数/分位数重算一遍,只报**汇总统计量的选择改变结论**的那些。
这是一次纯 Closure,但它保护的是本项目最老的一批数字。

---

## Entry 131, added by `E01·A14·R176` — 66.852% 保住了;而"跳过排不出来的对"值 12.9 个百分点

**【CLOSURE,明确标注】** `#130` 的 NEXT。不开新战线,保护本项目最老也最大的数字。

| # | 结论 | 判定 |
|---|---|---|
| 131a | **`#75` 站住,与汇总量的选择无关** | 同一批对上,**均值是最好的汇总量**(66.699%);唯一打平的 `trimmed20`(与均值顺序 Spearman = 0.993)只高 **+0.083 个百分点 = 0.3× 种子展布**。全部分位数汇总都更差(62.5–64.2%) |
| 131b | **跳过排序量并列的对,值 +12.9 个百分点** | `median` 跳过并列 76.233%,同一批对上 **63.298%**。分位数汇总把 31 个类别的 580 个有序对压成并列,于是**只在 18,001 个对上打分而均值在 35,486 个** —— 剩下那一半正是差距最大最容易的。**跳过模型排不出来的对是弃权,弃权调整后的准确率不是准确率**。这是 `#101b same_scale` 在排序任务上的形态,而 `same_scale` 必须是**第一个**门 |
| 131c | ⚠ **降级 `#130f`** | `#130f` 用 Spearman(稀有度, 中位数)=+0.437 把 `#128` 的范围陈述改成「人群时间表**确实**按常见→罕见排列」。**中位数排出的顺序确实与稀有度对齐,但它是一个更差的时间表**(63.30% vs 66.70%)。所以那句话是**输给均值的那个汇总量**的性质;表现最好的时间表(均值)与稀有度无关(+0.011)。README 已改回并加限定 |
| 131d | **`#130a` 的测量不变,机制回到 UNVERIFIED** | Δ = −0.2345(49×,种植对照,人内零)是稳的。但左尾解释在 `#130d` 已死,中位数解释在这里被削弱 —— **两个候选机制都不成立,而效应还在**。按 `P6`,这是 UNVERIFIED,不是 OVERTURNED |
| 131e | **我自己那条"随机顺序在 50%"的门,阈值是选的不是量的** | 第一次跑写的是 `abs(random-50) < 1.5`,而随机臂的种子 sd 是 **4.30**。52.87% 被判 FAIL,而它距 50 只有 0.65× 自身展布 —— **一个被选定的地板把一个正常的对照判死了**。改为按它自己的展布判(frontier §2:"a floor treated as measured when it was chosen") |

**按 §0.2 结账**:本项目最老的一个数字在一次针对性攻击下站住,并且现在知道它站住的**原因**
(均值序与 trimmed 序 Spearman 0.993,而分位数序是真的重排了却重排得更差)(**产出**);
一条值 12.9 个百分点的排序任务陷阱被量出来(**产出**);上一轮的一句范围陈述被降级
(**成本回收**);我自己一个选定阈值的门被修(**成本回收**)。

**NEXT**:`131d` 是现在最值钱的开口 —— Δ = −0.2345 是本弧最大的一个效应(49×),而它的
**两个候选机制都死了**。剩下的解释里最便宜的一个是**审查(censoring)**:一个人只能报告
他**已经获得**的东西,而罕见类别的获得率低,所以"在最早一格里看到它"的条件概率被压低。
做法:按人当前年龄分层重算 Δ —— 若 Δ 随年龄单调收缩,审查就是机制;若不随年龄变,审查死,
而我需要第四个候选。这是同一个问题的第二轮,若它也 UNVERIFIED,按 `#111c` 换方向,不追第三轮。

---

## Entry 132, added by `E01·A14·R177` — 审查解释不掉 Δ;而性版图在青春期结束时就基本定型了

`#131` 的 NEXT。Δ = −0.2345 是本弧最大的效应,两个候选机制已死,测最便宜的第三个:审查。

| # | 结论 | 判定 |
|---|---|---|
| 132a | **审查 UNVERIFIED** | Δ 按当前年龄 5 档:−0.2259 / −0.2274 / −0.2359 / −0.2380 / −0.2414,极差 **0.0156 = 1.6× 单档自助展布**,低于分辨率。每档的人内置换零都在 ±0.007 内 |
| 132b | ⚠ **我预注册的审查方向是错的,而模拟纠正了我** | 人为审查(删掉每人最晚获得的 f 比例类别)让 Δ **更不负**(−0.2345 → −0.2342 → −0.2270 → −0.2110),所以审查预测「年轻 → 更不负」。观测方向**一致**,幅度不够。**如果我没跑这个模拟,我会按 docstring 里写反的方向去读这张表,并宣布审查被排除** |
| 132c | **这是 Δ 机制上连续第二轮 UNVERIFIED —— 按 `#111c` 换方向** | `#130d` 左尾死,`#131c` 中位数时间表被削弱,`#132a` 审查 UNVERIFIED。不追第三轮。**Δ 作为测量继续成立(49×,种植对照,人内零),它的机制记为 OPEN** |

### 而分层仪器的正对照没开火 —— 那本身是这一轮最大的东西

审查的**前提**是"年轻人还没走完",可以直接测:

| 年龄档 | 类别数 | ≤17 岁获得的比例 | 最晚获得的那个 |
|---|---:|---:|---:|
| 14-17 | 12.4 | 91.7% | 17.3 岁 |
| 29-32 | **12.9** | **68.4%** | **22.6 岁** |

| # | 结论 | 判定 |
|---|---|---|
| 132d | **29–32 岁的人,68.4% 的性兴趣是 17 岁前获得的,最晚的那个平均在 22.6 岁** | **人内测量,不依赖任何横断面假设**。前几档的「≤17 比例」被年龄本身机械限制(14–17 岁的人不可能在 18 岁获得东西),所以**只有最年长档可读** |
| 132e | **15 年里曲目库只长 +4.5%(12.4 → 12.9)** | 横断面,队列混淆,标 WEAKER。但两条证据同向:**性版图在青春期结束时就基本定型** |
| 132f | **这给 Ivan 的模型 C 划了一条时间边界** | 模型 C(价值回流并重塑表征)要求版图持续被改写。若 68% 在 17 岁前就位、最晚的一个在 22.6 岁,那么**递归重塑若存在,它的作用窗口主要在青春期内,而不是成年后的持续过程**。这不是对 C 的反驳,是对它的**定域** |

**按 §0.2 结账**:一条关于人的新事实立住并且是人内测的(**产出**),它给三个模型之一划了
一条可检验的时间边界(**产出**),一个候选机制被判 UNVERIFIED 并触发换向(**成本回收**),
我预注册反了的方向被自己的模拟纠正(**成本回收**)。

**NEXT**(换方向,按 `#111c`):`132f` 是本项目第一次能对 Ivan 的**模型 C** 说话。
可测的版本是:**如果版图在 17 岁后基本不再增加新类别,那么 17 岁后发生的是什么?**
release 里有两个 17 岁后仍在变的量 —— 每个类别的**评分**(强度)和**广度**(勾选数)。
做法:在 29–32 岁档内,把「17 岁前获得的类别」与「17 岁后获得的类别」分开,比较它们的
评分分布与它们在坐标空间里的位置。若晚获得的那些**评分更低且更靠近已有兴趣的中心**,
那么成年后的变化是**深化而非扩张**,C 的作用被定域到"重排权重"而不是"新增表征"。

---

## Entry 133, added by `E01·A15·R178` — 成年后加进来的不是旁边那一块,是别处;而三分之一的人版图在 17 岁就关上了

`#132` 的 NEXT,新弧 A15。`#132f` 用「68.4% 在 17 岁前」给模型 C 划了时间边界;这一轮问那条
边界之后进来的东西是**深化**还是**扩张**。零是**人内置换早/晚标签**,所以对"他喜欢多少
东西"、"他整体早熟不早熟"、"他答了几个块"全部免疫。

| 证据线 | 晚 − 早 | 人内置换零 | 倍数 | 判定 |
|---|---:|---:|---:|---|
| 与早期集合的连通性 | −0.0047 | +0.0076 | 6.7× | 干净 |
| **去掉稀有度后的连通性** | **−0.0042** | +0.0076 | **6.5×** | **干净,保留 89%** |
| 评分 | −0.1815 | −0.0136 | 6.3× | ⚠ UNVERIFIED |
| 稀有度 | +0.2284 | +0.0016 | 20.3× | 复述 `#130a` |

| # | 结论 | 判定 |
|---|---|---|
| 133a | **EXPAND,不是 DEEPEN** | 晚获得的类别与这个人早期那些**更不相连**。`#114` 在两条连通性线上都**符号相反**,校正只会放大 |
| 133b | **不是稀有度的影子** | 在**配对层**把 `rar_a+rar_b`、`rar_a·rar_b`、`|rar_a−rar_b|` 回归掉后保留 **89%**,`check_residualized` 断言残差与稀有度和正交 |
| 133c | ⚠ **评分那条线是脏的,而且我在跑之前就写下了它会脏** | `#114` 的实际贡献是效应的 **248%**。预注册里写的是「评分这条证据线天生是脏的,而坐标位置那条不是:`#114` 说的是"多喜欢",不是"喜欢哪个"」。**这是本项目第一次预注册一条证据线会失败,然后它按预注册失败了** |
| 133d | **稀有度那条是复述,不是第三条证据** | 「晚获得的更罕见」与 `#130a`「最早一格更常见」是同一个原始年龄事实的两种测法。标 CONFIRMATORY,不进 README 当独立发现 |
| 133e | **标度写清楚:小,但准** | 连通性基线 +0.1398,效应 −0.0047 —— 比人内随机低约 **3.4%**。6.5–6.7× 自身展布。**不要说成大效应** |
| 133f | **被守卫拦下的 35% 本身是一句话** | `check_coverage` 因跳过 797/2,259 而报错。拆开:**17 岁后一个新兴趣都没有的 753 人 = 该档的 33.3%**,早期不足 3 个的 44 人 = 1.9%。**三分之一的人版图在 17 岁就关上了;剩下三分之二的人,它是向外开的** |
| 133g | **对模型 C 的含义** | C 的作用**不是**被定域为"重排既有权重"。成年后发生变化时,加进来的表征与既有的**更不相连** —— 这与"扩张"相容,与"深化"不相容 |

**残余缺口(本轮结构上测不了)**:还有一条**反向回忆**通路没排除 ——「这个东西跟我其余的
不搭,所以它一定是后来才有的」,即连通性低**导致**被报成晚。`#114` 的通道(越爱越早)已因
符号相反被排除,但这一条是**另一条通道**,横断面自报数据分不开。

**按 §0.2 结账**:一句关于人的新事实立住并带着它的标度(**产出**),一个被守卫拦下的跳过量
变成了第二句关于人的话(**产出**),模型 C 拿到一条方向性约束(**产出**),一条预注册会脏的
证据线按预注册脏掉(**成本回收,而且是廉价的那种 —— 它在跑之前就被标记了**)。

**NEXT**:`133` 的残余缺口是唯一还能便宜地打一下的地方。**反向回忆通路有一个可检验的印记:
若"不搭 → 记成晚"是真的,那么这个偏差应该随回忆本身的衰减而加深** —— 而 `#119` 已经独立
测到"记忆畸变随年龄加深,十五年里几乎翻倍"。做法:把连通性对比在**全部五个年龄档**上重算。
若 −0.0042 随年龄单调变大,反向回忆通路活;若各档平坦,它死,而 `133a` 升级。
`#119` 的年龄趋势(3.5×)给了这个检验一个现成的正对照。

---

## Entry 134, added by `E01·A15·R179` — 反向回忆的印记不存在;但那个效应绑在绝对年龄上,而它的量级要降级

`#133` 的 NEXT。`#133` 自己写下的残余缺口 ——「这个东西跟我其余的不搭,所以它一定是后来才
有的」,即连通性低**导致**被报成晚。这条通路是**回忆过程**的性质,所以印记是年龄梯度,
而 `#119` 的「畸变随年龄加深」给了现成的正对照。

| 年龄档 | n | 连通性差 | 倍数 | `#114` 斜率(正对照) |
|---|---:|---:|---:|---:|
| 14-17 | 418 | **−0.0088** | 5.8× | −0.1391 |
| 18-20 | 525 | −0.0072 | 6.1× | −0.1581 |
| 21-24 | 1,231 | −0.0045 | 5.3× | −0.1815 |
| 25-28 | 1,434 | −0.0057 | 7.7× | −0.1772 |
| 29-32 | 1,632 | −0.0069 | 10.9× | **−0.2920** |

| # | 结论 | 判定 |
|---|---|---|
| 134a | **效应在全部五个年龄档上都可分辨、同号** | 5.3–10.9×。**连 14–17 档都有** —— 所以它不是"成年后"独有的 |
| 134b | **反向回忆预测的单调年龄增长 ABSENT** | 趋势 **+0.000133/岁**,而**最年轻的一档效应最大**。年龄仪器同时在 `#114` 的斜率上从 −0.1391 走到 −0.2920(与 `#119` 独立测得的 3.5× 同向),**所以这不是没检出力的零**。`133a` 的**方向**升级:它是关于内容的,不是关于讲述的 |
| 134c | ⚠ **但 `133a` 的量级要降级** | 切点扫描:13.5 → −0.0165、15.5 → −0.0110、17.5 → −0.0069。**方向在每个绝对年龄切点上都一致,量级差 2.4 倍。可报的是方向,不是量级**(frontier §2:效应 X 只许可 ≤X 的断言) |
| 134d | **这个现象绑在绝对年龄上,不是"你自己序列里靠后的"** | **人内中位数分割只有 −0.0006(1.0×,不可分辨)**。中位数分割把绝对年龄打散,效应就消失。这是本弧第一次能说出这个区分 |
| 134e | ⚠ **`#133` 的"晚"实际是 19 岁以后,不是 18 岁以后** | 起始年龄分箱是 17.5 → 22 → 28,所以 `V>17.5`、`V>19.5`、`V>21.5` 选出**完全相同**的集合(n=1,632、晚 5.2 个、效应 −0.0069,三者逐位相同)。**一个我以为在扫参数的循环,其实有三格是同一格** |
| 134f | ⚠ **我预测的剂量方向又反了(本会话第二次,见 `#132b`)** | 写的是"切点越晚效应越强",实测切点越**早**效应越大。两次都是**模拟/扫描纠正了我的预注册方向**,而不是我读表读对了 |
| 134g | **我第一版把置换基线当"零"来断言 —— 问错了问题** | 它是两个随机半边之间的平均连通度,天然为正(+0.0047…+0.0075),零假设**不预测它为零**。所以它是 `offset_control` 的偏移量,不是 `negative_control` 的零。**"这个零应该是零吗"必须在选门之前问,而我选完门才问** |

**按 §0.2 结账**:`133a` 的方向在一个带正对照的年龄检验下升级(**产出**),现象被定位到
**绝对年龄**而不是相对顺序(**产出**),`133a` 的量级被降级为"只报方向"(**成本回收**),
一个假的参数扫描(三格同一格)被发现(**成本回收**),一次问错的零假设被纠正(**成本回收**)。

**NEXT**:`134d` 是新的开口,而且它很具体 —— 效应绑在**绝对年龄**上,那么**边界在哪里**?
切点扫描目前只有三个真正不同的格(13.5 / 15.5 / 17.5),因为分箱只有那么细。但**分箱本身
可以被绕过**:不用切点,直接对每个类别按它的**起始年龄分箱**求它与这个人其余类别的平均
连通度,得到一条「连通度 ~ 获得年龄」的曲线(每人一条,再平均)。若曲线在某个年龄有拐点,
那就是边界;若它是一条直线,那么"版图在 17 岁关上"这个说法要改成"连通度从一开始就单调
下降,没有边界"。**这个分叉直接决定 `#132`/`#133` 该不该继续用"青春期"这个词。**

---

## Entry 135, added by `E01·A15·R180` — 没有边界,但有两族:早来的东西是散的,晚来的东西是一整套

`#134` 的 NEXT。绕开被分箱卡死的切点扫描,直接画「连通度 ~ 获得年龄」的曲线。

| # | 结论 | 判定 |
|---|---|---|
| 135a | **没有拐点** | 增益 0.703 vs 置换零的 95 分位 0.754(**1.6×**)。两个对照都过:种植真拐点被检出(0.805),**种植纯直线不触发**(0.385)—— 所以这不是一个不能失败的检查。**曲线连续单调上升并在 15.5 附近穿零** |
| 135b | ⚠ **所以"版图在 17 岁关上"对连通度不成立** | `#132`/`#133` 的措辞要限定:那句话对**类别数**成立(`#132d`:68.4% 在 17 岁前),对**连通度**不成立 —— 连通度没有边界,它是一条连续的斜坡 |
| 135c | **曲线方向与 `#133a` 相反,而这解开了一个块结构** | `#133a`:晚→早更低;曲线:晚的类别连通度更高。同时为真 ⟺ **晚的彼此抱团**。测:早×早 **−0.01272(27.9×)** · 晚×晚 **+0.02103(23.1×)** · 晚×早 −0.00484(15.2×) |
| 135d | **早来的东西是散的,晚来的东西是一整套** | 最常晚到(19+):怀孕 32.1% · 感官 31.0% · 精神改变 28.1% · 施虐受虐 26.6% · 束缚 26.3% · 权力动态 26.0%。最常早到:外观 8.6% · 开始看色情 10.0% · 衣物 10.9% · 身体部位 11.1%。**题目层:最晚到的 8 个彼此 +0.06735,最早到的 8 个彼此 −0.00751,两组之间 −0.01902** |
| 135e | **这给了 `#75` 一个它没有的维度** | `#75` 只说了先后(外观 14.0 → 精神改变 17.0)。**新的是"成套"**:晚的那批彼此高度连通,早的那批彼此几乎不连通。**发育不是同一种东西按顺序到达,是两种不同组织度的东西先后到达** |
| 135f | **`#133a` 的表述要改** | 「晚获得的与早期那些更不相连」是真的但**是三个块里最小的那个**(−0.00484 vs 晚×晚 +0.02103)。**主效应是晚×晚的抱团,不是晚×早的断开。** README 已改 |

**本轮结构上做不到的**:「关系族**作为一个包**到达」需要更强的检验 —— 在一个人内部,
关系族的条目是不是比同样大小的随机子集在**时间上**更靠拢。本轮只证明了人的晚期集合
**富集**于这一族,没证明它们在**同一时间**到达。

**按 §0.2 结账**:一句关于人的新事实立住,而且它是本弧最锐利的一句(**产出**);
`#75` 拿到一个它原本没有的维度(**产出**);两个候选形状里"有边界"被排除,而排除它的
检验带着"纯直线不触发"这个对照(**产出**);`#133a` 的主次被纠正(**成本回收**)。

**NEXT**:`135` 自己写下的缺口就是下一个分离器,而且它是本项目第一次能问"**打包**"这件事。
做法:对每个人,取关系族(权力动态·束缚·施虐受虐·精神改变·感官)里他拥有的那些,量它们
获得年龄的**离散度**,对比同样大小的随机子集(人内置换,精确匹配)。
    PACKAGE  离散度显著更小 -> 它们确实是一起到的,「一整套」是字面意义的
    ORDERED  离散度与随机无异 -> 它们只是都靠后,而不是同时;「一整套」只能说结构上成套,
             不能说时间上成包
**这个分叉决定 `135d` 那句话能不能保留"一整套"这个词。**

---

## Entry 136, added by `E01·A15·R181` — 「一整套」是字面意义的:关系族一起到,具体族一个一个到

`#135` 的 NEXT。`#135d` 证明的是人的晚期集合**富集**于关系族,没证明它们**同时**到达。

设计上两件事先解决:族用**共现的谱分割**定义(与获得年龄不相交的仪器),离散度在
**题目去均值后的残差**上量(否则测到的是 `#75` 的时间表)。

| 族(共现谱分割) | n | 族内离散度 | 人内置换零 | 差 | 倍数 |
|---|---:|---:|---:|---:|---:|
| **A**(21 个:束缚·权力动态·施虐受虐·非自愿·性玩具·羞辱) | 9,904 | **2.2608** | 2.4954 | **−0.2347** | **17.7×** |
| **B**(10 个:身体部位·温柔·外观·衣物·vore·异常身体) | 8,921 | 2.4511 | 2.2488 | **+0.2023** | 12.1× |

| # | 结论 | 判定 |
|---|---|---|
| 136a | **关系性的那一族是一起到的;具体的那一族是一个一个到的** | 族 A 比同样多个随机挑的**靠拢 9.4%**(17.7×);族 B 比随机**分散 9.0%**(12.1×)。`#135d` 的"一整套"这个词**保留**,而且现在是字面意义的 |
| 136b | **「只是平均更晚」这个混淆在代数上不可能** | **sd 对平移不变**。数值验证:整族平移 5 年,离散度 2.2607600611 → 2.2607600611,差 **0.00e+00**。与救了 `#116`/`#128b` 的是同一种结构免疫 |
| 136c | **分箱地板要报出来** | 起始年龄的分箱是 **2.0 年**,这是离散度的地板。比这更紧的打包本设计看不见 |
| 136d | ⚠ **我在这一轮里写坏了两次对照,而正确答案不是对照** | 第一版「族内置换时间」——**sd 对置换不变**,逐位给同一个数(2.2608 vs 2.2608),一个不可能失败的检查(`#96a`)。第二版改成从这个人其余类别取捐赠残差 —— 撞上 `same_scale`(捐赠池 ≠ 零的池),给出无法解释的中间值 2.3160 |
| 136e | **一条可迁移的规则** | **当一个混淆对你的统计量在代数上不可能时,该做的是陈述并数值验证那个不变性,而不是造一个对照** —— 因为你造的任何对照都会在别的地方与零不同,于是你用一个新的混淆去查一个不存在的混淆。第一版和第二版都是这个错的两种形态:一个不变量,一个池不匹配 |

**按 §0.2 结账**:本弧最锐利的一句话拿到它的时间证据并保住了"一整套"这个词(**产出**);
一条关于"何时用不变性代替对照"的规则被写下来(**产出**);两个写坏的对照被抓住,而抓住
它们的是**逐位相同的输出**,不是我(**成本回收**)。

**NEXT**:A15 的决定已经安全,换弧。`136a` + `#135d` 合起来给出一个**新的**、更基础的开口:
两族的**组织度**不同(题目层彼此 +0.0674 vs −0.0075),而且它们的**到达方式**也不同
(一起 vs 分散)。这正好是 Ivan 三个模型的分界线 —— **模型 A(专用性内容系统)预测的是
一个统一的性内容检测器,而不是两个组织度不同的族**。可测的版本:两族对**同一个人身上的
非性变量**(#101/#102 测过的性别、五因素、成长环境)的响应是否不同?若族 A 有外部锚而
族 B 没有(或反之),那么"性"在这份数据里**不是一个东西**,而 `#69`「这是性内容在问卷里
问不出来」就有了一个建设性的补充:问不出**一个**,但也许能问出**两个**。

---

## Entry 137, added by `E01·A16·R182` — 两族不是两样东西:早/晚这条线在外部锚上是六个里最弱的

`#136` 的 NEXT,新弧 A16。`#135`/`#136` 的两族在**组织度**与**到达方式**上都不同 ——
那它们在**外部成因**上也不同吗?这压在模型 A(专用性内容系统预测一个统一检测器)上。

| # | 结论 | 判定 |
|---|---|---|
| 137a | **差是真的:开放性 +0.0593(4.1×)、无力感 −0.0441(2.8×),都过 11 变量的多重性门槛(2.68×)** | 族分数是人内中心化的评分,所以对广度与默许免疫;族由**共现**定义,与评分和非性变量都不相交 |
| 137b | **`#101` 找到的那个唯一外部锚 —— 性别 —— 对两族一视同仁** | +0.0706 vs +0.0816,**0.7×**。它锚住的是"性",不是其中某一族 |
| 137c | ⚠ **但差不是这两族特有的,而这才是答案** | 信度匹配的零(同一矩阵的 PC2–PC6,同样连贯的分割):最大外部锚差 0.0702 / 0.0869 / 0.1266 / 0.1459 / **0.2106**,而 PC1 是 **0.0593**。**早/晚这条线是六个同样连贯的分割里最弱的一个** |
| 137d | **UNIFIED。模型 A 的统一检测器没有被这条路线反驳** | 两族在组织度与到达方式上不同(`#135`/`#136` 不受影响),但在"各自挂到这个人身上什么非性属性"上,早/晚这条线什么也没切开 |
| 137e | **第一条零在算术上不可能赢,而且它与真实臂信度不匹配** | 随机等大小重分的阈值 **0.1072** > 本数据里存在过的最大外部相关(性别 +0.093)。**判别量的上界被一个本身就在噪声地板上的量卡住。** 而且随机分割的两个分数不连贯 → 噪声更大 → 零更宽(`same_scale`)。**换成正交特征向量的分割才是信度匹配的零** |
| 137f | ⚠ **我把一个门槛喂给了一个要展布的门** | `require_resolvable_first(z, alpha_z)` 会要求 `z > 2*alpha_z` —— 把多重性门槛又乘了 2。多重性门槛是**门槛**不是**展布**。与 `#129i`(喂了单位幅度伪影而不是实际贡献)同一族:**门要什么种类的量,必须先问** |
| 137g | **PC4 是一个具体的路标** | PC4 的分割给出开放性差 **+0.0893**、无力感差 **−0.0900**,都比 PC1 大。**若这张性版图上真有一条被心理学锚住的分界线,它不是早/晚那条** |

**按 §0.2 结账**:一条关于"性别锚住的是性本身而不是其中一族"的事实立住(**产出**);
一个诱人的分割被它自己的信度匹配零否掉,而否掉它的方法(用同一矩阵的正交方向做零)
是可复用的(**产出**);一条算术上不可能赢的路线被识别并命名(**产出**);一次喂错门的
参数被抓住(**成本回收**)。

**NEXT**:`137g` 是本轮最值钱的东西,而且它把问题从"我挑的这条线对不对"变成
"**哪条线是对的**"。做法:对 SIMR 的前 6 个特征向量各自的分割,系统地跑同一套外部锚检验,
把**每个分割 × 每个非性变量**的差做成一张完整的格子(6 × 11),对整格做多重性,
然后看**最强的那条线切开的是什么内容**。若 PC4 的两侧有一个可命名的心理学对立
(例如"施加 vs 承受"、"身体 vs 情境"),那么这张性版图的**主分界线**就第一次有了名字 ——
而它与获得时间无关,这本身就是对 `#75` 时间表叙事的一个限定。

---

## Entry 138, added by `E01·A16·R183` — 版图的主分界线不是"先来后到",是"物件 vs 叙事"

`#137` 的 NEXT。6 个特征向量分割 × 11 个非性变量的整格,多重性由**最大统计量零**
(把非性变量在人之间打乱,重算整格,取 |最大差|,200 次 → 95 分位 = 0.0560)一次性控制。

| # | 结论 | 判定 |
|---|---|---|
| 138a | **PC4 是那条线,而且它压倒性地强** | PC4 点亮 11 个变量里的 6 个;最强一格 PC4×性别 **+0.2529 = 9.0× 全族阈值**。合成正对照登顶整格。PC1(早/晚)只点亮 1 个(开放性 +0.0675) |
| 138b | **它切开的是"物件与装扮" ↔ "情境与叙事"** | A 侧:衣物 −0.444 · 性别扮演 −0.308 · 体液 −0.269 · 特定角色 −0.235 · 脏污 −0.227 · 变形 −0.196 · 性玩具 −0.144 · 羞辱 −0.114。B 侧:年龄相关 +0.321 · 惊悚 +0.245 · 非自愿 +0.222 · 兽/生物 +0.213 · 乱伦 +0.201 · 温柔 +0.187 · 残暴 +0.157 · 怀孕 +0.137 · 神话生物 +0.134。**命名是我对载荷的读法(D5),不是测量** |
| 138c | **五个锚在性别内部仍然成立** | 去掉性别后:开放性 **+0.0831**(5.0×)· 外向性 **+0.0883**(4.9×)· 神经质 **−0.0639**(3.7×)· 无力感 **−0.0590**(3.2×)· 年龄 +0.0327(2.0×)。**两个性别内的方向一致**,不是 Simpson 反转 |
| 138d | ⚠ **唯一一个纯粹是性别影子的,恰好是那个敏感变量** | 「成年后性侵受害」合并 −0.0713,**去掉性别后 −0.0148(0.8×)**。若不做这个控制,它会被当成一条关于创伤与性偏好的结论报出去 |
| 138e | **标度:小** | 去掉性别后的差在 **0.033–0.088**。可报的是**方向与排序**,不是量级 |
| 138f | **对 Ivan 三个模型的含义** | 版图最强的、被心理学锚住的分界线,是一个**表征格式**的区分(可操作的物件 vs 叙事情境),不是"是不是性内容"的区分。模型 A(专用性内容检测器)不预测这个;模型 B(对普通表征的情欲估值 `w(c,t)^T h`)**自然**预测它 —— 读出权重按被估值的**表征种类**分化。这是本项目第一次有一条证据在 A/B 之间给出方向 |
| 138g | **`#135`/`#136` 不受影响,但要重新定位** | 早/晚两族在**组织度**与**到达方式**上是真的(`#136a` 17.7×/12.1×),但它们**不是版图的主轴**。两条线是不同的东西,而 PC4 强一个数量级 |

**按 §0.2 结账**:这张性版图的主分界线第一次有了名字,而且带着它的外部锚与性别控制
(**产出**);一条会被误报成"创伤 → 性偏好"的相关被性别控制杀掉(**成本回收,而且是最该
花的那一笔**);`#135`/`#136` 被重新定位而不是撤回(**产出**)。

**NEXT**:`138f` 是本项目第一次有证据在模型 A 与 B 之间给出方向,而它现在只靠**一条**线的
内容读法(D5)。做法:把 PC4 的两侧当作两个**独立的估值读出**来检验它的可分离性 ——
若模型 B 对,两侧应当有**不同的内部结构**(例如两侧各自的分半信度、各自的维度数、
各自与广度的关系都应当不同);若模型 A 对,两侧只是同一个检测器输出的两个内容簇,
内部结构应当同构。**这是一个关于"结构"的检验,不再依赖我对载荷的命名。**

---

## Entry 139, added by `E01·A16·R184` — `#138f` 降级;真正结构不对称的是早/晚那条线

`#138` 的 NEXT。把命名拿掉,只看结构:四个量(分半信度 · 有效维度 · 第一因子占比 ·
与广度的相关),每个分割两侧各一份,**下采样到 k=9** 以拉平题目个数,评分用成对删除。

| # | 结论 | 判定 |
|---|---|---|
| 139a | ⚠ **`#138f` 降级** | PC4 两侧的结构不对称 **3.62,超出其余中位数 0.8×**,六个里排第 3。**不依赖命名的结构检验不支持模型 B**。外部锚的图样与两个模型都相容 —— 模型 A 的一个检测器,其输出簇完全可以在"哪些人得分高"上不同而内部结构相同。`#138a–e` 不受影响 |
| 139b | **真正结构不对称的是早/晚那条线** | PC1 = **5.73,超出其余中位数 3.0×**,六个里最大。拆开:晚族有效维度**多 2.42 个**,分半信度**低 0.35**,单因子占比**低 0.14**,与广度相关**高 0.19** |
| 139c | **合起来是一句关于人的话** | **关系族在时间上是一起到的(`#136a`,17.7×),但在偏好结构上不是一个东西 —— 它同时到达,却装着好几个独立的维度。具体族相反:分散地到达,却在评分上更像单一的一个东西** |
| 139d | **两条线各占一头,而且都不是"主轴"** | PC1(早/晚):外部锚最弱,结构不对称最大(3.0×)。PC4(物件/叙事):外部锚最强(9.0×),结构不对称 0.8×。**早/晚切开的是结构,物件/叙事切开的是人** |
| 139e | ⚠ **`hash()` 对 str 每进程加盐 —— 这一轮一开始不可复现** | 同一个脚本两次运行给出**不同排名**(PC4 一次第 3 一次第 2),因为种子是 `hash(k+nm)%10000`。改成 `zlib.crc32` 后两次运行逐位相同。**一个跨进程不可复现的轮次不是一个轮次**,而这个 bug 不会报错、不会警告,只会让排名轻轻晃动 |
| 139f | ⚠ **我写了一条不能失败的断言,而且它是假的** | 原文 `asserted(True, "…没有任何一条线的两侧在结构上可分辨地更不同")` —— 同一次运行里 PC1 就是 3.0×。**一个 `asserted(True, ...)` 把一句话变成装饰**,而装饰不会被自己的数据反驳。改成真检验后它当场开火,并直接产出了 `139b` |
| 139g | **范围** | PC1 的 B 侧恰好 9 个带评分类别,在 k=9 下**没有被下采样**;A 侧是 17 抽 9。不构成偏倚,但记下来 |

**按 §0.2 结账**:上一轮最兴奋的那条声明被它自己的后续检验降级(**成本回收,而且这笔花得
最值 —— 那个检验是我在 `#138` 的 NEXT 里为它专门设计的**);一条**新的**、可分辨的结构事实
立住并与 `#136a` 合成一句完整的话(**产出**);两个方法学缺陷被抓住,其中一个会让任何轮次
静默地不可复现(**产出:它是所有轮次共享的**)。

**NEXT**:`139c` 说关系族"同时到达却装着好几个独立维度"—— 那**是哪几个维度**?
`#75`/`A02` 已经命名过三条(谁服从 · 谁被看 · 谁接受,去衰减互相关 ≤ 0.362,有效维度 2.95/3),
而那是在**全部**类别上做的。做法:只在关系族内部重跑同一套命名(块分半 + 去衰减互相关),
看那 2.4 个额外维度是不是就是 A02 的三条 —— 若是,`139c` 就从"好几个"变成一个**具体的**
结构;若不是,关系族里有 A02 没找到的东西,而那是本项目下一个真正的开口。

---

## Entry 140, added by `E01·A16·R185` — 关系族装的就是 A02 那三条轴;而那个"零"是 0.349,不是 0

`#139` 的 NEXT。`#139c` 说关系族"同时到达却装着好几个独立维度"—— 是哪几个?
`A02` 已在**全部**类别上命名过三条(谁服从 · 谁被看 · 谁接受)。

| # | 结论 | 判定 |
|---|---|---|
| 140a | **关系族装的就是 `A02` 那三条轴** | 族内前 4 个成分与 A02 式坐标的**留出** \|相关\|:**0.813 / 0.813 / 0.910** / 0.323。正对照(成分 vs 它自己在另一半人上的重估)= 0.912 |
| 140b | ⚠ **而这个"零"不应该是零,它是 0.349,这改变了结论** | 随机置换 A02 坐标的载荷,留出 \|相关\| 仍是 **0.349** —— 所有成分都被一个一般因子主导。**按"零=0"读,成分4 的 0.323 会被读成"A02 没找到的新维度";按实测地板读,它在地板上,什么也不是。任何"这个成分与那个不同"的说法必须先越过 0.349** |
| 140c | **合起来:关系族是三条轴的交汇处** | 它"同时到达却不是一个东西",是因为它同时装着谁服从 · 谁被看 · 谁接受。而具体族在同一个 k 上更接近单一的一个东西 |
| 140d | ⚠ **有效维度必须在同一个 k 上比** | 未拉平时 12.35(k=17)vs 5.30(k=9);**参与比随 k 增长(纯噪声时 ≈ k)**,那个比较无效。下采样到 k=9 后 **7.64 vs 5.30**,结论方向不变但量级完全不同。这是 `#101b same_scale`,而我在 `#139` 刚为同一个理由做过下采样,这一轮又忘了 |
| 140e | **不解释在地板上的东西** | 成分4 的载荷(乱伦 −0.399 · 神话生物 −0.347 · 兽 −0.331 …)读起来很像一条"虚构/禁忌"轴 —— **不解释。一个在地板上的成分没有内容**,而这正是最容易被叙事捡走的地方 |
| 140f | **范围** | 关系族的类别**全部**也在全类别拟合里,所以"相同"永远无法与"由共享 item 强制"完全分开。报留出相关与地板,由读者定标(`P6` 安全侧) |

**按 §0.2 结账**:`#139c` 的"好几个维度"变成一个**具体且已命名**的结构(**产出**);一条关于
"随机地板不是零"的可迁移事实被量出来,而它当场改变了本轮的结论(**产出**);一次
`same_scale` 复发被抓住(**成本回收**);一条诱人的"新轴"叙事被地板挡住(**成本回收**)。

**NEXT**:A16 的决定已经安全,弧可以关。回到 `140b` 留下的那个洞 —— **0.349 的一般因子地板,
是这份 release 上所有"成分不同"类断言的共同门槛,而本项目此前从没量过它。**
做法:把这个地板量成一条曲线(随成分数、随题目数、随人数),并回头检查
`A02`「有效维度 2.95/3、去衰减互相关 ≤ 0.362」这条**现存声明** —— 0.362 与 0.349 几乎相等。
若 `A02` 的"三条近乎独立"是在一个从未被量过的地板上宣布的,那么本项目最老的一条结构声明
需要重新定价。**这是 Closure,但它保护的是 README 上的一行。**

---

## Entry 141, added by `E01·A17·R186` — A02 的结论站得住,但它引用错了自己的证据;而我的路标前提是错的

`#140` 的 NEXT。新弧 A17。

| # | 结论 | 判定 |
|---|---|---|
| 141a | ⚠ **先撤回我自己路标的前提** | `#140` 的 NEXT 写「0.362 与 0.349 几乎相等」。**那是两个不同的量**:0.349 是族内成分 vs 随机载荷坐标的**留出 \|相关\|**;0.362 是手工量表之间的 `r/sqrt(alpha·alpha)`。**数字接近是巧合被当成了等价** —— 这是"a label is not a description"的数值版 |
| 141b | **`A02` 的原始相关又小又紧,结论站得住** | P-G **0.080** [0.053,0.106] · P-C **0.112** [0.086,0.140] · G-C **0.085** [0.057,0.117]。三对全部 ≤ 0.112,**95% 上界全部 < 0.15** |
| 141c | ⚠ **但它写死的判定是刀尖** | max 去衰减 \|r\| 的自助区间 = **[0.259, 0.528]**,**跨过 `A02` 源码里写死的 0.4 的比例 = 30.2%**。三次重抽里有一次,同一个脚本会打印 "collapses — the axes may be one construct measured badly"。而那个 0.4 是**选的**,从未对过任何地板(frontier §2 第九条) |
| 141d | **不确定性全部来自一个 2 题量表** | `GAZE` alpha = **0.163**(自助 [0.114, 0.207]),去衰减因子 **4.27×**。`POWER` alpha 0.686、`coord4` 分半 0.337 都不是问题 |
| 141e | **README 重新定价** | 那一行原本引用「去衰减互相关 ≤ 0.362」。**改成引用原始相关 ≤ 0.112(95% 上界 < 0.15)**,并注明去衰减版本不稳。**结论不变,证据换成稳的那个** |
| 141f | **一条可迁移的方法学结论** | **去衰减是为保守而做的,但它把一个稳的量换成了一个不稳的量。** 在信度低的地方,校正因子自身的方差比被校正的量还大 —— 而**报出来的永远是校正后的那个,因为它看起来更严谨**。任何用短量表做的去衰减相关都有这个问题 |
| 141g | ⚠ **我第一版给量表的题目随机翻了符号,两个对照当场都塌了** | 正对照只到 0.567、负对照返回 `None`。原因同一个:**没定向的 alpha 不是信度**。`A02` 用的是有意的符号(`worship` 取 −1)。定向后正对照 **0.872**、负对照 **0.012** |

**按 §0.2 结账**:本项目最老的一条结构声明被重新定价 —— **结论保住,证据换成稳的那个**
(**产出**);一条关于"去衰减在低信度处不稳"的可迁移事实被量出来(**产出**);我自己上一轮
写下的路标前提被撤回(**成本回收**);一次"没定向的 alpha"被抓住(**成本回收**)。

**NEXT**:`141c` 打开的洞比它自己大 —— **本项目还有多少条现存声明,其判定阈值是写死在源码里
的一个选定数字?** `A02/R034` 的 `mx<0.4` 是被 `#141` 逮到的第一个,而 `tools/guard_lint.py`
(`#128a`)已经能扫全部 185 个 `run.py`。做法:给 linter 加一条规则 —— **任何与一个字面
常数比较并据此打印结论字符串的行**,都标出来;然后对**现存声明**那 20 来轮逐个看那个常数
是量出来的还是选出来的。这是 Closure,但 `#141c` 证明了它能翻结论。

---

## Entry 142, added by `E01·A17·R187` — `#118` 宣布的"加强"没有发生,而它的判定压在一个落进噪声带的常数上

`#141` 的 NEXT。`tools/guard_lint.py` 新增 `hardcoded_thresholds()`:找**与字面常数比较、
且该比较驱动结论字符串**的行;与自身展布比较(`2*spread`/`boot`/`sd`/`null`)白名单放行。

| # | 结论 | 判定 |
|---|---|---|
| 142a | **扫描落地:现存声明背后的 21 轮,12 处命中,5 处真的在用字面常数打判定** | `A03/R071` `ratio>1.3` · `A12/R168` `kr>0.5` · **`A12/R170` `r250>=3`** · `A02/R040` 与 `A02/R041` 各一个 `mx<50`(与 `#141` 的 `mx<0.4` 同一族) |
| 142b | ⚠ **`#118` 宣布的"加强"没有发生** | 加倍样本后强度 **3.0928 → 3.0484**,略微**下降**。而源码打印的是「加强成功,按新强度引用」。**那句话是假的** |
| 142c | ⚠ **而那个判定压在门槛上方 1.6%** | 写死 `r250>=3`,实测 **3.0484** |
| 142d | **更要紧的是:一位小数的精度本身就在噪声里** | 20 个不同的自助重抽种子给出 **3.023–3.163(sd 0.035)**。**一个设在 3 的门槛落在这个噪声带内** —— 所以这个预注册判定的输出由重抽种子决定,而不是由数据决定 |
| 142e | **`#118` 的效应本身不受影响** | +0.0339,**3.1×**(> 可分辨门槛 2×),46/68 为正,种植对照单调。README 那一行的**数字不变**,改的是它旁边那句关于"加强"的话 |
| 142f | ⚠ **我在这一轮里被自己的门抓到一次** | 我先判定「README 的 3.1× 是舍入错误,应作 3.0×」。**错的** —— 3.0928 确实舍入到 3.1,门当场标 FAIL。撤回。真正的问题不是舍入,是**一位小数的精度落在噪声之内**,而我差点把一个正确的数字改错 |

**按 §0.2 结账**:一条可复用的扫描规则落地并立刻找到 5 处(**产出**);一句写在源码里、
从未被检查过的"加强成功"被证伪(**成本回收**);一条关于"预注册阈值必须离开自己的噪声带"
的规则被量出来(**产出**);我自己一次错误的修正被门当场挡住(**成本回收 —— 而且它保护的是
一个本来正确的数字**)。

**NEXT**:`142d` 是这条链上最可迁移的东西,而且它可以直接变成一个守卫:
**任何预注册阈值,都必须先证明它离开了被比较量自身的噪声带。** `lib/gates.py` 现有的
`resolvable`/`require_resolvable_first` 检查的是"效应 vs 零",不检查"效应 vs 阈值"。
加一个 `threshold_outside_noise(name, value, threshold, spread)`,并回头把
`142a` 命中的另外 4 处(`A03/R071` 的 1.3 · `A12/R168` 的 0.5 · `A02/R040`/`R17` 的 50)
各量一次距离。**`#118` 已证明这类检查能翻掉一句已写下的结论。**

---

## Entry 143, added by `E01·A17·R188` — `#141` 整条撤回:我审计的是一个 118 条之前就已作废的版本

`#142` 的 NEXT 是给另外四处写死的常数各量距离。量到 `A02/R040` 的 `mx<50` 时看到实测是
**101% / 86%**,而它的预注册判定写着 `mx>=80 → ONE CONSTRUCT`。**顺着读下去,那不是发现,
是我上一轮的错误。**

| # | 结论 | 判定 |
|---|---|---|
| 143a | **`R16`/`R17` 的判定早已被 `R18` 作废,`Entry 24` 已经 settled 了整条线** | `R18` 诊断出 profile-profile 相关在 14 变量电池上零分布 sd ≈ 0.30 —— **仪器在构造上欠功效**,线被冻结。`Entry 24`(`A02·R19`)用验证过的信度阶梯(正对照 SUBSTANCE 对自身 **r_true = +1.018**,sham +0.023)测得 **POWER–SUBSTANCE r_true = +0.605**,共享 **37% 方差,不是 5%**;预注册判定落在中带 → **UNVERIFIED**;**「2.95 of 3」withdrawn** |
| 143b | ⚠ **`#141` 整条撤回** | 它重跑 `A02/R034`、自助它的去衰减值、发现 30.2% 跨过写死的 0.4 —— **这些计算本身没错**,但它审计的是一个**已经作废的版本**,而且它**把已撤回的「有效维度 2.95/3」写回了 README** |
| 143c | ⚠ **`#141e` 的方向是反的** | 它推荐"改引**原始**相关(≤0.112),因为那个稳"。**原始相关正是被衰减压平的那个。** `Entry 24` 的全部要点就是:一把**验证过的**去衰减给出 +0.605,而 `R10` 的 0.362 之所以错,是因为它**用了属于另一个测量的信度**(GCCA 坐标的块分半,去除一个手工复合量表) |
| 143d | **失败模式,写成一句可迁移的话** | **我用重跑那一轮的方式去审计一条声明,而没有读那条声明后来的账本条目。一个轮次自己的输出不是那条声明的当前状态;账本才是。** 这是 `feedback_false_self_retraction` 的镜像 —— **假的自我确认**:我"重新定价"到一个 118 条之前就被撤回的数 |
| 143e | **而真正的洞比 `#141` 大:README 从 `Entry 24` 起就没被改过** | 它带着「去衰减 ≤ 0.362,有效维度 2.95/3」又走了 **118 条**。**一条账本里的撤回,从来没有走到前页。** 这不是某一轮的错,是**两个文件之间没有任何机械连接** |
| 143f | **产出一:`tools/readme_ledger_audit.py`** | README 上每个可辨识数字,是否出现在账本的撤回语境里。127 个数字 → **29 个命中**(大多良性:那个数是撤回行里的**更正值**)。**输出是待读清单,不是判决**(`P6` 安全侧写在 docstring 里) |
| 143g | **产出二:第八个守卫 `threshold_outside_noise`** | `#142` 的 NEXT。预注册阈值必须离开被比较量**自身的噪声带**。回放:`#118` 的 `r250 = 3.0484` vs 写死的 `3`,噪声 0.035 → **1.4×,FAIL**;3.60 的例子 → 17.1×,PASS |

**按 §0.2 结账**:README 上一条挂了 118 条的错误声明被改正(**产出**);两个工具落地,其中
一个是**账本与前页之间此前完全不存在的机械连接**(**产出**);我上一轮一整轮的工作被撤回
(**成本回收**);一条关于"审计要读账本不要重跑轮次"的规则被写下来(**产出**)。

**NEXT**:`143e` 才是真正的开口 —— **29 个命中里,除了 A02 这一条,还有几条是真的?**
`tools/readme_ledger_audit.py` 已经把清单打出来了,但它是**待读清单**,必须人工分诊:
每一个命中都要回答「这个数字在那一行里,是**被撤的**还是**更正后的**」。
先看最像的三个:`0.012`(Entry 23 说 θ 的域封闭已被 `#17` 部分撤回)· `1.05`(Entry 85a
说「题目效应是交互的 11.7×/3.5×」已撤回,更正为 1.0×–2.4×)· `0.2515`(Entry 说
"账本里根本没有这个数")。**`#141` 已经证明:一条挂在前页的作废声明,会让整整一轮的工作
建在它上面。**

---

## Entry 144, added by `E01·A17·R189` — 一个不带出处的数字,是一个无法被撤回的数字

`#143` 的 NEXT:把 `readme_ledger_audit.py` 的 **29 个命中**逐个分诊。

| # | 结论 | 判定 |
|---|---|---|
| 144a | **29 个里只有 1 个是真的** | 其余是巧合子串、或那个数字是撤回行里的**更正值**、或 README 已在"被撤回的"表里正确标注。**一个高命中率的审计工具,其价值取决于分诊有没有真的做** |
| 144b | ⚠ **那一个是 `0.2515`,而账本原文是「账本里根本没有这个数」** | 它来自 A09 之前的 README 正文。`#26` 实测 **+0.2922** 未控制、**+0.2523** 加全指标,并把「其中 85% 是作答风格」降级为 **UNVERIFIED**(全部题项都是情欲内容且无反向计分)。**README 英文正文行 131 一直带着那个虚构的数与那句已被降级的话** |
| 144c | **一条保守方向的陈旧** | `0.385`:README 挂着 `#26` 的 UNVERIFIED 标记,而后来的条目说**强制单选设计让那条反对整个失效**。标记本身过期了,方向是保守的,记录但不改 |
| 144d | **而它暴露的东西比它本身大:README 有两套并行叙述** | 中文表(42 行)与英文正文(116 行),同一条声明各写一遍 —— 色情诱导 82.7% 在行 56 与 131;覆盖度定律 +0.815 在行 87 与 122;三条轴 2.95 在行 54 与 111。**`#143` 找到的是"撤回没走到前页";这一轮找到的是更坏的一种 —— 走到了前页的一半。** §P16「一个事实一个家」的直接违反 |
| 144e | ⚠ **而内部一致性检查没抓到它,原因就是答案** | 新加的 `internal_consistency()` 找到 5 个引用标记在两侧带不同数字,**但行 131 不在其中,因为行 131 不带任何引用标记**。**一个不带出处的数字,是一个无法被撤回的数字** |
| 144f | **审计盲区的大小被量了出来** | `uncited_numbers()`:README 有 **17 行带数字却不带任何出处**;去掉 Zenodo DOI 与 python 版本号,**前页有约 14 个数字没有回到账本的指针**。`#143` 的那个真错误就在这类行里 —— **它不是被漏掉的,它是结构上不可见的** |

**按 §0.2 结账**:一条挂在前页、账本明说"根本没有这个数"的虚构数字被改正(**产出**);
README 的双叙述结构被识别并量化(**产出**);审计的**盲区大小**第一次有了数字
(**产出 —— 一个工具报出自己看不见什么,比它报出看见了什么更有用**);29 个命中里
28 个是良性的这件事被写下来(**成本回收:一个不分诊的高命中率工具会制造下一个 `#141`**)。

**NEXT**:`144f` 给出了一个可以直接关掉的洞 —— **前页那 14 个不带出处的数字,每一个都要么
接上它的账本条目,要么被删。** 这不是审计,是修补,而且它有一个可验证的完成条件:
`uncited_numbers()` 的输出降到只剩 DOI 与版本号。做完之后,`readme_ledger_audit.py`
的三条规则(撤回语境 · 内部一致 · 无出处)就构成一个**前页不变量**,可以在每轮结束时跑一次。
**心理学上的收益是直接的**:那 14 个数字里包括「门控 P(enter|parent>0)=0.99」、
「人档案分半信度 0.727」、「集合比同样大小的基率集中 0.88%」—— 全都是关于人的话,
而它们现在无法被撤回。

---

## Entry 145, added by `E01·A17·R190` — 盲区被我高估了 14 倍;而剩下的那一个,挂在一句账本说是假的话上

`#144` 的 NEXT:把前页那 14 个不带出处的数字逐个接上账本或删掉。

| # | 结论 | 判定 |
|---|---|---|
| 145a | ⚠ **先撤回 `#144f`:那 14 个是粒度伪影** | `uncited_numbers()` 第一版**按行**判,而引用往往在同一**段**的别处。按段重判:**17 行 → 2 段**,其中一段是 `python 3.14`。**盲区被高估了 14 倍。检查的单位必须与"出处"的单位一致** —— 与 `#101b same_scale` 同一族,只是发生在**审计工具**上 |
| 145b | **而那唯一真的一个,比 14 个加起来更值** | 覆盖度定律那一段。数字(`+0.815`)是真的;**跟在它后面的那句话是假的** —— 「我没在任何已发表分析里见过有人提这一条」。`Entry 15` 早已杀掉它:数据发布者自己的 `analysis/swarm/14-missingness.md`,**早五个月**,写明门控列上的分组比较是在比不同的子人群,并附逐列缺失率与 18 列中 17 列 Cohen's d > 0.7 |
| 145c | **而这句假话在两套叙述里各有一份** | 中文 段 18 与英文 行 123,**两处都在,两处都已改**。这与 `#144d` 是同一个结构问题的第二个实例 |
| 145d | **剩下的是真的,而且很大** | `corr(一致性缺口, 覆盖度缺口) = +0.815`(9 次切分);pornhabit **0.2285 → 0.0871 → 0.0439**,**与性别的排序反转**(性别可比测量下 0.0778)。**危害是他们记录的;把它量化并因此翻掉自己一条结论(`#11`),是我做的** |
| 145e | **前页不变量现在是绿的** | 三条规则:撤回语境(29 命中,已分诊)· 内部一致(5 命中,均为英文正文是中文表的子集)· 无出处(**只剩 `python 3.14`**)。`#144` 的 NEXT 给出的可验证完成条件已达成 |
| 145f | **两轮下来的图样,值得单记** | `#144b` 是一个**虚构的数**,`#145b` 是一句**虚构的优先权**。两者都在前页,都被账本杀过,都活了下来。**前页最可引用的那些句子,恰恰是账本杀掉的那些** —— 因为可引用的句子是为效果写的,而效果与准确性不是同一个目标 |

**按 §0.2 结账**:一条方法学结论保住了它的**量化**部分并卸掉了它虚构的**优先权**部分
(**产出**);前页不变量第一次全绿,且完成条件是可验证的(**产出**);我上一轮对盲区的
14 倍高估被撤回(**成本回收**);一条关于"前页最可引用的句子最可能是假的"的观察被写下来
(**产出**)。

**NEXT**:A17 这条弧的决定已经安全 —— 前页与账本之间现在有三条机械连接,且完成条件可验证。
**换方向,回到人。** 本会话开出而未走完的最大一个心理学开口是 `#131d`:
`#130a` 的 Δ = −0.2345 是本弧最大的效应(49×),而它的**三个候选机制全部死了或未验证**
(左尾 `#130d` 死 · 中位数时间表 `#131c` 被削弱 · 审查 `#132a` UNVERIFIED)。
按 `#111c` 我当时换了方向,而现在有一个当时没有的工具:`#138`/`#140` 建立的
**共现谱分割**与**随机载荷地板**。做法:把 Δ 按 PC4 的两侧分开重算 ——
若"物件"侧与"叙事"侧的 Δ 差别很大,那么"最早那一格更常见"就不是一条关于时间的规律,
而是关于**哪种东西先被认成性的**。

---

## Entry 146, added by `E01·A18·R191` — Δ 的量级由"谁够格进入分析"决定,而那就是广度

`#145` 的 NEXT。把 `#130a` 的 Δ(−0.2345,49×)按 `#138` 的内容线切开。

| # | 结论 | 判定 |
|---|---|---|
| 146a | **内容线不解释它** | 六个分割的两侧都差很多(2.7–9.4×),而 **PC4 超出其余中位数只有 0.2×**。「最早那一格更常见」**不是**关于"哪种东西先被认成性的" |
| 146b | **第四个死掉的机制:稀有度离散度** | `corr(|Δ|, sd(该侧稀有度)) = +0.191`,归一化后离散系数纹丝不动(0.55 → 0.55)。前三个:左尾 `#130d` 死 · 中位数时间表 `#131c` 被削弱 · 审查 `#132a` UNVERIFIED |
| 146c | **而真正在动的东西找到了** | `corr(|Δ|, 合格者的平均类别数) = **+0.661**`;`corr(|Δ|, log 合格人数) = −0.570`。**\|Δ\| 在 12 个侧格上从 0.0154 到 0.4168 —— 27 倍**,而合格者的平均类别数从 17.5 到 27.0 |
| 146d | ⚠ **所以 `#130a` 的范围陈述要改** | **Δ = −0.2345 不是一个常数,是一个关于特定纳入规则所选出的人群的数。** 一句关于人的话:**「最早报告的那批兴趣更常见」不是普遍规律 —— 口味越广的人,他最早的那批越集中在大家都有的东西上** |
| 146e | ⚠ **我预注册的种植方向,第三次反了** | 把该侧**最罕见**的搬进最早一格,Δ 应当**更不负**。实测 −0.3632 → −0.2132 → −0.0049。种植是对的,错的是我写的期望(`#132b` 审查 · `#134f` 剂量 · 本轮)。**三次都是模拟/扫描纠正了我,没有一次是我读表读对的** |
| 146f | **一个必须记的缺陷** | `PC1` 的 B 侧只有 **45** 个合格者,人内置换零 **−0.037**(其余 11 格都在 ±0.023 内)。**小 n 的格子零不干净,那一行不可读** |

**按 §0.2 结账**:一条关于人的**新**事实立住(广度决定"最早的更常见"有多强)(**产出**);
`#130a` 的范围陈述被改正,从"一个常数"降为"一个关于选中人群的数"(**成本回收**);
第四个候选机制被证伪(**成本回收**);我第三次预注册反了的方向被自己的模拟纠正(**成本回收**)。

**NEXT**:`146c` 把 Δ 变成一个**关于广度的**量,而广度在本项目里已经被反复证明是**仪器量**
(`#5` 的覆盖度定律 `corr = +0.815`;`#129` 里坏掉的 S 就是 −z(勾选数))。
所以下一个分离器很具体:**把 Δ 在广度上做卡钳 1:1 匹配后重算**(`#129e` 的做法,
那一轮它把一个 −0.0502 的效应保住在 −0.0417)。
    survives  匹配后 Δ 仍显著为负 -> 「最早的更常见」是真的,只是强度随广度变
    dies      匹配后归零 -> 它一直是覆盖度,而本项目最大的效应之一要整条撤回
**这个分叉决定 `#130a` 能不能留在 README 上。**

---

## Entry 147, added by `E01·A18·R192` — `#130a` 得救而且是普遍的;而 `#146` 的广度联系是一条算术界

`#146` 的 NEXT 写的是「在广度上做卡钳 1:1 匹配」。**先纠正我自己的路标:那是错的工具** ——
匹配比较两个组,而 Δ 没有组(它的零是人内的,`E[Δ_i|null]` 恰好为 0)。广度不能让 Δ_i
偏离零,只能改变**幅度**。而幅度有一条纯算术的界:`Δ` 是 m 个的均值减 k 个的均值,
幅度被 **(k−m)/k** 卡住,**所以类别越多 |Δ| 越大是算术强制的**。
判据:用每个人**自己的**置换零展布归一化,那条界精确抵消。

| 广度层 | n | 平均类别数 | Δ | **归一化 z** | 倍数 |
|---|---:|---:|---:|---:|---:|
| 0 | 1,238 | 8.0 | −0.1783 | **−0.4685** | 17.5× |
| 2 | 1,977 | 11.5 | −0.2408 | −0.5790 | 26.9× |
| 4 | 2,090 | 19.2 | −0.2792 | **−0.5809** | 26.9× |

| # | 结论 | 判定 |
|---|---|---|
| 147a | **`#130a` 得救,而且比原来强:效应是普遍的** | **五个广度层每一层都可分辨为负(17.5–27.3×)**。不是口味广的人才有。整体 z = **−0.5515** |
| 147b | ⚠ **`#146c`/`#146d` 撤回** | 真实 z 趋势 **−0.00795/类别**,而正对照 B(一个**按构造与广度无关**的种植效应)产生的残余趋势是 **−0.00881/类别**,**比真实的还大**。**归一化后的残余梯度完全被一个与广度无关的效应解释掉** —— `#146c` 的 +0.661 是那条算术界。`#146d` 的「口味越广的人越明显」**撤回** |
| 147c | **归一化本身经过了正对照** | 正对照 A(人内打乱起始年龄的纯算术世界):整体 z **+0.0063**,趋势 −0.00334 —— 平坦且为零。正对照 B 的整体 z = **−0.7861**,种的效应被检出 |
| 147d | **换单位** | `Δ = −0.2345` 是**尺度依赖**的,不可跨人跨子集比较(`#146` 看到的 27 倍就是这个)。**可比的数是 `z = −0.5515`**,以每人自己的置换零展布为单位。README 已改 |
| 147e | **一条可迁移的话** | **一个"差值"统计量,如果它的两个成分共享分母的一部分,它的幅度就带着一条组合学的界。跨子集比较它之前,必须先用每个单位自己的零展布归一化 —— 否则测到的是集合大小。** `#139`/`#140` 的下采样是同一件事的另一种解法 |

**按 §0.2 结账**:本项目最大的效应之一被证明是**普遍的**,并且第一次有了可比的单位
(**产出**);上一轮那条诱人的心理学解释(「口味广的人更明显」)被它自己的正对照撤回
(**成本回收 —— 而且那个正对照是本轮为它专门造的**);一条关于差值统计量组合学界的规则被
写下来(**产出**)。

**NEXT**:`147a` 现在说的是一句很干净的话 —— **每个人最早报告的那批性兴趣,都比他自己
曲目库里随机抽的更常见,强度约 0.55 个自身展布,与他喜欢多少东西无关。**
而它的**机制**仍然空着(四个候选全死:左尾 · 中位数时间表 · 审查 · 稀有度离散度)。
既然效应与广度无关、与内容线无关、与年龄无关,剩下最有杠杆的一条是**共享的发育顺序本身**:
`#75` 的时间表是按题目均值排的,而 Δ 是一个**人内**量。做法:把每个人的起始年龄先减去
`#75` 的题目均值(`#128` 已经做过这一步并让符号翻转),然后在**归一化的 z** 上重算 Δ ——
若 z 归零,那么「最早的更常见」**就是**共享时间表在人内的投影,机制找到了;
若 z 不变,时间表被排除,而剩下的候选空间就非常小了。

---

## Entry 148, added by `E01·A18·R193` — 机制找到了:那不是关于个人的,是共享时间表投影到每个人身上的样子

`#147` 的 NEXT。跑之前先查:**`#130b` 当年就测过这一步**(扣掉时间表后 Δ 从 −0.2345 翻成
+0.0767)。这一轮把它钉在 `#147` 的可比单位上,并加一个当年没有的对照。

| 臂 | Δ | **z** | 倍数 | 解释掉 |
|---|---:|---:|---:|---:|
| 原始起始年龄 | −0.2303 | **−0.4462** | 45.3× | — |
| 减 `#75` 题目均值(共享时间表) | −0.0495 | **−0.0884** | 9.4× | **80%** |
| 减方差匹配的随机题目向量 | −0.1442 | −0.2851 | 27.4× | 36% |

| # | 结论 | 判定 |
|---|---|---|
| 148a | **机制是共享的发育时间表** | 它解释掉 **80%**,而方差相当的随机噪声只解释 **36%**;时间表特异的部分 **44 个百分点**,`offset_control` 给 **20.0×** 自身展布。**「最早报告的那批性兴趣更常见」不是关于个人的,是那张所有人共享的时间表投影到每个人身上的样子** |
| 148b | **剩下不到 20% 才是个人的** | 扣掉时间表后 z = −0.0884 仍有 **9.4×** —— 小,但不是零 |
| 148c | ⚠ **三个臂原本算的不是同一个统计量,而我自己的随机对照抓到了它** | 第一版"最早一格"= 全部并列者;原始年龄按 2 年分箱,`m` 可到 **24**,而减去任何**连续**题目向量后并列全消、`m` 恒为 1。随机对照当场开火(减随机向量也掉 48%)。全改 **m=1** 后三臂精确可比(`#101b same_scale`) |
| 148d | ⚠ **「这个零应该是零吗」—— 不应该,而我第一版问错了** | 改成 m=1 后随机臂**仍然**掉 36%。那是**稀释**不是伪影:方差相当的噪声同等程度打乱"谁最早"。**所以随机臂是基准不是零** —— `offset_control`,零的种类是"减去方差匹配的随机题目层向量后的同一个 z" |
| 148e | ⚠ **而账本五轮前就写着** | `#130b` 已经指认了时间表。其后 `#132a`(审查)与 `#146b`(稀有度离散度)各花一轮找机制 —— 找的是**同一本账里已经定位过**的东西。与 `#143` 同族:**用重跑的方式去问一个账本已经回答过的问题**。(`#130d`/`#131c` 打的是另一个问题「时间表**为什么**按稀有度排」,那条线仍开着) |

**按 §0.2 结账**:本项目追了五轮的一个机制被定位,并且带着它的稀释基准与百分比
(**产出**);一个 `same_scale` 缺陷被本轮自己的随机对照抓到(**成本回收**);
一次问错的零假设被纠正(**成本回收**);两轮找错方向的机制搜索被识别出来(**成本回收**)。

**NEXT**:`148b` 留下的 20% 是现在最干净的开口 —— **扣掉共享时间表之后,还有 9.4× 的
个人成分,而它是什么完全没查过。** 它有一个现成的候选:`#128b` 测到的
「扣掉时间表后,一个人自己那些罕见的兴趣来得**更早**」(−0.0328,8.8×)。
两者符号一致(都指向"个人偏离时间表的方向是朝罕见")。做法:把 `#148` 的残差 z
与 `#128b` 的人内 ρ 直接相关 —— 若它们是同一个量,那么本项目两条独立线索合并成一条;
若不相关,个人成分里有两件不同的事,而这是新东西。

---

## Entry 149, added by `E01·A19·R194` — 那 20% 是真的,但它不是一个特质

`#148` 的 NEXT。两个候选:`#148b` 的 **z_resid** 与 `#128b` 的 **ρ**,符号一致。

| 臂 | 值 | 倍数 |
|---|---:|---:|
| **同一批类别上** corr(z_resid, ρ) | **−0.5741** | **91.6×** |
| **换到这个人的另一半类别** | **−0.0910** | 7.8× |
| 人内置换地板 | −0.0034 | — |
| 种植正对照 | −0.6695 | 与种植量:z −0.428、ρ +0.559 |

| # | 结论 | 判定 |
|---|---|---|
| 149a | **同一批兴趣上,两个量读的是同一件事** | −0.5741,**91.6× 地板**;而地板本身是 −0.0034(人内置换把相关完全打掉) |
| 149b | **但换到他自己的另一半兴趣上,只剩 −0.0910** | 小 **6.3 倍**,仍 7.8× 可分辨、7.5× 高于地板。**真的,但弱** |
| 149c | ⚠ **我预注册的两个世界都不对** | 写的是 SAME(同一个量)与 TWO(两件事)。**答案是第三种:同一次测量上一致,作为特质不稳。** 而揭示它的那个臂,是我起错名字的那个 |
| 149d | **判定:那 20% 是真的,但不是一个稳定的人格特质** | **一个人偏离共享时间表的方向,更像一次次的具体情况,而不是一个跟着他走的性格** |
| 149e | **而这与 `#100` 形成一个干净的对照** | 稀有亲和特质本身**可靠**(跨不相交块分半 **+0.4611,23.1×**),而"什么时候得到它们"的个人偏离**不可靠**(跨半 −0.0910)。**喜欢什么罕见的东西是一个稳定的性质;什么时候得到它们,不是** |
| 149f | ⚠ **两个我问错的地方** | ① 我预期"两个量共享算术结构 → 地板不为零"——**错了**,人内置换把相关完全打掉。② **我把跨半臂当成"天花板"**,而它比同数据臂**小 6.3 倍** —— 一个天花板不可能低于它所限制的量。**它是另一个量,而我给它起错了名字。给一个臂起错名字,会让你读不出它在说什么** |

**按 §0.2 结账**:那 20% 从"完全没查过"变成"真的但不是特质",并与 `#100` 合成一句更大的话
(**产出**);两个预注册世界都被数据否掉,而第三种答案是它们没覆盖的(**产出:预注册的价值
不在于猜中,而在于让"都不对"这件事看得见**);两个我问错的地方被抓住(**成本回收**)。

**NEXT**:`149e` 把两条线并成了一句话 —— **喜欢什么罕见的东西是稳定的,什么时候得到不是。**
那么下一个可测的问题就很具体:**稳定的那一半(`#100` 的 S)与不稳的那一半(z_resid)之间,
还有关系吗?** `#129e` 测过 `corr(ρ, S) = −0.0417`(匹配后,3.1×),但那是用 ρ 测的;
`#149` 现在说 ρ 与 z_resid 在同一批数据上是同一件事,**所以 `corr(z_resid, S)` 应当也在
那个量级 —— 若不是,那么"同一件事"这个说法在跨到 S 时就破了**,而那会直接影响 `#129e`
这条现存声明的读法。这是一次便宜的一致性检验,而它检的是本项目两条主线的接缝。

---

## Entry 150, added by `E01·A19·R195` — 预注册的点预测命中了;而"一致"仍然不等于"证明了相同"

`#149` 的 NEXT。跑之前写死:既然 `corr(z, ρ) = −0.5718` 是负的,`corr(z_resid, S)` 应当与
`#129e` 的 `corr(ρ, S) = −0.0404` **符号相反、量级相当**。

| # | 结论 | 判定 |
|---|---|---|
| 150a | **预注册的点预测命中** | 预测 **+0.0404**,实测 **+0.0370**,判别量 `corr(z,S)+corr(ρ,S)` = **−0.0034**(**0.2×** 自身展布)。两臂的类别数残差都是 0.003 sd;正对照符号相反且都被拉动(−0.3442 / +0.4982);负对照(按人置换 S)+0.0078 / +0.0030 |
| 150b | ⚠ **但"一致"不等于"证明了相同"** | 判别量的 **95% 上界是 0.0317**,而效应本身只有 **0.0370** —— **这个设计只能排除大于效应 86% 的差异**。判定 **CONSISTENT,不是 PROVEN**。接缝没裂,也没被焊死 |
| 150c | ⚠ **我又用反了一个门,而这次错得有结构** | 假设本身是**「两者相同」**,而 `require_resolvable_first` 是为**「我要它非零」**设计的 —— 它把**想要的结果报成 FAIL**,并把整族标 MOOT。**一个为"证伪零"造的门,不能用来"证实等价"** |
| 150d | **第九个守卫:`equivalent_within`** | 等价界(TOST 式):`|diff| + 2*spread <= margin`,**margin 必须在跑之前指定**。它把 `150b` 那句话变成一个门能自己说出来的判决,而不是我事后加的一段散文 |
| 150e | **立住的那句话** | **「口味越罕见的人,他偏离共享时间表的方向越朝罕见」不依赖我用哪个统计量。** 用"他最早那个有多罕见"测,和用"他整体上罕见的来得多早"测,对同一个人格特质给出一致的答案。`#129e` 的读法不变,而它现在有了第二个独立测量 |

**按 §0.2 结账**:一条现存声明拿到第二个独立测量,并且是**预注册点预测**命中(**产出**);
第九个守卫落地,它把"一致 ≠ 证明相同"变成一个可自动开火的判决(**产出**);
我第二次把可分辨性门用在等价假设上被抓住(**成本回收**)。

**NEXT**:`150b` 的 86% 是这条线现在的分辨率上限,而它**不是数据不够** ——
判别量的展布 0.0142 由**两个相关各自的展布**决定,而它们都在 0.010 上下,n 已近万。
真正的限制是**效应小**(0.037)。所以要把接缝焊死,唯一的办法是**把效应放大**,
而本项目已经有一个现成的放大器:`#126` 的**强制选择**块 —— 它在构造上移除作答水平,
`#126` 用它把稀有亲和从 +0.1822 提到 +0.1449 的干净版本。做法:用强制选择块重建 S,
重跑这两个相关。**若效应放大到 0.10 以上,同样的 n 就能把等价边界压到有意义的水平。**

---

## Entry 151, added by `E01·A19·R196` — 强制选择没放大效应,它更吵;而"更干净"不等于"更好"

`#150` 的 NEXT:用 `#126` 的强制选择块重建 S 来放大效应,把接缝的等价边界压下去。
先处理重叠:10 个强制选择块里有 **8 个就是起始年龄的类别**,所以 z/ρ 在剔除后的 **23 个**
类别上重算(留块法,`#126c`/`#127` 的那个洞),`check_disjoint_items` 断言。

| S 的来源 | corr(z_resid, S) | corr(ρ, S) | 等价边界 |
|---|---:|---:|---:|
| 多选 S(32 块) | **+0.0383** | **−0.0485** | 只能排除 > 效应 **105%** |
| 强制选择 S_fc(10 块) | +0.0223 | **−0.0307** | 只能排除 > 效应 **185%** |

| # | 结论 | 判定 |
|---|---|---|
| 151a | **没放大,缩到 0.63×;等价边界从 105% 恶化到 185%** | `#150` 的 NEXT 建立在一个错的前提上 |
| 151b | **而歧义被分开了:它不是更干净,是更吵** | 我在跑之前写下了这个歧义。独立量:**S_fc 分半信度 0.2520**(5+5 块,Spearman-Brown)vs 多选 S 的 **0.4611**。**衰减预测 0.0485 × √(0.252/0.461) = 0.0359,实测 0.0307**,差 −0.0052,界在 0.0268(边界 0.0243)—— 擦边未过,记 CONSISTENT |
| 151c | **一条可迁移的规则** | **一个仪器"每个题目更干净"不等于"整体更好"。相关的量是信度,而信度是 题目数 × 干净度。** 强制选择在构造上移除了作答水平,但它只有 10 块 —— 净结果是**分辨率的损失**。用更干净但更短的量表去换分辨率,通常是亏的 |
| 151d | **能力边界** | **要焊死这个接缝,需要一个既移除作答水平、又有足够题目数的工具。这个 release 里没有。** 接缝的判定停在 `#150b` 的 CONSISTENT,不再追(`#111c`) |
| 151e | **一个顺带的观察,不解读** | 剔掉那 8 个重叠类别后,多选 S 那一臂的 `corr(ρ,S)` 从 `#150` 的 −0.0404 变成 **−0.0485**(类别集不同,不是同一个量)。**记录,不解读** |

**按 §0.2 结账**:一个放大方案被证伪,而**证伪它的方式给出了它失败的原因与数量**
(信度 0.252 vs 0.461,衰减预测命中)(**产出**);一条关于"干净 ≠ 好"的规则被写下来
(**产出**);一个能力边界被钉住,这条线可以停(**产出:知道在哪停下来,是省下来的成本**)。

**NEXT**:`151d` 关掉了这条线,**换方向**。本会话累计出的最大一句关于人的话现在是
`#148a`+`#149e`+`#150e` 合起来的那句:**「最早那批更常见」的 80% 是共享发育时间表,
剩下 20% 是真的但不是特质;而"口味越罕见的人偏离得越朝罕见"这句话不依赖用哪个统计量。**
它的**唯一薄弱处**是 `#148` 用的那张时间表 —— `#75` 的题目均值序,而 `#131c` 已经证明
**中位数序与均值序不同,且中位数序与稀有度对齐而均值序不对齐**。做法:把 `#148` 的三臂
分解用**中位数序**重跑一次。若"解释掉 80%"这个数随汇总量大幅变动,那么 `#148a` 的
**量级**要加限定;若它稳,`#148a` 就拿到了它现在缺的规格稳健性。

---

## Entry 152, added by `E01·A19·R197` — 方向立住,量级撤回:`#148a` 的分母错了,比值被高估 14 倍

`#151` 的 NEXT。用五个汇总量重跑 `#148` 的分解,并把稀释基准从一次抽样改成 20 次。

| 汇总量 | ρ(与稀有度) | 解释掉 | 稀释基准的 20 抽区间 |
|---|---:|---:|---|
| median | **+0.437** | **104%** | −90% .. +92% |
| q75 | +0.117 | 98% | −137% .. +184% |
| trim20 | −0.014 | 83% | −89% .. +90% |
| **mean(`#148` 用的)** | +0.011 | **80%** | **−61% .. +118%** |
| q25 | +0.329 | **61%** | −80% .. +66% |

| # | 结论 | 判定 |
|---|---|---|
| 152a | **预注册的方向预测命中** | 中位数序与稀有度对齐(+0.437)而均值序不(+0.011),预测中位数解释更多 —— 实测 **104% vs 80%** |
| 152b | **`#148a` 的方向立住,而且比原来强** | **五个汇总量全部解释掉一个大的正比例(61%–104%),符号一致。"共享时间表是机制"在五种规格下都稳** |
| 152c | ⚠ **但「80%」是规格依赖的** | 极差 **43 个百分点**。**必须写成区间 61%–104%** |
| 152d | ⚠ **而「44 个百分点,20.0×」用错了分母 —— 撤回** | `#148d` 的稀释基准是**一次**随机题目向量;重抽 20 次,均值臂的稀释解释比例跨 **−61% 到 +118%**,**`#148` 的 36% 只是其中一个点**。我判 20.0× 时用的分母是**按人自助**展布(0.0098),而这个差的不确定性由**基准自身实现**主导(展布 0.197 z 单位 = **44 个百分点**)。**正确的读法:净额 61 ± 44 个百分点 = 1.4×,不可分辨** |
| 152e | **一条可迁移的规则** | **当一个比较的基准是一次随机实现时,判它的分母是那个基准的实现展布,不是被比较量的抽样展布。** 用错它,把 1.4× 报成了 20.0× —— **高估 14 倍** |

**按 §0.2 结账**:`#148a` 的**方向**在五个规格下拿到稳健性,而且预注册的方向预测命中
(**产出**);它的**量级**被撤回并换成区间(**成本回收**);一个分母错误被定位并写成规则
(**产出 —— 它适用于本项目每一个"对照基准是一次抽样"的比较**)。

**NEXT**:`152e` 那条规则应当立刻回头扫一遍 —— **本项目还有多少个"净额/超出基准"的声明,
其基准是一次随机实现,而判它用的是被比较量自己的展布?** `tools/guard_lint.py` 已经有
扫描框架。可查的模式很具体:`offset_control(effect, offset, spread)` 里,`offset` 来自
一次抽样而 `spread` 来自 bootstrap。**`#148` 是第一个被抓到的;`#152e` 说这类比较在本项目
里到处都是。** 这是 Closure,但 `#148` 已经证明它能把一个比值改掉一个数量级。

---

## Entry 153, added by `E01·A20·R198` — 单次抽样的零不自动致命;伤害 = 零的实现 sd ÷ 效应

`#152` 的 NEXT。扫 139 个持久化结果文件:**53 个含零臂,17 个的零臂没有 seed 列**。
分诊后,两个支撑现存声明的零确认是单次抽样:**`#114`**(`default_rng(9)`,一次同类别内
打乱)与 **`#101`/`#102`**(`curveball(M, default_rng(8100))`,一次保边际实现)。

| # | 结论 | 判定 |
|---|---|---|
| 153a | **`#114` 的零重抽 20 次:原轮那一抽是最小的** | 效应 −0.2000 ± 0.0096(20.8×,逐位复现);原轮 seed 9 给 **−0.00035 = 效应的 0.2%**;20 抽范围 **0.2% .. 11.4%**,实现 sd = 效应的 **5.5%** |
| 153b | **但声明站住** | 最不利的一抽也只有 **11.4%**,`negative_control` 通过(门槛 50%)。**`#114` 不动,改的是它旁边那个数:「零 = 效应的 0%」→「0–11%,零的实现 sd 是效应的 5.5%」** |
| 153c | **一条比 `#152e` 更平衡的规则** | **单次抽样的零不自动致命,伤害 = 零的实现 sd ÷ 效应。** `#148`:基准 sd **44 个百分点** vs 净额 61 → **毁掉那个比值**;`#114`:零的 sd **5.5%** vs 效应 100% → **无害,只是让数字看起来比实际好**。**报单抽零时必须带上这个比值** |
| 153d | **扫描本身的范围** | 「没有 seed 列」**不等于**单次抽样(零可能在内部平均后才落盘)。17 个命中是**待读清单**,而分诊只走到了两个支撑现存声明的(`P6` 安全侧) |
| 153e | **未做** | `#101`/`#102` 的 curveball 零(23 块 × 20 抽,比本轮贵一个量级)留给下一轮 |

**按 §0.2 结账**:一条现存声明的零被重抽并保住(**产出**),而它旁边那个"0%"被换成诚实的
区间(**成本回收**);`#152e` 那条规则拿到它缺的另一半 —— **什么时候单抽是无害的**
(**产出:一条只会否定的规则会让人过度撤回**);扫描落地并明确了它自己看不见什么(**产出**)。

**NEXT**:`153e` —— `#101`/`#102` 的 curveball 零重抽 20 次。它比 `#114` 更值得查,原因有二:
① 它支撑的是 README 上「**它唯一挂得住的外部锚是性别**」这句话,而那句话的力量**全部来自
零**(「人格五因素全部 |r| ≤ 0.056」是一个**否定**,而否定的强度完全取决于零有多干净);
② curveball 是**保边际**的,它的实现方差比一次简单打乱**更难直觉**,而 `#105c` 已经证明
保边际零对某些统计量**结构上盲**。**若那 20 抽的展布把 0.056 抬到与 0.093 同量级,
那么「唯一的外部锚」这句话就要改。**

---

## Entry 154, added by `E01·A20·R199` — 我重抽的是一个已被取代的轮次;而那 10 抽仍然证明了 `R15` 的三次平均是必要的

`#153e` 的路标指向 `A11/R148` 的单抽 curveball 零。

| # | 结论 | 判定 |
|---|---|---|
| 154a | ⚠ **正对照当场开火:我复现不出 README 的数** | 我的真实臂给性别 **+0.0558**、agreeable **−0.0741**,而 README 写 +0.093、五因素 ≤0.056。**跑 `A11/R148` 本身:我的复现是精确的** —— 它自己就输出这些数,**而且它自己的判定是 UNVERIFIED,两个门 FAIL** |
| 154b | **读账本settled了它** | `Entry 101`:`R14` 的两个门失败**是我的错**(`101a` 线性残差化、`101b` 原始比去衰减);**`R15` 的修法**是 `affinity = S_real − mean(S_null)`,**三次** curveball 平均,**两个门都过**。README 的 +0.093 / ≤0.056 / 开放性 +0.023 全部出自 **R15**。**声明没有问题** |
| 154c | ⚠ **第三次「去看轮次而不是去读账本」** | `#143`(重跑 `A02/R034` 审计一个 118 条前已作废的声明)· `#148e`(两轮找一个账本已定位的机制)· 本轮。**而 `#153d` 刚写下「没有 seed 列不等于单次抽样」,我写了那句话,然后走过去了**(`feedback_confession_is_never_audited` 的同一形态) |
| 154d | **但那 10 抽仍然回答了本来的问题,而且答案是"平均是必要的"** | curveball 零臂的人格 max\|r\| 在 10 次独立实现上跨 **0.0320–0.0543**,实现 sd **0.0068**。**一次**零自己就能产生高到 **0.054** 的人格相关。`R15` 平均三次 → 残余 ≈ 0.004,而效应 +0.0930 → **23 倍余量** |
| 154e | **所以这条声明安全,而它安全的原因是设计而不是运气** | **用一次抽样的零,「人格五因素全部 ≤0.056」这句话会与零分不开**(零自己能到 0.054)。`R15` 的三次平均正是把它救出来的东西 |

**按 §0.2 结账**:一条现存声明的零被量化并确认足够(**产出:而且量出了"为什么三次是必要的"**);
我第三次走"轮次而非账本"的路被抓住,而抓住它的是**我自己上一轮写下却没遵守的告诫**
(**成本回收**);`#153` 那个扫描的假线索被定性(**产出:一个工具的盲区,要在它第一次
误导你时就记下来**)。

**NEXT**:`154c` 是本会话第三次同一个错误,而它有一个**机械的**修法 ——
`tools/readme_ledger_audit.py` 已经能把 README 的数字连回账本;缺的是**反方向**:
**给定一个轮次目录,它支撑的声明当前的账本状态是什么(是否已被后续条目取代)。**
做法:扫 `RETRACTIONS.md` 里所有 `added by \`E01·Axx·Ryy\`` 头,建 轮次 → 条目 的映射,
再对每个轮次报出**引用它的最后一条条目**。有了这张表,"我要重跑 R14"时一眼就能看到
"`Entry 101` 说 R15 取代了它"。**这个错误犯三次之后,应该由工具而不是由纪律来防。**

---

## Entry 155, added by `E01·A20·R200` — 犯三次之后,由工具防:轮次 → 账本的反向索引

`#154` 的 NEXT。本会话三次「去看轮次而不是读账本」(`#143` · `#148e` · `#154`)。

| # | 结论 | 判定 |
|---|---|---|
| 155a | **`tools/round_status.py` 落地** | 轮次 → 账本,顶出**最后一条提到它的条目**。关键是**兄弟轮次检测**,而它是照着骗过我的那一格调的:`Entry 101` 正文写的是**裸 `R15`**,按引用格式匹配抓不到 —— 改成扫**这一轮自己那条条目的正文**找兄弟轮次 + 取代类词汇。**现在查 `A11/R148` 直接打出「Entry 101 正文提到兄弟轮次 R15」** |
| 155b | **全量扫描:139 个轮次进过账本,11 个带取代提示** | 特异性够(不是一片红) |
| 155c | **压在带标记轮次上的现存声明只有 1 个** | `#101` → `A11/R148`,**而 `#154` 上一轮已经解决了它**(README 的数出自 `R15`)。**其余 16 条现存声明都不压在带标记的轮次上** |
| 155d | **盘上 195 个轮次,账本里从未出现过 59 个** | **但没有一个支撑现存声明**。它们是走过又没有留下信念更新的路 —— 记录,不处理 |
| 155e | **范围(`P6` 安全侧)** | 只在**命中**方向可读。**没有取代提示不等于这一轮是现行的**;输出是**必读清单**,不是判决 |

**按 §0.2 结账**:一个犯了三次的错误现在由工具防而不是由纪律防(**产出**);
"有没有现存声明压在被取代的轮次上"这个问题第一次被**机械地**回答,答案是**只有一个,而且
已解决**(**产出:一个闭合的决定**);59 个从未进账本的轮次被点出来(**产出:范围**)。

**NEXT**:三个工具现在覆盖三条边 —— `readme_ledger_audit`(README → 账本)·
`round_status`(轮次 → 账本)· `guard_lint`(轮次 → 守卫)。**缺的第四条是 账本 → 轮次**:
**一条账本条目声称的东西,那一轮的代码现在还跑得出来吗?** `#154a` 正是靠"跑一遍 `A11/R148`
本身"才发现我复现的是对的、README 才是错的 —— 那一步是手工的。做法:对**现存声明的 17 轮**
逐个重跑,把 stdout 与账本条目里引用的数字比对,报不一致。**`#141`、`#154` 都是在这一步上
被绊住的,而它现在完全没有自动化。**

---

## Entry 156, added by `E01·A20·R201` — `#117` 的轮次从写下那天起一直在崩,而崩的正是它自己记录的那个 bug

`#155` 的 NEXT:第四条边,**账本 → 轮次**。`tools/repro_audit.py` 重跑一轮,把账本条目里的
可辨识数字与新 stdout 比对(容差 2%)。**17 条现存声明,16 条跑通,1 条崩溃。**

| # | 结论 | 判定 |
|---|---|---|
| 156a | ⚠ **`#117` 的轮次 `rc=1`** | `A12/R169` 崩在 `sh=G.shift.values` → `AttributeError: 'function' object has no attribute 'values'`。**`shift` 是 pandas 的 DataFrame 方法** |
| 156b | ⚠ **而那正是 `#117e` 自己记录的 bug** | `#117e` 原文:「**第五次 pandas 访问器撞名,而这次撞的是我自己禁用名单上的第一个词**… **写下规则第三次没有阻止我违反它**」。**我诊断了它、为它写了 `check_columns` 守卫、把结论记进账本 —— 却没有把那一行改掉。** 轮次从那时起一直崩,而崩溃在**结论打印之后**,所以退出码没人看 |
| 156c | **两个工具在同一个具体缺陷上会合** | `guard_lint`(`#128a`)早就把 `A12/R169` 标为缺 `check_columns`;`repro_audit` 显示它**确实**崩在那一类 bug 上。**四个工具第一次相互印证** |
| 156d | **而声明本身完好** | 修好两处并**把守卫装上**后,那个从没跑到过的段落给出位移 **+0.0266 ± 0.0127**、零 **+0.0067**(25%)、**29/45**、**2.1×** —— 与 `Entry 117` **逐字相同**。**代码卫生的缺陷,不是声明的缺陷** |
| 156e | **"对不上"是必读清单,不是判决** | 账本条目会引用**前面条目**的数字、散文里的数字、以及本轮没跑到的臂。`#101` 的 21 个就是这样(条目覆盖 `R14`+`R15`,我只跑了 `R15`)。`#63` 是 **0/10**,`#116`/`#69`/`#100` 各 2 —— **多数轮次复现得很干净** |
| 156f | **顺带:59 个从未进账本的轮次,几乎全在早期弧** | A01 12/24 · A02 17/25 · A04 11/12 · A05 9/15,而 **A07 往后 100 来轮里只有 2 个**。**账本纪律是从 A07 开始的** —— 那不是"走过没留下更新",是**那时还没有记更新的做法** |

**按 §0.2 结账**:第四条边落地,而它第一次跑就抓到一个现存声明的轮次**从写下那天起一直在崩**
(**产出**);那个 bug 被修,守卫被**装上**而不只是存在(**成本回收**);声明本身被逐字复现确认
(**产出:一条声明第一次有了"它现在还跑得出来"的证据**);四个工具第一次在同一个缺陷上会合
(**产出**)。

**NEXT**:`156e` 里 `#88`(18/30 对不上)与 `#102`(14/19)是仅有的两个高不符率、
且**不能**用"条目覆盖多轮"解释的。`#88` 是 README 上一条现存声明
(**内容与个体化估值对认可概率的影响一样大:±22.6 / ±16.3 / ±23.7 pp**),
而它有 60% 的账本数字在新输出里找不到。**先查它** —— 这正是 `repro_audit` 造出来要抓的东西,
而它第一次跑就指到了本项目最"科学"的那条声明上。

---

## Entry 157, added by `E01·A20·R202` — 一条只写在账本散文里的校正,不会到达仪器(两轮之内第二例)

`#156` 的 NEXT 指向 `#88`。查下去,真正的问题在 **`#90`** —— README 头条
「±22.6 / ±16.3 / ±23.7 pp」的出处。

| # | 结论 | 判定 |
|---|---|---|
| 157a | ⚠ **`repro_audit` v1 看不见它要查的那种格式** | 数字正则要求 **≥2 位小数**,而 `22.6`/`16.3`/`23.7` 只有 1 位。它报 `#90` "5/6 对不上",而那 6 个里**根本不含三个头条数**。**一个检查看不见它被造出来要检查的格式** —— 已放宽 |
| 157b | ⚠ **这一轮打印的,正是它自己账本判为错的那一版** | `#90c` 把 clip 前的种植 sd 判为「第十七个 mis-specified statistic」,并记下校正后的 **±23.7 / CV 15.9%** —— **代码从没改**,至今打印 **±30.8 / CV 22.4%** |
| 157c | ⚠ **`#90d` 的噪声校正同样只活在散文里** | 人层 19.8 → 16.3 的二项噪声校正不在代码里,输出一直是未校正的 20.1 |
| 157d | **两条校正装进代码之后** | 交互 **±23.7 逐位命中**;人层 **±16.7**(噪声 11.2,账本记 10.4);族 CV 19.5%(账本 15.9%) |
| 157e | ⚠ **剩下的 5% 记录不猜** | 选项基率 **±22.6 → ±21.5**。`item_pp` 是**确定性**的,块数**仍是 23** 与账本一致,**这个漂移无法用块集变化解释**。**记录,不用推测解决** |
| 157f | **声明本身站住** | ±21.5 vs ±23.7 **仍然一样大**,而"一样大"正是这条声明说的话 |
| 157g | **一条规则,两轮之内两例** | **一条只写在账本散文里的校正,不会到达仪器。** `#117e` 是**崩溃**(退出码非零),`#90c` 是 **exit 0 的静默** —— 后者更难看。**诊断、记账、写守卫,都不等于改代码** |

**按 §0.2 结账**:README 头条三个数里两个被恢复(其中交互逐位命中),一个被标为未解释的
5% 漂移(**产出:一条声明第一次有了"代码现在能不能生成它"的答案**);两条只存在于散文里的
校正被装进仪器(**产出**);审计工具自身的一个格式盲区被修(**成本回收**);
`#117e`/`#90c` 合成一条规则(**产出**)。

**NEXT**:`157g` 说的病现在有两例,而它有一个**机械的**检测法 —— 账本条目里出现
「校正后 / corrected, the X falls from A to B」这类**成对数字**,而代码只能产生 A 不能产生 B。
`repro_audit` 已经能跑轮次并抓数字;缺的是**从账本条目里抽取"被判为错的旧值 A"**,
然后断言 **A 不应该出现在新输出里**。**现在的方向是反的**:我查的是"账本的数在不在新输出里",
而应该同时查"**账本判为错的数,还在不在新输出里**"。`#90` 的 ±30.8 会被这条立刻抓到。

---

## Entry 158, added by `E01·A20·R203` — 反向检验落地,并且它先证明了上一轮的修复

**【CLOSURE,明确标注】** `#157` 的 NEXT:`repro_audit` 的方向是反的 —— 只查「账本的数在不在
新输出里」,所以对 `#90` 的 ±30.8 无感,**而那个数本该消失**。

| # | 结论 | 判定 |
|---|---|---|
| 158a | ⚠ **抽取器连 `#90c` 自己那句话都抓不到** | v1 的正则被 markdown 的 `**`、`±`、`%` 打断(`falls from ±30.8 to **±23.7 pp**`)。放宽后全账本 **7 条**含校正对 |
| 158b | **反向检验先证明了上一轮的修复** | `#90` 的旧值 **30.8** 与 **22.4** 在新输出里各 **0 次**。**上一轮的修复被一个独立检查确认** |
| 158c | ⚠ **三类假阳性,每一类都是跑出来才知道的** | ① **参数→结果**(`Entry 69`:"at **g=0.15** moves ρ to −0.314")→ 排除赋值前缀;② **一条条目覆盖多轮,旧值在前一轮**(`Entry 126`:`R20` 打 +0.1822、`R21` 打 +0.1449)→ 判据改为「新值在该条目**任何一轮**里出现即良性」;③ **一轮本来就该报两个臂** → **不可机械区分,留在必读清单**(`P6`) |
| 158d | **结果** | **7 条含校正对,只有 `#90` 是现存声明,而它已修并已验证**。`#126` 在修好的判据下正确读作良性。**没有其它现存声明带着一个还活着的、被账本判为错的数** |
| 158e | **这一轮不产出关于人的新事实,而它该被这么记账** | 它保护的是 `#156`/`#157` 刚修好的两条声明。**Closure,标明** |

**按 §0.2 结账**:`#156`/`#157` 的两处修复现在由一个自动检查守着,而不是由我记得
(**产出**);"还有没有别的现存声明带着已判错的数"第一次被机械回答,答案是**没有**
(**产出:一个闭合的决定**);抽取器自身的格式盲区与三类假阳性被逐个定性(**成本回收**)。

**NEXT**:五条边现在齐了(README→账本 · 轮次→账本 · 轮次→守卫 · 账本→轮次 · 账本判错的数→轮次)。
**审计基建到此为止,换回人。** 本会话在人这一侧留下的最大开口是 `#149d`:
**扣掉共享时间表之后那 20% 的个人成分,是真的但不是一个特质**(同一批兴趣上两个读数
一致 −0.574,换到另一半兴趣只剩 −0.091)。而 `#149` 只证明了它**不稳**,没问**它为什么不稳**。
一个具体的候选:**它可能是"这个人在这一批兴趣上恰好怎么回答"** ——
即**块级**的东西,而不是人级的。可测版本:把 z_resid 按**块**而不是按人算,
看它在块之间的一致性是不是高于在人之间的。**若是,那 20% 就有了名字:它是块效应,不是人效应。**

---

## Entry 159, added by `E01·A21·R204` — 那 20% 几乎不是关于人的,它是关于东西的

`#158` 的 NEXT,换回人。`#149d` 只证明了个人成分**不稳**,没问**为什么**。
双向去均值后两个边际都为零,所以"是不是关于东西的"只能是:**交互是不是集中在特定类别上**。

| 边 | 真实 | 人内置换零 |
|---|---:|---:|
| ① **人侧**(一个人的类别劈两半,跨人相关) | **+0.0556**(SB **+0.1053**) | +0.0006 |
| ② **题目侧**(人劈两半,两套类别载荷相关) | **+0.9616** | −0.0867 |

| # | 结论 | 判定 |
|---|---|---|
| 159a | **相差 17 倍,答案是 ITEM** | 人侧 SB 信度 **0.105**;题目侧 **0.96** |
| 159b | **而它不是平凡的** | `corr(λ_j, 稀有度) = **−0.8948**` —— λ 大半是稀有度的函数,而稀有度在两半人里恒同。把 `rar` 与 `rar²` 都回归掉后,两半的残差 λ **仍相关 +0.8119** |
| 159c | **这给 `#149d` 改了名** | 它**不是**一条弱的人格特质,它是一条**强的、可复现的类别层结构**。**某些性兴趣就是会比共享时间表预期的更早或更晚到来,而这对完全不同的两组人是同一套图样;至于"哪个人更容易这样",几乎测不出来** |
| 159d | ⚠ **边②没有开火的正对照** | 种植类别特异信号只把它从 0.9616 推到 0.9753 —— **它已在天花板上**。边②的读数靠 `159b` 的平凡性检验,**不是靠一个通过的正对照**。**记录,不当作它通过了** |
| 159e | **内容读不出一条轴,所以不命名** | 去稀有度后载荷两极:一极 vore · 权力动态 · 急切 · 非自愿 · 羞辱(像关系族);另一极把**最常见的**「开始看色情」(89%)与**最罕见的**残暴(5%)、惊悚(4%)混在一起。**读不出一条轴,就不给它安一个**(`#140e` 的同一条纪律) |

**按 §0.2 结账**:一条挂了三轮的开放问题被关掉,而且答案换掉了它的名字(**产出**);
它经受住了一个会让它变平凡的检验(**产出**);边②缺一个开火的正对照被写下来而不是绕过
(**成本回收**);载荷的内容被拒绝命名(**成本回收:一个不被命名的图样,比一个被硬安上
名字的图样更可用**)。

**NEXT**:`159a` 的两个数(人侧 0.105 / 题目侧 0.96)现在是**同一个交互**的两种归属,
而本项目在**另一条线**上有一模一样的一对:`#100` 的稀有亲和特质**人侧**分半是 **+0.4611**。
**同一个人身上,"喜欢什么罕见的东西"人侧信度 0.46,"什么时候得到它们"人侧信度 0.105 ——
差 4.4 倍。** 可测的下一步:把两者放在**同一个设计**里(同一批人、同一种分半、同一个
Spearman-Brown),确认那 4.4 倍不是两套口径的产物。**若确认,本项目就有了一句很干净的话:
性偏好的"什么"是个人的,"何时"不是。**

---

## Entry 160, added by `E01·A21·R205` — 性偏好的「什么」是个人的,「何时」不是

`#159` 的 NEXT:把 `#100` 的 S(0.4611)与 `#159a` 的 z(0.105)放进**同一个设计**,
确认那 4.4 倍不是两套口径的产物。⚠ 先处理**单位数**:S 用 32 个块、z 用 ~13 个类别,
而信度随单位数涨,SB 只校正"劈成两半"。所以在**同一个 k** 上重量。

| 每半单位数 k | **S(喜欢什么)** | **z(何时得到)** | 比值 |
|---:|---:|---:|---:|
| 4 | **0.5742** | 0.0550 | **10.4×** |
| 5 | 0.6121 | 0.0840 | 7.3× |
| 6 | 0.6271 | 0.1075 | 5.8× |
| 7 | **0.6671** | 0.1056 | 6.3× |

| # | 结论 | 判定 |
|---|---|---|
| 160a | **4.4 倍是真的,而且被低估了** | 匹配 k 后是 **5.8–10.4 倍**;平均差 **+0.5321,26.6× 自身展布**,四个 k 上符号一致 |
| 160b | **而这不是仪器不行 —— 这一条才是关键** | 往 z 里种一个**人特异**信号,同一个仪器把信度从 0.055–0.108 推到 **0.67**;S 从 0.57–0.67 推到 0.99。**所以 z 的低信度是数据的性质,不是设计的局限** |
| 160c | **一句关于人的话** | **性偏好的「什么」是个人的,「何时」不是。** 你喜欢多罕见的东西,是一件跟着你走的性质(0.57–0.67);你的罕见兴趣比共享时间表早多少,几乎不是(0.055–0.108) |
| 160d | **范围,预注册过,不假装解决** | 匹配 k 拉平了**单位个数**,拉不平**每单位信息量**(一个块十几个选项 vs 一个类别一个年龄),**而这条差异的方向是让 S 占便宜**。反过来说:即使全算作 S 的优势,**z 的种植臂 0.67 仍然证明仪器有力** |

**按 §0.2 结账**:本会话最长的一条线(`#148` → `#149` → `#159` → `#160`)收束成**一句关于人的话**,
并且带着它的匹配设计、两个开火的种植对照与一条写明的范围(**产出**);而它比它的前身更强,
不是更弱(**产出**)。

**NEXT**:`160c` 现在是本项目最干净的一句话,而它有一个**可以立刻问的推论** ——
若「何时」几乎没有人层成分,那么 `#75` 的**发育时间表**(人群共享,66.852% 成对顺序准确率)
就应当**几乎榨干**了个人的顺序信息。`#63` 已经量过全体对上界 **[60.5%, 66.5%]**,
而 `#75` 达到 **66.852%** —— **它已经在上界之上**。
可测的下一步:把 `#63` 的上界**按人重算**(每个人自己的顺序有多少能被任何共享顺序解释),
看那个上界与 `#160` 的 0.105 是不是同一件事的两种说法。
**若是,本项目就有了一个闭合:「何时」的人层成分有多小,两条独立的路给出同一个数。**

---

## Entry 161, added by `E01·A21·R206` — 关于「何时」,两个完全不同的问法给出同一个小数字

`#160` 的 NEXT。账本里已有第三个可比的数:**`Entry 8`** 把「时间表贴合度是一种个人特质」
杀掉,SB = **0.214**。三个量放进**同一台机器**(同一批人、同一种劈分、同一个 k、同一个 SB)。

| 每半单位数 k | **S 喜欢什么** | **a 时间表贴合你** | **z 你偏离多远** |
|---:|---:|---:|---:|
| 4 | 0.5663 | 0.1182 | 0.1124 |
| 6 | 0.6330 | 0.1960 | 0.1065 |
| 7 | **0.6817** | 0.1935 | 0.1610 |
| **均值** | **0.6221** | **0.1634** | **0.1140** |

| # | 结论 | 判定 |
|---|---|---|
| 161a | **关于「何时」的两个完全不同的问法给出同一个小数字** | a 与 z 差 **+0.0494,1.3× 展布** —— 分不开。**「共享时间表贴合你多少」与「你偏离它多远」,作为个人性质是同一个量级的小** |
| 161b | **而它们与「什么」差一个数量级** | a 与 S 差 **+0.4587,12.2×**。**一个人的性偏好里,「什么」是他自己的;「何时」几乎都不是** |
| 161c | **`Entry 8` 的 0.214 立得住** | 这台机器上 k=6,7 给 **0.194–0.196**。那条早期撤回(「adherence 的相关是无物的相关」)的**基础数字**在新框架下复现 |
| 161d | **三个种植全部开火,三个零都在零附近** | S 0.57→0.99 · z 0.11→0.66 · a 0.12→0.55;置换零全部在 ±0.06 内。**所以三条曲线的相对位置可读** |
| 161e | **范围** | `a ≈ z` 是 **CONSISTENT 不是 PROVEN**:等价界只能把差界在 **0.124**,预设边界 0.10(与 `#150b` 同形)。单位信息量的差异继续携带(`#160d`) |

**按 §0.2 结账**:`#160c` 那句话拿到第二条独立的路,而两条路给出同一个小数字(**产出**);
一条早期撤回条目的基础数字在新框架下复现(**产出:一个五年前的撤回第一次被现在的机器确认**);
`a ≈ z` 的等价只报到设计能给的分辨率(**成本回收**)。

**NEXT**:`161b` 把本项目的两条主线并成一句话,而这句话有一个**尚未被问过的推论** ——
若「何时」几乎没有人层成分,那么 `#114`(**人把最爱的兴趣记得更早,约九个月**,19.8×)
就必须是一条**关于人群而不是关于个人**的规律。可测版本:把 `#114` 的回忆偏差斜率
**按人估计**,量它的分半信度,放进同一台机器。
    若 ~0.1  -> 回忆偏差也是人群层的,而"某些人比别人更会把最爱的记早"这句话不成立
    若 ~0.6  -> 回忆偏差是个人的,那么「何时」里**有**一个人层成分,只是不在顺序里而在**记忆**里
**这个分叉决定 `#114`/`#119` 那条线该被怎么读。**

---

## Entry 162, added by `E01·A21·R207` — 回忆偏差是所有人共有的;而它的效应比它自己零的噪声还小

`#161` 的 NEXT。`b_i` = 这个人**人内**的「起始年龄残差 ~ 该类别评分」斜率,量它的分半信度,
放进 `#160`/`#161` 的同一台机器。

| 每半单位数 k | **SB 信度** | 人内置换零(5 抽) | 人特异种植 |
|---:|---:|---:|---:|
| 4 | +0.0646 | −0.0060 ± 0.031 | **0.4121** |
| 5 | +0.0365 | +0.0265 ± 0.071 | **0.4347** |
| 6 | +0.0614 | −0.0162 ± 0.073 | **0.5471** |

| # | 结论 | 判定 |
|---|---|---|
| 162a | **人群层斜率复现 `#114`** | **−0.1975** vs `#114` 报的 −0.2000 |
| 162b | ⚠ **而「b = 0.054」不能当一个数来引用** | **零的实现 sd ÷ 效应 = 1.1×** —— **效应比它自己零的实现展布还小**。可以说的只有:**它在它自己零的噪声之内**。这是 `#153c` 那个比值在本项目至今最大的一次 |
| 162c | **但这是一个真的缺席,不是没检出力** | **人特异种植开到 0.41–0.55**,所以仪器**测得到**人层的回忆偏差 |
| 162d | **判定 POPULATION** | **「某些人比别人更会把最爱的记早」这句话不成立。回忆偏差是所有人共有的** |
| 162e | **`b` 与 `z` 的等价第一次真的过了预设边界** | \|差\| + 2sd = **0.0906 < 0.10**。两者是同一个量级的「无」。而 `b` 与 `S`(喜欢什么)差 **36.9×** |
| 162f | ⚠ **我在这一轮里犯了自己上一轮写下的错** | 第一版每个 k 只抽**一次**置换零,k=5 那一抽给 **+0.106**(k=4/6 给 −0.03/−0.05),门当场 FAIL。**`#153c` 是我上一轮写下的规则,一轮之后我自己踩了它** |

**按 §0.2 结账**:`#114` 那条声明被重新定位为**人群规律**,并且带着一个能证明"仪器有力"的
种植臂(**产出**);本项目第一次有一个等价检验**真的过了预设边界**(**产出**);
一个点估计被禁止引用,因为它比自己零的噪声还小(**成本回收**);
上一轮写下的规则一轮之后抓到了我自己(**成本回收**)。

**NEXT**:`162d` 让「何时」这一侧的四个量全部落位 —— S 0.62 · a 0.16 · z 0.11 · b 0.054,
**而后三个彼此分不开**。那么下一个问题不是"还有没有别的人层成分",而是
**「什么」那 0.62 里,有多少是它自己的?** `#100` 的 S 是"勾选项的平均冷门程度",
它同时含**两件事**:你**挑了哪些**(内容)与你**挑得多冷门**(位置)。可测版本:
把 S 拆成「内容成分」(与哪些具体选项被勾相关)与「位置成分」(与所勾选项的基率相关),
在同一台机器上各量一次分半信度。**若 0.62 几乎全在内容成分上,那么"稀有亲和"这个名字是错的
—— 稳的是"你喜欢哪些东西",而不是"你偏爱冷门"。这会直接改 `#100` 的措辞。**

---

## Entry 163, added by `E01·A22·R208` — 「偏爱冷门」是一个真的口味位置,不是"喜欢某一类恰好冷门的东西"

`#162` 的 NEXT:`#100` 的 S 同时含**内容**与**位置**,而它的名字从没被检验过。
分离器:**把块按内容相似度分到最不像的两端**,看倾向能不能传过去。

| | 随机分半 | **内容最不相似分半** | 保留 |
|---|---:|---:|---:|
| **真实** | **+0.5973** | **+0.5059** | **85%** |
| 纯**位置**种植 | +0.9859 | +0.9855 | **100%** |
| 纯**内容**种植 | +0.9268 | **−0.0212** | **−2%** |

| # | 结论 | 判定 |
|---|---|---|
| 163a | **POSITION。`#100` 的名字是对的,而这是它第一次被检验** | 真实保留 **85%**,而两个种植**完美分开**(位置 100% / 内容 −2%)—— 设计**能**分辨,真实数据落在位置那一侧。**「偏爱冷门」是一个真的、可迁移的口味位置,把话题换到最不像的另一端仍带得过去** |
| 163b | **但存在一个小的内容成分** | 掉的那 15% 是 **5.7× 自身展布**,不是噪声 |
| 163c | ⚠ **我的置换零是一次改标签,不是一次破坏** | 第一版 `x = Bc[随机另一个人]` —— 换的是**整条向量**,两半仍来自同一个人,**置换什么也没破坏**。零给到 **+0.38 = 效应的 71%**。正确的零是**每个块独立跨人置换**;改后 **+0.017**(5%) |
| 163d | ⚠ **种植建在真实信号之上,所以塌不下去** | 内容种植第一版只掉到 81%。改成建在**零背景**上,并把"不相似分半"改成**显式跨族**(第一版取"这个人自己块里 pc1 最低/最高的 k 个",块少的人两端落在同一族)。改后 **−2%** |
| 163e | **范围,预注册过** | "内容相似度"由共现给出,而共现与冷门程度相关(`#159b`:−0.89),所以最不相似的两半在稀有度构成上也不同。**由两个种植对照兜底** —— 若那个相关是主要驱动,纯内容种植不会塌到 −2% |

**按 §0.2 结账**:本项目唯一一个大的人层量,它的**名字**第一次被检验并通过(**产出**);
一个能分辨位置与内容的设计落地,两个种植 100% vs −2%(**产出**);
一个把零变成改标签的错误被抓住,而它曾让零虚报到效应的 71%(**成本回收**)。

**NEXT**:`163a` 让「什么」这一侧的图景完整了 —— **偏爱冷门是一个可迁移的位置倾向,
分半信度 0.60,跨内容最不像的两端仍保 85%**。而它现在有一个**没被问过的对手**:
`#101`/`#154` 说它唯一挂得住的外部锚是**性别**(+0.093),而性别本身也是一个"跨话题都成立"
的位置变量。可测版本:**把性别从每个块层分数里回归掉之后重跑这一轮**。
    若 85% 基本不变 -> 位置倾向不是性别的影子,它是自己的东西
    若塌下去       -> 「偏爱冷门」很大程度上是"男性在每个块上都挑更冷门的"的另一种说法,
                    而 `#100` 的 0.46 需要重新定价
**这是本项目最老的那条现存声明的最后一个未做的对照。**

---

## Entry 164, added by `E01·A22·R209` — 「偏爱冷门」是这个人自己的口味,不是他的性别

`#163` 的 NEXT。`#163a` 立住的是一个**位置**倾向,而性别本身也是一个"跨话题都成立"的位置变量
—— 男性若在**每个**块上都挑更冷门的,就会造出一个完美可迁移的位置倾向。
做法:把性别从**每一个块层分数**里回归掉,重跑 `#163` 的整台机器。

| | 去性别前 | **去性别后** | 保留 |
|---|---:|---:|---:|
| **内容最不相似分半** | +0.5059 | **+0.5166** | **102%** |
| 随机分半 | +0.5894 | +0.5914 | 100% |
| **纯性别驱动的位置种植** | +0.9439 | **−0.0009** | **−0%** |
| 纯**个人**位置种植(去性别后) | — | **+0.9858** | — |

| # | 结论 | 判定 |
|---|---|---|
| 164a | **INDEPENDENT。「偏爱冷门」是这个人自己的口味,不是他的性别** | 去性别后**保留 102%**(**4.9× 自身展布**,方向是**升**不是降);性别在块层分数里只解释掉中位 **1.32%** 的方差 |
| 164b | **两个对照必须一起看,而两个都开火** | 纯**性别驱动**的位置种植塌到 **−0%** → 残差化**确实有效**;纯**个人**位置种植存活到 **+0.9858** → 残差化**没有误伤**。**两个都开火,102% 才可读** |
| 164c | **`#100` 的 0.46 不动,而这是它最后一个未做的对照** | 本项目最老的那条现存声明,名字(`#163`)与最强的混淆(本轮)现在都被检验过 |
| 164d | **范围** | 性别在本 release 是二值自报,且与其它未测的位置变量(社会经验、暴露史)混在一起。**"不是性别的影子"不等于"不是任何人口学变量的影子"** —— 只是性别是 `#101`/`#154` 唯一量到过外部锚的那个,所以它是唯一一个有具体理由怀疑的 |

**按 §0.2 结账**:本项目最老的现存声明通过了它最后一个未做的对照(**产出**);
一个"残差化是否有效**且**是否误伤"的双对照设计落地,两端都开火(**产出**);
这一轮**没有**产生新的撤回 —— 而在连续多轮都以撤回收尾之后,这本身值得记一笔,
但它不是成绩,只是这一次我事先把两个对照都造对了。

**NEXT**:`#160`–`#164` 把「什么」与「何时」两侧都钉死了,而它们合起来指向一个
**还没被问过的第三件事**:`#163b` 说位置倾向里**有一个小的内容成分(15%,5.7×)**,
而 `#159` 说「何时」的那点结构**几乎全在类别上**(题目侧 0.96 vs 人侧 0.105)。
**那么"内容"这一侧到底有多大?** 可测版本:把"这个人勾了哪些具体选项"(而不是它们多冷门)
做成一个人层量 —— 例如跨不相交块的**选项身份**预测(用一半块上的勾选去预测另一半块上的勾选,
控制住冷门程度)—— 在同一台机器上量它的分半信度。
**若它明显高于 0.60,那么本项目一直在测的那个"位置"只是更大的一个"内容"的投影。**

---

## Entry 165, added by `E01·A22·R210` — 位置不是内容的投影;而「有没有内容」这个问题被泄漏底噪吃掉了

`#164` 的 NEXT 问的是**还没被问过的第三件事**:一直在测的「位置」(偏爱多冷门),会不会只是一个
更大的「内容」(勾了哪些具体选项)的投影?可测版本:在每半块的**选项矩阵**上取主成分给人打分,
量它的分半信度,并把同一半上的位置分数回归掉。**预注册阈值 0.60** —— 明显高于它,框架就要重排。

**答案是不。** k ∈ {6,8} × 3 个劈分种子,整格:

| | 内容 SB | 去位置后 | 位置 SB |
|---|---:|---:|---:|
| 真实 | +0.278 | **+0.262** ± 0.110 | **+0.598** |
| 每块独立跨人置换 | +0.072 | +0.079 | +0.048 |
| 纯内容种植 | +0.799 | +0.672 | — |
| 纯位置种植 | +0.570 | +0.362 | +0.916 |

位置比内容高 0.335,是自身展布的 **3.1×**;预注册的 0.60 也在噪声带外(3.1×)。
三个正对照全开火:内容种植抬高内容分(0.278→0.799),位置种植抬高位置分(0.598→0.916),
置换零塌到 0.072。**`#100` 起就一直在测的那个东西,是这套数据里最稳的人层量,不是别人的影子。**

**`#165a` 这个零不该是 0,而我差一点就把它当成 0。** 正对照二(位置种植在去位置后测不到)第一次
只是**擦边过**:纯位置种植去位置后仍读到 +0.362,**比真实的 +0.262 还高**。原因不是位置里有内容 ——
是**回归掉一个有噪声的协变量会留下衰减残留**,而主成分会把残留读成「内容」。
所以造了一个构造上内容为零的合成世界(置换背景 + 纯位置种植),按**位置信度匹配**到真实的 0.60:

| 种植强度 | 位置 SB | 去位置后内容 |
|---:|---:|---:|
| 0.00 | 0.048 | 0.079 |
| **0.05** | **0.686** | **0.076** |
| 0.10 | 0.845 | 0.043 |
| 0.25 | 0.905 | 0.263 |
| 0.60 | 0.915 | 0.255 |

**泄漏是非线性的**:强度 0.6 处漏 0.26,匹配真实信度的 0.05 处只漏 0.076。所以刚才那个正对照
是在一个**偏离真实数据很远的点**上做的 —— 它读到的 0.362 里,大部分是种植强度自己的产物。
在正确的点上:真实 +0.262 − 底噪 +0.076 = **+0.187,1.6× 展布 → UNVERIFIED**。
**「去掉位置之后还剩不剩内容」这个问题,本设计答不了。** 不是零,也不是非零。

**`#165b` 一个正对照的强度,是它自己的一个规格。** 我一直把种植强度当成"够大就行"的旋钮 ——
选 0.6 是因为它保证正对照会开火。但一个在 0.6 上成立的对照,不能用来判一个在 0.05 上的真实效应,
因为**泄漏/衰减这类量对强度是非线性的**。教训写成一条:**正对照的强度必须匹配被判对象的可观测量级**
(这里是位置信度 0.60),否则那个对照测的是另一个世界。这与 `#152`(基准是一次随机实现)同族:
**判一个比较,先问它的参照点是不是站在跟被判量同一个位置上**。

**`#165c` 我又一次把种植写成了噪声。** 第一版 `u` 在两半里各抽一次 —— 种下去的不是人层信号,
是两半各自独立的噪声,只会**稀释**真信号。正对照因此**反向开火**(+0.31 → +0.085),
而如果我把它读成"内容测不到",结论会完全反过来。同一形状的错在 `#163` 刚犯过一次
(置换写成了整向量换人,零读到 71%)。**共同点:一个作用在"人"上的操作,必须由人的全局身份索引,
不能由局部/半边的抽样索引。** 顺带修:四个臂原本各自劈块,n 从 682 跳到 1346,臂之间根本不可比;
劈分现在只由 (k, 种子) 决定,与臂无关。

**`#165d` 两个分数正交与否,答不了。** |r| = 0.319 ± 0.180(共享方差 10.2%):既不可分辨地非零
(1.8×),也不可分辨地小(等价界只能界到 0.678,边界 0.30)。**按 frontier §3 的盆地规则,
同一个问题两个方向都 UNVERIFIED 就不追第三轮。**

**一句关于人的话**:一个人的性偏好里,**「你偏爱多冷门的东西」比「你具体偏爱哪些东西」更是一个稳定的
个人特征**(0.60 对 0.26)。挑什么内容,很大程度上是共享的、随人群走的;而**在冷门—常见这条轴上
你站在哪儿,是你自己的**。这不是说内容不重要 —— 是说内容里**属于个人的那一份**,在这套数据里
小到测不清楚。

**NEXT**:`#165` 把「内容」这一侧钉在 0.26,但也暴露了一个**方法层面的空洞**:本项目所有
「去掉 X 之后还剩多少」的说法(`#104` 去掉勾选数、`#164` 去掉性别、本轮去掉位置),
用的都是**同一个有噪声的回归残差**,而只有本轮量过它的泄漏底噪。
`#164` 的 102% 与 `#104` 的 67%,**分母里可能都带着同一种衰减残留**。
可测版本:对这三处各造一个**构造上目标成分为零、而被去掉的那个协变量信度匹配真实值**的合成世界,
读它们各自的泄漏底噪,把三个"保留率"重新算一遍。
**若 `#164` 的 102% 在底噪校正后掉到 80% 以下,那条声明要降级 —— 而它现在挂在 README 上。**
