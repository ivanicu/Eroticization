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
