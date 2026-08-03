# The Eroticization Operator

**Is "sexual" a content category the brain detects, or a value the brain assigns to ordinary
world-representations?** Fifty-two self-attacking rounds against a public dataset, asking what
can be settled before anyone collects anything new.

Private curiosity project. Not written for publication and not seeking any.

Object: the **Big Kink Survey** public subsample — [Zenodo 10.5281/zenodo.18625141](https://zenodo.org/records/18625141),
15,503 respondents × 365 columns, aggressively binned, demographically stripped and noise-injected.
Data is not committed; `01_object_and_structure/R01_schema/run.py` re-downloads it.

---

## The three models being separated

Ivan's framing, kept verbatim because the whole repository is scoped by it:

| | |
|---|---|
| **A** dedicated sexual-content system | `Y_i(s) = f_i(h_sex(s))` — "sexual" is a content category like *face* in vision |
| **B** erotic valuation of ordinary representation | `v_i(s,c,t) = w_i(c,t)^T h(s)` — ordinary scene, a readout applied to it |
| **C** recursive mix | value re-enters and reshapes the representation it acted on |

---

## Headline

**Every surviving claim is sorted by whether it carries predictable variance ([RETRACTIONS #51](RETRACTIONS.md)):**
onset RSA **31%** of its ceiling · nestedness **25%** · sex deficit **9%** · cross-domain CCA **~0%** (a direction, not an amount) · maturational schedule **101%** — at the ceiling a global ordering can reach. **[#63]** accuracy is **[60.5%, 66.5%]**: 66.5% on the 63.7% of pairs the 2-year binning can order, 60.5% over all pairs with ties at half credit.


Everything measured is consistent with one expression, and it is B:

```
E_i(s)  ~  theta_i  x  population-shared ordering   +   a THIN shared subspace
           (scalar gain)                                (real, graded, predicts ~nothing per pair)
```

**[RESTATED — RETRACTIONS #49]** The shared cross-domain direction exists (CCA 0.198, tier-ordered)
but carries almost no predictable variance: pairwise block→block prediction R² is −0.002. It becomes
predictive only when aggregated over ~31 blocks (+0.017–0.037). The coordinate count is withdrawn (#18).

with a **third, orthogonal** organisation — *when* an interest arrives — that ignores the
coordinates entirely.

| what | number | scope |
|---|---:|---|
| shared cross-domain grammar, held-out CCA | **0.269 ± 0.014** vs 0.051 floor | 306 block pairs, ≥600 common people, 6 seeds. **[RESOLVABLE — RETRACTIONS #37]** ratio 19.6 |
| …surviving sex/age/personality/orientation removal | **0.200** | pairwise-block CCA, max over 3 components, median of 321 pairs. **[#35 withdrawn by #36]** the 0.149 figure was a different estimand (half-split factor CCA), not a discrepancy |
| coordinates surviving a block split-half | **WITHDRAWN — see RETRACTIONS #18** | count tracks K (4→4, 24→16, non-saturating); what survives is the profile's height above floor, 30–50× |
| correlation of those coordinates with the folk top/bottom axis | **≤0.159** | six coordinates, n=6,717 |
| POWER vs SUBSTANCE, ladder-disattenuated | **+0.605** ⚠ **[SCOPED — RETRACTIONS #39]** SUBSTANCE predicts nothing outside its 7 defining blocks (−0.0001 vs a role-free placebo's +0.0024) | **[RESTATED — RETRACTIONS #24]** 37% shared variance, not 5%; positive control recovers 1.018 for a measure against itself; GAZE unmeasurable here; "2.95 of 3" withdrawn |
| shared maturational schedule, share of people following it | **[0.747, 0.860]** | **[CONFIRMED — RETRACTIONS #31]** 72-cell specification curve; null stays 0.479–0.506 in every cell; positive control 0.993. The published 74.7% was the most conservative cell |
| coordinate-similarity predicting acquisition timing | **t=−0.46 → +1.26** | **[CONFIRMED — RETRACTIONS #43]** null survives rebuilding COORD from the validated person factors, at 344 pairs vs the original 153; TEMPO reproduces at t=+4.36 |
| breadth's coordinate diversity vs a size-matched base-rate set | **UNVERIFIED** | **[RETRACTIONS #46]** the measure cannot separate a synthetically concentrated population (−1.67%) from a base-rate one (−1.62%), and −0.88% is inside its own Jensen bias |
| breadth's nestedness | **24.0%** of chance→perfect | **[CONFIRMED — RETRACTIONS #47]** validated with known-answer populations: nested synthetic recovers 100%, base-rate synthetic −0.0001 |
| θ vs non-sexual variables, reliability-corrected | **r_true ≈ 0.10** | **[RESTATED ×3 — RETRACTIONS #17, #23]** adversity 0.105, mental illness 0.099, openness 0.090; null controls 0.000/0.010. ~1% of variance each |
| consumption → coordinates | **UNRESOLVABLE at triple matching** (0.036, ratio 1.40); reportable as 0.072–0.090 at looser matching. **[RETRACTIONS #45]** | **[RESTATED — RETRACTIONS #19]** modality is larger (+0.0546); both clear a known-null pool (p95 0.0216); the neuroticism reference was itself null |
| self-reported porn-induced fetish acquisition, timing | **0.5 yr, uniform** | 1 yr excluded; concentration null at 8.8× power |

---

## What each big round established (R = one iteration, one belief update, one commit)

**[R01 · the object](A01_object_and_structure)** — The item-level data is not the 68 category
ratings; it is 101 multiselect columns exploding to **1,332 options over 15,468 people**. Entry
to every block is gated on a parent rating: **P(enter | parent > 0) = 0.99**. This is
undocumented and it constrains every design downstream — naive cross-block holdout is
conditioned on liking the parent category. Person-profile split-half reliability **0.727**;
item base-rate reliability **0.999**. The noise injection did not destroy individual signal.

**[R02 · is there a shared grammar](A02_is_there_a_shared_grammar)** — Yes, and the first two
instruments could not see it. PC1-vs-PC1 gave 0.064 (retraction #2); held-out CCA gave 0.272 vs a
0.055 floor; leave-one-block-out gave **32/32 blocks positive**, median gain +0.0340 against a −0.0029 floor. **[BOUNDED — RETRACTIONS #25]** with all 32 coverage indicators in the baseline the gain is +0.0170 and 31/32 blocks stay positive; coverage is partly a mediator, so the honest statement is a bound of **[+0.017, +0.037]**. Factors learned from 31 domains predict a domain they never saw.

**[R03 · naming the coordinates](A03_naming_the_coordinates)** — Four coordinates survive a
**block** split-half (not a person split-half — the question is whether a coordinate is
recoverable from either half of the *domains*). Naming failed twice before succeeding. The
surviving three: *light restraint/toys ↔ insertive extremity and confinement* · *abjection and
filth ↔ blood, burning, weapons* (two **opposed** extremities, not one intensity axis) ·
*receiving a substance ↔ giving it*. None is the folk axis. **"Top/bottom" is one word over
three near-independent coordinates**: who submits, who is seen, who receives — disattenuated
mutual r ≤ 0.362, effective dimensionality 2.95 of 3.

**[R04 · acquisition and time](A04_acquisition_and_time)** — Interests arrive on a
population-shared schedule (content-like early: appearance 14.0, body parts 14.4, clothing 14.7;
relational late: power dynamics 16.8, bondage 16.9, mental alteration 17.0). **[SURVIVES FRAMING SWAP — RETRACTIONS #50]** Within-person "acquired together" tracks "liked together" at RSA **+0.599**, and unlike the cross-domain grammar (#49) it carries predictable variance: person-level onset→preference R² = +0.0136 against a −0.0043 null, 31% of the same-domain ceiling after stripping intensity leakage
(−0.126 — **[VERIFICATION — see RETRACTIONS #16]**, published as `01-age-onset` Finding 1) and recall anchoring, and after near-synonymy is excluded (0.594 for
280 pairs sharing no content word vs 0.646 for 64 that do). **[CONFIRMED — RETRACTIONS #22]** onset carries structure preference does not: observed top residual eigenvalue 0.959 against a purpose-built rival world at 0.441 ± 0.040 (95% upper 0.532), positive control fires at injection 0.3 and stays silent at 0.0. The attached '80% of the sd' is descriptive, not tested — interests a person does *not* like together were acquired at the
same time. But the organising variable is arrival time, not coordinate membership (retraction #9).

**[R05 · group differences and the instrument](A05_group_differences_and_the_instrument)** —
**Any group comparison on this release must be block-count matched or it partly measures survey
coverage**: corr(congruence deficit, coverage gap) = **+0.815** across nine splits. This is the
transferable methodological result and it is not documented anywhere I have read about this
dataset. Drawn-vs-live-action consumers — whose content contains no real bodies — differ by **0.0204 ± 0.0265**, bounded below sex-sized (0.093). **[UNVERIFIED — RETRACTIONS #34]** the written-vs-visual deficit is unresolvable at 5 seeds (0 of 9 corpus cuts have effect > 2× seed spread, while sex is resolvable in 9/9); ~44 seeds per cell would be needed. That is a real, bounded constraint on A.

**[R06 · induction](A06_induction)** — 82.7% of the 13,530 respondents with fetishes say porn
induced ones that would not otherwise exist. That claim carries **no acquisition timing
signature** (a uniform 0.5-year shift of *all* their onsets; a 1-year shift is excluded; the
concentration discriminator is null at 8.8× the power needed for a single-interest induction)
and **no structural signature** (misfit < 0.1 sd, extremity < 10%). It tracks **breadth**
(rho +0.2515, 85% surviving response-style control). Design consequence: retrospective
self-report of induction is unusable as an outcome measure. Phase 3 must be prospective.

**[R07 · breadth](A07_breadth)** — Quantity without shape. A person's set is **0.88%** more
concentrated in coordinate space than a size-matched base-rate set — real (t=−15.5) and
negligible. Sets are **24.2%** of the way from chance to perfect nesting. Breadth is moderately one trait (Spearman–Brown 0.557). **[UNVERIFIED — RETRACTIONS #26]** the '9–13% response style' figure comes from an index that is orthogonal to the POWER axis (+0.024) but correlates +0.385 with breadth itself; response style and erotic endorsement cannot be separated without balanced-keyed items, which this release lacks.

**[E01·R05 · is breadth the object](E01_sexual_as_a_value_not_a_category/A05_is_breadth_a_nuisance_or_the_object)** — Nearly domain-encapsulated. All 15 non-sexual
variables jointly: **R² = 0.012**. **[PARTLY RETRACTED — see RETRACTIONS #17]** That held only with survey progression controlled and a *gated* outcome. On the ungated `totalfetishcategory` with acquiescence controlled, childhood adversity → breadth is **r = +0.059, effect/floor 7.5** — real and small, matching the published d=0.151. Adult sexual assault, corporal punishment and upbringing remain near zero.

**[R09 · consumption](A09_consumption)** — Consumption touches both terms: it correlates with θ
(rho 0.17) *and* independently with the coordinates (**0.0439** after matching coverage, breadth
and sex — 3.6× the neuroticism reference). Direction is unavailable; this is exactly what a
prospective design would resolve.

**[R10 · additivity](A10_additivity)** — Ivan's model B requires it. The feature crosses every
substance boundary it meets (+0.24 to +0.58) **except source gender, where it is +0.017**.
Additivity is basis-dependent: in the folk basis (self/other) it fails; in a basis crossing role
with source gender it may hold. Small n on the decisive cells (3–4 pairs) — **PLAUSIBLE, not
CONFIRMED**.

---

## How the shared signal accumulates

**[RETRACTIONS #60]** `increment = 0.00723 × √(source domains)` — CV **6.4%** across n=1…31, against 11.2%
for log and 51.9% for linear, with a permuted-label null of exactly 0.0000 at every n. **[MECHANISM UNSUPPORTED — #61]** blocks are *not* interchangeable (subset variance 58× seed variance; block
size predicts subset quality at r=+0.816), so √(block count) and √(total sample) are confounded in that fit.
The curve stands as a description; what produces it is untested. Price list for new collection:
0.06 needs 69 domains, 0.08 needs 122, 0.10 needs 191 — extrapolated beyond the measured range.

## Within-block vs cross-domain

**[RETRACTIONS #59]** Two largely independent contributions, 87% additive, both against two exactly-null
controls: **cross-block factors +0.0409** (fitted on 31 *other* blocks) and **within-block structure
+0.0290**, combined **+0.0606**. The aggregate cross-domain signal is the *larger* — pairwise block→block
is ~0 (#49), and it accumulates across blocks.

## The item margin

**[RETRACTIONS #57]** Fifty-seven rounds measured structure over PEOPLE. The ITEM side, measured once:
item-neighbour structure recovers **+0.0206** held-out R² over a marginals base, against the person
factors' **+0.0289** and a random-neighbour control's **+0.0006** — **71% of the person margin**, on
578,989 held-out cells. **[#58]** Fitted jointly, the two are **67% the same structure**: combined increment +0.0359 against the person margin's +0.0293, so the item side's unique contribution is +0.0066 — ratio 1.65, **below resolution**.

## Governance

| file | what it holds |
|---|---|
| [`RETRACTIONS.md`](RETRACTIONS.md) | **14 entries.** Every claim killed, scoped or corrected, with what killed it. Twelve are a later round of mine destroying an earlier round of mine |
| [`FROZEN.md`](FROZEN.md) | Lines where further computation cannot identify what it is measuring, each with its unfreeze condition |
| [`PREREGISTRATION.md`](PREREGISTRATION.md) | What the next rounds will test, with thresholds fixed in advance. Rounds r01–r52 were **not** preregistered and are labelled exploratory throughout |
| [`ADVERSARY_FORECAST.md`](ADVERSARY_FORECAST.md) | What I expect an outside challenger to overturn, written before one arrives |
| [`STANDARD_AUDIT.md`](STANDARD_AUDIT.md) | These 52 rounds scored against the campaign standard, including everything they fail |

## Layout

**E · R · r** — epoch (the ontology shifted) contains big rounds (a decision became safe) contains
sub-rounds (one belief update). One epoch, six big rounds, 54 sub-rounds. Every count is discovered,
never chosen; see `~/.claude/CLAUDE.md` §P16.

Previously described as ten campaigns, 52 rounds. Each round is a directory with `run.py`, `README.md` and `results/`.
`lib/rounds.py` maps a round name to its path, because several rounds reuse an earlier round's
loaders and that dependency is made explicit rather than hidden. Environment is a self-contained
`.venv` (system python 3.14 has no `ensurepip`; pip was bootstrapped).

## Scope, stated once

Cross-sectional · one instrument · one population (18–32, US/Canada/Europe) · aggressively
anonymised, with correlations attenuated roughly 25% by design · every measure self-report.
No causal claim in this repository is identified, and none is made.
