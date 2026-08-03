# The Eroticization Operator

**Is "sexual" a content category the brain detects, or a value the brain assigns to ordinary
world-representations?** Fifty-two self-attacking rounds against a public dataset, asking what
can be settled before anyone collects anything new.

Private curiosity project. Not written for publication and not seeking any.

Object: the **Big Kink Survey** public subsample — [Zenodo 10.5281/zenodo.18625141](https://zenodo.org/records/18625141),
15,503 respondents × 365 columns, aggressively binned, demographically stripped and noise-injected.
Data is not committed; `01_object_and_structure/r01_schema/run.py` re-downloads it.

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

Everything measured is consistent with one expression, and it is B:

```
E_i(s)  ~  theta_i  x  population-shared ordering   +   ~4 individual coordinates
           (scalar gain)                                (role and relation)
```

with a **third, orthogonal** organisation — *when* an interest arrives — that ignores the
coordinates entirely.

| what | number | scope |
|---|---:|---|
| shared cross-domain grammar, held-out CCA | **0.272** vs 0.055 floor | 321 block pairs, ≥600 common people |
| …surviving sex/age/personality/orientation removal | **0.200** vs 0.055 | 73% of the transfer is not demographic |
| coordinates surviving a block split-half | **4** practically (7 above 3× floor) | held-out canonical r .357 / .290 / .193 / .130 |
| correlation of those coordinates with the folk top/bottom axis | **≤0.159** | six coordinates, n=6,717 |
| three role axes' mutual correlation, disattenuated | **≤0.362** | effective dimensionality 2.95 of 3, n=3,890 |
| shared maturational schedule, per-person rank agreement | **+0.232** | 74.7% of 9,691 positive vs 48.9% null, d=0.69 |
| coordinate-similarity predicting acquisition timing | **ΔR² +0.0008, t=−0.46** | absent, not hidden (predictors correlate −0.028) |
| breadth's coordinate diversity vs a size-matched base-rate set | **−0.88%** | out-of-sample loadings, 10 null draws each, n=1,189 |
| breadth's nestedness | **24.2%** of chance→perfect | excess +0.0660 over a 0.7278 base-rate null |
| θ vs all 15 non-sexual life-history and personality variables | **R² = 0.012** | MDE at n=15,000 is \|rho\|=0.023 |
| consumption → coordinates, blocks+θ+sex matched | **0.0439** | vs 0.0121 for neuroticism, 0.0778 for sex |
| self-reported porn-induced fetish acquisition, timing | **0.5 yr, uniform** | 1 yr excluded; concentration null at 8.8× power |

---

## What each big round established (R = one iteration, one belief update, one commit)

**[R01 · the object](R01_object_and_structure)** — The item-level data is not the 68 category
ratings; it is 101 multiselect columns exploding to **1,332 options over 15,468 people**. Entry
to every block is gated on a parent rating: **P(enter | parent > 0) = 0.99**. This is
undocumented and it constrains every design downstream — naive cross-block holdout is
conditioned on liking the parent category. Person-profile split-half reliability **0.727**;
item base-rate reliability **0.999**. The noise injection did not destroy individual signal.

**[R02 · is there a shared grammar](R02_is_there_a_shared_grammar)** — Yes, and the first two
instruments could not see it. PC1-vs-PC1 gave 0.064 (retraction #2); held-out CCA gave 0.272 vs a
0.055 floor; leave-one-block-out gave **32/32 blocks positive**, median gain +0.0340 against a
−0.0029 floor. Factors learned from 31 domains predict a domain they never saw.

**[R03 · naming the coordinates](R03_naming_the_coordinates)** — Four coordinates survive a
**block** split-half (not a person split-half — the question is whether a coordinate is
recoverable from either half of the *domains*). Naming failed twice before succeeding. The
surviving three: *light restraint/toys ↔ insertive extremity and confinement* · *abjection and
filth ↔ blood, burning, weapons* (two **opposed** extremities, not one intensity axis) ·
*receiving a substance ↔ giving it*. None is the folk axis. **"Top/bottom" is one word over
three near-independent coordinates**: who submits, who is seen, who receives — disattenuated
mutual r ≤ 0.362, effective dimensionality 2.95 of 3.

**[R04 · acquisition and time](R04_acquisition_and_time)** — Interests arrive on a
population-shared schedule (content-like early: appearance 14.0, body parts 14.4, clothing 14.7;
relational late: power dynamics 16.8, bondage 16.9, mental alteration 17.0). Within-person
"acquired together" tracks "liked together" at RSA **+0.599** after stripping intensity leakage
(a real confound, −0.126) and recall anchoring, and after near-synonymy is excluded (0.594 for
280 pairs sharing no content word vs 0.646 for 64 that do). **80% of the sd of onset structure is
not explained by preference** — interests a person does *not* like together were acquired at the
same time. But the organising variable is arrival time, not coordinate membership (retraction #9).

**[R05 · group differences and the instrument](R05_group_differences_and_the_instrument)** —
**Any group comparison on this release must be block-count matched or it partly measures survey
coverage**: corr(congruence deficit, coverage gap) = **+0.815** across nine splits. This is the
transferable methodological result and it is not documented anywhere I have read about this
dataset. Drawn-vs-live-action consumers — whose content contains no real bodies — differ by
**0.0204 ± 0.0265**, bounded below sex-sized (0.093). That is a real, bounded constraint on A.

**[R06 · induction](R06_induction)** — 82.7% of the 13,530 respondents with fetishes say porn
induced ones that would not otherwise exist. That claim carries **no acquisition timing
signature** (a uniform 0.5-year shift of *all* their onsets; a 1-year shift is excluded; the
concentration discriminator is null at 8.8× the power needed for a single-interest induction)
and **no structural signature** (misfit < 0.1 sd, extremity < 10%). It tracks **breadth**
(rho +0.2515, 85% surviving response-style control). Design consequence: retrospective
self-report of induction is unusable as an outcome measure. Phase 3 must be prospective.

**[R07 · breadth](R07_breadth)** — Quantity without shape. A person's set is **0.88%** more
concentrated in coordinate space than a size-matched base-rate set — real (t=−15.5) and
negligible. Sets are **24.2%** of the way from chance to perfect nesting. Breadth is moderately
one trait (Spearman–Brown 0.557) and **9–13% response style**, measured on the survey's own
non-category Likert items.

**[R08 · what theta is](R08_what_theta_is)** — Nearly domain-encapsulated. All 15 non-sexual
variables jointly: **R² = 0.012**. Childhood adversity, adult sexual assault, corporal
punishment and sexual repressiveness of upbringing are all bounded under **0.09 disattenuated**,
at a sample size where 0.023 would have been visible. The family of theories rooting broad or
unusual sexuality in adversity gets essentially no support here.

**[R09 · consumption](R09_consumption)** — Consumption touches both terms: it correlates with θ
(rho 0.17) *and* independently with the coordinates (**0.0439** after matching coverage, breadth
and sex — 3.6× the neuroticism reference). Direction is unavailable; this is exactly what a
prospective design would resolve.

**[R10 · additivity](R10_additivity)** — Ivan's model B requires it. The feature crosses every
substance boundary it meets (+0.24 to +0.58) **except source gender, where it is +0.017**.
Additivity is basis-dependent: in the folk basis (self/other) it fails; in a basis crossing role
with source gender it may hold. Small n on the decisive cells (3–4 pairs) — **PLAUSIBLE, not
CONFIRMED**.

---

## Governance

| file | what it holds |
|---|---|
| [`RETRACTIONS.md`](RETRACTIONS.md) | **14 entries.** Every claim killed, scoped or corrected, with what killed it. Twelve are a later round of mine destroying an earlier round of mine |
| [`FROZEN.md`](FROZEN.md) | Lines where further computation cannot identify what it is measuring, each with its unfreeze condition |
| [`PREREGISTRATION.md`](PREREGISTRATION.md) | What the next rounds will test, with thresholds fixed in advance. Rounds r01–r52 were **not** preregistered and are labelled exploratory throughout |
| [`ADVERSARY_FORECAST.md`](ADVERSARY_FORECAST.md) | What I expect an outside challenger to overturn, written before one arrives |
| [`STANDARD_AUDIT.md`](STANDARD_AUDIT.md) | These 52 rounds scored against the campaign standard, including everything they fail |

## Layout

Ten campaigns, 52 rounds. Each round is a directory with `run.py`, `README.md` and `results/`.
`lib/rounds.py` maps a round name to its path, because several rounds reuse an earlier round's
loaders and that dependency is made explicit rather than hidden. Environment is a self-contained
`.venv` (system python 3.14 has no `ensurepip`; pip was bootstrapped).

## Scope, stated once

Cross-sectional · one instrument · one population (18–32, US/Canada/Europe) · aggressively
anonymised, with correlations attenuated roughly 25% by design · every measure self-report.
No causal claim in this repository is identified, and none is made.
