# These 52 rounds, scored against the standard

Run under `realstat`. The standard is the 57-term campaign requirement, made failable by four
gates. **A standard nobody can fail is not a standard**, so this file leads with what fails.

## The one-line test

> Would this result probably have failed, had the claim been false?

Round by round, **31 of 52 yes, 14 partly, 7 no.** The seven are census/description rounds
(`r01`–`r06`, `r23`) which compute rather than test — correctly labelled DERIVATIONS here rather
than findings, per the arithmetic trap.

## The four gates

| gate | how these rounds do |
|---|---|
| **G1 estimand before method** | **Failed in 4 rounds, and it produced retraction #12.** `r37`–`r38` computed misfit and extremity and named them after; the estimand that mattered (breadth) was never stated, and matching deleted it. Also #2 and #3: the metric was chosen before asking what it could resolve |
| **G2 control saturation** | **Strongest area.** Positive controls in 14 rounds (`r13` instrument-not-blind, `r30` sex split, `r36` injected noise, `r43` out-of-sample, `r49` neuroticism reference). Permutation/degree-preserving/gating-preserving/binning nulls in 19. Matched placebo in `r16` (7/7 blocks). **MDEs stated in 6** (`r35`, `r36`, `r38`, `r46`, `r28`, `r33`) — every null that is load-bearing carries one |
| **G3 multiplicity over the whole grid** | **Weakest area. Done properly in 1 round of 52** (`r45`, family-wise permutation threshold over 34 columns). Elsewhere: uncorrected. The 321-pair and 418-pair sweeps report medians and counts above a floor rather than a corrected family, which is defensible for a floor comparison and is not a correction |
| **G4 specification curve** | **Done in 2** (`r47` three specifications published including the sign-flipping one; `r31` nine splits as a reference class). Everywhere else a single cell |

## Structurally impossible here, with what each would require

Not marked "planned" — that is an unavailability claim in the flattering direction.

| criterion | N/A because | what it would require |
|---|---|---|
| independently replicated | one public release, one team | a second release, or the triple-blind protocol below |
| causally identified · interventionally validated · counterfactually grounded | nothing was manipulated | Ivan's phase 2/3: meaning-flip and conditioning |
| temporally resolved | retrospective self-report in ~2-year bins | prospective collection |
| construct / criterion validated | no external gold standard for erotic value | psychophysiology, viewing time, behaviour |
| cross-model / cross-architecture | no model is used anywhere in this repository | — (this is a *strength* here: nothing routes through a judge) |
| cross-dataset · cross-domain · OOD | one dataset | SFBI, KOS, or a non-Western sample |
| position randomized · counterbalanced | no presentation-order field | a release that ships it |
| multi-seed · seed-robust | **not impossible — simply not done.** Most rounds run one seed | ≥3 seeds; this is the cheapest available upgrade and is scheduled |

## The one criterion this project meets unusually well

**Instrument-free.** Ten of eleven claims in the CoVal programme turned out to be a 2B model's
opinion stated as fact about a dataset, which is why `P14` has an `instrument` line. Nothing here
routes through a model. Every number is counted off the release or computed from it. The
`prior_art` line is a live risk instead: this dataset has a public explorer and a substantial
secondary literature, and **no round has yet checked whether a finding restates something already
published about it.** That is the highest-priority unrun check and it is `r53`.

## Honest score

Counting the 57 terms against what these rounds actually do: **19 met, 12 partly, 14 structurally
impossible with a stated requirement, 12 simply absent.** Met/(met+partly+absent) = **19/45 =
42%** — in the same band as the three projects audited on 2026-07-31 (6–49%), which is what an
honest detector should return. Near-100% would mean the detector is broken.

The 12 absent, in priority order: multi-seed · multiplicity over the grid · specification curve
as default · sham controls · dose-response beyond the two rounds that have it · hierarchical
modelling of the block structure · leakage audit of the option-text features · pre-registration ·
triple-blind replication · negative-control calibration in the group comparisons · noise floor by
replicate rather than by permutation in `r41`–`r44` · benchmark-degeneracy audit of the release.
