# Preregistration

**Rounds r01–r52 were not preregistered.** They are exploratory, and every claim in
[`README.md`](README.md) carries that status. Retrofitting a preregistration onto completed work
is the thing this file exists to prevent, so nothing below refers to a round that has run.

Thresholds are fixed here before the data is touched. A threshold chosen after seeing a result is
a narrative.

---

## r53 · Prior art, and it runs before anything else

**Why first.** The CoVal programme spent 140 rounds discovering an algorithm documented in its own
dataset card. This release has a public explorer, a Zenodo record and a secondary literature, and
**no round here has checked whether any finding restates something already published.**

**Procedure.** For each of the eleven headline claims in `README.md`: search the release
documentation, the explorer's own analysis pages, and the literature. Code each claim
`NOVEL / VERIFICATION / SUPERSEDED`.

**Pre-registered consequence.** Any claim coded VERIFICATION is relabelled in `README.md` as a
verification that the object does what it says — **not** as a finding. No exceptions, including
for the three I most want to keep (the three role axes, the +0.815 coverage artifact, the
theta encapsulation).

## r54 · Multi-seed, because it is the cheapest thing absent

**Claim under test.** The eleven headline numbers are seed-stable.

**Procedure.** Re-run every round contributing a headline number at ≥5 seeds. Report
`seed_spread / |effect|` for each.

**Pre-registered kill.** Any headline whose `seed_spread / |effect| > 0.5` is downgraded to a
direction and loses its number in `README.md`.

## r55 · The norm manipulation, the one separator this release cannot provide

**Claim under test.** Retraction #14 killed `n_i` as a *trait*. It did not test the prediction the
parameter was invented for: that normalising a feature moves erotic value in **opposite
directions** for transgression-driven and feature-driven people.

**Why it is the highest-leverage new collection.** A manipulation where two subgroups move in
opposite signs is far more diagnostic than any mean shift, and it is unfakeable by demand
characteristics — a participant guessing the socially correct answer produces a uniform sign.

**Design.** Same physical scenario; the only manipulated variable is stated prevalence
("most people find this ordinary" vs "this is rare"). Outcome measured as arousal / wanting /
pleasure separately, never as one score. Between-subjects on the framing, within on the scenario.

**Pre-registered thresholds.** Support requires (a) a subgroup × framing interaction with the two
subgroup slopes of opposite sign, and (b) the interaction surviving the split defined on a
*held-out* half of each participant's scenarios. Fixed now: interaction p < .01 and both simple
slopes' 95% CIs excluding zero **in opposite directions**. A main effect with no sign split is
recorded as a refutation, not as partial support.

## r56 · Triple-blind reanalysis of the shared-grammar claim

**Claim under test.** The 0.200 demographic-adjusted cross-domain CCA.

**Procedure.** Two clean-context agents receive the data, the schema, the question in plain words,
this file, distinct seeds, and a list of files they may not open. They receive **no** estimand, no
statistic, no controls, no numbers. Designing the statistic is the task.

**Reading rule, fixed now.** All three agree → design-independent. Agree on sign, differ on size →
report the spread as the finding. **Disagree on sign → the framing is the finding**, and the
assumption they differ on becomes r57. The three will not be averaged.
