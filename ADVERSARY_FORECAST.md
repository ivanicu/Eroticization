# Adversary forecast

No outside challenger has run against this repository. Every one of the fourteen retractions is
self-inflicted, which is exactly the condition under which self-review is void rather than weak —
a reviewer sampled from the weights that produced the material can only attack what was already
anticipated.

So this file is written **before** a challenger arrives, and its purpose is to be scored. When one
runs, what it finds that is not below is worth more than any individual verdict, because it
measures calibration about my own blind spots.

## What I expect to be overturned, ranked by my own probability

| p | prediction |
|---|---|
| 0.75 | **At least one headline is prior art.** The dataset has a public explorer with its own analyses. The three role axes are the most likely — "dominance, exhibitionism and receptivity are separable" is the kind of thing a factor analysis of this data would already have produced |
| 0.65 | **The four-coordinate count is unstable to the option floor.** `n≥20` was chosen once, in `r17`, and never swept. A challenger sweeping it finds 3 or 6 |
| 0.55 | **The +0.815 coverage correlation is inflated by two points.** Nine splits, and PORNHABIT sits far out on both axes. Without it the correlation is much weaker, and the *methodological* claim rests on a single influential observation |
| 0.50 | **`r24`–`r26`'s RSA survives its controls but not a hierarchical model.** Category pairs are not independent units; 344 pairs come from 27 categories, so the effective n is closer to 27 than to 344 and every p-value there is overstated |
| 0.40 | **The disattenuation in `r22` is doing too much work.** GAZE's alpha is 0.163; dividing by sqrt(0.163) is numerically unstable and I said so, but the three-axis claim still leans on it |
| 0.35 | **Breadth's 0.557 reliability makes θ's "no external correlate" partly a reliability ceiling**, and a better-measured θ would show the adversity correlations at 0.15 rather than under 0.09 |
| 0.25 | **`r16`'s 7/7 placebo result is one template, not seven replicates.** I said so in the commit; a challenger will say it louder and be right that the fluid family is n=1 in template |

## What I expect to survive

The leave-one-block-out result (`r13`), because it is the only claim here tested by prediction into
held-out *domains* rather than held-out people. And the gated-tree structure itself (`r04`), which
is a fact about the release that anyone can reproduce in one query.

## What I have no prediction about

Whether the whole world-decomposition is wrong — whether A, B and C are three answers to a
question that does not carve. That is the meta-separator and I cannot forecast my own blind spot;
naming it here is the only honest move available.

---

# Second forecast block — the `A09` / `A10` claims, written before any challenger runs

Added after `#75`. The block above predates arcs `A09` and `A10` entirely: it was written about a
world in which the item main effect had never been measured, so scoring `#1`, `#5`, `#7` now would
grade them against a question they were not asked. They stay open and unscored.

These are the new claims, and the same rule applies — **what a challenger finds that is not on this
list is worth more than any verdict on the ones that are.**

## What I expect to be overturned, ranked by my own probability

| p | prediction |
|---|---|
| 0.70 | **The cross-block transfer `C` is partly the gate, not the person.** Entry to every block is conditioned on the parent rating (`P(enter\|parent>0)=0.99`, `#04`). People present in many blocks are people who rated many parents highly — so "the same person deviates the same way across blocks" may be partly "the same person passed the same gate". The person-permutation null does **not** control for this: it destroys the person correspondence but every permuted person is still a gate-passer. The right control is a **gate-matched permutation** — permute only within strata of identical block-entry patterns — and I have not run it |
| 0.60 | **`#69`'s `R² = 0.994` is less decisive than it reads.** Across 23 blocks, `I` and `I_hat` both span 0.09–0.33, and any two quantities with that much shared range will correlate near 1. A challenger regresses the *residual* on something meaningful and finds the 0.006 that is left is not noise |
| 0.55 | **The item/interaction comparison is scale-dependent in a way I did not test.** Everything is squared error on binary cells. Under log-loss the interaction's share rises, because squared error on a 0/1 outcome under-weights confident-and-wrong predictions, which is exactly where a person-specific readout would differ from a base rate |
| 0.50 | **`#71`'s parameter-count argument is a bad accounting.** I called `C`'s person scores "free" because they are estimated elsewhere — but they are estimated on *the same people*, from 31 blocks, and that estimation is not free, it is just not billed to the target. A challenger charges it and the 22–800× advantage shrinks or inverts |
| 0.45 | **`#72`'s "no cliff" is a property of the estimator's regularisation, not the spectrum.** Soft-impute at rank K with 25 iterations has its own implicit shrinkage; the known-rank controls were generated with Gaussian factors and clean structure, which is the easiest possible case for finding a cliff. A challenger builds a known-rank world with **heavy-tailed loadings and unequal block sizes** and finds the cliff vanishes there too |
| 0.40 | **`#73`'s +0.90 is neighbour-set size, not similarity.** The random-neighbour control uses one random set per seed, so its variance across draws is under-sampled relative to the neighbour arm's. A challenger draws a fresh random set per *person* and the gap narrows |
| 0.30 | **The 23-block "identified" population is a selected sample.** Nine blocks were dropped for `\|dI\| > 0.01`, which correlates with block size — so every `A09`/`A10` number is conditioned on larger, better-measured blocks and does not generalise to the ones that were dropped |

## What I expect to survive

**`#69`'s gauge test**, because it is the one result here with a shape that is hard to fake: a
dispersion-matched surrogate reproduces `I` in **23 of 23** blocks. Whatever the item effect is, it
demonstrably cannot see *which* options are popular, and that alone kills the content reading
regardless of what happens to the `R² = 0.994`.

And **`#74`**, because it is a measurement of my own code against itself and either the cap bound or
it did not.

## What I have no prediction about

Whether the two-question split in `#70` is the right decomposition at all, or whether a challenger
will say both halves are asking something a survey of endorsements cannot address in either
direction. That is the meta-separator, and I cannot forecast my own ontology from inside it.
