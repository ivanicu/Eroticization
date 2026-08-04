# The Eroticization Operator

**"Sexual" — is it a content category the brain detects, or a value the brain assigns to ordinary
world-representations?** 206 self-attacking rounds against one public dataset, asking what can be
settled before anyone collects anything new. **164 numbered ledger entries** — [RETRACTIONS.md](RETRACTIONS.md)
carries the history; this page carries only what currently stands. **中文版:[`README_zh.md`](README_zh.md)** — English is canonical; corrections land here first.

Private curiosity project. Not written for publication and not seeking any.

Object: the **Big Kink Survey** public subsample — [Zenodo 10.5281/zenodo.18625141](https://zenodo.org/records/18625141),
15,503 respondents × 365 columns, aggressively binned, demographically stripped, noise-injected.
Data is not committed — neither the raw release nor the person-level derived tables; `E01…/A01…/R001_schema/run.py` re-downloads and rebuilds them. **[CORRECTED before publication]** 24 derived `data/derived/*.csv` files, two of them with one row per respondent, were tracked until now, so this sentence was false; they are untracked and gitignored as of this commit. They remain in the commit history.

## The three models being separated

Ivan's framing, kept verbatim because the whole repository is scoped by it:

| | |
|---|---|
| **A** dedicated sexual-content system | `Y_i(s) = f_i(h_sex(s))` — "sexual" is a content category like *face* in vision |
| **B** erotic valuation of ordinary representation | `v_i(s,c,t) = w_i(c,t)^T h(s)` — ordinary scene, a readout applied to it |
| **C** recursive mix | value re-enters and reshapes the representation it acted on |

---

---

# What this says about people

Each row is one sentence about people, at the strength it was actually measured. The ledger entry is
the receipt.

| | | receipt |
|---|---|---|
| **The "what" of a sexual preference is personal. The "when" is not.** | Same people, same split procedure, same number of units *k*: **"how uncommon the things you like are" split-half reliability 0.57–0.67**, **"how much earlier than the shared schedule your uncommon interests arrive" only 0.055–0.108** — a **5.8–10.4× gap** (mean difference +0.5321, **26.6×**, sign consistent across four values of *k*). **And this is not an instrument limit**: plant a person-specific signal into the "when" and the same instrument recovers it at **0.67**. **A second, independent way of asking "when" gives the same small number** (`#161`: "how well the shared schedule fits you" = 0.163, indistinguishable from "how far you deviate" = 0.114 at 1.3× spread; both **12.2×** below the "what") | `#159` `#160` `#161` |
| **Sexual interests arrive on a shared developmental schedule — concrete first, relational 2–3 years later** | appearance 14.0 · body parts 14.4 · clothing 14.7 → power dynamics 16.8 · bondage 16.9 · mental alteration 17.0. Using it to guess which of any two interests came first reaches **66.852 ± 0.191**, against a theoretical ceiling of 66.5% for "everyone shares one ordering" — the schedule has nearly exhausted what a shared ordering can provide | `#75` 8 seeds · 35,438 pairs · 6,230 people; all-pairs bound **[60.5%, 66.5%]** `#63` |
| **Content and individualised valuation move endorsement probability by the same amount** | option base rate **±21.5 pp** (ledger records ±22.6; **a 5% drift that is not explained**, `#157e`) · person overall rate **±16.7 pp** · person × option interaction **±23.7 pp** (`#157` moved two corrections that had lived only in ledger prose into the code; the interaction reproduces exactly). Content dominates every *predictability* measure (3.5–5.6×) **purely because it is estimated 179× better** (3,228 observations vs 18) | `#88` `#90` |
| **"Is this sexual content?" cannot be asked of this survey** | In a binary endorsement matrix an option's **base rate** and its **content** are the same number. The item main effect has a closed form: measured R² = **0.994**, slope 0.988, and completely insensitive to the *shape* of the base-rate distribution in 23 of 23 blocks | `#69` |
| **Taste for uncommon things is a continuum everyone is on, not a special minority** | The distribution widens symmetrically: the upper tail picks rarer options, the lower tail picks commoner ones, and the **median does not move**, each side 11–13× its own bootstrap sd | `#99` |
| **It is a reliable personality dimension, and it is not "how many boxes you ticked"** | Split-half across disjoint block sets **+0.432 ± 0.016** (person bootstrap, n=7,316, 400 draws — an error bar it went eight rounds without, `#167`; the −0.022 it used to carry is a curveball *null floor*, which answers a different question) against a planted ceiling of +0.832. **And the name itself is now tested** (`#163`): split the blocks into their two most content-dissimilar halves and **85% of the reliability survives**, while a pure *position* plant survives at 100% and a pure *content* plant collapses to −2% — **"prefers uncommon things" is a transferable position, not "likes one family that happens to be uncommon"**. **And it is not a shadow of sex** (`#164`): regressing sex out of *every* block-level score leaves the tendency at **102%** (4.9× its own spread, and the direction is *up*), while a purely sex-driven position plant collapses to **−0%** and a purely personal one survives at **+0.986** under the same treatment; **67%** survives removing the pick count, and its correlation with pick count (+0.608) is **lower** than the null's (+0.719). **And it is not a projection of a bigger "content" dimension** (`#165`): scoring people on the *option matrix* itself — which specific things they ticked — gives a split-half of only **+0.26**, against position's **+0.60** (3.1× its own spread, and the pre-registered 0.60 threshold sits outside the noise band at 3.1×). **Which content you prefer is largely shared; where you sit on the uncommon–common axis is yours** | `#100` `#104` `#165` |
| **A person's own uncommon interests sit ahead of their own common ones** | After removing the population-shared schedule and the person's overall precocity, the within-person «onset residual × category rarity» = **−0.0328** against a within-item cross-person permutation null of +0.0018, **8.8×**. **A discrete version (earliest bin vs the person's own repertoire) reproduces it independently: +0.0767 with the schedule removed, 14.5×**. `#114`'s recall bias contributes only **4%**, and with the opposite sign | `#128` `#130` |
| **The map's strongest psychologically-anchored dividing line is not "early vs late" — it is "objects vs scenarios"** | A 6 coherent splits × 11 non-sexual variables grid, with a max-statistic null giving a family-wise threshold of 0.0560. **PC4 lights up 6 variables, the strongest cell at 9.0×**: one side is clothing · genderplay · secretions · roles · dirtiness · transformation · toys, the other is age · horror · nonconsent · creatures · incest · gentleness · brutality · pregnancy · mythical. **Five anchors survive within-sex** (openness +0.083 · extroversion +0.088 · neuroticism −0.064 · powerlessness −0.059, both sexes agreeing in direction). ⚠ The naming is my reading of the loadings (D5); the scale is small (0.033–0.088). ⚠ **The model-level direction is withdrawn**: a naming-independent structural test (`#139`) shows PC4's two sides are **no more non-isomorphic than any other split** (0.8×), so this line is compatible with both model A and model B | `#138` `#139` |
| **Sex anchors "the sexual" itself, not one family within it** | early family +0.0706 vs late family +0.0816, difference **0.7×**. Openness (+0.0593, 4.1×) and powerlessness (−0.0441, 2.8×) really do anchor only the late family — **but that difference is not specific to these two families**: five other equally coherent splits of the same matrix give 0.0702–0.2106, so **the early/late line is the weakest of the six**. The two families differ in organisation and in arrival, but not in **external causes** | `#137` |
| **The relational family arrives together yet contains several independent dimensions** | The early/late line has the **largest structural asymmetry of six coherent splits** (5.73, 3.0×): the late family has **2.42 more effective dimensions**, **0.35 lower** split-half reliability, **0.14 lower** single-factor share. **It is one package in time and not one thing in preference structure** — the concrete family is the reverse: scattered in arrival, more like a single thing (effective dimensionality 7.64 vs 5.30 at matched item count). **And what it contains is the three axes `A02` already named** (who submits · who is seen · who receives; held-out correlations 0.81/0.81/0.91 against a random-loading floor of 0.349) — the relational family is the **intersection** of three near-independent axes | `#139` `#140` |
| **What arrives early is scattered; what arrives late is a set** | Three blocks against a within-person permutation null: early×early **−0.01272 (27.9×)** · **late×late +0.02103 (23.1×)** · late×early −0.00484 (15.2×). What arrives late is the relational family (power dynamics · bondage · sadomasochism · mental alteration · sensory, 26–31% after age 19); what arrives early is the concrete family (appearance 8.6% · clothing 10.9% · body parts 11.1%). At the item level: the late family's mutual connectivity is **+0.0674**, the early family's −0.0075. **And "a set" is literal**: the relational family's acquisition ages are **9.4% tighter** than an equally-sized random draw (17.7×, after removing the population schedule), while the concrete family is **9.0% more dispersed** (12.1×). **Development is not one kind of thing arriving in order; it is two kinds of thing, differently organised, arriving differently** — a dimension `#75`'s schedule does not have. ⚠ Connectivity has **no knee** (1.6×), so "closes at 17" holds only for the category count | `#135` `#136` |
| **For a third of people the map closed at 17; for the rest it opens outward** | In the 29–32 band, **33.3% report no new interest at all after 17**. Among those with later entries, the later-acquired categories are **less connected** to their earlier ones, retaining 89% after stripping rarity at the pair level, with `#114` opposite-signed, and **all five age bands agreeing in sign (5.3–10.9×) with no age gradient** — so it is not narrative reconstruction. **Expansion, not deepening** — model C cannot be localised to "reweighting what is already there". ⚠ **Direction only, not magnitude**: cuts from 13.5 to 17.5 give −0.0165 → −0.0069 (a factor of 2.4), and a **within-person median split is unresolvable (1.0×) — the phenomenon is tied to absolute age, not to "later in your own sequence"** | `#133` `#134` |
| **The sexual map is largely fixed by the end of adolescence** | Among 29–32 year olds, **68.4% of their own reported interests were acquired before 17**, and the latest one averages **22.6 years** (a within-person measurement, no cross-sectional assumption). Cross-sectionally the category count grows only from 12.4 to 12.9 over fifteen years (+4.5%; cohort-confounded, weaker). **This puts a time boundary on model C: if recursive reshaping exists, its window is mostly inside adolescence** | `#132` |
| **At the population level the map really does start from what everyone has** | The batch of interests a person reports earliest sits at the **33rd percentile** of their own repertoire by rarity (49.4×). **In units of each person's own permutation spread this is z = −0.5515**, and **every one of five breadth strata is resolvably negative (17.5–27.3×) — the effect is universal, independent of how many things a person likes** (`#147`). ⚠ `#146`'s "the broader your taste the stronger it is" is **withdrawn**: the raw Δ's 27-fold variation is the arithmetic bound `(k−m)/k`, and after normalisation the residual gradient is *smaller* than what a deliberately breadth-independent planted effect produces. The content line does not explain it either (PC4 only 0.2×), nor does rarity dispersion (+0.191). **The mechanism is located** (`#148` `#152`): **it is the shared developmental schedule projected onto each person** — **five different definitions of the schedule all explain most of it, sign consistent (61%–104%)**. ⚠ **The magnitude is specification-dependent**: `#148`'s "80%, exceeding a random baseline by 44 points at 20.0×" is **withdrawn** — that baseline was **one** random draw (20 draws span −61%..+118%) and the denominator was wrong; the correct reading is **61 ± 44 points = 1.4×, unresolvable**. The remaining **9.4×** was taken to be "personal" — and `#159` shows **it is almost not about people at all**: split a person's interests in half and their "earliness tendency" has a split-half reliability of only **0.105**; split the **people** in half and **the category-loading profile reproduces at 0.96 (0.81 after removing rarity)**. **Certain sexual interests simply arrive earlier or later than the shared schedule predicts, and it is the same profile for entirely different groups of people.** `#149` says **that part is real but not a trait**: two independent readings agree on the same set of interests (−0.5741, **91.6×**) but only **−0.0910** on the person's **other** half (6.3× smaller). **What uncommon things you like is a stable property (`#100`: +0.4611, 23.1×); when you got them is not.** ⚠ **Mechanism UNVERIFIED**: at the item level Spearman(rarity, **median** onset) = +0.437, but the median is a **worse** schedule (held-out pairwise order 63.30% vs the mean's 66.70%), and the mean ordering is unrelated to rarity (+0.011). The left-tail explanation is dead (−0.091). **Effect stable, mechanism open** | `#130` `#131` |
| **The rarer a person's eventual taste, the larger that head start** | −0.0459; **after 1:1 caliper matching on the number of categories answered, −0.0417 (3.1×, 91% retained)**, three matching specifications agreeing at 85–100%, permutation null 17%, planted positive control +0.5878. `#114`'s recall bias contributes 39%. **It reproduces on a discrete statistic: matched +0.0532, 5.1×**. **And `#150`'s pre-registered point prediction landed**: measured instead as "how rare is the single earliest one" gives **+0.0370** (predicted +0.0404, sign reversed by definition, discriminant 0.2×) — **the statement does not depend on which statistic I use**. ⚠ But CONSISTENT, not PROVEN: the equivalence bound can only exclude a discrepancy larger than 86% of the effect, and `#151` shows **this release has no instrument that can weld it shut** (forced choice is cleaner per item but has only 10 blocks; reliability 0.252 vs 0.461, and the net result is a loss of resolution) | `#128` `#129` `#130` `#150` `#151` |
| **Its only external anchor that holds is sex** | +0.093 (disattenuated +0.141). All Big Five \|r\| ≤ 0.056, **openness only +0.023** — it is not general novelty seeking. **The null's strength is now quantified** (`#154`): a single curveball null can by itself produce personality correlations as high as **0.054**, while `R15` averaged three draws, leaving a residual of ≈0.004 — a 23× margin. **This claim is safe by design, not by luck** | `#101` `#102` `#154` |
| **People remember their most-loved sexual interests as arriving about nine months earlier** | −0.2000 years per rating sd (**19.8× SE**); a within-category shuffle null is **0–11%** of the effect (`#153`: the original round drew it **once**, and that draw was the smallest of 20; the null's realisation sd is 5.5% of the effect); the planted ladder is monotone. **This is the first recall bias directly measured in this dataset, and the entire maturation schedule rests on top of it**. ⚠ **It is a population regularity, not a personal one** (`#162`): the per-person slope's split-half reliability is **at or below its own null's realisation noise** (0.054 vs a null sd of 0.058), while a planted person-specific bias is recovered at **0.41–0.55** — so this is a real absence, not a power failure. **"Some people are more prone than others to remembering their favourites as early" does not hold.** | `#114` `#162` |
| **And that memory distortion deepens with time, nearly doubling over fifteen years** | the 15-year-old band −0.0505 → the 30-year-old band −0.0917, an age trend of **3.5×**; the residual artefact is **opposite-signed**, so correcting only amplifies. **It is not generated at the moment of answering the survey** | `#119` |
| **Among interests of the same kind, which one came first predicts your whole remaining preference profile** | **+0.0159, 6.1× spread**, after subtracting each of the two ratings, their difference, and the person's mean rating. **The pairwise design is structurally immune to response level** — its sibling claim ("the earlier it arrives the more central it is") died on exactly that confound | `#107` `#110` `#116` `#117` |
| **And the remaining preferences are pulled toward whichever came first** | displacement **+0.0339, 3.1×**, positive in 46 of 68 pairs, with a generative positive control firing monotonically. ⚠ **Doubling the sample did not strengthen it** (3.093 → 3.048), and the source line claiming "strengthening succeeded" is withdrawn; the ratio's first decimal is itself inside the bootstrap noise (20 seeds: 3.023–3.163) | `#118` `#142`, labelled DESCRIPTION: direction does not discriminate causality |
| ~~**"Top/bottom" is one word over three near-unrelated things**~~ **downgraded to: the three role axes are strongly correlated but not the same thing** | **`Entry 24` withdrew "three near-independent axes" and "effective dimensionality 2.95 of 3" long ago.** With a validated reliability ladder (positive control: SUBSTANCE against itself r_true = **+1.018**; sham +0.023), **POWER–SUBSTANCE r_true = +0.605**, i.e. **37% shared variance, not 5%**. The pre-registered kill lands in the middle band (0.45 < 0.605 < 0.70) → **UNVERIFIED**: neither one construct nor confirmed distinct. ⚠ **`#141` is withdrawn in full** — it re-priced the claim by re-running `A02/R034` without reading `Entry 24`, put the withdrawn "2.95/3" back on this line, and recommended quoting the **observed** correlation (0.112), which is precisely the attenuated one | `Entry 24` · `#143` |
| **"Extreme" is not one axis; it is two mutually exclusive directions** | abjection and filth ↔ blood, burning, weapons — people who like one end systematically dislike the other | A02 |
| **82.7% say "porn gave me this interest", and that sentence leaves no trace in either their preference structure or their timeline** | No timing signature (a 1-year shift is excluded; the concentration discriminator is null at **8.8× the required power**) · no structural signature (misfit < 0.1 sd, extremity difference < 10%) · it **tracks how much the person ticked overall** (rho **+0.2922**). ⚠ But the step "85% of that is response style" was downgraded to **UNVERIFIED** by `#26` — every item is erotic content with no reverse-keying, so general agreeableness and general erotic endorsement are not separable in this release | A06 · `#26` |

---

# What was retracted (kept here because I once reported it)

| What I said | Now | Why |
|---|---|---|
| "Sexual is a value, not a category" (the epoch title) | **False, and the question was two questions** | Line 1 of the loader deleted the rival; and within a block base rate ≡ content, so that contrast can never adjudicate `#67` `#69` `#70` |
| "Domain-general beats domain-specific, 186/276" | **Inverted: 0 general / 201 specific** | The 186 was a −0.135 estimator handicap I had given the specific side `#82` |
| "The three components are the same size (in predictability), 1.05×" | **Withdrawn** | Subtracting a null does not fix an estimator; it charges the estimator's failure to the estimand `#86` `#87` `#88` |
| "A minority concentrates on rare options" | **Wording withdrawn, existence retained** | It is a symmetric widening, not a subgroup `#99` |
| "The earlier something enters the map, the more central it is" | **Withdrawn, 10% remains** | That was not time, it was "people who rate everything highly" `#115` |
| "The trait tracks earlier acquisition age" | **Withdrawn** | −0.030 against a null of −0.028; the null is 91% of the effect `#102` |
| "The √n accumulation across blocks reflects latent structure" | **Inverted** | At fixed budget, more blocks is worse; it is ordinary √N `#64` |

---

# What this data cannot do (each measured, not guessed)

| Cannot | The measured bound |
|---|---|
| See structure carried by ≤30% of people (variance-explained methods) | 5% carriers at **±50 pp** leave skill on the floor; ≈30% is required `#91` |
| Distinguish "a small group unusually intense about a few things" from "one heavy-tailed continuum" | Margin-preserving concentrated/diffuse plants all sit inside seed jitter at 3000 swaps per block; **the way out is an external anchor, not a stronger plant** `#122` |
| Distinguish memory becoming blurry from a story becoming settled | Real curvature **0.3×**; pushing it to 2× needs about **48× the sample** `#120` |
| Measure the interaction magnitude directly, without a model | **Algebraically impossible** — the residual's first moment is margin-determined and identically zero under a margin-preserving null `#105` |
| Count dimensions | An estimator with an honest floor cannot count; the one that can count has a floor that is an artefact `#89` |
| Recover presentation-order primacy | The release exports multiselect answers **alphabetically** (119 pairs, agreement 1.0000); display order is destroyed `#69` |
| Separate breadth from acquiescence | Needs reverse-keyed or forced-choice items; this release has neither |
| Causal direction | Cross-sectional, and a stripped release |

**And a transferable methodological result — but a much smaller one than I first published.** Block-count matching takes pornhabit's congruence deficit from 0.2285 → 0.0871 → 0.0439 **by direct measurement on that one split**, and **reverses its ordering against sex**; that correction stands and never depended on a correlation — **though it is a single-split direct measurement whose own precision has never been estimated (`#169`, scope stated rather than promised)**. What does **not** stand is the law I built on top of it: ~~`corr(deficit, coverage gap) = +0.815` over nine splits~~ was **killed by `#20`** — median r over K ∈ {3,5,8,12} is **+0.127**, the worst leave-one-out is **−0.294** (removing pornhabit, at 6 of 8 specifications), and the permutation null at this n is enormous: **|r| p95 = 0.820 at n=6 splits**, so +0.815 never cleared its own null. Sex and modality both carry large deficits at near-zero coverage gaps. **The unit was never people** — all 15,503 respondents appear in every split, but the estimand is a correlation over *splits*, and there are six. ⚠ And "nobody has pointed this hazard out" is **false, killed by `Entry 15`** — the explorer's own `analysis/swarm/14-missingness.md`, **five months earlier**, states that group comparisons on gated columns compare different subpopulations. **The hazard is theirs; one correction on one split is mine; the law was nobody's.** (`#11` `#15` `#20`)

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
**[DOWNGRADED — RETRACTIONS #24, re-affirmed #143]** "three near-independent axes / effective dimensionality 2.95 of 3" was withdrawn long ago. With a validated reliability ladder (positive control: SUBSTANCE against itself r_true = **+1.018**; sham +0.023), **POWER–SUBSTANCE r_true = +0.605** — they share **37% of variance, not 5%**. The pre-registered kill lands in the middle band (0.45 < 0.605 < 0.70): **UNVERIFIED**, neither one construct nor confirmed distinct. **RETRACTIONS #141 is withdrawn in full** — it re-priced the claim by re-running A02/R034 without reading the ledger entry that had already superseded it.

**[R04 · acquisition and time](A04_acquisition_and_time)** — Interests arrive on a
population-shared schedule (content-like early: appearance 14.0, body parts 14.4, clothing 14.7;
relational late: power dynamics 16.8, bondage 16.9, mental alteration 17.0). **[SURVIVES FRAMING SWAP — RETRACTIONS #50]** Within-person "acquired together" tracks "liked together" at RSA **+0.599**, and unlike the cross-domain grammar (#49) it carries predictable variance: person-level onset→preference R² = +0.0136 against a −0.0043 null, 31% of the same-domain ceiling after stripping intensity leakage
(−0.126 — **[VERIFICATION — see RETRACTIONS #16]**, published as `01-age-onset` Finding 1) and recall anchoring, and after near-synonymy is excluded (0.594 for
280 pairs sharing no content word vs 0.646 for 64 that do). **[CONFIRMED — RETRACTIONS #22]** onset carries structure preference does not: observed top residual eigenvalue 0.959 against a purpose-built rival world at 0.441 ± 0.040 (95% upper 0.532), positive control fires at injection 0.3 and stays silent at 0.0. The attached '80% of the sd' is descriptive, not tested — interests a person does *not* like together were acquired at the
same time. But the organising variable is arrival time, not coordinate membership (retraction #9).

**[R05 · group differences and the instrument](A05_group_differences_and_the_instrument)** —
**One split's group comparison on this release was partly measuring survey coverage**: block-count matching took pornhabit from 0.2285 → 0.0871 → 0.0439 and **reversed its ordering against sex** — measured directly on that split, and it stands. **[KILLED — RETRACTIONS #20]** the law I built on it — ~~corr(congruence deficit, coverage gap) = +0.815 across nine splits~~ — does not: median r over K ∈ {3,5,8,12} = +0.127, worst leave-one-out = −0.294 on removing pornhabit, and the permutation null has **|r| p95 = 0.820 at n=6 splits**, which +0.815 never cleared. The estimand was a correlation over *splits*, not over the 15,503 people who appear in every one of them.
**[CORRECTED — RETRACTIONS #15]** the sentence that used to follow — "not documented anywhere I have
read about this dataset" — **is false**: the explorer's own `analysis/swarm/14-missingness.md`, five
months earlier, states that group comparisons on gated columns compare different subpopulations.
**The hazard is theirs; one correction on one split is mine; the law was nobody's (`#20`).** Drawn-vs-live-action consumers — whose content contains no real bodies — differ by **0.0204 ± 0.0265**, bounded below sex-sized (0.093). **[UNVERIFIED — RETRACTIONS #34]** the written-vs-visual deficit is unresolvable at 5 seeds (0 of 9 corpus cuts have effect > 2× seed spread, while sex is resolvable in 9/9); ~44 seeds per cell would be needed. That is a real, bounded constraint on A.

**[R06 · induction](A06_induction)** — the retrospective induction claim, its two absent signatures and the breadth it actually tracks are stated once in the standing-claims table above; this section adds only what the table does not carry. **[CORRECTED — RETRACTIONS #26, caught by #144]** an earlier version of this page reported "rho +0.2515, 85% surviving response-style control" — that pair **appears nowhere in the ledger** and came from a pre-A09 draft of this page. The "85% is response style" step is **UNVERIFIED**: every item is erotic content with no reverse-keying, so general agreeableness and general erotic endorsement are not separable in this release. **Design consequence**: retrospective self-report of induction is unusable as an outcome measure — phase 3 must be prospective.

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
for log and 51.9% for linear, with a permuted-label null of exactly 0.0000 at every n. **[MECHANISM REVERSED — #61, #64]** blocks are *not* interchangeable (subset variance 58× seed variance), and at
**fixed total respondent-rows more blocks is worse**: n=8 gives 0.0042, n=16 gives 0.0028, n=31 gives 0.0025.
The accumulation tracks **total sample**, not block count — the √ shape is ordinary √N. Block boundaries are a
tax on the estimate, not a source of it. Price list for new collection:
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

## Tools: six edges

| tool | what it checks |
|---|---|
| `tools/readme_ledger_audit.py` | README → ledger (has a number here been withdrawn by a later entry; is any number uncited) |
| `tools/round_status.py` | round → ledger (current status of the claims this round supports; has a sibling round superseded it) |
| `tools/guard_lint.py` | round → guards (were the applicable guards called; is a verdict threshold a hard-coded literal) |
| `tools/repro_audit.py` | ledger → round (does the ledger's number still come out; **is a value the ledger condemned still being printed**) |
| `guard_lint.error_bar_scan()` | artifact → precision (**does the real arm have any source of jitter at all**, `#168`) |
| `lib/gates.py` | 10 guards + the `Gate` comparison rules every round attacks itself with |

Each tool carries its own `P6` proxy ledger: **readable in the hit direction only; the output is a
must-read list, not a verdict.**

## Governance

| file | what it holds |
|---|---|
| [`RETRACTIONS.md`](RETRACTIONS.md) | **167 entries.** Every claim killed, scoped or corrected, with what killed it. Twelve are a later round of mine destroying an earlier round of mine |
| [`FROZEN.md`](FROZEN.md) | Lines where further computation cannot identify what it is measuring, each with its unfreeze condition |
| [`PREREGISTRATION.md`](PREREGISTRATION.md) | What the next rounds will test, with thresholds fixed in advance. Rounds r01–r52 were **not** preregistered and are labelled exploratory throughout |
| [`ADVERSARY_FORECAST.md`](ADVERSARY_FORECAST.md) | What I expect an outside challenger to overturn, written before one arrives |
| [`STANDARD_AUDIT.md`](STANDARD_AUDIT.md) | These 52 rounds scored against the campaign standard, including everything they fail |

## Layout

**E · A · R** — an *epoch* closes when the ontology shifts (the object turned out to be a different
object), an *arc* closes when a decision becomes safe, a *round* closes on one belief update.
**One epoch, 25 arcs, 213 rounds.** `R` is numbered consecutively across the whole project, not
restarted per arc. Every count is discovered, never chosen; see `~/.claude/CLAUDE.md` §P16.

Each round is a directory with `run.py`, `README.md` and `results/`.
`lib/rounds.py` maps a round name to its path, because several rounds reuse an earlier round's
loaders and that dependency is made explicit rather than hidden. Environment is a self-contained
`.venv` (system python 3.14 has no `ensurepip`; pip was bootstrapped).

## Scope, stated once

Cross-sectional · one instrument · one population (18–32, US/Canada/Europe) · aggressively
anonymised, with correlations attenuated roughly 25% by design · every measure self-report.
No causal claim in this repository is identified, and none is made.
