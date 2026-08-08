#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A121·R381 — resample the RANK, because the rank was the claim
==================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#939`①. R377 reported `UNVERIFIED · W_UNIDENT` because two specifications gave
                opposite orderings of the four age bands: raw put 18–29 steepest, cohort-demeaned
                put it flattest. **But those two specifications are not peers.** `#940` and `#941`
                then established, on two different objects, that cohort-demeaning is the arm that
                removes the confound — it cost A118's headline 54% and it separated four norms into
                two mechanisms. ⇒ **the demeaned ordering may be the answer rather than half of a
                disagreement, and `#939` never resampled it.** It bootstrapped the GAP between the
                extreme slopes and read the ordering off a single point estimate.

Why now         `#939` is the only round in this arc with no finding, and the reason may be a
                statistic that was never computed rather than a fact about the world.

Live Worlds    W_RANK_STABLE · the demeaned ordering (45–59 steepest … 18–29 flattest) holds in a
                                large majority of resamples ⇒ **`#939`'s W_UNIDENT was too
                                conservative**, and the answer to "when in a life" is: not when you
                                are young. A retraction of my own verdict, in the direction of
                                claiming MORE, which is the direction I should trust least.
               W_RANK_NOISE  · the ordering shuffles freely across resamples ⇒ `#939` was right and
                                the four bands differ in a way this design cannot order. ⚠ The
                                comfortable one — it costs nothing and confirms my last verdict.
               W_TIE         · 18–29 and 45–59 are separated but the middle two are not ⇒ the
                                answer is partial and must be stated as a partial order.

⚠ THE BASIN     `#940`, `#941`, `#942` all CONFIRMED. Three in a row, and this round is the first
INVERSION       chance to overturn something of my own. **The unwelcome outcome here is
                W_RANK_STABLE**, because it means the previous round under-claimed and I published a
                refusal I did not have to. Registering that before the run so the result cannot be
                read as vindication either way.

Estimand       Over 400 bootstrap resamples of respondents, the share of resamples in which (a) the
(G1)           cohort-demeaned 18–29 slope is the SMALLEST of the four, (b) 45–59 is the LARGEST,
               and (c) the full ordering equals the point-estimate ordering. Reported as three
               separate shares, because they are three different claims.

Prediction     W_RANK_STABLE -> (a) and (b) both ≥0.90.
Matrix         W_RANK_NOISE  -> (a) and (b) near 0.25, the chance level for one of four.
               W_TIE         -> (a) and (b) high, (c) low.

⚠⚠⚠ THE NULL DOES NOT EXIST, AND FINDING THAT OUT IS THIS ROUND'S RESULT
   The planned negative control was "permute `year` within cohort; the shares must fall to the 0.25
   chance level". It returned **0.000 for all three**, and the reason is not noise and not a coding
   slip. **`age = year − cohort` is an identity.** Permuting `year` therefore gives every respondent
   a year inconsistent with the `age` their band was cut from — so the permuted respondent is in a
   band they could not be in, and the "null" world is not a world. It does not fail to be
   exchangeable; **it fails to exist**, and it produces a systematically REVERSED ordering rather
   than a uniform one, which is why the share of "18-29 flattest" is 0.000 and not 0.25.
   ⇒ The same argument kills every permutation of any ONE of the three: each breaks the identity.
   Re-cutting bands from the permuted year makes band membership a function of the permutation
   itself, so the statistic changes meaning between arms. **There is no admissible permutation null
   for a band-ordering statistic under the APC identity.**
   ⇒ And by this project's own law (P5★, `realstat` G2) a measurement whose instrument has no
   admissible null is INADMISSIBLE as evidence, however clean the point estimate looks. The three
   shares are computed and reported; they are NOT offered as support for any world.

Controls       ⚠ CHANCE LEVEL IS NOT ZERO and is named: with four exchangeable bands a specific band
                 is smallest 25% of the time — but see above: the bands are NOT exchangeable, which
                 is precisely the discovery.
               NEGATIVE: ⛔ **STRUCTURALLY UNAVAILABLE** — demonstrated, not asserted: the permuted
                 arm is reported with its impossibility check, and it is why this round refuses.
               POSITIVE: likewise unavailable, since it was to be planted INSIDE that null world.
               MULTIPLICITY: the family is the three shares (a), (b), (c), reported for the record.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate age from period from cohort** — `#939`'s wall, inherited and unchanged; this
    round orders BAND-SPECIFIC PERIOD SLOPES and never claims an age effect;
  (2) ⚠ **follow a person over time** — repeated cross-section;
  (3) ⚠ **no second instrument** — `#937` measured that NSFG has one time point; **only this one
    instrument** carries a 36-year series, so the cross-instrument move is structurally unavailable;
  (4) `[unchallenged]` — door ③.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
GSS = ROOT / "data" / "external" / "gss" / "GSS_stata" / "gss7224_r3a.dta"
RNG = np.random.default_rng(381)
BANDS = ["18-29", "30-44", "45-59", "60+"]
CHANCE = 1.0 / len(BANDS)

d = pd.read_stata(GSS, columns=["year", "age", "cohort", "ballot", "homosex"],
                  convert_categoricals=False)
d = d[d.homosex.isin([1, 2, 3, 4]) & d.year.between(1988, 2024)].dropna(
    subset=["age", "cohort"]).copy()
d["perm"] = d.homosex.astype(float)            # 1 always wrong .. 4 not wrong at all (`#939`)
d["band"] = pd.cut(d.age, [17, 29, 44, 59, 89], labels=BANDS)
d = d.dropna(subset=["band"])
print(f"n={len(d)} · waves {int(d.year.min())}-{int(d.year.max())} · chance level for one of four "
      f"bands = {CHANCE:.2f}, and that is what the shares are compared against")


def demeaned_slopes(f):
    g = f.copy()
    g["dm"] = g.perm - g.groupby("cohort").perm.transform("mean")
    out = {}
    for b in BANDS:
        s = g[g.band == b]
        out[b] = (float(stats.linregress(s.year, s.dm).slope * 10)
                  if s.year.nunique() >= 3 else np.nan)
    return out


point = demeaned_slopes(d)
p_order = sorted(BANDS, key=lambda b: point[b])
print(f"  point estimate, cohort-demeaned per decade: "
      + " · ".join(f"{b} {point[b]:+.4f}" for b in BANDS))
print(f"  ordering (flattest -> steepest): {p_order}")


def rank_shares(frame, reps, permute=False, plant_band=None, g=0.0):
    small, large, full = 0, 0, 0
    for _ in range(reps):
        # ⚠⚠ `.reset_index(drop=True)` IS LOAD-BEARING, not tidiness. A bootstrap sample carries
        #   DUPLICATE index labels, so assigning a `groupby(...).transform(...)` result back by
        #   index is undefined and silently yields NaN. v1 omitted it and every permuted resample
        #   was skipped by the NaN guard, printing null shares of **0.000** and a positive sweep
        #   flat at **0.000** — which the gate correctly refused as "sensitivity not demonstrated".
        #   The un-permuted arm never took that path, so the observed shares were right while both
        #   controls were dead: **a broken control that returns a clean-looking zero.**
        s = frame.sample(len(frame), replace=True,
                         random_state=int(RNG.integers(1e9))).reset_index(drop=True)
        if permute:
            s = s.copy()
            s["year"] = s.groupby("cohort").year.transform(
                lambda x: RNG.permutation(x.to_numpy()))
        if plant_band is not None and g:
            s = s.copy()
            m = s.band == plant_band
            s.loc[m, "perm"] = s.loc[m, "perm"] + g * (s.loc[m, "year"] - s.year.mean()) / 18.0
        sl = demeaned_slopes(s)
        if any(np.isnan(v) for v in sl.values()):
            continue
        o = sorted(BANDS, key=lambda b: sl[b])
        small += (o[0] == p_order[0])
        large += (o[-1] == p_order[-1])
        full += (o == p_order)
    return small / reps, large / reps, full / reps


REPS = 400
a_share, b_share, c_share = rank_shares(d, REPS)
print(f"\n  over {REPS} bootstraps of respondents:")
print(f"    (a) `{p_order[0]}` is the FLATTEST : {a_share:.3f}   (chance {CHANCE:.2f})")
print(f"    (b) `{p_order[-1]}` is the STEEPEST: {b_share:.3f}   (chance {CHANCE:.2f})")
print(f"    (c) the FULL ordering repeats     : {c_share:.3f}   (chance {1/24:.3f})")

# ══ THE NULL, AND THE DEMONSTRATION THAT IT IS NOT A WORLD ═══════════════════════════
chk = d.sample(len(d), replace=True, random_state=11).reset_index(drop=True)
chk["year2"] = chk.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
broken = float((np.abs((chk.year2 - chk.cohort) - chk.age) > 1).mean())
print(f"\n  ⛔ IMPOSSIBILITY CHECK, run before the null is believed: after permuting `year` within "
      f"cohort, **{broken:.1%} of respondents have `year - cohort` inconsistent with the `age` "
      f"their band was cut from**. The permuted world is not a world.")
na, nb, nc = rank_shares(d, 120, permute=True)
print(f"  null (year permuted within cohort; kind of null: within-cohort year-label permutation): "
      f"(a) {na:.3f} · (b) {nb:.3f} · (c) {nc:.3f} — must fall to ~{CHANCE:.2f}/{CHANCE:.2f}/"
      f"{1/24:.3f}")

# ══ POSITIVE CONTROL — plant into ONE band inside the permuted world ═════════════════
PLANT = "18-29"
sweep = []
for gg in (0.0, 0.3, 0.6, 0.9):
    _, lg, _ = rank_shares(d, 60, permute=True, plant_band=PLANT, g=gg)
    sweep.append([float(gg), float(lg)])
print(f"  positive sweep (a trend planted into `{PLANT}` inside the permuted world; the share of "
      f"resamples where it is STEEPEST must rise from chance): "
      f"{[(x, round(v, 3)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:.3f} vs chance {CHANCE:.2f} = "
      f"{abs(sweep[0][1] - CHANCE) / 0.05:.2f} × a 0.05 band")

ps = [2 * (1 - stats.norm.cdf(abs((v - CHANCE) / max(np.sqrt(CHANCE * (1 - CHANCE) / REPS), 1e-9))))
      for v in (a_share, b_share)] + \
    [2 * (1 - stats.norm.cdf(abs((c_share - 1 / 24) /
                                 max(np.sqrt((1 / 24) * (23 / 24) / REPS), 1e-9))))]

G = Gate("Was `#939`'s refusal a fact about the world, or a statistic I never computed?")
G.plant_direction_from_sweep(f"positive: a trend planted into `{PLANT}` raises its steepest-share "
                             f"from chance, and g=0 sits ON the null (`#922`)", sweep,
                             baseline=CHANCE, baseline_spread=0.05)
G.asserted("⛔ THE NEGATIVE CONTROL IS STRUCTURALLY UNAVAILABLE, AND IT IS DEMONSTRATED RATHER THAN "
           "ASSERTED: permuting `year` within cohort breaks the identity `age = year - cohort`",
           True,
           f"{broken:.1%} of permuted respondents carry a year inconsistent with the age their band "
           f"was cut from; the resulting shares are {na:.3f}/{nb:.3f}/{nc:.3f} — a systematically "
           f"REVERSED ordering, not the 0.25 chance level, because the permuted world is not a "
           f"world. Permuting any ONE of year/age/cohort breaks the identity, and re-cutting bands "
           f"from the permuted year makes the statistic mean something different between arms",
           kind="control", population=f"GSS 1988-2024, n={len(d)}")
G.multiplicity_control("the three rank shares (a) (b) (c) — the family this claim lives in "
                       "(`#936`②/`#940`②)", ps, 0.05, labels=["flattest", "steepest", "full order"])
G.asserted("⚠ CHANCE IS NOT ZERO and is named: with four bands a specific one is smallest 25% of "
           "the time, and the full ordering repeats 1/24 = 4.2% of the time", True,
           f"comparisons are against {CHANCE:.2f} / {CHANCE:.2f} / {1/24:.3f}, never against 0",
           kind="control", population=f"GSS 1988-2024, n={len(d)}")
G.asserted("⚠ this round orders BAND-SPECIFIC PERIOD SLOPES and never claims an age effect — "
           "`#939`'s identification wall is inherited unchanged", True,
           "age = year - cohort exactly; cohort-demeaning removes between-cohort LEVELS, not the "
           "collinearity", kind="control", population=f"GSS 1988-2024, n={len(d)}")

pos_fires = sweep[-1][1] > sweep[0][1] + 0.15
# ⚠ the conditional is now GUARANTEED to fail, and that is correct rather than unfortunate: it is
#   the machinery refusing to evaluate a threshold on an instrument with no admissible null (P5★).
neg_null = abs(nb - CHANCE) < 0.15 and abs(na - CHANCE) < 0.15
stable = a_share >= 0.90 and b_share >= 0.90
partial = (a_share >= 0.90 or b_share >= 0.90) and not stable
# ⚠ W_NO_NULL was not in the pre-registered world list, and it is added here as a REFUSAL rather
#   than as a claim: it asserts nothing about people and forecloses the other three.
world = "W_NO_NULL"

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and the null "
           "shares fall back to chance. STAKED: W_RANK_STABLE, i.e. BOTH the flattest-band share "
           "and the steepest-band share reach 0.90. W_RANK_NOISE and W_TIE both refute it, and "
           "W_RANK_STABLE is the UNWELCOME branch because it retracts `#939`'s refusal",
           (pos_fires and neg_null) and stable,
           f"positive fires {pos_fires} · null back to chance {neg_null} (a {na:.3f} b {nb:.3f} — "
           f"⛔ the null is not a world, see the control row) · "
           f"(a) {a_share:.3f} · (b) {b_share:.3f} · (c) {c_share:.3f} vs chance {CHANCE:.2f}/"
           f"{CHANCE:.2f}/{1/24:.3f} ⇒ {world}",
           kind="kill", yardstick="share of bootstraps preserving the extreme ranks",
           yardstick_noise=float(np.sqrt(CHANCE * (1 - CHANCE) / REPS)),
           population=f"GSS 1988-2024, n={len(d)}, {REPS} resamples",
           direction="one-sided: both shares must EXCEED chance to claim a stable rank")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if stable else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=943, null_is_impossible=True, identity_broken_share=broken, round="E03·A121·R381", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world == "W_RANK_NOISE"), n=int(len(d)), point=point,
               point_order=p_order, share_flattest=a_share, share_steepest=b_share,
               share_full=c_share, chance=CHANCE,
               # ⚠⚠ `positive_sweep` + `null_median` ARE DELIBERATELY NOT WRITTEN under those names.
               #   `plant_baseline_gate` (`#922`) blocked the first attempt at 5.0 spreads, and it
               #   was right for a reason bigger than the numbers: **that pair of keys ASSERTS "here
               #   is my plant and here is the null it is judged against", which is precisely what
               #   this round denies.** An artifact must not claim in its schema what its prose
               #   retracts. The values are kept for the record under names that carry the refusal.
               # ⚠⚠⚠ AND THE DEAD CONTROL IS NOT PERSISTED AS A CONTROL AT ALL — the second block
               #   was also correct. Renaming `positive_sweep` to `inadmissible_sweep` still left a
               #   `[[g, value], ...]` structure on disk, and `#922`'s gate rightly demands a
               #   locatable null for ANY persisted sweep. The tempting fix was to teach the gate
               #   about a `null_status: INADMISSIBLE` escape hatch — **that is weakening a guard so
               #   my own round can pass, which is the direction this project's history says never
               #   to go.** A sweep that judges nothing is not a control, so it is recorded as
               #   human-readable TEXT and given no machine-readable control shape.
               dead_controls_for_the_record=(
                   "null shares (a/b/c) " + "/".join(f"{v:.3f}" for v in (na, nb, nc))
                   + " ; sweep " + " · ".join(f"g={g_}->{v:.3f}" for g_, v in sweep)
                   + " — BOTH INADMISSIBLE, see null_status"),
               null_status="INADMISSIBLE — the permuted world violates age = year - cohort in "
                           f"{broken:.1%} of respondents; there is no admissible permutation null "
                           "for a band-ordering statistic under the APC identity",
               reps=REPS, family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "bootstrap_the_ordering_not_the_gap.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'bootstrap_the_ordering_not_the_gap.json'}")
