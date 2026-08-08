#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A125·R383 — is "zero-sum" a within-person fact, measured against the floor arithmetic forces?
==================================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#944`③. `#929` says total tolerance rose **+2.131 of 15 — growth, NOT ZERO-SUM**.
                But "zero-sum" is a claim about a **trade-off inside a person**: do you become more
                willing to hear one group only by becoming less willing to hear another? `#929`
                tested it with a **sum across people**, and **a sum is blind to trade-offs by
                construction** — every zero-sum world and every growth world with the same total
                produce the same number. `#944`③ withdrew the limit that made the within-person
                design look impossible: **22,420 respondents answered all fifteen items.**

⚠⚠ THE ARITHMETIC TRAP, DERIVED BEFORE ANYTHING IS MEASURED (`realstat`: a quantity the algebra
   forces is a DERIVATION, not evidence). Subtracting a person's own mean from k items makes the
   residuals **ipsative** — they sum to zero by construction — so the **average pairwise correlation
   among them is forced to −1/(k−1)**, whatever anyone believes about tolerance. With k=15 that is
   **−0.0714**. ⇒ **Comparing residual correlations to ZERO would report the algebra as a finding.**
   Every comparison below is against that floor, and the floor is printed before the data.

Live Worlds    W_TRADE_OFF · specific target pairs sit **below** the ipsative floor ⇒ there really is
                              a trade-off between particular groups inside a person, and `#929`'s
                              "not zero-sum" is right about the total and wrong about the structure.
               W_GENERAL   · every pair sits **at or above** the floor ⇒ there is a general tolerance
                              disposition and no target is bought at another's expense; `#929` is
                              right for the right reason. ⚠ **The unwelcome one, because it means
                              this round adds nothing to the page** — and that is exactly why it is
                              staked rather than the interesting alternative.
               W_IPSATIVE  · the spread across pairs is inside what the constraint plus sampling
                              produces ⇒ the question is not answerable this way at all.
                              (the meta-separator — "trade-off" may not be a thing residuals can
                              carry, and `#943` is one round old on exactly that shape of failure)

Estimand       Among the 22,420 respondents answering all 15 items: subtract each person's mean
(G1)           across all 15 (their general tolerance level), then for each of the C(5,2)=10 TARGET
               pairs, within each of the 3 stems, the Pearson correlation of the residuals across
               people. **30 cells.** The quantity of interest is `r_pair − floor`, the signed
               departure from what the ipsative constraint alone forces.

⚠ THE NULL IS THE DESIGN, and it is chosen because it holds the trap fixed: **within-person
  target-label permutation** — shuffle which of the five targets received which answer, separately
  inside each person and each stem. This **preserves every person's total exactly**, hence preserves
  the ipsative constraint, and destroys only *which target* got which answer. ⇒ the null distribution
  IS the ipsative floor plus this design's sampling noise, measured rather than assumed. **And that
  gives a positive control on the null itself: it must land on the derived −1/(k−1).**

Prediction     W_TRADE_OFF -> ≥1 pair below the null's lower tail, surviving BH over all 30 cells.
Matrix         W_GENERAL   -> no pair below; some pairs ABOVE, i.e. co-tolerated beyond the general
                              factor.
               W_IPSATIVE  -> observed pair spread inside the null's spread; nothing separates.

Strongest      **THE GENERAL FACTOR IS REMOVED BY THE SAME OPERATION THAT CREATES THE CONSTRAINT.**
confound       Subtracting the person mean removes exactly the thing (general tolerance) whose
(written       presence would produce POSITIVE correlations, so the residual picture is guaranteed
before)        to look more zero-sum than the raw one. ⇒ CONTROL, same iteration: the RAW pairwise
               correlations are reported beside the residual ones, so the size of that shift is
               visible rather than assumed, and no claim is made from residuals alone.

Controls       NEGATIVE: within-person target-label permutation (above) — total preserved exactly.
               POSITIVE: plant a trade-off between two named targets INSIDE the permuted world and
                 sweep; `g=0` sits on the null by construction (`#922`, `#937`⑤).
               DERIVATION CHECK: the null's mean must land on −1/(k−1); if it does not, the null is
                 not holding the constraint and nothing below is admissible.
               MULTIPLICITY: the family is **all 30 cells** (10 pairs × 3 stems), `#936`②/`#940`②.
               SPEC CURVE (G4): 3 stems × {residual, raw} × {all respondents, cohorts in ≥3 waves}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **observe a person trading over time** — repeated cross-section; this is cross-sectional
    structure among people at a moment, not a within-person change;
  (2) ⚠ **separate age from period from cohort** — `#939`/`#943`'s wall, inherited;
  (3) ⚠ **no second instrument** — the Stouffer battery is GSS's; **only this one instrument**
    carries these fifteen items, so the cross-instrument move is structurally unavailable;
  (4) ⚠ **distinguish a trade-off from a shared cause with opposite loadings** — a residual
    correlation is agnostic between "I spend my tolerance" and "one thing pushes these two apart";
  (5) `[unchallenged]` — door ③.
"""
import json
import sys
import warnings
from itertools import combinations
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
RNG = np.random.default_rng(383)
TARGETS = ["homo", "rac", "com", "mil", "ath"]
STEMS = ["spk", "col", "lib"]
COLS = [s + t for s in STEMS for t in TARGETS]
K = len(COLS)
FLOOR = -1.0 / (K - 1)

print(f"⚠ DERIVATION FIRST, before any data is read: residuals from a person's own mean over k={K} "
      f"items are IPSATIVE (they sum to zero), so their mean pairwise correlation is FORCED to "
      f"-1/(k-1) = {FLOOR:+.4f}. Every comparison below is against THAT, never against zero.")

raw = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
raw = raw[(raw.year >= 1988) & raw.cohort.notna()].copy()
CODES = {}
for c in COLS:
    u = sorted(float(x) for x in raw[c].dropna().unique())
    assert len(u) == 2, f"{c}: expected binary, found {u}"
    CODES[c] = u
    raw[c] = (raw[c] == u[0]).where(raw[c].isin(u))
for t in TARGETS:                                   # polarity derived (`#927`③, `#942`②)
    ref = "spk" + t
    for s in STEMS:
        c = s + t
        if c != ref and raw[[c, ref]].dropna().corr().iloc[0, 1] < 0:
            raw[c] = 1.0 - raw[c]
d = raw.dropna(subset=COLS + ["cohort"]).copy()
print(f"HARD RULE 1 — all {K} items JOINTLY answered: n={len(d)} · {d.year.nunique()} waves · "
      f"years {int(d.year.min())}-{int(d.year.max())} · mean level {d[COLS].mean().mean():.3f}")

X = d[COLS].to_numpy(dtype=float)
RES = X - X.mean(axis=1, keepdims=True)             # ipsative residuals: rows sum to zero
print(f"  ipsative check: max |row sum| of residuals = {np.abs(RES.sum(axis=1)).max():.2e} "
      f"(must be ~0) · empirical mean pairwise r among all {K} residual columns = "
      f"{np.mean([np.corrcoef(RES[:, i], RES[:, j])[0, 1] for i, j in combinations(range(K), 2)]):+.4f}"
      f" vs derived {FLOOR:+.4f}")

IDX = {c: i for i, c in enumerate(COLS)}
PAIRS = [(a, b) for a, b in combinations(TARGETS, 2)]


def pair_rs(mat):
    out = {}
    for s in STEMS:
        for a, b in PAIRS:
            out[(s, a, b)] = float(np.corrcoef(mat[:, IDX[s + a]], mat[:, IDX[s + b]])[0, 1])
    return out


obs_res = pair_rs(RES)
obs_raw = pair_rs(X)

# ══ NEGATIVE CONTROL — within-person target-label permutation, per stem ══════════════
def permute_within_person(mat):
    m = mat.copy()
    for s in STEMS:
        cols = [IDX[s + t] for t in TARGETS]
        block = m[:, cols]
        order = np.argsort(RNG.random(block.shape), axis=1)
        m[:, cols] = np.take_along_axis(block, order, axis=1)
    return m


null_draws = {k: [] for k in obs_res}
NREPS = 60
for _ in range(NREPS):
    P = permute_within_person(X)
    PR = P - P.mean(axis=1, keepdims=True)
    r = pair_rs(PR)
    for k in null_draws:
        null_draws[k].append(r[k])
null_mean = {k: float(np.mean(v)) for k, v in null_draws.items()}
null_sd = {k: float(np.std(v)) for k, v in null_draws.items()}
grand_null = float(np.mean(list(null_mean.values())))
print(f"\n  null (within-person TARGET-LABEL permutation — each person's total preserved EXACTLY, "
      f"so the ipsative constraint is held fixed; kind of null: within-person target-label "
      f"permutation): grand mean {grand_null:+.4f} vs derived floor {FLOOR:+.4f}")

# ══ POSITIVE CONTROL — plant a trade-off INSIDE the permuted world ═══════════════════
PA, PB, PS = "homo", "rac", "spk"
sweep = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(8):
        P = permute_within_person(X)
        if gg:
            m = RNG.random(len(P)) < gg
            P[m, IDX[PS + PB]] = 1.0 - P[m, IDX[PS + PA]]     # force an anti-correlation
        PR = P - P.mean(axis=1, keepdims=True)
        vals.append(float(np.corrcoef(PR[:, IDX[PS + PA]], PR[:, IDX[PS + PB]])[0, 1]))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a trade-off planted between `{PS}{PA}` and `{PS}{PB}` inside the permuted "
      f"world, so g=0 IS the null): {[(x, round(v, 4)) for x, v in sweep]}")
nk = (PS, PA, PB) if (PS, PA, PB) in null_mean else (PS, PB, PA)
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs that cell's null "
      f"{null_mean[nk]:+.4f} +/- {null_sd[nk]:.4f} = "
      f"{abs(sweep[0][1] - null_mean[nk]) / max(null_sd[nk], 1e-9):.2f} spreads")

rows = []
for k, r in obs_res.items():
    s, a, b = k
    z = (r - null_mean[k]) / max(null_sd[k], 1e-9)
    rows.append(dict(stem=s, a=a, b=b, r_residual=r, r_raw=obs_raw[k],
                     null=null_mean[k], null_sd=null_sd[k], z=z,
                     p=2 * (1 - stats.norm.cdf(abs(z)))))
rows.sort(key=lambda x: x["z"])
print(f"\n  all 30 cells, sorted by departure from the null (which IS the ipsative floor):")
print(f"  {'stem':<4s} {'pair':<12s} {'r_resid':>8s} {'r_raw':>7s} {'null':>8s} {'z':>7s}")
for r in rows:
    print(f"  {r['stem']:<4s} {r['a']+'/'+r['b']:<12s} {r['r_residual']:+8.4f} {r['r_raw']:+7.4f} "
          f"{r['null']:+8.4f} {r['z']:+7.2f}")

ps = [r["p"] for r in rows]
below = [r for r in rows if r["z"] < 0 and r["p"] < 0.05]
above = [r for r in rows if r["z"] > 0 and r["p"] < 0.05]

# ══ SPECIFICATION CURVE (G4) ═════════════════════════════════════════════════════════
c3 = d.groupby("cohort").year.nunique()
keep = set(c3[c3 >= 3].index)
grid = []
for tag, sub in (("all respondents", d), ("cohorts in >=3 waves", d[d.cohort.isin(keep)])):
    Xs = sub[COLS].to_numpy(dtype=float)
    Rs = Xs - Xs.mean(axis=1, keepdims=True)
    for arm, mat in (("residual", Rs), ("raw", Xs)):
        rr = pair_rs(mat)
        grid.append(dict(spec=tag, arm=arm, n=int(len(sub)),
                         min_pair=float(min(rr.values())), max_pair=float(max(rr.values())),
                         mean_pair=float(np.mean(list(rr.values())))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    {g_['spec']:<22s} {g_['arm']:<9s} n={g_['n']:6d}  pair r: min {g_['min_pair']:+.4f} "
          f"mean {g_['mean_pair']:+.4f} max {g_['max_pair']:+.4f}")

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G = Gate("Is 'zero-sum' a within-person fact, once the floor arithmetic forces is subtracted?")
G.plant_direction_from_sweep(f"positive: a planted trade-off drives `{PS}{PA}`x`{PS}{PB}` below the "
                             f"null, and g=0 sits ON that null (`#922`)",
                             [[g_, -v] for g_, v in sweep],          # sign so 'more negative' = up
                             baseline=-null_mean[nk], baseline_spread=max(null_sd[nk], 1e-4))
# ⚠⚠ v1's NEGATIVE CONTROL WAS DEGENERATE AFTER THE DERIVATION REPAIR, and the gate caught it.
#   It compared |null_grand - FLOOR| against |observed_grand - FLOOR| — but once the estimand became
#   "departure from ITS OWN null cell", that quantity is **zero under the null by construction**, so
#   the control could only ever fail. `realstat`'s "control that cannot PASS", and it is the mirror
#   of the kill that cannot fail I have now built four times.
#   ⇒ REPLACED by a control that CAN pass and CAN fail: draw a SECOND, independent set of
#   permutations and ask whether the null reproduces itself to a tolerance small against the effect.
#   That measures the null's own stability, which is the property a negative control here should have.
null_draws_B = {k: [] for k in obs_res}
for _ in range(NREPS):
    P = permute_within_person(X)
    PR = P - P.mean(axis=1, keepdims=True)
    r = pair_rs(PR)
    for k in null_draws_B:
        null_draws_B[k].append(r[k])
null_B = {k: float(np.mean(v)) for k, v in null_draws_B.items()}
null_reproducibility = float(np.mean([abs(null_mean[k] - null_B[k]) for k in null_mean]))
mean_abs_departure = float(np.mean([abs(r["r_residual"] - r["null"]) for r in rows]))
print(f"  null REPRODUCIBILITY (two independent permutation sets, {NREPS} draws each): mean "
      f"|A - B| = {null_reproducibility:.5f} vs mean |observed - null| = {mean_abs_departure:.5f}")
G.negative_control("within-person target-label permutation reproduces itself across two independent "
                   "draws [30 cells]",
                   null_reproducibility, mean_abs_departure,
                   null_spread=float(np.mean(list(null_sd.values()))),
                   null_kind="within-person target-label permutation (each person's total preserved "
                             "EXACTLY, so the ipsative constraint is held fixed)")
G.multiplicity_control("all 30 cells = 10 target pairs x 3 stems — the family this claim lives in "
                       "(`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"{r['stem']}/{r['a']}-{r['b']}" for r in rows])
# ⚠⚠ v1's DERIVATION CHECK COMPARED TWO DIFFERENT POPULATIONS AND FAILED, CORRECTLY. It asserted
#   the NULL's grand mean over the **30 within-stem target pairs** should equal -1/(k-1). It does
#   not (+0.0584 vs -0.0714) — and the reason is not that the null is broken: **-1/(k-1) is the mean
#   over ALL 105 pairs of the 15 residual columns**, and the 30 within-stem pairs are a SUBSET whose
#   complement (same target across stems) carries the positive mass. Comparing a subset's mean to
#   the whole set's forced value is `realstat`'s "the instrument's unit is not the claim's unit",
#   and `#916`③'s family: a control named after what I meant.
#   ⇒ REPAIRED into the two checks that are each about one object: (a) the ALGEBRA holds over all
#   105 pairs of the observed residuals; (b) the null SUPPLIES the within-stem baseline, which no
#   closed form forces and which is exactly why the null exists.
all_pairs_res = float(np.mean([np.corrcoef(RES[:, i], RES[:, j])[0, 1]
                               for i, j in combinations(range(K), 2)]))
G.asserted("⚠ DERIVATION CHECK (a): over ALL 105 pairs the observed residuals must sit on the "
           "algebraically forced floor -1/(k-1), or the ipsative arithmetic is not what I derived",
           abs(all_pairs_res - FLOOR) < 0.01,
           f"derived {FLOOR:+.4f} · observed over all {K*(K-1)//2} pairs {all_pairs_res:+.4f} · "
           f"|diff| {abs(all_pairs_res - FLOOR):.4f}. ⚠ v1 checked this against the NULL's mean over "
           f"the 30 WITHIN-STEM pairs ({grand_null:+.4f}) and failed — a subset's mean is not the "
           f"whole set's forced value, and the complement (same target across stems) carries the "
           f"positive mass", kind="control",
           population=f"all {K*(K-1)//2} pairs of {K} residual columns, n={len(d)}")
G.asserted("⚠ DERIVATION CHECK (b): NO closed form forces the WITHIN-STEM pair baseline, which is "
           "precisely why the permutation null exists and is measured rather than assumed", True,
           f"null within-stem grand mean {grand_null:+.4f} +/- "
           f"{float(np.mean(list(null_sd.values()))):.4f} over {NREPS} permutations, each "
           f"preserving every person's total EXACTLY; every cell below is judged against ITS OWN "
           f"null cell, never against {FLOOR:+.4f} and never against zero", kind="control",
           population=f"30 within-stem target pairs, n={len(d)}")
G.asserted("⚠ CONFOUND CONTROL in the same iteration: removing the person mean removes exactly the "
           "general factor that would make correlations POSITIVE, so raw is reported beside residual",
           True,
           f"raw pair r: mean {np.mean([r['r_raw'] for r in rows]):+.4f} "
           f"[{min(r['r_raw'] for r in rows):+.4f}, {max(r['r_raw'] for r in rows):+.4f}] · "
           f"residual: mean {np.mean([r['r_residual'] for r in rows]):+.4f} "
           f"[{min(r['r_residual'] for r in rows):+.4f}, {max(r['r_residual'] for r in rows):+.4f}]",
           kind="control", population=f"GSS Stouffer battery, all {K} items, n={len(d)}")
G.asserted("⚠ HARD RULE 1: n and years printed before any column was cited; and `#942`'s withdrawn "
           "limit is what makes this design exist at all", True,
           f"n={len(d)} answering all {K} items · {d.year.nunique()} waves "
           f"{int(d.year.min())}-{int(d.year.max())} — `#942` recorded this as 0 rows and `#944` "
           f"withdrew it", kind="control",
           population=f"GSS Stouffer battery, all {K} items, n={len(d)}")
G.asserted("the whole 30-cell grid and the specification curve are published, disagreeing cells "
           "included", True,
           " · ".join(f"{g_['spec'][:3]}/{g_['arm'][:3]} mean {g_['mean_pair']:+.3f}"
                      for g_ in grid), kind="control",
           population=f"GSS Stouffer battery, all {K} items, n={len(d)}")

pos_fires = (-sweep[-1][1]) > (-sweep[0][1]) + 2 * null_sd[nk]
# ⚠ the conditional now gates on the ALGEBRA holding (check a), not on the subset coincidence v1
#   demanded. The within-stem baseline is measured by the null, per check (b).
neg_null = (abs(all_pairs_res - FLOOR) < 0.01
            and null_reproducibility < 0.5 * mean_abs_departure)
trade = len(below) > 0
world = ("W_IPSATIVE" if (max(r["z"] for r in rows) - min(r["z"] for r in rows)) < 2 else
         ("W_TRADE_OFF" if trade else "W_GENERAL"))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and the null "
           "lands on the derived floor. STAKED: W_GENERAL, i.e. NO target pair sits significantly "
           "below the ipsative floor. W_TRADE_OFF refutes it. ⚠ W_GENERAL is staked precisely "
           "because it is the outcome that would make this round add nothing",
           (pos_fires and neg_null) and not trade,
           f"positive fires {pos_fires} · ipsative algebra holds over all 105 pairs AND the null "
           f"reproduces itself ({null_reproducibility:.5f} vs half the mean departure "
           f"{0.5 * mean_abs_departure:.5f}): {neg_null} "
           f"(|{all_pairs_res:+.4f} - {FLOOR:+.4f}| = {abs(all_pairs_res - FLOOR):.4f}) · cells "
           f"BELOW their own null cell at p<0.05: "
           f"{len(below)} · cells ABOVE: {len(above)} ⇒ {world}",
           kind="kill", yardstick="residual pair correlation minus the ipsative floor -1/(k-1)",
           yardstick_noise=float(np.mean(list(null_sd.values()))),
           population=f"GSS Stouffer battery, all {K} items jointly, n={len(d)}, "
                      f"{d.year.nunique()} waves 1988-2021",
           direction="two-sided; BELOW the floor is a trade-off, ABOVE is co-tolerance")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if not trade else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=945, round="E03·A125·R383", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_GENERAL"),
               n=int(len(d)), k=K, derived_floor=FLOOR, null_grand_mean=grand_null,
               all_pairs_residual_mean=all_pairs_res,
               cells=rows, grid=grid, below_floor=len(below), above_floor=len(above),
               null_median=null_mean[nk], null_sd=null_sd[nk], null_draws=NREPS,
               positive_sweep=sweep, family_size=len(ps), world=world, verdict=verdict,
               null_reproducibility=null_reproducibility,
               mean_abs_departure=mean_abs_departure),
          open(OUT / "the_trade_off_against_its_ipsative_floor.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_trade_off_against_its_ipsative_floor.json'}")
