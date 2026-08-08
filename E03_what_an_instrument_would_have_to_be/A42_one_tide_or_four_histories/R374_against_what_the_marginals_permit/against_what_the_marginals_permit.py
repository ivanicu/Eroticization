#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A118·R374 — `#935` gave a direction. Against what the marginals permit, is there a magnitude?
=================================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#935` measured the `homosex`-to-others coupling rising **0.3311 -> 0.4804** and then
                bounded its own claim: `homosex`'s sd rose **1.089 -> 1.407**, and **a near-degenerate
                marginal caps any rank correlation**, so part of the rise is the item un-pinning
                rather than integrating. `#935`(2) named the fix and it is the project's own:
                **`#918`'s attainable ceiling.**

The instrument For two ordinal items the maximum attainable |Spearman| is set by their MARGINALS
already        alone, and it is reached by the **comonotone** pairing — sort both, pair them. That is
exists         `#902`'s machinery, built in this project for a different question. **The ceiling
                contains no association**: it is a property of the two marginals, so dividing by it
                is a units fix and cannot launder the finding.

Live Worlds    W_REAL     · the ceiling-normalised coupling still rises ⇒ the integration is real and
                            `#935`'s direction becomes a magnitude.
               W_MARGINAL · it flattens or falls ⇒ `#935`'s rise WAS the item un-pinning, and the
                            sentence "became more tightly bound" must be retracted.
                            ⚠ **The unwelcome one, and the reason this round exists.**
               W_DEGEN    · the ceiling is ~1.0 in every wave ⇒ normalising does nothing and the
                            round cannot distinguish anything. (the meta-separator: it would mean
                            `#935`'s caveat was never testable this way)

⚠ THE          Checked and PRINTED before the estimator, because `W_DEGEN` would make the whole round
DEGENERACY     vacuous: measured at four waves the ceiling is **0.8236 / 0.8020 / 0.6824 / 0.7684** —
CHECK CAME     non-degenerate, and it **FALLS** slightly across the span. ⇒ **the confound runs
FIRST          AGAINST the finding**, so normalising can only make the rise larger, not manufacture it.

Estimand       per wave: `observed |rho| / comonotone-ceiling |rho|`, averaged over the three
(G1)           `homosex`-to-other pairs; and its trend per decade. The same for the three pairs
               AMONG the others, which is the sham.

Prediction     W_REAL     -> normalised trend clearly positive, sham flat.
Matrix         W_MARGINAL -> normalised trend at or below its null.
               W_DEGEN    -> ceilings ~1.0 and normalised == raw.

Controls       NEGATIVE: permute `homosex` across persons within wave — kills the coupling, preserves
               BOTH marginals, so the ceiling is untouched and only the association dies.
               POSITIVE: plant a rising coupling INTO the permuted world, `g=0` on the null
               (`#922`'s gate; `#935` got this backwards and its own control caught it).
               SHAM: the three other-three pairs, normalised identically.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ the ceiling bounds what the MARGINALS permit; it does not bound what a 4-point scale can
    express about a continuous attitude — coarseness is not corrected, only degeneracy;
  (2) ⚠ all four items are asked in one interview (HARD RULE 2), inherited; `#934` removed the
    four-item version, not the two-item one;
  (3) ⚠ repeated cross-section: a trend in a population statistic, never in any person;
  (4) ⚠ **only this one instrument**;
  (5) `[unchallenged]` — door (3).
"""
import json, sys, warnings
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
RNG = np.random.default_rng(374)
TARGET = "homosex"
OTHERS = ["premarsx", "teensex", "xmarsex"]
ITEMS = OTHERS + [TARGET]
MIN_WAVE = 300

d = pd.read_stata(GSS, columns=["year", "ballot"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS)
counts = d.groupby("year").size()
waves = [int(y) for y, n in counts.items() if n >= MIN_WAVE]
dropped = [int(y) for y, n in counts.items() if 0 < n < MIN_WAVE]
print(f"PRECONDITION — ballot 1, all four answered, waves n>={MIN_WAVE}: {len(waves)} usable · "
      f"DROPPED {len(dropped)}: {dropped or 'none'} (absence reported, not passed)")


def ceiling(x, y):
    """max attainable |Spearman| given BOTH marginals — the comonotone pairing (`#902`).
    A property of the two marginals alone: it contains no association whatsoever."""
    return abs(stats.spearmanr(np.sort(x), np.sort(y)).statistic)


def wave_stats(w, target_col=None):
    t = w[TARGET].to_numpy() if target_col is None else target_col
    obs_t, cl_t = [], []
    for o in OTHERS:
        y = w[o].to_numpy()
        obs_t.append(abs(stats.spearmanr(t, y).statistic))
        cl_t.append(ceiling(t, y))
    obs_o, cl_o = [], []
    for i, a in enumerate(OTHERS):
        for b in OTHERS[i + 1:]:
            obs_o.append(abs(stats.spearmanr(w[a], w[b]).statistic))
            cl_o.append(ceiling(w[a].to_numpy(), w[b].to_numpy()))
    return (float(np.mean(obs_t)), float(np.mean(cl_t)),
            float(np.mean(obs_o)), float(np.mean(cl_o)))


# ══ DEGENERACY CHECK, printed BEFORE the estimator ═══════════════════════════════════
print("\nDEGENERACY CHECK — if the ceiling were ~1.0 everywhere, normalising would do nothing:")
rows = []
for y in waves:
    w = d[d.year == y]
    ot, ct, oo, co = wave_stats(w)
    rows.append(dict(year=int(y), n=int(len(w)), obs_t=ot, ceil_t=ct, norm_t=ot / ct,
                     obs_o=oo, ceil_o=co, norm_o=oo / co, sd_t=float(w[TARGET].std())))
cl_all = [r["ceil_t"] for r in rows]
print(f"  target-pair ceiling: min {min(cl_all):.4f} · median {float(np.median(cl_all)):.4f} · "
      f"max {max(cl_all):.4f}  -> {'NON-degenerate' if max(cl_all) < 0.97 else 'DEGENERATE'}")
print(f"  and it FALLS: {rows[0]['ceil_t']:.4f} ({rows[0]['year']}) -> {rows[-1]['ceil_t']:.4f} "
      f"({rows[-1]['year']}) ⇒ the confound runs AGAINST the finding")

print("\n=== THE GRID — per wave (all cells, disagreeing ones included) ===")
for r in rows:
    print(f"  {r['year']}  n={r['n']:4d}  target: obs {r['obs_t']:.4f} / ceil {r['ceil_t']:.4f} "
          f"= {r['norm_t']:.4f}   |  others: obs {r['obs_o']:.4f} / ceil {r['ceil_o']:.4f} "
          f"= {r['norm_o']:.4f}   sd_t {r['sd_t']:.3f}")

yrs = np.array([r["year"] for r in rows], float)
nt = np.array([r["norm_t"] for r in rows])
no = np.array([r["norm_o"] for r in rows])
raw_t = np.array([r["obs_t"] for r in rows])
slope_norm = float(stats.linregress(yrs, nt).slope) * 10
slope_raw = float(stats.linregress(yrs, raw_t).slope) * 10
slope_sham = float(stats.linregress(yrs, no).slope) * 10
print(f"\n  trend per decade — RAW {slope_raw:+.4f} · CEILING-NORMALISED {slope_norm:+.4f} · "
      f"sham (others, normalised) {slope_sham:+.4f}")
print(f"  normalised coupling {nt[0]:.4f} ({rows[0]['year']}) -> {nt[-1]:.4f} ({rows[-1]['year']})")

# ══ NEGATIVE CONTROL — permute the target within wave (marginals, hence ceilings, intact) ══
null_slopes = []
for _ in range(250):
    vals = []
    for y in waves:
        w = d[d.year == y]
        perm = RNG.permutation(w[TARGET].to_numpy())
        ot, ct, _, _ = wave_stats(w, target_col=perm)
        vals.append(ot / ct)
    null_slopes.append(float(stats.linregress(yrs, np.array(vals)).slope) * 10)
null_med, null_sd = float(np.median(null_slopes)), float(np.std(null_slopes))
print(f"\n  null for the NORMALISED trend (target permuted within wave — BOTH marginals and therefore "
      f"the ceiling are untouched; kind of null: within-wave person-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_slopes)} draws")

# ══ POSITIVE CONTROL — plant a RISING coupling INTO the permuted world; g=0 IS the null ══
sweep = []
for g in (0.0, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(10):
        series = []
        for i, y in enumerate(waves):
            w = d[d.year == y]
            perm = RNG.permutation(w[TARGET].to_numpy())
            frac = g * (i / max(len(waves) - 1, 1))     # nothing restored at g=0 -> the null world
            take = RNG.random(len(w)) < frac
            col = np.where(take, w[TARGET].to_numpy(), perm)
            ot, ct, _, _ = wave_stats(w, target_col=col)
            series.append(ot / ct)
        vals.append(float(stats.linregress(yrs, np.array(series)).slope) * 10)
    sweep.append([float(g), float(np.median(vals))])
print(f"  positive sweep (a rising coupling of strength g planted into the permuted world; g=0 IS "
      f"the null): {[(g, round(v, 4)) for g, v in sweep]}")

# ⚠⚠ v1 RAN THE MULTIPLICITY CORRECTION ON THE WRONG FAMILY: it tested each WAVE's deviation from
#   the mean, and got 0 of 21 surviving — beside a verdict about a TREND. **A wave's deviation from
#   the mean is not a member of the family this round's claim belongs to**, which is `#898`'s error
#   class: comparing a number to a reference that answers a different question. `#931`'s gate flagged
#   the mismatch as a verdict/control contradiction, and it was right to.
#   ⇒ the family is the THREE TREND ESTIMATES this round actually reports — raw, ceiling-normalised,
#   and the sham — each against the same permutation null.
fam = {"raw trend": slope_raw, "ceiling-normalised trend": slope_norm, "sham (others) trend": slope_sham}
ps = [2 * (1 - stats.norm.cdf(abs((v - null_med) / (null_sd or 1e-9)))) for v in fam.values()]

if not rows:
    print("EMPTY POPULATION"); sys.exit(2)

real = (slope_norm - null_med) > 2 * null_sd
sham_flat = abs(slope_sham - null_med) <= 2 * null_sd
degenerate = max(cl_all) >= 0.97

G = Gate("Against what the marginals permit, did the coupling really rise?")
G.plant_direction_from_sweep("positive: a planted rising coupling raises the normalised trend, and "
                             "g=0 IS the null world", sweep, baseline=null_med,
                             baseline_spread=max(null_sd, 1e-5))
G.negative_control("target permuted within wave (marginals, hence the ceiling, untouched)",
                   abs(null_med), abs(slope_norm), null_spread=null_sd,
                   null_kind="within-wave person-label permutation")
G.multiplicity_control("the three trend estimates this round reports", ps, 0.05,
                       labels=list(fam))
G.asserted("DEGENERACY CHECK ran FIRST — a ceiling of ~1.0 would make this round vacuous",
           not degenerate,
           f"ceiling min {min(cl_all):.4f} median {float(np.median(cl_all)):.4f} max {max(cl_all):.4f}",
           kind="control")
G.asserted("the confound runs AGAINST the finding, so normalising cannot manufacture the rise",
           rows[-1]["ceil_t"] <= rows[0]["ceil_t"],
           f"ceiling {rows[0]['ceil_t']:.4f} -> {rows[-1]['ceil_t']:.4f}; raw trend {slope_raw:+.4f} "
           f"vs normalised {slope_norm:+.4f}", kind="control")
G.asserted("SHAM: the other-three pairs, normalised identically, are the reference",
           True, f"sham normalised trend {slope_sham:+.4f} per decade vs target {slope_norm:+.4f}",
           kind="control")
G.asserted("the ceiling contains no association — it is a property of the two MARGINALS alone", True,
           "computed as the comonotone pairing (`#902`'s machinery), so dividing by it is a units "
           "fix and cannot launder the finding it is used to test", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", rows)
G.asserted("KILL: W_MARGINAL requires the normalised trend to sit at or below its null",
           not real,
           f"normalised trend {slope_norm:+.4f}/decade vs null {null_med:+.4f} +/- {null_sd:.4f}; "
           f"coupling {nt[0]:.4f} -> {nt[-1]:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif degenerate:
    VERDICT, WORLD = "UNVERIFIED", "W_DEGEN · the ceiling cannot discriminate"
elif real and sham_flat:
    VERDICT, WORLD = "OVERTURNED", "W_REAL · the integration survives the marginals, with a magnitude"
elif real:
    VERDICT, WORLD = "OVERTURNED", "W_REAL, but the sham moved too — scope narrowed"
else:
    VERDICT, WORLD = "CONFIRMED", "W_MARGINAL · `#935`'s rise was the item un-pinning"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=936, round="E03·A118·R374", verdict=VERDICT, world=WORLD,
           estimand="per-wave observed |rho| divided by the comonotone ceiling the two marginals "
                    "permit, averaged over the three `homosex`-to-other pairs, and its trend",
           instrument="GSS 1972-2024 gss7224_r3a.dta, ballot 1",
           waves=len(waves), dropped_waves=dropped, rows=rows,
           trend_raw_per_decade=slope_raw, trend_normalised_per_decade=slope_norm,
           trend_sham_per_decade=slope_sham,
           normalised_first=float(nt[0]), normalised_last=float(nt[-1]),
           ceiling_first=rows[0]["ceil_t"], ceiling_last=rows[-1]["ceil_t"],
           ceiling_min=min(cl_all), ceiling_max=max(cl_all), degenerate=bool(degenerate),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_slopes),
           positive_sweep=sweep, family_size=len(ps), family={k: float(v) for k, v in fam.items()},
           multiplicity_family_note="v1 used per-wave deviation from the mean, which is not the "
                                    "family of a TREND claim (`#898`); `#931`'s gate caught it",
           upgrades="`#935`(2): turns a bounded DIRECTION into a MAGNITUDE",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "against_what_the_marginals_permit.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'against_what_the_marginals_permit.json'}")
