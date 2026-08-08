#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A121·R377 — the same decade, lived at four different ages
==============================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        A114 established that **at least 45% of the half-century change is conversion, not
                replacement** — living people changed their minds. A118 established **what** they
                changed it about, and that the change folded further INTO their sexual morality
                rather than out of it. **Nothing yet says WHEN IN A LIFE it happens.** That is the
                next question about people, and it is the one a reader asks immediately: *do minds
                change while they are still forming, or does the era carry everyone along at once?*

Why now         Two consecutive rounds (`#937`, `#938`) were about instruments, both UNVERIFIED.
                `#111c` says change direction rather than chase a third. This is the object.

Live Worlds    W_ERA      · every age band moves at the same rate. ⇒ the change is a property of the
                            SOCIETY and the age at which you meet it is irrelevant — nobody is too
                            old, and being young confers no special plasticity.
               W_PLASTIC  · movement concentrates in the young. ⇒ the change is generational
                            **within** living people too, and A114's "conversion" is mostly young
                            adults converting, which narrows it sharply.
               W_LATE     · the old move as much or more. ⚠ **The unwelcome one** — it contradicts
                            the impressionable-years account that every version of this literature
                            assumes, including the one I would write.
               W_UNIDENT  · the age-band difference is composition, not age. (the meta-separator:
                            the decomposition "age vs era" may not be a thing this data can carve)

⚠⚠ THE IDENTIFICATION WALL, WRITTEN BEFORE THE ESTIMAND AND NOT AFTER (`#924`'s lesson: reaching
   for a standard estimator without checking its precondition cost a whole column).
   **age = year − cohort, exactly.** So age, period and cohort main effects are collinear and **no
   amount of data separates them.** What this round therefore does NOT do is estimate an age effect.
   What IS estimable, and is a fact about people either way: **the period slope computed separately
   within each age band.** Whether a difference between those slopes is caused by age or by which
   cohorts happen to occupy that band is EXACTLY the confound, and it gets its control below, in the
   same iteration.

Estimand       For each of four age bands (18–29, 30–44, 45–59, 60+), the OLS slope of `homosex`
(G1)           (reversed, so higher = more permissive) on `year`, 1988–2024, GSS ballots carrying
               the item. Reported as change per decade, with the four slopes' spread as the
               quantity of interest. Population: 30,748 respondents with age, cohort and `homosex`.

Prediction     W_ERA     -> the four slopes agree within their resampling spread.
Matrix         W_PLASTIC -> 18–29 steepest, 60+ flattest, monotone.
               W_LATE    -> 60+ at least as steep as 18–29.
               W_UNIDENT -> the cohort-composition control removes the difference.

Strongest      **COHORT COMPOSITION.** The 18–29 band in 1988 holds people born ~1959–1970; in 2024
confound       it holds people born ~1995–2006. If cohorts differ in level AND the bands turn over,
(written       a slope difference appears with no age effect whatever. ⇒ CONTROL, same iteration:
before)        re-estimate every band's slope with **cohort fixed effects**, i.e. demean `homosex`
               within birth-year before taking the slope, so only WITHIN-cohort movement is left.
               ⚠ A second, subtler one: `homosex` is 4-point and its marginal moved a lot over 36
               years, so a slope in raw units is capped differently at different levels — the
               ceiling logic of `#936`. Reported beside the raw slope.

Controls       NEGATIVE: permute `year` WITHIN age band — bands, marginals and n untouched, only the
                 time ordering dies.
               POSITIVE: plant a differential trend into the PERMUTED world (so `g=0` sits on the
                 null, `#922`/`#937`⑤ — built backwards three times, so it is built from the null
                 here by default) and sweep.
               MULTIPLICITY: the family is **the four band slopes**, which is the family the claim
                 lives in (`#936`②).
               SPEC CURVE (G4): 4 age-band cuts × {raw, cohort-demeaned} × {all ballots, ballot 1
                 only} — every cell published, including disagreeing ones.

Stopping rule  Pre-registered CONDITIONAL kill. If positive fires and negative is null:
                 max−min slope spread > 2× its own bootstrap spread ⇒ the bands differ (W_PLASTIC
                 or W_LATE, read off the ordering); ≤ that ⇒ W_ERA. And if the cohort-demeaned
                 spread falls below half the raw spread ⇒ W_UNIDENT regardless of the above.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate age from period from cohort** — algebraically impossible, see the wall above;
    this round reports band-specific period slopes and never an age effect;
  (2) ⚠ **follow the same person over time** — GSS is repeated cross-section, so "changed their
    mind" is always an inference about groups, never an observation of a person (A116 is the
    nearest this project gets, and it established a sign only);
  (3) ⚠ **no second instrument for the TREND** — `#937` measured that NSFG has one time point;
    only this one instrument carries a 36-year series of this item;
  (4) `[unchallenged]` — door ③.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
GSS = ROOT / "data" / "external" / "gss" / "GSS_stata" / "gss7224_r3a.dta"
RNG = np.random.default_rng(377)

BANDS = [("18-29", 18, 29), ("30-44", 30, 44), ("45-59", 45, 59), ("60+", 60, 89)]

d = pd.read_stata(GSS, columns=["year", "age", "cohort", "ballot", "homosex"],
                  convert_categoricals=False)
d = d[d.homosex.isin([1, 2, 3, 4]) & d.year.between(1988, 2024)].dropna(
    subset=["age", "cohort"]).copy()
# ⚠ HARD RULE 1 — A VARIABLE NAME IS NOT A MEASUREMENT, and v1 got this backwards. GSS `homosex`
#   is ALREADY 1 "always wrong" -> 4 "not wrong at all", so higher is ALREADY more permissive; the
#   `5 - homosex` reversal v1 applied made every slope negative under a comment that said the
#   opposite. Nothing downstream would have caught it: the magnitudes and the controls are identical
#   under a sign flip, and only the sentence about people inverts. Used raw, matching A118.
d["perm"] = d.homosex.astype(float)            # 1 always wrong .. 4 not wrong at all
d["band"] = pd.cut(d.age, [17, 29, 44, 59, 89], labels=[b[0] for b in BANDS])
d = d.dropna(subset=["band"])
print(f"n={len(d)} · years {int(d.year.min())}-{int(d.year.max())} · "
      f"bands {dict(d.band.value_counts())}")


def slope_per_decade(f, col="perm"):
    if f[col].nunique() < 2 or f.year.nunique() < 3:
        return np.nan
    return float(stats.linregress(f.year, f[col]).slope * 10)


def demean_by_cohort(f):
    g = f.copy()
    g["perm"] = g.perm - g.groupby("cohort").perm.transform("mean")
    return g


def band_slopes(f, demean=False):
    g = demean_by_cohort(f) if demean else f
    return {b: slope_per_decade(g[g.band == b]) for b, _, _ in BANDS}


raw = band_slopes(d)
dem = band_slopes(d, demean=True)
spread_raw = float(np.nanmax(list(raw.values())) - np.nanmin(list(raw.values())))
spread_dem = float(np.nanmax(list(dem.values())) - np.nanmin(list(dem.values())))
print("\n  per decade, on a 4-point scale (higher = more permissive)")
for b, _, _ in BANDS:
    print(f"    {b:>6s}  raw {raw[b]:+.4f}   cohort-demeaned {dem[b]:+.4f}")
print(f"    spread max-min: raw {spread_raw:.4f} · cohort-demeaned {spread_dem:.4f}")

# ══ resampling spread of the SPREAD ══════════════════════════════════════════════════
boot = []
for _ in range(300):
    s = d.sample(len(d), replace=True, random_state=int(RNG.integers(1e9)))
    v = list(band_slopes(s).values())
    boot.append(np.nanmax(v) - np.nanmin(v))
boot_sd = float(np.std(boot))
print(f"  bootstrap spread OF that spread: {boot_sd:.4f}")

# ══ NEGATIVE CONTROL — permute year within band ══════════════════════════════════════
null = []
for _ in range(300):
    s = d.copy()
    s["year"] = s.groupby("band", observed=True).year.transform(
        lambda x: RNG.permutation(x.to_numpy()))
    v = list(band_slopes(s).values())
    null.append(np.nanmax(v) - np.nanmin(v))
null_med, null_sd = float(np.median(null)), float(np.std(null))
null_p95 = float(np.percentile(null, 95))
print(f"  null (year permuted WITHIN band — bands, marginals and n untouched; kind of null: "
      f"within-band year-label permutation): {null_med:+.4f} +/- {null_sd:.4f}")

# ══ POSITIVE CONTROL — plant a differential trend INTO the null; g=0 IS the null ═════
yc = d.year - d.year.mean()
sweep = []
for gg in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(30):
        s = d.copy()
        s["year"] = s.groupby("band", observed=True).year.transform(
            lambda x: RNG.permutation(x.to_numpy()))
        # plant: only the youngest band gains a trend, in the permuted (null) world
        bump = gg * (s.year - s.year.mean()) / 18.0
        s.loc[s.band == "18-29", "perm"] = s.loc[s.band == "18-29", "perm"] + \
            bump[s.band == "18-29"]
        v = list(band_slopes(s).values())
        vals.append(np.nanmax(v) - np.nanmin(v))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a differential trend planted into the YOUNGEST band, inside the permuted "
      f"world so g=0 IS the null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs null {null_med:+.4f} +/- {null_sd:.4f}"
      f" = {abs(sweep[0][1] - null_med) / max(null_sd, 1e-9):.2f} spreads")

# ══ SPECIFICATION CURVE (G4) — every cell, including the disagreeing ones ════════════
grid = []
for demean in (False, True):
    for tag, sub in (("all ballots", d), ("ballot 1 only", d[d.ballot == 1])):
        sl = band_slopes(sub, demean=demean)
        vals = [sl[b] for b, _, _ in BANDS]
        grid.append(dict(spec=f"{'cohort-demeaned' if demean else 'raw'} · {tag}",
                         n=int(len(sub)), **{b: sl[b] for b, _, _ in BANDS},
                         spread=float(np.nanmax(vals) - np.nanmin(vals)),
                         youngest_steepest=bool(sl["18-29"] == np.nanmax(vals)),
                         oldest_steepest=bool(sl["60+"] == np.nanmax(vals))))
print("\n  specification curve — all four cells, none dropped")
for r in grid:
    print(f"    {r['spec']:34s} n={r['n']:6d}  " +
          "  ".join(f"{b} {r[b]:+.3f}" for b, _, _ in BANDS) + f"   spread {r['spread']:.4f}")

ps = [2 * (1 - stats.norm.cdf(abs(raw[b] / max(boot_sd, 1e-9)))) for b, _, _ in BANDS]

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G = Gate("When in a life does a mind change — or does the era carry every age at once?")
G.plant_direction_from_sweep("positive: a trend planted into the youngest band widens the spread, "
                             "and g=0 sits ON the null this round judges against (`#922`)",
                             sweep, baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("year permuted within age band", abs(null_med), abs(spread_raw),
                   null_spread=null_sd, null_kind="within-band year-label permutation")
G.multiplicity_control("the four band slopes — the family this claim lives in (`#936`②)",
                       ps, 0.05, labels=[b for b, _, _ in BANDS])
G.asserted("⚠ the identification wall is registered, not worked around: age = year - cohort exactly, "
           "so NO age effect is estimated here", True,
           "what is reported is the period slope computed WITHIN each age band, and the cohort "
           "composition of those bands is the confound, controlled in the next row", kind="control",
           population="30,748 GSS respondents 1988-2024 with age, cohort and homosex")
G.asserted("⚠ CONFOUND CONTROL in the same iteration: cohort-demeaned slopes remove every "
           "between-cohort level difference", True,
           f"spread raw {spread_raw:.4f} -> cohort-demeaned {spread_dem:.4f} "
           f"({spread_dem / max(spread_raw, 1e-9):.2f}x); W_UNIDENT is declared below 0.50x",
           kind="control",
           population="30,748 GSS respondents 1988-2024 with age, cohort and homosex")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{r['spec']} spread {r['spread']:.4f}" for r in grid), kind="control",
           population="30,748 GSS respondents 1988-2024 with age, cohort and homosex")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * null_sd
# ⚠⚠ v1 WROTE `neg_null = abs(null_med) < 2*null_sd` AND THAT CONDITION CAN NEVER BE MET. `max-min`
#   of four noisy slopes is POSITIVE-DEFINITE: its permutation null is centred at +0.0275, not at 0,
#   so the conditional gated on it could never be evaluated and the kill could never fire. That is
#   `realstat`'s "control that cannot PASS" -- the mirror of `#938`'s kill that cannot fail, written
#   one round after measuring the corpus for exactly that family. **For a positive-definite
#   statistic the null is a DISTRIBUTION, not a zero**, and the precondition is that the observed
#   value clears its upper tail.
neg_null = spread_raw > null_p95
unident_by_ratio = spread_dem < 0.50 * spread_raw
differ = spread_raw > 2 * boot_sd and spread_raw > null_p95
order = sorted(BANDS, key=lambda b: -(raw[b[0]] if not np.isnan(raw[b[0]]) else -9))
order_dem = sorted(BANDS, key=lambda b: -(dem[b[0]] if not np.isnan(dem[b[0]]) else -9))
# ⚠⚠ AND THE PRE-REGISTERED QUANTITY WAS NOT THE QUANTITY THE CLAIM LIVES IN (`#936`(2), one level
#   up). I staked a RATIO OF SPREADS, but W_PLASTIC vs W_LATE is decided by the ORDERING, and the
#   ordering REVERSES under the confound control: raw puts 18-29 steepest, cohort-demeaned puts it
#   FLATTEST. A ratio of 0.56 sails past a 0.50 threshold while the entire psychological claim has
#   already changed sign. The ordering check is therefore evaluated too, and it is decisive.
order_stable = order[0][0] == order_dem[0][0]
world = ("W_UNIDENT" if (unident_by_ratio or not order_stable) else
         ("W_ERA" if not differ else
          ("W_PLASTIC" if order[0][0] == "18-29" else
           "W_LATE" if order[0][0] == "60+" else "W_MIXED")))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and the negative "
           "is null. STAKED: W_PLASTIC, i.e. the 18-29 band is steepest AND the spread exceeds 2x "
           "its own bootstrap spread. W_ERA and W_LATE both REFUTE it; and a cohort-demeaned spread "
           "below 0.50x the raw ⇒ W_UNIDENT; AND the steepest band must be the SAME under both "
           "specifications, because that ordering is where the claim actually lives",
           # ⚠⚠ v1 WROTE `(differ or not differ)` HERE — a tautology, i.e. the very defect R376
           #   (`#938`) had just spent a whole round measuring, written one round later by the
           #   author of that round. Caught before the first run, which is the only reason it is a
           #   note and not an entry. **The kill must name the outcome I EXPECT and be able to lose
           #   it**: the impressionable-years account predicts W_PLASTIC, so that is what is staked.
           (pos_fires and neg_null) and differ and order_stable and order[0][0] == "18-29",
           f"positive fires {pos_fires} · observed clears the null's 95th pct ({null_p95:.4f}) "
           f"{neg_null} · ORDERING raw {[b for b, _, _ in order][:1]} vs cohort-demeaned "
           f"{[b for b, _, _ in order_dem][:1]} -> {'stable' if order_stable else 'REVERSES'} · "
           f"spread {spread_raw:.4f} vs "
           f"2x bootstrap {2 * boot_sd:.4f} -> {'differ' if differ else 'indistinguishable'} · "
           f"demeaned/raw {spread_dem / max(spread_raw, 1e-9):.2f}x -> "
           f"{'W_UNIDENT' if unident_by_ratio else 'ratio ok'} · ordering "
           f"{[b for b, _, _ in order]} ⇒ {world}",
           kind="kill", yardstick="max-min of the four band slopes, per decade",
           # the staked claim is W_PLASTIC (young steepest AND separated); W_ERA and W_LATE both
           # refute it, and W_LATE is the one I would find hardest to write
           yardstick_noise=boot_sd,
           population="30,748 GSS respondents 1988-2024 with age, cohort and homosex",
           direction="two-sided; the ordering, not the sign, selects between W_PLASTIC and W_LATE")

print(G)
tv = str(G)
staked = (pos_fires and neg_null) and differ and order_stable and order[0][0] == "18-29"
verdict = (f"{'UNVERIFIED' if (not pos_fires or not neg_null) else ('CONFIRMED' if staked else 'OVERTURNED')}"
           f" · world {world} · staked W_PLASTIC {'held' if staked else 'FAILED'}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=939, round="E03·A121·R377", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_ERA"),
               n=int(len(d)), raw_slopes=raw, demeaned_slopes=dem,
               spread_raw=spread_raw, spread_demeaned=spread_dem, boot_sd=boot_sd,
               null_median=null_med, null_sd=null_sd, null_p95=null_p95, null_draws=len(null),
               order_raw=[b for b, _, _ in order], order_demeaned=[b for b, _, _ in order_dem],
               order_stable=bool(order_stable),
               positive_sweep=sweep, grid=grid, family_size=len(ps),
               world=world, verdict=verdict),
          open(OUT / "the_same_decade_lived_at_different_ages.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_same_decade_lived_at_different_ages.json'}")
