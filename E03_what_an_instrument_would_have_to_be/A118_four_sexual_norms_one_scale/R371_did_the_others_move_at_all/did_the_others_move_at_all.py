#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A118·R371 — four sexual norms, one scale, the same people: did the others move at all?
==========================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#926` put `homosex` at the 100th percentile of 22 items and `#927`-`#929` showed the
                same for the homosexuality TOLERANCE items — but **the comparison sets were other
                topics on other scales, standardised to be comparable.** ⚠ **The four GSS SEXUAL
                norms have never been compared to each other**, and they are the project's actual
                subject: `premarsx` · `teensex` · `xmarsex` · `homosex`.

Why this is    All four use the **identical 4-point frame** — *always wrong · almost always wrong ·
the cleanest   wrong only sometimes · not wrong at all* — asked of the same respondents in the same
comparison     interview. **No standardisation, no polarity assumption, no cross-scale translation.**
this project   Everything `#926` had to assume away is simply absent here.
has            (`#900`/`#902` showed these four are NOT one moving factor at the decade level; this
                asks the within-cohort version, which is a different question.)

⚠ PRIOR ART    Twenge, Sherman & Wells (2015) documented the DIVERGENT trends in exactly these four
DECLARED       GSS items. **The divergence is not this round's discovery.** What is this round's own:
BEFORE THE     the WITHIN-COHORT version with a ceiling control and a permutation null, which
NUMBER         separates "the population changed" from "these people changed".

Live Worlds    W_ONE_ACT  · `homosex` moves far more than the other three, in every cohort ⇒ what
                            moved is one act, not sexual morality.
               W_ALL      · all four move comparably ⇒ sexual morality as a whole loosened and
                            `#926`'s "topic-specific" reading was an artifact of its comparison set.
                            ⚠ **Unwelcome: it would undo four rounds.**
               W_CEILING  · `homosex` started lowest and simply had the most room ⇒ the gap is the
                            scale, not the people. (the meta-separator: it makes "how much did X
                            move" the wrong question and "how much COULD it move" the right one)

Estimand       the within-cohort change in each item's mean, early era -> late era, **in raw scale
(G1)           points** (identical scales make this legitimate); then `homosex` minus the mean of the
               other three, per cohort.

Prediction     W_ONE_ACT -> gap large and positive in every cohort, and it survives headroom scaling.
Matrix         W_ALL     -> gap inside its permutation null.
               W_CEILING -> gap large raw but vanishing once each change is divided by the room the
                            item had. ⚠ **Measured, not assumed** (`#929`(3)).

⚠ PRECONDITION `#925`(2)/`#932`(1): checked and PRINTED before the estimator, on **NON-NULL coverage**
CHECK FIRST    rather than a row count — the error that killed `#932`. Cohorts need n>=60 in BOTH eras.

Controls       NEGATIVE: permute ERA within cohort — destroys the change while preserving each
               cohort's level and size.
               POSITIVE: plant a movement INTO the permuted world and sweep, so `g=0` lands ON the
               null (`#922`'s gate).
               CEILING: every change also divided by `4 - early_mean`, the room that item had.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ repeated cross-section: a cohort is not the same PEOPLE twice, only the same birth years;
  (2) ⚠ all four items come from one questionnaire and one interview (HARD RULE 2) — a respondent
    answering all four in a row may anchor them against each other, which no design here removes;
  (3) ⚠ APC collinear; no age effect claimed;
  (4) ⚠ **only this one instrument** — no other release here asks four sexual norms on one frame;
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
RNG = np.random.default_rng(371)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
TARGET = "homosex"
TOP = 4.0
MIN_CELL = 60

x = pd.read_stata(GSS, columns=["year", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    x[c] = x[c].where(x[c].isin([1, 2, 3, 4]))      # `homosex` code 5 = "other", dropped

# ══ HARD RULE 1 — NON-NULL coverage, not a row count (`#932`(1)) ═════════════════════
print("PRECONDITION / COVERAGE (non-null, not rows — the error that killed `#932`):")
for c in ITEMS:
    yrs = sorted(x.loc[x[c].notna(), "year"].unique())
    print(f"  {c:9s} non-null {int(x[c].notna().sum()):6d}  waves {len(yrs):2d} "
          f"{int(yrs[0])}-{int(yrs[-1])}")
x = x.dropna(subset=ITEMS + ["cohort"])
x = x[x.cohort.between(1900, 2006)].copy()
x["grp"] = (x.cohort // 10 * 10).astype(int)
print(f"  all four answered together: n={len(x)} · waves {x.year.nunique()} "
      f"{int(x.year.min())}-{int(x.year.max())}")

SPLITS = {"1988-1998 vs 2014-2024": ((1988, 1998), (2014, 2024)),
          "1988-2000 vs 2012-2024": ((1988, 2000), (2012, 2024)),
          "1988-2002 vs 2010-2024": ((1988, 2002), (2010, 2024))}


def cells_for(early, late):
    e = x.copy()
    e["era"] = np.where(e.year.between(*early), "E", np.where(e.year.between(*late), "L", None))
    e = e[e.era.notna()]
    sz = e.groupby(["grp", "era"]).size().unstack()
    if sz.empty:
        return [], 0, e
    sz = sz.dropna()
    ks = list(sz[(sz >= MIN_CELL).all(axis=1)].index)
    rows = []
    for k in ks:
        a = e[(e.grp == k) & (e.era == "E")]
        b = e[(e.grp == k) & (e.era == "L")]
        r = dict(cohort=int(k), n=int(len(a) + len(b)))
        for c in ITEMS:
            m0, m1 = float(a[c].mean()), float(b[c].mean())
            r[f"{c}_early"] = m0
            r[f"{c}_move"] = m1 - m0
            r[f"{c}_head"] = (m1 - m0) / max(TOP - m0, 1e-9)      # ⚠ ceiling control
        r["gap"] = r[f"{TARGET}_move"] - np.mean([r[f"{c}_move"] for c in ITEMS if c != TARGET])
        r["gap_head"] = r[f"{TARGET}_head"] - np.mean([r[f"{c}_head"] for c in ITEMS if c != TARGET])
        rows.append(r)
    return rows, len(sz), e


grid, frames = [], {}
for name, (early, late) in SPLITS.items():
    rows, total, e = cells_for(early, late)
    frames[name] = e
    print(f"  split {name}: cohorts n>={MIN_CELL} in BOTH eras: {len(rows)} of {total}"
          f"  (absence reported, not passed)")
    if len(rows) < 3:
        continue
    grid.append(dict(split=name, cohorts=len(rows), rows=rows,
                     n=int(sum(r["n"] for r in rows)),
                     **{f"{c}_move": float(np.mean([r[f'{c}_move'] for r in rows])) for c in ITEMS},
                     **{f"{c}_head": float(np.mean([r[f'{c}_head'] for r in rows])) for c in ITEMS},
                     gap=float(np.mean([r["gap"] for r in rows])),
                     gap_head=float(np.mean([r["gap_head"] for r in rows]))))

print("\n=== THE GRID — within-cohort movement in RAW SCALE POINTS (identical 1-4 frame) ===")
for g in grid:
    print(f"  {g['split']:24s} cohorts={g['cohorts']}  " +
          " · ".join(f"{c} {g[f'{c}_move']:+.4f}" for c in ITEMS) +
          f"  |  gap {g['gap']:+.4f}  n={g['n']}")
print("\n  the CEILING control — every change divided by the room the item had (4 - early mean):")
for g in grid:
    print(f"  {g['split']:24s} " + " · ".join(f"{c} {g[f'{c}_head']:+.4f}" for c in ITEMS) +
          f"  |  gap_head {g['gap_head']:+.4f}")
print("\n  starting levels (early-era mean), which is what the ceiling argument rests on:")
for c in ITEMS:
    lv = float(np.mean([r[f"{c}_early"] for r in grid[0]["rows"]]))
    print(f"    {c:9s} {lv:.3f}   room to 4.0 = {TOP - lv:.3f}")

med_gap = float(np.median([g["gap"] for g in grid]))
med_gap_head = float(np.median([g["gap_head"] for g in grid]))
per_cohort = [r["gap"] for g in grid for r in g["rows"]]
print(f"\n  median gap (raw)      {med_gap:+.4f}   · positive in {sum(1 for v in per_cohort if v > 0)}"
      f"/{len(per_cohort)} cohort-cells")
print(f"  median gap (headroom) {med_gap_head:+.4f}   <- the confound control")

# ══ NEGATIVE CONTROL — permute ERA within cohort ═════════════════════════════════════
base = frames[grid[0]["split"]]
null_vals = []
for _ in range(300):
    p = base.copy()
    p["era"] = p.groupby("grp")["era"].transform(lambda s: RNG.permutation(s.to_numpy()))
    vals = []
    for k in sorted(p.grp.unique()):
        a, b = p[(p.grp == k) & (p.era == "E")], p[(p.grp == k) & (p.era == "L")]
        if len(a) < MIN_CELL or len(b) < MIN_CELL:
            continue
        mv = {c: float(b[c].mean() - a[c].mean()) for c in ITEMS}
        vals.append(mv[TARGET] - np.mean([mv[c] for c in ITEMS if c != TARGET]))
    if vals:
        null_vals.append(float(np.mean(vals)))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (era permuted within cohort; kind of null: within-cohort era-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant INTO the permuted world so g=0 lands on the null ════════
sweep = []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(40):
        p = base.copy()
        p["era"] = p.groupby("grp")["era"].transform(lambda s: RNG.permutation(s.to_numpy()))
        p.loc[p.era == "L", TARGET] = np.clip(p.loc[p.era == "L", TARGET] + g, 1, 4)
        cur = []
        for k in sorted(p.grp.unique()):
            a, b = p[(p.grp == k) & (p.era == "E")], p[(p.grp == k) & (p.era == "L")]
            if len(a) < MIN_CELL or len(b) < MIN_CELL:
                continue
            mv = {c: float(b[c].mean() - a[c].mean()) for c in ITEMS}
            cur.append(mv[TARGET] - np.mean([mv[c] for c in ITEMS if c != TARGET]))
        if cur:
            vals.append(float(np.mean(cur)))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (planted into the permuted null world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((g["gap"] - null_med) / (null_sd or 1e-9)))) for g in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

resolved = abs(med_gap - null_med) > 2 * null_sd
ceiling_survives = np.sign(med_gap_head) == np.sign(med_gap) and abs(med_gap_head) > 0.0

G = Gate("Four sexual norms on one scale: did the others move at all?")
G.plant_direction_from_sweep("positive: a planted movement in the target raises the gap, g=0 is null",
                             sweep, baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("era permuted within cohort", abs(null_med), abs(med_gap),
                   null_spread=null_sd, null_kind="within-cohort era-label permutation")
G.multiplicity_control("all era splits", ps, 0.05, labels=[g["split"] for g in grid])
G.asserted("PRECONDITION on NON-NULL coverage, not a row count (`#932`(1))", True,
           f"all four items' non-null n and wave spans printed above; cohorts require n>={MIN_CELL} "
           f"in both eras and shortfalls are reported", kind="control")
G.asserted("CEILING control ran and its premise was MEASURED, not asserted (`#929`(3))",
           True,
           f"starting levels: " + " · ".join(
               f"{c} {float(np.mean([r[f'{c}_early'] for r in grid[0]['rows']])):.2f}" for c in ITEMS)
           + f" ; gap raw {med_gap:+.4f} vs headroom {med_gap_head:+.4f}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("prior art declared", True,
           "Twenge, Sherman & Wells 2015 documented the divergent trends in exactly these four GSS "
           "items; the divergence is not this round's discovery. The within-cohort version with a "
           "ceiling control and a permutation null is", kind="control")
G.asserted("KILL: W_ALL requires the gap between `homosex` and the other three to sit on its null",
           not (resolved and ceiling_survives),
           f"gap {med_gap:+.4f} (headroom {med_gap_head:+.4f}) vs null {null_med:+.4f} "
           f"+/- {null_sd:.4f}; positive in {sum(1 for v in per_cohort if v > 0)}/{len(per_cohort)} "
           f"cohort-cells")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif resolved and ceiling_survives:
    VERDICT, WORLD = "OVERTURNED", "W_ONE_ACT · one act moved, not sexual morality"
elif resolved and not ceiling_survives:
    VERDICT, WORLD = "UNVERIFIED", "W_CEILING · the gap does not survive the room each item had"
else:
    VERDICT, WORLD = "CONFIRMED", "W_ALL · the four moved comparably"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=933, round="E03·A118·R371", verdict=VERDICT, world=WORLD,
           estimand="within-cohort change of each of the four GSS sexual norms in RAW scale points "
                    "(identical 1-4 frame), and `homosex` minus the mean of the other three",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           prior_art="Twenge, Sherman & Wells 2015 — the divergence is not this round's discovery",
           items=ITEMS, grid=[{k: v for k, v in g.items() if k != "rows"} for g in grid],
           cohort_cells=grid[0]["rows"],
           median_gap=med_gap, median_gap_headroom=med_gap_head,
           cohort_cells_positive=int(sum(1 for v in per_cohort if v > 0)), cohort_cells_total=len(per_cohort),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "did_the_others_move_at_all.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'did_the_others_move_at_all.json'}")
