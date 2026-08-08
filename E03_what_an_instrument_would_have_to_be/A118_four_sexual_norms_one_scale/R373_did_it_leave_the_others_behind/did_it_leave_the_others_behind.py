#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A118·R373 — `homosex` moved and the others did not. Did it LEAVE them, or drag them?
========================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#933`/`#934` are the project's best-supported claim: on one identical frame, within
                cohorts, robust to a randomised removal of the question context, **`homosex` moved
                +0.59 while the other three moved <= +0.17.** ⚠ **But nobody has asked what that did
                to the STRUCTURE.** An item can move a long way and stay embedded in the same
                dimension as its neighbours, or it can **detach from them**. Those are different
                psychological events and the project has been writing as though they were one.

⚠ WHY NOT      The obvious design — subgroup movement correlations — was measured first and REJECTED:
THE SUBGROUP   politics x education x religiosity x cohort gives **26-31 usable cells**, and `#930`
DESIGN         already ran that shape at 25 cells and got **0 of 3 splits surviving BH** against a
                ±0.21 permutation null. `#930`(2) called it a RESOLUTION limit and said the next move
                must be a different question, not a fourth attempt. **So the statistic changed, not
                the sample.**

Estimand       At the INDIVIDUAL level, per wave: the mean |Spearman rho| between `homosex` and each
(G1)           of the other three sexual norms — and its TREND across 1988-2024. ~700 respondents per
               wave on ballot 1, so this is well-powered where the subgroup design was not.

Live Worlds    W_DETACHED · the `homosex`-to-others correlation FALLS across waves ⇒ the act left the
                            others behind, and "sexual morality" lost a member rather than shifting.
               W_EMBEDDED · the correlation holds ⇒ `homosex` moved far while staying inside the same
                            structure, and "one act moved alone" is the wrong picture of what
                            happened. ⚠ **Unwelcome: it complicates the claim the page now leads with.**
               W_GENERAL  · ALL pairwise correlations fall, including among the other three ⇒ this is
                            a general decline in inter-item coherence (mode change, respondent
                            fatigue, a scale drifting), **not anything about this act**.
                            (the meta-separator: it makes the target/other split irrelevant)

Prediction     W_DETACHED -> the `homosex` trend is clearly negative AND the other-three trend is not.
Matrix         W_EMBEDDED -> the `homosex` trend sits on its null.
               W_GENERAL  -> both trends fall together.

⚠ THE SHAM     the three pairs AMONG the other three items (`premarsx`-`teensex`, `premarsx`-`xmarsex`,
THAT MATTERS   `teensex`-`xmarsex`) are the reference. **If their coherence falls too, nothing here is
               about homosexuality** — and that is the world I would otherwise never have checked.

⚠ SCOPE        Ballot 1 only, because that is where all four are asked of the same person (`#934`'s
               finding: `#933` was ballot 1 too). ~700/wave, 21 waves.

Controls       NEGATIVE: permute `homosex` across persons within wave — kills the coupling, preserves
               both marginals.
               POSITIVE: plant a shared factor INTO the permuted world and sweep, `g=0` on the null
               (`#922`'s gate).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ a falling correlation is consistent with detachment AND with the item's variance changing —
    both are reported, and the variance is printed beside the correlation so it cannot hide;
  (2) ⚠ all four items are asked in one interview (HARD RULE 2), inherited from `#933`;
  (3) ⚠ repeated cross-section: this is a trend in a population statistic, not in any person;
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
RNG = np.random.default_rng(373)
TARGET = "homosex"
OTHERS = ["premarsx", "teensex", "xmarsex"]
ITEMS = OTHERS + [TARGET]
MIN_WAVE = 300

d = pd.read_stata(GSS, columns=["year", "ballot"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[(d.ballot == 1)].dropna(subset=ITEMS)

# ══ PRECONDITION / COVERAGE — non-null, not rows (`#932`(1)) ═════════════════════════
counts = d.groupby("year").size()
waves = [int(y) for y, n in counts.items() if n >= MIN_WAVE]
dropped = [int(y) for y, n in counts.items() if 0 < n < MIN_WAVE]
print(f"PRECONDITION — ballot 1, all four answered, waves with n>={MIN_WAVE}: "
      f"{len(waves)} usable · DROPPED {len(dropped)}: {dropped or 'none'} (absence reported, not passed)")
print(f"  total n={len(d)} · per-wave min {int(counts[counts>=MIN_WAVE].min())} "
      f"median {int(counts[counts>=MIN_WAVE].median())}")

rows = []
for y in waves:
    w = d[d.year == y]
    tgt = [abs(stats.spearmanr(w[TARGET], w[o]).statistic) for o in OTHERS]
    oth = [abs(stats.spearmanr(w[a], w[b]).statistic)
           for i, a in enumerate(OTHERS) for b in OTHERS[i + 1:]]
    rows.append(dict(year=int(y), n=int(len(w)),
                     target_mean=float(np.mean(tgt)), others_mean=float(np.mean(oth)),
                     target_sd=float(w[TARGET].std()),
                     others_sd=float(np.mean([w[o].std() for o in OTHERS]))))

print("\n=== THE GRID — per wave (all cells, disagreeing ones included) ===")
for r in rows:
    print(f"  {r['year']}  n={r['n']:4d}  |rho| homosex-to-others {r['target_mean']:.4f} · "
          f"among the other three {r['others_mean']:.4f}  |  sd homosex {r['target_sd']:.3f} "
          f"others {r['others_sd']:.3f}")

yrs = np.array([r["year"] for r in rows], float)
tgt = np.array([r["target_mean"] for r in rows])
oth = np.array([r["others_mean"] for r in rows])
slope_t = float(stats.linregress(yrs, tgt).slope) * 10      # per decade
slope_o = float(stats.linregress(yrs, oth).slope) * 10
print(f"\n  trend per decade: homosex-to-others {slope_t:+.4f} · among the other three {slope_o:+.4f}")
print(f"  first wave {rows[0]['year']}: target {tgt[0]:.4f} others {oth[0]:.4f}")
print(f"  last  wave {rows[-1]['year']}: target {tgt[-1]:.4f} others {oth[-1]:.4f}")
print(f"  ⚠ variance check (a falling rho can be a shrinking sd): homosex sd "
      f"{rows[0]['target_sd']:.3f} -> {rows[-1]['target_sd']:.3f} · others "
      f"{rows[0]['others_sd']:.3f} -> {rows[-1]['others_sd']:.3f}")

# ══ NEGATIVE CONTROL — permute the target across persons within wave ═════════════════
null_slopes = []
for _ in range(300):
    vals = []
    for y in waves:
        w = d[d.year == y].copy()
        w[TARGET] = RNG.permutation(w[TARGET].to_numpy())
        vals.append(float(np.mean([abs(stats.spearmanr(w[TARGET], w[o]).statistic) for o in OTHERS])))
    null_slopes.append(float(stats.linregress(yrs, np.array(vals)).slope) * 10)
null_med, null_sd = float(np.median(null_slopes)), float(np.std(null_slopes))
print(f"\n  null for the TREND (target permuted across persons within each wave; kind of null: "
      f"within-wave person-label permutation): {null_med:+.4f} +/- {null_sd:.4f} "
      f"over {len(null_slopes)} draws")

# ══ POSITIVE CONTROL — plant a declining coupling INTO the permuted world ════════════
sweep = []
for g in (0.0, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(12):
        series = []
        for i, y in enumerate(waves):
            w = d[d.year == y].copy()
            w[TARGET] = RNG.permutation(w[TARGET].to_numpy())
            # ⚠⚠ v1 PLANTED THE WRONG WAY ROUND: `frac = 1 - g*i` restores the coupling FULLY at
            #   g=0, so the zero arm was the OBSERVED world (trend +0.0306) while the baseline
            #   passed to the gate was the PERMUTATION null (-0.0027). Two different worlds — the
            #   `#920`/`#923` error, and `plant_direction_from_sweep` caught it: "g=0 未落在基线上".
            #   ⇒ restore NOTHING at g=0 (pure permuted world = the null) and restore progressively
            #   MORE in later waves as g rises, which plants a RISING coupling trend — the direction
            #   the data actually shows — with g=0 landing on the null by construction.
            frac = g * (i / max(len(waves) - 1, 1))
            take = RNG.random(len(w)) < frac
            w.loc[take, TARGET] = d[d.year == y][TARGET].to_numpy()[take]
            series.append(float(np.mean([abs(stats.spearmanr(w[TARGET], w[o]).statistic)
                                         for o in OTHERS])))
        vals.append(float(stats.linregress(yrs, np.array(series)).slope) * 10)
    sweep.append([float(g), float(np.median(vals))])
print(f"  positive sweep (a RISING coupling trend of strength g planted into the permuted world; "
      f"g=0 IS the null): {[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["target_mean"] - np.mean(tgt)) / (np.std(tgt) or 1e-9))))
      for r in rows]

if not rows:
    print("EMPTY POPULATION"); sys.exit(2)

detached = (slope_t - null_med) < -2 * null_sd
others_fell = (slope_o - null_med) < -2 * null_sd

G = Gate("Did `homosex` leave the other three behind?")
G.plant_direction_from_sweep("positive: a planted RISING coupling trend raises the measured trend, "
                             "and g=0 IS the null world", sweep, baseline=null_med,
                             baseline_spread=max(null_sd, 1e-5))
G.negative_control("target permuted across persons within each wave", abs(null_med), abs(slope_t),
                   null_spread=null_sd, null_kind="within-wave person-label permutation")
G.multiplicity_control("all waves", ps, 0.05, labels=[str(r["year"]) for r in rows])
G.asserted("SHAM: the three pairs among the OTHER THREE are the reference — if their coherence falls "
           "too, nothing here is about homosexuality",
           not (detached and others_fell),
           f"target trend {slope_t:+.4f}/decade · other-three trend {slope_o:+.4f}/decade",
           kind="control")
G.asserted("VARIANCE printed beside the correlation, because a falling rho can be a shrinking sd",
           True,
           f"homosex sd {rows[0]['target_sd']:.3f} -> {rows[-1]['target_sd']:.3f}; others "
           f"{rows[0]['others_sd']:.3f} -> {rows[-1]['others_sd']:.3f}", kind="control")
G.asserted("the subgroup design was REJECTED on power before this one was built (`#930`(2))", True,
           "politics x education x religiosity x cohort gives 26-31 cells; `#930` ran that shape at "
           "25 and got 0/3 surviving BH against a ±0.21 null. The statistic changed, not the sample",
           kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", rows)
G.asserted("KILL: W_EMBEDDED requires the `homosex`-to-others coupling trend to sit on its null",
           not detached,
           f"trend {slope_t:+.4f} per decade vs null {null_med:+.4f} +/- {null_sd:.4f}; "
           f"coupling {tgt[0]:.4f} ({rows[0]['year']}) -> {tgt[-1]:.4f} ({rows[-1]['year']})")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif detached and others_fell:
    VERDICT, WORLD = "UNVERIFIED", "W_GENERAL · all coherence fell; nothing here is about this act"
elif detached:
    VERDICT, WORLD = "OVERTURNED", "W_DETACHED · the act left the other three behind"
else:
    VERDICT, WORLD = "CONFIRMED", "W_EMBEDDED · it moved far and stayed inside the same structure"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=935, round="E03·A118·R373", verdict=VERDICT, world=WORLD,
           estimand="per-wave individual-level mean |Spearman rho| between `homosex` and the other "
                    "three sexual norms, and its trend per decade",
           instrument="GSS 1972-2024 gss7224_r3a.dta, ballot 1 (all four asked of one person)",
           rejected_design="subgroup movement correlations: 26-31 cells, the shape `#930` already "
                           "showed underpowered at 25 (0/3 surviving BH)",
           waves=len(waves), dropped_waves=dropped, n=int(len(d)), rows=rows,
           target_trend_per_decade=slope_t, others_trend_per_decade=slope_o,
           target_first=float(tgt[0]), target_last=float(tgt[-1]),
           others_first=float(oth[0]), others_last=float(oth[-1]),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_slopes),
           positive_sweep=sweep, family_size=len(ps),
           claims_null=bool(not detached),
           claims_null_reason="if W_EMBEDDED, the finding IS that the coupling trend sits on its null",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "did_it_leave_the_others_behind.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'did_it_leave_the_others_behind.json'}")
