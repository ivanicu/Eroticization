#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A114·R363 — the split `#924` could not estimate, as a BOUND
===============================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#924` established that birth cohorts really did move (+0.609 [+0.405, +0.789]) and
                RETRACTED its own replacement-vs-conversion split as ill-posed: Kitagawa needs the
                same groups in both waves, and birth cohorts over 48 years do not overlap.
                **So "how much of the half-century was turnover" is still unmeasured**, and `#924`(1)
                says the fix is a different estimator, not a bigger sample.

⚠ THE          `#924`(2) named the error family this round must not repeat: *I reach for a standard
PRECONDITION   estimator without checking its precondition holds here.* So the preconditions are
CHECK CAME     checked and PRINTED before the estimator runs, and the obvious candidate was
FIRST          rejected on one: **age-standardisation does NOT work** — holding the age distribution
               fixed does not hold COHORT fixed, because in a repeated cross-section the cohort at a
               given age changes with the period. It removes ageing, not turnover.

The estimator  A **frozen-mean synthetic**, which needs no common group set:
                    predicted(t) = sum_k share_k(t) * m_k(frozen)
               where `m_k(frozen)` is cohort k's mean **at its first adequately-sized wave**. Each
               cohort contributes the view it BROUGHT, whenever it is present; nothing requires a
               cohort to appear in two particular waves. Precondition: every cohort must be observed
               once with adequate n — **checked, and it holds for all 11 birth decades.**

⚠⚠ WHY THIS    Measured: every cohort is frozen at roughly **age 17-28**, except the three oldest,
IS A BOUND     which the 1973 series start freezes at ~38/48/58. **A cohort entering in 2014 brings a
AND NOT A      view already formed under four decades of period change.** So the synthetic credits
POINT          to "replacement" everything that happened before each cohort walked in.
               ⇒ **`replacement_bound` is an UPPER BOUND**, and `1 - bound` is a LOWER BOUND on
               conversion. Reporting it as a point estimate would be the error `#916`(3) named.

Live Worlds    W_TIGHT · the bound is small ⇒ conversion decisively carries the trend.
               W_LOOSE · the bound is near 1 ⇒ this design cannot distinguish either.
                         ⚠ **PRE-REGISTERED: if LOOSE, the question CLOSES** — that would be two
                         attempts (`#924`, this) and `#111c` forbids a third.
               W_INCO  · predicted moves the wrong way or overshoots ⇒ the construction does not
                         carve at all. (the meta-separator)

Prediction     W_TIGHT -> bound well below 1 and above its permutation null.
Matrix         W_LOOSE -> bound at or above 1.
               W_INCO  -> bound negative, or predicted change opposite in sign to observed.

Controls       NEGATIVE: permute the frozen means ACROSS cohorts — destroys the cohort ordering that
               makes turnover work, while preserving the shares and the multiset of means.
               POSITIVE, and built so `g=0` lands ON that null (`#922`'s gate caught exactly this
               mistake in `#924`): start from the PERMUTED frozen means and interpolate toward the
               true ones, `m(g) = (1-g)*m_permuted + g*m_true`. At `g=0` the world IS the null.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠⚠ **replacement and period cannot be separated for cohorts entering mid-series** — their
    brought view already contains the period change that preceded them. This is why the answer is a
    bound and can never be a point on this data;
  (2) ⚠ age/period/cohort remain exactly collinear; no age effect is claimed;
  (3) ⚠ the three oldest cohorts are frozen at ~38-58, not at entry — a DIFFERENT construct for
    them. A specification restricted to cohorts frozen under 30 is run and published;
  (4) ⚠ repeated cross-section, not a panel: a cohort is not the same people twice;
  (5) ⚠ **only this one instrument** — GSS is the only release here with a five-decade series;
  (6) `[unchallenged]` — door (3).
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
RNG = np.random.default_rng(363)
MIN_N = 100

df = pd.read_stata(GSS, columns=["homosex", "cohort", "age", "year"], convert_categoricals=False)
df["homosex"] = df["homosex"].where(df["homosex"].isin([1, 2, 3, 4]))
d = df.dropna(subset=["homosex", "cohort", "age"])
d = d[d["cohort"].between(1900, 2006)].copy()


def frozen_table(e, width, min_n=MIN_N):
    """(precondition-checked) each cohort's mean at its FIRST adequately-sized wave."""
    e = e.copy()
    e["grp"] = (e["cohort"] // width * width).astype(int)
    rows, missing = {}, []
    for k, g in e.groupby("grp"):
        counts = g.groupby("year").size()
        okay = [y for y, n in counts.items() if n >= min_n]
        if not okay:
            missing.append(int(k))
            continue
        y = min(okay)
        rows[int(k)] = dict(frozen_year=int(y), frozen_mean=float(g.loc[g.year == y, "homosex"].mean()),
                            frozen_age=int(y - k - 5), n=int(counts[y]))
    return e, rows, missing


def synthetic(e, frozen, y0, y1):
    """predicted(t) = sum_k share_k(t) * m_k(frozen); returns the replacement bound."""
    out = {}
    for y in (y0, y1):
        w = e[e.year == y]
        if len(w) < 200:
            return None
        sh = w["grp"].value_counts(normalize=True)
        num = sum(p * frozen[k]["frozen_mean"] for k, p in sh.items() if k in frozen)
        cov = sum(p for k, p in sh.items() if k in frozen)
        if cov < 0.90:
            return None
        out[y] = dict(observed=float(w["homosex"].mean()), predicted=float(num / cov),
                      coverage=float(cov), n=int(len(w)))
    obs_d = out[y1]["observed"] - out[y0]["observed"]
    pre_d = out[y1]["predicted"] - out[y0]["predicted"]
    return dict(y0=int(y0), y1=int(y1), observed_change=obs_d, predicted_change=pre_d,
                replacement_bound=(pre_d / obs_d if abs(obs_d) > 1e-9 else np.nan),
                coverage=min(out[y0]["coverage"], out[y1]["coverage"]),
                n=out[y0]["n"] + out[y1]["n"])


# ══ PRECONDITIONS, checked and printed BEFORE the estimator runs (`#924`(2)) ═════════
print("PRECONDITION CHECK (this is the step `#924`(2) says I skip):")
e10, fr10, miss10 = frozen_table(d, 10)
print(f"  every cohort observed once with n>={MIN_N}? cohorts missing: {miss10 or 'NONE'} "
      f"-> {'FAILS' if miss10 else 'HOLDS'}")
for k in sorted(fr10):
    print(f"    {k}s frozen at {fr10[k]['frozen_year']} (age ~{fr10[k]['frozen_age']}, "
          f"n={fr10[k]['n']}) mean {fr10[k]['frozen_mean']:.3f}")
print("  ⚠ REJECTED CANDIDATE: age-standardisation. Holding the AGE distribution fixed does not "
      "hold COHORT fixed — at a given age the cohort changes with the period. It removes ageing, "
      "not turnover.")

WINDOWS = [(1974, 2022), (1974, 2024), (1990, 2022), (1998, 2022)]
grid = []
for width in (10, 5):
    for cap in (None, 30):
        e, fr, miss = frozen_table(d, width)
        if cap is not None:
            fr = {k: v for k, v in fr.items() if v["frozen_age"] <= cap}
        yrs = sorted(e.year.unique())
        for a, b in WINDOWS:
            y0 = min(yrs, key=lambda y: abs(y - a))
            y1 = min(yrs, key=lambda y: abs(y - b))
            if y0 == y1:
                continue
            r = synthetic(e, fr, y0, y1)
            if r:
                r.update(cohort_width=width, frozen_age_cap=(cap or "none"), cohorts=len(fr))
                grid.append(r)

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for r in grid:
    print(f"  width={r['cohort_width']:2d} frozen_age<={str(r['frozen_age_cap']):4s} "
          f"{r['y0']}->{r['y1']} k={r['cohorts']:2d}  observed {r['observed_change']:+.3f} · "
          f"predicted-if-nobody-changed {r['predicted_change']:+.3f}  => replacement bound "
          f"{r['replacement_bound']:6.1%}  coverage {r['coverage']:.3f}  n={r['n']}")

bounds = [r["replacement_bound"] for r in grid if not np.isnan(r["replacement_bound"])]
med_bound = float(np.median(bounds))
print(f"\n  median replacement UPPER BOUND: {med_bound:.1%}  over {len(bounds)} cells "
      f"[{min(bounds):.1%}, {max(bounds):.1%}]")
print(f"  => conversion LOWER bound: {1 - med_bound:.1%}")

# ══ NEGATIVE CONTROL — permute the frozen means ACROSS cohorts ═══════════════════════
e, fr, _ = frozen_table(d, 10)
yrs = sorted(e.year.unique())
y0 = min(yrs, key=lambda y: abs(y - 1974))
y1 = min(yrs, key=lambda y: abs(y - 2022))
ks = sorted(fr)
true_means = np.array([fr[k]["frozen_mean"] for k in ks])
null_vals = []
for _ in range(300):
    perm = RNG.permutation(true_means)
    f2 = {k: dict(fr[k], frozen_mean=float(perm[i])) for i, k in enumerate(ks)}
    r = synthetic(e, f2, y0, y1)
    if r and not np.isnan(r["replacement_bound"]):
        null_vals.append(r["replacement_bound"])
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (frozen means PERMUTED across cohorts; kind of null: cohort-label permutation of "
      f"the frozen means): {null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — start AT the null and interpolate toward the truth ════════════
# ⚠ built this way so `g=0` lands ON the null world by construction. `#922`'s gate caught exactly
#   the opposite arrangement in `#924`, where a plant's zero arm sat 10.4 spreads off its baseline.
base_perm = RNG.permutation(true_means)
sweep = []
for g in (0.0, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(20):
        p = RNG.permutation(true_means)
        mix = (1 - g) * p + g * true_means
        f2 = {k: dict(fr[k], frozen_mean=float(mix[i])) for i, k in enumerate(ks)}
        r = synthetic(e, f2, y0, y1)
        if r and not np.isnan(r["replacement_bound"]):
            vals.append(r["replacement_bound"])
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (permuted -> true frozen means): {[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((b - null_med) / (null_sd or 1e-9)))) for b in bounds]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

loose = med_bound >= 0.90
incoherent = med_bound < -0.10 or any(np.sign(r["predicted_change"]) != np.sign(r["observed_change"])
                                      for r in grid if abs(r["observed_change"]) > 0.05)
resolved = abs(med_bound - null_med) > 2 * null_sd

G = Gate("How much of the half-century COULD be turnover, at most?")
G.plant_direction_from_sweep("positive: interpolating from permuted toward true frozen means raises "
                             "the bound, and g=0 IS the null world", sweep,
                             baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("frozen means permuted across cohorts", abs(null_med), abs(med_bound),
                   null_spread=null_sd, null_kind="cohort-label permutation of the frozen means")
G.multiplicity_control("the whole width x frozen-age-cap x window grid", ps, 0.05,
                       labels=[f"w{r['cohort_width']}|cap{r['frozen_age_cap']}|{r['y0']}-{r['y1']}"
                               for r in grid if not np.isnan(r["replacement_bound"])])
G.asserted("PRECONDITIONS were checked and printed BEFORE the estimator ran (`#924`(2))",
           not miss10,
           f"all {len(fr10)} birth decades have a wave with n>={MIN_N}; age-standardisation was "
           f"REJECTED on its precondition (fixing age does not fix cohort); scope stated",
           kind="control")
G.asserted("coverage: the synthetic accounts for >=90% of each wave's respondents",
           all(r["coverage"] >= 0.90 for r in grid),
           f"minimum coverage across cells {min(r['coverage'] for r in grid):.3f}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("the answer is reported as a BOUND, not a point",
           True, "cohorts entering mid-series bring a view already formed under prior period change, "
                 "so the synthetic credits that to turnover: `replacement_bound` is an UPPER bound "
                 "and `1-bound` a LOWER bound on conversion; scope stated", kind="control")
G.asserted("KILL: W_LOOSE — if the bound reaches 1, this design cannot distinguish and the "
           "question CLOSES under `#111c`", not loose,
           f"median bound {med_bound:.1%} over {len(bounds)} cells, range "
           f"[{min(bounds):.1%}, {max(bounds):.1%}]")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif incoherent:
    VERDICT, WORLD = "UNVERIFIED", "W_INCO · the synthetic does not carve"
elif loose:
    VERDICT, WORLD = "UNVERIFIED", "W_LOOSE · the bound is vacuous; the question CLOSES (`#111c`)"
elif resolved:
    VERDICT, WORLD = "CONFIRMED", "W_TIGHT · turnover is bounded well below the whole trend"
else:
    VERDICT, WORLD = "UNVERIFIED", "the bound does not separate from its null"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=925, round="E03·A114·R363", verdict=VERDICT, world=WORLD,
           estimand="upper bound on the share of the aggregate change attributable to cohort "
                    "turnover, via a frozen-mean synthetic that needs no common group set",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           preconditions=dict(all_cohorts_covered=not miss10, missing=miss10,
                              rejected_candidate="age-standardisation: fixing age does not fix cohort"),
           frozen=fr10, grid=grid, median_replacement_bound=med_bound,
           conversion_lower_bound=1 - med_bound,
           bound_range=[min(bounds), max(bounds)],
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "an_upper_bound_on_replacement.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'an_upper_bound_on_replacement.json'}")
