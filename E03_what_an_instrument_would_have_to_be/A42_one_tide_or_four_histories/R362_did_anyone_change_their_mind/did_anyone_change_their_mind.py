#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A114·R362 — `#923` said Americans changed their minds. Did they, or did they get replaced?
==============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#923` put a sentence on the page: *"over fifty years Americans changed their minds
                about homosexuality almost entirely without changing who they were."* ⚠ **"Who they
                were" there meant only same-sex experience.** There is a far more basic sense in
                which the population changed: **the people who held the old view died, and were
                replaced by people born later.** If replacement carries the trend, nobody changed
                their mind — the electorate changed — and `#923`'s sentence is false as written.

Why Now         `#923`(1) named the residual as the frontier: ~96% of the change is unexplained by
                anything measured here. Cohort replacement is the single largest candidate, and it
                is the one that would embarrass the sentence I wrote last round.

⚠ PRIOR ART,    **This is a VERIFICATION on this release, not a discovery.** The cohort/period
DECLARED        structure of this trend is established in the literature (Treas 2002; Andersen &
BEFORE THE      Fetner 2008; Baunach 2012), which reports BOTH cohort succession and substantial
NUMBER          intracohort change, with period effects strengthening after ~1990. A round that
                "found" this and did not say so would be restating known work as its own — the
                failure `prior_art_in_card` exists to catch (`P15`). What this round adds is the
                measured split ON THIS RELEASE, with controls, to constrain `#923`'s residual.

Live Worlds     W_REPLACE · cohort turnover carries most of it. **Unwelcome: it guts `#923`.**
                W_CONVERT · within-cohort movement carries most of it: people did change their minds.
                W_BOTH    · comparable. `#923`'s sentence needs qualifying, not retracting.
                ⚠ meta-separator: the split presumes a cohort observed 40 years apart is the SAME
                population. **Differential mortality and participation break that** — the 1930
                cohort seen in 2020 is its surviving, participating remnant. That biases TOWARD
                apparent within-cohort liberalisation, i.e. toward the comfortable answer.

Estimand        (a) `comp_cohort` — the Kitagawa composition share of the aggregate change using
(G1)                BIRTH DECADE as the group. ⚠ DERIVATION, forced by the algebra given shares and
                    group means (`realstat`'s arithmetic trap), labelled as one and given a floor;
                (b) `within_mean` — the average movement WITHIN cohorts observed at both endpoints.
                    Identified without any APC assumption, and it is the testable half.

Prediction      W_REPLACE -> comp_cohort large, within_mean small.
Matrix          W_CONVERT -> the reverse.
                W_BOTH    -> both substantial.

Strongest       ⚠ **MORTALITY / PARTICIPATION SELECTION**, above. CONTROL, same iteration: the whole
Confound        analysis is repeated restricted to cohorts **under 60 at both endpoints**, where
                selective attrition is small. If within-cohort movement survives there, it is not
                an artifact of who was left alive to answer.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠⚠ **AGE, PERIOD AND COHORT ARE EXACTLY COLLINEAR** (cohort = period - age), so a full
    three-way decomposition is **UNIDENTIFIED, permanently, by arithmetic**. This round therefore
    estimates only the two quantities that ARE identified — within-cohort change between two
    waves, and the composition share given cohort shares — and claims nothing about "the age
    effect". No amount of data fixes this; it is not a power problem;
  (2) ⚠ mortality selection can be BOUNDED by the age restriction, never removed;
  (3) ⚠ prior art, above: the qualitative answer is known; the contribution here is the measured
    split on this release with controls;
  (4) ⚠ **only this one instrument** — GSS is the only release here with a five-decade series;
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
RNG = np.random.default_rng(362)

df = pd.read_stata(GSS, columns=["homosex", "cohort", "age", "year"], convert_categoricals=False)
df["homosex"] = df["homosex"].where(df["homosex"].isin([1, 2, 3, 4]))
d = df.dropna(subset=["homosex", "cohort", "age"])
d = d[d["cohort"].between(1900, 2006)].copy()

print(f"homosex x cohort x age: n={len(d)} · {d.year.nunique()} waves · "
      f"{int(d.year.min())}-{int(d.year.max())} · cohorts {int(d.cohort.min())}-{int(d.cohort.max())}")


def kitagawa(e, gcol, y0, y1):
    """DERIVATION (labelled): total change = composition + within, over the groups in `gcol`."""
    a, b = e[e.year == y0], e[e.year == y1]
    if len(a) < 200 or len(b) < 200:
        return None
    ks = sorted(set(a[gcol]) & set(b[gcol]))
    ks = [k for k in ks if (a[gcol] == k).sum() >= 25 and (b[gcol] == k).sum() >= 25]
    if len(ks) < 3:
        return None
    p0 = {k: float((a[gcol] == k).mean()) for k in ks}
    p1 = {k: float((b[gcol] == k).mean()) for k in ks}
    m0 = {k: float(a.loc[a[gcol] == k, "homosex"].mean()) for k in ks}
    m1 = {k: float(b.loc[b[gcol] == k, "homosex"].mean()) for k in ks}
    mbar = {k: (m0[k] + m1[k]) / 2 for k in ks}
    pbar = {k: (p0[k] + p1[k]) / 2 for k in ks}
    comp = sum((p1[k] - p0[k]) * mbar[k] for k in ks)
    within = sum(pbar[k] * (m1[k] - m0[k]) for k in ks)
    total = comp + within
    moves = [m1[k] - m0[k] for k in ks]
    return dict(y0=int(y0), y1=int(y1), groups=len(ks), total=total, comp=comp, within=within,
                comp_share=(comp / total if total else np.nan),
                within_mean=float(np.mean(moves)), within_min=float(np.min(moves)),
                within_max=float(np.max(moves)),
                n=int(len(a) + len(b)), n0=int(len(a)), n1=int(len(b)))


def prep(e, width):
    e = e.copy()
    e["grp"] = (e["cohort"] // width * width).astype(int)
    return e


WINDOWS = [(1974, 2022), (1974, 1998), (1998, 2022), (1990, 2022)]
grid = []
for width in (10, 5):
    for agecap in (None, 60):
        e = d if agecap is None else d[d["age"] < agecap]
        e = prep(e, width)
        for y0, y1 in WINDOWS:
            yrs = sorted(e.year.unique())
            a0 = min(yrs, key=lambda y: abs(y - y0))
            a1 = min(yrs, key=lambda y: abs(y - y1))
            if a0 == a1:
                continue
            r = kitagawa(e, "grp", a0, a1)
            if r:
                r.update(cohort_width=width, age_cap=(agecap or "none"))
                grid.append(r)

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for r in grid:
    print(f"  width={r['cohort_width']:2d} agecap={str(r['age_cap']):4s} {r['y0']}->{r['y1']} "
          f"grp={r['groups']:2d}  total {r['total']:+.3f} = comp {r['comp']:+.3f} + within "
          f"{r['within']:+.3f} | comp_share {r['comp_share']:6.1%} · within_mean {r['within_mean']:+.3f} "
          f"[{r['within_min']:+.2f},{r['within_max']:+.2f}]  n={r['n']}")

comp_shares = [r["comp_share"] for r in grid if not np.isnan(r["comp_share"])]
med_comp = float(np.median(comp_shares))
med_within = float(np.median([r["within_mean"] for r in grid]))
capped = [r for r in grid if r["age_cap"] == 60]
med_within_capped = float(np.median([r["within_mean"] for r in capped])) if capped else np.nan
print(f"\n  median comp_share  {med_comp:.1%}  over {len(comp_shares)} cells")
print(f"  median within_mean {med_within:+.3f}")
print(f"  median within_mean, AGE-CAPPED (<60 at both ends, mortality selection small): "
      f"{med_within_capped:+.3f}   <- the confound control")

# ⚠⚠⚠ THE DECOMPOSITION IS STRUCTURALLY ILL-POSED HERE, AND THE RUN SAID SO BEFORE I DID.
#   Kitagawa needs the SAME groups in both waves. **Birth cohorts over a 48-year window do not
#   overlap**: the 1900s cohort is absent in 2022 and the 2000s cohort is absent in 1974, so the
#   estimator runs on a small intersection and its `total` is NOT the population change. Measured:
#   totals of -0.200, -0.371, -0.203, +0.231 against a true 1973->2024 move of +1.251, and
#   `comp_share` values of 115%, 218%, 283%, -162% — a ratio dividing by a near-zero, wrong total.
#   ⇒ **`comp_share` is RETRACTED from this round, as ill-posed rather than imprecise.** No sample
#   size fixes it; it is the wrong estimator for a group that does not persist.
#   ⇒ What survives is the half that needs no common-group assumption: **the movement WITHIN
#   cohorts that ARE present at both endpoints.** That is reported, with its own null below.

# ══ NEGATIVE CONTROL — cohort shares held at their first-wave values ═════════════════
e10 = prep(d, 10)
yrs = sorted(e10.year.unique())
y0, y1 = min(yrs, key=lambda y: abs(y - 1974)), min(yrs, key=lambda y: abs(y - 2022))
base_p = e10.loc[e10.year == y0, "grp"].value_counts(normalize=True)
null_vals = []
for _ in range(120):
    parts = []
    for y in (y0, y1):
        w = e10[e10.year == y]
        n = len(w)
        picks = []
        for k, p in base_p.items():
            pool = w[w.grp == k]
            if len(pool) == 0:
                pool = e10[e10.grp == k]
            if len(pool) == 0:
                continue
            picks.append(pool.sample(max(int(round(p * n)), 1), replace=True,
                                     random_state=int(RNG.integers(1e9))))
        if not picks:
            parts = []
            break
        parts.append(pd.concat(picks))
    if len(parts) != 2:
        continue
    r = kitagawa(pd.concat(parts), "grp", y0, y1)
    if r and not np.isnan(r["comp_share"]):
        null_vals.append(r["comp"])
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null for the ABSOLUTE composition component (cohort shares held at their first-wave "
      f"values; kind of null: share-fixed resampling): {null_med:+.4f} +/- {null_sd:.4f} "
      f"({len(null_vals)} draws)")

# ══ POSITIVE CONTROL — plant a PURE-REPLACEMENT trend; sweep the ABSOLUTE component ══
# ⚠ `#923`: do NOT sweep `comp_share` — freezing the within-group means makes it 1.0 by
#   construction, so the plant would disturb the very quantity the statistic normalises by.
sweep = []
m_first = e10[e10.year == y0].groupby("grp")["homosex"].mean()
for g_amt in (0.0, 0.10, 0.20, 0.30, 0.40):
    vals = []
    for _ in range(25):
        parts = []
        for y, extra in ((y0, 0.0), (y1, g_amt)):
            w = e10[e10.year == y].copy()
            w["homosex"] = w["grp"].map(m_first)          # freeze every cohort's mean
            w = w.dropna(subset=["homosex"])
            if len(w) < 200:
                parts = []
                break
            if extra > 0:                                  # shift weight toward LATER cohorts
                late = w[w.grp >= w.grp.median()]
                if len(late) < 20:
                    parts = []
                    break
                k = int(round(extra * len(w)))
                w = pd.concat([w, late.sample(k, replace=True,
                                              random_state=int(RNG.integers(1e9)))])
            parts.append(w)
        if len(parts) != 2:
            continue
        r = kitagawa(pd.concat(parts), "grp", y0, y1)
        if r:
            vals.append(r["comp"])
    sweep.append([float(g_amt), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (planted replacement, median ABSOLUTE comp): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

# ══ THE TESTABLE HALF — is within-cohort movement above its own resampling floor? ════
boot = []
for _ in range(400):
    parts = [e10[e10.year == y].sample(len(e10[e10.year == y]), replace=True,
                                       random_state=int(RNG.integers(1e9))) for y in (y0, y1)]
    r = kitagawa(pd.concat(parts), "grp", y0, y1)
    if r:
        boot.append(r["within_mean"])
w_lo, w_hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
print(f"  within-cohort mean movement {int(y0)}->{int(y1)}: {float(np.median(boot)):+.3f} "
      f"[95% bootstrap {w_lo:+.3f}, {w_hi:+.3f}]  <- identified without any APC assumption")

# ══ THE NULL FOR THE SURVIVING QUANTITY — permute YEAR within cohort ═════════════════
# destroys within-cohort change over time while preserving every cohort's own level and size.
wnull = []
for _ in range(150):
    p = e10[e10.year.isin([y0, y1])].copy()
    p["year"] = p.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
    r = kitagawa(p, "grp", y0, y1)
    if r:
        wnull.append(r["within_mean"])
wn_med, wn_sd = float(np.median(wnull)), float(np.std(wnull))
print(f"  null for WITHIN-cohort movement (year permuted within cohort; kind of null: "
      f"within-cohort year-label permutation): {wn_med:+.4f} +/- {wn_sd:.4f} ({len(wnull)} draws)")

# ══ POSITIVE CONTROL FOR THE SURVIVING QUANTITY — and `#922`'s gate is why it exists ═══
# ⚠⚠ The first version of this round persisted the RETRACTED composition sweep beside the
#   WITHIN-cohort null, and `#922`'s new gate BLOCKED THE COMMIT: g=0 at +0.9392 against a null of
#   -0.0022 +/- 0.0905, 10.4 spreads apart. It was right, and for the exact reason it was built:
#   **those are two different quantities.** The surviving claim had a negative control and NO
#   positive control of its own. Planting a within-cohort trend into the year-permuted (null) world
#   and sweeping it is what the claim actually needed — the gate found a real hole, not a
#   bookkeeping one.
wsweep = []
for g_amt in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(20):
        p2 = e10[e10.year.isin([y0, y1])].copy()
        p2["year"] = p2.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
        p2.loc[p2.year == y1, "homosex"] = np.clip(p2.loc[p2.year == y1, "homosex"] + g_amt, 1, 4)
        r = kitagawa(p2, "grp", y0, y1)
        if r:
            vals.append(r["within_mean"])
    wsweep.append([float(g_amt), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep for WITHIN-cohort movement (planted into the year-permuted null world): "
      f"{[(g, round(v, 4)) for g, v in wsweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["within_mean"] - wn_med) / (wn_sd or 1e-9)))) for r in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

replace_wins = med_comp > 0.60
convert_alive = w_lo > 0 and med_within_capped > 0.15

G = Gate("Did anyone change their mind, or did the people change?")
G.plant_direction_from_sweep("positive: a planted within-cohort trend raises within_mean, and "
                             "g=0 is null", wsweep, baseline=wn_med, baseline_spread=wn_sd)
G.asserted("the Kitagawa split is RETRACTED as ill-posed, not reported as a result",
           True, "birth cohorts do not overlap across a 48-year window, so the estimator's `total` "
                 "is not the population change: measured totals -0.200/-0.371/-0.203/+0.231 against "
                 "a true +1.251, giving comp_share 115%/218%/283%/-162%. Wrong estimator for a "
                 "group that does not persist; no n fixes it; scope stated", kind="control")
G.negative_control("within-cohort movement vs year permuted within cohort",
                   abs(wn_med), abs(float(np.median(boot))),
                   null_spread=wn_sd, null_kind="within-cohort year-label permutation")
G.multiplicity_control("the whole width x agecap x window grid", ps, 0.05,
                       labels=[f"w{r['cohort_width']}|cap{r['age_cap']}|{r['y0']}-{r['y1']}"
                               for r in grid])
G.asserted("mortality-selection control ran: the analysis is repeated on cohorts under 60 at both ends",
           len(capped) >= 3 and not np.isnan(med_within_capped),
           f"{len(capped)} age-capped cells · within_mean capped {med_within_capped:+.3f} vs "
           f"uncapped {med_within:+.3f}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.has_error_bar("within-cohort movement carries an interval", float(np.median(boot)),
                (w_hi - w_lo) / 4, "bootstrap_人层")
G.asserted("APC collinearity is registered as unidentified, not estimated",
           True, "cohort = period - age exactly; this round estimates only within-cohort change "
                 "and the composition share, and claims nothing about an age effect; scope stated",
           kind="control")
G.asserted("prior art is declared: the qualitative answer is known from the literature",
           True, "Treas 2002 · Andersen & Fetner 2008 · Baunach 2012 report both cohort succession "
                 "and intracohort change; this round contributes the measured split on THIS "
                 "release with controls, and is a VERIFICATION", kind="control")
G.asserted("KILL: W_REPLACE requires within-cohort movement to be indistinguishable from zero",
           w_lo > 0 and abs(float(np.median(boot)) - wn_med) > 2 * wn_sd,
           f"within-cohort {float(np.median(boot)):+.3f} [{w_lo:+.3f},{w_hi:+.3f}] vs null "
           f"{wn_med:+.4f} +/- {wn_sd:.4f}; age-capped {med_within_capped:+.3f}")

tv = G.three_valued()
converted = w_lo > 0 and abs(float(np.median(boot)) - wn_med) > 2 * wn_sd
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif converted:
    VERDICT, WORLD = ("OVERTURNED",
                      "W_CONVERT-ALIVE · people DID change their minds; the SPLIT stays unmeasured")
else:
    VERDICT, WORLD = "UNVERIFIED", "within-cohort movement does not resolve"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=924, round="E03·A114·R362", verdict=VERDICT, world=WORLD,
           estimand="composition share of the aggregate change using BIRTH DECADE as the group "
                    "(DERIVATION) plus the within-cohort mean movement (identified, testable)",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           prior_art="Treas 2002 · Andersen & Fetner 2008 · Baunach 2012 — VERIFICATION, not discovery",
           unidentified="age/period/cohort are exactly collinear; no age effect is claimed",
           comp_share_RETRACTED="ill-posed: birth cohorts do not overlap across the window",
           grid=grid, median_comp_share_ill_posed=med_comp, median_within=med_within,
           median_within_agecapped=med_within_capped,
           within_ci=dict(median=float(np.median(boot)), lo=w_lo, hi=w_hi, draws=len(boot)),
           null_median=wn_med, null_sd=wn_sd, null_draws=len(wnull),
           comp_null_median_ill_posed=null_med, comp_null_sd_ill_posed=null_sd,
           positive_sweep=wsweep,
           positive_sweep_composition_RETRACTED=sweep,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "did_anyone_change_their_mind.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'did_anyone_change_their_mind.json'}")
