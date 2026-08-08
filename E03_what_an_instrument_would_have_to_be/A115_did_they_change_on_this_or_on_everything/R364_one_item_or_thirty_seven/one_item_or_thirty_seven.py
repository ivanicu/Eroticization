#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A115·R364 — did those cohorts change their minds about THIS, or about everything?
=====================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#925`(1) closed A114 with an IMPOSSIBILITY CLAIM: *nothing in `data/external/`
                records why anyone changed.* ⚠ **This project has been burned by exactly that twice**
                — `#912` declared five cells unreadable and `#913` opened four; `#918`(2) said the
                disclosure test needed an acquisition and `#919` found GSS ships `mode`. `#913`(3)
                named the class. So the claim gets checked, not inherited.

                And it is checkable, because "why" has one measurable form here: **`#924` measured
                that cohorts moved +0.609 on `homosex`. Did those same cohorts move on `homosex`
                SPECIFICALLY, or on everything at once?**

Why Now         If everything moved together, the project's subject dissolves at the decade unit the
                way `#905` dissolved it at the person unit — and A115 is where that gets decided.

⚠ THREE-WAY,   A ·  **`homosex`** — the moral wrongness of the act
NOT TWO        B ·  **same-topic TOLERANCE** — `spkhomo` `colhomo` `libhomo` (Stouffer civil-
                     liberties items about homosexuals: a different construct, same topic)
               C ·  **everything else** — the other long-series GSS attitude items (death penalty,
                     communists, racists, abortion, marijuana, suicide, school prayer, women in
                     politics, police force, ...)

Live Worlds    W_GENERAL · A ~ B ~ C. Everything moved together; there is no specific object, and
                           the decade unit repeats `#905`'s person-level verdict.
               W_TOPIC   · A ~ B >> C. Attitudes ABOUT HOMOSEXUALITY moved far more than everything
                           else — the object is the TOPIC, not "sexual morality".
               W_MORAL   · A >> B ~ C. The moral-wrongness judgement moved uniquely, more than even
                           same-topic tolerance.
               ⚠ Three distinct signatures, no flat row: no single outcome satisfies two worlds.

Estimand       For each item, the **within-cohort movement between two waves, in units of that
(G1)           item's own pooled SD** — scale-free, so binary and ordinal items are comparable.
               Then: `homosex`'s PERCENTILE in the distribution of |movement| across all items.
               (`#905` ran this design at the person level for one-factor structure; this is its
               decade/cohort analogue, and that lineage is declared.)

Prediction     W_GENERAL -> `homosex` near the middle of C.
Matrix         W_TOPIC   -> A and B both in the top decile of C.
               W_MORAL   -> A in the top decile, B not.

⚠ PRECONDITION `#925`(2): the check runs and PRINTS before the estimator, and items failing it are
CHECK FIRST    DROPPED WITH A COUNT, never silently included — each item needs adequate n in both
               endpoint waves within the cohorts used. **Absence is reported, not passed.**

Controls       NEGATIVE: permute YEAR within cohort, per item — destroys within-cohort change while
               preserving each cohort's level and size.
               POSITIVE: plant a within-cohort trend INTO the permuted world and sweep, so `g=0`
               lands ON the null (`#922`'s gate caught the opposite arrangement in `#924`).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ this measures WHETHER the change was specific, never WHY it happened — no mechanism is
    identified, and `#925`(1)'s claim is only being tested in its measurable form;
  (2) ⚠ items differ in wording, scale and salience; standardising by each item's own SD makes them
    comparable in units, NOT in meaning;
  (3) ⚠ the tolerance items end in 2021 and some others in 2021/2022 — the window is the common one,
    and items absent at an endpoint are dropped and counted;
  (4) ⚠ repeated cross-section, not a panel; age/period/cohort remain collinear, no age effect claimed;
  (5) ⚠ **only this one instrument** — no other release here carries 30+ waves of attitude items;
  (6) `[unchallenged]` — door (3).
"""
import json, re, sys, warnings
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
RNG = np.random.default_rng(364)

TARGET = "homosex"
TOPIC = ["spkhomo", "colhomo", "libhomo"]
ITEMS = ["cappun", "gunlaw", "spkath", "colath", "libath", "spkcom", "colcom", "libcom",
         "homosex", "polhitok", "polattak", "grass", "abany", "spkrac", "colrac", "librac",
         "spkhomo", "colhomo", "libhomo", "aged", "fepol", "eqwlth", "prayer", "letdie1",
         "suicide1", "suicide2", "suicide4", "spkmil", "colmil", "libmil"]

raw = pd.read_stata(GSS, columns=["year", "cohort", "age"] + ITEMS, convert_categoricals=False)
raw = raw[raw["cohort"].between(1900, 2006)].copy()
raw["grp"] = (raw["cohort"] // 10 * 10).astype(int)
for c in ITEMS:
    v = raw[c]
    raw[c] = v.where(v.between(1, 9) & (v != 8) & (v != 9))     # 8/9 = DK/NA on these scales
raw.loc[:, TARGET] = raw[TARGET].where(raw[TARGET].isin([1, 2, 3, 4]))

# ══ PRECONDITION CHECK, printed BEFORE the estimator (`#925`(2)) ═════════════════════
WINDOW = (1990, 2021)
yrs = sorted(raw.year.unique())
y0 = min(yrs, key=lambda y: abs(y - WINDOW[0]))
y1 = min(yrs, key=lambda y: abs(y - WINDOW[1]))
print(f"PRECONDITION CHECK — window {int(y0)}->{int(y1)}; an item is USABLE only if >=3 cohorts have "
      f"n>=40 in BOTH endpoint waves:")
usable, dropped = [], []
for c in ITEMS:
    e = raw.dropna(subset=[c])
    a, b = e[e.year == y0], e[e.year == y1]
    ks = [k for k in sorted(set(a.grp) & set(b.grp))
          if (a.grp == k).sum() >= 40 and (b.grp == k).sum() >= 40]
    (usable if len(ks) >= 3 else dropped).append((c, len(ks), int(len(a)), int(len(b))))
print(f"  usable {len(usable)} · DROPPED {len(dropped)}")
for c, k, na, nb in dropped:
    print(f"    ⚠ dropped {c:10s} cohorts={k} n({int(y0)})={na} n({int(y1)})={nb}")
if TARGET not in [c for c, *_ in usable]:
    print("  ⛔ the TARGET item itself failed the precondition — nothing to compare"); sys.exit(2)


def within_move(e, col, ya, yb, ks=None):
    """within-cohort movement in units of the item's own pooled SD (scale-free)."""
    a, b = e[e.year == ya], e[e.year == yb]
    if ks is None:
        ks = [k for k in sorted(set(a.grp) & set(b.grp))
              if (a.grp == k).sum() >= 40 and (b.grp == k).sum() >= 40]
    if len(ks) < 3:
        return np.nan, 0, 0
    sd = float(pd.concat([a[col], b[col]]).std())
    if not sd or np.isnan(sd) or sd <= 0:
        return np.nan, 0, 0
    moves = [float(b.loc[b.grp == k, col].mean() - a.loc[a.grp == k, col].mean()) / sd for k in ks]
    return float(np.mean(moves)), len(ks), int(len(a) + len(b))


rows = []
for c, _k, _na, _nb in usable:
    e = raw.dropna(subset=[c])
    m, k, n = within_move(e, c, y0, y1)
    if np.isnan(m):
        continue
    cat = "A target" if c == TARGET else ("B same-topic" if c in TOPIC else "C other")
    rows.append(dict(item=c, category=cat, move=m, absmove=abs(m), cohorts=k, n=n))
rows.sort(key=lambda r: -r["absmove"])

print(f"\n=== THE GRID — within-cohort movement, {int(y0)}->{int(y1)}, in item-SD units ===")
for i, r in enumerate(rows):
    star = "  <<<" if r["item"] == TARGET else ("  <<" if r["item"] in TOPIC else "")
    print(f"  {i+1:2d}. {r['item']:10s} [{r['category']:12s}] move {r['move']:+.4f} "
          f"|move| {r['absmove']:.4f}  cohorts={r['cohorts']:2d} n={r['n']:5d}{star}")

others = [r["absmove"] for r in rows if r["category"] == "C other"]
tgt = [r for r in rows if r["item"] == TARGET][0]
topic = [r for r in rows if r["category"] == "B same-topic"]
pct_target = float((np.array(others) < tgt["absmove"]).mean())
pct_topic = [float((np.array(others) < r["absmove"]).mean()) for r in topic]
print(f"\n  `{TARGET}` |move| {tgt['absmove']:.4f} -> percentile {pct_target:.1%} of {len(others)} "
      f"non-homosexuality items")
for r, p in zip(topic, pct_topic):
    print(f"  `{r['item']}` |move| {r['absmove']:.4f} -> percentile {p:.1%}")
print(f"  median |move| among the {len(others)} others: {float(np.median(others)):.4f}")

# ══ NEGATIVE CONTROL — permute YEAR within cohort, per item ══════════════════════════
null_vals = []
e_t = raw.dropna(subset=[TARGET])
for _ in range(200):
    p = e_t[e_t.year.isin([y0, y1])].copy()
    p["year"] = p.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
    m, _, _ = within_move(p, TARGET, y0, y1)
    if not np.isnan(m):
        null_vals.append(abs(m))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (year permuted within cohort; kind of null: within-cohort year-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant INTO the permuted world so g=0 lands on the null ════════
sweep = []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(20):
        p = e_t[e_t.year.isin([y0, y1])].copy()
        p["year"] = p.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
        p.loc[p.year == y1, TARGET] = np.clip(p.loc[p.year == y1, TARGET] + g, 1, 4)
        m, _, _ = within_move(p, TARGET, y0, y1)
        if not np.isnan(m):
            vals.append(abs(m))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (planted into the permuted null world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["absmove"] - null_med) / (null_sd or 1e-9)))) for r in rows]

if not rows:
    print("EMPTY POPULATION"); sys.exit(2)

top_decile = pct_target >= 0.90
topic_high = np.mean(pct_topic) >= 0.90 if pct_topic else False

G = Gate("Did they change on THIS, or on everything?")
G.plant_direction_from_sweep("positive: a planted within-cohort trend raises |move|, and g=0 is null",
                             sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("year permuted within cohort", null_med, tgt["absmove"],
                   null_spread=null_sd, null_kind="within-cohort year-label permutation")
G.multiplicity_control("all items in the comparison family", ps, 0.05,
                       labels=[r["item"] for r in rows])
G.asserted("PRECONDITIONS were checked and printed BEFORE the estimator, and failures were DROPPED "
           "WITH A COUNT (`#925`(2))", True,
           f"{len(usable)} usable · {len(dropped)} dropped: {[c for c, *_ in dropped] or 'none'}; "
           f"absence reported, not passed", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", rows)
G.asserted("the comparison set is large enough for a percentile to mean anything",
           len(others) >= 20, f"{len(others)} non-homosexuality items in the reference set",
           kind="control")
G.asserted("prior art / lineage declared", True,
           "`#905` ran this percentile design at the PERSON level for one-factor structure and put "
           "the sexual battery at the 10th percentile of forty; this is its decade/cohort analogue, "
           "not an independent idea; scope stated", kind="control")
G.asserted("KILL: W_GENERAL requires `homosex` to sit in the body of the other items",
           not top_decile,
           f"`{TARGET}` at the {pct_target:.1%} percentile of {len(others)} others; same-topic "
           f"tolerance at {[f'{p:.0%}' for p in pct_topic]}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif top_decile and topic_high:
    VERDICT, WORLD = "OVERTURNED", "W_TOPIC · attitudes about homosexuality moved, not sexual morality"
elif top_decile and not topic_high:
    VERDICT, WORLD = "OVERTURNED", "W_MORAL · the moral judgement moved uniquely"
else:
    VERDICT, WORLD = "CONFIRMED", "W_GENERAL · everything moved together"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=926, round="E03·A115·R364", verdict=VERDICT, world=WORLD,
           estimand="within-cohort movement per item in item-SD units, and `homosex`'s percentile "
                    "among the non-homosexuality items",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           window=[int(y0), int(y1)], usable=len(usable),
           dropped=[c for c, *_ in dropped], rows=rows,
           target_absmove=tgt["absmove"], target_percentile=pct_target,
           topic_percentiles=dict(zip([r["item"] for r in topic], pct_topic)),
           others_median=float(np.median(others)), n_others=len(others),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           checks="`#925`(1)'s impossibility claim tested in its measurable form",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "one_item_or_thirty_seven.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'one_item_or_thirty_seven.json'}")
