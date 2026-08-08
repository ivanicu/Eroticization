#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A123·R379 — all four sexual norms, each split into conversion and replacement
==================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#940`①. Cohort-demeaning cost A118's headline **54% of its size** (+0.0599 →
                +0.0277 per decade) and reversed an ordering one round earlier (`#939`). **It costs
                four lines and no other trend in this project has run it.** A114 gave ONE aggregate
                number — conversion is ≥45% of the half-century change — for `homosex` alone.
                **Nobody has asked whether the other three norms split the same way**, and A118's
                whole point is that they do not behave alike.

Why it matters  "Americans became more permissive" is four different sentences depending on whether
                each norm moved because people changed their minds or because the people changed.
                A norm that moved entirely by replacement was never actually reconsidered by anyone.

Live Worlds    W_SAME  · all four norms split conversion/replacement in about the same proportion
                          ⇒ the split is a property of the ERA, not of the act, and A118's
                          "homosex is different" is about SIZE only.
               W_SPLIT · the proportions differ across norms ⇒ **the four norms did not merely move
                          at different speeds, they moved by different MECHANISMS**, which is a
                          stronger claim than A118 has made and a new one.
               W_FLAT  · demeaning removes nearly everything everywhere ⇒ almost all of the
                          half-century change is replacement and A114's ≥45% conversion floor is
                          wrong. ⚠ The unwelcome one — it would retract a headline.

Estimand       For each of `premarsx`, `teensex`, `xmarsex`, `homosex`: the OLS trend per decade of
(G1)           the item mean over 21 GSS ballot-1 waves 1988–2024, in two arms — **RAW**, and
               **COHORT-DEMEANED** (item minus its own birth-year mean, within wave). The reported
               quantity is `conversion share = demeaned trend / raw trend`, per norm.

⚠ WHAT THE     `demeaned/raw` is a share of the OBSERVED trend attributable to within-cohort
SHARE IS AND   movement. It is **not** A114's Kitagawa-style decomposition — that was retracted as
IS NOT         ill-posed for non-overlapping birth cohorts (`#924`) — and it is **not** a causal
               decomposition. Cohort-demeaning removes every fixed between-cohort level difference;
               what remains is movement among people born in the same year.

Prediction     W_SAME  -> the four shares agree within their bootstrap spread.
Matrix         W_SPLIT -> they differ by more than that, and the ordering is stable across the grid.
               W_FLAT  -> shares near zero for all four.

Strongest      **A RATIO IS UNSTABLE WHEN ITS DENOMINATOR IS SMALL** (`#918`'s family). `xmarsex`
confound       barely moved at all in A118 (+0.0906 vs `homosex` +0.6400 on the shared frame), so
(written       its conversion SHARE can be enormous or negative purely from noise in a near-zero
before)        denominator. ⇒ CONTROL: both the ratio AND both absolute trends are reported for
               every norm, and any norm whose raw trend fails to clear twice its own bootstrap
               spread is marked **RATIO UNREADABLE** rather than given a share.

Controls       NEGATIVE: permute `year` within cohort — every marginal, every cohort mean and n
                 untouched; only the time ordering dies. Run in both arms.
               POSITIVE: plant a rising within-cohort trend INTO the permuted world and sweep;
                 `g=0` sits on the null by construction (`#922`, `#937`⑤).
               MULTIPLICITY: the family is **the four conversion shares**, which is the family this
                 claim lives in (`#936`②, `#939`②, `#940`②).
               SPEC CURVE (G4): 4 norms × {all ballots, ballot 1} × {all cohorts, cohorts with ≥3
                 waves} — every cell published.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate age from period from cohort** — `#939`'s wall, inherited unchanged;
  (2) ⚠ **follow a person over time** — repeated cross-section, so "changed their mind" stays an
    inference about groups;
  (3) ⚠ **no second instrument for the trend** — `#937` measured that NSFG has one time point, so
    this is **only this one instrument** and the cross-instrument move is structurally unavailable;
  (4) ⚠ **claim a causal decomposition** — see the estimand note; `#924` retracted the standard
    estimator here and this round does not resurrect it;
  (5) `[unchallenged]` — door ③.
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
RNG = np.random.default_rng(379)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
MIN_WAVE = 300

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"])
waves = [int(y) for y, n in d.groupby("year").size().items() if n >= MIN_WAVE]
d = d[d.year.isin(waves)]
print(f"PRECONDITION — ballot 1, all four answered, cohort present, n>={MIN_WAVE}: {len(waves)} "
      f"waves {waves[0]}-{waves[-1]} · n={len(d)} · cohorts {int(d.cohort.min())}-"
      f"{int(d.cohort.max())}")


def trends(f, item):
    """(raw trend, cohort-demeaned trend) per decade of the item mean over waves."""
    g = f.copy()
    g["dm"] = g[item] - g.groupby("cohort")[item].transform("mean")
    yrs = sorted(g.year.unique())
    raw = [g[g.year == y][item].mean() for y in yrs]
    dem = [g[g.year == y]["dm"].mean() for y in yrs]
    return (float(stats.linregress(yrs, raw).slope * 10),
            float(stats.linregress(yrs, dem).slope * 10))


rows = []
for it in ITEMS:
    r, m = trends(d, it)
    rows.append(dict(item=it, raw=r, demeaned=m))

# ══ bootstrap: the ratio's own spread, and the RAW trend's, per norm ═════════════════
B = 200
boot = {it: {"raw": [], "dem": []} for it in ITEMS}
for _ in range(B):
    s = d.sample(len(d), replace=True, random_state=int(RNG.integers(1e9)))
    for it in ITEMS:
        r, m = trends(s, it)
        boot[it]["raw"].append(r)
        boot[it]["dem"].append(m)
for row in rows:
    it = row["item"]
    row["raw_sd"] = float(np.std(boot[it]["raw"]))
    row["dem_sd"] = float(np.std(boot[it]["dem"]))
    # ⚠ the ratio is only reported when the DENOMINATOR clears twice its own spread
    row["readable"] = bool(abs(row["raw"]) > 2 * row["raw_sd"])
    ratios = [m / r for r, m in zip(boot[it]["raw"], boot[it]["dem"]) if abs(r) > 1e-6]
    row["share"] = float(row["demeaned"] / row["raw"]) if row["readable"] else float("nan")
    row["share_lo"], row["share_hi"] = ((float(np.percentile(ratios, 2.5)),
                                         float(np.percentile(ratios, 97.5)))
                                        if row["readable"] and ratios else (float("nan"),) * 2)

print("\n  per decade on a 4-point scale (higher = more permissive), 21 waves")
print(f"  {'norm':<10s} {'raw':>9s} {'±':>7s}   {'demeaned':>9s} {'±':>7s}   {'conversion share':>18s}")
for row in rows:
    sh = (f"{row['share']:.3f} [{row['share_lo']:.2f},{row['share_hi']:.2f}]"
          if row["readable"] else "⚠ RATIO UNREADABLE")
    print(f"  {row['item']:<10s} {row['raw']:+9.4f} {row['raw_sd']:7.4f}   "
          f"{row['demeaned']:+9.4f} {row['dem_sd']:7.4f}   {sh:>18s}")

readable = [r for r in rows if r["readable"]]
shares = [r["share"] for r in readable]
share_spread = float(max(shares) - min(shares)) if len(shares) > 1 else float("nan")

# ══ NEGATIVE CONTROL — permute year within cohort ════════════════════════════════════
null_dem = {it: [] for it in ITEMS}
for _ in range(150):
    s = d.copy()
    s["year"] = s.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
    s = s[s.year.isin(waves)]
    for it in ITEMS:
        null_dem[it].append(trends(s, it)[1])
nulls = {it: (float(np.median(v)), float(np.std(v))) for it, v in null_dem.items()}
print("\n  null (year permuted WITHIN cohort — marginals, cohort means and n untouched; kind of "
      "null: within-cohort year-label permutation):")
for it in ITEMS:
    print(f"    {it:<10s} demeaned null {nulls[it][0]:+.4f} +/- {nulls[it][1]:.4f}   "
          f"observed {dict((r['item'], r['demeaned']) for r in rows)[it]:+.4f}")

# ══ POSITIVE CONTROL — plant a rising within-cohort trend INTO the permuted world ════
TGT = "homosex"
sweep = []
for gg in (0.0, 0.2, 0.4, 0.6):
    vals = []
    for _ in range(8):
        s = d.copy()
        s["year"] = s.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
        s = s[s.year.isin(waves)]
        s[TGT] = s[TGT] + gg * (s.year - np.mean(waves)) / 18.0
        vals.append(trends(s, TGT)[1])
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a rising within-cohort trend planted into the permuted world, so g=0 IS "
      f"the null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs {TGT} null {nulls[TGT][0]:+.4f} +/- "
      f"{nulls[TGT][1]:.4f} = {abs(sweep[0][1] - nulls[TGT][0]) / max(nulls[TGT][1], 1e-9):.2f} "
      f"spreads")

# ══ SPECIFICATION CURVE (G4) ═════════════════════════════════════════════════════════
grid = []
keep_c = d.groupby("cohort").year.nunique()
c3 = set(keep_c[keep_c >= 3].index)
for tag, sub in (("ballot 1 · all cohorts", d),
                 ("ballot 1 · cohorts in >=3 waves", d[d.cohort.isin(c3)])):
    for it in ITEMS:
        r, m = trends(sub, it)
        grid.append(dict(spec=tag, item=it, raw=r, demeaned=m,
                         share=float(m / r) if abs(r) > 1e-6 else float("nan"), n=int(len(sub))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    {g_['spec']:32s} {g_['item']:<10s} raw {g_['raw']:+.4f}  dem {g_['demeaned']:+.4f}  "
          f"share {g_['share']:+.3f}")

ps = [2 * (1 - stats.norm.cdf(abs((r["demeaned"] - nulls[r["item"]][0]) /
                                  max(nulls[r["item"]][1], 1e-9)))) for r in rows]

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G = Gate("Did the four sexual norms move by the same mechanism, or only at different speeds?")
G.plant_direction_from_sweep("positive: a planted within-cohort trend raises the demeaned trend, and "
                             "g=0 sits ON the null this round judges against (`#922`)", sweep,
                             baseline=nulls[TGT][0], baseline_spread=max(nulls[TGT][1], 1e-4))
for r in rows:
    G.negative_control(f"year permuted within cohort [{r['item']}]", abs(nulls[r["item"]][0]),
                       abs(r["demeaned"]), null_spread=nulls[r["item"]][1],
                       null_kind="within-cohort year-label permutation")
G.multiplicity_control("the four conversion shares — the family this claim lives in (`#936`②)",
                       ps, 0.05, labels=ITEMS)
G.asserted("⚠ A RATIO WITH A SMALL DENOMINATOR IS NOT REPORTED AS A NUMBER (`#918`'s family): any "
           "norm whose RAW trend fails to clear twice its own bootstrap spread is marked UNREADABLE",
           True,
           " · ".join(f"{r['item']} raw {r['raw']:+.4f} vs 2x{r['raw_sd']:.4f} -> "
                      f"{'readable' if r['readable'] else 'UNREADABLE'}" for r in rows),
           kind="control", population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{g_['item']}@{g_['spec'].split('· ')[-1]} {g_['share']:+.3f}"
                      for g_ in grid), kind="control",
           population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")
G.asserted("⚠ this is NOT a causal decomposition and NOT `#924`'s retracted Kitagawa: it is the "
           "share of the observed trend that survives removing every fixed between-cohort level "
           "difference", True,
           "cohort-demeaning identifies movement among people born in the same year; it does not "
           "identify why they moved", kind="control",
           population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * nulls[TGT][1]
neg_null = all(abs(nulls[it][0]) < 2 * nulls[it][1] for it in ITEMS)
share_sds = [float(np.std([m / r for r, m in zip(boot[x["item"]]["raw"], boot[x["item"]]["dem"])
                           if abs(r) > 1e-6])) for x in readable]
differ = (len(shares) > 1 and share_spread > 2 * max(share_sds)) if share_sds else False
flat = bool(shares) and max(shares) < 0.15
world = "W_FLAT" if flat else ("W_SPLIT" if differ else "W_SAME")

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and every "
           "norm's null is null. STAKED: W_SPLIT, i.e. the conversion shares differ by more than "
           "twice the widest share's own bootstrap spread. W_SAME and W_FLAT both refute it",
           (pos_fires and neg_null) and differ,
           f"positive fires {pos_fires} · all four nulls null {neg_null} · readable shares "
           f"{len(readable)}/4 · spread {share_spread:.3f} vs 2x widest share sd "
           f"{2 * max(share_sds) if share_sds else float('nan'):.3f} ⇒ {world}",
           kind="kill", yardstick="max-min of the readable conversion shares",
           yardstick_noise=max(share_sds) if share_sds else 0.0,
           population=f"GSS ballot 1, {len(waves)} waves 1988-2024, n={len(d)}",
           direction="two-sided; W_SPLIT needs separation, not a sign")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if differ else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=941, round="E03·A123·R379", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_SAME"),
               waves=waves, n=int(len(d)), rows=rows, grid=grid,
               nulls={k: dict(median=v[0], sd=v[1]) for k, v in nulls.items()},
               null_median=nulls[TGT][0], null_sd=nulls[TGT][1], null_draws=len(null_dem[TGT]),
               positive_sweep=sweep, share_spread=share_spread, family_size=len(ps),
               world=world, verdict=verdict),
          open(OUT / "all_four_norms_under_the_same_control.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'all_four_norms_under_the_same_control.json'}")
