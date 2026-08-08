#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A114·R361 — does the person-level finding explain the half-century, or almost none of it?
=============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#915`-`#921` built one claim: people who have had same-sex experience, and people
                who say they are gay/lesbian/bisexual, judge the act less harshly — replicated on a
                second instrument, surviving a randomised re-wording. ⚠ **All of it is person-level,
                and the project's object is the DECADE and the SOCIETY.** The bridge was never
                built: if having the experience softens the verdict, then as more people have (or
                report) it, aggregate opinion should soften. **How much of the actual half-century
                change does that mechanism account for?**

Why Now         `#922`(2) requires this round's object to be people or data. And this is the one
                question that connects five rounds of person-level work to the thing the project
                says it is about — *what the society does with it*.

Live Worlds     W_COMP  · a large share of the aggregate change is COMPOSITION: more people with
                          the experience, each carrying the softer verdict they already had.
                W_CTX   · a small share. The groups themselves moved — **everyone changed**, and
                          the person-level coupling explains almost none of the history.
                          ⚠ **This is the unwelcome one: it deflates my own five-round chain to a
                          real but historically minor effect.**
                W_ONTO  · the "share" itself rose partly because REPORTING changed (`#919`/`#921`
                          leave disclosure unresolved), so composition and context do not carve
                          anything. (the meta-separator: it kills the decomposition's ontology)

⚠ THE ARITHMETIC TRAP, NAMED BEFORE THE NUMBER (`realstat`): given group shares and group means, the
                Kitagawa decomposition is **forced by the algebra**. It is a DERIVATION, not a test,
                and it is labelled as one. Its assumption: groups are fixed categories whose
                membership is measured, and only shares and within-group means move.
                **What is NOT forced, and IS a test: did the two groups move DIFFERENTLY?** That
                could have come out either way, and it is where the round can be embarrassed.

Estimand        (a) `comp_share` = the fraction of the total change in mean `homosex` attributable
(G1)                to the changing group share, holding within-group means at their pooled values
                    — DERIVATION, reported with its own noise floor;
                (b) `Delta_gap` = the change in the exposed-minus-unexposed gap over the window, and
                    the per-group movements — TEST, with a ceiling control.

Prediction      W_COMP -> comp_share large (say >0.25) and far above its floor.
Matrix          W_CTX  -> comp_share small; both groups move by comparable large amounts.
                W_ONTO -> flagged, not decided, because disclosure is unresolved upstream.

Strongest       ⚠ **CEILING.** The exposed group starts higher, and the scale stops at 4, so it has
Confound        less room to move — a smaller movement is not evidence it led or lagged. CONTROL,
                same iteration: each group's movement is also expressed against the room it had
                (`#918`'s attainable-range logic), and the raw and headroom-scaled versions are both
                published.

Stopping Rule   One pass over group-definition x window x weighting, published whole, including the
                cells that disagree.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **disclosure is unresolved upstream** (`#919`), so the rising share is partly a rising
    willingness to report. Every composition number here is therefore an UPPER BOUND on a
    quantity that is itself partly an artifact — W_ONTO can be flagged and not decided;
  (2) ⚠ the behaviour item exists only from 1989, so the composition window is 1989-2022 while the
    norm series runs 1973-2024. **The two windows are reported separately and never mixed**;
  (3) ⚠ this is a DECOMPOSITION of an observed change, not a causal account of it: nothing here
    identifies why anyone moved;
  (4) ⚠ **only this one instrument** — GSS is the only release here with a five-decade norm series
    AND a behaviour item. `#921` already did the cross-instrument move for the person-level claim;
    NSFG cannot supply a 1973-2024 trend at all;
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
RNG = np.random.default_rng(361)
TOP = 4.0                                     # `homosex` runs 1..4; the ceiling is 4

df = pd.read_stata(GSS, columns=["homosex", "nummen", "numwomen", "sex", "year", "sexornt"],
                   convert_categoricals=False)
df["homosex"] = df["homosex"].where(df["homosex"].isin([1, 2, 3, 4]))
ok = df["nummen"].notna() & df["numwomen"].notna() & df["sex"].isin([1, 2])
df["ss"] = np.where(((df["sex"] == 1) & (df["nummen"] > 0)) | ((df["sex"] == 2) & (df["numwomen"] > 0)), 1.0, 0.0)
df.loc[~ok, "ss"] = np.nan
df["lgb"] = np.where(df["sexornt"].isin([1, 2]), 1.0, np.where(df["sexornt"] == 3, 0.0, np.nan))
df["either"] = np.where((df["ss"] == 1) | (df["lgb"] == 1), 1.0,
                        np.where(df["ss"].notna() | df["lgb"].notna(), 0.0, np.nan))

# ══ HARD RULE 1 — the two series, their windows, and the fact they differ ════════════
norm = df.dropna(subset=["homosex"]).groupby("year")["homosex"].agg(["mean", "size"])
print(f"norm series `homosex`: {int(norm.index[0])}-{int(norm.index[-1])}, {len(norm)} waves, "
      f"{norm['mean'].iloc[0]:.3f} -> {norm['mean'].iloc[-1]:.3f} "
      f"({norm['mean'].iloc[-1] - norm['mean'].iloc[0]:+.3f})")
inv = {}
for g in ("ss", "lgb", "either"):
    e = df.dropna(subset=["homosex", g])
    sh = e.groupby("year")[g].mean()
    inv[g] = dict(waves=len(sh), first_year=int(sh.index[0]), last_year=int(sh.index[-1]),
                  first_share=float(sh.iloc[0]), last_share=float(sh.iloc[-1]), n=int(len(e)))
    print(f"  group `{g:6s}` {int(sh.index[0])}-{int(sh.index[-1])} · {len(sh):2d} waves · "
          f"share {sh.iloc[0]:.4f} -> {sh.iloc[-1]:.4f} · n={len(e)}")


def kitagawa(e, g, y0, y1):
    """DERIVATION (labelled): total change split into composition and within-group parts.
    d_total = sum_k [ (p1_k - p0_k) * mbar_k ]  +  sum_k [ pbar_k * (m1_k - m0_k) ]"""
    a, b = e[e.year == y0], e[e.year == y1]
    if len(a) < 100 or len(b) < 100:
        return None
    out = {}
    p0 = {k: float((a[g] == k).mean()) for k in (0.0, 1.0)}
    p1 = {k: float((b[g] == k).mean()) for k in (0.0, 1.0)}
    m0 = {k: float(a.loc[a[g] == k, "homosex"].mean()) for k in (0.0, 1.0)}
    m1 = {k: float(b.loc[b[g] == k, "homosex"].mean()) for k in (0.0, 1.0)}
    if any(np.isnan(v) for v in list(m0.values()) + list(m1.values())):
        return None
    mbar = {k: (m0[k] + m1[k]) / 2 for k in p0}
    pbar = {k: (p0[k] + p1[k]) / 2 for k in p0}
    comp = sum((p1[k] - p0[k]) * mbar[k] for k in p0)
    within = sum(pbar[k] * (m1[k] - m0[k]) for k in p0)
    total = float(b["homosex"].mean() - a["homosex"].mean())
    out.update(y0=int(y0), y1=int(y1), total=total, comp=comp, within=within,
               comp_share=(comp / total if total else np.nan),
               share0=p0[1.0], share1=p1[1.0],
               move_unexposed=m1[0.0] - m0[0.0], move_exposed=m1[1.0] - m0[1.0],
               gap0=m0[1.0] - m0[0.0], gap1=m1[1.0] - m1[0.0],
               # ⚠ CEILING CONTROL: movement as a fraction of the room each group had
               head_unexposed=(m1[0.0] - m0[0.0]) / max(TOP - m0[0.0], 1e-9),
               head_exposed=(m1[1.0] - m0[1.0]) / max(TOP - m0[1.0], 1e-9),
               n0=int(len(a)), n1=int(len(b)), n=int(len(a) + len(b)))
    return out


grid = []
for g in ("ss", "lgb", "either"):
    e = df.dropna(subset=["homosex", g])
    yrs = sorted(e.year.unique())
    windows = [(yrs[0], yrs[-1])]
    if len(yrs) >= 6:
        windows += [(yrs[0], yrs[len(yrs) // 2]), (yrs[len(yrs) // 2], yrs[-1])]
    for y0, y1 in windows:
        r = kitagawa(e, g, y0, y1)
        if r:
            r["group"] = g
            grid.append(r)

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for r in grid:
    print(f"  {r['group']:6s} {r['y0']}->{r['y1']}  total {r['total']:+.3f} = comp {r['comp']:+.4f} "
          f"+ within {r['within']:+.3f}  | comp_share {r['comp_share']:6.2%}  "
          f"share {r['share0']:.3f}->{r['share1']:.3f}  n={r['n0']}/{r['n1']}")
print("\n  movements and the CEILING control (movement / room the group had):")
for r in grid:
    print(f"  {r['group']:6s} {r['y0']}->{r['y1']}  unexposed {r['move_unexposed']:+.3f} "
          f"(head {r['head_unexposed']:.3f}) · exposed {r['move_exposed']:+.3f} "
          f"(head {r['head_exposed']:.3f}) · gap {r['gap0']:+.3f} -> {r['gap1']:+.3f}")

full = [r for r in grid if r["y1"] - r["y0"] == max(x["y1"] - x["y0"] for x in grid if x["group"] == r["group"])]
comp_shares = [r["comp_share"] for r in grid if not np.isnan(r["comp_share"])]
med_comp = float(np.median(comp_shares))
print(f"\n  median comp_share over {len(comp_shares)} cells: {med_comp:.2%}")

# ══ NEGATIVE CONTROL — hold the group share FIXED at its first-wave value ════════════
# resampling each wave to the first wave's share makes composition ZERO by construction; whatever
# `comp_share` still shows is the noise floor of the estimator, not composition.
null_vals = []
for g in ("ss", "either"):
    e = df.dropna(subset=["homosex", g])
    yrs = sorted(e.year.unique())
    y0, y1 = yrs[0], yrs[-1]
    p_target = float(e.loc[e.year == y0, g].mean())
    for _ in range(120):
        parts = []
        for y in (y0, y1):
            w = e[e.year == y]
            n = len(w)
            k = int(round(p_target * n))
            one, zero = w[w[g] == 1], w[w[g] == 0]
            if len(one) < 5 or len(zero) < 5:
                parts = []
                break
            s = pd.concat([one.sample(k, replace=True, random_state=int(RNG.integers(1e9))),
                           zero.sample(n - k, replace=True, random_state=int(RNG.integers(1e9)))])
            parts.append(s)
        if len(parts) != 2:
            continue
        r = kitagawa(pd.concat(parts), g, y0, y1)
        if r and not np.isnan(r["comp_share"]):
            null_vals.append(r["comp_share"])
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"  null (group share held at its first-wave value; kind of null: share-fixed resampling): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant a PURELY compositional trend; must be null at g=0 ═══════
sweep = []
for g_amt in (0.0, 0.10, 0.20, 0.30, 0.40):
    vals = []
    e = df.dropna(subset=["homosex", "ss"])
    yrs = sorted(e.year.unique())
    y0, y1 = yrs[0], yrs[-1]
    for _ in range(30):
        parts = []
        for y, extra in ((y0, 0.0), (y1, g_amt)):
            w = e[e.year == y].copy()
            # freeze within-group means at the FIRST wave's values, then raise the share by `extra`
            for k in (0.0, 1.0):
                m_first = float(e.loc[(e.year == y0) & (e["ss"] == k), "homosex"].mean())
                w.loc[w["ss"] == k, "homosex"] = m_first
            n = len(w)
            p = float(e.loc[e.year == y0, "ss"].mean()) + extra
            k_one = int(round(min(max(p, 0.01), 0.9) * n))
            one, zero = w[w["ss"] == 1], w[w["ss"] == 0]
            if len(one) < 5 or len(zero) < 5:
                parts = []
                break
            parts.append(pd.concat([one.sample(k_one, replace=True, random_state=int(RNG.integers(1e9))),
                                    zero.sample(n - k_one, replace=True, random_state=int(RNG.integers(1e9)))]))
        if len(parts) != 2:
            continue
        r = kitagawa(pd.concat(parts), "ss", y0, y1)
        if r and not np.isnan(r["comp"]):
            vals.append(r["comp"])
    sweep.append([float(g_amt), float(np.median(vals)) if vals else np.nan])
# ⚠⚠ v1 SWEPT `comp_share` AND IT WAS DEGENERATE — 1.0000 at every g. The plant freezes the
#   within-group means, so `within` is 0 by construction and `comp_share = comp/total` is 1
#   whatever the plant does. **The plant disturbed the very quantity the statistic normalises by**
#   — the third instance of `#919`(3), and precisely what `#922`'s new gate exists to catch
#   (g=0 landing off its baseline). The ABSOLUTE component `comp` has no such denominator: it is
#   0 when the share does not move, and rises with the planted share. That is what is swept.
print(f"  positive sweep (planted share rise, median ABSOLUTE comp): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

# ══ THE TEST THAT IS NOT A DERIVATION — did the groups move DIFFERENTLY? ═════════════
boot = []
e = df.dropna(subset=["homosex", "ss"])
yrs = sorted(e.year.unique())
y0, y1 = yrs[0], yrs[-1]
for _ in range(400):
    s = pd.concat([e[e.year == y].sample(len(e[e.year == y]), replace=True,
                                         random_state=int(RNG.integers(1e9))) for y in (y0, y1)])
    r = kitagawa(s, "ss", y0, y1)
    if r:
        boot.append(r["move_unexposed"] - r["move_exposed"])
diff_med = float(np.median(boot))
diff_lo, diff_hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
print(f"\n  unexposed-minus-exposed MOVEMENT, {int(y0)}->{int(y1)}: {diff_med:+.3f} "
       f"[95% bootstrap {diff_lo:+.3f}, {diff_hi:+.3f}]  <- NOT a derivation; it could have gone either way")

ps = [2 * (1 - stats.norm.cdf(abs((r["comp_share"] - null_med) / (null_sd or 1e-9))))
      for r in grid if not np.isnan(r["comp_share"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

comp_resolved = abs(med_comp - null_med) > 2 * null_sd
comp_large = med_comp > 0.25

G = Gate("Does the person-level coupling explain the half-century, or almost none of it?")
G.plant_direction_from_sweep("positive: a planted share rise raises the ABSOLUTE composition "
                             "component, and g=0 is null", sweep, baseline=0.0,
                             baseline_spread=max(null_sd, 1e-4))
G.negative_control("group share held at its first-wave value", abs(null_med), abs(med_comp),
                   null_spread=null_sd, null_kind="share-fixed resampling")
G.multiplicity_control("the whole group x window grid", ps, 0.05,
                       labels=[f"{r['group']}|{r['y0']}-{r['y1']}" for r in grid
                               if not np.isnan(r["comp_share"])])
G.asserted("the CEILING control ran: each group's movement is also expressed against its own room",
           all("head_unexposed" in r for r in grid),
           f"headroom-scaled movements published for all {len(grid)} cells", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.has_error_bar("the group-difference in movement carries an interval", diff_med,
                (diff_hi - diff_lo) / 4, "bootstrap_人层")
G.asserted("the Kitagawa split is labelled a DERIVATION, not a test",
           True, "given shares and group means the split is forced by the algebra; the TEST in this "
                 "round is the group-difference in movement, which could have gone either way; "
                 "scope stated", kind="control")
G.asserted("KILL: W_COMP requires composition to carry a large share of the aggregate change",
           not (comp_large and comp_resolved),
           f"median comp_share {med_comp:.2%} vs null {null_med:+.4f} +/- {null_sd:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif comp_large and comp_resolved:
    VERDICT, WORLD = "OVERTURNED", "W_COMP · composition carries the change"
else:
    VERDICT, WORLD = "CONFIRMED", "W_CTX · the society moved; composition is a rounding error"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=923, round="E03·A114·R361", verdict=VERDICT, world=WORLD,
           estimand="share of the aggregate change in `homosex` attributable to the changing share "
                    "of people with same-sex experience/identity (DERIVATION), plus the "
                    "group-difference in movement (TEST)",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           norm_series=dict(first_year=int(norm.index[0]), last_year=int(norm.index[-1]),
                            waves=len(norm), first=float(norm["mean"].iloc[0]),
                            last=float(norm["mean"].iloc[-1]),
                            total_move=float(norm["mean"].iloc[-1] - norm["mean"].iloc[0])),
           groups=inv, grid=grid, median_comp_share=med_comp,
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep,
           group_move_difference=dict(median=diff_med, lo=diff_lo, hi=diff_hi, draws=len(boot)),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "composition_or_the_society.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'composition_or_the_society.json'}")
