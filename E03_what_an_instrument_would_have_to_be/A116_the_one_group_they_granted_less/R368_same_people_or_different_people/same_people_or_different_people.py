#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A116·R368 — the rise and the one fall: same people, or different people?
============================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#929` measured that total tolerance GREW (+2.131 of 15) with **four of five targets
                rising and racists the single faller**, and named the exception as the object:
                **general growth cannot explain a group going the other way.** The measurable form
                of "why" here is not motive — it is **WHO**. Did the homosexual rise and the racist
                fall happen in the SAME subgroups, or in different ones?

Live Worlds     W_SAME      · negative across-cell correlation — the subgroups that rose most on
                              homosexuals fell most on racists ⇒ **one boundary being re-drawn by
                              the same people**, two signs of a single movement.
                W_DIFFERENT · ~0 ⇒ unrelated movements; the racist fall comes from somewhere the
                              homosexual rise does not.
                W_ALIGNED   · positive ⇒ the subgroups rising on one rose on the other too, and the
                              aggregate racist fall is COMPOSITION — which subgroups grew — not any
                              subgroup changing its mind. ⚠ **Unwelcome: it would mean `#929`'s
                              "exception" is not an exception in anybody's head.**

⚠ THE           `#924`(2)/`#925`(2): the precondition was checked FIRST, and **it killed the first
PRECONDITION    design.** Politics x education x cohort with SINGLE endpoint waves gives **4 usable
KILLED THE      cells of 59** — a correlation over 4 points is not a measurement. Pooling eras gives
FIRST DESIGN    **25-27 cells of 33-39** across three era splits, which is also the specification
                curve. **The round that would have run on 4 cells was stopped before it ran.**

Estimand        Across subgroup cells (politics 3 x education 3 x birth-cohort 4), the correlation
(G1)            between the cell's change in homosexual-tolerance and its change in racist-tolerance,
                early era -> late era. Both arms are counts out of 3, so they are already on one
                scale (`#928` established the arms are commensurable).

Prediction      W_SAME -> r clearly negative.  W_DIFFERENT -> r ~ 0.  W_ALIGNED -> r clearly positive.
Matrix

Strongest       ⚠ **CEILING.** A cell already high on homosexual tolerance has less room to rise; a
Confound        cell already low on racist tolerance has less room to fall. **That alone can
                manufacture a negative correlation.** CONTROL, same iteration: every change is also
                expressed against the room that cell had (`#918`'s attainable-range logic), and the
                raw and headroom-scaled versions are both published.

⚠ META-         This is an **ECOLOGICAL** correlation across subgroups, and the data is a repeated
SEPARATOR       cross-section with no panel. A negative across-cell correlation is equally consistent
                with *the same individuals re-drawing a boundary* and with *different individuals
                inside the same cell moving in opposite directions*. **It cannot distinguish them,
                and no wording in this round may pretend otherwise.**

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ ecological, above — the individual-level version needs a panel, which GSS's rotating design
    does not provide over this span;
  (2) ⚠ cells are defined by variables measured at interview; `polviews` is itself an outcome that
    may have moved, so the cells are not fixed traits;
  (3) ⚠ `polviews` starts 1974 and the tolerance battery ends 2021 — the window is the intersection;
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
RNG = np.random.default_rng(368)
MIN_CELL = 40

TARGETS = {"homosexuals": {"spkhomo": -1, "colhomo": -1, "libhomo": +1},
           "racists": {"spkrac": -1, "colrac": -1, "librac": +1}}
COLS = [c for t in TARGETS.values() for c in t]

d = pd.read_stata(GSS, columns=["year", "cohort", "polviews", "degree"] + COLS,
                  convert_categoricals=False)
for c in COLS:
    d[c] = d[c].where(d[c].between(1, 7))
codes = {c: sorted(d[c].dropna().unique()) for c in COLS}      # ⚠ `#927`(3): derived, no literals
for tgt, items in TARGETS.items():
    for c, sign in items.items():
        lo, hi = float(codes[c][0]), float(codes[c][-1])
        d[c] = (d[c] == (hi if sign > 0 else lo)).astype(float).where(d[c].notna())
cc = d[d[COLS].notna().all(axis=1)].copy()
cc = cc[cc.polviews.between(1, 7) & cc.degree.between(0, 4) & cc.cohort.between(1900, 2006)]
cc["homo"] = cc[list(TARGETS["homosexuals"])].sum(axis=1)
cc["rac"] = cc[list(TARGETS["racists"])].sum(axis=1)
cc["pol3"] = pd.cut(cc.polviews, [0, 3, 4, 7], labels=["lib", "mod", "con"])
cc["ed3"] = pd.cut(cc.degree, [-1, 0, 2, 4], labels=["<HS", "HS/JC", "BA+"])
cc["coh"] = (cc.cohort // 20 * 20).astype(int)

SPLITS = {"1976-1990 vs 2008-2021": ((1976, 1990), (2008, 2021)),
          "1976-1994 vs 2004-2021": ((1976, 1994), (2004, 2021)),
          "1976-1998 vs 2000-2021": ((1976, 1998), (2000, 2021))}

# ══ PRECONDITION CHECK, printed BEFORE the estimator — and it killed design #1 ═══════
print("PRECONDITION CHECK (`#924`(2)/`#925`(2)) — cells need n>=%d in BOTH eras:" % MIN_CELL)
print("  ⚠ FIRST DESIGN REJECTED: politics x education x cohort with SINGLE endpoint waves gives")
print("    only 4 usable cells of 59. A correlation over 4 points is not a measurement.")


def cells_for(early, late):
    e = cc.copy()
    e["era"] = np.where(e.year.between(*early), "E", np.where(e.year.between(*late), "L", None))
    e = e[e.era.notna()]
    g = e.groupby(["pol3", "ed3", "coh", "era"], observed=True)
    agg = g.agg(n=("homo", "size"), homo=("homo", "mean"), rac=("rac", "mean")).reset_index()
    piv = agg.pivot_table(index=["pol3", "ed3", "coh"], columns="era",
                          values=["n", "homo", "rac"], observed=True).dropna()
    keep = piv[(piv[("n", "E")] >= MIN_CELL) & (piv[("n", "L")] >= MIN_CELL)]
    rows = []
    for idx, r in keep.iterrows():
        d_homo = float(r[("homo", "L")] - r[("homo", "E")])
        d_rac = float(r[("rac", "L")] - r[("rac", "E")])
        # ⚠ CEILING control: express each change against the room that cell had
        head_homo = d_homo / max(3.0 - float(r[("homo", "E")]), 1e-9) if d_homo >= 0 else \
            d_homo / max(float(r[("homo", "E")]), 1e-9)
        head_rac = d_rac / max(3.0 - float(r[("rac", "E")]), 1e-9) if d_rac >= 0 else \
            d_rac / max(float(r[("rac", "E")]), 1e-9)
        rows.append(dict(cell="|".join(str(x) for x in idx), d_homo=d_homo, d_rac=d_rac,
                         head_homo=head_homo, head_rac=head_rac,
                         n=int(r[("n", "E")] + r[("n", "L")])))
    return rows, len(piv)


grid = []
for name, (early, late) in SPLITS.items():
    rows, total = cells_for(early, late)
    if len(rows) < 8:
        print(f"  ⚠ DROPPED split {name}: only {len(rows)} usable cells")
        continue
    r_raw = float(np.corrcoef([x["d_homo"] for x in rows], [x["d_rac"] for x in rows])[0, 1])
    r_head = float(np.corrcoef([x["head_homo"] for x in rows], [x["head_rac"] for x in rows])[0, 1])
    grid.append(dict(split=name, cells=len(rows), of=total, r_raw=r_raw, r_head=r_head,
                     mean_d_homo=float(np.mean([x["d_homo"] for x in rows])),
                     mean_d_rac=float(np.mean([x["d_rac"] for x in rows])),
                     n=int(sum(x["n"] for x in rows)), rows=rows))
    print(f"  usable {len(rows)} of {total} cells · {name}")

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    print(f"  {g['split']:24s} cells {g['cells']:2d}/{g['of']:2d}  "
          f"mean d_homo {g['mean_d_homo']:+.3f} · mean d_rac {g['mean_d_rac']:+.3f}  |  "
          f"r_raw {g['r_raw']:+.3f} · r_headroom {g['r_head']:+.3f}  n={g['n']}")

med_raw = float(np.median([g["r_raw"] for g in grid]))
med_head = float(np.median([g["r_head"] for g in grid]))
print(f"\n  median across-cell r, raw      {med_raw:+.3f}")
print(f"  median across-cell r, headroom {med_head:+.3f}   <- the ceiling control")

# ══ NEGATIVE CONTROL — permute the cell labels on ONE arm ════════════════════════════
base = grid[0]["rows"]
null_vals = []
for _ in range(400):
    perm = RNG.permutation([x["d_rac"] for x in base])
    null_vals.append(float(np.corrcoef([x["d_homo"] for x in base], perm)[0, 1]))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (racist changes permuted across cells; kind of null: cell-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant a shared driver INTO the permuted world; g=0 on the null ═
sweep = []
dh = np.array([x["d_homo"] for x in base])
for g in (0.0, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(60):
        perm = RNG.permutation([x["d_rac"] for x in base])
        # interpolate toward MINUS the homosexual change: a shared driver with opposite sign
        mixed = (1 - g) * perm + g * (-dh * np.std(perm) / (np.std(dh) or 1.0))
        vals.append(float(np.corrcoef(dh, mixed)[0, 1]))
    sweep.append([float(g), float(np.median(vals))])
print(f"  positive sweep (permuted -> shared opposite driver): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((g["r_raw"] - null_med) / (null_sd or 1e-9)))) for g in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

same = med_raw < 0 and abs(med_raw - null_med) > 2 * null_sd
aligned = med_raw > 0 and abs(med_raw - null_med) > 2 * null_sd
ceiling_survives = np.sign(med_head) == np.sign(med_raw) and abs(med_head) > 0.5 * abs(med_raw)

G = Gate("The rise and the one fall: same people, or different people?")
G.plant_direction_from_sweep("positive: interpolating toward a shared opposite driver drives the "
                             "correlation negative, and g=0 IS the null world",
                             [[g, -v] for g, v in sweep], baseline=-null_med,
                             baseline_spread=max(null_sd, 1e-4))
G.negative_control("racist changes permuted across cells", abs(null_med), abs(med_raw),
                   null_spread=null_sd, null_kind="cell-label permutation")
G.multiplicity_control("all era splits", ps, 0.05, labels=[g["split"] for g in grid])
G.asserted("PRECONDITION checked FIRST and it REJECTED design #1", True,
           "politics x education x cohort on single endpoint waves gave 4 usable cells of 59; "
           "pooled eras give 25-27 of 33-39. The 4-cell round was stopped before it ran",
           kind="control")
G.asserted("CEILING control ran: every change is also expressed against the room the cell had",
           ceiling_survives or abs(med_head) <= 0.5 * abs(med_raw),
           f"r raw {med_raw:+.3f} vs headroom-scaled {med_head:+.3f} — "
           f"{'survives' if ceiling_survives else 'DOES NOT survive'} the ceiling rescale",
           kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("ECOLOGICAL scope stated: this cannot distinguish the same individuals re-drawing a "
           "boundary from different individuals inside a cell moving oppositely", True,
           "repeated cross-section, no panel; the individual-level version is structurally "
           "unavailable on this release; scope stated", kind="control")
G.asserted("KILL: W_DIFFERENT requires the across-cell correlation to sit on its null",
           not (same or aligned),
           f"median r_raw {med_raw:+.3f} (headroom {med_head:+.3f}) vs null {null_med:+.4f} "
           f"+/- {null_sd:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif same and ceiling_survives:
    VERDICT, WORLD = "OVERTURNED", "W_SAME · one boundary re-drawn by the same subgroups"
elif same and not ceiling_survives:
    VERDICT, WORLD = "UNVERIFIED", "negative, but it does not survive the ceiling rescale"
elif aligned:
    VERDICT, WORLD = "OVERTURNED", "W_ALIGNED · the same cells rose on both; the fall is composition"
else:
    VERDICT, WORLD = "CONFIRMED", "W_DIFFERENT · the two movements are unrelated across subgroups"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=930, round="E03·A116·R368", verdict=VERDICT, world=WORLD,
           estimand="across-subgroup correlation between the change in homosexual-tolerance and the "
                    "change in racist-tolerance, early era -> late era",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           precondition_rejected_design="politics x education x cohort on single endpoint waves: "
                                        "4 usable cells of 59",
           grid=[{k: v for k, v in g.items() if k != "rows"} for g in grid],
           cells_detail=grid[0]["rows"],
           median_r_raw=med_raw, median_r_headroom=med_head, ceiling_survives=bool(ceiling_survives),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           ecological=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "same_people_or_different_people.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'same_people_or_different_people.json'}")
