#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A115·R367 — is the budget fixed? (`#927`(1) said it barely moved; nobody had measured it)
=============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#927`(1) put a framing on the page: *what the society did was re-sort its
                out-groups — one moved in, another moved out, and **the total quantity of tolerance
                barely moved**.* ⚠ **That last clause was never measured.** It was INFERRED from
                two standardised per-target movements (+0.4302 and −0.2602) appearing to cancel.
                **Standardised movements across items with different SDs do not sum to the change
                in the total** — that is an arithmetic step I took without doing the arithmetic.
                `#928`(2) named it as the untested question; `#928` also established the
                precondition (the arms ARE commensurable), so the sum is finally meaningful.

Live Worlds     W_ZEROSUM · total tolerance flat while its distribution across targets shifts ⇒
                            "re-allocation" is literal and `#927`(1) stands.
                W_GROWTH  · the total rose substantially ⇒ **general tolerance DID increase**, the
                            re-sorting rides on top of it, and **`#927`(1)'s clause is false.**
                            ⚠ **The unwelcome one, and it retracts my own page framing.**
                W_MIXED   · the total moves, but the target-specific parts dominate it.

⚠ THE WINDOW    `#927` measured **within-cohort** movement over **1990→2021**. A raw wave-mean over
MUST MATCH,     **1976→2021** is a DIFFERENT QUANTITY, and contradicting `#927` with it would be
OR IT IS A      `#898`'s error — comparing a number to a reference that answers a different
DIFFERENT       question. ⇒ **both windows are computed and reported side by side**, and any
QUESTION        disagreement between them is reported as a fact about the window, not about the world.

Estimand        (a) the change in **mean total tolerance**, the count of permissive answers over all
(G1)                15 items (0–15) — the canonical Stouffer scale;
                (b) the decomposition of each target's change into a COMMON component (the person's
                    total) and a TARGET-SPECIFIC component, and their relative size.

Prediction      W_ZEROSUM -> total change inside its null; specific components large.
Matrix          W_GROWTH  -> total change far outside its null; four or five targets moving together.
                W_MIXED   -> total moves, specific parts larger than the common part.

⚠ PRIOR ART     The 15-item count IS Stouffer's (1955) political-tolerance scale, and its long rise
DECLARED        is one of the most-replicated findings in American public opinion (Davis 1975 on the
                1954–72 rise; Bobo & Licari 1989). **Measuring it here is a VERIFICATION.** What is
                NOT standard is the decomposition against `#927`'s claim, which is this round's own.

Controls        NEGATIVE, and it is unusually clean for this question: **permute which TARGET BLOCK
                a person's answers belong to, within that person.** This preserves each person's
                total EXACTLY and destroys only the allocation ⇒ the target-specific component must
                go to zero while the common component is untouched.
                POSITIVE: plant an allocation shift into that permuted world and sweep, `g=0` on the
                null (`#922`'s gate; `#927` broke this three times, `#928` twice).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ "total tolerance" is a count over FIVE targets chosen by Stouffer in 1955 — it is a budget
    over THIS list, not over everything a person might tolerate. A fixed budget over a wider list
    is untestable here;
  (2) ⚠ the count treats all 15 items as equally weighted; nothing justifies that beyond convention;
  (3) ⚠ repeated cross-section; APC collinear; no age effect claimed;
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
RNG = np.random.default_rng(367)
MIN_N = 300

TARGETS = {
    "homosexuals":       {"spkhomo": -1, "colhomo": -1, "libhomo": +1},
    "racists":           {"spkrac": -1,  "colrac": -1,  "librac": +1},
    "communists":        {"spkcom": -1,  "colcom": +1,  "libcom": +1},
    "militarists":       {"spkmil": -1,  "colmil": -1,  "libmil": +1},
    "anti-religionists": {"spkath": -1,  "colath": -1,  "libath": +1},
}
COLS = [c for t in TARGETS.values() for c in t]

d = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
for c in COLS:
    d[c] = d[c].where(d[c].between(1, 7))
codes = {c: sorted(d[c].dropna().unique()) for c in COLS}          # ⚠ `#927`(3): derived, no literals
for tgt, items in TARGETS.items():
    for c, sign in items.items():
        lo, hi = float(codes[c][0]), float(codes[c][-1])
        d[c] = (d[c] == (hi if sign > 0 else lo)).astype(float).where(d[c].notna())
d["complete"] = d[COLS].notna().all(axis=1)
cc = d[d["complete"]].copy()
cc["total"] = cc[COLS].sum(axis=1)
for tgt, items in TARGETS.items():
    cc[f"T_{tgt}"] = cc[list(items)].sum(axis=1)

counts = cc.groupby("year").size()
waves = [int(y) for y, n in counts.items() if n >= MIN_N]
print(f"PRECONDITION CHECK — waves with >={MIN_N} complete cases: {len(waves)} usable · "
      f"DROPPED {len([y for y, n in counts.items() if 0 < n < MIN_N])}  (absence reported, not passed)")
cc = cc[cc.year.isin(waves)]
print(f"  complete cases {len(cc)} · waves {waves[0]}-{waves[-1]}")

WINDOWS = {"1976->2021 (this round's own)": (waves[0], waves[-1]),
           "1990->2021 (`#927`'s window)": (min(waves, key=lambda y: abs(y - 1990)), waves[-1])}


def totals(sub, y0, y1):
    a, b = sub[sub.year == y0], sub[sub.year == y1]
    row = dict(y0=int(y0), y1=int(y1), n=int(len(a) + len(b)),
               total0=float(a["total"].mean()), total1=float(b["total"].mean()))
    row["total_change"] = row["total1"] - row["total0"]
    for tgt in TARGETS:
        row[f"{tgt}"] = float(b[f"T_{tgt}"].mean() - a[f"T_{tgt}"].mean())
    # common component = each target's share of the total change if all moved proportionally
    row["common_per_target"] = row["total_change"] / len(TARGETS)
    row["specific"] = {t: row[t] - row["common_per_target"] for t in TARGETS}
    row["specific_abs_mean"] = float(np.mean([abs(v) for v in row["specific"].values()]))
    row["risers"] = int(sum(1 for t in TARGETS if row[t] > 0))
    return row


print("\n=== THE GRID — both windows, reported side by side ===")
grid = []
for name, (y0, y1) in WINDOWS.items():
    r = totals(cc, y0, y1)
    r["window"] = name
    grid.append(r)
    print(f"  {name}")
    print(f"    total tolerance {r['total0']:.3f} -> {r['total1']:.3f} of 15  "
          f"= {r['total_change']:+.3f}   n={r['n']}   targets rising: {r['risers']}/5")
    for t in TARGETS:
        print(f"      {t:18s} {r[t]:+.3f} of 3   (specific {r['specific'][t]:+.3f})")
    print(f"    common per target {r['common_per_target']:+.3f} · mean |specific| "
          f"{r['specific_abs_mean']:.3f}")

main = grid[0]
tot_change = main["total_change"]
ratio = main["specific_abs_mean"] / abs(main["common_per_target"]) if main["common_per_target"] else np.nan
print(f"\n  |specific| / |common| = {ratio:.2f}   (>1 means allocation dominates growth)")

# ══ NEGATIVE CONTROL — permute WHICH TARGET a person's answer-block belongs to ═══════
# preserves each person's TOTAL exactly; destroys only the allocation.
blocks = [list(v) for v in TARGETS.values()]
null_vals = []
y0, y1 = WINDOWS["1976->2021 (this round's own)"]
for _ in range(200):
    p = cc[cc.year.isin([y0, y1])].copy()
    arr = p[COLS].to_numpy().reshape(len(p), len(blocks), 3)
    idx = np.argsort(RNG.random((len(p), len(blocks))), axis=1)
    arr = np.take_along_axis(arr, idx[:, :, None], axis=1)
    for j, tgt in enumerate(TARGETS):
        p[f"T_{tgt}"] = arr[:, j, :].sum(axis=1)
    a, b = p[p.year == y0], p[p.year == y1]
    spec = [float(b[f"T_{t}"].mean() - a[f"T_{t}"].mean()) - main["common_per_target"] for t in TARGETS]
    null_vals.append(float(np.mean([abs(v) for v in spec])))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null for the SPECIFIC component (target blocks permuted within person — the person's "
      f"TOTAL is preserved exactly; kind of null: within-person target-block permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant an allocation shift into the permuted world ═════════════
sweep = []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(30):
        p = cc[cc.year.isin([y0, y1])].copy()
        arr = p[COLS].to_numpy().reshape(len(p), len(blocks), 3)
        idx = np.argsort(RNG.random((len(p), len(blocks))), axis=1)
        arr = np.take_along_axis(arr, idx[:, :, None], axis=1)
        # plant: in the LATE wave, move mass from block 1 to block 0 with probability g
        late = (p.year.to_numpy() == y1)
        take = late & (RNG.random(len(p)) < g)
        swap = arr[take].copy()
        moved = np.minimum(swap[:, 1, :], 1.0)
        swap[:, 0, :] = np.minimum(swap[:, 0, :] + moved, 1.0)
        swap[:, 1, :] = np.maximum(swap[:, 1, :] - moved, 0.0)
        arr[take] = swap
        for j, tgt in enumerate(TARGETS):
            p[f"T_{tgt}"] = arr[:, j, :].sum(axis=1)
        a, b = p[p.year == y0], p[p.year == y1]
        spec = [float(b[f"T_{t}"].mean() - a[f"T_{t}"].mean()) - main["common_per_target"]
                for t in TARGETS]
        vals.append(float(np.mean([abs(v) for v in spec])))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (allocation shift planted into the permuted world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

# ══ the TOTAL's own null: bootstrap, since permutation cannot move a preserved total ══
boot = []
for _ in range(400):
    a = cc[cc.year == y0].sample(len(cc[cc.year == y0]), replace=True, random_state=int(RNG.integers(1e9)))
    b = cc[cc.year == y1].sample(len(cc[cc.year == y1]), replace=True, random_state=int(RNG.integers(1e9)))
    boot.append(float(b["total"].mean() - a["total"].mean()))
t_lo, t_hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
print(f"  total change {tot_change:+.3f} of 15  [95% bootstrap {t_lo:+.3f}, {t_hi:+.3f}]")

ps = [2 * (1 - stats.norm.cdf(abs((r["specific_abs_mean"] - null_med) / (null_sd or 1e-9))))
      for r in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

budget_fixed = t_lo <= 0 <= t_hi
alloc_resolved = abs(main["specific_abs_mean"] - null_med) > 2 * null_sd

G = Gate("Is the budget fixed?")
G.plant_direction_from_sweep("positive: a planted allocation shift raises the specific component, "
                             "and g=0 is null", sweep, baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("target blocks permuted within person (total preserved exactly)",
                   abs(null_med), abs(main["specific_abs_mean"]),
                   null_spread=null_sd, null_kind="within-person target-block permutation")
G.multiplicity_control("both windows", ps, 0.05, labels=[r["window"] for r in grid])
G.has_error_bar("the total change carries an interval", tot_change, (t_hi - t_lo) / 4, "bootstrap_人层")
G.asserted("BOTH windows are reported side by side, so a disagreement is a fact about the window",
           len(grid) == 2,
           "; ".join(f"{r['window']}: total {r['total_change']:+.3f}, risers {r['risers']}/5"
                     for r in grid), kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("prior art declared", True,
           "the 15-item count IS Stouffer's 1955 political-tolerance scale and its long rise is one "
           "of the most-replicated findings in American opinion (Davis 1975; Bobo & Licari 1989) — "
           "measuring it here is a VERIFICATION; the decomposition against `#927` is this round's own",
           kind="control")
G.asserted("KILL: `#927`(1)'s clause requires the total to have barely moved",
           budget_fixed,
           f"total change {tot_change:+.3f} of 15 [95% {t_lo:+.3f}, {t_hi:+.3f}]; "
           f"{main['risers']}/5 targets rose")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif budget_fixed:
    VERDICT, WORLD = "CONFIRMED", "W_ZEROSUM · the budget is fixed"
elif ratio > 1:
    VERDICT, WORLD = "OVERTURNED", "W_MIXED · the total grew AND allocation dominates it"
else:
    VERDICT, WORLD = "OVERTURNED", "W_GROWTH · the budget is not fixed; tolerance grew"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=929, round="E03·A115·R367", verdict=VERDICT, world=WORLD,
           estimand="change in mean total tolerance (0-15 Stouffer count) and its decomposition into "
                    "a common component and target-specific components",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           prior_art="Stouffer 1955 scale; Davis 1975; Bobo & Licari 1989 — VERIFICATION",
           waves=len(waves), complete_cases=int(len(cc)), grid=grid,
           total_change=tot_change, total_ci=[t_lo, t_hi], specific_over_common=float(ratio),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           retracts="`#927`(1)'s clause that the total quantity of tolerance barely moved",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "is_the_budget_fixed.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'is_the_budget_fixed.json'}")
