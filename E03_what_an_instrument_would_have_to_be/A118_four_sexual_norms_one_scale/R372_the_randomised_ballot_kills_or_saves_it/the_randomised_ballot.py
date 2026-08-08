#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A118·R372 — the last alternative to `#933`, and GSS randomised it for me in 1988
====================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#933` is the project's strongest claim — on one identical 1-4 frame, within cohorts,
                `homosex` moved +0.6400 while `premarsx` +0.1453, `teensex` +0.1652, `xmarsex`
                +0.0906 — and it left **exactly one alternative standing**: all four are asked **in
                one interview, in a row**, so respondents may anchor them against each other
                (`#933`(2)). No design in `#933` removed it.

⚠ AND GSS      `#932`(1) says check COVERAGE before believing a test is available. Checking it found
RANDOMISED     something better than availability: **GSS randomly assigns respondents to BALLOTS, and
THE            the ballots split this very block.**
CONFOUND         **ballot 1** (n=17,355) — all four items: the full block, maximum anchoring
                 **ballot 2** (n=17,544) — `premarsx` + `teensex` only
                 **ballot 3** (n=17,407) — `xmarsex` + `homosex` only: **no `premarsx`/`teensex` to
                                            anchor against**
                ⇒ **`homosex` is asked WITH the rest of the block and WITHOUT it, in the same waves,
                by random assignment.** That is the confound `#933`(2) named, manipulated.

⚠ A SCOPE      **`#933` used respondents who answered ALL FOUR — which is ballot 1 and only ballot 1**
FACT ABOUT     (n = 15,056). **Ballot 3 is an independent, randomised, never-used replication arm of
`#933`         n = 15,598 for `homosex` and 16,152 for `xmarsex`.**

Live Worlds    W_ANCHORING · `homosex`'s within-cohort movement is markedly smaller on ballot 3 ⇒
                             `#933`'s gap was question context. ⚠ **Unwelcome: it deflates the
                             project's strongest claim one round after it was made.**
               W_ROBUST    · the movement is the same on both ⇒ the last alternative dies and
                             `#933` stands on a randomised replication.
               W_MIXED     · the level shifts but the movement does not, or vice versa — which would
                             mean context sets WHERE people answer without changing HOW THEY CHANGE.
                             (the meta-separator: "anchoring" would then be the wrong single word)

Estimand       within-cohort change in `homosex` and in `xmarsex`, early era -> late era, computed
(G1)           SEPARATELY BY BALLOT; then the ballot-1 minus ballot-3 difference of those changes.
               Same scale, same waves, same cohorts, **random assignment** — so a difference cannot
               be composition, era, scale or population.

Prediction     W_ANCHORING -> ballot1 movement >> ballot3 movement for `homosex`.
Matrix         W_ROBUST    -> the two agree within the ballot-permutation null.
               W_MIXED     -> levels differ, movements do not.

⚠ RANDOMISATION Checked and PRINTED before the estimator, because a "randomised" arm that is not
CHECK FIRST    balanced is just another observational contrast: age · education · sex · birth cohort
               compared across ballots.

Controls       NEGATIVE: permute the BALLOT label — it really is random, so this is the exact null.
               POSITIVE: plant a ballot-dependent movement INTO the permuted world and sweep, so
               `g=0` lands ON the null (`#922`'s gate).
               ⚠ PLACEBO: `xmarsex` is on BOTH ballot 1 and ballot 3 as well, and `#933` says it
               barely moved. If the ballot contrast shows up in `xmarsex` too, it is a ballot
               artifact rather than anything about anchoring the *target*.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ this manipulates WHICH ITEMS ACCOMPANY the target, not the order within a pair — `xmarsex`
    and `homosex` are still asked together on ballot 3, so a two-item anchoring effect survives
    untested and is NOT ruled out;
  (2) ⚠ repeated cross-section: a cohort is the same birth years, never the same people;
  (3) ⚠ APC collinear; no age effect claimed;
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
RNG = np.random.default_rng(372)
TARGET, PLACEBO = "homosex", "xmarsex"
MIN_CELL = 60

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort", "age", "educ", "sex", TARGET, PLACEBO],
                  convert_categoricals=False)
for c in (TARGET, PLACEBO):
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot.isin([1, 3]) & d.cohort.between(1900, 2006)].copy()
d["grp"] = (d.cohort // 10 * 10).astype(int)

# ══ COVERAGE, non-null, not rows (`#932`(1)) ═════════════════════════════════════════
print("COVERAGE (non-null, not rows):")
for b in (1, 3):
    s = d[d.ballot == b]
    yrs = sorted(s.year.unique())
    print(f"  ballot {b}: n={len(s):6d} · {TARGET} {int(s[TARGET].notna().sum()):5d} · "
          f"{PLACEBO} {int(s[PLACEBO].notna().sum()):5d} · waves {len(yrs)} "
          f"{int(yrs[0])}-{int(yrs[-1])}")
print("  ⚠ `#933` used respondents answering ALL FOUR — that is ballot 1 only. Ballot 3 is an "
      "independent randomised arm it never touched.")

# ══ RANDOMISATION CHECK, printed BEFORE the estimator ════════════════════════════════
print("\nRANDOMISATION CHECK — a randomised arm that is not balanced is just another observational "
      "contrast:")
bal = {}
for v in ("age", "educ", "sex", "cohort"):
    m = d.groupby("ballot")[v].mean()
    sd = d[v].std()
    bal[v] = dict(b1=float(m[1]), b3=float(m[3]), std_diff=float(abs(m[1] - m[3]) / (sd or 1)))
    print(f"  {v:7s} b1={m[1]:.3f} b3={m[3]:.3f}  standardised diff {bal[v]['std_diff']:.4f}")
worst = max(v["std_diff"] for v in bal.values())

SPLITS = {"1988-1998 vs 2014-2024": ((1988, 1998), (2014, 2024)),
          "1988-2000 vs 2012-2024": ((1988, 2000), (2012, 2024)),
          "1988-2002 vs 2010-2024": ((1988, 2002), (2010, 2024))}


def moves(frame, early, late, item):
    e = frame.copy()
    e["era"] = np.where(e.year.between(*early), "E", np.where(e.year.between(*late), "L", None))
    e = e[e.era.notna() & e[item].notna()]
    out = []
    for k in sorted(e.grp.unique()):
        a, b = e[(e.grp == k) & (e.era == "E")], e[(e.grp == k) & (e.era == "L")]
        if len(a) < MIN_CELL or len(b) < MIN_CELL:
            continue
        out.append(dict(cohort=int(k), move=float(b[item].mean() - a[item].mean()),
                        n=int(len(a) + len(b))))
    return out


grid = []
for name, (early, late) in SPLITS.items():
    row = dict(split=name)
    ok = True
    for item in (TARGET, PLACEBO):
        for b in (1, 3):
            mv = moves(d[d.ballot == b], early, late, item)
            if len(mv) < 3:
                ok = False
            row[f"{item}_b{b}"] = float(np.mean([m["move"] for m in mv])) if mv else np.nan
            row[f"{item}_b{b}_k"] = len(mv)
            row[f"{item}_b{b}_n"] = int(sum(m["n"] for m in mv))
    if not ok:
        print(f"  ⚠ DROPPED split {name}: fewer than 3 cohorts in some ballot arm")
        continue
    row["target_ballot_diff"] = row[f"{TARGET}_b1"] - row[f"{TARGET}_b3"]
    row["placebo_ballot_diff"] = row[f"{PLACEBO}_b1"] - row[f"{PLACEBO}_b3"]
    row["n"] = row[f"{TARGET}_b1_n"] + row[f"{TARGET}_b3_n"]
    grid.append(row)

print("\n=== THE GRID — within-cohort movement, by ballot (all cells, disagreeing ones included) ===")
for r in grid:
    print(f"  {r['split']:24s} {TARGET}: b1 {r[f'{TARGET}_b1']:+.4f} (k={r[f'{TARGET}_b1_k']}) · "
          f"b3 {r[f'{TARGET}_b3']:+.4f} (k={r[f'{TARGET}_b3_k']})  diff {r['target_ballot_diff']:+.4f}")
    print(f"  {'':24s} {PLACEBO}: b1 {r[f'{PLACEBO}_b1']:+.4f} · b3 {r[f'{PLACEBO}_b3']:+.4f}  "
          f"diff {r['placebo_ballot_diff']:+.4f}   <- placebo   n={r['n']}")

med_t = float(np.median([r["target_ballot_diff"] for r in grid]))
med_p = float(np.median([r["placebo_ballot_diff"] for r in grid]))
med_b1 = float(np.median([r[f"{TARGET}_b1"] for r in grid]))
med_b3 = float(np.median([r[f"{TARGET}_b3"] for r in grid]))
print(f"\n  {TARGET} movement: with the full block (b1) {med_b1:+.4f} · without it (b3) {med_b3:+.4f}")
print(f"  ballot difference  target {med_t:+.4f}   placebo {med_p:+.4f}")

# ══ NEGATIVE CONTROL — permute the BALLOT label (it really is random) ════════════════
early, late = SPLITS["1988-1998 vs 2014-2024"]
null_vals = []
for _ in range(300):
    p = d.copy()
    p["ballot"] = RNG.permutation(p["ballot"].to_numpy())
    a = moves(p[p.ballot == 1], early, late, TARGET)
    b = moves(p[p.ballot == 3], early, late, TARGET)
    if len(a) >= 3 and len(b) >= 3:
        null_vals.append(float(np.mean([m["move"] for m in a]) - np.mean([m["move"] for m in b])))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (BALLOT label permuted — the exact null, because assignment really is random; "
      f"kind of null: randomised-arm label permutation): {null_med:+.4f} +/- {null_sd:.4f} "
      f"over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant a ballot-dependent movement INTO the permuted world ═════
sweep = []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for _ in range(40):
        p = d.copy()
        p["ballot"] = RNG.permutation(p["ballot"].to_numpy())
        late_b1 = (p.ballot == 1) & p.year.between(*late)
        p.loc[late_b1, TARGET] = np.clip(p.loc[late_b1, TARGET] + g, 1, 4)
        a = moves(p[p.ballot == 1], early, late, TARGET)
        b = moves(p[p.ballot == 3], early, late, TARGET)
        if len(a) >= 3 and len(b) >= 3:
            vals.append(float(np.mean([m["move"] for m in a]) - np.mean([m["move"] for m in b])))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (ballot-dependent movement planted into the permuted world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["target_ballot_diff"] - null_med) / (null_sd or 1e-9))))
      for r in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

anchoring = abs(med_t - null_med) > 2 * null_sd and med_t > 0
placebo_clean = abs(med_p - null_med) <= 2 * null_sd

G = Gate("Does `homosex` move less when the rest of the block is not there?")
G.plant_direction_from_sweep("positive: a planted ballot-dependent movement raises the ballot "
                             "difference, and g=0 is null", sweep,
                             baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("ballot label permuted", abs(null_med), abs(med_t),
                   null_spread=null_sd, null_kind="randomised-arm label permutation")
G.multiplicity_control("all era splits", ps, 0.05, labels=[r["split"] for r in grid])
G.asserted("RANDOMISATION CHECK ran and the arms are balanced", worst < 0.05,
           "; ".join(f"{k} std-diff {v['std_diff']:.4f}" for k, v in bal.items()), kind="control")
G.asserted("PLACEBO: `xmarsex` is on both ballots too — a ballot artifact would show up there as well",
           placebo_clean,
           f"placebo ballot difference {med_p:+.4f} vs null {null_med:+.4f} +/- {null_sd:.4f}",
           kind="control")
G.asserted("COVERAGE checked non-null, not rows (`#932`(1))", True,
           f"ballot 1 {TARGET} n={int(d[(d.ballot==1)][TARGET].notna().sum())}, "
           f"ballot 3 n={int(d[(d.ballot==3)][TARGET].notna().sum())}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("KILL: W_ANCHORING requires `homosex` to move LESS without the rest of the block",
           not anchoring,
           f"{TARGET} movement b1 {med_b1:+.4f} vs b3 {med_b3:+.4f}, difference {med_t:+.4f}; "
           f"null {null_med:+.4f} +/- {null_sd:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif anchoring:
    VERDICT, WORLD = "OVERTURNED", "W_ANCHORING · `#933`'s gap depends on the question context"
else:
    VERDICT, WORLD = "CONFIRMED", "W_ROBUST · the movement survives removing the rest of the block"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=934, round="E03·A118·R372", verdict=VERDICT, world=WORLD,
           estimand="within-cohort movement of `homosex` and `xmarsex` by randomly-assigned GSS "
                    "BALLOT, and the ballot-1 minus ballot-3 difference",
           instrument="GSS 1972-2024 gss7224_r3a.dta, `ballot` split-questionnaire assignment",
           design="ballot 1 = all four items; ballot 3 = xmarsex+homosex only; randomly assigned",
           scope_of_933="`#933` used respondents answering ALL FOUR, i.e. ballot 1 ONLY (n=15,056); "
                        "ballot 3 is an independent randomised arm it never touched",
           randomisation=bal, worst_std_diff=worst, grid=grid,
           target_b1=med_b1, target_b3=med_b3, target_ballot_diff=med_t, placebo_ballot_diff=med_p,
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "the_randomised_ballot.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'the_randomised_ballot.json'}")
