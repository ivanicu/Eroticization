#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R355 — which ITEM carries the rank-1 departure, on a second instrument
===============================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#915` measured a rank-1 departure in the 2x2 [exposure x norm] matrix and read it
                as domain-specificity. `#916` measured that the departure pattern looks like ITEM
                BREADTH instead, returned `UNVERIFIED` because that site cannot run a placebo, and
                armed `#111c`: no third control on that file. ⚠ **But both rounds asked "what
                EXPLAINS the departure" while never asking "WHERE does it live".** If one item
                carries all of it, neither domain nor breadth is the right frame at all.

Why Now         `#916`(2) forbids another round on 2011-2013 female. This round changes BOTH the
                instrument (NSFG -> GSS, HARD RULE 4) and the question (what explains it -> which
                item carries it). And GSS supplies **the placebo NSFG structurally could not**: two
                norm items that neither exposure is tied to.

Live Worlds     W1 · **ONE ITEM** — removing the tautological item `homosex` collapses the
                     departure. Then `#915`/`#916` are about one item's self-endorsement, and the
                     {domain, breadth} decomposition was the wrong carve entirely.
                     **This is `#916`'s world C, it is the meta-separator, and it is unwelcome.**
                W2 · **DOMAIN** — removing either matched item leaves the departure high, because
                     the structure is two domains each with its own exposure.
                W3 · **DISTRIBUTED** — no single removal collapses it: a genuine second dimension.

Estimand        Generalises `#915`'s determinant rather than replacing it. Under ONE moral dial
(G1)            b(c -> i) = lambda_i * beta_c, so the 2xK coefficient matrix B[exposure, norm] is
                RANK 1 and its second singular value is **exactly zero, for any loadings**. The
                statistic is `sigma2/sigma1` (scale-free, so 2x4 and 2x3 are comparable), and
                LEAVE-ONE-ITEM-OUT identifies the carrier.
                ⚠ `#903`: a per-item MEAN OVER PAIRS is not per-item — one departing item lowers
                every pair it appears in. Leave-one-out on the MATRIX avoids exactly that, which is
                why the statistic is a singular value and not an average of determinants.

Prediction      W1 -> drop(`homosex`) collapses sigma2/sigma1; the other three drops do not.
Matrix          W2 -> drop(`homosex`) and drop(`xmarsex`) both collapse it; the placebo drops do not.
                W3 -> no drop collapses it.
                ⚠ No flat row: every world names a DIFFERENT drop as the collapsing one.

The instrument  GSS 1972-2024 (`gss7224_r3a.dta`), and HARD RULE 2 — **the shared-instrument threat
                is named**: norms and behaviours come from the SAME questionnaire, the SAME
                respondent, in the SAME interview, so a respondent who under-reports a behaviour may
                also report a stricter norm. That coupling is not removable here and is registered.
                  norms     `homosex` `premarsx` `teensex` `xmarsex`   (1 always wrong .. 4 not wrong)
                  exposures `ss` own same-sex partner since 18, built from `nummen`/`numwomen`
                            CROSSED WITH `sex` — ⚠ `sexsex` alone cannot say "same-sex" without
                            knowing the respondent's sex, so a name is not a measurement here;
                            `xm` `evstray == 1`, restricted to ever-married (`evstray` in {1,2};
                            code 3 is "never married" and is NOT a control, it is not at risk).
                  matched pairs  ss -> `homosex` · xm -> `xmarsex`
                  PLACEBO        `premarsx` + `teensex` — neither exposure's own domain

Stopping Rule   One pass over the drop x adjustment x exposure-definition x support grid, published
                whole. If two worlds' drops both collapse, that is W3 and it is reported as W3.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ shared instrument, above — norms and behaviours are one questionnaire;
  (2) cross-sectional per respondent; the arrow is not identified (inherited from `#915`);
  (3) `teensex` starts 1986 and `evstray` 1991, so the 2x4 lives on the intersection of waves and
    is NOT the full 1972-2024 series — the achieved wave set is printed, never assumed;
  (4) "tautological" is a claim about item wording, not a measured property;
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
RNG = np.random.default_rng(355)

NORMS = ["homosex", "premarsx", "teensex", "xmarsex"]
PLACEBO = ["premarsx", "teensex"]                 # neither exposure's own domain
COLS = NORMS + ["evstray", "nummen", "numwomen", "sexsex", "sex", "year", "age", "educ", "race", "attend"]

df = pd.read_stata(GSS, columns=COLS, convert_categoricals=False)
print(f"GSS rows {len(df)}")

# ── HARD RULE 1: build every derived column explicitly and print what it achieved ──
for c in NORMS:
    df[c] = df[c].where(df[c].isin([1, 2, 3, 4]))     # `homosex` also carries a code 5 (n=82) -> drop
# same-sex exposure needs the respondent's SEX; `sexsex` alone cannot express it
ss_men = (df["sex"] == 1) & (df["nummen"] > 0)
ss_wom = (df["sex"] == 2) & (df["numwomen"] > 0)
have_counts = df["nummen"].notna() & df["numwomen"].notna() & df["sex"].isin([1, 2])
df["ss"] = np.where(ss_men | ss_wom, 1.0, 0.0)
df.loc[~have_counts, "ss"] = np.nan
# `sexsex` alternative definition (spec-curve axis): 2 = both sexes, plus same-sex-only by `sex`
df["ss_alt"] = np.where(df["sexsex"] == 2, 1.0,
                        np.where((df["sex"] == 1) & (df["sexsex"] == 1), 1.0,
                                 np.where((df["sex"] == 2) & (df["sexsex"] == 3), 1.0, 0.0)))
df.loc[~df["sexsex"].isin([1, 2, 3]), "ss_alt"] = np.nan
# extramarital exposure: code 3 = never married = NOT AT RISK, excluded rather than coded 0
df["xm"] = np.where(df["evstray"] == 1, 1.0, np.where(df["evstray"] == 2, 0.0, np.nan))

inv = {}
for c in NORMS + ["ss", "ss_alt", "xm"]:
    s = df[c]
    yrs = sorted(df.loc[s.notna(), "year"].unique())
    inv[c] = dict(non_missing=int(s.notna().sum()), waves=len(yrs),
                  first=int(yrs[0]) if yrs else None, last=int(yrs[-1]) if yrs else None,
                  rate=(float(s.mean()) if c in ("ss", "ss_alt", "xm") else None))
    print(f"  {c:10s} n={int(s.notna().sum()):6d}  waves={len(yrs):3d} "
          f"{yrs[0] if yrs else '-'}..{yrs[-1] if yrs else '-'}"
          + (f"  rate={float(s.mean()):.4f}" if c in ("ss", "ss_alt", "xm") else ""))

core = df.dropna(subset=NORMS + ["ss", "xm"])
waves = sorted(core["year"].unique())
print(f"\n  the 2x4 lives on the INTERSECTION: n={len(core)}, {len(waves)} waves "
      f"{int(waves[0])}..{int(waves[-1])}  <- not the full 1972-2024 series")

ACHIEVED = {}


def match_rate(marked, target, tag):
    """`#916`(3): match by SUBSAMPLING and PRINT THE ACHIEVED RATE, never the requested one.
    Solving (M-k)/(N-k) = r gives k = (M - r*N)/(1 - r); thresholding a count cannot hit a rate."""
    idx = marked[marked == 1].index.to_numpy()
    N, M = float(marked.notna().sum()), float(len(idx))
    k = int(round(max((M - target * N) / (1.0 - target), 0.0)))
    out = marked.copy()
    if 0 < k <= len(idx):
        out.loc[RNG.choice(idx, size=k, replace=False)] = np.nan
    ACHIEVED[tag] = float(out.mean())
    return out


def beta(y, x, covs):
    ok = y.notna() & x.notna()
    for c in covs:
        ok &= c.notna()
    n = int(ok.sum())
    if n < 100:
        return np.nan, n
    Y = y[ok].to_numpy(float)
    X = np.column_stack([np.ones(n), x[ok].to_numpy(float)] + [c[ok].to_numpy(float) for c in covs])
    for j in range(1, X.shape[1]):
        sd = X[:, j].std()
        if sd > 0:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    Y = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        b, *_ = np.linalg.lstsq(X, Y, rcond=None)
        return float(b[1]), n
    except np.linalg.LinAlgError:
        return np.nan, n


def rank1_departure(frame, items, expos, adjust):
    """sigma2/sigma1 of the 2xK coefficient matrix. EXACTLY 0 under any one-factor model."""
    covs = [frame[c] for c in adjust if c in frame.columns]
    B, ns = [], []
    for e in expos:
        row = []
        for it in items:
            b, n = beta(frame[it], frame[e], covs)
            row.append(b)
            ns.append(n)
        B.append(row)
    B = np.array(B, float)
    if np.isnan(B).any():
        return np.nan, None, 0
    sv = np.linalg.svd(B, compute_uv=False)
    if sv[0] <= 0:
        return np.nan, B.tolist(), int(min(ns))
    return float(sv[1] / sv[0]), B.tolist(), int(min(ns))


ADJ_SETS = {"raw": [], "demog": ["age", "educ", "race"], "demog+relig": ["age", "educ", "race", "attend"]}
EXPO_DEFS = {"counts": ["ss", "xm"], "sexsex": ["ss_alt", "xm"]}

grid = []
for ename, expos in EXPO_DEFS.items():
    frame0 = df.dropna(subset=NORMS + expos)
    base = float(frame0[expos[0]].mean())
    for mname in ("unmatched", "prev-matched"):
        frame = frame0.copy()
        if mname == "prev-matched":
            frame[expos[1]] = match_rate(frame[expos[1]], base, f"{ename}|xm")
            frame = frame.dropna(subset=expos)
        # ⚠⚠ `#916`(3), COMMITTED AGAIN ONE ROUND AFTER REGISTERING IT. v1 recorded the achieved
        #   rate INSIDE `match_rate`, i.e. BEFORE the `dropna` on the next line — so the number in
        #   the artifact was not the rate the analysis actually saw. **A control's achieved value
        #   must be measured on the frame that is ANALYSED, not on the frame it was computed from.**
        #   The error class is unchanged and the instance is new: I measured the right quantity at
        #   the wrong POINT. ⇒ recorded here, on the final frame, for both arms.
        ACHIEVED[f"{ename}|{mname}|target_ss"] = float(frame[expos[0]].mean())
        ACHIEVED[f"{ename}|{mname}|achieved_xm"] = float(frame[expos[1]].mean())
        for aname, aset in ADJ_SETS.items():
            for drop in [None] + NORMS:
                items = [i for i in NORMS if i != drop]
                r, B, n = rank1_departure(frame, items, expos, aset)
                grid.append(dict(expo_def=ename, prevalence=mname, adjust=aname,
                                 dropped=drop or "none", k=len(items), ratio=r, n=n, coefs=B))
            # the PLACEBO: the two items neither exposure is tied to
            r, B, n = rank1_departure(frame, PLACEBO, expos, aset)
            grid.append(dict(expo_def=ename, prevalence=mname, adjust=aname,
                             dropped="PLACEBO(premarsx+teensex)", k=2, ratio=r, n=n, coefs=B))

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    r = "  nan " if g["ratio"] is None or np.isnan(g["ratio"]) else f"{g['ratio']:.4f}"
    print(f"  {g['expo_def']:7s} {g['prevalence']:12s} {g['adjust']:12s} "
          f"drop={g['dropped']:26s} k={g['k']} sigma2/sigma1={r}  n={g['n']:6d}")

full = [g["ratio"] for g in grid if g["dropped"] == "none" and not np.isnan(g["ratio"])]
med_full = float(np.median(full)) if full else np.nan
drops = {}
for it in NORMS:
    v = [g["ratio"] for g in grid if g["dropped"] == it and not np.isnan(g["ratio"])]
    drops[it] = float(np.median(v)) if v else np.nan
plac = [g["ratio"] for g in grid if g["dropped"].startswith("PLACEBO") and not np.isnan(g["ratio"])]
med_plac = float(np.median(plac)) if plac else np.nan

print(f"\n  full 2x4              sigma2/sigma1 = {med_full:.4f}")
for it in NORMS:
    tag = "  <- matched pair" if it in ("homosex", "xmarsex") else ""
    print(f"  drop {it:10s}      sigma2/sigma1 = {drops[it]:.4f}   "
          f"({drops[it] / med_full:.2f}x the full value){tag}")
print(f"  PLACEBO (premarsx+teensex) sigma2/sigma1 = {med_plac:.4f}")

# ══ NEGATIVE CONTROL — one-factor synthetic at the observed marginals ════════════════
n_syn = int(len(core))
lam = [0.8, 0.6, 0.5, 0.7]


def synth(n, g, reps):
    out = []
    for _ in range(reps):
        th = RNG.standard_normal(n)
        dom = {0: RNG.standard_normal(n), 1: RNG.standard_normal(n)}
        cols = {}
        for j, it in enumerate(NORMS):
            d = dom[0] if it == "homosex" else (dom[1] if it == "xmarsex" else 0.0)
            z = lam[j] * th + g * d + np.sqrt(max(1 - lam[j] ** 2, 1e-9)) * RNG.standard_normal(n)
            cols[it] = pd.Series(np.clip(np.digitize(z, np.quantile(z, [.25, .5, .75])) + 1, 1, 4))
        for j, e in enumerate(("ss", "xm")):
            z = 0.5 * th + g * dom[j] + RNG.standard_normal(n)
            cols[e] = pd.Series((z > np.quantile(z, 0.8)).astype(float))
        out.append(pd.DataFrame(cols))
    return out


NUL = [x for x, _, _ in (rank1_departure(d, NORMS, ["ss", "xm"], []) for d in synth(n_syn, 0.0, 60))
       if not np.isnan(x)]
null_med, null_sd = float(np.median(NUL)), float(np.std(NUL))
print(f"\n  one-factor null: median {null_med:.4f}  sd {null_sd:.4f}  (n={n_syn}, {len(NUL)} reps)")

sweep = []
for g in (0.0, 0.15, 0.3, 0.45, 0.6):
    v = [x for x, _, _ in (rank1_departure(d, NORMS, ["ss", "xm"], []) for d in synth(n_syn, g, 20))
         if not np.isnan(x)]
    sweep.append((g, float(np.median(v)) if v else np.nan))
print(f"  positive sweep (g, median sigma2/sigma1): {[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((g["ratio"] - null_med) / (null_sd or 1e-9))))
      for g in grid if not np.isnan(g["ratio"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

# collapse = the drop that brings the departure back toward its one-factor null
def collapsed(v):
    return (v - null_med) < 0.5 * (med_full - null_med)


carriers = [it for it in NORMS if not np.isnan(drops[it]) and collapsed(drops[it])]

G = Gate("Which ITEM carries the rank-1 departure?")
G.plant_direction_from_sweep("positive: planted domain coupling raises sigma2/sigma1, g=0 is null",
                             sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("synthetic ONE-FACTOR world at the observed n and marginals",
                   abs(null_med), abs(med_full), null_spread=null_sd,
                   null_kind="one-factor latent, 4 items, matched loadings")
G.multiplicity_control("the whole drop grid", ps, 0.05,
                       labels=[f"{g['expo_def']}|{g['prevalence']}|{g['adjust']}|{g['dropped']}"
                               for g in grid if not np.isnan(g["ratio"])])
# `#916`: an absent placebo must FAIL, never pass
G.asserted("placebo actually RAN (absence is not a pass)", not np.isnan(med_plac),
           f"placebo sigma2/sigma1 {med_plac} over {len(plac)} cells", kind="control")
G.asserted("placebo: the two items tied to NEITHER exposure stay near the one-factor null",
           not np.isnan(med_plac) and abs(med_plac - null_med) < abs(med_full - null_med),
           f"placebo {med_plac:.4f} vs full {med_full:.4f}, null {null_med:.4f}", kind="control")
# `#916`(3): print the ACHIEVED rate, never the requested one
# ⚠ v1 of this row asserted `bool(ACHIEVED)` — non-empty dict — which is a check that cannot fail.
#   It now compares the ACHIEVED rate to the TARGET on the analysed frame, in both exposure
#   definitions, and fails if either misses.
_gaps = {k.rsplit("|", 1)[0]: abs(ACHIEVED[f"{k.rsplit('|', 1)[0]}|achieved_xm"]
                                 - ACHIEVED[f"{k.rsplit('|', 1)[0]}|target_ss"])
         for k in ACHIEVED if k.endswith("|achieved_xm") and "prev-matched" in k}
G.asserted("prevalence control ACHIEVED its target on the ANALYSED frame, measured not intended",
           bool(_gaps) and all(v < 0.02 for v in _gaps.values()),
           f"|achieved - target| per arm {({k: round(v, 4) for k, v in _gaps.items()})}; "
           f"all rates {({k: round(v, 4) for k, v in ACHIEVED.items()})}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("KILL: W1 (one item) requires exactly ONE drop to collapse the departure",
           len(carriers) != 1,
           f"drops that collapse it: {carriers or 'none'}; full {med_full:.4f}, "
           f"null {null_med:.4f}, per-item {({k: round(v, 4) for k, v in drops.items()})}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif len(carriers) == 1:
    VERDICT, WORLD = "OVERTURNED", f"W1 · ONE ITEM (`{carriers[0]}`) — the decomposition was wrong"
elif set(carriers) >= {"homosex", "xmarsex"}:
    VERDICT, WORLD = "CONFIRMED", "W2 · DOMAIN — both matched items carry it"
elif not carriers:
    VERDICT, WORLD = "CONFIRMED", "W3 · DISTRIBUTED — no single item carries it"
else:
    VERDICT, WORLD = "UNVERIFIED", f"mixed carriers {carriers}"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=917, round="E03·A112·R355", verdict=VERDICT, world=WORLD,
           estimand="sigma2/sigma1 of the 2xK [exposure x norm] coefficient matrix, exactly 0 "
                    "under any one-factor model; leave-one-item-out identifies the carrier",
           instrument="GSS 1972-2024 gss7224_r3a.dta", shared_instrument_threat=True,
           inventory=inv, intersection_n=int(len(core)), intersection_waves=[int(w) for w in waves],
           grid=grid, full=med_full, drops=drops, placebo=med_plac, carriers=carriers,
           null_median=null_med, null_sd=null_sd, null_reps=len(NUL), positive_sweep=sweep,
           achieved_rates=ACHIEVED, family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "which_item_carries_it.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'which_item_carries_it.json'}")
