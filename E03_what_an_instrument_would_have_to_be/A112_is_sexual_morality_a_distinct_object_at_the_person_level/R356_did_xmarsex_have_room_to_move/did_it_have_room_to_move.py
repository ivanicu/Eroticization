#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R356 — did `xmarsex` carry none of the departure, or did it have no room to?
=====================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#917` localised the whole rank-1 departure to `homosex` and read it as: doing and
                approving travel together for same-sex behaviour but not for adultery. ⚠ **The
                comparison that carries the whole claim is `xmarsex`, and `xmarsex` is 77.4% on one
                endpoint with sd 0.7069 against `homosex`'s 1.4141 — exactly half.** An item pinned
                against its floor cannot express a graded association no matter what is true about
                the people answering it. **So the finding may be about how much ROOM each item had.**

Why Now         `#917` is one round old and already on the front page, where it replaced two
                readings that each lasted one round. **The cheapest thing that could destroy it goes
                first**, and this one is measurable rather than a reading of item text — which is
                what `#916`'s breadth rival was not.

⚠ WHY NOT      `#917`(2) proposed testing DISCLOSURE COUPLING by asking whether the departure shrinks
`#917`(2)       across 1991-2014 as disclosure stigma falls. **That design is dead on the gradient
                check, before any compute: a REAL coupling also attenuates over those decades**, as
                approval approaches its ceiling and the item's variance collapses. Both worlds
                predict shrinkage ⇒ **flat row** ⇒ the action separates nothing. `frontier` §1.2 says
                kill the action, not the world. **The disclosure world stays live and unaddressed**,
                and that is registered below rather than quietly dropped.

Live Worlds     W1 · **REAL** — the coupling is specific to `homosex` as an ACT/identity, and
                     survives once every item is scored against what it could attain.
                W2 · **CEILING ARTIFACT** — once each coefficient is expressed as a fraction of the
                     maximum its own marginal permits, `xmarsex` carries a departure too, or
                     `homosex` stops. `#917`'s psychology is then about item marginals.
                ⚠ BASIN NOTE: `#916` deflated `#915`, `#917` deflated `#916`. Three consecutive
                deflations is itself a basin, so the outcome I should find unwelcome HERE is **W1** —
                it would end the chain and force me to accept a psychological claim. The design is
                two-sided on purpose: normalisation can raise `xmarsex` OR leave it flat.

Estimand        Unchanged statistic, corrected UNITS. `sigma2/sigma1` of the 2xK coefficient matrix
(G1)            is exactly 0 under any one-factor model. ⚠ `beta()` already standardises the outcome,
                so raw VARIANCE is divided out — the residual threat is not variance but **ATTAINABLE
                RANGE**: a lumpy marginal caps the correlation a binary exposure can reach with it.
                So each cell is divided by its own **attainable ceiling**, computed as the
                correlation the item would show against a binary predictor of the SAME base rate
                that marks exactly the top of that item's own distribution. **That ceiling is a
                property of (marginal, base rate) alone — it contains no association.**

Prediction      W1 -> ceiling-normalised: `homosex` still the unique carrier.
Matrix          W2 -> ceiling-normalised: `xmarsex` joins it, or `homosex` drops out.

Stopping Rule   One pass over the ceiling-normalisation x adjustment x exposure-definition grid,
                published whole, with every item's ACHIEVED ceiling printed (`#916`(3): measure the
                control, never state its intent).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **the DISCLOSURE world is untested**, and this round does not touch it — `#917`(2)'s design
    was refuted on the gradient check and no replacement exists on this release. It would need a
    behaviour measure not self-reported to the same interviewer;
  (2) shared instrument (HARD RULE 2), inherited and unchanged;
  (3) the arrow is not identified; cross-sectional per respondent;
  (4) the ceiling is computed against a BINARY predictor because both exposures are binary — it does
    not bound what a continuous predictor could reach;
  (5) ⚠ **only this one instrument**, and here it is forced by the question rather than chosen: the
    object under test is whether `#917`'s GSS result is an artifact of THE GSS ITEMS' OWN MARGINALS.
    A second release has different items with different marginals, so it cannot answer this — it
    would be a new claim, not a replication. The cross-instrument move already happened at `#917`
    (NSFG -> GSS); this round is that result's own control.
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
RNG = np.random.default_rng(356)

NORMS = ["homosex", "premarsx", "teensex", "xmarsex"]
PLACEBO = ["premarsx", "teensex"]
COLS = NORMS + ["evstray", "nummen", "numwomen", "sexsex", "sex", "year", "age", "educ", "race", "attend"]

df = pd.read_stata(GSS, columns=COLS, convert_categoricals=False)
for c in NORMS:
    df[c] = df[c].where(df[c].isin([1, 2, 3, 4]))
have = df["nummen"].notna() & df["numwomen"].notna() & df["sex"].isin([1, 2])
df["ss"] = np.where(((df["sex"] == 1) & (df["nummen"] > 0)) | ((df["sex"] == 2) & (df["numwomen"] > 0)), 1.0, 0.0)
df.loc[~have, "ss"] = np.nan
df["ss_alt"] = np.where(df["sexsex"] == 2, 1.0,
                        np.where((df["sex"] == 1) & (df["sexsex"] == 1), 1.0,
                                 np.where((df["sex"] == 2) & (df["sexsex"] == 3), 1.0, 0.0)))
df.loc[~df["sexsex"].isin([1, 2, 3]), "ss_alt"] = np.nan
df["xm"] = np.where(df["evstray"] == 1, 1.0, np.where(df["evstray"] == 2, 0.0, np.nan))

# ══ HARD RULE 1 — the marginals that this whole round is about ════════════════════════
core = df.dropna(subset=NORMS + ["ss", "xm"])
marg = {}
print(f"the 2x4 intersection n={len(core)}")
for c in NORMS:
    sh = core[c].value_counts(normalize=True).sort_index()
    marg[c] = dict(sd=float(core[c].std()), max_cell=float(sh.max()),
                   shares={str(int(k)): float(v) for k, v in sh.items()})
    print(f"  {c:9s} sd={core[c].std():.4f}  max-cell={sh.max():.3f}  shares "
          f"{ {int(k): round(v, 3) for k, v in sh.items()} }")


def attainable_ceiling(y, rate):
    """Max |corr(y, x)| over binary x with P(x=1)=rate: mark exactly the top `rate` of y.
    A property of (marginal, base rate) ALONE — it contains no association whatsoever."""
    v = y.dropna().to_numpy(float)
    n = len(v)
    k = int(round(rate * n))
    if k < 1 or k >= n:
        return np.nan
    order = np.argsort(-v, kind="stable")
    x = np.zeros(n)
    x[order[:k]] = 1.0
    if x.std() == 0 or v.std() == 0:
        return np.nan
    return float(abs(np.corrcoef(v, x)[0, 1]))


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


CEIL = {}


def departure(frame, items, expos, adjust, normalise):
    covs = [frame[c] for c in adjust if c in frame.columns]
    B, ns = [], []
    for e in expos:
        rate = float(frame[e].mean())
        row = []
        for it in items:
            b, n = beta(frame[it], frame[e], covs)
            if normalise:
                c = attainable_ceiling(frame[it], rate)
                CEIL[f"{e}|{it}"] = c
                b = b / c if (c and not np.isnan(c) and c > 0) else np.nan
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
    frame = df.dropna(subset=NORMS + expos)
    for norm_on in (False, True):
        for aname, aset in ADJ_SETS.items():
            for drop in [None] + NORMS:
                items = [i for i in NORMS if i != drop]
                r, B, n = departure(frame, items, expos, aset, norm_on)
                grid.append(dict(expo_def=ename, normalised=norm_on, adjust=aname,
                                 dropped=drop or "none", ratio=r, n=n, coefs=B))
            r, B, n = departure(frame, PLACEBO, expos, aset, norm_on)
            grid.append(dict(expo_def=ename, normalised=norm_on, adjust=aname,
                             dropped="PLACEBO(premarsx+teensex)", ratio=r, n=n, coefs=B))

print("\n=== the attainable ceilings (a property of the MARGINAL, not of any association) ===")
for k, v in sorted(CEIL.items()):
    print(f"  {k:22s} ceiling |r|max = {v:.4f}")

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    r = "  nan " if g["ratio"] is None or np.isnan(g["ratio"]) else f"{g['ratio']:.4f}"
    print(f"  {g['expo_def']:7s} norm={str(g['normalised']):5s} {g['adjust']:12s} "
          f"drop={g['dropped']:26s} sigma2/sigma1={r}  n={g['n']:6d}")


def summarise(norm_on):
    full = [g["ratio"] for g in grid if g["dropped"] == "none" and g["normalised"] == norm_on
            and not np.isnan(g["ratio"])]
    med = float(np.median(full)) if full else np.nan
    d = {}
    for it in NORMS:
        v = [g["ratio"] for g in grid if g["dropped"] == it and g["normalised"] == norm_on
             and not np.isnan(g["ratio"])]
        d[it] = float(np.median(v)) if v else np.nan
    pl = [g["ratio"] for g in grid if g["dropped"].startswith("PLACEBO") and g["normalised"] == norm_on
          and not np.isnan(g["ratio"])]
    return med, d, (float(np.median(pl)) if pl else np.nan)


raw_full, raw_drops, raw_plac = summarise(False)
nrm_full, nrm_drops, nrm_plac = summarise(True)

# ══ NEGATIVE + POSITIVE controls on a one-factor synthetic ═══════════════════════════
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
            # ⚠ the synthetic reproduces the OBSERVED lumpy marginals, so the null lives in the
            #   same ceiling regime the data does — a null off that regime would not be its null.
            q = np.cumsum([marg[it]["shares"].get(str(i), 0.0) for i in (1, 2, 3)])
            cols[it] = pd.Series(np.clip(np.digitize(z, np.quantile(z, np.clip(q, 0, 1))) + 1, 1, 4))
        for j, e in enumerate(("ss", "xm")):
            rate = float(core[e].mean())
            z = 0.5 * th + g * dom[j] + RNG.standard_normal(n)
            cols[e] = pd.Series((z > np.quantile(z, 1 - rate)).astype(float))
        out.append(pd.DataFrame(cols))
    return out


def null_for(norm_on, reps=50):
    v = [x for x, _, _ in (departure(d, NORMS, ["ss", "xm"], [], norm_on) for d in synth(n_syn, 0.0, reps))
         if not np.isnan(x)]
    return float(np.median(v)), float(np.std(v))


null_raw = null_for(False)
null_nrm = null_for(True)
sweep = []
for g in (0.0, 0.15, 0.3, 0.45, 0.6):
    v = [x for x, _, _ in (departure(d, NORMS, ["ss", "xm"], [], True) for d in synth(n_syn, g, 18))
         if not np.isnan(x)]
    sweep.append((g, float(np.median(v)) if v else np.nan))

print(f"\n  RAW          full {raw_full:.4f} · null {null_raw[0]:.4f}+/-{null_raw[1]:.4f} · placebo {raw_plac:.4f}")
for it in NORMS:
    print(f"     drop {it:10s} {raw_drops[it]:.4f}  ({raw_drops[it]/raw_full:.2f}x)")
print(f"  CEILING-NORM full {nrm_full:.4f} · null {null_nrm[0]:.4f}+/-{null_nrm[1]:.4f} · placebo {nrm_plac:.4f}")
for it in NORMS:
    print(f"     drop {it:10s} {nrm_drops[it]:.4f}  ({nrm_drops[it]/nrm_full:.2f}x)")
print(f"  positive sweep (ceiling-normalised): {[(g, round(v, 4)) for g, v in sweep]}")


def carriers_of(full, drops, null_med):
    return [it for it in NORMS if not np.isnan(drops[it])
            and (drops[it] - null_med) < 0.5 * (full - null_med)]


car_raw = carriers_of(raw_full, raw_drops, null_raw[0])
car_nrm = carriers_of(nrm_full, nrm_drops, null_nrm[0])
print(f"\n  carriers RAW              : {car_raw}")
print(f"  carriers CEILING-NORMALISED: {car_nrm}")

ps = [2 * (1 - stats.norm.cdf(abs((g["ratio"] - (null_nrm[0] if g["normalised"] else null_raw[0]))
                                  / ((null_nrm[1] if g["normalised"] else null_raw[1]) or 1e-9))))
      for g in grid if not np.isnan(g["ratio"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

G = Gate("Did `xmarsex` carry none of the departure, or did it have no room to?")
G.plant_direction_from_sweep("positive: planted item coupling raises the normalised departure, g=0 null",
                             sweep, baseline=null_nrm[0], baseline_spread=null_nrm[1])
G.negative_control("one-factor synthetic reproducing the OBSERVED lumpy marginals",
                   abs(null_nrm[0]), abs(nrm_full), null_spread=null_nrm[1],
                   null_kind="one-factor latent at the observed marginals and base rates")
G.multiplicity_control("the whole raw+normalised drop grid", ps, 0.05,
                       labels=[f"{g['expo_def']}|{g['normalised']}|{g['adjust']}|{g['dropped']}"
                               for g in grid if not np.isnan(g["ratio"])])
G.asserted("placebo actually RAN in both scalings (absence is not a pass)",
           not np.isnan(raw_plac) and not np.isnan(nrm_plac),
           f"placebo raw {raw_plac:.4f} · normalised {nrm_plac:.4f}", kind="control")
G.asserted("the ceilings were MEASURED and they differ across items (else this round is vacuous)",
           bool(CEIL) and (max(CEIL.values()) - min(CEIL.values())) > 0.05,
           f"attainable |r|max spread {min(CEIL.values()):.4f}..{max(CEIL.values()):.4f} over "
           f"{len(CEIL)} (exposure,item) cells", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("KILL: W2 (ceiling artifact) requires the carrier set to CHANGE under normalisation",
           car_nrm == car_raw,
           f"carriers raw {car_raw} vs ceiling-normalised {car_nrm}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif car_nrm == car_raw == ["homosex"]:
    VERDICT, WORLD = "CONFIRMED", "W1 · REAL — `homosex` carries it at matched attainable range"
elif car_nrm != car_raw:
    VERDICT, WORLD = "OVERTURNED", f"W2 · CEILING ARTIFACT — carriers change to {car_nrm}"
else:
    VERDICT, WORLD = "UNVERIFIED", f"carriers stable but not homosex-only: {car_nrm}"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=918, round="E03·A112·R356", verdict=VERDICT, world=WORLD,
           estimand="sigma2/sigma1 of the 2xK coefficient matrix with each cell divided by the "
                    "attainable |r| ceiling implied by that item's marginal and the exposure's rate",
           instrument="GSS 1972-2024 gss7224_r3a.dta", marginals=marg, ceilings=CEIL,
           raw=dict(full=raw_full, drops=raw_drops, placebo=raw_plac,
                    null_median=null_raw[0], null_sd=null_raw[1], carriers=car_raw),
           normalised=dict(full=nrm_full, drops=nrm_drops, placebo=nrm_plac,
                           null_median=null_nrm[0], null_sd=null_nrm[1], carriers=car_nrm),
           positive_sweep=sweep, grid=grid, family_size=len(ps),
           disclosure_world="UNTESTED — `#917`(2)'s design died on the gradient check (both worlds "
                            "predict attenuation across 1991-2014); no replacement on this release",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "did_it_have_room_to_move.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'did_it_have_room_to_move.json'}")
