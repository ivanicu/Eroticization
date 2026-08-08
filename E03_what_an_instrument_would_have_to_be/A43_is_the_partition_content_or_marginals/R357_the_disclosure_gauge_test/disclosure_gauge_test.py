#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R357 — the disclosure rival, which `#918` said needed an acquisition, tested here
==========================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#915`-`#918` built to a single surviving claim — the coupling between doing and
                approving is specific to `homosex` — with exactly ONE live rival: **norms and
                behaviours are self-reported to the SAME interviewer**, so the coupling may be
                willingness-to-disclose rather than anything moral. `#918`(2) declared this needed an
                acquisition: "a behaviour measure not routed through the same interview".

⚠ WHY THAT     **`#918`(2) IS WRONG, and in the way this project has been burned by twice.** `#912`
CLAIM WAS       declared five cells UNREADABLE from a directory listing and `#913` opened four of
CHECKED         them; `#913`(3) named the class: *scoring a thing by its description when the thing
                is one command away.* So the impossibility claim was checked instead of inherited.
                **GSS ships `mode`** — `in-person` 21,436 · `web` 6,916 · `by phone` 3,129 ·
                `multimode` 485, over 11 waves 2004-2024 — and it overlaps this design at n=9,968.
                **`web` is self-administered: no interviewer.** That is a manipulation of disclosure
                pressure sitting inside the release.
                ⚠ *A fabricated impossibility beats a fabricated finding, because a wall makes
                stopping feel earned and is therefore never audited.*

The gauge test  `frontier` §1.3: name a transformation that should leave behaviour identical, then
                ask whether the MEASUREMENT is invariant under it and whether the PROPERTY is.
                **Mode is that transformation.** If the coupling is moral it is invariant to who is
                in the room; if it is disclosure it is not.

Live Worlds     W1 · **REAL** — the ceiling-normalised `ss -> homosex` coupling is invariant to mode.
                W2 · **DISCLOSURE** — it is LARGER where an interviewer is present (in-person) than
                     where none is (web), because the same reticence suppresses the behaviour report
                     and stiffens the norm report together.
                     **This is the unwelcome one: it would deflate the whole `#915`-`#918` chain.**
                W3 · **COMPOSITION** — mode is not randomised, so mode-groups differ in who they
                     are. Detected by the PLACEBO items moving too. (the meta-separator: the
                     {real, disclosure} split presumes mode changes only ANSWERING, not SAMPLING)

Estimand        The coefficient `b(ss -> homosex)` divided by its own **attainable ceiling** (`#918`:
(G1)            a property of `(marginal, base rate)` alone, containing no association), computed
                **WITHIN WAVE and BY MODE**, then contrasted across modes.
                ⚠ **WITHIN WAVE IS NOT OPTIONAL**: `web` is 2021+ almost entirely, so a naive
                mode contrast is an ERA contrast, and `#918` already established that a real
                coupling attenuates across these decades too. Only the within-wave contrast
                separates anything. Usable contrasts, measured before designing:
                  in-person vs phone  — 8 waves 2004-2018, in-person n 539..1387, phone n 46..182
                  in-person vs web    — 2022 only, n 408 vs 484

Prediction      W1 -> contrast ~ 0 in both families.
Matrix          W2 -> contrast > 0 (interviewer-present larger), in BOTH families.
                W3 -> the placebo items show the same contrast as `homosex`.

Strongest       ⚠ Mode is **not randomised**. CONTROL, same iteration: (a) demographic adjustment;
Confound        (b) the PLACEBO items `premarsx`/`teensex`, which the exposure is not tied to — under
                disclosure their coupling should NOT move with mode, under composition it should;
                (c) the norm's own marginal by mode is reported so a level shift cannot hide.

Stopping Rule   One pass over contrast-family x adjustment x wave, published whole. Under-powered
                cells are reported with their MDE, never dropped.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **mode is not randomised** — this is an observational gauge test, so W3 can only be
    DETECTED by the placebo, never excluded by design;
  (2) ⚠ the interviewer-absent arm is **one wave** (2022, n 408 vs 484); the 8-wave family contrasts
    in-person against phone, and **phone still has an interviewer** — it is a weaker manipulation,
    so a null there is a weaker acquittal and is reported as such;
  (3) `web` respondents in 2022 are a mode EXPERIMENT arm, not a random subsample of the frame;
  (4) ⚠ **only this one instrument**, forced by the question: the object under test is whether GSS's
    OWN `homosex` coupling survives GSS's OWN mode variation. A second release has different modes
    on different items and would be a new claim, not a replication — the cross-instrument move
    already happened at `#917` (NSFG -> GSS).
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
RNG = np.random.default_rng(357)

TARGET = "homosex"
PLACEBO = ["premarsx", "teensex"]
MODE = {1: "in-person", 2: "phone", 3: "multimode", 4: "web"}

df = pd.read_stata(GSS, columns=["mode", TARGET] + PLACEBO + ["nummen", "numwomen", "sex", "year",
                                                             "age", "educ", "race", "attend"],
                   convert_categoricals=False)
df["modename"] = df["mode"].map(MODE)
for c in [TARGET] + PLACEBO:
    df[c] = df[c].where(df[c].isin([1, 2, 3, 4]))
ok = df["nummen"].notna() & df["numwomen"].notna() & df["sex"].isin([1, 2])
df["ss"] = np.where(((df["sex"] == 1) & (df["nummen"] > 0)) | ((df["sex"] == 2) & (df["numwomen"] > 0)), 1.0, 0.0)
df.loc[~ok, "ss"] = np.nan

# ══ HARD RULE 1 — the mode variable this whole round rests on ════════════════════════
base = df.dropna(subset=["modename", TARGET, "ss"])
print(f"mode x {TARGET} x ss overlap: n={len(base)}")
inv = {}
for m, g in base.groupby("modename"):
    inv[m] = dict(n=int(len(g)), ss_rate=float(g["ss"].mean()), norm_mean=float(g[TARGET].mean()),
                  waves=[int(y) for y in sorted(g["year"].unique())])
    print(f"  {m:11s} n={len(g):5d}  ss_rate={g['ss'].mean():.4f}  {TARGET}_mean={g[TARGET].mean():.3f}"
          f"  waves {int(g['year'].min())}..{int(g['year'].max())}")


def attainable_ceiling(y, rate):
    """`#918`: max |corr| a binary predictor at this base rate can reach against this marginal.
    A property of (marginal, base rate) ALONE — it contains no association."""
    v = y.dropna().to_numpy(float)
    n = len(v)
    k = int(round(rate * n))
    if k < 1 or k >= n:
        return np.nan
    x = np.zeros(n)
    x[np.argsort(-v, kind="stable")[:k]] = 1.0
    if x.std() == 0 or v.std() == 0:
        return np.nan
    return float(abs(np.corrcoef(v, x)[0, 1]))


def norm_beta(frame, item, adjust):
    """Ceiling-normalised standardised slope of `item` on `ss`, plus n and the raw pieces."""
    covs = [frame[c] for c in adjust if c in frame.columns]
    y, x = frame[item], frame["ss"]
    m = y.notna() & x.notna()
    for c in covs:
        m &= c.notna()
    n = int(m.sum())
    if n < 60 or x[m].std() == 0:
        return np.nan, np.nan, n
    Y = y[m].to_numpy(float)
    X = np.column_stack([np.ones(n), x[m].to_numpy(float)] + [c[m].to_numpy(float) for c in covs])
    for j in range(1, X.shape[1]):
        sd = X[:, j].std()
        if sd > 0:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    Y = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n
    raw = float(b[1])
    ceil = attainable_ceiling(y[m], float(x[m].mean()))
    return (raw / ceil if (ceil and not np.isnan(ceil) and ceil > 0) else np.nan), raw, n


ADJ = {"raw": [], "demog": ["age", "educ", "race"], "demog+relig": ["age", "educ", "race", "attend"]}
FAMILIES = {"in-person vs phone (8 waves, BOTH have an interviewer)": ("in-person", "phone"),
            "in-person vs web (2022, web has NONE)": ("in-person", "web")}

grid = []
for fam, (a, b_) in FAMILIES.items():
    waves = sorted(set(base.loc[base.modename == a, "year"]) & set(base.loc[base.modename == b_, "year"]))
    for y in waves:
        w = base[base.year == y]
        ga, gb = w[w.modename == a], w[w.modename == b_]
        if len(ga) < 60 or len(gb) < 60:
            continue
        for aname, aset in ADJ.items():
            for item in [TARGET] + PLACEBO:
                na, ra, n_a = norm_beta(ga, item, aset)
                nb, rb, n_b = norm_beta(gb, item, aset)
                grid.append(dict(family=fam, year=int(y), adjust=aname, item=item,
                                 interviewer_present=na, interviewer_absent=nb,
                                 contrast=(na - nb) if not (np.isnan(na) or np.isnan(nb)) else np.nan,
                                 n=int(min(n_a, n_b)), n_a=n_a, n_b=n_b))

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    c = "  nan " if np.isnan(g["contrast"]) else f"{g['contrast']:+.4f}"
    print(f"  {g['family'][:34]:34s} {g['year']} {g['adjust']:12s} {g['item']:9s} "
          f"present={g['interviewer_present']:+.4f} absent={g['interviewer_absent']:+.4f} "
          f"contrast={c}  n={g['n']:5d}")


def summarise(fam, item):
    v = [g["contrast"] for g in grid if g["family"] == fam and g["item"] == item and not np.isnan(g["contrast"])]
    return (float(np.median(v)) if v else np.nan, len(v),
            float(np.std(v)) if len(v) > 1 else np.nan)


print("\n=== the contrast, by family (interviewer-PRESENT minus interviewer-ABSENT) ===")
summary = {}
for fam in FAMILIES:
    for item in [TARGET] + PLACEBO:
        med, k, sd = summarise(fam, item)
        summary[f"{fam}|{item}"] = dict(median=med, cells=k, sd=sd)
        tag = "  <- the carrier" if item == TARGET else "  (placebo)"
        print(f"  {fam[:34]:34s} {item:9s} median {med:+.4f} over {k:2d} cells (sd {sd:.4f}){tag}")

# ══ CONTROLS ═════════════════════════════════════════════════════════════════════════
# NEGATIVE: split each mode arm at random (no disclosure difference exists) -> contrast must be null
null_vals = []
for fam, (a, b_) in FAMILIES.items():
    waves = sorted(set(base.loc[base.modename == a, "year"]) & set(base.loc[base.modename == b_, "year"]))
    for y in waves:
        w = base[(base.year == y) & (base.modename.isin([a, b_]))]
        if len(w) < 200:
            continue
        for _ in range(40):
            perm = w.copy()
            perm["modename"] = RNG.permutation(perm["modename"].to_numpy())
            na, _, _ = norm_beta(perm[perm.modename == a], TARGET, [])
            nb, _, _ = norm_beta(perm[perm.modename == b_], TARGET, [])
            if not (np.isnan(na) or np.isnan(nb)):
                null_vals.append(na - nb)
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  label-permutation null (mode labels shuffled WITHIN wave, sizes preserved): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# POSITIVE: plant a disclosure effect of size g into the interviewer-present arm; must be null at g=0
sweep = []
for g in (0.0, 0.10, 0.20, 0.30, 0.40):
    vals = []
    for fam, (a, b_) in FAMILIES.items():
        waves = sorted(set(base.loc[base.modename == a, "year"]) & set(base.loc[base.modename == b_, "year"]))
        for y in waves[:4]:
            w = base[base.year == y]
            ga, gb = w[w.modename == a].copy(), w[w.modename == b_]
            if len(ga) < 60 or len(gb) < 60:
                continue
            # disclosure = reticent people BOTH hide the behaviour and stiffen the norm
            u = RNG.standard_normal(len(ga))
            ga[TARGET] = np.clip(np.round(ga[TARGET] - g * (u > 0) * 1.0), 1, 4)
            ga["ss"] = np.where((ga["ss"] == 1) & (u > 0) & (RNG.random(len(ga)) < g), 0.0, ga["ss"])
            na, _, _ = norm_beta(ga, TARGET, [])
            nb, _, _ = norm_beta(gb, TARGET, [])
            if not (np.isnan(na) or np.isnan(nb)):
                vals.append(na - nb)
    sweep.append((g, float(np.median(vals)) if vals else np.nan))
print(f"  positive sweep (planted disclosure g, median contrast): {[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((g["contrast"] - null_med) / (null_sd or 1e-9))))
      for g in grid if not np.isnan(g["contrast"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

tgt_web = summary[f"in-person vs web (2022, web has NONE)|{TARGET}"]["median"]
tgt_ph = summary[f"in-person vs phone (8 waves, BOTH have an interviewer)|{TARGET}"]["median"]
plac_meds = [summary[f"{f}|{p}"]["median"] for f in FAMILIES for p in PLACEBO
             if not np.isnan(summary[f"{f}|{p}"]["median"])]
plac_max = float(np.max(np.abs(plac_meds))) if plac_meds else np.nan

G = Gate("Is the `homosex` coupling invariant to who is in the room?")
G.plant_direction_from_sweep("positive: planted disclosure raises the contrast, and g=0 is null",
                             sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("mode labels permuted WITHIN wave, arm sizes preserved",
                   abs(null_med), abs(tgt_web) if not np.isnan(tgt_web) else 0.0,
                   null_spread=null_sd, null_kind="within-wave mode-label permutation")
G.multiplicity_control("the whole family x wave x adjustment x item grid", ps, 0.05,
                       labels=[f"{g['family'][:16]}|{g['year']}|{g['adjust']}|{g['item']}"
                               for g in grid if not np.isnan(g["contrast"])])
G.asserted("placebo actually RAN in both families (absence is not a pass)",
           len(plac_meds) >= 2, f"{len(plac_meds)} placebo family-medians computed", kind="control")
G.asserted("composition detector: the PLACEBO items must not move with mode as much as the carrier",
           not np.isnan(plac_max) and not np.isnan(tgt_web) and plac_max < abs(tgt_web) + 2 * null_sd,
           f"largest |placebo contrast| {plac_max:.4f} vs carrier(web family) {tgt_web:+.4f}",
           kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.resolvable("the interviewer-absent contrast is above the permutation floor",
             abs(tgt_web - null_med) if not np.isnan(tgt_web) else 0.0, null_sd)
G.asserted("KILL: W2 (disclosure) requires the contrast to be POSITIVE where the interviewer is absent",
           np.isnan(tgt_web) or abs(tgt_web - null_med) < 2 * null_sd,
           f"in-person minus web contrast {tgt_web:+.4f} vs null {null_med:+.4f} +/- {null_sd:.4f}; "
           f"in-person minus phone {tgt_ph:+.4f}")

tv = G.three_valued()
resolved_web = not np.isnan(tgt_web) and abs(tgt_web - null_med) > 2 * null_sd
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif not np.isnan(plac_max) and not np.isnan(tgt_web) and plac_max >= abs(tgt_web):
    VERDICT, WORLD = "UNVERIFIED", "W3 · COMPOSITION — the placebo moves as much as the carrier"
elif resolved_web and tgt_web > 0:
    VERDICT, WORLD = "OVERTURNED", "W2 · DISCLOSURE — the coupling depends on who is in the room"
elif not resolved_web:
    VERDICT, WORLD = "CONFIRMED", "W1 · REAL — invariant to mode at this design's resolution"
else:
    VERDICT, WORLD = "UNVERIFIED", f"contrast resolved but negative ({tgt_web:+.4f})"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=919, round="E03·A112·R357", verdict=VERDICT, world=WORLD,
           estimand="ceiling-normalised b(ss -> homosex) WITHIN WAVE, contrasted across interview "
                    "mode; mode is the gauge transformation",
           instrument="GSS 1972-2024 gss7224_r3a.dta, `mode` variable",
           mode_inventory=inv, grid=grid, summary=summary,
           carrier_contrast_web=tgt_web, carrier_contrast_phone=tgt_ph, placebo_max=plac_max,
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals), positive_sweep=sweep,
           corrects="`#918`(2) said this needed an acquisition; GSS ships `mode`",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "disclosure_gauge_test.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'disclosure_gauge_test.json'}")
