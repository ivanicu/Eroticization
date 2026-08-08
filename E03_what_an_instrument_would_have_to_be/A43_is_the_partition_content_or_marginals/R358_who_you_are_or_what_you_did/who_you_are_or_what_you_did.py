#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R358 — my own page sentence, tested: who you ARE, or what you DID?
===========================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#918` ends with a sentence now on the front page: *"People argue their own case
                about who they are, and do not argue it about what they have done."* ⚠ **That was
                never measured.** `#915`-`#918` measured a coupling between a BEHAVIOUR report and a
                NORM; the word "who they are" was supplied in the same breath as the pattern — which
                is precisely the failure `#911`(2) registered and I keep committing.

Why Now         GSS carries **both**: `sexornt` (self-labelled identity, 2008+) and `ss` (own
                same-sex partner). Their DISCORDANT cells are the discriminating cases and they are
                not small: **497 people report a same-sex partner while identifying heterosexual**,
                and 77 identify gay/lesbian/bisexual with no reported same-sex partner.
                ⚠ `#111c` closes the disclosure question after `#919`; this is a different question.

Live Worlds     W_ID   · identity survives conditioning on behaviour; behaviour does not survive
                         conditioning on identity. **My sentence stands.**
                W_ACT  · the reverse. **Unwelcome: it retracts the page sentence.**
                W_BOTH · both survive net of each other ⇒ the identity/act split is a FALSE
                         DICHOTOMY and the sentence was a rhetorical choice, not a finding.
                         (the meta-separator: an outcome killing my world-DECOMPOSITION)
                W_IDEO · identity's effect on `homosex` is matched by its effect on a NON-SEXUAL
                         moral item and on general ideology ⇒ this is liberalism, not a norm about
                         the act. **Unwelcome, and it threatens more than the sentence.**

Estimand        The contribution of identity to `homosex` NET of behaviour, and of behaviour NET of
(G1)            identity, each as a standardised slope divided by its own **attainable ceiling**
                (`#918`: a property of `(marginal, base rate)` alone, containing no association);
                plus the same two quantities computed against **`cappun`** (death penalty — moral,
                non-sexual, n=7,427) and **`polviews`** (n=7,547) as the specificity detectors.

Prediction      W_ID   -> id_net resolved, act_net not.
Matrix          W_ACT  -> act_net resolved, id_net not.
                W_BOTH -> both resolved on `homosex`.
                W_IDEO -> id_net on `cappun`/`polviews` comparable to id_net on `homosex`.

Strongest       ⚠ **GENERAL LIBERALISM.** People who identify as gay, lesbian or bisexual are more
Confound        liberal across the board, so identity may predict `homosex` through ideology and
                nothing sexual. CONTROL, same iteration: `cappun` and `polviews` are carried through
                the whole grid as outcomes, and `polviews` also enters the adjustment axis.

⚠ Inherited     `#919` established that the DISCLOSURE rival cannot be resolved on this release
limits          (positive control non-monotone, floor +/-0.3095 against a -0.1708 contrast, placebo
                moving more than the carrier). **Identity is also a report to the same interviewer.**
                So this round CANNOT separate "identity matters" from "identity-reporting and
                norm-reporting share a disclosure component". Registered, not re-litigated.

Stopping Rule   One pass over outcome x adjustment x identity-coding x era, published whole, with an
                MDE for every cell — the LGB x no-behaviour cells are n=10 and n=67 and a null there
                is silence unless its MDE is printed.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ disclosure, inherited from `#919` and NOT addressed here;
  (2) ⚠ the arrow is not identified — identity, behaviour and norm are simultaneous, and a person
    who accepts the norm may be likelier to adopt the label. This is decomposition, not causation;
  (3) ⚠ `sexornt` starts 2008, so this is **8 waves**, not the 1972-2024 series;
  (4) ⚠ **only this one instrument**, forced by the question: no other release here carries a
    self-labelled orientation AND a behaviour count AND the norm. The cross-instrument move already
    happened at `#917` (NSFG -> GSS);
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
RNG = np.random.default_rng(358)

OUTCOMES = {"homosex": (1, 4), "cappun": (1, 2), "polviews": (1, 7)}
COLS = ["sexornt", "homosex", "cappun", "polviews", "nummen", "numwomen", "sex", "year",
        "age", "educ", "race", "attend"]

df = pd.read_stata(GSS, columns=COLS, convert_categoricals=False)
df["homosex"] = df["homosex"].where(df["homosex"].isin([1, 2, 3, 4]))
df["cappun"] = df["cappun"].where(df["cappun"].isin([1, 2]))
df["polviews"] = df["polviews"].where(df["polviews"].between(1, 7))
ok = df["nummen"].notna() & df["numwomen"].notna() & df["sex"].isin([1, 2])
df["ss"] = np.where(((df["sex"] == 1) & (df["nummen"] > 0)) | ((df["sex"] == 2) & (df["numwomen"] > 0)), 1.0, 0.0)
df.loc[~ok, "ss"] = np.nan
df["lgb"] = np.where(df["sexornt"].isin([1, 2]), 1.0, np.where(df["sexornt"] == 3, 0.0, np.nan))

base = df.dropna(subset=["lgb", "ss", "homosex"])
print(f"identity x behaviour x norm: n={len(base)}  waves {[int(y) for y in sorted(base.year.unique())]}")
cells = base.groupby(["lgb", "ss"]).agg(n=("homosex", "size"), homosex=("homosex", "mean")).round(4)
print(cells)
inv = {f"lgb={int(k[0])},ss={int(k[1])}": dict(n=int(v["n"]), homosex_mean=float(v["homosex"]))
       for k, v in cells.iterrows()}
for c in OUTCOMES:
    print(f"  outcome {c:9s} non-missing within frame: {int(base[c].notna().sum())}")


def ceiling(y, rate):
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


def net_slopes(frame, item, adjust):
    """Standardised slopes of `lgb` and `ss` on `item`, EACH NET OF THE OTHER, ceiling-normalised.
    Returns (id_net, act_net, n)."""
    covs = [frame[c] for c in adjust if c in frame.columns]
    y = frame[item]
    m = y.notna() & frame["lgb"].notna() & frame["ss"].notna()
    for c in covs:
        m &= c.notna()
    n = int(m.sum())
    if n < 200:
        return np.nan, np.nan, n
    Y = y[m].to_numpy(float)
    cols = [frame["lgb"][m].to_numpy(float), frame["ss"][m].to_numpy(float)] + \
           [c[m].to_numpy(float) for c in covs]
    X = np.column_stack([np.ones(n)] + cols)
    for j in range(1, X.shape[1]):
        sd = X[:, j].std()
        if sd > 0:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    Yz = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        b, *_ = np.linalg.lstsq(X, Yz, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n
    c_id = ceiling(y[m], float(frame["lgb"][m].mean()))
    c_act = ceiling(y[m], float(frame["ss"][m].mean()))
    idn = float(b[1]) / c_id if (c_id and not np.isnan(c_id) and c_id > 0) else np.nan
    actn = float(b[2]) / c_act if (c_act and not np.isnan(c_act) and c_act > 0) else np.nan
    return idn, actn, n


ADJ = {"raw": [], "demog": ["age", "educ", "race"], "demog+relig": ["age", "educ", "race", "attend"],
       "demog+relig+ideology": ["age", "educ", "race", "attend", "polviews"]}
ERAS = {"all 2008-2022": None, "2008-2014": (2008, 2014), "2016-2022": (2016, 2022)}

grid = []
for item in OUTCOMES:
    for aname, aset in ADJ.items():
        if item == "polviews" and "polviews" in aset:
            continue                                    # never adjust an outcome for itself
        for ename, span in ERAS.items():
            fr = base if span is None else base[base.year.between(*span)]
            idn, actn, n = net_slopes(fr, item, aset)
            grid.append(dict(outcome=item, adjust=aname, era=ename, id_net=idn, act_net=actn, n=n))

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    f = lambda v: "  nan " if np.isnan(v) else f"{v:+.4f}"          # noqa: E731
    print(f"  {g['outcome']:9s} {g['adjust']:21s} {g['era']:13s} "
          f"id_net={f(g['id_net'])}  act_net={f(g['act_net'])}  n={g['n']:5d}")

# ══ NEGATIVE CONTROLS — and there must be TWO, one per quantity ══════════════════════
# ⚠⚠ v1 USED ONE NULL FOR BOTH, AND IT IS WRONG FOR act_net BY CONSTRUCTION. Permuting identity
#   within behaviour destroys identity's contribution and therefore hands the shared variance to
#   BEHAVIOUR: measured, that null put act_net at +0.4067 while the OBSERVED act_net is +0.1539, so
#   "act_net differs from its null" fired because the observed value is BELOW an inflated null.
#   **A null built by destroying the OTHER predictor is not this predictor's null.**
#   ⇒ each quantity gets the permutation that destroys ITSELF, holding the other fixed.
def _perm_null(which, reps=200):
    vals = []
    other = "ss" if which == "lgb" else "lgb"
    for _ in range(reps):
        p = base.copy()
        p[which] = p.groupby(other)[which].transform(lambda s: RNG.permutation(s.to_numpy()))
        i, a, _ = net_slopes(p, "homosex", [])
        v = i if which == "lgb" else a
        if not np.isnan(v):
            vals.append(v)
    return float(np.median(vals)), float(np.std(vals)), len(vals)


nid_m, nid_s, nid_k = _perm_null("lgb")
nact_m, nact_s, nact_k = _perm_null("ss")
print(f"\n  null for id_net  (identity permuted within behaviour): {nid_m:+.4f} +/- {nid_s:.4f} ({nid_k})")
print(f"  null for act_net (behaviour permuted within identity): {nact_m:+.4f} +/- {nact_s:.4f} ({nact_k})")

# ══ POSITIVE CONTROL — plant into a NULL WORLD, so g=0 genuinely has no effect ═══════
# ⚠⚠ v1 PLANTED INTO THE OBSERVED DATA and judged the sweep against the PERMUTATION baseline —
#   two different worlds, which is the `#905` failure, committed again. At g=0 the observed data
#   already carries id_net=+0.6894 while the baseline is -0.0015, so the control could only fail.
#   ⇒ plant into the permuted (identity-destroyed) world: g=0 then lands ON the baseline.
sweep_id, sweep_act = [], []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vi, va = [], []
    for _ in range(12):
        p = base.copy()
        p["lgb"] = p.groupby("ss")["lgb"].transform(lambda s: RNG.permutation(s.to_numpy()))
        p["homosex"] = np.clip(p["homosex"] + g * p["lgb"] * 1.0, 1, 4)   # identity ONLY
        i, a, _ = net_slopes(p, "homosex", [])
        if not np.isnan(i):
            vi.append(i)
        if not np.isnan(a):
            va.append(a)
    sweep_id.append((g, float(np.median(vi)) if vi else np.nan))
    sweep_act.append((g, float(np.median(va)) if va else np.nan))
print(f"  positive sweep (planted into the NULL world), IDENTITY planted:")
print(f"     id_net  {[(g, round(v, 4)) for g, v in sweep_id]}")
print(f"     act_net {[(g, round(v, 4)) for g, v in sweep_act]}  (must NOT rise)")

# ══ MDE for the small cells — a MEASURED bootstrap, not a rule of thumb ══════════════
# ⚠ v1 passed `1/sqrt(n)` as the spread source and the gate rejected it: not in the registry of
#   admissible sources. It was an analytic guess wearing a measurement's clothes.
n_small = int(min(v["n"] for k, v in inv.items() if k.startswith("lgb=1")))
lgb_cells = base[base.lgb == 1]
boot = []
for _ in range(400):
    s = lgb_cells.sample(len(lgb_cells), replace=True, random_state=int(RNG.integers(1e9)))
    a = s[s.ss == 1]["homosex"]
    b = s[s.ss == 0]["homosex"]
    if len(a) > 5 and len(b) > 5 and s["homosex"].std() > 0:
        boot.append((a.mean() - b.mean()) / s["homosex"].std())
mde = float(2 * np.std(boot)) if boot else np.nan
print(f"  smallest LGB cell n={n_small} -> bootstrap MDE on the standardised within-identity "
      f"behaviour gap ~{mde:.3f} ({len(boot)} resamples)")


def med(item, key, adj=None):
    v = [g[key] for g in grid if g["outcome"] == item and not np.isnan(g[key])
         and (adj is None or g["adjust"] == adj)]
    return float(np.median(v)) if v else np.nan


id_hom, act_hom = med("homosex", "id_net"), med("homosex", "act_net")
id_cap, act_cap = med("cappun", "id_net"), med("cappun", "act_net")
id_pol = med("polviews", "id_net")
id_hom_ideo = med("homosex", "id_net", "demog+relig+ideology")
act_hom_ideo = med("homosex", "act_net", "demog+relig+ideology")
print(f"\n  homosex : id_net {id_hom:+.4f}  act_net {act_hom:+.4f}")
print(f"  cappun  : id_net {id_cap:+.4f}  act_net {act_cap:+.4f}   <- non-sexual moral placebo")
print(f"  polviews: id_net {id_pol:+.4f}                            <- general ideology")
print(f"  homosex net of IDEOLOGY: id_net {id_hom_ideo:+.4f}  act_net {act_hom_ideo:+.4f}")

ps = [2 * (1 - stats.norm.cdf(abs((g["id_net"] - nid_m) / (nid_s or 1e-9)))) for g in grid
      if not np.isnan(g["id_net"])] + \
     [2 * (1 - stats.norm.cdf(abs((g["act_net"] - nact_m) / (nact_s or 1e-9)))) for g in grid
      if not np.isnan(g["act_net"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

id_res = abs(id_hom - nid_m) > 2 * nid_s
act_res = abs(act_hom - nact_m) > 2 * nact_s

G = Gate("Who you ARE, or what you DID?")
G.plant_direction_from_sweep("positive: an IDENTITY-only plant raises id_net, and g=0 is null",
                             sweep_id, baseline=nid_m, baseline_spread=nid_s)
G.asserted("the identity-only plant does NOT raise act_net (else the two are not separable)",
           abs(sweep_act[-1][1] - sweep_act[0][1]) < abs(sweep_id[-1][1] - sweep_id[0][1]),
           f"act_net moved {sweep_act[-1][1]-sweep_act[0][1]:+.4f} vs id_net "
           f"{sweep_id[-1][1]-sweep_id[0][1]:+.4f} under an identity-only plant", kind="control")
G.negative_control("identity permuted WITHIN behaviour strata", abs(nid_m), abs(id_hom),
                   null_spread=nid_s, null_kind="within-behaviour identity-label permutation")
G.multiplicity_control("the whole outcome x adjustment x era grid", ps, 0.05)
G.asserted("specificity: identity must predict `homosex` more than the NON-SEXUAL moral item",
           not np.isnan(id_cap) and abs(id_hom) > abs(id_cap),
           f"id_net homosex {id_hom:+.4f} vs cappun {id_cap:+.4f} (polviews {id_pol:+.4f})",
           kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.count_needs_interval("the LGB x no-behaviour cells are small and their MDE is stated",
                       n_small, n_small, mde, "bootstrap_人层", n_resamples=len(boot))
G.asserted("KILL: W_ID (my page sentence) requires identity to resolve and behaviour NOT to",
           not (id_res and not act_res),
           f"id_net resolved={id_res} ({id_hom:+.4f} vs {nid_m:+.4f}+/-{nid_s:.4f}) · "
           f"act_net resolved={act_res} ({act_hom:+.4f} vs {nact_m:+.4f}+/-{nact_s:.4f})")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif not np.isnan(id_cap) and abs(id_cap) >= abs(id_hom):
    VERDICT, WORLD = "OVERTURNED", "W_IDEO · general liberalism, not a norm about the act"
elif id_res and act_res:
    VERDICT, WORLD = "OVERTURNED", "W_BOTH · the identity/act split is a FALSE DICHOTOMY"
elif id_res:
    VERDICT, WORLD = "CONFIRMED", "W_ID · identity, net of behaviour"
elif act_res:
    VERDICT, WORLD = "OVERTURNED", "W_ACT · behaviour, net of identity"
else:
    VERDICT, WORLD = "UNVERIFIED", "neither resolves at this n"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=920, round="E03·A112·R358", verdict=VERDICT, world=WORLD,
           estimand="identity and behaviour contributions to `homosex`, each NET of the other, "
                    "ceiling-normalised, with `cappun`/`polviews` as specificity detectors",
           instrument="GSS 2008-2022 gss7224_r3a.dta", cells=inv, grid=grid,
           homosex=dict(id_net=id_hom, act_net=act_hom, id_resolved=bool(id_res),
                        act_resolved=bool(act_res)),
           cappun=dict(id_net=id_cap, act_net=act_cap), polviews=dict(id_net=id_pol),
           net_of_ideology=dict(id_net=id_hom_ideo, act_net=act_hom_ideo),
           null=dict(id_median=nid_m, id_sd=nid_s, id_draws=nid_k,
                     act_median=nact_m, act_sd=nact_s, act_draws=nact_k),
           positive_sweep_id=sweep_id, positive_sweep_act=sweep_act,
           smallest_lgb_cell=n_small, mde=mde, mde_resamples=len(boot), family_size=len(ps),
           inherited_limit="disclosure, from `#919`; identity is also a report to the interviewer",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "who_you_are_or_what_you_did.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'who_you_are_or_what_you_did.json'}")
