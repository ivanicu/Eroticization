#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R359 — `#920` on a second instrument, and identity under a RANDOMISED wording
======================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#920` decomposed the same-sex norm into identity and behaviour on GSS and found
                both resolved, identity ~2.9x behaviour. Two things were never done: it was never
                tried on a **second instrument** (HARD RULE 4), and its identity measure carries the
                shared-instrument threat (HARD RULE 2) that `#919` could not remove **because mode
                was not randomised and composition swamped it**.

Why Now         ⚠ **NSFG randomises the orientation question's wording.** `orient_a` and `orient_b`
                are disjoint half-samples (**13-14 respondents of ~6,000 answer both**) with
                DIFFERENT option sets — `orient_a` shows 3 substantive codes, `orient_b` shows 4.
                **Random assignment balances composition by design**, which is exactly the control
                `#919` could not build. So the instrument-dependence of IDENTITY is measurable here.

Live Worlds     W_STABLE  · `id_net` is the same under both randomised wordings, and the
                            identity/behaviour decomposition replicates on NSFG. `#920` travels.
                W_WORDING · `id_net` differs by randomised wording ⇒ identity's coupling with the
                            norm is partly a property of the QUESTION, not of the person.
                            **Unwelcome: it deflates `#920`'s headline.**
                W_NOREP   · the decomposition does not replicate on NSFG at all ⇒ `#920` was
                            GSS-specific. (also unwelcome, and a different kind of loss)
                ⚠ meta-separator: both wordings agreeing while NSFG's id/act RATIO differs sharply
                from GSS's would mean the decomposition is instrument-specific in a way none of the
                three worlds names.

Estimand        As `#920`: the contribution of identity to the same-sex norm NET of behaviour, and
(G1)            of behaviour NET of identity, each a standardised slope over its own attainable
                ceiling. Then the **randomised contrast** `id_net(orient_a) - id_net(orient_b)`.

⚠ THE CODING    **Value labels do not ship** (`#915`(2)), so which code means "gay/lesbian/bisexual"
IS ANCHORED,    is unknown from the file. It is NOT guessed: each identity code is scored by the
NOT GUESSED     share of its members reporting same-sex contact (`samesexany`), and the LGB codes are
                the ones whose members overwhelmingly do. ⚠ **That anchor uses the BEHAVIOUR
                variable and never the NORM**, so it cannot launder the finding it is used to test —
                the same logic as `#914`, which anchored a parse without reading a value.

Prediction      W_STABLE  -> |id_net(a) - id_net(b)| within its permutation null.
Matrix          W_WORDING -> outside it.
                W_NOREP   -> id_net on NSFG inside its own null in every specification.

Stopping Rule   One pass over file x wording-arm x adjustment, published whole. The A/B contrast has
                its own null by permuting the RANDOMISED label, which is the exact null because the
                assignment really is random.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ value labels do not ship — the coding is anchored on behaviour, which is evidence, not proof;
  (2) ⚠ the arrow is not identified (inherited from `#920`): identity, behaviour and norm are
    simultaneous;
  (3) ⚠ the randomisation balances WHO answers each wording; it does not make identity a
    non-self-report. Disclosure remains unresolved (inherited from `#919`);
  (4) the A/B split exists only in 2017-2019; 2011-2013 has a single `orient` and contributes to the
    replication but not to the wording contrast;
  (5) `[unchallenged]` — door (3).
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
D_ = ROOT / "data" / "external" / "nsfg"
RNG = np.random.default_rng(359)
PAT = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')

SITES = [("2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct", "2011-2013 female", ["orient"]),
         ("2017_2019_FemRespData.dat", "2017_2019_FemRespSetup.dct", "2017-2019 female", ["orient_a", "orient_b"]),
         ("2017_2019_MaleData.dat", "2017_2019_MaleSetup.dct", "2017-2019 male", ["orient_a", "orient_b"])]
NORM = "samesex"
COVS = ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]


def read(dat, dct, extra):
    lay = {m.group(2).lower(): (int(m.group(1)) - 1, int(m.group(3))) for m in
           (PAT.search(l) for l in (D_ / "setup" / dct).read_text(errors="replace").splitlines()) if m}
    want = [v for v in [NORM, "samesexany"] + extra + COVS if v in lay]
    df = pd.read_fwf(D_ / dat, colspecs=[(lay[v][0], lay[v][0] + lay[v][1]) for v in want],
                     names=want, dtype=str)
    for c in want:
        df[c] = pd.to_numeric(df[c].astype(str).str.strip().replace({"": None}), errors="coerce")
    df[NORM] = df[NORM].where(df[NORM].isin([1, 2, 3, 4, 5]))
    df["act"] = (df["samesexany"] == 1).astype(float).where(df["samesexany"].isin([1, 5]))
    return df


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


def net_slopes(fr, adjust):
    covs = [fr[c] for c in adjust if c in fr.columns]
    y = fr[NORM]
    m = y.notna() & fr["lgb"].notna() & fr["act"].notna()
    for c in covs:
        m &= c.notna()
    n = int(m.sum())
    if n < 200 or fr["lgb"][m].std() == 0 or fr["act"][m].std() == 0:
        return np.nan, np.nan, n
    Y = y[m].to_numpy(float)
    X = np.column_stack([np.ones(n), fr["lgb"][m].to_numpy(float), fr["act"][m].to_numpy(float)]
                        + [c[m].to_numpy(float) for c in covs])
    for j in range(1, X.shape[1]):
        sd = X[:, j].std()
        if sd > 0:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    Yz = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        b, *_ = np.linalg.lstsq(X, Yz, rcond=None)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n
    ci = ceiling(y[m], float(fr["lgb"][m].mean()))
    ca = ceiling(y[m], float(fr["act"][m].mean()))
    return (float(b[1]) / ci if (ci and ci > 0) else np.nan,
            float(b[2]) / ca if (ca and ca > 0) else np.nan, n)


# ══ HARD RULE 1 + the anchoring of the coding, on BEHAVIOUR and never on the norm ════
frames, anchor = {}, {}
for dat, dct, lab, cols in SITES:
    df = read(dat, dct, cols)
    for c in cols:
        sub = df[df[c].isin([1, 2, 3, 4])]
        rates = sub.groupby(c)["act"].agg(["mean", "size"])
        # LGB = every code whose members report same-sex contact far above the file's base rate
        basev = float(df["act"].mean())
        lgb_codes = [int(k) for k, r in rates.iterrows() if r["mean"] > max(3 * basev, 0.30)]
        anchor[f"{lab}|{c}"] = dict(base_rate=basev, lgb_codes=lgb_codes,
                                    by_code={int(k): dict(act_rate=round(float(r["mean"]), 4),
                                                          n=int(r["size"])) for k, r in rates.iterrows()})
        print(f"{lab:17s} {c:9s} base act rate {basev:.4f} · by code "
              f"{ {int(k): (round(float(r['mean']),3), int(r['size'])) for k,r in rates.iterrows()} } "
              f"=> LGB codes {lgb_codes}")
    frames[lab] = (df, cols)

ADJ = {"raw": [], "demog": ["age_r", "educat", "hisprace2"],
       "demog+relig": ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]}

grid = []
for lab, (df, cols) in frames.items():
    for c in cols:
        codes = anchor[f"{lab}|{c}"]["lgb_codes"]
        fr = df.copy()
        fr["lgb"] = np.where(fr[c].isin(codes), 1.0,
                             np.where(fr[c].isin([1, 2, 3, 4]), 0.0, np.nan))
        for aname, aset in ADJ.items():
            i, a, n = net_slopes(fr, aset)
            grid.append(dict(site=lab, arm=c, adjust=aname, id_net=i, act_net=a, n=n,
                             lgb_rate=float(fr["lgb"].mean())))

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    f = lambda v: "  nan " if np.isnan(v) else f"{v:+.4f}"          # noqa: E731
    print(f"  {g['site']:17s} {g['arm']:9s} {g['adjust']:12s} id_net={f(g['id_net'])} "
          f"act_net={f(g['act_net'])}  lgb_rate={g['lgb_rate']:.4f}  n={g['n']:5d}")


def med(key, site=None, arm=None):
    v = [g[key] for g in grid if not np.isnan(g[key])
         and (site is None or g["site"] == site) and (arm is None or g["arm"] == arm)]
    return float(np.median(v)) if v else np.nan


id_all, act_all = med("id_net"), med("act_net")
print(f"\n  NSFG overall: id_net {id_all:+.4f}  act_net {act_all:+.4f}  "
      f"(ratio {id_all/act_all:.2f}x)" if act_all else "")

# ══ THE RANDOMISED CONTRAST — and the LGB rate it produces ═══════════════════════════
ab = {}
for lab in ("2017-2019 female", "2017-2019 male"):
    a_i, b_i = med("id_net", lab, "orient_a"), med("id_net", lab, "orient_b")
    a_r = float(np.median([g["lgb_rate"] for g in grid if g["site"] == lab and g["arm"] == "orient_a"]))
    b_r = float(np.median([g["lgb_rate"] for g in grid if g["site"] == lab and g["arm"] == "orient_b"]))
    ab[lab] = dict(id_a=a_i, id_b=b_i, contrast=a_i - b_i, lgb_rate_a=a_r, lgb_rate_b=b_r)
    print(f"  {lab}: id_net A {a_i:+.4f} · B {b_i:+.4f} · contrast {a_i-b_i:+.4f} "
          f"| LGB rate A {a_r:.4f} vs B {b_r:.4f}  <- randomised wording")

# ══ NULLS — one per quantity, each destroying ITSELF (`#920`'s repair, carried) ══════
ref_lab = "2017-2019 female"
ref_df, ref_cols = frames[ref_lab]
ref = ref_df.copy()
ref["lgb"] = np.where(ref["orient_a"].isin(anchor[f"{ref_lab}|orient_a"]["lgb_codes"]), 1.0,
                      np.where(ref["orient_a"].isin([1, 2, 3, 4]), 0.0, np.nan))
ref = ref.dropna(subset=["lgb", "act", NORM])


def perm_null(which, reps=150):
    other = "act" if which == "lgb" else "lgb"
    vals = []
    for _ in range(reps):
        p = ref.copy()
        p[which] = p.groupby(other)[which].transform(lambda s: RNG.permutation(s.to_numpy()))
        i, a, _ = net_slopes(p, [])
        v = i if which == "lgb" else a
        if not np.isnan(v):
            vals.append(v)
    return float(np.median(vals)), float(np.std(vals)), len(vals)


nid_m, nid_s, nid_k = perm_null("lgb")
nact_m, nact_s, nact_k = perm_null("act")
print(f"\n  null for id_net  (identity permuted within behaviour): {nid_m:+.4f} +/- {nid_s:.4f} ({nid_k})")
print(f"  null for act_net (behaviour permuted within identity): {nact_m:+.4f} +/- {nact_s:.4f} ({nact_k})")

# ══ THE A/B CONTRAST'S OWN NULL — permute the RANDOMISED label, which IS its null ════
ab_null = []
for lab in ("2017-2019 female", "2017-2019 male"):
    df, cols = frames[lab]
    d2 = df[df["orient_a"].isin([1, 2, 3, 4]) | df["orient_b"].isin([1, 2, 3, 4])].copy()
    codes_a = anchor[f"{lab}|orient_a"]["lgb_codes"]
    codes_b = anchor[f"{lab}|orient_b"]["lgb_codes"]
    d2["armlab"] = np.where(d2["orient_a"].isin([1, 2, 3, 4]), "A", "B")
    d2["lgb"] = np.where(d2["armlab"] == "A", d2["orient_a"].isin(codes_a).astype(float),
                         d2["orient_b"].isin(codes_b).astype(float))
    d2 = d2.dropna(subset=["lgb", "act", NORM])
    for _ in range(120):
        p = d2.copy()
        p["armlab"] = RNG.permutation(p["armlab"].to_numpy())
        ia, _, _ = net_slopes(p[p.armlab == "A"], [])
        ib, _, _ = net_slopes(p[p.armlab == "B"], [])
        if not (np.isnan(ia) or np.isnan(ib)):
            ab_null.append(ia - ib)
abn_m, abn_s = float(np.median(ab_null)), float(np.std(ab_null))
print(f"  null for the A/B contrast (RANDOMISED label permuted): {abn_m:+.4f} +/- {abn_s:.4f} "
      f"({len(ab_null)} draws)  <- the exact null, because assignment really is random")

# ══ POSITIVE CONTROL — plant into the NULL world (`#920`'s repair, carried) ══════════
sweep_id, sweep_act = [], []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vi, va = [], []
    for _ in range(10):
        p = ref.copy()
        p["lgb"] = p.groupby("act")["lgb"].transform(lambda s: RNG.permutation(s.to_numpy()))
        p[NORM] = np.clip(p[NORM] - g * p["lgb"], 1, 5)      # identity ONLY; low = agrees
        i, a, _ = net_slopes(p, [])
        if not np.isnan(i):
            vi.append(-i)                                    # sign so "more agreement" is positive
        if not np.isnan(a):
            va.append(-a)
    sweep_id.append((g, float(np.median(vi)) if vi else np.nan))
    sweep_act.append((g, float(np.median(va)) if va else np.nan))
print(f"  positive sweep (planted into the NULL world), IDENTITY only:")
print(f"     id_net  {[(g, round(v, 4)) for g, v in sweep_id]}")
print(f"     act_net {[(g, round(v, 4)) for g, v in sweep_act]}  (must NOT rise)")

ps = [2 * (1 - stats.norm.cdf(abs((g["id_net"] - nid_m) / (nid_s or 1e-9)))) for g in grid
      if not np.isnan(g["id_net"])] + \
     [2 * (1 - stats.norm.cdf(abs((g["act_net"] - nact_m) / (nact_s or 1e-9)))) for g in grid
      if not np.isnan(g["act_net"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

contrasts = [v["contrast"] for v in ab.values() if not np.isnan(v["contrast"])]
max_contrast = float(np.max(np.abs(contrasts))) if contrasts else np.nan
id_res = abs(id_all - nid_m) > 2 * nid_s
act_res = abs(act_all - nact_m) > 2 * nact_s
wording_moves = not np.isnan(max_contrast) and abs(max_contrast - abn_m) > 2 * abn_s

G = Gate("Does `#920` survive a second instrument, and a randomised re-wording of identity?")
G.plant_direction_from_sweep("positive: an IDENTITY-only plant raises id_net, and g=0 is null",
                             sweep_id, baseline=-nid_m, baseline_spread=nid_s)
G.asserted("the identity-only plant does NOT raise act_net",
           abs(sweep_act[-1][1] - sweep_act[0][1]) < abs(sweep_id[-1][1] - sweep_id[0][1]),
           f"act_net moved {sweep_act[-1][1]-sweep_act[0][1]:+.4f} vs id_net "
           f"{sweep_id[-1][1]-sweep_id[0][1]:+.4f}", kind="control")
G.negative_control("identity permuted WITHIN behaviour strata", abs(nid_m), abs(id_all),
                   null_spread=nid_s, null_kind="within-behaviour identity-label permutation")
G.multiplicity_control("the whole site x arm x adjustment grid", ps, 0.05)
G.asserted("the identity CODING was anchored on behaviour, not guessed from labels that do not ship",
           all(v["lgb_codes"] for v in anchor.values()),
           f"LGB codes per arm { {k: v['lgb_codes'] for k, v in anchor.items()} }", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("the randomised wording DID change who counts as LGB (else the gauge test is vacuous)",
           any(abs(v["lgb_rate_a"] - v["lgb_rate_b"]) > 0.01 for v in ab.values()),
           f"LGB rates { {k: (round(v['lgb_rate_a'], 4), round(v['lgb_rate_b'], 4)) for k, v in ab.items()} }",
           kind="control")
G.asserted("KILL: W_WORDING requires the randomised A/B contrast to exceed its own permutation null",
           not wording_moves,
           f"max |id_net(A) - id_net(B)| {max_contrast:+.4f} vs null {abn_m:+.4f} +/- {abn_s:.4f} "
           f"(kind of null: randomised-arm label permutation)")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif wording_moves:
    VERDICT, WORLD = "OVERTURNED", "W_WORDING · identity's coupling moves with a randomised re-wording"
elif not id_res:
    VERDICT, WORLD = "OVERTURNED", "W_NOREP · the decomposition does not replicate on NSFG"
else:
    VERDICT, WORLD = "CONFIRMED", "W_STABLE · replicates, and survives a randomised re-wording"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=921, round="E03·A112·R359", verdict=VERDICT, world=WORLD,
           estimand="identity and behaviour contributions to the same-sex norm, each net of the "
                    "other, ceiling-normalised; plus the randomised orient_a vs orient_b contrast",
           instrument="NSFG 2011-2013 and 2017-2019 (a SECOND instrument; `#920` was GSS)",
           coding_anchor=anchor, grid=grid, ab=ab,
           nsfg=dict(id_net=id_all, act_net=act_all, id_resolved=bool(id_res),
                     act_resolved=bool(act_res)),
           nulls=dict(id_median=nid_m, id_sd=nid_s, id_draws=nid_k,
                      act_median=nact_m, act_sd=nact_s, act_draws=nact_k,
                      ab_median=abn_m, ab_sd=abn_s, ab_draws=len(ab_null)),
           positive_sweep_id=sweep_id, positive_sweep_act=sweep_act,
           max_ab_contrast=max_contrast, wording_moves=bool(wording_moves), family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "identity_under_randomised_wording.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'identity_under_randomised_wording.json'}")
