#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A129·R389 — bracket the content residual, because no single null holds both constraints
============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#953`①/`#952`②. `#951` reported that **13.9%** of the sexual-norm partition contrast
                survived its swap null — the residual on which A126–A127 now rests. But that swap
                held each item's column **SUM** on 4-point data and let column **SHAPE** drift
                (`#951` registered the amount itself). **`#952` then showed that on BINARY items,
                where shape is held EXACTLY, the same kind of effect vanishes and REVERSES.** So the
                13.9% is measured against a null that did not control the variable `#952` proved was
                carrying the effect.

⚠⚠ AND NO SINGLE NULL HOLDS BOTH CONSTRAINTS ON ORDINAL DATA — which is the round's real content
   **NULL-S (swap):** person totals EXACT · column sums EXACT · column shape only approximate.
   **NULL-C (per-column permutation):** each item's FULL marginal distribution EXACT · all
   person↔item association destroyed · but the DISTRIBUTION OF PERSON TOTALS is not preserved.
   ⇒ Each controls what the other leaks. **The honest output is a BRACKET over the two, not a
   point**, and reporting a single number here would be `realstat`'s "a bound replaced by a point".

⚠ G1, THE ADMISSIBILITY QUESTION, ANSWERED BEFORE THE RUN: is NULL-C admissible even though it
   destroys the general factor? **Yes, and for an exact reason.** Person-centring subtracts each
   person's own mean, so any person-level CONSTANT — which is precisely what a general permissiveness
   factor is — contributes **identically zero** to the residuals. Residual structure is specific by
   construction, so a null that destroys the general factor removes nothing the estimand uses.
   ⚠ What NULL-C DOES change is the spread of person totals: people who answer all four alike have
   all-zero residuals, and there are fewer such people under independent column draws. **That is
   measured below, not assumed**, and it is why NULL-C alone cannot carry the verdict.

Live Worlds    W_RESIDUAL_REAL · the departure survives under BOTH nulls ⇒ `#951`'s residual stands,
                                 now bracketed rather than a point, and A126–A127 keep a finding.
               W_ALL_MARGINAL  · it vanishes under at least one ⇒ **A126–A127 produce nothing but
                                 "it was the marginals"**, and `#949`/`#950`/`#951` collapse
                                 together. ⚠ **The unwelcome one, and `#952` says it is live.**
               W_NULLS_DISAGREE· the two nulls give very different answers ⇒ neither is admissible
                                 alone and the bracket IS the finding. (the meta-separator: "content
                                 vs marginal" may not be separable at all on 4-point ipsative data)

Estimand       The same PARTITION CONTRAST as `#950`/`#951`: Δ = mean residual r over the 2
(G1)           within-cluster pairs minus mean over the 4 crossing pairs, on person-centred
               z-scores. **The content residual is `Δ_obs − Δ_null` under EACH marginal-preserving
               null, and the claim is the MINIMUM of the two** — a residual is only as real as the
               weakest null it survives.

Prediction     W_RESIDUAL_REAL  -> both residuals > 2× their own null spread and same sign.
Matrix         W_ALL_MARGINAL   -> at least one inside its null spread.
               W_NULLS_DISAGREE -> both non-zero but differing by more than either's spread.

Strongest      **NULL-C CHANGES THE PERSON-TOTAL DISTRIBUTION**, and people with flat profiles have
confound       zero residuals. If NULL-C has systematically fewer flat people, its residual
(written       correlations are computed on a differently-shaped population and the comparison is
before)        not like-for-like. ⇒ CONTROL, same iteration: the person-total distribution and the
               share of all-equal profiles are measured under observed / NULL-S / NULL-C and
               reported, and NULL-C is never allowed to carry the verdict alone.

Controls       NEGATIVE-S / NEGATIVE-C as above, side by side, because the difference IS the round.
               POSITIVE: plant a content partition into EACH null world and sweep; `g=0` sits on
                 that null by construction (`#922`, `#937`⑤).
               ⚠ MARGINAL CHECK per null: column means, full column histograms, row totals.
               ⚠ PERSON-TOTAL CHECK: the confound above, measured.
               ⚠ NULL-REPRODUCIBILITY: two independent draw sets per null (`#948`①).
               MULTIPLICITY: 6 pairs × 2 nulls + 2 contrasts = 14.
               SPEC CURVE (G4): 2 nulls × {z-scored, raw-centred}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **hold person totals AND full column shape simultaneously on 4-point ipsative data** — that
    is the reason this round outputs a bracket; a Patefield/rc-sample fixes both margins of a
    CONTINGENCY TABLE, not of a person × item matrix with a per-person constraint;
  (2) ⚠ **observe a person over time** — repeated cross-section;
  (3) ⚠ **no second instrument** — the four norms are GSS's; **only this one instrument** carries
    them;
  (4) ⚠ **rescue `#951`'s residual if it falls** — a downgrade is not repairable by a better null;
  (5) `[unchallenged]` — door ③.
"""
import json
import sys
import warnings
from itertools import combinations
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
RNG = np.random.default_rng(389)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
K = len(ITEMS)
PAIRS = list(combinations(range(K), 2))
CLUSTER_A, CLUSTER_B = {0, 3}, {1, 2}
WITHIN = [p for p in PAIRS if set(p) <= CLUSTER_A or set(p) <= CLUSTER_B]
CROSS = [p for p in PAIRS if p not in WITHIN]

# ⚠ `#951`'s numbers READ FROM ITS ARTIFACT, never typed (`#840`; blocked me twice at `#952`)
REF951 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                        "A127_is_the_partition_content_or_marginal" /
                        "R387_a_null_that_keeps_the_marginals" / "results" /
                        "a_null_that_keeps_the_marginals.json"))
D951, S951 = REF951["delta_obs"], REF951["null_median"]
RESID951 = D951 - S951
print(f"⚠ `#951` from its ARTIFACT: Δ_obs {D951:+.4f} · swap null {S951:+.4f} · residual "
      f"{RESID951:+.4f} = {RESID951 / D951:.1%} · its registered shape drift "
      f"{REF951['shape_drift']:.4f}")

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"]).copy()
waves = sorted(int(y) for y in d.year.unique())
X = d[ITEMS].to_numpy(dtype=float)
N = len(X)
POOLED_SD = X.std(axis=0)
print(f"HARD RULE 1 — ballot 1, all four jointly: n={N} · {len(waves)} waves {waves[0]}-{waves[-1]}")
print("⚠ HARD RULE 2 — instrument: GSS ballot 1, four sexual-norm items, one block. Identical to "
      "`#951`; ONLY THE NULL CHANGES, and that is the round.")


def pair_rs(mat):
    Z = (mat - mat.mean(axis=0)) / POOLED_SD
    R = Z - Z.mean(axis=1, keepdims=True)
    return {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}


def contrast(rs):
    return float(np.mean([rs[p] for p in WITHIN]) - np.mean([rs[p] for p in CROSS]))


def swap_null(mat, n_swaps, rng):
    """NULL-S: person totals EXACT, column sums EXACT, column shape approximate."""
    m = mat.copy()
    B, done = 20000, 0
    while done < n_swaps:
        b = min(B, n_swaps - done)
        pp, qq = rng.integers(0, N, b), rng.integers(0, N, b)
        ii, jj = rng.integers(0, K, b), rng.integers(0, K, b)
        for p, q, i, j in zip(pp, qq, ii, jj):
            if p != q and i != j and m[p, i] > 1 and m[p, j] < 4 and m[q, i] < 4 and m[q, j] > 1:
                m[p, i] -= 1; m[p, j] += 1; m[q, i] += 1; m[q, j] -= 1
        done += b
    return m


def column_null(mat, rng):
    """NULL-C: each item's FULL marginal EXACT, all association destroyed."""
    m = mat.copy()
    for i in range(K):
        m[:, i] = rng.permutation(m[:, i])
    return m


D_OBS = contrast(pair_rs(X))
print(f"\n  observed partition contrast Δ = {D_OBS:+.4f}  (`#951` had {D951:+.4f})")

nulls = {}
for name, fn in (("NULL-S swap", lambda r: swap_null(X, 3_200_000, r)),
                 ("NULL-C column-perm", lambda r: column_null(X, r))):
    reps = 6 if "swap" in name else 40
    A = [contrast(pair_rs(fn(np.random.default_rng(389 + 100 * k)))) for k in range(reps)]
    B = [contrast(pair_rs(fn(np.random.default_rng(7000 + 100 * k)))) for k in range(reps)]
    nulls[name] = dict(mean=float(np.mean(A)), sd=float(np.std(A)),
                       repro=abs(float(np.mean(A)) - float(np.mean(B))), reps=reps)
    v = nulls[name]
    v["residual"] = D_OBS - v["mean"]
    v["z"] = v["residual"] / max(v["sd"], 1e-9)
    print(f"  {name:<20s} null {v['mean']:+.4f} +/- {v['sd']:.4f} · residual {v['residual']:+.4f} "
          f"({v['residual'] / D_OBS:.1%} of Δ) · z {v['z']:+.2f} · repro {v['repro']:.5f}")

# ══ MARGINAL + PERSON-TOTAL CHECKS ══════════════════════════════════════════════════
S_chk = swap_null(X, 3_200_000, np.random.default_rng(11))
C_chk = column_null(X, np.random.default_rng(11))
def hist(m):
    return np.array([[(m[:, i] == v).mean() for v in (1, 2, 3, 4)] for i in range(K)])
h0 = hist(X)
checks = {}
for nm, M in (("NULL-S", S_chk), ("NULL-C", C_chk)):
    checks[nm] = dict(col_mean_drift=float(np.abs(M.mean(axis=0) - X.mean(axis=0)).max()),
                      shape_drift=float(np.abs(hist(M) - h0).max()),
                      row_total_drift=float(np.abs(M.sum(axis=1) - X.sum(axis=1)).max()),
                      total_sd=float(M.sum(axis=1).std()),
                      flat_share=float((M.std(axis=1) == 0).mean()))
print(f"\n  ⚠ WHAT EACH NULL ACTUALLY HOLDS — measured, not asserted "
      f"(observed: person-total sd {X.sum(axis=1).std():.4f}, flat profiles "
      f"{(X.std(axis=1) == 0).mean():.4f})")
for nm, c in checks.items():
    print(f"    {nm}: col-mean drift {c['col_mean_drift']:.2e} · SHAPE drift {c['shape_drift']:.4f} "
          f"· row-total drift {c['row_total_drift']:.2e} · person-total sd {c['total_sd']:.4f} · "
          f"flat profiles {c['flat_share']:.4f}")
print("    ⇒ NULL-S holds totals and sums, leaks SHAPE · NULL-C holds shape exactly, leaks the "
      "PERSON-TOTAL distribution. Neither is complete, which is why the output is a bracket.")

# ══ POSITIVE CONTROL, one per null ══════════════════════════════════════════════════
sweeps = {}
for nm, base in (("NULL-S swap", S_chk), ("NULL-C column-perm", C_chk)):
    sw = []
    for gg in (0.0, 0.25, 0.50, 0.75):
        vals = []
        for _ in range(8):
            P = base.copy()
            if gg:
                m = RNG.random(N) < gg
                sh = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
                for i in CLUSTER_A:
                    P[m, i] += 0.5 * sh
                for i in CLUSTER_B:
                    P[m, i] -= 0.5 * sh
            vals.append(contrast(pair_rs(P)))
        sw.append([float(gg), float(np.median(vals))])
    sweeps[nm] = sw
    print(f"  positive sweep [{nm}] (g=0 IS that null): {[(x, round(v, 4)) for x, v in sw]} · "
          f"g=0 at {sw[0][1]:+.4f} vs null {nulls[nm]['mean']:+.4f} +/- {nulls[nm]['sd']:.4f} = "
          f"{abs(sw[0][1] - nulls[nm]['mean']) / max(nulls[nm]['sd'], 1e-9):.2f} spreads")

# ══ THE BRACKET ═════════════════════════════════════════════════════════════════════
res = {k: v["residual"] for k, v in nulls.items()}
lo_name = min(res, key=lambda k: res[k])
bracket = (res[lo_name], max(res.values()))
weakest = nulls[lo_name]
print(f"\n  ⇒ CONTENT RESIDUAL BRACKET: [{bracket[0]:+.4f}, {bracket[1]:+.4f}] "
      f"= [{bracket[0] / D_OBS:.1%}, {bracket[1] / D_OBS:.1%}] of Δ · the WEAKEST null is "
      f"{lo_name} at z {weakest['z']:+.2f}")

grid = []
for nm, M in (("NULL-S", S_chk), ("NULL-C", C_chk)):
    for arm, f in (("z-scored", pair_rs),
                   ("raw-centred", lambda m: {p: float(np.corrcoef(
                       (m - m.mean(axis=1, keepdims=True))[:, p[0]],
                       (m - m.mean(axis=1, keepdims=True))[:, p[1]])[0, 1]) for p in PAIRS})):
        grid.append(dict(null=nm, arm=arm, d_obs=contrast(f(X)), d_null=contrast(f(M))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    {g_['null']:<8s} {g_['arm']:<12s} Δ_obs {g_['d_obs']:+.4f}  Δ_null "
          f"{g_['d_null']:+.4f}  residual {g_['d_obs'] - g_['d_null']:+.4f}")

obs_rs = pair_rs(X)
ps = []
for nm, M in (("NULL-S", S_chk), ("NULL-C", C_chk)):
    nrs = pair_rs(M)
    key = "NULL-S swap" if nm == "NULL-S" else "NULL-C column-perm"
    for p in PAIRS:
        ps.append(2 * (1 - stats.norm.cdf(abs((obs_rs[p] - nrs[p]) /
                                              max(nulls[key]["sd"], 1e-9)))))
ps += [2 * (1 - stats.norm.cdf(abs(v["z"]))) for v in nulls.values()]

G = Gate("Does `#951`'s content residual survive a null that holds each item's FULL marginal?")
G.plant_direction_from_sweep(f"positive [{lo_name}]: a planted partition raises Δ, and g=0 sits ON "
                             f"that null (`#922`)", sweeps[lo_name],
                             baseline=weakest["mean"], baseline_spread=max(weakest["sd"], 1e-4))
G.negative_control(f"both nulls reproduce themselves across two independent draw sets ({lo_name} is "
                   f"the binding one)", max(v["repro"] for v in nulls.values()),
                   abs(weakest["residual"]),
                   null_spread=float(np.mean([v["sd"] for v in nulls.values()])),
                   null_kind="two marginal-preserving randomisations: checkerboard swap (totals + "
                             "column sums exact) and per-column permutation (full marginal exact)")
G.multiplicity_control("6 pairs x 2 nulls + 2 contrasts = 14 (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"{n}/{ITEMS[a][:5]}-{ITEMS[b][:5]}" for n in ("S", "C")
                               for a, b in PAIRS] + ["S/Δ", "C/Δ"])
G.asserted("⚠⚠ NEITHER NULL HOLDS BOTH CONSTRAINTS, AND BOTH LEAKS ARE MEASURED — which is why the "
           "output is a BRACKET and not a point", True,
           " · ".join(f"{nm}: shape drift {c['shape_drift']:.4f}, person-total sd {c['total_sd']:.4f} "
                      f"(observed {X.sum(axis=1).std():.4f}), flat profiles {c['flat_share']:.4f} "
                      f"(observed {(X.std(axis=1) == 0).mean():.4f})" for nm, c in checks.items()),
           kind="control", population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")
G.asserted("⚠ G1: person-centring removes any person-level CONSTANT exactly, so a general "
           "permissiveness factor contributes identically zero to the residuals — which is what "
           "makes NULL-C admissible despite destroying it", True,
           "residual structure is specific by construction; what NULL-C does change is the "
           "person-total spread, measured in the row above and the reason it cannot carry the "
           "verdict alone", kind="control",
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{g_['null']}/{g_['arm'][:3]} resid {g_['d_obs'] - g_['d_null']:+.4f}"
                      for g_ in grid), kind="control",
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")

pos_fires = sweeps[lo_name][-1][1] > sweeps[lo_name][0][1] + 2 * weakest["sd"]
repro_ok = max(v["repro"] for v in nulls.values()) < 0.5 * abs(weakest["residual"])
both_survive = all(v["z"] > 2 for v in nulls.values()) and \
    all(np.sign(v["residual"]) == np.sign(D_OBS) for v in nulls.values())
disagree = abs(res["NULL-S swap"] - res["NULL-C column-perm"]) > \
    2 * max(v["sd"] for v in nulls.values())
world = ("W_RESIDUAL_REAL" if both_survive and not disagree else
         ("W_NULLS_DISAGREE" if both_survive and disagree else "W_ALL_MARGINAL"))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and both nulls "
           "reproduce. STAKED: W_RESIDUAL_REAL, i.e. the residual clears 2x its own spread under "
           "BOTH nulls with the same sign. ⚠ W_ALL_MARGINAL is the unwelcome branch and it collapses "
           "`#949`/`#950`/`#951` together",
           (pos_fires and repro_ok) and both_survive and not disagree,
           f"positive fires {pos_fires} · nulls reproduce {repro_ok} · "
           + " · ".join(f"{k}: residual {v['residual']:+.4f} z {v['z']:+.2f}"
                        for k, v in nulls.items())
           + f" · bracket [{bracket[0]:+.4f}, {bracket[1]:+.4f}] = [{bracket[0] / D_OBS:.1%}, "
             f"{bracket[1] / D_OBS:.1%}] of Δ · nulls disagree {disagree} ⇒ {world}",
           kind="kill", yardstick="partition-contrast residual under the WEAKEST marginal-preserving "
                                  "null",
           yardstick_noise=weakest["sd"],
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024",
           direction="one-sided: W_RESIDUAL_REAL requires a POSITIVE residual under both")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and repro_ok) else ('CONFIRMED' if world == 'W_RESIDUAL_REAL' else 'OVERTURNED')}"
           f" · world {world} · content residual bracket [{bracket[0] / D_OBS:.1%}, "
           f"{bracket[1] / D_OBS:.1%}] of Δ (`#951` reported {RESID951 / D951:.1%})")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=954, round="E03·A129·R389", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world == "W_ALL_MARGINAL"),
               n=int(N), waves=waves, delta_obs=D_OBS, ref951=dict(delta=D951, null=S951,
                                                                   residual=RESID951),
               nulls=nulls, checks=checks, bracket=list(bracket),
               bracket_share=[bracket[0] / D_OBS, bracket[1] / D_OBS],
               weakest_null=lo_name, grid=grid,
               null_median=weakest["mean"], null_sd=weakest["sd"], null_draws=weakest["reps"],
               positive_sweep=sweeps[lo_name], family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "bracket_the_content_residual.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'bracket_the_content_residual.json'}")
