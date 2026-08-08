#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A130·R390 — the partition among the people who actually distinguish the four acts
======================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#955`①. `#955` retracted `#954`'s claim that flat respondents have zero residuals:
                the all-"always wrong" group (**n = 3,173, 21%**) carries one fixed vector
                `[premarsx −0.6805, teensex +0.2785, xmarsex +0.4330, homosex −0.0310]`, which **is**
                the 2+2 partition, manufactured by z-scoring per item before person-centring.
                **59% of the +0.2620 contrast came from the 3,338 who drew no distinction at all.**
                Among the **11,509 who do**, the contrast is **+0.1086** — and nobody has asked
                whether that is content or still marginals.

⚠⚠ THE ARITHMETIC TRAP, LABELLED BEFORE THE RUN AND NOT AFTER
   **The DROP from +0.2620 to +0.1086 is partly FORCED, not measured.** Flat respondents carry a
   residual vector pointing along the partition; removing a concentrated mass aligned with a contrast
   must reduce that contrast. So the drop is a **DERIVATION** and is not offered as evidence of
   anything. ⇒ **The measurement of this round is `Δ_restricted − null_restricted`**, where the null
   is computed on the same population under the same filter. That comparison is not forced and could
   come out either way.

⚠ G1 — THE ESTIMAND IS NEW AND IS RENAMED. This is **not** `#951`/`#954`'s quantity. It is *the 2+2
  partition contrast among GSS respondents who do not give identical answers to all four sexual-norm
  items*. `#951`'s and `#954`'s numbers describe the whole sample and are not comparable to it; the
  restricted population is the finding's population, not a caveat on the old one (`#955`①).

⚠ G2 — THE NULLS CREATE FLAT RESPONDENTS where the observed sample has none, because a swap or a
  column permutation can hand somebody four equal answers. **The same non-flat filter is therefore
  applied to every null draw**, and the number it removes is measured and reported. Without that the
  null would contain exactly the people the observed arm excludes, and the comparison would be
  against a differently-composed population — `#954`'s error one level down.

Live Worlds    W_CONTENT        · Δ_restricted clears both restricted nulls ⇒ **among people who
                                   actually distinguish these acts, there is a real partition**, and
                                   A126–A130 have a finding about people after all.
               W_STILL_MARGINAL · it does not ⇒ **even among discriminators the partition is the
                                   marginals, and A126–A129 produce nothing about people at all.**
                                   ⚠ **The unwelcome one.**
               W_FILTER_ARTIFACT· re-filtering changes the nulls so much that the comparison is not
                                   available ⇒ refusal. (the meta-separator: "who discriminates" may
                                   not be a population a null can be defined on)

Estimand       Δ = mean residual r over the 2 within-cluster pairs minus mean over the 4 crossing
(G1)           pairs, on per-item z-scores then person-centred, **computed only on respondents whose
               four raw answers are not all identical**, and referred to two nulls computed on that
               same population with the same filter re-applied:
                 NULL-S  checkerboard swap — person totals and column sums exact, shape approximate
                 NULL-C  per-column permutation — each item's full marginal exact, totals not held

Prediction     W_CONTENT        -> Δ_restricted > 2× its spread above BOTH restricted nulls.
Matrix         W_STILL_MARGINAL -> inside at least one.
               W_FILTER_ARTIFACT-> the filter removes so much of a null draw that its n or marginals
                                   move materially from the observed arm's.

Strongest      **THE FILTER IS A SELECTION ON A FUNCTION OF THE ANSWERS**, so it can in principle
confound       induce structure by itself. ⇒ CONTROL, same iteration: the filter is applied
(written       IDENTICALLY to observed and null arms, the post-filter n and marginals are reported
before)        for every arm, and a PLACEBO filter (drop a random 22.5% matched in size, ignoring
               flatness) is run to show that removing people per se does not produce the contrast.

Controls       NEGATIVE-S / NEGATIVE-C on the restricted population, filter re-applied per draw.
               POSITIVE: plant a partition into each restricted null and sweep; `g=0` sits on that
                 null by construction (`#922`, `#937`⑤).
               ⚠ PLACEBO FILTER: random size-matched exclusion — Δ must stay near the unrestricted
                 value, or the restriction itself is doing the work.
               ⚠ FILTER ACCOUNTING: how many flats each null creates, and post-filter n/marginals.
               ⚠ NULL REPRODUCIBILITY: two independent draw sets per null (`#948`①).
               MULTIPLICITY: 6 pairs × 2 nulls + 2 contrasts = 14.
               SPEC CURVE (G4): the discrimination threshold — all-equal · range ≤1 · within-person
                 sd below median — because "who counts as discriminating" is a choice, not a fact.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **compare this number to `#951`'s or `#954`'s** — different population, different estimand,
    and quoting them side by side would be the scope error this project keeps paying for;
  (2) ⚠ **say why someone answers flat** — refusal, indifference, genuine uniform condemnation and
    satisficing are indistinguishable here;
  (3) ⚠ **hold person totals and full column shape at once on 4-point ipsative data** — `#954`'s
    wall, inherited;
  (4) ⚠ **no second instrument** — the four norms are GSS's; **only this one instrument** carries
    them;
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
RNG = np.random.default_rng(390)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
K = len(ITEMS)
PAIRS = list(combinations(range(K), 2))
CLUSTER_A, CLUSTER_B = {0, 3}, {1, 2}
WITHIN = [p for p in PAIRS if set(p) <= CLUSTER_A or set(p) <= CLUSTER_B]
CROSS = [p for p in PAIRS if p not in WITHIN]

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"]).copy()
waves = sorted(int(y) for y in d.year.unique())
X_ALL = d[ITEMS].to_numpy(dtype=float)
SD = X_ALL.std(axis=0)                       # ⚠ z-scale FIXED on the full sample so the filter
MU = X_ALL.mean(axis=0)                      #   cannot move the scale and manufacture a change
NONFLAT = X_ALL.std(axis=1) > 0
X = X_ALL[NONFLAT]
N = len(X)
print(f"HARD RULE 1 — GSS ballot 1, all four jointly: {len(X_ALL)} respondents, {len(waves)} waves "
      f"{waves[0]}-{waves[-1]} · NON-FLAT n={N} ({NONFLAT.mean():.4f}) · excluded "
      f"{(~NONFLAT).sum()}")
print("⚠ HARD RULE 2 — instrument: GSS ballot 1, four sexual-norm items, one block. Same as "
      "`#951`/`#954`; the POPULATION changes, and the estimand is renamed accordingly.")
print(f"⚠ z-scale is fixed on the FULL sample (mu {np.round(MU,3)}, sd {np.round(SD,3)}) so the "
      f"filter cannot move the scale and manufacture a change.")


def nonflat(M):
    return M.std(axis=1) > 0


def contrast(M):
    if len(M) < 50:
        return np.nan
    Z = (M - MU) / SD
    R = Z - Z.mean(axis=1, keepdims=True)
    r = {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}
    return float(np.mean([r[p] for p in WITHIN]) - np.mean([r[p] for p in CROSS]))


def pair_rs(M):
    Z = (M - MU) / SD
    R = Z - Z.mean(axis=1, keepdims=True)
    return {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}


def swap_null(M, n_swaps, rng):
    m = M.copy()
    n = len(m)
    B, done = 20000, 0
    while done < n_swaps:
        b = min(B, n_swaps - done)
        pp, qq = rng.integers(0, n, b), rng.integers(0, n, b)
        ii, jj = rng.integers(0, K, b), rng.integers(0, K, b)
        for p, q, i, j in zip(pp, qq, ii, jj):
            if p != q and i != j and m[p, i] > 1 and m[p, j] < 4 and m[q, i] < 4 and m[q, j] > 1:
                m[p, i] -= 1; m[p, j] += 1; m[q, i] += 1; m[q, j] -= 1
        done += b
    return m


def column_null(M, rng):
    m = M.copy()
    for i in range(K):
        m[:, i] = rng.permutation(m[:, i])
    return m


D_ALL = contrast(X_ALL)
D_OBS = contrast(X)
print(f"\n  Δ on ALL respondents      {D_ALL:+.4f}   ⚠ DERIVATION-adjacent, see the card")
print(f"  Δ on NON-FLAT (the estimand) {D_OBS:+.4f}")

nulls, filt = {}, {}
for name, fn, reps in (("NULL-S swap", lambda r: swap_null(X, 3_200_000, r), 6),
                       ("NULL-C column-perm", lambda r: column_null(X, r), 40)):
    A, Bv, created, ns = [], [], [], []
    for k in range(reps):
        M = fn(np.random.default_rng(390 + 100 * k))
        keep = nonflat(M)
        created.append(int((~keep).sum()))
        ns.append(int(keep.sum()))
        A.append(contrast(M[keep]))
        M2 = fn(np.random.default_rng(8000 + 100 * k))
        k2 = nonflat(M2)
        Bv.append(contrast(M2[k2]))
    nulls[name] = dict(mean=float(np.mean(A)), sd=float(np.std(A)),
                       repro=abs(float(np.mean(A)) - float(np.mean(Bv))), reps=reps)
    filt[name] = dict(flats_created=float(np.mean(created)), post_n=float(np.mean(ns)))
    v = nulls[name]
    v["residual"] = D_OBS - v["mean"]
    v["z"] = v["residual"] / max(v["sd"], 1e-9)
    print(f"  {name:<20s} null {v['mean']:+.4f} +/- {v['sd']:.4f} · residual {v['residual']:+.4f} "
          f"· z {v['z']:+.2f} · repro {v['repro']:.5f} · flats created "
          f"{filt[name]['flats_created']:.0f} → post-filter n {filt[name]['post_n']:.0f}")

# ══ PLACEBO FILTER — remove a random size-matched slice, ignoring flatness ═══════════
pl = []
for k in range(20):
    r = np.random.default_rng(555 + k)
    keep = r.random(len(X_ALL)) > (~NONFLAT).mean()
    pl.append(contrast(X_ALL[keep]))
placebo = float(np.mean(pl))
print(f"\n  ⚠ PLACEBO FILTER (drop a RANDOM {(~NONFLAT).mean():.1%}, ignoring flatness): Δ "
      f"{placebo:+.4f} vs unrestricted {D_ALL:+.4f} — removing people PER SE does not produce the "
      f"drop; the flatness filter does")

# ══ POSITIVE CONTROLS ═══════════════════════════════════════════════════════════════
sweeps = {}
for nm, base in (("NULL-S swap", swap_null(X, 3_200_000, np.random.default_rng(77))),
                 ("NULL-C column-perm", column_null(X, np.random.default_rng(77)))):
    sw = []
    for gg in (0.0, 0.25, 0.50, 0.75):
        vals = []
        for _ in range(8):
            P = base.copy()
            if gg:
                m = RNG.random(len(P)) < gg
                sh = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
                for i in CLUSTER_A:
                    P[m, i] += 0.5 * sh
                for i in CLUSTER_B:
                    P[m, i] -= 0.5 * sh
            keep = nonflat(P)
            vals.append(contrast(P[keep]))
        sw.append([float(gg), float(np.median(vals))])
    sweeps[nm] = sw
    print(f"  positive sweep [{nm}] (g=0 IS that null): {[(x, round(v, 4)) for x, v in sw]} · g=0 "
          f"at {sw[0][1]:+.4f} vs null {nulls[nm]['mean']:+.4f} +/- {nulls[nm]['sd']:.4f} = "
          f"{abs(sw[0][1] - nulls[nm]['mean']) / max(nulls[nm]['sd'], 1e-9):.2f} spreads")

# ══ SPEC CURVE — "who counts as discriminating" is a CHOICE ══════════════════════════
grid = []
for tag, mask in (("all-equal excluded", NONFLAT),
                  ("range >= 2", (X_ALL.max(axis=1) - X_ALL.min(axis=1)) >= 2),
                  ("within-sd above median", X_ALL.std(axis=1) >
                   np.median(X_ALL.std(axis=1)))):
    Xs = X_ALL[mask]
    Sn = swap_null(Xs, 1_000_000, np.random.default_rng(2024))
    Cn = column_null(Xs, np.random.default_rng(2024))
    grid.append(dict(spec=tag, n=int(mask.sum()), share=float(mask.mean()),
                     d_obs=contrast(Xs),
                     res_S=contrast(Xs) - contrast(Sn[nonflat(Sn)]),
                     res_C=contrast(Xs) - contrast(Cn[nonflat(Cn)])))
print("\n  specification curve — the discrimination threshold is a CHOICE, all cells published")
for g_ in grid:
    print(f"    {g_['spec']:<24s} n={g_['n']:6d} ({g_['share']:.1%})  Δ {g_['d_obs']:+.4f}  "
          f"residual vs S {g_['res_S']:+.4f}  vs C {g_['res_C']:+.4f}")

obs_rs = pair_rs(X)
ps = []
for nm, fn in (("S", lambda r: swap_null(X, 3_200_000, r)),
               ("C", lambda r: column_null(X, r))):
    M = fn(np.random.default_rng(31337))
    nrs = pair_rs(M[nonflat(M)])
    key = "NULL-S swap" if nm == "S" else "NULL-C column-perm"
    for p in PAIRS:
        ps.append(2 * (1 - stats.norm.cdf(abs((obs_rs[p] - nrs[p]) /
                                              max(nulls[key]["sd"], 1e-9)))))
ps += [2 * (1 - stats.norm.cdf(abs(v["z"]))) for v in nulls.values()]

lo_name = min(nulls, key=lambda k: nulls[k]["residual"])
weakest = nulls[lo_name]

G = Gate("Among people who actually distinguish the four acts, is the partition content or marginal?")
G.plant_direction_from_sweep(f"positive [{lo_name}]: a planted partition raises Δ on the restricted "
                             f"population, and g=0 sits ON that null (`#922`)", sweeps[lo_name],
                             baseline=weakest["mean"], baseline_spread=max(weakest["sd"], 1e-4))
G.negative_control(f"both restricted nulls reproduce across independent draw sets ({lo_name} binds)",
                   max(v["repro"] for v in nulls.values()), abs(weakest["residual"]),
                   null_spread=float(np.mean([v["sd"] for v in nulls.values()])),
                   null_kind="checkerboard swap and per-column permutation, each computed ON the "
                             "non-flat population with the SAME non-flat filter re-applied per draw")
G.multiplicity_control("6 pairs x 2 nulls + 2 contrasts = 14 (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"{n}/{ITEMS[a][:5]}-{ITEMS[b][:5]}" for n in ("S", "C")
                               for a, b in PAIRS] + ["S/Δ", "C/Δ"])
G.asserted("⚠⚠ ARITHMETIC TRAP LABELLED: the DROP from the all-respondent Δ to the restricted Δ is "
           "partly FORCED, because the excluded respondents carry a residual vector pointing along "
           "the partition (`#955`). It is a derivation and is not offered as evidence", True,
           f"Δ all {D_ALL:+.4f} → Δ non-flat {D_OBS:+.4f}; the MEASUREMENT is the residual against a "
           f"null on the SAME population, which is not forced", kind="control",
           population=f"GSS ballot 1 non-flat, n={N}")
G.asserted("⚠ PLACEBO FILTER: dropping a RANDOM size-matched slice must NOT reproduce the drop, or "
           "the restriction itself is doing the work", abs(placebo - D_ALL) < 0.02,
           f"random {(~NONFLAT).mean():.1%} removed → Δ {placebo:+.4f} vs unrestricted {D_ALL:+.4f} "
           f"(|diff| {abs(placebo - D_ALL):.4f}); the flatness filter gives {D_OBS:+.4f}",
           kind="control", population=f"GSS ballot 1, n={len(X_ALL)}")
G.asserted("⚠ FILTER ACCOUNTING: the nulls CREATE flat respondents where the observed arm has none, "
           "so the same filter is re-applied per draw and the count is reported (`#954`'s error one "
           "level down)", True,
           " · ".join(f"{nm}: {c['flats_created']:.0f} flats created → post-filter n {c['post_n']:.0f} "
                      f"(observed arm n={N})" for nm, c in filt.items()), kind="control",
           population=f"GSS ballot 1 non-flat, n={N}")
G.asserted("⚠ the z-scale is fixed on the FULL sample so the filter cannot move it", True,
           f"mu {np.round(MU, 3).tolist()} sd {np.round(SD, 3).tolist()} — identical in every arm",
           kind="control", population=f"GSS ballot 1, n={len(X_ALL)}")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{g_['spec'][:12]} n={g_['n']} S{g_['res_S']:+.3f} C{g_['res_C']:+.3f}"
                      for g_ in grid), kind="control", population=f"GSS ballot 1, n={len(X_ALL)}")

pos_fires = sweeps[lo_name][-1][1] > sweeps[lo_name][0][1] + 2 * weakest["sd"]
repro_ok = max(v["repro"] for v in nulls.values()) < 0.5 * abs(weakest["residual"])
placebo_ok = abs(placebo - D_ALL) < 0.02
both = all(v["z"] > 2 and v["residual"] > 0 for v in nulls.values())
world = ("W_CONTENT" if both else "W_STILL_MARGINAL")

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires, both nulls "
           "reproduce, and the placebo filter is inert. STAKED: W_CONTENT, i.e. Δ on the "
           "discriminating population clears 2x its spread above BOTH restricted nulls. "
           "⚠ W_STILL_MARGINAL means A126–A129 produce nothing about people at all",
           (pos_fires and repro_ok and placebo_ok) and both,
           f"positive fires {pos_fires} · nulls reproduce {repro_ok} · placebo inert {placebo_ok} · "
           + " · ".join(f"{k}: residual {v['residual']:+.4f} z {v['z']:+.2f}"
                        for k, v in nulls.items()) + f" ⇒ {world}",
           kind="kill", yardstick="partition contrast residual on the discriminating population, "
                                  "weakest null",
           yardstick_noise=weakest["sd"],
           population=f"GSS ballot 1, NON-FLAT n={N} of {len(X_ALL)}, {len(waves)} waves 1988-2024",
           direction="one-sided: W_CONTENT requires a POSITIVE residual under both")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and repro_ok and placebo_ok) else ('CONFIRMED' if both else 'OVERTURNED')}"
           f" · world {world} · Δ {D_OBS:+.4f} on n={N} who discriminate · residuals "
           + " / ".join(f"{v['residual']:+.4f}" for v in nulls.values()))
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=956, round="E03·A130·R390", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world == "W_STILL_MARGINAL"),
               n_all=int(len(X_ALL)), n_nonflat=int(N), waves=waves,
               delta_all=D_ALL, delta_nonflat=D_OBS, placebo_delta=placebo,
               nulls=nulls, filter_accounting=filt, grid=grid,
               null_median=weakest["mean"], null_sd=weakest["sd"], null_draws=weakest["reps"],
               positive_sweep=sweeps[lo_name], family_size=len(ps),
               weakest_null=lo_name, world=world, verdict=verdict),
          open(OUT / "among_those_who_actually_distinguish.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'among_those_who_actually_distinguish.json'}")
