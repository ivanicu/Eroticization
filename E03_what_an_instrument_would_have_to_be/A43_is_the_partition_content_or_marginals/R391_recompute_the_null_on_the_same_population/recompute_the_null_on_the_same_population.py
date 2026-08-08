#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A131·R391 — the reversal, with each threshold's null recomputed on that threshold's population
===================================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#957`①/`#956`②. `#956` reported that the 2+2 partition contrast REVERSES as the
                discrimination threshold tightens: Δ +0.1086 (not-all-equal) → −0.0630 (range ≥ 2)
                → −0.2086 (within-person sd above median). That trajectory is now written onto four
                page blocks. **Nothing has established that it is about people.**

⚠⚠ THE GRADIENT CHECK FOUND THE DEFECT BEFORE THE RUN, IN MY OWN CODE, AND IT IS THE REASON THIS
   ROUND EXISTS. `#956`'s spec-curve rows filtered every NULL draw with `nonflat(...)` — the
   LOOSEST threshold — while filtering the OBSERVED arm with the row's own, tighter threshold:

       Xs = X_ALL[mask]                                    # observed: this row's threshold
       res_S = contrast(Xs) - contrast(Sn[nonflat(Sn)])    # null: the loosest threshold  ⚠

   So for `range ≥ 2` and `sd > median` the null was measured on a LESS SELECTED population than
   the observed arm. Δ falls as selection tightens, so the null sat too HIGH and the residual was
   pushed too NEGATIVE — **exactly the direction of the reported reversal.** That is `#954`'s error
   one level down, and it is the error `#956`'s own headline row was built to avoid. Second defect,
   same rows: 1,000,000 swaps where the headline used 3,200,000, i.e. an UNDER-MIXED null, which
   pushes the residual the other way. Two defects, opposite signs, neither controlled.

⚠ G1 — THE ESTIMAND. For a population P(t) defined by a discrimination rule t, the **content
  residual** ρ(t) = Δ_obs(P(t)) − E[Δ_null(P(t))], where the null is computed on P(t) and **t is
  re-applied to every null draw**. The question is whether ρ(t) changes SIGN as t tightens. This is
  not `#956`'s quantity: `#956`'s spec rows compared across two different populations.

⚠ THE ARITHMETIC TRAP, LABELLED BEFORE THE RUN. Person-centred residuals over k items are ipsative:
  Σᵢ Rᵢ = 0 per person ⇒ Σ_{i<j} Cov(Rᵢ,Rⱼ) = −½ Σᵢ Var(Rᵢ). With four items and equal residual
  variances the six pairwise correlations SUM to −2, hence

        Δ = W/2 − C/4  with  W + C = S   ⇒   Δ = (3/4)·W − S/4

  where W is the sum of the two within-cluster correlations. **Given S, Δ carries ONE degree of
  freedom.** So "the partition contrast" is an affine rescaling of "how strongly premarital–same-sex
  and teenage–extramarital are co-held", and any reversal is W falling. This is a DERIVATION; it is
  checked numerically at every threshold and it is not offered as evidence of anything.

Live Worlds    W_REAL_INVERSION · ρ(t) is significantly NEGATIVE at tight t and the placebo
                                  trajectory is flat ⇒ **among people who draw the sharpest
                                  distinctions, the crossing pairs are held TOGETHER** — a real and
                                  different structure, and A131 has an object.
               W_SELECTION      · the placebo trajectory reproduces the reversal ⇒ the trajectory is
                                  produced by SELECTING on a function of the answers, not by content.
               W_NULL_MISMATCH  · with t re-applied to the null and the chain converged, ρ(t) does
                                  NOT reverse ⇒ **`#956`'s spec curve was the two defects above**,
                                  and the sentence now on four page blocks must be corrected again.
                                  ⚠ **The unwelcome one, and it is the one I staked against.**

Prediction     W_REAL_INVERSION → ρ(tight) < −2 spreads; placebo-A range < 25% of observed range.
Matrix         W_SELECTION      → placebo-A reproduces the trajectory (range ≥ 25% of observed).
               W_NULL_MISMATCH  → |ρ(tight)| < 2 spreads once the null shares the population.

Strongest      **THE THRESHOLD SELECTS ON THE NORM OF THE RESIDUAL VECTOR.** Within-person sd is
confound       (up to the z-scaling) the norm of the very residuals whose correlations Δ is built
(written       from, so selecting on it is conditioning on a quantity made of the same numbers.
before)        ⇒ CONTROL, SAME ITERATION — **PLACEBO-A, and it is exact, not approximate**: permuting
               each person's four answers AMONG the items leaves that person's multiset unchanged,
               therefore leaves all-equal status, range and within-person sd EXACTLY unchanged. **The
               same individuals are selected at every threshold**, while which item got which answer
               is destroyed. If the trajectory survives that, it is the selection.

Controls       POSITIVE: plant a graded partition into the tightest population's null and sweep;
                 `g=0` must sit ON that null (`#922`).
               NEGATIVE: both nulls reproduce across independent draw sets; swap margin drift
                 measured (row and column), threshold re-applied per draw, post-filter n reported.
               PLACEBO-A (selection): within-person item-label permutation — selection identical by
                 construction, content destroyed.
               PLACEBO-B (removal): a RANDOM size-matched removal at EVERY threshold, not only the
                 loosest (`#956` ran this for one row).
               MIXING: swaps scaled per person, with a convergence curve on the tightest population.
               MULTIPLICITY: 6 thresholds × 2 nulls = 12.
               SPEC CURVE (G4): threshold × null × centring (z-scored / raw), all cells published.
               SEEDS: 3 independent seeds per null; the sign of ρ(tight) reported per seed.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **say WHY someone draws sharp distinctions** — discrimination, engagement, response style and
    genuine moral differentiation are indistinguishable on four items;
  (2) ⚠ **hold person totals and full column shape at once on 4-point ipsative data** — `#954`'s wall,
    inherited, and it is why two nulls are reported rather than one;
  (3) ⚠ **cross-instrument replication** — the four sexual-norm items are GSS's; there is
    **only this one instrument** carrying them, so no second site can be asked;
  (4) ⚠ **separate "sharp discriminators are different people" from "sharp discrimination is a
    different state of the same people"** — one observation per respondent, no panel;
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
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
K = len(ITEMS)
PAIRS = list(combinations(range(K), 2))
CLUSTER_A, CLUSTER_B = {0, 3}, {1, 2}
WITHIN = [p for p in PAIRS if set(p) <= CLUSTER_A or set(p) <= CLUSTER_B]
CROSS = [p for p in PAIRS if p not in WITHIN]
SWAPS_PER_PERSON = 300          # `#951` converged at 3.2M for n=14,847 = 215/person; 300 over-mixes
SEEDS = (391, 1391, 2391)

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"]).copy()
waves = sorted(int(y) for y in d.year.unique())
X_ALL = d[ITEMS].to_numpy(dtype=float)
N_ALL = len(X_ALL)
MU, SD = X_ALL.mean(axis=0), X_ALL.std(axis=0)   # z-scale FIXED on the full sample (`#956`)
sd_person = X_ALL.std(axis=1)
SD_MED, SD_P75 = float(np.median(sd_person)), float(np.quantile(sd_person, 0.75))

print(f"HARD RULE 1 — GSS ballot 1, all four sexual-norm items jointly: n={N_ALL} respondents, "
      f"{len(waves)} waves {waves[0]}-{waves[-1]} · item means {np.round(MU,3).tolist()}")
print("⚠ HARD RULE 2 — instrument: GSS ballot 1, four sexual-norm items, one block, the SAME "
      "instrument as `#949`-`#956`. What changes is the population rule and where the null is "
      "computed. No second instrument carries these items (registered, not planned).")
print(f"⚠ within-person sd thresholds fixed on the observed full sample: median {SD_MED:.4f} · "
      f"p75 {SD_P75:.4f} — applied as RULES to every null draw, never as counts")


# ── the population rules. Each is a function of a person's raw answer vector ──────────
def r_all(M):     return np.ones(len(M), bool)
def r_nonflat(M): return M.std(axis=1) > 0
def r_rng2(M):    return (M.max(axis=1) - M.min(axis=1)) >= 2
def r_rng3(M):    return (M.max(axis=1) - M.min(axis=1)) >= 3
def r_sdmed(M):   return M.std(axis=1) > SD_MED
def r_sdp75(M):   return M.std(axis=1) > SD_P75


RULES = [("all respondents", r_all), ("not all equal", r_nonflat), ("range >= 2", r_rng2),
         ("range >= 3", r_rng3), ("sd > median", r_sdmed), ("sd > p75", r_sdp75)]


def pair_rs(M, raw=False):
    Z = M if raw else (M - MU) / SD
    R = Z - Z.mean(axis=1, keepdims=True)
    return {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}


def contrast(M, raw=False):
    if len(M) < 200:
        return np.nan
    r = pair_rs(M, raw)
    return float(np.mean([r[p] for p in WITHIN]) - np.mean([r[p] for p in CROSS]))


def swap_null(M, rng, per_person=SWAPS_PER_PERSON):
    """checkerboard swap — person totals and item column SUMS exactly invariant."""
    m = M.copy()
    n = len(m)
    todo, B = int(per_person * n), 20000
    while todo > 0:
        b = min(B, todo)
        pp, qq = rng.integers(0, n, b), rng.integers(0, n, b)
        ii, jj = rng.integers(0, K, b), rng.integers(0, K, b)
        for p, q, i, j in zip(pp, qq, ii, jj):
            if p != q and i != j and m[p, i] > 1 and m[p, j] < 4 and m[q, i] < 4 and m[q, j] > 1:
                m[p, i] -= 1; m[p, j] += 1; m[q, i] += 1; m[q, j] -= 1
        todo -= b
    return m


def column_null(M, rng):
    """per-column permutation — each item's FULL marginal exactly invariant."""
    m = M.copy()
    for i in range(K):
        m[:, i] = rng.permutation(m[:, i])
    return m


def item_label_perm(M, rng):
    """PLACEBO-A — permute each person's four answers AMONG items.
    Leaves every person's multiset, hence all-equal status / range / sd, EXACTLY unchanged."""
    m = M.copy()
    for i in range(len(m)):
        m[i] = rng.permutation(m[i])
    return m


# ══ 0 · PLACEBO-A IS EXACT — verify the invariance before using it ═══════════════════
_p = item_label_perm(X_ALL, np.random.default_rng(7))
sel_ident = all(bool((rule(_p) == rule(X_ALL)).all()) for _, rule in RULES)
print(f"\n  ⚠ PLACEBO-A invariance check: item-label permutation selects the IDENTICAL individuals "
      f"under all {len(RULES)} rules → {sel_ident}  (max |Δsd| "
      f"{np.abs(_p.std(axis=1) - sd_person).max():.2e})")

PERM = [item_label_perm(X_ALL, np.random.default_rng(300 + j)) for j in range(6)]

# ══ 1 · THE MAIN GRID — the null is recomputed on EACH threshold's own population ════
rows, all_ps = [], []
for tag, rule in RULES:
    mask = rule(X_ALL)
    Xs = X_ALL[mask]
    d_obs = contrast(Xs)
    d_obs_raw = contrast(Xs, raw=True)
    rr = pair_rs(Xs)
    S6 = float(sum(rr.values()))
    W2 = float(rr[WITHIN[0]] + rr[WITHIN[1]])
    d_pred = 0.75 * W2 - S6 / 4.0                 # the DERIVATION, checked not assumed

    cells = {}
    for nname, fn in (("NULL-S", lambda M, r: swap_null(M, r)),
                      ("NULL-C", lambda M, r: column_null(M, r))):
        vals, per_seed, kept, drift_c, drift_r = [], [], [], [], []
        for s in SEEDS:
            M = fn(Xs, np.random.default_rng(s))
            drift_c.append(float(np.abs(M.sum(axis=0) - Xs.sum(axis=0)).max()))
            drift_r.append(float(np.abs(M.sum(axis=1) - Xs.sum(axis=1)).max()))
            keep = rule(M)                        # ⚠ THE REPAIR: this row's OWN rule, per draw
            kept.append(int(keep.sum()))
            v = contrast(M[keep])
            vals.append(v); per_seed.append(v)
        # independent draw set for reproducibility
        rep = []
        for s in SEEDS[:2]:
            Mr = fn(Xs, np.random.default_rng(50000 + s))
            rep.append(contrast(Mr[rule(Mr)]))
        mu_n, sd_n = float(np.mean(vals)), float(np.std(vals, ddof=1))
        res = d_obs - mu_n
        z = res / max(sd_n, 1e-9)
        cells[nname] = dict(null=mu_n, null_sd=sd_n, residual=res, z=z,
                            repro=abs(mu_n - float(np.mean(rep))),
                            post_n=float(np.mean(kept)), removed=int(mask.sum()) - int(np.mean(kept)),
                            col_drift=max(drift_c), row_drift=max(drift_r),
                            res_per_seed=[d_obs - v for v in per_seed])
        all_ps.append(2 * (1 - stats.norm.cdf(abs(z))))

    # PLACEBO-B — random size-matched removal reaching this row's n, from the full sample
    keepfrac = mask.mean()
    pb = [contrast(X_ALL[np.random.default_rng(900 + j).random(N_ALL) < keepfrac]) for j in range(12)]
    # PLACEBO-A — identical individuals, item→answer assignment destroyed
    pa = [contrast(P[mask]) for P in PERM]

    rows.append(dict(spec=tag, n=int(mask.sum()), share=float(mask.mean()),
                     d_obs=d_obs, d_obs_raw=d_obs_raw, S6=S6, W2=W2, d_pred=d_pred,
                     ident_err=abs(d_pred - d_obs),
                     placebo_A=float(np.mean(pa)), placebo_B=float(np.mean(pb)),
                     resS_per_seed=cells["NULL-S"]["res_per_seed"],
                     **{f"{k}_{f}": v[f] for k, v in cells.items()
                        for f in ("null", "null_sd", "residual", "z", "repro", "post_n",
                                  "removed", "col_drift", "row_drift")}))
    print(f"\n  {tag:<16s} n={int(mask.sum()):6d} ({mask.mean():.1%})  Δ {d_obs:+.4f}  "
          f"(raw-centred {d_obs_raw:+.4f})")
    print(f"      DERIVATION Σr₆ {S6:+.5f} · W {W2:+.5f} · (3/4)W−Σ/4 = {d_pred:+.5f} vs measured "
          f"{d_obs:+.5f} · |err| {abs(d_pred-d_obs):.2e}")
    for nname, c in cells.items():
        print(f"      {nname}  null {c['null']:+.4f} ± {c['null_sd']:.4f} · residual "
              f"{c['residual']:+.4f} · z {c['z']:+.2f} · repro {c['repro']:.5f} · post-filter n "
              f"{c['post_n']:.0f} (rule re-applied per draw, {c['removed']:+d} vs observed) · "
              f"drift col {c['col_drift']:.0f} row {c['row_drift']:.0f}")
    print(f"      PLACEBO-A (same people, item labels permuted) Δ {np.mean(pa):+.4f} · "
          f"PLACEBO-B (random size-matched) Δ {np.mean(pb):+.4f}")

# ══ 2 · MIXING CURVE on the tightest population ══════════════════════════════════════
Xt = X_ALL[RULES[-1][1](X_ALL)]
mixing = []
for pp in (50, 150, 300, 600):
    M = swap_null(Xt, np.random.default_rng(4242), per_person=pp)
    mixing.append([pp, contrast(M[RULES[-1][1](M)])])
print(f"\n  mixing curve on the tightest population (n={len(Xt)}): "
      f"{[(p, round(v,4)) for p, v in mixing]} — converged if the last two agree")

# ══ 3 · POSITIVE CONTROL on the tightest population, planted INTO its null ═══════════
base = swap_null(Xt, np.random.default_rng(77))
rng = np.random.default_rng(391)
sweep = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(8):
        P = base.copy()
        if gg:
            m = rng.random(len(P)) < gg
            sh = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
            for i in CLUSTER_A:
                P[m, i] += 0.5 * sh
            for i in CLUSTER_B:
                P[m, i] -= 0.5 * sh
        keep = RULES[-1][1](P)
        vals.append(contrast(P[keep]))
    sweep.append([float(gg), float(np.median(vals))])
sweep_neg = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(8):
        P = base.copy()
        if gg:
            m = rng.random(len(P)) < gg
            sh = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
            for i in CLUSTER_A:
                P[m, i] -= 0.5 * sh          # INVERTED plant
            for i in CLUSTER_B:
                P[m, i] += 0.5 * sh
        keep = RULES[-1][1](P)
        vals.append(contrast(P[keep]))
    sweep_neg.append([float(gg), float(np.median(vals))])
tight = rows[-1]
print(f"  INVERTED sweep on the tightest population (two-sided dose-response, `#950`'s error): "
      f"{[(g, round(v,4)) for g, v in sweep_neg]}")
print(f"  positive sweep on the tightest population (g=0 IS NULL-S): "
      f"{[(g, round(v,4)) for g, v in sweep]} · g=0 at {sweep[0][1]:+.4f} vs null "
      f"{tight['NULL-S_null']:+.4f} ± {tight['NULL-S_null_sd']:.4f} = "
      f"{abs(sweep[0][1]-tight['NULL-S_null'])/max(tight['NULL-S_null_sd'],1e-9):.2f} spreads")

# ══ 4 · BOOTSTRAP the observed Δ per threshold (persons resampled) ═══════════════════
boot = {}
for tag, rule in RULES:
    Xs = X_ALL[rule(X_ALL)]
    b = []
    for j in range(400):
        idx = np.random.default_rng(6000 + j).integers(0, len(Xs), len(Xs))
        b.append(contrast(Xs[idx]))
    boot[tag] = [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]
print("\n  bootstrap 95% interval on the observed Δ (persons resampled, 400 draws)")
for tag, ci in boot.items():
    print(f"    {tag:<16s} [{ci[0]:+.4f}, {ci[1]:+.4f}]")

# ══ 5 · THE VERDICT ═════════════════════════════════════════════════════════════════
obs_range = max(r["d_obs"] for r in rows) - min(r["d_obs"] for r in rows)
pa_range = max(r["placebo_A"] for r in rows) - min(r["placebo_A"] for r in rows)
pb_range = max(r["placebo_B"] for r in rows) - min(r["placebo_B"] for r in rows)
ident_ok = max(r["ident_err"] for r in rows) < 5e-3
mix_ok = abs(mixing[-1][1] - mixing[-2][1]) < 2 * tight["NULL-S_null_sd"]
pos_fires = sweep[-1][1] > sweep[0][1] + 2 * tight["NULL-S_null_sd"]
neg_fires = sweep_neg[-1][1] < sweep_neg[0][1] - 2 * tight["NULL-S_null_sd"]
pos_g0_ok = abs(sweep[0][1] - tight["NULL-S_null"]) < 2.0 * max(tight["NULL-S_null_sd"], 1e-4)
repro_ok = max(r[f"{n}_repro"] for r in rows for n in ("NULL-S", "NULL-C")) < 0.05
placebo_A_flat = pa_range < 0.25 * obs_range
controls_ok = pos_fires and neg_fires and pos_g0_ok and repro_ok and mix_ok

zS = tight["NULL-S_z"]
seed_signs = [float(np.sign(x)) for x in tight["resS_per_seed"]]
seed_agree = len(set(seed_signs)) == 1
if not controls_ok:
    world = "W_UNREADABLE"
elif not placebo_A_flat:
    world = "W_SELECTION"
elif zS < -2:
    world = "W_REAL_INVERSION"
else:
    world = "W_NULL_MISMATCH"

G = Gate("Does the partition contrast really reverse among sharp discriminators, once each "
         "threshold's null is computed on that threshold's own population?")
G.plant_direction_from_sweep(
    "positive: a planted partition raises Δ on the TIGHTEST population, and g=0 sits ON that "
    "population's own NULL-S (`#922`)", sweep,
    baseline=tight["NULL-S_null"], baseline_spread=max(tight["NULL-S_null_sd"], 1e-4))
G.negative_control(
    "both nulls reproduce across independent draw sets at every threshold; the binding cell is the "
    "tightest population's NULL-S",
    max(r[f"{n}_repro"] for r in rows for n in ("NULL-S", "NULL-C")),
    abs(tight["NULL-S_residual"]),
    null_spread=float(tight["NULL-S_null_sd"]),
    null_kind="checkerboard swap (person totals and item column sums exact) and per-column "
              "permutation (each item's full marginal exact), each computed ON that threshold's "
              "population with THAT threshold's rule re-applied to every draw")
G.asserted(
    "⚠ TWO-SIDED DOSE-RESPONSE, because the claim is a NEGATIVE residual and a one-sided visibility "
    "check would say nothing about the direction actually being claimed (`#950`'s error): the same "
    "plant applied with the partition INVERTED must drive Δ down",
    neg_fires,
    f"inverted sweep {[(g, round(v,4)) for g, v in sweep_neg]} · g=0 {sweep_neg[0][1]:+.4f} → "
    f"g=0.75 {sweep_neg[-1][1]:+.4f} against a null spread of {tight['NULL-S_null_sd']:.4f}",
    kind="control", population=f"GSS ballot 1, tightest population n={len(Xt)}")
G.asserted(
    "⚠ SEED AGREEMENT on the sign of the staked residual", seed_agree,
    f"per-seed residuals on the tightest population vs NULL-S: "
    f"{[round(x,4) for x in tight['resS_per_seed']]} · signs agree {seed_agree}",
    kind="control", population=f"GSS ballot 1, tightest population n={tight['n']}")
G.multiplicity_control("6 thresholds x 2 nulls = 12 (`#936`②)", all_ps, 0.05,
                       labels=[f"{t}/{n}" for t, _ in RULES for n in ("S", "C")])
G.asserted(
    "⚠⚠ THE REPAIR, and it is what this round is for: `#956` filtered every spec-curve NULL draw "
    "with the LOOSEST rule while filtering the observed arm with the row's own rule, so the null "
    "sat on a less-selected population and the residual was pushed negative — the direction of the "
    "reported reversal. Here each row's own rule is re-applied to every null draw",
    all(abs(r[f"{n}_removed"]) < 0.25 * r["n"] for r in rows for n in ("NULL-S", "NULL-C")),
    " · ".join(f"{r['spec']}: S n={r['NULL-S_post_n']:.0f} ({r['NULL-S_removed']:+d}) "
               f"C n={r['NULL-C_post_n']:.0f} ({r['NULL-C_removed']:+d}) vs observed {r['n']}"
               for r in rows), kind="control", population=f"GSS ballot 1, n={N_ALL}")
G.asserted(
    "⚠ PLACEBO-A, the control for the strongest confound: permuting each person's four answers "
    "among the items leaves the multiset — hence all-equal status, range and within-person sd — "
    "EXACTLY unchanged, so the identical individuals are selected at every threshold while the "
    "item->answer assignment is destroyed. If the trajectory survives that, it is the SELECTION",
    sel_ident and placebo_A_flat,
    f"selection identical under all rules: {sel_ident} · placebo-A trajectory range {pa_range:.4f} "
    f"vs observed range {obs_range:.4f} ({pa_range/max(obs_range,1e-9):.1%}); flat iff < 25%",
    kind="control", population=f"GSS ballot 1, n={N_ALL}")
G.asserted(
    "⚠ PLACEBO-B: a RANDOM size-matched removal at EVERY threshold (`#956` ran this only for the "
    "loosest row) must not reproduce the trajectory",
    pb_range < 0.25 * obs_range,
    f"placebo-B range {pb_range:.4f} vs observed {obs_range:.4f} "
    f"({pb_range/max(obs_range,1e-9):.1%}) · per row " +
    " · ".join(f"{r['spec'][:12]} {r['placebo_B']:+.4f}" for r in rows),
    kind="control", population=f"GSS ballot 1, n={N_ALL}")
G.asserted(
    "⚠⚠ ARITHMETIC TRAP LABELLED: under the ipsative constraint Δ = (3/4)W − Σr₆/4, so given the "
    "sum of the six pairwise correlations the contrast carries ONE degree of freedom and is an "
    "affine rescaling of the two within-cluster correlations. This is a DERIVATION, verified "
    "numerically at every threshold, and it is not offered as evidence",
    ident_ok,
    " · ".join(f"{r['spec'][:12]} Σ{r['S6']:+.4f} W{r['W2']:+.4f} pred{r['d_pred']:+.5f} "
               f"obs{r['d_obs']:+.5f} err{r['ident_err']:.1e}" for r in rows),
    kind="control", population=f"GSS ballot 1, n={N_ALL}")
G.asserted(
    "⚠ MIXING: swaps scale per person (300/person; `#951` converged at 215/person), with a "
    "convergence curve on the tightest population — `#956`'s spec rows used a FIXED 1M swaps, which "
    "is under-mixed for the loose rows and over-mixed for the tight ones", mix_ok,
    f"curve {[(p, round(v,4)) for p, v in mixing]} · last two differ by "
    f"{abs(mixing[-1][1]-mixing[-2][1]):.4f} against a null spread of {tight['NULL-S_null_sd']:.4f}",
    kind="control", population=f"GSS ballot 1, tightest population n={len(Xt)}")
G.asserted(
    "⚠ MARGIN DRIFT measured, not assumed: the swap must hold person totals and item column sums "
    "exactly", max(r["NULL-S_col_drift"] for r in rows) == 0
    and max(r["NULL-S_row_drift"] for r in rows) == 0,
    f"max column drift {max(r['NULL-S_col_drift'] for r in rows):.0f} · max row drift "
    f"{max(r['NULL-S_row_drift'] for r in rows):.0f} across all 6 populations",
    kind="control", population=f"GSS ballot 1, n={N_ALL}")
G.asserted("the whole specification grid is published, disagreeing cells included; the bootstrap "
           "interval on each observed Δ is reported beside it", True,
           " · ".join(f"{r['spec'][:12]} Δ{r['d_obs']:+.3f} CI[{boot[r['spec']][0]:+.3f},"
                      f"{boot[r['spec']][1]:+.3f}] zS{r['NULL-S_z']:+.2f} zC{r['NULL-C_z']:+.2f}"
                      for r in rows), kind="control", population=f"GSS ballot 1, n={N_ALL}")

G.asserted(
    "KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires with g=0 on its null, "
    "both nulls reproduce, and the chain is mixed. STAKED: W_REAL_INVERSION, i.e. on the tightest "
    "population the residual against NULL-S is more than 2 spreads BELOW zero, with the placebo-A "
    "trajectory flat. ⚠ I stake this because the reversal is already written onto four page blocks; "
    "W_NULL_MISMATCH means `#956`'s spec curve was my own two coding defects and every one of those "
    "blocks needs correcting again",
    controls_ok and placebo_A_flat and zS < -2,
    f"positive fires {pos_fires} · inverted fires {neg_fires} · g=0 on null {pos_g0_ok} · nulls reproduce {repro_ok} · mixed "
    f"{mix_ok} · placebo-A flat {placebo_A_flat} · tightest residual vs NULL-S "
    f"{tight['NULL-S_residual']:+.4f} (z {zS:+.2f}), vs NULL-C {tight['NULL-C_residual']:+.4f} "
    f"(z {tight['NULL-C_z']:+.2f}) ⇒ {world}",
    kind="kill", yardstick="partition contrast residual on the tightest discriminating population, "
                           "NULL-S",
    yardstick_noise=float(tight["NULL-S_null_sd"]),
    population=f"GSS ballot 1, within-person sd > p75, n={tight['n']} of {N_ALL}, "
               f"{len(waves)} waves {waves[0]}-{waves[-1]}",
    direction="one-sided: W_REAL_INVERSION requires a NEGATIVE residual on the tightest population")

print(G)
verdict = (("UNVERIFIED" if not controls_ok else
            ("CONFIRMED" if (placebo_A_flat and zS < -2) else "OVERTURNED"))
           + f" · world {world} · tightest population n={tight['n']} Δ {tight['d_obs']:+.4f} "
             f"residual vs NULL-S {tight['NULL-S_residual']:+.4f} (z {zS:+.2f})")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=958, round="E03·A131·R391", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world in ("W_NULL_MISMATCH", "W_SELECTION")),
               n_all=int(N_ALL), waves=waves, sd_median=SD_MED, sd_p75=SD_P75,
               grid=rows, bootstrap=boot, mixing=mixing,
               placebo_A_range=pa_range, placebo_B_range=pb_range, observed_range=obs_range,
               selection_identical_under_placebo_A=bool(sel_ident),
               null_median=tight["NULL-S_null"], null_sd=tight["NULL-S_null_sd"],
               null_draws=len(SEEDS), positive_sweep=sweep, inverted_sweep=sweep_neg, family_size=len(all_ps),
               seed_signs=seed_signs,
               seeds=list(SEEDS), world=world, verdict=verdict),
          open(OUT / "recompute_the_null_on_the_same_population.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'recompute_the_null_on_the_same_population.json'}")
