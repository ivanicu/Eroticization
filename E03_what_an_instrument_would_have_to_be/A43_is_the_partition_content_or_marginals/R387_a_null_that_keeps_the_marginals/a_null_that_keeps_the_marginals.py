#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A127·R387 — a null that keeps each item's marginal, because mine destroyed it
==================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#949` and `#950` both CONFIRMED the same 2+2 partition — `{premarsx, homosex}` held
                together, `{teensex, xmarsex}` held together, all four crossing pairs traded — and
                `#950` killed the secular-drift rival. **Two confirmations in a row on one structure
                is a basin, and the rule says design the step whose positive outcome is unwelcome.**
                Here it is, and it is not hypothetical.

⚠⚠ THE RIVAL, AND WHY MY OWN NULL IS BLIND TO IT — measured before anything else runs
   Item marginals: `premarsx` 23.5% "always wrong" / 53.0% "not wrong at all" · `homosex` 50.4% /
   38.4% · `teensex` **65.5%** / 7.4% · `xmarsex` **77.6%** / 2.2%.
   **The two items I called "a wronged party" are EXACTLY the two heavily floored items**, and the
   two I called "consenting adults" are exactly the two unfloored ones. The partition coincides
   perfectly with response-distribution SHAPE.
   ⚠⚠⚠ **And the null used at `#949`/`#950` — within-person ITEM-LABEL permutation — makes all four
   items have the SAME marginal by construction.** A null that erases the difference cannot test
   whether the difference produced the result. Every departure those rounds measured is a departure
   from a world in which the marginals were already equalised, so **the marginal rival was never in
   the comparison at all.**

Live Worlds    W_CONTENT  · Δ survives against a null that PRESERVES each item's marginal ⇒ the
                            partition is about what the items mean, and `#949`/`#950` stand.
               W_MARGINAL · Δ collapses ⇒ **the blocks are the two floored items against the two
                            unfloored ones, and both confirmed rounds die**, along with the sentence
                            now on the page in two languages. ⚠ **The unwelcome one, and the
                            marginals above make it the live one.**
               W_NULL_UNFIT· the marginal-preserving randomisation cannot mix — with 77.6% of
                            `xmarsex` at the floor there may be too few admissible moves — so the
                            question is not answerable this way and the honest output is a refusal.
                            (the meta-separator: "content vs marginal" may not be separable in
                            ordinal data this floored)

Estimand       The same PARTITION CONTRAST as `#950`: Δ = mean residual r over the 2 within-cluster
(G1)           pairs minus mean over the 4 crossing pairs, on person-centred z-scores. **What
               changes is only the null**, and that is the whole round: Δ is now referred to a
               **swap-randomised** null that holds BOTH each person's total AND each item's column
               sum fixed, so the only thing destroyed is which person's answer sits where.

⚠ THE NULL,    Checkerboard swap: pick persons p,q and items i,j; apply
STATED         `x[p,i]−1 · x[p,j]+1 · x[q,i]+1 · x[q,j]−1` when every value stays in [1,4].
EXACTLY        **Row sums and column sums are exactly invariant** under this move, so each person's
               total (hence the ipsative constraint) and each item's mean are preserved while the
               person↔item association is destroyed. ⚠ It preserves column SUMS exactly; the column
               SHAPE can drift, so the shape drift is measured and reported rather than assumed.

Prediction     W_CONTENT  -> Δ clears the swap null by a margin comparable to the label-permutation
Matrix         one.
               W_MARGINAL -> Δ sits inside the swap null.
               W_NULL_UNFIT -> the swap chain does not mix: the pair correlations barely move from
                            the observed values however many swaps are applied.

Strongest      **A SWAP CHAIN THAT HAS NOT MIXED LOOKS EXACTLY LIKE A CONFIRMED FINDING.** If too
confound       few moves are admissible, the "null" is still nearly the observed data, Δ_null ≈
(written       Δ_obs, and I would read that as W_MARGINAL — the unwelcome verdict — for a purely
before)        mechanical reason. ⇒ CONTROL, same iteration: a MIXING CURVE over swap counts, plus
               the fraction of cells that actually changed, plus a positive control planted INTO the
               swapped world. **A chain that has not mixed and a real marginal explanation are
               distinguished by whether the planted partition is still recoverable.**

Controls       NEGATIVE 1: the `#949`/`#950` null (within-person item-label permutation) — kept so
                 the two nulls can be compared side by side, since the difference between them IS
                 the round.
               NEGATIVE 2: the swap null above, at several chain lengths.
               POSITIVE: plant a content partition into the SWAPPED world and sweep; `g=0` sits on
                 the swap null by construction (`#922`, `#937`⑤).
               MIXING CURVE: Δ_null and the changed-cell fraction as a function of swap count; a
                 chain that has not moved is reported as W_NULL_UNFIT, never as W_MARGINAL.
               MARGINAL CHECK: item means and full histograms before vs after, to show the swap did
                 what it claims and to measure the shape drift it does not control.
               MULTIPLICITY: 6 pairs × 2 nulls + 2 contrasts = 14 cells.
               SPEC CURVE (G4): swap counts × {z-scored, raw-centred}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **preserve the full column DISTRIBUTION, only its sum** — the swap fixes column totals; any
    residual shape drift is measured and reported, not controlled;
  (2) ⚠ **observe a person over time** — repeated cross-section;
  (3) ⚠ **no second instrument** — the four norms are GSS's; **only this one instrument** carries
    them;
  (4) ⚠ **distinguish "content" from any third property that also tracks the marginals** — this
    round separates content from MARGINAL SHAPE and nothing finer;
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
RNG = np.random.default_rng(387)
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
X = d[ITEMS].to_numpy(dtype=float)
N = len(X)
print(f"HARD RULE 1 — ballot 1, all four jointly: n={N} · {len(waves)} waves {waves[0]}-{waves[-1]}")
print("⚠ HARD RULE 2 — instrument: GSS ballot 1, four sexual-norm items in one block. Identical to "
      "`#949`/`#950`; ONLY THE NULL CHANGES, and that is the round.")
print("\n⚠⚠ THE RIVAL, measured first — the blocks coincide with response-distribution SHAPE:")
for i, c in enumerate(ITEMS):
    print(f"   {c:<10s} mean {X[:, i].mean():.3f} · P(=1 'always wrong') {(X[:, i] == 1).mean():.3f}"
          f" · P(=4 'not wrong at all') {(X[:, i] == 4).mean():.3f}"
          f"   [{'cluster A' if i in CLUSTER_A else 'cluster B'}]")
print("   ⇒ cluster B = the two heavily floored items · cluster A = the two unfloored ones")

POOLED_SD = X.std(axis=0)


def pair_rs(mat):
    Z = (mat - mat.mean(axis=0)) / POOLED_SD
    R = Z - Z.mean(axis=1, keepdims=True)
    return {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}


def contrast(rs):
    return float(np.mean([rs[p] for p in WITHIN]) - np.mean([rs[p] for p in CROSS]))


def label_permute(mat):
    return np.take_along_axis(mat, np.argsort(RNG.random(mat.shape), axis=1), axis=1)


def swap_randomise(mat, n_swaps, rng):
    """Checkerboard swap preserving BOTH row sums and column sums exactly."""
    m = mat.copy()
    changed = np.zeros(m.shape, dtype=bool)
    accepted = 0
    B = 20000
    done = 0
    while done < n_swaps:
        b = min(B, n_swaps - done)
        pp = rng.integers(0, N, b)
        qq = rng.integers(0, N, b)
        ii = rng.integers(0, K, b)
        jj = rng.integers(0, K, b)
        for p, q, i, j in zip(pp, qq, ii, jj):
            if p == q or i == j:
                continue
            if m[p, i] > 1 and m[p, j] < 4 and m[q, i] < 4 and m[q, j] > 1:
                m[p, i] -= 1; m[p, j] += 1; m[q, i] += 1; m[q, j] -= 1
                changed[p, i] = changed[p, j] = changed[q, i] = changed[q, j] = True
                accepted += 1
        done += b
    return m, accepted, float(changed.mean())


obs = pair_rs(X)
D_OBS = contrast(obs)
print(f"\n  observed partition contrast Δ = {D_OBS:+.4f}")

# ══ NEGATIVE 1 — the `#949`/`#950` null, for side-by-side comparison ═════════════════
lab = [contrast(pair_rs(label_permute(X))) for _ in range(120)]
lab_m, lab_sd = float(np.mean(lab)), float(np.std(lab))
print(f"  NULL 1 (within-person item-label permutation — the `#949`/`#950` null; kind of null: "
      f"within-person item-label permutation): {lab_m:+.4f} +/- {lab_sd:.4f}  ⇒ z "
      f"{(D_OBS - lab_m) / lab_sd:+.2f}")
lab_means = label_permute(X).mean(axis=0)
print(f"    ⚠ and it EQUALISES the marginals: item means under it "
      f"{np.round(lab_means, 3).tolist()} vs observed {np.round(X.mean(axis=0), 3).tolist()}")

# ══ MIXING CURVE + NEGATIVE 2 — the swap null at several chain lengths ═══════════════
print(f"\n  ⚠ MIXING CURVE — a chain that has NOT mixed looks exactly like a confirmed finding, so "
      f"the chain length is swept and the changed-cell fraction reported:")
mix = []
for n_sw in (0, 50_000, 200_000, 800_000, 3_200_000):
    reps = 6 if n_sw else 1
    vals, accs, chgs = [], [], []
    for r in range(reps):
        S, acc, chg = swap_randomise(X, n_sw, np.random.default_rng(387 + 1000 * r + n_sw % 97))
        vals.append(contrast(pair_rs(S)))
        accs.append(acc)
        chgs.append(chg)
    mix.append(dict(swaps=n_sw, delta=float(np.mean(vals)), sd=float(np.std(vals)),
                    accepted=float(np.mean(accs)), changed=float(np.mean(chgs))))
    m = mix[-1]
    print(f"    swaps {n_sw:>9,d}  Δ_null {m['delta']:+.4f} +/- {m['sd']:.4f}  accepted "
          f"{m['accepted']:>9,.0f}  cells changed {m['changed']:.3f}")
SW = mix[-1]
swap_m, swap_sd = SW["delta"], max(SW["sd"], 1e-6)
z_swap = (D_OBS - swap_m) / swap_sd
print(f"  ⇒ against the LONGEST swap chain: Δ_obs {D_OBS:+.4f} vs null {swap_m:+.4f} +/- "
      f"{swap_sd:.4f} ⇒ z {z_swap:+.2f}")

# ══ MARGINAL CHECK — did the swap do what it claims? ════════════════════════════════
S_final, _, _ = swap_randomise(X, 3_200_000, np.random.default_rng(999))
mean_drift = float(np.abs(S_final.mean(axis=0) - X.mean(axis=0)).max())
hist_obs = np.array([[(X[:, i] == v).mean() for v in (1, 2, 3, 4)] for i in range(K)])
hist_swp = np.array([[(S_final[:, i] == v).mean() for v in (1, 2, 3, 4)] for i in range(K)])
shape_drift = float(np.abs(hist_swp - hist_obs).max())
print(f"\n  MARGINAL CHECK: max |item mean drift| {mean_drift:.2e} (must be ~0 — column sums are "
      f"invariant by construction) · max |histogram cell drift| {shape_drift:.4f} (NOT controlled, "
      f"reported)")
row_drift = float(np.abs(S_final.sum(axis=1) - X.sum(axis=1)).max())
print(f"  ROW CHECK: max |person total drift| {row_drift:.2e} — the ipsative constraint is held")

# ══ POSITIVE CONTROL — plant a content partition INTO the swapped world ═════════════
sweep = []
base_S, _, _ = swap_randomise(X, 3_200_000, np.random.default_rng(4242))
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(8):
        P = base_S.copy().astype(float)
        if gg:
            m = RNG.random(N) < gg
            sh = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
            for i in CLUSTER_A:
                P[m, i] += 0.5 * sh
            for i in CLUSTER_B:
                P[m, i] -= 0.5 * sh
        vals.append(contrast(pair_rs(P)))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a content partition planted into the SWAPPED world, so g=0 IS the swap "
      f"null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs swap null {swap_m:+.4f} +/- "
      f"{swap_sd:.4f} = {abs(sweep[0][1] - swap_m) / swap_sd:.2f} spreads")

# ══ SPECIFICATION CURVE (G4) ════════════════════════════════════════════════════════
grid = []
for n_sw in (200_000, 3_200_000):
    S, _, _ = swap_randomise(X, n_sw, np.random.default_rng(31337 + n_sw % 91))
    for arm, f in (("z-scored", lambda M: pair_rs(M)),
                   ("raw-centred", lambda M: {p: float(np.corrcoef(
                       (M - M.mean(axis=1, keepdims=True))[:, p[0]],
                       (M - M.mean(axis=1, keepdims=True))[:, p[1]])[0, 1]) for p in PAIRS})):
        grid.append(dict(swaps=n_sw, arm=arm, delta_obs=contrast(f(X)), delta_null=contrast(f(S))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    swaps {g_['swaps']:>9,d}  {g_['arm']:<12s} Δ_obs {g_['delta_obs']:+.4f}  Δ_null "
          f"{g_['delta_null']:+.4f}")

ps = []
for p in PAIRS:
    ps.append(2 * (1 - stats.norm.cdf(abs((obs[p] - lab_m) / max(lab_sd, 1e-9)))))
S_ps, _, _ = swap_randomise(X, 3_200_000, np.random.default_rng(5150))
sw_rs = pair_rs(S_ps)
for p in PAIRS:
    ps.append(2 * (1 - stats.norm.cdf(abs((obs[p] - sw_rs[p]) / max(swap_sd, 1e-9)))))
ps += [2 * (1 - stats.norm.cdf(abs((D_OBS - lab_m) / max(lab_sd, 1e-9)))),
       2 * (1 - stats.norm.cdf(abs(z_swap)))]

G = Gate("Is the 2+2 partition about what the items MEAN, or about their response distributions?")
G.plant_direction_from_sweep("positive: a planted content partition raises Δ, and g=0 sits ON the "
                             "SWAP null this round judges against (`#922`)", sweep,
                             baseline=swap_m, baseline_spread=swap_sd)
G.negative_control("swap randomisation holding BOTH each person's total and each item's column sum "
                   "fixed [longest chain]", abs(swap_m), abs(D_OBS), null_spread=swap_sd,
                   null_kind="checkerboard swap randomisation preserving row sums AND column sums "
                             "(the marginal-preserving null the earlier rounds lacked)")
G.multiplicity_control("6 pairs x 2 nulls + 2 contrasts = 14 cells (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"lab/{ITEMS[a][:5]}-{ITEMS[b][:5]}" for a, b in PAIRS]
                              + [f"swap/{ITEMS[a][:5]}-{ITEMS[b][:5]}" for a, b in PAIRS]
                              + ["lab/Δ", "swap/Δ"])
G.asserted("⚠⚠ THE EARLIER NULL EQUALISES THE MARGINALS AND SO COULD NOT TEST THE MARGINAL RIVAL — "
           "stated as the reason this round exists", True,
           f"item means observed {np.round(X.mean(axis=0), 3).tolist()} · under within-person "
           f"label permutation {np.round(lab_means, 3).tolist()} (all four equalised) · the swap "
           f"null holds them to {mean_drift:.2e}", kind="control",
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")
G.asserted("⚠ MIXING CURVE: a chain that has not mixed looks exactly like the unwelcome verdict, so "
           "chain length is swept and the changed-cell fraction reported", True,
           " · ".join(f"{m['swaps']:,}→Δ{m['delta']:+.4f}/chg{m['changed']:.2f}" for m in mix),
           kind="control", population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")
G.asserted("⚠ MARGINAL AND ROW CHECKS: the swap does what it claims (column sums and person totals "
           "invariant); the residual histogram drift is reported, not controlled", True,
           f"max |item mean drift| {mean_drift:.2e} · max |person total drift| {row_drift:.2e} · "
           f"max |histogram cell drift| {shape_drift:.4f}", kind="control",
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{g_['swaps']//1000}k/{g_['arm'][:3]} obs{g_['delta_obs']:+.3f} "
                      f"null{g_['delta_null']:+.3f}" for g_ in grid), kind="control",
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024")

mixed = mix[-1]["changed"] > 0.20 and abs(mix[-1]["delta"] - mix[0]["delta"]) > 2 * swap_sd
pos_fires = sweep[-1][1] > sweep[0][1] + 2 * swap_sd
marginal_share = swap_m / D_OBS
# ⚠⚠ v1 SET `content = z_swap > 2` AND ITS VERDICT STRING PRINTED **CONFIRMED** WHILE THE GATE'S OWN
#   `negative_control` ROW FAILED — the row asks whether the null is SMALL RELATIVE TO THE EFFECT,
#   and here the null is **86% of it**. Two facts, both true, and v1 reported only the flattering
#   one: the residual clears the null's SAMPLING SPREAD (z=+5.80), and the null's LEVEL accounts for
#   most of the observed value.
#   ⚠⚠⚠ The tempting fix was to reword the control so it compares departures instead of levels.
#   **That is weakening a guard so my own round can pass** — the move refused at `#943` — and it
#   would be `#916`③'s SIXTH instance dressed as a repair. The contract is pre-existing and binding,
#   so the VERDICT is what changes: `content` now requires BOTH that the residual clear the null's
#   spread AND that the null be under half the observed, which is exactly what the control demands.
content = (z_swap > 2) and (marginal_share < 0.5)
mostly_marginal = (z_swap > 2) and (marginal_share >= 0.5)
# W_MOSTLY_MARGINAL was NOT pre-registered; it is added as a REFUSAL that forecloses the other
# three rather than competing with them, in the shape `#943` used for W_NO_NULL.
world = ("W_NULL_UNFIT" if not mixed else
         ("W_CONTENT" if content else
          ("W_MOSTLY_MARGINAL" if mostly_marginal else "W_MARGINAL")))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires AND the swap "
           "chain demonstrably mixed. STAKED: W_CONTENT, i.e. Δ clears the MARGINAL-PRESERVING null "
           "by >2 spreads. ⚠ W_MARGINAL is the unwelcome branch and it kills `#949` AND `#950` "
           "together, along with a sentence now on the page in two languages",
           (pos_fires and mixed) and content,
           f"positive fires {pos_fires} · chain mixed {mixed} (cells changed "
           f"{mix[-1]['changed']:.3f}, Δ moved {mix[0]['delta']:+.4f}→{mix[-1]['delta']:+.4f}) · "
           f"Δ_obs {D_OBS:+.4f} vs swap null {swap_m:+.4f} +/- {swap_sd:.4f} ⇒ z {z_swap:+.2f} "
           f"(label-permutation null gave z {(D_OBS - lab_m) / lab_sd:+.2f}) · ⚠ the null is "
           f"{marginal_share:.1%} OF THE OBSERVED VALUE, so the pre-existing negative_control "
           f"contract (null < 50% of effect) FAILS and the verdict honours it ⇒ {world}",
           kind="kill", yardstick="partition contrast Δ against a marginal-preserving swap null",
           yardstick_noise=swap_sd,
           population=f"GSS ballot 1, n={N}, {len(waves)} waves 1988-2024",
           direction="one-sided: W_CONTENT requires Δ ABOVE the swap null")

print(G)
verdict = (f"{'UNVERIFIED' if (not (pos_fires and mixed)) or mostly_marginal else ('CONFIRMED' if content else 'OVERTURNED')}"
           f" · world {world}"
           + (f" · ⚠ {marginal_share:.1%} of Δ is the MARGINALS; the content residual is "
              f"{D_OBS - swap_m:+.4f} = {1 - marginal_share:.1%} of what `#949`/`#950` reported"
              if mostly_marginal else ""))
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=951, round="E03·A127·R387", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_MARGINAL"),
               n=int(N), waves=waves, delta_obs=D_OBS,
               item_means=X.mean(axis=0).tolist(),
               p_floor=[float((X[:, i] == 1).mean()) for i in range(K)],
               p_ceil=[float((X[:, i] == 4).mean()) for i in range(K)],
               label_null_mean=lab_m, label_null_sd=lab_sd,
               label_null_item_means=lab_means.tolist(),
               null_median=swap_m, null_sd=swap_sd, null_draws=6,
               z_label=(D_OBS - lab_m) / lab_sd, z_swap=z_swap,
               marginal_share=marginal_share, content_residual=D_OBS - swap_m,
               mixing=mix, mean_drift=mean_drift, row_drift=row_drift, shape_drift=shape_drift,
               positive_sweep=sweep, grid=grid, family_size=len(ps),
               world=world, verdict=verdict),
          open(OUT / "a_null_that_keeps_the_marginals.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'a_null_that_keeps_the_marginals.json'}")
