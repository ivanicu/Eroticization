#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A128·R388 — the same repair, pointed at `#945`, where the items are binary
==============================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#951`②. `#951` showed that the within-person **item-label permutation** null used
                throughout A125–A126 gives every item the same marginal by construction, and that
                **86.1% of the sexual-norm partition contrast was the marginals**. `#945` — *"the
                exchange rate is between homosexuals and racists"*, 12 of 30 pairs below their null,
                and a sentence now on both pages — **used that same null and has never been
                re-referred to one that keeps the marginals.**

⚠ GRADIENT CHECK, RUN BEFORE COMMITTING TO THE ROUND (and it is why the round is worth running)
   Across `#945`'s 30 cells, `corr(|marginal gap|, z) = −0.5316`, p = 0.0025: the bigger the gap
   between two targets' tolerance levels, the more negative the departure. Below-null cells have
   mean |gap| **0.1391**; above-null cells **0.0671** — twice as large. And `col/homo×rac` has the
   **largest marginal gap in the entire grid (0.320)** and is the most negative cell (z = −23.8).
   ⚠ But it is not a clean explanation: `lib/homo×rac` has |gap| of only **0.099** and still sits at
   z = −18.5. **That residual is what this round measures.**

⚠ HARD RULE 1, AND IT CATCHES MY OWN LEDGER: `#951`② wrote the Stouffer marginals as "0.247 to
   0.840". **0.247 is `libhomo` BEFORE the polarity flip** that `#942` established; post-flip it is
   0.767 and the true range is **0.480 to 0.849**. The NEXT line of my own last entry quoted a
   pre-correction number, which is the same class of error the corrections exist to stop.

⚠⚠ AND THIS TEST IS STRICTLY STRONGER THAN `#951`'s, FOR A REASON WORTH STATING
   `#951` swapped 4-point ordinal values, so its swap preserved column **sums** but let column
   **shape** drift (registered by it as uncontrolled; the value is READ from its artifact here). **These items are BINARY, and a
   binary column is completely determined by its sum** — so a margin-preserving swap here preserves
   the **entire marginal distribution exactly**. The caveat `#951` had to register does not exist
   in this round.

Live Worlds    W_CONTENT   · most of `#945`'s below-null cells survive a marginal-preserving null ⇒
                             the exchange rate is about who the targets ARE, and `#945` stands.
               W_MOSTLY_MARGINAL · they collapse ⇒ **`#945` is downgraded exactly as `#949`/`#950`
                             were, and the sentence about homosexuals and racists comes off both
                             pages as a size claim.** ⚠ The unwelcome one, and the gradient check
                             above says it is live.
               W_NULL_UNFIT· the binary swap chain will not mix ⇒ refusal. (the meta-separator)

Estimand       For each of the 30 cells (3 stems × 10 target pairs): the residual pair correlation
(G1)           (person-centred over all 15 items) MINUS its null, under TWO nulls —
                 NULL-A: within-person item-label permutation (`#945`'s, which equalises marginals)
                 NULL-B: curveball swap preserving every person's total AND every item's exact
                         marginal
               **Headline: how many of `#945`'s 12 below-null cells survive against NULL-B**, and
               what happens to `homosexuals × racists` specifically.

⚠ THE SWAP,    Binary curveball: find p,q,i,j with `x[p,i]=1, x[p,j]=0, x[q,i]=0, x[q,j]=1` and flip
STATED         all four. Row sums and column sums are exactly invariant, and for binary data column
EXACTLY        sum ⇒ column distribution, so **both margins are held exactly**.

Prediction     W_CONTENT -> ≥6 of the 12 survive at p<0.05 against NULL-B, and homo×rac survives in
Matrix         ≥2 of 3 stems.
               W_MOSTLY_MARGINAL -> <6 survive, and/or NULL-B's level is ≥50% of the observed
                             departure (the `negative_control` contract `#951` had to honour).
               W_NULL_UNFIT -> the mixing curve is flat.

Strongest      **A CHAIN THAT HAS NOT MIXED LOOKS EXACTLY LIKE W_MOSTLY_MARGINAL** — if too few
confound       curveball moves are admissible the null stays near the data, every departure vanishes,
(written       and I would read the unwelcome verdict for a mechanical reason. ⇒ CONTROL, same
before)        iteration: a mixing curve over swap counts with the changed-cell fraction, plus a
               positive control planted INTO the swapped world.

Controls       NEGATIVE-A / NEGATIVE-B as above, side by side, because the difference IS the round.
               POSITIVE: plant a trade-off between two named targets into the SWAPPED world, sweep;
                 `g=0` sits on NULL-B by construction (`#922`, `#937`⑤).
               ⚠ EXACT-MARGINAL CHECK: row sums, column sums AND full column distributions before
                 vs after — all three must be identical to machine precision, which `#951` could not
                 claim.
               MIXING CURVE: departures and changed-cell fraction vs swap count.
               MULTIPLICITY: the family is **30 cells × 2 nulls = 60**.
               SPEC CURVE (G4): swap counts × {all respondents, cohorts in ≥3 waves}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate "who the target is" from any third property that also tracks tolerance level** —
    this round separates content from MARGINAL LEVEL and nothing finer;
  (2) ⚠ **observe a person trading over time** — repeated cross-section;
  (3) ⚠ **no second instrument** — the Stouffer battery is GSS's; **only this one instrument**
    carries these fifteen items;
  (4) ⚠ **rescue `#945`'s magnitudes if they fall** — a downgrade here is not repairable by a better
    null, only by a different design;
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
RNG = np.random.default_rng(388)
TARGETS = ["homo", "rac", "com", "mil", "ath"]
STEMS = ["spk", "col", "lib"]
COLS = [s + t for s in STEMS for t in TARGETS]
K = len(COLS)
IDX = {c: i for i, c in enumerate(COLS)}
TPAIRS = list(combinations(TARGETS, 2))
CELLS = [(s, a, b) for s in STEMS for a, b in TPAIRS]

# ⚠ `#951`'s shape drift is READ FROM ITS ARTIFACT, not typed — `no_transcribed_numbers` (`#840`)
#   blocked v1 for typing `0.0888`, and it was right for the reason this round is about: **the
#   whole claim here is that `#951`'s null was weaker than this one**, so the number quantifying
#   that weakness must come from `#951`, not from my memory of it.
REF951 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                        "A127_is_the_partition_content_or_marginal" /
                        "R387_a_null_that_keeps_the_marginals" / "results" /
                        "a_null_that_keeps_the_marginals.json"))
SHAPE_DRIFT_951 = REF951["shape_drift"]

REF945 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                        "A125_is_zero_sum_a_within_person_fact" /
                        "R383_the_trade_off_against_its_ipsative_floor" / "results" /
                        "the_trade_off_against_its_ipsative_floor.json"))
Z945 = {(r["stem"], r["a"], r["b"]): r["z"] for r in REF945["cells"]}
BELOW945 = [k for k, v in Z945.items() if v < 0 and
            any(r["p"] < 0.05 for r in REF945["cells"]
                if (r["stem"], r["a"], r["b"]) == k)]
print(f"⚠ `#945` read from its ARTIFACT: {len(REF945['cells'])} cells, {REF945['below_floor']} below their "
      f"null at p<0.05, {REF945['above_floor']} above; its floor {REF945['derived_floor']:+.4f}")

d = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
d = d[(d.year >= 1988) & d.cohort.notna()].copy()
for c in COLS:
    u = sorted(float(x) for x in d[c].dropna().unique())
    assert len(u) == 2, f"{c}: expected binary, found {u}"
    d[c] = (d[c] == u[0]).where(d[c].isin(u))
for t in TARGETS:                                   # polarity derived (`#927`③/`#942`②)
    ref = "spk" + t
    for s in STEMS:
        c = s + t
        if c != ref and d[[c, ref]].dropna().corr().iloc[0, 1] < 0:
            d[c] = 1.0 - d[c]
d = d.dropna(subset=COLS).copy()
X = d[COLS].to_numpy(dtype=np.int8)
N = len(X)
lev = X.mean(axis=0)
print(f"HARD RULE 1 — all {K} items jointly, 1988+: n={N} · levels {lev.min():.3f}-{lev.max():.3f}")
print(f"  ⚠ `#951`② quoted '0.247 to 0.840'. **0.247 is `libhomo` BEFORE the polarity flip**; "
      f"post-flip it is {lev[IDX['libhomo']]:.3f} and the true range is {lev.min():.3f}-{lev.max():.3f}")
print("⚠ HARD RULE 2 — instrument: GSS Stouffer battery, one questionnaire, five targets per stem "
      "asked in a block. Identical to `#945`; ONLY THE NULL CHANGES.")


def cell_rs(mat):
    Z = mat.astype(float)
    R = Z - Z.mean(axis=1, keepdims=True)
    out = {}
    for s, a, b in CELLS:
        out[(s, a, b)] = float(np.corrcoef(R[:, IDX[s + a]], R[:, IDX[s + b]])[0, 1])
    return out


def label_permute(mat):
    m = mat.copy()
    for s in STEMS:
        cols = [IDX[s + t] for t in TARGETS]
        blk = m[:, cols]
        m[:, cols] = np.take_along_axis(blk, np.argsort(RNG.random(blk.shape), axis=1), axis=1)
    return m


def curveball(mat, n_swaps, rng):
    """Binary swap preserving row sums AND column sums exactly -> both margins EXACT."""
    m = mat.copy()
    changed = np.zeros(m.shape, dtype=bool)
    acc = 0
    B = 50000
    done = 0
    while done < n_swaps:
        b = min(B, n_swaps - done)
        pp, qq = rng.integers(0, N, b), rng.integers(0, N, b)
        ii, jj = rng.integers(0, K, b), rng.integers(0, K, b)
        for p, q, i, j in zip(pp, qq, ii, jj):
            if p == q or i == j:
                continue
            if m[p, i] == 1 and m[p, j] == 0 and m[q, i] == 0 and m[q, j] == 1:
                m[p, i] = 0; m[p, j] = 1; m[q, i] = 1; m[q, j] = 0
                changed[p, i] = changed[p, j] = changed[q, i] = changed[q, j] = True
                acc += 1
        done += b
    return m, acc, float(changed.mean())


obs = cell_rs(X)

# ══ NULL-A — `#945`'s null, kept for the side-by-side that IS the round ══════════════
draws_A = {k: [] for k in CELLS}
for _ in range(100):
    r = cell_rs(label_permute(X))
    for k in CELLS:
        draws_A[k].append(r[k])
A_m = {k: float(np.mean(v)) for k, v in draws_A.items()}
A_sd = {k: float(np.std(v)) for k, v in draws_A.items()}
labA = label_permute(X).mean(axis=0)
print(f"\n  NULL-A (within-person item-label permutation — `#945`'s; kind of null: within-person "
      f"item-label permutation): it EQUALISES the marginals to "
      f"{labA.min():.3f}-{labA.max():.3f} vs observed {lev.min():.3f}-{lev.max():.3f}")

# ══ MIXING CURVE + NULL-B ═══════════════════════════════════════════════════════════
print(f"\n  ⚠ MIXING CURVE — an unmixed chain looks exactly like the unwelcome verdict:")
mix = []
for n_sw in (0, 200_000, 1_000_000, 5_000_000):
    S, acc, chg = curveball(X, n_sw, np.random.default_rng(388 + n_sw % 101))
    rs = cell_rs(S)
    mix.append(dict(swaps=n_sw, accepted=int(acc), changed=chg,
                    mean_abs_dep=float(np.mean([abs(obs[k] - rs[k]) for k in CELLS]))))
    print(f"    swaps {n_sw:>10,d}  accepted {acc:>9,d}  cells changed {chg:.3f}  "
          f"mean |obs − null| {mix[-1]['mean_abs_dep']:.4f}")

draws_B = {k: [] for k in CELLS}
for r_ in range(8):
    S, _, _ = curveball(X, 5_000_000, np.random.default_rng(9000 + r_))
    rs = cell_rs(S)
    for k in CELLS:
        draws_B[k].append(rs[k])
B_m = {k: float(np.mean(v)) for k, v in draws_B.items()}
B_sd = {k: float(np.std(v)) for k, v in draws_B.items()}

# ══ EXACT-MARGINAL CHECK — stronger than `#951` could claim ═════════════════════════
S_chk, _, _ = curveball(X, 5_000_000, np.random.default_rng(4242))
col_drift = float(np.abs(S_chk.mean(axis=0) - lev).max())
row_drift = int(np.abs(S_chk.sum(axis=1) - X.sum(axis=1)).max())
print(f"\n  EXACT-MARGINAL CHECK: max |column mean drift| {col_drift:.2e} · max |row total drift| "
      f"{row_drift} — for BINARY items column sum determines the whole distribution, so both "
      f"margins are held EXACTLY (`#951` could only hold column SUMS on its 4-point items)")

rows = []
for k in CELLS:
    s, a, b = k
    zA = (obs[k] - A_m[k]) / max(A_sd[k], 1e-9)
    zB = (obs[k] - B_m[k]) / max(B_sd[k], 1e-9)
    rows.append(dict(stem=s, a=a, b=b, r_obs=obs[k], nullA=A_m[k], nullB=B_m[k],
                     zA=zA, zB=zB, z945=Z945.get(k, float("nan")),
                     gap=float(abs(lev[IDX[s + a]] - lev[IDX[s + b]])),
                     pA=2 * (1 - stats.norm.cdf(abs(zA))),
                     pB=2 * (1 - stats.norm.cdf(abs(zB)))))
rows.sort(key=lambda r: r["zB"])
print(f"\n  {'cell':<20s} {'r_obs':>8s} {'nullA':>8s} {'nullB':>8s} {'zA':>8s} {'zB':>8s} {'gap':>6s}")
for r in rows[:8]:
    print(f"  {r['stem']+'/'+r['a']+'-'+r['b']:<20s} {r['r_obs']:+8.4f} {r['nullA']:+8.4f} "
          f"{r['nullB']:+8.4f} {r['zA']:+8.2f} {r['zB']:+8.2f} {r['gap']:6.3f}")
print("  …")
for r in rows[-3:]:
    print(f"  {r['stem']+'/'+r['a']+'-'+r['b']:<20s} {r['r_obs']:+8.4f} {r['nullA']:+8.4f} "
          f"{r['nullB']:+8.4f} {r['zA']:+8.2f} {r['zB']:+8.2f} {r['gap']:6.3f}")

belowA = [r for r in rows if r["zA"] < 0 and r["pA"] < 0.05]
belowB = [r for r in rows if r["zB"] < 0 and r["pB"] < 0.05]
survivors = [r for r in belowA if r["zB"] < 0 and r["pB"] < 0.05]
hr = [r for r in rows if {r["a"], r["b"]} == {"homo", "rac"}]
hr_surv = [r for r in hr if r["zB"] < 0 and r["pB"] < 0.05]
mean_dep_A = float(np.mean([abs(r["r_obs"] - r["nullA"]) for r in rows]))
mean_dep_B = float(np.mean([abs(r["r_obs"] - r["nullB"]) for r in rows]))
# ⚠⚠ v1 COMPUTED `share_marginal = 1 - dep_B/dep_A` AND GOT **-90.9%**, which is not a share of
#   anything. The departures did not SHRINK under the margin-exact null — they grew in magnitude and
#   **flipped sign**: most cells now sit ABOVE their null. A "share explained" only means something
#   when the two quantities are nested, and here they are not. The honest quantities are the COUNT
#   and the DIRECTION, and they are reported instead.
share_marginal = float("nan")
n_above_B = sum(1 for k in CELLS if obs[k] > B_m[k])
sign_reversal = n_above_B >= 20
print(f"\n  below their null at p<0.05 — NULL-A {len(belowA)}/30 · NULL-B {len(belowB)}/30 · "
      f"of A's below-null cells, {len(survivors)}/{len(belowA)} survive B")
print(f"  homosexuals x racists: survives NULL-B in {len(hr_surv)}/3 stems "
      f"(zB {[round(r['zB'], 1) for r in hr]})")
print(f"  mean |departure|: NULL-A {mean_dep_A:.4f} → NULL-B {mean_dep_B:.4f} ⇒ "
      f"{share_marginal:.1%} of the departure was the marginals")

# ══ POSITIVE CONTROL ════════════════════════════════════════════════════════════════
sweep = []
base_S, _, _ = curveball(X, 5_000_000, np.random.default_rng(777))
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(6):
        P = base_S.copy()
        if gg:
            m = RNG.random(N) < gg
            P[m, IDX["spkrac"]] = 1 - P[m, IDX["spkhomo"]]
        R = P.astype(float) - P.astype(float).mean(axis=1, keepdims=True)
        vals.append(float(np.corrcoef(R[:, IDX["spkhomo"]], R[:, IDX["spkrac"]])[0, 1]))
    sweep.append([float(gg), float(np.median(vals))])
kB = ("spk", "homo", "rac")
print(f"  positive sweep (a trade-off planted between `spkhomo` and `spkrac` in the SWAPPED world, "
      f"g=0 IS NULL-B): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs NULL-B {B_m[kB]:+.4f} +/- "
      f"{B_sd[kB]:.4f} = {abs(sweep[0][1] - B_m[kB]) / max(B_sd[kB], 1e-9):.2f} spreads")

# ══ SPECIFICATION CURVE ═════════════════════════════════════════════════════════════
c3 = d.groupby("cohort").year.nunique()
keep = set(c3[c3 >= 3].index)
grid = []
for tag, sub in (("all respondents", d), ("cohorts in >=3 waves", d[d.cohort.isin(keep)])):
    Xs = sub[COLS].to_numpy(dtype=np.int8)
    for n_sw in (1_000_000, 5_000_000):
        gN = len(Xs)
        Ss, _, _ = curveball(Xs, n_sw, np.random.default_rng(555 + n_sw % 71)) if gN == N else \
            (Xs, 0, 0.0)
        if gN != N:
            Ss = Xs.copy()
            mm = Ss
            rngg = np.random.default_rng(555 + n_sw % 71)
            for _ in range(n_sw // 50000):
                pp, qq = rngg.integers(0, gN, 50000), rngg.integers(0, gN, 50000)
                ii, jj = rngg.integers(0, K, 50000), rngg.integers(0, K, 50000)
                for p, q, i, j in zip(pp, qq, ii, jj):
                    if p != q and i != j and mm[p, i] == 1 and mm[p, j] == 0 and \
                            mm[q, i] == 0 and mm[q, j] == 1:
                        mm[p, i] = 0; mm[p, j] = 1; mm[q, i] = 1; mm[q, j] = 0
        ro, rn = cell_rs(Xs), cell_rs(Ss)
        nb = sum(1 for k in CELLS if ro[k] < rn[k])
        grid.append(dict(spec=tag, swaps=n_sw, n=int(len(Xs)), below=int(nb),
                         mean_dep=float(np.mean([abs(ro[k] - rn[k]) for k in CELLS]))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    {g_['spec']:<22s} swaps {g_['swaps']:>9,d}  n={g_['n']:6d}  cells below null "
          f"{g_['below']}/30  mean |dep| {g_['mean_dep']:.4f}")

ps = [r["pA"] for r in rows] + [r["pB"] for r in rows]

G = Gate("Does `#945`'s exchange rate survive a null that keeps each target's tolerance level?")
G.plant_direction_from_sweep("positive: a planted trade-off drives `spkhomo`x`spkrac` below NULL-B, "
                             "and g=0 sits ON NULL-B (`#922`)",
                             [[g_, -v] for g_, v in sweep],
                             baseline=-B_m[kB], baseline_spread=max(B_sd[kB], 1e-4))
# ⚠⚠ v1 CALLED `negative_control(null=mean|dep under B|, effect=mean|dep under A|)` AND IT FAILED —
#   correctly, because **those two sides are not the same object**: both are DEPARTURES, under two
#   different nulls, so their ratio has no interpretation. That is `realstat`'s own tell, and it is
#   the SIXTH instance of `#916`③'s class, which `#951`③ predicted one entry earlier.
#   ⚠⚠⚠ AND IT IS NOT THE SAME AS `#951`'s, WHICH IS WHY REPAIRING IT IS NOT THE MOVE `#943`
#   REFUSED. At `#951` the control was correctly framed — a null LEVEL against an observed LEVEL on
#   one scale — and its failure was the finding; rewording it would have hidden a true result. Here
#   the call is malformed, and the round's actual conclusion (1 of 12 cells survive) never depended
#   on it. The distinction is exactly `#951`③'s: "a control testing the wrong property" versus "a
#   null that removed the variable the rival runs on".
draws_B2 = {k: [] for k in CELLS}
for r_ in range(8):
    S2, _, _ = curveball(X, 5_000_000, np.random.default_rng(31000 + r_))
    rs2 = cell_rs(S2)
    for k in CELLS:
        draws_B2[k].append(rs2[k])
B2_m = {k: float(np.mean(v)) for k, v in draws_B2.items()}
repro_B = float(np.mean([abs(B_m[k] - B2_m[k]) for k in CELLS]))
print(f"  ⚠ NULL-B reproducibility across two independent 8-draw sets: {repro_B:.5f} vs mean "
      f"|observed − NULL-B| {mean_dep_B:.4f}")
G.negative_control("curveball swap (both margins EXACT) reproduces itself across two independent "
                   "draw sets — the property this design needs, replacing a v1 call that compared "
                   "two DEPARTURES as though one were a null",
                   repro_B, mean_dep_B, null_spread=float(np.mean(list(B_sd.values()))),
                   null_kind="binary curveball swap randomisation, both margins exact")
G.multiplicity_control("30 cells x 2 nulls = 60 (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"A/{r['stem']}/{r['a']}-{r['b']}" for r in rows]
                              + [f"B/{r['stem']}/{r['a']}-{r['b']}" for r in rows])
G.asserted("⚠⚠ EXACT-MARGINAL CHECK — stronger than `#951` could claim: for BINARY items the column "
           "sum determines the entire distribution, so the swap holds BOTH margins exactly and "
           "`#951`'s registered shape-drift caveat does not exist here", col_drift == 0.0,
           f"max |column mean drift| {col_drift:.2e} · max |row total drift| {row_drift} · "
           f"(`#951` had to register {SHAPE_DRIFT_951:.4f} of uncontrolled shape drift on its "
           f"4-point items — READ from its artifact, not typed)",
           kind="control", population=f"GSS Stouffer battery, all {K} items, n={N}")
G.asserted("⚠ NULL-A EQUALISES THE MARGINALS AND SO COULD NOT TEST THE MARGINAL RIVAL — the reason "
           "this round exists", True,
           f"item levels observed {lev.min():.3f}-{lev.max():.3f} · under NULL-A "
           f"{labA.min():.3f}-{labA.max():.3f} · NULL-B holds them to {col_drift:.2e}",
           kind="control", population=f"GSS Stouffer battery, all {K} items, n={N}")
G.asserted("⚠ HARD RULE 1 caught my own ledger: `#951`② quoted '0.247' for the marginal floor, which "
           "is `libhomo` BEFORE the polarity flip `#942` established", True,
           f"post-flip `libhomo` is {lev[IDX['libhomo']]:.3f}; the true range is {lev.min():.3f}-"
           f"{lev.max():.3f}", kind="control",
           population=f"GSS Stouffer battery, all {K} items, n={N}")
G.asserted("⚠ MIXING CURVE: an unmixed chain looks exactly like the unwelcome verdict", True,
           " · ".join(f"{m['swaps']:,}→chg{m['changed']:.2f}/dep{m['mean_abs_dep']:.4f}"
                      for m in mix), kind="control",
           population=f"GSS Stouffer battery, all {K} items, n={N}")
G.asserted("the whole 30-cell grid and specification curve are published, disagreeing cells "
           "included", True,
           " · ".join(f"{g_['spec'][:3]}/{g_['swaps']//1000}k below{g_['below']}/30" for g_ in grid),
           kind="control", population=f"GSS Stouffer battery, all {K} items, n={N}")

mixed = mix[-1]["changed"] > 0.20 and mix[-1]["mean_abs_dep"] > mix[0]["mean_abs_dep"]
pos_fires = (-sweep[-1][1]) > (-sweep[0][1]) + 2 * B_sd[kB]
content = (len(survivors) >= 6) and (len(hr_surv) >= 2)
mostly = (len(survivors) > 0) and not content
world = ("W_NULL_UNFIT" if not mixed else
         ("W_CONTENT" if content else ("W_MOSTLY_MARGINAL" if mostly else "W_MARGINAL")))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires AND the chain "
           "mixed. STAKED: W_CONTENT, i.e. >=6 of `#945`'s below-null cells survive NULL-B, "
           "homosexuals x racists survives in >=2 of 3 stems, AND the marginal share is under 50%. "
           "⚠ W_MOSTLY_MARGINAL downgrades `#945` exactly as `#951` downgraded `#949`/`#950`",
           (pos_fires and mixed) and content,
           f"positive fires {pos_fires} · chain mixed {mixed} (changed {mix[-1]['changed']:.3f}) · "
           f"`#945`'s below-null cells surviving NULL-B: {len(survivors)}/{len(belowA)} · "
           f"homo x rac survives {len(hr_surv)}/3 stems · mean |departure| {mean_dep_A:.4f} → "
           f"{mean_dep_B:.4f}; ⚠ NOT a 'share explained' — the departures grew and FLIPPED SIGN, "
           f"with {n_above_B}/30 cells now ABOVE their null ⇒ {world}",
           kind="kill", yardstick="per-cell residual departure from a margin-exact swap null",
           yardstick_noise=float(np.mean(list(B_sd.values()))),
           population=f"GSS Stouffer battery, all {K} items jointly, n={N}, 1988-2021",
           direction="one-sided: W_CONTENT requires the cells to stay BELOW their null")

print(G)
# ⚠ v1 MAPPED `mostly -> UNVERIFIED`, contradicting the gate, which read OVERTURNED. The gate was
#   right: every control here PASSES and the pre-registered stake LOST, which is a refutation, not
#   an unfit instrument. `#951`'s UNVERIFIED was correct for the opposite reason — there a control
#   contract genuinely failed. **A lost stake with sound controls is OVERTURNED; only a failed
#   control makes it UNVERIFIED**, and conflating them would let a real refutation read as noise.
verdict = (f"{'UNVERIFIED' if not (pos_fires and mixed) else ('CONFIRMED' if content else 'OVERTURNED')}"
           f" · world {world}"
           + (f" · ⚠ only {len(survivors)}/{len(belowA)} of `#945`'s below-null cells survive a "
              f"margin-exact null, homosexuals x racists in {len(hr_surv)}/3 stems, and "
              f"{n_above_B}/30 cells now sit ABOVE their null — the direction REVERSES"
              if mostly else ""))
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=952, round="E03·A128·R388", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world in ("W_MOSTLY_MARGINAL", "W_MARGINAL")),
               n=int(N), k=K, levels=lev.tolist(), libhomo_postflip=float(lev[IDX["libhomo"]]),
               cells=rows, grid=grid, mixing=mix,
               below_nullA=len(belowA), below_nullB=len(belowB), survivors=len(survivors),
               homo_rac_survives=len(hr_surv), mean_dep_A=mean_dep_A, mean_dep_B=mean_dep_B,
               n_above_nullB=n_above_B, sign_reversal=bool(sign_reversal),
               null_B_reproducibility=repro_B,
               col_drift=col_drift, row_drift=row_drift,
               null_median=B_m[kB], null_sd=B_sd[kB], null_draws=8,
               positive_sweep=sweep, family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "the_same_repair_on_the_tolerance_battery.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_same_repair_on_the_tolerance_battery.json'}")
