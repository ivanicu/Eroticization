r"""#904 · E03·A108·R342 — leave-one-out: the one attempt `#111c` allows on "which act is special"

**COGNITIVE UPDATE CARD**
```
Core Gap        `#903` returned `UNVERIFIED` on "which act departs", and the reason was a design
                fault it named: **a per-item MEAN OVER PAIRS is not per-item** -- every pair
                involving a departing item sits in TWO items' averages, so one departing act lowers
                all four scores and "all four depart" is what a single departure looks like through
                that statistic. `#903`① specified the repair BEFORE this round existed, which is
                what makes running it a pre-registration rather than a threshold chosen after a
                result. `#903`② put `#111c` in force: **this is the ONE attempt.**
Why Now         `#902` established LINK-FREE that the four series are not one tide; `#901` said the
                residual dimension is same-sex-shaped; the raw pairwise ranks said extramarital.
                Nothing has yet been able to say WHICH ACT, and this is the last design available.
Live Worlds     EXTRA  only `xmarsex` departs -- adultery is judged on its own clock.
                SAME   only `homosex` departs, consistent with `#901`.
                BOTH   both, and the era has two special acts.
                NONE ⚠ THE UNWELCOME ONE -- no item departs once the contamination is removed, which
                       would mean `#901`'s shape and `#902`'s tension were both artifacts of
                       statistics that mix items together, and the arc ends with "not one tide" and
                       nothing about which act.
Discriminating  LEAVE-ONE-OUT. For held-out item i, build the common order from the OTHER THREE
Act             ONLY and ask how far i departs from it. Item i's own series never enters its own
                reference, which is precisely the contamination `#903` diagnosed.
                ⚠ AND THE REFERENCE IS THE **MEDIAN** RANK OF THE OTHER THREE, not the mean: with
                three items the median is robust to ONE deviant member, so a single departing act
                cannot drag the reference it is being compared against. Mean-rank is kept as a
                specification axis so the choice is visible rather than assumed.
Prediction      EXTRA -> only `xmarsex`'s |rho| below its own null
Matrix          SAME  -> only `homosex`'s
                BOTH  -> both, neither of the other two
                NONE  -> none
Confound        ⚠ WRITTEN BEFORE THE RUN: leave-one-out is LESS contaminated, not UNcontaminated --
                a departing act still sits in the other three items' references as 1 of 3. The
                median reference bounds that to "cannot move the median unless two of three move
                together", which is stated as a bound and not as elimination.
Controls        positive: plant a departure into ONE named item; only THAT item's statistic must
                move, and at g=0 nothing must · negative: the comonotone world, per item, where the
                held-out item IS comonotone with the other three so the deficit is exactly 0
Stopping Rule   `#111c`: this is the one attempt. If it returns UNVERIFIED, A108 CLOSES with
                "which act" unanswered and the page says so. Budget: one round.
Cost            21 waves x 4 items, ~2,500 resamples. CPU seconds.
Priority        It is the only design left, and #903 wrote it down before it could see the answer.
Expected        If NONE: two of my own recent readings shrink to "the system is not one tide".
Transform
```

⚠⚠ **`#901`①'s REMEDY, THIRD USE.** The outcome is **which subset of four items departs — 16
subsets.** Four are assigned to worlds above (`{extra}`·`{same}`·`{extra,same}`·`∅`). **The other
twelve are assigned here, before the run, to: report verbatim, return `UNVERIFIED` on the world
assignment, and do NOT adopt the nearest world.** It has changed the outcome on both previous uses.

`G1` **ESTIMAND**: for each item `i`, **`|ρ_i|` = |Spearman| between item `i`'s series across waves
and the common order built from the MEDIAN rank of the other three**, and its **deficit** against
the comonotone world's value for that same held-out item. **Population** GSS respondents on the 21
waves where all four items were asked, 1988–2024, per-wave n 868–2,680. **Instrument** GSS
`gss7224_r3a` — ⚠ **one instrument, `no second instrument`, `only this one instrument`** (`#897`:
SCCS's matching design resolves no effect at all; `#891`: only GSS ships question text); mode changed
mid-series, a specification axis. **Baseline** the comonotone world, per held-out item. **Regime**
T = 21.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES.** In a comonotone world the held-out item IS comonotone with
the other three, so its deficit has expectation exactly 0. ⇒ **`negative_control`**, **kind of null
named: a COMONOTONE BINOMIAL RESAMPLING NULL, PER HELD-OUT ITEM** — every item's own observed values
reassigned to waves in one common order (link-free by construction, marginals preserved), binomial
noise at that wave's actual n, and the leave-one-out statistic recomputed for the same held-out item.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a departure planted in ONE named item raises THAT item's deficit inside
                           the MEASURED monotone region, leaves the other three inside their nulls,
                           and does NOT fire at g=0):
       the Holm-corrected subset is one of the four assigned  -> that world
       it is any of the other twelve                          -> UNVERIFIED, subset reported verbatim
else:
       UNVERIFIED
```
⚠ **The positive control here is TWO-SIDED and that is the point**: it must show not only that the
planted item moves but that **the other three do not.** `#903` failed for exactly the second half.
`G3`: Holm over the family of four, all four reported. `G4`: {4 thresholds} × {mean · median
reference} × {all 21 waves · 1988–2018 only}.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **leave-one-out is LESS contaminated, not UNcontaminated** — a departing act is still 1 of 3 in
   the others' references, bounded by the median but not eliminated;
② **a departure says an act moves on its own clock, never WHY** — "betrayal", "consent", "visibility"
   are readings, not measurements;
③ **the comonotone null assumes a COMMON link across items** (`#902`①'s standing limit);
④ **T = 21**, MDE not computed (`#898`'s standing debt);
⑤ **mode is confounded with period**; **cohort is not separated from period**;
⑥ **cross-instrument N/A**; ⑦ ⚠ **`[unchallenged]`** — `door ③`; ⑧ no second coder, no second
   release, no test–retest.
"""
import itertools
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
from lib.gates import Gate

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(342)
NSIM = 1500
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
SHORT = {"premarsx": "pre", "teensex": "teen", "xmarsex": "extra", "homosex": "same"}
MODE_CHANGED = [2021, 2022, 2024]
PLANT_ITEM = "xmarsex"          # named BEFORE the run
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

print("=== (0) HARD RULE 1 — n, the years actually asked, and the VALUE SET ===")
d = pd.read_stata(F, columns=["year"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])} ({len(ys):2d} waves)  "
          f"codes={[int(v) for v in sorted(s[c].unique())]}")
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
WAVES = [int(y) for y in sorted(set.intersection(
    *[set(d[[c, "year"]].dropna().year.unique()) for c in ITEMS]))]
sub = d[d.year.isin(WAVES)]
print(f"  common waves {len(WAVES)}")
if len(WAVES) < 8:
    raise SystemExit("STOP: too few common waves; an empty design must never pass")

THRESH = {"≤1 always wrong": 1, "≤2 always/almost": 2, "≤3 any wrongness": 3, "mean 1–4": None}


def series(thr, keep):
    ps, ns = {}, {}
    for c in ITEMS:
        s = sub[sub.year.isin(keep)][[c, "year"]].dropna()
        g = s.groupby("year")[c]
        ns[c] = g.size()
        ps[c] = g.apply(lambda v: float((v <= thr).mean())) if thr else (4 - g.mean()) / 3.0
    return pd.DataFrame(ps).loc[keep], pd.DataFrame(ns).loc[keep]


def loo_abs_rho(P, agg="median"):
    """|Spearman| between each item and the common order of the OTHER THREE.

    ⚠ MEDIAN rank, not mean: with three items the median is robust to ONE deviant member, so a
    single departing act cannot drag the reference it is being compared against. That is the whole
    repair `#903`① specified."""
    A = P.to_numpy(float)
    k = A.shape[1]
    rk = np.apply_along_axis(stats.rankdata, 0, A)
    out = np.empty(k)
    for i in range(k):
        others = [j for j in range(k) if j != i]
        # orient each other item with the group before aggregating, so a negatively-loading item
        # does not cancel the reference it belongs to
        ref0 = rk[:, others].mean(1)
        oriented = np.column_stack([
            rk[:, j] if (stats.spearmanr(rk[:, j], ref0).statistic or 1) >= 0 else -rk[:, j]
            for j in others])
        ref = np.median(oriented, axis=1) if agg == "median" else oriented.mean(1)
        r = stats.spearmanr(A[:, i], ref).statistic
        out[i] = abs(float(r)) if np.isfinite(r) else np.nan
    return out


def comonotone_probs(P):
    A = P.to_numpy(float)
    rk = np.apply_along_axis(stats.rankdata, 0, A)
    common = rk.mean(1)
    order = np.argsort(np.argsort(common))
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        sgn = np.sign(stats.spearmanr(A[:, j], common).statistic) or 1.0
        srt = np.sort(A[:, j]) if sgn > 0 else np.sort(A[:, j])[::-1]
        out[:, j] = srt[order]
    return out


def resample_loo(P, N, probs, nsim, agg="median"):
    nn = N.to_numpy(float)
    out = []
    for _ in range(nsim):
        k = RNG.binomial(nn.astype(int), np.clip(probs, 1e-4, 1 - 1e-4))
        v = loo_abs_rho(pd.DataFrame(k / nn, index=P.index, columns=P.columns), agg)
        if np.isfinite(v).all():
            out.append(v)
    return np.vstack(out)


print("\n=== (1) LEAVE-ONE-OUT: each item against the MEDIAN order of the other three ===")
P0, N0 = series(2, WAVES)
OBSV = loo_abs_rho(P0)
CP = comonotone_probs(P0)
NULLV = resample_loo(P0, N0, CP, NSIM)
DEF = NULLV.mean(0) - OBSV
P95 = np.percentile(NULLV.mean(0)[None, :] - NULLV, 95, axis=0)
PVAL = np.array([float((((NULLV.mean(0)[i] - NULLV[:, i])) >= DEF[i]).mean()) for i in range(4)])
print(f"  {'held-out':10s} {'|rho| LOO':>10s} {'comonotone':>11s} {'deficit':>9s} {'null 95th':>10s} {'p':>7s}")
for i, c in enumerate(ITEMS):
    print(f"  {SHORT[c]:10s} {OBSV[i]:10.4f} {NULLV.mean(0)[i]:11.4f} {DEF[i]:+9.4f} "
          f"{P95[i]:10.4f} {PVAL[i]:7.4f}")
order = np.argsort(PVAL)
holm = np.zeros(4, bool)
for rank, i in enumerate(order):
    if PVAL[i] <= 0.05 / (4 - rank):
        holm[i] = True
    else:
        break
SUBSET = frozenset(SHORT[c] for i, c in enumerate(ITEMS) if holm[i])
print(f"\n  Holm over the family of FOUR (all four reported whatever they do): "
      f"**{{{', '.join(sorted(SUBSET)) or '∅'}}}**")

print("\n=== (2) POSITIVE CONTROL — TWO-SIDED: the planted item must move AND the others must not ===")
j = ITEMS.index(PLANT_ITEM)
w2 = RNG.standard_normal(len(WAVES)); w2 -= w2.mean(); w2 /= w2.std()
zc = stats.norm.ppf(np.clip(CP, 1e-4, 1 - 1e-4))
sweep, others_max = [], []
for g in (0.0, 0.05, 0.10, 0.20, 0.35, 0.60):
    z = zc.copy()
    z[:, j] += g * float(zc.std()) * w2
    V = resample_loo(P0, N0, stats.norm.cdf(z), 250)
    dd = NULLV.mean(0)[None, :] - V
    sweep.append((g, float(np.median(dd[:, j]))))
    om = max(float(np.median(dd[:, i])) for i in range(4) if i != j)
    others_max.append(om)
    print(f"  g={g:<5.2f} deficit({SHORT[PLANT_ITEM]}) {sweep[-1][1]:+.4f}   "
          f"largest OTHER deficit {om:+.4f}  (its null 95th ≈ {P95[[i for i in range(4) if i!=j]].max():.4f})")
TURN = int(np.argmax([v for _, v in sweep]))
mono = sweep[:TURN + 1]
OTHERS_CLEAN = all(others_max[k] <= P95[[i for i in range(4) if i != j]].max()
                   for k in range(TURN + 1))
PC_OK = (abs(mono[0][1]) < 3 * float(NULLV.std(0)[j])) and (mono[-1][1] > P95[j]) and (TURN >= 2) \
        and all(mono[i][1] <= mono[i + 1][1] + 1e-9 for i in range(len(mono) - 1)) and OTHERS_CLEAN
print(f"  ⚠ turning point MEASURED at g={sweep[TURN][0]:g} · monotone in [0, {sweep[TURN][0]:g}] · "
      f"g=0 {mono[0][1]:+.4f} lands on zero ⇒ CAN fail")
print(f"  ⚠ **THE SECOND HALF, which `#903` failed: do the OTHER THREE stay inside their nulls "
      f"while `{SHORT[PLANT_ITEM]}` is planted? -> {OTHERS_CLEAN}**")
print(f"  ⇒ positive control fires: {PC_OK}")

print("\n=== (3) G3/G4 — 4 thresholds × {median · mean reference} × {all 21 waves · 1988–2018} ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for tname, thr in THRESH.items():
    for agg in ("median", "mean"):
        for sname, keep in (("all", WAVES), ("pre-2021", PRE)):
            Pk, Nk = series(thr, keep)
            ov = loo_abs_rho(Pk, agg)
            nv = resample_loo(Pk, Nk, comonotone_probs(Pk), 300, agg)
            dk = nv.mean(0) - ov
            pk = np.array([float(((nv.mean(0)[i] - nv[:, i]) >= dk[i]).mean()) for i in range(4)])
            o = np.argsort(pk); hh = np.zeros(4, bool)
            for r_, i_ in enumerate(o):
                if pk[i_] <= 0.05 / (4 - r_):
                    hh[i_] = True
                else:
                    break
            ss = tuple(sorted(SHORT[c] for i_, c in enumerate(ITEMS) if hh[i_]))
            rows.append((tname, agg, sname, [float(x) for x in dk], ss))
            print(f"  {tname:18s} {agg:6s} {sname:8s} " +
                  " ".join(f"{SHORT[c]}={dk[i]:+.3f}" for i, c in enumerate(ITEMS)) +
                  f"  ⇒ {{{', '.join(ss) or '∅'}}}")
subs = {}
for r in rows:
    subs[r[4]] = subs.get(r[4], 0) + 1
TOP, TOPN = max(subs.items(), key=lambda x: x[1])
print(f"\n  **grid: {TOPN}/{len(rows)} cells return {{{', '.join(TOP) or '∅'}}}** · "
      f"{len(subs)} distinct subsets across the grid")

print("\n=== (4) THE CONDITIONAL KILL — 16 subsets, all assigned before the run ===")
ASSIGNED = {frozenset({"extra"}): "EXTRA", frozenset({"same"}): "SAME",
            frozenset({"extra", "same"}): "BOTH", frozenset(): "NONE"}
G = Gate("Leave-one-out: which act departs from the common order?")
G.plant_direction_from_sweep(f"positive: a planted departure raises `{SHORT[PLANT_ITEM]}`'s LOO "
                             f"deficit (inside the MEASURED monotone region)", mono, baseline=0.0,
                             baseline_spread=float(NULLV.std(0)[j]), half_of=max(P95[j], 1e-4))
G.negative_control("comonotone world, per held-out item", 0.0, float(np.max(np.abs(DEF))),
                   null_spread=float(NULLV.std(0).max()),
                   null_kind="COMONOTONE BINOMIAL RESAMPLING NULL, PER HELD-OUT ITEM — every item's "
                             "own observed values reassigned to waves in one common order "
                             "(link-free, marginals preserved), binomial noise at that wave's "
                             "actual n, leave-one-out statistic recomputed for the same held-out item")
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", ("the positive control did not license a reading"
                                    + ("" if OTHERS_CLEAN else " — and it failed on its SECOND half: "
                                       "planting one item moved the others too, so leave-one-out did "
                                       "NOT remove the contamination `#903` diagnosed"))
elif SUBSET in ASSIGNED:
    W = ASSIGNED[SUBSET]
    VERDICT = "OVERTURNED" if W == "NONE" else "CONFIRMED"
    WORLD = f"{W} — subset {{{', '.join(sorted(SUBSET)) or '∅'}}}"
else:
    VERDICT, WORLD = "UNVERIFIED", (f"UNLISTED SUBSET {{{', '.join(sorted(SUBSET))}}} — no world "
                                    f"predicted it; reported verbatim, nearest world NOT adopted")
print(G)
print(f"\n  subset {{{', '.join(sorted(SUBSET)) or '∅'}}} · grid {TOPN}/{len(rows)} · "
      f"others-stay-clean under the plant: {OTHERS_CLEAN}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ `#111c`: `#903` was the first UNVERIFIED on 'which act is special'. If this is the")
print("     second, A108 CLOSES with the question unanswered rather than buying a third round.")

art = dict(entry=904, round="E03·A108·R342", verdict=VERDICT, world=WORLD, waves=WAVES,
           loo_abs_rho={c: float(OBSV[i]) for i, c in enumerate(ITEMS)},
           comonotone={c: float(NULLV.mean(0)[i]) for i, c in enumerate(ITEMS)},
           deficit={c: float(DEF[i]) for i, c in enumerate(ITEMS)},
           null_p95={c: float(P95[i]) for i, c in enumerate(ITEMS)},
           pvalues={c: float(PVAL[i]) for i, c in enumerate(ITEMS)},
           holm_subset=sorted(SUBSET), plant_item=PLANT_ITEM,
           positive_sweep=sweep, others_max=others_max, others_clean=bool(OTHERS_CLEAN),
           turning_point=sweep[TURN][0], positive_ok=bool(PC_OK),
           grid_rows=rows, grid_top=list(TOP), grid_top_n=TOPN, grid_distinct=len(subs),
           unchallenged=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "leave_one_out.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'leave_one_out.json'}")
