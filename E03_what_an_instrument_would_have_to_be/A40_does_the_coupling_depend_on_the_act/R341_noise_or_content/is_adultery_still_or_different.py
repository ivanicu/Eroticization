r"""#903 · E03·A108·R341 — is adultery a different KIND of judgement, or just one that has not moved?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#902` refuted the one-factor family LINK-FREE (min pairwise |Spearman| 0.4143 vs a
                comonotone null of 0.7883, 6.1x its spread) and then recorded a tension it REFUSED
                to claim: the raw pairwise ranks make EXTRAMARITAL the odd one out (its two lowest
                pairs) while `#901`'s residual dimension isolated SAME-SEX. `#902`① says the
                reconciliation was "an ARGUMENT, not a measurement", and that inventing the estimand
                inside the round that noticed the tension would be a threshold chosen after the
                result. **This is that estimand, pre-registered in a round of its own.**
Why Now         It is the only measurable question A107 left, and it is the one with the psychology:
                extramarital sex is the single act Americans did NOT liberalise on (92.3% -> 85.2%
                over 36 years). Whether that is a different KIND of moral judgement or merely a
                SLOWER one is a question about people, and the two answers imply different objects.
Live Worlds     NOISE-ONLY ⚠ THE UNWELCOME ONE -- every item's departure sits inside its own null;
                           `xmarsex` simply barely moved, and `#901`'s same-sex pattern does not
                           reappear in this more direct statistic, which WEAKENS `#901`.
                EXTRA      `xmarsex` is the outlier: adultery is judged on a different axis --
                           betrayal of a person rather than disapproval of an act.
                SAME       `homosex` is the outlier, consistent with `#901`.
                BOTH       both depart, and the era has two special acts rather than one.
Discriminating  Per item, its MEAN |Spearman| with the other three, minus what the COMONOTONE world
Act             gives THAT ITEM. The comonotone null preserves each item's own observed values, so
                it automatically prices "small movement -> noisy rank order" -- which is exactly the
                deflationary rival, priced rather than argued.
Prediction      NOISE-ONLY -> no item's deficit beyond its own null
Matrix          EXTRA      -> only `xmarsex`
                SAME       -> only `homosex`
                BOTH       -> both, and neither of the other two
Confound        ⚠ WRITTEN BEFORE THE RUN. The deficit is a difference of two |Spearman| averages and
                is bounded, so a maximal plant SATURATES and the positive control will be
                non-monotone past a turning point -- the `#885`/`#887`/`#902` family. The turning
                point is MEASURED and monotonicity required only inside it, as `#902` did.
Controls        positive: plant a graded item-specific departure into ONE named item; its deficit
                must rise inside the measured monotone region and must NOT fire at g=0 · offset: the
                comonotone binomial resampling null, per item
Stopping Rule   Whatever fires, A108 is one round. If NOISE-ONLY, `#901` is downgraded and the arc
                closes with the tension resolved AGAINST my newer entry.
Cost            21 waves x 4 items, ~3,000 resamples. CPU seconds.
Priority        `#902`② closed A107; this is the one question it left that is not Closure.
Expected        If EXTRA: the project gains a second CONTENT claim -- some acts are judged as harms
Transform       to a person and move on their own clock.
```

⚠⚠ **`#901`①'s REMEDY, SECOND USE, AND THE OUTCOME SPACE IS BIGGER THAN THE WORLD LIST.** The
outcome is **which subset of the four items has a deficit beyond its own null — 16 subsets.** Four
are assigned to worlds above (`∅`→NOISE-ONLY · `{extra}`→EXTRA · `{same}`→SAME · `{extra,same}`→BOTH).
**The other TWELVE are UNLISTED and are assigned here, before the run, to a single explicit
outcome:** *no world predicted this subset; report it verbatim and return `UNVERIFIED` on the world
assignment rather than adopting whichever world is nearest.* **Last time reality returned a fifth
world out of a space I had not written down; this time all sixteen cells have a destination.**

`G1` **ESTIMAND**: for each item, **`deficit_i` = (its mean |Spearman| with the other three under the
comonotone null) − (its observed mean |Spearman| with the other three)**. **Population** GSS
respondents on the 21 waves where all four items were asked, 1988–2024, per-wave n 868–2,680.
**Instrument** GSS `gss7224_r3a` — ⚠ **one instrument, `no second instrument`, `only this one
instrument`** (`#897`: SCCS's matching design resolves no effect at all; `#891`: only GSS ships
question text) — mode changed mid-series, a specification axis. **Baseline** the comonotone world's
own per-item distribution. **Regime** T = 21.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES.** Under NOISE-ONLY each item's observed mean |Spearman| is
drawn from the comonotone world's distribution *for that item*, so the deficit has expectation
exactly 0. ⇒ **`negative_control`**, **kind of null named: a COMONOTONE BINOMIAL RESAMPLING NULL,
PER ITEM — each item's own observed values reassigned to waves in one common order (link-free by
construction, every marginal preserved), then binomial noise at that wave's actual n, and the
statistic recomputed for the same item.**

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a planted item-specific departure raises that item's deficit inside the
                           MEASURED monotone region, and does NOT fire at g=0):
       the subset of items with deficit > its own null's 95th percentile, family-wise corrected,
       is one of the four assigned subsets  -> that world
       it is any of the other twelve        -> UNVERIFIED on the world assignment, subset reported
else:
       UNVERIFIED
```
`G3`: the family is 4 items; Holm correction over the four, and **all four reported whatever they do**.
`G4`: {4 condemnation thresholds} × {all 21 waves · 1988–2018 only}.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **a deficit says an item departs from comonotonicity, never WHY** — "betrayal" is a reading, not a
   measurement, and this design cannot distinguish it from any other reason `xmarsex` might move on
   its own clock;
② **the comonotone null assumes a COMMON link across items** (`#902`①'s standing limit);
③ **T = 21** and the MDE of the per-item deficit is **not computed** (`#898`'s standing debt);
④ **mode is confounded with period**; **cohort is not separated from period**;
⑤ **cross-instrument N/A**; ⑥ no second coder, no second release, no test–retest;
⑦ ⚠ **`[unchallenged]`** — `door ③`, and `#899`'s pre-registration table is what a real adversary
   should be scored against.
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
RNG = np.random.default_rng(341)
NSIM = 1500
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
SHORT = {"premarsx": "pre", "teensex": "teen", "xmarsex": "extra", "homosex": "same"}
MODE_CHANGED = [2021, 2022, 2024]
PLANT_ITEM = "xmarsex"          # named BEFORE the run; the control tests the instrument, not the world
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R340 = (ROOT / "E03_what_an_instrument_would_have_to_be/A42_one_tide_or_four_histories/"
        "R340_was_the_second_dimension_just_the_wrong_link/results/link_free_attack.json")

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
PRIOR = json.loads(R340.read_text())
print(f"  common waves {len(WAVES)} · `#902` read from its artifact (`#840`'s RULE; its own scope "
      f"was the `homosex` item alone, so only the practice transfers): "
      f"min|ρ| {PRIOR['min_abs_spearman']:.4f} vs null 5th {PRIOR['null_p05']:.4f}")
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


def mean_abs_rho(P):
    """Per item: its mean |Spearman| with the other three."""
    A = P.to_numpy(float)
    k = A.shape[1]
    R = np.eye(k)
    for i, j in itertools.combinations(range(k), 2):
        r = stats.spearmanr(A[:, i], A[:, j]).statistic
        if not np.isfinite(r):
            return None
        R[i, j] = R[j, i] = abs(float(r))
    return np.array([R[i, [x for x in range(k) if x != i]].mean() for i in range(k)])


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


def resample_mean_rho(P, N, probs, nsim):
    nn = N.to_numpy(float)
    out = []
    for _ in range(nsim):
        k = RNG.binomial(nn.astype(int), np.clip(probs, 1e-4, 1 - 1e-4))
        v = mean_abs_rho(pd.DataFrame(k / nn, index=P.index, columns=P.columns))
        if v is not None:
            out.append(v)
    return np.vstack(out)


print("\n=== (1) PER-ITEM MEAN |Spearman| WITH THE OTHER THREE, and its comonotone expectation ===")
P0, N0 = series(2, WAVES)
OBSV = mean_abs_rho(P0)
CP = comonotone_probs(P0)
NULLV = resample_mean_rho(P0, N0, CP, NSIM)
DEF = NULLV.mean(0) - OBSV
P95 = np.percentile(NULLV.mean(0)[None, :] - NULLV, 95, axis=0)
PVAL = np.array([float((( NULLV.mean(0)[i] - NULLV[:, i]) >= DEF[i]).mean()) for i in range(4)])
print(f"  {'item':10s} {'observed':>9s} {'comonotone':>11s} {'deficit':>9s} {'null 95th':>10s} "
      f"{'p':>7s}")
for i, c in enumerate(ITEMS):
    print(f"  {SHORT[c]:10s} {OBSV[i]:9.4f} {NULLV.mean(0)[i]:11.4f} {DEF[i]:+9.4f} "
          f"{P95[i]:10.4f} {PVAL[i]:7.4f}")
order = np.argsort(PVAL)
holm = np.zeros(4, bool)
for rank, i in enumerate(order):
    if PVAL[i] <= 0.05 / (4 - rank):
        holm[i] = True
    else:
        break
SUBSET = frozenset(SHORT[c] for i, c in enumerate(ITEMS) if holm[i])
print(f"\n  ⚠ Holm over the family of FOUR (all four reported whatever they do): "
      f"**{{{', '.join(sorted(SUBSET)) or '∅'}}}**")

print("\n=== (2) POSITIVE CONTROL — plant an item-specific departure into `xmarsex`, named before ===")
j = ITEMS.index(PLANT_ITEM)
w2 = RNG.standard_normal(len(WAVES)); w2 -= w2.mean(); w2 /= w2.std()
sweep, sw_sd = [], []
for g in (0.0, 0.05, 0.10, 0.20, 0.35, 0.60):
    pg = CP.copy()
    z = stats.norm.ppf(np.clip(CP, 1e-4, 1 - 1e-4))
    z[:, j] = z[:, j] + g * float(z.std()) * w2
    pg = stats.norm.cdf(z)
    V = resample_mean_rho(P0, N0, pg, 250)
    dd = NULLV.mean(0)[j] - V[:, j]
    sweep.append((g, float(np.median(dd))))
    sw_sd.append(float(dd.std(ddof=1)))
for (g, v), s in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} deficit({SHORT[PLANT_ITEM]}) {v:+.4f} ± {s:.4f}")
TURN = int(np.argmax([v for _, v in sweep]))
mono = sweep[:TURN + 1]
PC_OK = (abs(mono[0][1]) < 3 * sw_sd[0]) and (mono[-1][1] > P95[j]) and (TURN >= 2) and \
        all(mono[i][1] <= mono[i + 1][1] + 1e-9 for i in range(len(mono) - 1))
print(f"  ⚠ turning point MEASURED at g={sweep[TURN][0]:g} (the deficit is BOUNDED, so a maximal "
      f"plant saturates — the `#885`/`#887`/`#902` family, remedied by measuring rather than excusing)")
print(f"  monotone region g ∈ [0, {sweep[TURN][0]:g}] · g=0 {mono[0][1]:+.4f} lands on zero ⇒ the "
      f"control CAN fail · top {mono[-1][1]:+.4f} · threshold {P95[j]:.4f} strictly between: {PC_OK}")

print("\n=== (3) G3/G4 — 4 thresholds × {all 21 waves · 1988–2018 only}; every item reported ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for tname, thr in THRESH.items():
    for sname, keep in (("all waves", WAVES), ("pre-2021 only", PRE)):
        Pk, Nk = series(thr, keep)
        ov = mean_abs_rho(Pk)
        nv = resample_mean_rho(Pk, Nk, comonotone_probs(Pk), 400)
        dk = nv.mean(0) - ov
        pk = np.array([float(((nv.mean(0)[i] - nv[:, i]) >= dk[i]).mean()) for i in range(4)])
        o = np.argsort(pk); hh = np.zeros(4, bool)
        for r_, i_ in enumerate(o):
            if pk[i_] <= 0.05 / (4 - r_):
                hh[i_] = True
            else:
                break
        ss = frozenset(SHORT[c] for i_, c in enumerate(ITEMS) if hh[i_])
        rows.append((tname, sname, len(keep), [float(x) for x in dk], [float(x) for x in pk],
                     sorted(ss)))
        print(f"  {tname:18s} {sname:14s} deficits " +
              " ".join(f"{SHORT[c]}={dk[i]:+.3f}" for i, c in enumerate(ITEMS)) +
              f"  ⇒ {{{', '.join(sorted(ss)) or '∅'}}}")
subs = {}
for r in rows:
    subs[tuple(r[5])] = subs.get(tuple(r[5]), 0) + 1
TOP, TOPN = max(subs.items(), key=lambda x: x[1])
print(f"\n  **grid: {TOPN}/{len(rows)} cells return the same subset {{{', '.join(TOP) or '∅'}}}**")

print("\n=== (4) THE CONDITIONAL KILL — all SIXTEEN subsets had a destination before the run ===")
ASSIGNED = {frozenset(): "NOISE-ONLY", frozenset({"extra"}): "EXTRA",
            frozenset({"same"}): "SAME", frozenset({"extra", "same"}): "BOTH"}
G = Gate("Is adultery a different KIND of judgement, or just one that has not moved?")
G.plant_direction_from_sweep(f"positive: a planted departure raises `{SHORT[PLANT_ITEM]}`'s deficit "
                             f"(inside the MEASURED monotone region)", mono, baseline=0.0,
                             baseline_spread=sw_sd[0], half_of=max(P95[j], 1e-4))
G.negative_control("comonotone world, per item", 0.0, float(np.max(np.abs(DEF))),
                   null_spread=float(NULLV.std(0).max()),
                   null_kind="COMONOTONE BINOMIAL RESAMPLING NULL, PER ITEM — each item's own "
                             "observed values reassigned to waves in one common order (link-free by "
                             "construction, every marginal preserved), binomial noise at that "
                             "wave's actual n, statistic recomputed for the same item")
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif SUBSET in ASSIGNED:
    W = ASSIGNED[SUBSET]
    VERDICT = "CONFIRMED" if W != "NOISE-ONLY" else "OVERTURNED"
    WORLD = (f"{W} — subset {{{', '.join(sorted(SUBSET)) or '∅'}}}"
             + ("; `#901`'s same-sex pattern does NOT reappear in this more direct statistic, and "
                "`#901` is DOWNGRADED" if W == "NOISE-ONLY" else ""))
else:
    VERDICT = "UNVERIFIED"
    WORLD = (f"UNLISTED SUBSET {{{', '.join(sorted(SUBSET))}}} — no world predicted it; reported "
             f"verbatim, and the nearest world is NOT adopted (`#901`①, second use)")
print(G)
print(f"\n  subset {{{', '.join(sorted(SUBSET)) or '∅'}}} · grid agreement {TOPN}/{len(rows)}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ A deficit says an item DEPARTS from comonotonicity, never WHY. 'Betrayal' is a reading,")
print("     not a measurement, and this design cannot distinguish it from any other reason an act")
print("     might move on its own clock. `[unchallenged]` — `door ③`.")

art = dict(entry=903, round="E03·A108·R341", verdict=VERDICT, world=WORLD, waves=WAVES,
           observed_mean_abs_rho={c: float(OBSV[i]) for i, c in enumerate(ITEMS)},
           comonotone_mean={c: float(NULLV.mean(0)[i]) for i, c in enumerate(ITEMS)},
           deficit={c: float(DEF[i]) for i, c in enumerate(ITEMS)},
           null_p95={c: float(P95[i]) for i, c in enumerate(ITEMS)},
           pvalues={c: float(PVAL[i]) for i, c in enumerate(ITEMS)},
           holm_subset=sorted(SUBSET), plant_item=PLANT_ITEM,
           positive_sweep=sweep, positive_sd=sw_sd, turning_point=sweep[TURN][0],
           positive_ok=bool(PC_OK), grid_rows=rows, grid_top=list(TOP), grid_top_n=TOPN,
           outcome_space="16 subsets; 4 assigned to worlds, the other 12 assigned BEFORE the run to "
                         "'report verbatim and return UNVERIFIED on the world assignment'",
           unchallenged=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "noise_or_content.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'noise_or_content.json'}")
