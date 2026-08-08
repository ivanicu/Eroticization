r"""#900 · E03·A107·R338 — one tide, or four histories? The four sexual norms, 21 waves, rank-tested

**COGNITIVE UPDATE CARD**
```
Core Gap        A105 established that moral judgement is CASE-indexed WITHIN one act (abortion):
                `D = +0.6688` at 5.6x its one-factor null. Whether that generalises to the project's
                actual subject -- sex -- and to the DECADE unit has never been asked. Evidence
                cannot currently distinguish "America has ONE sexual tide that four acts ride at
                different rates" from "four acts with four separate histories".
Why Now         It is the first question in this project that uses the decade as the unit AND the
                subject matter the project is named for, and `#899`② cleared `#887`①'s ban for it
                by CHECKING the ban's scope rather than assuming it.
Live Worlds     TIDE  one moving permissiveness factor with FIXED loadings. Then
                      z_i(t) = b_i - lambda_i * F(t), the item-centred T x 4 probit matrix has rank
                      EXACTLY 1, and the four series are one number per era plus four constants.
                FOUR  the acts have separate histories; a second dimension is required.
                MODE  ⚠ META -- the second dimension is the INSTRUMENT, not the object. GSS changed
                      mode in 2021 (web/mail push) and 2022/2024 are mixed-mode; a mode shift
                      injects exactly the structure FOUR predicts.
                ⚠ UNWELCOME BRANCH IS **TIDE**, and that is why this step was chosen: my own A105
                      result predicts FOUR, so a rank-1 world would say case-indexing is a
                      within-person within-act phenomenon that does NOT aggregate, narrowing the
                      only surviving finding in E03.
Discriminating  PC1's variance share of the item-centred probit matrix, against a null in which the
Act             truth IS rank 1 and only binomial sampling at the observed per-wave n's is added.
Prediction      TIDE -> PC1 share inside the rank-1 null
Matrix          FOUR -> PC1 share below the null's 5th percentile, AND the deficit survives dropping
                        the mode-changed waves
                MODE -> PC1 share below the null on all 21 waves, INSIDE it on 1988-2018
Confound        ⚠ WRITTEN BEFORE THE RUN, and it is not the mode change -- it is POWER IN THE OTHER
                DIRECTION. Per-wave n is 868-2,680, so a probit cell has SE ~0.037 against a signal
                range ~1.2: signal/noise ~32. **At this precision the rank test will reject rank 1
                for a second dimension far too small to matter to anyone.** ⇒ a rejection answers
                "is there ANY second dimension", NEVER "does it matter" -- `#898`'s error class,
                pre-empted. So the round reports BOTH the test AND the SIZE, in percentage points.
Controls        positive: a graded second factor planted into the rank-1 world; PC1 share must fall
                monotonically and must NOT fire at g=0 · offset: the rank-1 binomial resampling null
Stopping Rule   Whatever the verdict, the size is reported in points and the sentence about people
                is written from the SIZE, not from the test. Budget: one round.
Cost            21 waves x 4 items, ~4,000 simulated matrices. CPU seconds.
Priority        It is the only open question that is about people, about sex, and at a unit this
                project has never used.
Expected        If TIDE: E03's finding shrinks to within-person and the epoch's object was never
Transform       the era. If FOUR: the era has more than one dimension and "permissiveness" as a
                single summary is wrong at the population level.
```

⚠⚠ **PRIOR ART, DECLARED BEFORE THE RESULT.** **That these four series DIVERGED is textbook** —
Twenge, Sherman & Wells (2015, *Archives of Sexual Behavior*) documented premarital and same-sex
attitudes liberalising while extramarital hardened. `D6`, from reading, not from a search run in
this session. **Therefore: reporting divergence would be a VERIFICATION, not a discovery, and this
round does not claim it.** What the literature does **not** settle is the question here: **a
one-factor model with different loadings PREDICTS different rates of change**, so "they moved at
different rates" is *compatible with one factor* and is not evidence against it. **The rank test is
the part that is not in the descriptive literature, and it is the only part this round may claim.**

`G1` **ESTIMAND**: **PC1's share of the variance of the item-mean-centred `T x 4` probit matrix**,
plus **the reconstruction error in percentage points** when each series is rebuilt from PC1 alone.
**Population** GSS respondents on the **21 waves where all four items were asked** (1988–2024).
**Instrument** GSS `gss7224_r3a`, one questionnaire — ⚠ **and it CHANGED MID-SERIES**: 2021 was a
web/mail push and 2022/2024 are mixed-mode, which is a HARD RULE 2 instrument threat and is a
specification axis below, not a footnote. **Baseline** a rank-1 world at the observed per-wave n's.
**Regime** per-wave n 868–2,680 per item.

⚠ **"SHOULD THIS ZERO BE ZERO?" — NO.** Under a TRUE rank-1 world PC1's share is **not** 1: binomial
sampling at finite n pushes it down by a systematic, positive amount. ⇒ **`offset_control`**, and
**the kind of null is named: a RANK-1 BINOMIAL RESAMPLING NULL at the observed per-wave n's** — the
rank-1 approximation is fitted, every cell is redrawn as a binomial at that wave's actual n, and
PC1's share is recomputed. *The deficit that null produces is the part of any observed deficit that
costs nothing.*

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a planted second factor drives PC1's share DOWN monotonically, with
                           floor and ceiling MEASURED, and it does NOT fire at g=0):
       PC1 share >= the rank-1 null's 5th percentile                       -> TIDE
       PC1 share <  it on all 21 waves AND on 1988-2018 alone              -> FOUR
       PC1 share <  it on all 21 waves but INSIDE it on 1988-2018          -> MODE
else:
       UNVERIFIED
```
**AND THE SIZE IS REPORTED WHATEVER THE TEST SAYS**, because a rank test at n~1,500 answers a
different question than the sentence anyone wants to write.

`G3`/`G4`: {4 condemnation thresholds: `<=1` · `<=2` · `<=3` · the 1–4 mean} × {all 21 waves ·
1988–2018 only} × {PC1 share · reconstruction error}. Whole grid published.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **21 waves is 21 rows** — a rank test on a `21 x 4` matrix has limited power against a *weak*
   second dimension, and its MDE is not computed here (`#898`'s debt, one level over);
② **mode is confounded with period** — 2021+ is both a mode change and the most recent era, and
   nothing in this release separates them; the axis is reported, not resolved;
③ **`homosex` carries a fifth code** ("other") the other three do not; it is dropped and counted,
   and the asymmetry between items is registered rather than smoothed;
④ **cross-instrument N/A — `no second instrument` and `only this one instrument`**: `#897` measured
   that SCCS's matching design cannot resolve any effect at all, and `#891` that only GSS ships
   question text;
⑤ **causally identified N/A** — repeated cross-sections, no intervention on an era;
⑥ **cohort and period are not separated** — this is a period design and a cohort-replacement world
   would produce the same series; that is a different round and is not run here;
⑦ no second coder, no second release, no test–retest.
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
RNG = np.random.default_rng(338)
NSIM = 1000
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
LABEL = {"premarsx": "sex before marriage", "teensex": "sex at 14–16",
         "xmarsex": "sex outside marriage", "homosex": "same-sex relations"}
MODE_CHANGED = [2021, 2022, 2024]          # web/mail push 2021; mixed-mode after
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

print("=== (0) HARD RULE 1 — n, the years actually asked, and the VALUE SET, before any citation ===")
d = pd.read_stata(F, columns=["year"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])} ({len(ys):2d} waves)  "
          f"codes={[int(v) for v in sorted(s[c].unique())]}  ({LABEL[c]})")
n5 = int((d["homosex"] == 5).sum())
print(f"  ⚠ **`homosex` carries a FIFTH code ('other', n={n5}) the other three do not** — dropped "
      f"below, counted here, and the asymmetry is registered rather than smoothed")
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
waves = sorted(set.intersection(*[set(d[[c, "year"]].dropna().year.unique()) for c in ITEMS]))
WAVES = [int(y) for y in waves]
print(f"\n  common waves (all four asked): **{len(WAVES)}** — {WAVES}")
sub = d[d.year.isin(WAVES)]
NPW = {c: sub[[c, "year"]].dropna().groupby("year").size() for c in ITEMS}
print("  per-wave n: " + " · ".join(f"{c}={NPW[c].min()}–{NPW[c].max()}" for c in ITEMS))
if len(WAVES) < 8:
    raise SystemExit("STOP: too few common waves; an empty design must never pass")

THRESH = {"≤1 always wrong": 1, "≤2 always/almost": 2, "≤3 any wrongness": 3, "mean 1–4": None}


def series(thr, keep):
    """p_i(t) and n_i(t) on the kept waves. `thr=None` ⇒ the 1–4 mean, rescaled to (0,1)."""
    ps, ns = {}, {}
    for c in ITEMS:
        s = sub[sub.year.isin(keep)][[c, "year"]].dropna()
        g = s.groupby("year")[c]
        ns[c] = g.size()
        ps[c] = g.apply(lambda v: float((v <= thr).mean())) if thr else (4 - g.mean()) / 3.0
    return pd.DataFrame(ps).loc[keep], pd.DataFrame(ns).loc[keep]


def pc1_share(P):
    """PC1's variance share of the ITEM-MEAN-CENTRED probit matrix.

    Under one moving factor with fixed loadings, z_i(t) = b_i − λ_i·F(t) ⇒ the centred matrix has
    rank EXACTLY 1 for ANY loadings, so different rates of change are NOT evidence against it."""
    z = stats.norm.ppf(np.clip(P.to_numpy(float), 1e-4, 1 - 1e-4))
    z = z - z.mean(0, keepdims=True)
    if not np.isfinite(z).all() or z.std() == 0:
        return np.nan, None
    sv = np.linalg.svd(z, compute_uv=False)
    return float(sv[0] ** 2 / (sv ** 2).sum()), z


def rank1_fit(z):
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    return u[:, :1] * s[0] @ vt[:1, :]


print("\n=== (1) THE FOUR SERIES, and PC1's share of their centred probit matrix ===")
P0, N0 = series(2, WAVES)
for c in ITEMS:
    print(f"  {LABEL[c]:22s} {100*P0[c].iloc[0]:5.1f}% → {100*P0[c].iloc[-1]:5.1f}%  "
          f"({WAVES[0]}→{WAVES[-1]}) say always/almost-always wrong")
SHARE, Z = pc1_share(P0)
print(f"\n  **PC1 share (all {len(WAVES)} waves, threshold ≤2) = {SHARE:.5f}**  "
      f"⇒ PC2+ carries {100*(1-SHARE):.3f}% of the centred variance")

print("\n=== (2) OFFSET CONTROL — the RANK-1 BINOMIAL RESAMPLING NULL at the observed per-wave n ===")


def null_shares(P, N, nsim=NSIM):
    z = stats.norm.ppf(np.clip(P.to_numpy(float), 1e-4, 1 - 1e-4))
    z = z - z.mean(0, keepdims=True)
    fit = rank1_fit(z) + stats.norm.ppf(np.clip(P.to_numpy(float), 1e-4, 1 - 1e-4)).mean(0, keepdims=True)
    p1 = np.clip(stats.norm.cdf(fit), 1e-4, 1 - 1e-4)
    nn = N.to_numpy(float)
    out = []
    for _ in range(nsim):
        k = RNG.binomial(nn.astype(int), p1)
        Q = pd.DataFrame(k / nn, index=P.index, columns=P.columns)
        v, _ = pc1_share(Q)
        if np.isfinite(v):
            out.append(v)
    return np.asarray(out)


NULL = null_shares(P0, N0)
OFF, OSD, O05 = float(np.median(NULL)), float(NULL.std(ddof=1)), float(np.percentile(NULL, 5))
print(f"  rank-1 null: median {OFF:.5f} · sd {OSD:.5f} · **5th percentile {O05:.5f}**")
print(f"  ⇒ a TRUE rank-1 world does NOT return 1.0 here — it returns {OFF:.5f}, and that deficit "
      f"is free")
print(f"  observed {SHARE:.5f} vs the null's 5th percentile {O05:.5f} ⇒ "
      f"{'BELOW — rank 1 rejected' if SHARE < O05 else 'INSIDE — rank 1 survives'}")

print("\n=== (3) POSITIVE CONTROL — a graded SECOND factor; at g=0 it must NOT fire ===")
zc = stats.norm.ppf(np.clip(P0.to_numpy(float), 1e-4, 1 - 1e-4))
base = zc.mean(0, keepdims=True)
r1 = rank1_fit(zc - base)
w2 = RNG.standard_normal(len(WAVES))
w2 -= w2.mean()
l2 = np.array([1.0, -1.0, 1.0, -1.0])          # a contrast the first factor cannot express
scale = np.abs(r1).std()
sweep, sw_sd = [], []
for g in (0.0, 0.05, 0.10, 0.20, 0.35, 0.50):
    vals = []
    for _ in range(120):
        z2 = r1 + g * scale * np.outer(w2 / w2.std(), l2) + base
        p2 = np.clip(stats.norm.cdf(z2), 1e-4, 1 - 1e-4)
        k = RNG.binomial(N0.to_numpy(float).astype(int), p2)
        Q = pd.DataFrame(k / N0.to_numpy(float), index=P0.index, columns=P0.columns)
        v, _ = pc1_share(Q)
        if np.isfinite(v):
            vals.append(v)
    sweep.append((g, float(np.mean(vals))))
    sw_sd.append(float(np.std(vals, ddof=1)))
for (g, v), s in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} PC1 share {v:.5f} ± {s:.5f}")
PC_OK = (abs(sweep[0][1] - OFF) < 3 * OSD) and (sweep[-1][1] < O05) and \
        all(sweep[i][1] >= sweep[i + 1][1] - 1e-9 for i in range(len(sweep) - 1))
print(f"  floor(g=0) {sweep[0][1]:.5f} lands on the null median {OFF:.5f} ⇒ the control CAN fail · "
      f"ceiling(g=0.5) {sweep[-1][1]:.5f} · threshold {O05:.5f} strictly between: {PC_OK}")

print("\n=== (4) THE SIZE, IN PERCENTAGE POINTS — because a rank test answers a different question ===")
print("  ⚠ per-wave n is 868–2,680, so a probit cell has SE ≈ 0.037 against a signal range ≈ 1.2.")
print("     At that precision the test rejects rank 1 for a second dimension far too small to")
print("     matter. **`#898`'s error class, pre-empted: the test says IS THERE, the sentence wants")
print("     DOES IT MATTER — so the size is reported whatever the test says.**")
zfull = stats.norm.ppf(np.clip(P0.to_numpy(float), 1e-4, 1 - 1e-4))
rec = stats.norm.cdf(rank1_fit(zfull - base) + base)
err = 100 * np.abs(rec - P0.to_numpy(float))
print(f"  rebuilding all four series from PC1 alone, the largest error is "
      f"**{err.max():.2f} percentage points** (mean {err.mean():.2f})")
worst = np.unravel_index(err.argmax(), err.shape)
print(f"     worst cell: {LABEL[ITEMS[worst[1]]]} in {WAVES[worst[0]]} — "
      f"observed {100*P0.to_numpy()[worst]:.1f}% vs one-tide {100*rec[worst]:.1f}%")
for j, c in enumerate(ITEMS):
    print(f"     {LABEL[c]:22s} max error {err[:, j].max():5.2f} pts · mean {err[:, j].mean():4.2f}")

print("\n=== (5) G3/G4 — 4 thresholds × {all 21 waves · 1988–2018 only, the MODE axis} ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for tname, thr in THRESH.items():
    for sname, keep in (("all waves", WAVES), ("pre-2021 only", PRE)):
        Pk, Nk = series(thr, keep)
        sh, _ = pc1_share(Pk)
        nl = null_shares(Pk, Nk, 300)
        p05 = float(np.percentile(nl, 5))
        zk = stats.norm.ppf(np.clip(Pk.to_numpy(float), 1e-4, 1 - 1e-4))
        bk = zk.mean(0, keepdims=True)
        e = 100 * np.abs(stats.norm.cdf(rank1_fit(zk - bk) + bk) - Pk.to_numpy(float))
        rows.append((tname, sname, len(keep), float(sh), p05, float(e.max())))
        print(f"  {tname:18s} {sname:14s} T={len(keep):2d}  PC1 {sh:.5f}  null5th {p05:.5f}  "
              f"⇒ {'BELOW' if sh < p05 else 'inside':6s}  max error {e.max():5.2f} pts")
below = sum(1 for r in rows if r[3] < r[4])
below_pre = sum(1 for r in rows if r[1] == "pre-2021 only" and r[3] < r[4])
print(f"\n  **grid: {below}/{len(rows)} cells reject rank 1** · of the pre-2021 cells, "
      f"{below_pre}/{sum(1 for r in rows if r[1]=='pre-2021 only')} reject")
print(f"  max reconstruction error across the whole grid: {max(r[5] for r in rows):.2f} points")

print("\n=== (6) THE CONDITIONAL KILL ===")
G = Gate("One tide the four acts ride, or four separate histories?")
G.plant_direction_from_sweep("positive: a planted second factor drives PC1's share down",
                             sweep, baseline=OFF, baseline_spread=OSD, half_of=max(OFF - O05, 1e-5))
G.offset_control("PC1 share vs the rank-1 null", SHARE, OFF, OSD,
                 null_kind="RANK-1 BINOMIAL RESAMPLING NULL at the observed per-wave n — the rank-1 "
                           "approximation is fitted, every cell redrawn as a binomial at that "
                           "wave's actual n, and PC1's share recomputed")
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif SHARE >= O05:
    VERDICT, WORLD = "CONFIRMED", "TIDE — one moving factor reproduces all four series"
elif below_pre >= 3:
    VERDICT, WORLD = "OVERTURNED", ("FOUR — a second dimension is required, and it survives dropping "
                                    "every mode-changed wave")
else:
    VERDICT, WORLD = "OVERTURNED", ("MODE — rank 1 is rejected on the full series but survives on "
                                    "1988–2018; the second dimension is the INSTRUMENT")
print(G)
print(f"\n  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print(f"  ⇒ and the SIZE, which is the part a person can use: rebuilding every series from one tide "
      f"is wrong by at most **{max(r[5] for r in rows):.2f} percentage points**.")

art = dict(entry=900, round="E03·A107·R338", verdict=VERDICT, world=WORLD,
           waves=WAVES, n_per_wave={c: [int(NPW[c].min()), int(NPW[c].max())] for c in ITEMS},
           homosex_code5=n5, pc1_share=SHARE, null_median=OFF, null_sd=OSD, null_p05=O05,
           positive_sweep=sweep, positive_sd=sw_sd, positive_ok=bool(PC_OK),
           max_reconstruction_error_pts=float(err.max()),
           per_item_max_error={c: float(err[:, j].max()) for j, c in enumerate(ITEMS)},
           series_first_last={c: [float(P0[c].iloc[0]), float(P0[c].iloc[-1])] for c in ITEMS},
           grid_rows=rows, grid_below=below, grid_below_pre2021=below_pre,
           prior_art="the DIVERGENCE of these series is textbook (Twenge, Sherman & Wells 2015); "
                     "a one-factor model with different loadings PREDICTS different rates, so "
                     "divergence is not evidence against it — the RANK test is the new part",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "one_tide_or_four.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'one_tide_or_four.json'}")
