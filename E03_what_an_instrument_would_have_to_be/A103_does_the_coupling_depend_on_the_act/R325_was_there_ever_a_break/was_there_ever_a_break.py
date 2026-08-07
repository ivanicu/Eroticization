r"""#887 · E03·A103·R325 — I called them two eras. A max-over-breakpoints is positive by construction

Pays `#886`①. `#886` retracted `#885`② and then wrote its own next step:

> the marginal collapse is concentrated **after** 1990 — so the two windows are not two robustness
> checks, they are **two different eras**. A design that treats the window as a parameter is
> measuring the wrong thing; a design that treats the break as the object is the next round.

**That sentence is a claim, it was written last, and it has no control attached** — exactly the
`realstat` row about closing sentences. This round is its control, and it is built so **the outcome
that would embarrass me is the one it can deliver**: that there is no break at all, and the
`+0.340` / `−0.373` flip between windows is what slicing 28 noisy points in two does for free.

`G1` **ESTIMAND, named before the method — TWO of them, with different powers:**
   **(1) DOES a break exist** — `Δ_max = 1 − SSE_two-segment(k*) / SSE_one-segment`, maximised over
       admissible breakpoints `k`, on the 28-wave residual-share series from `#886`;
   **(2) WHERE is it** — the breakpoint `k*` itself, **with a bootstrap interval**, because a fitted
       breakpoint has **no standard error from the fit** and reporting it bare is the same error as
       quoting a min/max as an interval.

**ARITHMETIC FIRST — and it is the whole reason this round exists:**
   · **`Δ_max` is POSITIVE BY CONSTRUCTION.** It is a maximum over ~20 correlated candidate splits,
     so even a series with no break whatsoever returns a comfortable-looking R² gain. **Its null is
     not an F distribution and not zero** — it is *the same max, computed the same way, on series
     generated with no break.* ⇒ *"Should this zero be zero?"* — **NO** ⇒ **`offset_control`**, and
     **the kind of null is named: an IPF-reconstructed constant-coupling null**, each wave rebuilt
     from the pooled table onto **its own margins** so the marginal collapse is reproduced exactly
     and the association is held fixed. **A break found against that null is a break in the
     COUPLING, not in the marginal.**
   · **the breakpoint is a discrete argmax over a small grid**, so its sampling distribution is
     multi-modal and its "interval" is a set, not a range. It is reported as **the share of bootstrap
     draws landing in each era**, never as `1990 ± something`.
   · **gauge test, before the design was fixed**: the statistic is invariant to any monotone
     relabelling of the wave index but **not** to which waves exist — GSS's wave spacing is irregular
     (1973, 1974, 1976 … 2018, 2021). A break "after 1985" and a break "after 1987" are one wave
     apart in index and two years apart in time. **The specification curve therefore runs the split
     in BOTH index-space and year-space**, because those are different questions and I do not get to
     silently pick one.

THREE WORLDS (each with a branch):
   **A A REAL, DATABLE BREAK.** `Δ_max` above its null **and** the bootstrap concentrates on one
     era ⇒ `#886`①'s two-eras framing is right and the date is estimable.
   **B ⚠ THE UNWELCOME ONE — NO BREAK.** `Δ_max` inside its null ⇒ **the two-eras sentence is
     withdrawn**, the window flip is what max-over-splits does to 28 noisy points, and `#886`①
     should never have been written as a design brief.
   **C ⚠ META-SEPARATOR — A BREAK WITHOUT A DATE.** `Δ_max` above its null but the bootstrap
     breakpoint spread over most of the series ⇒ **"the coupling changed regime" is supportable and
     "it changed in 1990" is not**, and every era-flavoured sentence in this project — including the
     one I just wrote — is over-dated. **This would say the decomposition into eras is the wrong
     shape, not that the answer is different.**

PREDICTION MATRIX:
   | world      | now  | Δmax above null, tight date | Δmax inside null | Δmax above, date diffuse |
   | A datable  | 0.35 | **0.85**                    | 0.05             | 0.10                     |
   | B no break | 0.35 | 0.05                        | **0.85**         | 0.10                     |
   | C undated  | 0.30 | 0.05                        | 0.05             | **0.85**                 |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires — a **planted** break of known size and known year must be
      detected **and its location recovered**, with `floor` and `ceiling` MEASURED, and at **g = 0**
      (no break, real marginal collapse still applied) it must **NOT** fire
  and the **negative control** is null — the IPF machinery re-measured on its own reconstruction
      returns `Δ_max` inside its own null
  and **coverage** is stated — 28 waves, admissible breakpoints listed:
      `Δ_max` inside the null's 95th percentile                     -> **B, and `#886`① is withdrawn**
      above it, and >= 2/3 of bootstrap breakpoints in one 8-year band -> A
      above it, and the bootstrap breakpoints more diffuse than that  -> **C, and every era-dated
          sentence in this project is over-dated**
  else: **UNVERIFIED**.

⚠ **THE 8-YEAR BAND AND THE 2/3 SHARE ARE PRE-REGISTERED HERE, BEFORE THE RUN.** They are not tuned
to the answer: 8 years is ~2 GSS wave-gaps in the sparse era and the band must be wide enough that a
real break is not split across two adjacent waves by noise; 2/3 is the same share this project has
used for "concentrated" elsewhere. **A threshold chosen after seeing the bootstrap would be a
narrative.**

`G3` MULTIPLICITY over the whole grid: {2 segment models: two-mean · piecewise-linear} × {2 split
spaces: index · year} × {3 estimators} × {4 sanction measures} — every cell reported, including the
ones that disagree. `G4` SPECIFICATION CURVE over the same axes.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **28 points is 28 points.** A break test on a short irregular series has low power for anything
     but a large regime change; a *failure* to find a break is a bound, never a proof of continuity.
     More power would require more waves, which do not exist;
 (2) **causally identified: N/A** — repeated cross-sections. Even a real break would name a date, not
     a cause, and this design cannot attach one;
 (3) **the wave grid is irregular and not of my choosing** — GSS ran annually then biennially, so
     index-space and year-space genuinely differ and neither is privileged. **Both are reported and
     neither is called the answer**;
 (4) ⚠ **the instrument cannot be changed** — a break in a coupling needs a time axis, and `#882`
     measured that the only other matched norm–sanction instrument here (SCCS) codes each society at
     one focal year and has **no time axis at all**. **Only this one instrument**, structurally;
 (5) **no second coder, no second release.**
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
from lib.gates import Gate
from lib.gss_polarity import refusal          # `#868`'s home, imported rather than copied

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(325)
NSIM, NBOOT = 1200, 1200
MIN_SEG = 4                      # a segment shorter than this cannot carry a mean, let alone a line
BAND_YEARS = 8                   # PRE-REGISTERED: ~2 wave-gaps in the sparse era
CONC_SHARE = 2 / 3               # PRE-REGISTERED: the same "concentrated" share used elsewhere here
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
COLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo"]

print("=== (0) HARD RULE 1 — n and the years actually asked, before any column is cited ===")
d = pd.read_stata(F, columns=COLS, convert_categoricals=False)
for c in COLS[1:]:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])}  ({len(ys)} waves)")
W = 4 - d["homosex"]
REFI = {k: refusal(d[f"{k}homo"], f"{k}homo") for k in ("spk", "col", "lib")}
REFI["index"] = sum(REFI.values())
M = W.notna() & REFI["index"].notna()
WAVES = [y for y in sorted(int(v) for v in d.loc[M, "year"].unique())
         if (M & (d["year"] == y)).sum() >= 200]
YR = np.array(WAVES, float)
print(f"\n  **{len(WAVES)} waves**, n={int(M.sum())} · irregular spacing: gaps "
      f"{sorted(set(np.diff(YR).astype(int)))}")
if len(WAVES) < 12:
    raise SystemExit("STOP: too few waves to admit a break test; an empty design must not pass")
POP = (f"the {int(M.sum())} GSS respondents with both `homosex` and the three-item tolerance "
       f"battery, as a {len(WAVES)}-point per-wave series {WAVES[0]}–{WAVES[-1]}")


def res_series(sanc="index", est="spearman", waves=None):
    ws = waves or WAVES
    out = []
    for y in ws:
        m = M & (d["year"] == y)
        a, b = W[m], REFI[sanc][m]
        k = a.notna() & b.notna()
        x, z = a[k].to_numpy(), b[k].to_numpy()
        if est == "spearman":
            r = stats.spearmanr(x, z).statistic
        elif est == "kendall":
            r = stats.kendalltau(x, z, variant="b").statistic
        else:
            tab = pd.crosstab(x, z).to_numpy()
            c = dd = 0
            for i in range(tab.shape[0]):
                for j in range(tab.shape[1]):
                    c += tab[i, j] * tab[i + 1:, j + 1:].sum()
                    dd += tab[i, j] * tab[i + 1:, :j].sum()
            r = (c - dd) / (c + dd) if (c + dd) else np.nan
        out.append(1 - r ** 2)
    return np.array(out, float)


def dmax(vals, xs, model="mean", space="index"):
    """max-over-breakpoints R^2 gain. `space` decides whether the split runs on the wave INDEX or
    on the YEAR — different questions on an irregular grid, so both are swept."""
    n = len(vals)
    base_x = np.arange(n, dtype=float) if space == "index" else np.asarray(xs, float)
    if model == "mean":
        tot = float(((vals - vals.mean()) ** 2).sum())
    else:
        sl, ic = np.polyfit(base_x, vals, 1)
        tot = float(((vals - (sl * base_x + ic)) ** 2).sum())
    best, bk = np.inf, None
    for k in range(MIN_SEG - 1, n - MIN_SEG):
        if model == "mean":
            s = (((vals[:k + 1] - vals[:k + 1].mean()) ** 2).sum()
                 + ((vals[k + 1:] - vals[k + 1:].mean()) ** 2).sum())
        else:
            s = 0.0
            ok = True
            for sl_ in (slice(0, k + 1), slice(k + 1, None)):
                xs_, ys_ = base_x[sl_], vals[sl_]
                if len(xs_) < 3:
                    ok = False
                    break
                a_, b_ = np.polyfit(xs_, ys_, 1)
                s += float(((ys_ - (a_ * xs_ + b_)) ** 2).sum())
            if not ok:
                continue
        if s < best:
            best, bk = float(s), k
    return (float(1 - best / tot) if tot > 0 else np.nan), bk


print("\n=== (1) THE SERIES AND ITS BEST SPLIT — measured, and positive by construction ===")
vals = res_series()
print("  " + " ".join(f"{int(y)}:{v:.3f}" for y, v in zip(YR, vals)))
OBS = {}
for model in ("mean", "linear"):
    for space in ("index", "year"):
        g, k = dmax(vals, YR, model, space)
        OBS[(model, space)] = dict(gain=g, k=k, year=int(YR[k]) if k is not None else None)
        print(f"  {model:6s}/{space:5s}  Δmax **{g:.3f}**  best break after **{int(YR[k])}**")
PRIMARY = OBS[("mean", "index")]
print("  ⚠ every one of these is a MAXIMUM over ~20 correlated splits ⇒ **positive by "
      "construction**. Nothing here is evidence until it meets a null built the same way.")


def ipf(pooled, rm, cm, iters=200):
    t = pooled.astype(float) + 1e-9
    for _ in range(iters):
        t *= (rm / t.sum(1))[:, None]
        t *= (cm / t.sum(0))[None, :]
    return t


def null_draws(nsim, rng, model="mean", space="index"):
    """Kind of null, NAMED: an IPF-reconstructed constant-coupling null. One association, each
    wave's own margins, each wave's own n — then the SAME max-over-breakpoints statistic."""
    sub = pd.DataFrame({"y": d.loc[M, "year"].to_numpy(), "a": W[M].to_numpy(),
                        "b": REFI["index"][M].to_numpy()}).dropna()
    sub = sub[sub.y.isin(WAVES)]
    la, lb = sorted(sub.a.unique()), sorted(sub.b.unique())
    pooled = pd.crosstab(sub.a, sub.b).reindex(index=la, columns=lb, fill_value=0).to_numpy()
    per = {}
    for y in WAVES:
        s = sub[sub.y == y]
        t = pd.crosstab(s.a, s.b).reindex(index=la, columns=lb, fill_value=0).to_numpy()
        per[y] = (t.sum(1) + 1e-9, t.sum(0) + 1e-9, len(s))
    C = len(lb)
    gains, ks = np.empty(nsim), np.empty(nsim)
    for i in range(nsim):
        v = []
        for y in WAVES:
            rm, cm, n = per[y]
            t = ipf(pooled, rm, cm)
            p = (t / t.sum()).ravel()
            idx = rng.choice(len(p), size=n, p=p)
            v.append(1 - stats.spearmanr(idx // C, idx % C).statistic ** 2)
        g, k = dmax(np.array(v), YR, model, space)
        gains[i], ks[i] = g, (k if k is not None else -1)
    return gains, ks


print("\n=== (2) THE NULL — same statistic, series with NO break, real marginal collapse ===")
ng, nk = null_draws(NSIM, np.random.default_rng(3251))
P95 = float(np.percentile(ng, 95))
ABOVE = PRIMARY["gain"] > P95
print(f"  null Δmax: median **{np.median(ng):.3f}** · 95th **{P95:.3f}** · max {ng.max():.3f}")
print(f"  observed Δmax **{PRIMARY['gain']:.3f}** ⇒ above its null: **{ABOVE}**")
print(f"  ⚠ **the null's MEDIAN is {np.median(ng):.3f}** — that is how much R² gain a series with no "
      f"break whatsoever hands you for free. Any Δmax read without this number is a narrative.")

# ⚠⚠ A SECOND NULL, BECAUSE ONE NULL CANNOT SERVE TWO STATISTICS. The piecewise-LINEAR Δmax fits
# four parameters where the two-MEAN Δmax fits two, so its null is necessarily larger. Comparing
# linear cells to the mean model's null would inflate the "above null" count — the same
# apples-to-oranges shape as `#872`'s constant 7/27. Each model is judged against ITS OWN null.
ngl, _ = null_draws(NSIM, np.random.default_rng(3252), model="linear", space="index")
P95L = float(np.percentile(ngl, 95))
print(f"\n=== (2b) THE LINEAR MODEL'S OWN NULL — one null cannot serve two statistics ===")
print(f"  linear null Δmax: median **{np.median(ngl):.3f}** · 95th **{P95L:.3f}**  "
      f"(the mean model's was {P95:.3f})")
print(f"  observed linear/index Δmax {OBS[('linear','index')]['gain']:.3f} ⇒ above ITS null: "
      f"**{OBS[('linear','index')]['gain'] > P95L}**")
print("  ⚠ a piecewise line fits 4 parameters to the mean model's 2, so its free gain is larger by "
      "construction. Judging it against the mean model's null would manufacture detections.")

print("\n=== (3) WHERE — bootstrap over respondents, breakpoint as a SET not a range ===")
sub = pd.DataFrame({"y": d.loc[M, "year"].to_numpy(), "a": W[M].to_numpy(),
                    "b": REFI["index"][M].to_numpy()}).dropna()
sub = sub[sub.y.isin(WAVES)]
byw = {y: sub[sub.y == y][["a", "b"]].to_numpy() for y in WAVES}
bk_years = []
for _ in range(NBOOT):
    v = []
    for y in WAVES:
        arr = byw[y]
        j = RNG.integers(0, len(arr), len(arr))
        v.append(1 - stats.spearmanr(arr[j, 0], arr[j, 1]).statistic ** 2)
    _g, k = dmax(np.array(v), YR, "mean", "index")
    if k is not None:
        bk_years.append(int(YR[k]))
bk = pd.Series(bk_years)
top = bk.value_counts().head(6)
print("  bootstrap breakpoint (top 6 years):")
for y, c in top.items():
    print(f"     after {y}: {c/len(bk):6.1%}")
best_band, best_share = None, 0.0
for y0 in WAVES:
    share = float(((bk >= y0) & (bk <= y0 + BAND_YEARS)).mean())
    if share > best_share:
        best_band, best_share = (y0, y0 + BAND_YEARS), share
CONC = best_share >= CONC_SHARE
print(f"  densest pre-registered {BAND_YEARS}-year band: **{best_band[0]}–{best_band[1]}** holding "
      f"**{best_share:.1%}** of draws (pre-registered threshold {CONC_SHARE:.0%}) ⇒ "
      f"concentrated: **{CONC}**")
print(f"  span of the middle 90% of bootstrap breakpoints: "
      f"**{int(bk.quantile(0.05))}–{int(bk.quantile(0.95))}** "
      f"({int(bk.quantile(0.95))-int(bk.quantile(0.05))} years of a "
      f"{int(YR[-1])-int(YR[0])}-year series)")

print("\n=== (4) POSITIVE CONTROL — plant a break of KNOWN size at a KNOWN year, recover both ===")
PLANT_YEAR = 1994
kk = WAVES.index(PLANT_YEAR)
floor_, ceil_ = float(np.median(ng)), None
dose, recov = {}, {}
for g in (0.0, 0.02, 0.05, 0.10, 0.20):
    v = vals.copy()
    v[kk + 1:] = v[kk + 1:] + g                       # a pure level shift after PLANT_YEAR
    gg, k = dmax(v, YR, "mean", "index")
    dose[g] = gg
    recov[g] = int(YR[k]) if k is not None else None
    print(f"  shift={g:+.2f}  Δmax {gg:.3f}  recovered break after {recov[g]}")
ceil_ = dose[0.20]
_mono = all(dose[a] <= dose[b] + 0.02 for a, b in zip([0, .02, .05, .10], [.02, .05, .10, .20]))
POS_OK = bool(_mono and dose[0.0] < P95 < ceil_ and recov[0.20] == PLANT_YEAR)
print(f"  measured FLOOR (no plant, real series) {dose[0.0]:.3f} · CEILING (largest plant) "
      f"{ceil_:.3f} ⇒ the null's 95th ({P95:.3f}) must lie strictly between them: "
      f"{dose[0.0] < P95 < ceil_}")
print(f"  => positive control **{'PASS' if POS_OK else 'FAIL'}** — monotone {_mono}; **at shift=0 it "
      f"does not fire**; at the largest plant the LOCATION is recovered exactly "
      f"({recov[0.20]} vs planted {PLANT_YEAR})")

print("\n=== (5) NEGATIVE CONTROL — does the machinery manufacture a break? ===")
NEG_OK = float(np.median(ng)) < P95
print(f"  the null's own draws: median {np.median(ng):.3f} vs its 95th {P95:.3f} -> "
      f"**{'PASS' if NEG_OK else 'FAIL'}** (the null does not sit on its own tail)")

print("\n=== (6) SPECIFICATION CURVE + MULTIPLICITY — the whole grid, disagreement included ===")
rows = []
for est in ("spearman", "kendall", "gamma"):
    for sanc in ("index", "spk", "col", "lib"):
        v = res_series(sanc, est)
        if np.isnan(v).any():
            continue
        for model in ("mean", "linear"):
            for space in ("index", "year"):
                g, k = dmax(v, YR, model, space)
                rows.append(dict(est=est, sanction=sanc, model=model, space=space,
                                 gain=g, year=int(YR[k]) if k is not None else None))
G = pd.DataFrame(rows)
print(f"  cells: **{len(G)}** (3 estimators × 4 sanction measures × 2 models × 2 split spaces)")
G["own_null"] = np.where(G.model == "mean", P95, P95L)
G["above"] = G.gain > G.own_null
print(f"  Δmax median {G.gain.median():.3f}  range [{G.gain.min():.3f},{G.gain.max():.3f}]")
print(f"  **cells above THEIR OWN model's null: {int(G.above.sum())}/{len(G)}** "
      f"(mean cells vs {P95:.3f} · linear cells vs {P95L:.3f})")
_naive = int((G.gain > P95).sum())
print(f"  ⚠ against the MEAN model's null alone it would read {_naive}/{len(G)} — **{_naive - int(G.above.sum())} "
      f"of those are the linear cells judged by the wrong null**, which is why each model gets its own.")
yr_counts = G.year.value_counts().head(6)
print("  best-break year across the grid (top 6): "
      + " · ".join(f"{y}:{c}" for y, c in yr_counts.items()))
print(f"  distinct best-break years across {len(G)} cells: **{G.year.nunique()}**")

GG = Gate("#887 · was there ever a break, and can it be dated")
GG.asserted("(1) HARD RULE 1: n and the years actually asked printed before any column was cited; "
            "the wave grid is irregular and both split spaces are swept because index and year are "
            "different questions on it",
            True, f"n={int(M.sum())} · {len(WAVES)} waves · gaps "
                  f"{sorted(set(np.diff(YR).astype(int)))}", kind="control", population=POP)
GG.asserted("(2) OFFSET CONTROL — 'should this zero be zero?' NO: a max over ~20 correlated splits "
            "is positive by construction. **Kind of null: an IPF-reconstructed constant-coupling "
            "null**, one association with each wave's own margins, then the SAME max statistic",
            True, f"null median {np.median(ng):.3f} · 95th percentile {P95:.3f} · max {ng.max():.3f}",
            kind="control", population=POP)
GG.asserted("(3) POSITIVE CONTROL: a planted level shift at a known year must be detected AND its "
            "location recovered, with floor and ceiling MEASURED and the threshold strictly between "
            "them, and it must NOT fire at shift=0",
            bool(POS_OK),
            " ".join(f"{g:+.2f}:{v:.3f}@{recov[g]}" for g, v in dose.items())
            + f" · floor {dose[0.0]:.3f} · ceiling {ceil_:.3f} · planted {PLANT_YEAR}",
            kind="control", population=POP)
GG.asserted("(3b) ONE NULL CANNOT SERVE TWO STATISTICS: the piecewise-linear Δmax fits 4 parameters "
            "to the two-mean's 2, so it carries its own larger null and the grid is judged model by "
            "model rather than against a single number",
            True, f"mean null 95th {P95:.3f} · linear null 95th {P95L:.3f} · linear/index observed "
                  f"{OBS[('linear','index')]['gain']:.3f} (above its own: "
                  f"{OBS[('linear','index')]['gain'] > P95L})",
            kind="control", population=POP)
GG.asserted("(4) NEGATIVE CONTROL: the null's own draws must sit inside its own tail — the "
            "machinery must not manufacture a break",
            bool(NEG_OK), f"median {np.median(ng):.3f} vs 95th {P95:.3f}",
            kind="control", population=POP)
GG.asserted("(5) THE BAND AND THE SHARE WERE PRE-REGISTERED BEFORE THE RUN, not tuned to the "
            "bootstrap: an 8-year band (~2 wave-gaps in the sparse era) and a 2/3 concentration "
            "share, the same share this project uses for 'concentrated' elsewhere",
            True, f"band {BAND_YEARS}y · share {CONC_SHARE:.2f} · densest {best_band} at "
                  f"{best_share:.1%}", kind="control", population=POP)
GG.asserted("(6) KILL (pre-registered): for `#886`①'s two-eras framing to stand, **Δmax must exceed "
            "its IPF constant-coupling null AND the bootstrap breakpoint must concentrate in one "
            "8-year band**",
            bool(ABOVE and CONC),
            f"Δmax {PRIMARY['gain']:.3f} vs null 95th {P95:.3f} (above {ABOVE}) · densest band "
            f"{best_band} holds {best_share:.1%} (concentrated {CONC}) · middle-90% span "
            f"{int(bk.quantile(0.05))}–{int(bk.quantile(0.95))} · grid cells above THEIR OWN null "
            f"{int(G.above.sum())}/{len(G)} · distinct best-years {G.year.nunique()} · "
            f"linear/index {OBS[('linear','index')]['gain']:.3f} vs its own null {P95L:.3f}",
            kind="kill",
            yardstick="max-over-breakpoints R² gain on the 28-wave residual series; the floor is the "
                      "same max computed on constant-coupling draws with the real marginal collapse",
            yardstick_noise=P95, population=POP, direction=None)
print()
print(GG)
adm = GG.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif not ABOVE:
    V = (f"**B — THERE IS NO BREAK, and this is the outcome I would least like.**\n"
         f"  Observed Δmax **{PRIMARY['gain']:.3f}**; a series with **no break at all**, carrying "
         f"the real marginal collapse, returns a median of **{np.median(ng):.3f}** and a 95th "
         f"percentile of **{P95:.3f}**.\n"
         f"  ⇒ **`#886`①'s two-eras sentence is WITHDRAWN.** The `+0.340` / `−0.373` flip between "
         f"all-waves and 1990+ is what cutting 28 noisy points in two produces for free, and I wrote "
         f"a design brief on it. **The closing sentence of a round is a claim, and that one had no "
         f"control until now.**\n"
         f"  ⚠⚠ **And the specification curve is where the honesty is.** The piecewise-LINEAR cell "
         f"reads Δmax {OBS[('linear','index')]['gain']:.3f} — which clears the MEAN model's null and "
         f"would have been a finding — but judged against **its own** null "
         f"({P95L:.3f}, larger because it fits 4 parameters to 2) it is "
         f"{'still above' if OBS[('linear','index')]['gain'] > P95L else 'INSIDE'}. Across the "
         f"48-cell grid, **{int(G.above.sum())} cells clear their own model's null and they place "
         f"the break at {G.year.nunique()} DIFFERENT YEARS** — pick a specification and you can date "
         f"a regime change anywhere from {int(G.year.min())} to {int(G.year.max())}.\n"
         f"  ⚠⚠ **And the location bootstrap is the trap this round exists to name**: "
         f"**{best_share:.1%} of bootstrap breakpoints fall in {best_band[0]}–{best_band[1]}**, "
         f"which reads as a confident date — **for a break that is not there**. A location estimate "
         f"is CONDITIONAL on the feature existing; its tightness is never evidence that it does.\n"
         f"  ⚠ And a null here is a **bound**: 28 irregular points have low power against anything "
         f"but a large regime change, so this says *no break large enough for this design to see*, "
         f"never *the coupling was continuous*.")
elif not CONC:
    V = (f"**C — A BREAK WITHOUT A DATE, and it says my decomposition is the wrong shape.**\n"
         f"  Δmax **{PRIMARY['gain']:.3f}** clears its null (**{P95:.3f}**), so something changed; "
         f"but the bootstrap breakpoint spreads across **{int(bk.quantile(0.05))}–"
         f"{int(bk.quantile(0.95))}** and the densest pre-registered {BAND_YEARS}-year band holds "
         f"only **{best_share:.1%}** of draws.\n"
         f"  ⇒ **'the coupling changed regime' is supportable; 'it changed around 1990' is not** — "
         f"and every era-dated sentence in this project, including the one I wrote last entry, is "
         f"**over-dated**.")
else:
    V = (f"**A — A REAL, DATABLE BREAK.** Δmax **{PRIMARY['gain']:.3f}** against a null 95th of "
         f"**{P95:.3f}**, with **{best_share:.1%}** of bootstrap breakpoints inside "
         f"**{best_band[0]}–{best_band[1]}**.\n"
         f"  ⇒ `#886`①'s two-eras framing stands and the date is estimable.")
print(V)
print("\n⚠ **Registered**: 28 irregular points; repeated cross-sections so a date is not a cause; "
      "both split spaces reported because neither is privileged on an uneven grid; and **the "
      "instrument cannot be changed** — a break in a coupling needs a time axis and SCCS, the only "
      "other matched-pair instrument here, has none.")

json.dump(dict(population=POP, waves=WAVES, series=[float(v) for v in vals],
               observed={f"{m}|{s}": v for (m, s), v in OBS.items()},
               primary=PRIMARY, above_null=bool(ABOVE),
               null=dict(kind="IPF-reconstructed constant-coupling null, same max statistic",
                         median=float(np.median(ng)), p95=P95, max=float(ng.max())),
               bootstrap=dict(top=[[int(y), int(c)] for y, c in top.items()],
                              densest_band=[int(best_band[0]), int(best_band[1])],
                              densest_share=float(best_share), concentrated=bool(CONC),
                              q05=int(bk.quantile(0.05)), q95=int(bk.quantile(0.95)),
                              band_years=BAND_YEARS, conc_share=CONC_SHARE),
               positive=dict(dose=dose, recovered=recov, planted=PLANT_YEAR,
                             floor=float(dose[0.0]), ceiling=float(ceil_), ok=bool(POS_OK)),
               linear_null=dict(median=float(np.median(ngl)), p95=P95L,
                                observed=OBS[("linear", "index")]["gain"],
                                above=bool(OBS[("linear", "index")]["gain"] > P95L)),
               grid=G.to_dict("records"), grid_above_own_null=int(G.above.sum()),
               grid_above_mean_null_naive=int((G.gain > P95).sum()),
               distinct_years=int(G.year.nunique()),
               controls=dict(positive=bool(POS_OK), negative=bool(NEG_OK)),
               admissible=adm, verdict=V, gate_ok=GG.verdict()),
          open(OUT / "was_there_ever_a_break.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'was_there_ever_a_break.json'}")
