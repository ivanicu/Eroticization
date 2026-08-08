r"""#902 · E03·A107·R340 — I attack `#900`/`#901` at the hole I built into them: the LINK

**COGNITIVE UPDATE CARD**
```
Core Gap        `#901`② invited an adversary: *"if an adversary can show a fixed-loading one-factor
                model reproducing these series, `#900` and `#901` both fall."* I cannot dispatch one
                (`#899` records that), so I am the adversary, and I have found the hole I built in.
                **`#900`'s rank test presupposes the PROBIT LINK.** `z_i(t) = b_i − λ_i·F(t)` is
                rank-1 *on the probit scale*. A one-factor world with a LOGIT link, or any other
                monotone link, is NOT rank-1 after a probit transform -- so **"PC1 share 0.90 < its
                null" may be detecting LINK MISSPECIFICATION and not a second dimension at all.**
Why Now         It is the only attack that would take down both of the last two entries at once,
                and I built the vulnerability myself one round ago. `door ③` says self-review is
                void -- so this round does not claim to be the adversary, it claims to be the
                cheapest thing I can do while one is unavailable.
Live Worlds     TWO   there really are two dimensions; the finding survives a LINK-FREE test.
                LINK  ⚠ THE UNWELCOME ONE -- one factor with some monotone link reproduces the four
                      series, and `#900` + `#901` are both artifacts of choosing probit.
                NEITHER ⚠ META -- the link-free statistic has no power here, so the question is
                      unanswerable at 21 waves and both entries must be downgraded to
                      probit-conditional rather than retracted.
Discriminating  A LINK-FREE INVARIANT. Under one factor with ANY common monotone link and fixed
Act             loadings, p_i(t) = g(b_i − λ_i·F(t)) is a monotone function of F(t) for every item
                ⇒ **all four series are comonotone in t, so every pairwise Spearman across waves is
                ±1 up to sampling noise.** No link appears anywhere in that statement.
Prediction      TWO  -> min |pairwise Spearman| BELOW its comonotone null
Matrix          LINK -> min |pairwise Spearman| INSIDE that null, AND some link makes PC1's share
                        inside its own rank-1 null
                NEITHER -> the positive control cannot separate a planted second factor from the
                        null at this T, so the statistic has no power
Confound        ⚠ WRITTEN BEFORE THE RUN: `xmarsex` moved only 7 points in 36 years, so its
                wave-to-wave ORDER is mostly noise and its pairwise Spearman will be far below 1
                **even under a true one-factor world**. ⇒ the null cannot be "1"; it must be the
                comonotone world's own distribution at the observed per-wave n.
Controls        positive: plant a graded second factor; min|ρ| must fall and must NOT fire at g=0 ·
                offset: a COMONOTONE BINOMIAL RESAMPLING NULL -- each item's OWN observed values
                reassigned to waves in one common order, then binomial noise at the observed n
Stopping Rule   If LINK, `#900` and `#901` are retracted the same day they were written. If TWO,
                they stand and are now link-free. Either way A107 closes. Budget: one round.
Cost            21 waves x 4 items, ~4,000 resamples x 2 statistics. CPU seconds.
Priority        Attacking my own newest result beats extending it, and `#901`② named this exact
                attack before I knew where the hole was.
Expected        If LINK: two entries fall and the epoch ends with nothing about eras.
Transform       If TWO: the claim stops depending on a link I chose for convenience.
```

⚠ **`#901`①'s REMEDY IS APPLIED HERE, ITS FIRST USE.** *When the outcome space is a finite set,
enumerate it and assign every member to a world or to "unlisted", before the run.* The outcome space
here is `(min|ρ| inside / below its null) × (any link puts PC1 inside its null: yes / no)` — **four
cells**, and all four are assigned above: `below × no` → TWO · `inside × yes` → LINK ·
`inside × no` → **LINK-partial**, the link family survives on the weaker statistic only, and the
honest verdict is a downgrade rather than a retraction · `below × yes` → **contradiction between two
statistics, which is `UNVERIFIED` and not a choice between them.** **No cell is unlisted.**

`G1` **ESTIMAND**: **the minimum absolute pairwise Spearman correlation among the four series across
waves**, and, secondarily, **PC1's variance share under five links**. **Population** GSS respondents
on the 21 waves where all four items were asked, 1988–2024, per-wave n 868–2,680. **Instrument** GSS
`gss7224_r3a` — ⚠ **one instrument, `no second instrument`, `only this one instrument`** (`#897`:
SCCS's matching design resolves no effect at all; `#891`: only GSS ships question text) — and its
**mode changed mid-series**, a specification axis. **Baseline** the comonotone world at the observed
n. **Regime** T = 21.

⚠ **"SHOULD THIS ZERO BE ZERO?" — NO.** Under a true one-factor world the pairwise Spearmans are
**1**, not 0, and the observed ones fall below 1 purely from sampling — by an amount that depends on
how far each item actually moved. ⇒ **`offset_control`**, **kind of null named: a COMONOTONE
BINOMIAL RESAMPLING NULL — each item's OWN observed values reassigned to waves in one common order
(so the world is perfectly comonotone and link-free by construction, with every item's marginal
preserved), then binomial noise at that wave's actual n.**

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a planted second factor drives min|ρ| DOWN monotonically and does NOT
                           fire at g=0):
       min|ρ| < the comonotone null's 5th percentile AND no link puts PC1 inside its rank-1 null
                                                                       -> TWO
       min|ρ| inside AND some link puts PC1 inside                      -> LINK, retract #900/#901
       min|ρ| inside AND no link puts PC1 inside                        -> LINK-partial, downgrade
       min|ρ| below  AND some link puts PC1 inside                      -> UNVERIFIED (contradiction)
else:
       UNVERIFIED
```
`G3`/`G4`: {5 links} × {4 condemnation thresholds} × {all 21 waves · 1988–2018 only}. Published whole.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **comonotonicity assumes a COMMON link across items** — four items each with its own monotone link
   is a different and weaker family, and **is not tested**, because with item-specific links the
   model has more freedom than the data has structure;
② **`#900`②'s drifting-loading rival remains UNIDENTIFIABLE** — more parameters than cells;
③ **T = 21** — the MDE of the link-free statistic is not computed (`#898`'s standing debt);
④ **mode is confounded with period**; **cohort is not separated from period**;
⑤ **cross-instrument N/A**; ⑥ no second coder, no second release, no test–retest;
⑦ ⚠ **I am not an adversary.** `door ③`: a reviewer sampled from the weights that wrote `#900` can
   only attack what those weights already anticipated. **This round is marked `[unchallenged]`, and
   `#899`'s pre-registration table is the thing a real adversary should be scored against.**
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
RNG = np.random.default_rng(340)
NSIM = 1200
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
SHORT = {"premarsx": "pre", "teensex": "teen", "xmarsex": "extra", "homosex": "same"}
MODE_CHANGED = [2021, 2022, 2024]
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R338 = (ROOT / "E03_what_an_instrument_would_have_to_be/A42_one_tide_or_four_histories/"
        "R338_does_one_moving_factor_reproduce_all_four/results/one_tide_or_four.json")

LINKS = {
    "probit": (lambda p: stats.norm.ppf(np.clip(p, 1e-4, 1 - 1e-4)), stats.norm.cdf),
    "logit": (lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4))),
              lambda z: 1 / (1 + np.exp(-z))),
    "cloglog": (lambda p: np.log(-np.log(1 - np.clip(p, 1e-4, 1 - 1e-4))),
                lambda z: 1 - np.exp(-np.exp(np.clip(z, -30, 30)))),
    "identity": (lambda p: np.asarray(p, float), lambda z: np.clip(z, 1e-4, 1 - 1e-4)),
    "arcsine": (lambda p: np.arcsin(np.sqrt(np.clip(p, 0, 1))),
                lambda z: np.sin(np.clip(z, 0, np.pi / 2)) ** 2),
}

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
PRIOR = json.loads(R338.read_text())
print(f"  common waves {len(WAVES)} · `#900` (read from its artifact, `#840`'s RULE, whose own "
      f"scope was `homosex` alone so only the practice transfers): PC1 {PRIOR['pc1_share']:.5f}")
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


def min_abs_spearman(P):
    """The LINK-FREE statistic. Under one factor with ANY common monotone link and fixed loadings
    every pair is comonotone or anti-comonotone, so every |Spearman| is 1 up to sampling noise."""
    A = P.to_numpy(float)
    vals = []
    for i, j in itertools.combinations(range(A.shape[1]), 2):
        r = stats.spearmanr(A[:, i], A[:, j]).statistic
        if not np.isfinite(r):
            return np.nan, {}
        vals.append((f"{SHORT[ITEMS[i]]}×{SHORT[ITEMS[j]]}", abs(float(r))))
    return min(v for _, v in vals), dict(vals)


def comonotone_probs(P):
    """A perfectly comonotone world with EVERY item's own observed values preserved.

    The common order is the average within-item rank; each item's own sorted values are reassigned
    along it (reversed if that item runs against the common order). No link is used anywhere."""
    A = P.to_numpy(float)
    rk = np.apply_along_axis(stats.rankdata, 0, A)
    common = rk.mean(1)
    order = np.argsort(np.argsort(common))          # rank of each wave in the common order
    out = np.empty_like(A)
    for j in range(A.shape[1]):
        sgn = np.sign(stats.spearmanr(A[:, j], common).statistic) or 1.0
        srt = np.sort(A[:, j]) if sgn > 0 else np.sort(A[:, j])[::-1]
        out[:, j] = srt[order]
    return out


def resample(P, N, probs, nsim, stat):
    nn = N.to_numpy(float)
    out = []
    for _ in range(nsim):
        k = RNG.binomial(nn.astype(int), np.clip(probs, 1e-4, 1 - 1e-4))
        Q = pd.DataFrame(k / nn, index=P.index, columns=P.columns)
        v = stat(Q)
        v = v[0] if isinstance(v, tuple) else v
        if np.isfinite(v):
            out.append(v)
    return np.asarray(out)


print("\n=== (1) THE LINK-FREE STATISTIC — six pairwise Spearmans across 21 waves ===")
P0, N0 = series(2, WAVES)
MINR, PAIRS = min_abs_spearman(P0)
for k, v in sorted(PAIRS.items(), key=lambda x: x[1]):
    print(f"  |ρ| {v:.4f}   {k}")
print(f"\n  **min |pairwise Spearman| = {MINR:.4f}**  — under ANY one-factor model with a common "
      f"monotone link this is 1 up to sampling noise")

print("\n=== (2) OFFSET CONTROL — the COMONOTONE BINOMIAL RESAMPLING NULL (no link anywhere) ===")
CP = comonotone_probs(P0)
NULL = resample(P0, N0, CP, NSIM, min_abs_spearman)
OFF, OSD, O05 = float(np.median(NULL)), float(NULL.std(ddof=1)), float(np.percentile(NULL, 5))
print(f"  comonotone null: median {OFF:.4f} · sd {OSD:.4f} · **5th percentile {O05:.4f}**")
print(f"  ⚠ **the null is NOT 1.0** — `xmarsex` moved only ~7 points in 36 years, so its wave order "
      f"is mostly noise even in a perfectly comonotone world. Written before the run.")
print(f"  observed {MINR:.4f} vs 5th percentile {O05:.4f} ⇒ "
      f"{'BELOW — the one-factor family is refuted LINK-FREE' if MINR < O05 else 'INSIDE — it survives'}")

print("\n=== (3) POSITIVE CONTROL — plant a second factor; at g=0 it must NOT fire ===")
l2 = np.array([1.0, 1.0, 1.0, -1.0])       # `#901`'s own recovered pattern, used as the plant
w2 = RNG.standard_normal(len(WAVES)); w2 -= w2.mean(); w2 /= w2.std()
zc = stats.norm.ppf(np.clip(CP, 1e-4, 1 - 1e-4))
scale = float(zc.std())
sweep, sw_sd = [], []
for g in (0.0, 0.05, 0.10, 0.20, 0.35, 0.50):
    pg = stats.norm.cdf(zc + g * scale * np.outer(w2, l2))
    v = resample(P0, N0, pg, 300, min_abs_spearman)
    sweep.append((g, float(np.median(v))))
    sw_sd.append(float(v.std(ddof=1)))
for (g, v), s in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} min|ρ| {v:.4f} ± {s:.4f}")
# ⚠ `min|ρ|` is BOUNDED BELOW BY 0, so past a turning point a bigger plant makes the pairs
#   comonotone again with respect to the NEW dominant structure and the statistic RISES. That is
#   the `#885`/`#887` family verbatim -- *a control that cannot PASS because the statistic is even
#   in the thing being planted* -- and this project's recorded remedy is to MEASURE THE TURNING
#   POINT and require monotonicity only inside it, rather than to excuse the non-monotonicity.
TURN = int(np.argmin([v for _, v in sweep]))
mono = sweep[:TURN + 1]
SAT = sweep[TURN][0]
PC_OK = (abs(mono[0][1] - OFF) < 3 * OSD) and (mono[-1][1] < O05) and (TURN >= 2) and \
        all(mono[i][1] >= mono[i + 1][1] - 1e-9 for i in range(len(mono) - 1))
print(f"  ⚠ **turning point MEASURED at g={SAT:g}** (min|ρ| bottoms at {sweep[TURN][1]:.4f} and rises "
      f"to {sweep[-1][1]:.4f} at g=1 — the statistic SATURATES at its 0 bound)")
print(f"  monotone region g ∈ [0, {SAT:g}] · floor(g=0) {mono[0][1]:.4f} lands on the null median "
      f"{OFF:.4f} ⇒ the control CAN fail · bottom {mono[-1][1]:.4f} · threshold {O05:.4f} strictly "
      f"between: {PC_OK}")

print("\n=== (4) THE LINK SWEEP — `#901`②'s attack, run directly: does ANY link save one factor? ===")


def pc1_share_link(P, link):
    f, _ = LINKS[link]
    z = f(P.to_numpy(float))
    z = z - z.mean(0, keepdims=True)
    if not np.isfinite(z).all():
        return np.nan
    sv = np.linalg.svd(z, compute_uv=False)
    return float(sv[0] ** 2 / (sv ** 2).sum())


def rank1_probs_link(P, link):
    f, inv = LINKS[link]
    z = f(P.to_numpy(float))
    b = z.mean(0, keepdims=True)
    u, s, vt = np.linalg.svd(z - b, full_matrices=False)
    return np.clip(inv(u[:, :1] * s[0] @ vt[:1, :] + b), 1e-4, 1 - 1e-4)


link_rows = []
for lname in LINKS:
    obs = pc1_share_link(P0, lname)
    nl = resample(P0, N0, rank1_probs_link(P0, lname), 400, lambda Q: pc1_share_link(Q, lname))
    p05 = float(np.percentile(nl, 5))
    link_rows.append((lname, obs, p05, bool(obs >= p05)))
    print(f"  {lname:9s} PC1 {obs:.5f}  its own rank-1 null 5th pct {p05:.5f}  ⇒ "
          f"{'INSIDE — this link SAVES one factor' if obs >= p05 else 'BELOW'}")
ANY_LINK_SAVES = any(r[3] for r in link_rows)
print(f"\n  **any link puts PC1 inside its own null: {ANY_LINK_SAVES}**")

print("\n=== (5) G3/G4 — 4 thresholds × {all 21 waves · 1988–2018} on the LINK-FREE statistic ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for tname, thr in THRESH.items():
    for sname, keep in (("all waves", WAVES), ("pre-2021 only", PRE)):
        Pk, Nk = series(thr, keep)
        mk, _ = min_abs_spearman(Pk)
        nl = resample(Pk, Nk, comonotone_probs(Pk), 400, min_abs_spearman)
        p05 = float(np.percentile(nl, 5))
        rows.append((tname, sname, len(keep), float(mk), p05, bool(mk < p05)))
        print(f"  {tname:18s} {sname:14s} T={len(keep):2d}  min|ρ| {mk:.4f}  null5th {p05:.4f}  "
              f"⇒ {'BELOW' if mk < p05 else 'inside'}")
below = sum(1 for r in rows if r[5])
below_pre = sum(1 for r in rows if r[1] == "pre-2021 only" and r[5])
print(f"\n  **grid: {below}/{len(rows)} cells refute the one-factor family LINK-FREE** · "
      f"pre-2021 cells: {below_pre}/{sum(1 for r in rows if r[1]=='pre-2021 only')}")

print("\n=== (6) THE CONDITIONAL KILL — all four cells of the outcome space were assigned ===")
G = Gate("Was `#900`/`#901`'s second dimension just the wrong link?")
G.plant_direction_from_sweep("positive: a planted second factor drives min|ρ| down (inside the "
                             "MEASURED monotone region)", mono,
                             baseline=OFF, baseline_spread=OSD, half_of=max(OFF - O05, 1e-4))
G.offset_control("min|ρ| vs the comonotone null", MINR, OFF, OSD,
                 null_kind="COMONOTONE BINOMIAL RESAMPLING NULL — each item's OWN observed values "
                           "reassigned to waves in one common order (perfectly comonotone and "
                           "LINK-FREE by construction, every marginal preserved), then binomial "
                           "noise at that wave's actual n")
BELOW = MINR < O05
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif BELOW and not ANY_LINK_SAVES:
    VERDICT, WORLD = "CONFIRMED", ("TWO — the second dimension survives a LINK-FREE test and no link "
                                   "saves one factor; `#900`/`#901` stand and no longer depend on probit")
elif not BELOW and ANY_LINK_SAVES:
    VERDICT, WORLD = "OVERTURNED", "LINK — one factor with some monotone link reproduces the series; `#900`/`#901` RETRACTED"
elif not BELOW and not ANY_LINK_SAVES:
    VERDICT, WORLD = "OVERTURNED", ("LINK-partial — the one-factor family survives the link-free "
                                    "statistic; `#900`/`#901` DOWNGRADED to probit-conditional")
else:
    VERDICT, WORLD = "UNVERIFIED", ("contradiction — the link-free statistic refutes one factor while "
                                    "some link saves it; that is not a choice between them")
print(G)
print(f"\n  min|ρ| {MINR:.4f} vs null 5th {O05:.4f} ⇒ below: {BELOW} · any link saves one factor: "
      f"{ANY_LINK_SAVES} · grid {below}/{len(rows)}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ `[unchallenged]` — I am not an adversary. `door ③`: a reviewer sampled from the weights")
print("     that wrote `#900` can only attack what those weights already anticipated. `#899`'s")
print("     pre-registration table remains the thing a real adversary should be scored against.")

art = dict(entry=902, round="E03·A107·R340", verdict=VERDICT, world=WORLD, waves=WAVES,
           pairwise=PAIRS, min_abs_spearman=MINR,
           null_median=OFF, null_sd=OSD, null_p05=O05, below=bool(BELOW),
           positive_sweep=sweep, positive_sd=sw_sd, positive_ok=bool(PC_OK),
           turning_point_g=SAT, monotone_region=[m[0] for m in mono],
           saturation_note='min|rho| is bounded below by 0, so past the turning point a bigger plant makes the pairs comonotone again and the statistic RISES — the #885/#887 family, remedied by MEASURING the turning point rather than excusing it',
           link_rows=link_rows, any_link_saves=bool(ANY_LINK_SAVES),
           grid_rows=rows, grid_below=below, grid_below_pre2021=below_pre,
           outcome_space_enumerated="#901①'s remedy, first use: 4 cells (below/inside × link saves/"
                                    "not), all four assigned to a world before the run, none unlisted",
           unchallenged=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "link_free_attack.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'link_free_attack.json'}")
