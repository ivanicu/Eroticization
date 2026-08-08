r"""#905 · E03·A109·R343 — is "not one tide" a fact about SEX, or about attitude series in general?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#900`/`#902` established LINK-FREE that four GSS sexual-norm series are not one
                moving factor (min pairwise |Spearman| 0.4143 vs a comonotone null of 0.7883, 6.1x
                its spread; five links, none saves one factor; 8/8 grid cells). **Every word of
                that was written as if it were a fact about SEXUAL morality. It has never been
                asked whether ANY four GSS attitude series would fail the same test.** If they all
                do, the statistic is detecting a property of survey series and the psychological
                reading collapses -- while every number stays exactly as measured.
Why Now         `#904`① closed A108 and forbade a third statistic on "which act is special", so the
                direction must change. This changes the OBJECT rather than the statistic: same test,
                same instrument, same waves, DIFFERENT CONTENT.
Live Worlds     SPECIFIC   sexual norms depart unusually; most non-sexual batteries are comonotone.
                RANK-ONLY  everything departs and sex departs most -- the claim shrinks to a ranking
                           inside a generic phenomenon.
                ⚠ GENERIC  ⚠ THE UNWELCOME ONE -- sex sits in the middle of the placebo distribution
                           and most batteries depart: `#900`/`#902` are about GSS attitude series,
                           not about sex, and the psychology comes out of them entirely.
                WEAK       sex departs LESS than typical -- `#902` is a mild instance of a generic
                           pattern, which would be the most humiliating outcome available.
Discriminating  Run the IDENTICAL link-free comonotonicity test on many 4-item batteries drawn from
Act             three non-sexual GSS domains (confidence in institutions, spending priorities,
                Stouffer tolerance toward non-sexual targets), restricted to the SAME 21 waves, and
                place the sexual battery in that distribution.
Prediction      SPECIFIC  -> sex above the 90th percentile of placebos AND placebo median z < 2
Matrix          RANK-ONLY -> sex above the 90th AND placebo median z >= 2
                GENERIC   -> sex inside the middle 80% AND placebo median z >= 2
                WEAK      -> sex below the 10th percentile
Confound        ⚠ WRITTEN BEFORE THE RUN: the placebo items are 3-point or binary while the sexual
                items are 4-point, and coarser items have noisier wave ORDER, which INFLATES their
                deficits and would make sex look unusually comonotone for a reason that is not
                about sex. Each battery is judged against ITS OWN comonotone null at ITS OWN n and
                marginals, which prices exactly that -- and the z is reported per battery so the
                inflation is visible rather than assumed away.
Controls        positive: plant a departure into ONE item of a placebo battery; that battery's z
                must rise inside the MEASURED monotone region and must NOT fire at g=0 · negative:
                each battery's own comonotone binomial resampling null, where z has expectation 0
Stopping Rule   One round. If GENERIC or WEAK, `#900`/`#902` keep their numbers and lose their
                subject, and the page says so in those words.
Cost            ~40 batteries x 300 resamples x 6 Spearmans on 21 points. CPU minutes.
Priority        It is the only attack left that can take the epoch's surviving result away, and it
                attacks the SUBJECT rather than the statistic -- which nothing has done yet.
Expected        If GENERIC: E03 ends with a methods result and no psychology, and that is the
Transform       honest end.
```

⚠⚠ **`#901`①'s REMEDY, FOURTH USE.** The outcome space is `(sex percentile: high / middle / low) ×
(placebo median departs: yes / no)` — **six cells, all six assigned before the run**: `high×no`→
SPECIFIC · `high×yes`→RANK-ONLY · `middle×yes`→GENERIC · `low×yes`→WEAK · `middle×no` and `low×no`
→ **contradiction with `#902`'s own measurement, which is `UNVERIFIED` and not a choice between
them.** **No cell is unlisted.** The remedy has changed the outcome on three of its three uses.

`G1` **ESTIMAND**: for each 4-item battery, **`z = (comonotone-null median of min pairwise
|Spearman| − observed min) / null sd`**, and then **the sexual battery's percentile within the
distribution of placebo `z`s**. **Population** GSS respondents on the **21 waves where all four
sexual items were asked, 1988–2024**; every placebo battery is restricted to those same waves so `T`
and the wave spacing are held fixed. **Instrument** GSS `gss7224_r3a` — ⚠ **one instrument, and here
that is the DESIGN and not a limitation**: holding the questionnaire, the waves and the sample
fixed while varying only the CONTENT is what makes the comparison mean anything. Mode changed
mid-series, a specification axis. **Baseline** each battery's own comonotone world. **Regime**
`T = 21`, per-wave n in the high hundreds to low thousands per item.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES.** Under comonotonicity `z` has expectation exactly 0 for
every battery, sexual or placebo. ⇒ **`negative_control`**, **kind of null named: a COMONOTONE
BINOMIAL RESAMPLING NULL, PER BATTERY** — each item's own observed values reassigned to waves in one
common order (link-free by construction, marginals preserved), binomial noise at that wave's actual
n, statistic recomputed for the same battery.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a departure planted into one item of a placebo battery raises THAT
                           battery's z inside the MEASURED monotone region, and does NOT fire at g=0):
       evaluate the six-cell table above
else:
       UNVERIFIED
```
`G3`/`G4`: {2 item reductions: lowest code · at-or-below the median code} × {all 21 waves ·
1988–2018 only} × {within-domain · cross-domain batteries}. Every cell published.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **a placebo domain is a CHOICE** — three were used and they do not exhaust GSS; a domain I did not
   pick could sit anywhere, and the percentile is with respect to the batteries actually drawn;
② **the placebo items are 3-point or binary and the sexual items 4-point** — priced by each
   battery's own null, **not** eliminated;
③ **`#902`①'s common-link assumption stands**, and `#900`②'s drifting-loading rival is still
   unidentifiable;
④ **T = 21**, MDE not computed (`#898`'s standing debt);
⑤ **mode is confounded with period**; **cohort is not separated from period**;
⑥ **cross-INSTRUMENT N/A — `no second instrument`, `only this one instrument`** (`#897`: SCCS's
   matching design resolves no effect at all; `#891`: only GSS ships question text). ⚠ **This round
   is cross-CONTENT, which is a different axis and is not a substitute for it**;
⑦ ⚠ **`[unchallenged]`** — `door ③`; `#899`'s pre-registration table is what a real adversary should
   be scored against;
⑧ no second coder, no second release, no test–retest.
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
RNG = np.random.default_rng(343)
NNULL, NBAT = 300, 40
SEX = ["premarsx", "teensex", "xmarsex", "homosex"]
POOL = {
    "confidence": ["confinan", "conbus", "conclerg", "coneduc", "confed", "conlabor", "conpress",
                   "conmedic", "contv", "conjudge", "consci", "conlegis", "conarmy"],
    "spending": ["natspac", "natenvir", "natheal", "natcity", "natcrime", "natdrug", "nateduc",
                 "natrace", "natarms", "nataid", "natfare"],
    "tolerance": ["spkath", "spkrac", "spkcom", "colath", "colrac", "colcom", "libath", "librac",
                  "libcom"],
}
MODE_CHANGED = [2021, 2022, 2024]
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R340 = (ROOT / "E03_what_an_instrument_would_have_to_be/A42_one_tide_or_four_histories/"
        "R340_was_the_second_dimension_just_the_wrong_link/results/link_free_attack.json")

ALL = SEX + [c for v in POOL.values() for c in v]
print("=== (0) HARD RULE 1 — n, the years actually asked, and the VALUE SET, for EVERY column ===")
d = pd.read_stata(F, columns=["year"] + ALL, convert_categoricals=False)
CODES = {}
for c in ALL:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    CODES[c] = [int(v) for v in sorted(s[c].unique())]
    tag = "SEX" if c in SEX else next(k for k, v in POOL.items() if c in v)
    print(f"  {c:9s} [{tag:10s}] n={len(s):6d}  {int(ys[0])}–{int(ys[-1])} ({len(ys):2d} waves)  "
          f"codes={CODES[c]}")
for c in SEX:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))          # `homosex`'s 5th code, per `#900`
WAVES = [int(y) for y in sorted(set.intersection(
    *[set(d[[c, "year"]].dropna().year.unique()) for c in SEX]))]
print(f"\n  the sexual battery's 21 common waves: {WAVES}")
PRIOR = json.loads(R340.read_text())
print(f"  `#902` read from its artifact (`#840`'s RULE; its own scope was the `homosex` item alone, "
      f"so only the practice transfers): min|ρ| {PRIOR['min_abs_spearman']:.4f} vs null 5th "
      f"{PRIOR['null_p05']:.4f}")
sub = d[d.year.isin(WAVES)]
COV = {c: int(sub[[c, "year"]].dropna().year.nunique()) for c in ALL}
USABLE = [c for c in ALL if COV[c] == len(WAVES)]
print(f"  columns covering ALL {len(WAVES)} of those waves: **{len(USABLE)}/{len(ALL)}** — "
      f"dropped: {[c for c in ALL if c not in USABLE]}")
if len(USABLE) < 12:
    raise SystemExit("STOP: too few usable columns; an empty design must never pass")

REDUCE = {"lowest code": "low", "≤ median code": "mid"}


def series(items, keep, red):
    ps, ns = {}, {}
    for c in items:
        s = sub[sub.year.isin(keep)][[c, "year"]].dropna()
        g = s.groupby("year")[c]
        ns[c] = g.size()
        thr = min(CODES[c]) if red == "low" else int(np.median(CODES[c]))
        ps[c] = g.apply(lambda v: float((v <= thr).mean()))
    return pd.DataFrame(ps).loc[keep], pd.DataFrame(ns).loc[keep]


def min_abs_spearman(P):
    A = P.to_numpy(float)
    vals = []
    for i, j in itertools.combinations(range(A.shape[1]), 2):
        r = stats.spearmanr(A[:, i], A[:, j]).statistic
        if not np.isfinite(r):
            return np.nan
        vals.append(abs(float(r)))
    return min(vals)


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


def zdev(items, keep, red, nnull=NNULL, probs=None):
    P, N = series(items, keep, red)
    if P.isna().any().any() or P.std().min() == 0:
        return np.nan, np.nan
    obs = min_abs_spearman(P)
    src = comonotone_probs(P) if probs is None else probs
    nn = N.to_numpy(float)
    vals = []
    for _ in range(nnull):
        k = RNG.binomial(nn.astype(int), np.clip(src, 1e-4, 1 - 1e-4))
        v = min_abs_spearman(pd.DataFrame(k / nn, index=P.index, columns=P.columns))
        if np.isfinite(v):
            vals.append(v)
    a = np.asarray(vals)
    if len(a) < 30 or a.std(ddof=1) == 0:
        return np.nan, obs
    return float((np.median(a) - obs) / a.std(ddof=1)), obs


print("\n=== (1) THE SEXUAL BATTERY, and 40 PLACEBO batteries on the SAME 21 waves ===")
Z_SEX, OBS_SEX = zdev(SEX, WAVES, "low")
print(f"  sexual battery: observed min|ρ| {OBS_SEX:.4f} ⇒ **z = {Z_SEX:.2f}**")
pool_usable = {k: [c for c in v if c in USABLE] for k, v in POOL.items()}
bats = []
for k, v in pool_usable.items():
    for combo in itertools.islice(itertools.combinations(v, 4), 0, 10):
        bats.append((k, list(combo)))
cross = [c for v in pool_usable.values() for c in v]
for _ in range(NBAT - len(bats)):
    bats.append(("cross-domain", list(RNG.choice(cross, 4, replace=False))))
PL = []
for dom, items in bats:
    z, o = zdev(items, WAVES, "low")
    if np.isfinite(z):
        PL.append((dom, items, float(z), float(o)))
PZ = np.array([p[2] for p in PL])
print(f"  placebo batteries computed: **{len(PL)}** · z median {np.median(PZ):.2f} · "
      f"IQR [{np.percentile(PZ,25):.2f}, {np.percentile(PZ,75):.2f}] · "
      f"10th/90th [{np.percentile(PZ,10):.2f}, {np.percentile(PZ,90):.2f}]")
PCT = float((PZ < Z_SEX).mean() * 100)
print(f"  **the sexual battery sits at the {PCT:.1f}th percentile of the placebo distribution**")
print("  the five most-departing placebos:")
for dom, items, z, o in sorted(PL, key=lambda x: -x[2])[:5]:
    print(f"     z={z:6.2f}  min|ρ|={o:.4f}  [{dom}] {'+'.join(items)}")
print("  the five least-departing placebos:")
for dom, items, z, o in sorted(PL, key=lambda x: x[2])[:5]:
    print(f"     z={z:6.2f}  min|ρ|={o:.4f}  [{dom}] {'+'.join(items)}")

print("\n=== (2) POSITIVE CONTROL — plant a departure into ONE item of a placebo battery ===")
# ⚠ v2. v1 passed the planted probabilities as the NULL generator while computing the OBSERVED
#   statistic from the REAL data — so at g=0 it returned the real battery's z (12.44) instead of 0.
#   **The control's two sides were drawn from two different worlds**, which is the same family
#   `#898` named and R335/R336/R339 each hit once. **v2 generates the OBSERVED series from the
#   planted world too**, and fits that world's own comonotone null to it: at g=0 the synthetic
#   world IS comonotone, so z must land on 0 and the control CAN fail.
pdom, pitems = next((d_, i_) for d_, i_ in bats if d_ == "confidence")
P0, N0 = series(pitems, WAVES, "low")
CP0 = comonotone_probs(P0)
z0 = stats.norm.ppf(np.clip(CP0, 1e-4, 1 - 1e-4))
w2 = RNG.standard_normal(len(WAVES)); w2 -= w2.mean(); w2 /= w2.std()
nn0 = N0.to_numpy(float)


def z_from_world(probs, reps=70):   # ⚠ 25 reps left the g=0.10 cell at −0.07 while its
    # neighbours sat at +0.20 and +0.79 — all inside the estimator's own spread. A control
    # that fails on SIMULATION NOISE is fixed by measuring better, not by excusing it
    # (`#898`'s remedy, reused).
    """Draw a synthetic OBSERVED dataset from `probs`, then score it against ITS OWN comonotone
    null — both sides from the same world, which is the whole repair."""
    out = []
    for _ in range(reps):
        k = RNG.binomial(nn0.astype(int), np.clip(probs, 1e-4, 1 - 1e-4))
        Q = pd.DataFrame(k / nn0, index=P0.index, columns=P0.columns)
        if Q.std().min() == 0:
            continue
        obs = min_abs_spearman(Q)
        src = comonotone_probs(Q)
        vals = []
        for _ in range(150):
            kk = RNG.binomial(nn0.astype(int), np.clip(src, 1e-4, 1 - 1e-4))
            v = min_abs_spearman(pd.DataFrame(kk / nn0, index=P0.index, columns=P0.columns))
            if np.isfinite(v):
                vals.append(v)
        aa = np.asarray(vals)
        if len(aa) > 30 and aa.std(ddof=1) > 0 and np.isfinite(obs):
            out.append((np.median(aa) - obs) / aa.std(ddof=1))
    return float(np.median(out)) if out else np.nan


sweep = []
for g in (0.0, 0.10, 0.20, 0.40, 0.60, 0.90):
    z = z0.copy()
    z[:, 0] += g * float(z0.std()) * w2
    sweep.append((g, z_from_world(stats.norm.cdf(z))))
    print(f"  g={g:<5.2f} z(synthetic world with `{pitems[0]}` planted) {sweep[-1][1]:7.2f}")
TURN = int(np.nanargmax([v for _, v in sweep]))
mono = sweep[:TURN + 1]
PC_OK = (abs(mono[0][1]) < 3.0) and (mono[-1][1] > 3.0) and (TURN >= 2) and \
        all(mono[i][1] <= mono[i + 1][1] + 1e-9 for i in range(len(mono) - 1))
print(f"  ⚠ turning point MEASURED at g={sweep[TURN][0]:g} · monotone in [0, {sweep[TURN][0]:g}] · "
      f"g=0 {mono[0][1]:.2f} lands on ZERO (both sides now drawn from the SAME world) ⇒ the control "
      f"CAN fail ⇒ fires: {PC_OK}")

print("\n=== (3) G3/G4 — 2 reductions × {all 21 waves · 1988–2018} ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for rname, red in REDUCE.items():
    for sname, keep in (("all", WAVES), ("pre-2021", PRE)):
        zs, _ = zdev(SEX, keep, red, 200)
        pz = []
        for dom, items in bats[:24]:
            z, _ = zdev(items, keep, red, 120)
            if np.isfinite(z):
                pz.append(z)
        pz = np.asarray(pz)
        pct = float((pz < zs).mean() * 100) if len(pz) else np.nan
        rows.append((rname, sname, float(zs), float(np.median(pz)), pct, len(pz)))
        print(f"  {rname:14s} {sname:8s} z(sex) {zs:6.2f} · placebo median {np.median(pz):6.2f} "
              f"· sex at the {pct:5.1f}th pct  (n_placebo={len(pz)})")
HIGH = sum(1 for r in rows if r[4] >= 90)
GEN = sum(1 for r in rows if r[3] >= 2)
print(f"\n  **grid: sex above the 90th percentile in {HIGH}/{len(rows)} cells · placebo median z ≥ 2 "
      f"in {GEN}/{len(rows)} cells**")

print("\n=== (4) THE CONDITIONAL KILL — six cells, all assigned before the run ===")
G = Gate("Is 'not one tide' about SEX, or about attitude series in general?")
G.plant_direction_from_sweep("positive: a planted departure raises that battery's z (inside the "
                             "MEASURED monotone region)", mono, baseline=0.0, baseline_spread=1.0,
                             half_of=2.0)
G.negative_control("comonotone world, per battery", 0.0, float(Z_SEX), null_spread=1.0,
                   null_kind="COMONOTONE BINOMIAL RESAMPLING NULL, PER BATTERY — each item's own "
                             "observed values reassigned to waves in one common order (link-free, "
                             "marginals preserved), binomial noise at that wave's actual n; z has "
                             "expectation 0 there by construction")
PLACEBO_DEPARTS = bool(np.median(PZ) >= 2)
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif PCT >= 90 and not PLACEBO_DEPARTS:
    VERDICT, WORLD = "CONFIRMED", "SPECIFIC — sexual norms depart unusually and most placebos do not"
elif PCT >= 90 and PLACEBO_DEPARTS:
    VERDICT, WORLD = "CONFIRMED", ("RANK-ONLY — everything departs and sex departs most; the claim "
                                   "shrinks to a ranking inside a generic phenomenon")
elif PCT <= 10:
    VERDICT, WORLD = "OVERTURNED", ("WEAK — sexual norms depart LESS than a typical GSS battery; "
                                    "`#900`/`#902` are a mild instance of a generic pattern")
elif PLACEBO_DEPARTS:
    VERDICT, WORLD = "OVERTURNED", ("GENERIC — sex sits in the middle and most batteries depart; "
                                    "`#900`/`#902` keep their numbers and lose their subject")
else:
    VERDICT, WORLD = "UNVERIFIED", ("contradiction — placebos do not depart yet sex is only middling "
                                    "among them; that is not a choice between the two measurements")
print(G)
print(f"\n  z(sex) {Z_SEX:.2f} · placebo median {np.median(PZ):.2f} · sex at the {PCT:.1f}th "
      f"percentile · placebo median departs (z≥2): {PLACEBO_DEPARTS}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ A placebo domain is a CHOICE: three were used and they do not exhaust GSS, so the")
print("     percentile is with respect to the batteries actually drawn. `[unchallenged]` — `door ③`.")

art = dict(entry=905, round="E03·A109·R343", verdict=VERDICT, world=WORLD, waves=WAVES,
           z_sex=float(Z_SEX), obs_sex=float(OBS_SEX), percentile=PCT,
           placebo=[dict(domain=p[0], items=p[1], z=p[2], min_rho=p[3]) for p in PL],
           placebo_median=float(np.median(PZ)), placebo_iqr=[float(np.percentile(PZ, 25)),
                                                             float(np.percentile(PZ, 75))],
           placebo_departs=PLACEBO_DEPARTS, positive_sweep=sweep, positive_ok=bool(PC_OK),
           turning_point=sweep[TURN][0], grid_rows=rows, grid_high=HIGH, grid_generic=GEN,
           usable_columns=USABLE, unchallenged=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "sex_or_series.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'sex_or_series.json'}")
