r"""#901 · E03·A107·R339 — the second dimension has a magnitude; does it have a SHAPE?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#900` established that no fixed-loading one-factor model reproduces the four GSS
                sexual-norm series (PC1 share 0.90014 against a rank-1 binomial null of 0.98307,
                24.4x its own spread) and that one tide misses by up to 6.23 points. **It did not
                ask WHICH ACTS the second dimension separates**, and `#900`① registered why: a
                loading vector read off a REJECTED model is `#893`③'s trap one level up -- **a
                component of a rejected model is not itself a measured object** until it has a null.
Why Now         It is the only remaining question in this arc with psychology in it. "There is a
                second thing moving" is a fact about a matrix; "it separates sex with a minor from
                everything else" would be a fact about people.
Live Worlds     MINOR    the second dimension separates `teensex` from the other three -- an
                         age/consent axis that moved on its own schedule.
                BETRAYAL it separates `xmarsex` -- breach of a promise to a specific person, not
                         disapproval of an act.
                SPLIT    it separates the two liberalising items (`premarsx`, `homosex`) from the
                         two that barely moved -- i.e. it is a RATE axis, not a content axis.
                ⚠ NOISE  ⚠ THE UNWELCOME ONE -- the sign pattern is not stable under resampling, so
                         `#900`'s second dimension is real in MAGNITUDE and has no readable SHAPE,
                         and nothing may be said about which acts. This is what `#900`① feared.
Discriminating  Bootstrap the whole 21x4 design at the observed per-wave n, extract PC2 each time,
Act             ANCHOR its arbitrary sign against a fixed reference, and measure how often the same
                loading sign pattern returns. A pattern that survives resampling is a measured
                object; one that does not is the trap.
Prediction      MINOR    -> the modal pattern isolates `teensex`, and it is stable
Matrix          BETRAYAL -> the modal pattern isolates `xmarsex`, stable
                SPLIT    -> {premarsx, homosex} vs {teensex, xmarsex}, stable
                NOISE    -> modal frequency near the 1/8 a rank-1 world produces
Confound        ⚠ WRITTEN BEFORE THE RUN. **A PCA component's SIGN IS ARBITRARY** (the constitution's
                own guard: an eigenvector's sign is a coin flip and a label written off an unanchored
                one is exactly backwards half the time). Anchoring is therefore not cosmetic: without
                it the "modal pattern" is a modal pattern of a coin. Anchored to a fixed reference
                item, and the anchoring is reported as part of the statistic.
Controls        positive: plant a KNOWN second factor with a KNOWN sign pattern; the procedure must
                recover it, and at g=0 it must NOT · offset: a rank-1 world's PC2 is pure noise, so
                its modal-pattern frequency is the floor -- MEASURED, not assumed to be 1/8
Stopping Rule   If NOISE, the arc closes with a magnitude and no shape, and the page says so.
                Budget: one round. `#111c` -- A107 carries no UNVERIFIED yet.
Cost            21 waves x 4 items, ~3,000 bootstraps. CPU seconds.
Priority        `#900`② (drifting loadings) is UNIDENTIFIABLE here and is registered as such rather
                than run -- with T x 4 free loadings a drifting-loading model fits anything. This
                question is answerable; that one is not.
Expected        If a stable shape: the era has a named second axis and the project can say what it
Transform       is. If NOISE: `#900` is the end of what this instrument can say about eras.
```

⚠ **PRIOR ART.** The two-factor structure of GSS sexual attitudes at a POINT in time is studied; the
question here is the shape of the residual TIME structure after one moving factor is removed, which
is downstream of `#900`'s rank test and is not a restatement of a cross-sectional factor solution.
`D5` on that boundary — **stated as a boundary I am not certain of**, which is the honest grade.

⚠ **AND THE RIVAL `#900` COULD NOT KILL IS STILL ALIVE AND IS NOT KILLED HERE EITHER.** `#900`②: a
one-factor world whose **loadings themselves drift** reproduces anything, because `T x 4` free
loadings have more parameters than the matrix has cells. **It is UNIDENTIFIABLE on this release, and
"planned" is forbidden, so it is registered as impossible rather than deferred.** Everything below
is therefore **conditional on the two-factor reading**, and says so.

`G1` **ESTIMAND**: **the modal sign pattern of PC2's loadings, and the share of bootstrap draws that
return it.** **Population** GSS respondents on the 21 waves where all four items were asked,
1988–2024. **Instrument** GSS `gss7224_r3a` — ⚠ mode changed in 2021 (web/mail push; 2022/2024
mixed-mode), a specification axis. **Baseline** a rank-1 world's own modal frequency. **Regime**
per-wave n 868–2,680 per item.

⚠ **"SHOULD THIS ZERO BE ZERO?" — NO.** Under a rank-1 world PC2 is pure sampling noise, and its
modal sign-pattern frequency is **not** 1/8: the four items have different n and different
variances, so some patterns are commoner than others for free. ⇒ **`offset_control`**, **kind of
null named: a RANK-1 BINOMIAL RESAMPLING NULL at the observed per-wave n, with the modal frequency
MEASURED from it rather than assumed to be 1/8.**

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (a planted pattern is recovered, and at g=0 the modal frequency sits on
                           the rank-1 floor):
       modal frequency > the rank-1 null's 95th percentile AND the same pattern is modal in a
         majority of the specification grid            -> that world, named by the pattern
       otherwise                                        -> NOISE, and the arc closes with a
                                                          magnitude and no shape
else:
       UNVERIFIED
```
`G3`/`G4`: {4 condemnation thresholds} × {all 21 waves · 1988–2018 only}. Whole grid published.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **`#900`②'s drifting-loading rival is UNIDENTIFIABLE here** — more parameters than cells — so
   every shape below is conditional on the two-factor reading, never a refutation of that rival;
② **a sign pattern is a direction, not a magnitude** — it says which acts sit on opposite sides,
   never how far;
③ **21 rows** — the MDE of this stability test against a weak-but-real shape is not computed
   (`#898`'s debt again, and named rather than marked planned);
④ **mode is confounded with period**; **cohort and period are not separated**;
⑤ **cross-instrument N/A — `no second instrument`, `only this one instrument`** (`#897`, `#891`);
⑥ no second coder, no second release, no test–retest.
"""
import json
import pathlib
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
from lib.gates import Gate

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(339)
NBOOT = 3000
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
SHORT = {"premarsx": "pre", "teensex": "teen", "xmarsex": "extra", "homosex": "same"}
MODE_CHANGED = [2021, 2022, 2024]
ANCHOR = "teensex"          # the sign reference, FIXED before the run — see the confound note
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R338 = (ROOT / "E03_what_an_instrument_would_have_to_be/A42_one_tide_or_four_histories/"
        "R338_does_one_moving_factor_reproduce_all_four/results/one_tide_or_four.json")

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
print(f"  common waves {len(WAVES)} · `#900` read from its artifact (never retyped, `#840`'s RULE, "
      f"whose own scope was the `homosex` item alone so only the practice transfers): "
      f"PC1 share {PRIOR['pc1_share']:.5f}, rank-1 null median {PRIOR['null_median']:.5f}")
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


def pc2_pattern(P):
    """Anchored sign pattern of PC2's loadings.

    ⚠ A PCA component's sign is ARBITRARY — the constitution's own guard records a label written off
    an unanchored eigenvector being exactly backwards. Anchoring to a FIXED reference item makes the
    pattern a statement about which acts sit on OPPOSITE sides, which is sign-invariant."""
    z = stats.norm.ppf(np.clip(P.to_numpy(float), 1e-4, 1 - 1e-4))
    z = z - z.mean(0, keepdims=True)
    if not np.isfinite(z).all():
        return None
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    v = vt[1]
    j = ITEMS.index(ANCHOR)
    if v[j] < 0:
        v = -v
    return tuple(int(np.sign(x)) for x in v)


def pretty(pat):
    if pat is None:
        return "n/a"
    pos = [SHORT[c] for c, s in zip(ITEMS, pat) if s > 0]
    neg = [SHORT[c] for c, s in zip(ITEMS, pat) if s < 0]
    return "{" + "+".join(pos) + "} vs {" + "+".join(neg) + "}"


def boot_patterns(P, N, nsim, mat=None):
    """Binomial resampling at the observed per-wave n. `mat` overrides the generating probabilities
    (used for the rank-1 null and for the plant)."""
    nn = N.to_numpy(float)
    src = P.to_numpy(float) if mat is None else mat
    c = Counter()
    for _ in range(nsim):
        k = RNG.binomial(nn.astype(int), np.clip(src, 1e-4, 1 - 1e-4))
        pat = pc2_pattern(pd.DataFrame(k / nn, index=P.index, columns=P.columns))
        if pat is not None:
            c[pat] += 1
    tot = sum(c.values())
    return c, tot


def rank1_probs(P):
    z = stats.norm.ppf(np.clip(P.to_numpy(float), 1e-4, 1 - 1e-4))
    b = z.mean(0, keepdims=True)
    u, s, vt = np.linalg.svd(z - b, full_matrices=False)
    return stats.norm.cdf(u[:, :1] * s[0] @ vt[:1, :] + b), (u, s, vt, b)


print("\n=== (1) THE OBSERVED PC2 PATTERN, and how often resampling returns it ===")
P0, N0 = series(2, WAVES)
OBS = pc2_pattern(P0)
print(f"  observed PC2 sign pattern (anchored on `{ANCHOR}`): **{pretty(OBS)}**  {OBS}")
cnt, tot = boot_patterns(P0, N0, NBOOT)
modal, modal_n = cnt.most_common(1)[0]
MODAL_F = modal_n / tot
print(f"  bootstrap over {tot} draws: modal pattern **{pretty(modal)}** at **{100*MODAL_F:.1f}%**")
for pat, k in cnt.most_common(4):
    print(f"     {pretty(pat):34s} {100*k/tot:5.1f}%")

print("\n=== (2) OFFSET CONTROL — a RANK-1 world's PC2 is pure noise; its modal frequency MEASURED ===")
r1p, _ = rank1_probs(P0)
ncnt, ntot = boot_patterns(P0, N0, NBOOT, mat=r1p)
nmodal, nmodal_n = ncnt.most_common(1)[0]
NULL_F = nmodal_n / ntot
null_fs = np.array(sorted((v / ntot for v in ncnt.values()), reverse=True))
print(f"  rank-1 null: modal pattern {pretty(nmodal)} at **{100*NULL_F:.1f}%** "
      f"({len(ncnt)} distinct patterns seen)")
print(f"  ⚠ **the floor is NOT 1/8 = 12.5%** — the four items differ in n and variance, so some "
      f"patterns are commoner for free. Measured floor: **{100*NULL_F:.1f}%**")
boot_null = []
for _ in range(200):
    c2, t2 = boot_patterns(P0, N0, 60, mat=r1p)
    boot_null.append(c2.most_common(1)[0][1] / t2)
N95 = float(np.percentile(boot_null, 95))
print(f"  rank-1 modal-frequency 95th percentile over 200 replications of 60 draws: **{100*N95:.1f}%**")

print("\n=== (3) POSITIVE CONTROL — plant a KNOWN pattern; at g=0 it must sit on the rank-1 floor ===")
PLANT = np.array([1.0, -1.0, 1.0, -1.0])       # {pre, extra} vs {teen, same} — chosen before the run
z0 = stats.norm.ppf(np.clip(P0.to_numpy(float), 1e-4, 1 - 1e-4))
b0 = z0.mean(0, keepdims=True)
u, s, vt, _ = rank1_probs(P0)[1]
w2 = RNG.standard_normal(len(WAVES)); w2 -= w2.mean(); w2 /= w2.std()
scale = float(np.abs((u[:, :1] * s[0] @ vt[:1, :])).std())
# ⚠ v2. v1 tracked the frequency of WHATEVER PATTERN WAS MODAL, and at low g the modal pattern is
#   still the NOISE pattern while at high g it is the planted one — so the series compared two
#   different objects across the crossover and read as non-monotone (26.5% → 25.8% → 52.0% → …).
#   `realstat` §4, "the control fails for its own reasons": its two sides were not the same object.
#   **v2 tracks the frequency of the PLANTED pattern specifically**, whose floor is that pattern's
#   own rank-1 rate — not the modal rate, which belongs to a different pattern.
want = tuple(int(np.sign(x)) for x in (PLANT if PLANT[ITEMS.index(ANCHOR)] > 0 else -PLANT))
FLOOR_WANT = ncnt.get(want, 0) / ntot
fw = []
for _ in range(200):
    c2, t2 = boot_patterns(P0, N0, 60, mat=r1p)
    fw.append(c2.get(want, 0) / t2)
FW95 = float(np.percentile(fw, 95))
print(f"  the PLANTED pattern {pretty(want)} occurs in a rank-1 world at {100*FLOOR_WANT:.1f}% "
      f"(95th over 200 replications of 60 draws: {100*FW95:.1f}%) — that is ITS floor, not the modal one")
sweep, sw_pat = [], []
for g in (0.0, 0.05, 0.10, 0.20, 0.35, 0.50):
    zz = stats.norm.ppf(np.clip(r1p, 1e-4, 1 - 1e-4)) + g * scale * np.outer(w2, PLANT)
    c3, t3 = boot_patterns(P0, N0, 400, mat=stats.norm.cdf(zz))
    m3, _ = c3.most_common(1)[0]
    sweep.append((g, c3.get(want, 0) / t3))
    sw_pat.append(m3)
    print(f"  g={g:<5.2f} planted-pattern rate {100*sweep[-1][1]:5.1f}%   (modal was {pretty(m3)})")
PC_OK = (sweep[0][1] <= FW95) and (sweep[-1][1] > FW95) and (sw_pat[-1] == want) and \
        all(sweep[i][1] <= sweep[i + 1][1] + 1e-9 for i in range(len(sweep) - 1))
print(f"  planted pattern {pretty(want)} modal at g=0.5: {sw_pat[-1] == want} · "
      f"g=0 sits on ITS OWN floor ({100*sweep[0][1]:.1f}% vs 95th {100*FW95:.1f}%) ⇒ fires: {PC_OK}")

print("\n=== (4) G3/G4 — 4 thresholds × {all 21 waves · 1988–2018 only} ===")
PRE = [y for y in WAVES if y not in MODE_CHANGED]
rows = []
for tname, thr in THRESH.items():
    for sname, keep in (("all waves", WAVES), ("pre-2021 only", PRE)):
        Pk, Nk = series(thr, keep)
        ck, tk = boot_patterns(Pk, Nk, 800)
        mk, kk = ck.most_common(1)[0]
        r1k, _ = rank1_probs(Pk)
        nk, ntk = boot_patterns(Pk, Nk, 800, mat=r1k)
        rows.append((tname, sname, pretty(mk), kk / tk, nk.most_common(1)[0][1] / ntk))
        print(f"  {tname:18s} {sname:14s} modal {pretty(mk):34s} {100*kk/tk:5.1f}%  "
              f"rank-1 floor {100*rows[-1][4]:5.1f}%")
pats = Counter(r[2] for r in rows)
TOP, TOPN = pats.most_common(1)[0]
above = sum(1 for r in rows if r[3] > r[4])
print(f"\n  **grid: {TOPN}/{len(rows)} cells share the modal pattern `{TOP}`** · "
      f"{above}/{len(rows)} cells exceed their own rank-1 floor")

print("\n=== (5) THE CONDITIONAL KILL ===")
G = Gate("Does `#900`'s second dimension have a readable shape, or only a magnitude?")
G.plant_direction_from_sweep("positive: the PLANTED pattern's own recovery rate rises with g",
                             sweep, baseline=FLOOR_WANT, baseline_spread=max(FW95 - FLOOR_WANT, 0.02),
                             half_of=0.15)
G.offset_control("modal frequency vs the rank-1 floor", MODAL_F, NULL_F, max(N95 - NULL_F, 1e-3),
                 null_kind="RANK-1 BINOMIAL RESAMPLING NULL at the observed per-wave n, with the "
                           "modal frequency MEASURED from it rather than assumed to be 1/8")
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif MODAL_F > N95 and TOPN > len(rows) / 2:
    VERDICT, WORLD = "CONFIRMED", f"the second dimension has a SHAPE: {pretty(OBS)} (modal in {TOPN}/{len(rows)} cells)"
else:
    VERDICT, WORLD = "OVERTURNED", ("NOISE — the sign pattern does not survive resampling above its "
                                    "own rank-1 floor; `#900`'s second dimension is real in "
                                    "MAGNITUDE and has no readable SHAPE")
print(G)
print(f"\n  observed modal {100*MODAL_F:.1f}% vs rank-1 floor {100*NULL_F:.1f}% "
      f"(95th {100*N95:.1f}%) · grid agreement {TOPN}/{len(rows)}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ CONDITIONAL ON THE TWO-FACTOR READING. `#900`②'s drifting-loading rival has more")
print("     parameters than this matrix has cells and is UNIDENTIFIABLE here — registered as")
print("     impossible, not deferred, because marking an impossible criterion 'planned' is forbidden.")

art = dict(entry=901, round="E03·A107·R339", verdict=VERDICT, world=WORLD, waves=WAVES,
           anchor=ANCHOR, observed_pattern=list(OBS) if OBS else None, observed_pretty=pretty(OBS),
           modal_pattern=pretty(modal), modal_freq=MODAL_F,
           null_modal=pretty(nmodal), null_freq=NULL_F, null_p95=N95,
           top_patterns={pretty(p): k / tot for p, k in cnt.most_common(5)},
           positive_sweep=sweep, positive_patterns=[pretty(p) for p in sw_pat], positive_ok=bool(PC_OK),
           planted_pattern_floor=FLOOR_WANT, planted_pattern_floor_p95=FW95,
           planted=pretty(want), grid_rows=rows, grid_top=TOP, grid_top_n=TOPN,
           grid_above_floor=above,
           unidentifiable="#900②'s drifting-loading rival has more parameters than cells and cannot "
                          "be excluded on this release; every shape here is conditional on the "
                          "two-factor reading",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "shape_of_the_second_dimension.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'shape_of_the_second_dimension.json'}")
