r"""#906 · E03·A109·R344 — the same attack, pointed at the project's OTHER surviving result

**COGNITIVE UPDATE CARD**
```
Core Gap        `#905` took the SUBJECT away from `#900`/`#902`: "not one tide" is a property of GSS
                attitude series, not of sexual morality, and the sexual battery is among the tamest.
                `#905`② then said the obvious thing out loud -- **every earlier "it is about sex"
                claim in this project is suspect by the same argument, and A105's case-indexing
                (`#892`, `D = +0.6688` at 5.6x its one-factor null) has NEVER been run against a
                placebo domain.** That is the project's other surviving result and the only one left
                that says anything about people.
Why Now         If it is generic, E03 ends with NOTHING about sex, and the honest page says so. That
                is worth knowing before another round is spent extending it.
Live Worlds     SPECIFIC   case-indexing is distinctive to the abortion domain; non-sexual crossed
                           designs sit near their one-factor nulls.
                RANK-ONLY  everything shows case-indexing and abortion shows the most.
                ⚠ GENERIC  ⚠ THE UNWELCOME ONE -- abortion sits in the middle: `#892` keeps its
                           number and loses its subject exactly as `#900`/`#902` just did, and E03
                           has no surviving claim about sexual morality at all.
                WEAK       abortion shows LESS case-indexing than a typical GSS crossed design.
Discriminating  GSS's Stouffer battery is a TARGET x ACTION matrix -- for each of {atheist, racist,
Act             communist, militarist} it asks whether that person may SPEAK, TEACH, and keep a BOOK
                in the library. **That is structurally identical to A105's design**: two "norm" items
                indexed by case (`spk` for two targets) against two "sanction" classes indexed by
                the same cases ({col,lib} for each target). Six non-sexual target pairs, n = 34,092
                to 39,991 each.
Prediction      SPECIFIC  -> abortion above the 90th percentile of the six AND placebo median D
Matrix                       inside its own one-factor null
                RANK-ONLY -> abortion above the 90th AND placebo median beyond its null
                GENERIC   -> abortion inside the middle AND placebo median beyond its null
                WEAK      -> abortion below the 10th percentile
Confound        ⚠ TWO, BOTH WRITTEN BEFORE THE RUN.
                (1) **n differs by 20x** (1,960 vs ~37,000). The null's spread scales with 1/sqrt(n),
                    so a z-comparison is rigged. **The headline comparison is therefore at MATCHED
                    n**: every placebo is subsampled to A105's 1,960 and D is compared as a SIZE.
                (2) **coarseness** -- `#905`'s lesson. A105 was 4-point x binary; Stouffer is binary
                    x binary. The log-interaction is one-factor-proof by algebra so coarseness moves
                    the SPREAD, not the centre, and each design is judged against its own null.
                ⚠ (3) AND ONE THAT CUTS AGAINST ME: A105 DELETED its identical-proposition twins to
                    kill the wording world. Stouffer has no identical twin to delete (speak/teach/
                    library are three different actions), so the placebo is already in the
                    distance-matched regime and is if anything FAVOURED to show case-indexing.
Controls        positive: plant case-indexing into a placebo pair, both sides drawn from the SAME
                world (`#905`'s repair), monotone inside the MEASURED turning point, not firing at
                g=0 · negative: the synthetic one-factor world, where D = 0 by `#892`'s derivation
Stopping Rule   One round. If GENERIC or WEAK, E03 closes with no surviving claim about sex and the
                page says that in those words. Budget: one round.
Cost            6 pairs x subsampling x ~300 sims. CPU minutes.
Priority        It is the last attack that can empty the epoch, and `#905`② named it before I could
                see the answer.
Expected        If GENERIC: the whole of E03 becomes a methods result, and the honest deliverable
Transform       is a set of instruments plus a register of what this data cannot say.
```

⚠⚠ **`#901`①'s REMEDY, FIFTH USE.** Outcome space = `(abortion percentile: high / middle / low) ×
(placebo median beyond its own null: yes / no)` — **six cells, all six assigned before the run**,
with `middle×no` and `low×no` going to **`UNVERIFIED` (a contradiction with `#892`'s own
measurement), not to whichever world is nearest.** The remedy has changed the outcome on four of its
four uses.

`G1` **ESTIMAND**: the **case-indexing interaction on `log|ρ|`** — for a pair of cases `(1,2)`,
`D = [m(norm₁, class₁) − m(norm₁, class₂)] − [m(norm₂, class₁) − m(norm₂, class₂)]` — which is
**identically zero under any one-factor model with any loadings** (`#892`'s derivation, and `#893`③'s
correction that the invariance belongs to the INTERACTION and never to a bracket). **Population**
GSS respondents carrying all six items of a target pair. **Instrument** GSS `gss7224_r3a` — ⚠ **one
instrument, and here that is the DESIGN**: holding the questionnaire fixed while varying only the
CONTENT is what makes the comparison mean anything; **cross-INSTRUMENT remains N/A — `no second
instrument`, `only this one instrument`** (`#897`, `#891`). **Baseline** each design's own
one-factor world. **Regime** matched at `n = 1,960`, with full-n reported beside it.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES**, by the derivation above. ⇒ **`negative_control`**, and the
**kind of null is named: a SYNTHETIC ONE-FACTOR WORLD fitted to each design's own loadings and each
item's own marginal, at the matched n** — where `D = 0` for any loadings, so the resampling spread
is the whole of what the null contributes.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (planted case-indexing raises D inside the MEASURED monotone region, both
                           sides from the same world, and it does NOT fire at g=0):
       evaluate the six-cell table above, on MATCHED n
else:
       UNVERIFIED
```
`G3`/`G4`: {6 target pairs} × {matched n · full n} × {3 estimators}. Every cell published.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **six placebo pairs is six** — the Stouffer battery has four non-sexual targets and that is all;
   the percentile is with respect to those six and a domain I cannot build could sit anywhere;
② **the placebo's "norm" is `spk` and its "sanction" is `col`/`lib`** — all three are tolerance
   items, whereas A105's norm was a MORAL judgement and its sanction a LEGAL one. **The placebo is
   structurally matched and semantically not**, and that is a real difference this design cannot
   remove;
③ **coarseness and n are matched or priced, never eliminated**;
④ **causally identified N/A** — repeated cross-sections;
⑤ ⚠ **`[unchallenged]`** — `door ③`; `#899`'s pre-registration table is what a real adversary should
   be scored against;
⑥ no second coder, no second release, no test–retest.
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
from lib.gss_polarity import refusal

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(344)
NSIM, NSUB = 300, 200
TARGETS = ["ath", "rac", "com", "mil"]          # non-sexual only; `homo` is the sexual target
ACTIONS = ["spk", "col", "lib"]
NAME = {"ath": "atheist", "rac": "racist", "com": "communist", "mil": "militarist"}
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R330 = (ROOT / "E03_what_an_instrument_would_have_to_be/A105_is_it_the_act_or_what_the_two_"
        "questions_ask/R330_is_the_distance_a_wording_fact_or_a_moral_one/results/"
        "wording_or_cases.json")

COLS = [f"{a}{t}" for t in TARGETS for a in ACTIONS]
print("=== (0) HARD RULE 1 — n, the years actually asked, and the VALUE SET, for EVERY column ===")
d = pd.read_stata(F, columns=["year"] + COLS, convert_categoricals=False)
for c in COLS:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])} ({len(ys):2d} waves)  "
          f"codes={[int(v) for v in sorted(s[c].unique())]}")
R = {c: refusal(d[c], c) for c in COLS}         # 1 = would NOT allow; `#868`'s polarity table
PRIOR = json.loads(R330.read_text())
D_ABORTION = float(PRIOR["D"])
N_ABORTION = int(PRIOR["n"])
print(f"\n  `#892`'s abortion D, READ from `R330`'s artifact and not retyped (`#840`'s RULE; its own")
print(f"  scope was the `homosex` item alone, so only the practice transfers): **{D_ABORTION:+.4f}** "
      f"at n = {N_ABORTION}")


def _rho(x, y, est):
    if len(x) < 40 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    if est == "spearman":
        return float(stats.spearmanr(x, y).statistic)
    if est == "kendall":
        return float(stats.kendalltau(x, y, variant="b").statistic)
    tab = pd.crosstab(x, y).to_numpy()
    c = dd = 0
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            c += tab[i, j] * tab[i + 1:, j + 1:].sum()
            dd += tab[i, j] * tab[i + 1:, :j].sum()
    return float((c - dd) / (c + dd)) if (c + dd) else np.nan


def interaction(M, t1, t2, est="spearman"):
    """log|ρ| case-indexing interaction. Zero under ANY one-factor model, for ANY loadings."""
    lg = {}
    for nt in (t1, t2):
        for ct in (t1, t2):
            vals = []
            for a in ("col", "lib"):
                r = abs(_rho(M[f"spk{nt}"], M[f"{a}{ct}"], est))
                if not np.isfinite(r) or r <= 1e-6:
                    return np.nan
                vals.append(np.log(r))
            lg[(nt, ct)] = float(np.mean(vals))
    return (lg[(t1, t1)] - lg[(t1, t2)]) - (lg[(t2, t1)] - lg[(t2, t2)])


def frame(t1, t2):
    need = [f"{a}{t}" for t in (t1, t2) for a in ACTIONS]
    M = pd.DataFrame({c: R[c] for c in need}).dropna()
    return M


def one_factor_null(M, t1, t2, n, nsim=NSIM):
    """A synthetic one-factor world at the matched n, fitted to this design's own loadings and
    each item's own marginal. D = 0 there for any loadings, so the spread is the whole null."""
    cols = list(M.columns)
    z = (M.to_numpy(float) - M.to_numpy(float).mean(0)) / np.where(M.to_numpy(float).std(0) > 0,
                                                                   M.to_numpy(float).std(0), 1)
    gl = z.mean(1)
    lam = {c: float(np.corrcoef(z[:, i], gl)[0, 1]) for i, c in enumerate(cols)}
    mg = {c: float(M[c].mean()) for c in cols}
    out = []
    for _ in range(nsim):
        g = RNG.standard_normal(n)
        sim = {}
        for c in cols:
            u = lam[c] * g + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
            sim[c] = (u > np.quantile(u, 1 - min(max(mg[c], 0.02), 0.98))).astype(float)
        v = interaction(pd.DataFrame(sim), t1, t2)
        if np.isfinite(v):
            out.append(v)
    return np.asarray(out)


print("\n=== (1) SIX NON-SEXUAL TARGET PAIRS, matched to A105's n = 1,960 ===")
rows = []
for t1, t2 in itertools.combinations(TARGETS, 2):
    M = frame(t1, t2)
    d_full = interaction(M, t1, t2)
    subs = []
    for _ in range(NSUB):
        idx = RNG.choice(len(M), N_ABORTION, replace=False)
        v = interaction(M.iloc[idx], t1, t2)
        if np.isfinite(v):
            subs.append(v)
    d_match = float(np.median(subs))
    nul = one_factor_null(M, t1, t2, N_ABORTION)
    p95 = float(np.percentile(np.abs(nul), 95))
    rows.append(dict(pair=f"{NAME[t1]}×{NAME[t2]}", t1=t1, t2=t2, n=int(len(M)),
                     D_full=float(d_full), D_matched=d_match, null_p95=p95,
                     null_med=float(np.median(nul)), beyond=bool(abs(d_match) > p95)))
    print(f"  {NAME[t1]:11s} × {NAME[t2]:11s} n={len(M):6d}  D_full {d_full:+.4f}  "
          f"D@1960 {d_match:+.4f}  one-factor null 95th {p95:.4f}  ⇒ "
          f"{'BEYOND' if abs(d_match) > p95 else 'inside'}")
PD = np.array([r["D_matched"] for r in rows])
PCT = float((PD < D_ABORTION).mean() * 100)
DEPARTS = bool(np.median([abs(r["D_matched"]) for r in rows]) >
               np.median([r["null_p95"] for r in rows]))
print(f"\n  placebo D@1960: median {np.median(PD):+.4f} · range [{PD.min():+.4f}, {PD.max():+.4f}]")
print(f"  **abortion D = {D_ABORTION:+.4f} sits at the {PCT:.1f}th percentile of the six**")
print(f"  placebo median |D| beyond the median null 95th: {DEPARTS}")

print("\n=== (2) POSITIVE CONTROL — plant case-indexing; both sides from the SAME world ===")
t1, t2 = "ath", "com"
M0 = frame(t1, t2)
cols0 = list(M0.columns)
zz = (M0.to_numpy(float) - M0.to_numpy(float).mean(0)) / np.where(M0.to_numpy(float).std(0) > 0,
                                                                  M0.to_numpy(float).std(0), 1)
gl0 = zz.mean(1)
lam0 = {c: float(np.corrcoef(zz[:, i], gl0)[0, 1]) for i, c in enumerate(cols0)}
mg0 = {c: float(M0[c].mean()) for c in cols0}


def world(g, n=N_ABORTION):
    """One synthetic dataset: one-factor plus a planted CASE-INDEXING component of size g."""
    lat = RNG.standard_normal(n)
    case = RNG.standard_normal(n)          # the case-specific dimension
    sim = {}
    for c in cols0:
        own = c[3:]
        s = (1.0 if own == t1 else -1.0) * g
        u = lam0[c] * lat + s * case + np.sqrt(max(1 - lam0[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        sim[c] = (u > np.quantile(u, 1 - min(max(mg0[c], 0.02), 0.98))).astype(float)
    return pd.DataFrame(sim)


sweep = []
for g in (0.0, 0.10, 0.20, 0.35, 0.50, 0.80):
    vals = [interaction(world(g), t1, t2) for _ in range(120)]
    vals = [v for v in vals if np.isfinite(v)]
    sweep.append((g, float(np.median(vals))))
    print(f"  g={g:<5.2f} D {sweep[-1][1]:+.4f}   (synthetic world, observed AND null from the same "
          f"generator)")
TURN = int(np.nanargmax([v for _, v in sweep]))
mono = sweep[:TURN + 1]
NUL0 = one_factor_null(M0, t1, t2, N_ABORTION)
P95_0 = float(np.percentile(np.abs(NUL0), 95))
PC_OK = (abs(mono[0][1]) < P95_0) and (mono[-1][1] > P95_0) and (TURN >= 2) and \
        all(mono[i][1] <= mono[i + 1][1] + 1e-9 for i in range(len(mono) - 1))
print(f"  ⚠ turning point MEASURED at g={sweep[TURN][0]:g} · monotone in [0, {sweep[TURN][0]:g}] · "
      f"g=0 {mono[0][1]:+.4f} sits inside the null 95th {P95_0:.4f} ⇒ the control CAN fail ⇒ "
      f"fires: {PC_OK}")

print("\n=== (3) G3/G4 — 6 pairs × {matched n · full n} × 3 estimators ===")
grid = []
for r in rows:
    M = frame(r["t1"], r["t2"])
    for est in ("spearman", "kendall", "gamma"):
        df = interaction(M, r["t1"], r["t2"], est)
        subs = []
        for _ in range(60):
            idx = RNG.choice(len(M), N_ABORTION, replace=False)
            v = interaction(M.iloc[idx], r["t1"], r["t2"], est)
            if np.isfinite(v):
                subs.append(v)
        dm = float(np.median(subs)) if subs else np.nan
        grid.append((r["pair"], est, float(df), dm))
        print(f"  {r['pair']:25s} {est:9s} D_full {df:+.4f}  D@1960 {dm:+.4f}")
gm = [g[3] for g in grid if np.isfinite(g[3])]
above = sum(1 for v in gm if v > D_ABORTION)
print(f"\n  **grid: {above}/{len(gm)} matched-n cells exceed abortion's {D_ABORTION:+.4f}**")

print("\n=== (4) THE CONDITIONAL KILL — six cells, all assigned before the run ===")
G = Gate("Is case-indexing about SEX either, or about GSS crossed designs?")
G.plant_direction_from_sweep("positive: planted case-indexing raises D (inside the MEASURED "
                             "monotone region)", mono, baseline=0.0,
                             baseline_spread=float(NUL0.std(ddof=1)), half_of=max(P95_0, 1e-4))
G.negative_control("synthetic one-factor world at the matched n", float(np.median(NUL0)),
                   D_ABORTION, null_spread=float(NUL0.std(ddof=1)),
                   null_kind="SYNTHETIC ONE-FACTOR WORLD fitted to each design's own loadings and "
                             "each item's own marginal, at the matched n — D = 0 there for ANY "
                             "loadings, so the resampling spread is the whole of the null")
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif PCT >= 90 and not DEPARTS:
    VERDICT, WORLD = "CONFIRMED", "SPECIFIC — case-indexing is distinctive to the abortion domain"
elif PCT >= 90 and DEPARTS:
    VERDICT, WORLD = "CONFIRMED", ("RANK-ONLY — every crossed design shows case-indexing and "
                                   "abortion shows the most")
elif PCT <= 10:
    VERDICT, WORLD = "OVERTURNED", ("WEAK — abortion shows LESS case-indexing than a typical GSS "
                                    "crossed design")
elif DEPARTS:
    VERDICT, WORLD = "OVERTURNED", ("GENERIC — abortion sits in the middle; `#892` keeps its number "
                                    "and loses its subject, and E03 has no surviving claim about sex")
else:
    VERDICT, WORLD = "UNVERIFIED", ("contradiction — the placebos do not depart yet abortion is only "
                                    "middling among them; not a choice between the two measurements")
print(G)
print(f"\n  abortion D {D_ABORTION:+.4f} · placebo median {np.median(PD):+.4f} · abortion at the "
      f"{PCT:.1f}th percentile · placebo median departs: {DEPARTS}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ THE PLACEBO IS STRUCTURALLY MATCHED AND SEMANTICALLY NOT: its 'norm' is a tolerance")
print("     item and its 'sanction' is another tolerance item, while A105's norm was a MORAL")
print("     judgement and its sanction a LEGAL one. That difference is real and is not removable")
print("     here. ⚠ And A105 DELETED its identical twins; Stouffer has none to delete, so the")
print("     placebo is if anything FAVOURED to show case-indexing — which cuts against me, not for.")

art = dict(entry=906, round="E03·A109·R344", verdict=VERDICT, world=WORLD,
           D_abortion_read_from_R330=D_ABORTION, n_abortion=N_ABORTION,
           pairs=rows, placebo_median=float(np.median(PD)),
           placebo_range=[float(PD.min()), float(PD.max())], percentile=PCT, departs=DEPARTS,
           positive_sweep=sweep, positive_ok=bool(PC_OK), turning_point=sweep[TURN][0],
           null_p95_ath_com=P95_0, grid=grid, grid_above=above, grid_n=len(gm),
           unchallenged=True,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "case_indexing_placebo.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'case_indexing_placebo.json'}")
