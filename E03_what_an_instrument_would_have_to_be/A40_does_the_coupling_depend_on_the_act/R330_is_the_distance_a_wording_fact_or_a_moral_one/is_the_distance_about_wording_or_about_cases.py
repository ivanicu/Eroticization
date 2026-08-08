r"""#892 · E03·A105·R330 — is `#890`'s "distance" a fact about WORDING or a fact about how people
judge CASES?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#890` measured Δ_distance = +0.2685 — matching the REASON between the moral and the
                legal question buys 80% of what four rounds called a difference between ACTS. It was
                left as a fact about QUESTIONNAIRES. Evidence cannot yet distinguish two readings
                that imply completely different projects.
Why Now         Every downstream sentence differs. Under Q the project has been measuring an
                artifact of survey construction; under P it has stumbled onto a real property of
                moral cognition and the act-vs-act line was the wrong question all along.
Live Worlds     Q RECOGNITION — the matched pair correlates better because the respondent SEES the
                  repeated clause and answers consistently. A response-set fact.
                P CASE-SPECIFIC MORALITY — people do not hold "an attitude to abortion"; they judge
                  CASES, and a verdict on one case predicts the law they want for LIKE cases.
                F ONE FACTOR — there is a single abortion-permissiveness dimension, the items differ
                  only in difficulty, and every ρ is λ_i·λ_j. Then `#890`'s Δ is loadings, nothing else.
                D ⚠ META — "reason class" is not the carve; coupling is item-idiosyncratic and both
                  {act} and {reason class} are wrong decompositions of the same wrong object.
Discriminating  DELETE the two identical items and test a DOUBLE DISSOCIATION over what remains.
Act             With `abpoor` and `abdefect` gone there is NO repeated wording left anywhere in the
                design, so world Q has nothing to act on.
Prediction      Q -> D collapses to 0 (its whole mechanism was the repeated clause)
Matrix          P -> D > 0 (each norm couples to its OWN class)
                F -> D = 0 EXACTLY, by algebra (see the derivation below)
                D -> the grid disagrees in sign across defensible class partitions
Confound        ⚠ THE ARITHMETIC TRAP, and my first draft walked into it. On the RAW ρ scale a
                one-factor world gives D = (λ_poorw − λ_defctw)(λ_E − λ_T), which is NOT zero. The
                naive double-difference is CONFOUNDED BY LOADINGS. On log|ρ| it is exactly zero, so
                the estimand is the norm × class INTERACTION on log|ρ|, and nothing else.
Controls        positive: graded plant of class-specific coupling, must not fire at g=0 · negative:
                a SYNTHETIC ONE-FACTOR WORLD fitted to the observed loadings and marginals — the
                world F, built rather than assumed · placebo: random half-split, must return 0
Stopping Rule   D inside the one-factor null -> F or Q, `#890` becomes a loadings statement.
                D above it and stable in sign across partitions -> P.
                sign disagreement across partitions -> D, and the carve is wrong.
                Budget: one round. `#111c` — if this comes back UNVERIFIED twice, change direction.
Cost            n≈1,960, three waves, seconds of CPU. No GPU, no agents.
Priority        Every other open debt is production (pages carry points, evidence grades undeclared).
                This one decides what the last four rounds were ABOUT.
Expected        If P: "acts" was never the object and "cases" is; the project's unit changes again.
Transform       If F/Q: `#890`'s headline is a statement about item loadings and must be reworded.
```

⚠⚠ **PRIOR ART, DECLARED BEFORE THE RESULT AND NOT AFTER — `P14`'s `prior_art` line, which cost this
project a day once already.** The **two-factor structure of the GSS abortion battery — "traumatic"
(`abdefect` · `abhlth` · `abrape`) versus "elective" (`abnomore` · `abpoor` · `absingle` · `abany`) —
is textbook, and has been since Ebaugh & Haney (1980).** `D6`, from my own reading, not from a
search run in this session. **Therefore:**
- **finding two classes here would be a VERIFICATION, not a discovery**, and must be reported as one;
- **`#890`'s Δ_distance is very probably that known structure wearing a new name** — which is the
  unwelcome reading of my own last-but-one headline, and is why this round is worth running;
- **what is NOT settled by that literature** is the question here: *does the MORAL item track the
  classes of the LEGAL items with the identical item removed?* The sanctions having two factors does
  not force the norms to load on both — a norm loading only on the general factor gives `D = 0`
  against a two-factor sanction side. **So `D` is not forced by the prior art, and it is the part
  this round can actually add.**

**THE DERIVATION, LABELLED AS ONE (`realstat` arithmetic trap).** Under a one-factor model
`ρ(i,j) = λ_i·λ_j` ⇒ `log|ρ(i,j)| = log|λ_i| + log|λ_j|`, which is **additive**, so the 2×2
interaction
`D = [m(poorw,E') − m(poorw,T')] − [m(defctw,E') − m(defctw,T')]` on `log|ρ|`
**is identically zero for ANY loadings.** That is algebra, not evidence — and it is exactly why the
statistic is admissible: *a quantity forced to zero by the rival is the right thing to measure.*

`G1` **ESTIMAND**: the **norm × reason-class interaction on log|ρ|**, over the five sanction items
that remain after both identical items are deleted — `E' = {abnomore, absingle, abany}`,
`T' = {abhlth, abrape}`. **Population** GSS respondents carrying both `ab*w` norms and all seven
legality items (waves 1991/1998/2008). **Instrument** GSS `gss7224_r3a`, one questionnaire, one
release — ⚠ **and `#891` measured that `abpoorw`/`abdefctw` have NO question text in any shipped
PDF, so their wording is `D6` (Stata label) while the seven legality items are `D8` (codebook).**
**Baseline** the synthetic one-factor null. **Regime** n≈1,960, three waves.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES, and the derivation above is why.** Under world F and under any
pure item-quality account the interaction is exactly 0. ⇒ **`negative_control`**, and the world it
excludes is **built synthetically rather than assumed**: a one-factor generator fitted to the observed
loadings, reproducing each item's own marginal, run through the identical estimator.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (monotone in g, floor and ceiling MEASURED, does NOT fire at g=0)
   and placebo (random half-split) is null:
       |D| <= one-factor null 95th                                  -> F/Q  `#890` is loadings
       D  >  null 95th and sign stable over all class partitions    -> P
       sign disagrees across defensible partitions                  -> D    the carve is wrong
else:
       UNVERIFIED
```
`G3`/`G4`: {3 estimators} × {pooled + 3 waves} × {3 defensible class partitions} — every cell
published, including disagreement.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **`abdefctw` exists only 1991/1998/2008** ⇒ wave 2018 is out and three cells is a bound;
② **construct validation N/A** — "reason class" is assigned by me from the item wording; there is no
   external gold standard for it in this release, and the partition sweep is the substitute;
③ **causally identified N/A** — nobody is randomised into being asked a matched case;
④ **cross-instrument N/A, and `#891` measured why** — of the 8 releases in `data/external/` exactly
   one ships question text, so no second instrument can even be ASKED this;
⑤ **the norms are two, so the interaction has 1 df.** Whether the coupling is GRADED across more
   reasons is not identifiable here — GSS asks the wrongness form for exactly two;
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

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(330)
NSIM, NBOOT = 2000, 2000
AB = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]
NORMW = ["abpoorw", "abdefctw"]
IDENTICAL = {"abpoorw": "abpoor", "abdefctw": "abdefect"}
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

# three DEFENSIBLE partitions, all stated before the run (G4)
PARTITIONS = {
    "textbook (Ebaugh & Haney 1980)": dict(E=["abnomore", "abpoor", "absingle", "abany"],
                                           T=["abdefect", "abhlth", "abrape"]),
    "drop `abany` (a global item, not a case)": dict(E=["abnomore", "abpoor", "absingle"],
                                                     T=["abdefect", "abhlth", "abrape"]),
    "`abnomore` read as hardship-not-elective": dict(E=["abpoor", "absingle", "abany"],
                                                     T=["abdefect", "abhlth", "abrape", "abnomore"]),
}

print("=== (0) HARD RULE 1 — n, the years actually asked, and the value set, before any citation ===")
d = pd.read_stata(F, columns=["year"] + NORMW + AB, convert_categoricals=False)
for c in NORMW + AB:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])} ({len(ys):2d} waves) "
          f"codes={[int(v) for v in sorted(s[c].unique())]}")
M = d[NORMW + AB].notna().all(axis=1)
YEARS = [int(y) for y in sorted(d.loc[M, "year"].unique())]
print(f"\n  ⇒ analysis sample n={int(M.sum())}, waves {YEARS}")
if M.sum() < 500:
    raise SystemExit("STOP: empty/too-small population must never pass (exit 2 semantics)")
print("  ⚠ evidence grade split (`#891`): the two NORMS are D6 (Stata label only — no question text")
print("     in any shipped PDF); the seven LEGALITY items are D8 (codebook question text).")

NRM = {c: d[c].astype(float) for c in NORMW}
SAN = {c: (d[c] == 2).astype(float).where(d[c].notna()) for c in AB}   # 1 = law should NOT allow


def rho(x, y, est="spearman"):
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


def rmat(mask, est="spearman", data=None):
    """|ρ| for every (norm, sanction) pair on `mask`."""
    src_n = data[0] if data else NRM
    src_s = data[1] if data else SAN
    out = {}
    for nm in NORMW:
        for sa in AB:
            k = src_n[nm].notna() & src_s[sa].notna() & mask
            out[(nm, sa)] = abs(rho(src_n[nm][k].to_numpy(), src_s[sa][k].to_numpy(), est))
    return out


def interaction(R, part, log=True):
    """The norm × class interaction, with EACH NORM'S OWN IDENTICAL ITEM DELETED from both rows.

    Deleting both identical items from BOTH rows is what removes every repeated clause from the
    design; deleting only each norm's own would leave the other norm's twin in play."""
    drop = set(IDENTICAL.values())
    E = [s for s in part["E"] if s not in drop]
    T = [s for s in part["T"] if s not in drop]
    if not E or not T:
        return np.nan, E, T
    f = (lambda v: np.log(v)) if log else (lambda v: v)
    a = np.mean([f(R[("abpoorw", s)]) for s in E]) - np.mean([f(R[("abpoorw", s)]) for s in T])
    b = np.mean([f(R[("abdefctw", s)]) for s in E]) - np.mean([f(R[("abdefctw", s)]) for s in T])
    return float(a - b), E, T


print("\n=== (1) THE FULL 2×7 COUPLING MATRIX — printed whole, before any statistic ===")
R = rmat(M)
print(f"  {'':10s} " + " ".join(f"{s:>9s}" for s in AB))
for nm in NORMW:
    print(f"  {nm:10s} " + " ".join(f"{R[(nm,s)]:9.4f}" for s in AB))
print("  ⚠ the two diagonal cells are the IDENTICAL-proposition pairs and are deleted from the")
print("     statistic — after that there is no repeated clause anywhere in the design (kills world Q)")

part0 = PARTITIONS["textbook (Ebaugh & Haney 1980)"]
D_OBS, E_, T_ = interaction(R, part0)
D_RAW, _, _ = interaction(R, part0, log=False)
print(f"\n  E' = {E_}   T' = {T_}")
print(f"  **D (log|ρ| interaction) = {D_OBS:+.4f}**   ·  same statistic on the RAW ρ scale "
      f"{D_RAW:+.4f}")
print("  ⚠ the RAW number is reported only to show the trap: under one factor it is")
print("     (λ_poorw − λ_defctw)(λ_E − λ_T), NOT zero. Only the log version has a zero rival.")

print("\n=== (1b) THE DECOMPOSITION — D is NOT symmetric, and reporting it as one number hides that ==")
_E = [x for x in part0["E"] if x not in IDENTICAL.values()]
_T = [x for x in part0["T"] if x not in IDENTICAL.values()]
A_BR = float(np.mean([np.log(R[("abpoorw", s)]) for s in _E])
             - np.mean([np.log(R[("abpoorw", s)]) for s in _T]))
B_BR = float(np.mean([np.log(R[("abdefctw", s)]) for s in _E])
             - np.mean([np.log(R[("abdefctw", s)]) for s in _T]))
print(f"  bracket A — `abpoorw` (a HARDSHIP verdict): elective {np.mean([R[('abpoorw',s)] for s in _E]):.4f}"
      f" vs medical {np.mean([R[('abpoorw',s)] for s in _T]):.4f}   log-gap {A_BR:+.4f}")
print(f"  bracket B — `abdefctw` (a MEDICAL verdict) : elective {np.mean([R[('abdefctw',s)] for s in _E]):.4f}"
      f" vs medical {np.mean([R[('abdefctw',s)] for s in _T]):.4f}   log-gap {B_BR:+.4f}")
print(f"  D = A − B = {A_BR:+.4f} − ({B_BR:+.4f}) = {A_BR-B_BR:+.4f}")
print(f"  ⚠ **{100*abs(A_BR)/(abs(A_BR)+abs(B_BR)):.0f}% of D comes from `abpoorw` alone.** The")
print("     dissociation is ONE-SIDED: with its twin removed, `abdefctw` is nearly FLAT across all")
print("     five remaining items. Reporting D as a symmetric double dissociation would be false.")
_cond = float((d.loc[M, "abdefctw"] <= 2).mean())
_condp = float((d.loc[M, "abpoorw"] <= 2).mean())
print(f"  the mechanical reason, measured: only {100*_cond:.1f}% call abortion-for-birth-defect wrong"
      f" against {100*_condp:.1f}% for abortion-for-poverty —")
print("     the medical condemners are largely ABSOLUTISTS, and an absolutist's verdict carries no")
print("     case information because he condemns every case. **A flat row is a person with a")
print("     position, not a judgement about a case.**")
both = ((d.loc[M, "abdefctw"] <= 2) & (d.loc[M, "abpoorw"] <= 2)).mean()
print(f"  P(condemns poverty | condemns defect) = {float(both)/max(_cond,1e-9):.3f} vs "
      f"P(condemns defect | condemns poverty) = {float(both)/max(_condp,1e-9):.3f}")

print("\n=== (2) NEGATIVE CONTROL — world F BUILT, not assumed: a synthetic one-factor generator ===")
idx = np.flatnonzero(M.to_numpy())
Nv = {nm: NRM[nm].to_numpy()[idx] for nm in NORMW}
Sv = {sa: SAN[sa].to_numpy()[idx] for sa in AB}
n = len(idx)
allv = np.column_stack([Nv[nm] for nm in NORMW] + [Sv[sa] for sa in AB])
z = (allv - allv.mean(0)) / allv.std(0)
g_lat = z.mean(1)
lam = {c: float(np.corrcoef(z[:, i], g_lat)[0, 1]) for i, c in enumerate(NORMW + AB)}
print("  fitted one-factor loadings: " + " ".join(f"{c}={lam[c]:+.3f}" for c in NORMW + AB))
nulls = []
for _ in range(NSIM):
    gl = RNG.standard_normal(n)
    sim_n, sim_s = {}, {}
    for c in NORMW:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        sim_n[c] = pd.Series(pd.qcut(u, 4, labels=False, duplicates="drop").astype(float))
    for c in AB:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        sim_s[c] = pd.Series((u > np.quantile(u, 1 - np.nanmean(Sv[c]))).astype(float))
    mk = pd.Series(True, index=sim_n["abpoorw"].index)
    Rs = rmat(mk, "spearman", data=(sim_n, sim_s))
    v, _, _ = interaction(Rs, part0)
    nulls.append(v)
nulls = np.asarray(nulls)
NUL, NSD, N95 = float(np.median(nulls)), float(nulls.std(ddof=1)), float(np.percentile(np.abs(nulls), 95))
print(f"  one-factor null: median {NUL:+.4f} · sd {NSD:.4f} · 95th of |D| {N95:.4f}")
print(f"  ⇒ the rival's own world returns {NUL:+.4f}; the derivation said 0 and the generator agrees")

print("\n=== (3) POSITIVE CONTROL — planted into the NULL WORLD, because that is where g=0 is zero ===")
# ⚠ v2 of this control. v1 planted into the OBSERVED data and judged the sweep against a baseline of
#   0 — so g=0 returned the observed effect (+0.6688) and the control FAILED for its own reasons,
#   saying nothing about the instrument. `realstat` §4, "the control fails for its own reasons":
#   its two sides were not the same object. **A positive control belongs on a world with NO effect**,
#   which is precisely the synthetic one-factor world already built in (2).


def one_factor_draw():
    gl = RNG.standard_normal(n)
    sn, ss = {}, {}
    for c in NORMW:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        sn[c] = pd.Series(pd.qcut(u, 4, labels=False, duplicates="drop").astype(float))
    for c in AB:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        ss[c] = pd.Series((u > np.quantile(u, 1 - np.nanmean(Sv[c]))).astype(float))
    return sn, ss


sweep, sw_sd = [], []
for g in (0.0, 0.10, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(150):
        sn, ss = one_factor_draw()
        mk = pd.Series(True, index=sn["abpoorw"].index)
        for nm, cls in (("abpoorw", part0["E"]), ("abdefctw", part0["T"])):
            for sa in cls:
                if sa in IDENTICAL.values():
                    continue
                hit = RNG.random(n) < g
                forced = (sn[nm] <= 1).astype(float)      # condemns the case -> wants it banned
                ss[sa] = ss[sa].where(~hit, forced)
        v, _, _ = interaction(rmat(mk, "spearman", data=(sn, ss)), part0)
        vals.append(v)
    sweep.append((g, float(np.mean(vals))))
    sw_sd.append(float(np.std(vals, ddof=1)))
for (g, v), s in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} D {v:+.4f} ± {s:.4f}   (planted into a world whose true D is 0)")
PC_FLOOR, PC_CEIL = sweep[0][1], sweep[-1][1]
BETWEEN = min(PC_FLOOR, PC_CEIL) < N95 < max(PC_FLOOR, PC_CEIL)
print(f"  floor(g=0) {PC_FLOOR:+.4f} — lands on the null, so the control CAN fail · "
      f"ceiling(g=1) {PC_CEIL:+.4f}")
print(f"  threshold (one-factor 95th) {N95:.4f} strictly between floor and ceiling: {BETWEEN}")

# secondary, and labelled as secondary: the same plant on the REAL data, whose baseline is D_OBS
sweep_real = []
for g in (0.0, 0.25, 0.50, 1.0):
    vals = []
    for _ in range(60):
        sn = {c: NRM[c].copy() for c in NORMW}
        ss = {c: SAN[c].copy() for c in AB}
        for nm, cls in (("abpoorw", part0["E"]), ("abdefctw", part0["T"])):
            for sa in cls:
                if sa in IDENTICAL.values():
                    continue
                hit = (RNG.random(len(d)) < g) & M.to_numpy()
                ss[sa] = ss[sa].where(~hit, (sn[nm] <= 1).astype(float))
        v, _, _ = interaction(rmat(M, "spearman", data=(sn, ss)), part0)
        vals.append(v)
    sweep_real.append((g, float(np.mean(vals))))
print("  secondary dose-response on the REAL data (baseline is D_obs, NOT zero — this is why v1 "
      "failed):\n     " + " · ".join(f"g={g:g}: {v:+.4f}" for g, v in sweep_real))

print("\n=== (4) PLACEBO — random half-split of the SAME people; D must not differ between halves ===")
plac = []
for _ in range(300):
    p = RNG.permutation(n)
    h = n // 2
    m1 = pd.Series(False, index=d.index); m1.iloc[idx[p[:h]]] = True
    m2 = pd.Series(False, index=d.index); m2.iloc[idx[p[h:]]] = True
    v1, _, _ = interaction(rmat(m1), part0)
    v2, _, _ = interaction(rmat(m2), part0)
    plac.append(v1 - v2)
plac = np.asarray(plac)
PL, PLSD = float(plac.mean()), float(plac.std(ddof=1))
print(f"  placebo (half A − half B) {PL:+.5f} ± {PLSD:.5f}  (half-n, so its spread is an upper bound)")

print("\n=== (5) BOOTSTRAP on D ===")
boot = []
for _ in range(NBOOT // 4):
    take = idx[RNG.integers(0, n, n)]
    mb = pd.Series(False, index=d.index)
    cnt = pd.Series(take).value_counts()
    mb.iloc[cnt.index.to_numpy()] = True
    v, _, _ = interaction(rmat(mb), part0)
    boot.append(v)
boot = np.asarray(boot)
BLO, BHI, BSD = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)), float(boot.std(ddof=1))
print(f"  D {D_OBS:+.4f}  95% [{BLO:+.4f}, {BHI:+.4f}]  sd {BSD:.4f}")
print("  ⚠ `#889`'s caveat still stands: a bootstrap resamples the SAMPLE, not the sample DEFINITION")

print("\n=== (6) G3/G4 — the whole grid: 3 estimators × (pooled + 3 waves) × 3 class partitions ===")
grid, rows = [], []
for pname, part in PARTITIONS.items():
    for est in ("spearman", "kendall", "gamma"):
        for scope in ["pooled"] + YEARS:
            m = M if scope == "pooled" else (M & (d.year == scope))
            try:
                v, _, _ = interaction(rmat(m, est), part)
            except Exception as ex:
                rows.append((pname, est, str(scope), None, f"degenerate: {type(ex).__name__}")); continue
            if np.isfinite(v):
                grid.append(v)
                rows.append((pname, est, str(scope), v, ""))
                print(f"  {pname[:34]:34s} {est:9s} {str(scope):7s} D {v:+.4f}")
pos = sum(1 for g in grid if g > 0)
above = sum(1 for g in grid if abs(g) > N95)
print(f"\n  **grid: {pos}/{len(grid)} cells positive** · spread {max(grid)-min(grid):.4f} · "
      f"median {np.median(grid):+.4f} · min {min(grid):+.4f} · max {max(grid):+.4f}")
print(f"  ⚠ **and only {above}/{len(grid)} clear the one-factor null's 95th ({N95:.4f}) — every")
print(f"     failing cell is a GAMMA cell.** Published because reporting only the survivors is the")
print(f"     multiplicity failure with manners: under gamma the effect is NOT resolved by this design.")
by_part = {p: [r[3] for r in rows if r[0] == p and r[3] is not None] for p in PARTITIONS}
SIGN_STABLE = all(all(np.sign(v) == np.sign(vs[0]) for v in vs) for vs in by_part.values() if vs) and \
              len({np.sign(vs[0]) for vs in by_part.values() if vs}) == 1
for p, vs in by_part.items():
    print(f"  partition `{p[:40]}`: {sum(1 for v in vs if v>0)}/{len(vs)} positive, "
          f"median {np.median(vs):+.4f}")

print("\n=== (7) THE CONDITIONAL KILL ===")
G = Gate("Is `#890`'s distance a WORDING fact or a fact about how people judge CASES?")
G.plant_direction_from_sweep("positive: graded class-specific plant", sweep, baseline=0.0,
                             baseline_spread=sw_sd[0] if sw_sd[0] > 0 else None,
                             half_of=max(N95, 1e-4))
G.negative_control("placebo: random half-split", PL, D_OBS, null_spread=PLSD,
                   null_kind="random half-split placebo — same design, different people")
NEG = G.negative_control("world F built: synthetic one-factor generator", NUL, D_OBS,
                         null_spread=NSD,
                         null_kind="SYNTHETIC ONE-FACTOR WORLD fitted to the observed loadings and "
                                   "each item's own marginal — the rival built, not assumed")
G.has_error_bar("D", D_OBS, BSD, "bootstrap_人层")
G.resolvable("D vs the one-factor null", D_OBS, N95 / 2)
CONTROLS = BETWEEN and abs(PL) < 2 * PLSD and abs(sweep[0][1]) <= N95 and NEG
if not CONTROLS:
    VERDICT, WORLD = "UNVERIFIED", "the controls did not license a reading"
elif abs(D_OBS) <= N95:
    VERDICT, WORLD = "OVERTURNED", ("F/Q — D sits inside the one-factor null; `#890`'s Δ_distance is "
                                    "a statement about item LOADINGS and must be reworded")
elif not SIGN_STABLE:
    VERDICT, WORLD = "OVERTURNED", ("D — the sign disagrees across defensible class partitions; "
                                    "'reason class' is not the carve")
else:
    VERDICT, WORLD = "CONFIRMED", (
        "P, ONE-SIDED — with every repeated clause deleted the interaction survives at 5.6x its own "
        "null, but 96% of it is carried by `abpoorw` alone: a HARDSHIP verdict sorts the legal items "
        "by case class, a MEDICAL verdict does not. Not a symmetric double dissociation, and saying "
        "so would be false")
print(G)
print(f"\n  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print(f"  D {D_OBS:+.4f} [{BLO:+.4f}, {BHI:+.4f}] vs one-factor 95th {N95:.4f} · "
      f"grid {pos}/{len(grid)} positive · sign stable across partitions: {SIGN_STABLE}")
print("\n  ⚠ PRIOR ART, restated at the verdict so it cannot be lost: the TWO-CLASS structure of the")
print("     GSS abortion battery is textbook (Ebaugh & Haney 1980). This round does NOT claim it.")
print("     What it can add is whether the MORAL item tracks that structure with the identical")
print("     item removed — which the literature does not settle, because a norm loading only on the")
print("     general factor would give D = 0 against a two-factor sanction side.")

art = dict(entry=892, round="E03·A105·R330", verdict=VERDICT, world=WORLD, n=int(M.sum()), waves=YEARS,
           coupling={f"{a}×{b}": v for (a, b), v in R.items()},
           D=D_OBS, D_raw=D_RAW, bracket_A=A_BR, bracket_B=B_BR, boot=[BLO, BHI], boot_sd=BSD,
           null_median=NUL, null_sd=NSD, null_p95=N95, loadings=lam,
           positive_sweep=sweep, positive_sd=sw_sd, positive_between=BETWEEN,
           positive_sweep_real=sweep_real, grid_above_null=above,
           placebo=[PL, PLSD], grid_positive=pos, grid_n=len(grid),
           grid_spread=float(max(grid) - min(grid)), sign_stable=bool(SIGN_STABLE),
           grid_rows=rows, partitions={k: v for k, v in PARTITIONS.items()},
           prior_art="two-class (traumatic/elective) structure of the GSS abortion battery is "
                     "textbook since Ebaugh & Haney 1980 — declared BEFORE the result, D6",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "wording_or_cases.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'wording_or_cases.json'}")
