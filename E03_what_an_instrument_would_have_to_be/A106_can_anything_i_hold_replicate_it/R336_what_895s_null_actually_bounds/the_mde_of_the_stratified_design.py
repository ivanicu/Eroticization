r"""#898 · E03·A106·R336 — what `#895`'s null actually bounds, and it is not "no effect"

**TYPE: CLOSURE, labelled before the run.** It opens no world and needs none. **What it protects:**
`#895`'s verdict — *the person-level gradient is inside its stratum-label permutation null on 2 of 3
stratifiers* — which `#895` itself registered as **structurally incomplete**:

> `#895` ⑥ · *"POWER IS NOT MEASURED HERE AND THE NULL IS NOT A DEMONSTRATION OF HOMOGENEITY. The
> design's MDE at these stratum sizes was not computed, so 'inside the null' is a BOUND on the
> gradient, not proof there is none. A real person-level effect smaller than ~0.24 in `T` would be
> invisible to this round."*

**That `~0.24` was an eyeball off the null's 95th percentile, not a power computation** — and the
two are different numbers, because detection needs the effect to clear the threshold *on most
draws*, not on one. **`#897`② registered the gap; this pays it with the machinery `#897` just
built.** *A null without an MDE is a shrug; a null with one is a claim.*

`G1` **ESTIMAND**: the **MDE of `R333`'s stratified design** — the smallest planted person-level
gradient whose `|T|` exceeds that stratifier's own null 95th percentile in **≥80%** of draws.
**Population** the same 1,960 GSS respondents, waves 1991/1998/2008. **Instrument** GSS
`gss7224_r3a`; ⚠ the two `ab*w` norms are `D6` (Stata label only — `#891`), the seven legality items
`D8`; ⚠ **cross-instrument N/A — `no second instrument` and `only this one instrument`**, and
`#897` measured why: SCCS's matching design cannot resolve any effect at all. **Baseline** `R333`'s
own per-stratifier null, **READ from `R333`'s artifact and not retyped** (`#840`'s RULE, and ⚠
naming `#840`'s own scope so nothing is inherited: `#840`/`#838`/`#839`/`#846` were measured on the
`homosex` item ALONE, so only the practice transfers, never their findings). **Regime** three strata
per stratifier, n ≈ 516–902.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES**, for the null this MDE is measured against: under a
population-constant tilt a between-stratum contrast is exactly 0. ⇒ **`negative_control`**, **kind of
null named: `R333`'s STRATUM-LABEL PERMUTATION NULL with stratum sizes preserved**, reused rather
than recomputed so that the MDE is on the same scale as the verdict it qualifies.
⚠ **And the MDE is a DERIVATION from a simulated design, not a measurement of the world** — it can
never say the gradient is absent, only how large one would have to be to have been seen.

**PRE-REGISTERED READING — no kill, because Closure has no branch to kill:**
```
MDE < min(observed T)   -> #895's null is a DEMONSTRATION for every stratifier: any real gradient
                           would have shown, and the ones below it are genuinely small
MDE inside the observed  -> #895's null is a BOUND, and each stratifier must be labelled separately
   T range                 as RESOLVED / UNRESOLVED rather than reported together
MDE > max(observed T)   -> #895 measured NOTHING about people; its "inside the null" is silence,
                           and the sentence on the page must say so
```
`G3`: all 3 stratifiers × the whole g-grid, published including the ones that never reach 80%.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① an MDE is a property of the DESIGN, never of the world;
② **the plant is one particular shape** — a graded case-specific coupling applied to one stratum —
   and a real person-level effect of a different shape would have a different MDE. **This is a
   bound for THIS plant**, and no sweep over plant shapes is run;
③ `attend`/`polviews`/`educ` are proxies, mutually correlated, not three independent tests;
④ **cross-instrument N/A** — `no second instrument`, `only this one instrument`;
⑤ no second coder, no second release, no test–retest.
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

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(336)
AB = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]
NORMW = ["abpoorw", "abdefctw"]
TWINS = {"abpoor", "abdefect"}
E_, T_ = ["abnomore", "absingle", "abany"], ["abhlth", "abrape"]
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
R333 = (ROOT / "E03_what_an_instrument_would_have_to_be/A105_is_it_the_act_or_what_the_two_"
        "questions_ask/R333_does_case_judgement_depend_on_who_you_are/results/"
        "person_or_population.json")

print("=== (0) HARD RULE 1 — n, years, codes, before any column is cited ===")
d = pd.read_stata(F, columns=["year", "attend", "polviews", "educ"] + NORMW + AB,
                  convert_categoricals=False)
M = d[NORMW + AB].notna().all(axis=1)
YEARS = [int(y) for y in sorted(d.loc[M, "year"].unique())]
for c in NORMW + AB + ["attend", "polviews", "educ"]:
    s = d.loc[M, [c]].dropna()
    print(f"  {c:9s} n_in_sample={len(s):5d}/{int(M.sum())} codes={[int(v) for v in sorted(s[c].unique())][:9]}")
print(f"  ⇒ base sample n={int(M.sum())}, waves {YEARS}")
if M.sum() < 500:
    raise SystemExit("STOP: an empty/too-small population must never pass")

PRIOR = json.loads(R333.read_text())
NULL95 = {k: float(v["p95"]) for k, v in PRIOR["null"].items()}
OBS_T = {k: float(v["T"]) for k, v in PRIOR["observed"].items()}
SIZES = {k: [int(x) for x in v["sizes"]] for k, v in PRIOR["observed"].items()}
print("\n  read from `R333`'s artifact (never retyped — `#840`'s RULE):")
for k in NULL95:
    print(f"    {k:9s} observed T {OBS_T[k]:+.4f} · its own null 95th pct {NULL95[k]:.4f} · "
          f"strata n={SIZES[k]}")

NRM = {c: d[c].astype(float) for c in NORMW}
SAN = {c: (d[c] == 2).astype(float).where(d[c].notna()) for c in AB}


def bracket(mask, ss=None):
    ss = ss or SAN
    v = {}
    for s in E_ + T_:
        k = NRM["abdefctw"].notna() & ss[s].notna() & mask
        if k.sum() < 60:
            return np.nan
        r = abs(stats.spearmanr(NRM["abdefctw"][k].to_numpy(), ss[s][k].to_numpy()).statistic)
        if not np.isfinite(r) or r <= 1e-6:
            return np.nan
        v[s] = np.log(r)
    return float(np.mean([v[s] for s in E_]) - np.mean([v[s] for s in T_]))


print("\n=== (1) THE POWER CURVE — a graded person-level plant, detection MEASURED per stratifier ===")
print("  ⚠ the plant is applied to a RANDOM third redrawn each replicate (`#895`'s v3 fix: a")
print("     covariate drawn once gives a control whose spread is exactly zero, and a constant")
print("     cannot be a null). Its size is then read against EACH stratifier's own null.")
idx = np.flatnonzero(M.to_numpy())
n = len(idx)
GRID = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.0)
curve = {}
for g in GRID:
    ts = []
    for _ in range(200):   # ⚠ 60 gave a non-monotone wobble at g=0.10→0.15 (0.0944→0.0813);
        # a control that fails on SIMULATION noise is fixed by measuring better, not by excusing it
        cov = RNG.random(len(d))
        G1 = pd.Series(cov < 1 / 3, index=d.index)
        G3 = pd.Series(cov >= 2 / 3, index=d.index)
        ss = {c: SAN[c].copy() for c in AB}
        for sx in T_:
            hit = pd.Series(RNG.random(len(d)) < g, index=d.index) & G1 & M
            ss[sx] = ss[sx].where(~hit, (NRM["abdefctw"] <= 1).astype(float))
        v1, v3 = bracket(M & G1, ss), bracket(M & G3, ss)
        if np.isfinite(v1) and np.isfinite(v3):
            ts.append(v3 - v1)
    a = np.asarray(ts)
    curve[g] = dict(median=float(np.median(a)), sd=float(a.std(ddof=1)),
                    det={k: float(np.mean(np.abs(a) > NULL95[k])) for k in NULL95})
    print(f"  g={g:<5.2f} |T| median {abs(curve[g]['median']):.4f} ± {curve[g]['sd']:.4f}   detection "
          + " ".join(f"{k}={100*curve[g]['det'][k]:5.1f}%" for k in NULL95))

print("\n=== (2) THE MDE PER STRATIFIER, and what it says about each observed contrast ===")
MDE = {}
for k in NULL95:
    hit = [g for g in GRID if curve[g]["det"][k] >= 0.80]
    MDE[k] = abs(curve[hit[0]]["median"]) if hit else None
    obs = abs(OBS_T[k])
    if MDE[k] is None:
        verdict = "80% power NEVER REACHED on this grid — the design cannot see this stratifier"
    elif obs >= MDE[k]:
        verdict = f"observed |T| {obs:.4f} ≥ MDE {MDE[k]:.4f} ⇒ **RESOLVED** (and `#895` read it as fired)"
    else:
        verdict = (f"observed |T| {obs:.4f} < MDE {MDE[k]:.4f} ⇒ **UNRESOLVED, not absent** — "
                   f"`#895`'s null here is SILENCE, not homogeneity")
    print(f"  {k:9s} MDE " + ("not reached" if MDE[k] is None else f"{MDE[k]:.4f}") + f"  ⇒ {verdict}")

vals = [v for v in MDE.values() if v is not None]
MDE_MED = float(np.median(vals)) if vals else None
obs_all = [abs(v) for v in OBS_T.values()]
print(f"\n  median MDE over stratifiers: " + ("n/a" if MDE_MED is None else f"{MDE_MED:.4f}")
      + f"  ·  observed |T| range [{min(obs_all):.4f}, {max(obs_all):.4f}]")

print("\n=== (3) THE READING — pre-registered, and Closure has no kill ===")
if MDE_MED is None:
    READ = "the design never reaches 80% power on this grid ⇒ `#895` measured NOTHING about people"
elif MDE_MED > max(obs_all):
    READ = ("`#895`'s null is SILENCE for every stratifier — the MDE exceeds every observed contrast, "
            "so no gradient of the observed size could have been seen and the page must say so")
elif MDE_MED < min(obs_all):
    READ = ("`#895`'s null is a DEMONSTRATION — any gradient as large as the smallest observed would "
            "have shown, so the small ones are genuinely small")
else:
    READ = ("`#895`'s null is a BOUND and the stratifiers must be labelled SEPARATELY: the ones above "
            "the MDE are resolved, the ones below are UNRESOLVED and not absent")
print(f"  ⇒ **{READ}**")

G = Gate("What does `#895`'s null actually bound?")
G.plant_direction_from_sweep("positive: the statistic responds to the person-level plant",
                             [(g, curve[g]["median"]) for g in GRID], baseline=curve[0.0]["median"],
                             baseline_spread=curve[0.0]["sd"], half_of=0.15)
G.negative_control("g=0 lands on zero", curve[0.0]["median"], max(obs_all),
                   null_spread=curve[0.0]["sd"],
                   null_kind="`R333`'s STRATUM-LABEL PERMUTATION NULL with sizes preserved, REUSED "
                             "rather than recomputed so the MDE is on the same scale as the verdict "
                             "it qualifies")
print(G)
print(f"  gate three-valued : {G.three_valued()}")
print("\n  ⚠ AN MDE IS A PROPERTY OF THE DESIGN, NOT OF THE WORLD, and this one is a bound for ONE")
print("     PLANT SHAPE — a graded case-specific coupling applied to one stratum. A real effect of a")
print("     different shape would have a different MDE, and no sweep over plant shapes was run.")

art = dict(entry=898, round="E03·A106·R336", type="CLOSURE",
           protects="#895's 'the gradient is inside its null' — its own register item ⑥ said the MDE "
                    "was never computed, so 'inside the null' was a shrug",
           n=int(M.sum()), waves=YEARS, read_from_R333=dict(null95=NULL95, observed_T=OBS_T, sizes=SIZES),
           curve={str(g): curve[g] for g in GRID}, mde=MDE, mde_median=MDE_MED,
           observed_abs_range=[min(obs_all), max(obs_all)], reading=READ,
           caveat="an MDE is a property of the DESIGN, not of the world, and this is a bound for ONE "
                  "plant shape; no sweep over plant shapes was run",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "mde_of_the_stratified_design.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'mde_of_the_stratified_design.json'}")
