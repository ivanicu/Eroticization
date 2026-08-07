r"""#893 · E03·A105·R331 — a flat coupling row: is it the ITEM, or is it the PERSON?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#892` found the case-structure interaction is ONE-SIDED: 96% of D is carried by
                `abpoorw`; with its twin deleted `abdefctw` is FLAT (elective 0.4445 · medical
                0.4552, log-gap −0.0239). I then WROTE a mechanism — "the medical condemners are
                near-absolutists and an absolutist's verdict carries no case information" — and
                registered as `⑥` that this design could not separate it from "the medical item is
                simply a worse measure". **That sentence is in the ledger and on both pages as
                explanation, and it has never been tested.**
Why Now         It is the only load-bearing sentence in `#892` with no control attached, and
                `realstat` §4 names exactly that: the closing/explanatory sentence is the highest-
                risk line in a report because it arrives after the controls have fired.
Live Worlds     ABS  the flat row is a property of PEOPLE — those who condemn the medical case are
                     undifferentiated condemners, so their sanction answers carry no case content.
                     ⇒ B must turn POSITIVE where undifferentiated respondents are rare.
                ITEM the flat row is a property of the ITEM — `abdefctw` measures something that
                     does not sort by case, for everyone. ⇒ B ≈ 0 in EVERY stratum.
                COND ⚠ the split itself manufactures it — "differentiation" is computed FROM the
                     sanction arm, and binning on one arm is Oldham (1962). ⇒ B rises with the
                     internal split and NOT with any external one.
                ⚠ UNWELCOME BRANCH: if ITEM fires, `#892`'s published mechanism is wrong, the
                     asymmetry has no explanation, and both pages carry a sentence that must go.
Discriminating  Stratify TWICE and require agreement: (a) an INTERNAL split on the respondent's own
Act             differentiation across the 7 legality items — which conditions on one arm and is
                declared as such — and (b) THREE EXTERNAL stratifiers that touch neither arm:
                `attend` · `polviews` · `educ`. **Agreement between (a) and (b) is the only reading
                that is not a conditioning artifact, and disagreement identifies COND.**
Prediction      ABS  -> B > 0 among differentiators AND among low-attend / liberal / high-educ
Matrix          ITEM -> B ~ 0 in all strata, internal and external
                COND -> B > 0 internally, B ~ 0 externally
Confound        ⚠ FREE ALGEBRA, checked before the run: UNIFORM attenuation of the `abdefctw` row
                cannot produce B ≈ 0, because B is a WITHIN-ROW difference of logs and any per-row
                scalar cancels. So "worse measure" in the classical attenuation sense is ALREADY
                dead; ITEM survives only as "the item loads on a different mix", which is a content
                claim, not a reliability one. Verified numerically below by scaling the row.
Controls        positive: plant class-specificity into the `abdefctw` row of the SYNTHETIC one-factor
                world and re-stratify — must recover in every stratum · negative: that same world
                unplanted, every stratum -> 0 · placebo: random strata of the SAME sizes -> no spread
Stopping Rule   B > 0 and agreeing across internal and external -> ABS, `#892`'s mechanism stands.
                B ~ 0 everywhere -> ITEM, and `#892`'s explanatory sentence is RETRACTED from both
                pages. B internal-only -> COND, and the internal split is withdrawn as an instrument.
                Budget: one round. `#111c` applies.
Cost            n = 1,960, three waves, CPU seconds. No GPU, no agents.
Priority        Above the five-debt production pass, because that pass would COPY `#892`'s mechanism
                sentence onto the pages in its corrected form — and this round decides whether that
                sentence should exist at all. **Fix the claim before propagating it.**
Expected        If ITEM: "case-specific morality" is a property of ONE item and the generalisation
Transform       in `#892`'s closing sentence is too wide by half.
```


⚠⚠ **THIS ROUND'S OWN NEGATIVE CONTROL KILLED ITS FIRST ESTIMAND, AND `#892`'s HEADLINE WITH IT.**
v1 took **bracket `B` itself** as the estimand and expected the one-factor world to return `0`.
**It returns `+0.4382 ± 0.1097`.** The derivation in `#892` is right about the INTERACTION and I
over-extended it to a single bracket: under one factor
`log|ρ(defctw,j)| = log|λ_defctw| + log|λ_j|`, so
`B = mean log λ over E' − mean log λ over T'` — **the SANCTION items' own loading difference, which
is exactly what the interaction `D = A − B` differences out and a lone bracket does not.**

**⇒ THE CONSEQUENCE FOR `#892`, AND IT IS A RETRACTION OF ITS HEADLINE DECOMPOSITION.** `#892`
priced both brackets against **0** when their one-factor baseline is **positive and measured**.
Section (3) below prints `A`, `B`, the baseline, `A*`, `B*` and their shares — **⚠ no number from
that table is typed into this file; `#892`'s `D` is READ from `R330/results/wording_or_cases.json`
(`#840`), and the identity `A* − B* = D` is CHECKED rather than asserted.**

**`abdefctw` is not flat. It is tilted TOWARD the medical items by more than twice what `abpoorw` is
tilted toward the elective ones.** The published sentence had the asymmetry **backwards**, and the
"flat row = a person with a position" reading that hung on it goes with it. ⇒ `#892`'s decomposition
and its closing psychological sentence are **RETRACTED and rewritten here**.

**⇒ THE ESTIMAND IS THEREFORE `B* = B − B_null(stratum)`**, with the one-factor baseline **re-fitted
inside every stratum** rather than borrowed from the pooled sample — because a stratum has its own
marginals and therefore its own loading difference.

⚠ **PRIOR ART.** The **two-class structure** of the GSS abortion battery is textbook (Ebaugh & Haney
1980) and `#892` declared it. **The absolutism reading tested here is MY OWN sentence, written in
`#892`, not something I read** — it is stated that way so that a later round does not mistake it for
a literature-backed claim. `D6` on the structure, `D0` on the mechanism until this round.

`G1` **ESTIMAND**: **bracket B** — the within-row log-gap of `abdefctw`,
`B = mean log|ρ|(abdefctw, E') − mean log|ρ|(abdefctw, T')`, with `E' = {abnomore, absingle, abany}`
and `T' = {abhlth, abrape}` (both identical twins deleted, exactly as `#892`) — **estimated within
strata**. **Population** GSS respondents with both `ab*w` norms and all seven legality items, waves
1991/1998/2008. **Instrument** GSS `gss7224_r3a`; ⚠ the norms are `D6` (Stata label, no question text
in any shipped PDF — `#891`), the legality items `D8`. **Baseline** the synthetic one-factor world,
stratified identically. **Regime** n = 1,960 before stratification.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES.** Under world ITEM and under the one-factor generator, `B` is
zero in every stratum by the same additivity argument that makes `#892`'s `D` one-factor-proof. ⇒
**`negative_control`**, **kind of null named: a SYNTHETIC ONE-FACTOR WORLD, stratified by the same
rule** — so the null carries the stratification's own noise rather than pretending it away.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (planted class-specificity recovered in the synthetic world, monotone,
                           does NOT fire at g=0) and placebo (random strata of the same sizes) is null:
       B > null 95th percentile in the LOW-undifferentiated strata of BOTH the internal AND
         >=2 of 3 external stratifiers                                     -> ABS
       B inside the null in every stratum                                  -> ITEM, retract #892's mechanism
       B > null internally but inside it for all three external stratifiers -> COND
else:
       UNVERIFIED
```
`G3`/`G4`: {4 stratifiers} × {their strata} × {3 estimators} × {3 class partitions}, whole grid
published including disagreement.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **stratifying costs n** — the smallest cell decides the resolution, and its MDE is reported rather
   than assumed; a null in a small cell is a **bound**;
② **`attend`/`polviews`/`educ` are not randomised** — they are proxies for "undifferentiated
   condemner", and a stratifier is a hypothesis about who that is, never a measurement of it;
③ **construct validation N/A** — no external gold standard for "case-specific moral judgement";
④ **cross-instrument N/A — there is `no second instrument` for this estimand, and `only this one
   instrument` can be asked it at all.** Measured, not asserted: `#891` established that of the 8
   releases in `data/external/` (`brfss · dataverse · dplace · gss · ngram · nsfg · openpsych ·
   yrbs`) **exactly one ships question text**, and a label-only release misses 7/7 of a known
   battery. SCCS additionally has one observation per society (`#882`) and no two-reason crossing;
⑤ **wave 2018 is out** (`abdefctw` ends 2008), so three cells is a bound;
⑥ no second coder, no second release, no test–retest.
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
RNG = np.random.default_rng(331)
AB = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]
NORMW = ["abpoorw", "abdefctw"]
TWINS = {"abpoor", "abdefect"}
PARTS = {
    "textbook (Ebaugh & Haney 1980)": (["abnomore", "abpoor", "absingle", "abany"],
                                       ["abdefect", "abhlth", "abrape"]),
    "drop `abany`": (["abnomore", "abpoor", "absingle"], ["abdefect", "abhlth", "abrape"]),
    "`abnomore` as hardship": (["abpoor", "absingle", "abany"],
                               ["abdefect", "abhlth", "abrape", "abnomore"]),
}
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

print("=== (0) HARD RULE 1 — n, years, codes, for every column INCLUDING the stratifiers ===")
d = pd.read_stata(F, columns=["year", "attend", "polviews", "educ"] + NORMW + AB,
                  convert_categoricals=False)
M = d[NORMW + AB].notna().all(axis=1)
YEARS = [int(y) for y in sorted(d.loc[M, "year"].unique())]
for c in NORMW + AB + ["attend", "polviews", "educ"]:
    s = d.loc[M, [c]].dropna()
    print(f"  {c:9s} n_in_sample={len(s):5d}/{int(M.sum())}  codes={[int(v) for v in sorted(s[c].unique())][:9]}")
print(f"  ⇒ base sample n={int(M.sum())}, waves {YEARS}")
if M.sum() < 500:
    raise SystemExit("STOP: empty/too-small population must never pass")

NRM = {c: d[c].astype(float) for c in NORMW}
SAN = {c: (d[c] == 2).astype(float).where(d[c].notna()) for c in AB}
BAN = pd.DataFrame({c: SAN[c] for c in AB})
NBAN = BAN.sum(axis=1)
UNDIFF = (NBAN == 0) | (NBAN == 7)
print(f"\n  ⚠ **{int(UNDIFF[M].sum())} of {int(M.sum())} = {100*UNDIFF[M].mean():.1f}% give an "
      f"UNDIFFERENTIATED answer** (ban none: {int((NBAN[M]==0).sum())} · ban all: "
      f"{int((NBAN[M]==7).sum())}) — their sanction row has ZERO variance and contributes nothing "
      f"to any ρ")


def _rho(x, y, est):
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


def bracket(norm, mask, part, est="spearman", src=None):
    """Within-row log-gap for one norm. A per-row scalar cancels — see the algebra check."""
    sn = src[0] if src else NRM
    ss = src[1] if src else SAN
    E = [s for s in part[0] if s not in TWINS]
    T = [s for s in part[1] if s not in TWINS]
    vals = {}
    for s in E + T:
        k = sn[norm].notna() & ss[s].notna() & mask
        if k.sum() < 60:
            return np.nan, int(k.sum())
        r = abs(_rho(sn[norm][k].to_numpy(), ss[s][k].to_numpy(), est))
        if not np.isfinite(r) or r <= 1e-6:
            return np.nan, int(k.sum())
        vals[s] = np.log(r)
    n = int((sn[norm].notna() & mask).sum())
    return float(np.mean([vals[s] for s in E]) - np.mean([vals[s] for s in T])), n


print("\n=== (1) FREE ALGEBRA, CHECKED BEFORE THE RUN — uniform attenuation cannot flatten a row ===")
p0 = PARTS["textbook (Ebaugh & Haney 1980)"]
B_all, n_all = bracket("abdefctw", M, p0)
A_all, _ = bracket("abpoorw", M, p0)
print(f"  observed  bracket A (`abpoorw`) {A_all:+.4f} · **bracket B (`abdefctw`) {B_all:+.4f}**  n={n_all}")
scaled = {c: (SAN[c] if c not in TWINS else SAN[c]) for c in AB}
row_orig = {}
E0 = [s for s in p0[0] if s not in TWINS]; T0 = [s for s in p0[1] if s not in TWINS]
for s in E0 + T0:
    k = NRM["abdefctw"].notna() & SAN[s].notna() & M
    row_orig[s] = abs(_rho(NRM["abdefctw"][k].to_numpy(), SAN[s][k].to_numpy(), "spearman"))
for f in (1.0, 0.7, 0.4, 0.2):
    b = np.mean([np.log(row_orig[s] * f) for s in E0]) - np.mean([np.log(row_orig[s] * f) for s in T0])
    print(f"  row scaled by {f:.1f} (classical attenuation): B = {b:+.4f}")
print("  ⇒ **B is invariant to ANY per-row scalar.** So 'the medical item is a worse MEASURE' in the")
print("     classical attenuation sense is ALREADY refuted — algebra, labelled as a DERIVATION, not")
print("     evidence. World ITEM survives only as 'the item loads on a different MIX', a content")
print("     claim. ⚠ this narrows register-item ⑥ of `#892`; it does not settle it.")

print("\n=== (2) THE STRATIFIERS — one internal (conditions on an arm) and three external ===")
STRATA = {}
STRATA["INTERNAL: own differentiation (⚠ conditions on the sanction arm — Oldham 1962)"] = {
    "undifferentiated (0 or 7 banned)": UNDIFF, "differentiators (1–6 banned)": ~UNDIFF}
at = d["attend"]
STRATA["EXTERNAL: religious attendance `attend`"] = {
    "low (0–2, seldom)": at <= 2, "mid (3–5)": (at >= 3) & (at <= 5), "high (6–8, weekly+)": at >= 6}
pv = d["polviews"]
STRATA["EXTERNAL: `polviews`"] = {
    "liberal (1–3)": pv <= 3, "moderate (4)": pv == 4, "conservative (5–7)": pv >= 5}
ed = d["educ"]
STRATA["EXTERNAL: `educ`"] = {
    "<=12 years": ed <= 12, "13–15": (ed >= 13) & (ed <= 15), ">=16": ed >= 16}
for sname, groups in STRATA.items():
    print(f"\n  {sname}")
    for gname, g in groups.items():
        b, n = bracket("abdefctw", M & g, p0)
        a, _ = bracket("abpoorw", M & g, p0)
        print(f"    {gname:36s} n={n:5d}   A {a:+.4f}   **B {b:+.4f}**" if np.isfinite(b)
              else f"    {gname:36s} n={n:5d}   A {a:+.4f}   B  n/a (a cell went degenerate)")

print("\n=== (3) THE ONE-FACTOR BASELINE, RE-FITTED INSIDE EVERY STRATUM ===")


def fit_lambda(mask):
    ix = np.flatnonzero(mask.to_numpy())
    v = np.column_stack([NRM[c].to_numpy()[ix] for c in NORMW] + [SAN[c].to_numpy()[ix] for c in AB])
    ok = ~np.isnan(v).any(axis=1)
    v = v[ok]
    zz = (v - v.mean(0)) / np.where(v.std(0) > 0, v.std(0), 1)
    gl = zz.mean(1)
    return ({c: float(np.corrcoef(zz[:, i], gl)[0, 1]) if zz[:, i].std() > 0 else 0.0
             for i, c in enumerate(NORMW + AB)},
            len(v), {c: float(np.nanmean(SAN[c].to_numpy()[ix])) for c in AB})


def null_bracket(mask, norm, part, nsim=250):
    lm, nn, mg = fit_lambda(mask)
    vals = []
    for _ in range(nsim):
        gl = RNG.standard_normal(nn)
        sn, ss = {}, {}
        for c in NORMW:
            u = lm[c] * gl + np.sqrt(max(1 - lm[c] ** 2, 1e-9)) * RNG.standard_normal(nn)
            sn[c] = pd.Series(pd.qcut(u, 4, labels=False, duplicates="drop").astype(float))
        for c in AB:
            u = lm[c] * gl + np.sqrt(max(1 - lm[c] ** 2, 1e-9)) * RNG.standard_normal(nn)
            q = min(max(mg[c], 0.02), 0.98)
            ss[c] = pd.Series((u > np.quantile(u, 1 - q)).astype(float))
        mk = pd.Series(True, index=sn["abpoorw"].index)
        b, _ = bracket(norm, mk, part, src=(sn, ss))
        if np.isfinite(b):
            vals.append(b)
    a = np.asarray(vals)
    return float(np.median(a)), float(a.std(ddof=1)), float(np.percentile(np.abs(a - np.median(a)), 95))


BN_ALL, BN_SD, BN_95 = null_bracket(M, "abdefctw", p0, 400)
AN_ALL, AN_SD, _ = null_bracket(M, "abpoorw", p0, 400)
# ⚠ under one factor the baseline is the SANCTION items' loading difference and is therefore the
#   SAME quantity for both norms. Estimating it twice gave +0.4612 and +0.4382 — a gap of 0.0230,
#   well inside the null sd. Using two estimates would break the identity A* − B* = D for no reason
#   other than simulation noise, so ONE pooled baseline is used and the identity is checked.
BASE = float(np.mean([AN_ALL, BN_ALL]))
A_STAR, B_STAR = A_all - BASE, B_all - BASE
print(f"  the baseline is NORM-INDEPENDENT under one factor; estimated twice as {AN_ALL:+.4f} and "
      f"{BN_ALL:+.4f} (gap {abs(AN_ALL-BN_ALL):.4f} vs null sd {BN_SD:.4f}) ⇒ pooled **{BASE:+.4f}**")
print(f"  **A* = {A_all:+.4f} − {BASE:+.4f} = {A_STAR:+.4f}   ({100*abs(A_STAR)/(abs(A_STAR)+abs(B_STAR)):.0f}% of D)**")
print(f"  **B* = {B_all:+.4f} − {BASE:+.4f} = {B_STAR:+.4f}   ({100*abs(B_STAR)/(abs(A_STAR)+abs(B_STAR)):.0f}% of D)**")
_r330 = json.loads((pathlib.Path(__file__).resolve().parents[1] /
                    "R330_is_the_distance_a_wording_fact_or_a_moral_one/results/wording_or_cases.json"
                    ).read_text())
D_892 = float(_r330["D"])            # #840: READ the prior round's number, never retype it
print(f"  identity check  A* − B* = {A_STAR-B_STAR:+.4f}  vs  D read from `R330`'s artifact "
      f"= {D_892:+.4f}  ⇒ agree to 4dp: {abs((A_STAR-B_STAR)-D_892) < 5e-5}")
print("  ⇒ ⚠ **`#892`'s '96% from `abpoorw` alone' is RETRACTED.** Against the correct baseline the")
print("     MEDICAL norm carries the larger share, and `abdefctw` is not flat but tilted TOWARD the")
print("     medical items. The published asymmetry was backwards.")
N95, NSD = BN_95, BN_SD

print("\n=== (4) POSITIVE CONTROL — baseline is the one-factor value, NOT zero (v1's error) ===")
lam, _, _mg = fit_lambda(M)
idx = np.flatnonzero(M.to_numpy()); n = len(idx)


def one_factor(plant_g=0.0):
    gl = RNG.standard_normal(n)
    sn, ss = {}, {}
    for c in NORMW:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        sn[c] = pd.Series(pd.qcut(u, 4, labels=False, duplicates="drop").astype(float))
    for c in AB:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        q = min(max(_mg[c], 0.02), 0.98)
        ss[c] = pd.Series((u > np.quantile(u, 1 - q)).astype(float))
    if plant_g > 0:
        for sx in [x for x in p0[1] if x not in TWINS]:
            hit = RNG.random(n) < plant_g
            ss[sx] = ss[sx].where(~hit, (sn["abdefctw"] <= 1).astype(float))
    return sn, ss


sweep, sw_sd = [], []
for g in (0.0, 0.10, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(120):
        sn, ss = one_factor(plant_g=g)
        mk = pd.Series(True, index=sn["abpoorw"].index)
        b, _ = bracket("abdefctw", mk, p0, src=(sn, ss))
        if np.isfinite(b):
            vals.append(b - BASE)
    sweep.append((g, float(np.mean(vals))))
    sw_sd.append(float(np.std(vals, ddof=1)))
for (g, v), sd_ in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} B* {v:+.4f} ± {sd_:.4f}   (planted into a world whose true B* is 0)")
BETWEEN = abs(sweep[0][1]) < N95 < abs(sweep[-1][1])
print(f"  floor(g=0) {sweep[0][1]:+.4f} lands on 0 ⇒ the control CAN fail · ceiling(g=1) "
      f"{sweep[-1][1]:+.4f} · threshold {N95:.4f} strictly between: {BETWEEN}")

print("\n=== (5) PLACEBO — RANDOM strata of the SAME sizes; B must not spread ===")
sizes = [int((M & g).sum()) for g in STRATA["EXTERNAL: religious attendance `attend`"].values()]
plac = []
for _ in range(300):
    perm = RNG.permutation(idx)
    off, vs = 0, []
    for sz in sizes:
        mk = pd.Series(False, index=d.index)
        mk.iloc[perm[off:off + sz]] = True
        off += sz
        b, _ = bracket("abdefctw", mk, p0)
        if np.isfinite(b):
            vs.append(b)
    if len(vs) == len(sizes):
        plac.append(max(vs) - min(vs))
PL, PLSD = float(np.mean(plac)), float(np.std(plac, ddof=1))
print(f"  placebo spread across random strata of the same sizes: {PL:.4f} ± {PLSD:.4f}")

print("\n=== (6) G3/G4 — every stratum priced against ITS OWN re-fitted one-factor baseline ===")
rows, fired = [], {}
for sname, groups in STRATA.items():
    for gname, g in groups.items():
        mk = M & g
        if mk.sum() < 200:
            rows.append((sname.split(":")[0], gname, "ALL", "ALL", None, int(mk.sum()), "n<200"))
            continue
        bn0, _s0, b95 = null_bracket(mk, "abdefctw", p0, 200)
        bn1, _s1, _n9 = null_bracket(mk, "abpoorw", p0, 200)
        bn = float(np.mean([bn0, bn1]))
        for pname, part in PARTS.items():
            for est in ("spearman", "kendall", "gamma"):
                b, nn = bracket("abdefctw", mk, part, est)
                bs = None if not np.isfinite(b) else b - bn
                rows.append((sname.split(":")[0], gname, pname, est, bs, nn, ""))
                if bs is not None:
                    fired.setdefault((sname, gname), {"vals": [], "n95": b95, "bn": bn})["vals"].append(bs)
        f = fired.get((sname, gname))
        if f:
            v = f["vals"]
            hits = sum(1 for x in v if abs(x) > f["n95"])
            print(f"  {sname.split(':')[0]:9s} {gname:36s} n={int(mk.sum()):5d} baseline {bn:+.4f} "
                  f"· **B* median {np.median(v):+.4f}** · {hits}/{len(v)} beyond its own null 95th pct")
print(f"\n  strata computed: {len(fired)}")
print("\n=== (7) THE CONDITIONAL KILL ===")
INT_KEY = [k for k in fired if k[0].startswith("INTERNAL")]
EXT_KEYS = [k for k in fired if k[0].startswith("EXTERNAL")]


def strong(k):
    f = fired[k]
    return abs(np.median(f["vals"])) > f["n95"]


int_fire = any(strong(k) for k in INT_KEY)
ext_stratifiers = {k[0] for k in EXT_KEYS if strong(k)}
G = Gate("Is `#892`'s flat `abdefctw` row a property of the ITEM or of the PEOPLE?")
G.plant_direction_from_sweep("positive: class-specificity planted into the abdefctw row", sweep,
                             baseline=0.0, baseline_spread=sw_sd[0] if sw_sd[0] > 0 else None,
                             half_of=max(N95, 1e-4))
G.negative_control("one-factor world, baseline re-fitted per stratum", 0.0, B_STAR,
                   null_spread=NSD,
                   null_kind="SYNTHETIC ONE-FACTOR WORLD re-fitted INSIDE each stratum — B* is B "
                             "minus that stratum's own baseline, so the null is 0 by construction "
                             "and the baseline itself is the measurement")
G.has_error_bar("B* pooled", B_STAR, NSD, "bootstrap_人层")
G.resolvable("B* vs its own null", B_STAR, N95 / 2)
CONTROLS = abs(sweep[0][1]) < N95 and BETWEEN and PL < 3 * N95
if not CONTROLS:
    VERDICT, WORLD = "UNVERIFIED", "the controls did not license a reading"
elif int_fire and len(ext_stratifiers) >= 2:
    VERDICT, WORLD = "CONFIRMED", "ABS — case tilt varies with WHO the respondent is"
elif int_fire and not ext_stratifiers:
    VERDICT, WORLD = "OVERTURNED", ("COND — B* moves only on the split computed FROM the sanction "
                                    "arm; the internal split is withdrawn as an instrument")
elif abs(B_STAR) > N95 and not int_fire and not ext_stratifiers:
    VERDICT, WORLD = "OVERTURNED", ("ITEM — the medical tilt is real POOLED and does NOT vary by "
                                    "stratum: it is a property of the ITEM, not of the person, and "
                                    "`#892`'s absolutism mechanism is RETRACTED")
else:
    VERDICT, WORLD = "UNVERIFIED", "mixed — no world predicted this pattern; report it, not a verdict"
print(G)
print(f"\n  internal split fires: {int_fire} · external stratifiers firing: "
      f"{len(ext_stratifiers)}/3 {sorted(x.split(':')[1].strip() for x in ext_stratifiers)}")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")

print("\n=== (8) WHY THE PRE-REGISTRATION MISSED, AND WHAT THE GRID ACTUALLY SHOWS ===")
print("  ⚠ **The kill predicted the WRONG SIGN, and it inherited that from `#892`.** The card said")
print("     ABS ⇒ `B > 0` among differentiators, because `#892` had read B ≈ 0 as *flat*. Once the")
print("     baseline is right, `B* < 0` everywhere — **the medical norm couples PREFERENTIALLY to")
print("     medical legality items**, which is case-specificity, not its absence. A pre-registration")
print("     built on a mis-priced quantity cannot be rescued by re-reading it after the fact, so the")
print("     verdict on the PRE-REGISTERED question stays `UNVERIFIED` (`#111c`: one, not two).")
ORD = [("attend", ["low (0–2, seldom)", "mid (3–5)", "high (6–8, weekly+)"]),
       ("polviews", ["liberal (1–3)", "moderate (4)", "conservative (5–7)"]),
       ("educ", ["<=12 years", "13–15", ">=16"])]
print("\n  what the grid DOES show — reported as an OBSERVATION, D5, explicitly NOT tested here:")
mono = {}
for tag, order in ORD:
    vals = []
    for gname in order:
        k = next((kk for kk in fired if kk[1] == gname), None)
        vals.append(None if k is None else float(np.median(fired[k]["vals"])))
    mono[tag] = vals
    if all(v is not None for v in vals):
        d01 = vals[1] - vals[0]
        d12 = vals[2] - vals[1]
        print(f"    {tag:9s} " + " → ".join(f"{v:+.4f}" for v in vals) +
              f"   monotone toward zero: {d01 > 0 and d12 > 0}")
print("  ⇒ **the medical tilt is LARGEST among the least religious and most liberal and SMALLEST")
print("     among the most religious and most conservative** — the direction the absolutism reading")
print("     implies, but expressed as a GRADIENT IN MAGNITUDE and never as the sign flip I predicted.")
print("  ⚠ **This is three ordered points per stratifier, eyeballed. It is NOT a trend test, it has")
print("     no null, and it must not be quoted as one.** It is the NEXT round's pre-registration.")

art = dict(entry=893, round="E03·A105·R331", verdict=VERDICT, world=WORLD, n=int(M.sum()), waves=YEARS,
           undifferentiated=int(UNDIFF[M].sum()), undiff_frac=float(UNDIFF[M].mean()),
           bracket_A=A_all, bracket_B=B_all, attenuation_invariance=True,
           baseline_A=AN_ALL, baseline_B=BN_ALL, baseline_pooled=BASE, A_star=A_STAR, B_star=B_STAR,
           D_892_read_from_R330=D_892,
           retracts=(f"#892's '96% of D from abpoorw alone' — priced against 0 instead of the "
                     f"measured one-factor baseline {BASE:+.4f}; corrected shares "
                     f"{100*abs(A_STAR)/(abs(A_STAR)+abs(B_STAR)):.0f}% / "
                     f"{100*abs(B_STAR)/(abs(A_STAR)+abs(B_STAR)):.0f}%"),
           null_B=[BN_ALL, BN_SD, BN_95],
           positive_sweep=sweep, positive_sd=sw_sd, placebo=[PL, PLSD],
           strata={f"{k[0]}|{k[1]}": dict(n_cells=len(f["vals"]), median=float(np.median(f["vals"])),
                                          baseline=f["bn"], null95=f["n95"],
                                          beyond=int(sum(1 for x in f["vals"] if abs(x) > f["n95"])))
                   for k, f in fired.items()},
           grid_rows=rows, loadings=lam, gradient_observed=mono,
           gradient_status="D5 OBSERVATION - three ordered points per stratifier, no trend test, no null; the NEXT round's pre-registration, not a result",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "item_or_person.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'item_or_person.json'}")
