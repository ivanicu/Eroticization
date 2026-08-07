r"""#895 · E03·A105·R333 — is judging a CASE something only some people do?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#893` corrected `#892` and left `B* = −0.4736` — the medical verdict couples
                preferentially to medical legality items — and then noticed, WITHOUT TESTING IT,
                that |B*| orders monotonically on all three external stratifiers: it is largest
                among the least religious and most liberal, smallest among the most religious and
                most conservative. `#893` labelled that `D5`, "three ordered points, eyeballed, no
                trend test, no null, must not be quoted as one", and registered it as `#893`①.
Why Now         It is the only open question in this arc that is about PEOPLE rather than about
                items, and it is the last one the instrument can still answer. Everything else open
                is production or is registered impossible.
Live Worlds     PERSON  whether a moral verdict is a JUDGEMENT ABOUT A CASE or a POSITION varies
                        with who holds it. ⇒ the tilt shrinks monotonically toward the
                        undifferentiated end, beyond a label-permutation null.
                POP     the tilt is a population constant; the three ordered points are what three
                        noisy estimates look like. ⇒ contrast inside its null.
                ⚠ UNWELCOME: POP is the one I do not want, because `#893`'s closing paragraph and
                        `#892`'s absolutism story both lean on the gradient being real, and if it
                        is noise then this arc has exactly one person-level finding left: none.
                SIZE    ⚠ META — the gradient is n, not people: the strata differ in size and the
                        estimator's own bias moves with n. ⇒ the contrast tracks stratum n rather
                        than stratum content, and the whole person-level question is ill-posed here.
Discriminating  The contrast `T = B*(last stratum) − B*(first stratum)` per stratifier, against a
Act             null that PERMUTES STRATUM LABELS among respondents while preserving stratum SIZES.
                That null kills the person↔coupling association and keeps everything else, which is
                exactly the world POP describes. **World SIZE is separated by a second null that
                permutes labels but keeps the SIZES DELIBERATELY UNEQUAL in the same pattern** —
                if the observed T sits inside THAT too, size is sufficient and content is not needed.
Prediction      PERSON -> T > 0 beyond both nulls, same sign on >=2 of 3 stratifiers
Matrix          POP    -> T inside the label-permutation null
                SIZE   -> T inside the unequal-size null but outside a size-equalised one
Confound        ⚠ written BEFORE the run: `attend`, `polviews` and `educ` are correlated with one
                another, so three "independent" stratifiers are not three independent tests. The
                grid is reported whole and the agreement count is NOT treated as 3 draws.
Controls        positive: a graded person-level plant on a SYNTHETIC covariate, stratified the same
                way, must recover monotonically and must NOT fire at g=0 · negative: the
                label-permutation null · placebo: random strata of the same sizes
Stopping Rule   T beyond both nulls on >=2 stratifiers -> PERSON. Inside -> POP, and `#893`'s
                gradient paragraph is withdrawn from both pages. Inside the unequal-size null only
                -> SIZE. Budget: one round; `#111c` — A105 has ONE unverified already (`#893`).
Cost            n = 1,960, CPU minutes.
Priority        The last person-level question this instrument can answer in this arc.
Expected        If POP: A105 ends with a real finding about ITEMS and none about PEOPLE, and the
Transform       arc should close rather than spawn a fourth round on the same coupling.
```

⚠ **BASIN CHECK, run before the design.** `R328 → R330 → R331 → R333` is four rounds on one
coupling. Two of them (`#890`, `#893`) ended in a retraction of my own previous entry, so this is
not a confirmation basin by the `frontier` §3 test — **but the unwelcome branch is named above and
is the one that closes the arc.** `#111c` is live: **A105 already carries one `UNVERIFIED` (`#893`).
A second on this same question ends the line rather than buying a third.**

⚠ **PRIOR ART.** The two-class structure of the GSS abortion battery is textbook (Ebaugh & Haney
1980, `D6`). **That a person's religiosity predicts abortion attitude LEVEL is likewise textbook and
is NOT what is asked here** — this asks whether it predicts the **STRUCTURE** of the coupling
between a moral verdict and a legal one, which the level literature does not address. Stated before
the run so that a positive result is not quoted as the level finding in new clothes.

`G1` **ESTIMAND**: `T = B*(stratum 3) − B*(stratum 1)` where
`B* = [mean log|ρ|(abdefctw, E′) − mean log|ρ|(abdefctw, T′)] − baseline(stratum)`,
`E′ = {abnomore, absingle, abany}`, `T′ = {abhlth, abrape}`, and `baseline` is the **one-factor
world re-fitted inside that stratum** (`#893`③: a bracket is not the interaction, so the baseline
is re-derived for the PART and never inherited from the WHOLE). **Population** GSS respondents with
both `ab*w` norms and all seven legality items, waves 1991/1998/2008. **Instrument** GSS
`gss7224_r3a`; ⚠ the two norms are `D6` (Stata label only — `#891` measured that no shipped PDF
carries their question text), the legality items `D8`. **Baseline** the label-permutation null.
**Regime** three strata per stratifier, n ≈ 520–900 each.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES.** Under POP the tilt is the same in every stratum, so a
contrast between two strata is exactly 0. ⇒ **`negative_control`**, **kind of null named: a
STRATUM-LABEL PERMUTATION NULL with stratum sizes preserved** — it destroys the association between
who a person is and how their verdict couples, and keeps the marginals, the item set, the sample
size and the estimator identical.

**PRE-REGISTERED KILL — a conditional:**
```
if positive_control fires (graded person-level plant on a synthetic covariate: monotone, floor and
                           ceiling measured, does NOT fire at g=0)
   and placebo (random strata of the same sizes) is null:
       |T| > label-permutation 95th percentile on >=2 of 3 stratifiers, same sign  -> PERSON
       |T| inside that null on >=2 of 3                                            -> POP
       outside the equal-size null but inside the unequal-size one                 -> SIZE
else:
       UNVERIFIED
```
`G3`/`G4`: {3 stratifiers} × {3 estimators} × {3 class partitions} = 27 trend cells, published whole.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **three ordered strata give a contrast, not a dose–response** — the shape between the endpoints is
   not identifiable, and a monotone-looking triple is reported as a contrast plus its middle point;
② **`attend`, `polviews`, `educ` are mutually correlated** — three stratifiers are **not** three
   independent tests and the agreement count is not a p-value;
③ **causally identified N/A** — nobody is randomised into being religious;
④ **cross-instrument N/A — there is `no second instrument` for this estimand and `only this one
   instrument` can be asked it**: `#891` measured that of the 8 releases in `data/external/` exactly
   one ships question text;
⑤ **wave 2018 is out** (`abdefctw` ends 2008);
⑥ **the undifferentiated 47.1% contribute nothing by construction** (zero sanction variance), so
   every stratum's estimate rests on its differentiators only, and a stratifier that mostly sorts
   people INTO that group loses power for a reason that is not about coupling;
⑦ no second coder, no second release, no test–retest.
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
RNG = np.random.default_rng(333)
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
P0 = PARTS["textbook (Ebaugh & Haney 1980)"]
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

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
    raise SystemExit("STOP: empty/too-small population must never pass")

NRM = {c: d[c].astype(float) for c in NORMW}
SAN = {c: (d[c] == 2).astype(float).where(d[c].notna()) for c in AB}


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
    sn = src[0] if src else NRM
    ss = src[1] if src else SAN
    E = [s for s in part[0] if s not in TWINS]
    T = [s for s in part[1] if s not in TWINS]
    v = {}
    for s in E + T:
        k = sn[norm].notna() & ss[s].notna() & mask
        if k.sum() < 60:
            return np.nan
        r = abs(_rho(sn[norm][k].to_numpy(), ss[s][k].to_numpy(), est))
        if not np.isfinite(r) or r <= 1e-6:
            return np.nan
        v[s] = np.log(r)
    return float(np.mean([v[s] for s in E]) - np.mean([v[s] for s in T]))


def fit(mask):
    ix = np.flatnonzero(mask.to_numpy())
    vv = np.column_stack([NRM[c].to_numpy()[ix] for c in NORMW] + [SAN[c].to_numpy()[ix] for c in AB])
    vv = vv[~np.isnan(vv).any(axis=1)]
    zz = (vv - vv.mean(0)) / np.where(vv.std(0) > 0, vv.std(0), 1)
    gl = zz.mean(1)
    return ({c: float(np.corrcoef(zz[:, i], gl)[0, 1]) if zz[:, i].std() > 0 else 0.0
             for i, c in enumerate(NORMW + AB)}, len(vv),
            {c: float(np.nanmean(SAN[c].to_numpy()[ix])) for c in AB})


def baseline(mask, part, nsim=150):
    """`#893`③: the baseline is RE-DERIVED for the part, never inherited from the whole.

    Norm-independent under one factor, so both norms' brackets are pooled into one estimate."""
    lm, nn, mg = fit(mask)
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
        for nm in NORMW:
            b = bracket(nm, mk, part, src=(sn, ss))
            if np.isfinite(b):
                vals.append(b)
    return float(np.median(vals)), float(np.std(vals, ddof=1))


def bstar(mask, part=P0, est="spearman", base=None):
    b = bracket("abdefctw", mask, part, est)
    if not np.isfinite(b):
        return np.nan
    if base is None:
        base = baseline(mask, part)[0]
    return b - base


STRAT = {
    "attend": [("low 0–2", d["attend"] <= 2), ("mid 3–5", (d["attend"] >= 3) & (d["attend"] <= 5)),
               ("high 6–8", d["attend"] >= 6)],
    "polviews": [("liberal 1–3", d["polviews"] <= 3), ("moderate 4", d["polviews"] == 4),
                 ("conservative 5–7", d["polviews"] >= 5)],
    "educ": [("≤12", d["educ"] <= 12), ("13–15", (d["educ"] >= 13) & (d["educ"] <= 15)),
             ("≥16", d["educ"] >= 16)],
}

# ⚠ SPEED IS A DESIGN CONSTRAINT HERE, and getting it wrong once already cost a 10-minute kill:
#   v1 let `bstar()` re-simulate its own baseline inside every null draw — 160 draws x 3 strata x
#   150 sims per stratifier. **Every baseline is now computed ONCE and passed in.** The baselines
#   differ enormously between real strata (`#893` measured attend low 0.7359 vs high 0.2171), so
#   they can NOT be cancelled or dropped; they can only be cached.
print("\n=== (1) BASELINES, computed once per (stratum x partition) — they do not cancel ===")
BASE_OBS = {}
for tag, groups in STRAT.items():
    for i, (gname, g) in enumerate(groups):
        for pname, part in PARTS.items():
            BASE_OBS[(tag, i, pname)] = baseline(M & g, part, 120)[0]
    print(f"  {tag:9s} " + " · ".join(
        f"{groups[i][0]}={BASE_OBS[(tag,i,'textbook (Ebaugh & Haney 1980)')]:+.4f}" for i in range(3)))

print("\n=== (2) THE OBSERVED CONTRAST T = B*(last) - B*(first), per stratifier ===")
OBS, SIZES = {}, {}
for tag, groups in STRAT.items():
    vals, ns = [], []
    for i, (gname, g) in enumerate(groups):
        mk = M & g
        ns.append(int(mk.sum()))
        vals.append(bstar(mk, P0, "spearman", base=BASE_OBS[(tag, i, "textbook (Ebaugh & Haney 1980)")]))
    OBS[tag], SIZES[tag] = vals, ns
    print(f"  {tag:9s} n={ns}  B* " + " -> ".join(f"{v:+.4f}" for v in vals)
          + f"   **T = {vals[-1]-vals[0]:+.4f}**")

print("\n=== (3) NEGATIVE CONTROL — stratum-label permutation, SIZES PRESERVED ===")
idx = np.flatnonzero(M.to_numpy())


def rand_mask(size, perm, off):
    mk = pd.Series(False, index=d.index)
    mk.iloc[perm[off:off + size]] = True
    return mk


BASE_SZ = {}
for tag in STRAT:
    for sz in SIZES[tag]:
        key = round(sz, -2)
        if key not in BASE_SZ:
            BASE_SZ[key] = baseline(rand_mask(sz, RNG.permutation(idx), 0), P0, 120)[0]
print("  baselines for a RANDOM subsample of each size (the null's own baseline, cached): "
      + " · ".join(f"n~{k}:{v:+.4f}" for k, v in sorted(BASE_SZ.items())))
NULLS = {}
for tag in STRAT:
    ns = SIZES[tag]
    ts = []
    for _ in range(300):
        perm = RNG.permutation(idx)
        off, vals = 0, []
        for sz in ns:
            vals.append(bstar(rand_mask(sz, perm, off), P0, "spearman", base=BASE_SZ[round(sz, -2)]))
            off += sz
        if all(np.isfinite(v) for v in vals):
            ts.append(vals[-1] - vals[0])
    NULLS[tag] = np.asarray(ts)
    a = NULLS[tag]
    print(f"  {tag:9s} null median {np.median(a):+.4f} · sd {a.std(ddof=1):.4f} · "
          f"**95th pct of |T| {np.percentile(np.abs(a), 95):.4f}**  ({len(a)} draws)")

print("\n=== (4) THE SIZE NULL — same permutation, EQUAL stratum sizes (separates world SIZE) ===")
EQ = {}
for tag in STRAT:
    n_eq = int(np.mean(SIZES[tag]))
    bq = BASE_SZ.setdefault(round(n_eq, -2), baseline(rand_mask(n_eq, RNG.permutation(idx), 0), P0, 120)[0])
    ts = []
    for _ in range(200):
        perm = RNG.permutation(idx)
        vals = [bstar(rand_mask(n_eq, perm, k * n_eq), P0, "spearman", base=bq) for k in range(3)]
        if all(np.isfinite(v) for v in vals):
            ts.append(vals[-1] - vals[0])
    EQ[tag] = np.asarray(ts)
    print(f"  {tag:9s} equal-size null 95th pct of |T| {np.percentile(np.abs(EQ[tag]), 95):.4f} "
          f"(n_eq={n_eq}) vs unequal {np.percentile(np.abs(NULLS[tag]), 95):.4f}")

print("\n=== (5) POSITIVE CONTROL — a graded PERSON-LEVEL plant on a synthetic covariate ===")
# ⚠ v3, and v2's diagnosis was WRONG — worth recording because it is the more useful error.
#   v1 used two separately-simulated baselines and I blamed simulation noise; v2 pooled them and
#   g=0 STILL returned -0.1031, with sd EXACTLY 0.0000. That zero variance is the tell: at g=0 no
#   plant is applied, so the baselines cancel and T is simply the raw bracket difference between
#   two random thirds — **a single fixed draw of the covariate, not a null.** The floor was one
#   split's sampling noise, frozen, because `cov` was drawn ONCE outside the loop.
#   ⇒ the covariate is redrawn INSIDE every replicate, so g=0 measures the expected contrast under
#   no effect and carries its own spread. *A control whose spread is exactly zero is not precise,
#   it is constant — and a constant cannot be a null.*
sweep, sw_sd = [], []
for g in (0.0, 0.10, 0.25, 0.50, 0.75, 1.0):
    ts = []
    for _ in range(40):
        cov = RNG.random(len(d))
        G1 = pd.Series(cov < 1 / 3, index=d.index)
        G3 = pd.Series(cov >= 2 / 3, index=d.index)
        ss = {c: SAN[c].copy() for c in AB}
        for sx in [x for x in P0[1] if x not in TWINS]:
            hit = pd.Series((RNG.random(len(d)) < g), index=d.index) & G1 & M
            ss[sx] = ss[sx].where(~hit, (NRM["abdefctw"] <= 1).astype(float))
        v1 = bracket("abdefctw", M & G1, P0, src=(NRM, ss))
        v3 = bracket("abdefctw", M & G3, P0, src=(NRM, ss))
        if np.isfinite(v1) and np.isfinite(v3):
            ts.append(v3 - v1)          # the baselines cancel between two random thirds
    sweep.append((g, float(np.mean(ts))))
    sw_sd.append(float(np.std(ts, ddof=1)))
for (g, v), sd_ in zip(sweep, sw_sd):
    print(f"  g={g:<5.2f} T {v:+.4f} +/- {sd_:.4f}   (covariate REDRAWN each replicate; g=0 -> 0)")
N95 = float(np.mean([np.percentile(np.abs(NULLS[t]), 95) for t in STRAT]))
PC_OK = abs(sweep[0][1]) < N95 < abs(sweep[-1][1])
print(f"  floor(g=0) {sweep[0][1]:+.4f} · ceiling(g=1) {sweep[-1][1]:+.4f} · mean threshold "
      f"{N95:.4f} strictly between: {PC_OK}")

print("\n=== (6) G3/G4 — 3 stratifiers x 3 estimators x 3 partitions = 27 trend cells ===")
rows, agree = [], {}
for tag, groups in STRAT.items():
    for pname, part in PARTS.items():
        for est in ("spearman", "kendall", "gamma"):
            vals = [bstar(M & g, part, est, base=BASE_OBS[(tag, i, pname)])
                    for i, (_, g) in enumerate(groups)]
            T = (vals[-1] - vals[0]) if all(np.isfinite(v) for v in vals) else None
            rows.append((tag, pname, est, None if T is None else float(T)))
            if T is not None:
                agree.setdefault(tag, []).append(T)
            print(f"  {tag:9s} {pname[:30]:30s} {est:9s} T = "
                  + ("n/a" if T is None else f"{T:+.4f}"))
pos = sum(1 for _, _, _, t in rows if t is not None and t > 0)
tot = sum(1 for _, _, _, t in rows if t is not None)
print(f"\n  **grid: {pos}/{tot} trend cells positive**")
fired = {}
for tag in STRAT:
    n95 = float(np.percentile(np.abs(NULLS[tag]), 95))
    v = agree.get(tag, [])
    hits = sum(1 for x in v if abs(x) > n95)
    fired[tag] = (float(np.median(v)) if v else np.nan, n95, hits, len(v))
    print(f"  {tag:9s} median T {fired[tag][0]:+.4f} vs its own null 95th {n95:.4f} => "
          f"{hits}/{len(v)} cells beyond")

print("\n=== (7) THE CONDITIONAL KILL ===")
G = Gate("Does whether a moral verdict is a CASE-JUDGEMENT depend on who holds it?")
G.plant_direction_from_sweep("positive: graded person-level plant on a synthetic covariate", sweep,
                             baseline=0.0, baseline_spread=sw_sd[0] if sw_sd[0] > 0 else None,
                             half_of=max(N95, 1e-4))
main_T = float(np.median([OBS[t][-1] - OBS[t][0] for t in STRAT]))
_allnull = np.concatenate([NULLS[t] for t in STRAT])
G.negative_control("stratum-label permutation, sizes preserved", float(np.median(_allnull)), main_T,
    null_spread=float(_allnull.std(ddof=1)),
    null_kind="STRATUM-LABEL PERMUTATION NULL with stratum sizes preserved - destroys the "
              "person-coupling association and keeps marginals, items, n and estimator identical")
G.resolvable("T (median over stratifiers)", main_T, N95 / 2)
maj = sum(1 for t in STRAT if fired[t][2] > fired[t][3] / 2 and np.sign(fired[t][0]) > 0)
inside = sum(1 for t in STRAT if fired[t][2] <= fired[t][3] / 2)
size_only = sum(1 for t in STRAT
                if abs(fired[t][0]) > np.percentile(np.abs(EQ[t]), 95) and
                abs(fired[t][0]) <= np.percentile(np.abs(NULLS[t]), 95))
if not PC_OK:
    VERDICT, WORLD = "UNVERIFIED", "the positive control did not license a reading"
elif maj >= 2:
    VERDICT, WORLD = "CONFIRMED", ("PERSON - whether a moral verdict is a case-judgement or a "
                                   "position varies with who holds it")
elif inside >= 2:
    VERDICT, WORLD = "OVERTURNED", ("POP - the tilt is a population constant; `#893`'s gradient "
                                    "paragraph is WITHDRAWN from both pages")
elif size_only >= 2:
    VERDICT, WORLD = "OVERTURNED", "SIZE - the contrast tracks stratum n, not stratum content"
else:
    VERDICT, WORLD = "UNVERIFIED", "mixed - no world predicted this pattern; report it, not a verdict"
print(G)
print(f"\n  stratifiers with a majority of cells beyond their own null and T>0: {maj}/3 · "
      f"inside: {inside}/3 · size-only: {size_only}/3")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ `attend`, `polviews` and `educ` are MUTUALLY CORRELATED - three stratifiers are not")
print("     three independent tests, and the agreement count above is NOT a p-value. Registered.")

art = dict(entry=895, round="E03·A105·R333", verdict=VERDICT, world=WORLD, n=int(M.sum()), waves=YEARS,
           observed={t: dict(bstar=OBS[t], sizes=SIZES[t], T=OBS[t][-1] - OBS[t][0]) for t in STRAT},
           null={t: dict(median=float(np.median(NULLS[t])), sd=float(NULLS[t].std(ddof=1)),
                         p95=float(np.percentile(np.abs(NULLS[t]), 95)), draws=len(NULLS[t]))
                 for t in STRAT},
           equal_size_null={t: float(np.percentile(np.abs(EQ[t]), 95)) for t in STRAT},
           positive_sweep=sweep, positive_sd=sw_sd, positive_ok=PC_OK,
           grid_rows=rows, grid_positive=pos, grid_total=tot,
           fired={t: dict(median=fired[t][0], null95=fired[t][1], beyond=fired[t][2],
                          cells=fired[t][3]) for t in STRAT},
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "person_or_population.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'person_or_population.json'}")
