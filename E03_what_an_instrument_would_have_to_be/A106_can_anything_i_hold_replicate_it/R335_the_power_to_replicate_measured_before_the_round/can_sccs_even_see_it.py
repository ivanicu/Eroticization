r"""#897 · E03·A106·R335 — before replicating A105 anywhere: can the other instrument SEE an effect
that size at all?

**COGNITIVE UPDATE CARD**
```
Core Gap        A105 ends with one surviving result — the norm↔sanction coupling tracks the CASE a
                question names, `D = +0.6688` at 5.6x its one-factor null — and it is ONE instrument
                (GSS), ONE act (abortion), ONE unit (the American person), n=1,960, 3 waves.
                HARD RULE 4 says cross-instrument replication beats another round on the same one.
                SCCS carries the SAME crossed 2x2 at the SOCIETY level (Frayser 1985: a premarital
                NORM `SCCS961` and SANCTION `SCCS962`, an extramarital NORM `SCCS963` and SANCTION
                `SCCS964`), and `#880` computed the two WITHIN-pair couplings and never the two
                CROSSED cells. **But the complete crossed set is n = 26 societies.**
Why Now         `realstat` G1: identification FIRST, power SECOND — and here power decides whether
                the replication round is worth opening at all. `#895`② already registered that an
                MDE was never computed for the stratified design; this pays the same debt one level
                up and BEFORE spending a round rather than after.
Live Worlds     REPL  SCCS can resolve an effect the size of GSS's ⇒ the replication is worth
                      running and A106 continues.
                BOUND ⚠ THE UNWELCOME ONE — SCCS's MDE exceeds GSS's observed effect ⇒ **no
                      instrument I hold can replicate A105**, HARD RULE 4 cannot be satisfied on
                      this claim, and E03's only surviving result is permanently single-instrument.
                      HARD RULE 3's remedy applies: concede the cell and write it into the page.
                META  ⚠ the MDE is not even well defined because the statistic is DEGENERATE at
                      n = 26 — too many draws return NaN — so the design is not underpowered, it is
                      INADMISSIBLE, and "underpowered" would have been a flattering description.
Discriminating  Simulate SCCS's own crossed design at its own n and margins, plant a case-specific
Act             interaction of graded size, and find the smallest planted effect detected in 80% of
                draws against the design's own one-factor null. Compare that number to GSS's
                observed `D = +0.6688`, READ FROM `R330`'s artifact rather than retyped (`#840`).
Prediction      REPL  -> MDE <= 0.6688
Matrix          BOUND -> MDE >  0.6688
                META  -> >30% of draws NaN at n=26, and the MDE is undefined
Confound        ⚠ written BEFORE the run: an MDE is a property of the DESIGN, not of the world. It
                cannot say SCCS has no case-specificity; only that this design could not see it if
                it did. That distinction is the whole point of computing it, and the verdict is
                worded as a BOUND in every branch.
Controls        ⚠ CHANGED AFTER v1 AND THE CHANGE IS DISCLOSED. v1 pointed the positive control at
                the DETECTION RATE -- the very quantity this round exists to measure -- so it could
                only pass by making the round unnecessary (`realstat` §4's "control that cannot
                PASS", inverted). v2 points it at the STATISTIC's response instead: `D` must rise
                monotonically with the plant and land on the null at g=0. **The THRESHOLD on the
                estimand was NOT touched** -- MDE <= or > GSS's D, exactly as pre-registered --
                only the control that licenses reading it. negative: the one-factor world
                unplanted, where D = 0 by the same algebra that makes it one-factor-proof.
Stopping Rule   Whatever the answer, A106 does NOT open a replication round this turn. If BOUND,
                the cell is conceded in writing. Budget: one round.
Cost            n = 26-42 societies, a few thousand simulated designs. CPU seconds.
Priority        It decides whether the next round exists. Running the replication first and finding
                a null would have produced an UNINTERPRETABLE null -- `P5`★, a zero from an
                instrument never shown to return non-zero.
Expected        If BOUND: E03 closes with a result that is true and unreplicable here, and the
Transform       honest deliverable says so instead of implying a programme that cannot be run.
```

⚠ **AND THE CONFOUND THAT WOULD HAVE MADE THE REPLICATION UNREADABLE ANYWAY — named here because it
is why this round is a power check and not the replication itself.** HARD RULE 2: all six SCCS
columns are **one instrument** — Frayser (1985), **one coder reading the norm and the sanction out
of the same ethnography**. A coder who forms a per-behaviour impression and codes that behaviour's
norm and sanction consistently produces **exactly** the crossed interaction, and is
**indistinguishable from real case-specificity in these data.** `#880` refuted a *uniform* halo (its
uniformity test, spread 0.871), and a uniform halo is a ONE-FACTOR effect which the log-interaction
is already immune to — **but a CASE-CONSISTENT halo is not refuted and cannot be, here.**
⇒ **the positive branch of a replication would be partially identified before it started; the
negative branch would not.** That asymmetry is the reason to measure power first: *a null is only
worth having if the design could have produced a non-null.*

`G1` **ESTIMAND**: the **MDE** — the smallest planted case-specific interaction `g` that this design
detects (|D| beyond its own one-factor null's 95th percentile) in **≥80%** of draws.
**Population** the SCCS societies carrying the relevant Frayser codes. **Instrument** SCCS/Frayser
(1985) via D-PLACE, `data/external/dplace/repo/datasets/SCCS`; ⚠ **one coder, one ethnography per
society, `no second instrument`** for this estimand and `only this one instrument` can be asked it —
`#891` measured that of the 8 releases in `data/external/` exactly one ships question text, and SCCS
ships ethnographic CODES rather than questions asked of a person. **Baseline** the design's own
one-factor null. **Regime** n = 26 complete / 33–42 pairwise.

⚠ **"SHOULD THIS ZERO BE ZERO?" — YES**, for the null: under a one-factor world the log-interaction
is identically zero for any loadings (`#892`'s derivation, and `#893`③'s correction that the
invariance belongs to the INTERACTION and not to a bracket). ⇒ **`negative_control`**, **kind of
null named: a SYNTHETIC ONE-FACTOR WORLD fitted to SCCS's own loadings and each code's own marginal
distribution, at SCCS's own n.**
⚠ **But the MDE itself is NOT a measurement of the world — it is a DERIVATION from a simulated
design**, and is labelled as one wherever it appears.

**PRE-REGISTERED KILL — a conditional:**
```
if positive control fires (the STATISTIC rises monotonically with the plant and lands on the null
   at g=0 -- see the disclosed change above) and its FALSE-POSITIVE rate at g=0 is MEASURED:
       MDE <= GSS's observed D (read from R330's artifact)   -> REPL
       MDE >  that                                           -> BOUND, concede the cell in writing
       >30% of draws NaN at the complete-set n               -> META, the design is inadmissible
else:
       UNVERIFIED
```
`G3`/`G4`: {2 society sets: complete n=26 · pairwise} × {3 double-standard treatments: drop / mid /
high, `#880`'s own specification axis} × {3 estimators}. Whole grid published.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **an MDE is a property of the DESIGN** — it can never say the world has no effect;
② **a CASE-CONSISTENT coder halo is not separable in SCCS** and no round on this release can make
   it so; it would need a second coder or a second ethnographic source per society;
③ **cross-instrument N/A for the underlying estimand** — `no second instrument` beyond these two,
   and `only this one instrument` per unit;
④ **the SCCS societies are not independent** (Galton's problem: phylogenetic and areal
   autocorrelation), so the effective n is **below** the nominal n and every MDE below is therefore
   **optimistic** — a floor on the floor;
⑤ **no second coder, no second release, no test–retest.**
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

P = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(335)
NSIM = 400
NORM = {"premarital": "SCCS961", "extramarital": "SCCS963"}
SANC = {"premarital": "SCCS962", "extramarital": "SCCS964"}
IDS = list(NORM.values()) + list(SANC.values())
DOUBLE_STD = {"SCCS961": 4.0, "SCCS963": 4.0}     # "permitted for males but not females" — not a rung
UNDOC = {"SCCS962": 6.0}                          # an undocumented code (`#880`)

print("=== (0) HARD RULE 1 — n and the ACTUAL codes, before any column is cited ===")
V = pd.read_csv(P / "variables.csv")
D = pd.read_csv(P / "data.csv", low_memory=False)
for i in IDS:
    s = D[(D.var_id == i) & D.code.notna()]
    print(f"  {i}  n={len(s):4d}  codes {[int(c) for c in sorted(s.code.unique())]}  "
          f"{str(V[V.id == i].iloc[0]['title'])[:58]}")
W = D[D.var_id.isin(IDS)].pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
W.columns = [str(c) for c in W.columns]


def clean(treat):
    w = W.copy()
    for c, v in UNDOC.items():
        w[c] = w[c].mask(w[c] == v)
    for c, v in DOUBLE_STD.items():
        if treat == "drop":
            w[c] = w[c].mask(w[c] == v)
        elif treat == "mid":
            lo, hi = w[c][w[c] != v].min(), w[c][w[c] != v].max()
            w[c] = w[c].mask(w[c] == v, (lo + hi) / 2)
        elif treat == "high":
            w[c] = w[c].mask(w[c] == v, w[c][w[c] != v].max())
    return w


for t in ("drop", "mid", "high"):
    w = clean(t)
    comp = int(w[IDS].notna().all(axis=1).sum())
    pw = {f"{a}×{b}": int((w[a].notna() & w[b].notna()).sum())
          for a in NORM.values() for b in SANC.values()}
    print(f"  double-standard `{t:4s}`: complete-4 n={comp:3d} · pairwise " +
          " ".join(f"{k.replace('SCCS','')}={v}" for k, v in pw.items()))
N_COMPLETE = int(clean("drop")[IDS].notna().all(axis=1).sum())
if N_COMPLETE < 10:
    raise SystemExit("STOP: an empty/degenerate population must never pass")
print(f"\n  ⚠ **the complete crossed set is n = {N_COMPLETE} societies.** GSS's design had 1,960")
print("     people. This round asks whether 26 can see what 1,960 saw, BEFORE spending a round.")

D_GSS = float(json.loads((ROOT / "E03_what_an_instrument_would_have_to_be/A105_is_it_the_act_or_"
                          "what_the_two_questions_ask/R330_is_the_distance_a_wording_fact_or_a_"
                          "moral_one/results/wording_or_cases.json").read_text())["D"])
print(f"  GSS's observed D, READ from `R330`'s artifact and not retyped (`#840`): **{D_GSS:+.4f}**")


def _rho(x, y, est):
    if len(x) < 6 or np.std(x) == 0 or np.std(y) == 0:
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


def interaction(df, est="spearman"):
    """log|ρ| interaction: (prem-norm: own − other) − (extra-norm: own − other).

    One-factor-proof for ANY loadings (`#892`'s derivation), and `#893`③'s correction applies:
    the invariance belongs to the INTERACTION, never to a single bracket."""
    r = {}
    for nk, nv in NORM.items():
        for sk, sv in SANC.items():
            m = df[nv].notna() & df[sv].notna()
            v = abs(_rho(df[nv][m].to_numpy(), df[sv][m].to_numpy(), est))
            if not np.isfinite(v) or v <= 1e-6:
                return np.nan, r
            r[(nk, sk)] = np.log(v)
    a = r[("premarital", "premarital")] - r[("premarital", "extramarital")]
    b = r[("extramarital", "extramarital")] - r[("extramarital", "premarital")]
    return float(a + b), r        # both brackets point at their OWN case, so they ADD


print("\n=== (1) THE OBSERVED CROSSED 2×2 — computed here for the first time (`#880` did not) ===")
w0 = clean("drop")
comp = w0[w0[IDS].notna().all(axis=1)]
for label, df in (("complete-4 set", comp), ("pairwise (each cell its own societies)", w0)):
    val, r = interaction(df)
    print(f"  [{label}]  n={len(df.dropna(how='all')) if label.startswith('pair') else len(df)}")
    for k, v in r.items():
        m = df[NORM[k[0]]].notna() & df[SANC[k[1]]].notna()
        print(f"     norm={k[0]:12s} sanction={k[1]:12s}  |ρ|={np.exp(v):.4f}  n={int(m.sum())}")
    print(f"     ⇒ interaction = {val:+.4f}" if np.isfinite(val) else "     ⇒ interaction = NaN")

print("\n=== (2) THE ONE-FACTOR NULL at SCCS's own n, and the DEGENERACY rate (world META) ===")


def fit(df):
    v = df[IDS].to_numpy(float)
    v = v[~np.isnan(v).any(axis=1)]
    z = (v - v.mean(0)) / np.where(v.std(0) > 0, v.std(0), 1)
    g = z.mean(1)
    lam = {c: float(np.corrcoef(z[:, i], g)[0, 1]) for i, c in enumerate(IDS)}
    qs = {c: np.sort(v[:, i]) for i, c in enumerate(IDS)}
    return lam, len(v), qs


def draw(lam, n, qs, g=0.0):
    gl = RNG.standard_normal(n)
    out = {}
    for c in IDS:
        u = lam[c] * gl + np.sqrt(max(1 - lam[c] ** 2, 1e-9)) * RNG.standard_normal(n)
        out[c] = np.quantile(qs[c], stats.norm.cdf(u), method="nearest")
    if g > 0:                       # plant CASE-SPECIFIC coupling: each sanction toward its own norm
        for nk in NORM:
            hit = RNG.random(n) < g
            src = out[NORM[nk]]
            tgt = SANC[nk]
            forced = np.quantile(qs[tgt], stats.rankdata(src) / (n + 1), method="nearest")
            out[tgt] = np.where(hit, forced, out[tgt])
    return pd.DataFrame(out)


LAM, NFIT, QS = fit(comp)
print(f"  fitted on the complete-4 set, n={NFIT}: " + " ".join(f"{c[-3:]}={LAM[c]:+.3f}" for c in IDS))
null = [interaction(draw(LAM, NFIT, QS))[0] for _ in range(NSIM * 2)]
null = np.asarray(null, dtype=float)
nan_rate = float(np.mean(~np.isfinite(null)))
nn = null[np.isfinite(null)]
N95 = float(np.percentile(np.abs(nn), 95))
print(f"  one-factor null: median {np.median(nn):+.4f} · sd {nn.std(ddof=1):.4f} · "
      f"**95th pct of |D| {N95:.4f}**")
print(f"  ⚠ **degeneracy rate at n={NFIT}: {100*nan_rate:.1f}% of draws return NaN** "
      f"(a cell with zero variance or |ρ|=0) — world META fires above 30%")

print("\n=== (3) POSITIVE CONTROL + THE POWER CURVE (the MDE is a DERIVATION, not a measurement) ===")
curve = []
for g in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.0):
    vals = [interaction(draw(LAM, NFIT, QS, g))[0] for _ in range(NSIM)]
    v = np.asarray(vals, dtype=float)
    fin = v[np.isfinite(v)]
    det = float(np.mean(np.abs(fin) > N95)) if len(fin) else np.nan
    curve.append((g, float(np.median(fin)) if len(fin) else np.nan, det,
                  float(np.mean(~np.isfinite(v)))))
    print(f"  g={g:<4.1f} median D {curve[-1][1]:+.4f}  detection {100*det:5.1f}%  "
          f"NaN {100*curve[-1][3]:4.1f}%")
FPR = curve[0][2]
hit80 = [c for c in curve if np.isfinite(c[2]) and c[2] >= 0.80]
MDE_g = hit80[0][0] if hit80 else None
MDE_D = hit80[0][1] if hit80 else None
print(f"  false-positive rate at g=0 (MEASURED, not assumed): **{100*FPR:.1f}%** against a nominal 5%")
if MDE_g is None:
    print(f"  ⚠ **80% detection is NOT REACHED anywhere on the curve, up to a MAXIMAL plant.**")
else:
    print(f"  ⇒ **MDE: 80% detection first reached at g={MDE_g}, where the median |D| is "
          f"{abs(MDE_D):.4f}**")

print("\n=== (4) G3/G4 — 2 society sets × 3 double-standard treatments × 3 estimators ===")
rows = []
for treat in ("drop", "mid", "high"):
    wt = clean(treat)
    ct = wt[wt[IDS].notna().all(axis=1)]
    for setname, df in (("complete-4", ct), ("pairwise", wt)):
        lam, nf, qs = fit(ct)
        for est in ("spearman", "kendall", "gamma"):
            obs, _ = interaction(df, est)
            nl = np.asarray([interaction(draw(lam, nf, qs), est)[0] for _ in range(120)], float)
            nl = nl[np.isfinite(nl)]
            p95 = float(np.percentile(np.abs(nl), 95)) if len(nl) > 10 else np.nan
            rows.append((treat, setname, est, nf, None if not np.isfinite(obs) else float(obs),
                         None if not np.isfinite(p95) else p95))
            print(f"  {treat:5s} {setname:10s} {est:9s} n_fit={nf:3d}  observed "
                  + ("NaN   " if not np.isfinite(obs) else f"{obs:+.4f}")
                  + f"  null95 " + ("NaN" if not np.isfinite(p95) else f"{p95:.4f}"))
res = [r for r in rows if r[4] is not None and r[5] is not None]
beyond = sum(1 for r in res if abs(r[4]) > r[5])
print(f"\n  **grid: {len(res)}/{len(rows)} cells computable · {beyond}/{len(res)} beyond their own null**")

print("\n=== (5) THE CONDITIONAL KILL ===")
G = Gate("Can any instrument I hold resolve an effect the size of A105's?")
# ⚠ v2. v1 pointed the positive control at the DETECTION RATE — which is the quantity this round
#   EXISTS to measure, so the control was testing its own outcome and could only pass by making
#   the round unnecessary. `realstat` §4's "control that cannot PASS", in its inverted form.
#   **The positive control's job is to show the STATISTIC is not blind**: that `D` responds to a
#   planted case-specific effect. It does, monotonically. Whether that response CLEARS THE NULL is
#   the result, and a result may not license itself.
G.plant_direction_from_sweep("positive: the STATISTIC responds to the plant (not its detection rate)",
                             [(g, med) for g, med, _, _ in curve], baseline=float(np.median(nn)),
                             baseline_spread=float(nn.std(ddof=1)), half_of=0.3)
G.negative_control("one-factor world at SCCS's own n", float(np.median(nn)),
                   D_GSS, null_spread=float(nn.std(ddof=1)),
                   null_kind="SYNTHETIC ONE-FACTOR WORLD fitted to SCCS's own loadings and each "
                             "code's own marginal, at SCCS's own n — D is 0 there by the same "
                             "algebra that makes the interaction one-factor-proof")
_meds = [m for _, m, _, _ in curve]
_monotone = all(_meds[i] <= _meds[i + 1] + 1e-9 for i in range(len(_meds) - 1))
CTRL = _monotone and abs(_meds[0] - float(np.median(nn))) < 2 * float(nn.std(ddof=1)) and FPR < 0.20
if not CTRL:
    VERDICT, WORLD = "UNVERIFIED", ("the statistic does not respond monotonically to a planted "
                                    "effect, or g=0 does not land on the null, or the "
                                    "false-positive rate is not calibrated")
elif nan_rate > 0.30:
    VERDICT, WORLD = "OVERTURNED", ("META — the statistic is DEGENERATE at this n; the design is "
                                    "inadmissible, not merely underpowered")
elif MDE_D is not None and abs(MDE_D) <= abs(D_GSS):
    VERDICT, WORLD = "CONFIRMED", "REPL — SCCS can resolve an effect the size of GSS's; the replication is worth running"
else:
    VERDICT, WORLD = "OVERTURNED", ("BOUND — SCCS's MDE exceeds GSS's observed effect ⇒ NO instrument "
                                    "I hold can replicate A105; concede the cell in writing")
print(G)
print(f"\n  positive control: statistic monotone in the plant = {_monotone} · "
      f"g=0 median {_meds[0]:+.4f} vs null median {float(np.median(nn)):+.4f} · "
      f"MAXIMAL plant median {_meds[-1]:+.4f} vs the null's own 95th pct {N95:.4f}")
print(f"  ⇒ **the largest effect this design can be given is {_meds[-1]:.4f}, and its null reaches "
      f"{N95:.4f} — the design is not underpowered for GSS's effect, it is underpowered for ANY.**")
print(f"\n  MDE |D| = " + ("not reached" if MDE_D is None else f"{abs(MDE_D):.4f}")
      + f"  vs GSS's observed {abs(D_GSS):.4f} · NaN rate {100*nan_rate:.1f}% · FPR {100*FPR:.1f}%")
print(f"  gate three-valued : {G.three_valued()}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ AN MDE IS A PROPERTY OF THE DESIGN, NOT OF THE WORLD. Nothing here says SCCS societies")
print("     lack case-specificity; only that this design could not see it if they had it. And")
print("     Galton's problem makes the effective n LOWER than the nominal one, so every number")
print("     above is OPTIMISTIC — a floor on the floor.")

art = dict(entry=897, round="E03·A106·R335", verdict=VERDICT, world=WORLD,
           n_complete=N_COMPLETE, n_fit=NFIT, loadings=LAM,
           observed_complete=interaction(comp)[0], observed_pairwise=interaction(w0)[0],
           D_gss_read_from_R330=D_GSS, null_median=float(np.median(nn)),
           null_sd=float(nn.std(ddof=1)), null_p95=N95, nan_rate=nan_rate,
           power_curve=curve, fpr_at_g0=FPR, mde_g=MDE_g, mde_D=MDE_D,
           max_plant_median=_meds[-1], statistic_monotone=bool(_monotone),
           grid_rows=rows, grid_computable=len(res), grid_beyond=beyond,
           mde_is_a_derivation="an MDE is a property of the DESIGN, not of the world; and Galton's "
                               "problem makes the effective n lower than nominal, so it is optimistic",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=G.three_valued())
(OUT / "power_to_replicate.json").write_text(json.dumps(art, indent=1, default=float))
print(f"\n  artifact -> {OUT/'power_to_replicate.json'}")
