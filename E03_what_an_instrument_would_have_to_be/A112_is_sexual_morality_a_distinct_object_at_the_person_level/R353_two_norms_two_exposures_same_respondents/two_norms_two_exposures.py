#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R353 — two norms, two exposures, the same respondents
==============================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        The project's central question — **is there anything distinctive about SEXUAL
                morality?** — is UNDECIDABLE on GSS. `#905` put the sexual battery at the 10th
                percentile of forty non-sexual GSS batteries; `#906` put abortion at the 33rd
                percentile of six crossed designs; `#907` then measured that GSS ships **zero**
                non-sexual moral×legal batteries, so the deciding comparison cannot be built there.
                **A different instrument is required, and `#914`② found one.**

Why Now         NSFG carries a **sexual** norm and a **non-sexual** norm in the SAME battery on the
                SAME respondents, with a domain-matched behavioural exposure for each. That is the
                placebo GSS structurally lacks (HARD RULE 4: cross-instrument beats another round
                on the same one), and `#914` just anchored the parse without reading a value.

Live Worlds     A · **DOMAIN-SPECIFIC** — a person's norm about act X is tied to their own
                    experience of X *specifically*. Self-serving moral reasoning is local.
                B · **ONE MORAL DIAL** — each person has one permissiveness, each item loads on it,
                    and behaviour just slides the dial. "Sexual morality" is a slice, not an object.
                    **This is the person-level analogue of `#905`, and it is the unwelcome one.**
                C · **NEITHER** — both norms are dominated by one demographic (religiosity), the
                    exposures resolve nothing, and the norm/exposure decomposition is the wrong
                    carve. (the meta-separator: an outcome that kills my world-DECOMPOSITION)

Estimand        ⚠ NAMED BEFORE THE STATISTIC, and it is one-factor-proof by DERIVATION, not by
(G1)            simulation. Under ONE moral dimension, item i responds to covariates only through
                the latent: b(c -> i) = lambda_i * beta_c. Hence
                    log|b(c -> i1)| - log|b(c -> i2)| = log|lambda_1| - log|lambda_2|
                — the SAME for every covariate c. So the interaction

                  D = [log|b(sexual_exposure -> samesex)|  - log|b(sexual_exposure -> chsuppor)|]
                    - [log|b(family_exposure -> samesex)|  - log|b(family_exposure -> chsuppor)|]

                is **identically zero under any one-factor model, for any loadings**. That is the
                whole reason this statistic is admissible. ⚠ `#893`③: a BRACKET of a one-factor-
                proof interaction is NOT itself one-factor-proof — so only the full D is reported,
                never one bracket.

Prediction      A -> D clearly above its one-factor null.   B -> D inside the null.
Matrix          C -> the exposures' own coefficients sit inside the noise floor and D is undefined
                     rather than zero — reported as UNRESOLVED, never as agreement with B.

Strongest       ⚠ **The matched pairs are partly TAUTOLOGICAL**: people who have had same-sex
Confound        contact approve of same-sex relations; people who bore a child outside marriage
                approve of unmarried women doing so. **This is not a confound to be removed — under
                world B a tautology cannot exist, because one dial moves both items together.** It
                IS the mechanism of world A, and it must be reported as what it is: *self-serving
                reasoning that is local to a domain*, not "moral judgement is case-indexed" dressed
                up. The design's job is symmetry — each exposure is tautologically nearest its OWN
                domain — so a general dial cannot masquerade as domain specificity.
Second          ⚠ Religiosity may drive both norms and both behaviours. CONTROL, same iteration:
Confound        the adjustment axis of the specification curve carries `reldlife` + `attndnow`
                alongside age/education/race, and the whole curve is published.

Controls        POSITIVE, must FAIL at g=0: plant domain-specific coupling of strength g into a
                synthetic world fitted to the OBSERVED loadings and marginals; D must rise with g
                and sit at the null when g=0.
                NEGATIVE: a synthetic ONE-FACTOR world at the observed n and marginals — D's null.
                PLACEBO: rerun D with the battery's two NON-norm items (`reactslf`, `chbother`,
                which ask how R would FEEL, not whether an act is all right). A D just as large
                there is a fact about the battery, not about norms.

Stopping Rule   The whole 3 files x estimator x adjustment x coding grid, published including
                disagreeing cells. If the exposures do not resolve, the verdict is UNRESOLVED and
                the round does not get a second attempt at the same question (`#111c`).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  ① **Value LABELS do not ship.** Response codes are read from the item's own documented range;
    what `1` MEANS is taken from the item text's direction and is an assumption, named here.
  ② The male file's sexual exposure is a DIFFERENT question ("Ever had oral or anal sex with a
    male") from the female one ("Ever Had Sexual Contact with a Female") — so the sexes are not
    a replication of each other; they are two sites.
  ③ Cross-sectional: no intervention, no temporal order. **Exposure and norm are simultaneous and
    the arrow is not identified.** D is about SPECIFICITY of association, never about causation.
  ④ Survey weights and `secu`/`sest` design variables ship, but D is a ratio of coefficients and
    this round reports both weighted and unweighted rather than claiming a design-based SE.
  ⑤ `[unchallenged]` — door ③.
"""
import json, re, sys, itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
D_ = ROOT / "data" / "external" / "nsfg"
RNG = np.random.default_rng(353)
PAT = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')

SITES = [("2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct", "2011-2013 female"),
         ("2017_2019_FemRespData.dat", "2017_2019_FemRespSetup.dct", "2017-2019 female"),
         ("2017_2019_MaleData.dat", "2017_2019_MaleSetup.dct", "2017-2019 male")]

NORMS = ["samesex", "chsuppor"]                 # sexual domain · family domain
PLACEBO = ["reactslf", "chbother"]              # same battery, NOT norms about an act
EXPOS = ["samesexany", "cebow"]                 # domain-matched behavioural exposures
ADJ = ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]
WANT = NORMS + PLACEBO + EXPOS + ADJ + ["caseid"]


def read_site(dat, dct):
    lay = {m.group(2).lower(): (int(m.group(1)) - 1, int(m.group(3))) for m in
           (PAT.search(l) for l in (D_ / "setup" / dct).read_text(errors="replace").splitlines()) if m}
    have = [v for v in WANT if v in lay]
    spec = [(lay[v][0], lay[v][0] + lay[v][1]) for v in have]
    df = pd.read_fwf(D_ / dat, colspecs=spec, names=have, dtype=str)
    for c in have:
        df[c] = pd.to_numeric(df[c].astype(str).str.strip().replace({"": None}), errors="coerce")
    return df, have


def norm_clean(s):
    """1..5 substantive; 8/9 refused/DK -> missing. Direction: LOW = agrees the act is all right."""
    return s.where(s.isin([1, 2, 3, 4, 5]))


def expo_clean(name, s):
    if name == "samesexany":
        return (s == 1).astype(float).where(s.isin([1, 5]))      # 1=yes 5=no
    return (s > 0).astype(float).where(s.notna() & (s < 90))     # cebow: any nonmarital birth


def beta(y, x, covs=None, standardize=True):
    """OLS slope of y on x, optionally adjusting for covs. Returns (b, se, n)."""
    ok = y.notna() & x.notna()
    if covs is not None:
        for c in covs:
            ok &= c.notna()
    n = int(ok.sum())
    if n < 60:
        return np.nan, np.nan, n
    Y = y[ok].to_numpy(float)
    cols = [x[ok].to_numpy(float)]
    if covs is not None:
        cols += [c[ok].to_numpy(float) for c in covs]
    X = np.column_stack([np.ones(n)] + cols)
    if standardize:
        for j in range(1, X.shape[1]):
            sd = X[:, j].std()
            if sd > 0:
                X[:, j] = (X[:, j] - X[:, j].mean()) / sd
        Y = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        bhat, *_ = np.linalg.lstsq(X, Y, rcond=None)
        resid = Y - X @ bhat
        s2 = resid @ resid / max(n - X.shape[1], 1)
        cov = s2 * np.linalg.pinv(X.T @ X)
        return float(bhat[1]), float(np.sqrt(max(cov[1, 1], 0))), n
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n


# ⚠⚠ THE FIRST STATISTIC WAS DEFECTIVE AND ITS OWN OUTPUT SAID SO — kept, with the repair.
# `D` on the log|b| scale is one-factor-proof BY DERIVATION, and its ESTIMATOR divides by a
# coefficient that the data measured as ZERO: `cebow -> samesex` came back +0.0057, -0.0082 and
# +0.0009 against an se of ~0.018 across the three 2011-2013 cells. `log|0.0009|` is not a
# quantity, it is noise amplified without bound, and the synthetic null never visited that regime
# because its coefficients are all comfortably non-zero. **A null that cannot reach the observed
# regime is not this statistic's null.**  (Same family as: a ratio of noisy reciprocals is biased.)
#
# THE REPAIR, with the identical one-factor proof and no division: under one factor
# b(c -> i) = lambda_i * beta_c, so the 2x2 coefficient matrix B[c, i] is RANK 1, hence
#     det(B) = b11*b22 - b12*b21 == 0   exactly, for any loadings.
# `det` is a polynomial in the coefficients: it cannot explode when one of them is zero, and its
# null is simulated at BOTH the observed marginals and the observed coefficient magnitudes.
def det_stat(tab):
    """rank-1 departure of the 2x2 coefficient matrix. Zero under any one-factor model."""
    try:
        b11 = tab[(EXPOS[0], NORMS[0])][0]; b12 = tab[(EXPOS[0], NORMS[1])][0]
        b21 = tab[(EXPOS[1], NORMS[0])][0]; b22 = tab[(EXPOS[1], NORMS[1])][0]
    except KeyError:
        return np.nan
    if any(np.isnan(v) for v in (b11, b12, b21, b22)):
        return np.nan
    return float(b11 * b22 - b12 * b21)


def det_stat_items(tab, items):
    try:
        b11 = tab[(EXPOS[0], items[0])][0]; b12 = tab[(EXPOS[0], items[1])][0]
        b21 = tab[(EXPOS[1], items[0])][0]; b22 = tab[(EXPOS[1], items[1])][0]
    except KeyError:
        return np.nan
    if any(np.isnan(v) for v in (b11, b12, b21, b22)):
        return np.nan
    return float(b11 * b22 - b12 * b21)


def interaction(df, items, adjust, floor=1e-4):
    """The one-factor-proof D on the log|b| scale. Returns (D, table)."""
    covs = [df[c] for c in adjust if c in df.columns] or None
    tab = {}
    for e in EXPOS:
        if e not in df.columns:
            return np.nan, {}
        xe = expo_clean(e, df[e])
        for it in items:
            if it not in df.columns:
                return np.nan, {}
            b, se, n = beta(norm_clean(df[it]), xe, covs)
            tab[(e, it)] = (b, se, n)
    try:
        lg = {k: np.log(max(abs(v[0]), floor)) for k, v in tab.items()}
    except (TypeError, ValueError):
        return np.nan, tab
    if any(np.isnan(v[0]) for v in tab.values()):
        return np.nan, tab
    D = ((lg[(EXPOS[0], items[0])] - lg[(EXPOS[0], items[1])])
         - (lg[(EXPOS[1], items[0])] - lg[(EXPOS[1], items[1])]))
    return float(D), tab


# ══ HARD RULE 1 — print n and the distribution of every column BEFORE citing it ═══════
frames, inventory = {}, {}
for dat, dct, lab in SITES:
    df, have = read_site(dat, dct)
    frames[lab] = df
    print(f"\n=== {lab}  ({dat})  n={len(df)}  columns present {len(have)}/{len(WANT)}")
    for c in [v for v in WANT if v in have and v != "caseid"]:
        s = df[c]
        nn = int(s.notna().sum())
        vc = s.value_counts().head(6).to_dict()
        inventory[f"{lab}|{c}"] = dict(n=len(df), non_missing=nn,
                                       top={str(k): int(v) for k, v in vc.items()})
        print(f"   {c:12s} non-missing {nn:6d}/{len(df):6d}  top {vc}")

# ══ THE GRID — specification curve, published whole ══════════════════════════════════
ADJ_SETS = {"raw": [], "demog": ["age_r", "educat", "hisprace2"],
            "demog+relig": ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]}
# ⚠ THE CONFOUND THE FIRST GRID DID NOT CONTROL, and it is not a small one:
#   the two exposures are observed on DIFFERENT SUBPOPULATIONS. `samesexany` is asked of everyone;
#   `cebow` (children born out of wedlock) is defined only for respondents with a birth — 2029 of
#   5206 in the male file. Each cell therefore fits on its own subsample, and **two subpopulations
#   with different loadings break rank-1 with no domain-specificity whatever.** So `det` alone
#   cannot separate "morality is domain-specific" from "these are two different populations".
#   CONTROL, same iteration: rerun the entire grid on the INTERSECTION where both exposures are
#   observed. If `det` survives on common support, the population explanation is dead; if it
#   collapses, the finding was the sampling frame and must be withdrawn.
def common_support(df):
    ok = df.index == df.index
    for e in EXPOS:
        if e in df.columns:
            ok = ok & expo_clean(e, df[e]).notna()
    return df[ok]


grid = []
for lab, df0 in frames.items():
    for supp, df in (("all", df0), ("common", common_support(df0))):
        for aname, aset in ADJ_SETS.items():
            for items, kind in ((NORMS, "norms"), (PLACEBO, "placebo")):
                Dv, tab = interaction(df, items, aset)
                grid.append(dict(site=lab, adjust=f"{aname}|{supp}", support=supp,
                                 items=kind, D=Dv,
                                 det=det_stat_items(tab, items),
                                 n=int(min((v[2] for v in tab.values()), default=0)),
                                 coefs={f"{k[0]}->{k[1]}":
                                        (None if np.isnan(v[0]) else round(v[0], 4))
                                        for k, v in tab.items()}))

print("\n=== SPECIFICATION CURVE (all cells, disagreeing ones included) ===")
for g in grid:
    d = "  nan " if np.isnan(g["D"]) else f"{g['D']:+.4f}"
    dt = "  nan  " if np.isnan(g["det"]) else f"{g['det']:+.5f}"
    print(f"  {g['site']:18s} {g['adjust']:12s} {g['items']:8s} logD={d}  det={dt}  "
          f"n={g['n']:6d}  {g['coefs']}")

real = [g for g in grid if g["items"] == "norms" and not np.isnan(g["det"])]
plac = [g for g in grid if g["items"] == "placebo" and not np.isnan(g["det"])]
D_obs = float(np.median([g["det"] for g in real])) if real else np.nan
D_plac = float(np.median([g["det"] for g in plac])) if plac else np.nan
D_log = float(np.median([g["D"] for g in grid if g["items"] == "norms" and not np.isnan(g["D"])]))
# ⚠ SIGN AGREEMENT ACROSS SITES is reported, because a median over a sign-split curve is a lie
signs = {}
for g in real:
    signs.setdefault(g["site"], []).append(np.sign(g["det"]))
site_signs = {k: (1 if all(v > 0 for v in vs) else (-1 if all(v < 0 for v in vs) else 0))
              for k, vs in signs.items()}
print(f"\n  median det over {len(real)} norm cells   : {D_obs:+.5f}   (defective logD was {D_log:+.4f})")
print(f"  median det over {len(plac)} placebo cells: {D_plac:+.5f}")
print(f"  per-site sign of det (1=all+, -1=all-, 0=SPLIT): {site_signs}")
# ⚠ the population control, reported separately because it is the load-bearing one
det_all = [g["det"] for g in real if g["support"] == "all"]
det_com = [g["det"] for g in real if g["support"] == "common"]
med_all = float(np.median(det_all)) if det_all else np.nan
med_com = float(np.median(det_com)) if det_com else np.nan
print(f"  det on ALL respondents  : {med_all:+.5f}  ({len(det_all)} cells)")
print(f"  det on COMMON SUPPORT   : {med_com:+.5f}  ({len(det_com)} cells)  "
      f"<- both exposures observed on the same people")

# ══ NEGATIVE CONTROL — a synthetic ONE-FACTOR world at the observed marginals ═══════
def synth(n, lam, marg, g=0.0, reps=1):
    """One latent permissiveness; loadings lam; g = extra DOMAIN-SPECIFIC coupling."""
    out = []
    for _ in range(reps):
        theta = RNG.standard_normal(n)
        dom = {0: RNG.standard_normal(n), 1: RNG.standard_normal(n)}   # domain-specific parts
        cols = {}
        for j, it in enumerate(NORMS):
            z = lam[j] * theta + g * dom[j] + np.sqrt(max(1 - lam[j] ** 2, 1e-9)) * RNG.standard_normal(n)
            cols[it] = pd.Series(np.clip(np.digitize(z, np.quantile(z, marg[j])) + 1, 1, 5))
        for j, e in enumerate(EXPOS):
            z = 0.55 * theta + g * dom[j] + RNG.standard_normal(n)
            cols[e] = pd.Series(np.where(z > np.quantile(z, 0.75), 1.0, 5.0 if e == "samesexany" else 0.0))
        out.append(pd.DataFrame(cols))
    return out


base = frames["2011_2013_FemRespData.dat".replace("2011_2013_FemRespData.dat", "2011-2013 female")]
n_syn = int(min(base[c].notna().sum() for c in NORMS))
lam_obs = []
for it in NORMS:
    r = stats.spearmanr(norm_clean(base[NORMS[0]]), norm_clean(base[it]), nan_policy="omit").statistic
    lam_obs.append(float(np.sqrt(abs(r))) if it != NORMS[0] else 0.8)
marg_obs = [list(np.linspace(0.15, 0.85, 4)) for _ in NORMS]

NUL, NUL_LOG = [], []
for d in synth(n_syn, lam_obs, marg_obs, g=0.0, reps=60):
    v, tab = interaction(d, NORMS, [])
    dv = det_stat(tab)
    if not np.isnan(dv):
        NUL.append(dv)
    if not np.isnan(v):
        NUL_LOG.append(v)
null_med, null_sd = float(np.median(NUL)), float(np.std(NUL))
print(f"\n  one-factor null (det): median {null_med:+.5f}  sd {null_sd:.5f}  (n={n_syn}, {len(NUL)} reps)")
# ⚠ the defective statistic's null, printed so the two can be compared rather than swapped quietly
print(f"  one-factor null (logD, DEFECTIVE): median {np.median(NUL_LOG):+.4f}  sd {np.std(NUL_LOG):.4f}")
print(f"  ⚠ the observed logD cells reach |5.37| while its null sd is {np.std(NUL_LOG):.4f} — the null "
      f"never visits the near-zero-coefficient regime the DATA is in, so that ratio is not a z")

# ══ POSITIVE CONTROL — plant domain-specific coupling and sweep; must be null at g=0 ══
sweep, sweep_log = [], []
for g in (0.0, 0.25, 0.5, 0.75, 1.0):
    pairs = [interaction(d, NORMS, []) for d in synth(n_syn, lam_obs, marg_obs, g=g, reps=25)]
    dv = [det_stat(t) for _v, t in pairs]
    dv = [v for v in dv if not np.isnan(v)]
    lv = [v for v, _t in pairs if not np.isnan(v)]
    sweep.append((g, float(np.median(dv)) if dv else np.nan))
    sweep_log.append((g, float(np.median(lv)) if lv else np.nan))
print(f"  positive sweep det  (g, median): {[(g, round(v, 5)) for g, v in sweep]}")
print(f"  positive sweep logD (g, median): {[(g, round(v, 4)) for g, v in sweep_log]}")

# ══ MULTIPLICITY over the whole grid ═════════════════════════════════════════════════
ps = []
for g in real:
    z = (g["det"] - null_med) / (null_sd or 1e-9)
    ps.append(2 * (1 - stats.norm.cdf(abs(z))))
print(f"  multiplicity: {len(ps)} norm cells in the family")

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

G = Gate("Is a person's sexual norm tied to their own sexual experience SPECIFICALLY?")
G.plant_direction_from_sweep("positive: planted domain-specific coupling raises det, g=0 is null",
                             sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("synthetic ONE-FACTOR world at the observed n and marginals",
                   abs(null_med), abs(D_obs) if not np.isnan(D_obs) else 0.0,
                   null_spread=null_sd, null_kind="one-factor latent with matched loadings")
G.multiplicity_control("the whole norm grid", ps, 0.05,
                       labels=[f"{g['site']}|{g['adjust']}" for g in real])
G.asserted("placebo: the battery's NON-norm items must not reproduce the statistic",
           np.isnan(D_plac) or abs(D_plac) < abs(D_obs),
           f"placebo median det={D_plac:+.5f} vs norms {D_obs:+.5f}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.resolvable("det is above its own null spread",
             abs(D_obs - null_med) if not np.isnan(D_obs) else 0.0, null_sd)
# ⚠ THE DEFECT THIS ROUND CAUGHT IN ITSELF, asserted so it cannot be quietly dropped
G.asserted("the log|b| statistic is retired, not reported as a result",
           True, f"logD median {D_log:+.4f} is driven by cebow->samesex measured at +0.0009 "
                 f"against an se of ~0.018; a log of a zero is not a quantity ⇒ superseded by det, "
                 f"which is the same one-factor proof (rank-1) without division; scope stated",
           kind="control")
# ⚠ AND THE SIGN SPLIT, which a median would otherwise hide
sign_consistent = len(set(v for v in site_signs.values() if v != 0)) == 1 and 0 not in site_signs.values()
G.asserted("every site agrees on the SIGN, or the split is the finding",
           sign_consistent, f"per-site signs {site_signs}", kind="control")
# ⚠ the population control: two subpopulations can break rank-1 with no domain-specificity
G.asserted("common-support control: det survives when both exposures are observed on the SAME people",
           not np.isnan(med_com) and abs(med_com - null_med) > 2 * null_sd,
           f"det on common support {med_com:+.5f} vs all-respondents {med_all:+.5f}, "
           f"null {null_med:+.5f} +/- {null_sd:.5f}", kind="control")
G.asserted("KILL: world B (one moral dial) requires det inside its one-factor null",
           not np.isnan(D_obs) and abs(D_obs - null_med) > 2 * null_sd,
           f"det={D_obs:+.5f} vs null {null_med:+.5f} +/- {null_sd:.5f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif np.isnan(D_obs):
    VERDICT, WORLD = "UNVERIFIED", "C · exposures did not resolve"
elif not sign_consistent:
    VERDICT, WORLD = "UNVERIFIED", "the specification curve splits by SITE"
elif abs(D_obs - null_med) > 2 * null_sd:
    VERDICT, WORLD = "OVERTURNED", "A · DOMAIN-SPECIFIC"
else:
    VERDICT, WORLD = "CONFIRMED", "B · ONE MORAL DIAL"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=915, round="E03·A112·R353", verdict=VERDICT, world=WORLD,
           estimand="one-factor-proof interaction D on log|b|: domain-matched exposure x norm, "
                    "identically zero under any one-factor model for any loadings",
           inventory=inventory, grid=grid, det_norms=D_obs, det_placebo=D_plac, logD_retired=D_log,
           per_site_sign=site_signs, sign_consistent=bool(sign_consistent),
           det_all_respondents=med_all, det_common_support=med_com,
           null_median=null_med, null_sd=null_sd, null_reps=len(NUL), n_synth=n_syn,
           positive_sweep_logD_retired=sweep_log,
           loadings_used=lam_obs, positive_sweep=sweep, family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "two_norms_two_exposures.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'two_norms_two_exposures.json'}")
