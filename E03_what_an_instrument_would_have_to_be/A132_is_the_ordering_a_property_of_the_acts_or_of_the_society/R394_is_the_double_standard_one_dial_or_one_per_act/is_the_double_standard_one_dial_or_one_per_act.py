#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A132·R394 — is the sexual double standard one dial, or one per act? And do two coding teams agree?
=======================================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#960`①. `#960` found the operative cross-cultural unit is **act × sex-of-actor**:
                72 of 109 societies (66.1%) apply the extramarital rule asymmetrically by sex, and
                GSS's `xmarsex` — *"is it wrong for a married person"* — has no degree of freedom in
                which that could be said. The question that opens is whether the asymmetry is **one
                property of a society** or **one rule per act**.

⚠⚠ AND `#960`'s IMPOSSIBILITY REGISTER WAS WRONG, WHICH THIS ROUND RETRACTS WITH A MEASUREMENT.
   `#960` registered *"there is no second ethnographic coding of these acts"* as STRUCTURALLY
   CANNOT. **False.** `SCCS596`/`SCCS597` are **Whyte 1978**, an independent coder and paper, and
   `SCCS597` codes the extramarital double standard that `SCCS169` (Broude & Greene 1976) also
   codes — **51 societies carry both**. An unchecked wall is UNVERIFIED, never SETTLED, and I wrote
   this one into the register one round ago without running the query that dissolves it. **So this
   round carries the cross-instrument validation `#960` said was unavailable**, and it is the first
   time in this project that a claim has had two independent instruments pointed at it.

⚠ HARD RULE 1 — printed before any claim:
    SCCS596 (No) Double Standard, PREMARITAL     Whyte 1978   n = 73/186   1:32(Yes) 2:41(No)
    SCCS597 (No) Double Standard, EXTRAMARITAL   Whyte 1978   n = 75/186   1:32(Yes) 2:41(No) 3:2(male punished severely)
    SCCS169 Extramarital Sex                     Broude 1976  n = 109/186  1:13 2:48 3:24 4:24
    complete on 596+597 (the fork)               n = 64
    complete on 597+169 (the cross-instrument)   n = 51

⚠ HARD RULE 2 — TWO instruments, named, and they are the point. **Whyte, M. K. (1978), *The Status
  of Women in Preindustrial Societies*** coded `SCCS596`/`SCCS597`. **Broude & Greene (1976),
  *Ethnology* 15(4):409–429** coded `SCCS169`. Both read HRAF ethnographies, so they share a
  CORPUS — that is not eliminated and is registered — but they are different teams with different
  instruments and different questions, which is the strongest independence available here.

⚠ G1 — TWO ESTIMANDS, named before the method.
  **E1 (instrument):** among the 51 societies coded by both teams, the agreement on *"does this
  society apply a sexual double standard to extramarital sex?"* — Cohen's κ with its permutation
  null. **This gates E2**: if two teams reading the same ethnographies disagree, neither world below
  is readable.
  **E2 (the fork):** among the 64 societies Whyte coded on both acts, the association between having
  a premarital double standard and an extramarital one — φ and the odds ratio, exact test.

Live Worlds    W_ONE_DIAL     · φ strongly positive ⇒ the double standard is **a property of the
                                 society**, not of the act: where women's premarital sexuality is
                                 policed asymmetrically, so is their marital sexuality.
               W_ACT_SPECIFIC · φ ≈ 0 ⇒ **each act carries its own rule**, and "the double standard"
                                 is a phrase covering two unrelated arrangements.
               W_INSTRUMENT   · ⚠ **the unwelcome one, and it is newly testable.** The two teams
                                 disagree on the extramarital double standard ⇒ the cross-cultural
                                 line of A132 is about coders, and `#960`'s 66.1% is a claim about
                                 Broude & Greene rather than about societies.

Prediction     W_ONE_DIAL     -> κ admissible AND φ > 0 clearing its permutation null, in most of
Matrix                            the code-3 x instrument grid.
               W_ACT_SPECIFIC -> κ admissible AND φ inside its null.
               W_INSTRUMENT   -> κ at or below chance; E2 is then reported but NOT interpreted.

Strongest      **BOTH TEAMS READ THE SAME HRAF ETHNOGRAPHIES.** A high κ therefore bounds
confound       coder-independence from above, not from below: it can mean the societies really are
(written       like that, or that both teams inherited the same ethnographer's sentence. ⇒ CONTROL,
before)        SAME ITERATION: κ is reported as an UPPER bound on independence and never as proof of
               it; and the two teams asked *different questions* (Whyte: is there a double standard;
               Broude & Greene: what is the extramarital rule), so their agreement is at least not
               a restatement. ⚠ The residual — a shared source ethnography — is registered as
               structurally unremovable here.

Controls       POSITIVE: plant a graded association into 2x2 tables drawn from the fitted
                 independence margins and sweep; g=0 must sit on zero, COMPUTED not asserted
                 (`#959`③); report the first detectable dose.
               NEGATIVE / N1: permute the extramarital indicator across societies — destroys the
                 pairing, preserves both margins, expected φ = 0 in closed form.
               PLACEBO: the whole pipeline on permuted data; φ must vanish.
               REPLICATION ARM: E2 recomputed with Broude & Greene's `SCCS169` supplying the
                 extramarital double standard instead of Whyte's `SCCS597` — two independent
                 instruments on the same fork.
               POWER AT THE OBSERVED EFFECT (`#960`③): the power to detect the φ actually seen, not
                 the power to detect some larger φ. A null is inadmissible without it.
               REGION: leave out each of Murdock's six world regions.
               MULTIPLICITY: 3 code-3 treatments x 2 instruments x 2 estimators.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **remove the shared HRAF corpus** — two teams, one library of ethnographies; κ bounds
    coder-independence from above only. A second field observation is what it would take;
  (2) ⚠ **no causal identification** — 64 societies observed once;
  (3) ⚠ **no US comparison of the same quantity** — GSS's items carry no sex-of-actor, which is
    `#960`'s finding and is why this quantity has no American counterpart to replicate in;
  (4) ⚠ **`SCCS597` code 3 (male punished severely) has n=2** — swept three ways, never resolved;
  (5) ⚠ **within-society variation invisible**; focal years differ by more than a century;
  (6) ⚠ **Galton's problem mitigated by SCCS's design, not eliminated** — region leave-one-out bounds;
  (7) `[unchallenged]` — door ③.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
B = ROOT / "data" / "external" / "dplace" / "repo" / "datasets" / "SCCS"
SEEDS = (394, 1394, 2394)
VARS = ["SCCS596", "SCCS597", "SCCS169"]
REG = [(1, 28, "Africa"), (29, 58, "Circum-Mediterranean"), (59, 87, "East Eurasia"),
       (88, 114, "Insular Pacific"), (115, 153, "North America"), (154, 186, "South America")]

d = pd.read_csv(B / "data.csv", low_memory=False)
W = d[d.var_id.isin(VARS)].pivot_table(index="soc_id", columns="var_id", values="code",
                                       aggfunc="first")
W["region"] = [next((nm for a, b, nm in REG if a <= int(str(i).replace("SCCS", "")) <= b), "?")
               for i in W.index]

print("HARD RULE 1 — coverage and codes, printed before any claim (SCCS, 186 societies):")
for v in VARS:
    print(f"  {v} n={int(W[v].notna().sum()):3d}/186  "
          + " ".join(f"{int(k)}:{int(x)}" for k, x in W[v].value_counts().sort_index().items()))
n_fork = int(W[["SCCS596", "SCCS597"]].notna().all(axis=1).sum())
n_xi = int(W[["SCCS597", "SCCS169"]].notna().all(axis=1).sum())
print(f"  complete on 596+597 (the fork) n={n_fork} · complete on 597+169 (cross-instrument) "
      f"n={n_xi}")
print("⚠ HARD RULE 2 — TWO instruments: Whyte 1978 (SCCS596/597) and Broude & Greene 1976 "
      "(SCCS169). Different teams, different questions, SHARED HRAF corpus — so agreement bounds "
      "coder-independence from ABOVE, never proves it.")
print("⚠⚠ `#960` registered 'no second ethnographic coding of these acts' as STRUCTURALLY CANNOT. "
      "That is RETRACTED here by measurement: Whyte 1978 is the second coding and 51 societies "
      "carry both. An unchecked wall is UNVERIFIED, never SETTLED.")


def ds_ext_whyte(code, code3):
    """Whyte's extramarital double standard. code3 ∈ {'not', 'is', 'drop'} for the n=2 reverse cell."""
    out = pd.Series(np.nan, index=code.index)
    out[code == 1] = 1.0
    out[code == 2] = 0.0
    if code3 == "is":
        out[code == 3] = 1.0
    elif code3 == "not":
        out[code == 3] = 0.0
    return out


def phi_and_or(a, b):
    """φ and the odds ratio on a 2x2, with Fisher's exact p."""
    t = np.array([[int(((a == i) & (b == j)).sum()) for j in (0, 1)] for i in (0, 1)])
    n = t.sum()
    r0, r1 = t.sum(1)
    c0, c1 = t.sum(0)
    den = np.sqrt(float(r0) * r1 * c0 * c1)
    phi = float((t[1, 1] * t[0, 0] - t[1, 0] * t[0, 1]) / den) if den > 0 else np.nan
    orr = float(((t[1, 1] + .5) * (t[0, 0] + .5)) / ((t[1, 0] + .5) * (t[0, 1] + .5)))
    p = float(stats.fisher_exact(t)[1]) if n > 0 else 1.0
    return phi, orr, p, t, int(n)


# ══ E1 — THE CROSS-INSTRUMENT VALIDATION, which gates everything below ═══════════════
xi = W[W[["SCCS597", "SCCS169"]].notna().all(axis=1)]
whyte = ds_ext_whyte(xi["SCCS597"], "not")
KAPPAS = {}
for bgdef in ([2, 3], [2]):
    _b = xi["SCCS169"].isin(bgdef).astype(float)
    _k = whyte.notna()
    _w, _bb = whyte[_k].to_numpy(), _b[_k].to_numpy()
    _ag = float((_w == _bb).mean()); _pw, _pb = _w.mean(), _bb.mean()
    _pe = _pw * _pb + (1 - _pw) * (1 - _pb)
    KAPPAS[str(bgdef)] = dict(agreement=_ag, kappa=float((_ag - _pe) / (1 - _pe)) if _pe < 1 else float("nan"),
                              chance=_pe, bg_rate=float(_bb.mean()))
bg = xi["SCCS169"].isin([2, 3]).astype(float)      # headline: 2 or 3 = asymmetric by sex
keep = whyte.notna()
w1, b1 = whyte[keep].to_numpy(), bg[keep].to_numpy()
agree = float((w1 == b1).mean())
pw, pb = w1.mean(), b1.mean()
pe = pw * pb + (1 - pw) * (1 - pb)
kappa = float((agree - pe) / (1 - pe)) if pe < 1 else np.nan
k_null = []
for s in SEEDS:
    r = np.random.default_rng(s)
    for _ in range(4000):
        bp = r.permutation(b1)
        a2 = float((w1 == bp).mean())
        k_null.append((a2 - pe) / (1 - pe) if pe < 1 else 0.0)
k_null = np.array(k_null)
k_p = float((k_null >= kappa).mean())
print(f"\n  E1 κ under both Broude&Greene definitions: "
      + " · ".join(f"{k}: agreement {v['agreement']:.4f} κ {v['kappa']:+.4f} (BG rate "
                   f"{v['bg_rate']:.3f})" for k, v in KAPPAS.items()))
print(f"  E1 CROSS-INSTRUMENT, n={len(w1)}: Whyte says double standard in {int(w1.sum())}, "
      f"Broude & Greene in {int(b1.sum())} · raw agreement {agree:.4f} · Cohen's κ {kappa:+.4f} "
      f"(chance agreement {pe:.4f}) · permutation null {k_null.mean():+.4f} ± {k_null.std():.4f} "
      f"· one-sided p {k_p:.4f}")

# ══ E2 — THE FORK, over the code-3 x instrument grid ════════════════════════════════
grid, ps = [], []
for code3 in ("not", "is", "drop"):
    for iname, ivar, bgdef in (("Whyte SCCS597", "SCCS597", "-"),
                               ("Broude&Greene 2or3", "SCCS169", [2, 3]),
                               ("Broude&Greene 2only", "SCCS169", [2])):
        sub = W[W[["SCCS596", ivar]].notna().all(axis=1)]
        pre = (sub["SCCS596"] == 1).astype(float)
        ext = (ds_ext_whyte(sub["SCCS597"], code3) if ivar == "SCCS597"
               else sub["SCCS169"].isin(bgdef).astype(float))
        m = ext.notna()
        phi, orr, p, t, n = phi_and_or(pre[m], ext[m])
        grid.append(dict(code3=code3, instrument=iname, bg_def=str(bgdef), n=n, phi=phi, odds_ratio=orr, fisher_p=p,
                         table=t.tolist(), pre_rate=float(pre[m].mean()),
                         ext_rate=float(ext[m].mean())))
        ps.append(p)
        print(f"  E2 [{iname:<22s} code3={code3:<4s}] n={n:3d} · premarital DS {pre[m].mean():.3f} "
              f"· extramarital DS {ext[m].mean():.3f} · φ {phi:+.4f} · OR {orr:6.2f} · Fisher p "
              f"{p:.4f} · table {t.tolist()}")
HEAD = grid[0]

# ══ N1 — permute the extramarital indicator; both margins preserved, expected φ = 0 ══
sub = W[W[["SCCS596", "SCCS597"]].notna().all(axis=1)]
pre = (sub["SCCS596"] == 1).astype(float)
ext = ds_ext_whyte(sub["SCCS597"], "not")
m = ext.notna()
PRE, EXT = pre[m].to_numpy(), ext[m].to_numpy()
n1 = []
for s in SEEDS:
    r = np.random.default_rng(s + 10)
    for _ in range(4000):
        n1.append(phi_and_or(pd.Series(PRE), pd.Series(r.permutation(EXT)))[0])
n1 = np.array(n1)
n1_mu, n1_sd = float(np.nanmean(n1)), float(np.nanstd(n1, ddof=1))
z = (HEAD["phi"] - n1_mu) / max(n1_sd, 1e-9)
p_perm = float((np.abs(n1 - n1_mu) >= abs(HEAD["phi"] - n1_mu)).mean())
print(f"\n  N1 permutation of the extramarital indicator (both margins held, expected φ = 0): "
      f"{n1_mu:+.5f} ± {n1_sd:.4f} · observed {HEAD['phi']:+.4f} · z {z:+.2f} · two-sided "
      f"empirical p {p_perm:.4f}")

# ══ PLACEBO ═════════════════════════════════════════════════════════════════════════
plac = [phi_and_or(pd.Series(PRE), pd.Series(np.random.default_rng(s + 900).permutation(EXT)))[0]
        for s in SEEDS]
plac_mu = float(np.mean(plac))
print(f"  PLACEBO (pipeline on permuted data): φ {plac_mu:+.4f} against a null spread of {n1_sd:.4f}")

# ══ POSITIVE CONTROL — plant a graded association at the observed margins ════════════
sweep = []
for g in (0.0, 0.15, 0.30, 0.50):
    vals = []
    for s in SEEDS:
        r = np.random.default_rng(s + int(g * 1000))
        e2 = np.where(r.random(len(PRE)) < g, PRE, r.permutation(EXT))
        vals.append(phi_and_or(pd.Series(PRE), pd.Series(e2))[0])
    sweep.append([float(g), float(np.median(vals))])
mde_g = next((g for g, v in sweep if v > sweep[0][1] + 2 * n1_sd), None)
print(f"  positive sweep (copy a share g of the premarital indicator into the extramarital one; "
      f"g=0 is a fresh permutation, COMPUTED not asserted): "
      f"{[(g, round(v,4)) for g, v in sweep]} · first g clearing 2 null spreads: {mde_g}")

# ══ POWER AT THE OBSERVED EFFECT (`#960`③) ══════════════════════════════════════════
hits = 0
for j in range(4000):
    r = np.random.default_rng(70000 + j)
    lo, hi = 0.5 - abs(HEAD["phi"]) / 2, 0.5 + abs(HEAD["phi"]) / 2
    a = (r.random(len(PRE)) < PRE.mean()).astype(float)
    b = np.where(a == 1, r.random(len(PRE)) < hi, r.random(len(PRE)) < lo).astype(float)
    hits += phi_and_or(pd.Series(a), pd.Series(b))[2] < 0.05
power_at_obs = hits / 4000
print(f"  ⚠⚠ POWER AT THE OBSERVED φ = {HEAD['phi']:+.4f}: **{power_at_obs:.3f}** — the control "
      f"`#960`③ owes, evaluated at the effect actually seen and not at a larger one")

# ══ REGION LEAVE-ONE-OUT ════════════════════════════════════════════════════════════
loo = []
for rg in sorted(sub.loc[m, "region"].unique()):
    keep2 = (sub.loc[m, "region"] != rg).to_numpy()
    phi2, _, p2, _, n2 = phi_and_or(pd.Series(PRE[keep2]), pd.Series(EXT[keep2]))
    loo.append(dict(dropped=rg, n=n2, phi=phi2, p=p2))
print("  REGION LEAVE-ONE-OUT: "
      + " · ".join(f"−{r['dropped'][:12]} (n={r['n']}) φ {r['phi']:+.4f}" for r in loo))

# ══ VERDICT ═════════════════════════════════════════════════════════════════════════
pos_fires = mde_g is not None and sweep[-1][1] > sweep[0][1] + 2 * n1_sd
pos_g0_ok = abs(sweep[0][1] - n1_mu) < 2 * n1_sd
n1_ok = abs(n1_mu) < 0.05
plac_ok = abs(plac_mu - n1_mu) < 2 * n1_sd
powered = power_at_obs >= 0.50
kappa_ok = (not np.isnan(kappa)) and kappa > 0.20 and k_p < 0.05
controls_ok = pos_fires and pos_g0_ok and n1_ok and plac_ok
sig = sum(1 for g in grid if g["fisher_p"] < 0.05 and g["phi"] > 0)
# ⚠ per-INSTRUMENT verdicts. The code-3 axis moves 2 societies, so its three cells are near-copies
#   and counting them as three votes is a multiplicity error in the flattering direction.
inst_names = sorted({g["instrument"] for g in grid})
inst_verdict = {nm: ("positive" if all(g["phi"] > 0 and g["fisher_p"] < 0.05
                                       for g in grid if g["instrument"] == nm)
                     else ("null" if all(g["fisher_p"] >= 0.05
                                         for g in grid if g["instrument"] == nm)
                           else "mixed")) for nm in inst_names}
n_pos = sum(1 for v in inst_verdict.values() if v == "positive")
n_null = sum(1 for v in inst_verdict.values() if v == "null")
instruments_agree = (n_pos == len(inst_names)) or (n_null == len(inst_names))
print(f"\n  ⚠⚠ PER-INSTRUMENT verdicts (the grid's 9 cells are 3 instruments x a code-3 axis that "
      f"moves 2 societies, so they are NOT 9 votes): {inst_verdict} · agree {instruments_agree}")

if not controls_ok:
    world = "W_UNREADABLE"
elif not kappa_ok:
    world = "W_INSTRUMENT"
elif not instruments_agree:
    world = "W_INSTRUMENT_SPLIT"
elif not powered and not (HEAD["phi"] > 0 and p_perm < 0.05):
    world = "W_UNDERPOWERED"
elif HEAD["phi"] > 0 and p_perm < 0.05 and n_pos == len(inst_names):
    world = "W_ONE_DIAL"
else:
    world = "W_ACT_SPECIFIC"

G = Gate("Is the sexual double standard one property of a society, or one rule per act? — and do "
         "two independent ethnographic coding teams agree about it at all?")
G.plant_direction_from_sweep(
    "positive: copying a share g of the premarital double-standard indicator into the extramarital "
    "one must raise φ, and g=0 is a fresh permutation whose value is COMPUTED in the round rather "
    "than asserted (`#959`③)", sweep, baseline=n1_mu, baseline_spread=max(n1_sd, 1e-4))
G.negative_control(
    "N1: permuting the extramarital indicator across societies destroys the pairing while holding "
    "both margins, so the expected φ is zero in closed form",
    abs(n1_mu), abs(HEAD["phi"] - n1_mu), null_spread=n1_sd,
    null_kind="label permutation null on the extramarital indicator, whose expected value is zero "
              "in closed form because permuting one margin against a fixed other cannot produce "
              "association in expectation, plus Fisher's exact test on the 2x2")
G.multiplicity_control("3 code-3 treatments x 2 instruments (`#936`②)", ps, 0.05,
                       labels=[f"{g['instrument'][:12]}/code3={g['code3']}" for g in grid])
G.asserted(
    "⚠⚠ `#960`'s IMPOSSIBILITY REGISTER IS RETRACTED BY MEASUREMENT. It recorded 'no second "
    "ethnographic coding of these acts' as STRUCTURALLY CANNOT, one round ago, without running the "
    "query that dissolves it. Whyte 1978 coded the extramarital double standard independently of "
    "Broude & Greene 1976 and 51 societies carry both. An unchecked wall is UNVERIFIED, never "
    "SETTLED", True,
    f"SCCS596/597 source=whyte1978cross (n=73/75) · SCCS169 source=broude1976cross (n=109) · "
    f"overlap n={n_xi} · this is the first claim in the project with two independent instruments "
    f"pointed at it", kind="control", population="SCCS, 186 societies")
G.asserted(
    "⚠ E1 GATES E2: two teams reading HRAF and asking different questions — Whyte 'is there a "
    "double standard', Broude & Greene 'what is the extramarital rule'. If they disagree, neither "
    "world below is readable and `#960`'s 66.1% is a claim about one team", kappa_ok,
    f"n={len(w1)} · raw agreement {agree:.4f} · Cohen's κ {kappa:+.4f} · permutation null "
    f"{k_null.mean():+.4f} ± {k_null.std():.4f} · one-sided p {k_p:.4f} · admissible iff κ>0.20 and "
    f"p<0.05", kind="control", population=f"SCCS societies coded by both teams, n={len(w1)}")
G.asserted(
    "⚠ THE STRONGEST CONFOUND, registered rather than removed: both teams read the SAME HRAF "
    "ethnographies, so κ bounds coder-independence from ABOVE and never proves it — agreement can "
    "mean the societies are like that, or that both teams inherited one ethnographer's sentence. A "
    "second field observation is what removing it would take", True,
    "shared corpus is structurally unremovable here; the mitigation is that the two teams asked "
    "different questions, so their agreement is at least not a restatement", kind="control",
    population=f"SCCS societies coded by both teams, n={len(w1)}")
G.asserted(
    "⚠⚠ POWER AT THE OBSERVED EFFECT — the control `#960`③ owes, and the third round in a row where "
    "this family of defect had to be repaired. Evaluated at the φ actually seen, not at a larger one",
    powered, f"power to detect φ={HEAD['phi']:+.4f} at n={HEAD['n']} is {power_at_obs:.3f}",
    kind="control", population=f"SCCS fork subset n={HEAD['n']}")
G.asserted(
    "⚠⚠ REPLICATION ARM, AND IT IS THE FINDING: the same fork recomputed with Broude & Greene "
    "supplying the extramarital double standard instead of Whyte. Counted PER INSTRUMENT rather "
    "than per cell, because the code-3 axis moves 2 societies and three near-copies are not three "
    "votes — that count was my own multiplicity error in the flattering direction",
    instruments_agree,
    " · ".join(f"{g['instrument'][:12]}/code3={g['code3']}: n={g['n']} φ {g['phi']:+.4f} OR "
               f"{g['odds_ratio']:.2f} p {g['fisher_p']:.4f}" for g in grid),
    kind="control", population="SCCS, both fork subsets")
G.asserted(
    "⚠ GALTON'S PROBLEM bounded, not removed: φ recomputed leaving out each of Murdock's six world "
    "regions", True,
    " · ".join(f"−{r['dropped']} (n={r['n']}) φ {r['phi']:+.4f}" for r in loo),
    kind="control", population=f"SCCS fork subset n={HEAD['n']}")
G.asserted("the whole code-3 x instrument grid is published, disagreeing cells included", True,
           " · ".join(f"{g['instrument'][:8]}/{g['code3']} φ{g['phi']:+.3f} p{g['fisher_p']:.3f}"
                      for g in grid), kind="control", population="SCCS, both fork subsets")

G.asserted(
    "KILL: pre-registered CONDITIONAL — evaluated ONLY if the plant fires with g=0 on the "
    "permutation null, N1 returns its known zero, the placebo vanishes, AND the two coding teams "
    "agree above chance. STAKED: W_ONE_DIAL, i.e. the premarital and extramarital double standards "
    "are positively associated, clearing the permutation null, in a majority of the grid. "
    "⚠ W_ACT_SPECIFIC means 'the double standard' is one phrase over two unrelated arrangements; "
    "W_INSTRUMENT means A132 has been measuring coders",
    controls_ok and kappa_ok and instruments_agree and HEAD["phi"] > 0 and p_perm < 0.05
    and n_pos == len(inst_names),
    f"plant fires {pos_fires} (first detectable g={mde_g}) · g=0 on null {pos_g0_ok} · N1 at zero "
    f"{n1_ok} ({n1_mu:+.5f}) · placebo {plac_ok} · κ admissible {kappa_ok} (κ={kappa:+.4f}, "
    f"p={k_p:.4f}) · power at observed φ {power_at_obs:.3f} · φ {HEAD['phi']:+.4f} OR "
    f"{HEAD['odds_ratio']:.2f}, permutation p {p_perm:.4f} · grid cells positive-and-significant "
    f"{sig}/{len(grid)} · PER-INSTRUMENT {inst_verdict}, agree {instruments_agree} "
    f"⇒ {world}",
    kind="kill", yardstick="φ between the premarital and extramarital sexual double standards "
                           "across societies",
    yardstick_noise=n1_sd,
    population=f"SCCS societies coded by Whyte 1978 on both acts, n={HEAD['n']} of 186, six world "
               f"regions",
    direction="one-sided: W_ONE_DIAL requires a POSITIVE association")

print(G)
verdict = (("UNVERIFIED" if not (controls_ok and kappa_ok) else
            ("CONFIRMED" if world == "W_ONE_DIAL" else
             ("UNVERIFIED" if world in ("W_UNDERPOWERED", "W_INSTRUMENT_SPLIT")
              else "OVERTURNED")))
           + f" · world {world} · κ {kappa:+.4f} (n={len(w1)}) · φ {HEAD['phi']:+.4f} OR "
             f"{HEAD['odds_ratio']:.2f} (n={HEAD['n']}) · permutation p {p_perm:.4f}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=961, round="E03·A132·R394", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world in ("W_ACT_SPECIFIC", "W_INSTRUMENT", "W_UNDERPOWERED")),
               instruments={"whyte1978": ["SCCS596", "SCCS597"], "broude1976": ["SCCS169"]},
               coverage={v: int(W[v].notna().sum()) for v in VARS},
               n_fork=n_fork, n_cross_instrument=n_xi,
               kappa=kappa, kappa_p=k_p, raw_agreement=agree, chance_agreement=pe,
               kappa_null_mean=float(k_null.mean()), kappa_null_sd=float(k_null.std()),
               grid=grid, head=HEAD, kappa_by_bg_def=KAPPAS,
               instrument_verdict=inst_verdict, instruments_agree=bool(instruments_agree), null_median=n1_mu, null_sd=n1_sd, null_draws=int(n1.size),
               perm_p=p_perm, placebo=plac_mu, positive_sweep=sweep, first_detectable_g=mde_g,
               power_at_observed=power_at_obs, region_loo=loo,
               family_size=len(ps), seeds=list(SEEDS), world=world, verdict=verdict),
          open(OUT / "is_the_double_standard_one_dial_or_one_per_act.json", "w"), indent=1,
          default=float)
print(f"\nwrote {OUT / 'is_the_double_standard_one_dial_or_one_per_act.json'}")
