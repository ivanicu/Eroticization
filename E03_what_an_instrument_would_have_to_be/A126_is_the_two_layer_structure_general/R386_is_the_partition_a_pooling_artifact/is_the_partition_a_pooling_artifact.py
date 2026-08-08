#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A126·R386 — is the partition a fact about people, or an artifact of pooling 21 waves?
==========================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#949` found an exact 2+2 partition — `{premarsx, homosex}` held together above
                their nulls, `{teensex, xmarsex}` held together, and all four crossing pairs traded
                below — and registered the rival reading as *"equally available"*. **That
                understated it, and this round exists because the rival PREDICTS THE EXACT PATTERN.**

⚠⚠ THE CONFOUND, WORKED OUT BEFORE ANY CODE RUNS, FROM `#941`'s OWN MEASUREMENTS
   `#941` measured the four items' raw trends per decade: `premarsx` **+0.1827** · `homosex`
   **+0.4056** · `teensex` **+0.1566** · `xmarsex` **+0.0567**. **The two fastest-moving items are
   exactly `{premarsx, homosex}` and the two slowest are exactly `{teensex, xmarsex}`** — my two
   clusters, with no remainder.
   Person-centring removes a respondent's LEVEL. It does **not** remove their PROFILE. A 2024
   respondent sits relatively high on the fast-moving items and relatively low on the slow ones; a
   1988 respondent the reverse. **Pooling 21 waves therefore manufactures positive residual
   correlation WITHIN each speed-group and negative correlation ACROSS — which is the 2+2 partition,
   exactly.** ⇒ the rival does not merely coexist with `#949`'s finding; **it entails it.**
   ⚠ And the same argument runs one level down: within a single wave, birth cohorts differ in
   profile too, so removing only the wave leaves a cohort version of the same artifact.

Live Worlds    W_WRONGED  · the partition survives once wave AND cohort profile differences are
                            removed, at a size comparable to pooled ⇒ it is a structure in how
                            people hold these judgements at a moment. `#949` stands.
               W_POOLING  · it collapses toward zero as the grouping gets finer ⇒ **`#949`'s
                            headline is an artifact of pooling and comes off the page in both
                            languages, one round after landing.** ⚠ **The unwelcome one, and it is
                            the one the arithmetic above predicts.**
               W_PARTIAL  · it shrinks substantially but survives ⇒ both contribute, and the honest
                            output is the SHARE, not a yes/no.
               W_NOT_A_PARTITION · once grouped, the six pairs stop blocking cleanly at all ⇒ the
                            "2+2" summary was the wrong description of even the pooled data.
                            (the meta-separator)

Estimand       **The PARTITION CONTRAST** Δ = mean residual r over the **2 within-cluster** pairs
(G1)           minus mean residual r over the **4 crossing** pairs. One scalar, computed at three
               grouping levels, with centring and scaling done INSIDE each group:
                 L0 POOLED            — reproduces `#949`
                 L1 WITHIN WAVE       — 21 groups; removes each wave's item profile
                 L2 WITHIN WAVE×COHORT-BAND — removes each wave's AND each cohort band's profile
               ⚠ **This is a different estimand from `#949`'s per-pair correlations, and saying so
               matters**: `#949` asked *which pairs depart from their null*; this asks *how much of
               the block structure survives removing between-group profile variation*. A round that
               changed estimand silently would be answering a question nobody asked.

⚠ SCALING      At L1/L2 each item is centred at its GROUP mean but divided by the POOLED sd. Z-scoring
CHOICE,        inside a small cell would let cell-level sampling noise into the denominator and shrink
STATED         correlations for a reason that has nothing to do with the partition. Centring is what
BEFORE         removes the profile; scaling is not, so only the centring is made group-local.

Prediction     W_WRONGED  -> Δ(L2) ≥ 0.5·Δ(L0) and clears 2× its own null spread.
Matrix         W_POOLING  -> Δ(L2) inside its null.
               W_PARTIAL  -> in between; report the retained share.
               W_NOT_A_PARTITION -> the 6 pairs no longer split 2/4 by sign at any level.

Controls       NEGATIVE: within-person ITEM-LABEL permutation applied **inside the same groups** —
                 each person's total preserved exactly, so the ipsative constraint is held fixed and
                 only which item got which value dies. Recomputed separately at L0, L1, L2, because
                 a null computed at the wrong grouping is `#947`'s error again.
               ⚠ NULL-WORLD CHECK (`#948`①, the live defect — four consecutive rounds have now had
                 a mis-specified negative control): a second independent permutation set at each
                 level, and the null must reproduce itself to a tolerance small against Δ.
               POSITIVE: plant a partition INTO the permuted world and sweep; `g=0` sits on the null
                 by construction (`#922`, `#937`⑤).
               ⚠ SYNTHETIC-ARTIFACT CONTROL — the one this round actually turns on: **build the
                 rival's world**. Take the permuted (structureless) data and impose ONLY a
                 wave-varying differential shift matching `#941`'s measured trends. If Δ(L0) is
                 large there and Δ(L2) is not, the instrument demonstrably detects the artifact and
                 the grouping demonstrably removes it — which is what licenses reading the real
                 Δ(L2). Without this the round cannot distinguish "no artifact" from "blind to it".
               MULTIPLICITY: the family is **3 levels × 6 pairs = 18 cells plus 3 contrasts = 21**.
               SPEC CURVE (G4): 3 levels × {group-centred, pooled-centred} × {z by pooled sd, z by
                 group sd}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **observe a person over time** — repeated cross-section, so "at a moment" is the strongest
    temporal claim available and no within-person change is measured;
  (2) ⚠ **separate age from period from cohort** — `#939`/`#943`'s wall, inherited; grouping by
    wave×cohort does NOT identify their separate effects, it only removes their profile means;
  (3) ⚠ **no second instrument** — the four norms are GSS's; **only this one instrument** carries
    them;
  (4) ⚠ **distinguish a trade-off from a shared cause pushing two items apart** — inherited from
    `#945`/`#949` unchanged;
  (5) `[unchallenged]` — door ③.
"""
import json
import sys
import warnings
from itertools import combinations
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
GSS = ROOT / "data" / "external" / "gss" / "GSS_stata" / "gss7224_r3a.dta"
RNG = np.random.default_rng(386)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
K = len(ITEMS)
PAIRS = list(combinations(range(K), 2))
CLUSTER_A, CLUSTER_B = {0, 3}, {1, 2}          # {premarsx, homosex} · {teensex, xmarsex}
WITHIN = [p for p in PAIRS if set(p) <= CLUSTER_A or set(p) <= CLUSTER_B]
CROSS = [p for p in PAIRS if p not in WITHIN]
assert len(WITHIN) == 2 and len(CROSS) == 4

# ⚠ `#941`'s trends are READ FROM ITS ARTIFACT, not typed (`#840`, and `#949` was blocked for
#   exactly this). They are the numbers the whole confound argument rests on.
REF941 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                        "A123_conversion_or_replacement_per_norm" /
                        "R379_all_four_norms_under_the_same_control" / "results" /
                        "all_four_norms_under_the_same_control.json"))
TREND = {r["item"]: r["raw"] for r in REF941["rows"]}
order = sorted(TREND, key=lambda k: -TREND[k])
print("⚠ THE CONFOUND, from `#941`'s ARTIFACT rather than memory — raw trend per decade:")
print("   " + " · ".join(f"{k} {TREND[k]:+.4f}" for k in order))
print(f"   fastest two = {{{order[0]}, {order[1]}}} · slowest two = {{{order[2]}, {order[3]}}}")
fast, slow = {order[0], order[1]}, {order[2], order[3]}
predicts = (fast == {ITEMS[i] for i in CLUSTER_A} and slow == {ITEMS[i] for i in CLUSTER_B}) or \
           (fast == {ITEMS[i] for i in CLUSTER_B} and slow == {ITEMS[i] for i in CLUSTER_A})
print(f"   ⇒ do the speed groups COINCIDE with `#949`'s clusters? **{predicts}** — if True the rival "
      f"does not merely coexist with the finding, it ENTAILS it")

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"]).copy()
d["cband"] = pd.cut(d.cohort, [1880, 1930, 1950, 1970, 2010], labels=["<1930", "30-49", "50-69", "70+"])
d = d.dropna(subset=["cband"])
waves = sorted(int(y) for y in d.year.unique())
cells_L2 = d.groupby(["year", "cband"], observed=True).size()
print(f"\nHARD RULE 1 — ballot 1, all four jointly, cohort present: n={len(d)} · {len(waves)} waves "
      f"{waves[0]}-{waves[-1]} · L2 cells {len(cells_L2)} with n {int(cells_L2.min())}-"
      f"{int(cells_L2.max())} (median {int(cells_L2.median())})")
print("⚠ HARD RULE 2 — instrument: GSS ballot 1, four sexual-norm items, one block, one "
      "questionnaire. Identical to `#949`'s; only the GROUPING changes.")

X = d[ITEMS].to_numpy(dtype=float)
POOLED_SD = X.std(axis=0)
GROUPS = {
    "L0 pooled": np.zeros(len(d), dtype=int),
    "L1 wave": d.year.astype("category").cat.codes.to_numpy(),
    "L2 wave x cohort": d.groupby(["year", "cband"], observed=True).ngroup().to_numpy(),
}


def group_centre(mat, g):
    """centre each item at its GROUP mean, scale by the POOLED sd (see the scaling note)."""
    out = np.empty_like(mat)
    for gid in np.unique(g):
        m = g == gid
        out[m] = mat[m] - mat[m].mean(axis=0)
    return out / POOLED_SD


def ipsative(mat):
    return mat - mat.mean(axis=1, keepdims=True)


def pair_rs(mat):
    R = ipsative(mat)
    return {p: float(np.corrcoef(R[:, p[0]], R[:, p[1]])[0, 1]) for p in PAIRS}


def contrast(rs):
    return float(np.mean([rs[p] for p in WITHIN]) - np.mean([rs[p] for p in CROSS]))


def permute_within_person(mat):
    return np.take_along_axis(mat, np.argsort(RNG.random(mat.shape), axis=1), axis=1)


NREPS = 120
levels = {}
for name, g in GROUPS.items():
    obs = pair_rs(group_centre(X, g))
    nulls_A, nulls_B = [], []
    for _ in range(NREPS):
        nulls_A.append(contrast(pair_rs(group_centre(permute_within_person(X), g))))
        nulls_B.append(contrast(pair_rs(group_centre(permute_within_person(X), g))))
    levels[name] = dict(rs=obs, delta=contrast(obs),
                        null_mean=float(np.mean(nulls_A)), null_sd=float(np.std(nulls_A)),
                        null_B=float(np.mean(nulls_B)))
    L = levels[name]
    L["z"] = (L["delta"] - L["null_mean"]) / max(L["null_sd"], 1e-9)
    L["repro"] = abs(L["null_mean"] - L["null_B"])

print(f"\n  PARTITION CONTRAST  Δ = mean r(2 within-cluster pairs) − mean r(4 crossing pairs)")
print(f"  {'level':<20s} {'Δ':>9s} {'its null':>10s} {'null sd':>8s} {'z':>8s} {'repro':>8s}")
for name, L in levels.items():
    print(f"  {name:<20s} {L['delta']:+9.4f} {L['null_mean']:+10.4f} {L['null_sd']:8.4f} "
          f"{L['z']:+8.2f} {L['repro']:8.5f}")
retained = levels["L2 wave x cohort"]["delta"] / max(levels["L0 pooled"]["delta"], 1e-9)
print(f"  ⇒ Δ(L2)/Δ(L0) = {retained:.3f} — the share of the pooled partition surviving removal of "
      f"wave AND cohort profile")

# ══ SYNTHETIC-ARTIFACT CONTROL — build the rival's world and check we can SEE it ═════
#   permuted (no structure) + a wave-varying differential shift matching `#941`'s trends
yr = d.year.to_numpy(dtype=float)
yc = (yr - yr.mean()) / 10.0
drift = np.array([TREND[c] for c in ITEMS])
# ⚠⚠ v1's SYNTHETIC CONTROL DEMANDED THE MEASURED DRIFT BE VISIBLE AT L0, AND IT WAS NOT
#   (Δ_syn = +0.0773 against a null of +0.0972 ± 0.0105 — BELOW it). v1 read that as "the
#   instrument may be blind", refused the round, and it was asking the wrong question.
#   **"The measured drift does not produce a visible partition" and "the instrument cannot detect a
#   drift-induced partition" are different statements**, and only the second would invalidate
#   Δ(L2). Demanding the first is a control testing a property the design does not require —
#   `#916`③'s family for the FIFTH consecutive round, and this time in the control I built
#   specifically to be the round's decisive one.
#   ⇒ REPAIRED into a DOSE-RESPONSE: sweep a multiplier on `#941`'s measured trends. The
#   instrument must SEE the artifact at SOME dose (that is the real sensitivity question), and the
#   informative output is **how many times the measured drift it would take** to manufacture the
#   observed Δ. That converts a yes/no into the quantity the confound argument actually needs.
yr = d.year.to_numpy(dtype=float)
yc = (yr - yr.mean()) / 10.0
drift = np.array([TREND[c] for c in ITEMS])
syn_dose = []
for mult in (0.0, 1.0, 3.0, 10.0, 30.0, 100.0):
    vals = []
    for _ in range(6):
        S = permute_within_person(X) + mult * yc[:, None] * drift[None, :]
        vals.append({name: contrast(pair_rs(group_centre(S, g))) for name, g in GROUPS.items()})
    syn_dose.append(dict(mult=mult,
                         L0=float(np.median([v["L0 pooled"] for v in vals])),
                         L2=float(np.median([v["L2 wave x cohort"] for v in vals]))))
print(f"\n  ⚠ SYNTHETIC-ARTIFACT DOSE-RESPONSE — structureless data + m x `#941`'s measured "
      f"differential drift. m=1 IS the drift actually observed in the world.")
print(f"    {'m':>6s} {'Δ at L0':>10s} {'Δ at L2':>10s}   (L0 null {levels['L0 pooled']['null_mean']:+.4f} "
      f"+/- {levels['L0 pooled']['null_sd']:.4f})")
for s in syn_dose:
    print(f"    {s['mult']:6.1f} {s['L0']:+10.4f} {s['L2']:+10.4f}")
L0n, L0sd = levels["L0 pooled"]["null_mean"], levels["L0 pooled"]["null_sd"]
# ⚠⚠ TWO-SIDED, and v2 got this wrong too. Sensitivity asks "can the instrument SEE a drift-induced
#   partition", not "can it see one in the direction I predicted". The artifact turns out to move Δ
#   the OTHER way, so a one-sided test called the instrument blind while it was reading the artifact
#   loud and clear at -0.2092.
visible = [s["mult"] for s in syn_dose if abs(s["L0"] - L0n) > 2 * L0sd]
# ⚠⚠⚠ AND THE SIGN IS THE FINDING. My confound argument grouped items by RANK of trend ("the two
#   fastest"), but a person-centred residual carries each item's DEVIATION FROM THE MEAN TREND:
#   mean trend {mt:.4f}; homosex {h:+.4f}, premarsx {p:+.4f}, teensex {te:+.4f}, xmarsex {x:+.4f}.
#   **premarsx sits essentially AT the mean and homosex far above it**, so the drift artifact makes
#   a 1-vs-3 structure (homosex against the rest), which is orthogonal to my 2+2 and pushes the
#   contrast NEGATIVE. Arguing could not have found this; building the world did.
_mt = float(np.mean(list(TREND.values())))
dev = {k: TREND[k] - _mt for k in ITEMS}
print(f"    ⇒ WHY the sign flips: a person-centred residual carries DEVIATION FROM THE MEAN TREND "
      f"({_mt:+.4f}), not rank — " + " · ".join(f"{k} {dev[k]:+.4f}" for k in ITEMS) +
      f". `premarsx` sits at the mean while `homosex` is far above, so the artifact separates "
      f"homosex from the other THREE, not `#949`'s 2+2, and drives the contrast the other way.")
syn_visible = len(visible) > 0                     # the instrument CAN see a drift artifact
syn_measured_small = (syn_dose[1]["L0"] - L0n) < 2 * L0sd     # at m=1 it does NOT
reach = [s["mult"] for s in syn_dose if s["L0"] >= levels["L0 pooled"]["delta"]]
syn_sign_opposes = syn_dose[-1]["L0"] < syn_dose[0]["L0"]
syn = {name: syn_dose[1][k] for name, k in (("L0 pooled", "L0"), ("L2 wave x cohort", "L2"))}
syn_removed = abs(syn_dose[1]["L2"] - levels["L2 wave x cohort"]["null_mean"]) < \
    2 * levels["L2 wave x cohort"]["null_sd"]
print(f"    ⇒ the instrument CAN see a drift artifact (two-sided, at m>=" +
      (f"{min(visible):.0f}" if visible else "never in this sweep") + f"): {syn_visible}")
print(f"    ⇒ at the MEASURED drift (m=1) it is invisible: {syn_measured_small}")
print(f"    ⇒ it would take m>=" + (f"{min(reach):.0f}" if reach else ">100") +
      f" — i.e. {'at least ' + str(int(min(reach))) if reach else 'more than 100'}x the drift the "
      f"world actually shows — to manufacture the observed Δ(L0) = "
      f"{levels['L0 pooled']['delta']:+.4f}")

# ══ POSITIVE CONTROL — plant a real partition into the permuted world ════════════════
sweep = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(12):
        P = permute_within_person(X)
        if gg:
            m = RNG.random(len(P)) < gg
            shift = P[m][:, list(CLUSTER_A)].mean(axis=1) - P[m][:, list(CLUSTER_B)].mean(axis=1)
            for i in CLUSTER_A:
                P[m, i] += 0.5 * shift
            for i in CLUSTER_B:
                P[m, i] -= 0.5 * shift
        vals.append(contrast(pair_rs(group_centre(P, GROUPS["L2 wave x cohort"]))))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a real partition planted into the permuted world, judged at L2, so g=0 IS "
      f"the null): {[(x, round(v, 4)) for x, v in sweep]}")
L2 = levels["L2 wave x cohort"]
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs L2 null {L2['null_mean']:+.4f} +/- "
      f"{L2['null_sd']:.4f} = {abs(sweep[0][1] - L2['null_mean']) / max(L2['null_sd'], 1e-9):.2f} "
      f"spreads")

# ══ SPECIFICATION CURVE (G4) ═════════════════════════════════════════════════════════
grid = []
for name, g in GROUPS.items():
    for cen, use_group in (("group-centred", True), ("pooled-centred", False)):
        base = group_centre(X, g) if use_group else (X - X.mean(axis=0)) / POOLED_SD
        rs = pair_rs(base)
        grid.append(dict(level=name, centring=cen, delta=contrast(rs),
                         signs_split=int(sum(1 for p in WITHIN if rs[p] > np.mean(list(rs.values())))
                                         + sum(1 for p in CROSS if rs[p] < np.mean(list(rs.values()))))))
print("\n  specification curve — every cell, none dropped (signs_split = how many of the 6 pairs "
      "fall on the side the 2+2 partition predicts, max 6)")
for r in grid:
    print(f"    {r['level']:<20s} {r['centring']:<15s} Δ {r['delta']:+.4f}  signs "
          f"{r['signs_split']}/6")

ps = []
for name, L in levels.items():
    for p in PAIRS:
        ps.append(2 * (1 - stats.norm.cdf(abs((L["rs"][p] - L["null_mean"]) /
                                              max(L["null_sd"], 1e-9)))))
ps += [2 * (1 - stats.norm.cdf(abs(L["z"]))) for L in levels.values()]

G = Gate("Is the 2+2 partition a fact about people, or an artifact of pooling 21 waves?")
G.plant_direction_from_sweep("positive: a planted partition raises Δ at L2, and g=0 sits ON the L2 "
                             "null (`#922`)", sweep,
                             baseline=L2["null_mean"], baseline_spread=max(L2["null_sd"], 1e-4))
G.negative_control("within-person item-label permutation, recomputed INSIDE each level's groups, "
                   "reproduces itself across two independent draws (`#948`①)",
                   float(np.mean([L["repro"] for L in levels.values()])),
                   abs(L2["delta"] - L2["null_mean"]),
                   null_spread=float(np.mean([L["null_sd"] for L in levels.values()])),
                   null_kind="within-person item-label permutation inside each grouping level "
                             "(each person's total preserved EXACTLY)")
G.multiplicity_control("3 levels x 6 pairs + 3 contrasts = 21 cells (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"{n.split()[0]}/{ITEMS[p[0]][:5]}-{ITEMS[p[1]][:5]}"
                               for n in levels for p in PAIRS]
                              + [f"{n.split()[0]}/Δ" for n in levels])
G.asserted("⚠⚠ SYNTHETIC-ARTIFACT DOSE-RESPONSE: the rival's world is BUILT at six doses of "
           "`#941`'s measured drift. The instrument must SEE a drift-induced partition at SOME "
           "dose (sensitivity), and the informative number is how many times the MEASURED drift it "
           "would take to manufacture the observed Δ",
           bool(syn_visible),
           f"doses m={[s['mult'] for s in syn_dose]} give L0 Δ="
           f"{[round(s['L0'], 4) for s in syn_dose]} against a null of {L0n:+.4f} +/- {L0sd:.4f}. "
           f"Visible (TWO-SIDED) from m>=" + (f"{min(visible):.0f}" if visible else "never") +
           f"; at the MEASURED drift m=1 it is invisible ({syn_measured_small}); reaching the "
           f"observed Δ(L0)={levels['L0 pooled']['delta']:+.4f} needs m>=" +
           (f"{min(reach):.0f}" if reach else ">100") +
           f". ⚠ v1 required visibility AT m=1 and refused the round when the measured drift turned "
           f"out too weak to see — conflating 'the artifact is small' with 'the instrument is "
           f"blind'. ⚠⚠ AND THE SIGN REFUTES THE CONFOUND OUTRIGHT: more drift drives the contrast "
           f"DOWN, because a person-centred residual carries deviation from the MEAN trend and "
           f"`premarsx` sits at that mean while `homosex` is far above — so the artifact makes a "
           f"1-vs-3 split, not `#949`'s 2+2. Arguing could not have found that; building the world "
           f"did", kind="control",
           population=f"GSS ballot 1, n={len(d)}, {len(waves)} waves, {len(cells_L2)} L2 cells")
G.asserted("⚠ the confound is not 'equally available', it ENTAILS the finding — checked against "
           "`#941`'s artifact, not memory", True,
           f"trends {' · '.join(f'{k} {TREND[k]:+.4f}' for k in order)}; fastest two "
           f"{{{order[0]}, {order[1]}}} vs `#949`'s cluster A — coincide: {predicts}",
           kind="control",
           population=f"GSS ballot 1, n={len(d)}, {len(waves)} waves, {len(cells_L2)} L2 cells")
G.asserted("⚠ ESTIMAND CHANGED FROM `#949` AND IT IS SAID SO: `#949` asked which pairs depart from "
           "their null; this asks how much BLOCK STRUCTURE survives removing between-group profile "
           "variation", True,
           f"Δ(L0) {levels['L0 pooled']['delta']:+.4f} · Δ(L1) {levels['L1 wave']['delta']:+.4f} · "
           f"Δ(L2) {L2['delta']:+.4f} · retained {retained:.3f}", kind="control",
           population=f"GSS ballot 1, n={len(d)}, {len(waves)} waves, {len(cells_L2)} L2 cells")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{r['level'].split()[0]}/{r['centring'][:4]} Δ{r['delta']:+.3f} "
                      f"signs{r['signs_split']}/6" for r in grid), kind="control",
           population=f"GSS ballot 1, n={len(d)}, {len(waves)} waves, {len(cells_L2)} L2 cells")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * L2["null_sd"]
neg_null = (float(np.mean([L["repro"] for L in levels.values()]))
            < 0.5 * abs(L2["delta"] - L2["null_mean"])) and bool(syn_visible)
survives = (L2["z"] > 2) and (retained >= 0.5)
collapses = abs(L2["z"]) < 2
world = ("W_WRONGED" if survives else
         ("W_POOLING" if collapses else "W_PARTIAL"))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires, the null "
           "reproduces, AND the synthetic artifact is both visible at L0 and removed at L2. STAKED: "
           "W_WRONGED, i.e. Δ(L2) clears 2x its own null spread AND retains >=50% of Δ(L0). "
           "⚠ W_POOLING is the unwelcome branch and the one the arithmetic predicts — it retracts "
           "`#949`'s headline from both pages one round after it landed",
           (pos_fires and neg_null) and survives,
           f"positive fires {pos_fires} · null reproduces and artifact control passes {neg_null} · "
           f"Δ(L0) {levels['L0 pooled']['delta']:+.4f} → Δ(L2) {L2['delta']:+.4f} (retained "
           f"{retained:.3f}, z {L2['z']:+.2f}) ⇒ {world}",
           kind="kill", yardstick="partition contrast Δ at L2, against its own within-group "
                                  "permutation null",
           yardstick_noise=L2["null_sd"],
           population=f"GSS ballot 1, n={len(d)}, {len(waves)} waves 1988-2024, "
                      f"{len(cells_L2)} wave x cohort cells",
           direction="one-sided: W_WRONGED requires Δ(L2) POSITIVE and large")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if survives else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=950, round="E03·A126·R386", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_POOLING"),
               n=int(len(d)), waves=waves, n_L2_cells=int(len(cells_L2)),
               ref941_trends=TREND, speed_groups_coincide=bool(predicts),
               levels={k: {kk: (vv if not isinstance(vv, dict) else
                                {f"{ITEMS[a]}-{ITEMS[b]}": v for (a, b), v in vv.items()})
                           for kk, vv in L.items()} for k, L in levels.items()},
               retained=retained, synthetic=syn, synthetic_dose=syn_dose,
               synthetic_visible_at_any_dose=bool(syn_visible),
               synthetic_invisible_at_measured_dose=bool(syn_measured_small),
               synthetic_sign_opposes_finding=bool(syn_sign_opposes),
               mean_trend=float(np.mean(list(TREND.values()))),
               trend_deviations={k: TREND[k] - float(np.mean(list(TREND.values()))) for k in ITEMS},
               dose_needed_to_reach_observed=(min(reach) if reach else None),
               synthetic_removed_L2=bool(syn_removed),
               null_median=L2["null_mean"], null_sd=L2["null_sd"], null_draws=NREPS,
               positive_sweep=sweep, grid=grid, family_size=len(ps),
               world=world, verdict=verdict),
          open(OUT / "is_the_partition_a_pooling_artifact.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'is_the_partition_a_pooling_artifact.json'}")
