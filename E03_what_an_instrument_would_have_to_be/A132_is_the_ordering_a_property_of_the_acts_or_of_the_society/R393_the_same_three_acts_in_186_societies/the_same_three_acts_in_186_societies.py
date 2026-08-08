#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A132·R393 — the same three acts, rated in 186 societies instead of asked of 15,056 Americans
=================================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#959`②. `#959` closed A131 by showing the 2+2 pairing carries nothing beyond the
                four item permission rates — so **the rates ARE the structure**: `premarsx` .5305 >
                `homosex` .3859 > `teensex` .0748 > `xmarsex` .0221. But four items give four
                numbers, and at least three readings fit them exactly: **harm to a wronged party**
                (`#949`'s reading, still on the page as untested) · **how far US liberalisation got**
                (`#941`: the two fastest-moving items are the two most permitted) · **consent and
                capacity**. GSS cannot separate them — four points, three curves.

⚠ HARD RULE 4, AND IT IS WHY THIS ROUND LEAVES THE INSTRUMENT. A fifth GSS round cannot separate
  readings that all fit the same four numbers. The separator has to come from a site where the
  *society* varies: if the ordering is about the ACTS, it should survive being measured by
  ethnographers in 186 mostly preindustrial societies; if it is about how far the United States got
  between 1988 and 2024, it should not.

⚠ HARD RULE 1 — MEASURED BEFORE ANY CLAIM, AND IT REWROTE THIS DESIGN.
    SCCS165 Premarital Sex Attitudes (Female)  n = 130/186   codes 1:30 2:28 3:22 4:11 5:4 6:35
    SCCS169 Extramarital Sex                   n = 109/186   codes 1:13 2:48 3:24 4:24
    SCCS176 Homosexuality                      n =  40/186   codes 1:9 2:4 3:6 4:4 5:17
    complete on all three                      n =  27
  **Twenty-seven societies.** The marginal rates rest on three DIFFERENT society sets, which is
  `#956`'s population mismatch one level up, so the identified estimand is the PAIRED one and the
  marginal one is a specification with its selection named.

⚠ HARD RULE 2 — THE INSTRUMENT, NAMED. All three variables are **Broude & Greene 1976, *Ethnology*
  15(4):409–429, "Cross-cultural codes on twenty sexual attitudes and practices"**, coded from HRAF
  ethnographies. **One team, one paper, one corpus.** Every number below is a claim about what that
  team read and recorded. There is **no second ethnographic coding of these acts** in this release.

⚠ G1 — THE ESTIMAND. Among the 27 SCCS societies coded on all three, **the share of societies in
  which each act sits at its variable's most permissive code**, and the three pairwise orderings.
  The decisive cell is **homosexuality vs extramarital sex**: in GSS the two are 0.3859 and 0.0221,
  a 17× gap in favour of same-sex relations.

Live Worlds    W_ACT_UNIVERSAL   · the SCCS ordering matches GSS — premarital most permitted,
                                    homosexuality next, extramarital least ⇒ the ordering is a
                                    property of the acts and survives an instrument change of
                                    maximal distance. `#949`'s wronged-party reading lives.
               W_SOCIETY_SPECIFIC· homosexuality sits at or BELOW extramarital across societies ⇒
                                    GSS's ordering is where the United States got to, not a fact
                                    about acts, and the wronged-party reading dies by evidence
                                    rather than by caution.
               W_RECORD          · ⚠ **the meta-separator, and it says my decomposition may be
                                    wrong.** Ethnographers recorded premarital norms for 130
                                    societies and homosexuality for 40. If the ordering tracks
                                    COVERAGE rather than either acts or societies, the object is
                                    neither — it is the ethnographic record, and both worlds above
                                    are asking about an artefact of what got written down.

Prediction     W_ACT_UNIVERSAL    -> homosexuality permitted in MORE of the 27 than extramarital,
Matrix                               McNemar discordance favouring it, in most of the cut grid.
               W_SOCIETY_SPECIFIC -> equal or fewer, in most of the cut grid.
               W_RECORD           -> the societies missing SCCS176 differ systematically on
                                     SCCS165/SCCS169 from those that have it.

Strongest      **THE BINARISATION CHOOSES THE ANSWER.** Three ordinal scales of different lengths
confound       (6, 4, 5 levels) with no common metric; a lenient cut on one act and a strict cut on
(written       another manufactures any ordering. ⇒ CONTROL, SAME ITERATION: every admissible cut
before)        for each variable is swept (3 x 2 x 3 = 18 cells) and the whole grid is published,
               with the headline at the least arbitrary choice — each variable's single most
               permissive code. ⚠ And **SCCS176 code 2 carries NO LABEL in the codebook**; it is
               swept both ways and never resolved.

Controls       POSITIVE: shift SCCS176 toward permissive by g scale levels; the permitted rate must
                 rise monotonically and the ordering must flip at a measurable dose; g=0 must sit on
                 the observed value. Reports the MDE in societies.
               NEGATIVE / N1: permute the three binary indicators WITHIN each society — preserves
                 how many acts that society permits, equalises the per-act rates, so the ordering
                 statistic has a KNOWN expected value of zero.
               PLACEBO: the whole pipeline on N1 data; the ordering must vanish.
               SHAM: all three pairwise orderings are reported, not only the decisive one.
               COVERAGE ARM: are societies missing SCCS176 different on the other two?
               REGION: the ordering recomputed leaving out each of Murdock's six world regions,
                 because 27 societies from six regions is where Galton's problem would show.
               MULTIPLICITY: 3 pairs x 18 cut cells.
               SEEDS: 3.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **only this one instrument** — Broude & Greene 1976 coded all three variables from HRAF;
    separating the society from the ethnographer would require a second, independent coding team;
  (2) ⚠ **no teenage-sex analogue from the same team** — SCCS827–830 are Barry & Schlegel 1984 on a
    different 1–10 scale, so the four-way GSS ordering CANNOT be reproduced; this is a THREE-way
    test and the fourth act is conceded in writing;
  (3) ⚠ **`SCCS176` code 2 has no label in `codes.csv`** — 4 societies carry a code whose meaning is
    not in the release; swept, not resolved;
  (4) ⚠ **no causal identification, no intervention** — 27 societies observed once;
  (5) ⚠ **within-society individual variation is invisible**, and the ethnographies' focal years
    differ by more than a century, so "the society" is one coder's summary of one period;
  (6) ⚠ **Galton's problem is mitigated by SCCS's design, not eliminated** — the region leave-one-out
    bounds it, it does not remove it;
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
SEEDS = (393, 1393, 2393)
ACTS = {"SCCS165": "premarital", "SCCS169": "extramarital", "SCCS176": "homosexuality"}
CORE = list(ACTS)
# GSS anchors, READ from `#959`'s artifact rather than typed (`#840`)
G = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                   "A131_why_the_contrast_reverses_among_sharp_discriminators" /
                   "R392_the_six_ways_to_split_four_acts" / "results" /
                   "the_six_ways_to_split_four_acts.json"))
GSS_RATE = dict(zip(["premarital", "teenage", "extramarital", "homosexuality"], G["p_top_full"]))

d = pd.read_csv(B / "data.csv", low_memory=False)
soc = pd.read_csv(B / "societies.csv")
W = d[d.var_id.isin(CORE)].pivot_table(index="soc_id", columns="var_id", values="code",
                                       aggfunc="first")
REG = [(1, 28, "Africa"), (29, 58, "Circum-Mediterranean"), (59, 87, "East Eurasia"),
       (88, 114, "Insular Pacific"), (115, 153, "North America"), (154, 186, "South America")]


def region_of(sid):
    n = int(str(sid).replace("SCCS", ""))
    for a, b, nm in REG:
        if a <= n <= b:
            return nm
    return "?"


print(f"HARD RULE 1 — SCCS, {d.soc_id.nunique()} societies in the release; per-variable coverage and "
      f"code distribution, printed before any claim:")
for vid, name in ACTS.items():
    vc = W[vid].value_counts().sort_index()
    print(f"  {vid} {name:<14s} n={int(W[vid].notna().sum()):3d}/186  "
          + " ".join(f"{int(k)}:{int(v)}" for k, v in vc.items()))
P = W[W[CORE].notna().all(axis=1)].copy()
P["region"] = [region_of(i) for i in P.index]
N = len(P)
print(f"  complete on all three (the identified population): n = {N}, regions "
      f"{P.region.value_counts().to_dict()}")
print("⚠ HARD RULE 2 — instrument: Broude & Greene 1976, Ethnology 15(4):409-429, all three "
      "variables coded by ONE team from the HRAF corpus. Only this one instrument exists here.")
print(f"⚠ GSS anchors read from `#959`'s artifact, not typed: "
      f"{ {k: round(v, 4) for k, v in GSS_RATE.items()} }")

# ── the cut grid. Each entry: permitted iff code <= cut. Pre-registered, all published ──
CUTS = {"SCCS165": [1, 2, 3],       # Expected | +Tolerated | +Mildly disapproved
        "SCCS169": [1, 2],          # both allowed | +husband allowed
        "SCCS176": [1, 2, 3]}       # Accepted/ignored | +the UNLABELLED code 2 | +Ridiculed
HEAD = {"SCCS165": 1, "SCCS169": 1, "SCCS176": 1}   # least arbitrary: each variable's TOP code


def binarise(df, cuts):
    return pd.DataFrame({v: (df[v] <= cuts[v]).astype(float) for v in CORE}, index=df.index)


def mcnemar(a, b):
    """exact paired test on discordant societies: b=1,a=0 vs a=1,b=0."""
    n01 = int(((a == 0) & (b == 1)).sum())
    n10 = int(((a == 1) & (b == 0)).sum())
    if n01 + n10 == 0:
        return n10, n01, 1.0
    return n10, n01, float(stats.binomtest(n10, n01 + n10, 0.5).pvalue)


def ordering_stat(Bn):
    """the decisive quantity: rate(homosexuality) - rate(extramarital). GSS says +0.3638."""
    return float(Bn["SCCS176"].mean() - Bn["SCCS169"].mean())


BH = binarise(P, HEAD)
rates = {ACTS[v]: float(BH[v].mean()) for v in CORE}
D_OBS = ordering_stat(BH)
D_GSS = GSS_RATE["homosexuality"] - GSS_RATE["extramarital"]
print(f"\n  HEADLINE, paired n={N}, each variable at its MOST PERMISSIVE code:")
for v in CORE:
    k = int(BH[v].sum())
    lo, hi = stats.binomtest(k, N).proportion_ci(0.95)
    print(f"    {ACTS[v]:<14s} permitted in {k:2d}/{N} = {k/N:.4f}  [{lo:.4f}, {hi:.4f}]   "
          f"GSS {GSS_RATE[ACTS[v]]:.4f}")
print(f"    DECISIVE: rate(homosexuality) - rate(extramarital) = {D_OBS:+.4f}   "
      f"GSS says {D_GSS:+.4f}")
pairs = {}
for a, b in (("SCCS176", "SCCS169"), ("SCCS165", "SCCS169"), ("SCCS165", "SCCS176")):
    n10, n01, p = mcnemar(BH[b], BH[a])
    pairs[f"{ACTS[a]}>{ACTS[b]}"] = dict(n_a_only=n10, n_b_only=n01, p=p)
    print(f"    McNemar {ACTS[a]:>14s} vs {ACTS[b]:<14s}: {n10} societies permit only the first, "
          f"{n01} only the second, exact p = {p:.4f}")

# ══ THE CUT GRID (G4) — the strongest confound, controlled ═══════════════════════════
grid, ps = [], []
for c1 in CUTS["SCCS165"]:
    for c2 in CUTS["SCCS169"]:
        for c3 in CUTS["SCCS176"]:
            cuts = {"SCCS165": c1, "SCCS169": c2, "SCCS176": c3}
            Bn = binarise(P, cuts)
            n10, n01, p = mcnemar(Bn["SCCS169"], Bn["SCCS176"])
            grid.append(dict(cut165=c1, cut169=c2, cut176=c3,
                             r_pre=float(Bn["SCCS165"].mean()), r_ext=float(Bn["SCCS169"].mean()),
                             r_hom=float(Bn["SCCS176"].mean()), delta=ordering_stat(Bn),
                             n_hom_only=n10, n_ext_only=n01, p=p))
            ps.append(p)
same_sign = sum(1 for g in grid if np.sign(g["delta"]) == np.sign(D_GSS))
gss_dir = sum(1 for g in grid if g["delta"] > 0)
print(f"\n  CUT GRID — {len(grid)} cells, all published. Cells where homosexuality is MORE permitted "
      f"than extramarital (GSS's direction): {gss_dir}/{len(grid)}")
for g in grid:
    print(f"    165<={g['cut165']} 169<={g['cut169']} 176<={g['cut176']} | pre {g['r_pre']:.3f} "
          f"ext {g['r_ext']:.3f} hom {g['r_hom']:.3f} | Δ(hom−ext) {g['delta']:+.4f} | "
          f"discordant {g['n_hom_only']}v{g['n_ext_only']} p={g['p']:.4f}")

# ══ N1 — a null with a KNOWN value: permute the three indicators WITHIN society ══════
n1 = []
for s in SEEDS:
    r = np.random.default_rng(s)
    vals = []
    for _ in range(4000):
        M = BH.to_numpy().copy()
        for i in range(len(M)):
            M[i] = r.permutation(M[i])
        vals.append(M[:, CORE.index("SCCS176")].mean() - M[:, CORE.index("SCCS169")].mean())
    n1.append(vals)
n1 = np.array(n1)
n1_mu, n1_sd = float(n1.mean()), float(n1.std(ddof=1))
z_obs = (D_OBS - n1_mu) / max(n1_sd, 1e-9)
p_perm = float((np.abs(n1.ravel() - n1_mu) >= abs(D_OBS - n1_mu)).mean())
print(f"\n  N1 within-society label permutation (known expected value 0): {n1_mu:+.5f} ± {n1_sd:.4f} "
      f"· observed {D_OBS:+.4f} · z {z_obs:+.2f} · two-sided empirical p {p_perm:.4f}")

# ══ PLACEBO — the whole pipeline on N1 data, the ordering must vanish ════════════════
plac = []
for s in SEEDS:
    r = np.random.default_rng(s + 500)
    M = BH.to_numpy().copy()
    for i in range(len(M)):
        M[i] = r.permutation(M[i])
    plac.append(M[:, CORE.index("SCCS176")].mean() - M[:, CORE.index("SCCS169")].mean())
plac_mu = float(np.mean(plac))
print(f"  PLACEBO (pipeline on permuted data): Δ {plac_mu:+.4f} against a null spread of {n1_sd:.4f}")

# ══ POSITIVE CONTROL — shift homosexuality toward permissive by g scale levels ═══════
sweep, mde_g = [], None
for g in (0, 1, 2, 3, 4):
    Q = P.copy()
    Q["SCCS176"] = np.maximum(Q["SCCS176"] - g, 1)
    Bn = binarise(Q, HEAD)
    v = ordering_stat(Bn)
    sweep.append([float(g), float(v)])
    if mde_g is None and v - sweep[0][1] > 2 * n1_sd:
        mde_g = g
print(f"  positive sweep (shift SCCS176 down by g levels; g=0 IS the observed): "
      f"{[(g, round(v,4)) for g, v in sweep]} · g=0 at {sweep[0][1]:+.4f} vs observed "
      f"{D_OBS:+.4f} · first g clearing 2 null spreads: {mde_g}")

# ══ MDE of the paired design, measured by simulation ═════════════════════════════════
mde_rows = []
for true_gap in (0.10, 0.20, 0.30, 0.40):
    hits = 0
    for j in range(2000):
        r = np.random.default_rng(9000 + j)
        base = float(BH["SCCS169"].mean())
        a = r.random(N) < base
        b = r.random(N) < min(base + true_gap, 1.0)
        _, _, p = mcnemar(pd.Series(a.astype(float)), pd.Series(b.astype(float)))
        hits += p < 0.05
    mde_rows.append([true_gap, hits / 2000])
print(f"  MDE of the paired design at n={N} (power to detect a true rate gap, McNemar exact): "
      + " · ".join(f"gap {g:.2f} -> power {p:.2f}" for g, p in mde_rows))

# ⚠⚠ POWER AT THE OBSERVED EFFECT, and this is a repair of my own control. The check above asks
#    "can this design detect ANYTHING" and passes at gap 0.40. The question that decides whether a
#    non-significant result is evidence is "can it detect THIS" -- the observed gap. A power check
#    that passes because some larger effect would be visible is a check that cannot fail.
hits = 0
for j in range(4000):
    r = np.random.default_rng(40000 + j)
    base = float(BH["SCCS169"].mean())
    a = r.random(N) < base
    b = r.random(N) < min(base + abs(D_OBS), 1.0)
    _, _, pv = mcnemar(pd.Series(a.astype(float)), pd.Series(b.astype(float)))
    hits += pv < 0.05
power_at_obs = hits / 4000
print(f"  ⚠⚠ POWER AT THE OBSERVED GAP |Δ| = {abs(D_OBS):.4f}: **{power_at_obs:.3f}**. A "
      f"non-significant result from a design with this power is SILENCE, not evidence for the null")

# ── the gendered asymmetry, a descriptive count that needs no null ───────────────────
e = W["SCCS169"].dropna()
asym = int(((e == 2) | (e == 3)).sum()); n_e = int(len(e))
lo_a, hi_a = stats.binomtest(asym, n_e).proportion_ci(0.95)
print(f"\n  ⚠ THE CODEBOOK ITSELF: of {n_e} societies with a coded extramarital norm, **{asym} "
      f"({asym/n_e:.1%}) [{lo_a:.3f}, {hi_a:.3f}] apply it ASYMMETRICALLY BY SEX** — code 2 "
      f"'husband only' n={int((e==2).sum())}, code 3 'both forbidden, women punished' "
      f"n={int((e==3).sum())}. Symmetric rules: code 1 both allowed n={int((e==1).sum())}, code 4 "
      f"both condemned n={int((e==4).sum())}. And SCCS165 is titled 'Premarital Sex Attitudes - "
      f"FEMALE'. GSS's `xmarsex` asks about 'a married person' and cannot express any of this.")

# ══ COVERAGE ARM — W_RECORD. Do societies missing SCCS176 differ on the other two? ═══
have = W["SCCS176"].notna()
cov = {}
for v in ("SCCS165", "SCCS169"):
    a = W.loc[have & W[v].notna(), v]
    b = W.loc[(~have) & W[v].notna(), v]
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    cov[v] = dict(n_with=len(a), n_without=len(b), mean_with=float(a.mean()),
                  mean_without=float(b.mean()), p=float(u.pvalue))
    print(f"\n  COVERAGE ARM {v}: societies WITH a homosexuality code n={len(a)} mean {a.mean():.3f} "
          f"· WITHOUT n={len(b)} mean {b.mean():.3f} · Mann-Whitney p {u.pvalue:.4f}")
print(f"  ⚠ coverage itself: premarital {int(W['SCCS165'].notna().sum())}/186 · extramarital "
      f"{int(W['SCCS169'].notna().sum())}/186 · homosexuality "
      f"{int(W['SCCS176'].notna().sum())}/186 — the record is 3.25x thinner on homosexuality")

# ══ REGION LEAVE-ONE-OUT — Galton's problem, bounded not removed ════════════════════
loo = []
for rg in sorted(P.region.unique()):
    sub = P[P.region != rg]
    Bn = binarise(sub, HEAD)
    loo.append(dict(dropped=rg, n=len(sub), delta=ordering_stat(Bn)))
print(f"\n  REGION LEAVE-ONE-OUT: "
      + " · ".join(f"−{r['dropped'][:12]} (n={r['n']}) Δ {r['delta']:+.4f}" for r in loo))

# ══ MARGINAL SPECIFICATION — each variable on its OWN society set (selection named) ══
marg = {ACTS[v]: float((W[v].dropna() <= HEAD[v]).mean()) for v in CORE}
marg_n = {ACTS[v]: int(W[v].notna().sum()) for v in CORE}
print(f"  MARGINAL spec (different society sets, `#956`'s mismatch one level up): "
      + " · ".join(f"{k} {marg[k]:.4f} (n={marg_n[k]})" for k in marg)
      + f" · Δ(hom−ext) {marg['homosexuality'] - marg['extramarital']:+.4f}")

# ══ VERDICT ═════════════════════════════════════════════════════════════════════════
pos_fires = mde_g is not None and sweep[-1][1] > sweep[0][1] + 2 * n1_sd
pos_g0_ok = abs(sweep[0][1] - D_OBS) < 1e-9
plac_ok = abs(plac_mu - n1_mu) < 2 * n1_sd
n1_ok = abs(n1_mu) < 0.02
powered = power_at_obs >= 0.50   # repaired: power AT THE OBSERVED effect, not at any effect
coverage_biased = any(c["p"] < 0.05 for c in cov.values())
controls_ok = pos_fires and pos_g0_ok and plac_ok and n1_ok

if not controls_ok:
    world = "W_UNREADABLE"
elif not powered:
    world = "W_UNDERPOWERED"
elif coverage_biased:
    world = "W_RECORD"
elif D_OBS > 0 and p_perm < 0.05:
    world = "W_ACT_UNIVERSAL"
else:
    world = "W_SOCIETY_SPECIFIC"

Gt = Gate("Is the permissiveness ordering across sexual acts a property of the ACTS or of the "
          "SOCIETY? — the same three acts, rated in 186 societies instead of asked of Americans")
Gt.plant_direction_from_sweep(
    "positive: shifting the homosexuality codes toward permissive by g scale levels must raise "
    "rate(homosexuality) − rate(extramarital); g=0 IS the observed value by construction, computed "
    "in the round and not asserted (`#959`③)", sweep, baseline=D_OBS,
    baseline_spread=max(n1_sd, 1e-4))
Gt.negative_control(
    "N1, a null with a KNOWN value of zero: permuting the three permitted-indicators WITHIN each "
    "society preserves how many acts that society permits and equalises the per-act rates",
    abs(n1_mu), abs(D_OBS - n1_mu), null_spread=n1_sd,
    null_kind="within-society act-label permutation null, whose expected value is zero in closed "
              "form because permuting labels within a row cannot change the column means in "
              "expectation, plus an exact McNemar binomial on discordant societies")
Gt.multiplicity_control("3 pairs x 18 cut cells (`#936`②)", ps, 0.05,
                        labels=[f"165<={g['cut165']}/169<={g['cut169']}/176<={g['cut176']}"
                                for g in grid])
Gt.asserted(
    "⚠ HARD RULE 1: coverage and code distribution printed before any claim, and it rewrote the "
    "design — the identified population is the PAIRED one, because the three marginal rates rest on "
    "three different society sets", True,
    f"SCCS165 n={int(W['SCCS165'].notna().sum())} · SCCS169 n={int(W['SCCS169'].notna().sum())} · "
    f"SCCS176 n={int(W['SCCS176'].notna().sum())} · complete on all three n={N}",
    kind="control", population=f"SCCS, 186 societies, paired subset n={N}")
Gt.asserted(
    "⚠ HARD RULE 2: the instrument is Broude & Greene 1976 (Ethnology 15(4):409-429), one team "
    "coding all three variables from the HRAF corpus. Every number here is a claim about what that "
    "team read and recorded, and there is no second ethnographic coding — **only this one "
    "instrument** carries these acts cross-culturally", True,
    "SCCS165 / SCCS169 / SCCS176 all carry source=broude1976cross; the adolescent-sex variables "
    "SCCS827-830 are Barry & Schlegel 1984 on a different 1-10 scale and are NOT comparable, so the "
    "fourth GSS act is conceded rather than approximated", kind="control",
    population=f"SCCS paired subset n={N}")
Gt.asserted(
    "⚠⚠ THE STRONGEST CONFOUND, controlled in the same iteration: three ordinal scales of different "
    "lengths (6, 4, 5) with no common metric, so the binarisation can manufacture any ordering. "
    "Every admissible cut is swept and the whole grid published; the headline is the least arbitrary "
    "choice, each variable's single most permissive code",
    True, f"{len(grid)} cells · cells in GSS's direction (homosexuality more permitted than "
          f"extramarital) {gss_dir}/{len(grid)} · Δ range "
          f"[{min(g['delta'] for g in grid):+.4f}, {max(g['delta'] for g in grid):+.4f}]",
    kind="control", population=f"SCCS paired subset n={N}")
Gt.asserted(
    "⚠ THE META-SEPARATOR, and it is a live world rather than a caveat: the ethnographic record is "
    "3.25x thinner on homosexuality than on premarital norms, so the object may be neither the acts "
    "nor the societies but WHAT GOT WRITTEN DOWN. Tested: do societies missing a homosexuality code "
    "differ on the other two?", not coverage_biased,
    " · ".join(f"{v}: with n={c['n_with']} mean {c['mean_with']:.3f} vs without n={c['n_without']} "
               f"mean {c['mean_without']:.3f}, p={c['p']:.4f}" for v, c in cov.items()),
    kind="control", population="SCCS, all 186 societies")
Gt.asserted(
    "⚠⚠ POWER AT THE OBSERVED EFFECT, and this row is a REPAIR OF MY OWN CONTROL. The first version "
    "asked whether the design detects ANY gap and passed at 0.40 — a check that cannot fail. What "
    "decides whether a non-significant result is evidence is the power at the gap actually observed",
    powered,
    " · ".join(f"true gap {g:.2f} -> power {p:.2f}" for g, p in mde_rows)
    + f" · **power at the OBSERVED gap |Δ|={abs(D_OBS):.4f} is {power_at_obs:.3f}**", kind="control",
    population=f"SCCS paired subset n={N}")
Gt.asserted(
    "⚠ GALTON'S PROBLEM bounded, not removed: the decisive quantity recomputed leaving out each of "
    "Murdock's six world regions", True,
    " · ".join(f"−{r['dropped']} (n={r['n']}) Δ {r['delta']:+.4f}" for r in loo),
    kind="control", population=f"SCCS paired subset n={N}")
Gt.asserted(
    "⚠ MARGINAL SPECIFICATION reported with its selection named: each variable on its own society "
    "set gives more power and a different population, which is `#956`'s mismatch one level up", True,
    " · ".join(f"{k} {marg[k]:.4f} (n={marg_n[k]})" for k in marg)
    + f" · Δ(hom−ext) {marg['homosexuality'] - marg['extramarital']:+.4f} vs paired {D_OBS:+.4f}",
    kind="control", population="SCCS, per-variable society sets")
Gt.asserted("the whole cut grid is published, disagreeing cells included", True,
            " · ".join(f"[{g['cut165']},{g['cut169']},{g['cut176']}] Δ{g['delta']:+.3f}"
                       for g in grid), kind="control", population=f"SCCS paired subset n={N}")

Gt.asserted(
    "KILL: pre-registered CONDITIONAL — evaluated ONLY if the plant fires with g=0 on the observed "
    "value, N1 returns its known zero, and the placebo vanishes. STAKED: W_ACT_UNIVERSAL, i.e. "
    "homosexuality is permitted in MORE of the 27 societies than extramarital sex, in GSS's "
    "direction, clearing the within-society permutation null. ⚠ W_SOCIETY_SPECIFIC means GSS's "
    "ordering is where the United States got to and `#949`'s wronged-party reading dies",
    controls_ok and powered and (not coverage_biased) and D_OBS > 0 and p_perm < 0.05,
    f"plant fires {pos_fires} (first detectable shift g={mde_g}) · POWER AT THE OBSERVED GAP "
    f"{power_at_obs:.3f} · g=0 on observed {pos_g0_ok} · N1 "
    f"at zero {n1_ok} ({n1_mu:+.5f}) · placebo {plac_ok} · coverage unbiased "
    f"{not coverage_biased} · Δ(hom−ext) {D_OBS:+.4f} vs GSS {D_GSS:+.4f}, permutation p "
    f"{p_perm:.4f}, grid cells in GSS's direction {gss_dir}/{len(grid)} ⇒ {world}",
    kind="kill", yardstick="rate(homosexuality permitted) − rate(extramarital permitted) across "
                           "societies",
    yardstick_noise=n1_sd,
    population=f"SCCS societies coded by Broude & Greene 1976 on all three acts, n={N} of 186, "
               f"six world regions",
    direction="one-sided: W_ACT_UNIVERSAL requires the same POSITIVE sign as GSS")

print(Gt)
# ⚠ `#916`③(c): the verdict STRING must be computed from the same conditions the gate reads, never
#   typed alongside them. The first version of this line omitted `powered` and printed OVERTURNED
#   while the gate printed UNVERIFIED -- a false refutation, which is as permanent as a false
#   acquittal because nobody re-examines a claim its own author withdrew.
verdict = (("UNVERIFIED" if not (controls_ok and powered) else
            ("CONFIRMED" if (not coverage_biased and D_OBS > 0 and p_perm < 0.05)
             else "OVERTURNED"))
           + f" · world {world} · paired n={N} · rates pre {rates['premarital']:.4f} hom "
             f"{rates['homosexuality']:.4f} ext {rates['extramarital']:.4f} · Δ(hom−ext) "
             f"{D_OBS:+.4f} against GSS {D_GSS:+.4f}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=960, round="E03·A132·R393", gate_verdict=str(Gt).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in Gt.rows],
               claims_null=(world in ("W_SOCIETY_SPECIFIC", "W_RECORD")),
               instrument="Broude & Greene 1976, Ethnology 15(4):409-429, HRAF coding",
               coverage={v: int(W[v].notna().sum()) for v in CORE}, n_paired=int(N),
               regions=P.region.value_counts().to_dict(),
               rates_paired=rates, rates_marginal=marg, marginal_n=marg_n,
               gss_rate=GSS_RATE, delta_obs=D_OBS, delta_gss=D_GSS,
               mcnemar=pairs, grid=grid, cells_in_gss_direction=int(gss_dir),
               null_median=n1_mu, null_sd=n1_sd, null_draws=int(n1.size),
               perm_p=p_perm, placebo=plac_mu, positive_sweep=sweep, first_detectable_shift=mde_g,
               power=mde_rows, power_at_observed=power_at_obs,
               extramarital_asymmetric=asym, extramarital_coded=n_e,
               coverage_arm=cov, region_loo=loo,
               family_size=len(ps), seeds=list(SEEDS), world=world, verdict=verdict),
          open(OUT / "the_same_three_acts_in_186_societies.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_same_three_acts_in_186_societies.json'}")
