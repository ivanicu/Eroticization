#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A126·R385 — is the two-layer structure a fact about out-groups, or about moral attitudes?
==============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#945` found that tolerance has **two layers**: at the raw level every pair of
                targets is positively correlated (a general disposition — tolerant people are
                tolerant of nearly everyone), and **net of a person's own level, 12 of 30 pairs sit
                BELOW their ipsative null** — a specific exchange rate, sharpest between homosexuals
                and racists. **I do not know whether that is a fact about OUT-GROUPS or a fact about
                how people hold MORAL ATTITUDES at all.** Nothing in `#945` distinguishes them, and
                the distinction changes what the finding means.

⚠ BASIN, and    A123, A124 and A125 are three consecutive arcs on the **same** Stouffer battery.
why THIS step   That is a basin of OBJECT even though the last three verdicts were self-refutations,
                and `#111c` plus HARD RULE 4 both say the move is a different instrument. This round
                takes `#945`'s exact design to **the four sexual norms** — the project's actual
                subject, a different battery, the same respondents.
                **The unwelcome outcome is W_OUTGROUP_ONLY**: it would narrow a claim I made two
                rounds ago from "how people hold moral attitudes" to "how people rank out-groups".

Live Worlds    W_GENERAL_LAW · the sexual norms show the same two layers — positive raw pairs, and
                                ≥1 pair below its ipsative null ⇒ **the two-layer structure is how
                                moral attitudes are held**, not a property of out-groups, and
                                `#945` generalises.
               W_OUTGROUP_ONLY· raw pairs positive, but NO pair below its null ⇒ sexual norms are
                                held as a single graded permissiveness with no trades. `#945` is
                                about out-groups specifically. ⚠ **The unwelcome one.**
               W_INVERTED    · sexual norms show trades where tolerance did not, or vice versa in
                                a way no general story covers ⇒ the "two layers" framing is the
                                wrong decomposition for both. (the meta-separator)

⚠⚠ THE ARITHMETIC, DERIVED BEFORE THE DATA (a quantity the algebra forces is a DERIVATION):
   person-centred residuals over k items are ipsative, so their **mean pairwise correlation is
   forced to −1/(k−1)**. Here **k=4 ⇒ −1/3 = −0.3333**, a far harsher constraint than `#945`'s
   k=15 ⇒ −0.0714. **With only four items, most pairs MUST be strongly negative**, and comparing
   any of them to zero — or to `#945`'s floor — would be reporting the algebra. Every cell is
   judged against **its own** within-person label-permutation null.

⚠ AND A SCALE   `#945`'s items were BINARY with similar marginals. These four are **4-point with very
PROBLEM THAT    different means** — `premarsx` 2.982, `teensex` 1.605, `xmarsex` 1.339, `homosex`
`#945` DID NOT  2.336 — so subtracting a person's raw mean would conflate *"this person is permissive"*
HAVE           with *"this item is permissive"*, and the pair structure would be driven by item
                difficulty rather than by anybody's trade-off. ⇒ **each item is z-scored across
                people before person-centring**, and the permutation null shuffles the z-scores, so
                a label swap is meaningful. The raw-centred arm is reported beside it as the
                specification that shows how much this choice matters.

Estimand       Among the 14,847 respondents answering all four (GSS ballot 1, 21 waves 1988–2024):
(G1)           z-score each item across people, subtract each person's mean across the four, then
               for each of the C(4,2)=6 pairs the Pearson correlation of residuals across people.
               **The quantity is `r_pair − null_pair`**, the departure from what the ipsative
               constraint plus these marginals alone produce.

Prediction     W_GENERAL_LAW  -> ≥1 of 6 pairs below its null, surviving BH over all 6.
Matrix         W_OUTGROUP_ONLY-> none below; some above.
               W_INVERTED     -> a pattern neither the general-factor nor the trade-off story covers.

Strongest      **THE SAME OPERATION THAT CREATES THE CONSTRAINT REMOVES THE GENERAL FACTOR**
confound       (inherited from `#945` unchanged): person-centring removes exactly the general
(written       permissiveness that produces the positive raw correlations, so the residual picture
before)        is guaranteed to look more zero-sum than the raw one. ⇒ CONTROL, same iteration: the
               raw pair correlations are reported beside the residual ones, and no claim is made
               from residuals alone.

⚠ GAUGE TEST   A pair correlation is symmetric in the two items, so this design says *which pairs*
(zero compute) trade and **cannot say which item is being spent for which**. With only four items the
               pairs are nameable, unlike `#945`'s aggregate, so the loss is smaller here — but it
               is registered rather than discovered.

Controls       NEGATIVE: within-person item-label permutation — each person's total across the four
                 z-scores preserved EXACTLY, so the ipsative constraint is held fixed and only
                 *which item got which value* dies.
               POSITIVE: plant a trade-off between two named items INSIDE the permuted world and
                 sweep; `g=0` sits on the null by construction (`#922`, `#937`⑤).
               DERIVATION CHECK: the observed residuals' mean over all 6 pairs must land on −1/3.
               ⚠ NULL-WORLD CHECK (`#948`①, the live defect): a second independent permutation set,
                 to show the null reproduces itself — because three consecutive rounds have now had
                 a negative control that demanded a property the design does not require.
               MULTIPLICITY: the family is all **6 pairs × 2 centring arms = 12 cells**.
               SPEC CURVE (G4): {z-scored, raw-centred} × {all, cohorts in ≥3 waves}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **observe a person trading over time** — repeated cross-section;
  (2) ⚠ **distinguish a trade-off from a shared cause pushing two items apart** — inherited from
    `#945` unchanged, and a residual correlation is agnostic between them;
  (3) ⚠ **separate age from period from cohort** — `#939`/`#943`'s wall, inherited;
  (4) ⚠ **no second instrument for THIS construct** — the four norms are GSS's; **only this one
    instrument** carries them, so the cross-instrument move is unavailable *within* the round. The
    round IS itself the cross-battery test `#945` lacked, which is a different thing and is not
    claimed as the same;
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
RNG = np.random.default_rng(385)
ITEMS = ["premarsx", "teensex", "xmarsex", "homosex"]
K = len(ITEMS)
FLOOR = -1.0 / (K - 1)
PAIRS = list(combinations(range(K), 2))

# ⚠ `#945`'s floor is READ FROM ITS ARTIFACT, never typed. v1 wrote `-0.0714` into the derivation
#   print and `no_transcribed_numbers` (`#840`) blocked the commit — right, and for the reason that
#   matters here: **this round's whole argument is that the two floors are not comparable**, so the
#   number I am contrasting against must come from the round I am contrasting with.
REF945 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                        "A125_is_zero_sum_a_within_person_fact" /
                        "R383_the_trade_off_against_its_ipsative_floor" / "results" /
                        "the_trade_off_against_its_ipsative_floor.json"))
FLOOR_945, K_945 = REF945["derived_floor"], REF945["k"]

print(f"⚠ DERIVATION FIRST, before the data is read: person-centred residuals over k={K} items are "
      f"IPSATIVE, so their mean pairwise correlation is FORCED to -1/(k-1) = {FLOOR:+.4f}. That is "
      f"{abs(FLOOR / FLOOR_945):.1f}x harsher than `#945`'s k={K_945} floor of {FLOOR_945:+.4f} "
      f"(READ from its artifact, not typed), so **most pairs here MUST be strongly "
      f"negative** and comparing them to zero — or to `#945`'s floor — would report the algebra.")

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))      # 1 always wrong .. 4 not wrong at all (`#939`)
d = d[(d.ballot == 1)].dropna(subset=ITEMS + ["cohort"]).copy()
waves = sorted(int(y) for y in d.year.unique())
print(f"HARD RULE 1 — all four jointly, ballot 1, cohort present: n={len(d)} · {len(waves)} waves "
      f"{waves[0]}-{waves[-1]}")
print("  item means (1=always wrong .. 4=not wrong at all): "
      + " · ".join(f"{c} {d[c].mean():.3f} (sd {d[c].std():.3f})" for c in ITEMS))
print("⚠ HARD RULE 2 — the instrument: GSS ballot 1, four sexual-norm items asked in one block, one "
      "questionnaire. Same respondents as `#941`, a DIFFERENT battery from `#945`'s Stouffer items.")

X_raw = d[ITEMS].to_numpy(dtype=float)
Z = (X_raw - X_raw.mean(axis=0)) / X_raw.std(axis=0)     # z per item, across people
print(f"  ⚠ z-scored per item BEFORE person-centring, because the item means span "
      f"{X_raw.mean(axis=0).min():.3f}-{X_raw.mean(axis=0).max():.3f}; without it the pair "
      f"structure would be item difficulty, not anybody's trade-off")


def ipsative(mat):
    return mat - mat.mean(axis=1, keepdims=True)


def pair_rs(mat):
    R = ipsative(mat)
    return {(i, j): float(np.corrcoef(R[:, i], R[:, j])[0, 1]) for i, j in PAIRS}


def permute_within_person(mat):
    order = np.argsort(RNG.random(mat.shape), axis=1)
    return np.take_along_axis(mat, order, axis=1)


obs_z = pair_rs(Z)
obs_rawc = pair_rs(X_raw)
all_pairs_mean = float(np.mean(list(obs_z.values())))
print(f"  ipsative check: mean over all {len(PAIRS)} pairs of the z-arm residuals = "
      f"{all_pairs_mean:+.4f} vs derived {FLOOR:+.4f}")

# ══ NEGATIVE CONTROL — within-person item-label permutation, total preserved exactly ══
NREPS = 200
draws_A = {k: [] for k in PAIRS}
draws_B = {k: [] for k in PAIRS}
for store in (draws_A, draws_B):
    for _ in range(NREPS):
        r = pair_rs(permute_within_person(Z))
        for k in store:
            store[k].append(r[k])
null_mean = {k: float(np.mean(v)) for k, v in draws_A.items()}
null_sd = {k: float(np.std(v)) for k, v in draws_A.items()}
null_B = {k: float(np.mean(v)) for k, v in draws_B.items()}
reproducibility = float(np.mean([abs(null_mean[k] - null_B[k]) for k in PAIRS]))
print(f"\n  null (within-person ITEM-LABEL permutation — each person's total across the four "
      f"z-scores preserved EXACTLY; kind of null: within-person item-label permutation): grand mean "
      f"{np.mean(list(null_mean.values())):+.4f}")
print(f"  ⚠ NULL-WORLD CHECK (`#948`①): two independent permutation sets reproduce to "
      f"{reproducibility:.5f}")

rows = []
for (i, j) in PAIRS:
    r = obs_z[(i, j)]
    z = (r - null_mean[(i, j)]) / max(null_sd[(i, j)], 1e-9)
    rows.append(dict(a=ITEMS[i], b=ITEMS[j], r_resid=r, r_raw_pearson=float(
        np.corrcoef(X_raw[:, i], X_raw[:, j])[0, 1]), r_rawcentred=obs_rawc[(i, j)],
        null=null_mean[(i, j)], null_sd=null_sd[(i, j)], z=z,
        p=2 * (1 - stats.norm.cdf(abs(z)))))
rows.sort(key=lambda x: x["z"])
mean_abs_dep = float(np.mean([abs(r["r_resid"] - r["null"]) for r in rows]))
print(f"\n  all {len(PAIRS)} pairs, sorted by departure from their own null "
      f"(the null IS the ipsative floor for these marginals)")
print(f"  {'pair':<22s} {'r_resid':>8s} {'r_raw':>7s} {'null':>8s} {'z':>8s}")
for r in rows:
    print(f"  {r['a']+'/'+r['b']:<22s} {r['r_resid']:+8.4f} {r['r_raw_pearson']:+7.4f} "
          f"{r['null']:+8.4f} {r['z']:+8.2f}")

# ══ POSITIVE CONTROL — plant a trade-off INSIDE the permuted world ═══════════════════
PA, PB = 0, 3                                        # premarsx x homosex
sweep = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(20):
        P = permute_within_person(Z)
        if gg:
            m = RNG.random(len(P)) < gg
            P[m, PB] = -P[m, PA]                     # force an anti-correlation
        vals.append(pair_rs(P)[(PA, PB)])
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a trade-off planted between `{ITEMS[PA]}` and `{ITEMS[PB]}` inside the "
      f"permuted world, so g=0 IS the null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs that cell's null "
      f"{null_mean[(PA, PB)]:+.4f} +/- {null_sd[(PA, PB)]:.4f} = "
      f"{abs(sweep[0][1] - null_mean[(PA, PB)]) / max(null_sd[(PA, PB)], 1e-9):.2f} spreads")

# ══ SPECIFICATION CURVE (G4) ═════════════════════════════════════════════════════════
c3 = d.groupby("cohort").year.nunique()
keep = set(c3[c3 >= 3].index)
grid = []
for tag, sub in (("all respondents", d), ("cohorts in >=3 waves", d[d.cohort.isin(keep)])):
    Xs = sub[ITEMS].to_numpy(dtype=float)
    Zs = (Xs - Xs.mean(axis=0)) / Xs.std(axis=0)
    for arm, mat in (("z-scored", Zs), ("raw-centred", Xs)):
        rr = pair_rs(mat)
        grid.append(dict(spec=tag, arm=arm, n=int(len(sub)),
                         min_pair=float(min(rr.values())), max_pair=float(max(rr.values())),
                         mean_pair=float(np.mean(list(rr.values())))))
print("\n  specification curve — every cell, none dropped")
for g_ in grid:
    print(f"    {g_['spec']:<22s} {g_['arm']:<12s} n={g_['n']:6d}  pair r: min {g_['min_pair']:+.4f} "
          f"mean {g_['mean_pair']:+.4f} max {g_['max_pair']:+.4f}")

ps = [r["p"] for r in rows] + [
    2 * (1 - stats.norm.cdf(abs((obs_rawc[(i, j)] - null_mean[(i, j)]) /
                                max(null_sd[(i, j)], 1e-9)))) for (i, j) in PAIRS]
below = [r for r in rows if r["z"] < 0 and r["p"] < 0.05]
above = [r for r in rows if r["z"] > 0 and r["p"] < 0.05]

G = Gate("Is the two-layer structure a fact about out-groups, or about moral attitudes at all?")
G.plant_direction_from_sweep(f"positive: a planted trade-off drives `{ITEMS[PA]}`x`{ITEMS[PB]}` "
                             f"below its null, and g=0 sits ON that null (`#922`)",
                             [[g_, -v] for g_, v in sweep],
                             baseline=-null_mean[(PA, PB)],
                             baseline_spread=max(null_sd[(PA, PB)], 1e-4))
G.negative_control("within-person item-label permutation reproduces itself across two independent "
                   "draws — the check `#948`① asks for after three mis-specified controls",
                   reproducibility, mean_abs_dep,
                   null_spread=float(np.mean(list(null_sd.values()))),
                   null_kind="within-person item-label permutation (each person's total across the "
                             "four z-scores preserved EXACTLY)")
G.multiplicity_control("all 12 cells = 6 pairs x 2 centring arms (`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"z/{r['a'][:5]}-{r['b'][:5]}" for r in rows]
                              + [f"raw/{ITEMS[i][:5]}-{ITEMS[j][:5]}" for (i, j) in PAIRS])
G.asserted("⚠ DERIVATION CHECK: the observed residuals must sit on the algebraically forced floor "
           "-1/(k-1), or the ipsative arithmetic is not what I derived", abs(all_pairs_mean - FLOOR) < 0.01,
           f"derived {FLOOR:+.4f} · observed mean over all {len(PAIRS)} pairs {all_pairs_mean:+.4f}. "
           f"⚠ k={K} makes this floor {abs(FLOOR/FLOOR_945):.1f}x harsher than `#945`'s "
           f"{FLOOR_945:+.4f} (k={K_945}, read from its artifact), which is why no number here is "
           f"comparable to `#945`'s in absolute terms — only the DEPARTURE is",
           kind="control", population=f"GSS ballot 1, all four norms jointly, n={len(d)}")
G.asserted("⚠ SCALE CONTROL `#945` did not need: items are z-scored BEFORE person-centring, because "
           "their means span 1.339-2.982 and raw centring would make the pair structure item "
           "DIFFICULTY rather than anybody's trade-off", True,
           " · ".join(f"{r['a'][:5]}-{r['b'][:5]} z {r['r_resid']:+.3f} vs raw-centred "
                      f"{r['r_rawcentred']:+.3f}" for r in rows[:3]) + " · …",
           kind="control", population=f"GSS ballot 1, all four norms jointly, n={len(d)}")
G.asserted("⚠ CONFOUND CONTROL in the same iteration: person-centring removes exactly the general "
           "permissiveness that makes the raw pairs positive, so raw is reported beside residual",
           True,
           f"raw Pearson pair r: mean {np.mean([r['r_raw_pearson'] for r in rows]):+.4f} "
           f"[{min(r['r_raw_pearson'] for r in rows):+.4f}, "
           f"{max(r['r_raw_pearson'] for r in rows):+.4f}] · residual: mean "
           f"{np.mean([r['r_resid'] for r in rows]):+.4f}", kind="control",
           population=f"GSS ballot 1, all four norms jointly, n={len(d)}")
G.asserted("the whole grid and specification curve are published, disagreeing cells included", True,
           " · ".join(f"{g_['spec'][:3]}/{g_['arm'][:3]} mean {g_['mean_pair']:+.3f}"
                      for g_ in grid), kind="control",
           population=f"GSS ballot 1, all four norms jointly, n={len(d)}")

pos_fires = (-sweep[-1][1]) > (-sweep[0][1]) + 2 * null_sd[(PA, PB)]
neg_null = (abs(all_pairs_mean - FLOOR) < 0.01) and (reproducibility < 0.5 * mean_abs_dep)
general = len(below) > 0
world = ("W_GENERAL_LAW" if general else
         ("W_OUTGROUP_ONLY" if len(above) > 0 else "W_INVERTED"))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and the null "
           "sits on the derived floor and reproduces itself. STAKED: W_GENERAL_LAW, i.e. at least "
           "one sexual-norm pair sits significantly BELOW its own ipsative null, as 12 of 30 "
           "tolerance pairs did. ⚠ W_OUTGROUP_ONLY is the unwelcome branch — it narrows `#945` "
           "from 'how moral attitudes are held' to 'how out-groups are ranked'",
           (pos_fires and neg_null) and general,
           f"positive fires {pos_fires} · null on floor and reproducible {neg_null} "
           f"(|{all_pairs_mean:+.4f} - {FLOOR:+.4f}| = {abs(all_pairs_mean - FLOOR):.4f}; repro "
           f"{reproducibility:.5f} vs half the mean departure {0.5*mean_abs_dep:.5f}) · pairs BELOW "
           f"their null at p<0.05: {len(below)}/6 · ABOVE: {len(above)}/6 ⇒ {world}",
           kind="kill", yardstick="residual pair correlation minus its own permutation null",
           yardstick_noise=float(np.mean(list(null_sd.values()))),
           population=f"GSS ballot 1, all four norms jointly, n={len(d)}, {len(waves)} waves "
                      f"1988-2024",
           direction="two-sided; BELOW its null is a trade-off, ABOVE is co-holding")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if general else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=949, round="E03·A126·R385", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world == "W_OUTGROUP_ONLY"),
               n=int(len(d)), k=K, waves=waves, derived_floor=FLOOR,
               ref_945_floor=FLOOR_945, ref_945_k=K_945,
               observed_all_pairs_mean=all_pairs_mean, cells=rows, grid=grid,
               below=len(below), above=len(above),
               null_median=null_mean[(PA, PB)], null_sd=null_sd[(PA, PB)], null_draws=NREPS,
               null_reproducibility=reproducibility, mean_abs_departure=mean_abs_dep,
               positive_sweep=sweep, family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "the_same_test_on_the_sexual_norms.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_same_test_on_the_sexual_norms.json'}")
