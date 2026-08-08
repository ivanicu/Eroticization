#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A123·R380 — five targets of tolerance, each split into conversion and replacement
======================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#941`①. R379 measured a **conversion share** per sexual norm and found the four did
                not move by the same mechanism: `premarsx` 0.192, `homosex` 0.487. A115 has a second
                object entirely — GSS's tolerance battery, five TARGETS × three STEMS — where `#929`
                established that total tolerance rose **+2.131 of 15** with 4 of 5 targets rising,
                i.e. **growth, not zero-sum**. That is a SIGN and a SIZE. **It has no mechanism.**
                If tolerance of homosexuals is high-conversion while tolerance of racists is
                high-replacement, then "growth not zero-sum" means something quite specific: the
                same society did two different things to two of its out-groups.

Why now         The machinery exists and is one round old, the object is different (tolerance, not
                sexual permissiveness), and the answer changes how `#929` must be worded.

Live Worlds    W_TARGET · the conversion share varies across the five TARGETS ⇒ Americans changed
                           their own minds about some out-groups and merely outlived their opinion
                           of others. **A mechanism for `#929`.**
               W_ERA    · all five share about the same conversion share ⇒ tolerance rose as one
                           era-wide process and the target is irrelevant to HOW, only to how much.
               W_STEM   · the share depends more on which STEM is asked (speak / teach / library)
                           than on the target ⇒ **the quantity is a property of the QUESTION, not of
                           the attitude**, and R379's per-norm shares deserve the same suspicion.
                           ⚠ The unwelcome one, and it is the meta-separator: it would say the
                           conversion share is an artifact of item format.

Estimand       For each of 5 targets × 3 stems: trend per decade of the proportion answering
(G1)           "allowed", over GSS waves 1988–2024, in two arms — RAW and COHORT-DEMEANED (item
               minus its own birth-year mean, within wave). Conversion share = demeaned ÷ raw.
               ⚠ The **stem is a specification axis, not a nuisance to average away** — W_STEM is a
               live world and averaging first would make it unobservable.

Prediction     W_TARGET -> between-target spread of shares exceeds the between-stem spread.
Matrix         W_ERA    -> both spreads small.
               W_STEM   -> between-stem spread equals or exceeds between-target spread.

Strongest      **A RATIO WITH A SMALL DENOMINATOR** (`#918`, and `#941`② is an open debt about
confound       exactly this). Any cell whose RAW trend fails to clear twice its own bootstrap spread
(written       is marked UNREADABLE and given no share. ⚠ Second: these are BINARY items, so a
before)        proportion near 0 or 1 has compressed variance and a floor/ceiling on its trend —
               reported as the raw level beside every share.

Controls       NEGATIVE: permute `year` within cohort — marginals, cohort means and n untouched.
               POSITIVE: plant a rising within-cohort trend INTO the permuted world, `g=0` on the
                 null by construction (`#922`, `#937`⑤).
               MULTIPLICITY: the family is **all 15 cells**, not the five targets — `#936`②/`#940`②
                 is this project's most frequent live defect and the grid IS the family here.
               SPEC CURVE (G4): all 15 cells published, disagreeing ones included.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate age from period from cohort** — `#939`'s wall, inherited;
  (2) ⚠ **follow a person over time** — repeated cross-section;
  (3) ⚠ **no second instrument** — this battery is GSS's; `#937` measured that the cross-instrument
    move needs a comparable quantity and there is **only this one instrument** here;
  (4) ⚠ **claim a causal decomposition** — cohort-demeaning identifies within-birth-year movement,
    never why;
  (5) `[unchallenged]` — door ③.
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
GSS = ROOT / "data" / "external" / "gss" / "GSS_stata" / "gss7224_r3a.dta"
RNG = np.random.default_rng(380)
TARGETS = [("homo", "homosexuals"), ("rac", "racists"), ("com", "communists"),
           ("mil", "militarists"), ("ath", "anti-religionists")]
STEMS = [("spk", "speak"), ("col", "teach"), ("lib", "library book")]
COLS = [s + t for s, _ in STEMS for t, _ in TARGETS]

d = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
d = d[(d.year >= 1988) & d.cohort.notna()]

# ⚠⚠ `#927`③ REPEATED VERBATIM AND CAUGHT BY THE PRECONDITION, NOT BY ME. v1 wrote `(d[c] == 1)` for
#   all fifteen items because every one of them is LABELLED "allowed / not allowed". **The `col`
#   stem's CODES ARE 4 AND 5**, not 1 and 2 — which is the exact fact `#927`③ recorded ("a label
#   list is not a code list") after `colhomo` cost a round there. v1 therefore silently produced
#   FIVE all-false columns and zero usable rows. The codes are now DERIVED from the data per column
#   and asserted, so the next stem GSS codes differently cannot repeat this.
CODES = {}
for c in COLS:
    u = sorted(float(x) for x in d[c].dropna().unique())
    assert len(u) == 2, f"{c}: expected a binary item, found codes {u}"
    CODES[c] = u
    d[c] = (d[c] == u[0]).where(d[c].isin(u))            # provisional: lower code -> 1
print("PRECONDITION — codes DERIVED per item, never assumed (`#927`③): "
      + " · ".join(f"{c}{tuple(int(x) for x in v)}" for c, v in list(CODES.items())[::5]))

# ⚠⚠⚠ AND DERIVING THE CODES IS NOT DERIVING THE DIRECTION — `#927`③ one level deeper, caught by
#   the data and not by me. After the code fix, tolerance of homosexuals read speak 0.840, teach
#   0.797, **library 0.247**. That 0.247 cannot be a tolerance rate for the same target on the same
#   people: the `lib` stem asks whether to **REMOVE** the book, so its low code is INTOLERANCE and
#   the whole stem is sign-flipped. Assuming "lower = allowed" was a second label-vs-code error
#   inside the fix for the first one.
#   ⇒ POLARITY IS DERIVED TOO: `spk` is the reference, and any item whose correlation with its own
#   target's `spk` item is NEGATIVE is flipped. This uses only the data, needs no codebook, and it
#   is its own positive control — a stem that agrees needs no flip and must not get one.
POLARITY = {}
for tcode, _ in TARGETS:
    ref = "spk" + tcode
    for scode, _ in STEMS:
        c = scode + tcode
        r = d[[c, ref]].dropna().corr().iloc[0, 1] if c != ref else 1.0
        POLARITY[c] = bool(r < 0)
        if POLARITY[c]:
            d[c] = 1.0 - d[c]
flipped = [c for c, f in POLARITY.items() if f]
print(f"PRECONDITION — polarity DERIVED against each target's `spk` item: flipped {flipped or 'none'}")
assert all(not POLARITY["spk" + tc] for tc, _ in TARGETS), "the reference stem must never flip"

# ⚠ AND THE 15 ITEMS ARE NOT JOINTLY ANSWERED: requiring all fifteen leaves ZERO rows. Each cell
#   uses its own complete cases, so n differs per cell and is reported per cell rather than once.
waves = sorted(int(y) for y in d.year.unique())
n_per = {c: int(d[c].notna().sum()) for c in COLS}
print(f"PRECONDITION — {len(waves)} waves {waves[0]}-{waves[-1]} · per-item n "
      f"{min(n_per.values())}-{max(n_per.values())} (⚠ all 15 jointly = 0, so cells are NOT the "
      f"same people and no within-person comparison is made)")


def trends(f, col):
    g = f[f[col].notna()].copy()          # per-cell complete cases; see the precondition note
    g["dm"] = g[col] - g.groupby("cohort")[col].transform("mean")
    yrs = sorted(g.year.unique())
    raw = [g[g.year == y][col].mean() for y in yrs]
    dem = [g[g.year == y]["dm"].mean() for y in yrs]
    return (float(stats.linregress(yrs, raw).slope * 10),
            float(stats.linregress(yrs, dem).slope * 10),
            float(g[col].mean()))


cells = []
for s, sl in STEMS:
    for t, tl in TARGETS:
        r, m, lvl = trends(d, s + t)
        cells.append(dict(stem=sl, target=tl, col=s + t, raw=r, demeaned=m, level=lvl,
                          n=int(d[s + t].notna().sum()), codes=CODES[s + t]))

B = 150
bootc = {c["col"]: {"raw": [], "dem": []} for c in cells}
for _ in range(B):
    smp = d.sample(len(d), replace=True, random_state=int(RNG.integers(1e9)))
    for c in cells:
        r, m, _ = trends(smp, c["col"])
        bootc[c["col"]]["raw"].append(r)
        bootc[c["col"]]["dem"].append(m)
for c in cells:
    bc = bootc[c["col"]]
    c["raw_sd"] = float(np.std(bc["raw"]))
    c["readable"] = bool(abs(c["raw"]) > 2 * c["raw_sd"])
    c["share"] = float(c["demeaned"] / c["raw"]) if c["readable"] else float("nan")
    rr = [m / r for r, m in zip(bc["raw"], bc["dem"]) if abs(r) > 1e-6]
    c["share_sd"] = float(np.std(rr)) if rr else float("nan")

print("\n  trend per decade in proportion 'allowed'; conversion share = demeaned / raw")
print(f"  {'stem':<13s} {'target':<18s} {'level':>6s} {'raw':>9s} {'demeaned':>9s} {'share':>16s}")
for c in cells:
    sh = f"{c['share']:.3f} ±{c['share_sd']:.3f}" if c["readable"] else "⚠ UNREADABLE"
    print(f"  {c['stem']:<13s} {c['target']:<18s} {c['level']:6.3f} {c['raw']:+9.4f} "
          f"{c['demeaned']:+9.4f} {sh:>16s}")

ok = [c for c in cells if c["readable"]]
by_t = {tl: [c["share"] for c in ok if c["target"] == tl] for _, tl in TARGETS}
by_s = {sl: [c["share"] for c in ok if c["stem"] == sl] for _, sl in STEMS}
t_means = {k: float(np.mean(v)) for k, v in by_t.items() if v}
s_means = {k: float(np.mean(v)) for k, v in by_s.items() if v}
spread_t = float(max(t_means.values()) - min(t_means.values())) if len(t_means) > 1 else np.nan
spread_s = float(max(s_means.values()) - min(s_means.values())) if len(s_means) > 1 else np.nan
print(f"\n  by TARGET (averaged over stems): "
      + " · ".join(f"{k} {v:.3f}" for k, v in sorted(t_means.items(), key=lambda x: -x[1])))
print(f"  by STEM   (averaged over targets): "
      + " · ".join(f"{k} {v:.3f}" for k, v in sorted(s_means.items(), key=lambda x: -x[1])))
print(f"  between-TARGET spread {spread_t:.3f} · between-STEM spread {spread_s:.3f}")

# ══ NEGATIVE CONTROL ═════════════════════════════════════════════════════════════════
REF = "spkhomo"
nulls = {}
for c in cells:
    nulls[c["col"]] = []
for _ in range(80):
    s = d.copy()
    s["year"] = s.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
    for c in cells:
        nulls[c["col"]].append(trends(s, c["col"])[1])
nullstat = {k: (float(np.median(v)), float(np.std(v))) for k, v in nulls.items()}
n_null_ok = sum(1 for k, (m_, sd_) in nullstat.items() if abs(m_) < 2 * sd_)
print(f"\n  null (year permuted WITHIN cohort; kind of null: within-cohort year-label permutation): "
      f"{n_null_ok}/15 cells null · {REF} {nullstat[REF][0]:+.4f} +/- {nullstat[REF][1]:.4f}")

# ══ POSITIVE CONTROL ═════════════════════════════════════════════════════════════════
sweep = []
for gg in (0.0, 0.15, 0.30, 0.45):
    vals = []
    for _ in range(6):
        s = d.copy()
        s["year"] = s.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
        s[REF] = s[REF] + gg * (s.year - np.mean(waves)) / 18.0
        vals.append(trends(s, REF)[1])
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (rising within-cohort trend planted into the permuted world, g=0 IS the "
      f"null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs null {nullstat[REF][0]:+.4f} +/- "
      f"{nullstat[REF][1]:.4f} = "
      f"{abs(sweep[0][1] - nullstat[REF][0]) / max(nullstat[REF][1], 1e-9):.2f} spreads")

ps = [2 * (1 - stats.norm.cdf(abs((c["demeaned"] - nullstat[c["col"]][0]) /
                                  max(nullstat[c["col"]][1], 1e-9)))) for c in cells]

G = Gate("Did Americans change their minds about some out-groups and merely outlive others?")
G.plant_direction_from_sweep("positive: a planted within-cohort trend raises the demeaned trend, "
                             "g=0 sits ON the null this round judges against (`#922`)", sweep,
                             baseline=nullstat[REF][0], baseline_spread=max(nullstat[REF][1], 1e-4))
G.negative_control(f"year permuted within cohort [{REF}]", abs(nullstat[REF][0]),
                   abs([c for c in cells if c["col"] == REF][0]["demeaned"]),
                   null_spread=nullstat[REF][1], null_kind="within-cohort year-label permutation")
G.multiplicity_control("all 15 target x stem cells — the grid IS the family (`#936`②/`#940`②, this "
                       "project's most frequent live defect)", ps, 0.05,
                       labels=[f"{c['stem'][:4]}/{c['target'][:5]}" for c in cells])
G.asserted("⚠ a ratio with a small denominator is not given a share (`#918`, open debt `#941`②)",
           True, f"{len(ok)}/15 cells readable; unreadable: "
                 f"{[c['col'] for c in cells if not c['readable']] or 'none'}", kind="control",
           population=f"GSS {len(waves)} waves 1988-2024, n={len(d)}")
G.asserted("⚠ the STEM is a specification axis, not averaged away — W_STEM is a live world and "
           "averaging first would make it unobservable", True,
           f"between-TARGET spread {spread_t:.3f} vs between-STEM spread {spread_s:.3f}",
           kind="control", population=f"GSS {len(waves)} waves 1988-2024, n={len(d)}")
G.asserted("the whole 15-cell grid is published, disagreeing cells included", True,
           " · ".join(f"{c['stem'][:3]}/{c['target'][:4]} {c['share']:+.2f}" for c in ok),
           kind="control", population=f"GSS {len(waves)} waves 1988-2024, n={len(d)}")
G.asserted("⚠ codes AND polarity both DERIVED per item, never assumed from the label (`#927`③, "
           "twice in one round): `col` is 4/5 while `spk`/`lib` are 1/2, and the `lib` stem asks "
           "about REMOVING the book so its low code is INTOLERANCE", True,
           f"codes {[(c, tuple(int(x) for x in v)) for c, v in list(CODES.items())[::5]]} · "
           f"flipped against each target's own `spk` item: {flipped or 'none'} · the tell was "
           f"tolerance of homosexuals reading 0.840/0.797/0.247 across three stems",
           kind="control", population=f"GSS {len(waves)} waves 1988-2024")
G.asserted("⚠ NOT a causal decomposition: cohort-demeaning identifies within-birth-year movement, "
           "never why anyone moved; and the 15 cells are NOT the same people (all 15 jointly = 0 "
           "rows), so no within-person comparison is made", True,
           f"per-item n {min(n_per.values())}-{max(n_per.values())}", kind="control",
           population=f"GSS {len(waves)} waves 1988-2024, n={len(d)}")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * nullstat[REF][1]
neg_null = n_null_ok >= 13
target_wins = spread_t > spread_s
world = ("W_STEM" if (not np.isnan(spread_s) and spread_s >= spread_t) else
         ("W_TARGET" if spread_t > 0.10 else "W_ERA"))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and >=13/15 "
           "nulls are null. STAKED: W_TARGET, i.e. the between-TARGET spread of conversion shares "
           "exceeds the between-STEM spread AND exceeds 0.10. W_STEM and W_ERA both refute it",
           (pos_fires and neg_null) and target_wins and spread_t > 0.10,
           f"positive fires {pos_fires} · nulls null {n_null_ok}/15 · between-TARGET "
           f"{spread_t:.3f} vs between-STEM {spread_s:.3f} -> "
           f"{'target' if target_wins else 'STEM'} dominates ⇒ {world}",
           kind="kill", yardstick="between-target spread of conversion shares, vs between-stem",
           yardstick_noise=float(np.nanmean([c["share_sd"] for c in ok])),
           population=f"GSS {len(waves)} waves 1988-2024, n={len(d)}",
           direction="one-sided: W_TARGET requires target > stem")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if (target_wins and spread_t > 0.10) else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=942, round="E03·A123·R380", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_ERA"),
               waves=waves, n=int(len(d)), n_per_item=n_per, codes=CODES,
               polarity_flipped=flipped, cells=cells, by_target=t_means, by_stem=s_means,
               spread_target=spread_t, spread_stem=spread_s,
               nulls={k: dict(median=v[0], sd=v[1]) for k, v in nullstat.items()},
               null_median=nullstat[REF][0], null_sd=nullstat[REF][1], null_draws=len(nulls[REF]),
               positive_sweep=sweep, family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "five_targets_of_tolerance_same_control.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'five_targets_of_tolerance_same_control.json'}")
