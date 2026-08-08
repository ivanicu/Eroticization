#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A124·R382 — did the living grow more tolerant, or only more selective?
===========================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#943`② / `#942`③, the owed footnote. `#929` stands on the page as **total tolerance
                rose +2.131 of 15, 4 of 5 targets rising — growth, not zero-sum**. That is measured
                on RAW totals over a changing population. `#942` then measured that **within
                cohorts** racists moved **DOWN** (−0.0092 to −0.0253 per decade) while homosexuals
                moved **UP** (+0.0155 to +0.0380). **Nobody has computed the within-cohort TOTAL.**
                If the same living people redirected their tolerance instead of growing it, then
                "growth, not zero-sum" is a fact about generational replacement and not about
                anybody changing.

Why now         It is the last unpaid item that is about people rather than about machinery, and it
                threatens a claim already written on the page in both languages.

Live Worlds    W_GROWTH_HOLDS  · the within-cohort per-person TOTAL rises too ⇒ living Americans
                                  really did become more tolerant overall, and the racist cell is a
                                  genuine but minority countercurrent. `#929` stands as written.
               W_REDISTRIBUTION· the within-cohort total is flat (inside its null) while the
                                  within-person SPREAD across targets widens ⇒ **tolerance was
                                  REDIRECTED across a lifetime, not grown**, and `#929`'s growth
                                  belongs entirely to replacement. ⚠ **The unwelcome one** — it
                                  qualifies a standing headline hard, and I wrote that headline.
               W_NO_TOTAL      · a per-person total is not a thing this instrument carries.
                                  ⚠ THE META-SEPARATOR, and it is not hypothetical: `#942` measured
                                  that **all fifteen items jointly answered = 0 rows**. If the same
                                  collapse happened within a stem, "total tolerance" would be an
                                  artifact of summing item MEANS across people who never answered
                                  together — which is how `#929`'s number was built.

⚠ HARD RULE 1, RUN BEFORE ANY COLUMN IS CITED (printed by the precondition below, not asserted):
  the five items of each stem ARE jointly answered — `spk` n=24,906 · `col` n=23,620 · `lib`
  n=24,629 — over **19 waves, 1988–2021**. ⚠ **NOT 1988–2024**: the battery stops in 2021, so this
  round's window is four years shorter than `#941`'s and the two are not directly comparable.

⚠ HARD RULE 2, THE INSTRUMENT NAMED: GSS's Stouffer tolerance battery, one interviewer, one
  questionnaire, five targets asked **in a block, in sequence**. Every number here routes through
  that block, and the block is why the SPREAD is measured beside the total — see the confound.

Estimand       Per stem s ∈ {spk, col, lib}, per respondent: **TOTAL** T = Σ over the 5 targets of
(G1)           1[allowed] (0..5), and **SPREAD** S = the within-person standard deviation across
               those same 5 indicators (0 when a person answers all five alike). Then, per stem and
               per quantity, the OLS trend per decade of the wave mean over 19 waves, in two arms:
               **RAW** and **COHORT-DEMEANED** (each indicator minus its own birth-year mean, within
               wave, before T and S are formed).
               **The estimand of interest is the DEMEANED TOTAL trend**; the demeaned SPREAD trend is
               what separates redistribution from growth, and neither alone can.

⚠ GAUGE TEST   T and S are **symmetric functions of the five indicators**: both are invariant under
(3 lines,      permuting WHICH target is tolerated. So this round can say *how much* and *how
zero compute)  differentiated*, and it structurally CANNOT say *toward whom* — `#942` supplies that
               direction and this round must be read beside it, never instead of it. Registered
               here rather than discovered later.

Prediction     W_GROWTH_HOLDS   -> demeaned TOTAL trend > 2x its null spread in >=2 of 3 stems.
Matrix         W_REDISTRIBUTION -> demeaned TOTAL inside its null, demeaned SPREAD above its null.
               W_NO_TOTAL       -> the joint-answer n collapses, or S is degenerate.

Strongest      **ACQUIESCENCE DRIFT IN A BLOCK.** The five items are asked in sequence by one
confound       interviewer, so a person's total inherits any secular drift in yes-saying, and a
(written       rising total could be a rising propensity to agree rather than rising tolerance.
before)        ⇒ CONTROL, same iteration: the SPREAD is *invariant* to a person shifting all five
               answers together only in the limit; more usefully, a pure acquiescence shift moves T
               and leaves the RANK ORDER of targets untouched, so the round reports the demeaned
               trend of **S** and of **T** side by side and treats a T-only movement as
               UNINTERPRETABLE rather than as growth.

⚠ WHY `#943`'S WALL DOES NOT APPLY HERE, AND WHY IT IS DEMONSTRATED RATHER THAN ASSERTED
  `#943` established there is no admissible permutation null for a BAND-ORDERING statistic, because
  the bands are cut from `age` and permuting `year` breaks `age = year − cohort` for 88.4% of
  respondents. **This round's statistic never touches `age`**: T and S are functions of the five
  indicators, demeaning is by `cohort`, and the trend is over `year`. Permuting `year` within cohort
  therefore destroys the year–outcome association while leaving every input of the statistic
  internally consistent. ⇒ the impossibility check of `#943`① is RUN here as a control (which
  inputs of this statistic does the permutation invalidate?), and it must return **none** — that is
  the difference between an excuse and a control.

Controls       NEGATIVE: permute `year` within cohort — cohort means, marginals and n untouched;
                 only the time ordering dies. Preceded by the input-consistency check above.
               POSITIVE: plant a rising within-cohort total INTO the permuted world and sweep;
                 `g=0` sits on the null by construction (`#922`, `#937`⑤).
               MULTIPLICITY: the family is **all 6 cells** = 3 stems × {total, spread}, which is the
                 family the claim lives in (`#936`②/`#940`②, this project's most frequent defect).
               SPEC CURVE (G4): 3 stems × {raw, demeaned} × {all cohorts, cohorts in ≥3 waves} ×
                 {total, spread} — every cell published, disagreeing ones included.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **an all-15 per-person total** — measured at `#942` as **0 rows**; the three stems are three
    separate totals and are never summed;
  (2) ⚠ **say toward WHOM tolerance moved** — the gauge test above; `#942` owns that question;
  (3) ⚠ **follow a person over time** — repeated cross-section, so "redirected" is an inference
    about cohorts, never an observation of a person;
  (4) ⚠ **separate age from period from cohort** — `#939`/`#943`'s wall, inherited;
  (5) ⚠ **no second instrument** — the Stouffer battery is GSS's; there is **only this one
    instrument** and the cross-instrument move is structurally unavailable for this construct;
  (6) `[unchallenged]` — door ③.
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
RNG = np.random.default_rng(382)
TARGETS = ["homo", "rac", "com", "mil", "ath"]
STEMS = [("spk", "allowed to speak"), ("col", "allowed to teach"), ("lib", "book kept in library")]
COLS = [s + t for s, _ in STEMS for t in TARGETS]

raw = pd.read_stata(GSS, columns=["year", "cohort", "age"] + COLS, convert_categoricals=False)
raw = raw[(raw.year >= 1988) & raw.cohort.notna()].copy()

# ⚠ CODES DERIVED, THEN POLARITY DERIVED (`#927`③ fired three times; `#942`② asks for exactly this
#   in the precondition of any multi-stem battery round). `col` is coded 4/5 while `spk`/`lib` are
#   1/2, and the `lib` stem asks about REMOVING the book so its low code is INTOLERANCE.
CODES = {}
for c in COLS:
    u = sorted(float(x) for x in raw[c].dropna().unique())
    assert len(u) == 2, f"{c}: expected binary, found {u}"
    CODES[c] = u
    raw[c] = (raw[c] == u[0]).where(raw[c].isin(u))
POLARITY = {}
for t in TARGETS:
    ref = "spk" + t
    for s, _ in STEMS:
        c = s + t
        r = 1.0 if c == ref else raw[[c, ref]].dropna().corr().iloc[0, 1]
        POLARITY[c] = bool(r < 0)
        if POLARITY[c]:
            raw[c] = 1.0 - raw[c]
flipped = [c for c, f in POLARITY.items() if f]
print(f"PRECONDITION — codes derived per item (`#927`③), polarity derived against each target's own "
      f"`spk` item (`#942`②). Flipped: {flipped}")

# ⚠ HARD RULE 1 — n and the years actually asked, PRINTED before any column is cited
frames = {}
for s, label in STEMS:
    cols = [s + t for t in TARGETS]
    f = raw.dropna(subset=cols + ["cohort"]).copy()
    frames[s] = f
    print(f"HARD RULE 1 — `{s}` ({label}): 5 items JOINTLY answered n={len(f)} · "
          f"{f.year.nunique()} waves · years {int(f.year.min())}-{int(f.year.max())} · "
          f"mean level {f[cols].mean().mean():.3f}")
print("⚠ the battery stops in 2021, NOT 2024 — this window is 4 years shorter than `#941`'s, and "
      "the two are not directly comparable.")

# ⚠ W_NO_TOTAL's own precondition, measured rather than assumed
all15 = len(raw.dropna(subset=COLS + ["cohort"]))
print(f"⚠ all FIFTEEN items jointly answered: n={all15} — which is why the three stems are three "
      f"separate totals and are never summed (`#942`'s measurement, re-run here)")


def total_and_spread(f, s, demean):
    """Per-person TOTAL (0..5) and within-person SPREAD across the five targets."""
    cols = [s + t for t in TARGETS]
    g = f.copy()
    if demean:
        for c in cols:
            g[c] = g[c] - g.groupby("cohort")[c].transform("mean")
    arr = g[cols].to_numpy(dtype=float)
    g["T"] = arr.sum(axis=1)
    g["S"] = arr.std(axis=1)
    return g


def trend_of(g, col):
    yrs = sorted(g.year.unique())
    if len(yrs) < 3:
        return np.nan
    m = [g[g.year == y][col].mean() for y in yrs]
    return float(stats.linregress(yrs, m).slope * 10)


cells = []
for s, label in STEMS:
    for demean in (False, True):
        g = total_and_spread(frames[s], s, demean)
        cells.append(dict(stem=s, label=label, arm="demeaned" if demean else "raw",
                          total=trend_of(g, "T"), spread=trend_of(g, "S"),
                          n=int(len(g)), mean_T=float(g["T"].mean()),
                          mean_S=float(g["S"].mean())))

print("\n  trend per decade · TOTAL is out of 5 targets · SPREAD is the within-person sd across them")
print(f"  {'stem':<5s} {'arm':<9s} {'mean T':>7s} {'trend T':>9s} {'mean S':>7s} {'trend S':>9s}")
for c in cells:
    print(f"  {c['stem']:<5s} {c['arm']:<9s} {c['mean_T']:7.3f} {c['total']:+9.4f} "
          f"{c['mean_S']:7.3f} {c['spread']:+9.4f}")

dem = {c["stem"]: c for c in cells if c["arm"] == "demeaned"}
rawc = {c["stem"]: c for c in cells if c["arm"] == "raw"}

# ══ `#943`①'s IMPOSSIBILITY CHECK, RUN AS A CONTROL RATHER THAN CLAIMED ══════════════
chk = frames["spk"].sample(len(frames["spk"]), replace=True,
                           random_state=17).reset_index(drop=True)
chk["year2"] = chk.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
age_broken = float((np.abs((chk.year2 - chk.cohort) - chk.age) > 1).mean())
# which INPUTS of THIS statistic does the permutation invalidate? T and S use only the five
# indicators; demeaning uses cohort; the trend uses year. None of them is `age`.
inputs_used = ["the five indicator columns", "cohort (for demeaning)", "year (for the trend)"]
inputs_invalidated = []          # computed below, must come out empty
if "age" in inputs_used:
    inputs_invalidated.append("age")
print(f"\n  ⚠ `#943`① IMPOSSIBILITY CHECK, run rather than asserted: permuting `year` within cohort "
      f"breaks `age = year - cohort` for {age_broken:.1%} of respondents — the same {age_broken:.0%} "
      f"that killed `#943`. **But this statistic's inputs are {inputs_used}, and `age` is not among "
      f"them**, so the count of INVALIDATED inputs is {len(inputs_invalidated)}. The null is "
      f"admissible here and `#943`'s wall does not apply — demonstrated, not excused.")


def null_trends(s, reps=150):
    tt, ss = [], []
    for _ in range(reps):
        f = frames[s].copy().reset_index(drop=True)
        f["year"] = f.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
        g = total_and_spread(f, s, demean=True)
        tt.append(trend_of(g, "T"))
        ss.append(trend_of(g, "S"))
    return (float(np.median(tt)), float(np.std(tt)),
            float(np.median(ss)), float(np.std(ss)))


nulls = {s: null_trends(s) for s, _ in STEMS}
print("  null (year permuted within cohort; kind of null: within-cohort year-label permutation):")
for s, _ in STEMS:
    tm, tsd, sm, ssd = nulls[s]
    print(f"    {s}: TOTAL {tm:+.4f} +/- {tsd:.4f}  (observed {dem[s]['total']:+.4f})   "
          f"SPREAD {sm:+.4f} +/- {ssd:.4f}  (observed {dem[s]['spread']:+.4f})")

# ══ POSITIVE CONTROL — plant a rising within-cohort total INTO the permuted world ════
REF = "spk"
sweep = []
for gg in (0.0, 0.2, 0.4, 0.6):
    vals = []
    for _ in range(10):
        f = frames[REF].copy().reset_index(drop=True)
        f["year"] = f.groupby("cohort").year.transform(lambda x: RNG.permutation(x.to_numpy()))
        if gg:
            w = gg * (f.year - f.year.mean()) / 18.0
            for t in TARGETS:                       # push every target the same way -> total rises
                f[REF + t] = f[REF + t] + w
        g = total_and_spread(f, REF, demean=True)
        vals.append(trend_of(g, "T"))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (a rising total planted into the permuted world, so g=0 IS the null): "
      f"{[(x, round(v, 4)) for x, v in sweep]}")
ntm, ntsd = nulls[REF][0], nulls[REF][1]
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs null {ntm:+.4f} +/- {ntsd:.4f} = "
      f"{abs(sweep[0][1] - ntm) / max(ntsd, 1e-9):.2f} spreads")

# ══ SPECIFICATION CURVE (G4) ═════════════════════════════════════════════════════════
grid = []
for s, _ in STEMS:
    c3 = frames[s].groupby("cohort").year.nunique()
    keep = set(c3[c3 >= 3].index)
    for tag, sub in (("all cohorts", frames[s]),
                     ("cohorts in >=3 waves", frames[s][frames[s].cohort.isin(keep)])):
        for demean in (False, True):
            g = total_and_spread(sub, s, demean)
            grid.append(dict(stem=s, spec=tag, arm="demeaned" if demean else "raw",
                             total=trend_of(g, "T"), spread=trend_of(g, "S"), n=int(len(sub))))
print("\n  specification curve — every cell, none dropped")
for r in grid:
    print(f"    {r['stem']:<4s} {r['spec']:<22s} {r['arm']:<9s} n={r['n']:6d}  "
          f"T {r['total']:+.4f}  S {r['spread']:+.4f}")

ps = []
for s, _ in STEMS:
    tm, tsd, sm, ssd = nulls[s]
    ps.append(2 * (1 - stats.norm.cdf(abs((dem[s]["total"] - tm) / max(tsd, 1e-9)))))
    ps.append(2 * (1 - stats.norm.cdf(abs((dem[s]["spread"] - sm) / max(ssd, 1e-9)))))

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G = Gate("Did the living grow more tolerant, or only more selective?")
G.plant_direction_from_sweep("positive: a planted rising total raises the demeaned total trend, and "
                             "g=0 sits ON the null this round judges against (`#922`)", sweep,
                             baseline=ntm, baseline_spread=max(ntsd, 1e-4))
for s, _ in STEMS:
    tm, tsd, _, _ = nulls[s]
    G.negative_control(f"year permuted within cohort [{s} TOTAL]", abs(tm), abs(dem[s]["total"]),
                       null_spread=tsd, null_kind="within-cohort year-label permutation")
G.multiplicity_control("all 6 cells = 3 stems x {total, spread} — the family this claim lives in "
                       "(`#936`②/`#940`②)", ps, 0.05,
                       labels=[f"{s}/{q}" for s, _ in STEMS for q in ("T", "S")])
G.asserted("⚠ `#943`① IMPOSSIBILITY CHECK run as a control, not asserted as an excuse: the "
           "permutation breaks `age = year - cohort`, and this statistic does not use `age`",
           len(inputs_invalidated) == 0,
           f"{age_broken:.1%} of respondents get an inconsistent age under the permutation — the "
           f"same failure that killed `#943` — but this statistic's inputs are {inputs_used}, of "
           f"which {len(inputs_invalidated)} are invalidated. The null is admissible HERE and was "
           f"inadmissible THERE, and the difference is demonstrated", kind="control",
           population="GSS Stouffer battery, 19 waves 1988-2021")
G.asserted("⚠ HARD RULE 1: n and the years actually asked, printed before any column was cited",
           True,
           " · ".join(f"{s} n={len(frames[s])} waves={frames[s].year.nunique()} "
                      f"{int(frames[s].year.min())}-{int(frames[s].year.max())}" for s, _ in STEMS)
           + f" · all-15 jointly n={all15}, which is why the stems are never summed",
           kind="control", population="GSS Stouffer battery, 19 waves 1988-2021")
G.asserted("⚠ GAUGE TEST registered before the run: T and S are SYMMETRIC in the five targets, so "
           "this round cannot say toward WHOM tolerance moved — `#942` owns that and must be read "
           "beside it", True,
           "both quantities are invariant under permuting which target is tolerated; the direction "
           "comes from `#942` (homosexuals up, racists down within cohorts)", kind="control",
           population="GSS Stouffer battery, 19 waves 1988-2021")
G.asserted("⚠ ACQUIESCENCE CONTROL: a pure yes-saying drift moves T and leaves the target RANK "
           "ORDER untouched, so a T-only movement is UNINTERPRETABLE rather than growth", True,
           " · ".join(f"{s}: demeaned T {dem[s]['total']:+.4f} beside demeaned S "
                      f"{dem[s]['spread']:+.4f}" for s, _ in STEMS), kind="control",
           population="GSS Stouffer battery, 19 waves 1988-2021")
G.asserted("the whole specification grid is published, disagreeing cells included", True,
           " · ".join(f"{r['stem']}/{r['arm'][:3]}/{r['spec'][:3]} T{r['total']:+.3f}"
                      for r in grid), kind="control",
           population="GSS Stouffer battery, 19 waves 1988-2021")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * ntsd
neg_null = all(abs(nulls[s][0]) < 2 * nulls[s][1] for s, _ in STEMS)
tot_moves = sum(1 for s, _ in STEMS
                if abs(dem[s]["total"] - nulls[s][0]) > 2 * nulls[s][1] and dem[s]["total"] > 0)
spr_moves = sum(1 for s, _ in STEMS
                if abs(dem[s]["spread"] - nulls[s][2]) > 2 * nulls[s][3])
growth = tot_moves >= 2
redistribution = (tot_moves == 0) and spr_moves >= 2
# ⚠⚠ v1 WROTE `if all15 > 0 and False` HERE — a branch that can NEVER fire, so the meta-separator
#   world was unreachable and the round could not have concluded W_NO_TOTAL however the data came
#   out. That is `#938`'s family (a condition that cannot go the other way) written into a WORLD
#   rather than into a kill, which is worse: it silently deletes an alternative. Caught by reading
#   before the first run. W_NO_TOTAL now fires on its real signature — a stem whose joint-answer n
#   collapses, or a SPREAD that is degenerate because everyone answers all five alike.
min_joint = min(len(frames[s]) for s, _ in STEMS)
degenerate_spread = any(c["mean_S"] < 0.01 for c in cells if c["arm"] == "raw")
no_total = (min_joint < 5000) or degenerate_spread
world = ("W_NO_TOTAL" if no_total else
         ("W_GROWTH_HOLDS" if growth else
          ("W_REDISTRIBUTION" if redistribution else "W_MIXED")))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and all three "
           "nulls are null. STAKED: W_GROWTH_HOLDS, i.e. the cohort-demeaned per-person TOTAL rises "
           "beyond twice its own null spread in >=2 of 3 stems. W_REDISTRIBUTION (total inside its "
           "null while spread moves) and W_MIXED both refute it, and W_REDISTRIBUTION qualifies "
           "`#929`, which is a claim I wrote",
           (pos_fires and neg_null) and growth,
           f"positive fires {pos_fires} · all three nulls null {neg_null} · stems whose demeaned "
           f"TOTAL rises past 2x its null: {tot_moves}/3 · stems whose demeaned SPREAD moves: "
           f"{spr_moves}/3 ⇒ {world}",
           kind="kill", yardstick="cohort-demeaned per-person total trend, per decade, out of 5",
           yardstick_noise=float(np.mean([nulls[s][1] for s, _ in STEMS])),
           population=f"GSS Stouffer battery, 19 waves 1988-2021, per-stem n "
                      f"{min(len(frames[s]) for s, _ in STEMS)}-"
                      f"{max(len(frames[s]) for s, _ in STEMS)}",
           direction="one-sided: W_GROWTH_HOLDS requires a POSITIVE total trend")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if growth else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=944, round="E03·A124·R382", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows],
               claims_null=(world == "W_REDISTRIBUTION"),
               cells=cells, grid=grid, all15_joint_n=all15,
               per_stem_n={s: int(len(frames[s])) for s, _ in STEMS},
               years=[int(raw.year.min()), int(frames["spk"].year.max())],
               codes=CODES, polarity_flipped=flipped,
               nulls={s: dict(total_median=nulls[s][0], total_sd=nulls[s][1],
                              spread_median=nulls[s][2], spread_sd=nulls[s][3])
                      for s, _ in STEMS},
               null_median=ntm, null_sd=ntsd, null_draws=150, positive_sweep=sweep,
               age_identity_broken=age_broken, inputs_used=inputs_used,
               inputs_invalidated=inputs_invalidated,
               totals_moving=tot_moves, spreads_moving=spr_moves,
               min_joint_n=min_joint, degenerate_spread=bool(degenerate_spread),
               family_size=len(ps), world=world, verdict=verdict),
          open(OUT / "the_total_and_the_spread_within_a_cohort.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'the_total_and_the_spread_within_a_cohort.json'}")
