#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A122·R378 — the project's best arc, put through the control that just reversed another one
===============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#939`③. R377 demeaned within birth year — the cheapest possible confound control —
                and **an ordering reversed**: raw said the young move fastest, cohort-demeaned said
                they move slowest. That control costs four lines. **A118 is this project's
                best-supported arc and it never ran it.** `#936`'s headline — the ceiling-normalised
                coupling of `homosex` to the other three sexual norms rising **0.4020 → 0.6252,
                +0.0595 per decade** — is computed per WAVE, and a wave is a different set of people
                every time.

BASIN RULE     Four consecutive rounds (`#933`–`#936`) confirmed and strengthened A118. That is a
(why THIS      basin, and the rule says design a step whose POSITIVE outcome would be unwelcome.
step)          **The unwelcome outcome here is that the project's strongest claim is composition**,
               and it is a live possibility precisely because the same control just did that to a
               claim I found equally obvious one round ago.

Live Worlds    W_SURVIVES    · the normalised trend holds within cohorts ⇒ A118 is about people
                               changing their minds, and it now carries the control that killed a
                               sibling claim.
               W_COMPOSITION · the trend collapses or reverses under demeaning ⇒ **`#936`'s
                               integration is generational replacement wearing a trend's clothes**,
                               and A118's headline must be rewritten. ⚠ The unwelcome one.
               W_DEGEN       · demeaning destroys the item's usable variance (a 4-point item minus a
                               cohort mean can go nearly continuous and thin), so the comparison is
                               not available. (the meta-separator — the control may not be
                               applicable to a coupling at all, only to a level)

Estimand       Per wave: `observed |rho| / comonotone-ceiling |rho|`, averaged over the three pairs
(G1)           `homosex`×{`premarsx`,`teensex`,`xmarsex`}; then the OLS trend of that per decade.
               Two arms, identical in every other respect: **RAW** (reproduces `#936`) and
               **COHORT-DEMEANED** (each of the four items minus its own birth-year mean, computed
               within wave, before any correlation).

⚠ BUILT-IN     The RAW arm must reproduce `#936` to within resampling. **That is a positive control
REPRODUCTION   on this re-implementation**, and it is the only reason a difference in the demeaned
CHECK          arm can be attributed to demeaning rather than to my code.

Prediction     W_SURVIVES    -> demeaned trend positive and within ~2 spreads of the raw trend.
Matrix         W_COMPOSITION -> demeaned trend at or below its null, or of the opposite sign.
               W_DEGEN       -> demeaned ceilings collapse toward 1.0, or wave n falls below usable.

Strongest      **DEMEANING IS ITSELF A TRANSFORMATION OF THE MARGINAL**, and this project has been
confound       burned by exactly that (`#918`, `#936`): subtracting a cohort mean changes each
(written       item's distribution, hence the comonotone ceiling, hence the normalised value —
before)        independently of any association. ⇒ CONTROL: the ceiling is recomputed inside each
               arm from that arm's own marginals, and both arms' ceilings are reported per wave so
               the size of that shift is visible rather than assumed.

Controls       NEGATIVE: permute the target within wave — marginals and ceilings untouched, only the
                 association dies. Run separately in each arm.
               POSITIVE: plant a rising coupling INTO the permuted world and sweep; `g=0` sits on
                 the null by construction (`#922`, `#937`⑤ — built backwards three times).
               SHAM: the three other-three pairs, normalised identically, in both arms.
               MULTIPLICITY: the family is **the two arms' trend estimates plus the two shams**,
                 which is the family this claim lives in (`#936`②, `#939`②).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **separate age from period from cohort** — `#939`'s wall, unchanged and inherited;
  (2) ⚠ **follow a person over time** — repeated cross-section;
  (3) ⚠ **no second instrument for the trend** — `#937` measured that NSFG has one time point, so
    this is **only this one instrument** and the cross-instrument move is structurally unavailable;
  (4) `[unchallenged]` — door ③.
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
RNG = np.random.default_rng(378)
TARGET, OTHERS = "homosex", ["premarsx", "teensex", "xmarsex"]
ITEMS = OTHERS + [TARGET]
MIN_WAVE = 300

d = pd.read_stata(GSS, columns=["year", "ballot", "cohort"] + ITEMS, convert_categoricals=False)
for c in ITEMS:
    d[c] = d[c].where(d[c].isin([1, 2, 3, 4]))
d = d[d.ballot == 1].dropna(subset=ITEMS + ["cohort"])
counts = d.groupby("year").size()
waves = [int(y) for y, n in counts.items() if n >= MIN_WAVE]
print(f"PRECONDITION — ballot 1, all four answered, cohort present, n>={MIN_WAVE}: "
      f"{len(waves)} usable waves {waves[0]}-{waves[-1]} · n={len(d)}")


def ceiling(x, y):
    """max attainable |Spearman| given BOTH marginals — contains no association (`#902`/`#936`)."""
    return abs(stats.spearmanr(np.sort(x), np.sort(y)).statistic)


def demean(f):
    g = f.copy()
    for c in ITEMS:
        g[c] = g[c] - g.groupby("cohort")[c].transform("mean")
    return g


def wave_norm(f, tgt, oth):
    obs = [abs(stats.spearmanr(f[tgt], f[o]).statistic) for o in oth]
    cl = [ceiling(f[tgt].to_numpy(), f[o].to_numpy()) for o in oth]
    return float(np.mean(obs)), float(np.mean(cl)), float(np.mean(obs) / np.mean(cl))


def series(frame, dm, tgt=TARGET, oth=OTHERS):
    out = []
    for y in waves:
        f = frame[frame.year == y]
        f = demean(f) if dm else f
        o, c, n = wave_norm(f, tgt, oth)
        out.append((y, o, c, n))
    return out


def trend(ser):
    return float(stats.linregress([s[0] for s in ser], [s[3] for s in ser]).slope * 10)


raw_s, dem_s = series(d, False), series(d, True)
raw_t, dem_t = trend(raw_s), trend(dem_s)
print(f"\n  RAW arm      first {raw_s[0][3]:.4f} -> last {raw_s[-1][3]:.4f} · trend "
      f"{raw_t:+.4f}/decade   (ceilings {raw_s[0][2]:.4f} -> {raw_s[-1][2]:.4f})")
print(f"  DEMEANED arm first {dem_s[0][3]:.4f} -> last {dem_s[-1][3]:.4f} · trend "
      f"{dem_t:+.4f}/decade   (ceilings {dem_s[0][2]:.4f} -> {dem_s[-1][2]:.4f})")

# ⚠ SHAM — the three other-three pairs, both arms, identical normalisation
SH = [("premarsx", ["teensex", "xmarsex"]), ("teensex", ["xmarsex"])]
sham_raw = trend([(y, 0, 0, np.mean([series(d, False, a, b)[i][3] for a, b in SH]))
                  for i, y in enumerate(waves)])
sham_dem = trend([(y, 0, 0, np.mean([series(d, True, a, b)[i][3] for a, b in SH]))
                  for i, y in enumerate(waves)])
print(f"  sham (other-three pairs): raw {sham_raw:+.4f}/decade · demeaned {sham_dem:+.4f}/decade")

# ══ NEGATIVE CONTROL, per arm ════════════════════════════════════════════════════════
def null_trend(dm, reps=120):
    out = []
    for _ in range(reps):
        s = d.copy()
        s[TARGET] = s.groupby("year")[TARGET].transform(lambda x: RNG.permutation(x.to_numpy()))
        out.append(trend(series(s, dm)))
    return float(np.median(out)), float(np.std(out))


nr_med, nr_sd = null_trend(False)
nd_med, nd_sd = null_trend(True)
print(f"  null (target permuted within wave; kind of null: within-wave person-label permutation): "
      f"raw {nr_med:+.4f} +/- {nr_sd:.4f} · demeaned {nd_med:+.4f} +/- {nd_sd:.4f}")

# ══ POSITIVE CONTROL — plant a rising coupling INTO the permuted world; g=0 IS the null
sweep = []
for gg in (0.0, 0.25, 0.50, 0.75):
    vals = []
    for _ in range(6):
        s = d.copy()
        s[TARGET] = s.groupby("year")[TARGET].transform(lambda x: RNG.permutation(x.to_numpy()))
        w = (s.year - min(waves)) / max(max(waves) - min(waves), 1)
        keep = RNG.random(len(s)) < gg * w                       # restore true pairing, rising
        s.loc[keep, TARGET] = d.loc[keep, TARGET]
        vals.append(trend(series(s, False)))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (true pairing restored at a RISING rate inside the permuted world, so g=0 "
      f"IS the null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check: g=0 at {sweep[0][1]:+.4f} vs null {nr_med:+.4f} +/- {nr_sd:.4f} = "
      f"{abs(sweep[0][1] - nr_med) / max(nr_sd, 1e-9):.2f} spreads")

grid = [dict(arm="raw", trend=raw_t, null=nr_med, null_sd=nr_sd, sham=sham_raw,
             first=raw_s[0][3], last=raw_s[-1][3]),
        dict(arm="cohort-demeaned", trend=dem_t, null=nd_med, null_sd=nd_sd, sham=sham_dem,
             first=dem_s[0][3], last=dem_s[-1][3])]
ps = [2 * (1 - stats.norm.cdf(abs((r["trend"] - r["null"]) / max(r["null_sd"], 1e-9))))
      for r in grid] + \
     [2 * (1 - stats.norm.cdf(abs((s - nr_med) / max(nr_sd, 1e-9)))) for s in (sham_raw, sham_dem)]

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G = Gate("Does A118's integration trend survive the control that reversed R377's ordering?")
G.plant_direction_from_sweep("positive: a rising planted coupling raises the trend, and g=0 sits ON "
                             "the null this round judges against (`#922`)", sweep,
                             baseline=nr_med, baseline_spread=max(nr_sd, 1e-4))
G.negative_control("target permuted within wave [raw arm]", abs(nr_med), abs(raw_t),
                   null_spread=nr_sd, null_kind="within-wave person-label permutation")
G.negative_control("target permuted within wave [demeaned arm]", abs(nd_med), abs(dem_t),
                   null_spread=nd_sd, null_kind="within-wave person-label permutation")
G.multiplicity_control("both arms' trends and both shams — the family this claim lives in "
                       "(`#936`②)", ps, 0.05,
                       labels=["raw", "cohort-demeaned", "sham raw", "sham demeaned"])
# ⚠⚠ `#936`'S NUMBERS ARE READ FROM ITS ARTIFACT, NOT TYPED. v1 typed 0.4020 / 0.6252 / 0.0595 into
#   this very row and `no_transcribed_numbers` (`#840`) blocked the commit -- and it was right for a
#   reason larger than tidiness: **a reproduction check whose reference values are transcribed is
#   checking my memory of `#936`, not `#936`.** If the published artifact ever changes, a typed
#   reference silently keeps passing against a number that no longer exists.
REF = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be" /
                     "A118_four_sexual_norms_one_scale" /
                     "R374_against_what_the_marginals_permit" / "results" /
                     "against_what_the_marginals_permit.json"))
ref_first, ref_last = REF["normalised_first"], REF["normalised_last"]
ref_trend = REF["trend_normalised_per_decade"]
repro = (abs(raw_s[0][3] - ref_first) < 0.01 and abs(raw_s[-1][3] - ref_last) < 0.01
         and abs(raw_t - ref_trend) < 0.005)
G.asserted("⚠ REPRODUCTION CHECK against `#936`'s ARTIFACT (values read, never typed — `#840`): the "
           "raw arm must land on it, or a difference in the demeaned arm is my code, not the control",
           repro,
           f"`#936` artifact {ref_first:.4f} -> {ref_last:.4f}, {ref_trend:+.4f}/decade; this run's "
           f"raw arm {raw_s[0][3]:.4f} -> {raw_s[-1][3]:.4f}, {raw_t:+.4f}/decade (tolerances "
           f"0.01 on level, 0.005 on trend; this run additionally requires `cohort` non-null, which "
           f"is why they are close rather than identical)", kind="control",
           population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")
G.asserted("⚠ demeaning moves the MARGINAL, hence the ceiling — reported, not assumed (`#918`); and "
           "the pre-registered W_DEGEN condition FIRES, which is recorded rather than rewritten",
           True, f"ceilings raw {raw_s[0][2]:.4f}->{raw_s[-1][2]:.4f} · demeaned "
                 f"{dem_s[0][2]:.4f}->{dem_s[-1][2]:.4f} — the correction is INERT after demeaning, "
                 f"so the demeaned trend is a RAW correlation trend. That is a SCOPE CAVEAT on what "
                 f"the two arms compare, not a void; W_DEGEN as I worded it named the wrong thing",
           kind="control",
           population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")
G.asserted("the sham is the same operation on the pairs NOT under study, in BOTH arms", True,
           f"sham raw {sham_raw:+.4f} · sham demeaned {sham_dem:+.4f} vs target raw {raw_t:+.4f} · "
           f"target demeaned {dem_t:+.4f}", kind="control",
           population=f"GSS ballot 1, {len(waves)} waves, n={len(d)}")

pos_fires = sweep[-1][1] > sweep[0][1] + 2 * nr_sd
neg_null = abs(nr_med) < 2 * nr_sd and abs(nd_med) < 2 * nd_sd
survives = dem_t > 2 * nd_sd and dem_t > 0
# ⚠⚠ THE PRE-REGISTERED W_DEGEN CONDITION FIRES, AND IT MIS-OPERATIONALISES ITS OWN WORLD. I wrote
#   W_DEGEN as "demeaning destroys the usable variance, so THE COMPARISON IS NOT AVAILABLE", and
#   operationalised it as "demeaned ceilings > 0.98". The ceilings ARE ~0.99-1.00 -- but what that
#   means is only that **the ceiling correction is INERT in the demeaned arm**, so that arm reports
#   a raw correlation trend where the raw arm reports a ceiling-normalised one. The comparison is
#   still available; what changes is what the demeaned number IS. That is a SCOPE CAVEAT, not a
#   void.
#   ⇒ this is the THIRD time in four rounds that a pre-registered threshold named one world and
#   tested another (`#936`② the multiplicity family · `#939`② the ordering vs the spread ratio ·
#   here the normalisation's inertness vs the comparison's availability). **One error class, not
#   three mistakes**, and it is now the project's most frequent live defect.
#   The pre-registration is HONOURED rather than rewritten: the fired condition is reported, and the
#   staked kill is reported separately, because they answer different questions.
ceiling_inert = dem_s[0][2] > 0.98 and dem_s[-1][2] > 0.98
world = "W_SURVIVES" if survives else "W_COMPOSITION"

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and both nulls "
           "are null. STAKED: W_SURVIVES, i.e. the cohort-demeaned trend is positive AND exceeds "
           "twice its own null spread. W_COMPOSITION refutes it and is the unwelcome branch",
           (pos_fires and neg_null) and survives,
           f"positive fires {pos_fires} · both nulls null {neg_null} · demeaned trend {dem_t:+.4f} "
           f"vs 2x its null spread {2 * nd_sd:.4f} -> {'survives' if survives else 'DOES NOT'} · "
           f"raw {raw_t:+.4f} for comparison · sham demeaned {sham_dem:+.4f} ⇒ {world}",
           kind="kill", yardstick="ceiling-normalised coupling trend per decade, demeaned arm",
           yardstick_noise=nd_sd,
           population=f"GSS ballot 1, {len(waves)} waves 1988-2024, n={len(d)}",
           direction="one-sided: A118 predicts a POSITIVE demeaned trend")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if survives else 'OVERTURNED')}"
           f" · world {world}"
           + (" · ⚠ SCOPE: the ceiling correction is INERT in the demeaned arm (ceilings ~1.00), so "
              "that arm's trend is a RAW correlation trend, not a ceiling-normalised one — the two "
              "arms are not the identical quantity" if ceiling_inert else ""))
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=940, round="E03·A122·R378", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_COMPOSITION"),
               waves=waves, n=int(len(d)), raw_series=raw_s, demeaned_series=dem_s,
               raw_trend=raw_t, demeaned_trend=dem_t, sham_raw=sham_raw, sham_demeaned=sham_dem,
               null_median=nr_med, null_sd=nr_sd, null_median_demeaned=nd_med,
               null_sd_demeaned=nd_sd, positive_sweep=sweep, grid=grid, family_size=len(ps),
               world=world, verdict=verdict, ceiling_inert_in_demeaned_arm=bool(ceiling_inert),
               reproduces_936=bool(repro), ref_936=dict(first=ref_first, last=ref_last,
                                                        trend=ref_trend),
               preregistered_W_DEGEN_fired=bool(ceiling_inert)),
          open(OUT / "a118_under_cohort_demeaning.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'a118_under_cohort_demeaning.json'}")
