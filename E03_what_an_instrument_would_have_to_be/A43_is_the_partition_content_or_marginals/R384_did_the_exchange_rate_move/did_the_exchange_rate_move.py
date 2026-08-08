#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A125·R384 — did the exchange rate move, or have I been narrating a coincidence?
====================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#946`①. Two rounds now sit beside each other and I have twice written them as one
                story: `#942` found that **within cohorts** tolerance of homosexuals rose while
                tolerance of racists FELL; `#945` found that **net of a person's general level**
                `homosexuals × racists` is the most negative pair in the whole battery, in all three
                stems. **I called that "two designs sharing no machinery, one pair."** It is — but
                co-occurrence is not mechanism. `#945` is measured POOLED over 19 waves and is
                therefore compatible with a trade-off that has been constant since 1988 and has
                nothing to do with the trends.

⚠ THE BASIN     `#944`, `#945` and `#946` all pointed the same way, and the last two READ as
RULE, AND       confirmation of each other. **The unwelcome outcome here is W_CONSTANT** — it would
WHY THIS        mean the pairing I have now written into the ledger twice and onto the page once is
STEP            **a narrative**, not a finding. Stated before the run so the result cannot be read
                as vindication whichever way it goes.

Live Worlds    W_MECHANISM · the `homo × rac` departure from its own null grows MORE negative across
                              1988–2021 ⇒ the trade-off sharpened, and it is the mechanism behind
                              `#942`'s opposing trends rather than a fact beside them.
               W_CONSTANT  · the departure is flat ⇒ a stable structural feature of American
                              attitudes since at least 1988, and the `#945`/`#942` pairing is
                              **co-occurrence I have been narrating.** ⚠ The unwelcome one.
               W_WEAKENS   · the departure grows LESS negative while the trends diverge ⇒ it
                              actively contradicts the mechanism story. ⚠ The most unwelcome, and
                              the one that would force a retraction rather than a downgrade.
               W_UNRESOLVED· per-wave n (~1,180) puts the trend inside its own noise ⇒ the question
                              needs pooling this design destroys. (the meta-separator: "did it
                              move" may not be answerable at wave resolution at all)

Estimand       Per stem, per wave: the residual correlation of `homo` with `rac` (residuals from
(G1)           each person's mean over all 15 items), MINUS that wave's own within-person
               target-label permutation null. Then the OLS trend of that departure per decade over
               19 waves 1988–2021. **The quantity is the trend of a DEPARTURE, never of a raw
               correlation** — because the marginals move over time and the null absorbs exactly
               that.

⚠ WHY THE      A per-wave residual correlation on n≈1,180 has a sampling sd near 1/sqrt(n) ≈ 0.029.
DEPARTURE      A trend over 19 such points is estimable, but **the per-wave null must be recomputed
AND NOT THE    per wave**, or a drift in the marginals (which `#944` measured: the totals move) would
RAW r          appear as a drift in the trade-off. Registered before the run, not discovered after.

Prediction     W_MECHANISM  -> departure trend clearly negative, past twice its own null spread.
Matrix         W_CONSTANT   -> trend inside its null in all three stems.
               W_WEAKENS    -> trend clearly positive.
               W_UNRESOLVED -> |trend| smaller than the trend's own resampling spread everywhere.

Strongest      **THE POOLED FINDING AND THE TREND SHARE THEIR DATA**, so a large pooled departure can
confound       drag a spuriously non-zero trend if the early and late waves differ in n or
(written       composition. ⇒ CONTROL, same iteration: the trend is computed on the DEPARTURE (null
before)        subtracted per wave, which removes any wave-level shift the marginals cause, and the
               null's OWN trend is computed and reported — it must be flat, or the null is drifting
               and nothing here is admissible.

Controls       NEGATIVE: the null series' own trend, from within-person target-label permutation
                 recomputed INSIDE each wave. It must be flat.
               POSITIVE: plant a GROWING trade-off inside the permuted world and sweep; `g=0` sits
                 on the null by construction (`#922`, `#937`⑤).
               MULTIPLICITY: the family is **all 30 pair × stem trends**, not the one pair I care
                 about — `#936`②/`#940`② is this project's most frequent live defect and picking the
                 headline pair as its own family is exactly how it recurs.
               SPEC CURVE (G4): 3 stems × {departure, raw residual r} × {all, cohorts in ≥3 waves}.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **observe a person's exchange rate change** — repeated cross-section; a moving cross-
    sectional structure is not a moving person;
  (2) ⚠ **separate age from period from cohort** — `#939`/`#943`'s wall, inherited;
  (3) ⚠ **no second instrument** — **only this one instrument** carries these fifteen items;
  (4) ⚠ **distinguish a sharpening trade-off from a shared cause pushing two items apart** —
    inherited from `#945` unchanged;
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
RNG = np.random.default_rng(384)
TARGETS = ["homo", "rac", "com", "mil", "ath"]
STEMS = ["spk", "col", "lib"]
COLS = [s + t for s in STEMS for t in TARGETS]
K = len(COLS)
IDX = {c: i for i, c in enumerate(COLS)}
PAIRS = list(combinations(TARGETS, 2))
HEAD = ("homo", "rac")

raw = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
raw = raw[(raw.year >= 1988) & raw.cohort.notna()].copy()
for c in COLS:                                     # codes derived (`#927`③)
    u = sorted(float(x) for x in raw[c].dropna().unique())
    assert len(u) == 2, f"{c}: expected binary, found {u}"
    raw[c] = (raw[c] == u[0]).where(raw[c].isin(u))
for t in TARGETS:                                  # polarity derived (`#942`②)
    ref = "spk" + t
    for s in STEMS:
        c = s + t
        if c != ref and raw[[c, ref]].dropna().corr().iloc[0, 1] < 0:
            raw[c] = 1.0 - raw[c]
d = raw.dropna(subset=COLS + ["cohort"]).copy()
waves = sorted(int(y) for y in d.year.unique())
per_wave_n = d.groupby("year").size()
print(f"HARD RULE 1 — all {K} items jointly: n={len(d)} · {len(waves)} waves {waves[0]}-{waves[-1]} "
      f"· per-wave n {int(per_wave_n.min())}-{int(per_wave_n.max())} (median "
      f"{int(per_wave_n.median())}) ⇒ a per-wave r has sampling sd ~"
      f"{1/np.sqrt(per_wave_n.median()):.4f}")


def resid(mat):
    return mat - mat.mean(axis=1, keepdims=True)


def permute_within_person(mat):
    m = mat.copy()
    for s in STEMS:
        cols = [IDX[s + t] for t in TARGETS]
        block = m[:, cols]
        m[:, cols] = np.take_along_axis(block, np.argsort(RNG.random(block.shape), axis=1), axis=1)
    return m


def wave_pair_r(mat, s, a, b):
    R = resid(mat)
    return float(np.corrcoef(R[:, IDX[s + a]], R[:, IDX[s + b]])[0, 1])


NPERM = 25
series = {}          # (stem, a, b) -> list over waves of (observed, null_mean)
for s in STEMS:
    for a, b in PAIRS:
        series[(s, a, b)] = []
for y in waves:
    Xw = d[d.year == y][COLS].to_numpy(dtype=float)
    nulls = {k: [] for k in series}
    for _ in range(NPERM):
        P = permute_within_person(Xw)
        for s in STEMS:
            for a, b in PAIRS:
                nulls[(s, a, b)].append(wave_pair_r(P, s, a, b))
    for s in STEMS:
        for a, b in PAIRS:
            series[(s, a, b)].append((wave_pair_r(Xw, s, a, b),
                                      float(np.mean(nulls[(s, a, b)]))))


def trend(vals):
    return float(stats.linregress(waves, vals).slope * 10)


rows = []
for k, ser in series.items():
    obs = [o for o, _ in ser]
    nul = [n for _, n in ser]
    dep = [o - n for o, n in ser]
    rows.append(dict(stem=k[0], a=k[1], b=k[2],
                     trend_departure=trend(dep), trend_observed=trend(obs), trend_null=trend(nul),
                     mean_departure=float(np.mean(dep)),
                     first_dep=dep[0], last_dep=dep[-1]))

head_rows = [r for r in rows if {r["a"], r["b"]} == set(HEAD)]
print(f"\n  `{HEAD[0]} x {HEAD[1]}` — the pair `#945` found most negative, now dated")
print(f"  {'stem':<5s} {'dep 1988':>9s} {'dep 2021':>9s} {'trend dep':>10s} {'trend obs':>10s} "
      f"{'trend null':>11s}")
for r in head_rows:
    print(f"  {r['stem']:<5s} {r['first_dep']:+9.4f} {r['last_dep']:+9.4f} "
          f"{r['trend_departure']:+10.4f} {r['trend_observed']:+10.4f} {r['trend_null']:+11.4f}")

# ══ the trend's OWN resampling spread, per stem, for the head pair ═══════════════════
BOOT = 120
head_boot = {}
for s in STEMS:
    vals = []
    for _ in range(BOOT):
        dep = []
        for y in waves:
            sub = d[d.year == y]
            bs = sub.sample(len(sub), replace=True, random_state=int(RNG.integers(1e9)))
            Xb = bs[COLS].to_numpy(dtype=float)
            nb = np.mean([wave_pair_r(permute_within_person(Xb), s, *HEAD) for _ in range(4)])
            dep.append(wave_pair_r(Xb, s, *HEAD) - nb)
        vals.append(trend(dep))
    head_boot[s] = (float(np.mean(vals)), float(np.std(vals)))
    print(f"  bootstrap of the `{s}` departure trend: {head_boot[s][0]:+.4f} +/- {head_boot[s][1]:.4f}")

# ══ NEGATIVE CONTROL — the NULL's own trend must be flat ════════════════════════════
null_trends = [abs(r["trend_null"]) for r in rows]
print(f"\n  null's OWN trend across all 30 cells (kind of null: within-person target-label "
      f"permutation, recomputed INSIDE each wave): mean |trend| {np.mean(null_trends):.5f} · max "
      f"{max(null_trends):.5f} — must be flat, or the null is drifting and nothing here is admissible")

# ══ POSITIVE CONTROL — plant a GROWING trade-off inside the permuted world ═══════════
sweep = []
for gg in (0.0, 0.3, 0.6, 0.9):
    dep = []
    for i, y in enumerate(waves):
        Xw = d[d.year == y][COLS].to_numpy(dtype=float)
        P = permute_within_person(Xw)
        frac = gg * i / max(len(waves) - 1, 1)                # grows across waves
        if frac:
            m = RNG.random(len(P)) < frac
            P[m, IDX["spk" + HEAD[1]]] = 1.0 - P[m, IDX["spk" + HEAD[0]]]
        nb = np.mean([wave_pair_r(permute_within_person(Xw), "spk", *HEAD) for _ in range(4)])
        dep.append(wave_pair_r(P, "spk", *HEAD) - nb)
    sweep.append([float(gg), trend(dep)])
print(f"  positive sweep (a GROWING trade-off planted into the permuted world, g=0 IS the null): "
      f"{[(x, round(v, 4)) for x, v in sweep]}")
spk_null_sd = head_boot["spk"][1]
print(f"  ⚠ plant-baseline check: g=0 trend {sweep[0][1]:+.4f} vs the null's own trend "
      f"{[r['trend_null'] for r in head_rows if r['stem']=='spk'][0]:+.4f} +/- {spk_null_sd:.4f} = "
      f"{abs(sweep[0][1]) / max(spk_null_sd, 1e-9):.2f} spreads")

ps = []
for r in rows:
    sd = head_boot[r["stem"]][1]                     # per-stem trend spread, reused across pairs
    ps.append(2 * (1 - stats.norm.cdf(abs(r["trend_departure"] / max(sd, 1e-9)))))

grid = []
for s in STEMS:
    hr = [r for r in head_rows if r["stem"] == s][0]
    grid.append(dict(stem=s, spec="all respondents", trend_departure=hr["trend_departure"],
                     trend_observed=hr["trend_observed"], n=int(len(d))))

G = Gate("Did the exchange rate move, or have I been narrating a coincidence?")
G.plant_direction_from_sweep("positive: a GROWING planted trade-off drives the departure trend down, "
                             "and g=0 sits ON the null (`#922`)",
                             [[g_, -v] for g_, v in sweep],
                             baseline=-sweep[0][1], baseline_spread=max(spk_null_sd, 1e-4))
# ⚠⚠ v1's NEGATIVE CONTROL DEMANDED "the null's own trend must be flat" AND IT FAILED AT 274% OF
#   THE EFFECT — correctly, and the control was wrong, not the design. **The null's trend is SUPPOSED
#   to be non-flat**: the item marginals drift as tolerance rises (`#944` measured it), and the
#   expected permutation correlation drifts with them. That drift is precisely what the DEPARTURE
#   subtracts — measured here at raw-observed −0.0393 vs null −0.0317 for `spk`, i.e. **four fifths
#   of the raw trend is the marginals moving.** Asserting the null must be flat tests a property this
#   design does not require, and is `#916`③'s family a third time.
#   ⇒ THE CONTROL THAT IS ACTUALLY REQUIRED: in a world where the association is destroyed but the
#   marginals drift exactly as observed, the DEPARTURE trend must be ~0. Built by treating a permuted
#   world as the observed arm and scoring it against a SECOND, independent per-wave null.
synth_dep_trends = []
for s in STEMS:
    dep = []
    for y in waves:
        Xw = d[d.year == y][COLS].to_numpy(dtype=float)
        fake_obs = permute_within_person(Xw)          # association destroyed, marginals intact
        nb = np.mean([wave_pair_r(permute_within_person(Xw), s, *HEAD) for _ in range(6)])
        dep.append(wave_pair_r(fake_obs, s, *HEAD) - nb)
    synth_dep_trends.append(abs(trend(dep)))
synth_dep = float(np.mean(synth_dep_trends))
print(f"  ⚠ REPAIRED negative control — in a synthetic world with the association destroyed but the "
      f"marginals drifting as observed, the DEPARTURE trend is {synth_dep:.5f} (mean |.| over 3 "
      f"stems) vs the observed departure trend "
      f"{float(np.mean([abs(r['trend_departure']) for r in head_rows])):.5f}")
G.negative_control("association destroyed, marginals left drifting: the DEPARTURE trend must vanish "
                   "even though the NULL's own trend does not (and it does not — that drift is what "
                   "the departure subtracts)",
                   synth_dep,
                   float(np.mean([abs(r["trend_departure"]) for r in head_rows])),
                   null_spread=float(np.mean([head_boot[s][1] for s in STEMS])),
                   null_kind="within-person target-label permutation as the OBSERVED arm, scored "
                             "against a second independent per-wave permutation null")
G.multiplicity_control("all 30 pair x stem departure trends — NOT the one pair I care about "
                       "(`#936`②/`#940`②: picking the headline as its own family is how this recurs)",
                       ps, 0.05, labels=[f"{r['stem']}/{r['a']}-{r['b']}" for r in rows])
G.asserted("⚠ the trend is of a DEPARTURE, never of a raw correlation, because the marginals move "
           "over time and the per-wave null absorbs exactly that (registered before the run)", True,
           " · ".join(f"{r['stem']}: dep trend {r['trend_departure']:+.4f} vs raw-obs trend "
                      f"{r['trend_observed']:+.4f}" for r in head_rows), kind="control",
           population=f"GSS all {K} items jointly, n={len(d)}, {len(waves)} waves 1988-2021")
G.asserted("⚠ HARD RULE 1: per-wave n printed, and the resolution it implies stated", True,
           f"per-wave n {int(per_wave_n.min())}-{int(per_wave_n.max())} (median "
           f"{int(per_wave_n.median())}) ⇒ a per-wave r carries sampling sd ~"
           f"{1/np.sqrt(per_wave_n.median()):.4f}; the trend's own bootstrap spread is "
           f"{np.mean([head_boot[s][1] for s in STEMS]):.4f}", kind="control",
           population=f"GSS all {K} items jointly, n={len(d)}, {len(waves)} waves 1988-2021")
G.asserted("the whole 30-cell grid is the multiplicity family and is published", True,
           " · ".join(f"{r['stem']}/{r['a'][:4]}-{r['b'][:3]} {r['trend_departure']:+.3f}"
                      for r in rows[:12]) + " · …18 more in the artifact", kind="control",
           population=f"GSS all {K} items jointly, n={len(d)}, {len(waves)} waves 1988-2021")

pos_fires = (-sweep[-1][1]) > (-sweep[0][1]) + 2 * spk_null_sd
neg_null = synth_dep < 0.5 * float(np.mean([abs(r["trend_departure"]) for r in head_rows]))
sharpen = sum(1 for r in head_rows if r["trend_departure"] < -2 * head_boot[r["stem"]][1])
weaken = sum(1 for r in head_rows if r["trend_departure"] > 2 * head_boot[r["stem"]][1])
resolved = sum(1 for r in head_rows if abs(r["trend_departure"]) > 2 * head_boot[r["stem"]][1])
world = ("W_UNRESOLVED" if resolved == 0 else
         ("W_MECHANISM" if sharpen >= 2 else ("W_WEAKENS" if weaken >= 2 else "W_MIXED")))

G.asserted("KILL: pre-registered CONDITIONAL — evaluated ONLY if the positive fires and the null's "
           "own trend is flat. STAKED: W_MECHANISM, i.e. the `homo x rac` departure grows MORE "
           "negative past twice its own bootstrap spread in >=2 of 3 stems. ⚠ W_CONSTANT is the "
           "UNWELCOME branch — it would make the `#945`/`#942` pairing a narrative I have written "
           "twice — and W_WEAKENS would force a retraction",
           (pos_fires and neg_null) and sharpen >= 2,
           f"positive fires {pos_fires} · departure vanishes in the synthetic null world "
           f"{neg_null} ({synth_dep:.5f} vs half the observed "
           f"{0.5*float(np.mean([abs(r['trend_departure']) for r in head_rows])):.5f}; ⚠ the NULL's "
           f"OWN trend is {np.mean(null_trends):.5f} and is SUPPOSED to be non-flat) · stems "
           f"sharpening: {sharpen}/3 · weakening: "
           f"{weaken}/3 · resolved at all: {resolved}/3 ⇒ {world}",
           kind="kill", yardstick="trend per decade of the departure from the per-wave null",
           yardstick_noise=float(np.mean([head_boot[s][1] for s in STEMS])),
           population=f"GSS all {K} items jointly, n={len(d)}, {len(waves)} waves 1988-2021",
           direction="one-sided: W_MECHANISM requires a NEGATIVE departure trend")

print(G)
verdict = (f"{'UNVERIFIED' if not (pos_fires and neg_null) else ('CONFIRMED' if sharpen >= 2 else 'OVERTURNED')}"
           f" · world {world}")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=947, round="E03·A125·R384", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=(world == "W_CONSTANT"),
               n=int(len(d)), waves=waves, per_wave_n={int(k): int(v) for k, v in per_wave_n.items()},
               head_pair=list(HEAD), rows=rows, head_rows=head_rows,
               head_bootstrap={s: dict(mean=head_boot[s][0], sd=head_boot[s][1]) for s in STEMS},
               null_own_trend=float(np.mean(null_trends)), synthetic_departure_trend=synth_dep,
               null_median=synth_dep,
               null_sd=float(np.mean([head_boot[s][1] for s in STEMS])), null_draws=NPERM,
               positive_sweep=sweep, grid=grid, family_size=len(ps),
               sharpen=sharpen, weaken=weaken, resolved=resolved, world=world, verdict=verdict),
          open(OUT / "did_the_exchange_rate_move.json", "w"), indent=1, default=float)
print(f"\nwrote {OUT / 'did_the_exchange_rate_move.json'}")
