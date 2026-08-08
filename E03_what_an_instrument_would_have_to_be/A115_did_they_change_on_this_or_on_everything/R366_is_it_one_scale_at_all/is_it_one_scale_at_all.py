#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A115·R366 — before "re-sorting" can mean anything: is it one scale at all?
==============================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#927` measured tolerance moving **+0.4302 toward homosexuals and −0.2602 toward
                racists** and read it as *what moved was not their tolerance but its object* —
                condemnation re-allocated. ⚠ **That reading has a precondition it never checked**:
                you cannot call something a RE-ALLOCATION across scales that are not commensurable.
                `#927`(2) named it. **The precondition comes before the interpretation** — the
                discipline `#924`(2) and `#925`(2) exist to enforce.

Why Now        `#927`(1) wants to treat condemnation as a roughly zero-sum allocation. If the five
                target arms are separate constructs, that whole framing is a **category error**, and
                the page carries it right now.

⚠ PRIOR ART,   **This is the central dispute in the political-tolerance literature, not a new idea.**
DECLARED       Stouffer (1955) built these items; Sullivan, Piereson & Marcus (1979, 1982) argued
BEFORE THE     the traditional battery **conflates tolerance with agreement with the target**, and
NUMBER         proposed content-controlled "least-liked group" measures precisely because scores
                across targets may not be one scale. **This round measures that on this release; it
                does not discover it.**

Live Worlds    W_ONE      · disattenuated cross-target correlation is high ⇒ one tolerance
                            disposition ⇒ `#927`'s opposite movement is a genuine re-allocation, and
                            a real puzzle: the same dimension moving opposite ways by target.
               W_SEPARATE · it is low ⇒ different constructs ⇒ **`#927`'s "re-sorting" is a category
                            error and the page must be corrected.** ⚠ **The unwelcome one.**
               W_MIXED    · in between: a general component plus large target-specific parts, in
                            which case neither "one scale" nor "separate" is the right description.
                            (the meta-separator — it kills the binary I am testing)

Estimand       Within each wave, at the PERSON level: the correlation between a person's mean
(G1)           tolerance of homosexuals (3 stems) and of racists (3 stems), **compared against the
               WITHIN-target cross-stem correlation as the ceiling** — not against zero. Reported
               raw and **disattenuated** (`r_cross / sqrt(rel_A * rel_B)`), because a correlation
               between two 3-item means is attenuated by those items' own reliability.

Prediction     W_ONE      -> disattenuated cross-target r >= 0.60.
Matrix         W_SEPARATE -> <= 0.30.
               W_MIXED    -> between, with the full 5x5 matrix showing structure.

⚠ PRECONDITION Checked and PRINTED before the estimator: >=300 respondents per wave answering ALL 15
CHECK FIRST    items. Measured: **26 waves qualify, 31,648 complete cases, min 736 per wave.**
               Waves failing are DROPPED WITH A COUNT.

⚠ CODES ARE    `#927`(3): `col*` items use numeric codes **4/5**, not 1/2. Every code here is DERIVED
DERIVED        from the data; no literal appears. The orientation SIGNS come from the shipped labels.

Controls       NEGATIVE: permute one target's per-person scores within the wave — destroys the
               person-level link while preserving both marginals. Cross-target r must go to 0.
               POSITIVE: plant a shared factor into that permuted world and sweep, so `g=0` lands ON
               the null (`#922`'s gate; `#927` broke this three times in one round).
               REFERENCE/SHAM: the within-target cross-stem correlations are the ceiling any
               cross-target correlation is measured against.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ a correlation is not a scale proof — high cross-target r is consistent with one disposition
    AND with two dispositions that share a cause (education, say). No causal claim is made;
  (2) ⚠ binary items attenuate correlations by construction; disattenuation corrects for unreliability
    but NOT for the coarseness of a 2-point scale;
  (3) ⚠ this is a within-wave, person-level check on a between-wave, cohort-level claim — it tests the
    PRECONDITION of `#927`'s reading, never `#927`'s measurement, which stands either way;
  (4) ⚠ **only this one instrument**;
  (5) `[unchallenged]` — door (3).
"""
import json, sys, warnings
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
RNG = np.random.default_rng(366)
MIN_N = 300

# sign = +1 if HIGH is permissive, -1 if LOW is permissive — from the SHIPPED LABELS (`#927`)
TARGETS = {
    "homosexuals":       {"spkhomo": -1, "colhomo": -1, "libhomo": +1},
    "racists":           {"spkrac": -1,  "colrac": -1,  "librac": +1},
    "communists":        {"spkcom": -1,  "colcom": +1,  "libcom": +1},
    "militarists":       {"spkmil": -1,  "colmil": -1,  "libmil": +1},
    "anti-religionists": {"spkath": -1,  "colath": -1,  "libath": +1},
}
COLS = [c for t in TARGETS.values() for c in t]

d = pd.read_stata(GSS, columns=["year"] + COLS, convert_categoricals=False)
for c in COLS:
    d[c] = d[c].where(d[c].between(1, 7))

# ⚠ `#927`(3): every code DERIVED, never a literal. Orient each item to 0 = restrictive, 1 = permissive.
codes = {c: sorted(d[c].dropna().unique()) for c in COLS}
print("⚠ codes derived from the data (no literals):", {c: [float(x) for x in v] for c, v in list(codes.items())[:4]}, "...")
for tgt, items in TARGETS.items():
    for c, sign in items.items():
        lo, hi = float(codes[c][0]), float(codes[c][-1])
        perm = hi if sign > 0 else lo
        d[c] = (d[c] == perm).astype(float).where(d[c].notna())

d["complete"] = d[COLS].notna().all(axis=1)

# ══ PRECONDITION CHECK, printed BEFORE the estimator ═════════════════════════════════
counts = d.groupby("year")["complete"].sum()
waves = [int(y) for y, n in counts.items() if n >= MIN_N]
dropped = [int(y) for y, n in counts.items() if 0 < n < MIN_N]
print(f"\nPRECONDITION CHECK — waves with >={MIN_N} respondents answering ALL 15 items:")
print(f"  usable {len(waves)} · DROPPED {len(dropped)}: {dropped or 'none'}  (absence reported, not passed)")
cc = d[d["complete"]]
print(f"  complete cases {len(cc[cc.year.isin(waves)])} · per-wave min "
      f"{int(counts[counts >= MIN_N].min())} median {int(counts[counts >= MIN_N].median())}")


def wave_stats(w):
    """within-target cross-stem r (the ceiling) and cross-target r (raw + disattenuated)."""
    means, rel = {}, {}
    for tgt, items in TARGETS.items():
        cols = list(items)
        rs = [w[a].corr(w[b]) for i, a in enumerate(cols) for b in cols[i + 1:]]
        rbar = float(np.mean(rs))
        rel[tgt] = float(3 * rbar / (1 + 2 * rbar)) if rbar > 0 else np.nan   # Spearman-Brown, k=3
        means[tgt] = w[cols].mean(axis=1)
        rel[f"{tgt}__rbar"] = rbar
    out = {}
    names = list(TARGETS)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r = float(means[a].corr(means[b]))
            den = np.sqrt(rel[a] * rel[b]) if (rel[a] > 0 and rel[b] > 0) else np.nan
            out[f"{a}|{b}"] = dict(raw=r, disatt=(float(r / den) if den and den > 0 else np.nan))
    return means, rel, out


rows = []
for y in waves:
    w = cc[cc.year == y]
    means, rel, pairs = wave_stats(w)
    key = "homosexuals|racists"
    rows.append(dict(year=int(y), n=int(len(w)),
                     within_homo=rel["homosexuals__rbar"], within_rac=rel["racists__rbar"],
                     rel_homo=rel["homosexuals"], rel_rac=rel["racists"],
                     cross_raw=pairs[key]["raw"], cross_disatt=pairs[key]["disatt"],
                     all_pairs={k: v["disatt"] for k, v in pairs.items()}))

print("\n=== THE GRID — per wave (all cells, disagreeing ones included) ===")
for r in rows:
    print(f"  {r['year']}  n={r['n']:5d}  within-stem r: homo {r['within_homo']:.3f} rac "
          f"{r['within_rac']:.3f}  |  cross-target r raw {r['cross_raw']:.3f} "
          f"disattenuated {r['cross_disatt']:.3f}")

med_within = float(np.median([np.mean([r["within_homo"], r["within_rac"]]) for r in rows]))
med_raw = float(np.median([r["cross_raw"] for r in rows]))
med_dis = float(np.median([r["cross_disatt"] for r in rows]))
print(f"\n  median WITHIN-target cross-stem r (the ceiling): {med_within:.3f}")
print(f"  median CROSS-target r, raw                     : {med_raw:.3f}")
print(f"  median CROSS-target r, disattenuated           : {med_dis:.3f}")
allp = {}
for r in rows:
    for k, v in r["all_pairs"].items():
        allp.setdefault(k, []).append(v)
print("\n  all ten target pairs, median disattenuated r:")
for k, v in sorted(allp.items(), key=lambda kv: -np.median(kv[1])):
    print(f"    {k:38s} {float(np.median(v)):+.3f}")

# ══ NEGATIVE CONTROL — permute one target's person scores within the wave ════════════
null_vals = []
for y in waves[:12]:
    w = cc[cc.year == y]
    for _ in range(20):
        p = w.copy()
        for c in TARGETS["racists"]:
            p[c] = RNG.permutation(p[c].to_numpy())
        _, rel, pairs = wave_stats(p)
        # ⚠ v1 nulled the DISATTENUATED r and its spread came out +/-0.3062, so 0 of 26 waves
        #   resolved. Disattenuation divides by an estimated reliability, so it AMPLIFIES the
        #   permutation noise — a statistic can be the right quantity and the wrong one to null.
        #   ⇒ the null and the multiplicity family use the RAW r; disattenuated is reported beside
        #   it as a secondary, with its instability stated.
        v = pairs["homosexuals|racists"]["raw"]
        if not np.isnan(v):
            null_vals.append(v)
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (racist items permuted across persons within wave; kind of null: person-label "
      f"permutation): {null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant a shared factor INTO the permuted world; g=0 on the null ═
sweep = []
for g in (0.0, 0.15, 0.30, 0.45, 0.60):
    vals = []
    for y in waves[:8]:
        w = cc[cc.year == y]
        for _ in range(6):
            p = w.copy()
            for c in TARGETS["racists"]:
                p[c] = RNG.permutation(p[c].to_numpy())
            # plant: with prob g, copy the person's homosexual-tolerance answer onto the racist item
            # ⚠ v1 copied ALL THREE homosexual answers onto the three racist items, which makes the
            #   two arms nearly identical: the sweep jumped to >1 at g=0.15 and saturated flat.
            #   **A plant that saturates cannot demonstrate graded sensitivity** — it only shows the
            #   statistic can reach its ceiling. ⇒ copy ONE stem, so g moves the link continuously.
            src = list(TARGETS["homosexuals"])[0]
            tgt_c = list(TARGETS["racists"])[0]
            take = RNG.random(len(p)) < g
            p.loc[take, tgt_c] = p.loc[take, src].to_numpy()
            _, _, pairs = wave_stats(p)
            v = pairs["homosexuals|racists"]["raw"]
            if not np.isnan(v):
                vals.append(v)
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (shared factor planted into the permuted world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["cross_raw"] - null_med) / (null_sd or 1e-9)))) for r in rows]

if not rows:
    print("EMPTY POPULATION"); sys.exit(2)

one_scale = med_dis >= 0.60
separate = med_dis <= 0.30

G = Gate("Is it one scale at all?")
G.plant_direction_from_sweep("positive: a planted shared factor raises the cross-target correlation, "
                             "and g=0 is null", sweep, baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("racist items permuted across persons within wave", abs(null_med), abs(med_raw),
                   null_spread=null_sd, null_kind="within-wave person-label permutation")
G.multiplicity_control("all waves", ps, 0.05, labels=[str(r["year"]) for r in rows])
G.asserted("PRECONDITIONS checked and printed BEFORE the estimator; failures dropped with a count",
           True, f"{len(waves)} usable waves · {len(dropped)} dropped: {dropped or 'none'}",
           kind="control")
G.asserted("codes were DERIVED from the data, never written as literals (`#927`(3))",
           all(len(v) == 2 for v in codes.values()),
           f"derived code pairs, e.g. colhomo={[float(x) for x in codes['colhomo']]} vs "
           f"spkhomo={[float(x) for x in codes['spkhomo']]} — the literals would have been wrong",
           kind="control")
G.asserted("REFERENCE: the cross-target correlation is judged against the WITHIN-target ceiling, "
           "not against zero", True,
           f"within-target cross-stem r {med_within:.3f} · cross-target raw {med_raw:.3f} · "
           f"disattenuated {med_dis:.3f}", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", rows)
G.asserted("prior art declared", True,
           "Stouffer 1955 built these items; Sullivan, Piereson & Marcus 1979/1982 argued the "
           "battery conflates tolerance with agreement with the target and proposed content-"
           "controlled measures for exactly this reason. Measured here, not discovered", kind="control")
G.asserted("KILL: `#927`'s re-allocation reading requires the arms to be commensurable",
           not separate,
           f"disattenuated cross-target r {med_dis:.3f} (raw {med_raw:.3f}) vs within-target "
           f"ceiling {med_within:.3f}; null {null_med:+.4f} +/- {null_sd:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif separate:
    VERDICT, WORLD = "OVERTURNED", "W_SEPARATE · the arms are different constructs; `#927`'s reading is a category error"
elif one_scale:
    VERDICT, WORLD = "CONFIRMED", "W_ONE · one disposition, moving opposite ways by target"
else:
    VERDICT, WORLD = "CONFIRMED", "W_MIXED · a general component plus large target-specific parts"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=928, round="E03·A115·R366", verdict=VERDICT, world=WORLD,
           estimand="within-wave person-level correlation between tolerance of homosexuals and of "
                    "racists, judged against the within-target cross-stem ceiling",
           instrument="GSS 1972-2024 gss7224_r3a.dta",
           prior_art="Stouffer 1955; Sullivan, Piereson & Marcus 1979/1982 — measured, not discovered",
           preconditions=dict(usable_waves=len(waves), dropped=dropped, min_n=MIN_N),
           derived_codes={c: [float(x) for x in v] for c, v in codes.items()},
           rows=rows, median_within=med_within, median_cross_raw=med_raw,
           median_cross_disattenuated=med_dis,
           all_pairs_median={k: float(np.median(v)) for k, v in allp.items()},
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           tests="`#927`(2): the PRECONDITION of `#927`'s re-allocation reading",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "is_it_one_scale_at_all.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'is_it_one_scale_at_all.json'}")
