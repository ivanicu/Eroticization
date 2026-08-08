#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A112·R354 — is `#915`'s asymmetry about the DOMAIN, or about how broad the item is?
=======================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#915` measured that own same-sex experience moves BOTH norms (-0.21, -0.13) while
                own non-marital birth moves only its own (-0.10 vs +0.005). I wrote that as
                "crossing a sexual line loosens the whole moral frame". ⚠ **There is a deflationary
                rival that predicts the identical numbers**: `samesex` is a BROAD item (a whole
                category of relationship) and `chsuppor` is a NARROW one (one specific situation).
                Generalisation may track ITEM BREADTH, and have nothing to do with sex.
                **The rival goes first. If it wins, `#915`'s psychology is a measurement artifact.**

Why Now         `#915` is on the front page as the project's only live claim about people. It has
                been there for one round. **The cheapest thing that could destroy it should be run
                before anything is built on top of it.**

Live Worlds     A · **DOMAIN** — generalisation tracks the domain an act belongs to.
                B · **BREADTH** — generalisation tracks how broad the item is; sex is incidental.
                    **This is the unwelcome one, and it is the reason for this round.**
                C · **NEITHER** — the driver is something the decomposition does not name (e.g.
                    these two items are about TEENAGERS, an age effect wearing a breadth costume).
                    (the meta-separator: an outcome that kills my world-DECOMPOSITION itself)

The instrument  2011-2013 female is the ONLY file that ships four items in one battery:
                  `samesex`   "Sexual relations between two same-sex adults is all right"  SEXUAL, BROAD
                  `sxok18`    "All right for unmarried 18 year olds have sex..."           SEXUAL, NARROW
                  `sxok16`    "All right for unmarried 16 year olds have sex..."           SEXUAL, NARROW
                  `chsuppor`  "Okay for unmarried woman to have and raise a child"         FAMILY, NARROW
                Two exposures as in `#915`: `samesexany` (own same-sex contact), `cebow` (own child
                born out of wedlock). n = 5,595 / 5,598 / 5,586 / 5,598 substantive of 5,601.

Estimand        Unchanged from `#915` and one-factor-proof by DERIVATION: under one moral dial
(G1)            b(c -> i) = lambda_i * beta_c, so every 2x2 [exposure x item] coefficient matrix is
                RANK 1 and det(B) == 0 exactly, for any loadings. The round computes det for FOUR
                item pairs, and the worlds differ on which pairs depart:

                  pair                       breadth   domain    A predicts   B predicts
                  (samesex , chsuppor)       differs   differs   != 0         != 0     <- `#915`, uninformative alone
                  (sxok18  , chsuppor)       MATCHED   differs   != 0         ~ 0
                  (sxok16  , chsuppor)       MATCHED   differs   != 0         ~ 0
                  (samesex , sxok18 )        differs   SAME      ~ 0          != 0

                ⚠ **The last three rows are the round.** Rows 2-3 and row 4 carry OPPOSITE
                signatures, so no single outcome can satisfy both worlds. There is no flat row and
                no flat `+` column.

Third rival     ⚠ **EXPOSURE PREVALENCE**, named before the run and not covered by the item design:
                `samesexany` is rare, `cebow > 0` is common within the birth-having subsample. A
                rarer, more extreme exposure marks a more distinctive group, which could produce a
                rank-1 departure with no domain and no breadth involved. CONTROL, same iteration: a
                specification in which `cebow` is dichotomised at `samesexany`'s OWN base rate.

Stopping Rule   One pass over the four pairs x adjustment x support x prevalence-matching grid,
                published whole. Rows 2-4 disagreeing with each other = world C, reported as such
                and NOT adjudicated by picking the pair I prefer.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **SINGLE SITE.** `sxok18`/`sxok16` ship in 2011-2013 female and nowhere else, so this round
    **cannot be replicated across cycles or sexes**. It is one file, and `only this one instrument`
    is a fact about the release rather than a shortcut: no second release asks these four items.
  (2) ⚠ "BROAD" and "NARROW" are MY reading of the item text, not a measured property. The
    assignment is published above so it can be disputed; it is not derived from anything.
  (3) Cross-sectional; the arrow is not identified (inherited from `#915`).
  (4) Value labels do not ship (inherited from `#915`(2)); `det` is invariant to the direction flip.
  (5) `[unchallenged]` — door (3).
"""
import json, re, sys
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
RNG = np.random.default_rng(354)
PAT = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')

DAT, DCT = "2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct"
ITEMS = ["samesex", "sxok18", "sxok16", "chsuppor"]
PLACEBO = ["reactslf", "chbother"]
EXPOS = ["samesexany", "cebow"]
ADJ = ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]
WANT = ITEMS + PLACEBO + EXPOS + ADJ

# (first item, second item, breadth_relation, domain_relation, A_predicts, B_predicts)
PAIRS = [("samesex", "chsuppor", "differs", "differs", "nonzero", "nonzero"),
         ("sxok18", "chsuppor", "MATCHED", "differs", "nonzero", "zero"),
         ("sxok16", "chsuppor", "MATCHED", "differs", "nonzero", "zero"),
         ("samesex", "sxok18", "differs", "SAME", "zero", "nonzero")]

lay = {m.group(2).lower(): (int(m.group(1)) - 1, int(m.group(3)), m.group(4)) for m in
       (PAT.search(l) for l in (D_ / "setup" / DCT).read_text(errors="replace").splitlines()) if m}
have = [v for v in WANT if v in lay]
df = pd.read_fwf(D_ / DAT, colspecs=[(lay[v][0], lay[v][0] + lay[v][1]) for v in have],
                 names=have, dtype=str)
for c in have:
    df[c] = pd.to_numeric(df[c].astype(str).str.strip().replace({"": None}), errors="coerce")

# ══ HARD RULE 1 — n and distribution of every column BEFORE it is cited ═══════════════
inventory = {}
print(f"{DAT}  n={len(df)}  columns present {len(have)}/{len(WANT)}")
for c in have:
    s = df[c]
    inventory[c] = dict(non_missing=int(s.notna().sum()),
                        top={str(k): int(v) for k, v in s.value_counts().head(6).items()},
                        label=lay[c][2])
    print(f"  {c:12s} non-miss {int(s.notna().sum()):5d}/{len(df)}  "
          f"top {dict(list(s.value_counts().head(5).items()))}  {lay[c][2][:52]}")


def norm_clean(s):
    return s.where(s.isin([1, 2, 3, 4, 5]))


ACHIEVED = {}


def expo_clean(name, s, match_rate=None):
    if name == "samesexany":
        return (s == 1).astype(float).where(s.isin([1, 5]))
    v = s.where(s.notna() & (s < 90))
    marked = (v > 0).astype(float).where(v.notna())
    if match_rate is None:
        ACHIEVED["unmatched"] = float(marked.mean())
        return marked
    # ⚠⚠ VERSION 1 OF THIS CONTROL DID NOT DO WHAT ITS NAME SAID, and it nearly retracted `#915`.
    #   It thresholded at `v.quantile(1 - match_rate)`. `cebow` is a DISCRETE COUNT with atoms at
    #   0 (1,166) and 1 (988), so quantile(0.818) = 2.0 and `v > 2` marks **0.0770**, not 0.1820 —
    #   **less than half the target**, and "more than two non-marital births" is a RARER, MORE
    #   EXTREME CONSTRUCT, not a rescaling of the same one. Its smaller `det` was therefore
    #   evidence about a different variable, and reading it as "prevalence explains `#915`" would
    #   have been a FALSE RETRACTION — as permanent as a false acquittal (`#782`).
    #   ⇒ A base rate cannot be hit by thresholding a discrete count. **Keep the construct
    #   (`cebow > 0`) and randomly DROP marked cases until the marked rate equals the target.**
    #   ⚠ AND VERSION 2 GOT THE ARITHMETIC WRONG IN THE OTHER DIRECTION — measured, target 0.1820,
    #   achieved 0.3291 while the UNMATCHED rate was 0.2501, i.e. matching made the rate WORSE.
    #   Dropping a marked case removes it from the NUMERATOR **and the DENOMINATOR**, so
    #   `n_keep = r*N` overshoots. Solving `(M-k)/(N-k) = r` gives `k = (M - r*N)/(1 - r)`.
    #   *Three versions of one control in one round, and the first two both printed a number.*
    idx = marked[marked == 1].index.to_numpy()
    N, M = float(marked.notna().sum()), float(len(idx))
    k = int(round(max((M - match_rate * N) / (1.0 - match_rate), 0.0)))
    if 0 < k <= len(idx):
        drop = RNG.choice(idx, size=k, replace=False)
        marked = marked.copy()
        marked.loc[drop] = np.nan          # dropped, not recoded as 0: they are not controls
    ACHIEVED["prev-matched"] = float(marked.mean())
    return marked


def beta(y, x, covs):
    ok = y.notna() & x.notna()
    for c in covs:
        ok &= c.notna()
    n = int(ok.sum())
    if n < 60:
        return np.nan, n
    Y = y[ok].to_numpy(float)
    X = np.column_stack([np.ones(n), x[ok].to_numpy(float)] + [c[ok].to_numpy(float) for c in covs])
    for j in range(1, X.shape[1]):
        sd = X[:, j].std()
        if sd > 0:
            X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    Y = (Y - Y.mean()) / (Y.std() or 1.0)
    try:
        b, *_ = np.linalg.lstsq(X, Y, rcond=None)
        return float(b[1]), n
    except np.linalg.LinAlgError:
        return np.nan, n


def det_for(frame, i1, i2, adjust, match_rate=None):
    covs = [frame[c] for c in adjust if c in frame.columns]
    B, ns = {}, []
    for e in EXPOS:
        xe = expo_clean(e, frame[e], match_rate)
        for it in (i1, i2):
            b, n = beta(norm_clean(frame[it]), xe, covs)
            B[(e, it)] = b
            ns.append(n)
    if any(np.isnan(v) for v in B.values()):
        return np.nan, {}, 0
    d = B[(EXPOS[0], i1)] * B[(EXPOS[1], i2)] - B[(EXPOS[0], i2)] * B[(EXPOS[1], i1)]
    return float(d), {f"{k[0]}->{k[1]}": round(v, 4) for k, v in B.items()}, int(min(ns))


base_rate = float(expo_clean("samesexany", df["samesexany"]).mean())
print(f"\n  `samesexany` base rate = {base_rate:.4f}  <- the prevalence control matches `cebow` to this")

ADJ_SETS = {"raw": [], "demog": ["age_r", "educat", "hisprace2"],
            "demog+relig": ["age_r", "educat", "hisprace2", "reldlife", "attndnow"]}

grid = []
for i1, i2, br, dom, predA, predB in PAIRS + [(PLACEBO[0], PLACEBO[1], "n/a", "n/a", "zero", "zero")]:
    if i1 not in df.columns or i2 not in df.columns:
        continue
    for aname, aset in ADJ_SETS.items():
        for mname, mrate in (("unmatched", None), ("prev-matched", base_rate)):
            d, B, n = det_for(df, i1, i2, aset, mrate)
            grid.append(dict(pair=f"{i1} vs {i2}", breadth=br, domain=dom,
                             A_predicts=predA, B_predicts=predB, adjust=aname,
                             prevalence=mname, det=d, n=n, coefs=B))

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for g in grid:
    d = "   nan  " if np.isnan(g["det"]) else f"{g['det']:+.5f}"
    print(f"  {g['pair']:22s} breadth={g['breadth']:8s} domain={g['domain']:8s} "
          f"{g['adjust']:12s} {g['prevalence']:12s} det={d}  n={g['n']:5d}")

# ══ NEGATIVE CONTROL — one-factor synthetic at the observed marginals ════════════════
n_syn = int(min(norm_clean(df[i]).notna().sum() for i in ITEMS))
marg = [list(np.linspace(0.15, 0.85, 4)) for _ in range(2)]


def synth(n, g, reps):
    out = []
    for _ in range(reps):
        theta = RNG.standard_normal(n)
        dom = {0: RNG.standard_normal(n), 1: RNG.standard_normal(n)}
        cols = {}
        for j, it in enumerate(("A", "B")):
            z = 0.8 * theta + g * dom[j] + np.sqrt(max(1 - 0.64, 1e-9)) * RNG.standard_normal(n)
            cols[it] = pd.Series(np.clip(np.digitize(z, np.quantile(z, marg[j])) + 1, 1, 5))
        for j, e in enumerate(EXPOS):
            z = 0.55 * theta + g * dom[j] + RNG.standard_normal(n)
            cols[e] = pd.Series(np.where(z > np.quantile(z, 0.75), 1.0, 5.0 if e == "samesexany" else 0.0))
        out.append(pd.DataFrame(cols))
    return out


NUL = [x for x in (det_for(d, "A", "B", [])[0] for d in synth(n_syn, 0.0, 80)) if not np.isnan(x)]
null_med, null_sd = float(np.median(NUL)), float(np.std(NUL))
print(f"\n  one-factor null: median {null_med:+.5f}  sd {null_sd:.5f}  (n={n_syn}, {len(NUL)} reps)")

sweep = []
for g in (0.0, 0.25, 0.5, 0.75, 1.0):
    v = [x for x in (det_for(d, "A", "B", [])[0] for d in synth(n_syn, g, 25)) if not np.isnan(x)]
    sweep.append((g, float(np.median(v)) if v else np.nan))
print(f"  positive sweep (g, median det): {[(g, round(v, 5)) for g, v in sweep]}")

# ══ READ THE TWO SIGNATURES ══════════════════════════════════════════════════════════
def cells(pred_pair):
    return [g for g in grid if g["pair"] == pred_pair and not np.isnan(g["det"])]


def resolved(g):
    return abs(g["det"] - null_med) > 2 * null_sd


matched_cross = [g for g in grid if g["breadth"] == "MATCHED" and not np.isnan(g["det"])]
same_domain = [g for g in grid if g["domain"] == "SAME" and not np.isnan(g["det"])]
plac = [g for g in grid if g["domain"] == "n/a" and not np.isnan(g["det"])]

med_matched = float(np.median([g["det"] for g in matched_cross])) if matched_cross else np.nan
med_same = float(np.median([g["det"] for g in same_domain])) if same_domain else np.nan
med_plac = float(np.median([g["det"] for g in plac])) if plac else np.nan
frac_matched_res = np.mean([resolved(g) for g in matched_cross]) if matched_cross else np.nan
frac_same_res = np.mean([resolved(g) for g in same_domain]) if same_domain else np.nan

print(f"\n  breadth-MATCHED cross-domain (sxok18/sxok16 vs chsuppor): median det {med_matched:+.5f}  "
      f"resolved in {frac_matched_res:.0%} of {len(matched_cross)} cells")
print(f"  SAME-domain breadth-differing (samesex vs sxok18)       : median det {med_same:+.5f}  "
      f"resolved in {frac_same_res:.0%} of {len(same_domain)} cells")
print(f"  placebo (reactslf vs chbother)                          : median det {med_plac:+.5f}")

ps = [2 * (1 - stats.norm.cdf(abs((g["det"] - null_med) / (null_sd or 1e-9))))
      for g in grid if not np.isnan(g["det"])]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

G = Gate("Is `#915`'s asymmetry about the DOMAIN, or about how broad the item is?")
G.plant_direction_from_sweep("positive: planted domain coupling raises det, and g=0 is null",
                             sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("synthetic ONE-FACTOR world at the observed n and marginals",
                   abs(null_med), abs(med_matched) if not np.isnan(med_matched) else 0.0,
                   null_spread=null_sd, null_kind="one-factor latent, matched loadings")
G.multiplicity_control("the whole four-pair grid", ps, 0.05,
                       labels=[f"{g['pair']}|{g['adjust']}|{g['prevalence']}"
                               for g in grid if not np.isnan(g["det"])])
# ⚠⚠ VERSION 1 OF THIS ASSERTION COULD NOT FAIL: it read `np.isnan(med_plac) or ...`, so a placebo
#   that never ran PASSED. It never ran — `chbother` (2,446 non-missing) and `cebow` (3,141, and
#   defined only for respondents with a birth) are **near-disjoint by questionnaire design**, so the
#   intersection is under the n>=60 guard and every placebo cell returned n=0. An absent control
#   scoring as a passed control is the "empty population passes" failure, in my own gate.
#   ⇒ absence now FAILS, and the unavailability is registered instead of being enjoyed.
# ⚠ `#840`: read the prior round's number from ITS ARTIFACT, never transcribe it into this file.
_p915 = (ROOT / "E03_what_an_instrument_would_have_to_be"
         / "A112_is_sexual_morality_a_distinct_object_at_the_person_level"
         / "R353_two_norms_two_exposures_same_respondents" / "results"
         / "two_norms_two_exposures.json")
PLAC_915 = json.loads(_p915.read_text())["det_placebo"] if _p915.exists() else None
G.asserted("placebo actually RAN at this site (absence is not a pass)",
           not np.isnan(med_plac),
           f"placebo det {med_plac} over {len(plac)} cells; `chbother` x `cebow` are near-disjoint "
           f"by design, so this site STRUCTURALLY cannot supply the placebo `#915` had at its two "
           f"2017-2019 sites ({PLAC_915}); scope stated", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("prevalence control ACHIEVED the target rate, not merely attempted it",
           abs(ACHIEVED.get("prev-matched", 0.0) - base_rate) < 0.02,
           f"target {base_rate:.4f} · achieved {ACHIEVED.get('prev-matched', float('nan')):.4f} "
           f"· unmatched {ACHIEVED.get('unmatched', float('nan')):.4f}; v1 thresholded a discrete "
           f"count and achieved 0.0770 while claiming to match 0.1820; scope stated", kind="control")
# ── THE KILL, and it is a CONDITIONAL on two cells with OPPOSITE signatures ──
G.asserted("KILL: world B (breadth) requires the breadth-MATCHED cross-domain pairs to be null",
           not np.isnan(med_matched) and abs(med_matched - null_med) > 2 * null_sd,
           f"breadth-matched median det {med_matched:+.5f} vs null {null_med:+.5f} "
           f"+/- {null_sd:.5f}; resolved in {frac_matched_res:.0%} of cells")

tv = G.three_valued()
matched_nonzero = not np.isnan(med_matched) and abs(med_matched - null_med) > 2 * null_sd
same_nonzero = not np.isnan(med_same) and abs(med_same - null_med) > 2 * null_sd
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif matched_nonzero and not same_nonzero:
    VERDICT, WORLD = "OVERTURNED", "A · DOMAIN (breadth refuted)"
elif same_nonzero and not matched_nonzero:
    VERDICT, WORLD = "CONFIRMED", "B · BREADTH — `#915`'s psychology is an item artifact"
elif matched_nonzero and same_nonzero:
    VERDICT, WORLD = "UNVERIFIED", "C · both signatures fire — the decomposition is wrong"
else:
    VERDICT, WORLD = "UNVERIFIED", "C · neither fires — nothing resolves at this n"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=916, round="E03·A112·R354", verdict=VERDICT, world=WORLD,
           estimand="det of the 2x2 [exposure x item] coefficient matrix, == 0 under any "
                    "one-factor model; computed for four item pairs whose worlds predict opposite",
           site=DAT, inventory=inventory, grid=grid, pairs=[list(p) for p in PAIRS],
           det_breadth_matched=med_matched, det_same_domain=med_same, det_placebo=med_plac,
           frac_resolved_matched=float(frac_matched_res), frac_resolved_same=float(frac_same_res),
           null_median=null_med, null_sd=null_sd, null_reps=len(NUL), n_synth=n_syn,
           positive_sweep=sweep, base_rate=base_rate, achieved_rates=ACHIEVED,
           family_size=len(ps),
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "domain_or_breadth.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'domain_or_breadth.json'}")
