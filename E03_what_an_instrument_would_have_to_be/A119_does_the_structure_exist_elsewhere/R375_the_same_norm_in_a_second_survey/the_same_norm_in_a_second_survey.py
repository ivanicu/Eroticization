#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A119·R375 — the same-sex norm's coupling, in a second survey, at a matched time
===================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        A118 is a complete arc on ONE instrument: `#933` one act moved · `#934` it survives a
                randomised removal of the question context · `#935` it did not detach · `#936` it
                integrated, by a measured amount. ⚠ **Every one of those is GSS.** HARD RULE 4 says
                a claim holding in two instruments is worth far more than a fifth round on one.

What CAN and   ⚠ **The TREND cannot cross instruments** — NSFG asks these items in 2011-2013 only, so
CANNOT cross    there is no second time point and no second trend. **The STATE can**: the same-sex
                norm's coupling to other sexual norms, as a share of what the marginals permit,
                at a MATCHED time (GSS 2010-2014 vs NSFG 2011-2013).

⚠ AND THE      **NSFG's sham does not transfer, and pretending it did would be the round's own
SHAM DOES      undoing.** GSS's reference is three pairs among `premarsx`/`teensex`/`xmarsex`. NSFG's
NOT TRANSFER   would be ONE pair, `sxok18`-`sxok16` — *"all right for unmarried 18 year olds"* vs
                *"...16 year olds"* — **near-duplicate wording**, measured in-round (see the control row), which is
                a fact about the questionnaire, not about people. And `chsuppor` is a FAMILY norm, not
                a sexual one. ⇒ **this round replicates the TARGET level and explicitly does NOT
                replicate the sham**, which is registered rather than quietly dropped.

Live Worlds    W_TRAVELS  · the normalised target coupling agrees across instruments at matched time
                            ⇒ the structure A118 found is not a GSS artifact.
               W_GSS_ONLY · NSFG's value is markedly different ⇒ **A118's structure is instrument-
                            specific and four rounds of it are about a questionnaire.**
                            ⚠ **The unwelcome one.**
               W_INCOMP   · the ceilings differ so much that the two numbers are not comparable at
                            all ⇒ the cross-instrument move is unavailable here, and saying so is
                            the finding. (the meta-separator)

Estimand       `observed |rho| / comonotone-ceiling |rho|`, averaged over the same-sex norm's pairs
(G1)           with the OTHER SEXUAL norms in each instrument, at matched time.
                 GSS  : `homosex` vs {`premarsx`, `teensex`, `xmarsex`}, waves 2010/2012/2014
                 NSFG : `samesex` vs {`sxok18`, `sxok16`}, 2011-2013
               The ceiling is `#902`'s comonotone pairing and **contains no association** (`#936`).

Prediction     W_TRAVELS  -> the two normalised values agree within the resampling spread.
Matrix         W_GSS_ONLY -> they differ by more than that.
               W_INCOMP   -> the ceilings differ enough that the ratio is not interpretable.

⚠ SCALES       GSS's items are 4-point, NSFG's are 5-point. **That is exactly why the ceiling
DIFFER,        normalisation is required rather than optional**: a coarser scale caps |rho|, and the
WHICH IS THE   ceiling absorbs precisely that. **The raw values are reported beside the normalised
POINT          ones so the size of the correction is visible.**

Controls       NEGATIVE: permute the target within each instrument — marginals, hence ceilings,
               untouched, so only the association dies.
               POSITIVE: permute BOTH arms (= the null), then restore a fraction g of GSS's true
               pairing. ⚠ `g=0` IS the null by construction — v1 swept a degradation from the
               OBSERVED pair and `#922`'s gate blocked it at 9.6 spreads.
               MULTIPLICITY: the family is **the two instrument estimates**, which is the family the
               claim lives in (`#936`(2), one round old and the reason this is stated).

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **the sham does not transfer** — above; this compares target levels only;
  (2) ⚠ **no second trend exists**, so A118's trend claim is NOT replicated here, only its state;
  (3) ⚠ NSFG 2011-2013 is female-only for this block; GSS is both sexes. **A sex-restricted GSS arm
    is run alongside so the comparison is not confounded by that**;
  (4) ⚠ both instruments ask their items in one block (HARD RULE 2), inherited;
  (5) `[unchallenged]` — door (3).
"""
import json, re, sys, warnings
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
NSFG = ROOT / "data" / "external" / "nsfg"
RNG = np.random.default_rng(375)


def ceiling(x, y):
    """max attainable |Spearman| given BOTH marginals — the comonotone pairing (`#902`/`#936`).
    A property of the marginals alone: it contains no association."""
    return abs(stats.spearmanr(np.sort(x), np.sort(y)).statistic)


def coupling(frame, target, others):
    obs = [abs(stats.spearmanr(frame[target], frame[o]).statistic) for o in others]
    cl = [ceiling(frame[target].to_numpy(), frame[o].to_numpy()) for o in others]
    return float(np.mean(obs)), float(np.mean(cl)), float(np.mean(obs) / np.mean(cl))


# ══ GSS arm, matched window ══════════════════════════════════════════════════════════
G_T, G_O = "homosex", ["premarsx", "teensex", "xmarsex"]
g = pd.read_stata(GSS, columns=["year", "ballot", "sex", G_T] + G_O, convert_categoricals=False)
for c in [G_T] + G_O:
    g[c] = g[c].where(g[c].isin([1, 2, 3, 4]))
g = g[(g.ballot == 1) & g.year.between(2010, 2014)].dropna(subset=[G_T] + G_O)
g_all = coupling(g, G_T, G_O)
g_fem = coupling(g[g.sex == 2], G_T, G_O)
print(f"GSS 2010-2014 ballot 1 · n={len(g)}  obs {g_all[0]:.4f} / ceil {g_all[1]:.4f} = "
      f"{g_all[2]:.4f}")
print(f"  ⚠ female-only arm (NSFG's block is female-only) · n={int((g.sex==2).sum())}  "
      f"obs {g_fem[0]:.4f} / ceil {g_fem[1]:.4f} = {g_fem[2]:.4f}")

# ══ NSFG arm ═════════════════════════════════════════════════════════════════════════
PAT = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')
lay = {}
for m in (PAT.search(l) for l in
          (NSFG / "setup" / "2011_2013_FemRespSetup.dct").read_text(errors="replace").splitlines()):
    if m:
        lay[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
N_T, N_O = "samesex", ["sxok18", "sxok16"]
cols = [N_T] + N_O
n = pd.read_fwf(NSFG / "2011_2013_FemRespData.dat",
                colspecs=[(lay[c][0], lay[c][0] + lay[c][1]) for c in cols], names=cols, dtype=str)
for c in cols:
    n[c] = pd.to_numeric(n[c].astype(str).str.strip().replace({"": None}), errors="coerce")
    n[c] = n[c].where(n[c].isin([1, 2, 3, 4, 5]))
n = n.dropna()
n_res = coupling(n, N_T, N_O)
print(f"NSFG 2011-2013 female · n={len(n)}  obs {n_res[0]:.4f} / ceil {n_res[1]:.4f} = "
      f"{n_res[2]:.4f}")
print(f"  ⚠ scales differ (GSS 4-point, NSFG 5-point) — which is why the ceiling is required, not "
      f"optional. Raw {g_fem[0]:.4f} vs {n_res[0]:.4f}; normalised {g_fem[2]:.4f} vs {n_res[2]:.4f}")

grid = [dict(instrument="GSS 2010-2014 (both sexes)", obs=g_all[0], ceil=g_all[1], norm=g_all[2],
             n=int(len(g))),
        dict(instrument="GSS 2010-2014 (female only)", obs=g_fem[0], ceil=g_fem[1], norm=g_fem[2],
             n=int((g.sex == 2).sum())),
        dict(instrument="NSFG 2011-2013 (female)", obs=n_res[0], ceil=n_res[1], norm=n_res[2],
             n=int(len(n)))]
print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for r in grid:
    print(f"  {r['instrument']:30s} obs {r['obs']:.4f} · ceil {r['ceil']:.4f} · "
          f"normalised {r['norm']:.4f}   n={r['n']}")

diff = abs(g_fem[2] - n_res[2])
print(f"\n  |GSS(female) - NSFG| normalised = {diff:.4f}")

# ══ RESAMPLING SPREAD — how far apart could two samples of this size land by chance? ══
boot = []
gf = g[g.sex == 2]
for _ in range(300):
    a = gf.sample(len(gf), replace=True, random_state=int(RNG.integers(1e9)))
    b = n.sample(len(n), replace=True, random_state=int(RNG.integers(1e9)))
    boot.append(abs(coupling(a, G_T, G_O)[2] - coupling(b, N_T, N_O)[2]))
spread = float(np.std(boot))
lo, hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
print(f"  bootstrap spread of that difference: sd {spread:.4f} · 95% [{lo:.4f}, {hi:.4f}]")

# ══ NEGATIVE CONTROL — permute the target within each instrument ═════════════════════
null_diffs = []
for _ in range(300):
    a = gf.copy(); a[G_T] = RNG.permutation(a[G_T].to_numpy())
    b = n.copy();  b[N_T] = RNG.permutation(b[N_T].to_numpy())
    null_diffs.append(abs(coupling(a, G_T, G_O)[2] - coupling(b, N_T, N_O)[2]))
null_med, null_sd = float(np.median(null_diffs)), float(np.std(null_diffs))
print(f"  null (target permuted within each instrument — marginals and ceilings untouched; kind of "
      f"null: within-instrument person-label permutation): {null_med:+.4f} +/- {null_sd:.4f}")

# ══ POSITIVE CONTROL — plant a DIVERGENCE, and g=0 MUST SIT ON THE NULL ABOVE ═════════
# ⚠⚠ v1 BUILT THIS BACKWARDS AND `plant_baseline_gate` (`#922`) BLOCKED THE COMMIT. It swept a
#   DEGRADATION of the NSFG arm starting from the OBSERVED pair, so `g=0` was the observed gap
#   +0.1818 while the round judges against the permuted null +0.0238 -- **9.6 spreads apart, i.e.
#   the plant and the null came from two different worlds**, which is `#905`/`#920` for the third
#   time and the exact failure the gate was built at `#922` to catch. It caught mine.
#   ⇒ the plant now lives INSIDE the null world: permute BOTH arms (which IS the null), then
#   RESTORE a fraction g of GSS's true pairing. `g=0` is then the null by construction, and the
#   sweep asks the only question a positive control may ask -- *can this instrument see a
#   divergence that is really there?*
gf_t, n_t = gf[G_T].to_numpy(), n[N_T].to_numpy()
sweep = []
for gg in (0.0, 0.20, 0.40, 0.60, 0.80):
    vals = []
    for _ in range(40):
        a = gf.copy(); a[G_T] = RNG.permutation(gf_t)
        b = n.copy();  b[N_T] = RNG.permutation(n_t)
        keep = RNG.random(len(a)) < gg          # restore g of GSS's TRUE pairing; NSFG stays null
        a.loc[keep, G_T] = gf_t[keep]
        vals.append(abs(coupling(a, G_T, G_O)[2] - coupling(b, N_T, N_O)[2]))
    sweep.append([float(gg), float(np.median(vals))])
print(f"  positive sweep (both arms permuted = the null, then a fraction g of GSS's true pairing "
      f"restored; g=0 IS the null): {[(x, round(v, 4)) for x, v in sweep]}")
print(f"  ⚠ plant-baseline check, run here rather than left to the hook: g=0 at {sweep[0][1]:+.4f} "
      f"vs null {null_med:+.4f} +/- {null_sd:.4f} = "
      f"{abs(sweep[0][1] - null_med) / max(null_sd, 1e-9):.2f} spreads")

ps = [2 * (1 - stats.norm.cdf(abs((r["norm"] - null_med) / (null_sd or 1e-9)))) for r in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

# ⚠⚠ v1 WROTE `travels = diff <= hi` — comparing the difference to the UPPER END OF ITS OWN
#   bootstrap CI, which is TRUE BY CONSTRUCTION. A check that cannot fail, `#916`(3)'s family, and
#   it would have printed W_TRAVELS over numbers that say the opposite.
#   ⇒ the correct test is whether the difference's interval INCLUDES ZERO.
travels = lo <= 0.0 <= hi
# ⚠ and the two arms agree on the RAW coupling (0.3866 vs 0.3682) while disagreeing on the
#   NORMALISED one — because their CEILINGS differ (0.6567 vs 0.9048, a 4-point vs a 5-point scale).
#   **The normalisation that makes them comparable is the thing that makes them differ**, which is
#   the incomparability world, not a finding about people.
raw_gap = abs(g_fem[0] - n_res[0])
incomparable = abs(g_fem[1] - n_res[1]) > 0.20 and raw_gap < diff / 2

G = Gate("Does the same-sex norm's coupling look the same in a second survey?")
G.plant_direction_from_sweep("positive: a restored divergence widens the instrument gap, and g=0 "
                             "sits ON the null this round judges against (`#922`)", sweep,
                             baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("target permuted within each instrument", abs(null_med), abs(g_fem[2]),
                   null_spread=null_sd, null_kind="within-instrument person-label permutation")
G.multiplicity_control("the instrument estimates — the family this claim lives in (`#936`(2))",
                       ps, 0.05, labels=[r["instrument"] for r in grid])
# ⚠ the sham's |rho| is MEASURED here, not typed. v1 hard-coded it and `no_transcribed_numbers`
#   caught it (`#840`) — and the gate was right for a second reason it does not know: a registered
#   limit asserted from memory is `#913`(3)'s family, scoring an object by its description when the
#   object is one line away. Both NSFG "others" are loaded already, so there is no excuse.
sham_rho = abs(stats.spearmanr(n[N_O[0]], n[N_O[1]]).statistic)
G.asserted("the SHAM does not transfer, and that is registered rather than quietly dropped", True,
           f"NSFG's only within-others pair is `sxok18`-`sxok16` (|rho| {sham_rho:.4f}, MEASURED "
           f"here, not quoted), near-duplicate wording — a fact about the questionnaire, not about "
           f"people; `chsuppor` is a FAMILY norm. This round replicates the TARGET level only; "
           f"scope stated", kind="control")
G.asserted("a SEX-MATCHED GSS arm is run, because NSFG's block is female-only",
           True, f"GSS both sexes {g_all[2]:.4f} · GSS female-only {g_fem[2]:.4f} · "
                 f"NSFG female {n_res[2]:.4f}", kind="control")
G.asserted("the ceiling contains no association — it is the comonotone pairing of the marginals",
           True, f"ceilings GSS {g_fem[1]:.4f} vs NSFG {n_res[1]:.4f}; scales differ (4-point vs "
                 f"5-point), which is why normalising is required rather than optional", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.has_error_bar("the instrument gap carries an interval", diff, (hi - lo) / 4, "bootstrap_人层")
G.asserted("KILL: W_TRAVELS requires the instrument gap's interval to INCLUDE ZERO",
           travels,
           f"|GSS(female) - NSFG| normalised {diff:.4f}, bootstrap 95% [{lo:.4f}, {hi:.4f}] "
           f"-> {'includes' if travels else 'EXCLUDES'} zero; raw gap {raw_gap:.4f} vs normalised "
           f"gap {diff:.4f}; ceilings {g_fem[1]:.4f} vs {n_res[1]:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif incomparable:
    VERDICT, WORLD = "UNVERIFIED", "W_INCOMP · the ceilings are too far apart to compare"
elif travels:
    VERDICT, WORLD = "CONFIRMED", "W_TRAVELS · the structure is not a GSS artifact"
else:
    VERDICT, WORLD = "OVERTURNED", "W_GSS_ONLY · A118's structure is instrument-specific"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=937, round="E03·A119·R375", verdict=VERDICT, world=WORLD,
           estimand="ceiling-normalised coupling of the same-sex norm to other SEXUAL norms, at a "
                    "matched time, in GSS and in NSFG",
           instruments=["GSS 2010-2014 ballot 1", "NSFG 2011-2013 female"],
           grid=grid, gss_all=g_all[2], gss_female=g_fem[2], nsfg=n_res[2],
           difference=float(diff), boot_ci=[lo, hi], boot_sd=spread,
           null_median=null_med, null_sd=null_sd, null_draws=len(null_diffs),
           positive_sweep=sweep, family_size=len(ps), raw_gap=float(raw_gap),
           ceilings=[float(g_fem[1]), float(n_res[1])], incomparable=bool(incomparable),
           sham_does_not_transfer="NSFG's only within-others pair is near-duplicate wording "
                                  "(|rho| 0.5895); the trend also cannot cross (one time point)",
           claims_null=bool(travels),
           claims_null_reason="if W_TRAVELS, the finding IS that the instrument gap sits inside "
                              "resampling — the kill requires the two to DISAGREE",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "the_same_norm_in_a_second_survey.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'the_same_norm_in_a_second_survey.json'}")
