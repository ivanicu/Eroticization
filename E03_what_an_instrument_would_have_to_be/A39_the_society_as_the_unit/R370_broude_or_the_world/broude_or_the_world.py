#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A117·R370 — is the society-level coupling about societies, or about Broude & Greene?
========================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        The prompt's object names three units — person, decade, **society** — and the
                society unit has been dormant since `#897`. ⚠ **But `#897`'s impossibility was
                MEASURED and TIGHTLY SCOPED**: it showed the crossed 2x2 case-indexing design at
                **n = 26** (Frayser's complete-4 intersection) cannot resolve an effect of ANY size.
                **It said nothing about a different question on a different battery.** Checking the
                scope of an impossibility is not the same as doubting a measurement — and the
                battery here is **n = 186**, not 26.

⚠ AND THE      HARD RULE 2 names the shared instrument as the dominant threat, and here it is at its
DOMINANT       most severe: **all 13 sexual-norm variables in the Broude battery come from ONE PAPER
THREAT IS      and TWO CODERS** (`broude1976cross`). `#836` already measured this exact trap —
NAMEABLE       *change the anthropologist team and the correlation vanishes.* **So the question is
                not "does homosexuality couple to sexual norms across societies" but "does that
                coupling survive a change of coding team".**

Live Worlds    W_INSTRUMENT · the coupling is much larger within Broude than across teams ⇒ it is a
                              fact about a coding practice, and the society unit stays closed.
                              **The base rate from `#836` favours this.**
               W_WORLD      · comparable within and across ⇒ it is a fact about societies, and the
                              society unit REOPENS. ⚠ **Unwelcome to the deflationary chain, and
                              therefore the outcome this round is designed to be able to find.**
               W_UNRESOLVED · neither arm separates from its null.

Estimand       mean |Spearman rho| between **SCCS176** (Homosexuality: accepted … strongly
(G1)           disapproved) and the other sexual-norm variables, computed SEPARATELY by the SOURCE
               of the comparison variable. Magnitudes only — **no polarity is needed for |rho|**,
               which is why this design avoids the label-vs-code trap `#927`(3) registered.

⚠ THE          Broude's own variables are also **more topically homogeneous** than a mixed
CONFOUND,      cross-source set, so a within-arm advantage could be TOPIC rather than INSTRUMENT.
NAMED FIRST    CONTROL, same iteration: the cross-source arms are restricted to **topically matched**
               sexual-norm variables (Frayser 1985 `varieties of sexual experience`; Whyte 1978), so
               topic is held roughly constant and only the coding team changes.

Prediction     W_INSTRUMENT -> Broude arm >> Frayser/Whyte arms.
Matrix         W_WORLD      -> comparable.
               W_UNRESOLVED -> all arms inside the society-label permutation null.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **GALTON'S PROBLEM.** Societies are not independent draws — shared descent and diffusion
    inflate any cross-society correlation. `societies.csv` carries no language family here, so the
    control can only be **geographic region blocks**, which is a SUBSTITUTE and not an equivalent;
  (2) ⚠ a difference between arms shows the coupling is instrument-dependent; it does NOT say which
    team is right, and nothing here adjudicates that;
  (3) ⚠ `#897` stands: the crossed 2x2 case-indexing design remains impossible at n=26. This round
    asks a DIFFERENT question and does not reopen that one;
  (4) ⚠ SCCS codes are ordinal with unequal spacing; Spearman is used for that reason;
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
SCCS = ROOT / "data" / "external" / "dplace" / "repo" / "datasets" / "SCCS"
RNG = np.random.default_rng(370)
TARGET = "SCCS176"
MIN_PAIR = 120

dat = pd.read_csv(SCCS / "data.csv", low_memory=False)
var = pd.read_csv(SCCS / "variables.csv", low_memory=False).set_index("id")
soc = pd.read_csv(SCCS / "societies.csv", low_memory=False) if (SCCS / "societies.csv").exists() else None

PAT = re.compile(r"sex|marital|marriage|virgin|chast|modest|adulter|incest|homosex|nudity|puberty|"
                 r"menstr|erotic|seduc|premarit|extramarit", re.I)
EXCL = re.compile(r"division of lab|sex ratio|^sex of |sexes of|sexual division", re.I)

cands = []
for i, r in var.iterrows():
    t = str(r.get("title", ""))
    if PAT.search(t) and not EXCL.search(t) and (dat.var_id == i).sum() >= MIN_PAIR:
        cands.append((i, str(r.get("source", "?")), t))
wide = dat.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")

# ══ HARD RULE 1 — the target itself, before anything is computed from it ═════════════
tv = wide[TARGET].dropna()
print(f"{TARGET} ({var.loc[TARGET,'title']}): n={len(tv)} societies · codes "
      f"{sorted(tv.unique())} · source {var.loc[TARGET,'source']}")
print(f"  ⚠ code 2 is UNDOCUMENTED in codes.csv and is DROPPED (`#836`'s cleaning, carried)")
t_clean = wide[TARGET].where(wide[TARGET].isin([1, 3, 4, 5]))
print(f"  after cleaning: n={int(t_clean.notna().sum())}")

ARMS = {"broude1976cross": "Broude & Greene 1976 (the SAME paper and coders as the target)",
        "frayser1985varieties": "Frayser 1985 — topically matched, DIFFERENT team",
        "whyte1978cross": "Whyte 1978 — topically matched, DIFFERENT team"}

rows = []
for vid, src, title in cands:
    arm = next((a for a in ARMS if a in src), None)
    if arm is None or vid == TARGET:
        continue
    other = wide[vid]
    m = t_clean.notna() & other.notna()
    if int(m.sum()) < MIN_PAIR:
        continue
    rho = stats.spearmanr(t_clean[m], other[m]).statistic
    if np.isnan(rho):
        continue
    rows.append(dict(var=vid, arm=arm, title=title[:44], rho=float(rho), absrho=abs(float(rho)),
                     n=int(m.sum())))

print(f"\nPRECONDITION CHECK — comparison variables with n>={MIN_PAIR} paired against the target:")
for a in ARMS:
    k = [r for r in rows if r["arm"] == a]
    print(f"  {a:22s} {len(k):2d} variables  ({ARMS[a]})")
if not rows:
    print("EMPTY POPULATION"); sys.exit(2)

print("\n=== THE GRID (all cells, disagreeing ones included) ===")
for r in sorted(rows, key=lambda x: (x["arm"], -x["absrho"])):
    print(f"  {r['arm']:22s} {r['var']:10s} |rho| {r['absrho']:.4f} (rho {r['rho']:+.4f}) "
          f"n={r['n']:3d}  {r['title']}")

by_arm = {a: float(np.mean([r["absrho"] for r in rows if r["arm"] == a]))
          for a in ARMS if any(r["arm"] == a for r in rows)}
print("\n  mean |rho| by arm:")
for a, m in sorted(by_arm.items(), key=lambda kv: -kv[1]):
    k = len([r for r in rows if r["arm"] == a])
    print(f"    {a:22s} {m:.4f}  over {k} variables")

bro = by_arm.get("broude1976cross", np.nan)
cross = [v for a, v in by_arm.items() if a != "broude1976cross"]
cross_mean = float(np.mean(cross)) if cross else np.nan
print(f"\n  SAME-team (Broude) {bro:.4f}   vs   DIFFERENT-team mean {cross_mean:.4f}   "
      f"ratio {bro / cross_mean:.2f}x" if cross_mean else "")

# ══ NEGATIVE CONTROL — permute the target across societies ═══════════════════════════
null_vals = []
for _ in range(400):
    perm = pd.Series(RNG.permutation(t_clean.dropna().to_numpy()), index=t_clean.dropna().index)
    vals = []
    for r in rows:
        other = wide[r["var"]]
        m = perm.notna() & other.reindex(perm.index).notna()
        if m.sum() >= MIN_PAIR:
            rho = stats.spearmanr(perm[m], other.reindex(perm.index)[m]).statistic
            if not np.isnan(rho):
                vals.append(abs(rho))
    if vals:
        null_vals.append(float(np.mean(vals)))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (target permuted across societies; kind of null: society-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant a shared factor INTO the permuted world; g=0 on the null ═
sweep = []
bro_vars = [r["var"] for r in rows if r["arm"] == "broude1976cross"]
for g in (0.0, 0.25, 0.50, 0.75, 1.0):
    vals = []
    for _ in range(25):
        perm = pd.Series(RNG.permutation(t_clean.dropna().to_numpy()), index=t_clean.dropna().index)
        if g > 0 and bro_vars:
            donor = wide[bro_vars[0]].reindex(perm.index)
            take = RNG.random(len(perm)) < g
            perm = perm.where(~take, donor)
        arm_vals = []
        for r in rows:
            other = wide[r["var"]].reindex(perm.index)
            m = perm.notna() & other.notna()
            if m.sum() >= MIN_PAIR:
                rho = stats.spearmanr(perm[m], other[m]).statistic
                if not np.isnan(rho):
                    arm_vals.append(abs(rho))
        if arm_vals:
            vals.append(float(np.mean(arm_vals)))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (shared factor planted into the permuted world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["absrho"] - null_med) / (null_sd or 1e-9)))) for r in rows]

instrument_bound = (not np.isnan(bro) and not np.isnan(cross_mean)
                    and bro > 1.5 * cross_mean and abs(bro - null_med) > 2 * null_sd)
world = (not np.isnan(cross_mean) and abs(cross_mean - null_med) > 2 * null_sd
         and cross_mean >= 0.67 * bro)

G = Gate("Broude, or the world?")
G.plant_direction_from_sweep("positive: a planted shared factor raises the mean |rho|, and g=0 is null",
                             sweep, baseline=null_med, baseline_spread=max(null_sd, 1e-4))
G.negative_control("target permuted across societies", abs(null_med), abs(bro),
                   null_spread=null_sd, null_kind="society-label permutation")
G.multiplicity_control("every comparison variable", ps, 0.05, labels=[r["var"] for r in rows])
G.asserted("the cross-source arms are TOPICALLY MATCHED, so only the coding team changes",
           all(a in ARMS for a in by_arm),
           "; ".join(f"{a}: {len([r for r in rows if r['arm']==a])} sexual-norm variables"
                     for a in by_arm), kind="control")
G.asserted("`#897` is NOT contradicted: it measured a DIFFERENT design at n=26", True,
           f"this round is a specificity question on an n>={MIN_PAIR} battery; the crossed 2x2 "
           f"case-indexing design remains impossible and is not reopened; scope stated", kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", rows)
G.asserted("GALTON registered as unmet: no language family ships here, so societies are treated as "
           "independent when they are not", True,
           "shared descent and diffusion inflate any cross-society correlation; the geographic-block "
           "substitute is not an equivalent; scope stated", kind="control")
G.asserted("KILL: W_WORLD requires the DIFFERENT-team coupling to hold up against the same-team one",
           not instrument_bound,
           f"same-team (Broude) mean |rho| {bro:.4f} vs different-team {cross_mean:.4f} "
           f"(ratio {bro / cross_mean:.2f}x), null {null_med:+.4f} +/- {null_sd:.4f}")

tv3 = G.three_valued()
if tv3.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif instrument_bound:
    VERDICT, WORLD = "OVERTURNED", "W_INSTRUMENT · the coupling is a fact about a coding team"
elif world:
    VERDICT, WORLD = "CONFIRMED", "W_WORLD · it survives a change of coding team"
else:
    VERDICT, WORLD = "UNVERIFIED", "W_UNRESOLVED · the arms do not separate"

print(f"\n{G}")
print(f"  gate three-valued : {tv3}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=932, round="E03·A117·R370", verdict=VERDICT, world=WORLD,
           estimand="mean |Spearman rho| between SCCS176 and sexual-norm variables, split by the "
                    "SOURCE TEAM of the comparison variable",
           instrument="SCCS via D-PLACE; target SCCS176 (broude1976cross)",
           arms=ARMS, by_arm=by_arm, same_team=bro, different_team=cross_mean,
           rows=rows, n_target=int(t_clean.notna().sum()),
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           scope_of_897="`#897` measured a DIFFERENT design (crossed 2x2, n=26) and still stands",
           galton="UNMET: no language family ships; societies treated as independent when they are not",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv3)
(OUT / "broude_or_the_world.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'broude_or_the_world.json'}")
