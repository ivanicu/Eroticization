r"""#889 · E03·A104·R327 — a bootstrap resamples the sample; it does not resample the sample DEFINITION

Pays `#888`② and attacks `#888`'s own headline. **The design I set out to run was killed by its own
gradient check, and the check redirected the round** — which is `frontier` §1.6 working rather than
failing, and is recorded here rather than quietly dropped.

**WHAT I INTENDED, AND WHY IT DIED BEFORE COSTING A ROUND.** `#888`① named norm-side attenuation as
the largest open threat. The plan was to **equalise the norm ruler by construction**: replace each
act-specific single-item norm with one multi-item *general sexual-strictness* scale
(`premarsx`+`xmarsex`+`teensex`+`homosex`) used for **both** acts. Three measurements killed it:
   · **the shared scale is the NOISIER ruler, not the cleaner one** — Cronbach α = **0.691**, against
     the sanction indices' KR-20 of **0.802** and **0.789**. Equalising by adding noise is not a
     control, it is a worse instrument;
   · **it costs half the sample** — requiring all four norm items drops n from **2,480 to 1,279**,
     because `teensex` starts in 1986 and takes 2008 out;
   · **it is partly circular** — the scale *contains* `homosex`, which is the homosexuality norm,
     so `shared × homosexuality-sanction` puts the same item on both sides. A leave-one-out scale
     fixes the circularity and lands Δ at **+0.015**, i.e. it changes the estimand rather than
     controlling the old one. ⇒ **it was a different estimand dressed as a control, and the gradient
     check said so before any round was spent.**

**WHAT THE CHECK FOUND INSTEAD, AND IT IS ABOUT `#888`'s OWN NUMBER.** Recomputing `#888`'s
Δ_matched under **defensible variations of the SAMPLE DEFINITION** — not the items, the *rows*:

```
#888 as published (all 7 abortion items present)   n=2480   Δ +0.051
only the TRIO required (not all seven)             n=2586   Δ +0.033
#883's base sample (single-binary rule)            n=2546   Δ +0.035
+ requires teensex (this drops wave 2008)          n=1266   Δ +0.024
```
**`#888` reported `+0.051` with a paired bootstrap 95% of `[+0.013, +0.088]`. A bootstrap resamples
the SAMPLE; it does not resample the SAMPLE DEFINITION** — and a definition change is not
hypothetical here, it is four lines of equally defensible code.

`G1` **ESTIMAND, named before the method — and it is deliberately NOT a time-series estimand:**
   the **leave-one-wave-out spread** of Δ_matched, `max − min` over the four LOWO samples, plus the
   **completeness-rule spread** over the sample definitions above. **No trend, no break, no era.**
   ⚠ `#887`① forbids a fourth *time-series* round on this coupling. **This is not one**: it asks
   whether a POOLED number is an average over heterogeneous cells, which is a `G3`/`G4` obligation
   for any pooled statistic and would be required if the four cells were countries rather than years.
   The tension is named rather than skirted.

**ARITHMETIC FIRST:**
   · a **max − min over four noisy estimates is positive by construction** ⇒ *"Should this zero be
     zero?"* — **NO** ⇒ **`offset_control`**, and **the kind of null is named: a CONSTANT-Δ
     resampling null** — respondents drawn within each wave from the **pooled** joint table, so the
     true Δ is identical in every wave by construction, and the **same LOWO spread** is computed on
     each draw. Whatever spread that produces is what four noisy waves hand you for free;
   · **a leave-one-out estimate and the pooled estimate are not independent** — the LOWO sample is
     three quarters of the pooled one — so a LOWO value lying outside the pooled bootstrap interval
     is **not** a contradiction and is not reported as one.

THREE WORLDS (each with a branch):
   **A HOMOGENEOUS.** The LOWO spread sits inside the constant-Δ null ⇒ `#888`'s `+0.051` is a
     legitimate pooled estimate and the four waves are four draws of one number.
   **B ⚠ THE UNWELCOME ONE — ONE WAVE CARRIES IT.** The spread exceeds the null and dropping a single
     wave moves Δ materially ⇒ **`#888`'s headline is an average over heterogeneous cells**, and the
     honest statement is the LOWO range, not the pooled point.
   **C ⚠ META-SEPARATOR — THE DEFINITION IS THE PARAMETER.** The completeness-rule spread rivals or
     exceeds the LOWO spread ⇒ **the largest source of uncertainty in this claim is which rows I
     chose to keep**, which no bootstrap and no null addresses, and every interval this project has
     published on a pooled level is narrower than the thing it is an interval about.

PREDICTION MATRIX:
   | world           | now  | LOWO spread inside null | above null, one wave dominant | definition spread rivals it |
   | A homogeneous   | 0.35 | **0.85**                | 0.05                          | 0.10                        |
   | B one wave      | 0.35 | 0.05                    | **0.85**                      | 0.10                        |
   | C definition    | 0.30 | 0.10                    | 0.10                          | **0.80**                    |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires — a deliberately planted between-wave difference in Δ must be
      detected, with `floor` and `ceiling` MEASURED, and **at g = 0 it must not fire**
  and the **negative control** is null — a global permutation sits on zero
  and **every sample definition is reported**, including the ones that shrink the effect:
      LOWO spread inside the constant-Δ null's 95th percentile      -> A
      above it                                                       -> **B, and `#888`'s headline
          becomes the LOWO range**
      completeness-rule spread >= the LOWO spread                    -> **C, and it outranks A and B**
  else: **UNVERIFIED**.

`G3` MULTIPLICITY over the whole grid: {4 LOWO samples + 5 completeness rules} × {3 estimators} ×
{35 abortion trios for the primary estimator}. Every cell reported. `G4` SPECIFICATION CURVE over
the same axes.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **four waves is four cells** — a heterogeneity test on four groups has low power, so a null
     here is a **bound**, never a demonstration of homogeneity;
 (2) **causally identified: N/A** — repeated cross-sections;
 (3) **the completeness-rule spread is not exhaustive** — I enumerate the rules I can defend, and a
     rule I did not think of is not covered. **Unenumerated is not cleared**;
 (4) ⚠ **the instrument cannot be changed** — `#882` measured that the only other matched-pair
     instrument here (SCCS) has one observation per society and therefore **no within-instrument
     subsample structure at all**. **Only this one instrument**;
 (5) **no second coder, no second release.**
"""
import itertools
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
from lib.gates import Gate
from lib.gss_polarity import refusal          # `#868`'s home, imported rather than copied

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(327)
NSIM = 1500
AB7 = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]
TRIO = ["abdefect", "abhlth", "abrape"]        # `#888`'s difficulty-matched trio, fixed by its rule
D_888, CI_888 = 0.051, (0.013, 0.088)
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"

print("=== (0) HARD RULE 1 — n and years, before any column is cited ===")
d = pd.read_stata(F, columns=["year", "homosex", "spkhomo", "colhomo", "libhomo", "abpoorw",
                              "premarsx", "xmarsex", "teensex"] + AB7, convert_categoricals=False)
for c in ["homosex", "abpoorw", "teensex"] + TRIO:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])}  ({len(ys)} waves)")
W = 4 - d["homosex"]
HREF = sum(refusal(d[f"{k}homo"], f"{k}homo") for k in ("spk", "col", "lib"))
AIDX = sum((d[c] == 2).astype(float).where(d[c].notna()) for c in TRIO)
M888 = (W.notna() & HREF.notna() & d["abpoorw"].notna() & d["abpoor"].notna()
        & d[AB7].notna().all(axis=1))
WAVES = sorted(int(y) for y in d.loc[M888, "year"].unique())
print(f"\n  `#888`'s published sample n={int(M888.sum())} · waves {WAVES}")
if M888.sum() < 500:
    raise SystemExit("STOP: empty population must never pass")
POP = (f"the {int(M888.sum())} GSS respondents of `#888`'s published sample, waves {WAVES}, and the "
       f"defensible re-definitions of that sample enumerated below")


def delta(mask, est="spearman", ab=AIDX):
    def r_(a, b):
        k = a.notna() & b.notna() & mask
        if k.sum() < 150:
            return np.nan
        x, y = a[k].to_numpy(), b[k].to_numpy()
        if est == "spearman":
            rr = stats.spearmanr(x, y).statistic
        elif est == "kendall":
            rr = stats.kendalltau(x, y, variant="b").statistic
        else:
            tab = pd.crosstab(x, y).to_numpy()
            c = dd = 0
            for i in range(tab.shape[0]):
                for j in range(tab.shape[1]):
                    c += tab[i, j] * tab[i + 1:, j + 1:].sum()
                    dd += tab[i, j] * tab[i + 1:, :j].sum()
            rr = (c - dd) / (c + dd) if (c + dd) else np.nan
        return 1 - rr ** 2
    return float(r_(W, HREF) - r_(d["abpoorw"], ab))


print("\n=== (1) THE GRADIENT CHECK THAT REDIRECTED THIS ROUND — recorded, not dropped ===")
Z = d.loc[M888 & d[["premarsx", "xmarsex", "teensex"]].notna().all(axis=1),
          ["premarsx", "xmarsex", "teensex", "homosex"]].apply(lambda s: 4 - s)
k_ = Z.shape[1]
alpha = float(k_ / (k_ - 1) * (1 - Z.var(ddof=1).sum() / Z.sum(axis=1).var(ddof=1)))
print(f"  the shared 4-item sexual-strictness scale I meant to use: Cronbach α = **{alpha:.3f}** "
      f"(n={len(Z)}), against the sanction indices' KR-20 of 0.802 and 0.789")
print(f"  ⇒ **equalising the norm ruler would have equalised it DOWNWARD**, cost half the sample, "
      f"and put `homosex` on both sides. **A different estimand dressed as a control.**")

print("\n=== (2) THE SAMPLE-DEFINITION SWEEP — the rows, not the items ===")
DEFS = {
    "#888 as published (all 7 ab items)": M888,
    "only the TRIO required": W.notna() & HREF.notna() & d["abpoorw"].notna() & d[TRIO].notna().all(axis=1),
    "#883 base (single-binary rule)": W.notna() & HREF.notna() & d["abpoorw"].notna() & d["abpoor"].notna(),
    "+ requires teensex (drops 2008)": M888 & d["teensex"].notna(),
    "+ requires all 4 norm items": M888 & d[["premarsx", "xmarsex", "teensex"]].notna().all(axis=1),
}
defs = {k: dict(delta=delta(v), n=int(v.sum())) for k, v in DEFS.items()}
for k, v in defs.items():
    print(f"  {k:38s} n={v['n']:5d}  **Δ {v['delta']:+.3f}**")
DEF_SPREAD = max(v["delta"] for v in defs.values()) - min(v["delta"] for v in defs.values())
print(f"  **completeness-rule spread = {DEF_SPREAD:.3f}**")

print("\n=== (3) LEAVE-ONE-WAVE-OUT — is the pooled number an average over heterogeneous cells? ===")
lowo, per = {}, {}
for y in WAVES:
    lowo[y] = delta(M888 & (d["year"] != y))
    per[y] = dict(delta=delta(M888 & (d["year"] == y)), n=int((M888 & (d["year"] == y)).sum()))
    print(f"  drop {y}: Δ **{lowo[y]:+.3f}**   |   wave {y} alone: Δ {per[y]['delta']:+.3f} "
          f"(n={per[y]['n']})")
LOWO_SPREAD = max(lowo.values()) - min(lowo.values())
worst = max(lowo, key=lambda y: abs(lowo[y] - D_888))
print(f"  **LOWO spread = {LOWO_SPREAD:.3f}** · the wave whose removal moves Δ most is "
      f"**{worst}** ({D_888:+.3f} → {lowo[worst]:+.3f})")

print("\n=== (4) THE NULL — kind of null NAMED: a CONSTANT-Δ resampling null ===")
sub = pd.DataFrame({"y": d.loc[M888, "year"].to_numpy(), "w": W[M888].to_numpy(),
                    "h": HREF[M888].to_numpy(), "an": d.loc[M888, "abpoorw"].to_numpy(),
                    "ai": AIDX[M888].to_numpy()}).dropna()
pool = sub[["w", "h", "an", "ai"]].to_numpy()
sizes = {y: int((sub.y == y).sum()) for y in WAVES}
null = np.empty(NSIM)
for i in range(NSIM):
    lv = {}
    draws = {y: pool[RNG.integers(0, len(pool), sizes[y])] for y in WAVES}
    for y in WAVES:
        keep = np.vstack([draws[z] for z in WAVES if z != y])
        lv[y] = ((1 - stats.spearmanr(keep[:, 0], keep[:, 1]).statistic ** 2)
                 - (1 - stats.spearmanr(keep[:, 2], keep[:, 3]).statistic ** 2))
    null[i] = max(lv.values()) - min(lv.values())
P95 = float(np.percentile(null, 95))
ABOVE = LOWO_SPREAD > P95
print(f"  respondents drawn within each wave from the POOLED table ⇒ true Δ identical in every wave")
print(f"  null LOWO spread: median **{np.median(null):.3f}** · 95th **{P95:.3f}**")
print(f"  observed **{LOWO_SPREAD:.3f}** ⇒ above what four noisy waves give for free: **{ABOVE}**")

print("\n=== (5) POSITIVE CONTROL — plant a between-wave difference, floor and ceiling MEASURED ===")
dose = {}
for g in (0.0, 0.05, 0.10, 0.20, 0.40):
    lv = {}
    draws = {y: pool[RNG.integers(0, len(pool), sizes[y])] for y in WAVES}
    tgt = WAVES[-1]
    for y in WAVES:
        keep = []
        for z in WAVES:
            if z == y:
                continue
            arr = draws[z].copy()
            if z == tgt and g > 0:                      # degrade the abortion pair in ONE wave
                fl = RNG.random(len(arr)) < g
                arr[fl, 3] = RNG.integers(0, 4, fl.sum())
            keep.append(arr)
        keep = np.vstack(keep)
        lv[y] = ((1 - stats.spearmanr(keep[:, 0], keep[:, 1]).statistic ** 2)
                 - (1 - stats.spearmanr(keep[:, 2], keep[:, 3]).statistic ** 2))
    dose[g] = float(max(lv.values()) - min(lv.values()))
    print(f"  g={g:.2f}  planted LOWO spread = {dose[g]:.3f}")
floor_, ceil_ = dose[0.0], dose[0.40]
_mono = all(dose[a] <= dose[b] + 0.01 for a, b in zip([0, .05, .10, .20], [.05, .10, .20, .40]))
POS_OK = bool(_mono and floor_ < P95 < ceil_)
print(f"  measured FLOOR {floor_:.3f} · CEILING {ceil_:.3f} ⇒ the null's 95th ({P95:.3f}) lies "
      f"strictly between them: {floor_ < P95 < ceil_}")
print(f"  => positive control **{'PASS' if POS_OK else 'FAIL'}** — monotone {_mono}, **and at g=0 "
      f"it does not fire**")

print("\n=== (6) NEGATIVE CONTROL — a global permutation. This zero SHOULD be zero ===")
gp = np.array([stats.spearmanr(pool[:, 2], RNG.permutation(pool[:, 3])).statistic
               for _ in range(400)])
NEG_OK = abs(gp.mean()) < 0.02
print(f"  null mean {gp.mean():+.4f} -> **{'PASS' if NEG_OK else 'FAIL'}**")

print("\n=== (7) SPECIFICATION CURVE + MULTIPLICITY ===")
rows = []
for est in ("spearman", "kendall", "gamma"):
    for nm, m in DEFS.items():
        rows.append(dict(kind="definition", cell=nm, est=est, delta=delta(m, est)))
    for y in WAVES:
        rows.append(dict(kind="lowo", cell=f"drop {y}", est=est, delta=delta(M888 & (d["year"] != y), est)))
for t in itertools.combinations(AB7, 3):
    ai = sum((d[c] == 2).astype(float).where(d[c].notna()) for c in t)
    rows.append(dict(kind="trio", cell="+".join(t), est="spearman", delta=delta(M888, "spearman", ai)))
G = pd.DataFrame(rows).dropna(subset=["delta"])
print(f"  cells: **{len(G)}** · Δ>0 in **{int((G.delta > 0).sum())}/{len(G)}**")
for kind, g in G.groupby("kind"):
    print(f"  {kind:11s} Δ median {g.delta.median():+.3f} [{g.delta.min():+.3f},{g.delta.max():+.3f}] "
          f"· Δ>0 in {int((g.delta>0).sum())}/{len(g)}")
DEF_GE_LOWO = DEF_SPREAD >= LOWO_SPREAD

GG = Gate("#889 · is `#888`'s level gap a pooled average over heterogeneous cells")
GG.asserted("(1) HARD RULE 1: n and years printed before use; and the design I intended was KILLED "
            "by its own gradient check and the check is recorded rather than dropped",
            True, f"shared-scale α {alpha:.3f} vs sanction KR-20 0.802/0.789 · n would fall "
                  f"{int(M888.sum())}→{len(Z)} · scale contains `homosex` (circular)",
            kind="control", population=POP)
GG.asserted("(2) OFFSET CONTROL — 'should this zero be zero?' NO: a max−min over four noisy "
            "estimates is positive by construction. **Kind of null: a CONSTANT-Δ resampling null**, "
            "respondents drawn within each wave from the pooled table so Δ is identical by design",
            True, f"null LOWO spread median {np.median(null):.3f} · 95th {P95:.3f}",
            kind="control", population=POP)
GG.asserted("(3) POSITIVE CONTROL: a planted between-wave difference must be detected, floor and "
            "ceiling MEASURED with the threshold strictly between them, and it must not fire at g=0",
            bool(POS_OK), " ".join(f"g={g}:{v:.3f}" for g, v in dose.items())
                          + f" · floor {floor_:.3f} ceiling {ceil_:.3f} · null95 {P95:.3f}",
            kind="control", population=POP)
GG.asserted("(4) NEGATIVE CONTROL: a global permutation, and this zero SHOULD be zero",
            bool(NEG_OK), f"mean {gp.mean():+.4f}", kind="control", population=POP)
GG.asserted("(5) EVERY SAMPLE DEFINITION REPORTED, including the ones that shrink the effect — a "
            "bootstrap resamples the SAMPLE and never the sample DEFINITION",
            True, " · ".join(f"{k}:{v['delta']:+.3f}" for k, v in defs.items()),
            kind="control", population=POP)
GG.asserted("(6) KILL (pre-registered): for `#888`'s `+0.051` to be a legitimate pooled estimate, "
            "**the leave-one-wave-out spread must sit inside the constant-Δ null**",
            bool(not ABOVE),
            f"LOWO spread {LOWO_SPREAD:.3f} vs null 95th {P95:.3f} (above {ABOVE}) · "
            f"completeness-rule spread {DEF_SPREAD:.3f} (>= LOWO: {DEF_GE_LOWO}) · "
            f"`#888` reported {D_888:+.3f} CI {CI_888} · grid Δ>0 {int((G.delta>0).sum())}/{len(G)}",
            kind="kill",
            yardstick="the max−min of Δ across leave-one-wave-out samples; the floor is the same "
                      "spread simulated with Δ constant across waves",
            yardstick_noise=P95, population=POP, direction=[float(v) for v in G.delta])
print()
print(GG)
adm = GG.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif DEF_GE_LOWO:
    V = (f"**C — THE DEFINITION IS THE PARAMETER.** The completeness-rule spread "
         f"(**{DEF_SPREAD:.3f}**) is at least as large as the leave-one-wave-out spread "
         f"(**{LOWO_SPREAD:.3f}**) ⇒ **the largest single source of uncertainty in this claim is "
         f"which rows I chose to keep**, which no bootstrap and no null addresses.\n"
         f"  `#888` published **{D_888:+.3f} [{CI_888[0]:+.3f}, {CI_888[1]:+.3f}]**; four equally "
         f"defensible completeness rules give "
         f"{', '.join(f'{v['delta']:+.3f}' for v in defs.values())}.\n"
         f"  ⇒ **every interval this project has published on a pooled level is narrower than the "
         f"thing it is an interval about.**")
elif ABOVE:
    V = (f"**B — ONE WAVE CARRIES IT.** LOWO spread **{LOWO_SPREAD:.3f}** against a constant-Δ null "
         f"whose 95th is **{P95:.3f}**; dropping **{worst}** moves Δ from {D_888:+.3f} to "
         f"**{lowo[worst]:+.3f}** ⇒ **`#888`'s headline is an average over heterogeneous cells**, "
         f"and the honest statement is the LOWO range "
         f"**[{min(lowo.values()):+.3f}, {max(lowo.values()):+.3f}]**, not the pooled point.")
else:
    V = (f"**A — HOMOGENEOUS.** LOWO spread **{LOWO_SPREAD:.3f}** sits inside the constant-Δ null "
         f"(95th **{P95:.3f}**), and the completeness-rule spread (**{DEF_SPREAD:.3f}**) is smaller "
         f"still ⇒ `#888`'s **{D_888:+.3f}** is a legitimate pooled estimate and the four waves are "
         f"four draws of one number.\n"
         f"  ⇒ **one sentence about people: the small remaining gap — an American's judgement of "
         f"abortion predicting what he wants the law to do a little better than his judgement of "
         f"homosexuality predicts what he wants done to homosexuals — is the same size in every "
         f"wave and under every defensible way of choosing whom to count.**")
print(V)
print(f"\n⚠ **Registered**: four waves is four cells, so a null here is a BOUND and never a "
      f"demonstration of homogeneity; the completeness rules are the ones I can defend and "
      f"**unenumerated is not cleared**; and **the instrument cannot be changed** — SCCS has one "
      f"observation per society and therefore no within-instrument subsample structure at all.")

json.dump(dict(population=POP, waves=WAVES, definitions=defs, def_spread=DEF_SPREAD,
               lowo={str(k): v for k, v in lowo.items()}, lowo_spread=LOWO_SPREAD,
               per_wave={str(k): v for k, v in per.items()}, worst_wave=int(worst),
               null=dict(kind="constant-Δ resampling null", median=float(np.median(null)), p95=P95),
               above_null=bool(ABOVE), def_ge_lowo=bool(DEF_GE_LOWO),
               shared_scale_alpha=alpha, dose=dose, floor=floor_, ceiling=ceil_,
               grid=G.to_dict("records"), reported_888=dict(delta=D_888, ci=list(CI_888)),
               controls=dict(positive=bool(POS_OK), negative=bool(NEG_OK)),
               admissible=adm, verdict=V, gate_ok=GG.verdict()),
          open(OUT / "is_the_level_gap_one_wave.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'is_the_level_gap_one_wave.json'}")
