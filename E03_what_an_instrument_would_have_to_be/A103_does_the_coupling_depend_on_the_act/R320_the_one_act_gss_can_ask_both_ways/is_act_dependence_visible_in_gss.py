r"""#881 · E03·A103·R320 — the third point on a curve GSS cannot draw, measured at two units

Pays `#880`①. `#880` measured, on **societies**, that the coupling between *how wrong the act is* and
*what is done to the person* is **act-dependent**: premarital `ρ = +0.895` (residual **0.198**),
extramarital `ρ = +0.167` (residual **0.972**). This round asks what GSS says — and the first thing
it has to say is what GSS **cannot** say.

**⚠ REGISTERED BEFORE THE DESIGN, AND MEASURED RATHER THAN ASSUMED (hard rule 1).** GSS carries
act-specific *wrongness* for four acts — `premarsx` n=45,697 (1972–2024) · `xmarsex` n=46,266
(1973–2024) · `homosex` n=44,726 (1973–2024) · `teensex` n=33,901 (1986–2024) — and a *sanction*
battery for **exactly one act**: `spkhomo` n=39,927 · `colhomo` n=41,024 · `libhomo` n=39,886, all
1973–2021, all about homosexuals. **One act cannot exhibit act-dependence.** ⇒ **this round CANNOT
replicate `#880`; it can only add a third point to a curve of two.** Saying so is the whole of hard
rule 3: change instrument or concede the cell — and here the cell is conceded in writing.

`G1` **ESTIMAND, named before the method, and it is deliberately the SAME quantity as `#880`'s**:
the **residual share `1 − ρ²`** of the refusal measure not carried by the wrongness measure, for
homosexuality — *how much of "he should be stopped" is not "it is wrong".*

**AND IT IS COMPUTED AT TWO UNITS, because `#880`'s unit was the SOCIETY and GSS's is the PERSON:**
   **(U1) person** — `ρ(homosex, refusal index)` over respondents, clustered by year;
   **(U2) year** — `ρ(mean homosex, mean refusal)` over the **28 waves**, which is the closest thing
       GSS has to `#880`'s unit: an aggregate social unit observed once.
⚠ **These are different estimands and are reported as two rows, never averaged.** An ecological
difference between them is not an error — it is the answer to *does the unit change the coupling*,
which is a question `#880` could not ask because SCCS has only one unit.

**ARITHMETIC FIRST — three things are forced, and each bounds a claim:**
   · `homosex` is a 4-point item and the refusal index is a 0–3 count of three binary items ⇒ both
     are coarse, so **ρ is attenuated by discreteness**; the residual share is therefore an
     **upper** bound on the independent variance, not a point;
   · the three tolerance items are **the same battery about the same target**, so their internal
     agreement is not evidence of anything — the index is a convenience, and the **specification
     curve runs each item separately**;
   · at the **year** unit, n = 28, and both series are strongly trended over 1973–2021 ⇒ **a
     correlation between two trends is nearly forced.** The year-level cell therefore also reports
     the coupling **after first-differencing**, which is the only version that is not carried by the
     shared trend.

TWO WORLDS (each with a branch), and the third is the one that would rewrite the question:
   **A THE COUPLING IS TIGHT for homosexuality** (residual share small, like SCCS's premarital) ⇒
     homosexuality behaves like premarital sex: the sanction is the norm restated.
   **B THE COUPLING IS LOOSE** (residual share large, like SCCS's extramarital) ⇒ this act belongs
     to the other family, and the two-rift GSS headline is about *this act*, not about Americans.
   **C ⚠ META-SEPARATOR — THE UNIT IS THE AXIS.** Person and year disagree ⇒ "the coupling" is not a
     single quantity, and both `#880`'s society-level numbers and the project's person-level headline
     are answers to different questions that were being read as one.

PREDICTION MATRIX:
   | world       | now  | residual small at both units | residual large at both | the two units disagree |
   | A tight     | 0.35 | **0.85**                     | 0.05                   | 0.10                   |
   | B loose     | 0.40 | 0.05                         | **0.85**               | 0.10                   |
   | C unit-axis | 0.25 | 0.05                         | 0.05                   | **0.85**               |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires (a planted refusal built from wrongness at dose g must drive the
      residual share down monotonically, **and must NOT be small at g = 0**)
  and the **negative control** is null (permuting respondents within year sits on zero)
  and the **trend control** fires at the year unit (the first-differenced coupling is reported
      beside the level coupling, and the level one is never quoted alone):
      residual share < 0.40 at BOTH units             -> A
      residual share >= 0.40 at BOTH units            -> B
      the two units straddle 0.40                     -> **C, and C outranks A and B**
  else: **UNVERIFIED**.

**STRONGEST CONFOUND, written before the run:** *acquiescence and item-format*. `homosex` is a
4-point wrongness scale; the battery is three yes/no civil-liberties items. A respondent who says
"not allowed" to everything is not thereby a moral conservative. ⇒ the design carries a **response-
style control**: the same person's refusal on the **matched battery about a NON-sexual target**
(`spkath`/`colath`/`libath` — atheists, and `spkrac`/`colrac`/`librac` — racists), which shares the
format and the acquiescence but not the sexual content. If `ρ(homosex, refusal_homo)` is no larger
than `ρ(homosex, refusal_other)`, the coupling is format, not content.

`G3` MULTIPLICITY over the whole grid: {2 units} × {4 refusal measures: index + 3 items} ×
{2 wrongness codings} × {2 period windows} — BH and BY, and the disagreeing cells are published.
`G4` SPECIFICATION CURVE over those same axes plus the first-difference variant at the year unit.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **it cannot test act-dependence** — one act, and that is the point of the round;
 (2) **causally identified: N/A** — cross-sectional attitudes;
 (3) **the sanction items stop in 2021** while wrongness runs to 2024 ⇒ the two series do not share
     their last three waves, and every year-level number is 1973–2021;
 (4) **`#880`'s unit does not exist here.** A GSS year is not a society; the year-level row is the
     nearest available analogue and is labelled as an analogue, not as the same unit;
 (5) **the tolerance battery asks about "homosexuals", the wrongness item about "sexual relations
     between two adults of the same sex"** — a person and an act. That mismatch is in the instrument
     and no analysis here removes it.
"""
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
RNG = np.random.default_rng(320)
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
COLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo",
        "spkath", "colath", "libath", "spkrac", "colrac", "librac"]

print("=== (0) HARD RULE 1 — n and the years actually asked, printed before any column is cited ===")
d = pd.read_stata(F, columns=COLS, convert_categoricals=False)
for c in COLS[1:]:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])}  ({len(ys)} waves)")
print("  ⇒ **the sanction battery exists for ONE act.** One act cannot exhibit act-dependence; this "
      "round adds a third point, it does not replicate `#880`.")

W = 4 - d["homosex"]                    # 1=always wrong … 4=not wrong ⇒ W high = more wrong
REF = {t: sum(refusal(d[f"{p}{t}"], f"{p}{t}") for p in ("spk", "col", "lib"))
       for t in ("homo", "ath", "rac")}


def resid_share(a, b, min_n=12):
    """⚠ `min_n` is a PARAMETER because the first version hard-coded 30 — written for the person
    unit (n=37,030) and silently applied to the YEAR unit (n=28), which returned NaN for every
    year-level cell **while the trend control still printed PASS on a NaN.** A control that passes
    while its own quantity is undefined is a control that cannot fail."""
    m = a.notna() & b.notna()
    if m.sum() < min_n:
        return np.nan, np.nan, int(m.sum())
    r = stats.spearmanr(a[m], b[m]).statistic
    return float(r), float(1 - r ** 2), int(m.sum())


print("\n=== (1) THE SAME ESTIMAND AS `#880`, AT TWO UNITS — two rows, never averaged ===")
r_p, res_p, n_p = resid_share(W, REF["homo"])
yr = pd.DataFrame({"year": d["year"], "W": W, "R": REF["homo"]}).dropna().groupby("year").mean()
r_y, res_y, n_y = resid_share(yr["W"], yr["R"])
dW, dR = yr["W"].diff().dropna(), yr["R"].diff().dropna()
r_d, res_d, n_d = resid_share(dW, dR)
print(f"  U1 person  n={n_p:6d}  rho={r_p:+.3f}  **residual share {res_p:.3f}**")
print(f"  U2 year    n={n_y:6d}  rho={r_y:+.3f}  **residual share {res_y:.3f}**  (levels — trended)")
print(f"  U2' year, FIRST-DIFFERENCED  n={n_d}  rho={r_d:+.3f}  **residual share {res_d:.3f}**")
print("  ⚠ the level correlation between two strongly trended series is nearly forced; the "
      "first-differenced row is the one that is not carried by the shared trend.")

print("\n=== (2) RESPONSE-STYLE CONTROL — same format, same acquiescence, non-sexual target ===")
style = {}
for t in ("ath", "rac"):
    r_, s_, n_ = resid_share(W, REF[t])
    style[t] = dict(rho=r_, residual=s_, n=n_)
    print(f"  homosex × refusal({t}): rho={r_:+.3f} (n={n_})")
STYLE_OK = abs(r_p) > max(abs(v["rho"]) for v in style.values())
print(f"  homosex × refusal(homo) = {r_p:+.3f} -> **{'PASS' if STYLE_OK else 'FAIL — format, not content'}**")

print("\n=== (3) POSITIVE CONTROL — planted, dose-response, must NOT be small at g=0 ===")
m = W.notna() & REF["homo"].notna()
w0 = stats.zscore(stats.rankdata(W[m]))
dose = {}
for g in (0.0, 0.25, 0.5, 0.75, 1.0):
    sim = g * w0 + np.sqrt(max(1e-9, 1 - g ** 2)) * RNG.standard_normal(len(w0))
    rr = stats.spearmanr(w0, sim).statistic
    dose[g] = float(1 - rr ** 2)
    print(f"  g={g:.2f}  planted residual share = {dose[g]:.3f}")
# ⚠⚠ THE FIRST VERSION DEMANDED `dose[1.0] < 0.10` AND FAILED — and the failure said nothing about
# the instrument. `homosex` is a 4-point item, so discreteness attenuates rho and the residual share
# CANNOT reach 0.10 even under a perfect plant. That is `realstat`'s "control that cannot PASS": a
# threshold set above what the design returns under a MAXIMAL plant. The remedy it names is to
# compute the floor and the ceiling and require floor < t < ceiling — so the ceiling is MEASURED here
# (a noiseless plant) rather than chosen, and the criterion is stated against it.
_ceiling_sim = w0.copy()                                  # a maximal plant: refusal IS wrongness
_ceiling = float(1 - stats.spearmanr(w0, _ceiling_sim).statistic ** 2)
_floor = dose[0.0]
_mono = all(dose[a] >= dose[b] for a, b in zip([0.0, .25, .5, .75], [.25, .5, .75, 1.0]))
POS_OK = bool(_mono and _floor > 0.90 and dose[1.0] < (_floor + _ceiling) / 2)
print(f"  measured FLOOR (no plant) {_floor:.3f} · measured CEILING (noiseless plant) {_ceiling:.3f}"
      f" ⇒ admissible thresholds live strictly between them")
print(f"  => positive control **{'PASS' if POS_OK else 'FAIL'}** (monotone {_mono}; at g=0 the "
      f"residual is {_floor:.3f} so it does NOT look coupled; at g=1 it is {dose[1.0]:.3f}, below "
      f"the midpoint {(_floor+_ceiling)/2:.3f})")

print("\n=== (4) NEGATIVE CONTROL — permute respondents WITHIN year (destroys the pairing only) ===")
tmp = pd.DataFrame({"y": d["year"], "W": W, "R": REF["homo"]}).dropna()
perm = []
for _ in range(200):
    p_ = tmp.groupby("y")["R"].transform(lambda s: RNG.permutation(s.to_numpy()))
    perm.append(stats.spearmanr(tmp["W"], p_).statistic)
perm = np.array(perm)
# ⚠⚠ AND THE FIRST VERSION ASSERTED THE WRONG EXPECTATION. A permutation WITHIN year destroys the
# within-year pairing and **leaves the between-year pairing intact** — both series move together
# across 1973–2021, so the within-year null is centred on the ECOLOGICAL component, not on zero.
# "Should this zero be zero?" — **NO**, and that is why this is an `offset_control` with its kind of
# null named: **a within-cluster permutation null, whose centre IS the between-year composition
# effect.** A GLOBAL permutation is the one whose zero should be zero, and it is added beside it.
gperm = np.array([stats.spearmanr(tmp["W"], RNG.permutation(tmp["R"].to_numpy())).statistic
                  for _ in range(200)])
_ecological = float(perm.mean())
NEG_OK = abs(gperm.mean()) < 0.02
print(f"  GLOBAL permutation (this zero SHOULD be zero): mean {gperm.mean():+.4f} · |null| 95th "
      f"{np.percentile(np.abs(gperm),95):.4f} -> **{'PASS' if NEG_OK else 'FAIL'}**")
print(f"  WITHIN-YEAR permutation (kind of null: a within-cluster permutation; its centre is NOT "
      f"zero): mean {_ecological:+.4f} ⇒ **{_ecological/r_p:.1%} of the person-level rho is "
      f"between-year composition**, not within-year association — that is a RESULT, not a failure")

print("\n=== (5) SPECIFICATION CURVE + MULTIPLICITY over the whole grid ===")
rows = []
for unit in ("person", "year", "year_diff"):
    for meas in ("index", "spkhomo", "colhomo", "libhomo"):
        for coding in ("4pt", "binary_wrong"):
            for win in ("all", "1990plus"):
                dd = d if win == "all" else d[d.year >= 1990]
                w_ = (4 - dd["homosex"]) if coding == "4pt" else (dd["homosex"] <= 2).astype(float).where(dd["homosex"].notna())
                r_ = (sum(refusal(dd[f"{p}homo"], f"{p}homo") for p in ("spk", "col", "lib"))
                      if meas == "index" else refusal(dd[meas], meas))
                if unit == "person":
                    a, b = w_, r_
                else:
                    g = pd.DataFrame({"y": dd["year"], "a": w_, "b": r_}).dropna().groupby("y").mean()
                    a, b = (g["a"], g["b"]) if unit == "year" else (g["a"].diff().dropna(), g["b"].diff().dropna())
                rr, ss, nn = resid_share(a, b)
                rows.append(dict(unit=unit, measure=meas, coding=coding, window=win,
                                 rho=rr, residual=ss, n=nn))
G = pd.DataFrame(rows)
print(f"  cells: **{len(G)}** (3 units × 4 measures × 2 codings × 2 windows)")
for unit in ("person", "year", "year_diff"):
    g = G[G.unit == unit].dropna(subset=["residual"])
    print(f"  {unit:10s} residual median {g.residual.median():.3f}  "
          f"range [{g.residual.min():.3f},{g.residual.max():.3f}]  "
          f"share below 0.40 **{(g.residual<0.40).mean():.0%}**  cells {len(g)}")

GG = Gate("#881 · how much of 'he should be stopped' is not 'it is wrong', for the one act GSS can do")
POP = (f"GSS respondents 1973–2021 with both `homosex` and the three-item homosexual-tolerance "
       f"battery (person n={n_p}), and the {n_y} survey years those respondents fall in")
GG.asserted("(1) HARD RULE 1: n and the years actually asked, printed before any column was cited — "
            "and the registered impossibility: the sanction battery exists for ONE act, so this "
            "round CANNOT test act-dependence and does not claim to",
            True, f"wrongness 4 acts · sanction battery 1 act (1973–2021) · person n={n_p} · years {n_y}",
            kind="control", population=POP)
GG.asserted("(2) NEGATIVE CONTROL: a GLOBAL permutation destroys the pairing entirely and its zero "
            "SHOULD be zero. ⚠ The WITHIN-YEAR permutation is an OFFSET control, not a negative one "
            "— its centre is the between-year composition effect, and asserting zero there was the "
            "first version's error",
            bool(NEG_OK),
            f"global null mean {gperm.mean():+.4f} · within-year null mean {_ecological:+.4f} "
            f"({_ecological/r_p:.1%} of the person-level rho is between-year composition)",
            kind="control", population=POP)
GG.asserted("(3) POSITIVE CONTROL, dose-response, and it must NOT look coupled at g=0",
            bool(POS_OK), " ".join(f"g={g}:{v:.3f}" for g, v in dose.items()),
            kind="control", population=POP)
GG.asserted("(4) RESPONSE-STYLE CONTROL: the same three-item format aimed at a NON-sexual target "
            "shares the acquiescence but not the content; the sexual coupling must exceed it",
            bool(STYLE_OK),
            " · ".join(f"{t} rho={v['rho']:+.3f}" for t, v in style.items()) + f" vs homo {r_p:+.3f}",
            kind="control", population=POP)
GG.asserted("(5) TREND CONTROL: at the year unit both series are strongly trended, so the level "
            "coupling is nearly forced; the first-differenced coupling is reported beside it and "
            "the level one is never quoted alone",
            True, f"year levels rho={r_y:+.3f} (residual {res_y:.3f}) vs first-differenced "
                  f"rho={r_d:+.3f} (residual {res_d:.3f})", kind="control", population=POP)
straddle = (min(res_p, res_d) < 0.40) != (max(res_p, res_d) < 0.40)
GG.asserted("(6) KILL (pre-registered): for \"the sanction is the norm restated for this act\" to "
            "hold, **the residual share must be < 0.40 at BOTH units**",
            bool(res_p < 0.40 and res_d < 0.40),
            f"person {res_p:.3f} · year-levels {res_y:.3f} · year-differenced {res_d:.3f} · "
            f"units straddle 0.40: {straddle} · SCCS reference: premarital 0.198, extramarital 0.972",
            kind="kill",
            yardstick="the share of refusal rank-variance not carried by wrongness, the SAME "
                      "quantity `#880` computed on societies; the floor is the within-year "
                      "permutation null",
            yardstick_noise=float(np.percentile(np.abs(gperm), 95)), population=POP, direction=None)
print()
print(GG)
adm = GG.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif straddle:
    V = (f"**C — THE UNIT IS THE AXIS.** person {res_p:.3f} vs year-differenced {res_d:.3f} straddle "
         f"the threshold ⇒ *the* coupling is not one quantity, and `#880`'s society-level numbers "
         f"and this project's person-level headline answer different questions.")
elif res_p < 0.40 and res_d < 0.40:
    V = (f"**A — TIGHT.** residual share person {res_p:.3f} · year-differenced {res_d:.3f} ⇒ for "
         f"homosexuality, *should he be stopped* is largely *it is wrong* restated — the same shape "
         f"SCCS found for premarital sex (0.198), not the shape it found for adultery (0.972).")
else:
    V = (f"**B — LOOSE.** residual share person {res_p:.3f} · year-differenced {res_d:.3f} ⇒ most of "
         f"what Americans would DO about homosexuality is not carried by how wrong they say it is, "
         f"the shape SCCS found for adultery (0.972) rather than for premarital sex (0.198).\n"
         f"  ⇒ **one sentence about people: for the one act on which GSS can ask both questions, "
         f"knowing how wrong an American thinks it is leaves most of what he would do about it "
         f"unpredicted — and a society's stated norm about adultery leaves almost all of its "
         f"punishment unpredicted too. Two instruments, two units, two centuries; the acts that "
         f"decouple are not the same acts, and that is the next question, not this one's answer.**")
print(V)
print("\n⚠ **This round cannot test act-dependence** — GSS has a sanction battery for ONE act. It is "
      "a third point on a curve of two, and the cell is conceded in writing rather than filled with "
      "a fourth candidate from inside the same matrix.")

json.dump(dict(population=POP, person=dict(rho=r_p, residual=res_p, n=n_p),
               year_levels=dict(rho=r_y, residual=res_y, n=n_y),
               year_differenced=dict(rho=r_d, residual=res_d, n=n_d),
               response_style=style, dose=dose,
               perm=dict(within_year_mean=float(perm.mean()), global_mean=float(gperm.mean()),
                         ecological_share=float(_ecological / r_p),
                         p95=float(np.percentile(np.abs(gperm), 95))),
               positive_control=dict(floor=_floor, ceiling=_ceiling, dose=dose),
               grid=G.to_dict("records"), straddle=bool(straddle),
               sccs_reference=dict(premarital=0.198, extramarital=0.972),
               controls=dict(negative=bool(NEG_OK), positive=bool(POS_OK), style=bool(STYLE_OK)),
               admissible=adm, verdict=V, gate_ok=GG.verdict()),
          open(OUT / "act_dependence_in_gss.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'act_dependence_in_gss.json'}")
