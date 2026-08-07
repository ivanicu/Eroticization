r"""#883 · E03·A103·R322 — two acts, one instrument, the same 2,694 people: the first version that is not an analogy

Pays `#882`①. Every cross-instrument sentence this project has written about act-dependence has been
an **analogy between different acts**: SCCS decouples adultery and couples premarital sex, GSS could
only speak about homosexuality. `#882` found the pair that removes the analogy — GSS asks, about
**abortion**, both *is it **wrong*** and *should the **law allow** it*, of the same respondents.

**SO THIS ROUND IS THE PAIRED DESIGN.** Two acts, one questionnaire, one field period, and — this is
the part that matters — **the same 2,694 people answered both pairs.** Person composition, cohort,
period, response style and acquiescence are removed **by construction**, not by adjustment.

`G1` **ESTIMAND, named before the method**: `Δ = residual_share(homosexuality) −
residual_share(abortion)`, where `residual_share = 1 − ρ²` is the share of the sanction not carried
by the norm — **the same quantity `#880` computed on societies and `#881` on Americans**, so the
three numbers are commensurable by construction rather than by argument.

**ARITHMETIC FIRST — three things are forced, and each one bounds or strengthens the claim:**
   · **the SIGN is a coding artifact and is not a finding.** `abpoorw` runs 1 = always wrong → 4 =
     not wrong; `abpoor` runs 1 = yes, the law should allow → 2 = no. So the abortion coupling comes
     out NEGATIVE for the same reason the homosexuality one comes out positive. **Only |ρ| enters
     the estimand**, and the sign is printed so a reader can see it was handled;
   · **the format bias runs AGAINST the finding.** Homosexuality's sanction is a 0–3 count of three
     items; abortion's is a single **binary**. A binary is the *more* attenuated instrument, so if
     the abortion coupling comes out *tighter*, format cannot be why. Reported as a direction, and
     the format control below measures it anyway rather than resting on the argument;
   · `residual_share` is bounded in [0,1], so `Δ` is bounded in [−1,1] and **its null is not zero**:
     two attenuated couplings measured on one sample share their sampling error. *"Should this zero
     be zero?"* — **NO** ⇒ the null for `Δ` is an **`offset_control`**, and the **kind of null is
     named: a paired bootstrap over the 2,694 respondents**, which preserves the sharing.

FOUR WORLDS (each with a branch):
   **A ACT-DEPENDENCE, inside one instrument.** `Δ` clears its paired-bootstrap interval and survives
     the format control ⇒ what a person would DO about an act is tied to how wrong he thinks it is
     **to a degree that depends on the act**, in the same head, at the same moment.
   **B FORMAT.** The gap collapses when the homosexuality sanction is reduced to a single binary item
     ⇒ the finding was a 3-item count against a 1-item flag.
   **C QUESTION-MATCHING.** ⚠ **The strongest confound, written before the run.** Abortion's two
     items share a qualifier (*"if the family has a low income"*) while homosexuality's ask about
     **an act** and about **a person** (*"should a homosexual be allowed to speak"*). The tight
     abortion coupling could be *semantic matching*, not psychology. **The cross-reason sham tests
     it**: `abpoorw` (wrong if poor) × `abdefect` (law should allow for birth defects) — same act,
     **different reason**. If cross ≈ within, the coupling is general abortion attitude and the
     qualifier is doing nothing.
   **D ⚠ META-SEPARATOR — THE AXIS IS THE SANCTION, NOT THE ACT.** *Make it illegal* and *let him
     keep his library book* are not the same severity of thing to want. If the two couplings differ
     because one sanction is **the law** and the other is **a civil liberty**, then "act-dependence"
     is the wrong name for the phenomenon and the ordering of this project's whole question changes:
     it would be **sanction-severity dependence**, and SCCS's premarital/extramarital split would
     have to be re-read the same way.

PREDICTION MATRIX:
   | world        | now  | Δ clears its null | Δ dies under 1-item format | cross-reason ≈ within |
   | A act        | 0.35 | **0.85**          | 0.05                       | 0.10                  |
   | B format     | 0.15 | 0.10              | **0.85**                   | 0.10                  |
   | C matching   | 0.25 | 0.30              | 0.10                       | **0.85**              |
   | D sanction   | 0.25 | **0.85**          | 0.05                       | 0.10                  |
⚠ **A and D predict the SAME observable here** — that is stated rather than hidden, and it is why
world D cannot be selected by this round. Separating them needs a second sanction *format* for the
same act, which GSS does not carry. **Registered, not planned.**

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires (a planted sanction built from the norm at dose g drives the
      residual share down monotonically, with **floor and ceiling MEASURED** and the criterion
      strictly between them, and **it does not look coupled at g = 0**)
  and the **negative control** is null (a global permutation of respondents sits on zero)
  and the **format control** is reported (the homosexuality sanction reduced to each single item):
      Δ's paired-bootstrap 95% interval excludes 0, and Δ stays >0 in every single-item cell -> A
      Δ's interval includes 0, or Δ flips under any single-item cell                         -> B
      cross-reason sham >= within-reason coupling                                            -> **C, and C outranks A and B**
  else: **UNVERIFIED**.

`G3` MULTIPLICITY over the whole grid: {2 acts} × {4 sanction measures for homosexuality: index +
3 single items} × {2 abortion reasons} × {4 waves} — BH and BY, and the disagreeing cells published.
`G4` SPECIFICATION CURVE over those axes plus wave-by-wave, because a coupling measured on four
waves spanning 27 years is four measurements, not one.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **A and D are not separable here.** It would require a second sanction format for the same act —
     a legality item about homosexuality, or a civil-liberties item about abortion. GSS carries
     neither;
 (2) **causally identified: N/A** — cross-sectional attitudes, no intervention;
 (3) **abortion is not a member of the family SCCS coded** — it is a reproductive decision, not a
     sexual act. This round compares two acts *within GSS*; it does not thereby compare them to
     SCCS's premarital/extramarital pair;
 (4) **the abortion norm items exist on 4 waves only** (1991 · 1998 · 2008 · 2018) against a
     homosexuality battery running 1973–2021 ⇒ the paired sample is the abortion item's sample, and
     every number here is 1991–2018;
 (5) **no second coder / second release** — one questionnaire, one field house. Independent
     replication would need a different survey programme asking both questions about one act.
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
from lib.gss_polarity import refusal        # `#868`'s home, imported rather than copied

OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(322)
NBOOT = 4000
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
COLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo",
        "abpoorw", "abpoor", "abdefctw", "abdefect"]

print("=== (0) HARD RULE 1 — n and the years actually asked, before any column is cited ===")
d = pd.read_stata(F, columns=COLS, convert_categoricals=False)
for c in COLS[1:]:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])}  ({len(ys)} waves)")

W_HOMO = 4 - d["homosex"]                                  # high = more wrong
REF = {k: refusal(d[f"{k}homo"], f"{k}homo") for k in ("spk", "col", "lib")}
REF["index"] = sum(REF.values())
ACTS = {
    "homosexuality": dict(norm=W_HOMO, sancs={k: v for k, v in REF.items()}),
    "abortion_poor": dict(norm=d["abpoorw"], sancs={"abpoor": d["abpoor"]}),
    "abortion_defect": dict(norm=d["abdefctw"], sancs={"abdefect": d["abdefect"]}),
}
PAIRED = (W_HOMO.notna() & REF["index"].notna()
          & d["abpoorw"].notna() & d["abpoor"].notna())
print(f"\n  **the paired sample — people who answered BOTH pairs: n = {int(PAIRED.sum())}** "
      f"· waves {sorted(int(y) for y in d.loc[PAIRED, 'year'].unique())}")
print("  ⇒ person composition, cohort, period, response style and acquiescence are removed BY "
      "CONSTRUCTION, not by adjustment.")
if PAIRED.sum() < 200:
    raise SystemExit("STOP: an empty/underpowered paired sample must never be counted as a pass")
POP = (f"the {int(PAIRED.sum())} GSS respondents who answered both the homosexuality norm+battery "
       f"and the abortion(low-income) norm+legality pair, waves 1991/1998/2008/2018")


def rs(a, b, mask):
    m = a.notna() & b.notna() & mask
    r = stats.spearmanr(a[m], b[m]).statistic
    return float(r), float(1 - r ** 2), int(m.sum())


print("\n=== (1) THE TWO COUPLINGS, ON THE SAME PEOPLE ===")
r_h, s_h, n_h = rs(W_HOMO, REF["index"], PAIRED)
r_a, s_a, n_a = rs(d["abpoorw"], d["abpoor"], PAIRED)
print(f"  homosexuality  rho={r_h:+.3f} (|rho|={abs(r_h):.3f})  **residual share {s_h:.3f}**  n={n_h}")
print(f"  abortion(poor) rho={r_a:+.3f} (|rho|={abs(r_a):.3f})  **residual share {s_a:.3f}**  n={n_a}")
DELTA = s_h - s_a
print(f"  **Δ = {DELTA:+.3f}** — the sign of each rho is a CODING artifact (abpoor: 1=law should "
      f"allow, 2=no), so only |rho| enters the estimand.")

print("\n=== (2) OFFSET CONTROL — the null for Δ is NOT zero. Kind of null: a PAIRED BOOTSTRAP ===")
idx = np.flatnonzero(PAIRED.to_numpy())
wh, rh_, ap, aq = (W_HOMO.to_numpy()[idx], REF["index"].to_numpy()[idx],
                   d["abpoorw"].to_numpy()[idx], d["abpoor"].to_numpy()[idx])
boot = np.empty(NBOOT)
for i in range(NBOOT):
    j = RNG.integers(0, len(idx), len(idx))
    a_ = 1 - stats.spearmanr(wh[j], rh_[j]).statistic ** 2
    b_ = 1 - stats.spearmanr(ap[j], aq[j]).statistic ** 2
    boot[i] = a_ - b_
lo, hi = np.percentile(boot, [2.5, 97.5])
DELTA_OK = (lo > 0) or (hi < 0)
print(f"  Δ = {DELTA:+.3f}  paired-bootstrap 95% [{lo:+.3f}, {hi:+.3f}]  "
      f"**{'EXCLUDES 0' if DELTA_OK else 'includes 0'}**")
print("  ⚠ the bootstrap resamples PEOPLE, so both couplings move together — which is exactly the "
      "sharing that makes a zero-centred null wrong here.")

print("\n=== (3) FORMAT CONTROL — the homosexuality sanction reduced to each SINGLE binary item ===")
fmt = {}
for k in ("spk", "col", "lib", "index"):
    rr, ss, nn = rs(W_HOMO, REF[k], PAIRED)
    fmt[k] = dict(rho=rr, residual=ss, n=nn, delta=float(ss - s_a))
    print(f"  {k:6s} |rho|={abs(rr):.3f} residual={ss:.3f}  Δ vs abortion = {ss - s_a:+.3f}")
FMT_OK = all(v["delta"] > 0 for v in fmt.values())
print(f"  => Δ stays positive in every single-item cell: **{FMT_OK}** "
      f"(a binary is the MORE attenuated instrument, so format runs AGAINST the finding)")

print("\n=== (4) THE STRONGEST CONFOUND — question-matching. Cross-REASON sham, same act ===")
mask_ab = d["abpoorw"].notna() & d["abdefect"].notna() & d["abdefctw"].notna() & d["abpoor"].notna()
r_within_p, s_within_p, n_wp = rs(d["abpoorw"], d["abpoor"], mask_ab)
r_cross_p, s_cross_p, n_cp = rs(d["abpoorw"], d["abdefect"], mask_ab)
r_within_d, s_within_d, n_wd = rs(d["abdefctw"], d["abdefect"], mask_ab)
r_cross_d, s_cross_d, n_cd = rs(d["abdefctw"], d["abpoor"], mask_ab)
print(f"  within-reason  wrong(poor)  × law(poor)   |rho|={abs(r_within_p):.3f}  n={n_wp}")
print(f"  CROSS-reason   wrong(poor)  × law(defect) |rho|={abs(r_cross_p):.3f}  n={n_cp}")
print(f"  within-reason  wrong(defect)× law(defect) |rho|={abs(r_within_d):.3f}  n={n_wd}")
print(f"  CROSS-reason   wrong(defect)× law(poor)   |rho|={abs(r_cross_d):.3f}  n={n_cd}")
SHAM_OK = (abs(r_cross_p) < abs(r_within_p)) and (abs(r_cross_d) < abs(r_within_d))
print(f"  => cross < within in BOTH directions: **{'PASS' if SHAM_OK else 'FAIL — world C'}**")
print("  ⚠ if the qualifier did nothing, cross would equal within and the tight coupling would be "
      "general abortion attitude rather than reason-specific.")

print("\n=== (5) POSITIVE CONTROL — dose-response, floor and ceiling MEASURED, fails at g=0 ===")
z = stats.zscore(stats.rankdata(W_HOMO.to_numpy()[idx]))
dose = {}
for g in (0.0, 0.25, 0.5, 0.75, 1.0):
    sim = g * z + np.sqrt(max(1e-9, 1 - g ** 2)) * RNG.standard_normal(len(z))
    dose[g] = float(1 - stats.spearmanr(z, sim).statistic ** 2)
    print(f"  g={g:.2f}  planted residual share = {dose[g]:.3f}")
floor_, ceil_ = dose[0.0], float(1 - stats.spearmanr(z, z).statistic ** 2)
mono = all(dose[a] >= dose[b] for a, b in zip([0, .25, .5, .75], [.25, .5, .75, 1.0]))
POS_OK = bool(mono and floor_ > 0.90 and dose[1.0] < (floor_ + ceil_) / 2)
print(f"  measured FLOOR {floor_:.3f} · measured CEILING {ceil_:.3f} ⇒ thresholds live between them")
print(f"  => **{'PASS' if POS_OK else 'FAIL'}** (monotone {mono}; at g=0 it does NOT look coupled)")

print("\n=== (6) NEGATIVE CONTROL — a global permutation. This zero SHOULD be zero ===")
gp = np.array([stats.spearmanr(wh, RNG.permutation(rh_)).statistic for _ in range(300)])
NEG_OK = abs(gp.mean()) < 0.02
print(f"  null mean {gp.mean():+.4f} · |null| 95th {np.percentile(np.abs(gp),95):.4f} -> "
      f"**{'PASS' if NEG_OK else 'FAIL'}**")

print("\n=== (7) SPECIFICATION CURVE + MULTIPLICITY — wave by wave, because 4 waves is 4 measurements ===")
rows = []
for wave in sorted(int(y) for y in d.loc[PAIRED, "year"].unique()):
    mw = PAIRED & (d["year"] == wave)
    for k in ("index", "spk", "col", "lib"):
        rr, ss, nn = rs(W_HOMO, REF[k], mw)
        ra, sa, na = rs(d["abpoorw"], d["abpoor"], mw)
        rows.append(dict(wave=wave, sanction=k, res_homo=ss, res_abort=sa,
                         delta=float(ss - sa), n=min(nn, na)))
G = pd.DataFrame(rows)
print(f"  cells: **{len(G)}** (4 waves × 4 sanction measures)")
for wave, g in G.groupby("wave"):
    print(f"  {wave}  n≈{int(g.n.median()):4d}  res_homo {g.res_homo.median():.3f}  "
          f"res_abort {g.res_abort.median():.3f}  **Δ median {g.delta.median():+.3f}**  "
          f"share Δ>0 {(g.delta>0).mean():.0%}")
print(f"  whole grid: **{(G.delta>0).sum()}/{len(G)} cells with Δ>0**, "
      f"median {G.delta.median():+.3f}, range [{G.delta.min():+.3f},{G.delta.max():+.3f}]")

GG = Gate("#883 · two acts, one instrument, the same people — is the coupling act-dependent")
GG.asserted("(1) HARD RULE 1: n and the years actually asked printed before any column was cited; "
            "the paired sample is the abortion item's sample and every number is 1991–2018",
            True, f"paired n={int(PAIRED.sum())} · waves 1991/1998/2008/2018",
            kind="control", population=POP)
GG.asserted("(2) NEGATIVE CONTROL: a global permutation of respondents destroys the pairing "
            "entirely and its zero SHOULD be zero",
            bool(NEG_OK), f"null mean {gp.mean():+.4f}", kind="control", population=POP)
GG.asserted("(3) POSITIVE CONTROL, dose-response, floor and ceiling MEASURED rather than chosen, "
            "and it must not look coupled at g=0",
            bool(POS_OK), " ".join(f"g={g}:{v:.3f}" for g, v in dose.items())
                          + f" · floor {floor_:.3f} ceiling {ceil_:.3f}",
            kind="control", population=POP)
GG.asserted("(4) OFFSET CONTROL for Δ — 'should this zero be zero?' NO: two attenuated couplings on "
            "one sample share their sampling error. **Kind of null: a paired bootstrap over the "
            "respondents**, which preserves that sharing",
            bool(DELTA_OK), f"Δ={DELTA:+.3f} 95% [{lo:+.3f},{hi:+.3f}]",
            kind="control", population=POP)
GG.asserted("(5) FORMAT CONTROL: the homosexuality sanction reduced to each single binary item — a "
            "binary is the MORE attenuated instrument, so format runs against the finding",
            bool(FMT_OK), " · ".join(f"{k}:Δ{v['delta']:+.3f}" for k, v in fmt.items()),
            kind="control", population=POP)
GG.asserted("(6) SHAM (selects world C): cross-REASON, same act — wrong(poor)×law(defect) and "
            "wrong(defect)×law(poor) must both be weaker than their within-reason partners, or the "
            "coupling is general abortion attitude and the qualifier does nothing",
            bool(SHAM_OK),
            f"within poor {abs(r_within_p):.3f} vs cross {abs(r_cross_p):.3f} · "
            f"within defect {abs(r_within_d):.3f} vs cross {abs(r_cross_d):.3f}",
            kind="control", population=POP)
GG.asserted("(7) KILL (pre-registered): for the coupling to be act-dependent inside one instrument, "
            "**Δ's paired-bootstrap interval must exclude 0 AND Δ must stay positive in every "
            "single-item format cell**",
            bool(DELTA_OK and FMT_OK and DELTA > 0),
            f"Δ={DELTA:+.3f} [{lo:+.3f},{hi:+.3f}] · format cells all positive {FMT_OK} · "
            f"grid {(G.delta>0).sum()}/{len(G)} · SCCS reference 0.198/0.972 · GSS-homo 0.771",
            kind="kill",
            yardstick="the difference of two residual shares (1−ρ²) measured on the SAME people; "
                      "the floor is the paired bootstrap of that difference",
            yardstick_noise=float(hi - lo), population=POP,
            direction=[v["delta"] for v in fmt.values()])
print()
print(GG)
adm = GG.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif not SHAM_OK:
    V = ("**C — QUESTION-MATCHING.** The cross-reason sham is not weaker than the within-reason "
         "coupling, so the tight abortion pairing is general abortion attitude and the shared "
         "qualifier is doing nothing.")
elif DELTA_OK and FMT_OK and DELTA > 0:
    V = (f"**A/D — THE COUPLING IS ACT-DEPENDENT INSIDE ONE INSTRUMENT, and A and D cannot be "
         f"separated here.** On the **same {int(PAIRED.sum())} people**: homosexuality residual "
         f"**{s_h:.3f}**, abortion residual **{s_a:.3f}**, **Δ = {DELTA:+.3f}** with a paired "
         f"bootstrap 95% of [{lo:+.3f}, {hi:+.3f}], positive in **{(G.delta>0).sum()}/{len(G)}** "
         f"grid cells and in every single-item format cell.\n"
         f"  ⇒ **one sentence about people: ask the same person how wrong two different things are "
         f"and what should be done about them, and for one of those things his answer to the first "
         f"question nearly tells you the second, while for the other it barely tells you anything. "
         f"The gap is not in the person and not in the questionnaire — it is in the act.**\n"
         f"  ⚠ **and the name may still be wrong**: *make it illegal* and *let him keep his library "
         f"book* are not the same severity of demand, so this could be **sanction-severity "
         f"dependence** rather than act-dependence. GSS carries no second sanction format for either "
         f"act, so **world D is registered as unseparable here, not dismissed.**")
else:
    V = (f"**B — FORMAT, or unresolved.** Δ={DELTA:+.3f} with 95% [{lo:+.3f},{hi:+.3f}]; format "
         f"cells all-positive = {FMT_OK} ⇒ the gap does not survive its own controls.")
print(V)
print("\n⚠ **Registered, and it is the load-bearing limitation**: worlds A and D predict the same "
      "observable here. Separating *the act* from *what the sanction costs* needs a second sanction "
      "format for the same act — a legality item about homosexuality, or a civil-liberties item "
      "about abortion. **GSS carries neither.**")

json.dump(dict(population=POP, paired_n=int(PAIRED.sum()),
               homosexuality=dict(rho=r_h, residual=s_h, n=n_h),
               abortion_poor=dict(rho=r_a, residual=s_a, n=n_a),
               delta=DELTA, delta_ci=[float(lo), float(hi)], delta_excludes_zero=bool(DELTA_OK),
               format_control=fmt, sham=dict(within_poor=r_within_p, cross_poor=r_cross_p,
                                             within_defect=r_within_d, cross_defect=r_cross_d,
                                             ok=bool(SHAM_OK)),
               dose=dose, floor=floor_, ceiling=ceil_,
               perm=dict(mean=float(gp.mean()), p95=float(np.percentile(np.abs(gp), 95))),
               grid=G.to_dict("records"),
               reference=dict(sccs_premarital=0.198, sccs_extramarital=0.972, gss_homo_all=0.771),
               controls=dict(negative=bool(NEG_OK), positive=bool(POS_OK), offset=bool(DELTA_OK),
                             format=bool(FMT_OK), sham=bool(SHAM_OK)),
               admissible=adm, verdict=V, gate_ok=GG.verdict()),
          open(OUT / "two_acts_one_instrument.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'two_acts_one_instrument.json'}")
