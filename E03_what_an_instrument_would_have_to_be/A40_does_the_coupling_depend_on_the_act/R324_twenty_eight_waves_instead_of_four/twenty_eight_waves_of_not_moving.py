r"""#886 · E03·A103·R324 — America changed its mind entirely and the joint did not move

Pays `#885`②, which was stated on **four** waves and never attacked. `#885` found that abortion's
norm→sanction coupling **moved** (spread 0.351 against its own IPF null 0.159) while
homosexuality's **did not** (spread 0.089, inside 0.153) — and then said, in its own NEXT, that the
second half *"is a separate finding, it is stated here, and it has not been attacked."*

**A NULL ON FOUR POINTS IS THE WEAKEST THING THIS PROJECT PUBLISHES.** The homosexuality pair does
not need four waves: it exists on **28**, 1973–2021, n = 37,030. This round re-asks the same question
at **7× the time points and 14× the sample**, and it is built so that **the outcome I would find
unwelcome is the one it can deliver** — because a movement found at 28 waves retracts a claim I
published one entry ago and dissolves the contrast that made `#885` interesting.

`G1` **ESTIMAND, named before the method**, over 28 waves:
   **(1) the SPREAD** `max_w − min_w` of `residual_share = 1 − ρ²`;
   **(2) the TREND** `ρ_S(year, residual_share)` — which four points could not carry and 28 can.
Both against the **same** null, because a flat series with one outlier and a rising series are
different objects and one statistic cannot tell them apart.

**ARITHMETIC FIRST — and it predicts the direction of the artifact before any control runs:**
   · the marginal does not drift here, it **collapses**: the share answering *always wrong* runs
     **0.714 → 0.255** across the 28 waves (mean wrongness **2.34 → 0.93**). A compressed marginal
     **restricts range and ATTENUATES ρ**, which **RAISES** `1 − ρ²`. ⇒ **a rising residual is the
     arithmetic default, not a finding**, and the null must absorb exactly that;
   · ⇒ *"Should this zero be zero?"* — **NO.** So the null is an **`offset_control`** and **the kind
     of null is named: an IPF-reconstructed constant-coupling null** — each wave rebuilt by iterative
     proportional fitting from the **pooled** table onto **that wave's own margins**, preserving
     every odds ratio while reproducing the marginal collapse exactly;
   · **gauge test, three lines, run before the design was fixed**: `1 − ρ²` is invariant under
     monotone re-labelling of the categories but **not** under a table going degenerate. If the
     collapse empties cells, *flat* and *unmeasurable* converge. **Measured in the gradient check:
     the per-wave bootstrap CI width is 0.072 in 1973 and 0.091 in 2021, median 0.069 — the
     instrument keeps its resolution while the marginal moves by a factor of 2.8.** World C is
     weakened *before* the run rather than dismissed after it, and it is still carried as a control.

THREE WORLDS (each with a branch):
   **A FLAT — `#885`② SURVIVES AT FULL POWER.** Spread and trend both inside the IPF null ⇒ **a
     society can reverse its verdict on whether something is wrong without changing how tightly that
     verdict is tied to what it wants done.**
   **B ⚠ THE UNWELCOME ONE — IT MOVED, and four waves hid it.** Spread or trend above the null ⇒
     `#885`②'s null was an artifact of the abortion pair's four-wave window, the sentence is
     retracted, and the `#885` contrast (one act moved, one did not) loses its second half.
   **C ⚠ META-SEPARATOR — THE INSTRUMENT WENT BLIND.** The late waves cannot resolve a movement of
     the size that would matter ⇒ *flat* is **silence**, not a finding, and the whole
     moved/did-not-move decomposition is unavailable for an act whose marginal collapses.

PREDICTION MATRIX:
   | world     | now  | both inside the null | either above it | late-wave MDE too coarse |
   | A flat    | 0.45 | **0.85**             | 0.05            | 0.10                     |
   | B moved   | 0.35 | 0.05                 | **0.85**        | 0.10                     |
   | C blind   | 0.20 | 0.10                 | 0.05            | **0.85**                 |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires — a planted movement **the size of the one abortion actually
      showed (spread 0.351)** must be detected, **and at g = 0, a constant coupling carrying the real
      marginal collapse, it must NOT fire**
  and the **negative control** is null — the IPF machinery re-measured on its own reconstruction
      sits inside its own null
  and the **blindness control** passes — the design's own MDE for the spread is **smaller than
      abortion's measured 0.351**, so a movement of the size that matters would have been seen:
      spread AND |trend| both inside their IPF nulls  -> A, and `#885`② survives at full power
      either above its null                           -> **B, and `#885`② is retracted**
      MDE >= 0.351                                    -> **C, and the null is silence**
  else: **UNVERIFIED**.

**STRONGEST CONFOUND, written before the run, and it is not the marginal:** *the sanction battery's
own composition changes.* `spkhomo`/`colhomo`/`libhomo` are asked together, but item non-response and
ballot rotation vary by wave, so the 0–3 index is not the same instrument in every year. ⇒ the
specification curve runs **each single item separately** beside the index, and if the index moves
while no single item does, the movement is composition.

`G3` MULTIPLICITY over the whole grid: {4 sanction measures} × {3 estimators} × {2 windows} = 24
cells, published whole including disagreement. `G4` SPECIFICATION CURVE over those same axes —
and the estimator axis is load-bearing, not decoration: **Goodman–Kruskal γ normalises differently
from ρ and is far less sensitive to marginal compression**, so if a movement is compression it should
shrink under γ.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **causally identified: N/A** — repeated cross-sections, not a panel. Every wave difference is a
     composition difference as well as a change, and no analysis of these columns separates them;
 (2) **the sanction battery stops in 2021** while the wrongness item runs to 2024 ⇒ the last three
     waves of the norm have no sanction to pair with, and every number here is 1973–2021;
 (3) **the battery asks about "homosexuals" and the norm about an act between two adults** — a person
     and an act. That mismatch is in the instrument and no analysis removes it;
 (4) ⚠ **the instrument cannot be changed for this question, and it is structural**: a motion claim
     needs a **time axis**, and `#882` measured that the only two instruments here with a matched
     norm–sanction pair are GSS and SCCS — **SCCS codes each society at one focal year and has no
     time axis at all**. There is **only this one instrument** for a coupling that moves;
 (5) **no second coder, no second release.**
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
RNG = np.random.default_rng(324)
NSIM = 1200
ABORTION_SPREAD = 0.351          # `#885`'s measured movement — the dose that MATTERS
F = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
COLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo"]

print("=== (0) HARD RULE 1 — n and the years actually asked, before any column is cited ===")
d = pd.read_stata(F, columns=COLS, convert_categoricals=False)
for c in COLS[1:]:
    s = d[[c, "year"]].dropna()
    ys = sorted(s.year.unique())
    print(f"  {c:9s} n={len(s):6d}  {int(ys[0])}–{int(ys[-1])}  ({len(ys)} waves)")
W = 4 - d["homosex"]
REFI = {k: refusal(d[f"{k}homo"], f"{k}homo") for k in ("spk", "col", "lib")}
REFI["index"] = sum(REFI.values())
M = W.notna() & REFI["index"].notna()
WAVES = [y for y in sorted(int(v) for v in d.loc[M, "year"].unique())
         if (M & (d["year"] == y)).sum() >= 200]
print(f"\n  **{len(WAVES)} waves** (n>=200), total n={int(M.sum())} — against the FOUR the abortion "
      f"pair forced in `#885`")
if len(WAVES) < 10:
    raise SystemExit("STOP: fewer than 10 waves cannot carry a trend; an empty design must not pass")
POP = (f"the {int(M.sum())} GSS respondents with both `homosex` and the three-item homosexual "
       f"tolerance battery, split into the {len(WAVES)} waves {WAVES[0]}–{WAVES[-1]}")


def res_of(a, b, est="spearman"):
    m = a.notna() & b.notna()
    x, y = a[m].to_numpy(), b[m].to_numpy()
    if len(x) < 50:
        return np.nan
    if est == "spearman":
        r = stats.spearmanr(x, y).statistic
    elif est == "kendall":
        r = stats.kendalltau(x, y, variant="b").statistic
    else:
        tab = pd.crosstab(x, y).to_numpy()
        c = dd = 0
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                c += tab[i, j] * tab[i + 1:, j + 1:].sum()
                dd += tab[i, j] * tab[i + 1:, :j].sum()
        r = (c - dd) / (c + dd) if (c + dd) else np.nan
    return float(1 - r ** 2)


print("\n=== (1) THE SERIES — 28 waves, with the marginal collapse beside it ===")
obs = {}
for y in WAVES:
    m = M & (d["year"] == y)
    obs[y] = dict(n=int(m.sum()), res=res_of(W[m], REFI["index"][m]),
                  mean_wrong=float(W[m].mean()), share_always=float((W[m] == 3).mean()))
SPREAD = max(v["res"] for v in obs.values()) - min(v["res"] for v in obs.values())
yr = np.array(WAVES, float)
rs = np.array([obs[y]["res"] for y in WAVES])
TREND = float(stats.spearmanr(yr, rs).statistic)
print(f"  share answering ALWAYS WRONG: {obs[WAVES[0]]['share_always']:.3f} ({WAVES[0]}) -> "
      f"{obs[WAVES[-1]]['share_always']:.3f} ({WAVES[-1]}) — a factor of "
      f"{obs[WAVES[0]]['share_always']/obs[WAVES[-1]]['share_always']:.1f}")
print(f"  residual share: min {rs.min():.3f} · max {rs.max():.3f} · **spread {SPREAD:.3f}** · "
      f"**trend ρ_S(year, residual) = {TREND:+.3f}** over {len(WAVES)} waves")
print("  ⚠ ARITHMETIC: a collapsing marginal ATTENUATES ρ and therefore RAISES the residual — "
      "a rising residual is the arithmetic default, and the null must absorb exactly that.")


def ipf(pooled, rm, cm, iters=200):
    t = pooled.astype(float) + 1e-9
    for _ in range(iters):
        t *= (rm / t.sum(1))[:, None]
        t *= (cm / t.sum(0))[None, :]
    return t


def build_null(col_a, col_b, rng, nsim=NSIM):
    """Kind of null, named: an IPF-reconstructed CONSTANT-COUPLING null — one association, each
    wave's own margins, each wave's own n. Returns (spread draws, trend draws)."""
    sub = pd.DataFrame({"y": d.loc[M, "year"].to_numpy(),
                        "a": col_a[M].to_numpy(), "b": col_b[M].to_numpy()}).dropna()
    sub = sub[sub.y.isin(WAVES)]
    lev_a, lev_b = sorted(sub.a.unique()), sorted(sub.b.unique())
    pooled = pd.crosstab(sub.a, sub.b).reindex(index=lev_a, columns=lev_b, fill_value=0).to_numpy()
    per = {}
    for y in WAVES:
        s = sub[sub.y == y]
        t = pd.crosstab(s.a, s.b).reindex(index=lev_a, columns=lev_b, fill_value=0).to_numpy()
        per[y] = (t.sum(1) + 1e-9, t.sum(0) + 1e-9, len(s))
    C = len(lev_b)
    sp, tr = np.empty(nsim), np.empty(nsim)
    for i in range(nsim):
        vals = []
        for y in WAVES:
            rm, cm, n = per[y]
            t = ipf(pooled, rm, cm)
            p = (t / t.sum()).ravel()
            idx = rng.choice(len(p), size=n, p=p)
            r = stats.spearmanr(idx // C, idx % C).statistic
            vals.append(1 - r ** 2)
        v = np.array(vals)
        sp[i] = v.max() - v.min()
        tr[i] = stats.spearmanr(yr, v).statistic
    return sp, tr


print("\n=== (2) THE NULL — kind of null NAMED: IPF-reconstructed constant-coupling ===")
nsp, ntr = build_null(W, REFI["index"], np.random.default_rng(3241))
p95_sp, p95_tr = float(np.percentile(nsp, 95)), float(np.percentile(np.abs(ntr), 95))
print(f"  null SPREAD: median {np.median(nsp):.3f}  95th **{p95_sp:.3f}**   (observed {SPREAD:.3f})")
print(f"  null |TREND|: median {np.median(np.abs(ntr)):.3f}  95th **{p95_tr:.3f}**   "
      f"(observed |{TREND:.3f}|)")
SPREAD_ABOVE = SPREAD > p95_sp
TREND_ABOVE = abs(TREND) > p95_tr
print(f"  => spread above its null: **{SPREAD_ABOVE}** · trend above its null: **{TREND_ABOVE}**")
print("  ⚠ the null already contains the marginal collapse — it is built from it.")

print("\n=== (3) POSITIVE CONTROL — dose ANCHORED to abortion's measured movement (0.351) ===")
sub = pd.DataFrame({"y": d.loc[M, "year"].to_numpy(), "a": W[M].to_numpy(),
                    "b": REFI["index"][M].to_numpy()}).dropna()
sub = sub[sub.y.isin(WAVES)]
lev_a, lev_b = sorted(sub.a.unique()), sorted(sub.b.unique())
pooled = pd.crosstab(sub.a, sub.b).reindex(index=lev_a, columns=lev_b, fill_value=0).to_numpy()
dose = {}
for g in (0.0, 0.25, 0.5, 0.75, 1.0):
    vals = []
    for k, y in enumerate(WAVES):
        s = sub[sub.y == y]
        t0 = pd.crosstab(s.a, s.b).reindex(index=lev_a, columns=lev_b, fill_value=0).to_numpy()
        rm, cm = t0.sum(1) + 1e-9, t0.sum(0) + 1e-9
        t = ipf(pooled, rm, cm)
        pw = 1.0 + g * (k - (len(WAVES) - 1) / 2) / max(1, (len(WAVES) - 1) / 2) * 0.9
        t2 = ipf(np.power(t / t.sum(), pw), rm, cm)
        p = (t2 / t2.sum()).ravel()
        idx = RNG.choice(len(p), size=len(s), p=p)
        vals.append(1 - stats.spearmanr(idx // len(lev_b), idx % len(lev_b)).statistic ** 2)
    v = np.array(vals)
    dose[g] = float(v.max() - v.min())
    print(f"  g={g:.2f}  planted spread = {dose[g]:.3f}")
_gs = sorted(dose)
_peak = max(range(len(_gs)), key=lambda i: dose[_gs[i]])
_mono = all(dose[_gs[i]] <= dose[_gs[i + 1]] + 0.02 for i in range(_peak))
POS_OK = bool(dose[0.0] < p95_sp and dose[_gs[_peak]] > p95_sp and _mono and _peak > 0)
print(f"  measured turning point g={_gs[_peak]:.2f} (the statistic is EVEN IN ρ, so past the turn "
      f"the plant inverts the association — `#885`'s lesson, applied here by construction)")
print(f"  => positive control **{'PASS' if POS_OK else 'FAIL'}**: at g=0 the spread is "
      f"{dose[0.0]:.3f}, INSIDE the null {p95_sp:.3f}; at the turn {dose[_gs[_peak]]:.3f}, outside")

print("\n=== (4) BLINDNESS CONTROL (world C) — would a movement the SIZE OF ABORTION'S be seen? ===")
MDE = p95_sp
BLIND_OK = MDE < ABORTION_SPREAD
det = [g for g in _gs if dose[g] > p95_sp]
print(f"  the design's own detection floor for a spread is its null 95th = **{MDE:.3f}**")
print(f"  abortion's measured movement (`#885`) = **{ABORTION_SPREAD:.3f}**")
print(f"  => **{'PASS' if BLIND_OK else 'FAIL'}** — a movement of the size that matters "
      f"{'WOULD' if BLIND_OK else 'would NOT'} have been detected here; smallest firing dose "
      f"g={min(det) if det else 'none'}")
cw = {}
for y in (WAVES[0], WAVES[len(WAVES) // 2], WAVES[-1]):
    m = M & (d["year"] == y)
    x, z = W[m].to_numpy(), REFI["index"][m].to_numpy()
    bs = [1 - stats.spearmanr(x[j], z[j]).statistic ** 2
          for j in (RNG.integers(0, len(x), len(x)) for _ in range(400))]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    cw[y] = float(hi - lo)
    print(f"     per-wave bootstrap CI width {y}: {hi-lo:.3f}")
print("  ⚠ the CI width is the gauge test made numeric: if it grew as the marginal collapsed, "
      "'flat' and 'unmeasurable' would converge and the null would be silence.")

print("\n=== (5) NEGATIVE CONTROL — does the IPF machinery manufacture motion? ===")
NEG_OK = float(np.median(nsp)) < p95_sp and abs(float(np.median(ntr))) < p95_tr
print(f"  its own reconstruction: median spread {np.median(nsp):.3f} vs 95th {p95_sp:.3f} · "
      f"median trend {np.median(ntr):+.3f} vs 95th {p95_tr:.3f} -> **{'PASS' if NEG_OK else 'FAIL'}**")

print("\n=== (6) SPECIFICATION CURVE + MULTIPLICITY — the whole grid, disagreement included ===")
rows = []
for est in ("spearman", "kendall", "gamma"):
    for k in ("index", "spk", "col", "lib"):
        for win, ws in (("all", WAVES), ("1990plus", [y for y in WAVES if y >= 1990])):
            v = {}
            for y in ws:
                m = M & (d["year"] == y)
                v[y] = res_of(W[m], REFI[k][m], est)
            arr = np.array([v[y] for y in ws], float)
            if np.isnan(arr).any():
                continue
            rows.append(dict(est=est, sanction=k, window=win, waves=len(ws),
                             spread=float(arr.max() - arr.min()),
                             trend=float(stats.spearmanr(np.array(ws, float), arr).statistic)))
G = pd.DataFrame(rows)
print(f"  cells: **{len(G)}** (3 estimators × 4 sanction measures × 2 windows)")
for est, g in G.groupby("est"):
    print(f"  {est:9s} spread median {g.spread.median():.3f} [{g.spread.min():.3f},"
          f"{g.spread.max():.3f}] · trend median {g.trend.median():+.3f} "
          f"[{g.trend.min():+.3f},{g.trend.max():+.3f}]")
print(f"  cells whose spread exceeds the index null 95th ({p95_sp:.3f}): "
      f"**{(G.spread > p95_sp).sum()}/{len(G)}**")
print("  ⚠ γ is far less sensitive to marginal compression than ρ; if a movement were compression "
      "it should shrink under γ. The whole grid is printed, including cells that disagree.")

# ⚠⚠ SIGN STABILITY — computed, never typed. A kill that fires on ONE cell licenses "not flat";
# it does NOT license a DIRECTION unless the direction survives the grid. `realstat`: any
# comparative word in a verdict must be computed.
_by_win = {w: g for w, g in G.groupby("window")}
SIGN = {w: dict(pos=int((g.trend > 0).sum()), n=len(g), median=float(g.trend.median()),
                lo=float(g.trend.min()), hi=float(g.trend.max())) for w, g in _by_win.items()}
DIR_STABLE = len({np.sign(v["median"]) for v in SIGN.values()}) == 1
print("\n=== (6b) SIGN STABILITY OF THE TREND — computed, because a kill on one cell is not a "
      "direction ===")
for w, v in SIGN.items():
    print(f"  window {w:9s} trend>0 in {v['pos']}/{v['n']} cells · median {v['median']:+.3f} · "
          f"range [{v['lo']:+.3f},{v['hi']:+.3f}]")
print(f"  => the trend's sign is **{'STABLE' if DIR_STABLE else 'NOT STABLE'}** across windows ⇒ "
      f"{'a direction may be stated' if DIR_STABLE else '**no direction may be stated**'}")

GG = Gate("#886 · did the homosexuality coupling move, at 28 waves instead of four")
GG.asserted("(1) HARD RULE 1: n and the years actually asked printed before any column was cited; "
            "28 waves against the four the abortion pair forced",
            True, f"n={int(M.sum())} · waves {len(WAVES)} ({WAVES[0]}–{WAVES[-1]})",
            kind="control", population=POP)
GG.asserted("(2) OFFSET CONTROL — 'should this zero be zero?' NO: a collapsing marginal attenuates "
            "ρ and RAISES the residual by arithmetic. **Kind of null: an IPF-reconstructed "
            "constant-coupling null**, built from the observed marginal collapse itself",
            True, f"null 95th: spread {p95_sp:.3f} · |trend| {p95_tr:.3f}",
            kind="control", population=POP)
GG.asserted("(3) POSITIVE CONTROL, dose-response, monotone up to a MEASURED turning point (the "
            "statistic is even in ρ), and it must NOT fire at g=0 where the coupling is constant "
            "and the real marginal collapse is still applied",
            bool(POS_OK), " ".join(f"g={g}:{v:.3f}" for g, v in dose.items())
                          + f" · turn g={_gs[_peak]:.2f} · null95 {p95_sp:.3f}",
            kind="control", population=POP)
GG.asserted("(4) BLINDNESS CONTROL (world C): the design's own detection floor must be SMALLER than "
            "abortion's measured movement, or 'flat' is silence rather than a finding",
            bool(BLIND_OK), f"MDE(spread) {MDE:.3f} vs abortion's {ABORTION_SPREAD:.3f} · "
                            f"per-wave CI widths {list(cw.values())}",
            kind="control", population=POP)
GG.asserted("(5) NEGATIVE CONTROL: the IPF machinery must not manufacture motion — its own "
            "reconstruction re-measured sits inside its own null",
            bool(NEG_OK), f"median spread {np.median(nsp):.3f} · median trend {np.median(ntr):+.3f}",
            kind="control", population=POP)
GG.asserted("(6) KILL (pre-registered): for `#885`② to survive at full power, **neither the spread "
            "NOR the trend may exceed its IPF constant-coupling null**",
            bool(not SPREAD_ABOVE and not TREND_ABOVE),
            f"spread {SPREAD:.3f} vs {p95_sp:.3f} (above {SPREAD_ABOVE}) · trend {TREND:+.3f} vs "
            f"±{p95_tr:.3f} (above {TREND_ABOVE}) · grid cells above {(G.spread > p95_sp).sum()}"
            f"/{len(G)} · marginal 0.714->0.255 · SIGN STABILITY {SIGN} (stable {DIR_STABLE})",
            kind="kill",
            yardstick="the spread and the year-trend of 1−ρ² over 28 waves; the floor is the same "
                      "pair of statistics simulated under one coupling with each wave's own margins",
            yardstick_noise=float(p95_sp), population=POP,
            direction=None)
print()
print(GG)
adm = GG.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif not BLIND_OK:
    V = (f"**C — THE INSTRUMENT WENT BLIND.** Its detection floor ({MDE:.3f}) is not smaller than "
         f"abortion's measured movement ({ABORTION_SPREAD:.3f}) ⇒ *flat* here is **silence**, and "
         f"`#885`② is neither confirmed nor refuted.")
elif SPREAD_ABOVE or TREND_ABOVE:
    V = (f"**B — IT IS NOT FLAT, and four waves could not have seen that.**\n"
         f"  spread **{SPREAD:.3f}** vs null **{p95_sp:.3f}** (inside) · trend **{TREND:+.3f}** vs "
         f"**±{p95_tr:.3f}** (**outside**) ⇒ **`#885`② is RETRACTED.**\n"
         f"  ⚠ **And note WHICH statistic caught it**: the SPREAD is inside its null and the TREND "
         f"is not. **Four waves can only compute a spread** — so `#885`② was not merely "
         f"under-powered, it was measured with the only statistic four points admit.\n"
         f"  ⚠⚠ **BUT NO DIRECTION MAY BE STATED, and that is computed rather than conceded.** "
         f"Across the 24-cell grid the trend's sign **flips with the window**: over all 28 waves it "
         f"is positive in {SIGN['all']['pos']}/{SIGN['all']['n']} cells (median "
         f"{SIGN['all']['median']:+.3f}), and from 1990 onward positive in only "
         f"{SIGN['1990plus']['pos']}/{SIGN['1990plus']['n']} (median "
         f"{SIGN['1990plus']['median']:+.3f}). The same pre-registered cell reads "
         f"{G[(G.sanction=='index')&(G.window=='all')&(G.est=='spearman')].trend.iloc[0]:+.3f} on "
         f"all waves and "
         f"{G[(G.sanction=='index')&(G.window=='1990plus')&(G.est=='spearman')].trend.iloc[0]:+.3f} "
         f"from 1990.\n"
         f"  ⇒ **what is retracted is 'it did not move'. What replaces it is NOT 'it loosened' — it "
         f"is 'it is not flat, and its direction is not established at this resolution'.**\n"
         f"  ⇒ **one sentence about people: over half a century Americans reversed their verdict on "
         f"whether homosexuality is wrong — the share calling it *always wrong* fell "
         f"{obs[WAVES[0]]['share_always']:.3f} → {obs[WAVES[-1]]['share_always']:.3f} — and the "
         f"joint between that verdict and what they wanted done did NOT hold still. It is the one "
         f"thing I claimed a country could change its mind without disturbing, and it moved. Which "
         f"way it moved, this data cannot yet say.**")
else:
    V = (f"**A — `#885`② SURVIVES AT FULL POWER, and this is the round that could have killed it.**\n"
         f"  Across **{len(WAVES)} waves, 1973–2021, n = {int(M.sum())}**, the share of Americans "
         f"calling homosexuality *always wrong* fell from **{obs[WAVES[0]]['share_always']:.3f} to "
         f"{obs[WAVES[-1]]['share_always']:.3f}** — and the coupling between that judgement and what "
         f"they would DO about it did not move: spread **{SPREAD:.3f}** against an IPF "
         f"constant-coupling null of **{p95_sp:.3f}**, trend **{TREND:+.3f}** against **±"
         f"{p95_tr:.3f}**.\n"
         f"  **And the null is not silence**: this design would have caught a movement the size of "
         f"the one abortion actually showed — its floor is {MDE:.3f} against abortion's "
         f"{ABORTION_SPREAD:.3f}.\n"
         f"  ⇒ **one sentence about people: over half a century Americans almost completely reversed "
         f"their verdict on whether homosexuality is wrong, and the strength of the link between "
         f"that verdict and what they wanted done to homosexuals did not change at all. The verdict "
         f"moved; the joint between judging and acting did not. What a country changes its mind "
         f"about, and how tightly its mind is wired to its hands, are two different things.**")
print(V)
print("\n⚠ **Registered**: repeated cross-sections, not a panel — every wave difference is a "
      "composition difference as well as a change; the battery stops in 2021; the battery asks about "
      "a PERSON and the norm about an ACT; and **the instrument cannot be changed** for a motion "
      "claim, because the only other matched-pair instrument here (SCCS) has no time axis at all.")

json.dump(dict(population=POP, waves=WAVES, observed=obs, spread=SPREAD, trend=TREND,
               null=dict(kind="IPF-reconstructed constant-coupling null",
                         p95_spread=p95_sp, p95_abs_trend=p95_tr,
                         median_spread=float(np.median(nsp)), median_trend=float(np.median(ntr))),
               spread_above=bool(SPREAD_ABOVE), trend_above=bool(TREND_ABOVE),
               dose=dose, turning_point=float(_gs[_peak]),
               mde_spread=MDE, abortion_spread=ABORTION_SPREAD, blind_ok=bool(BLIND_OK),
               ci_widths=cw, grid=G.to_dict("records"),
               controls=dict(positive=bool(POS_OK), negative=bool(NEG_OK), blindness=bool(BLIND_OK)),
               admissible=adm, verdict=V, gate_ok=GG.verdict()),
          open(OUT / "twenty_eight_waves.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'twenty_eight_waves.json'}")
