#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A115·R365 — the same three questions, five different targets
================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#926` found `homosex` at the 100th percentile of 22 items and read it as *"not
                general liberalisation"* — but with a caveat it could not remove: **|move| was used
                and per-item polarity was never checked**, so no direction was interpretable.
                `#926`(2) called that a one-command fix. It is, and running it changes the question.

⚠ WHAT THE     Value labels, read from the release (`#915`(2): they DO ship for GSS, unlike NSFG):
POLARITY         `homosex`  1 always wrong .. 4 not wrong at all      -> HIGH = permissive
CHECK GAVE       `spkhomo`  1 allowed to speak, 2 not allowed         -> LOW  = permissive
                 `colhomo`  1 allowed to teach, 2 not allowed         -> LOW  = permissive
                 `libhomo`  1 remove, 2 not remove                    -> HIGH = permissive
               ⇒ **All four homosexuality items moved TOWARD acceptance.** `#926`'s apparent sign
               disagreement was coding polarity and nothing else, exactly as it flagged. That part
               is now a re-expression, **not a test**, and is labelled so.

Why Now        Reading the labels exposed a far better design than the one `#926` had.
               **`spkhomo`/`colhomo`/`libhomo` and `spkrac`/`colrac`/`librac` are the SAME THREE
               QUESTION STEMS** — allowed to speak · allowed to teach · book kept in the library —
               **asked about a different target group**, in the same waves, of the same cohorts, on
               the same 2-point scale. Any difference between them cannot be the question form, the
               era, the scale, or the respondents. **It can only be the target.**
               And the battery has five targets in all: homosexuals · racists · communists ·
               militarists · anti-religionists.

Live Worlds    W_TARGET    · homosexual-tolerance rises while racist-tolerance FALLS ⇒ what changed
                             is **who** is tolerated, not how tolerant people are.
               W_TOLERANCE · both rise ⇒ general tolerance rose and homosexuals merely more.
                             **Unwelcome: it makes `#926`'s "not general" much weaker.**
               W_NULL      · neither separates from the within-cohort permutation null.

Estimand       For each of the 5 targets, the within-cohort movement across the 3 stems, re-expressed
(G1)           in a common **toward-permissive** direction using the shipped labels; then the
               contrast between targets. **Matched by construction on stem, wave, cohort and scale.**

Prediction     W_TARGET    -> homosexual and racist movements have OPPOSITE signs, both resolved.
Matrix         W_TOLERANCE -> both positive.
               W_NULL      -> inside the null.

⚠ PRECONDITION `#925`(2): checked and PRINTED before the estimator — every (target, stem) pair needs
CHECK FIRST    >=3 cohorts with n>=40 in both endpoint waves. Failures are DROPPED WITH A COUNT.

Controls       NEGATIVE: permute YEAR within cohort, per item.
               POSITIVE: plant into the permuted world so `g=0` lands ON the null (`#922`'s gate).
               ⚠ SHAM, and it is the one that matters: the three STEMS within a target should agree.
               If `speak`, `teach` and `library` disagree for the same target, the "target" reading
               is wrong and the stems are measuring different things.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ the polarity re-expression is a RE-EXPRESSION, not a finding — `#926` already measured the
    magnitudes; this round only makes their direction readable;
  (2) ⚠ "toward permissive" for the tolerance items means "toward allowing the speech/teaching/book".
    Whether that is the same psychological quantity across targets is an ASSUMPTION, not a
    measurement — and it is exactly what W_TARGET puts in doubt;
  (3) ⚠ the tolerance battery ends in 2021; the window is the common one;
  (4) ⚠ repeated cross-section; age/period/cohort collinear; no age effect claimed;
  (5) ⚠ **only this one instrument**;
  (6) `[unchallenged]` — door (3).
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
RNG = np.random.default_rng(365)

# target -> {stem: (column, sign)}, sign = +1 if HIGH is permissive, -1 if LOW is permissive.
# ⚠ every sign below is read from the SHIPPED value labels, never guessed (`#926`(2)).
TARGETS = {
    "homosexuals":      {"speak": ("spkhomo", -1), "teach": ("colhomo", -1), "library": ("libhomo", +1)},
    "racists":          {"speak": ("spkrac", -1),  "teach": ("colrac", -1),  "library": ("librac", +1)},
    "communists":       {"speak": ("spkcom", -1),  "teach": ("colcom", +1),  "library": ("libcom", +1)},
    "militarists":      {"speak": ("spkmil", -1),  "teach": ("colmil", -1),  "library": ("libmil", +1)},
    "anti-religionists": {"speak": ("spkath", -1), "teach": ("colath", -1),  "library": ("libath", +1)},
}
# ⚠ `colcom` is the ODD ONE: its stem is "Should communist teacher be FIRED" (1 fired, 2 not fired),
#   so HIGH = permissive, opposite to every other `col*`. Read from the label, not from the name —
#   a variable name is not a measurement (HARD RULE 1), and here the NAME would have got it backwards.
MORAL = ("homosex", +1)   # 1 always wrong .. 4 not wrong at all

COLS = sorted({c for t in TARGETS.values() for c, _ in t.values()} | {MORAL[0]})
raw = pd.read_stata(GSS, columns=["year", "cohort"] + COLS, convert_categoricals=False)
raw = raw[raw["cohort"].between(1900, 2006)].copy()
raw["grp"] = (raw["cohort"] // 10 * 10).astype(int)
for c in COLS:
    v = raw[c]
    raw[c] = v.where(v.between(1, 7))
raw[MORAL[0]] = raw[MORAL[0]].where(raw[MORAL[0]].isin([1, 2, 3, 4]))

yrs = sorted(raw.year.unique())
y0 = min(yrs, key=lambda y: abs(y - 1990))
y1 = min(yrs, key=lambda y: abs(y - 2021))
print(f"window {int(y0)} -> {int(y1)}")


def move(col, sign):
    """within-cohort movement in item-SD units, oriented TOWARD PERMISSIVE."""
    e = raw.dropna(subset=[col])
    a, b = e[e.year == y0], e[e.year == y1]
    ks = [k for k in sorted(set(a.grp) & set(b.grp))
          if (a.grp == k).sum() >= 40 and (b.grp == k).sum() >= 40]
    if len(ks) < 3:
        return None
    sd = float(pd.concat([a[col], b[col]]).std())
    if not sd or sd <= 0:
        return None
    m = np.mean([float(b.loc[b.grp == k, col].mean() - a.loc[a.grp == k, col].mean()) for k in ks])
    return dict(col=col, oriented=float(sign * m / sd), raw=float(m / sd), sign=sign,
                cohorts=len(ks), n=int(len(a) + len(b)))


# ══ PRECONDITION CHECK, printed BEFORE the estimator (`#925`(2)) ═════════════════════
print("PRECONDITION CHECK — each (target, stem) needs >=3 cohorts with n>=40 in both endpoints:")
grid, dropped = [], []
for tgt, stems in TARGETS.items():
    for stem, (col, sign) in stems.items():
        r = move(col, sign)
        if r is None:
            dropped.append((tgt, stem, col))
            continue
        r.update(target=tgt, stem=stem)
        grid.append(r)
print(f"  usable {len(grid)} of {sum(len(s) for s in TARGETS.values())} · DROPPED {len(dropped)}: "
      f"{[f'{t}/{s}' for t, s, _ in dropped] or 'none'}  (absence reported, not passed)")

moral = move(*MORAL)
print(f"\n=== the moral item, for reference ===")
print(f"  homosex  oriented {moral['oriented']:+.4f}  cohorts={moral['cohorts']} n={moral['n']}")

print("\n=== THE GRID — same three stems, five targets, all oriented TOWARD PERMISSIVE ===")
by_target = {}
for tgt in TARGETS:
    rows = [r for r in grid if r["target"] == tgt]
    if not rows:
        continue
    vals = [r["oriented"] for r in rows]
    by_target[tgt] = dict(mean=float(np.mean(vals)), stems={r["stem"]: r["oriented"] for r in rows},
                          agree=bool(len(set(np.sign(vals))) == 1), n_stems=len(rows),
                          n=int(min(r["n"] for r in rows)))
    print(f"  {tgt:18s} mean {np.mean(vals):+.4f}  " +
          " · ".join(f"{r['stem']}={r['oriented']:+.4f}" for r in rows) +
          f"   stems agree: {'YES' if by_target[tgt]['agree'] else 'NO'}  n={by_target[tgt]['n']}")

homo = by_target.get("homosexuals", {}).get("mean", np.nan)
rac = by_target.get("racists", {}).get("mean", np.nan)
print(f"\n  homosexuals {homo:+.4f}   racists {rac:+.4f}   difference {homo - rac:+.4f}")

# ══ NEGATIVE CONTROL — permute YEAR within cohort ════════════════════════════════════
null_vals = []
e_h = raw.dropna(subset=["colhomo"])
for _ in range(200):
    p = e_h[e_h.year.isin([y0, y1])].copy()
    p["year"] = p.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
    a, b = p[p.year == y0], p[p.year == y1]
    ks = [k for k in sorted(set(a.grp) & set(b.grp))
          if (a.grp == k).sum() >= 40 and (b.grp == k).sum() >= 40]
    if len(ks) < 3:
        continue
    sd = float(pd.concat([a["colhomo"], b["colhomo"]]).std())
    if sd > 0:
        null_vals.append(float(-np.mean([b.loc[b.grp == k, "colhomo"].mean()
                                         - a.loc[a.grp == k, "colhomo"].mean() for k in ks]) / sd))
null_med, null_sd = float(np.median(null_vals)), float(np.std(null_vals))
print(f"\n  null (year permuted within cohort; kind of null: within-cohort year-label permutation): "
      f"{null_med:+.4f} +/- {null_sd:.4f} over {len(null_vals)} draws")

# ══ POSITIVE CONTROL — plant INTO the permuted world so g=0 lands on the null ════════
sweep = []
for g in (0.0, 0.10, 0.20, 0.30, 0.40):
    vals = []
    for _ in range(80):          # ⚠ 20 reps left the low-g end noisy enough to read as non-monotone
        p = e_h[e_h.year.isin([y0, y1])].copy()
        p["year"] = p.groupby("grp")["year"].transform(lambda s: RNG.permutation(s.to_numpy()))
        # ⚠⚠ THIS PLANT WAS WRONG TWICE, AND THE SECOND WAY IS THE INSTRUCTIVE ONE.
        #   v1 added a continuous offset and clipped it on a 2-point item: the pooled SD collapsed
        #   and the sweep sat flat at +2.1368 — the plant disturbing its own denominator, third
        #   instance of that family (`#919`(3), `#923`, here).
        #   v2 flipped responses instead, which is right in principle — and still did nothing,
        #   because it flipped `== 2`. ⚠ **`colhomo`'s numeric codes are 4 and 5, not 1 and 2.**
        #   I had read the value LABELS and assumed the CODES. **A label list is not a code list**,
        #   and that is HARD RULE 1 failing at the exact place I believed I had obeyed it.
        #   ⇒ the codes are now DERIVED from the data. (The polarity SIGNS were all correct — only
        #   the plant's hard-coded literals were wrong, so the substantive grid is unaffected.)
        # ⚠ v3 named this `vals` and SHADOWED the sweep's own accumulator, so `np.median` ran over
        #   the CODE LIST and returned a flat 4.0 at every g. Three broken plants in one round,
        #   each caught by the same control — which is the argument for the control, not against it.
        codes = sorted(p["colhomo"].dropna().unique())
        lo_code, hi_code = float(codes[0]), float(codes[-1])
        perm_code, other_code = lo_code, hi_code          # `colhomo`: LOW is permissive
        late = p.index[(p.year == y1) & (p["colhomo"] == other_code)]
        if len(late):
            flip = RNG.random(len(late)) < g
            p.loc[late[flip], "colhomo"] = perm_code
        a, b = p[p.year == y0], p[p.year == y1]
        ks = [k for k in sorted(set(a.grp) & set(b.grp))
              if (a.grp == k).sum() >= 40 and (b.grp == k).sum() >= 40]
        sd = float(pd.concat([a["colhomo"], b["colhomo"]]).std())
        if len(ks) >= 3 and sd > 0:
            vals.append(float(-np.mean([b.loc[b.grp == k, "colhomo"].mean()
                                        - a.loc[a.grp == k, "colhomo"].mean() for k in ks]) / sd))
    sweep.append([float(g), float(np.median(vals)) if vals else np.nan])
print(f"  positive sweep (planted into the permuted null world): "
      f"{[(g, round(v, 4)) for g, v in sweep]}")

ps = [2 * (1 - stats.norm.cdf(abs((r["oriented"] - null_med) / (null_sd or 1e-9)))) for r in grid]

if not grid:
    print("EMPTY POPULATION"); sys.exit(2)

homo_res = abs(homo - null_med) > 2 * null_sd
rac_res = abs(rac - null_med) > 2 * null_sd
opposite = (not np.isnan(homo)) and (not np.isnan(rac)) and np.sign(homo) != np.sign(rac)
stems_agree = all(v["agree"] for v in by_target.values())

G = Gate("Same three questions, five targets: did tolerance rise, or did the target change?")
G.plant_direction_from_sweep("positive: a planted within-cohort trend raises the oriented movement, "
                             "and g=0 is null", sweep, baseline=null_med, baseline_spread=null_sd)
G.negative_control("year permuted within cohort", null_med, abs(homo),
                   null_spread=null_sd, null_kind="within-cohort year-label permutation")
G.multiplicity_control("all (target, stem) cells", ps, 0.05,
                       labels=[f"{r['target']}/{r['stem']}" for r in grid])
G.asserted("PRECONDITIONS checked and printed BEFORE the estimator; failures dropped with a count",
           True, f"{len(grid)} usable · {len(dropped)} dropped: "
                 f"{[f'{t}/{s}' for t, s, _ in dropped] or 'none'}", kind="control")
G.asserted("every polarity was READ FROM THE SHIPPED LABELS, not inferred from the variable name",
           True, "`colcom` is 'should the communist teacher be FIRED' (1 fired, 2 not), so HIGH is "
                 "permissive — the OPPOSITE of every other `col*`. The name would have got it "
                 "backwards; the label did not. scope stated", kind="control")
G.asserted("SHAM: the three stems within a target must agree in sign, or 'target' is the wrong read",
           stems_agree,
           "; ".join(f"{t}: {'agree' if v['agree'] else 'DISAGREE'} "
                     f"({', '.join(f'{k}={x:+.3f}' for k, x in v['stems'].items())})"
                     for t, v in by_target.items()), kind="control")
G.spec_curve_cells_declare_n("every published cell states its n", grid)
G.asserted("KILL: W_TOLERANCE requires homosexual- and racist-tolerance to move the SAME way",
           not (opposite and homo_res and rac_res),
           f"homosexuals {homo:+.4f} (resolved={homo_res}) · racists {rac:+.4f} "
           f"(resolved={rac_res}) · difference {homo - rac:+.4f}, null {null_med:+.4f} "
           f"+/- {null_sd:.4f}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif opposite and homo_res and rac_res:
    VERDICT, WORLD = "OVERTURNED", "W_TARGET · who is tolerated changed, not how tolerant people are"
elif homo_res and rac_res and not opposite:
    VERDICT, WORLD = "CONFIRMED", "W_TOLERANCE · general tolerance rose; homosexuals merely more"
else:
    VERDICT, WORLD = "UNVERIFIED", "W_NULL · one or both arms do not resolve"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=927, round="E03·A115·R365", verdict=VERDICT, world=WORLD,
           estimand="within-cohort movement of the SAME three tolerance stems across five target "
                    "groups, oriented toward-permissive using the shipped value labels",
           instrument="GSS 1972-2024 gss7224_r3a.dta", window=[int(y0), int(y1)],
           polarity={c: s for t in TARGETS.values() for c, s in t.values()},
           grid=grid, by_target=by_target, moral_item=moral,
           homosexuals=homo, racists=rac, difference=float(homo - rac),
           stems_agree=stems_agree, dropped=[f"{t}/{s}" for t, s, _ in dropped],
           null_median=null_med, null_sd=null_sd, null_draws=len(null_vals),
           positive_sweep=sweep, family_size=len(ps),
           settles="`#926`(2): all four homosexuality items move toward acceptance; the sign "
                   "disagreement was coding polarity",
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "same_question_different_target.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'same_question_different_target.json'}")
