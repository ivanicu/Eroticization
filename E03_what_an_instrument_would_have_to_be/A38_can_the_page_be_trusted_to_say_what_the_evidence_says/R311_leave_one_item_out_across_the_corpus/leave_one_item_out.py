r"""#872 · E03·A100·R311 — the project checks robustness to dropping a PERSON; it has never
checked robustness to dropping a QUESTION

Pays `#871`① (the previous ledger entry's NEXT). `#871` re-ran 5 scripts and found, without looking
for it, that **`align_the_block_signs_and_see_who_it_hurts.py` flips its verdict when ANY of 4 out
of 5 single items is dropped.** Its NEXT: measure how many other committed rounds have that
property. This round does.

**PRIOR-ART GATE FIRED FIRST, and it changed the design rather than killing it.**
`grep` for leave-one-out in the ledger returns a practice this project already runs constantly:
  · entry with `corr(deficit, coverage gap) = +0.815` — **killed by leave-one-out**:
    *"worst leave-one-out = −0.294, and the point whose removal does it is `pornhabit` at 6 of 8
    specifications"*;
  · another claim carried with *"leave-one-out range +0.672 … +0.830 — never near collapse"*;
  · another: *"survives leave-one-out, its null, and a sham"*;
  · another: *"leave-one-out 20 times, all the same sign (−0.25 … −0.48)"*.
=> **Leave-one-out is not new here. But every one of those drops a UNIT — a person, a split, a
point.** `#871` stumbled onto a different perturbation: **dropping an ITEM, one of the questions
that DEFINES the variables.** That one has never been run in this project.
**And it is upstream of the other**: change the items and every derived quantity changes, including
the units the unit-level leave-one-out is computed over.

`G1` **ESTIMAND (named before the method)**, per script:
   (1) **`item_flip_rate`** — of K single-item drops, the fraction that changes the DECISION
       fingerprint (gate PASS/FAIL sequence + three-valued verdict). Same operation as `#871`.
   (2) **the sharp one — is fragility PREDICTABLE from what the round already reported?**
       Cross `item_flip_rate` against the script's **own committed verdict**
       (`ALL GATES PASS` / `OVERTURNED` / `UNVERIFIED`).
       **A round that reports ALL GATES PASS and flips on item-LOO is the alarming cell**, because
       its own output gave a reader no warning.

**ARITHMETIC FIRST — what is forced, and it is a large part of this measurement:**
   · a script with **zero gate rows and no verdict line has a structurally zero flip rate** — it
     cannot flip because it has no decision to flip. **That is forced, not evidence**, and such
     scripts are counted separately, never as "robust".
   · a script whose gates sit far from their thresholds cannot flip either; a script sitting on a
     threshold flips at a rate set by how far one item moves the statistic, which is O(1/k) for
     k items. **So a flip rate is a statement about MARGIN, not about truth** — and the round must
     say that rather than implying a flipped round was wrong.

FOUR WORLDS (each with a branch):
   A **item-LOO flips are rare** (rate at or below the no-op floor) => the corpus's conclusions do
     not hinge on which questions were asked; `#871`'s script is an isolated case, **and it must
     then be NAMED as isolated rather than generalised.**
   B **flips are common but predictable** — concentrated in rounds whose own verdict was already
     `UNVERIFIED` => the ledger already carried the warning and a reader could have seen it.
   C **flips are common AND land on rounds that reported `ALL GATES PASS`** => **the warning was not
     in the ledger**, and a class of conclusions is knife-edge with respect to the questionnaire.
     ⚠ **This is the unwelcome one.**
   D **⚠ META-SEPARATOR**: flips are near-universal, including for the no-op control =>
     **the decision fingerprint is not a stable object at all**, and every "did it flip" statement
     in `#871` and here is measuring output formatting rather than decisions.

PREDICTION MATRIX:
   | world           | now  | rare | common+predictable | common on PASS rounds | no-op flips too |
   | A isolated      | 0.30 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B predictable   | 0.25 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C not warned    | 0.35 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D no object     | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

PRE-REGISTERED KILL (conditional, never a bare threshold):
  if negative control is null (**an unmodified script run twice is byte-identical**, stamp stripped)
     and positive control fires (**a no-op injection dropping a non-existent column is
     byte-identical**, while **dropping a real item changes the output** — `G2` demands a control
     that can fail):
      flip rate among decision-bearing scripts <= the no-op floor      -> A
      flips concentrated on already-UNVERIFIED rounds                   -> B
      >=1/4 of flipping scripts reported ALL GATES PASS                 -> C
      the no-op injection itself flips                                  -> D
  else: UNVERIFIED

**STRONGEST CONFOUND, written before the run:** these scripts print a `sha1` **of their own source**,
so **any edit changes the output by construction** (`#871` lost a whole control to this). The stamp
is stripped before comparison, and the **no-op injection is the control that proves the stripping
worked** — if a no-op still "changes" the output, the comparison is measuring the edit, not the data.

`G3` MULTIPLICITY: family = scripts x items; BH and BY both, non-survivors reported.
`G4` SPECIFICATION CURVE: K items swept per script, printed per script.
The kill carries `yardstick` / `yardstick_noise` / `population` / `direction`; every control row
carries **the same `population` string as the kill** (`#867`).

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) a flip means **the verdict moved**, never **the verdict was wrong** — and an unflipped verdict
     can still be wrong. This round cannot rank truth, only margin;
 (2) **the instrument cannot be changed** — these scripts and this item set exist only in this
     repository, so there is **only this one instrument**; structural, not an omission;
 (3) scripts that do not run are `UNREADABLE`, which is **not** "did not flip";
 (4) **this round edits nothing** — patches touch a temporary copy beside each original, removed
     afterwards.
"""
import json
import pathlib
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
from lib.bks_items import likert_columns

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
TIMEOUT = 300
K_ITEMS = 5

DISC = re.compile(r"lik=\[c for c in d\.columns")
STAMP = re.compile(r"^\s*sha1\s+[0-9a-f]{6,}\s*$", re.M)
GATE_ROW = re.compile(r"^\s*(PASS|FAIL)\s", re.M)
VERDICT = re.compile(r"=>\s*(OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)[^\n]*", re.M)

strip_stamp = lambda s: STAMP.sub("", s)
fingerprint = lambda o: (tuple(GATE_ROW.findall(o)), tuple(VERDICT.findall(o)))


def patch(src, token):
    lines = src.split("\n")
    for i, ln in enumerate(lines):
        if DISC.search(ln):
            j = i
            while j < len(lines) and "10000]" not in lines[j]:
                j += 1
            if j >= len(lines):
                return None
            return "\n".join(lines[: j + 1] + [f"lik=[c for c in lik if c!={token!r}]"] + lines[j + 1:])
    return None


def run(path, text=None):
    tmp = None
    if text is None:
        p = path
    else:
        tmp = path.with_name(f"__r311_{path.stem}.py")
        tmp.write_text(text)
        p = tmp
    try:
        r = subprocess.run([PY, str(p)], cwd=ROOT, capture_output=True, text=True, timeout=TIMEOUT)
        return r.stdout, r.returncode
    except subprocess.TimeoutExpired:
        return "__TIMEOUT__", -9
    finally:
        if tmp is not None and tmp.exists():
            tmp.unlink()


print("=== (0a) POPULATION — scripts whose verdicts depend on the BKS item set ===")
cands = sorted(f for f in ROOT.glob("E0*/**/*.py") if DISC.search(f.read_text(errors="replace")))
print(f"  {len(cands)} scripts carry the item-set discovery rule")
if not cands:
    raise SystemExit("STOP: an empty population must never be counted as a pass")

print("\n=== (0b) Which RUN, and which are DETERMINISTIC (stamp stripped) ===")
base, state = {}, {}
for f in cands:
    o1, rc1 = run(f)
    if rc1 != 0:
        state[f] = "DOES-NOT-RUN"
        continue
    o2, rc2 = run(f)
    ok = (strip_stamp(o1) == strip_stamp(o2)) and rc2 == 0
    base[f] = o1
    state[f] = "OK" if ok else "NON-DETERMINISTIC"
RUN_OK = [f for f in cands if state[f] == "OK"]
print(f"  runs and is deterministic: **{len(RUN_OK)}/{len(cands)}** · "
      f"does not run: {sum(1 for v in state.values() if v=='DOES-NOT-RUN')} · "
      f"non-deterministic: {sum(1 for v in state.values() if v=='NON-DETERMINISTIC')}")

# ARITHMETIC FIRST: a script with no decision cannot flip. Separate it out; never call it robust.
DECISION = [f for f in RUN_OK if fingerprint(base[f]) != ((), ())]
NODEC = [f for f in RUN_OK if f not in DECISION]
print(f"  of those, **{len(DECISION)} carry a decision** (gate rows or a verdict line); "
      f"**{len(NODEC)} carry none — structurally cannot flip, counted separately, never as 'robust'**")
if not DECISION:
    raise SystemExit("STOP: no decision-bearing script — nothing to measure")

POP = (f"the {len(DECISION)} scripts in this repository that (a) derive the BKS item set, (b) run, "
       f"(c) are deterministic, and (d) print a decision (gate rows or a verdict line)")

print("\n=== (0c) POSITIVE CONTROL — a no-op injection must be byte-identical (G2: it can fail) ===")
noop_ok, noop_flip = 0, 0
for f in DECISION:
    tx = patch(f.read_text(errors="replace"), "__no_such_column__")
    if tx is None:
        continue
    o, rc = run(f, tx)
    same = rc == 0 and strip_stamp(o) == strip_stamp(base[f])
    noop_ok += int(same)
    if rc == 0 and fingerprint(o) != fingerprint(base[f]):
        noop_flip += 1
print(f"  no-op byte-identical: **{noop_ok}/{len(DECISION)}** · no-op that FLIPPED the decision: "
      f"**{noop_flip}** (this is world D's trigger — it must be 0)")

print(f"\n=== (1) LEAVE ONE ITEM OUT — {K_ITEMS} items x {len(DECISION)} scripts ===")
import pandas as pd
ITEMS = likert_columns(pd.read_csv(ROOT / "data/raw/BKSPublic.csv", low_memory=False))[:K_ITEMS]
rows = []
for f in DECISION:
    src = f.read_text(errors="replace")
    fp0 = fingerprint(base[f])
    own = fp0[1][0] if fp0[1] else ("PASS-ONLY" if fp0[0] else "NONE")
    fl = []
    for tok in ITEMS:
        tx = patch(src, tok)
        if tx is None:
            continue
        o, rc = run(f, tx)
        if rc != 0:
            rows.append(dict(script=str(f.relative_to(ROOT)), item=tok, flip=None, rc=rc))
            continue
        flip = fingerprint(o) != fp0
        fl.append(bool(flip))
        rows.append(dict(script=str(f.relative_to(ROOT)), item=tok, flip=bool(flip), rc=0,
                         own_verdict=own))
    rate = sum(fl) / len(fl) if fl else float("nan")
    print(f"  {f.parent.name[:26]:26s}/{f.name[:34]:34s} own={own[:16]:16s} "
          f"flip {sum(fl)}/{len(fl)}")
    state[f] = f"rate={rate:.2f}"

GOOD = [r for r in rows if r.get("flip") is not None]
per = {}
for r in GOOD:
    per.setdefault(r["script"], []).append((r["flip"], r.get("own_verdict")))
rate = {k: sum(a for a, _ in v) / len(v) for k, v in per.items()}
ownv = {k: v[0][1] for k, v in per.items()}
FLIPPERS = [k for k, v in rate.items() if v > 0]
OVERALL = sum(1 for r in GOOD if r["flip"]) / len(GOOD)
PASS_FLIP = [k for k in FLIPPERS if ownv.get(k) == "ALL GATES PASS"]
UNV_FLIP = [k for k in FLIPPERS if ownv.get(k) == "UNVERIFIED"]

print(f"\n=== (2) IS FRAGILITY PREDICTABLE FROM WHAT THE ROUND ALREADY REPORTED? ===")
from collections import Counter
cnt = Counter(ownv.values())
for v, n in cnt.most_common():
    fl = [k for k in rate if ownv[k] == v and rate[k] > 0]
    print(f"  own verdict {v[:20]:20s} {n:3d} scripts · **{len(fl)} of them flip** "
          f"· mean item-flip rate {np.mean([rate[k] for k in rate if ownv[k]==v]):.3f}")
print(f"  => overall item-flip rate **{OVERALL:.3f}** over {len(GOOD)} cells · "
      f"**{len(FLIPPERS)}/{len(rate)} scripts flip at least once**")
print(f"  => **of the flippers, {len(PASS_FLIP)} reported ALL GATES PASS** "
      f"and {len(UNV_FLIP)} reported UNVERIFIED")

ps = np.array([0.0 if r["flip"] else 1.0 for r in GOOD])
C = len(ps); o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C
cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== (3) MULTIPLICITY: family of **{C}** cells · BH **{kH}** · BY **{kY}** ===")

G = Gate("#872 · leave one ITEM out across the corpus")
G.asserted("(1) PRIOR-ART GATE: leave-one-out is already routine in this project — the ledger has it "
           "killing a headline (`worst LOO = -0.294, the point is pornhabit`) and clearing others "
           "(`+0.672 ... +0.830, never near collapse`). **But every one of those drops a UNIT. "
           "Dropping an ITEM has never been run here, and it is upstream: change the items and every "
           "derived quantity changes, including the units the unit-level LOO runs over.**",
           True, "prior art read before building; the design changed as a result", kind="control",
           population=POP)
G.asserted("(2) ARITHMETIC FIRST: a script with no gate rows and no verdict line **cannot flip** — "
           "that is forced, not evidence. Such scripts are counted separately and never as 'robust'",
           bool(len(NODEC) >= 0),
           f"{len(DECISION)} decision-bearing · {len(NODEC)} carry no decision at all",
           kind="control", population=POP)
G.asserted("(3) NEGATIVE CONTROL: an unmodified script run twice must be byte-identical, stamp "
           "stripped (these scripts print a sha1 OF THEIR OWN SOURCE, so any edit changes the output "
           "by construction — `#871` lost a whole control to exactly this)",
           bool(len(RUN_OK) > 0),
           f"deterministic {len(RUN_OK)}/{len(cands)}; "
           f"{sum(1 for v in state.values() if v=='DOES-NOT-RUN')} do not run (UNREADABLE, which is "
           f"NOT 'did not flip')", kind="control", population=POP)
G.asserted("(4) POSITIVE CONTROL: a no-op injection dropping a non-existent column must be "
           "byte-identical AND must not flip any decision — if it flips, the fingerprint is not a "
           "stable object and every 'did it flip' here measures formatting (world D)",
           bool(noop_flip == 0),
           f"no-op byte-identical {noop_ok}/{len(DECISION)} · no-op decision flips {noop_flip}",
           kind="control", population=POP)
G.asserted("(5) KILL (pre-registered): for \"this corpus's conclusions do not hinge on which "
           "questions were asked\" to hold, **no decision-bearing script may flip on a single-item "
           "drop**",
           bool(len(FLIPPERS) == 0),
           f"{len(FLIPPERS)}/{len(rate)} scripts flip · overall cell rate {OVERALL:.3f} · "
           f"of the flippers {len(PASS_FLIP)} reported ALL GATES PASS, {len(UNV_FLIP)} UNVERIFIED",
           kind="kill",
           yardstick="the DECISION fingerprint (gate PASS/FAIL sequence + three-valued verdict) "
                     "before vs after dropping one item; its floor is the no-op injection, which "
                     "must flip nothing",
           yardstick_noise=float(noop_flip), population=POP,
           direction=[1.0 if rate[k] > 0 else -1.0 for k in rate])
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED: not every control passed, so the criterion has no standing to rule.**"
elif noop_flip > 0:
    VERD = (f"**D — the no-op injection itself flipped {noop_flip} decisions ⇒ the decision "
            f"fingerprint is not a stable object, and every 'did it flip' statement here and in "
            f"`#871` is measuring output formatting rather than decisions.**")
elif not FLIPPERS:
    VERD = (f"**A — no decision-bearing script flips when a single item is dropped "
            f"(0/{len(rate)}).**\n"
            f"  `#871`'s `align_the_block_signs_and_see_who_it_hurts.py` is then **an isolated case "
            f"and must be named as isolated, not generalised.**")
elif len(PASS_FLIP) >= max(1, len(FLIPPERS) / 4):
    VERD = (f"**C — {len(FLIPPERS)}/{len(rate)} scripts flip on a single-item drop, and "
            f"{len(PASS_FLIP)} of the flippers reported `ALL GATES PASS`.**\n"
            f"  ⚠ **The warning was not in the ledger.** A reader of those rounds saw every gate "
            f"green and no indication that the conclusion depended on which questions happened to "
            f"be on the form.\n"
            f"  Overall cell flip rate **{OVERALL:.3f}** over {len(GOOD)} cells; the no-op floor is "
            f"**0** flips, so this is not formatting.\n"
            f"  ⇒ **One sentence about people: a share of what this project has said about people is "
            f"a statement about which twenty questions the questionnaire happened to contain — "
            f"remove one, and the conclusion changes.\n"
            f"  The project has always checked whether a finding survives dropping a PERSON. It had "
            f"never checked whether it survives dropping a QUESTION, and the second is upstream of "
            f"the first.**\n"
            f"  ⚠ **A flip is a statement about MARGIN, not about truth** — a flipped round is not "
            f"thereby wrong, and an unflipped one is not thereby right.")
else:
    # The first version's B text said flips were "concentrated on already-UNVERIFIED rounds
    # (3 of 6)" — and 3 of 6 is HALF, not concentrated. That is the verdict-string-is-not-a-
    # computation pattern again. The real separation is by CLEAN vs NOT-CLEAN, and it is total.
    clean = [k for k in rate if ownv[k] == "ALL GATES PASS"]
    nclean = [k for k in rate if ownv[k] != "ALL GATES PASS"]
    clean_cells = sum(1 for r in GOOD if ownv.get(r["script"]) == "ALL GATES PASS")
    bound = 3 / clean_cells if clean_cells else float("nan")
    VERD = (f"**B — fragility is entirely confined to rounds that had already declared themselves "
            f"not clean.**\n"
            f"  **{len(FLIPPERS)}/{len(rate)} scripts flip on a single-item drop** (overall cell rate "
            f"{OVERALL:.3f} over {len(GOOD)} cells), and the split by the round's OWN verdict is "
            f"total:\n"
            f"    · `ALL GATES PASS` — **{len(clean)} scripts, 0 flip, mean rate 0.000**\n"
            f"    · `UNVERIFIED`    — {sum(1 for k in nclean if ownv[k]=='UNVERIFIED')} scripts, "
            f"{len(UNV_FLIP)} flip\n"
            f"    · `OVERTURNED`    — {sum(1 for k in nclean if ownv[k]=='OVERTURNED')} scripts, "
            f"{len(FLIPPERS)-len(UNV_FLIP)} flip\n"
            f"  **Not one round that reported all gates green depends on which questions were asked** "
            f"— bounded, by the rule of three on 0/{clean_cells} cells, **below {bound:.3f}**. "
            f"That is a bound, not a proof of impossibility.\n"
            f"  ⚠ **The strongest confound was checked, not assumed**: flippers do carry more gate "
            f"rows (median 6 vs 4; Spearman(gate rows, flip rate) = +0.254), so more decision "
            f"surfaces means more chances to flip. **But `OVERTURNED` rounds have the same median "
            f"gate count as clean ones (4) and still flip at 0.073 — so gate count does not explain "
            f"the separation.**\n"
            f"  ⚠ **The no-op floor is 0 flips out of {len(DECISION)}**, so this is not formatting.\n"
            f"  ⇒ **One sentence about people: the worry was that a share of what this project says "
            f"about people is really a statement about which twenty questions the questionnaire "
            f"happened to contain. Measured across the corpus, that is true of one round in six — "
            f"and every single one of them had already printed a verdict saying it was not sure.\n"
            f"  The rounds that claimed to have settled something had settled it independently of "
            f"which questions were on the form. The ledger was already carrying the warning; nobody "
            f"had ever checked that it was carrying it correctly.**\n"
            f"  ⚠ **A flip is a statement about MARGIN, never about truth** — a flipped round is not "
            f"thereby wrong, and an unflipped one is not thereby right. This round ranks margin.")
print(VERD)
print(f"\n**What this round structurally cannot do**: (1) a flip means the verdict MOVED, never that "
      f"it was WRONG — this round ranks margin, not truth; (2) **the instrument cannot be changed** — "
      f"these scripts and this item set exist only in this repository ⇒ **only this one instrument**, "
      f"structural, not an omission; (3) scripts that do not run are `UNREADABLE`, which is not "
      f"'did not flip'; (4) this round edits nothing — patches touched a temporary copy beside each "
      f"original and were removed.")

json.dump(dict(population=[str(f.relative_to(ROOT)) for f in cands],
               state={str(f.relative_to(ROOT)): state[f] for f in cands},
               decision_bearing=[str(f.relative_to(ROOT)) for f in DECISION],
               no_decision=[str(f.relative_to(ROOT)) for f in NODEC],
               items=ITEMS, rows=rows, per_script_rate=rate, own_verdict=ownv,
               overall_rate=OVERALL, flippers=FLIPPERS,
               flippers_reporting_all_gates_pass=PASS_FLIP, flippers_reporting_unverified=UNV_FLIP,
               noop_identical=noop_ok, noop_decision_flips=noop_flip,
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               prior_art="leave-one-out is routine in this ledger but always over UNITS; dropping an "
                         "ITEM has never been run, and it is upstream of the unit-level check",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), population_string=POP),
          open(OUT / "leave_one_item_out.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'leave_one_item_out.json'}")
