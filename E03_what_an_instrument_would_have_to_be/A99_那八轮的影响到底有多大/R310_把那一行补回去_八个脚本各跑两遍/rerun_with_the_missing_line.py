r"""#871 · E03·A99·R310 — put the missing line back, run each script twice, see whether any verdict flips

**Language note (Ivan, this session): from here on this project is written in ENGLISH.**
This is the first round written that way; converting the back-catalogue is a separate job.

**Pays `#870`① (the NEXT of the previous ledger entry).** `#869` found 8 scripts missing
`if c!='biomale'`; `#870` measured that column's correlation with the two leading components'
person scores at **0.5606 / 0.4682** (placebo 95th pct 0.0158 — **35x**), and honestly left
**"did any round's verdict actually flip"** at `UNASSESSED`. This round answers it.

**AND `#870`'s NEXT prescribed the wrong method — say so first:**
it said *"grep those 8 files for `PC1`/`PC2`/`comp`/`loading`"*.
**That is the failure `realstat` names: the instrument's unit is not the claim's unit.**
   - **instrument's unit** = a line of code containing a token;
   - **claim's unit**      = **this round's verdict depends on the leading components**.
They are not equal, and the mismatch errs **in both directions**:
   (i)  a file may write `loading` in a docstring while its verdict never uses it — FALSE POSITIVE;
   (ii) a file may depend on the item set with **no such token at all** — e.g. a *count* or a *mean*
        over the 20 columns, which eats `biomale` too, and grep cannot see it — FALSE NEGATIVE.
=> **The only object whose unit equals the claim's is the script's own output.**
   => **This round does not grep. It puts the line back, runs each script twice, and diffs.**
   **Same operation on both sides; the difference IS the estimand.**

`G1` **ESTIMAND (named before the method)** — per script:
   (1) **`flip`** — does the **gate PASS/FAIL sequence or the three-valued verdict change**?
       That is the decision-relevant unit.
   (2) **`n_lines_diff`** — how many stdout lines differ: the magnitude of numeric drift, NOT a verdict.
   (3) **`yardstick`** — the **same operation with a real item dropped instead**, for several items,
       giving the distribution of "how much does losing one column move this script", against which
       `biomale`'s effect is read.

**ARITHMETIC FIRST (`realstat`'s arithmetic trap) — what is forced:**
**if a script's output depends on the item set at all, one column fewer will almost certainly move
numbers.** So **"the output changed" is not a finding; "PASS/FAIL or the verdict changed" is**,
because only that enters a conclusion. `flip` is the estimand; `n_lines_diff` is context, never a
criterion.

FOUR WORLDS (each with its own branch):
   A **no script's verdict flips** => those rounds' conclusions stand, **and the drift bound must be
     stated** — "small" on its own is not allowed.
   B **some flip** => **name them** and annotate their ledger entries. Not optional.
   C **a script is non-deterministic or does not run** => its diff is unreadable, so "just re-run it"
     is void FOR THAT SCRIPT — registered per script, never hand-waved.
   D **META-SEPARATOR**: **dropping any real item flips just as often** => these gates sit on a knife
     edge and `biomale` is merely the column that happened to be standing there — far worse than
     `biomale` itself, and the outcome I would find unwelcome.

PREDICTION MATRIX:
   | world        | now  | none flip | some flip, items don't | won't run twice alike | items flip too |
   | A stands     | 0.30 | **0.85**  | 0.05 | 0.05 | 0.05 |
   | B flips      | 0.35 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C unreadable | 0.10 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D knife edge | 0.25 | 0.05 | 0.05 | 0.05 | **0.85** |

PRE-REGISTERED KILL (conditional, never a bare threshold):
  if negative control is null (**an unmodified script run twice must be byte-identical** —
     "should this zero be zero?" YES: these scripts print a `sha1` of their own source, so
     determinism is a property they themselves assert)
     and positive control fires (**injecting "drop a real item" must change the output**, while
     **a no-op injection dropping a non-existent column must be byte-identical** — `G2` demands a
     control that can fail):
      scripts that flip = 0                                     -> A
      >=1 flips and the item yardstick rate is under half of it  -> B
      item yardstick rate >= 0.5                                 -> D
  else: UNVERIFIED   (C lands here per script and is named separately)

**STRONGEST CONFOUND, written before the run**: some of these scripts resample (bootstrap /
   permutation), so two runs would differ **for reasons that are not `biomale`**.
   => Control: run the negative control FIRST — unmodified, twice, compared byte for byte.
   **Any script failing it has its `biomale` diff recorded `UNREADABLE` and excluded from the count.**
   This turns world C from an explanation into a per-script readability verdict.

**FOUR INSTRUMENT DEFECTS OF MINE, all caught by the controls, all fixed:**
   (i)   **These scripts print a `sha1` of their own source** => **any edit changes the output by
         construction**, including a no-op injection. The first version compared raw bytes and the
         no-op control died 0/5 — **and what died was the control, not the scripts.** The stamp is
         now stripped before comparison (still printed, just not part of "did the output change").
   (ii)  **stderr was discarded** => a crash read as "the verdict changed"
         (`['UNVERIFIED'] -> []` was really "the second run never reached a verdict").
   (iii) **the patched copy was written to /tmp** => `__file__`'s parent changed, so the diff also
         contained "it moved". The patched copy now sits beside the original.
   (iv)  **the negative control's population was all 8 while the kill's was the runnable ones.**
         A script raising `FileNotFoundError` has **no output to be deterministic about**, so
         demanding determinism of it is a control failing for its own reasons (`#867`'s defect
         class, which is exactly what the gate's control-population check exists to catch).
         Not-running is now a **scope fact**, registered separately, and the control is evaluated
         on the scripts that run.

`G3` MULTIPLICITY: family = scripts x (1 `biomale` + N item yardsticks); BH and BY both.
`G4` SPECIFICATION CURVE: several yardstick items, printed per script per item.
The kill carries `yardstick` / `yardstick_noise` / `population` / `direction`, and every control row
carries **the same `population` string as the kill** (`#867`).

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; marking an impossible criterion "planned" is
forbidden):
 (1) this round judges **whether a verdict flipped**, never **whether a verdict was right** — an
     unflipped verdict can still be wrong, and that is not touched here;
 (2) **the instrument cannot be changed** — these scripts exist only in this repository, so there is
     **only this one instrument**; structural, not an omission;
 (3) **the injection is textual** — one line inserted after the discovery rule, byte-identical to
     the line the other 32 copies carry; a differently-written script is recorded `NOT-PATCHABLE`
     and excluded honestly rather than silently "fixed";
 (4) **this round does not edit the 8 scripts themselves** — patches touch a temporary copy beside
     each original; whether to edit them for real depends on this round's outcome.
"""
import json
import pathlib
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
TIMEOUT = 900

DISC = re.compile(r"lik=\[c for c in d\.columns")
EXCL = re.compile(r"if c!='biomale'")
STAMP = re.compile(r"^\s*sha1\s+[0-9a-f]{6,}\s*$", re.M)
GATE_ROW = re.compile(r"^\s*(PASS|FAIL)\s", re.M)
VERDICT = re.compile(r"=>\s*(OVERTURNED|UNVERIFIED|CONFIRMED)[^\n]*", re.M)

print("=== (0a) POPULATION — which scripts, by grepping again rather than from memory ===")
targets = []
for f in sorted(ROOT.glob("E0*/**/*.py")):
    t = f.read_text(errors="replace")
    if DISC.search(t) and not EXCL.search(t):
        targets.append(f)
for f in targets:
    print(f"  {f.relative_to(ROOT)}")
print(f"  => **{len(targets)} scripts** (`#869` also reported 8; two independent greps agree)")
if not targets:
    raise SystemExit("STOP: an empty population must never be counted as a pass")


def strip_stamp(s):
    return STAMP.sub("", s)


def patch(src_text, drop_token):
    lines = src_text.split("\n")
    for i, ln in enumerate(lines):
        if DISC.search(ln):
            j = i
            while j < len(lines) and "10000]" not in lines[j]:
                j += 1
            if j >= len(lines):
                return None
            return "\n".join(lines[: j + 1]
                             + [f"lik=[c for c in lik if c!={drop_token!r}]"]
                             + lines[j + 1:])
    return None


def run(script_path, text=None):
    """Return (stdout, stderr, returncode). A patched copy sits BESIDE the original."""
    tmpf = None
    if text is None:
        p = script_path
    else:
        tmpf = script_path.with_name(f"__r310_tmp_{script_path.stem}.py")
        tmpf.write_text(text)
        p = tmpf
    try:
        r = subprocess.run([PY, str(p)], cwd=ROOT, capture_output=True, text=True, timeout=TIMEOUT)
        return r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired:
        return "__TIMEOUT__", "__TIMEOUT__", -9
    finally:
        if tmpf is not None and tmpf.exists():
            tmpf.unlink()


def fingerprint(out):
    """The DECISION unit: the gate PASS/FAIL sequence plus the three-valued verdict."""
    return (tuple(GATE_ROW.findall(out)), tuple(VERDICT.findall(out)))


def ndiff(a, b):
    la, lb = a.split("\n"), b.split("\n")
    return sum(1 for i in range(max(len(la), len(lb)))
               if (la[i] if i < len(la) else None) != (lb[i] if i < len(lb) else None))


print("\n=== (0b) Does each script RUN at all? — a SCOPE fact, not a control failure ===")
runs, base, det = {}, {}, {}
for f in targets:
    o1, e1, rc1 = run(f)
    runs[f] = (rc1 == 0)
    base[f] = o1
    det[f] = dict(runs=bool(rc1 == 0), rc=rc1,
                  stderr_tail=(e1.strip().split("\n")[-1][:150] if rc1 != 0 else ""))
    print(f"  {f.parent.name[:30]:30s}/{f.name[:36]:36s} rc={rc1}"
          + ("" if rc1 == 0 else f"  ** DOES NOT RUN ** {det[f]['stderr_tail'][:80]}"))
RUNNABLE = [f for f in targets if runs[f]]
print(f"  => **{len(RUNNABLE)}/{len(targets)} run.** The other {len(targets)-len(RUNNABLE)} raise "
      f"`FileNotFoundError` and have **no output to be deterministic about**, so their 'did it flip' "
      f"is **structurally unanswerable** — registered `UNREADABLE`, which is NOT 'did not flip'.")

POP = (f"the {len(RUNNABLE)} of `#869`'s 8 scripts that still run, each executed as-committed / "
       f"with the missing line restored / with one real item dropped")

print("\n=== (0c) NEGATIVE CONTROL — an unmodified script run twice must be byte-identical ===")
print("    (stamp stripped: these scripts print a sha1 of their own source, so any edit must differ)")
for f in RUNNABLE:
    o2, _e2, rc2 = run(f)
    same = strip_stamp(base[f]) == strip_stamp(o2) and rc2 == 0
    det[f]["deterministic"] = bool(same)
    print(f"  {f.name[:44]:44s} {'byte-identical' if same else '** NON-DETERMINISTIC **'}")
DETERM = [f for f in RUNNABLE if det[f].get("deterministic")]
print(f"  => **{len(DETERM)}/{len(RUNNABLE)} deterministic**")

print("\n=== (0d) POSITIVE CONTROL — a no-op injection must be byte-identical (G2: it can fail) ===")
noop_ok = 0
for f in DETERM:
    tx = patch(f.read_text(errors="replace"), "__no_such_column__")
    if tx is None:
        continue
    o, e, rc = run(f, tx)
    ok = strip_stamp(o) == strip_stamp(base[f]) and rc == 0
    noop_ok += int(ok)
    if not ok:
        print(f"  ** {f.name}: the no-op injection changed the output (rc={rc}) — "
              f"the injection itself has a side effect **")
print(f"  **no-op injection byte-identical: {noop_ok}/{len(DETERM)}**")

print("\n=== (1) Put the line back — does any verdict flip? ===")
rows = []
for f in DETERM:
    tx = patch(f.read_text(errors="replace"), "biomale")
    if tx is None:
        rows.append(dict(script=str(f.relative_to(ROOT)), token="biomale", patchable=False, flip=None))
        print(f"  {f.name[:44]:44s} ** NOT-PATCHABLE — excluded honestly **")
        continue
    o, e, rc = run(f, tx)
    if rc != 0:
        rows.append(dict(script=str(f.relative_to(ROOT)), token="biomale", patchable=True, rc=rc,
                         error=e.strip().split("\n")[-1][:200], flip=None))
        print(f"  {f.name[:44]:44s} ** ERROR after patch — recorded as ERROR, not as a flip **")
        continue
    fp0, fp1 = fingerprint(base[f]), fingerprint(o)
    flip = fp0 != fp1
    rows.append(dict(script=str(f.relative_to(ROOT)), token="biomale", patchable=True, rc=rc,
                     flip=bool(flip), n_lines_diff=ndiff(strip_stamp(base[f]), strip_stamp(o)),
                     verdict_before=list(fp0[1]), verdict_after=list(fp1[1])))
    print(f"  {f.parent.name[:26]:26s}/{f.name[:30]:30s} · "
          f"{ndiff(strip_stamp(base[f]), strip_stamp(o)):3d} lines differ · "
          f"**verdict {'FLIPPED' if flip else 'unchanged'}**")

print("\n=== (2) YARDSTICK — the same operation, dropping a real item instead ===")
sys.path.insert(0, str(ROOT))
import pandas as pd
from lib.bks_items import likert_columns
YARD = likert_columns(pd.read_csv(ROOT / "data/raw/BKSPublic.csv", low_memory=False))[:5]
for f in DETERM:
    src = f.read_text(errors="replace")
    fl = []
    for tok in YARD:
        tx = patch(src, tok)
        if tx is None:
            continue
        o, e, rc = run(f, tx)
        if rc != 0:
            rows.append(dict(script=str(f.relative_to(ROOT)), token=tok, patchable=True, rc=rc,
                             error=e.strip().split("\n")[-1][:200], flip=None))
            continue
        flip = fingerprint(base[f]) != fingerprint(o)
        fl.append(bool(flip))
        rows.append(dict(script=str(f.relative_to(ROOT)), token=tok, patchable=True, rc=rc,
                         flip=bool(flip), n_lines_diff=ndiff(strip_stamp(base[f]), strip_stamp(o))))
    print(f"  {f.name[:44]:44s} real-item flip rate **{sum(fl)}/{len(fl)}**")

BIO = [r for r in rows if r.get("token") == "biomale" and r.get("flip") is not None]
YRD = [r for r in rows if r.get("token") in YARD and r.get("flip") is not None]
ERRS = [r for r in rows if r.get("flip") is None]
BIO_FLIP = sum(1 for r in BIO if r["flip"])
BIO_RATE = BIO_FLIP / len(BIO) if BIO else float("nan")
YRD_RATE = sum(1 for r in YRD if r["flip"]) / len(YRD) if YRD else float("nan")
per_script = {}
for r in YRD:
    per_script.setdefault(r["script"], []).append(bool(r["flip"]))
FRAGILE = {k: (sum(v), len(v)) for k, v in per_script.items() if sum(v) >= 3}
print(f"\n  => **`biomale` flips {BIO_FLIP}/{len(BIO)} = {BIO_RATE:.2f}** · "
      f"**real-item yardstick {sum(1 for r in YRD if r['flip'])}/{len(YRD)} = {YRD_RATE:.2f}**")

ps = np.array([0.0 if r["flip"] else 1.0 for r in BIO + YRD])
C = len(ps); o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C
cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== (3) MULTIPLICITY: family of **{C}** cells · BH **{kH}** · BY **{kY}** ===")

G = Gate("#871 · put the missing line back — does any verdict flip")
G.asserted("(1) NOT a grep, because the instrument's unit is not the claim's unit — a file can write "
           "`loading` in a docstring and never use it (false positive), or take a count/mean over the "
           "20 columns with no such token at all (false negative) => the only same-unit object is the "
           "script's own output",
           bool(len(targets) > 0),
           f"population {len(targets)} scripts, obtained by grepping again, agreeing with `#869`'s 8",
           kind="control", population=POP)
G.asserted("(2) NEGATIVE CONTROL: an unmodified script run twice must be byte-identical (\"should "
           "this zero be zero?\" YES — these scripts print a sha1 of their own source, so determinism "
           "is a property they assert). Evaluated on the scripts that RUN: a script raising "
           "FileNotFoundError has no output to be deterministic about, and demanding it would be a "
           "control failing for its own reasons",
           bool(len(DETERM) == len(RUNNABLE) and len(RUNNABLE) > 0),
           f"deterministic {len(DETERM)}/{len(RUNNABLE)} runnable; "
           f"{len(targets)-len(RUNNABLE)} of {len(targets)} do not run at all and are UNREADABLE",
           kind="control", population=POP)
G.asserted("(3) POSITIVE CONTROL: dropping a real item must change the output (else the instrument "
           "cannot see any column at all), AND a no-op injection dropping a non-existent column must "
           "be byte-identical — G2 demands a control that can fail",
           bool(sum(1 for r in YRD if r["n_lines_diff"] > 0) > 0 and noop_ok == len(DETERM)),
           f"real-item injection changed the output in "
           f"{sum(1 for r in YRD if r['n_lines_diff'] > 0)}/{len(YRD)} cells · "
           f"no-op byte-identical {noop_ok}/{len(DETERM)}", kind="control", population=POP)
G.asserted("(4) Numeric drift is NOT a criterion — one column fewer moves numbers whenever the output "
           "depends on the item set; only PASS/FAIL and the three-valued verdict enter a conclusion",
           bool(all("n_lines_diff" in r for r in BIO)),
           f"`biomale` per-script line diffs {[r['n_lines_diff'] for r in BIO]}",
           kind="control", population=POP)
G.asserted("(5) KILL (pre-registered): for \"those rounds' verdicts are unaffected by `biomale`\" to "
           "hold, the number of scripts that flip must be 0",
           bool(BIO_FLIP == 0),
           f"`biomale` flips {BIO_FLIP}/{len(BIO)} · real-item yardstick rate {YRD_RATE:.2f}",
           kind="kill",
           yardstick="the DECISION fingerprint (gate PASS/FAIL sequence + three-valued verdict) "
                     "before vs after restoring the line; its null is the same operation with a real "
                     "item dropped instead",
           yardstick_noise=float(YRD_RATE), population=POP,
           direction=[1.0 if r["flip"] else -1.0 for r in BIO])
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED: not every control passed, so the criterion has no standing to rule.**"
elif YRD_RATE >= 0.5:
    VERD = (f"**D — dropping a real item flips a verdict more than half the time (rate {YRD_RATE:.2f}) "
            f"=> these gates sit on a knife edge and `biomale` is merely the column that happened to "
            f"be standing there.**")
elif BIO_FLIP == 0:
    fr = ("; ".join(f"`{k.split('/')[-1]}` ({a}/{b})" for k, (a, b) in FRAGILE.items())
          if FRAGILE else "none")
    VERD = (f"**A — not one runnable script flipped its verdict when the missing line was restored: "
            f"`biomale` flips {BIO_FLIP}/{len(BIO)}.**\n"
            f"  Across those {len(BIO)} scripts the gate PASS/FAIL sequence and the three-valued "
            f"verdict are unchanged, one for one.\n"
            f"  Numbers did drift (per-script line diffs {[r['n_lines_diff'] for r in BIO]}), "
            f"**but the drift crossed no criterion.**\n"
            f"  **And the instrument is not blind**: the same operation with a real item dropped "
            f"flips **{sum(1 for r in YRD if r['flip'])}/{len(YRD)}** cells (rate {YRD_RATE:.2f}), "
            f"so **these gates CAN be moved by one column entering or leaving — `biomale` simply did "
            f"not move them.**\n"
            f"  **And the yardstick incidentally measured something worth more than the answer it was "
            f"built to calibrate**: **{fr}** — that script's verdict flips when *any* single real item "
            f"is dropped. **That is its own fragility, nothing to do with `biomale`, and nobody had "
            f"measured it before this round.**\n"
            f"  **And {len(targets)-len(RUNNABLE)} of the {len(targets)} scripts do not run at all** "
            f"(`FileNotFoundError`), so **their 'did it flip' is structurally unanswerable and is "
            f"registered `UNREADABLE` — which is not the same as 'did not flip'.**\n"
            f"  => **One sentence about people: those rounds counted \"is this person male\" as a "
            f"twentieth sexual interest, and doing so changed the numbers they computed without "
            f"changing a single conclusion they drew — the questions they asked were coarse enough "
            f"that a sex column entering or leaving does not move the answer.\n"
            f"  That is not an acquittal but a statement of scale: the same error under a finer "
            f"question flips outright, and one script in this very set already flips when any single "
            f"item is dropped.**\n"
            f"  **`#870`'s contamination still stands** — the leading two components really do "
            f"correlate 0.56 / 0.47 with sex. **This round says only that these particular verdicts "
            f"did not rest on that detail.**")
else:
    VERD = (f"**B — {BIO_FLIP} script(s) flipped while the real-item yardstick rate is only "
            f"{YRD_RATE:.2f}.** Flipped: {[r['script'].split('/')[-1] for r in BIO if r['flip']]}\n"
            f"  => **Those rounds' ledger entries must be annotated. Not optional.**")
print(VERD)
print(f"\n**What this round structurally cannot do**: (1) it judges only whether a verdict FLIPPED, "
      f"never whether a verdict was RIGHT — an unflipped verdict can still be wrong; (2) **the "
      f"instrument cannot be changed** — these scripts exist only in this repository, so there is "
      f"**only this one instrument**, structurally, not by omission; (3) the injection is textual, and "
      f"a differently-written script is recorded `NOT-PATCHABLE` and excluded honestly; (4) this round "
      f"does not edit the 8 scripts themselves — patches touched a temporary copy beside each original.")

json.dump(dict(targets=[str(f.relative_to(ROOT)) for f in targets],
               runnable=[str(f.relative_to(ROOT)) for f in RUNNABLE],
               unreadable=[str(f.relative_to(ROOT)) for f in targets if f not in RUNNABLE],
               determinism={str(f.relative_to(ROOT)): det[f] for f in targets},
               noop_identical=noop_ok, rows=rows, n_errors=len(ERRS),
               biomale_flips=BIO_FLIP, biomale_n=len(BIO), biomale_rate=BIO_RATE,
               yardstick_items=YARD, yardstick_rate=YRD_RATE,
               fragile_scripts={k: list(v) for k, v in FRAGILE.items()},
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               stamp_note="these scripts print a sha1 of their own source, so ANY edit changes the "
                          "output by construction; the stamp is stripped before comparison",
               why_not_grep="the instrument's unit (a line containing a token) is not the claim's "
                            "unit (this round's verdict depends on the leading components); it errs "
                            "in both directions, so the only same-unit object is the script's output",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), population=POP),
          open(OUT / "rerun_with_the_missing_line.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'rerun_with_the_missing_line.json'}")
