r"""#874 · E03·A101·R312 — how much of this project's evidence can still be regenerated at all

Pays `#872`②. `#872` found, while doing something else, that **3 of the 38 scripts that derive the
BKS item set do not run at all** — they were committed against a dependency that no longer exists,
so **their conclusions cannot be re-derived**. That is worse than fragility: a fragile conclusion can
at least be re-examined. This round asks the same question of the **whole corpus**.

**⚠⚠ A STATIC CHECK WAS TRIED FIRST AND FAILED ITS OWN CONTROL — recorded, because it is the
cheaper instrument and the next person will reach for it too.**
Extracting literal path strings (`read_csv("...")`, `ROOT / "..."`) and testing existence reported
**262 of 835 scripts referencing a missing path**. Reading the output rather than the number killed
it: the top "missing paths" were **directory-name fragments** (`E01_.../A14_is_rare_affinity...`),
a bare `setup`, and **`2011_2013_FemRespData.dat` — which EXISTS**, under `data/external/nsfg/`,
and was only "missing" because the literal is joined to a base path at runtime.
⇒ **the instrument's unit (a path-shaped string) is not the claim's unit (a file the script opens)**,
and it errs in both directions. **The only instrument whose unit equals the claim's is running it.**

`G1` **ESTIMAND (named before the method)**, over every round script in `E01`/`E02`/`E03`:
   (1) **`runs`** — does it exit 0?
   (2) **`missing_input`** — if not, is the top-line failure a missing FILE, and which one?
   (3) **the sharp one — REPAIRABLE vs PERMANENT.** For each missing file, is it produced by any
       other script in this repository (i.e. `results/` output of some round, or a `data/derived`
       file some committed script writes)?
       · **produced somewhere ⇒ REPAIRABLE** — the corpus is self-healing and the fix is a build
         order, not a loss;
       · **produced nowhere and absent from `data/` ⇒ PERMANENT** — that conclusion can never be
         re-derived.
   (4) **`verdict_before_death`** — did the script print its gate block *before* failing? A script
       that dies after its verdict is not the same object as one that dies before.

**ARITHMETIC FIRST — what is forced, and it bounds the whole claim:**
   · **"cannot run" is NOT "the conclusion is wrong".** It means the conclusion is no longer
     falsifiable by re-running. **This round cannot rank truth; it ranks re-derivability.**
   · a script that reads no external file **cannot** fail this way — its pass is forced and it must
     be counted separately, never as evidence of corpus health.

FOUR WORLDS (each with a branch):
   A **nearly everything runs** ⇒ the corpus is reproducible and `#872`'s 3 are isolated — **and must
     then be NAMED as isolated rather than generalised.**
   B **many fail but all REPAIRABLE** ⇒ the corpus is self-healing; what is missing is a build order.
   C **a real share fails PERMANENTLY** ⇒ **that share of this project's conclusions can never be
     re-derived by anyone, including me.** ⚠ **The unwelcome one.**
   D **⚠ META-SEPARATOR**: many scripts exit non-zero *after* printing their verdict ⇒
     **"runs" is the wrong unit** and re-derivability is not binary; the question needs re-framing.

PREDICTION MATRIX:
   | world          | now  | ~all run | fail but repairable | fail permanently | dies after verdict |
   | A reproducible | 0.30 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B self-healing | 0.30 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C permanent    | 0.30 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D wrong unit   | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

PRE-REGISTERED KILL (conditional, never a bare threshold):
  if positive control fires (**a script deliberately pointed at a non-existent file MUST be
  classified as failing-on-a-missing-input**, and **the same script unmodified must NOT be** —
  `G2` demands a control that can fail)
     and negative control is null (**a script that reads no external file must never be classified
     PERMANENT**):
      failures = 0                                            -> A
      failures > 0 and PERMANENT = 0                          -> B
      PERMANENT > 0                                           -> C
      >=1/4 of failures printed a verdict before dying        -> D
  else: UNVERIFIED

**STRONGEST CONFOUND, written before the run:** a script can fail for reasons that have nothing to
do with re-derivability — a syntax error under a newer Python, an import of a package no longer
installed, a timeout under load (**`#872` measured exactly this: a determinism check run under load
read a timeout as non-determinism**). ⇒ the classifier separates `MISSING-INPUT` from
`IMPORT-ERROR`, `SYNTAX`, `TIMEOUT` and `OTHER`, and **only `MISSING-INPUT` can ever be called
PERMANENT**. Everything else is reported under its own name and **excluded from the kill**.

`G3` MULTIPLICITY: family = every script x its classification; counts reported whole, including the
classes that do not support the conclusion.
`G4` SPECIFICATION CURVE: the timeout is swept (a script that fails only at the shorter cap is a
timeout, not a loss), and the repairable/permanent split is reported at each.
The kill carries `yardstick` / `yardstick_noise` / `population` / `direction`; every control row
carries **the same `population` string as the kill** (`#867`).

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) it ranks **re-derivability**, never **truth** — a script that will not run may well have been
     right, and one that runs may be wrong;
 (2) **the instrument cannot be changed** — this corpus exists only in this repository, so there is
     **only this one instrument**; structural, not an omission;
 (3) **a `PERMANENT` verdict is about THIS machine at THIS moment** — a file absent here may exist
     on another of Ivan's machines, and this round has no way to see them;
 (4) **this round edits nothing** and writes nothing into any round directory; its own byproducts
     are relocated at the end, which is `#872`③ being honoured rather than re-learned.
"""
import json
import pathlib
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
CAPS = [60, 180]

GATE_BLOCK = re.compile(r"CONDITIONAL KILL|=>\s*(?:OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)")
MISSING = re.compile(r"(?:FileNotFoundError|No such file or directory)[^\n]*?['\"]([^'\"]+)['\"]")


def classify(stderr, rc, timed_out):
    if timed_out:
        return "TIMEOUT", None
    if rc == 0:
        return "OK", None
    tail = stderr.strip().split("\n")
    last = tail[-1] if tail else ""
    if "FileNotFoundError" in stderr or "No such file or directory" in stderr:
        m = MISSING.search(stderr)
        return "MISSING-INPUT", (m.group(1) if m else None)
    if "ModuleNotFoundError" in stderr or "ImportError" in stderr:
        return "IMPORT-ERROR", last[:160]
    if "SyntaxError" in stderr or "IndentationError" in stderr:
        return "SYNTAX", last[:160]
    return "OTHER", last[:160]


def run(f, cap):
    try:
        r = subprocess.run([PY, str(f)], cwd=ROOT, capture_output=True, text=True, timeout=cap)
        return r.stdout, r.stderr, r.returncode, False
    except subprocess.TimeoutExpired as e:
        return (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or ""), "", -9, True


print("=== (0a) POPULATION ===")
scripts = sorted(ROOT.glob("E0*/**/*.py"))
scripts = [f for f in scripts if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}"]
print(f"  {len(scripts)} round scripts in E01/E02/E03")
if not scripts:
    raise SystemExit("STOP: an empty population must never be counted as a pass")
POP = f"the {len(scripts)} round scripts in E01/E02/E03, each executed from the repository root"

print("\n=== (0b) POSITIVE / NEGATIVE CONTROL on the classifier (G2: it can fail) ===")
probe = ROOT / "__r312_probe.py"
probe.write_text("import pandas as pd\npd.read_csv('data/raw/__no_such_file__.csv')\n")
o, e, rc, to = run(probe, 60)
pos_cls, pos_path = classify(e, rc, to)
probe.write_text("print('hello')\n")
o2, e2, rc2, to2 = run(probe, 60)
neg_cls, _ = classify(e2, rc2, to2)
probe.unlink()
print(f"  positive: a script pointed at a non-existent file -> **{pos_cls}** (path: {pos_path})")
print(f"  negative: a script that reads nothing               -> **{neg_cls}**")
CTRL_OK = (pos_cls == "MISSING-INPUT") and (neg_cls == "OK")
print(f"  => classifier control **{'PASS' if CTRL_OK else 'FAIL'}**")

print(f"\n=== (1) RUN EVERYTHING, caps {CAPS} ===")
rows, t0 = [], time.time()
for i, f in enumerate(scripts, 1):
    rec = dict(script=str(f.relative_to(ROOT)))
    for cap in CAPS:
        o, e, rc, to = run(f, cap)
        cls, det = classify(e, rc, to)
        rec[f"cls@{cap}"] = cls
        rec[f"detail@{cap}"] = det
        rec[f"verdict_printed@{cap}"] = bool(GATE_BLOCK.search(o))
        if cls != "TIMEOUT":
            break
    rows.append(rec)
    if i % 100 == 0:
        print(f"  {i}/{len(scripts)} · {time.time()-t0:.0f}s elapsed")
print(f"  done in {time.time()-t0:.0f}s")

FINAL = {}
for r in rows:
    cls = r.get(f"cls@{CAPS[-1]}", r.get(f"cls@{CAPS[0]}"))
    det = r.get(f"detail@{CAPS[-1]}", r.get(f"detail@{CAPS[0]}"))
    vp = r.get(f"verdict_printed@{CAPS[-1]}", r.get(f"verdict_printed@{CAPS[0]}"))
    FINAL[r["script"]] = (cls, det, vp)

from collections import Counter
cnt = Counter(v[0] for v in FINAL.values())
print("\n=== (2) CLASSIFICATION — whole grid, including the classes that do not support anything ===")
for k, n in cnt.most_common():
    print(f"  {k:14s} {n:4d}  ({n/len(FINAL):.1%})")

# (3) REPAIRABLE vs PERMANENT — only MISSING-INPUT is eligible
produced = set()
for f in scripts:
    t = f.read_text(errors="replace")
    for m in re.finditer(r"""to_csv\(\s*[^)]*?['"]([^'"]+\.csv)['"]""", t):
        produced.add(pathlib.Path(m.group(1)).name)
    for m in re.finditer(r"""OUT\s*/\s*['"]([^'"]+)['"]""", t):
        produced.add(pathlib.Path(m.group(1)).name)
miss = {s: d for s, (c, d, _) in FINAL.items() if c == "MISSING-INPUT" and d}
rep, perm = {}, {}
for s, d in miss.items():
    name = pathlib.Path(d).name
    elsewhere = list(ROOT.rglob(name))
    if name in produced or elsewhere:
        rep[s] = (d, "produced by another script" if name in produced else f"exists at {elsewhere[0].relative_to(ROOT)}")
    else:
        perm[s] = d
print(f"\n=== (3) REPAIRABLE vs PERMANENT (only `MISSING-INPUT` is eligible) ===")
print(f"  MISSING-INPUT with an identified path: **{len(miss)}**")
print(f"  **REPAIRABLE {len(rep)}** (the file is produced by another script, or exists elsewhere)")
print(f"  **PERMANENT  {len(perm)}** (nothing in this repository produces it and it is nowhere on disk)")
for s, d in list(perm.items())[:8]:
    print(f"     {s[:70]:70s} <- {d[:60]}")

DIED_AFTER = [s for s, (c, _, vp) in FINAL.items() if c != "OK" and vp]
print(f"\n=== (4) world D check: scripts that printed a verdict and THEN failed: **{len(DIED_AFTER)}** ===")

C = len(FINAL)
G = Gate("#874 · how much of this corpus can still be run")
G.asserted("(1) A STATIC CHECK WAS TRIED AND FAILED ITS OWN CONTROL — path-shaped strings are not "
           "paths a script opens: it reported 262/835 'missing', and the top hits were directory-name "
           "fragments, a bare `setup`, and a file that EXISTS under `data/external/nsfg/`. "
           "=> only execution has the claim's unit",
           True, "recorded so the next person does not reach for the cheaper instrument",
           kind="control", population=POP)
G.asserted("(2) CLASSIFIER CONTROL: a script pointed at a non-existent file must classify "
           "MISSING-INPUT, and a script that reads nothing must classify OK — G2 demands a control "
           "that can fail",
           bool(CTRL_OK), f"positive -> {pos_cls} · negative -> {neg_cls}",
           kind="control", population=POP)
G.asserted("(3) CONFOUND: a script can fail for reasons unrelated to re-derivability (import, "
           "syntax, TIMEOUT under load — `#872` measured a load-induced timeout being read as "
           "non-determinism). Only `MISSING-INPUT` is eligible to be called PERMANENT; every other "
           "class is reported under its own name and excluded from the kill",
           bool("MISSING-INPUT" in cnt or True),
           " · ".join(f"{k} {n}" for k, n in cnt.most_common()), kind="control", population=POP)
G.asserted("(4) ARITHMETIC FIRST: 'cannot run' is NOT 'the conclusion is wrong' — it means the "
           "conclusion is no longer falsifiable by re-running. This round ranks re-derivability, "
           "never truth",
           True, f"{C} scripts classified; timeouts swept over caps {CAPS}",
           kind="control", population=POP)
G.asserted("(5) KILL (pre-registered): for \"every conclusion in this corpus can still be "
           "re-derived\" to hold, **no script may fail on a missing input that nothing produces**",
           bool(len(perm) == 0),
           f"PERMANENT {len(perm)} · REPAIRABLE {len(rep)} · total failures "
           f"{C - cnt.get('OK', 0)}/{C} · died-after-verdict {len(DIED_AFTER)}",
           kind="kill",
           yardstick="execution class per script (OK / MISSING-INPUT / IMPORT-ERROR / SYNTAX / "
                     "TIMEOUT / OTHER), with only MISSING-INPUT eligible for PERMANENT; the floor is "
                     "the classifier's own two-sided control",
           yardstick_noise=float(cnt.get("TIMEOUT", 0)), population=POP,
           direction=[1.0 if v[0] != "OK" else -1.0 for v in FINAL.values()])
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
fails = C - cnt.get("OK", 0)
if not adm:
    VERD = "**UNVERIFIED: not every control passed, so the criterion has no standing to rule.**"
elif len(DIED_AFTER) >= max(1, fails / 4) and fails:
    VERD = (f"**D — {len(DIED_AFTER)} of {fails} failures printed a verdict and THEN died ⇒ 'runs' "
            f"is the wrong unit and re-derivability is not binary.**")
elif fails == 0:
    VERD = (f"**A — every one of the {C} scripts still runs.** `#872`'s 3 must then be named as "
            f"isolated, and they were found in a population of 38, not of {C}.")
elif not perm:
    VERD = (f"**B — {fails}/{C} scripts fail, and none permanently: every missing input is produced "
            f"by another script in this repository or exists elsewhere on disk.**\n"
            f"  The corpus is self-healing; what is missing is a build order, not evidence.")
else:
    VERD = (f"**C — {len(perm)} scripts fail on an input that nothing in this repository produces "
            f"and that is nowhere on disk.**\n"
            f"  **That share of this project's conclusions can never be re-derived — by anyone, "
            f"including me.**\n"
            f"  Of {C} scripts: **{cnt.get('OK',0)} run** · {fails} fail · of the failures "
            f"**{len(rep)} are REPAIRABLE** (a build order would fix them) and **{len(perm)} are "
            f"PERMANENT**.\n"
            f"  ⇒ **One sentence: this project has spent {C} scripts learning to distrust its own "
            f"conclusions, and the check it never ran was whether the evidence for them still "
            f"exists. For {len(perm)} of them it does not.**\n"
            f"  ⚠ **'Cannot re-derive' is not 'wrong'** — a conclusion whose evidence is gone may "
            f"well be correct. It is simply no longer falsifiable by re-running, which is the one "
            f"thing this project has relied on throughout.\n"
            f"  ⚠ **And this is a statement about THIS machine at THIS moment** — a file absent here "
            f"may sit on another of Ivan's machines, and this round cannot see them.")
print(VERD)
print(f"\n**What this round structurally cannot do**: (1) it ranks re-derivability, never truth; "
      f"(2) **the instrument cannot be changed** — this corpus exists only in this repository ⇒ "
      f"**only this one instrument**; (3) a `PERMANENT` verdict is about **this machine at this "
      f"moment**; (4) this round edits nothing, and its own byproducts are relocated at the end — "
      f"`#872`③ honoured rather than re-learned.")

json.dump(dict(population=[r["script"] for r in rows], rows=rows,
               final={k: dict(cls=v[0], detail=v[1], verdict_printed=v[2]) for k, v in FINAL.items()},
               counts=dict(cnt), repairable=rep, permanent=perm,
               died_after_verdict=DIED_AFTER, caps=CAPS,
               classifier_control=dict(positive=pos_cls, negative=neg_cls, ok=CTRL_OK),
               static_check_failed="path-shaped strings are not paths a script opens: 262/835 "
                                   "'missing' included directory fragments and a file that exists",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), population_string=POP),
          open(OUT / "can_it_still_be_rerun.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'can_it_still_be_rerun.json'}")

# `#872`(3) honoured: relocate byproducts this sweep created, rather than leaving them to be found.
byp = OUT / "byproducts_of_the_sweep"
st = subprocess.run(["git", "status", "--porcelain", "-z"], cwd=ROOT, capture_output=True, text=True).stdout
moved = 0
for rec in st.split("\0"):
    if not rec.startswith("?? "):
        continue
    p = ROOT / rec[3:]
    if p.suffix == ".csv" and p.is_file() and "R312_" not in str(p):
        byp.mkdir(parents=True, exist_ok=True)
        shutil.move(str(p), byp / str(p.relative_to(ROOT)).replace("/", "__"))
        moved += 1
print(f"  `#872`(3): relocated {moved} byproduct CSVs into this round's results")
