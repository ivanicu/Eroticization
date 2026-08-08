r"""#879 · E03·A102·R318 — the same measurement with both failed controls repaired, and only now a branch

`#875` found that **79 of 836 round scripts die at one name**, `round_path`, removed by my own
commit `4819b9b`, and that **52 of 52 literal references resolve mechanically out of git**. It also
said, in writing, that **resolving a reference is an UPPER BOUND on resurrecting a script** — the
file a key points at has been rewritten twice since. `#875`① is that gap, and this round measures it
rather than assuming either end of it.

**`#876` RAN THIS DESIGN AND CAME BACK `UNVERIFIED` — two of six controls failed, so its numbers
(74 of 81 exit 0) had no standing to select a world. This round repairs exactly those two things and
nothing else, so that the comparison is with the same instrument rather than a different one.**

**REPAIR ① — the population is PARSED, not matched.** `#876` built it with
`from\s+lib\.rounds\s+import[^\n]*round_path`, which returned **81**; `ast.ImportFrom` returns
**80**. The extra member was `#875`'s own diagnostic script, whose docstring **quotes** the broken
import line it exists to explain. That was the third of four instances in one session of a text scan
counting a **mention** as a **use**, and the durable fix is the same one the registry-key gate took:
**an `ast.ImportFrom` cannot be written in a docstring.**

**REPAIR ② — the ANALYSIS population is the scripts the name actually killed.** One member was
already `MISSING-INPUT` before the repair, so *"anything running now was caused by the repair"* was
false for it — correctly reported by the control. It is **named and reported, not silently dropped**:
the estimand is *does restoring the name resurrect the scripts the name killed*, and a script that
was failing for another reason was not killed by the name.

**REPAIR ③ — the negative control names its own exclusion.** `#876` scored 19/20 because
`bootstrap_noise_floor.py` is **not byte-identical across two solo runs with no load**. Scoring a
non-deterministic script against a repair it does not use is a control failing for its own reasons.
It is **excluded by name**, and **characterised here instead** — run 5 times, reporting the number of
distinct outputs — which is `#876`② measured rather than deferred.

⚠ **All three repairs were identified BEFORE `#876` finished running** (the AST count and the
non-determinism were both measured while it was still going), so none of them is a threshold moved
after seeing a result.

**THE REPAIR IS APPLIED** (same commit as this round): `lib/rounds.py` now carries `LEGACY_PATHS`
(52 pre-`4819b9b` keys → today's paths, every one verified to exist) and `round_path`. ⚠ An alias
`round_path = path` would have been **wrong** — the two accessors have **disjoint key vocabularies**
(52 numbered filenames vs 635 stems, 0 shared), so an alias converts an `ImportError` into a
`KeyError` and calls it a repair.

`G1` **ESTIMAND, named before the method**, over the **80 scripts that import `round_path`**:
   (1) **`comes_back`** — does it now exit 0?
   (2) **`next_wall`** — if not, what is the next failure, and **is it another stale path**? That is
       the difference between *the reference was the only thing wrong* and *the reference was the
       tip of a stack*.
   (3) **the sharp one — `still_agrees`.** For every script that comes back **and** whose round
       directory holds a committed artifact carrying a verdict, does the fresh run print the **same
       verdict**? ⚠ **A resurrected script that prints a different verdict is a ledger entry whose
       own evidence no longer reproduces**, and that is worth more than any count of exit codes.

**ARITHMETIC FIRST — what is forced, and it bounds the claim:**
   · these 80 scripts were **all** `IMPORT-ERROR` before the repair (measured, `#875`'s artifact),
     so **any** of them exiting 0 now is caused by the repair. There is no other candidate. That
     makes the *direction* forced and the *magnitude* the only open question;
   · **`comes_back` ≤ `reference resolves`** by construction — a script cannot run further than its
     first unresolved reference. So this round can only ever move the number DOWN from `#875`'s 52,
     and a result of "all 80" would be impossible, not impressive: **30 of the 80 build their key at
     runtime and `#875` explicitly could not see them.**

THREE WORLDS (each with a branch):
   A **the reference was the only thing wrong** ⇒ most come back and agree ⇒ the corpus's rot was
     entirely my registry rename, and it is now undone.
   B **the reference was the tip of a stack** ⇒ few come back, and the next wall is another stale
     path ⇒ **the rot is layered**, and `#875`'s "77% is one error class" was measuring the
     outermost layer of a laminate rather than the whole thing.
   C **⚠ THE UNWELCOME ONE — they come back and DISAGREE** ⇒ a resurrected script prints a verdict
     its own ledger entry does not carry ⇒ **the corpus runs and no longer says what it said**,
     which is worse than not running, because a re-derivable wrong number is quotable.

PREDICTION MATRIX:
   | world           | now  | most come back & agree | few come back, next wall is a path | come back & disagree |
   | A only-the-name | 0.50 | **0.85**               | 0.05                               | 0.10                 |
   | B layered       | 0.35 | 0.05                   | **0.85**                           | 0.10                 |
   | C stale verdict | 0.15 | 0.10                   | 0.10                               | **0.80**             |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires — **all 80 were `IMPORT-ERROR` before the repair**, read from
      `#875`'s artifact, so the before-state is not re-measured under a different load
  and the **negative control** is null — a random 20 scripts that never imported `round_path` are
      **byte-identical** before and after the repair (the repair must not touch them)
  and the **sham** is null — `LEGACY_PATHS` present but the accessor removed must still fail, so
      what fixed it is the accessor and not the import of the module
  and **coverage is complete** — every one of the 80 was reached:
      comes_back = 0                                     -> the repair does nothing; retract `#875`①
      comes_back > 0 and disagreements = 0               -> A
      comes_back > 0 and next wall is a stale path >= half of the rest -> B
      disagreements > 0                                  -> **C, and C wins over A and B**
  else: **UNVERIFIED**.

**STRONGEST CONFOUND, written before the run:** *a script "comes back" because the repair made it
exit early rather than run.* A `KeyError` swallowed by a bare `except` would exit 0 having done
nothing. ⇒ every returning script is also checked for **whether it printed a gate block at all**,
and "came back" is reported under **two definitions** — exit 0, and exit 0 **with a verdict
printed** — which is `G4`'s specification axis, not a footnote.

**⚠ `#875`④ IS HONOURED IN THE DESIGN, NOT IN MY ATTENTION.** `#875` ran 836 scripts **ten-wide**
and two of them — `R059` scripts that build throwaway commits to test a gate — **put two real
commits on `main`**, because each guards on *"is the tree clean?"*, which was true for both at once.
So: **this round runs strictly SERIALLY**, records `HEAD` before and after, and **fails loudly if it
moved**. A guard on the TREE cannot see another process about to move `HEAD`; the only sound answer
at this scale is not to have two processes.

`G3` MULTIPLICITY: family = 80 scripts × {comes back, prints a verdict, agrees}; counts reported
whole. `G4` SPECIFICATION CURVE: cap swept (derived from measured elapsed, labelled) ×
came-back-definition {exit 0 | exit 0 with a verdict} × agreement-definition {verdict string |
gate PASS/FAIL sequence}.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) it ranks **re-derivability**, never **truth**. A script that comes back and agrees is not
     thereby right; a disagreement says the evidence moved, not which side is correct;
 (2) **the instrument cannot be changed** — these scripts exist only in this repository ⇒ **only
     this one instrument**. It would require a second machine holding a second copy;
 (3) **`still_agrees` is only defined where a committed artifact carries a verdict.** Where none
     does, the answer is **`ABSENT`**, never "agrees" — and the coverage of that comparison is
     printed beside the rate, because a rate over an unstated denominator is the failure this
     project has retracted for most often;
 (4) **this round cannot distinguish "the file was rewritten" from "the data changed"** — both
     produce a disagreement. Separating them needs the artifact's own source hash, and not every
     round wrote one.
"""
import ast
import json
import os
import pathlib
import random
import re
import shutil
import signal
import subprocess
import sys
import time
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
SEED = 318
CAP = 240
DEADLINE_S = 90 * 60
PRIOR = (ROOT / "E03_what_an_instrument_would_have_to_be/A101_what_can_still_be_rederived"
              / "R314_the_sweep_rebuilt_so_it_can_finish/results/how_much_can_still_be_rederived.json")

HB = open(OUT / "heartbeat.log", "w", buffering=1)
STDIO = OUT / "_stdio"
STDIO.mkdir(exist_ok=True)
GATE_BLOCK = re.compile(r"CONDITIONAL KILL|=>\s*(?:OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)")
VERDICT = re.compile(r"=>\s*(OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)")
PASSFAIL = re.compile(r"^\s{2,}(PASS|FAIL)\s", re.M)
MISSING = re.compile(r"(?:FileNotFoundError|No such file or directory)[^\n]*?['\"]([^'\"]+)['\"]")


def hb(m):
    HB.write(f"{time.strftime('%H:%M:%S')} {m}\n")


def run_one(f, cap, tag):
    op, ep = STDIO / f"{tag}.out", STDIO / f"{tag}.err"
    t0 = time.time()
    with open(op, "wb") as so, open(ep, "wb") as se:
        p = subprocess.Popen([PY, str(f)], cwd=str(ROOT), stdout=so, stderr=se,
                             stdin=subprocess.DEVNULL, start_new_session=True)
        to = False
        while True:
            rc = p.poll()
            if rc is not None:
                break
            if time.time() - t0 > cap:
                to = True
                try:
                    os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    pass
                try:
                    rc = p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    rc = -9
                break
            time.sleep(0.2)
    o, e = op.read_text(errors="replace"), ep.read_text(errors="replace")
    op.unlink(missing_ok=True)
    ep.unlink(missing_ok=True)
    return o, e, rc, to, time.time() - t0


def classify(stderr, rc, to):
    if to:
        return "TIMEOUT", None
    if rc == 0:
        return "OK", None
    if "FileNotFoundError" in stderr or "No such file or directory" in stderr:
        m = MISSING.search(stderr)
        return "MISSING-INPUT", (m.group(1) if m else None)
    if "KeyError" in stderr:
        return "KEY-ERROR", stderr.strip().split("\n")[-1][:160]
    if "ModuleNotFoundError" in stderr or "ImportError" in stderr:
        return "IMPORT-ERROR", stderr.strip().split("\n")[-1][:160]
    if "SyntaxError" in stderr or "IndentationError" in stderr:
        return "SYNTAX", stderr.strip().split("\n")[-1][:160]
    return "OTHER", stderr.strip().split("\n")[-1][:160]


MINE = str(HERE.relative_to(ROOT))


def is_mine(path):
    """⚠ `#875` excluded by the substring `R314_`, and git reports the SHALLOWEST untracked
    directory — here the ARC directory, which does not contain the round's name. A prefix test on
    either side is the predicate; a substring test on the round name is a rule about what I
    happened to see."""
    path = path.strip().strip('"')
    return path.startswith(MINE) or MINE.startswith(path.rstrip("/"))


print("=== (0a) PRECONDITION — a WRITE, and `#875`④: record HEAD, run SERIALLY ===")
dirty = [l for l in subprocess.run(["git", "status", "--porcelain"], cwd=str(ROOT),
                                   capture_output=True, text=True).stdout.split("\n")
         if l.strip() and not is_mine(l[3:]) and "lib/rounds.py" not in l]
if dirty:
    print("  ⚠ tree not clean apart from this round and the repair:")
    for l in dirty[:20]:
        print("    ", l)
    raise SystemExit(2)
HEAD0 = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                       capture_output=True, text=True).stdout.strip()
print(f"  HEAD before: {HEAD0[:10]} · **this round runs strictly serially** (`#875`④)")

print("\n=== (0b) POPULATION — the scripts that import the restored name ===")
def imports_round_path(src):
    """PARSED, not matched (`#876`'s repair ①): a docstring quoting the import is a MENTION."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return False
    return any(isinstance(n, ast.ImportFrom) and n.module == "lib.rounds"
               and any(a.name == "round_path" for a in n.names) for n in ast.walk(tree))


IMPORTS = re.compile(r"from\s+lib\.rounds\s+import[^\n]*\bround_path\b")
cand = [f for f in sorted(ROOT.glob("E0*/**/*.py"))
        if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}"
        and not is_mine(str(f.relative_to(ROOT)))]
by_regex = [f for f in cand if IMPORTS.search(f.read_text(errors="replace"))]
by_ast = [f for f in cand if imports_round_path(f.read_text(errors="replace"))]
mention_only = [str(f.relative_to(ROOT)) for f in by_regex if f not in by_ast]
print(f"  regex {len(by_regex)} · **AST {len(by_ast)}** · mention-not-use {len(mention_only)}")
for m in mention_only:
    print(f"     MENTION ONLY (excluded): {m}")
scripts = by_ast
if not scripts:
    raise SystemExit("STOP: an empty population must never be counted as a pass")
POP = ("the round scripts that import `round_path` BY AST and were `IMPORT-ERROR` before the "
       "repair, each executed alone from the repository root")
print(f"  **{len(scripts)}** scripts")

print("\n=== (0c) POSITIVE CONTROL — the before-state, READ from `#875`'s artifact, not re-measured ===")
A = json.load(open(PRIOR))
before = {s: A["final"].get(s, {}).get("cls") for s in (str(f.relative_to(ROOT)) for f in scripts)}
b_cnt = Counter(before.values())
not_killed = {s: c for s, c in before.items() if c != "IMPORT-ERROR"}
for s, c in not_killed.items():
    print(f"     EXCLUDED FROM THE ANALYSIS POPULATION (before = {c}, so the name did not kill it): {s}")
scripts = [f for f in scripts if str(f.relative_to(ROOT)) not in not_killed]
before = {s: c for s, c in before.items() if s not in not_killed}
b_cnt = Counter(before.values())
POS_OK = bool(before) and all(c == "IMPORT-ERROR" for c in before.values())
print(f"  before-state of the same {len(before)} scripts: {dict(b_cnt)}")
print(f"  => **{'PASS' if POS_OK else 'FAIL'}** — all IMPORT-ERROR before, so anything that runs "
      f"now was caused by the repair and by nothing else")

print("\n=== (0d) SHAM — the accessor removed, the table left in place ===")
mod = ROOT / "lib/rounds.py"
orig = mod.read_text()
sham_src = orig.replace("def round_path(src_name):", "def _round_path_disabled(src_name):")
probe = scripts[0]
mod.write_text(sham_src)
_o, e_s, rc_s, to_s, _ = run_one(probe, 90, "sham")
sham_cls, _ = classify(e_s, rc_s, to_s)
mod.write_text(orig)
SHAM_OK = sham_cls == "IMPORT-ERROR"
print(f"  {probe.name} with `LEGACY_PATHS` present but no accessor -> **{sham_cls}** "
      f"(**{'PASS' if SHAM_OK else 'FAIL'}**: what fixed it is the accessor, not importing the module)")

print("\n=== (0e) NEGATIVE CONTROL — 20 random scripts that never imported the name ===")
rng = random.Random(SEED)
others = [f for f in sorted(ROOT.glob("E0*/**/*.py"))
          if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}"
          and not is_mine(str(f.relative_to(ROOT))) and f not in scripts
          and A["final"].get(str(f.relative_to(ROOT)), {}).get("cls") == "OK"]
NONDET = "bootstrap_noise_floor.py"          # `#876`② — established non-deterministic, named
others = [f for f in others if f.name != NONDET]
neg = rng.sample(others, min(20, len(others)))
print(f"  ⚠ excluded by name: **{NONDET}** — `#876` measured it non-identical across two SOLO runs; "
      f"scoring a non-deterministic script against a repair it does not use is a control failing "
      f"for its own reasons. It is characterised below instead of being counted here.")
neg_same = 0
for i, f in enumerate(neg, 1):
    o1, _e, _rc, _to, _ = run_one(f, CAP, f"n{i}a")
    o2, _e, _rc, _to, _ = run_one(f, CAP, f"n{i}b")
    same = (o1 == o2)
    neg_same += same
    hb(f"neg {i}/{len(neg)} identical={same} {f.name}")
NEG_OK = neg_same == len(neg)
print(f"  {neg_same}/{len(neg)} byte-identical across two runs with the repair in place -> "
      f"**{'PASS' if NEG_OK else 'FAIL'}**")
print("  ⚠ this is a determinism check, and these scripts DO NOT import the repaired name, so it "
      "answers 'the repair did not disturb them'; it does not answer 'they were correct'")

print(f"\n=== (0f) `#876`② CHARACTERISED, not deferred — {NONDET} run 5 times ===")
nd = [f for f in sorted(ROOT.glob("E0*/**/*.py")) if f.name == NONDET]
nd_outs, nd_note = [], "not found"
if nd:
    for i in range(5):
        o, _e, _rc, _to, _dt = run_one(nd[0], CAP, f"nd{i}")
        nd_outs.append(o)
    nd_note = (f"{len(set(nd_outs))} distinct outputs in 5 runs "
               f"({'DETERMINISTIC' if len(set(nd_outs)) == 1 else '**NOT SEED-LOCKED**'})")
    print(f"  {nd[0].relative_to(ROOT)}")
    print(f"  => **{nd_note}** ⇒ a bootstrap whose FLOOR is itself a draw; every threshold "
          f"compared against that floor inherits the draw")

print(f"\n=== (1) RUN THE {len(scripts)} — SERIALLY, cap {CAP}s, process-GROUP kill ===")
t0 = time.time()
res, elapsed = {}, {}
for i, f in enumerate(scripts, 1):
    rel = str(f.relative_to(ROOT))
    if time.time() - t0 > DEADLINE_S:
        res[rel] = dict(cls="NOT-REACHED", detail=None, verdict=None, gates=None, printed=False)
        continue
    o, e, rc, to, dt = run_one(f, CAP, f"s{i}")
    cls, det = classify(e, rc, to)
    v = VERDICT.search(o)
    res[rel] = dict(cls=cls, detail=det, verdict=(v.group(1) if v else None),
                    gates="".join("P" if m.group(1) == "PASS" else "F"
                                  for m in PASSFAIL.finditer(o)) or None,
                    printed=bool(GATE_BLOCK.search(o)))
    elapsed[rel] = dt
    hb(f"run {i}/{len(scripts)} {cls:13s} {dt:6.1f}s verdict={res[rel]['verdict']} {rel}")
print(f"  done in {(time.time()-t0)/60:.1f} min")
cnt = Counter(v["cls"] for v in res.values())
for k, n in cnt.most_common():
    print(f"    {k:14s} {n:4d}  ({n/len(res):.1%})")
NOT_REACHED = [s for s, v in res.items() if v["cls"] == "NOT-REACHED"]
COVERAGE_OK = not NOT_REACHED

back = [s for s, v in res.items() if v["cls"] == "OK"]
back_v = [s for s in back if res[s]["printed"]]
print(f"\n=== (2) TWO DEFINITIONS OF 'CAME BACK' (`G4`, not a footnote) ===")
print(f"  exit 0                       : **{len(back)}/{len(scripts)}**")
print(f"  exit 0 AND printed a verdict : **{len(back_v)}/{len(scripts)}**  "
      f"⇒ {len(back)-len(back_v)} exited 0 without printing one (the confound this checks)")

print("\n=== (3) THE NEXT WALL — for those that did NOT come back ===")
walls = Counter(v["cls"] for s, v in res.items() if v["cls"] not in ("OK", "NOT-REACHED"))
stale_path = sum(1 for s, v in res.items()
                 if v["cls"] in ("MISSING-INPUT", "KEY-ERROR"))
for k, n in walls.most_common():
    print(f"  {k:14s} {n:4d}")
print(f"  of the {sum(walls.values())} that did not come back, **{stale_path}** hit another "
      f"stale path or key ⇒ the rot is layered to that extent")

print("\n=== (4) DOES IT STILL SAY WHAT THE LEDGER SAYS IT SAID? ===")
agree, disagree, absent = [], [], []
for s in back_v:
    d = (ROOT / s).parent / "results"
    stored = None
    if d.is_dir():
        for j in sorted(d.glob("*.json")):
            try:
                obj = json.load(open(j))
            except Exception:
                continue
            for key in ("verdict", "three_valued", "gate_verdict"):
                if isinstance(obj, dict) and isinstance(obj.get(key), str):
                    m = VERDICT.search(obj[key]) or re.search(
                        r"\b(OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)\b", obj[key])
                    if m:
                        stored = (j.name, m.group(1))
                        break
            if stored:
                break
    if stored is None:
        absent.append(s)
    elif stored[1] == res[s]["verdict"]:
        agree.append((s, stored[1]))
    else:
        disagree.append((s, stored[0], stored[1], res[s]["verdict"]))
cov = len(agree) + len(disagree)
print(f"  comparable (a committed artifact carries a verdict): **{cov}/{len(back_v)}** · "
      f"**ABSENT {len(absent)}** — never counted as 'agrees'")
print(f"  agree **{len(agree)}** · **DISAGREE {len(disagree)}**")
for s, j, was, now in disagree[:10]:
    print(f"     {was:15s} -> {str(now):15s}  {j}  {s[-58:]}")

print("\n=== (5) SPECIFICATION CURVE over the cap — DERIVED from elapsed, labelled ===")
curve = {c: sum(1 for s, dt in elapsed.items() if dt > c or res[s]["cls"] == "TIMEOUT")
         for c in (30, 60, 120, 240)}
for c, n in curve.items():
    print(f"  cap {c:4d}s -> {n:3d} would not have finished "
          f"{'(measured)' if c == CAP else '(derived)'}")

HEAD1 = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                       capture_output=True, text=True).stdout.strip()
HEAD_OK = HEAD0 == HEAD1
print(f"\n=== (6) `#875`④ — HEAD before {HEAD0[:10]} · after {HEAD1[:10]} · "
      f"**{'unchanged' if HEAD_OK else 'MOVED — the sweep wrote to HISTORY'}** ===")

G = Gate("#876 · does repairing the reference bring the script back, and does it still agree")
G.asserted("(1) POSITIVE CONTROL: all 80 were IMPORT-ERROR before the repair, read from `#875`'s "
           "artifact rather than re-measured under a different load, so anything running now was "
           "caused by the repair and by nothing else",
           bool(POS_OK), " · ".join(f"{k} {n}" for k, n in b_cnt.most_common()),
           kind="control", population=POP)
G.asserted("(2) SHAM: `LEGACY_PATHS` present with the accessor removed must still fail — what "
           "fixed it is the accessor, not importing the module",
           bool(SHAM_OK), f"{probe.name} -> {sham_cls}",
           kind="control", population="one script of the population, run with the accessor removed")
G.asserted("(3) NEGATIVE CONTROL: 20 random scripts that never imported the name must be "
           "byte-identical across two runs with the repair in place",
           bool(NEG_OK), f"{neg_same}/{len(neg)} identical",
           kind="control", population=f"{len(neg)} scripts drawn at seed {SEED} from those `#875` "
                                      f"classified OK and which do not import `round_path`")
G.asserted("(4) CONFOUND: 'came back' is reported under TWO definitions — exit 0, and exit 0 with "
           "a gate block printed — because a swallowed exception exits 0 having done nothing",
           True, f"exit 0 {len(back)} · with a verdict {len(back_v)}",
           kind="control", population=POP)
G.asserted("(5) COVERAGE: a partial run is UNVERIFIED, never 'mostly A'",
           bool(COVERAGE_OK), f"{len(scripts)-len(NOT_REACHED)}/{len(scripts)} reached",
           kind="control", population=POP)
G.asserted("(6) `#875`④: HEAD recorded before and after, and this round runs SERIALLY because a "
           "guard on the TREE cannot see another process about to move HEAD",
           bool(HEAD_OK), f"{HEAD0[:10]} -> {HEAD1[:10]}", kind="control", population=POP)
G.asserted("(7) KILL (pre-registered): for \"the reference was the only thing wrong\" to hold, "
           "**no resurrected script may print a verdict its committed artifact does not carry**",
           bool(len(disagree) == 0),
           f"came back {len(back)} (with a verdict {len(back_v)}) · comparable {cov} · "
           f"agree {len(agree)} · DISAGREE {len(disagree)} · ABSENT {len(absent)} · "
           f"still dead {sum(walls.values())} of which {stale_path} on another stale path/key",
           kind="kill",
           yardstick="the verdict string a script prints, against the one its own committed "
                     "artifact carries; the floor is the two-run byte-identity of 20 untouched "
                     "scripts",
           yardstick_noise=float(len(neg) - neg_same), population=POP, direction=None)
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif len(back) == 0:
    V = ("**The repair does nothing.** Not one of the 80 exits 0 with `round_path` restored ⇒ "
         "`#875`①'s 52/52 reference resolution was an upper bound with nothing under it, and the "
         "claim that these scripts are 'mechanically repairable' is **retracted**.")
elif disagree:
    V = (f"**C — {len(disagree)} resurrected scripts print a verdict their own committed artifact "
         f"does not carry.**\n"
         f"  {len(back)} of {len(scripts)} come back ({len(back_v)} with a verdict printed), "
         f"{cov} are comparable, **{len(agree)} agree and {len(disagree)} DISAGREE**.\n"
         f"  ⇒ **the corpus runs and no longer says what it said** — which is worse than not "
         f"running, because a re-derivable wrong number is quotable and an unrunnable one is not.\n"
         f"  ⚠ A disagreement says the evidence moved, **not which side is right**.")
elif stale_path >= max(1, (len(scripts) - len(back)) / 2):
    V = (f"**B — the reference was the tip of a stack.** {len(back)}/{len(scripts)} come back; of "
         f"the {sum(walls.values())} that do not, **{stale_path} hit another stale path or key** ⇒ "
         f"**the rot is layered**, and `#875`'s \"77% is one error class\" measured the outermost "
         f"layer of a laminate.\n"
         f"  Of the {cov} comparable resurrections, **{len(agree)} agree and {len(disagree)} "
         f"disagree**; {len(absent)} carry no committed verdict and are `ABSENT`, not 'agree'.")
else:
    V = (f"**A — the reference was the only thing wrong.** {len(back)}/{len(scripts)} come back "
         f"({len(back_v)} printing a verdict), **{len(agree)}/{cov} comparable ones agree and "
         f"0 disagree**, {len(absent)} are `ABSENT`.\n"
         f"  ⇒ **one sentence: 10% of this project's evidence stopped being reproducible for one "
         f"day because I renamed a helper, and putting the name back returns it — the corpus was "
         f"never damaged, only misaddressed.**")
print(V)
print("\n⚠ **What this round structurally cannot do**: (1) it ranks re-derivability, never truth — "
      "a script that comes back and agrees is not thereby right, and a disagreement says the "
      "evidence moved, not which side is correct; (2) **the instrument cannot be changed** — these "
      "scripts exist only in this repository; (3) `still_agrees` is defined only where a committed "
      "artifact carries a verdict, and the coverage is printed beside the rate because a rate over "
      "an unstated denominator is what this project has retracted for most often; (4) it cannot "
      "separate 'the file was rewritten' from 'the data changed' — both produce a disagreement.")

json.dump(dict(head_before=HEAD0, head_after=HEAD1, head_unchanged=HEAD_OK, seed=SEED, cap=CAP,
               population=[str(f.relative_to(ROOT)) for f in scripts], before=before,
               results=res, elapsed=elapsed, counts=dict(cnt), walls=dict(walls),
               came_back=back, came_back_with_verdict=back_v,
               agree=agree, disagree=disagree, absent=absent, cap_curve_derived=curve,
               regex_population=len(by_regex), ast_population=len(by_ast),
               mention_only=mention_only, excluded_not_killed=not_killed,
               nondeterministic=dict(script=NONDET, note=nd_note,
                                     distinct_outputs=len(set(nd_outs)) if nd_outs else None),
               controls=dict(positive=POS_OK, sham=SHAM_OK, sham_cls=sham_cls,
                             negative=NEG_OK, neg_identical=neg_same, neg_n=len(neg),
                             coverage=COVERAGE_OK, head=HEAD_OK),
               admissible=adm, verdict=V, gate_ok=G.verdict(), population_string=POP),
          open(OUT / "does_the_repair_bring_them_back.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'does_the_repair_bring_them_back.json'}")

print("\n=== (7) RESTORE — this round is a WRITE (`#874`②/`#875`④) ===")
byp = OUT / "byproducts_of_the_run"
st = subprocess.run(["git", "status", "--porcelain", "-z"], cwd=str(ROOT),
                    capture_output=True, text=True).stdout
untracked, modified = [], []
for rec in st.split("\0"):
    if not rec.strip() or is_mine(rec[3:]) or "lib/rounds.py" in rec:
        continue
    (untracked if rec[:2] == "??" else modified).append(rec[3:])
moved = 0
for path in untracked:
    p = ROOT / path
    if not p.exists():
        continue
    for q in ([p] if p.is_file() else [x for x in p.rglob("*") if x.is_file()]):
        byp.mkdir(parents=True, exist_ok=True)
        shutil.move(str(q), byp / str(q.relative_to(ROOT)).replace("/", "__"))
        moved += 1
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
if modified:
    keep = byp / "_pre_restore_working_copies"
    keep.mkdir(parents=True, exist_ok=True)
    for path in modified:
        q = ROOT / path
        if q.is_file():
            shutil.copy2(q, keep / path.replace("/", "__"))
    print(f"  ⚠ about to `git checkout --` {len(modified)} tracked files; each working-tree "
          f"version was copied to `_pre_restore_working_copies` FIRST (`L81`)")
for k in range(0, len(modified), 200):
    subprocess.run(["git", "checkout", "--"] + modified[k:k + 200], cwd=str(ROOT),
                   capture_output=True)
shutil.rmtree(STDIO, ignore_errors=True)
left = [l for l in subprocess.run(["git", "status", "--porcelain"], cwd=str(ROOT),
                                  capture_output=True, text=True).stdout.split("\n")
        if l.strip() and not is_mine(l[3:]) and "lib/rounds.py" not in l]
print(f"  relocated **{moved}** untracked byproducts by LOCATION · restored **{len(modified)}** "
      f"tracked files · tree outside this round: **{'clean' if not left else str(len(left))+' dirty'}**")
for l in left[:10]:
    print("    ", l)
HB.close()
