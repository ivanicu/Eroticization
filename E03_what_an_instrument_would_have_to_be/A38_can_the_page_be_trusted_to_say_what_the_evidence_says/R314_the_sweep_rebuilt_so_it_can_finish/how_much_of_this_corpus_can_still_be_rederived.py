r"""#875 · E03·A101·R314 — the sweep rebuilt so it can finish, and what it says about the claims

Pays `#874`①② and answers `#872`②, which has been open for three rounds.

**THE PSYCHOLOGICAL SENTENCE THIS ROUND IS FOR.** Every claim this project makes about people —
that Americans are a gradient when asked *is it wrong* and not when asked *should he be stopped*,
that refusal is one thing inside a person and not three, that religion and education split two
different rifts — rests on a script that produced it. **A claim whose evidence can no longer be
regenerated is not false; it is no longer falsifiable**, and this project's entire method has been
to keep re-attacking its own conclusions. So the question is not academic housekeeping: it is
**how many of my sentences about people are still attackable at all.**

**WHY THIS ROUND EXISTS AND THE PREVIOUS ONE DOES NOT HAVE A RESULT.** `#874` ran this design and
**wedged at 109 minutes**, for two instrument reasons, both mine:
  · `subprocess.run(capture_output=True, timeout=...)` **does not time out when a GRANDCHILD holds
    the pipe** — the parent sat at 0.0% CPU with nothing written for 30 minutes;
  · the cost meter had been run on **38 scripts selected for another purpose** (the fast BKS-item
    ones), so the projection was built on a biased sample and was wrong by more than an order.
**Fixes, all three of them in the design rather than in my attention:**
  (a) `Popen(start_new_session=True)` + **`os.killpg` on the whole process GROUP**, and
      **stdout/stderr go to FILES, not pipes** — with no pipe, there is nothing for a grandchild to
      hold; `stdin=DEVNULL` so a script that waits for input cannot wait forever;
  (b) a **heartbeat**, one line per script, `flush=True`, so *wedged* and *slow* are distinguishable
      **from outside** without waiting;
  (c) the cost meter runs on a **RANDOM sample of the actual population**, seeded, and its projection
      is printed **before** the sweep commits.
And `#874`②: **a sweep that re-runs committed scripts is a WRITE.** So:
  (d) the round **refuses to start on a dirty tree** (exit 2), snapshots nothing because the clean
      tree IS the snapshot, and **restores tracked files with `git checkout` and relocates untracked
      byproducts BY LOCATION** — anything untracked under a round directory, whatever its extension.
      `#874` relocated by extension and missed everything that was not `.csv`.
  (e) ⚠ **caught while this round was already running, and it is the reason it was restarted:** a
      restore that reverts *anything modified* **cannot tell the sweep's byproduct from a file a
      human edited during the window**, and `git checkout --` is an `rm` of the working-tree
      version. So every discarded version is **copied into the byproducts first** (`L81`: never
      `rm`, always `mv`), and the constraint — *do not edit tracked files while this runs* — is
      **printed by the round**, not left as something I am expected to remember.
  (f) ⚠ **and the precondition earned itself on the first launch**: it refused to start, because a
      `curve.csv` timestamped **two minutes AFTER I had verified `#874`'s tree was clean** was
      sitting in an eight-year-old round directory. `#874` killed the parent and not the group, so
      its grandchildren went on writing after the cleanup. **"The tree is clean" is a measurement
      with an expiry date whenever a runaway job has just been killed.**

`G1` **ESTIMAND, named before the method**, over every round script in `E01`/`E02`/`E03`:
   (1) **`class`** — OK / MISSING-INPUT / IMPORT-ERROR / SYNTAX / TIMEOUT / OTHER;
   (2) **REPAIRABLE vs PERMANENT** — for a MISSING-INPUT, is the file produced by another script in
       this repository, or present anywhere on disk? produced ⇒ REPAIRABLE (a build order, not a
       loss); nowhere ⇒ **PERMANENT**, that conclusion can never be re-derived;
   (3) **the decision-relevant cut, which `#874` did not have** — is the script's round **cited in
       the README**? A dead script behind no live claim costs nothing. A dead script behind a README
       row is **a sentence about people with no recoverable evidence**, and that is the number.
   (4) **`verdict_before_death`** — did it print its gate block before dying? A script that dies
       after its verdict is not the same object as one that dies before.

**ARITHMETIC FIRST — what is forced, and it bounds everything below:**
   · **"cannot run" is NOT "the conclusion is wrong."** It means the conclusion is no longer
     falsifiable by re-running. **This round ranks re-derivability, never truth.**
   · the specification curve over caps is **DERIVED from the measured elapsed time**, not measured:
     a script that finished in 45 s "would have finished at cap 60" is algebra, not a test. It is
     labelled as a derivation. It is conservative in one direction only — elapsed under parallel
     load is INFLATED, so the derived curve **over**-states how many scripts are slow.

FOUR WORLDS (each with a branch):
   A **nearly everything runs** ⇒ the corpus is reproducible and `#872`'s 3 are isolated — and must
     then be NAMED as isolated in a population of 38, never generalised to the corpus.
   B **many fail, all REPAIRABLE** ⇒ the corpus is self-healing; what is missing is a build order.
   C **a real share fails PERMANENTLY** ⇒ that share of this project's conclusions can never be
     re-derived by anyone, including me. ⚠ **The unwelcome one — and this round is designed so that
     its POSITIVE outcome is the one I would rather not have** (§3 BASIN RULE: the last five rounds
     all confirmed that the instrument is repairable; this one can say it is not).
   D **⚠ META-SEPARATOR**: many scripts exit non-zero *after* printing their verdict ⇒ **"runs" is
     the wrong unit**, re-derivability is not binary, and the world-decomposition A/B/C is itself
     mis-cut.

PREDICTION MATRIX:
   | world          | now  | ~all run | fail but repairable | fail permanently | dies after verdict |
   | A reproducible | 0.25 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B self-healing | 0.35 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C permanent    | 0.30 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D wrong unit   | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the classifier's two-sided instrument control fires
  and the **population-matched positive control** fires — the **3 scripts `#872` independently
      established as dead must be classified non-OK**; this is the realstat remedy of running the
      instrument where the answer is already known, and it can fail
  and the **negative control** is null — a script containing **no file-reading call at all** must
      never be classified MISSING-INPUT
  and **coverage is complete** — every script in the population was actually reached:
      failures = 0                                            -> A
      failures > 0 and PERMANENT = 0                          -> B
      PERMANENT > 0                                           -> C
      >=1/4 of failures printed a verdict before dying        -> D
  else: **UNVERIFIED** — and a partial sweep is UNVERIFIED, never "mostly A".

**STRONGEST CONFOUND, written before the run:** *a script fails because the machine was loaded, not
because its evidence is gone.* `#872` measured exactly this — a determinism check run under load
read a timeout as non-determinism. ⇒ **every non-OK script is re-run SERIALLY, alone**, and only the
serial verdict counts. The parallel screen is a screen; the serial pass is the measurement. The
disagreement between them is reported as its own number, because it is the size of the confound.
Second confound: a failure unrelated to re-derivability (import, syntax). ⇒ only `MISSING-INPUT` is
ever eligible for PERMANENT; every other class is reported under its own name and excluded from the
kill.

`G3` MULTIPLICITY: family = every script x its class; counts reported whole, including the classes
that support nothing. `G4` SPECIFICATION CURVE: caps {30,60,120,240} (derived, labelled) x
parallel-vs-serial x REPAIRABLE-definition {produced-by-a-script | anywhere-on-disk}.
Every control row carries its **own** population string, and where it differs from the kill's the
tool says so (`#867`) — **the difference is the information, and writing them identical to silence
the warning is the failure that rule exists to catch.**

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) it ranks **re-derivability**, never **truth** — a script that will not run may well have been
     right, and one that runs may be wrong. No design here can change that;
 (2) **the instrument cannot be changed** — this corpus exists only in this repository, so there is
     **only this one instrument**. Cross-instrument replication, which this project otherwise
     demands, is structurally unavailable here. It would require a second machine holding a second
     copy of this corpus, which does not exist;
 (3) a `PERMANENT` verdict is about **THIS machine at THIS moment** — a file absent here may sit on
     another of Ivan's machines, and this round has no way to see them. It would require access to
     those machines;
 (4) **elapsed time under parallel load is not a clean runtime measurement** — the derived cap curve
     is an upper bound on slowness, not a runtime profile. A clean profile would require running the
     whole population serially, which is the cost this design exists to avoid.
"""
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
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)
PY = sys.executable
SEED = 314
CAP = 240                      # one pass; the shorter caps are DERIVED from elapsed, not re-run
WORKERS = 10                   # 43 GB available, 24 threads; each script may hold a GSS frame
DEADLINE_S = 100 * 60          # pre-registered stopping rule: partial coverage is UNVERIFIED
SAMPLE_N = 20                  # cost meter, RANDOM sample of the actual population (`#874`(1))
COST_BUDGET_S = 12 * 60        # and the cost meter has its own bound, announced when it bites

HB = open(OUT / "heartbeat.log", "w", buffering=1)          # (b) flushed, one line per script
STDIO = OUT / "_stdio"
STDIO.mkdir(exist_ok=True)

GATE_BLOCK = re.compile(r"CONDITIONAL KILL|=>\s*(?:OVERTURNED|UNVERIFIED|CONFIRMED|ALL GATES PASS)")
MISSING = re.compile(r"(?:FileNotFoundError|No such file or directory)[^\n]*?['\"]([^'\"]+)['\"]")
READS_A_FILE = re.compile(r"read_csv|read_stata|read_sas|read_spss|read_excel|read_parquet|"
                          r"read_table|read_fwf|\bopen\s*\(|\.read_text\s*\(|\.read_bytes\s*\(|"
                          r"json\.load\s*\(|np\.load")


def hb(msg):
    HB.write(f"{time.strftime('%H:%M:%S')} {msg}\n")        # buffering=1 => line-flushed


def run_one(f, cap, tag):
    """(a) process GROUP + files instead of pipes: nothing for a grandchild to hold."""
    op, ep = STDIO / f"{tag}.out", STDIO / f"{tag}.err"
    t0 = time.time()
    with open(op, "wb") as so, open(ep, "wb") as se:
        p = subprocess.Popen([PY, str(f)], cwd=str(ROOT), stdout=so, stderr=se,
                             stdin=subprocess.DEVNULL, start_new_session=True)
        timed_out = False
        while True:
            rc = p.poll()
            if rc is not None:
                break
            if time.time() - t0 > cap:
                timed_out = True
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
    o = op.read_text(errors="replace")
    e = ep.read_text(errors="replace")
    op.unlink(missing_ok=True)
    ep.unlink(missing_ok=True)
    return o, e, rc, timed_out, time.time() - t0


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


# ---------------------------------------------------------------- (0a) PRECONDITION: this is a WRITE
print("=== (0a) PRECONDITION — this round is a WRITE, so it refuses a dirty tree (`#874`(2)) ===")
dirty = subprocess.run(["git", "status", "--porcelain"], cwd=str(ROOT),
                       capture_output=True, text=True).stdout
pre = [l for l in dirty.split("\n") if l.strip() and "R314_" not in l]
if pre:
    print("  ⚠ the tree is NOT clean; re-running committed scripts would overwrite live work:")
    for l in pre[:20]:
        print("    ", l)
    raise SystemExit(2)
HEAD = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT),
                      capture_output=True, text=True).stdout.strip()
print(f"  clean tree at {HEAD[:10]} — the clean tree IS the snapshot; restore is `git checkout`")

# ---------------------------------------------------------------- (0b) POPULATION
print("\n=== (0b) POPULATION ===")
scripts = sorted(ROOT.glob("E0*/**/*.py"))
scripts = [f for f in scripts if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}"
           and "R314_" not in str(f)]
if not scripts:
    raise SystemExit("STOP: an empty population must never be counted as a pass")
POP = (f"the {len(scripts)} round scripts under E01/E02/E03, each executed from the repository "
       f"root, non-OK ones re-run serially")
print(f"  {len(scripts)} round scripts (E01/E02/E03, excluding _archive and this round)")

# ---------------------------------------------------------------- (0c) COST METER — RANDOM sample
print(f"\n=== (0c) COST METER on a RANDOM sample of {SAMPLE_N} (`#874`(1): never a subset picked "
      f"for another purpose) ===")
rng = random.Random(SEED)
sample = rng.sample(scripts, min(SAMPLE_N, len(scripts)))
hb(f"cost-meter start n={len(sample)}")
times, t_cost = [], time.time()
for i, f in enumerate(sample, 1):
    if time.time() - t_cost > COST_BUDGET_S:
        print(f"  ⚠ cost meter hit its own {COST_BUDGET_S//60}-min budget after {len(times)} of "
              f"{len(sample)} — **{len(sample)-len(times)} dropped, said out loud rather than "
              f"silently truncated**; the projection below rests on {len(times)} draws")
        break
    _o, _e, _rc, _to, dt = run_one(f, CAP, f"cost{i}")
    times.append(dt)
    hb(f"cost {i}/{len(sample)} {dt:6.1f}s {f.relative_to(ROOT)}")
times_sorted = sorted(times)
med = times_sorted[len(times_sorted) // 2]
p90 = times_sorted[int(0.9 * (len(times_sorted) - 1))]
mean = sum(times) / len(times)
proj_serial = mean * len(scripts)
proj_par = proj_serial / WORKERS
print(f"  median {med:.1f}s · mean {mean:.1f}s · p90 {p90:.1f}s · max {max(times):.1f}s")
print(f"  projected SERIAL {proj_serial/60:.0f} min · projected at {WORKERS} workers "
      f"**{proj_par/60:.0f} min** · pre-registered deadline {DEADLINE_S/60:.0f} min")
_verdict_cost = ("proceeding" if proj_par < DEADLINE_S else
                 "⚠ projection EXCEEDS the deadline — the sweep will run to the deadline and "
                 "report partial coverage as UNVERIFIED, never as A")
print(f"  ⇒ {_verdict_cost}")
hb(f"cost-meter done mean={mean:.1f}s projected_parallel={proj_par/60:.0f}min")

# ---------------------------------------------------------------- (0d) INSTRUMENT CONTROL (probe)
print("\n=== (0d) INSTRUMENT CONTROL — two-sided, on a synthetic probe ===")
probe = ROOT / "__r314_probe.py"
probe.write_text("import pandas as pd\npd.read_csv('data/raw/__no_such_file_r314__.csv')\n")
_o, e1, rc1, to1, _ = run_one(probe, 60, "probe_pos")
pos_cls, pos_path = classify(e1, rc1, to1)
probe.write_text("print('hello')\n")
_o, e2, rc2, to2, _ = run_one(probe, 60, "probe_neg")
neg_cls, _ = classify(e2, rc2, to2)
probe.unlink()
PROBE_POP = ("a synthetic 2-line probe — ⚠ NOT the round population; this row validates the "
             "CLASSIFIER, not the corpus")
print(f"  positive: pointed at a non-existent file -> **{pos_cls}** ({pos_path})")
print(f"  negative: reads nothing                  -> **{neg_cls}**")
CTRL_OK = (pos_cls == "MISSING-INPUT") and (neg_cls == "OK")
print(f"  => instrument control **{'PASS' if CTRL_OK else 'FAIL'}**")

# ---------------------------------------------------------------- (1) PARALLEL SCREEN
print(f"\n=== (1) PARALLEL SCREEN — {WORKERS} workers, cap {CAP}s, process-GROUP kill ===")
t0 = time.time()
done = [0]
res = {}


def work(idx_f):
    i, f = idx_f
    if time.time() - t0 > DEADLINE_S:
        return str(f.relative_to(ROOT)), ("NOT-REACHED", None, False, 0.0)
    o, e, rc, to, dt = run_one(f, CAP, f"s{i}")
    cls, det = classify(e, rc, to)
    done[0] += 1
    hb(f"screen {done[0]}/{len(scripts)} {cls:13s} {dt:6.1f}s {f.relative_to(ROOT)}")
    return str(f.relative_to(ROOT)), (cls, det, bool(GATE_BLOCK.search(o)), dt)


with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    for k, v in ex.map(work, list(enumerate(scripts))):
        res[k] = v
screen_s = time.time() - t0
print(f"  screen done in {screen_s/60:.1f} min")
scr_cnt = Counter(v[0] for v in res.values())
for k, n in scr_cnt.most_common():
    print(f"    {k:14s} {n:4d}")
NOT_REACHED = [s for s, v in res.items() if v[0] == "NOT-REACHED"]
COVERAGE_OK = not NOT_REACHED
print(f"  coverage: {len(scripts)-len(NOT_REACHED)}/{len(scripts)} reached "
      f"**{'complete' if COVERAGE_OK else 'PARTIAL — the kill is UNVERIFIED'}**")

# ---------------------------------------------------------------- (2) SERIAL CONFIRMATION
print(f"\n=== (2) SERIAL CONFIRMATION — the load confound `#872` measured; only this pass counts ===")
suspect = [s for s, v in res.items() if v[0] not in ("OK", "NOT-REACHED")]
print(f"  {len(suspect)} non-OK scripts, re-run ALONE at cap {CAP}s")
FINAL, flipped = {}, []
hb(f"serial start n={len(suspect)}")
for i, s in enumerate(suspect, 1):
    o, e, rc, to, dt = run_one(ROOT / s, CAP, f"c{i}")
    cls, det = classify(e, rc, to)
    hb(f"serial {i}/{len(suspect)} {res[s][0]}->{cls} {dt:6.1f}s {s}")
    if cls != res[s][0]:
        flipped.append((s, res[s][0], cls))
    FINAL[s] = (cls, det, bool(GATE_BLOCK.search(o)), dt)
for s, v in res.items():
    FINAL.setdefault(s, v)
print(f"  **{len(flipped)}** changed class between the loaded screen and the solo re-run "
      f"— that number IS the size of the confound")
for s, a, b in flipped[:10]:
    print(f"     {a:13s} -> {b:13s}  {s[:78]}")

cnt = Counter(v[0] for v in FINAL.values())
print("\n=== (3) CLASSIFICATION — the whole grid, including the classes that support nothing ===")
for k, n in cnt.most_common():
    print(f"  {k:14s} {n:4d}  ({n/len(FINAL):.1%})")

# ---------------------------------------------------------------- (4) NEGATIVE CONTROL, on the pop
no_read = [str(f.relative_to(ROOT)) for f in scripts
           if not READS_A_FILE.search(f.read_text(errors="replace"))]
nr_missing = [s for s in no_read if FINAL.get(s, ("",))[0] == "MISSING-INPUT"]
NEG_POP = (f"the {len(no_read)} scripts of the population containing no file-reading call — a "
           f"NAMED SUBSET of the kill's population, not the same string, and that is the honest "
           f"statement (`#867`)")
print(f"\n=== (4) NEGATIVE CONTROL — {len(no_read)} scripts open no file; MISSING-INPUT among "
      f"them: **{len(nr_missing)}** (must be 0) ===")
NEG_OK = (len(no_read) > 0) and (len(nr_missing) == 0)

# ---------------------------------------------------------------- (5) POSITIVE CONTROL, on the pop
KNOWN_DEAD = [
    "E01_sexual_as_a_value_not_a_category/A05_what_the_second_dimension_is/R005_shame_and_the_position/rho_panel_controlling_age.py",
    "E01_sexual_as_a_value_not_a_category/A06_can_component_three_earn_a_name/R006_尺子装对之后_第三个维度的正名是_breadth_而_189_的三个/is_the_advantage_only_in_large_samples.py",
    "E01_sexual_as_a_value_not_a_category/A06_can_component_three_earn_a_name/R006_尺子装对之后_第三个维度的正名是_breadth_而_189_的三个/the_als_loading_as_a_person_variable.py",
]
kd_present = [s for s in KNOWN_DEAD if s in FINAL]
kd_cls = {s: FINAL[s][0] for s in kd_present}
POS_OK = bool(kd_present) and all(c != "OK" for c in kd_cls.values())
KD_POP = (f"the {len(kd_present)} scripts `#872` independently established as DOES-NOT-RUN — a "
          f"3-element subset of the kill's population whose answer was known BEFORE this "
          f"instrument existed")
print(f"\n=== (5) POSITIVE CONTROL on the POPULATION — the 3 scripts `#872` already knew were dead "
      f"(realstat: run the instrument where the answer is known) ===")
for s, c in kd_cls.items():
    print(f"  {c:14s}  {s[-72:]}")
print(f"  => **{'PASS' if POS_OK else 'FAIL'}** — this control CAN fail: if the sweep called them "
      f"OK, the instrument would be blind to exactly the thing it was built for")

# ---------------------------------------------------------------- (6) REPAIRABLE vs PERMANENT
produced = set()
for f in scripts:
    t = f.read_text(errors="replace")
    for m in re.finditer(r"""to_csv\(\s*[^)]*?['"]([^'"]+\.csv)['"]""", t):
        produced.add(pathlib.Path(m.group(1)).name)
    for m in re.finditer(r"""OUT\s*/\s*['"]([^'"]+)['"]""", t):
        produced.add(pathlib.Path(m.group(1)).name)
miss = {s: d for s, (c, d, _v, _t) in FINAL.items() if c == "MISSING-INPUT" and d}
rep, perm, rep_strict = {}, {}, {}
for s, d in miss.items():
    name = pathlib.Path(d).name
    elsewhere = list(ROOT.rglob(name))
    if name in produced:
        rep[s] = (d, "produced by another script in this repository")
        rep_strict[s] = rep[s]
    elif elsewhere:
        rep[s] = (d, f"not produced by any script, but present at "
                     f"{elsewhere[0].relative_to(ROOT)}")
    else:
        perm[s] = d
print(f"\n=== (6) REPAIRABLE vs PERMANENT (only MISSING-INPUT is eligible) ===")
print(f"  MISSING-INPUT with an identified path: **{len(miss)}**")
print(f"  **REPAIRABLE {len(rep)}** — of which {len(rep_strict)} under the STRICT definition "
      f"(produced by a script) and {len(rep)-len(rep_strict)} only under the LOOSE one (exists "
      f"somewhere on disk)  <- G4 specification axis")
print(f"  **PERMANENT  {len(perm)}** — nothing produces it and it is nowhere on disk")
for s, d in list(perm.items())[:10]:
    print(f"     {s[-70:]:70s} <- {d[:56]}")

# ------------------------------------------------- (7) THE DECISION-RELEVANT CUT: is it a live claim?
led = (ROOT / "RETRACTIONS.md").read_text(errors="replace")
readme = (ROOT / "README.md").read_text(errors="replace") + (ROOT / "README_zh.md").read_text(errors="replace")
ear = {}
for m in re.finditer(r'^## Entr(?:y|ies) (\d+) · `([^`]+)`', led, re.M):
    for r_id in re.findall(r'R(\d+)', m.group(2)):
        ear.setdefault(f"R{int(r_id):03d}", set()).add(int(m.group(1)))
cited = {n for n in re.findall(r'Entry (\d+)', readme)}
cited = {int(n) for n in cited}


def round_of(script_path):
    for part in script_path.split("/"):
        m = re.match(r'^R(\d+)_', part)
        if m:
            return f"R{int(m.group(1)):03d}"
    return None


def is_live(script_path):
    r_id = round_of(script_path)
    if r_id is None:
        return None
    return bool(ear.get(r_id, set()) & cited)


# positive control on THIS instrument too: a round whose entry is cited must come back live
probe_entry = 528
probe_rounds = [r for r, es in ear.items() if probe_entry in es]
CHAIN_OK = bool(probe_rounds) and (probe_entry in cited)
print(f"\n=== (7) THE CUT THAT MATTERS — is the script's round behind a claim the README still "
      f"makes? ===")
print(f"  chain control: Entry {probe_entry} maps to {probe_rounds} and is cited in the README: "
      f"**{'PASS' if CHAIN_OK else 'FAIL'}** (the chain script->round->entry->README is an "
      f"instrument and gets its own known answer)")
failed_all = [s for s, v in FINAL.items() if v[0] not in ("OK",)]
live_dead = [s for s in failed_all if is_live(s)]
live_perm = [s for s in perm if is_live(s)]
print(f"  scripts that do not run: {len(failed_all)} · of them behind a LIVE README claim: "
      f"**{len(live_dead)}**")
print(f"  PERMANENT: {len(perm)} · of them behind a LIVE README claim: **{len(live_perm)}**")
for s in live_perm[:10]:
    print(f"     {s[-90:]}")

# ---------------------------------------------------------------- (8) DERIVED CAP CURVE (labelled)
elapsed = {s: v[3] for s, v in FINAL.items() if v[0] != "NOT-REACHED"}
curve = {}
for c in (30, 60, 120, 240):
    curve[c] = sum(1 for s, dt in elapsed.items() if dt > c or FINAL[s][0] == "TIMEOUT")
print(f"\n=== (8) SPECIFICATION CURVE over the cap — ⚠ **DERIVED from measured elapsed, not "
      f"re-measured** (the arithmetic trap, labelled) ===")
for c, n in curve.items():
    print(f"  cap {c:4d}s -> {n:4d} scripts would not have finished "
          f"{'(measured: this is the cap that ran)' if c == CAP else '(derived)'}")
print(f"  ⚠ elapsed was recorded under {WORKERS}-way load ⇒ inflated ⇒ the derived curve "
      f"OVER-states slowness. It bounds one direction only.")

DIED_AFTER = [s for s, v in FINAL.items() if v[0] not in ("OK", "NOT-REACHED") and v[2]]
print(f"\n=== (9) world D — scripts that printed a verdict and THEN failed: **{len(DIED_AFTER)}** ===")

# ---------------------------------------------------------------- (10) THE GATE
C = len([v for v in FINAL.values() if v[0] != "NOT-REACHED"])
fails = C - cnt.get("OK", 0)
G = Gate("#875 · how much of this corpus can still be re-derived, on a sweep that can finish")
G.asserted("(1) INSTRUMENT CONTROL: a script pointed at a non-existent file must classify "
           "MISSING-INPUT and a script that reads nothing must classify OK",
           bool(CTRL_OK), f"positive -> {pos_cls} · negative -> {neg_cls}",
           kind="control", population=PROBE_POP)
G.asserted("(2) POSITIVE CONTROL ON THE POPULATION: the 3 scripts `#872` independently established "
           "as dead must be classified non-OK — the instrument run where the answer is already "
           "known, and it can fail",
           bool(POS_OK), " · ".join(f"{pathlib.Path(s).name}->{c}" for s, c in kd_cls.items()),
           kind="control", population=KD_POP)
G.asserted("(3) NEGATIVE CONTROL: a script containing no file-reading call must never be "
           "classified MISSING-INPUT",
           bool(NEG_OK), f"{len(no_read)} such scripts · MISSING-INPUT among them {len(nr_missing)}",
           kind="control", population=NEG_POP)
G.asserted("(4) LOAD CONFOUND (`#872`): every non-OK script was re-run SERIALLY and only the solo "
           "verdict counts; the screen-vs-solo disagreement is reported as the size of the confound",
           True, f"{len(suspect)} re-run alone · {len(flipped)} changed class",
           kind="control", population=POP)
G.asserted("(5) COVERAGE: a partial sweep is UNVERIFIED, never 'mostly A' — the pre-registered "
           "deadline stops the run and says what it did not reach",
           bool(COVERAGE_OK), f"{C}/{len(scripts)} reached · not reached {len(NOT_REACHED)} · "
                              f"screen {screen_s/60:.1f} min vs deadline {DEADLINE_S/60:.0f} min",
           kind="control", population=POP)
G.asserted("(6) CHAIN CONTROL: script -> round -> ledger entry -> README citation is an instrument; "
           "an entry known to be cited must come back cited",
           bool(CHAIN_OK), f"Entry {probe_entry} -> rounds {probe_rounds} -> cited {probe_entry in cited}",
           kind="control", population="the ledger's 875 entries and the two README pages")
G.asserted("(7) ARITHMETIC: 'cannot run' is NOT 'the conclusion is wrong'; and the cap curve is "
           "DERIVED from elapsed, not re-measured — labelled as a derivation, conservative in one "
           "direction only",
           True, f"caps derived {list(curve)} · elapsed under {WORKERS}-way load ⇒ inflated",
           kind="control", population=POP)
G.asserted("(8) KILL (pre-registered): for \"every conclusion in this corpus can still be "
           "re-derived\" to hold, **no script may fail on a missing input that nothing produces "
           "and that is nowhere on disk**",
           bool(len(perm) == 0),
           f"PERMANENT {len(perm)} (live-claim {len(live_perm)}) · REPAIRABLE {len(rep)} "
           f"(strict {len(rep_strict)}) · failures {fails}/{C} · died-after-verdict {len(DIED_AFTER)}",
           kind="kill",
           yardstick="execution class per script, measured SOLO, with only MISSING-INPUT eligible "
                     "for PERMANENT; the floor is the classifier's two-sided control plus the 3 "
                     "known-dead scripts",
           yardstick_noise=float(len(flipped)), population=POP,
           direction=None)
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif len(DIED_AFTER) >= max(1, fails / 4) and fails:
    VERD = (f"**D — {len(DIED_AFTER)} of {fails} failures printed a verdict and THEN died ⇒ 'runs' "
            f"is the wrong unit; the A/B/C decomposition is itself mis-cut.**")
elif fails == 0:
    VERD = (f"**A — every one of the {C} scripts still runs.** `#872`'s 3 must then be named as "
            f"isolated **in a population of 38**, never generalised to {C}.")
elif not perm:
    VERD = (f"**B — {fails}/{C} scripts fail and NONE permanently.** Every missing input is either "
            f"produced by another script here ({len(rep_strict)}) or present elsewhere on disk "
            f"({len(rep)-len(rep_strict)}).\n"
            f"  The corpus is self-healing: what is missing is a build order, not evidence.\n"
            f"  ⇒ **One sentence about people: every claim this project makes about how Americans "
            f"condemn is still attackable — the evidence behind it can be regenerated.**")
else:
    VERD = (f"**C — {len(perm)} of {C} scripts fail on an input that nothing in this repository "
            f"produces and that is nowhere on disk.**\n"
            f"  Of {C}: **{cnt.get('OK',0)} run** · {fails} fail · of the failures **{len(rep)} "
            f"REPAIRABLE** and **{len(perm)} PERMANENT**.\n"
            f"  **{len(live_perm)} of the PERMANENT ones sit behind a claim the README still "
            f"makes.** That is the number that costs something.\n"
            f"  ⇒ **One sentence about people: {len(live_perm)} of this project's live sentences "
            f"about how Americans condemn can no longer be attacked, because the evidence that "
            f"produced them can no longer be produced.**\n"
            f"  ⚠ **'Cannot re-derive' is not 'wrong'** — those conclusions may well be correct. "
            f"They are simply no longer falsifiable by re-running, which is the one thing this "
            f"project has relied on from round one.")
print(VERD)
print(f"\n**What this round structurally cannot do**: (1) it ranks re-derivability, never truth; "
      f"(2) **the instrument cannot be changed** — this corpus exists only in this repository, so "
      f"cross-instrument replication is structurally unavailable and would require a second machine "
      f"holding a second copy; (3) a PERMANENT verdict is about **this machine at this moment**; "
      f"(4) elapsed under parallel load is not a runtime profile — the cap curve bounds slowness "
      f"from one side only.")

json.dump(dict(head=HEAD, seed=SEED, cap=CAP, workers=WORKERS, deadline_s=DEADLINE_S,
               population=[str(f.relative_to(ROOT)) for f in scripts],
               cost_meter=dict(n_drawn=len(sample), n_timed=len(times), median=med, mean=mean,
                               p90=p90, max=max(times), budget_s=COST_BUDGET_S,
                               projected_parallel_min=proj_par / 60,
                               sample=[str(f.relative_to(ROOT)) for f in sample]),
               screen=dict(seconds=screen_s, counts=dict(scr_cnt), not_reached=NOT_REACHED),
               serial_flips=flipped,
               final={k: dict(cls=v[0], detail=v[1], verdict_printed=v[2], elapsed=v[3])
                      for k, v in FINAL.items()},
               counts=dict(cnt), repairable=rep, repairable_strict=list(rep_strict),
               permanent=perm, live_dead=live_dead, live_permanent=live_perm,
               died_after_verdict=DIED_AFTER, cap_curve_derived=curve,
               controls=dict(instrument=dict(positive=pos_cls, negative=neg_cls, ok=CTRL_OK),
                             known_dead=kd_cls, known_dead_ok=POS_OK,
                             no_read_n=len(no_read), no_read_missing=nr_missing, neg_ok=NEG_OK,
                             coverage_ok=COVERAGE_OK, chain_ok=CHAIN_OK),
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), population_string=POP),
          open(OUT / "how_much_can_still_be_rederived.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'how_much_can_still_be_rederived.json'}")

# ------------------------------------------------- (11) RESTORE — `#874`(2): the sweep is a WRITE
print("\n=== (11) RESTORE — a sweep that re-runs committed scripts is a WRITE, and the restore is "
      "part of the design, not something I notice afterwards (`#874`(2)) ===")
byp = OUT / "byproducts_of_the_sweep"
st = subprocess.run(["git", "status", "--porcelain", "-z"], cwd=str(ROOT),
                    capture_output=True, text=True).stdout
moved, restored = 0, 0
untracked, modified = [], []
for rec in st.split("\0"):
    if not rec.strip():
        continue
    code, path = rec[:2], rec[3:]
    if "R314_" in path:
        continue
    if code == "??":
        untracked.append(path)
    else:
        modified.append(path)
# relocate untracked BY LOCATION, not by extension (`#874` moved only .csv and missed the rest)
for path in untracked:
    p = ROOT / path
    if not p.exists():
        continue
    targets = [p] if p.is_file() else [q for q in p.rglob("*") if q.is_file()]
    for q in targets:
        byp.mkdir(parents=True, exist_ok=True)
        dest = byp / str(q.relative_to(ROOT)).replace("/", "__")
        shutil.move(str(q), dest)
        moved += 1
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
# ⚠ `git checkout --` is an `rm` of the working-tree version, and this restore cannot tell the
# sweep's byproduct from a file a HUMAN edited during the window. L81 says never rm, always mv:
# so every version about to be discarded is COPIED into the byproducts first, and the constraint
# ("do not edit tracked files while this runs") is printed rather than assumed.
if modified:
    keep = byp / "_pre_restore_working_copies"
    keep.mkdir(parents=True, exist_ok=True)
    for path in modified:
        p = ROOT / path
        if p.is_file():
            shutil.copy2(p, keep / path.replace("/", "__"))
    print(f"  ⚠ about to `git checkout --` {len(modified)} tracked files; the working-tree "
          f"version of each was copied to {keep.relative_to(ROOT)} FIRST (L81), because this "
          f"restore cannot distinguish a sweep byproduct from a human edit made during the run")
for k in range(0, len(modified), 200):          # chunked: a long path list can exceed ARG_MAX
    subprocess.run(["git", "checkout", "--"] + modified[k:k + 200], cwd=str(ROOT),
                   capture_output=True)
    restored += len(modified[k:k + 200])
shutil.rmtree(STDIO, ignore_errors=True)
after = subprocess.run(["git", "status", "--porcelain"], cwd=str(ROOT),
                       capture_output=True, text=True).stdout
left = [l for l in after.split("\n") if l.strip() and "R314_" not in l]
print(f"  relocated **{moved}** untracked byproducts (by LOCATION, any extension) · restored "
      f"**{restored}** tracked files with `git checkout`")
print(f"  tree outside this round after restore: **{'clean' if not left else str(len(left)) + ' still dirty'}**")
for l in left[:10]:
    print("    ", l)
HB.close()
