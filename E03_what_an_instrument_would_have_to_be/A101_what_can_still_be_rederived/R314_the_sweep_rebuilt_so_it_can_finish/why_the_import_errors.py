r"""#875 · the dominant failure class, named — and whether it is repairable, measured not asserted

The sweep's design **excludes `IMPORT-ERROR` from the kill** as "unrelated to re-derivability".
**That exclusion is only correct if the cause is the ENVIRONMENT.** So the cause has to be
identified before the exclusion can stand — otherwise the round discards a whole class on an
assumption.
⚠ **And the number that made me look was itself a scope error worth recording**: the first 200
scripts screened came back **~40% `IMPORT-ERROR`**, and I read that as the corpus rate. It is not —
the screen walks `E01` first, and that is exactly where these scripts live. Over all 836 the class
is **79, i.e. 9.4%**. *A rate measured on a prefix is a rate about the prefix*, and the ordering was
alphabetical, which is not random with respect to anything.

**What it is** (diagnosed by executing each import statement of a failing script in isolation):
   `from lib.rounds import round_path` -> `ImportError: cannot import name 'round_path'`.
**Not a missing package. A missing name in this repository's own shared module.**

**Where it came from** (`git log -S`, two commits, both mine):
   · `5807d87` flattened `src/00..51.py` into one directory per round;
   · `4819b9b` (**2026-08-06**, one day before this round) flattened again — Ivan: 「多个 run 才能够
     算一个 round」 — and in the same commit renamed `PATHS` -> `ROUNDS` and `round_path` -> `path`.
   The registry was rebuilt with a **new key vocabulary** each time: old keys are numbered filenames
   (`'24_attack_rsa.py'`), new keys are stems (`'accumulation_curve'`), and the two sets are
   **disjoint — 52 old, 635 new, 0 shared**. So a one-line alias `round_path = path` does NOT fix
   it; the callers pass keys that no longer exist.

**THE QUESTION THIS SCRIPT ANSWERS, and it is the only one that changes a decision:**
   for each caller, is the key it passes **resolvable at all** — old key -> old path (from the
   registry as it stood before the rename, recoverable with `git show`) -> current path (via git's
   rename detection)? A resolvable key is a **mechanical repair**; an unresolvable one is a real
   loss.

`G1` **ESTIMAND**: the share of `round_path(...)` call sites whose key can be carried forward to a
file that exists today. Population = every script under E01/E02/E03 containing a `round_path(` call.
**Instrument**: the pre-rename registry read out of git, plus `git log --follow` / `--name-status`
rename detection. **Baseline**: a key that resolves to a path that exists today.

**POSITIVE CONTROL, on the population, with a known answer**: at least one key must resolve, and the
resolution must land on a file that EXISTS. A resolver that resolves nothing is silence.
**NEGATIVE CONTROL**: an invented key (`__no_such_key__.py`) must resolve to nothing. A resolver
that resolves everything is a resolver that is not looking anything up.

⚠ **What this cannot do**: it measures whether the PATH can be carried forward, never whether the
script would then RUN — the file it points at has itself been rewritten twice. Repairability of the
reference is an upper bound on repairability of the script.
"""
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
RENAME_COMMIT = "4819b9b"
CALL = re.compile(r"round_path\(\s*['\"]([^'\"]+)['\"]\s*\)")


def git(*a):
    return subprocess.run(["git", *a], cwd=str(ROOT), capture_output=True, text=True).stdout


print("=== (1) POPULATION — and it is TWO populations, which is the scope this round nearly lost ===")
scripts = [f for f in sorted(ROOT.glob("E0*/**/*.py"))
           if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}" and "R314_" not in str(f)]
IMPORTS = re.compile(r"from\s+lib\.rounds\s+import[^\n]*\bround_path\b")
importers = [str(f.relative_to(ROOT)) for f in scripts
             if IMPORTS.search(f.read_text(errors="replace"))]
callers = {}
for f in scripts:
    keys = CALL.findall(f.read_text(errors="replace"))
    if keys:
        callers[str(f.relative_to(ROOT))] = sorted(set(keys))
n_calls = sum(len(v) for v in callers.values())
print(f"  **{len(importers)}** scripts IMPORT `round_path` — every one of them dies at import, and "
      f"a separate static pass over every `from lib.* import *` in the corpus found this is the "
      f"**only** broken name in the whole shared library")
print(f"  **{len(callers)}** of them have a call site with a LITERAL key ({n_calls} references) — "
      f"the rest pass a variable or build the key, and **this script can only speak about the "
      f"literal ones.** ⚠ The repairability number below is scoped to those, not to all "
      f"{len(importers)}: a reference this instrument cannot see is not a reference it has cleared.")
if not callers:
    raise SystemExit("STOP: an empty population must never be counted as a pass")

print(f"\n=== (2) THE PRE-RENAME REGISTRY, read out of git at {RENAME_COMMIT}^ ===")
old_src = git("show", f"{RENAME_COMMIT}^:lib/rounds.py")
ns = {}
exec(compile(old_src.split("def round_path")[0], "old_rounds", "exec"), ns)
OLD = ns.get("PATHS", {})
sys.path.insert(0, str(ROOT))
from lib.rounds import ROUNDS as NEW
print(f"  old `PATHS` {len(OLD)} keys · new `ROUNDS` {len(NEW)} keys · shared "
      f"**{len(set(OLD) & set(NEW))}**")

print("\n=== (3) RENAME MAP — every path git has ever renamed, followed to today ===")
raw = git("log", "--all", "--name-status", "-M", "--diff-filter=R", "--format=%H")
ren = {}
for line in raw.split("\n"):
    if line.startswith("R"):
        parts = line.split("\t")
        if len(parts) == 3:
            ren[parts[1]] = parts[2]


def follow(p, depth=0):
    seen = set()
    while p in ren and p not in seen and depth < 40:
        seen.add(p)
        p = ren[p]
        depth += 1
    return p


print(f"  {len(ren)} rename edges recorded in history")


def resolve(key):
    old = OLD.get(key)
    if old is None:
        return None, "key not in the pre-rename registry either"
    now = follow(old)
    if (ROOT / now).exists():
        return now, "resolved" if now != old else "unchanged and still present"
    return None, f"followed to {now}, which does not exist"


print("\n=== (4) CONTROLS ===")
sample_key = next(iter(OLD))
p_ok, p_why = resolve(sample_key)
POS = p_ok is not None
n_ok, n_why = resolve("__no_such_key__.py")
NEG = n_ok is None
print(f"  positive (a key that IS in the old registry): {sample_key!r} -> {p_ok} [{p_why}] "
      f"**{'PASS' if POS else 'FAIL'}**")
print(f"  negative (an invented key): '__no_such_key__.py' -> {n_ok} [{n_why}] "
      f"**{'PASS' if NEG else 'FAIL'}**")

print("\n=== (5) RESOLUTION over the whole population ===")
res, unresolved = {}, {}
for s, keys in callers.items():
    for k in keys:
        tgt, why = resolve(k)
        (res if tgt else unresolved)[f"{s}::{k}"] = tgt or why
print(f"  resolvable **{len(res)}/{n_calls}** ({len(res)/n_calls:.1%})")
print(f"  unresolvable **{len(unresolved)}/{n_calls}**")
reasons = {}
for v in unresolved.values():
    reasons[v.split(",")[0][:60]] = reasons.get(v.split(",")[0][:60], 0) + 1
for r, n in sorted(reasons.items(), key=lambda x: -x[1])[:6]:
    print(f"     {n:4d}  {r}")
scripts_fully = [s for s, keys in callers.items()
                 if all(f"{s}::{k}" in res for k in keys)]
print(f"  scripts ALL of whose keys resolve: **{len(scripts_fully)}/{len(callers)}**")

print("\n" + "=" * 100)
if not (POS and NEG):
    V = "**UNVERIFIED — the resolver's own two-sided control did not pass.**"
elif len(res) == 0:
    V = ("**The references cannot be carried forward at all.** The pre-rename registry no longer "
         "describes anything that exists; these scripts are dead by path, not by data.")
elif len(unresolved) == 0:
    V = (f"**Every one of the {n_calls} LITERAL references, in {len(callers)} of the "
         f"{len(importers)} scripts that import the name, is mechanically repairable** — old key -> "
         f"old path -> today's path, entirely out of git. **The working tree alone does not contain "
         f"what is needed; the version history does.**\n"
         f"  ⚠ The other {len(importers)-len(callers)} importers build their key at runtime and "
         f"**this instrument cannot see them** — unseen is not cleared.")
else:
    V = (f"**{len(res)} of {n_calls} literal references ({len(res)/n_calls:.0%}) are mechanically "
         f"repairable out of git; {len(unresolved)} are not.**\n"
         f"  **{len(scripts_fully)} of {len(callers)} scripts** would have every reference "
         f"restored.\n"
         f"  ⇒ the dominant failure in this corpus is **not lost evidence — it is my own "
         f"refactors**: a shared registry rebuilt twice with a new key vocabulary, and the callers "
         f"left behind. The repair exists only because git remembers what the working tree forgot.")
print(V)
print("\n⚠ **Upper bound only**: this measures whether the REFERENCE can be carried forward, never "
      "whether the script would then run — the file it points at has itself been rewritten twice.")

json.dump(dict(importers=importers, callers=callers, n_calls=n_calls, old_keys=len(OLD), new_keys=len(NEW),
               shared_keys=len(set(OLD) & set(NEW)), rename_edges=len(ren),
               resolved=res, unresolved=unresolved,
               scripts_fully_resolvable=scripts_fully,
               controls=dict(positive=POS, positive_detail=[sample_key, p_ok, p_why],
                             negative=NEG, negative_detail=n_why),
               verdict=V),
          open(OUT / "why_the_import_errors.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'why_the_import_errors.json'}")
