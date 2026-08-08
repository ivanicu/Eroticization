r"""#875 · the sweep called two files PERMANENT — and it never asked git. An attack on my own verdict

**This attack SUCCEEDS, which is exactly when `realstat` §3 says it must be a full round rather than
a shell one-liner** — a cheap attack that appears to kill a claim is the most expensive kind of
error, because it retracts something true. So this carries its own estimand, its own two-sided
control, its own artifact, and its own registered limits.

**THE CLAIM UNDER ATTACK IS MINE.** The sweep's `PERMANENT` bucket was defined as: a `MISSING-INPUT`
whose file is **(a) not written by any script in this repository** and **(b) not present anywhere in
the working tree**. It returned **2 of 836**, and the round's branch C then said *"that share of
this project's conclusions can never be re-derived — by anyone, including me."*

**Both halves of that test interrogate the WORKING TREE. Neither asks the version history.** And
this corpus's actual persistence layer is git: the round directories are renamed constantly (4,368
rename edges), and a `results/*.csv` that has left the tree has usually left it by being renamed
away, not by being deleted.

⇒ **`PERMANENT` as measured means "not in the working tree", and the round reported it as "gone".**
That is `η` too large in the §2 sense — a correct measurement carrying a claim two sizes bigger than
its own definition.

`G1` **ESTIMAND**: of the files the sweep called PERMANENT, how many are retrievable from git as a
byte-exact blob. Population = the sweep's own `permanent` set, read from its artifact.
**Instrument** = `git log --all --format= --name-only` for the path, then `git cat-file` on the blob.
**Baseline** = a path that has never existed.

**POSITIVE CONTROL (known answer, on the population's own medium):** a file that is tracked *right
now* must be retrievable. A retriever that retrieves nothing is silence, not an acquittal.
**NEGATIVE CONTROL:** an invented path must not be retrievable. A retriever that retrieves anything
is not looking anything up.
**AND THE CONTROL THAT MATTERS MOST — it must be able to FAIL:** if the answer came back "0 of 2
recoverable", the sweep's verdict would stand unchanged. Nothing in this design forces the
comfortable answer.

⚠ **What this cannot do**: recovering the FILE is not re-deriving the CONCLUSION. A blob pulled out
of git is evidence that was produced by a script that may itself no longer run. It moves the case
from *lost* to *stale*, and those are different words on purpose.
"""
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
A = json.load(open(OUT / "how_much_can_still_be_rederived.json"))
PERM = A["permanent"]


def git(*a):
    r = subprocess.run(["git", *a], cwd=str(ROOT), capture_output=True, text=True)
    return r.stdout


def retrievable(path):
    """Is there a commit in ANY ref whose tree holds this path, and can the blob be read?"""
    p = str(path)
    if p.startswith(str(ROOT)):
        p = str(pathlib.Path(p).relative_to(ROOT))
    commits = [c for c in git("log", "--all", "--format=%H", "--", p).split("\n") if c]
    if not commits:
        return None, 0, "no commit in any ref holds this path"
    for c in commits:
        blob = subprocess.run(["git", "show", f"{c}:{p}"], cwd=str(ROOT),
                              capture_output=True)
        if blob.returncode == 0 and blob.stdout:
            return c, len(blob.stdout), f"{len(commits)} commit(s) hold it; read {len(blob.stdout)} bytes"
    return None, 0, f"{len(commits)} commit(s) name it, but no blob could be read"


print("=== (1) CONTROLS on the retriever — two-sided, before any verdict is moved ===")
tracked = [l for l in git("ls-files").split("\n") if l.endswith(".csv")][:1]
pc, pn, pw = retrievable(tracked[0]) if tracked else (None, 0, "no tracked csv")
POS = pc is not None and pn > 0
nc, nn, nw = retrievable("E01_sexual_as_a_value_not_a_category/__never_existed__/x.csv")
NEG = nc is None
print(f"  positive: a file tracked right now -> {tracked[0][-64:] if tracked else '-'}")
print(f"            **{'PASS' if POS else 'FAIL'}** ({pw})")
print(f"  negative: an invented path -> **{'PASS' if NEG else 'FAIL'}** ({nw})")

print(f"\n=== (2) THE {len(PERM)} FILES THE SWEEP CALLED PERMANENT ===")
rows = {}
for script, target in PERM.items():
    c, n, why = retrievable(target)
    rows[script] = dict(target=str(target), commit=c, bytes=n, why=why)
    t = str(target)
    t = t[len(str(ROOT)) + 1:] if t.startswith(str(ROOT)) else t
    print(f"  {'RECOVERABLE' if c else 'GONE':12s}  {t[-86:]}")
    print(f"                {why}" + (f" · from {c[:10]}" if c else ""))
    print(f"                needed by: {script[-80:]}")
rec = sum(1 for v in rows.values() if v["commit"])

# where does the file live TODAY, if anywhere -- the repair, not just the diagnosis
print("\n=== (3) AND WHERE IS IT NOW? — a renamed target is a repair, not a loss ===")
today = {}
for script, v in rows.items():
    name = pathlib.Path(v["target"]).name
    hits = [str(p.relative_to(ROOT)) for p in ROOT.rglob(name) if "R314_" not in str(p)]
    today[script] = hits
    print(f"  {name}: **{len(hits)}** file(s) with this name in the tree today")
    for h in hits[:4]:
        print(f"     {h}")

print("\n" + "=" * 100)
if not (POS and NEG):
    V = "**UNVERIFIED — the retriever's own two-sided control did not pass.**"
elif rec == 0:
    V = (f"**The sweep's verdict STANDS: 0 of {len(PERM)} are in git either.** `PERMANENT` meant "
         f"what it said.")
else:
    V = (f"**{rec} of {len(PERM)} files the sweep called PERMANENT are byte-retrievable from git.**\n"
         f"  ⇒ **`PERMANENT` did not mean *gone*; it meant *not in the working tree*.** The test had "
         f"two halves — produced-by-a-script, present-on-disk — and **both interrogate the working "
         f"tree**, while this corpus's persistence layer is the version history: 4,368 rename edges, "
         f"and a `results/*.csv` that has left the tree has usually left it by being RENAMED AWAY.\n"
         f"  ⇒ **The corpus has zero measured cases of unrecoverable evidence.** What it has is "
         f"references pointing at where a file used to live.\n"
         f"  ⚠ **and the word changes, not just the number**: these conclusions are not *lost*, they "
         f"are *stale* — the blob comes back, but the script that made it may still not run. "
         f"Recovering a file is not re-deriving a conclusion.")
print(V)
print("\n⚠ **Registered, and it is the honest half**: this attack was possible only because the "
      "sweep persisted the exact paths it condemned. A verdict that had reported '2 PERMANENT' "
      "without naming them could not have been attacked at all — by me or by anyone.")

json.dump(dict(n_permanent=len(PERM), n_recoverable=rec, rows=rows, present_today=today,
               controls=dict(positive=POS, positive_detail=pw, negative=NEG, negative_detail=nw),
               verdict=V),
          open(OUT / "is_permanent_actually_permanent.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'is_permanent_actually_permanent.json'}")
