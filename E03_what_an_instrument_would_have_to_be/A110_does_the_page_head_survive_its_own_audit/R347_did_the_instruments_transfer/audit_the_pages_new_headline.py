r"""#909 · E03·A110·R347 — the page's new headline, audited against the artifacts that would refute it

**COGNITIVE UPDATE CARD**
```
Core Gap        `#908` rewrote the page head to lead with SIX INSTRUMENTS as the project's residue,
                and in the same entry registered `#908`②: **that claim is the page's new headline
                and it has NO CONTROL.** The failable version was written there: *can each be
                pointed at an object it was NOT built for, and does it still carry its positive
                control there?* `#908` asserted **"two have been, four have not"** — from memory,
                on the same day, about my own work. **That number has never been measured.**
Why Now         It is the newest load-bearing sentence on the deliverable and the only one nobody
                has attacked. `§0.2` says the residue is the product; a residue claim that is wrong
                in either direction is worse than none.
Live Worlds     CORRECT     exactly two instruments have a controlled off-domain application.
                ⚠ UNDERSTATED  more than two do — the page understates its own residue, which is
                            the error `#899` pre-registered as *"possibly over-modest"* and which
                            `§0.2` says has its optimum at doing nothing.
                OVERSTATED  fewer than two do — "six instruments" is inflated and the head is the
                            kind of sentence this project exists to prevent.
                BLIND       ⚠ META — the artifacts do not record enough to tell, in which case the
                            page's headline is unauditable and must say so.
Discriminating  For each instrument, find every round whose SCRIPT uses it, split into the round
Act             that BUILT it and later rounds that applied it to a different object, and read that
                round's PERSISTED ARTIFACT for whether its positive control fired there. **The
                evidence is the artifact, not my memory of the round.**
Prediction      CORRECT -> exactly 2 · UNDERSTATED -> >2 · OVERSTATED -> <2 · BLIND -> artifacts silent
Matrix
Confound        ⚠ WRITTEN BEFORE THE RUN: **"uses it" is a text match on a function name, and a
                text scan counts MENTIONS as USES** — this project has a memory file about exactly
                that. So the scan parses the script with `ast` and counts a DEFINITION or a CALL,
                never a string or a comment. And **a shared function name is not a shared
                instrument**, so each signature is checked against the round that built it.
Controls        positive: `#906` (R344) MUST come back as a controlled transfer of the case-indexing
                interaction — it is the one I am certain of · negative: an invented signature must
                return zero rounds · empty population -> exit 2, never 0
Stopping Rule   One round. Whatever the count, the page is corrected to it in the same commit,
                because a headline that survives its own audit unchanged is the only kind worth
                keeping. Budget: one round.
Cost            ~900 scripts parsed, 217 artifacts read. CPU seconds.
Priority        It is the page's newest claim and the only one with no control.
Expected        If UNDERSTATED: I have been shrinking my own residue while writing rounds that warn
Transform       against exactly that, which is a finding about the loop and not about GSS.
```

⚠⚠ **`#901`①'s REMEDY, SEVENTH USE.** Outcome space = the integer count of instruments with ≥1
controlled off-domain application, `0..6`, **partitioned before the run**: `<2`→OVERSTATED ·
`==2`→CORRECT · `>2`→UNDERSTATED · **and "artifacts silent for ≥half the transfers"→BLIND, which
overrides all three.** No value is unlisted.

`G1` **ESTIMAND**: **per instrument, the number of distinct rounds that (a) use it by AST definition
or call, (b) are LATER than the round that built it, and (c) persist an artifact recording a positive
control.** **Population** every `*.py` under `E0*/A*/R*` (868) and every `results/*.json` (217).
**Instrument** this repository's own source and artifacts — ⚠ **`no second instrument`,
`only this one instrument`**, because the object under test IS the project's own claim about itself.
**Baseline** `#908`②'s asserted "two". **Regime** whole repository.

⚠ **"SHOULD THIS ZERO BE ZERO?" — the answer is a COUNT, not an effect**, so the failable object is
the **detector**: it gets a positive control (`#906` must be recovered), a negative control (an
invented signature returns zero), and an empty-population exit. **Saying which control shape applies
rather than forcing a null distribution is the honest reading of the rule** (`#907` did the same).

**PRE-REGISTERED KILL — a conditional:**
```
if `#906` is recovered as a controlled transfer AND an invented signature returns zero rounds:
       count the instruments with >=1 controlled off-domain application and read the partition above
else:
       UNVERIFIED — the detector cannot see, and the page's headline stays unaudited and says so
```
`G3`: every instrument reported with its builder round and every transfer, including the ones with
no artifact.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **"a different object" is judged by ME** — the scan finds reuse, and whether the object differs is
   a reading; `door ③` says self-review is void and this round is **`[unchallenged]`**;
② **`ast` sees literal calls** — a call built at runtime is **UNSEEN, not absent** (`#875` measured
   that as 30 of 80 importers);
③ **two of the six are GATES, not statistics** — their "transfer" is that they run on every commit
   over the whole repository, which is a different kind of evidence and is reported as such, not
   averaged in;
④ **an artifact recording `positive_ok` does not prove the control was WELL-SPECIFIED** — this
   session alone rebuilt five that were not;
⑤ no second coder, no second release.
"""
import ast
import json
import pathlib
import sys

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

# instrument -> (signature names, the entry that BUILT it, what it was built for)
INSTR = {
    "link-free comonotonicity test": (["min_abs_spearman"], 902, "the four GSS sexual-norm series"),
    "one-factor-proof case-indexing interaction": (["interaction"], 892, "GSS abortion norm×legality"),
    "rank-1 binomial resampling null": (["rank1_fit", "rank1_probs", "rank1_probs_link"], 900,
                                        "the four sexual series' probit matrix"),
    "comonotone binomial resampling null": (["comonotone_probs"], 902, "the four sexual series"),
}
GATES = {"tools/supersede_gate.py": 894, "tools/registry_keys_gate.py": 877}
FAKE = "flurbish_null"          # the negative control's invented signature


def used_names(path):
    """DEFINITIONS and CALLS only — never a string, never a comment.

    ⚠ `feedback_a_text_scan_counts_mentions_as_uses`: the better the ledger, the noisier the grep.
    This project has a memory file about counting a report that QUOTES a defect as the defect."""
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except Exception:
        return set()
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.add(n.name)
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Name):
            out.add(n.func.id)
    return out


print("=== (0) POPULATION — parsed, not globbed, and counted before anything is claimed ===")
scripts = sorted(ROOT.glob("E0*/A*/R*/*.py"))
arts = sorted(ROOT.glob("E0*/A*/R*/results/*.json"))
print(f"  round scripts: {len(scripts)} · persisted artifacts: {len(arts)}")
if not scripts or not arts:
    print("EMPTY POPULATION — exit 2, never 0.")
    sys.exit(2)

ART = {}
for a in arts:
    try:
        d = json.loads(a.read_text())
    except Exception:
        continue
    e = d.get("entry")
    if isinstance(e, int):
        ART.setdefault(e, {}).update({k: d.get(k) for k in
                                      ("positive_ok", "gates", "gate_verdict", "verdict", "round")})
print(f"  artifacts carrying an entry number: {len(ART)}  "
      f"(the rest predate the convention and are UNSEEN, not absent)")


def entry_of(script):
    """The entry number a round script declares in its own docstring header."""
    txt = script.read_text(errors="replace")[:400]
    import re
    m = re.search(r"#(\d{3,4})\s*·", txt)
    return int(m.group(1)) if m else None


USES = {}
for s in scripts:
    e = entry_of(s)
    if e is None:
        continue
    USES.setdefault(e, set()).update(used_names(s))
print(f"  scripts declaring an entry number: {sum(1 for e in USES)}")


def transfers(sigs, built):
    """Rounds LATER than the builder that define or call one of the signatures."""
    return sorted(e for e, names in USES.items() if e > built and (set(sigs) & names))


print("\n=== (1) POSITIVE CONTROL — `#906` must come back as a transfer of the case-indexing tool ===")
ci = transfers(INSTR["one-factor-proof case-indexing interaction"][0], 892)
pos_ok = 906 in ci and bool(ART.get(906, {}).get("positive_ok"))
print(f"  case-indexing transfers found: {ci}")
print(f"  `#906` present: {906 in ci} · its artifact records positive_ok="
      f"{ART.get(906, {}).get('positive_ok')} ⇒ control fires: **{pos_ok}**")
neg = transfers([FAKE], 0)
neg_ok = not neg
print(f"  negative control — invented signature `{FAKE}` returns {len(neg)} rounds ⇒ **{neg_ok}**")

print("\n=== (2) EVERY INSTRUMENT, ITS BUILDER, AND EVERY LATER ROUND THAT USES IT ===")
rows, silent = [], 0
for name, (sigs, built, forwhat) in INSTR.items():
    tr = transfers(sigs, built)
    controlled = [e for e in tr if ART.get(e, {}).get("positive_ok") is True
                  or (ART.get(e, {}).get("gate_verdict") or "").startswith("ALL GATES PASS")]
    nosig = [e for e in tr if e not in ART]
    silent += len(nosig)
    rows.append((name, built, forwhat, tr, controlled, nosig))
    print(f"\n  {name}")
    print(f"     built in `#{built}` for {forwhat}")
    print(f"     later rounds using it : {tr}")
    print(f"     ...with a control recorded in their artifact: **{controlled}**")
    if nosig:
        print(f"     ...with NO artifact entry (UNSEEN, not absent): {nosig}")

print("\n=== (3) THE TWO GATES — a different kind of evidence, reported and not averaged in ===")
hook = pathlib.Path.home() / ".claude/neural-commit/hooks/pre-commit"
hooktxt = hook.read_text(errors="replace") if hook.exists() else ""
for g, built in GATES.items():
    wired = g.split("/")[-1] in hooktxt
    print(f"  {g:34s} built in `#{built}` · wired into the pre-commit hook: **{wired}** "
          f"⇒ it runs over the WHOLE repository on every commit, which is transfer by construction")

print("\n=== (4) THE COUNT, AND THE PARTITION WRITTEN BEFORE THE RUN ===")
n_stat = sum(1 for r in rows if r[4])
n_gate = sum(1 for g in GATES if g.split("/")[-1] in hooktxt)
TOTAL = n_stat + n_gate
print(f"  statistical instruments with >=1 CONTROLLED off-domain application: **{n_stat}/4**")
print(f"  gates running repository-wide on every commit:                      **{n_gate}/2**")
print(f"  ⇒ **{TOTAL}/6 instruments have a demonstrated application beyond the round that built "
      f"them**, against `#908`②'s asserted **2**")
print(f"  transfers whose round has no artifact entry (UNSEEN, not absent): {silent}")

CTRL = pos_ok and neg_ok
if not CTRL:
    VERDICT, WORLD = "UNVERIFIED", ("the detector cannot see — the page's headline stays unaudited "
                                    "and the page must say so")
elif silent >= sum(len(r[3]) for r in rows) / 2:
    VERDICT, WORLD = "UNVERIFIED", ("BLIND — half or more of the transfers have no artifact, so the "
                                    "headline is unauditable")
elif TOTAL > 2:
    VERDICT, WORLD = "OVERTURNED", (f"UNDERSTATED — {TOTAL} of 6, not 2. The page understates its own "
                                    f"residue, which is the error `#899` pre-registered as "
                                    f"'possibly over-modest'")
elif TOTAL == 2:
    VERDICT, WORLD = "CONFIRMED", "CORRECT — `#908`②'s two is the measured number"
else:
    VERDICT, WORLD = "OVERTURNED", f"OVERSTATED — only {TOTAL} of 6; the head is inflated"
print(f"\n  positive control {pos_ok} · negative control {neg_ok}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ AN ARTIFACT RECORDING `positive_ok` DOES NOT PROVE THE CONTROL WAS WELL-SPECIFIED —")
print("     this session alone rebuilt five that were not. The count is of controls that FIRED,")
print("     not of controls that were RIGHT. `[unchallenged]` — `door ③`.")

art = dict(entry=909, round="E03·A110·R347", verdict=VERDICT, world=WORLD,
           n_scripts=len(scripts), n_artifacts=len(arts), n_with_entry=len(ART),
           positive_ok=bool(pos_ok), negative_ok=bool(neg_ok),
           claimed_by_908=2, measured_total=TOTAL, statistical=n_stat, gates=n_gate,
           silent_transfers=silent,
           instruments=[dict(name=r[0], built_in=r[1], built_for=r[2], transfers=r[3],
                             controlled=r[4], no_artifact=r[5]) for r in rows],
           caveat="an artifact recording positive_ok does not prove the control was WELL-SPECIFIED; "
                  "this session alone rebuilt five that were not",
           unchallenged=True)
(OUT / "did_the_instruments_transfer.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'did_the_instruments_transfer.json'}")
