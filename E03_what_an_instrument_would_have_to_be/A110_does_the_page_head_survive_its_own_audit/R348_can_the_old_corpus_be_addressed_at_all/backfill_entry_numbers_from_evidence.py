r"""#910 · E03·A110·R348 — "back-filling is mechanical" was a guess, and it is wrong for two fifths

**COGNITIVE UPDATE CARD**
```
Core Gap        `#909` measured 6/6 instruments transferred and then registered `#909`①: the count
                is a FLOOR, because **442 of 457 persisted artifacts carry no entry number**, so
                anything transferred in E01/E02 is invisible. It also said, in the same breath,
                that **"back-filling entry numbers into old artifacts is MECHANICAL"** — asserted,
                not measured, one sentence after a round whose entire finding was that I assert
                counts about my own work without measuring them.
Why Now         Every count of this project's own residue is a lower bound until the corpus is
                addressable, and `#909`② (a gate for under-claiming) cannot see 97% of the corpus
                without it. **This is the enabling step, and its cost has never been checked.**
Live Worlds     MECHANICAL   >=90% of artifacts resolve to exactly one entry from evidence, and
                             `#909`①'s word stands.
                ⚠ PARTIAL    ⚠ THE UNWELCOME ONE — a large residue cannot be addressed at all, so
                             "mechanical" was wrong and the corpus is permanently only partly
                             auditable. Every residue count stays a floor, forever.
                CONFLICTED   ⚠ META — the routes disagree, which would mean the project's own
                             addressing scheme is inconsistent and neither number can be trusted.
Discriminating  Resolve each artifact by THREE independent routes and compare them: (1) an `entry`
Act             key it already carries, (2) the `#NNN ·` header its own round script declares,
                (3) the ledger heading `## Entry NNN · \`E0x·Ayy·RZZZ\`` whose R matches its
                directory. **Three routes that agree is evidence; one route is a guess.**
Prediction      MECHANICAL -> union >=90%, zero conflicts
Matrix          PARTIAL    -> union <90%
                CONFLICTED -> conflicts on >=5% of the resolvable set
Confound        ⚠ MEASURED BEFORE THE RUN AND IT IS WHY THIS ROUND EXISTS: **66 R-numbers map to
                MORE THAN ONE ledger entry** — `R6` alone maps to thirteen (226, 227, 230, 231,
                235...). That is not a bug: `P16`'s ORIGINAL scheme had a big `R` theme containing
                `rNN` sub-rounds, and the "R = one belief update" convention came later. **The
                numbering scheme CHANGED MID-PROJECT**, so route (3) is ambiguous exactly where the
                corpus is oldest — which is where `#909`① wanted to look.
Controls        positive: `R348`/`R347` must resolve to `#910`/`#909` by every route that has them,
                and agree · negative: an invented R-number resolves to nothing · empty population
                exits 2, never 0
Stopping Rule   One round. Whatever the residue, `#909`①'s word is corrected to it and the artifacts
                that CAN be addressed are annotated. Budget: one round.
Cost            458 artifacts, 868 scripts, one ledger. CPU seconds.
Priority        It is the prerequisite for `#909`②, and its own cost claim is unmeasured.
Expected        If PARTIAL: this project's corpus is permanently only partly auditable, and every
Transform       "what stands" count carries a floor marker rather than a number.
```

⚠⚠ **`#901`①'s REMEDY, EIGHTH USE.** Outcome space = `(union share: >=90% / <90%) × (conflict share:
<5% / >=5%)` — **four cells, all assigned before the run**: `>=90% × <5%`→MECHANICAL ·
`<90% × <5%`→PARTIAL · `anything × >=5%`→CONFLICTED (**both sub-cells**, because inconsistent
addressing makes the union meaningless whatever its size). No cell is unlisted.

`G1` **ESTIMAND**: **the share of persisted artifacts that resolve to exactly one ledger entry from
evidence, by route, and the share that resolve to none.** **Population** every `results/*.json` under
`E0*/A*/R*`. **Instrument** this repository's own ledger, scripts and directory names — ⚠ **`no
second instrument`, `only this one instrument`**, because the object under test is the project's own
addressability. **Baseline** `#909`①'s asserted "mechanical". **Regime** whole repository.

⚠ **"SHOULD THIS ZERO BE ZERO?" — the answer is a SHARE, not an effect**, so the failable object is
the **resolver**: positive control, negative control, empty-population exit. Forcing a null
distribution onto a census would be the wrong shape and saying so is the honest reading (`#907`,
`#909` did the same).

**PRE-REGISTERED KILL — a conditional:**
```
if `R347`->`#909` resolves and agrees across the routes that have it
   AND an invented R-number resolves to nothing:
       read the four-cell table above
else:
       UNVERIFIED — the resolver cannot see, and no share below means anything
```
`G3`: every route reported separately, with its own coverage and its disagreements.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① ⚠ **an ambiguous `R` cannot be disambiguated from the filesystem** — the old scheme's sub-rounds
   were `rNN` INSIDE a big `R`, and the artifacts do not record which; **a machine may not invent a
   WHY (`L80`), so those are left UNRESOLVED rather than assigned to the first entry**;
② **210 of 894 ledger headings carry no round address at all** — governance and production entries —
   so route (3) has a ceiling it cannot exceed;
③ **an artifact's directory is evidence of WHERE it was written, not of WHICH entry claimed it**;
④ no second coder, no second release;
⑤ ⚠ **`[unchallenged]`** — `door ③`.
"""
import ast
import collections
import json
import pathlib
import re
import sys

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

print("=== (0) POPULATION — counted before anything is claimed ===")
led = (ROOT / "RETRACTIONS.md").read_text()
heads = re.findall(r"^## Entry (\d+) · `([^`]+)`", led, re.M)
all_heads = len(re.findall(r"^## Entry ", led, re.M))
arts = sorted(ROOT.glob("E0*/A*/R*/results/*.json"))
scripts = sorted(ROOT.glob("E0*/A*/R*/*.py"))
print(f"  ledger entries {all_heads} · of which carrying a round address {len(heads)} "
      f"({100*len(heads)/all_heads:.1f}%)")
print(f"  persisted artifacts {len(arts)} · round scripts {len(scripts)}")
if not arts or not heads:
    print("EMPTY POPULATION — exit 2, never 0.")
    sys.exit(2)

R2E = collections.defaultdict(set)
for e, addr in heads:
    m = re.search(r"R(\d+)", addr)
    if m:
        R2E[int(m.group(1))].add(int(e))
AMB = {k: sorted(v) for k, v in R2E.items() if len(v) > 1}
print(f"\n  ⚠ **{len(AMB)} of {len(R2E)} R-numbers map to MORE THAN ONE entry** — the confound "
      f"measured before the run")
for k in sorted(AMB)[:4]:
    print(f"     R{k} -> {AMB[k][:6]}{' …' if len(AMB[k]) > 6 else ''}")
print("  ⇒ `P16`'s ORIGINAL scheme had a big `R` theme containing `rNN` sub-rounds; "
      "'R = one belief update' came later. **The numbering changed mid-project.**")

SCRIPT_ENTRY = {}
for s in scripts:
    m = re.search(r"#(\d{3,4})\s*·", s.read_text(errors="replace")[:400])
    if m:
        SCRIPT_ENTRY.setdefault(s.parent, int(m.group(1)))
print(f"  round dirs whose script declares an entry: {len(SCRIPT_ENTRY)}")


def routes(a):
    """Three independent routes to an entry number. Returns {route: entry-or-None}."""
    r = {}
    try:
        d = json.loads(a.read_text())
        r["artifact_key"] = d.get("entry") if isinstance(d.get("entry"), int) else None
    except Exception:
        r["artifact_key"] = None
    r["script_header"] = SCRIPT_ENTRY.get(a.parent.parent)
    m = re.match(r"R(\d+)", a.parent.parent.name)
    n = int(m.group(1)) if m else None
    r["ledger_address"] = (sorted(R2E[n])[0] if n in R2E and len(R2E[n]) == 1 else None)
    return r


print("\n=== (1) POSITIVE CONTROL — `R347` must resolve to `#909` and the routes must agree ===")
ctrl = [a for a in arts if a.parent.parent.name.startswith("R347")]
pos_ok = False
for a in ctrl:
    r = routes(a)
    vals = {v for v in r.values() if v is not None}
    pos_ok = (909 in vals) and len(vals) == 1
    print(f"  {a.parent.parent.name[:36]:38s} {r} ⇒ agree on one value: {len(vals) == 1}")
neg_ok = 99999 not in R2E
print(f"  negative control — invented R99999 resolves to nothing: **{neg_ok}**")

print("\n=== (2) COVERAGE BY ROUTE, AND WHERE THE ROUTES DISAGREE ===")
res, conflicts, unresolved = {}, [], []
cov = collections.Counter()
for a in arts:
    r = routes(a)
    vals = {v for v in r.values() if v is not None}
    for k, v in r.items():
        if v is not None:
            cov[k] += 1
    if not vals:
        unresolved.append(a)
    elif len(vals) > 1:
        conflicts.append((a, r))
    else:
        res[a] = (vals.pop(), [k for k, v in r.items() if v is not None])
for k in ("artifact_key", "script_header", "ledger_address"):
    print(f"  route {k:16s} resolves {cov[k]:4d}/{len(arts)}  ({100*cov[k]/len(arts):5.1f}%)")
UNION = len(res) + len(conflicts)
print(f"\n  **union (>=1 route): {UNION}/{len(arts)} = {100*UNION/len(arts):.1f}%** · "
      f"agreeing on one value: {len(res)} · **conflicting: {len(conflicts)}** · "
      f"**unresolved: {len(unresolved)}**")
for a, r in conflicts[:6]:
    print(f"     ⚠ CONFLICT {a.parent.parent.name[:34]:36s} {r}")
amb_arts = sum(1 for a in unresolved
               if (m := re.match(r"R(\d+)", a.parent.parent.name)) and int(m.group(1)) in AMB)
print(f"  of the unresolved, **{amb_arts}** sit under an AMBIGUOUS R (the old sub-round scheme) — "
      f"⚠ **left unassigned, because a machine may not invent a WHY (`L80`)**")

print("\n=== (3) THE ANNOTATION — a NEW key, never overwriting an existing `entry` ===")
written = 0
for a, (e, via) in res.items():
    try:
        d = json.loads(a.read_text())
    except Exception:
        continue
    if isinstance(d.get("entry"), int):
        continue                      # already addressed; nothing to add
    d["entry_backfilled"] = e
    d["entry_source"] = "+".join(via) + " (`#910`; three routes compared, agreeing)"
    a.write_text(json.dumps(d, indent=1, default=str))
    written += 1
print(f"  artifacts annotated with `entry_backfilled`: **{written}** "
      f"(existing `entry` keys untouched — `L81`)")

print("\n=== (4) THE CONDITIONAL KILL — four cells, all assigned before the run ===")
share = UNION / len(arts)
confl = len(conflicts) / max(UNION, 1)
if not (pos_ok and neg_ok):
    VERDICT, WORLD = "UNVERIFIED", "the resolver cannot see; no share above means anything"
elif confl >= 0.05:
    VERDICT, WORLD = "OVERTURNED", (f"CONFLICTED — the routes disagree on {100*confl:.1f}% of the "
                                    f"resolvable set; the addressing scheme is inconsistent")
elif share >= 0.90:
    VERDICT, WORLD = "CONFIRMED", f"MECHANICAL — {100*share:.1f}% resolve; `#909`①'s word stands"
else:
    VERDICT, WORLD = "OVERTURNED", (f"PARTIAL — only {100*share:.1f}% resolve, so 'back-filling is "
                                    f"MECHANICAL' is WRONG for {100*(1-share):.1f}% of the corpus, "
                                    f"and every residue count stays a floor")
print(f"  union {100*share:.1f}% · conflicts {100*confl:.1f}% of resolvable · "
      f"positive {pos_ok} · negative {neg_ok}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ AN ARTIFACT'S DIRECTORY IS EVIDENCE OF WHERE IT WAS WRITTEN, NOT OF WHICH ENTRY")
print("     CLAIMED IT. The three routes agreeing is what makes this evidence rather than a guess,")
print("     and the ambiguous-R residue is left UNASSIGNED. `[unchallenged]` — `door ③`.")

art = dict(entry=910, round="E03·A110·R348", verdict=VERDICT, world=WORLD,
           n_artifacts=len(arts), n_scripts=len(scripts), n_ledger=all_heads,
           n_headings_with_address=len(heads), ambiguous_R=len(AMB), total_R=len(R2E),
           ambiguous_examples={str(k): AMB[k] for k in sorted(AMB)[:6]},
           coverage={k: cov[k] for k in ("artifact_key", "script_header", "ledger_address")},
           union=UNION, union_share=share, agreeing=len(res), conflicts=len(conflicts),
           unresolved=len(unresolved), unresolved_under_ambiguous_R=amb_arts,
           annotated=written, positive_ok=bool(pos_ok), negative_ok=bool(neg_ok),
           corrects="#909①'s 'back-filling entry numbers into old artifacts is MECHANICAL'",
           unchallenged=True)
(OUT / "addressability.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'addressability.json'}")
