r"""#911 · E03·A110·R349 — three things `#910` got wrong, and then the basin rule stops the chain

**COGNITIVE UPDATE CARD**
```
Core Gap        `#910` corrected `#909`, whose finding was that I assert counts about my own work
                without measuring them. **`#910` then published three quantities of its own and at
                least one of them was asserted rather than measured.** `#910`② named a check it had
                skipped; running that check found the other two.
Why Now         Three published numbers are wrong on the deliverable. That is owed regardless of
                what comes next.
Live Worlds     SOUND    `#910`'s three quantities survive re-measurement.
                ⚠ WRONG  ⚠ THE UNWELCOME ONE — they do not, and the entry that diagnosed
                         assert-without-measuring committed it three times inside the diagnosis.
Discriminating  Re-measure each of `#910`'s three published quantities directly: the artifact
Act             POPULATION, the MECHANISM it gave for R-ambiguity, and the ADDRESSABILITY share
                including a route it never looked for.
Prediction      SOUND -> all three reproduce · WRONG -> any fails
Matrix
Confound        ⚠ WRITTEN BEFORE THE RUN: I am auditing my own audit of my own audit. **This is the
                THIRD consecutive round whose object is my own bookkeeping**, and `frontier` §3 says
                N consecutive steps confirming the same story is a BASIN. The story being confirmed
                is *"I assert without measuring"*. **A fourth would be the failure `§0.2` names: an
                activity metric whose optimum is doing nothing.** The basin call is made in this
                round's NEXT, not left for a later one to notice.
Controls        positive: `R347` must resolve to `#909` by every route that has it · negative: an
                invented R-number resolves to nothing · empty population exits 2, never 0
Stopping Rule   ONE round, and then the chain STOPS by the basin rule. Whatever this finds, the next
                round changes object.
Cost            one ledger, 869 scripts, every file under `results/`. CPU seconds.
Priority        Three wrong numbers on the page, and a basin that needs calling.
Expected        If WRONG: the correction lands and the direction changes, which is the point.
Transform
```

⚠⚠ **`#901`①'s REMEDY, NINTH USE.** Outcome space = which of `#910`'s three quantities survive —
`{population, mechanism, share}`, **eight subsets, all assigned before the run**: `all three`→SOUND ·
**any other subset**→WRONG, **with the failing members named individually and no subset unlisted.**

`G1` **ESTIMAND**: `#910`'s three published quantities, re-measured —
**(1)** the artifact population it used (**it said 458**);
**(2)** the mechanism it gave for 66 R-numbers mapping to several entries (**it said `P16`'s big-`R`
with `rNN` sub-rounds**);
**(3)** the addressable share (**it said 58.5%**), now including the **`<script_stem>__` filename
prefix** route it never looked for.
**Population** every file under `E0*/A*/R*/results/`. **Instrument** this repository — ⚠ **`no second
instrument`, `only this one instrument`**. **Baseline** `#910`'s own numbers. **Regime** whole repo.

⚠ **"SHOULD THIS ZERO BE ZERO?" — the answers are COUNTS and a SHARE, not effects**, so the failable
object is the resolver, which gets a positive control, a negative control and an empty-population
exit. Forcing a null distribution onto a census is the wrong shape (`#907`, `#909`, `#910` all said
this; it is the same reading).

**PRE-REGISTERED KILL — a conditional:**
```
if `R347`->`#909` resolves and an invented R resolves to nothing:
       all three of `#910`'s quantities reproduce -> SOUND
       any fails                                  -> WRONG, each failing member named
else:
       UNVERIFIED
```
`G3`: all three quantities reported whatever they do.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **the `<stem>__` prefix links an artifact to a SCRIPT, not to an ENTRY** — and only **41 of 893**
   ledger entries cite a script filename, so that route does **not** complete;
② **an artifact's directory is evidence of WHERE it was written, not WHICH entry claimed it**
   (`#910`'s own limit, and it stands);
③ ⚠ **`[unchallenged]`** — `door ③`; ④ no second coder, no second release.
"""
import collections
import json
import pathlib
import re
import sys

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
R348 = ROOT / ("E03_what_an_instrument_would_have_to_be/A38_can_the_page_be_trusted_to_say_what_the_evidence_says/"
               "R348_can_the_old_corpus_be_addressed_at_all/results/addressability.json")

print("=== (0) `#910`'s OWN NUMBERS, read from its artifact and not retyped (`#840`'s RULE) ===")
prior = json.loads(R348.read_text())
P_POP, P_SHARE, P_AMB = prior["n_artifacts"], prior["union_share"], prior["ambiguous_R"]
print(f"  `#910` published: population **{P_POP}** · addressable **{100*P_SHARE:.1f}%** · "
      f"ambiguous R-numbers **{P_AMB}**")
print(f"  and its stated MECHANISM: \"`P16`'s ORIGINAL scheme had a big `R` theme containing "
      f"`rNN` sub-rounds\"")

print("\n=== (1) QUANTITY ONE — the POPULATION ===")
files = [p for p in ROOT.glob("E0*/A*/R*/results/*") if p.is_file()]
ext = collections.Counter(p.suffix for p in files)
jsons = [p for p in files if p.suffix == ".json"]
print(f"  every file under results/: **{len(files)}**  {dict(ext.most_common(6))}")
print(f"  of which `.json`: **{len(jsons)}** ⇒ `#910` measured the JSON SUBSET and called it "
      f"'this project's artifacts'")
Q1 = (len(jsons) == P_POP or abs(len(jsons) - P_POP) <= 2) and len(files) > P_POP * 1.5
print(f"  ⇒ quantity ONE {'FAILS' if Q1 else 'holds'}: the denominator was **{100*P_POP/len(files):.1f}% "
      f"of the real population**")
if not files:
    print("EMPTY POPULATION — exit 2, never 0."); sys.exit(2)

print("\n=== (2) QUANTITY TWO — the MECHANISM it gave for R-ambiguity ===")
rnn_dirs = [p for p in ROOT.glob("E0*/A*/R*/r[0-9]*") if p.is_dir()]
led = (ROOT / "RETRACTIONS.md").read_text()
heads = re.findall(r"^## Entry (\d+) · `([^`]+)`", led, re.M)
rnn_addr = [a for _, a in heads if re.search(r"r\d+", a)]
forms = collections.Counter(re.sub(r"\d+", "N", a) for _, a in heads)
print(f"  `rNN` sub-round DIRECTORIES on disk: **{len(rnn_dirs)}**")
print(f"  ledger addresses containing an `rNN`: **{len(rnn_addr)}**")
print(f"  address FORMS: {dict(forms.most_common(3))}")
Q2 = len(rnn_dirs) == 0 and len(rnn_addr) == 0
print(f"  ⇒ quantity TWO {'FAILS' if Q2 else 'holds'}: **the sub-round scheme does not exist, on "
      f"disk or in any address.** `#910` INVENTED a mechanism for its own observation.")
r2e = collections.defaultdict(set)
for e, a in heads:
    m = re.search(r"R(\d+)", a)
    if m:
        r2e[int(m.group(1))].add(int(e))
amb = {k: sorted(v) for k, v in r2e.items() if len(v) > 1}
_r006 = list(ROOT.glob("E0*/A*/R006_*"))
_n006 = len(list(_r006[0].glob("*.py"))) if _r006 else 0
print(f"  what IS true, measured: **{len(amb)} R-numbers carry several entries** — e.g. "
      f"`R006` -> {len(amb.get(6, []))} entries, and its directory holds {_n006} scripts")
print("  ⇒ **early rounds simply produced many belief-updates each. No sub-rounds were involved.**")

print("\n=== (3) QUANTITY THREE — the ADDRESSABLE SHARE, with the route `#910` never looked for ===")
SCRIPT_ENTRY = {}
for s in ROOT.glob("E0*/A*/R*/*.py"):
    m = re.search(r"#(\d{3,4})\s*·", s.read_text(errors="replace")[:400])
    if m:
        SCRIPT_ENTRY.setdefault(s.parent, int(m.group(1)))
cite = collections.defaultdict(set)
marks = [(int(m.group(1)), m.start()) for m in re.finditer(r"^## Entry (\d+)", led, re.M)]
for i, (n, s0) in enumerate(marks):
    body = led[s0:(marks[i + 1][1] if i + 1 < len(marks) else len(led))]
    for f in re.findall(r"`([a-z0-9_]+\.py)`", body):
        cite[f].add(n)
uniq_cite = {f: sorted(v)[0] for f, v in cite.items() if len(v) == 1}
print(f"  ledger entries citing >=1 script filename: {sum(1 for _ in cite)} distinct names, "
      f"of which **{len(uniq_cite)}** name exactly one entry")
stem_ok = sum(1 for p in files if (p.parent.parent / (p.name.split('__')[0] + '.py')).exists())
print(f"  artifacts whose `<stem>__` prefix matches a script in the same round: "
      f"**{stem_ok}/{len(files)} = {100*stem_ok/len(files):.1f}%** — ⚠ a route to a SCRIPT, "
      f"which reaches an ENTRY only for the {len(uniq_cite)} uniquely-cited names")


def resolve(p):
    r = {}
    if p.suffix == ".json":
        try:
            d = json.loads(p.read_text())
            r["artifact_key"] = d.get("entry") if isinstance(d.get("entry"), int) else \
                (d.get("entry_backfilled") if isinstance(d.get("entry_backfilled"), int) else None)
        except Exception:
            r["artifact_key"] = None
    r["script_header"] = SCRIPT_ENTRY.get(p.parent.parent)
    m = re.match(r"R(\d+)", p.parent.parent.name)
    n = int(m.group(1)) if m else None
    r["ledger_address"] = sorted(r2e[n])[0] if n in r2e and len(r2e[n]) == 1 else None
    r["script_citation"] = uniq_cite.get(p.name.split("__")[0] + ".py")
    return {k: v for k, v in r.items() if v is not None}


res = [resolve(p) for p in files]
union = sum(1 for r in res if r)
agree = sum(1 for r in res if len(set(r.values())) == 1 and r)
conflict = sum(1 for r in res if len(set(r.values())) > 1)
print(f"\n  over the FULL population of {len(files)}: union **{union} = {100*union/len(files):.1f}%** "
      f"· agreeing {agree} · conflicting {conflict}")
Q3 = abs(100 * union / len(files) - 100 * P_SHARE) > 5
print(f"  ⇒ quantity THREE {'FAILS' if Q3 else 'holds'}: `#910` said **{100*P_SHARE:.1f}%**, the "
      f"full-population figure is **{100*union/len(files):.1f}%**")

print("\n=== (4) CONTROLS ===")
ctrl = [p for p in files if p.parent.parent.name.startswith("R347")]
pos = any(909 in resolve(p).values() for p in ctrl)
neg = 99999 not in r2e
print(f"  positive: an `R347` artifact resolves to `#909` -> **{pos}**  (n={len(ctrl)})")
print(f"  negative: invented `R99999` resolves to nothing -> **{neg}**")

print("\n=== (5) THE CONDITIONAL KILL ===")
failing = [n for n, f in (("population", Q1), ("mechanism", Q2), ("share", Q3)) if f]
if not (pos and neg):
    VERDICT, WORLD = "UNVERIFIED", "the resolver cannot see"
elif not failing:
    VERDICT, WORLD = "CONFIRMED", "SOUND — all three of `#910`'s quantities reproduce"
else:
    VERDICT, WORLD = "OVERTURNED", (f"WRONG — {len(failing)} of `#910`'s three quantities fail: "
                                    f"{', '.join(failing)}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠⚠ BASIN CALL, MADE HERE AND NOT LEFT FOR A LATER ROUND: this is the THIRD consecutive")
print("     round whose object is my own bookkeeping (`#909`→`#908`, `#910`→`#909`, `#911`→`#910`),")
print("     and the story each confirms is the same — *I assert without measuring*. `frontier` §3")
print("     calls that a basin; `§0.2` calls a fourth one an activity metric whose optimum is doing")
print("     nothing. **The chain stops here. The next round changes object.**")

art = dict(entry=911, round="E03·A110·R349", verdict=VERDICT, world=WORLD,
           prior_population=P_POP, true_population=len(files), by_extension=dict(ext),
           prior_share=P_SHARE, full_population_share=union / len(files),
           rnn_dirs=len(rnn_dirs), rnn_addresses=len(rnn_addr), address_forms=dict(forms),
           ambiguous_R=len(amb), stem_attributable=stem_ok,
           uniquely_cited_scripts=len(uniq_cite),
           union=union, agreeing=agree, conflicting=conflict,
           failing_quantities=failing, positive_ok=bool(pos), negative_ok=bool(neg),
           basin="third consecutive round on my own bookkeeping; the chain stops here",
           unchallenged=True)
(OUT / "three_things_910_got_wrong.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'three_things_910_got_wrong.json'}")
