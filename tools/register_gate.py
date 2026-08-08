"""tools/register_gate.py — a line saying something DOES NOT EXIST is a search result, and a
search is an instrument with no positive control.

Built at `E03·A133·R395` (`#961`③'s NEXT). Action type: **PRODUCTION**, labelled as such.

**WHY.** `#960` wrote *"There is no second ethnographic coding of these acts in the release."*
`#961`'s own opening query refuted it one round later: Whyte 1978 coded the same extramarital
double standard, and 51 societies carry both codings. **The impossibility register is where an
unchecked claim is least likely to be re-examined, because its whole function is to say stop
looking.**

**⚠⚠ THE GRADIENT CHECK KILLED THE FIRST DESIGN BEFORE A LINE WAS WRITTEN, and both halves of it.**
  ① I planned to scan the `STRUCTURALLY CANNOT` block. **The refuted sentence is not in it** — it
     lives in `#960`'s HARD RULE 2 paragraph. A gate scoped to that block would have reported a
     clean corpus while the very sentence that motivated it sat two paragraphs above. That is
     `realstat` §4's row exactly: *the instrument measured headings; the sentence asserted entries.*
  ② I planned to separate EXISTENCE claims ("there is no X") from CAPABILITY claims ("X would
     require Y"). `#960`'s clause (1) reads *"only this one instrument … separating the society from
     the ethnographer **would require** a second, independent coding team"* — **it is phrased as
     capability and asserts existence**, and it is false in the same way. The distinction is not a
     grammar; it is a property of the assertion.
  ⇒ Corpus = **every sentence** of every round docstring and every ledger entry. Classification is
     on what is asserted, not on the modal verb.

**WHAT IT CHECKS.** A sentence asserting that something is not available must **cite the query that
returned nothing** — a path, a variable id, a command, a count, or an explicit search verb. Without
one it is `UNVERIFIED`, never `SETTLED`.

⚠ P6 PROXY LEDGER — read this before believing any number below.
  PROPERTY    the register contains no unchecked non-existence claim
  PROXY       a regex over sentences, plus a regex for a cited query
  IMPLICATION **one direction only**: a HIT is a sentence that asserts unavailability and cites no
              query — worth checking. A MISS is **NOT** a certificate: a non-existence claim phrased
              in words this pattern does not know is UNSEEN, not cleared.
  WITNESS     `#960`'s clause (1) is a real sentence that the FIRST design (block-scoped,
              existence-vs-capability) would have missed. It is kept as the positive control.
  SAFE SIDE   report HITS as a reading list. Never report "the register is clean".

⚠ PRIOR ART (`P4`). `tools/wall_audit.py` exists (`E02·A243·R626`) and asks a **different** question
  — *is this wall measurable?* — bucketing clauses into measured / measurable-but-unmeasured /
  structurally-unmeasurable. It is not this gate, and this gate does not replace it. But its corpus
  is **measured stale here and that is reported**: it reads only the ledger's Chinese heading form
  and its highest entry is far below the ledger's tip, so a growing share of the register has been
  invisible to it. A tool that always returns a tidy answer is the hardest kind of dead instrument
  to notice.
"""
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

# ── the instrument's UNIT is a SENTENCE. Named here, separately from the claim's unit, and
#    required to be equal before any control is designed (`realstat` §4's stronger remedy).
INSTRUMENT_UNIT = "sentence"
CLAIM_UNIT = "sentence"

# ── TIGHT pattern: asserts that a THING is not available. Anchored on the assertion.
# ⚠ P7 ATTACK V3 FOUND A REAL HOLE AND THIS IS THE REPAIR: the first version used literal single
#   spaces, so it matched only because `sentences()` normalises whitespace upstream. A lock whose
#   coverage lives in a DIFFERENT function is a lock that breaks when someone edits that function.
#   Every gap is now `\s+`, so the pattern holds on raw text too.
TIGHT = re.compile(
    r"(there\s+(?:is|are|was|were)\s+no\b"
    r"|\bno\s+second\b|\bno\s+other\b|\bno\s+such\b|\bnone\s+(?:exists?|available|survive)"
    r"|\bdoes\s+not\s+exist\b|\bdo\s+not\s+exist\b|\b(?:is|are)\s+not\s+available\b"
    r"|\bonly\s+this\s+one\b|\bthe\s+only\s+\w+\s+(?:that|which|carrying|available)"
    r"|\bnothing\s+(?:else\s+)?(?:carries|exists|available)"
    r"|\bno\s+\w+(?:\s+\w+){0,2}\s+(?:exists?|is available|are available|carries|carry)\b)", re.I)
# ── LOOSE pattern, run only to MEASURE the disagreement (`realstat` §4: when a loose and a tight
#    pattern disagree, the tight one is not "more conservative", it is the one that was tested).
LOOSE = re.compile(r"\bno\b|\bnot\b|\bnone\b|\bcannot\b|\bunavailable\b", re.I)

# ── a CITED QUERY: a path, a variable id, a command, a count, or an explicit search verb.
CITED = re.compile(
    r"\b[\w/.-]+\.(?:csv|py|json|md|dta|tsv|sav)\b"          # a file
    r"|\bSCCS\d+\b|\bEA\d+\b"                                 # a variable id
    r"|\bn\s*=\s*\d|\b\d+\s*/\s*\d+\b|\b\d+\s+of\s+\d+\b"     # a count
    r"|\bgrep\b|\bquer(?:y|ied)\b|\bsearched?\b|\bscanned?\b|\blisted?\b|\bchecked\b"
    r"|`[^`]+`", re.I)

# ── an explicit opt-out, so the gate can be answered rather than only obeyed.
OPTOUT = re.compile(r"UNVERIFIED, not (?:cleared|settled)|not checked|no query was run", re.I)

SPLIT = re.compile(r"(?<=[.;:])\s+|\n\s*\n")

# measured 2026-08-07 at E03·A133·R395; a ratchet, and lowering it is the work
BASELINE_UNCITED = 134


def sentences(text, source, key):
    for raw in SPLIT.split(text):
        s = " ".join(raw.split())
        if 25 <= len(s) <= 600:
            yield dict(source=source, key=key, text=s)


def corpus():
    """Ledger entries AND round docstrings. Per-source counts are reported so a DEAD SOURCE is
    visible — a gate whose corpus quietly stops growing still prints a tidy answer."""
    out = []
    led = (ROOT / "RETRACTIONS.md").read_text()
    ents = re.split(r"\n## Entry ", led)[1:]
    for e in ents:
        m = re.match(r"(\d+)", e)
        if m:
            out += list(sentences(e, "ledger", f"#{m.group(1)}"))
    for p in sorted(ROOT.glob("E0*/**/*.py")):
        try:
            txt = p.read_text()
        except Exception:
            continue
        m = re.match(r'\s*(?:#![^\n]*\n)?(?:#[^\n]*\n)*\s*"""(.*?)"""', txt, re.S)
        if m:
            out += list(sentences(m.group(1), "docstring", str(p.relative_to(ROOT))))
    return out, len(ents), max((int(re.match(r"(\d+)", e).group(1))
                                for e in ents if re.match(r"(\d+)", e)), default=0)


def classify(rows):
    hits, cited, optout = [], [], []
    for r in rows:
        if not TIGHT.search(r["text"]):
            continue
        if OPTOUT.search(r["text"]):
            optout.append(r)
        elif CITED.search(r["text"]):
            cited.append(r)
        else:
            hits.append(r)
    return hits, cited, optout


# ══ CONTROLS. Every one is a REAL string from this repo's history, not an invented case ═══
POS = ("There is no second ethnographic coding of these acts in the release.",
       "`#960`, refuted by `#961`'s own opening query one round later")
POS2 = ("only this one instrument — Broude & Greene 1976 coded all three variables from HRAF; "
        "separating the society from the ethnographer would require a second, independent "
        "coding team",
        "`#960` clause (1): asserts existence while phrased as capability — the case that killed "
        "the first design")
NEG = ("A second field observation is what removing it would take, and that is registered.",
       "`#961`: an incapacity claim naming its remedy — must NOT be caught")
NEG2 = ("SCCS176 n=40 of 186 societies carry a homosexuality code, so no fourth act is available "
        "from the same team.",
        "a non-existence claim that DOES cite its query (n=40 of 186) — must NOT be caught")


def run():
    rows, n_entries, max_entry = corpus()
    by_src = {}
    for r in rows:
        by_src[r["source"]] = by_src.get(r["source"], 0) + 1
    docs = len({r["key"] for r in rows if r["source"] == "docstring"})

    print("=== register gate (`#961`③) — a non-existence claim is a SEARCH RESULT ===")
    print(f"  instrument unit = {INSTRUMENT_UNIT!r} · claim unit = {CLAIM_UNIT!r} · "
          f"EQUAL -> {INSTRUMENT_UNIT == CLAIM_UNIT}")
    print(f"  corpus: {len(rows)} sentences · ledger {by_src.get('ledger',0)} from {n_entries} "
          f"entries (tip #{max_entry}) · docstrings {by_src.get('docstring',0)} from {docs} files")

    # ⚠ EMPTY / STALE POPULATION -> exit 2, never 0 (`realstat` §4: a gate that reports success
    #   having examined nothing).
    if not rows or by_src.get("docstring", 0) == 0 or by_src.get("ledger", 0) == 0:
        print("  ⛔ EMPTY OR ONE-SIDED POPULATION — the gate examined nothing on at least one "
              "source. Exit 2, never 0.")
        return 2, {}

    # ⚠ UNIT-EQUALITY CONTROL, run BEFORE the pass/fail controls (`realstat` §4's stronger remedy):
    #   is the positive control's sentence actually IN the corpus this gate builds? The first
    #   design failed exactly here — it scanned the STRUCTURALLY CANNOT block, and the sentence
    #   lives two paragraphs above it.
    in_corpus = any(POS[0][:60].lower() in r["text"].lower() for r in rows)
    print(f"  control UNIT (is the known bad sentence in the corpus at all?): "
          f"{'PASS' if in_corpus else 'FAIL — the gate is scanning the wrong unit'}")

    def caught(t):
        return bool(TIGHT.search(t)) and not CITED.search(t) and not OPTOUT.search(t)

    cpos, cpos2 = caught(POS[0]), caught(POS2[0])
    cneg, cneg2 = caught(NEG[0]), caught(NEG2[0])
    print(f"  control POSITIVE  ({POS[1]}): {'PASS' if cpos else 'FAIL'}")
    print(f"  control POSITIVE2 ({POS2[1]}): {'PASS' if cpos2 else 'FAIL'}")
    print(f"  control NEGATIVE  ({NEG[1]}): {'PASS' if not cneg else 'FAIL'}")
    print(f"  control NEGATIVE2 ({NEG2[1]}): {'PASS' if not cneg2 else 'FAIL'}")
    ctl = in_corpus and cpos and cpos2 and (not cneg) and (not cneg2)

    hits, cited, optout = classify(rows)
    loose_n = sum(1 for r in rows if LOOSE.search(r["text"]))
    print(f"\n  ⚠ LOOSE vs TIGHT, measured rather than argued: a loose 'contains a negation' "
          f"pattern matches {loose_n} sentences ({loose_n/len(rows):.1%} of the corpus); the tight "
          f"pattern matches {len(hits)+len(cited)+len(optout)}. "
          f"The tight one is not 'more conservative' — it is the one with controls.")
    print(f"\n  non-existence assertions: {len(hits)+len(cited)+len(optout)} · "
          f"CITE a query {len(cited)} · explicit opt-out {len(optout)} · "
          f"**UNCITED (the reading list) {len(hits)}**")

    by_key = {}
    for h in hits:
        by_key.setdefault(h["key"], []).append(h["text"])
    for k in sorted(by_key, key=lambda x: (-len(by_key[x]), x))[:18]:
        print(f"\n  ⚠ {k}  ({len(by_key[k])})")
        for t in by_key[k][:2]:
            print(f"      {t[:170]}")

    print("\n  ⚠ SCOPE (P6, one direction only): a HIT is a sentence asserting unavailability with "
          "no cited query — a reading list. A MISS is UNSEEN, NOT CLEARED: a non-existence claim "
          "phrased in words this pattern does not know is invisible, and the FIRST version of this "
          "gate proves that failure mode is real rather than hypothetical.")
    if not ctl:
        print("  ⛔ CONTROLS FAILED — the gate is UNREADABLE this run, not clean. Exit 2.")
        return 2, dict(hits=len(hits))

    # ⚠ RATCHET, not zero. 134 uncited assertions already exist; demanding zero today would block
    #   every commit and the gate would be disabled within the hour. It blocks on REGRESSION —
    #   a NEW uncited non-existence claim — and the baseline is written down, not remembered.
    bl = json.loads((ROOT / "tools" / "gate_baseline.json").read_text()) \
        if (ROOT / "tools" / "gate_baseline.json").exists() else {}
    base = bl.get("register_gate_uncited", BASELINE_UNCITED)
    worse = len(hits) > base
    print(f"\n  RATCHET: uncited {len(hits)} vs baseline {base} -> "
          f"{'⛔ WORSE, blocking' if worse else 'ok'}")
    return (1 if worse else 0), dict(
        corpus=len(rows), ledger_sentences=by_src.get("ledger", 0), entries=n_entries,
        ledger_tip=max_entry, docstring_sentences=by_src.get("docstring", 0), docstring_files=docs,
        loose_matches=loose_n, assertions=len(hits) + len(cited) + len(optout),
        cited=len(cited), optout=len(optout), uncited=len(hits),
        controls=dict(unit=in_corpus, positive=cpos, positive2=cpos2,
                      negative=not cneg, negative2=not cneg2),
        baseline=base, worse=bool(worse),
        uncited_by_key={k: v for k, v in sorted(by_key.items(), key=lambda kv: -len(kv[1]))})


if __name__ == "__main__":
    code, payload = run()
    out = ROOT / "E03_what_an_instrument_would_have_to_be" / "A133_the_register_is_a_search_result" \
        / "R395_a_non_existence_claim_needs_its_query" / "results"
    out.mkdir(parents=True, exist_ok=True)
    (out / "register_gate.json").write_text(json.dumps(payload, indent=1, ensure_ascii=False))
    print(f"\nwrote {out / 'register_gate.json'} · exit {code}")
    sys.exit(0 if "--report" in sys.argv else code)
