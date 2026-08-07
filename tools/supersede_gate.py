r"""`#894` · the ledger retracts a claim; the PAGE keeps stating it.

**The defect this exists for, measured 2026-08-07.** `#893` retracted `#892`'s headline
decomposition — *"96% of D comes from `abpoorw` alone; with its twin removed `abdefctw` is FLAT"* —
and its closing psychological sentence, as **backwards**. Both pages still carried `#892`'s row
stating exactly that, with **no marker of any kind**. A reader landing on that row has no way to
learn it is dead; they would have to read forward to a row they have no reason to open.

**And this is not the first time.** `#885`② was retracted by `#886`, `#886`① by `#887`, `#876`② by
`#884`, `#890`'s framing by `#892` — **and the page carries every one of those rows unmarked.** The
ledger is honest and the deliverable is not, which is the worst of the two orders to get right:
*the ledger is read by me, the page is read by everyone else.*

**Why a gate and not a fix.** `#835`-family lesson, and `P7`'s: the same defect three times means
build infrastructure, not a third patch. A marker applied by hand today is a marker nobody applies
in three rounds' time.

⚠ **P6 PROXY LEDGER — this instrument is a SEARCH, and a search has no positive control unless you
give it one.**
```
PROPERTY     the page states a claim the ledger has since retracted
PROXY        entry M's body contains a retraction verb bound to a `#N` reference (N < M), and the
             page row anchored `(Entry N` carries no supersede marker
IMPLICATION  one direction only: **verb+ref present and marker absent -> that row is unmarked**
             (reliable). The converse fails: a marker does NOT prove the row now reads correctly,
             and a missing verb does NOT prove nothing was retracted -- a retraction written
             without any of the listed verbs is INVISIBLE to this gate.
SAFE SIDE    report only "unmarked"; never certify a row as "correctly superseded".
CONTROLS     positive: `#893` retracting `#892` is a known true pair and MUST be recovered.
             negative: an invented entry number MUST NOT be found.
             empty population -> exit 2, never 0.
```
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
LEDGER = ROOT / "RETRACTIONS.md"
PAGES = ("README.md", "README_zh.md")
MARK = "⚠ **SUPERSEDED"

VERB = r"(?:RETRACTED|OVERTURNED|WITHDRAWN|DOWNGRADED|is\s+backwards|was\s+BACKWARDS)"
# a retraction sentence binds a `#N` reference to a verb within a short window, in either order
PAT = re.compile(r"`#(\d{2,4})`[^\n]{0,160}?" + VERB + r"|" + VERB + r"[^\n]{0,160}?`#(\d{2,4})`")


def entries(text):
    marks = [(int(m.group(1)), m.start()) for m in re.finditer(r"^## Entry (\d+)", text, re.M)]
    out = []
    for i, (n, s0) in enumerate(marks):
        s1 = marks[i + 1][1] if i + 1 < len(marks) else len(text)
        out.append((n, text[s0:s1]))
    return out


def retraction_pairs(text):
    """(target, by) for every retraction sentence. A target must be an EARLIER entry."""
    pairs = set()
    for n, body in entries(text):
        for m in PAT.finditer(body):
            tgt = m.group(1) or m.group(2)
            t = int(tgt)
            if t < n:
                pairs.add((t, n))
    return pairs


ANCHOR = "(Entry {}"


def is_row(line):
    """A TABLE ROW, not a prose mention. `#894`'s first run crashed on `(Entry 716).` in prose."""
    t = line.strip()
    return t.startswith("| ") and t.endswith("|")


def row_index(lines, entry):
    """Index of the table row anchored `(Entry N` where N is followed by a NON-DIGIT.

    ⚠ `#894`'s first run matched `(Entry 68` inside `(Entry 686 · 689 · 690)` and reported a row
    for entry 68 that does not exist. **A substring match on a NUMBER is a measuring instrument
    with no boundary** — the same failure `realstat` §4 records three times in one hour. The
    anchor must terminate."""
    import re as _re
    pat = _re.compile(r"\(Entry " + str(entry) + r"(?![0-9])")
    for k, line in enumerate(lines):
        # ⚠ USE vs MENTION, third defect in this instrument and found by running it: `#894`'s own
        #   page row DESCRIBES the entry-68 bug and writes the anchor `(Entry 68` inside a code
        #   span. The gate then matched its own description and reported a row for entry 68.
        #   An anchor inside backticks is a MENTION. Strip code spans before matching.
        if is_row(line) and pat.search(_re.sub(r"`[^`]*`", "", line)):
            return k
    return None


def row_span(page, entry):
    lines = page.split("\n")
    k = row_index(lines, entry)
    return None if k is None else lines[k]


def main():
    led = LEDGER.read_text()
    pairs = retraction_pairs(led)

    # ---- controls, before any verdict -------------------------------------------------
    pos_ok = any(t == 892 and b == 893 for t, b in pairs)
    neg_ok = not any(t == 9999 or b == 9999 for t, b in pairs)
    print("=== controls ===")
    print(f"  positive: `#893` retracting `#892` recovered -> {'PASS' if pos_ok else 'FAIL'}")
    print(f"  negative: invented entry 9999 absent         -> {'PASS' if neg_ok else 'FAIL'}")
    if not pairs:
        print("EMPTY POPULATION — no retraction pair found at all. Exit 2, never 0.")
        return 2
    if not (pos_ok and neg_ok):
        print("The instrument failed its own control ⇒ every 'unmarked' below is UNVERIFIED.")
        return 2

    pagetexts = {f: (ROOT / f).read_text() for f in PAGES if (ROOT / f).exists()}
    print(f"\n=== {len(pairs)} retraction pair(s) declared in the ledger ===")
    unmarked, absent = [], []
    for tgt, by in sorted(pairs):
        for f, txt in pagetexts.items():
            row = row_span(txt, tgt)
            if row is None:
                absent.append((tgt, by, f))
            elif MARK not in row:
                unmarked.append((tgt, by, f))
    for tgt, by, f in sorted(unmarked):
        print(f"  ⚠ UNMARKED  entry {tgt} (retracted by {by})  in {f}")
    for tgt, by, f in sorted(absent):
        print(f"  · no row     entry {tgt} (retracted by {by})  in {f} — cannot be marked, not a hit")
    print(f"\n  pairs {len(pairs)} · rows unmarked {len(unmarked)} · targets with no page row {len(absent)}")
    print("  ⚠ Scope: only retractions written with one of the listed verbs near a `#N` reference are")
    print("     visible. A retraction phrased otherwise is UNSEEN, not cleared. And a marker present")
    print("     does NOT certify the row now reads correctly — this gate reports absence only.")
    return 1 if unmarked else 0


if __name__ == "__main__":
    sys.exit(main())
