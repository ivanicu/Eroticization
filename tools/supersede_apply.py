r"""`#894` · apply the supersede marker, quoting the LEDGER rather than paraphrasing it.

**L80 — a machine may not invent a WHY.** So this does not write a reason. For each unmarked pair it
lifts the retracting **sentence the ledger already wrote**, trims it, and prepends it to the row as a
quotation with its entry number. If the sentence cannot be extracted, the marker degrades to a bare
pointer — which asserts only what the gate established — rather than to a guess.

**Idempotent**: a row that already carries the marker is skipped, so re-running is safe.
"""
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from supersede_gate import (LEDGER, MARK, PAGES, ROOT, entries, retraction_pairs,
                            row_index)

VERB = r"(?:RETRACTED|OVERTURNED|WITHDRAWN|DOWNGRADED|is\s+backwards|was\s+BACKWARDS)"


def retraction_sentence(led, target, by):
    """The ledger's OWN sentence binding `#target` to a retraction verb, inside entry `by`."""
    body = next((b for n, b in entries(led) if n == by), "")
    flat = re.sub(r"\s+", " ", body)
    # ⚠ a sentence boundary is NOT every `.`: `+0.051` is one token and splitting inside it
    #   produced the marker *“`#888`'s `+0”* on the first run — `realstat` §4, "truncated string
    #   read as data". A terminator is `.` `·` `|` NOT flanked by digits.
    END = r"(?<![0-9])[.·|](?![0-9])"
    parts = [x.strip() for x in re.split(END, flat) if x.strip()]
    cand = [x for x in parts if f"`#{target}`" in x and re.search(VERB, x)]
    if not cand:
        cand = [x for x in parts if f"#{target}" in x and re.search(VERB, x)]
    if not cand:
        return None
    best = min(cand, key=len)
    best = re.sub(r"^[^A-Za-z`⚠*]+", "", best).strip()
    return (best[:300].rstrip() + "…") if len(best) > 300 else best


def main():
    led = LEDGER.read_text()
    pairs = sorted(retraction_pairs(led))
    changed = 0
    for f in PAGES:
        p = ROOT / f
        if not p.exists():
            continue
        txt = p.read_text()
        lines = txt.split("\n")
        for tgt, by in pairs:
            i = row_index(lines, tgt)
            if i is None or MARK in lines[i]:
                continue
            sent = retraction_sentence(led, tgt, by)
            if sent:
                marker = (f"{MARK} IN PART BY `#{by}`** — the ledger's own words: *“{sent}”* "
                          f"⟶ read `#{by}` before using anything below. ")
            else:
                marker = (f"{MARK} IN PART BY `#{by}`** — entry `#{by}` retracts part of this row; "
                          f"the sentence could not be lifted mechanically, so **no reason is invented "
                          f"here (L80)** — read `#{by}`. ")
            body = lines[i]
            body = body.lstrip()
            assert body.startswith("| "), body[:60]
            lines[i] = "| " + marker + body[2:]
            changed += 1
        p.write_text("\n".join(lines))
    print(f"markers applied: {changed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
