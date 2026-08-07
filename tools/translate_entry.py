"""Translate ONE ledger entry and every anchor that points into it, atomically — or refuse.

`#872`① established the constraint by measuring rather than by starting:
  · the back-catalogue is **~1,410,942 CJK characters across 719 files**;
  · **517 of 519 phrase anchors are Chinese**, and `readme_gate.dangling_anchors` **BLOCKS** when an
    anchor's phrase is not found in the entry it names;
  ⇒ **translating a heading without its README rows in the same commit cannot even be committed.**
  ⇒ the conversion must go **entry by entry, atomically** — never file by file.

`#873`'s prep added two facts the plan needs, both measured here rather than assumed:
  · **643 of 858 entries have a parseable `## Entry N · \\`EAR\\` — title` heading.** The rest use a
    different shape (e.g. a plural `## Entries 43–47`). **The tool must refuse those, not guess**;
  · **10 of 585 anchors quote the entry BODY, not the heading** ⇒ translating headings alone breaks
    those ten. The tool therefore checks every anchor against the WHOLE entry, not the title.

DESIGN, and the part that matters: **this tool's default is a dry run, and it refuses rather than
half-applies.** Before writing anything it verifies that every anchor's NEW phrase is a substring of
the NEW entry text; if any would dangle, nothing is written. That makes the failure mode "no change"
rather than "a broken page", which is the only safe direction for an operation that will be repeated
858 times.

⚠ It does not translate. It **moves a translation into place consistently.** Producing the English
text is a separate act of judgement and is not automated here — a machine may not invent a WHY, and
a ledger entry is mostly WHY.
"""
import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
LEDGER = ROOT / "RETRACTIONS.md"
PAGES = (ROOT / "README.md", ROOT / "README_zh.md")
HEAD = re.compile(r'^## Entry (\d+) · `([^`]+)` — (.+)$', re.M)
ANCHOR = re.compile(r'\[#(\d+)「([^」]+)」\]')


def entry_span(text, n):
    """(start, end) of entry `n` in the ledger, or None. End is the next `## Entry`/`## Entries`."""
    m = re.search(rf'^## Entr(?:y|ies) {n}\b.*$', text, re.M)
    if not m:
        return None
    nxt = re.search(r'^## Entr(?:y|ies) \d+', text[m.end():], re.M)
    return (m.start(), m.end() + nxt.start() if nxt else len(text))


def parse_heading(text, n):
    span = entry_span(text, n)
    if span is None:
        return None
    head_line = text[span[0]:text.index("\n", span[0])]
    m = HEAD.match(head_line)
    return m.groups() if m else None


def anchors_for(n):
    out = []
    for p in PAGES:
        t = p.read_text()
        for m in ANCHOR.finditer(t):
            if m.group(1) == str(n):
                out.append((p, m.group(2)))
    return out


def convert(n, new_title, phrase_map=None, apply=False):
    """Rewrite entry `n`'s heading and every anchor into it. Returns a report; writes only if apply."""
    led = LEDGER.read_text()
    parsed = parse_heading(led, n)
    rep = {"entry": n, "ok": False, "reasons": []}
    if parsed is None:
        rep["reasons"].append("heading does not parse as `## Entry N · `EAR` — title` — REFUSED, "
                              "not guessed (215 of 858 entries have another shape)")
        return rep
    num, ear, old_title = parsed
    span = entry_span(led, n)
    old_entry = led[span[0]:span[1]]
    anchors = anchors_for(n)
    rep.update(ear=ear, old_title=old_title, n_anchors=len(anchors))

    # every anchor's OLD phrase must currently resolve inside this entry, or the ledger is already broken
    for p, phr in anchors:
        if phr not in old_entry:
            rep["reasons"].append(f"anchor 「{phr}」 in {p.name} is ALREADY dangling — fix that first")
    if rep["reasons"]:
        return rep

    pm = dict(phrase_map or {})
    for _p, phr in anchors:
        if phr in pm:
            continue
        if phr == old_title:
            pm[phr] = new_title          # the anchor quotes the WHOLE heading: it follows the heading
        else:
            # ⚠ Caught by attacking the tool: an anchor quoting a FRAGMENT of the heading was being
            # silently remapped to the WHOLE new title. That never dangles — so the safety check
            # passed — but it changes what the anchor means, from "this phrase" to "this entry".
            # A check that cannot fail is not a check. A fragment now REQUIRES an explicit mapping.
            pm[phr] = None
    unmapped = [k for k, v in pm.items() if v is None]
    if unmapped:
        rep["reasons"].append(f"{len(unmapped)} anchor(s) quote a FRAGMENT (of the heading or of the "
                              f"body), not the whole heading, and need an explicit phrase_map — "
                              f"defaulting them to the new title would silently change what the "
                              f"anchor points at: {unmapped}")
        return rep

    new_entry = old_entry.replace(old_title, new_title, 1)
    for old_p, new_p in pm.items():
        if old_p != new_title and old_p in new_entry:
            new_entry = new_entry.replace(old_p, new_p)

    # the check that makes this safe: every NEW phrase must resolve in the NEW entry, or refuse
    for _p, phr in anchors:
        if pm[phr] not in new_entry:
            rep["reasons"].append(f"new phrase 「{pm[phr]}」 would NOT be found in the rewritten "
                                  f"entry — refusing, so the failure mode is 'no change'")
    if rep["reasons"]:
        return rep

    rep["ok"] = True
    rep["plan"] = {"ledger": f"{old_title!r} -> {new_title!r}",
                   "anchors": {k: v for k, v in pm.items()}}
    if apply:
        LEDGER.write_text(led[:span[0]] + new_entry + led[span[1]:])
        for p in PAGES:
            t = p.read_text()
            for old_p, new_p in pm.items():
                t = t.replace(f"[#{n}「{old_p}」]", f"[#{n}「{new_p}」]")
            p.write_text(t)
        rep["applied"] = True
    return rep


def audit():
    """How many entries this tool could handle at all, and where the rest fail. Measured, not guessed."""
    led = LEDGER.read_text()
    nums = sorted({int(m.group(1)) for m in re.finditer(r'^## Entr(?:y|ies) (\d+)', led, re.M)})
    parses = sum(1 for n in nums if parse_heading(led, n))
    anch = {}
    for p in PAGES:
        for m in ANCHOR.finditer(p.read_text()):
            anch.setdefault(int(m.group(1)), []).append(m.group(2))
    body_only = 0
    for n, phrs in anch.items():
        h = parse_heading(led, n)
        if not h:
            continue
        for phr in phrs:
            if phr not in h[2]:
                body_only += 1
    return dict(entries=len(nums), parseable=parses, unparseable=len(nums) - parses,
                anchored_entries=len(anch), anchors=sum(len(v) for v in anch.values()),
                anchors_quoting_body=body_only)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", action="store_true")
    ap.add_argument("--entry", type=int)
    ap.add_argument("--title")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    if a.audit or not a.entry:
        r = audit()
        print("scope of the entry-by-entry conversion, measured:")
        for k, v in r.items():
            print(f"  {k:22s} {v}")
        print("\n  => an entry whose heading does not parse is REFUSED, never guessed;")
        print("     an anchor quoting the body needs an explicit phrase map.")
        sys.exit(0)
    rep = convert(a.entry, a.title, apply=a.apply)
    print(rep)
    sys.exit(0 if rep["ok"] else 2)
