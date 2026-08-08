#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verdict_contradiction_gate — a round's verdict may not be contradicted by its own control rows.

Built at `#931` to pay `#930`(3). In `#930` the gate printed
    OVERTURNED -- controls sound, and the pre-registered threshold fired AGAINST the expectation
while, three lines above it, its OWN multiplicity row read **cells surviving 0** and its OWN ceiling
row read **DOES NOT survive the ceiling rescale**. I noticed by reading. Nothing in the machinery did.

⚠ WHY `three_valued()` CANNOT SEE THIS, read from the source rather than guessed: it consults only
the pass/fail BOOLEAN of each row. A control can PASS in the sense that "the correction ran" while
its DETAIL TEXT reports evidence against the very verdict being printed. `multiplicity_control`
passes whenever BH was applied — including when BH leaves nothing standing.

⚠⚠ AND HALF THE DEFECT WAS MINE, NOT THE LIBRARY'S. `#930`'s ceiling row asserted
    `ceiling_survives or abs(med_head) <= 0.5 * abs(med_raw)`
whose two branches cover the magnitude condition completely: it can only fail on a SIGN mismatch, so
for the thing it claimed to test it **could not fail**. That is `#916`(3)'s family — naming a control
after what I meant rather than what it does — and no library change prevents it. **A tool that reads
the detail text catches both halves, which is why this is a detail-text gate and not a patch to
`three_valued()`.**

⚠ P6 PROXY LEDGER
  PROPERTY    the verdict is consistent with what the round's own controls found
  PROXY       a control row's detail matches a CONTRADICTION phrase while the verdict asserts a
              positive finding (OVERTURNED / CONFIRMED / ALL GATES PASS)
  IMPLICATION one direction only: **a match ⇒ the round really does contain a contradiction**
              (reliable). No match is NOT proof of consistency — a round can contradict itself in
              words this pattern does not know.
  WITNESS     `#930`: verdict OVERTURNED, multiplicity `cells surviving 0`, ceiling `DOES NOT survive`.
  SAFE SIDE   report contradiction only. **Never certify a verdict as sound.**

Exit codes: 0 clean · 1 at least one contradiction · 2 empty population or a control failed.
"""
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CUTOFF = 931          # blocks from this entry on; earlier rounds are named, never retro-blocked

POSITIVE_VERDICT = re.compile(r"OVERTURNED|CONFIRMED|ALL GATES PASS")
CONTRADICTION = [
    (re.compile(r"cells?\s+surviving\s+0\b|存活\s*0\b|存活\s*无"), "the multiplicity correction left NOTHING standing"),
    (re.compile(r"DOES NOT survive|does not survive"), "a robustness rescale was reported as not survived"),
    (re.compile(r"未存活.*\]\s*$"), None),        # informational only; not a contradiction on its own
]
REQUIRED = "gates=[(name, ok, ...)] rows AND a `gate_verdict` string"


def audit(paths):
    hits, unreadable, ok = [], [], []
    for p in paths:
        try:
            a = json.load(open(p))
        except Exception:                                   # noqa: BLE001
            unreadable.append((p, "unparseable"))
            continue
        if not isinstance(a, dict):
            unreadable.append((p, "not an object"))
            continue
        verdict = a.get("gate_verdict") or a.get("verdict")
        rows = a.get("gates")
        if not verdict or not isinstance(rows, list) or not rows:
            unreadable.append((p, "no `gate_verdict` and/or no `gates` rows"))
            continue
        if not POSITIVE_VERDICT.search(str(verdict)):
            ok.append(dict(path=str(p), entry=a.get("entry"), verdict=str(verdict)[:40]))
            continue
        found = []
        for r in rows:
            detail = " ".join(str(x) for x in r) if isinstance(r, (list, tuple)) else str(r)
            for pat, why in CONTRADICTION:
                if why and pat.search(detail):
                    name = r[0] if isinstance(r, (list, tuple)) and r else "?"
                    found.append((str(name)[:60], why))
        if found:
            hits.append(dict(path=str(p), entry=a.get("entry"), verdict=str(verdict)[:60], found=found))
        else:
            ok.append(dict(path=str(p), entry=a.get("entry"), verdict=str(verdict)[:40]))
    return hits, unreadable, ok


def _unit_tests():
    """Positive: `#930`'s own shape must be CAUGHT. Negative: the same shape with survivors must
    not be. Blind: an artifact with no gate rows is UNREADABLE, never PASS."""
    import tempfile
    bad = dict(entry=9101, gate_verdict="OVERTURNED -- controls sound",
               gates=[["the whole grid", True, True, "method=bh q=0.05 · cells tested 3 · cells surviving 0"],
                      ["ceiling", True, True, "r raw +0.460 vs +0.226 — DOES NOT survive the ceiling rescale"]])
    good = dict(entry=9102, gate_verdict="OVERTURNED -- controls sound",
                gates=[["the whole grid", True, True, "method=bh q=0.05 · cells tested 3 · cells surviving 3"]])
    blind = dict(entry=9103, gate_verdict="ALL GATES PASS")
    out = []
    with tempfile.TemporaryDirectory() as td:
        for nm, obj in (("bad", bad), ("good", good), ("blind", blind)):
            q = pathlib.Path(td) / f"{nm}.json"
            q.write_text(json.dumps(obj))
            out.append(q)
        h, u, o = audit(out)
    return (len(h) == 1 and h[0]["entry"] == 9101 and len(h[0]["found"]) == 2,
            any(x["entry"] == 9102 for x in o),
            len(u) == 1 and u[0][1].startswith("no `gate_verdict`"))


def main():
    print("=== verdict-contradiction gate (`#931`, pays `#930`(3)) ===")
    pc, nc, bc = _unit_tests()
    print(f"  control positive (`#930`'s shape: OVERTURNED beside 'surviving 0' and 'DOES NOT "
          f"survive') must be CAUGHT: {'PASS' if pc else 'FAIL'}")
    print(f"  control negative (same verdict, survivors present) must NOT be caught:            "
          f"{'PASS' if nc else 'FAIL'}")
    print(f"  control blind    (no gate rows -> UNREADABLE, never PASS):                        "
          f"{'PASS' if bc else 'FAIL'}")
    if not (pc and nc and bc):
        print("  ⛔ a control failed — this gate's verdict on the corpus is INADMISSIBLE")
        return 2

    paths = sorted(ROOT.glob("E0*/A*/R*/results/*.json"))
    if not paths:
        print("  ⛔ EMPTY POPULATION — exit 2, never 0")
        return 2
    hits, unreadable, ok = audit(paths)
    print(f"\n  artifacts {len(paths)} · with a verdict AND gate rows {len(hits) + len(ok)} · "
          f"⚠ UNREADABLE {len(unreadable)} · CONTRADICTED {len(hits)}")
    blocking = []
    for h in sorted(hits, key=lambda x: -(x["entry"] or 0)):
        e = h["entry"] or 0
        tag = "⛔ BLOCKS" if e >= CUTOFF else "⚠ named, pre-cutoff"
        if e >= CUTOFF:
            blocking.append(h)
        print(f"    {tag} #{e}  verdict says `{h['verdict'][:34]}…` but:")
        for name, why in h["found"]:
            print(f"        · `{name}` — {why}")
    if unreadable:
        print(f"    ⚠ {len(unreadable)} artifact(s) carry no `gate_verdict`/`gates` and are "
              f"UNREADABLE, NOT cleared")
    print(f"\n  ⚠ Scope: a match proves a contradiction is present; NO match is not proof of "
          f"consistency — a round can contradict itself in words this pattern does not know "
          f"(P6 safe side). From #{CUTOFF} on, a contradiction BLOCKS. Required: {REQUIRED}")
    return 1 if blocking else 0


if __name__ == "__main__":
    sys.exit(main())
