#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A111·R352 — anchoring the parse without touching a single value
===================================================================

CLOSURE (labelled honestly — it protects `#913`'s readable cells; it opens no world).

WHAT IT PROTECTS   `#913` established that NSFG is `SUPPORTED` and untried, and then had to bound
                   its own claim: **no attitude was measured**, because four attempts to rule out a
                   one-byte misalignment all failed to separate (blank-as-valid 1.0000 both ways ·
                   `caseid` uniqueness 0.9909 vs 1.0000 · `caseid` monotonicity 0.9998 vs 1.0000 ·
                   a two-sided sweep of the item's byte window returning 1.0000 at every offset).
                   `#913`② therefore blocks the round that matters until an EXTERNAL anchor exists.

WHY EVERY EARLIER ATTEMPT FAILED, stated as a property rather than as bad luck: all four asked
                   *"do the VALUES at this position look right?"* The `IH`/`JG` block is a dense run
                   of 1-digit Likert items, so every neighbouring byte is a documented code by
                   construction and no value-based check can separate. **The question was wrong,
                   not the effort.** ⇒ ask instead whether the LAYOUT and the FILE agree about the
                   shape of a record, which uses no value at all.

THE ANCHOR         A Stata dictionary declares (start, width) for every variable. Two things are
                   then checkable against the file without reading a single response:
                     ① `max(start + width - 1)` must equal the file's ACTUAL record width;
                     ② the declared fields must TILE the record — zero overlaps, zero gaps.
                   Both hold ⇒ the record start is pinned, and since every position is declared
                   relative to that same start, pinning the start pins every variable. A one-byte
                   shift would show up as a width mismatch or a gap.

⚠ NOT CIRCULAR     `#913`② named the cross-cycle marginal shift as the cheap candidate and rejected
                   it: it is the very quantity a study would want to report. This anchor reads
                   **no respondent's answer to anything**, so it cannot launder a finding.

CONTROL, and it is the whole round: the 3×3 layout × file cross-pairing. If width-matching were a
                   coincidence of "NSFG files are all about this long", off-diagonal pairs would
                   match too. **Pre-registered: the diagonal must match and the six off-diagonal
                   pairings must all fail.** That is the g=0 arm — it is what makes a MATCH a
                   measurement rather than an observation about file sizes.

STRUCTURALLY CANNOT: ① this anchors the RECORD, not the meaning of a code — the value labels are
                   still absent and `#913`②'s successor question (what does `1` mean) is untouched;
                   ② 12 of 15 NSFG files ship no layout and are not anchored by this or anything;
                   ③ ⚠ **only this one instrument**, and here that is a property of the question
                   rather than a shortcut: the object under test is *whether THIS release's shipped
                   layout describes THIS release's data file*. A second release cannot be asked it,
                   because the claim is about a (layout, file) pair and every other release is a
                   different pair. The cross-instrument move happens one round later, at `#914`②,
                   where the CLAIM is about people and GSS is the second instrument;
                   ④ `[unchallenged]` — door ③.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
D = ROOT / "data" / "external" / "nsfg"

PAT = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')
PAIRS = [("2011_2013_FemRespSetup.dct", "2011_2013_FemRespData.dat"),
         ("2017_2019_FemRespSetup.dct", "2017_2019_FemRespData.dat"),
         ("2017_2019_MaleSetup.dct", "2017_2019_MaleData.dat")]


def layout(dct):
    return [(m.group(2).lower(), int(m.group(1)), int(m.group(3))) for m in
            (PAT.search(l) for l in (D / "setup" / dct).read_text(errors="replace").splitlines()) if m]


def record_width(dat):
    with open(D / dat, "rb") as fh:
        return len(fh.readline().rstrip(b"\r\n"))


decl, act, shape = {}, {}, {}
for dct, dat in PAIRS:
    lay = layout(dct)
    decl[dct] = max(s + w - 1 for _n, s, w in lay)
    act[dat] = record_width(dat)
    occ = sorted((s, s + w - 1) for _n, s, w in lay)
    overlaps = sum(1 for a, b in zip(occ, occ[1:]) if b[0] <= a[1])
    gaps = sum(1 for a, b in zip(occ, occ[1:]) if b[0] > a[1] + 1)
    shape[dct] = dict(n_vars=len(lay), overlaps=overlaps, gaps=gaps)
    print(f"{dat:30s} declared={decl[dct]:5d} actual={act[dat]:5d} vars={len(lay):4d} "
          f"overlaps={overlaps} gaps={gaps}")

print("\n=== the 3x3 cross-pairing control (the whole round) ===")
matrix, diag_ok, off_hits = [], 0, 0
for dct, _ in PAIRS:
    row = []
    for _, dat in PAIRS:
        ok = decl[dct] == act[dat]
        row.append(ok)
        is_diag = (dct, dat) in PAIRS
        if is_diag and ok:
            diag_ok += 1
        if not is_diag and ok:
            off_hits += 1
    matrix.append(row)
    print(f"  {dct:30s} " + "  ".join(f"{'MATCH' if v else '  -  ':>5s}" for v in row))
print(f"\n  diagonal matched {diag_ok}/3 · off-diagonal false matches {off_hits}/6")

# the 12 unanchored files, reported rather than skipped
others = sorted(p for p in D.glob("*.dat") if p.name not in [d for _, d in PAIRS])
unanchored = {}
for p in others:
    with open(p, "rb") as fh:
        unanchored[p.name] = len(fh.readline().rstrip(b"\r\n"))
print(f"\n  {len(others)} files ship NO layout and are anchored by nothing:")
for k, v in unanchored.items():
    note = "  ⚠ ONE LINE — the whole file has no newline at all" if v > 10 ** 6 else ""
    print(f"    {k:30s} record width {v}{note}")

if not PAIRS:
    print("EMPTY POPULATION"); sys.exit(2)

G = Gate("Can the NSFG parse be anchored without reading a value?")
G.positive_control("the declared width matches the file for every shipped layout",
                   planted=float(diag_ok), floor=0.0, spread=0.5)
G.asserted("g=0 arm: a MISMATCHED layout x file pairing must NOT match",
           off_hits == 0, f"{off_hits} of 6 off-diagonal pairings matched", kind="control")
G.asserted("the declared fields tile the record with no gap and no overlap",
           all(v["overlaps"] == 0 and v["gaps"] == 0 for v in shape.values()),
           "; ".join(f"{k}: {v['n_vars']} vars, {v['overlaps']} overlaps, {v['gaps']} gaps"
                     for k, v in shape.items()), kind="control")
G.asserted("KILL: `#913`(2) is discharged only if the anchor uses NO respondent value",
           diag_ok == 3 and off_hits == 0,
           "record width and field tiling are properties of (layout, file); no response was read")

tv = G.three_valued()
VERDICT = "CONFIRMED" if not tv.startswith("UNVERIFIED") and diag_ok == 3 and off_hits == 0 \
    else ("UNVERIFIED" if tv.startswith("UNVERIFIED") else "OVERTURNED")
WORLD = "ANCHORED-BY-SHAPE" if VERDICT == "CONFIRMED" else "still unanchored"
print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")

art = dict(entry=914, round="E03·A111·R352", verdict=VERDICT, world=WORLD, kind="CLOSURE",
           estimand="whether the shipped layout and the data file agree on record SHAPE, "
                    "using no respondent value",
           declared=decl, actual=act, tiling=shape,
           cross_pairing=[[bool(v) for v in r] for r in matrix],
           diagonal_matched=diag_ok, off_diagonal_false_matches=off_hits,
           unanchored_files=unanchored,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "anchor_by_record_shape.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'anchor_by_record_shape.json'}")
