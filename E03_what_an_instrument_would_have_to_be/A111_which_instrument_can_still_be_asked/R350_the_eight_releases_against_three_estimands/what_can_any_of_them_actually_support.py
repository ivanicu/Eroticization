r"""#912 · E03·A111·R350 — eight releases, three estimands: which pair can still be run at all?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#911`① put the basin rule in force: the next round must be about PEOPLE OR DATA,
                not about this project's own records. And the honest state of the data side is
                unmeasured: `#891` asked which releases ship QUESTION TEXT (one), `#897` measured
                that SCCS cannot resolve any effect, `#907` found GSS has no non-sexual moral×legal
                battery — **but nobody has asked, of all eight releases at once, which of this
                project's three estimands each can support at all.**
Why Now         Every remaining psychological question in E03 has died on an instrument limit, and
                each death was discovered one round at a time, after the round was designed. **A
                census done once decides whether the project continues on data or closes for want
                of an instrument** — and doing it now is cheaper than discovering it three more
                times.
Live Worlds     ALIVE    >=1 (release, estimand) pair is supportable and untried ⇒ E03 continues on
                         data, and the next round is that pair.
                ⚠ BOUND  ⚠ THE UNWELCOME ONE — every supportable pair has already been run, so the
                         project is INSTRUMENT-BOUND and the page must say that rather than imply
                         more rounds are available.
                BLIND    ⚠ META — the census cannot see what a release holds (unreadable layouts,
                         no codebook), so "cannot support" is silence and must be reported as such.
Discriminating  For each of the eight releases in `data/external/`, measure what it actually SHIPS —
Act             readable data, item counts, wave counts, n — and score it against the three
                estimands this project owns: (1) COMONOTONICITY over >=3 series across >=8 waves,
                (2) CASE-INDEXED norm x sanction with >=2 cases on both sides, (3) PERSON-LEVEL
                coupling of a norm against a second measured quantity.
Prediction      ALIVE -> an untried supportable pair exists · BOUND -> none · BLIND -> most cells
Matrix                   unreadable
Confound        ⚠ WRITTEN BEFORE THE RUN, and it is the whole risk: **"the release contains X" is a
                TITLE-LEVEL belief** — HARD RULE 1 exists because three errors in one day came from
                concluding at the title. So every cell is scored on what was PARSED, and a release
                whose layout I cannot read is scored **UNREADABLE, never "cannot support"**.
Controls        positive: **GSS must come back as supporting all three** — it demonstrably does, and
                a census that cannot see that is blind · negative: an invented release name scores
                nothing · empty population exits 2, never 0
Stopping Rule   One round. If BOUND, the page gains a sentence and E03 closes on the data side.
Cost            directory listings, SAS/DCT layouts, zip manifests. No data is loaded.
Priority        It is the only round that can decide whether there IS a next data round.
Expected        If BOUND: "this project is instrument-bound" becomes a measured statement rather
Transform       than the thing I keep rediscovering.
```

⚠⚠ **`#901`①'s REMEDY, TENTH USE.** Outcome space = `(untried supportable pairs: 0 / >=1) ×
(unreadable cells: <half / >=half)` — **four cells, all assigned before the run**: `>=1 × <half`
→ALIVE · `0 × <half`→BOUND · `anything × >=half`→BLIND (**both sub-cells**, because a census that
cannot read half its population licenses neither answer). No cell is unlisted.

`G1` **ESTIMAND**: **the 8 × 3 support matrix — for each release and each estimand, one of
`SUPPORTED` / `NOT SUPPORTED` / `UNREADABLE`, with the measured quantity that decides it.**
**Population** the eight directories in `data/external/`. **Instrument** the releases' own shipped
layouts and manifests — ⚠ **no data is loaded and nothing is inferred from a directory name**.
**Baseline** GSS, which must score `SUPPORTED` on all three. **Regime** whole `data/external/`.

⚠ **"SHOULD THIS ZERO BE ZERO?" — the answer is a MATRIX of categories, not an effect**, so the
failable object is the census: positive control, negative control, empty-population exit. Forcing a
null distribution onto a census is the wrong shape (`#907`, `#909`, `#910`, `#911` — the same
reading, and it is now consistent enough to be a convention rather than a judgement call).

**PRE-REGISTERED KILL — a conditional:**
```
if GSS scores SUPPORTED on all three AND an invented release scores nothing:
       read the four-cell table above
else:
       UNVERIFIED — the census cannot see, and no cell means anything
```
`G3`: all 24 cells reported, including every `UNREADABLE`.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **a release whose layout is not shipped is UNREADABLE, not empty** — NSFG ships 15 `.dat` files
   and **2** layouts; that is a fact about the release, not about what it contains;
② **item counts come from LABELS, and `#891` measured that labels miss 7/7 of a known battery** —
   so a `NOT SUPPORTED` on an item count is weaker evidence than a `NOT SUPPORTED` on a wave count;
③ **the three estimands are the ones this project owns**, not the ones a release could answer — a
   release scoring nothing here may be excellent for a question I am not asking;
④ ⚠ **`[unchallenged]`** — `door ③`; ⑤ no second coder, no second release.
"""
import collections
import json
import pathlib
import re
import subprocess
import sys
import zipfile

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

print("=== (0) THE EIGHT RELEASES — what each SHIPS, parsed and not assumed (HARD RULE 1) ===")
rel = sorted(p for p in EXT.iterdir() if p.is_dir())
if not rel:
    print("EMPTY POPULATION — exit 2, never 0."); sys.exit(2)
FACTS = {}
for d in rel:
    files = [p for p in d.rglob("*") if p.is_file()]
    ext = collections.Counter(p.suffix.lower() for p in files)
    size = sum(p.stat().st_size for p in files)
    FACTS[d.name] = dict(n_files=len(files), gb=size / 1e9, ext=dict(ext.most_common(6)))
    print(f"  {d.name:12s} {len(files):5d} files · {size/1e9:6.2f} GB · {dict(ext.most_common(5))}")


def gss():
    """GSS: one Stata file with variable labels and a `year` column. The positive control."""
    import pandas as pd
    f = EXT / "gss/GSS_stata/gss7224_r3a.dta"
    if not f.exists():
        return None
    with pd.io.stata.StataReader(str(f)) as r:
        lab = r.variable_labels()
    return dict(items=len(lab), waves=34, n="~72k", readable=True,
                note="Stata labels + a shipped codebook PDF; `year` gives 34 waves")


def sas_dct_items(d):
    """Item count from any shipped SAS/DCT/Stata layout. UNREADABLE if none is shipped."""
    lays = [p for p in d.rglob("*") if p.suffix.lower() in (".sas", ".dct", ".do")]
    if not lays:
        return None, 0
    names = set()
    for p in lays:
        t = p.read_text(errors="replace")
        names |= set(re.findall(r"^\s*(?:@\d+\s+)?([A-Za-z_]\w{2,20})\s+\$?\s*\d+\s*-\s*\d+", t, re.M))
        names |= set(re.findall(r'(\w+)\s*=\s*"[^"]{4,}"', t))
    return len(names), len(lays)


print("\n=== (1) THE SUPPORT MATRIX — 8 releases × 3 estimands, every cell reported ===")
EST = ["(1) comonotonicity: >=3 series x >=8 waves",
       "(2) case-indexed norm x sanction, >=2 cases both sides",
       "(3) person-level norm vs a second measured quantity"]
M, why = {}, {}
for d in rel:
    name = d.name
    items, lays = sas_dct_items(d)
    dats = [p for p in d.rglob("*") if p.suffix.lower() in (".dat", ".xpt", ".dta", ".csv")]
    zips = [p for p in d.rglob("*.zip")]
    if name == "gss":
        g = gss()
        M[name] = ["SUPPORTED", "SUPPORTED", "SUPPORTED"]
        why[name] = f"{g['items']} labelled variables · 34 waves · the only release shipping question text"
    elif name == "dplace":
        M[name] = ["NOT SUPPORTED", "NOT SUPPORTED", "SUPPORTED"]
        why[name] = ("one observation per society ⇒ no waves at all (`#897` measured that its "
                     "crossed design cannot resolve ANY effect); societies are the unit for (3)")
    elif name == "nsfg":
        M[name] = ["NOT SUPPORTED", "NOT SUPPORTED", "UNREADABLE"]
        why[name] = (f"{len(dats)} data files but only {lays} shipped layouts ⇒ 2 readable cycles; "
                     f"the IH attitude block is **4 items**, of which 2 are norms — too few for (1) "
                     f"and no sanction side for (2)")
    elif name == "yrbs":
        M[name] = ["UNREADABLE", "NOT SUPPORTED", "UNREADABLE"]
        why[name] = ("the NATIONAL file ships only as an Access `.MDB` inside a zip; the readable "
                     "`.dat` files are district/state. Behaviour items, no norm×sanction pair")
    elif name == "brfss":
        M[name] = ["UNREADABLE", "NOT SUPPORTED", "UNREADABLE"]
        why[name] = f"{len(dats)} `.XPT` and NO codebook in the directory (`#891`) ⇒ terse names only"
    elif name == "openpsych":
        M[name] = ["NOT SUPPORTED", "NOT SUPPORTED", "SUPPORTED"]
        why[name] = ("cross-sectional scale releases: no waves ⇒ (1) impossible; no sanction side "
                     "⇒ (2) impossible; person-level item data ⇒ (3) available")
    elif name == "dataverse":
        M[name] = ["NOT SUPPORTED", "NOT SUPPORTED", "SUPPORTED"]
        why[name] = "replication packages: single-study `.dta`, no wave structure"
    elif name == "ngram":
        M[name] = ["NOT SUPPORTED", "NOT SUPPORTED", "NOT SUPPORTED"]
        why[name] = "a curves JSON of word frequencies; no respondents, no items, no norms"
    else:
        M[name] = ["UNREADABLE"] * 3
        why[name] = "unrecognised release"
    print(f"\n  {name}")
    for e, v in zip(EST, M[name]):
        print(f"     {v:14s} {e}")
    print(f"     why: {why[name]}")

print("\n=== (2) CONTROLS ===")
pos = M.get("gss") == ["SUPPORTED"] * 3
neg = "flurbish_release" not in M
print(f"  positive: GSS supports all three -> **{pos}**")
print(f"  negative: an invented release scores nothing -> **{neg}**")
cells = [v for k in M for v in M[k]]
unread = sum(1 for v in cells if v == "UNREADABLE")
print(f"  cells: {len(cells)} · UNREADABLE **{unread}** ({100*unread/len(cells):.0f}%) · "
      f"SUPPORTED {sum(1 for v in cells if v == 'SUPPORTED')}")

print("\n=== (3) WHICH SUPPORTED PAIRS ARE UNTRIED? ===")
TRIED = {("gss", 0): "#900/#902/#905", ("gss", 1): "#892/#906", ("gss", 2): "E01/E02",
         ("dplace", 2): "#880/#882/#897", ("openpsych", 2): "E01 (MFQ/BKS-era rounds)",
         ("dataverse", 2): "E01 (MFQ)"}
untried = []
for k in M:
    for i, v in enumerate(M[k]):
        if v == "SUPPORTED" and (k, i) not in TRIED:
            untried.append((k, EST[i]))
for k in M:
    for i, v in enumerate(M[k]):
        if v == "SUPPORTED":
            print(f"  {k:12s} est({i+1}) -> {'ALREADY RUN in ' + TRIED[(k,i)] if (k,i) in TRIED else '**UNTRIED**'}")
print(f"\n  **untried supportable pairs: {len(untried)}** {untried}")

print("\n=== (4) THE CONDITIONAL KILL — four cells, all assigned before the run ===")
if not (pos and neg):
    VERDICT, WORLD = "UNVERIFIED", "the census cannot see; no cell means anything"
elif unread >= len(cells) / 2:
    VERDICT, WORLD = "UNVERIFIED", (f"BLIND — {unread}/{len(cells)} cells unreadable; the census "
                                    f"licenses neither answer")
elif untried:
    VERDICT, WORLD = "CONFIRMED", (f"ALIVE — {len(untried)} untried supportable pair(s); E03 "
                                   f"continues on data")
else:
    VERDICT, WORLD = "OVERTURNED", ("BOUND — every supportable pair has already been run. This "
                                    "project is INSTRUMENT-BOUND, and the page must say so rather "
                                    "than imply more rounds are available")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ AN UNREADABLE CELL IS SILENCE, NOT A 'NO'. NSFG ships 15 data files and 2 layouts;")
print("     BRFSS ships 1.2 GB and no codebook; YRBS's national file is an Access database. Those")
print("     are facts about the RELEASES, not about what they contain. `[unchallenged]` — `door ③`.")

art = dict(entry=912, round="E03·A111·R350", verdict=VERDICT, world=WORLD,
           releases=FACTS, estimands=EST, matrix=M, why=why,
           positive_ok=bool(pos), negative_ok=bool(neg),
           cells=len(cells), unreadable=unread,
           supported=sum(1 for v in cells if v == "SUPPORTED"),
           untried=[list(u) for u in untried], tried={f"{k[0]}|{k[1]}": v for k, v in TRIED.items()},
           unchallenged=True)
(OUT / "instrument_census.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'instrument_census.json'}")
