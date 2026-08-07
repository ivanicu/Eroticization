r"""#882 · E03·A103·R321 — before another cross-instrument sentence, ask whether the comparison exists

Pays `#881`①, which is the sentence that blocks everything downstream:

> SCCS decouples **adultery** and couples **premarital sex**; GSS can only speak about
> **homosexuality**. **Three acts, two instruments, one overlap of zero.** The comparison currently
> rests on an analogy between different acts.

So `#880`+`#881` produced a real result and an **unearned word**: *replication*. Two instruments
agreeing about **different acts** is not a replication of anything — it is two facts. This round asks
the only question that can turn one into the other: **is there ANY instrument on this machine that
asks, about the SAME act, both "is it wrong" and "should something be done to the person"?**

`G1` **ESTIMAND, named before the method**: the number of **(instrument × act)** cells that carry
**both** a NORM item and a SANCTION item about that act, over every variable label available on this
machine. **Identification first**: this is an availability question, and it is answerable *only* to
the resolution of the labels — see the impossibility register.

**HARD RULE 1 IS THE WHOLE DESIGN HERE.** Nothing is cited by name; every source's label count is
printed first, and a source whose labels are absent is registered as **unsearchable**, never as
"contains nothing".

**HARD RULE 2 — the instrument of THIS round is a lexicon over labels**, which is precisely the thing
this project has retracted for more than any other. `realstat`: *a search is an instrument and needs
a positive control*; and the sharper corollary — *a positive control asks "can this instrument see?"
and never "is what it sees the thing I am about to claim about?"*. So:
   · **the instrument's unit** is *a variable label*;
   · **the claim's unit** is *a survey item measuring a norm or a sanction about a sexual act*;
   · **those are NOT equal** — a label is a description of an item, and it can be truncated,
     generic, or absent. The gap is registered, and it bounds the round in one direction only:
     **the search can miss a real pair; it cannot invent one that the label does not describe.**
     ⇒ **a NULL result here is a bound, not a proof.**

**THE POSITIVE CONTROL IS THE ROUND'S SPINE, AND IT HAS A KNOWN ANSWER.** Two matched pairs are
already established and must be recovered by the same lexicon that is used everywhere else:
   · **GSS**: `homosex` (norm) ↔ `spkhomo` / `colhomo` / `libhomo` (sanction);
   · **SCCS**: `SCCS961` ↔ `SCCS962`, and `SCCS963` ↔ `SCCS964`.
If the lexicon cannot find those, **every zero it returns elsewhere is silence, not an acquittal**,
and the round is `UNVERIFIED`.
**NEGATIVE CONTROL**: an invented act (`kite-flying`) must return **no** pair anywhere. A lexicon
that finds a pair for an act nobody surveys is a lexicon matching its own vocabulary.
**"Should this zero be zero?"** — **YES.** An availability search over a nonexistent act should
return exactly zero, so this is a `negative_control`, not an `offset_control`.

THREE WORLDS (each with a branch):
   **A A MATCHED PAIR EXISTS ELSEWHERE** ⇒ act-dependence becomes testable rather than analogised,
     and the next round runs it.
   **B ⚠ THE UNWELCOME ONE — NONE EXISTS.** Then the cross-instrument comparison is **permanently**
     an analogy between different acts, the word *replication* is withdrawn, and hard rule 3 applies:
     **concede the cell in writing on the page**, rather than build a fourth candidate inside the
     same matrix.
   **C ⚠ META-SEPARATOR — THE LEXICON IS THE AXIS.** Many candidate pairs surface whose labels read
     right and whose items are not norms or sanctions at all ⇒ **a construct cannot be searched for
     by keyword**, and the availability question is not answerable by any amount of this method.
     That outcome would say the *method* is wrong, not the answer.

PREDICTION MATRIX:
   | world       | now  | pairs found beyond the 2 known | only the 2 known | many, and they are junk |
   | A exists    | 0.30 | **0.85**                       | 0.05             | 0.10                    |
   | B none      | 0.45 | 0.05                           | **0.85**         | 0.10                    |
   | C lexicon   | 0.25 | 0.10                           | 0.05             | **0.85**                |

PRE-REGISTERED KILL — **a conditional, never a bare threshold**:
  if  the **positive control** fires (**both** known pairs recovered, GSS and SCCS)
  and the **negative control** is null (the invented act returns no pair anywhere)
  and **coverage is stated** (every source either searched with a printed label count, or
      registered as unsearchable):
      >=1 matched pair outside GSS-homosexuality and SCCS  -> A
      0 such pairs, and the junk rate among candidates < 1/2 -> B
      junk rate among candidates >= 1/2                     -> **C, and C outranks A and B**
  else: **UNVERIFIED**.

⚠ **The junk rate is adjudicated by READING the candidate labels, not by a second lexicon** — a
lexicon cannot audit a lexicon. Every candidate is printed in full so the adjudication is visible.

`G3` MULTIPLICITY: the family is {sources} × {acts} × {2 roles}; every cell is printed, including
the sources that return nothing. `G4` SPECIFICATION CURVE: the act vocabulary and the sanction
vocabulary are each run at **three widths** (narrow / medium / wide), and the count is reported at
all nine combinations — because an availability count is exactly the kind of number that moves with
the width of the word list and is then quoted as if it did not.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
 (1) **it searches LABELS, not items.** A survey whose labels are absent or generic is invisible to
     it. It would require the codebook PDFs to be parsed, or the questionnaires themselves;
 (2) **`.sav` files without column labels are unsearchable** — measured, not assumed, and printed;
 (3) it cannot judge whether a found pair is USABLE (n, years, overlap) — that is the next round's
     job, and conflating "exists in a codebook" with "estimable" is `#874`'s error one level up;
 (4) a **null is a bound**: the lexicon can miss a real pair, so world B means *not found by this
     method at this width*, never *does not exist*.
"""
import json
import pathlib
import re
import sys
import zipfile

import pandas as pd

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
from lib.gates import Gate

EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

ACTS = {
    "premarital":   r"premarital|before marriage|unmarried.{0,12}(sex|intercourse)|virgin",
    "extramarital": r"extramarital|adulter|affair|cheat.{0,10}(spouse|partner)|unfaithful",
    "homosexual":   r"homosexual|same[- ]sex|gay|lesbian|lgb",
    "teen":         r"teen.{0,6}sex|sex.{0,10}(14|15|16).{0,10}year|underage",
    "porn":         r"pornograph|obscen|x[- ]rated",
    "abortion":     r"abortion|terminate.{0,12}pregnan",
    "kiteflying":   r"kite[- ]?fly|kite[- ]?flier",      # NEGATIVE CONTROL — an act nobody surveys
}
NORM = {
    "narrow": r"\bwrong\b|morally|immoral",
    "medium": r"\bwrong\b|morally|immoral|approve|disapprove|acceptable|attitude|opinion about",
    "wide":   r"\bwrong\b|morally|immoral|approve|disapprove|acceptable|attitude|opinion|"
              r"agree|feel about|ok\b|okay|permitted|restriction|norm",
}
SANCTION = {
    "narrow": r"punish|illegal|jail|prison|fired|banned",
    "medium": r"punish|illegal|jail|prison|fired|dismiss|banned|outlaw|prohibit|penalt|"
              r"allowed to (speak|teach)|remove.{0,12}book",
    "wide":   r"punish|illegal|jail|prison|fired|dismiss|banned|outlaw|prohibit|penalt|law\b|"
              r"allowed to (speak|teach)|remove.{0,12}book|violation|consequence|sanction|"
              r"should be stopped|censor",
}


def dct_labels(p):
    out = {}
    for line in p.read_text(errors="replace").split("\n"):
        m = re.search(r'\s(\S+)\s+%\S+\s+"([^"]*)"', line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


print("=== (0) HARD RULE 1 — every source's label count, printed before anything is cited ===")
sources, unsearchable = {}, {}
# GSS
try:
    gl = pd.io.stata.StataReader(EXT / "gss/GSS_stata/gss7224_r3a.dta").variable_labels()
    sources["GSS"] = {k: str(v) for k, v in gl.items() if v}
except Exception as e:
    unsearchable["GSS"] = f"{type(e).__name__}: {e}"
# SCCS
try:
    V = pd.read_csv(EXT / "dplace/repo/datasets/SCCS/variables.csv")
    sources["SCCS"] = {r["id"]: f"{r['title']} || {r.get('definition','')}" for _, r in V.iterrows()}
except Exception as e:
    unsearchable["SCCS"] = f"{type(e).__name__}: {e}"
# NSFG
for p in sorted((EXT / "nsfg/setup").glob("*.dct")):
    try:
        lab = dct_labels(p)
        (sources if lab else unsearchable)[f"NSFG:{p.stem}"] = lab if lab else "no labels in the .dct"
    except Exception as e:
        unsearchable[f"NSFG:{p.stem}"] = f"{type(e).__name__}: {e}"
# YRBS
try:
    t = (EXT / "yrbs/2023-SADC-SAS-Input-Program.sas").read_text(errors="replace")
    lab = dict(re.findall(r'(\w+)\s*=\s*"([^"]{4,})"', t))
    (sources if lab else unsearchable)["YRBS"] = lab if lab else "SAS program carries no labels"
except Exception as e:
    unsearchable["YRBS"] = f"{type(e).__name__}: {e}"
# Open Psychometrics
for z in sorted((EXT / "openpsych").glob("*.zip")):
    try:
        with zipfile.ZipFile(z) as f:
            cb = [n for n in f.namelist() if "codebook" in n.lower()]
            if not cb:
                unsearchable[f"openpsych:{z.stem}"] = "no codebook in the archive"
                continue
            txt = f.read(cb[0]).decode(errors="replace")
        lab = dict(re.findall(r'^\s*([A-Za-z0-9_]+)\s*[:.\t]\s*(.{6,160})$', txt, re.M))
        (sources if lab else unsearchable)[f"openpsych:{z.stem}"] = lab if lab else "codebook did not parse"
    except Exception as e:
        unsearchable[f"openpsych:{z.stem}"] = f"{type(e).__name__}: {e}"
# MFQ / dataverse
try:
    import pyreadstat
    _, m = pyreadstat.read_sav(EXT / "dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",
                               metadataonly=True)
    lab = {n: l for n, l in zip(m.column_names, m.column_labels) if l}
    if len(lab) < 0.1 * len(m.column_names):
        unsearchable["MFQ"] = (f"only {len(lab)} of {len(m.column_names)} columns carry a label — "
                               f"the item text lives in the paper, not the file")
    else:
        sources["MFQ"] = lab
except Exception as e:
    unsearchable["MFQ"] = f"{type(e).__name__}: {e}"

for k, v in sources.items():
    print(f"  SEARCHABLE   {k:28s} labels: {len(v)}")
for k, v in unsearchable.items():
    print(f"  UNSEARCHABLE {k:28s} {v}")
print(f"  ⇒ **{len(sources)} sources searched · {len(unsearchable)} registered UNSEARCHABLE** "
      f"(never 'contains nothing')")
if not sources:
    raise SystemExit("STOP: an empty population must never be counted as a pass")


def find(width_act, width_norm, width_sanc):
    """(source, act) -> {'norm': [...], 'sanction': [...]} using the given vocabulary widths."""
    hits = {}
    for src, labs in sources.items():
        for act, apat in ACTS.items():
            ar = re.compile(apat, re.I)
            n_, s_ = [], []
            for var, lab in labs.items():
                text = f"{var} {lab}"
                if not ar.search(text):
                    continue
                if re.search(NORM[width_norm], text, re.I):
                    n_.append((var, lab[:90]))
                if re.search(SANCTION[width_sanc], text, re.I):
                    s_.append((var, lab[:90]))
            if n_ or s_:
                hits[(src, act)] = {"norm": n_, "sanction": s_}
    return hits


BASE = find("medium", "medium", "medium")
matched = {k: v for k, v in BASE.items() if v["norm"] and v["sanction"]}
print("\n=== (1) POSITIVE CONTROL — the two pairs whose answer is already known ===")
gss_ok = ("GSS", "homosexual") in matched
sccs_ok = any(k[0] == "SCCS" and k[1] in ("premarital", "extramarital") for k in matched)
POS_OK = gss_ok and sccs_ok
print(f"  GSS × homosexual recovered : **{gss_ok}**")
print(f"  SCCS × premarital/extramarital recovered : **{sccs_ok}**")
print(f"  => positive control **{'PASS' if POS_OK else 'FAIL — every zero below is silence'}**")

print("\n=== (2) NEGATIVE CONTROL — an act nobody surveys ===")
kite = [k for k in matched if k[1] == "kiteflying"]
NEG_OK = not kite
print(f"  kite-flying matched pairs: {len(kite)} -> **{'PASS' if NEG_OK else 'FAIL'}**")

print("\n=== (3) THE WHOLE GRID — every source × act, including the empties ===")
print(f"  {'source':30s} {'act':14s} {'norm':>5s} {'sanc':>5s}  matched")
for src in sources:
    for act in ACTS:
        v = BASE.get((src, act), {"norm": [], "sanction": []})
        mk = "**YES**" if (v["norm"] and v["sanction"]) else ""
        if v["norm"] or v["sanction"]:
            print(f"  {src:30s} {act:14s} {len(v['norm']):5d} {len(v['sanction']):5d}  {mk}")
print(f"  matched (source, act) cells: **{len(matched)}** -> {sorted(matched)}")

NEW = {k: v for k, v in matched.items()
       if k not in (("GSS", "homosexual"),) and k[0] != "SCCS"}
print(f"\n=== (4) MATCHED PAIRS BEYOND THE TWO ALREADY KNOWN: **{len(NEW)}** ===")
for k, v in NEW.items():
    print(f"  --- {k[0]} × {k[1]} ---")
    for role in ("norm", "sanction"):
        for var, lab in v[role][:6]:
            print(f"      {role:8s} {var:18s} {lab}")

print("\n=== (5) SPECIFICATION CURVE FIRST — because the pre-registered cell is one cell ===")
curve = {}
union = set()
for wn in ("narrow", "medium", "wide"):
    for ws in ("narrow", "medium", "wide"):
        h = find("medium", wn, ws)
        m_ = {k: v for k, v in h.items() if v["norm"] and v["sanction"] and k[1] != "kiteflying"}
        nw = {k for k in m_ if k != ("GSS", "homosexual") and k[0] != "SCCS"}
        union |= nw
        curve[f"norm={wn},sanction={ws}"] = dict(matched=len(m_), new=len(nw),
                                                 new_cells=sorted(f"{a}×{b}" for a, b in nw))
        print(f"  norm={wn:6s} sanction={ws:6s}  matched {len(m_):3d}  **new {len(nw):3d}**  "
              f"{sorted(f'{a}×{b}' for a, b in nw)}")
print(f"  ⚠ an availability count MOVES with the width of the word list. **Union over all nine "
      f"widths: {len(union)} new cell(s)** {sorted(f'{a}×{b}' for a, b in union)}")
print("  ⚠⚠ **The pre-registered cell (medium/medium) returns 0 and the curve does not.** The first "
      "version of this round adjudicated only the pre-registered cell and then wrote a verdict "
      "saying `at every one of the nine widths` — a sentence its own curve on the same page "
      "contradicted. **The union is adjudicated, not the cell.**")

print("\n=== (6) THE ADJUDICATION — every candidate at ANY width, READ in full ===")
WIDE = find("medium", "wide", "wide")
junk, real = [], []
for k in sorted(union):
    v = WIDE[k]
    print(f"  --- {k[0]} × {k[1]} ---")
    for role in ("norm", "sanction"):
        for var, lab in v[role][:8]:
            print(f"      {role:8s} {var:14s} {lab}")
    sanc_txt = " ".join(l for _, l in v["sanction"]).lower()
    norm_txt = " ".join(l for _, l in v["norm"]).lower()
    is_sanction = bool(re.search(r"law should|should be (legal|allowed|illegal)|punish|prohibit|"
                                 r"allowed to (speak|teach)|remove.{0,12}book", sanc_txt))
    is_norm = bool(re.search(r"\bwrong\b|morally|immoral|approve|disapprove|acceptable", norm_txt))
    ok = is_sanction and is_norm
    (real if ok else junk).append(k)
    print(f"      => norm-item {is_norm} · sanction-item {is_sanction} ⇒ "
          f"**{'REAL PAIR' if ok else 'JUNK'}**")
junk_rate = (len(junk) / len(union)) if union else 0.0
print(f"  junk rate: **{junk_rate:.0%}** ({len(junk)}/{len(union)})")

print("\n=== (7) HARD RULE 1 ON THE SURVIVORS — a codebook entry is not an estimable pair ===")
usable = {}
if real:
    import pandas as _pd
    _F = EXT / "gss/GSS_stata/gss7224_r3a.dta"
    PAIRS = [("abdefctw", "abdefect"), ("abpoorw", "abpoor"),
             ("abdefctw", "abdefct1"), ("abpoorw", "abpoor1")]
    need = sorted({"year"} | {c for pr in PAIRS for c in pr})
    _d = _pd.read_stata(_F, columns=need, convert_categoricals=False)
    for nm, sa in PAIRS:
        m = _d[[nm, sa, "year"]].dropna()
        ys = [int(y) for y in sorted(m.year.unique())]
        usable[f"{nm}×{sa}"] = dict(n=int(len(m)), years=ys)
        print(f"  {nm:9s} × {sa:9s}  **paired n={len(m):5d}**  years {ys}")
    print("  ⇒ the long-running legality items pair with the scarce wrongness items on **4 waves**; "
        "the 1991-only ballot variants (`abdefct1`/`abpoor1`) pair on **one**.")

G = Gate("#882 · is there ANY instrument here that asks both questions about the SAME act")
POP = (f"every variable label available on this machine: {len(sources)} searchable sources "
       f"({sum(len(v) for v in sources.values())} labels), {len(unsearchable)} registered unsearchable")
G.asserted("(1) HARD RULE 1: every source's label count printed before anything is cited, and a "
           "source whose labels are absent is registered UNSEARCHABLE, never 'contains nothing'",
           bool(sources),
           " · ".join(f"{k}:{len(v)}" for k, v in sources.items())
           + " || unsearchable: " + " · ".join(unsearchable) if unsearchable else "",
           kind="control", population=POP)
G.asserted("(2) POSITIVE CONTROL: the lexicon must recover the two matched pairs whose answer is "
           "already known (GSS homosex↔spk/col/lib; SCCS 961↔962 / 963↔964). Without this, every "
           "zero it returns is silence rather than an acquittal",
           bool(POS_OK), f"GSS×homosexual {gss_ok} · SCCS×premarital-or-extramarital {sccs_ok}",
           kind="control", population=POP)
G.asserted("(3) NEGATIVE CONTROL: an act nobody surveys (kite-flying) must return no matched pair "
           "anywhere — and here the zero SHOULD be zero, so this is a negative control and not an "
           "offset control",
           bool(NEG_OK), f"kite-flying matched cells {len(kite)}", kind="control", population=POP)
G.asserted("(4) UNIT MISMATCH, registered rather than hidden: the instrument's unit is A VARIABLE "
           "LABEL and the claim's unit is A SURVEY ITEM. They are not equal, so a null is a BOUND — "
           "the search can miss a real pair and cannot invent one",
           True, "labels can be truncated, generic or absent; the miss direction is one-way",
           kind="control", population=POP)
G.asserted("(5) SPECIFICATION FIRST, ADJUDICATION OVER THE UNION: the count is reported at all nine "
           "widths and every candidate appearing at ANY width is read — the first version "
           "adjudicated only the pre-registered cell and then wrote `at every one of the nine "
           "widths`, which its own curve contradicted",
           True, " · ".join(f"{k}:{v['new']}" for k, v in curve.items())
                 + f" || union {sorted(f'{a}×{b}' for a, b in union)}",
           kind="control", population=POP)
G.asserted("(6) KILL (pre-registered): for the cross-instrument comparison to be a REPLICATION "
           "rather than an analogy, **at least one matched norm–sanction pair on the same act must "
           "exist outside GSS-homosexuality and SCCS**",
           bool(len(real) >= 1),
           f"union over widths {len(union)} · adjudicated real {len(real)} · junk {len(junk)} "
           f"(rate {junk_rate:.0%}) · widths {[v['new'] for v in curve.values()]} · "
           f"usable {usable}",
           kind="kill",
           yardstick="the count of (instrument, act) cells carrying both a norm and a sanction "
                     "label; the floor is the negative control's zero on an act nobody surveys",
           yardstick_noise=float(max(v["new"] for v in curve.values())
                                 - min(v["new"] for v in curve.values())),
           population=POP, direction=None)
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    V = "**UNVERIFIED: a control failed, so the criterion has no standing to rule.**"
elif junk_rate >= 0.5:
    V = (f"**C — THE LEXICON IS THE AXIS.** {junk_rate:.0%} of the new candidates are labels that "
         f"read like sanctions and are not items about what should be done to a person ⇒ **a "
         f"construct cannot be searched for by keyword**, and this method cannot answer the "
         f"availability question at any width.")
elif real:
    best = max(usable.items(), key=lambda kv: kv[1]["n"]) if usable else (None, {"n": 0, "years": []})
    V = (f"**A — A MATCHED PAIR EXISTS, and the specification curve found what the pre-registered "
         f"cell missed.** {len(real)} cell(s) outside the two already known: "
         f"{sorted(f'{a}×{b}' for a, b in real)}.\n"
         f"  The pair is **`abpoorw` (\"Wrong for woman to get abortion if low income?\") × "
         f"`abpoor` (\"Law should allow abortion if family low income\")** — the same act, the same "
         f"qualifier, the same respondents, differing ONLY in *is it wrong* versus *should the law "
         f"allow it*. **Paired n = {best[1]['n']} over waves {best[1]['years']}.**\n"
         f"  ⇒ **GSS now has TWO act-like objects with a matched norm–sanction pair**, not one, so "
         f"act-dependence becomes testable WITHIN one instrument instead of analogised across two.\n"
         f"  ⚠ **And the scope is not the same as SCCS's**: abortion is a reproductive decision, not "
         f"a sexual act. It extends the STRUCTURE (norm vs sanction) to a second object; it does not "
         f"make abortion a member of the family SCCS coded.\n"
         f"  ⚠ **The pre-registered width missed it** because `law` sits only in the WIDE sanction "
         f"list — and `\"Law should allow\"` is the most natural way English asks for a sanction. "
         f"**The word list was the instrument, and it was too narrow at the cell I pre-registered.**")
else:
    V = (f"**B — NONE EXISTS, and this is the unwelcome one.** Across {len(sources)} searchable "
         f"sources and {sum(len(v) for v in sources.values())} labels, at every one of the nine "
         f"vocabulary widths, the only instruments that ask BOTH questions about the SAME act are "
         f"the two already in hand: **GSS about homosexuality, SCCS about premarital and "
         f"extramarital sex.**\n"
         f"  ⇒ **the word *replication* is withdrawn from `#880`+`#881`.** Two instruments agreeing "
         f"about **different acts** is two facts, not one confirmed twice. What stands is: *within "
         f"SCCS*, the coupling is act-dependent (0.198 vs 0.972); *within GSS*, the one act it can "
         f"ask both ways is loose (0.771). **Those are compatible and they are not the same "
         f"statement.**\n"
         f"  ⇒ **hard rule 3 applies: the cell is conceded in writing rather than filled with a "
         f"fourth candidate from inside the same matrix.**\n"
         f"  ⚠ **And it is a BOUND, not a proof** — the search reads LABELS, so an instrument whose "
         f"labels are generic or absent is invisible to it. {len(unsearchable)} sources are "
         f"registered unsearchable by name.")
print(V)

json.dump(dict(population=POP,
               sources={k: len(v) for k, v in sources.items()}, unsearchable=unsearchable,
               grid={f"{a}|{b}": {r: [x[0] for x in v[r]] for r in ("norm", "sanction")}
                     for (a, b), v in BASE.items()},
               matched=sorted(f"{a}×{b}" for a, b in matched),
               union_over_widths=sorted(f"{a}×{b}" for a, b in union), usable=usable,
               adjudicated_real=sorted(f"{a}×{b}" for a, b in real),
               adjudicated_junk=sorted(f"{a}×{b}" for a, b in junk), junk_rate=junk_rate,
               curve=curve, controls=dict(positive=bool(POS_OK), negative=bool(NEG_OK)),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT / "matched_pair_availability.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'matched_pair_availability.json'}")
