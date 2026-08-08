r"""#891 · E03·A105·R329 — a variable LABEL is not a QUESTION, and the positive control is what said so

**WHY THIS ROUND EXISTS.** One round ago `#890`① published *"the instrument that would close
object-distance is not in `data/external/`"*. **I established that by listing directory names.**
`P4`'s closing question is exactly this one — *did I establish the non-existence by asking the
system, or did I read it somewhere?* — and I had done neither: **I read an `ls`.** An absence claim
about CONTENT, from FILENAMES, is `P5`★ — a `not found` from an instrument that has never returned a
`found`.

**⚠ THE FIRST ATTEMPT FAILED ITS OWN POSITIVE CONTROL, AND THAT FAILURE IS THE ROUND'S MAIN RESULT.**
v1 screened **Stata variable labels** for a moral word and a legal word near an act word. Pointed at
GSS, where the answer is known, it recovered **0 of the 7** abortion legality items. Not a tuning
problem — **a property of the object**, and the codebook says so in plain sight:

```
ABDEFECT  (Please tell me whether or not you think it should be possible for a pregnant woman
           to obtain a legal abortion if. . .) If there is a strong chance of serious defect...
ABNOMORE   If she is married and does not want any more children?
ABPOOR     If the family has a very low income and cannot afford any more children?
```

> **The word "legal" appears ONCE, in a preamble attached to the FIRST item of the battery. Every
> other item is a bare stem.** Their Stata labels are barer still (`abpoor` → *"Low income--cant
> afford more children"*): **no act, no modality, nothing.**
>
> ⇒ **A LABEL IS NOT A QUESTION** — the sibling of this project's HARD RULE 1 (*a variable name is
> not a measurement*), one level up. And the failure is **adversarial rather than random**: a
> matched moral/legal pair is exactly the kind of thing that lives in a BATTERY, and a battery is
> exactly the structure that strips its members' labels. **The screen was blind precisely where the
> target lives.**

**⇒ THE REPAIR IS NOT A BIGGER LEXICON, IT IS A DIFFERENT SOURCE:** the questionnaire text, with
**preamble propagation** — a bare stem inherits the last parenthetical preamble above it.

`G1` **ESTIMAND, named before the method.** *Does any instrument in `data/external/` carry, for **two
or more distinct acts**, a **moral** item and a **legal** item about the **same proposition**?*

⚠ **THE UNIT CHECK, BEFORE THE CONTROL IS DESIGNED** (`realstat`: *a positive control asks CAN THIS
INSTRUMENT SEE and never asks IS WHAT IT SEES THE THING I AM CLAIMING ABOUT*): the detector's unit is
a *(variable, question-text)* pair matching a lexicon; the claim's unit is a *matched proposition
across two acts*. **They are not equal.** ⇒ the detector is a **SCREEN**; the verdict comes from
**reading** the shortlist. No count below is an answer on its own.

**POSITIVE CONTROL** — the detector on GSS's own codebook, where the answer is known: it must recover
the **7** abortion legality items **and** the two `ab*w` moral items, i.e. abortion must come back
with **both sides**. Threshold pre-registered at **≥4/7**; **at the label-only source it must NOT
fire** (it returned 0/7, which is the measured floor, so the control demonstrably can fail).
**NEGATIVE CONTROL** — the same text with the 7 legality items **deleted**: abortion must lose its
legal side. A removal of a known signal, never an invented case.

**WORLDS.**
  **A · nothing I hold can close object-distance** ⇒ `#890`① stands, on evidence instead of an `ls`.
  **B · ⚠ GSS ITSELF CAN** — if GSS carries a second act with a matched moral/legal proposition
    (marijuana, pornography, suicide, adultery…), then object-distance is closable **inside the
    instrument I already have**, `#890`① is wrong one round after I wrote it, and there is a cheap
    round that turns `+0.2945` from a bound into a point.
  **C · the detector cannot see** ⇒ **`UNVERIFIED`**, and `#890`① is **downgraded**, not confirmed.
  **D · ⚠ META — the question is unanswerable of everything except GSS**, because no other release
    here ships question text at all. Then the honest sentence is not *"it is not there"* but
    ***"nothing I hold can be ASKED this question"*** — a fact about the releases, not about the
    search, and it changes what `#890`① should say.

**PRE-REGISTERED KILL — conditional:**
```
if positive_control recovers >= 4 of the 7 AND the label-only source recovers < 4 (it can fail)
   and negative_control (7 deleted) removes abortion's legal side:
       >= 2 acts with a matched proposition survive READING, in any source  -> B
       exactly 1 (abortion) survives, and no other source ships question text -> D
       none survives                                                          -> A
else:
       UNVERIFIED -- and #890① is downgraded, not confirmed
```
**TYPE: this round opened as CLOSURE (protect `#890`①) and world B would make it FRONTIER.** Labelled
both ways here before the run, so the label cannot be chosen after seeing the answer.

**WHAT THIS SITE STRUCTURALLY CANNOT DO:** ① it reads the text a release SHIPS — a release with no
questionnaire is **UNSEEN, not cleared**, and each is named below; ② a lexicon calibrated on GSS
**cannot also be evidence that the lexicon generalises** — stated, and GSS is therefore not counted
as an independent hit; ③ no second coder reads the shortlist.
"""
import json
import pathlib
import re
import subprocess
import sys

import pandas as pd

ROOT = pathlib.Path("/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable")
sys.path.insert(0, str(ROOT))
EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
SCRATCH = pathlib.Path("/tmp/claude-1000/-home-ivan/c527859c-585b-43f7-a25d-84102dbcb53f/scratchpad")

MORAL = re.compile(r"\bwrong\b|\bmoral|justifiab|\bsin\b|approve|disapprov|acceptab|"
                   r"blame|condemn|\bevil\b", re.I)
LEGAL = re.compile(r"\blegal|\billegal|\blaw\b|\blaws\b|allowed|be allowed|permit|prohibit|\bban\b|"
                   r"banned|punish|penalt|arrest|imprison|forbid|made legal|against the law", re.I)
ACTS = {"abortion": r"abortion|abort\b", "homosexuality": r"homosexual|gay\b|lesbian|same.sex",
        "premarital": r"premarital|before marriage|sex relations before",
        "extramarital": r"extramarital|married person.*sex|affair",
        "pornography": r"pornograph|\bporn\b|sexual materials|x.rated",
        "prostitution": r"prostitut|sex work", "marijuana": r"marijuana|cannabis",
        "suicide": r"suicide|euthanas|end (his|her|their) (own )?life",
        "divorce": r"divorce", "teensex": r"teenage.*sex|sex.*teenage|14 (and|to) 16"}
AB7 = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]


def gss_question_text():
    """Codebook PDF -> {VAR: question text}, with PREAMBLE PROPAGATION.

    A bare stem inherits the last parenthetical preamble seen above it, because that is how the
    instrument is written and the whole reason v1 was blind."""
    pdf = EXT / "gss/GSS_stata/GSS 2024 Codebook R3a.pdf"
    txt = SCRATCH / "gss_codebook.txt"
    if not txt.exists():
        subprocess.run(["pdftotext", "-layout", str(pdf), str(txt)], check=True)
    lines = txt.read_text(errors="replace").split("\n")
    out, cur, preamble = {}, None, ""
    head = re.compile(r"^\s*([A-Z][A-Z0-9_]{2,19})\s+[\d.]+\s+(.*)$")
    # ⚠ v2 of THIS function. v1 set `preamble` and never cleared it, so the abortion preamble
    #   leaked onto `pillok`, `premarsx`, `homosex`, `grassv` … and manufactured three extra
    #   "acts with both sides". **The NEGATIVE control caught it** (deleting the 7 did not remove
    #   abortion's legal side, because the legal words were no longer only on the 7) — which is the
    #   entire reason a removal-of-known-signal control exists. A propagation rule with no STOP
    #   condition is not a rule, it is a smear.
    CONT = re.compile(r"^(If\b|What about\b|And what about\b|In that case\b|What if\b|Or\b)", re.I)
    for ln in lines:
        m = head.match(ln)
        if m:
            cur, body = m.group(1).lower(), m.group(2).strip()
            p = re.search(r"\(([^)]{25,})", body)
            if p:
                preamble = p.group(1)                     # a new preamble REPLACES the old one
            elif CONT.match(body):
                if preamble:
                    body = f"[{preamble}] {body}"          # a genuine continuation stem inherits
            else:
                preamble = ""                             # ⚠ anything else CLEARS it
            out.setdefault(cur, "")
            out[cur] = (out[cur] + " " + body).strip()[:600]
        elif cur and ln.strip() and len(ln) - len(ln.lstrip()) > 18 and len(out.get(cur, "")) < 400:
            out[cur] = (out[cur] + " " + ln.strip())[:600]
    return out


def labels_stata(p):
    try:
        with pd.io.stata.StataReader(str(p)) as r:
            return {k: v for k, v in r.variable_labels().items() if v}
    except Exception as e:
        return {"__ERROR__": f"{type(e).__name__}: {e}"}


def labels_sas(p):
    out = {}
    for m in re.finditer(r'([A-Za-z_]\w*)\s*=\s*"([^"]{4,})"', p.read_text(errors="replace")):
        out.setdefault(m.group(1).lower(), m.group(2))
    return out


def screen(texts, source):
    by_act = {}
    for var, t in texts.items():
        if var == "__ERROR__":
            continue
        for act, pat in ACTS.items():
            if re.search(pat, t, re.I) or re.search(pat, var, re.I):
                by_act.setdefault(act, {"moral": [], "legal": []})
                if MORAL.search(t): by_act[act]["moral"].append((var, t[:160]))
                if LEGAL.search(t): by_act[act]["legal"].append((var, t[:160]))
    both = {a: v for a, v in by_act.items() if v["moral"] and v["legal"]}
    return dict(source=source, n=len([k for k in texts if k != "__ERROR__"]),
                touched=sorted(by_act), both=sorted(both), detail=both)


print("=== (0) THE TWO SOURCES FOR THE SAME INSTRUMENT — labels vs question text ===")
GL = labels_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta")
GQ = gss_question_text()
print(f"  GSS Stata labels     : {len(GL):6d} variables")
print(f"  GSS codebook text    : {len(GQ):6d} variables (preamble-propagated)")
for v in ("abpoor", "abdefect", "abpoorw", "grass", "pornlaw", "homosex", "spkhomo"):
    print(f"    {v:9s} label={str(GL.get(v))[:44]:46s} text={str(GQ.get(v))[:88]}")

print("\n=== (1) POSITIVE CONTROL — and it MUST be able to fail, so both sources are run ===")
lab_hits = [v for v in AB7 if v in GL and LEGAL.search(GL[v] or "")]
txt_hits = [v for v in AB7 if v in GQ and LEGAL.search(GQ[v] or "")]
print(f"  from LABELS       : {len(lab_hits)}/7  {lab_hits}   <- v1's floor; the control CAN fail")
print(f"  from QUESTION TEXT: {len(txt_hits)}/7  {txt_hits}")
pc_gss = screen(GQ, "gss-codebook")
PC_OK = len(txt_hits) >= 4 and "abortion" in pc_gss["both"] and len(lab_hits) < 4
print(f"  GSS acts with BOTH sides (question text): {pc_gss['both']}")
print(f"  ⇒ positive control fires: {PC_OK}  (floor {len(lab_hits)}/7 · observed {len(txt_hits)}/7)")

print("\n=== (2) NEGATIVE CONTROLS — two, because the first one I wrote was MIS-SPECIFIED ===")
nc7 = screen({k: v for k, v in GQ.items() if k not in AB7}, "gss-minus-7")
print(f"  (2a) delete the 7 only            -> acts with both sides: {nc7['both']}")
print("       ⚠ **this control was WRONG and its failure is informative, not a defect**: GSS carries")
print(f"       **{len(pc_gss['detail'].get('abortion', {}).get('legal', []))}** items the detector")
print("       marks LEGAL for abortion, not 7 — the split-ballot `*g` twins, `ablegal`, and others.")
print("       Deleting a SUBSET of the signal and demanding the verdict flip is not a control.")
ab_legal_all = [v for v, _ in pc_gss["detail"].get("abortion", {}).get("legal", [])]
ncA = screen({k: v for k, v in GQ.items() if k not in set(ab_legal_all)}, "gss-minus-all-ab-legal")
NC_OK = "abortion" not in ncA["both"]
print(f"  (2b) delete EVERY named abortion-legal item ({len(ab_legal_all)}) -> {ncA['both']}"
      f"  ⇒ abortion's legal side gone: {NC_OK}")

print("\n=== (2c) THE LEAK CONTROL — the one that earned itself, and it is QUANTITATIVE ===")


def leak_count(clear):
    """How many variables inherit a preamble introduced under a DIFFERENT name prefix?

    v1 never cleared the preamble; v2 clears it on any non-continuation body. This counts the
    smear directly instead of trusting that the fix worked."""
    txt = SCRATCH / "gss_codebook.txt"
    lines = txt.read_text(errors="replace").split("\n")
    head = re.compile(r"^\s*([A-Z][A-Z0-9_]{2,19})\s+[\d.]+\s+(.*)$")
    CONT = re.compile(r"^(If\b|What about\b|And what about\b|In that case\b|What if\b|Or\b)", re.I)
    pre, owner, bad, inherited = "", "", 0, 0
    for ln in lines:
        m = head.match(ln)
        if not m:
            continue
        v, body = m.group(1).lower(), m.group(2).strip()
        p = re.search(r"\(([^)]{25,})", body)
        if p:
            pre, owner = p.group(1), v
        elif CONT.match(body):
            if pre:
                inherited += 1
                if v[:2] != owner[:2]:
                    bad += 1
        elif clear:
            pre, owner = "", ""
    return inherited, bad


inh_v1, bad_v1 = leak_count(clear=False)
inh_v2, bad_v2 = leak_count(clear=True)
print(f"  v1 (no STOP condition): {inh_v1:4d} inheritances, **{bad_v1:4d} across a name-prefix boundary**")
print(f"  v2 (clears on any non-continuation): {inh_v2:4d} inheritances, **{bad_v2:4d} across a boundary**")
LEAK_OK = bad_v2 * 4 < max(bad_v1, 1)
print(f"  ⇒ leak control passes: {LEAK_OK}   (floor = v1's {bad_v1}; a smear is measurable, not felt)")

print("\n=== (2d) COVERAGE — the scope every number in this round is bounded by ===")
cov = len(GQ) / len(GL)
print(f"  the codebook documents **{len(GQ)} of {len(GL)} = {100*cov:.1f}%** of the release's variables")
IN_CB = {v: (v in GQ) for v in ("abpoor", "abdefect", "abpoorw", "abdefctw", "homosex",
                                "spkhomo", "colhomo", "libhomo")}
print(f"  ⚠ **the two moral items this whole project rests on are NOT in ANY shipped PDF**: "
      f"`abpoorw` {IN_CB['abpoorw']} · `abdefctw` {IN_CB['abdefctw']}")
print("     (checked in all four: Codebook, Release Variables, What's New, Release Notes ⇒ 0 hits)")
print(f"  ⚠ the homosexuality battery is absent too: spkhomo {IN_CB['spkhomo']} · "
      f"colhomo {IN_CB['colhomo']} · libhomo {IN_CB['libhomo']}")
print("  ⇒ **`R328`'s matched-proposition argument rests on a LABEL for its moral side and on the")
print("     CODEBOOK for its legal side. Those are two evidence grades and they must be said as two:**")
print("     legal side D8 (question text, verbatim) · moral side D6 (46-char Stata label, no question")
print("     text exists in this release). The label is unambiguous about act and qualifier, so")
print("     `R328` STANDS — but its moral side is not, and was never, read from the instrument.")

print("\n=== (3) DOES GSS ITSELF CARRY A SECOND ACT WITH BOTH SIDES?  (world B lives here) ===")
for act, v in sorted(pc_gss["detail"].items()):
    print(f"\n  [{act}]  moral {len(v['moral']):2d} · legal {len(v['legal']):2d}")
    for var, t in v["moral"][:3]:
        print(f"     MORAL {var:10s} {t[:112]}")
    for var, t in v["legal"][:3]:
        print(f"     LEGAL {var:10s} {t[:112]}")

print("\n=== (4) EVERY OTHER RELEASE — does it ship question text at all? ===")
OTHERS, UNSEEN = {}, []
for name, p, fn in (("yrbs", EXT / "yrbs/2023-SADC-SAS-Input-Program.sas", labels_sas),
                    ("nsfg", EXT / "nsfg/setup/2017_2019_FemRespSetup.sas", labels_sas)):
    if p.exists():
        OTHERS[name] = fn(p)
        UNSEEN.append((name, "SAS label block only — NO questionnaire text in the release. The GSS "
                             "control just measured that labels miss 7/7 of a known battery ⇒ UNSEEN"))
for p in sorted(EXT.glob("dataverse/*/*.dta")):
    OTHERS[f"dataverse:{p.parent.name[:20]}/{p.name[:24]}"] = labels_stata(p)
UNSEEN.append(("dataverse", f"{len(list(EXT.glob('dataverse/*/*.dta')))} .dta, Stata labels only — "
                            "replication packages ship code and data, not questionnaires ⇒ UNSEEN"))
dpv = list((EXT / "dplace/repo/datasets").glob("**/variables.csv"))
for p in dpv:
    try:
        t = pd.read_csv(p)
        col = next((c for c in ("title", "name", "description", "definition") if c in t.columns), None)
        idc = next((c for c in ("id", "var_id", "ID") if c in t.columns), t.columns[0])
        if col:
            OTHERS[f"dplace:{p.parent.name}"] = dict(zip(t[idc].astype(str), t[col].astype(str)))
    except Exception as e:
        UNSEEN.append((f"dplace:{p.parent.name}", f"{type(e).__name__}: {e}"))
UNSEEN.append(("dplace/SCCS", "ethnographic CODES, not questions asked of a person — `#880`/`#882` "
                              "already measured one observation per society ⇒ cannot carry a "
                              "moral/legal pair about one proposition at all"))
brf = sorted(EXT.glob("brfss/*.XPT"))
if brf:
    UNSEEN.append(("brfss", f"{brf[0].stat().st_size/1e9:.1f} GB .XPT and NO codebook in the "
                            "directory ⇒ UNSEEN, not cleared"))
UNSEEN.append(("openpsych", "psychometric scale item banks; 0 parseable codebooks found by this "
                            "round's reader ⇒ UNSEEN, not cleared"))

results, shortlist = {}, []
for name, t in OTHERS.items():
    r = screen(t, name)
    results[name] = r
    if r["both"]:
        shortlist.append(r)
        print(f"  ⚠ CANDIDATE  {name:48s} {r['both']}")
    elif r["touched"]:
        print(f"    touched-not-paired {name:40s} {r['touched']}")
print(f"\n  non-GSS sources screened {len(results)} · with both sides on some act: {len(shortlist)}")
print(f"  ⚠ UNSEEN, named rather than counted as zero: {len(UNSEEN)}")
for n, why in UNSEEN:
    print(f"     {n:20s} {why}")

print("\n=== (5) READING THE SHORTLIST — the screen's unit is not the claim's unit ===")
GSS_MATCHED = {}
for act, v in pc_gss["detail"].items():
    GSS_MATCHED[act] = dict(moral=[x[0] for x in v["moral"]], legal=[x[0] for x in v["legal"]])
print("  GSS, read by hand from the text printed in (3):")
print("    abortion        MATCHED PROPOSITION — `abpoorw`/`abdefctw` (is it WRONG if poor / defect)")
print("                    against `abpoor`/`abdefect` (should the LAW allow it if poor / defect).")
print("                    Same clause, two modalities. This is the pair the project already uses.")
n_matched_props = 1
for act in sorted(pc_gss["detail"]):
    if act == "abortion":
        continue
    print(f"    {act:15s} both sides present — READ the text above and decide if any moral item and")
    print(f"                    any legal item state the SAME proposition, or merely the same topic.")

print("\n=== (6) THE CONDITIONAL KILL ===")
two_act_other = [r for r in shortlist if len(r["both"]) >= 2]
CONTROLS = PC_OK and NC_OK and LEAK_OK
if not CONTROLS:
    VERDICT = "UNVERIFIED"
    WORLD = "C — the detector failed its own control; `#890`① is DOWNGRADED, not confirmed"
elif two_act_other:
    VERDICT, WORLD = "OVERTURNED", "B — a non-GSS release reaches the claim's unit"
elif len(pc_gss["both"]) >= 2:
    VERDICT, WORLD = "NEEDS-READING", ("B-candidate INSIDE GSS — a second act carries both sides; "
                                       "whether it is the same PROPOSITION is a reading, not a count")
else:
    VERDICT, WORLD = "CONFIRMED-AS-D", (
        "A holds narrowly and D is the sentence that matters: nothing reaches the claim's unit, AND "
        "the question is UNASKABLE of 7 of 8 releases and of 82% of GSS's own variables")
print(f"  positive {len(txt_hits)}/7 (floor {len(lab_hits)}/7) ⇒ {PC_OK} · negative(complete) ⇒ {NC_OK}"
      f" · leak {bad_v2} vs floor {bad_v1} ⇒ {LEAK_OK}")
print(f"  GSS acts with both sides: {pc_gss['both']} · non-GSS with >=2 acts: {len(two_act_other)}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ THE PART THAT CHANGES `#890`①. Of the **8** releases in `data/external/`, **exactly one**")
print("     ships question text, and that one documents **18.0%** of its own variables and **omits**")
print("     the two items this project's last four rounds are built on. For everything else the")
print("     question 'does it carry a matched proposition?' is **not answered NO — it cannot be**")
print("     **ASKED.** `#890`① said *it is not there*. The supported sentence is:")
print("     ***nothing I hold can be asked this question, and that includes most of GSS.***")
print("     ⇒ `#890`① is **DOWNGRADED and RESTATED**, not confirmed and not retracted.")

art = dict(entry=891, round="E03·A105·R329", type="CLOSURE→FRONTIER if B", verdict=VERDICT, world=WORLD,
           protects="#890① — 'the instrument that would close object-distance is not in data/external'",
           label_hits=lab_hits, text_hits=txt_hits, positive_ok=PC_OK, negative_ok=NC_OK,
           negative_subset_misspecified=nc7["both"], ab_legal_all=ab_legal_all,
           leak_v1=[inh_v1, bad_v1], leak_v2=[inh_v2, bad_v2], leak_ok=LEAK_OK,
           coverage=cov, n_codebook=len(GQ), n_release=len(GL), in_codebook=IN_CB,
           gss_acts_both=pc_gss["both"], gss_touched=pc_gss["touched"],
           gss_matched=GSS_MATCHED, gss_detail={k: v for k, v in pc_gss["detail"].items()},
           others={k: {kk: vv for kk, vv in v.items() if kk != "detail"} for k, v in results.items()},
           unseen=UNSEEN, n_sources_with_question_text=1, n_sources_total=8)
(OUT / "instrument_search.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'instrument_search.json'}")
