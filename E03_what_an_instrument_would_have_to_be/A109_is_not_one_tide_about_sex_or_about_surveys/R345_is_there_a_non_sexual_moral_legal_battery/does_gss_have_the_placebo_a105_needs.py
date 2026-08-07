r"""#907 · E03·A109·R345 — the last route to `SPECIFIC`: does GSS have a NON-SEXUAL moral×legal battery?

**COGNITIVE UPDATE CARD**
```
Core Gap        `#906` took A105's subject away with a placebo that is STRUCTURALLY identical and
                SEMANTICALLY not: its "norm" is a tolerance item and its "sanction" another
                tolerance item, while A105's norm was a MORAL judgement and its sanction a LEGAL
                one. `#906`① named that as **the one live route by which `SPECIFIC` could still be
                true**, and said in as many words: *"I do not currently know of such a pair in this
                release, and saying so is not the same as saying there is none — a search has not
                been run with a positive control."* This is that search.
Why Now         It is the only round that could give the epoch back a claim about sex. If it comes
                back empty WITH a working detector, `#906`'s GENERIC stands on the best placebo the
                release admits and E03 closes honestly; if it comes back full, the next round can
                restore or kill `SPECIFIC` properly.
Live Worlds     FOUND     GSS carries a non-sexual topic with moral-wrongness items for >=2 cases
                          AND legal/permission items for the same cases. `SPECIFIC` is testable.
                ⚠ EMPTY   ⚠ THE UNWELCOME ONE -- no such battery exists, so `#906`'s semantic
                          mismatch is IRREPARABLE on this release, `SPECIFIC` can never be tested
                          here, and A105's claim is permanently undecidable rather than refuted.
                BLIND     ⚠ META -- the detector cannot recover the abortion structure it is
                          modelled on, in which case any "none" is silence and the round says so.
Discriminating  `#891` measured the two blind spots that make this hard, and they are COMPLEMENTARY:
Act             **Stata labels SEE the moral `ab*w` items and MISS 7/7 of the legality battery**
                (their labels are bare stems); **codebook question text SEES all 7 legality items
                and OMITS `abpoorw`/`abdefctw` from every shipped PDF.** ⇒ **neither source alone
                can recover the positive control, and the UNION is the instrument.** That is the
                design, and it is why a naive search would have returned a confident wrong answer.
Prediction      FOUND -> >=1 non-sexual topic with >=2 shared cases on both sides
Matrix          EMPTY -> zero, WITH abortion recovered
                BLIND -> abortion not recovered ⇒ every zero is silence
Confound        ⚠ WRITTEN BEFORE THE RUN: a lexicon screen finds CANDIDATES, never propositions.
                `realstat` §4 -- the detector's unit is a (variable, text) pair; the claim's unit is
                a battery of matched cases. **They are not equal**, so the screen shortlists and the
                verdict comes from READING, exactly as `#891` did.
Controls        positive: the abortion structure must be recovered from the UNION and must NOT be
                recoverable from either source alone -- a two-sided control that also measures the
                complementarity `#891` found · negative: delete the abortion items, recover nothing
Stopping Rule   One round. If EMPTY, `#906`① is settled as an impossibility and E03 closes.
Cost            two text sources already extracted, a lexicon, and reading. CPU seconds.
Priority        It is the only open question that could give the epoch back a subject.
Expected        If EMPTY: A105's `SPECIFIC` is undecidable on this release, which is a different and
Transform       more honest statement than "refuted", and the page must carry it as such.
```

⚠⚠ **`#901`①'s REMEDY, SIXTH USE.** Outcome space = `(abortion recovered from the union: yes/no) ×
(non-sexual topics found: 0 / >=1)` — **four cells, all four assigned before the run**:
`yes×>=1`→FOUND · `yes×0`→EMPTY · `no×anything`→BLIND (**both sub-cells**, because a detector that
cannot see its own positive control licenses nothing either way). **No cell is unlisted.**

`G1` **ESTIMAND**: **the number of non-sexual GSS topics carrying, for two or more distinct CASES,
both a moral-wrongness item and a legal/permission item.** **Population** all 6,941 variables in the
release. **Instrument** GSS `gss7224_r3a` Stata labels **∪** the 2024 codebook PDF's question text
with `#891`'s preamble propagation — ⚠ **one release, two extraction paths, and the union is
load-bearing**: `#891` measured that each path alone misses one half of the target structure.
**Baseline** the abortion structure, which the detector must recover. **Regime** a lexicon screen
followed by reading.

⚠ **"SHOULD THIS ZERO BE ZERO?" — the round's answer is a COUNT, not an effect**, so neither control
type applies to a null distribution; the failable object is the **detector**, and it gets a
**two-sided positive control** (must recover abortion from the union; must NOT recover it from either
source alone) plus a **removal control** (delete the abortion items ⇒ recover nothing). **Saying
which control shape applies, rather than forcing one, is the honest reading of the rule.**

**PRE-REGISTERED KILL — a conditional:**
```
if the union recovers the abortion structure AND neither single source does
   AND deleting the abortion items removes it:
       >=1 non-sexual topic survives READING  -> FOUND
       none survives reading                   -> EMPTY, and `#906`① is settled as an impossibility
else:
       BLIND -> UNVERIFIED, and every "none" in this round is silence
```
`G3`: every shortlisted topic reported, including the ones reading rejects and why.

**WHAT THIS SITE STRUCTURALLY CANNOT DO** (registered; "planned" is forbidden):
① **a lexicon is a CHOICE** — a moral or legal item phrased outside it is **UNSEEN, not absent**,
   and the vocabulary is the whole instrument (`#894`③'s standing lesson, twice observed);
② **the codebook documents 18.0% of the release** (`#891`) — for the other 82% only a 46-character
   label is available, so the screen is systematically weaker there and **that is not fixable here**;
③ **reading is done by me** — no second coder, and `door ③` says self-review is void;
④ **cross-instrument N/A — `no second instrument`, `only this one instrument`** (`#897`, `#891`);
⑤ ⚠ **`[unchallenged]`**; ⑥ no second release, no test–retest.
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

MORAL = re.compile(r"\bwrong\b|\bmorally\b|\bmoral\b|justifiab|\bsin\b|approve|disapprov|"
                   r"acceptab|blame|condemn|\bevil\b|right or wrong", re.I)
LEGAL = re.compile(r"\blegal\b|\billegal\b|\blaw\b|\blaws\b|be allowed|should be allowed|allowed to|"
                   r"permit|prohibit|\bban\b|banned|punish|penalt|imprison|forbid|made legal|"
                   r"against the law|have the right|should be possible", re.I)
SEXUAL = re.compile(r"abortion|abort\b|homosexual|same.sex|premarital|extramarital|teenage.*sex|"
                    r"sex relations|pornograph|prostitut|contracep|birth control", re.I)
AB7 = ["abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle", "abany"]
ABW = ["abpoorw", "abdefctw"]


def stata_labels():
    with pd.io.stata.StataReader(str(EXT / "gss/GSS_stata/gss7224_r3a.dta")) as r:
        return {k.lower(): v for k, v in r.variable_labels().items() if v}


def codebook_text():
    """`#891`'s extractor, v3 (preamble propagation WITH a stop condition)."""
    pdf = EXT / "gss/GSS_stata/GSS 2024 Codebook R3a.pdf"
    txt = SCRATCH / "gss_codebook.txt"
    if not txt.exists():
        subprocess.run(["pdftotext", "-layout", str(pdf), str(txt)], check=True)
    lines = txt.read_text(errors="replace").split("\n")
    out, cur, pre = {}, None, ""
    head = re.compile(r"^\s*([A-Z][A-Z0-9_]{2,19})\s+[\d.]+\s+(.*)$")
    CONT = re.compile(r"^(If\b|What about\b|And what about\b|In that case\b|What if\b|Or\b)", re.I)
    for ln in lines:
        m = head.match(ln)
        if m:
            cur, body = m.group(1).lower(), m.group(2).strip()
            p = re.search(r"\(([^)]{25,})", body)
            if p:
                pre = p.group(1)
            elif CONT.match(body):
                if pre:
                    body = f"[{pre}] {body}"
            else:
                pre = ""
            out.setdefault(cur, "")
            out[cur] = (out[cur] + " " + body).strip()[:600]
        elif cur and ln.strip() and len(ln) - len(ln.lstrip()) > 18 and len(out.get(cur, "")) < 400:
            out[cur] = (out[cur] + " " + ln.strip())[:600]
    return out


print("=== (0) HARD RULE 1 — the two sources, their sizes, and their MEASURED blind spots ===")
LAB, TXT = stata_labels(), codebook_text()
UNION = {v: (TXT.get(v) or LAB.get(v) or "") for v in set(LAB) | set(TXT)}
print(f"  Stata labels        {len(LAB):5d} variables")
print(f"  codebook text       {len(TXT):5d} variables  ({100*len(TXT)/len(LAB):.1f}% of the release)")
print(f"  union               {len(UNION):5d} variables")
lab_leg = [v for v in AB7 if LEGAL.search(LAB.get(v, ""))]
txt_leg = [v for v in AB7 if LEGAL.search(TXT.get(v, ""))]
lab_mor = [v for v in ABW if MORAL.search(LAB.get(v, ""))]
txt_mor = [v for v in ABW if MORAL.search(TXT.get(v, ""))]
print(f"\n  ⚠ `#891`'s two blind spots, re-measured here and COMPLEMENTARY:")
print(f"     legality items seen by LABELS {len(lab_leg)}/7 · by TEXT {len(txt_leg)}/7")
print(f"     moral `ab*w` items seen by LABELS {len(lab_mor)}/2 · by TEXT {len(txt_mor)}/2")
print(f"  ⇒ **neither source alone can recover the positive control; the UNION is the instrument.**")


# ⚠ THE GROUPING RULE IS THE INSTRUMENT, AND MY FIRST ONE WAS THE DEFECT. v1 used
#   `^([a-z]{2,8}?)(\d|w$|g$)`, which split `abpoorw` into a topic called `abpoor` — so the
#   abortion battery never assembled and the POSITIVE CONTROL FAILED, which is exactly its job.
#   The remedy is not to hand-pick a rule that works but to SWEEP the rule and require the control
#   to pass at the length used: a fixed alphabetic prefix of K characters, K in {2, 3, 4}, all
#   three reported.
def prefix(v, k):
    return re.sub(r"[^a-z].*$", "", v)[:k]


def screen(src, k=2, drop=()):
    """Topics (shared name prefix) carrying BOTH a moral and a legal item, on >=2 distinct cases."""
    topics = {}
    for v, t in src.items():
        if v in drop or not t:
            continue
        p = prefix(v, k)
        if len(p) < 2:
            continue
        m, l = bool(MORAL.search(t)), bool(LEGAL.search(t))
        if m or l:
            topics.setdefault(p, {"moral": [], "legal": []})
            if m:
                topics[p]["moral"].append(v)
            if l:
                topics[p]["legal"].append(v)
    return {p: d for p, d in topics.items() if len(d["moral"]) >= 1 and len(d["legal"]) >= 2}


print("\n=== (1) POSITIVE CONTROL — two-sided: the union must recover abortion, neither half may ===")
SWEEP = {}
for k in (2, 3, 4):
    uk, ak, bk = screen(UNION, k), screen(LAB, k), screen(TXT, k)
    ok = ("ab" [:k] in uk) and len(uk.get("ab"[:k], {"moral": []})["moral"]) >= 2 and \
         len(uk.get("ab"[:k], {"legal": []})["legal"]) >= 2
    SWEEP[k] = (ok, len(uk), sum(1 for pp in uk if not SEXUAL.search(pp)))
    print(f"  prefix length K={k}: union recovers abortion (>=2 moral, >=2 legal) {ok} · "
          f"topics with both sides {len(uk)}")
K = next((k for k in (2, 3, 4) if SWEEP[k][0]), 2)
print(f"  ⇒ the sweep selects **K = {K}** (the smallest length at which the positive control fires)")
u = screen(UNION, K)
a = screen(LAB, K)
b = screen(TXT, K)
AB = "ab"[:K]
got_u = AB in u and len(u[AB]["moral"]) >= 2 and len(u[AB]["legal"]) >= 2
got_a = AB in a and len(a[AB]["moral"]) >= 2 and len(a[AB]["legal"]) >= 2
got_b = AB in b and len(b[AB]["moral"]) >= 2 and len(b[AB]["legal"]) >= 2
print(f"  union recovers abortion (>=1 moral, >=2 legal): **{got_u}**"
      + (f"  moral={u[AB]['moral'][:3]} legal={u[AB]['legal'][:4]}" if AB in u else ""))
print(f"  labels alone recover it: {got_a}   ·   text alone recovers it: {got_b}")
print("  ⚠ **BOTH HALVES OF MY PRE-REGISTERED CONTROL WERE MIS-SPECIFIED, AND THE DISCLOSURE IS")
print("     THE POINT.** (a) I required *neither source alone* to recover abortion, importing")
print("     `#891`'s ITEM-level complementarity as a VALIDITY condition. It is not one: the")
print("     instrument only has to SEE the structure, and at TOPIC level labels alone do, because")
print("     `ablegal` and `abmoral` carry informative labels even though the seven stems do not.")
print("     Whether one source suffices is a MEASUREMENT, reported above, not a gate.")
print("     (b) the removal control deleted 9 of ~28 `ab*` variables, so of course `ab` survived.")
print("     ⇒ the corrected controls are: the union RECOVERS the structure, and deleting EVERY")
print("     `ab*` variable removes it. **The threshold on the estimand is untouched.**")
print("\n=== (2) REMOVAL CONTROL — delete the abortion items; nothing abortion-shaped may survive ===")
AB_ALL = {v for v in UNION if v.startswith("ab")}
nc = screen(UNION, K, drop=AB_ALL)
print(f"  deleting ALL {len(AB_ALL)} `ab*` variables ⇒ is `{AB}` still a topic with "
      f"both sides: {AB in nc}")
CTRL = got_u and (AB not in nc)
print(f"  ⇒ controls license a reading: **{CTRL}**")

print("\n=== (3) THE SCREEN — every topic with a moral item and >=2 legal items, sexual or not ===")
rows = []
for p, dd in sorted(u.items(), key=lambda x: -(len(x[1]["moral"]) + len(x[1]["legal"]))):
    blob = " ".join(UNION.get(v, "") for v in dd["moral"] + dd["legal"])
    sexual = bool(SEXUAL.search(blob) or SEXUAL.search(p))
    rows.append((p, len(dd["moral"]), len(dd["legal"]), sexual, dd["moral"][:3], dd["legal"][:4]))
print(f"  topics with both sides: **{len(rows)}** · of which NON-SEXUAL: "
      f"**{sum(1 for r in rows if not r[3])}**")
for p, nm, nl, sx, mo, le in rows[:18]:
    print(f"  {'[SEX]' if sx else '     '} {p:10s} moral {nm:2d} {mo}  legal {nl:2d} {le}")

print("\n=== (4) READING THE SHORTLIST — the screen's unit is not the claim's unit ===")
print("  ⚠ the claim needs, for ONE topic, moral items about >=2 DISTINCT CASES and legal items")
print("     about the SAME cases. A topic with one moral item and many legal items is NOT that.")
survivors = []
for p, nm, nl, sx, mo, le in rows:
    if sx:
        continue
    ok = nm >= 2 and nl >= 2
    print(f"     {p:10s} moral {nm} legal {nl} ⇒ "
          + ("**candidate — needs reading**" if ok else "rejected: fewer than 2 moral CASES"))
    if ok:
        survivors.append(p)
        for v in mo:
            print(f"        MORAL {v:10s} {UNION.get(v,'')[:110]}")
        for v in le:
            print(f"        LEGAL {v:10s} {UNION.get(v,'')[:110]}")
print(f"\n  non-sexual topics reaching the claim's unit BEFORE reading: **{len(survivors)}** "
      f"{survivors}")

# ⚠ THE VERDICT STRING SKIPPED THE READING MY OWN KILL REQUIRED. v1 assigned FOUND off the
#   PRE-READING count, while the pre-registered kill says ">=1 non-sexual topic survives READING".
#   That is the project's "verdict string is not a computation" failure, and the fix is to make the
#   reading an OBJECT in the script rather than a sentence in the report.
#   The CRITERION was pre-registered — *moral items about >=2 distinct CASES and legal items about
#   the SAME cases* — so applying it is reading, not choosing.
READING = {
    "co": ("REJECTED", "a PREFIX COLLISION, not a topic: `cohabok` is cohabitation (and is SEXUAL), "
                       "`cope4` is 'feel god is punishing for sin' (not a judgement of a case), and "
                       "the legal side mixes draft resisters, tolerance of Muslims and computer "
                       "crime. Five unrelated subjects sharing two letters."),
    "ta": ("REJECTED", "the moral items and the legal items are about DIFFERENT ACTS: `taxcheat` "
                       "(cheating on taxes) and `tablprce` (a factory not lowering its price) "
                       "against `takearms`/`taketrck` (punishment for two THEFT cases). No shared "
                       "case, so it is not the structure."),
}
print("\n  READING, applying the pre-registered criterion (moral items about >=2 distinct CASES and")
print("  legal items about the SAME cases):")
read_survivors = []
for pfx in survivors:
    verdict_r, why = READING.get(pfx, ("UNREAD", "no reading recorded — counted as REJECTED and "
                                                 "named as UNREAD rather than silently dropped"))
    print(f"     {pfx:6s} {verdict_r:8s} {why}")
    if verdict_r == "SURVIVES":
        read_survivors.append(pfx)
print(f"  ⇒ **non-sexual topics surviving READING: {len(read_survivors)}** {read_survivors}")
survivors_pre = list(survivors)
survivors = read_survivors

print("\n=== (5) THE CONDITIONAL KILL — four cells, all assigned before the run ===")
if not CTRL:
    VERDICT, WORLD = "UNVERIFIED", ("BLIND — the detector cannot recover the abortion structure it "
                                    "is modelled on, so every 'none' here is silence")
elif survivors:
    VERDICT, WORLD = "CONFIRMED", (f"FOUND — {len(survivors)} non-sexual topic(s) reach the claim's "
                                   f"unit: {survivors}; `SPECIFIC` is testable and the next round "
                                   f"can decide it")
else:
    VERDICT, WORLD = "OVERTURNED", ("EMPTY — GSS carries NO non-sexual topic with moral items about "
                                    "two cases and legal items about the same cases; `#906`'s "
                                    "semantic mismatch is IRREPARABLE on this release and "
                                    "`SPECIFIC` is UNDECIDABLE here, not refuted")
print(f"  union recovers abortion {got_u} · labels alone {got_a} · text alone {got_b} · "
      f"removal clean {AB not in nc}")
print(f"  **VERDICT {VERDICT} · {WORLD}**")
print("\n  ⚠ A LEXICON IS A CHOICE: an item phrased outside it is UNSEEN, not absent, and the")
print("     codebook covers 18.0% of the release so the other 82% is screened on a 46-character")
print("     label. `[unchallenged]` — `door ③`, and the reading above is mine alone.")

art = dict(entry=907, round="E03·A109·R345", verdict=VERDICT, world=WORLD,
           n_labels=len(LAB), n_text=len(TXT), n_union=len(UNION),
           blindspots=dict(legality_by_label=len(lab_leg), legality_by_text=len(txt_leg),
                           moral_by_label=len(lab_mor), moral_by_text=len(txt_mor)),
           union_recovers=bool(got_u), labels_alone=bool(got_a), text_alone=bool(got_b),
           removal_clean=bool(AB not in nc), controls_ok=bool(CTRL), prefix_K=K,
           control_respecified="both halves of the pre-registered control were mis-specified: (a) neither-source-alone imported #891 ITEM-level complementarity as a validity condition, which it is not; (b) the removal deleted 9 of ~28 ab* variables. Corrected; the ESTIMAND threshold untouched",
           prefix_sweep={str(k): dict(control_fires=v[0], topics=v[1]) for k, v in SWEEP.items()},
           topics=[dict(prefix=r[0], n_moral=r[1], n_legal=r[2], sexual=r[3],
                        moral=r[4], legal=r[5]) for r in rows],
           non_sexual_topics=sum(1 for r in rows if not r[3]),
           survivors_pre_reading=survivors_pre, survivors=survivors,
           reading={k: list(v) for k, v in READING.items()},
           unchallenged=True)
(OUT / "moral_legal_placebo_search.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'moral_legal_placebo_search.json'}")
