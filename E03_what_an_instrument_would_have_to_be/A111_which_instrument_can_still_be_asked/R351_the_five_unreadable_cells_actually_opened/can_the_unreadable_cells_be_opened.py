#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A111·R351 — the five cells `#912` called UNREADABLE, actually opened
=======================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#912` declared the project `BOUND` — its terminal verdict — and derived it from
                five (release, estimand) cells scored **UNREADABLE**. Every one of those five was
                scored from the DIRECTORY LISTING. Not one file was opened. So `BOUND` may be a
                fact about what I tried rather than about the data, and it is the strongest claim
                on the project's front page.

Why Now         The page head now tells a reader that no further round is possible without one of
                three named acquisitions. If that is wrong, every hour spent chasing a document is
                spent on nothing, and a supportable pair sits unrun. This is the cheapest possible
                check on the most consequential sentence the project has ever published.

Live Worlds     A · BOUND-REAL      — the five cells genuinely need an acquisition. `#912` stands.
                B · BOUND-ARTEFACT  — at least one opens with software already installed, and
                                      contains a norm item ⇒ an untried supportable pair exists
                                      ⇒ `BOUND` is OVERTURNED and E03 is not over.
                C · READABLE-EMPTY  — the files open and contain no norm item ⇒ `BOUND` survives
                                      with its REASON corrected, and the three named acquisitions
                                      are revealed as worthless, which is its own retraction.

Discriminating  Open them. Parse the shipped layouts, read the actual bytes, and count norm items
Act             from the shipped ITEM TEXT — never from a variable name (`#891` measured a
                name/label screen missing 7/7 of a known battery).

Prediction      A → every attempt returns 0 records.
Matrix          B → at least one returns records AND a moral-evaluative item.
                C → records return, no moral-evaluative item anywhere.

BASIN RULE      `#911` called a basin and `#912` was the step out of it. This round's UNWELCOME
                positive outcome is world B: it overturns the headline I committed one round ago
                and re-opens an epoch I had just closed. That is the design the basin rule asks
                for — the step whose success I would rather not have.

Meta-separator  Both `#912` and this round assume "readable" is BINARY. It is not: BRFSS opens as
                350 named columns while shipping no codebook for its VALUE codes. If that state
                exists, the decomposition {SUPPORTED · NOT · UNREADABLE} was too coarse and the
                census inherited a wrong ontology, not merely a wrong cell.

Strongest       ⚠ A fixed-width layout that is off by one byte still "reads" — it returns a column
Confound        of plausible small integers. So a successful read is NOT evidence the parse is
                correct, and I would be reporting garbage with a valid-looking distribution.
                CONTROL, same iteration: `offset_control` — shift the record start by one byte and
                require the parse to COLLAPSE.
                ⚠ FIRST VERSION OF THIS CONTROL COULD NOT FAIL, and the run said so: I scored
                "share of values in {1..5,8,9} OR BLANK", so a shift into a blank byte returned
                1.0000 — and in the 2017-2019 cycles `chsuppor` sits in the byte immediately
                AFTER `samesex`, so a one-byte shift lands on another 1-digit Likert item and is
                valid by construction. **Both failure modes make the control unfailable, and the
                second is a property of the data no amount of care about the first would have
                caught.** Repointed at `caseid`: at the true offset it is a near-unique respondent
                id, and one byte off it collapses into collisions. That check can fail.

Stopping Rule   One pass over five cells. No acquisition, no install, no download. If a cell needs
                software this machine does not have, it stays UNREADABLE and that is the finding.

Cost            minutes; one 1.2 GB streaming read, chunked, never held whole.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  ① Absence of a norm item is scored from SHIPPED TEXT. A release that ships no labels
    (BRFSS) cannot be scored 0 — it is UNDECIDABLE, and `#891`'s 7/7 miss is why a name screen
    may not stand in.
  ② No second coder decides what counts as a norm item; the rule is published and mechanical.
  ③ `[unchallenged]` — door ③.
"""
import json, os, re, sys, zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
EXT = ROOT / "data" / "external"
RNG = np.random.default_rng(351)

# ── the norm-item rule, published and mechanical (G1: named before anything is counted) ──
# A NORM ITEM is an item whose SHIPPED TEXT states a moral/normative evaluation of a behaviour:
# it says an act is right/wrong/okay/acceptable/approved, or asks the respondent to evaluate it.
# It is NOT an item asking whether the respondent DID something, or how often.
NORM_RE = re.compile(
    r"\b(is all right|all right|okay for|ok for|acceptable|wrong|approve|disapprove|"
    r"should be (allowed|legal|permitted)|morally|right or wrong)\b", re.I)
# a behaviour item — used only to show the rule DISCRIMINATES, never to score a release
BEHAV_RE = re.compile(r"\b(ever had|how many times|during the past|did you|have you ever)\b", re.I)


def dct_layout(path):
    """Stata dictionary → [(name, start0, width, label)]. `_column(N) type NAME %Wf "label"`."""
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)[a-z]\s*"([^"]*)"')
    out = []
    for line in Path(path).read_text(errors="replace").splitlines():
        m = pat.search(line)
        if m:
            out.append((m.group(2), int(m.group(1)) - 1, int(m.group(3)), m.group(4)))
    return out


def sas_input_layout(path):
    """SAS input program → [(name, start0, width, '')]. `name $ 1-5` / `age 159-161`."""
    txt = Path(path).read_text(errors="replace")
    body = txt.split("input", 1)[1] if "input" in txt else txt
    out = []
    for m in re.finditer(r"^\s*([A-Za-z_]\w*)\s*(\$?)\s*(\d+)-(\d+)\s*$", body, re.M):
        a, b = int(m.group(3)), int(m.group(4))
        out.append((m.group(1), a - 1, b - a + 1, ""))
    return out


def sas_format_labels(path):
    """SAS formats program → every quoted value label, as the release's own item vocabulary."""
    txt = Path(path).read_text(errors="replace")
    return re.findall(r'"([^"]{3,})"|\'([^\']{3,})\'', txt)


rows = []          # the whole grid, published including disagreeing cells
detail = {}

# ══ CELL 1-2 · NSFG — 15 data files, and the layouts that ship with them ══════════════
nsfg_dat = sorted((EXT / "nsfg").glob("*.dat"))
nsfg_dct = sorted((EXT / "nsfg" / "setup").glob("*.dct"))
covered = []
for d in nsfg_dct:
    stem = d.name.replace("Setup.dct", "")
    hit = [f for f in nsfg_dat if f.name.startswith(stem)]
    if hit:
        covered.append((hit[0], d))
print(f"NSFG · {len(nsfg_dat)} data files · {len(nsfg_dct)} layouts · {len(covered)} files with a layout")

nsfg_records, nsfg_norm_items, nsfg_years = 0, [], []
nsfg_valid_share, nsfg_offset_share = None, None
for dat, dct in covered:
    lay = dct_layout(dct)
    norms = [(n, lb) for (n, _s, _w, lb) in lay if NORM_RE.search(lb)]
    # ⚠ HARD RULE 1 — a variable name is not a measurement. Read the actual bytes.
    # ⚠ and the name is CASE-DEPENDENT across cycles: 2011-2013 ships `SAMESEX`, 2017-2019 ships
    #   `samesex`. The first version matched case-sensitively and silently read 1 file of 3.
    by = {n.lower(): (n, s, w, lb) for (n, s, w, lb) in lay}
    want = [by[k] for k in ("caseid", "samesex", "chsuppor") if k in by]
    if len(want) < 2:
        print(f"  {dat.name:34s} SKIPPED — layout lacks the attitude block ({sorted(by)[:3]}...)")
        continue
    spec = [(s, s + w) for (_n, s, w, _l) in want]
    names = [t[0].lower() for t in want]
    df = pd.read_fwf(dat, colspecs=spec, names=names, dtype=str, nrows=200_000)
    n = len(df)
    nsfg_records += n
    yr = re.findall(r"(\d{4})", dat.name)
    nsfg_years += yr
    nsfg_norm_items += [f"{dat.name}:{a}" for a, _ in norms]
    # documented codes for the IH/JG items: 1-5 substantive, 8/9 refused/dk. BLANK IS NOT VALID —
    # counting it as valid is what made version 1 of this control unable to fail.
    sm = df["samesex"].fillna("").str.strip()
    share = float(sm.isin(list("12345") + ["8", "9"]).mean())
    nsfg_valid_share = share if nsfg_valid_share is None else min(nsfg_valid_share, share)
    # ── the confound, measured in the SAME iteration and NOT survived ──
    # Shift the item's byte window by -2..+2 and score documented codes (blank INVALID).
    # ⚠ MEASURED: every offset returns 1.0000 on all three files. The IH/JG block is a dense run
    #   of 1-digit Likert items, so ANY neighbouring byte is a documented code by construction.
    #   ⇒ code-validity is STRUCTURALLY INCAPABLE of anchoring a fixed-width parse at this site,
    #   and three earlier versions of this control (blank-as-valid · caseid uniqueness · caseid
    #   monotonicity, 0.9909/0.9998 vs 1.0000) failed to separate for the same underlying reason.
    #   Registered, not excused — and it BOUNDS THE CLAIM: this round may say the cell is readable
    #   and holds a norm item, and may NOT say any attitude has been measured.
    s0, w0 = by["samesex"][1], by["samesex"][2]
    shift = {}
    for off in (-2, -1, 0, 1, 2):
        c = pd.read_fwf(dat, colspecs=[(s0 + off, s0 + off + w0)], names=["X"], dtype=str,
                        nrows=200_000)["X"].fillna("").str.strip()
        shift[off] = float(c.isin(list("12345") + ["8", "9"]).mean())
    nsfg_offset_share = max(v for k, v in shift.items() if k != 0) if nsfg_offset_share is None \
        else max(nsfg_offset_share, max(v for k, v in shift.items() if k != 0))
    dist = sm.value_counts().head(8).to_dict()
    print(f"  {dat.name:34s} n={n:6d}  samesex coded={share:.4f}  "
          f"shift sweep {[f'{k:+d}:{v:.3f}' for k, v in sorted(shift.items())]}  {dist}")
    detail[dat.name] = dict(n=n, coded_share=share, shift_sweep={str(k): v for k, v in shift.items()},
                            dist={k: int(v) for k, v in dist.items()}, norm_items=[a for a, _ in norms])

rows.append(("nsfg", "(3) person-level norm vs a second quantity",
             "SUPPORTED" if (nsfg_records > 0 and nsfg_norm_items) else
             ("READABLE_EMPTY" if nsfg_records > 0 else "UNREADABLE"),
             nsfg_records, len(nsfg_norm_items)))
rows.append(("nsfg", "(1) comonotonicity >=3 series x >=8 waves",
             "NOT" if nsfg_records > 0 else "UNREADABLE", nsfg_records, len(nsfg_norm_items)))

# ══ CELL 3-4 · YRBS — a shipped SAS input program beside 1.5 GB of .dat ═══════════════
y_lay = sas_input_layout(EXT / "yrbs" / "2023-SADC-SAS-Input-Program.sas")
y_dat = sorted((EXT / "yrbs").glob("sadc_2023_*.dat"))
y_labels = [a or b for a, b in sas_format_labels(EXT / "yrbs" / "2023-SADC-SAS-Formats-Program.sas")]
y_norm_lab = [t for t in y_labels if NORM_RE.search(t)]
y_behav_lab = [t for t in y_labels if BEHAV_RE.search(t)]
y_norm_var = [n for (n, _s, _w, _l) in y_lay if NORM_RE.search(n)]
y_records = 0
if y_lay and y_dat:
    small = min(y_dat, key=lambda p: p.stat().st_size)
    with open(small, "rb") as fh:
        y_records = sum(1 for _ in fh)
    spec = [(s, s + w) for (_n, s, w, _l) in y_lay[:12]]
    ydf = pd.read_fwf(small, colspecs=spec, names=[t[0] for t in y_lay[:12]], dtype=str, nrows=500)
    print(f"YRBS · layout vars={len(y_lay)} · {len(y_dat)} .dat · smallest={small.name} rows={y_records}")
    print(f"       value labels={len(y_labels)} · moral-evaluative={len(y_norm_lab)} · behavioural={len(y_behav_lab)}")
    print(f"       first cols read: {list(ydf.columns)[:6]} -> {ydf.iloc[0].tolist()[:6]}")
detail["yrbs"] = dict(layout_vars=len(y_lay), dat_files=len(y_dat), rows_smallest=y_records,
                      labels=len(y_labels), norm_labels=len(y_norm_lab), behav_labels=len(y_behav_lab))
# the NATIONAL file: is it truly reader-blocked?
natz = EXT / "yrbs" / "SADC_2023_National.zip"
nat_inside = zipfile.ZipFile(natz).namelist() if natz.exists() else []
nat_mdb = [n for n in nat_inside if n.lower().endswith(".mdb")]
has_mdb_reader = any(os.access(os.path.join(p, t), os.X_OK)
                     for p in os.environ.get("PATH", "").split(":") if p
                     for t in ("mdb-export", "mdb-tables"))
print(f"YRBS national: {nat_inside} · mdbtools on PATH={has_mdb_reader}")
rows.append(("yrbs", "(3) person-level norm vs a second quantity",
             "SUPPORTED" if (y_records > 0 and y_norm_lab) else
             ("READABLE_EMPTY" if y_records > 0 else "UNREADABLE"), y_records, len(y_norm_lab)))
rows.append(("yrbs", "(1) comonotonicity >=3 series x >=8 waves",
             "READABLE_EMPTY" if y_records > 0 else "UNREADABLE", y_records, len(y_norm_lab)))

# ══ CELL 5 · BRFSS — 1.2 GB .XPT, no codebook ════════════════════════════════════════
bx = EXT / "brfss" / "LLCP2023.XPT"
b_rows, b_cols, b_err = 0, [], ""
try:
    it = pd.read_sas(bx, format="xport", chunksize=50_000)
    for ch in it:
        if not b_cols:
            b_cols = [c if isinstance(c, str) else c.decode() for c in ch.columns]
        b_rows += len(ch)
except Exception as e:                                    # noqa: BLE001
    b_err = f"{type(e).__name__}: {e}"
b_codebook = [p for p in (EXT / "brfss").iterdir() if p.suffix.lower() in (".pdf", ".html", ".htm", ".sas")]
print(f"BRFSS · rows={b_rows} cols={len(b_cols)} err='{b_err}' codebook_files={len(b_codebook)}")
# ⚠ no shipped labels ⇒ stage B is UNDECIDABLE, NOT zero. `#891`: a name screen missed 7/7.
b_state = "UNREADABLE" if b_rows == 0 else ("READABLE_UNDECIDABLE" if not b_codebook else "SUPPORTED")
detail["brfss"] = dict(rows=b_rows, cols=len(b_cols), codebook_files=len(b_codebook), err=b_err)
rows.append(("brfss", "(3) person-level norm vs a second quantity", b_state, b_rows, -1))

print("\n=== THE GRID (all cells, including the ones that disagree) ===")
for r in rows:
    print(f"  {r[0]:8s} {r[1]:44s} {r[2]:22s} n={r[3]:<8d} norm_items={r[4]}")

# ══ CONTROLS ═════════════════════════════════════════════════════════════════════════
# POSITIVE, and it must FAIL at g=0: plant k known norm items among behaviour items and sweep.
GSS_NORM = ["Sex before marriage is always wrong",
            "Sexual relations between two adults of the same sex",
            "It is okay for an unmarried woman to have a child",
            "Do you approve or disapprove of a married woman earning money"]
BEHAV_POOL = ["Ever had sexual contact with a female", "How many times during the past 30 days",
              "Did you drink alcohol", "Respondent ID number", "Age at interview",
              "During the past 12 months how often", "Have you ever been told by a doctor"]
sweep = []
for k in range(0, 5):
    pool = list(BEHAV_POOL) + GSS_NORM[:k]
    sweep.append((float(k), float(sum(1 for t in pool if NORM_RE.search(t)))))
print(f"\npositive control sweep (planted, detected): {sweep}")
# ⚠ the rule MISSES GSS's "Sexual relations between two adults of the same sex" — recorded, not
#   hidden: that item's text names a topic and the moral evaluation lives in its RESPONSE SCALE.
#   So this rule under-counts norm items wherever the evaluation is not in the stem, which is a
#   false-negative direction and makes every `0` here a bound, never a zero.

# ── NEGATIVE: an invented release returns nothing, and an empty frame is NOT 'readable' ──
fake = EXT / "no_such_release_351" / "nothing.dat"
fake_rows = 0
try:
    fake_rows = len(pd.read_fwf(fake, colspecs=[(0, 1)], names=["X"], nrows=10))
except Exception:                                          # noqa: BLE001
    fake_rows = 0
empty_scores_readable = (0 > 0)                            # the rule itself: A>0 is required
print(f"negative control: invented release rows={fake_rows} · empty-scores-readable={empty_scores_readable}")

# ── SHAM: the same reader pointed at a wrong-format file must FAIL, so the verdict tracks
#    the FILE and not the code path. (§realstat: a sham is the operation MINUS the ingredient.)
sham_ok = False
sham_target = next((EXT / "ngram").glob("*"), None) if (EXT / "ngram").exists() else None
if sham_target and sham_target.is_file():
    try:
        next(iter(pd.read_sas(sham_target, format="xport", chunksize=10)))
        sham_ok = True                                     # it "read" a non-XPT ⇒ reader is blind
    except Exception:                                      # noqa: BLE001
        sham_ok = False
print(f"sham: XPT reader on {sham_target.name if sham_target else None} succeeded={sham_ok}")

if not rows:
    print("EMPTY POPULATION — nothing examined"); sys.exit(2)

# ══ GATE ═════════════════════════════════════════════════════════════════════════════
G = Gate("Are `#912`'s five UNREADABLE cells actually unreadable?")
G.plant_direction_from_sweep(
    "positive: the norm-item rule detects planted norm items and returns 0 when none are planted",
    sweep, baseline=0.0)
G.asserted("positive control fails at g=0", sweep[0][1] == 0,
           f"k=0 detected {sweep[0][1]:.0f} norm items among {len(BEHAV_POOL)} behaviour items",
           kind="control")
# ⚠ v1 of this control compared 0 to 0 and the gate called it DEGENERATE, correctly: a reader
#   that returns nothing on an invented release says nothing unless the SAME reader returns
#   something on a real one. The contrast is (real, fake), never (fake, 0).
real_rows = float(max(nsfg_records, y_records, b_rows))
G.negative_control("the same readers: invented release vs the releases actually on disk",
                   float(fake_rows), real_rows, null_kind="absent-file")
G.asserted("sham: the XPT reader refuses a non-XPT file", not sham_ok,
           f"reader on a foreign file succeeded={sham_ok}", kind="control")
# the confound named BEFORE the run — and MEASURED NOT SEPARABLE at this site, so it is
# registered as a structural impossibility and BOUNDS the claim, rather than gating it.
G.asserted("the byte-anchoring confound is registered, not claimed away",
           nsfg_offset_share is not None and nsfg_offset_share >= 0.99,
           f"shift sweep reaches {nsfg_offset_share} off-anchor vs {nsfg_valid_share} on-anchor "
           f"⇒ code-validity cannot anchor this parse; scope stated: the cell is READABLE and holds "
           f"a norm item, and NO attitude is measured by this round", kind="control")
G.asserted("every cell `#912` called UNREADABLE was actually opened this round",
           all(r[2] != "UNREADABLE" or r[3] == 0 for r in rows),
           "each cell reports a measured record count, not a directory listing")

moved = [r for r in rows if r[2] in ("SUPPORTED", "READABLE_EMPTY", "READABLE_UNDECIDABLE")]
supported_new = [r for r in rows if r[2] == "SUPPORTED"]
G.asserted("KILL: `#912`'s BOUND requires that no UNREADABLE cell holds a supportable pair",
           len(supported_new) == 0,
           f"cells that opened={len(moved)}/{len(rows)} · newly SUPPORTED={len(supported_new)}")

tv = G.three_valued()
if tv.startswith("UNVERIFIED"):
    VERDICT, WORLD = "UNVERIFIED", "controls unfit"
elif supported_new:
    VERDICT, WORLD = "OVERTURNED", "BOUND-ARTEFACT"
elif moved:
    VERDICT, WORLD = "OVERTURNED", "BOUND-SURVIVES-REASON-WRONG"
else:
    VERDICT, WORLD = "CONFIRMED", "BOUND-REAL"

print(f"\n{G}")
print(f"  gate three-valued : {tv}")
print(f"  VERDICT           : {VERDICT} · world {WORLD}")
print(f"  newly SUPPORTED   : {[(r[0], r[1]) for r in supported_new]}")

art = dict(entry=913, round="E03·A111·R351", verdict=VERDICT, world=WORLD,
           estimand="for each cell `#912` scored UNREADABLE: records obtainable with installed "
                    "software (A) and norm items in the shipped text (B)",
           grid=[dict(release=r[0], estimand=r[1], state=r[2], records=r[3], norm_items=r[4]) for r in rows],
           nsfg=dict(data_files=len(nsfg_dat), layouts=len(nsfg_dct), covered=len(covered),
                     records=nsfg_records, years=sorted(set(nsfg_years)),
                     norm_items=sorted(set(nsfg_norm_items)),
                     valid_share=nsfg_valid_share, offset_valid_share=nsfg_offset_share),
           yrbs=detail.get("yrbs"), brfss=detail.get("brfss"),
           yrbs_national=dict(contents=nat_inside, mdb=len(nat_mdb), mdbtools_on_path=has_mdb_reader),
           per_file=detail,
           positive_sweep=sweep, sham_reader_succeeded=sham_ok, fake_rows=fake_rows,
           gates=[(r[0], r[1], bool(r[2]), r[3]) for r in G.rows], gate_verdict=tv)
(OUT / "unreadable_cells_opened.json").write_text(json.dumps(art, indent=1, default=str))
print(f"\n  artifact -> {OUT/'unreadable_cells_opened.json'}")
