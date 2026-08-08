r"""#875 · the sweep's own citation regex was too tight — recount, and say which direction it erred

The sweep decides whether a dead script sits behind **a claim the README still makes** by walking
`script -> round -> ledger entry -> "Entry N" in a README`. The last step used `Entry (\d+)`.

**That rule cannot see a continuation.** The pages cite runs of entries — `(Entry 529 · 530)`,
`(Entry 640 · 641 · 642)` — and the tight rule keeps only the first number of each run. Measured
against the two pages: **177 cited entries by the tight rule, 205 with continuations, 28 missed
(13.7%).**

**The direction is the point.** A missed citation makes a round look *not live*, so the sweep's
"dead script behind a live claim" figure is an **UNDERCOUNT** — it errs toward the comfortable
answer. That is the direction one never notices, because the number one wanted came out.

This script recomputes the figure from the sweep's own artifact with the loose rule, and reports
**both**, because replacing a number quietly is how a corrected instrument stops being auditable.

`G1` **ESTIMAND**: the number of scripts that do not run whose round is cited in a README, under two
citation rules. Population = the sweep's own non-OK set, read from its artifact.
**POSITIVE CONTROL, known answer:** `Entry 530` is cited on the page as the second half of
`(Entry 529 · 530)`; the loose rule must find it and the tight rule must not. A recount that
returned the same number under both rules would mean the loose rule is not looser — a check that
cannot fail.
**NEGATIVE CONTROL:** an entry number that appears nowhere (`Entry 99999`) must be absent under both.
"""
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
A = json.load(open(OUT / "how_much_can_still_be_rederived.json"))
FINAL = A["final"]

led = (ROOT / "RETRACTIONS.md").read_text(errors="replace")
pages = (ROOT / "README.md").read_text(errors="replace") + (ROOT / "README_zh.md").read_text(errors="replace")

ear = {}
for m in re.finditer(r'^## Entr(?:y|ies) (\d+) · `([^`]+)`', led, re.M):
    for r in re.findall(r'R(\d+)', m.group(2)):
        ear.setdefault(f"R{int(r):03d}", set()).add(int(m.group(1)))

TIGHT = {int(n) for n in re.findall(r'Entry (\d+)', pages)}
LOOSE = set(TIGHT)
for m in re.finditer(r'Entry\s+(\d+(?:\s*[·,、]\s*\d+)*)', pages):
    LOOSE |= {int(n) for n in re.findall(r'\d+', m.group(1))}

print("=== (1) CONTROLS on the citation rule itself ===")
POS = (530 in LOOSE) and (530 not in TIGHT)
NEG = (99999 not in LOOSE) and (99999 not in TIGHT)
print(f"  positive: Entry 530 is the second half of `(Entry 529 · 530)` — loose sees it "
      f"({530 in LOOSE}), tight does not ({530 not in TIGHT}) -> **{'PASS' if POS else 'FAIL'}**")
print(f"  negative: Entry 99999 appears on neither page -> **{'PASS' if NEG else 'FAIL'}**")
print(f"  cited entries: tight **{len(TIGHT)}** · loose **{len(LOOSE)}** · missed "
      f"**{len(LOOSE - TIGHT)}**")


def round_of(p):
    for part in p.split("/"):
        m = re.match(r'^R(\d+)_', part)
        if m:
            return f"R{int(m.group(1)):03d}"
    return None


def live(p, cited):
    r = round_of(p)
    return bool(r and (ear.get(r, set()) & cited))


dead = [s for s, v in FINAL.items() if v["cls"] not in ("OK", "NOT-REACHED")]
perm = set(A.get("permanent", {}))
rows = {}
for label, cited in (("tight", TIGHT), ("loose", LOOSE)):
    rows[label] = dict(
        dead_live=sum(1 for s in dead if live(s, cited)),
        perm_live=sum(1 for s in perm if live(s, cited)),
        rounds_live=len({round_of(s) for s in dead if live(s, cited)}),
    )
print("\n=== (2) RECOUNT ===")
print(f"  {'rule':6s} {'dead & cited':>13s} {'PERMANENT & cited':>19s} {'distinct rounds':>16s}")
for label in ("tight", "loose"):
    r = rows[label]
    print(f"  {label:6s} {r['dead_live']:13d} {r['perm_live']:19d} {r['rounds_live']:16d}")
d = rows["loose"]["dead_live"] - rows["tight"]["dead_live"]

print("\n" + "=" * 100)
if not (POS and NEG):
    V = "**UNVERIFIED — the citation rule's own two-sided control did not pass.**"
else:
    V = (f"**The sweep undercounted by {d} scripts** ({rows['tight']['dead_live']} -> "
         f"{rows['loose']['dead_live']}), across "
         f"{rows['loose']['rounds_live']} distinct rounds. **PERMANENT & cited: "
         f"{rows['tight']['perm_live']} -> {rows['loose']['perm_live']}.**\n"
         f"  Both numbers are reported because a corrected instrument that quietly replaces its "
         f"predecessor stops being auditable — and because the correction moved the figure in the "
         f"**unflattering** direction, which is the direction a tight regex never moves it by "
         f"accident.")
print(V)
print("\n⚠ **Still an undercount**: a page can point at a finding without writing `Entry N` at all "
      "(prose reference, a number quoted with no citation). This rule can only see explicit "
      "citations, so `not cited` never means `not relied upon`.")

json.dump(dict(tight_cited=len(TIGHT), loose_cited=len(LOOSE), missed=sorted(LOOSE - TIGHT),
               rows=rows, controls=dict(positive=POS, negative=NEG), verdict=V),
          open(OUT / "recount_the_live_claims.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'recount_the_live_claims.json'}")
