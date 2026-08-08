r"""#875 · are these 300-odd failures one error or three hundred? — the question `§0.2` forces

Ivan's standing demand of any report carrying a large failure count: **「你是不是都在犯同一种错误?
这个错误是什么?」** A count is not an answer. A ledger of 330 dead scripts is worth nothing until it
says whether they are 330 independent losses or **one loss, propagated**.

**THE HYPOTHESIS THIS TESTS** (written before it was run, and it has a losing branch):
   **CASCADE** — one broken name kills ~80 *producer* scripts; those producers were the only thing
   that wrote certain derived files; every *consumer* of those files then dies as `MISSING-INPUT`.
   ⇒ the corpus has **one** re-derivability defect wearing 330 costumes.
   **INDEPENDENT** — the missing inputs are unrelated to the dead producers: raw data that was
   never committed, one-off scratch files, external downloads.
   ⇒ the corpus has **many** small losses, and there is no single repair.

**PREDICTION MATRIX** (coarse, the shape is the point):

   | world       | now  | most missing files are produced by a DEAD script | by a LIVE script | by nothing |
   | CASCADE     | 0.45 | **0.85**                                         | 0.10             | 0.05       |
   | INDEPENDENT | 0.45 | 0.05                                             | 0.15             | **0.80**   |
   | BUILD-ORDER | 0.10 | 0.10                                             | **0.75**         | 0.15       |

**BUILD-ORDER** is the third world and it is the one I would find *most* welcome, so it is named
explicitly rather than allowed to win by default: the producer is alive and simply was never re-run.

`G1` **ESTIMAND**: of the scripts classified `MISSING-INPUT` by the sweep, the share whose missing
file is written by another script in this corpus, split by **whether that producer itself runs**.
Population = the sweep's own `MISSING-INPUT` set, read from its artifact. Instrument = a static scan
for write-calls naming that file. Baseline = a file no script writes.

⚠ **THE INSTRUMENT IS THE STATIC SCAN THAT ALREADY FAILED ONCE.** `#874` recorded that extracting
path-shaped strings and testing existence reported 262/835 "missing" and the top hits were
directory-name fragments and a file that existed. **That failure was about READS.** This scan is
about **WRITES**, which is a different and narrower shape (`to_csv(...)`, `np.save(...)`,
`OUT / "..."`), but the lesson transfers: it gets a **two-sided control** before it is believed —
a file known to be written must be found, and an invented filename must not be.
"""
import json
import pathlib
import re
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
ART = OUT / "how_much_can_still_be_rederived.json"
A = json.load(open(ART))
FINAL = A["final"]

scripts = [f for f in sorted(ROOT.glob("E0*/**/*.py"))
           if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}" and "R314_" not in str(f)]
WRITE = re.compile(r"""(?:to_csv|to_parquet|np\.save|np\.savez|savefig|write_text|write_bytes)\s*\(\s*[^)]*?['"]([^'"]+)['"]"""
                   r"""|OUT\s*/\s*['"]([^'"]+)['"]"""
                   r"""|DERIVED\s*/\s*['"]([^'"]+)['"]""")
writers = {}
for f in scripts:
    t = f.read_text(errors="replace")
    for m in WRITE.finditer(t):
        nm = pathlib.Path(next(g for g in m.groups() if g)).name
        writers.setdefault(nm, set()).add(str(f.relative_to(ROOT)))

print("=== (1) CONTROL on the write-scan — two-sided, before any count is believed ===")
known = None
for nm, ws in writers.items():
    if nm.endswith(".csv") and (ROOT / "data" / "derived" / nm).exists():
        known = nm
        break
POS = known is not None
NEG = "__no_such_output_file__.csv" not in writers
print(f"  positive: a file that exists in data/derived AND is named by a write-call: {known!r} "
      f"-> **{'PASS' if POS else 'FAIL'}** (written by {len(writers.get(known, [])) if known else 0} script(s))")
print(f"  negative: an invented filename must not be found -> **{'PASS' if NEG else 'FAIL'}**")
print(f"  the scan knows {len(writers)} distinct output filenames across {len(scripts)} scripts")

print("\n=== (2) THE CASCADE TEST ===")
miss = {s: v["detail"] for s, v in FINAL.items() if v["cls"] == "MISSING-INPUT" and v.get("detail")}
print(f"  MISSING-INPUT with an identified path: **{len(miss)}**")
buckets = Counter()
by_producer_state = Counter()
detail_rows = []
for s, d in miss.items():
    nm = pathlib.Path(d).name
    ws = writers.get(nm, set())
    if not ws:
        buckets["written by NO script in this corpus"] += 1
        detail_rows.append((s, nm, "none", "-"))
        continue
    states = {FINAL.get(w, {}).get("cls", "NOT-IN-POPULATION") for w in ws}
    if states == {"OK"}:
        buckets["written by a script that RUNS (build order)"] += 1
        st = "OK"
    elif "OK" in states:
        buckets["written by several, at least one of which RUNS"] += 1
        st = "mixed"
    else:
        buckets["written ONLY by scripts that are themselves DEAD"] += 1
        st = "/".join(sorted(states))
    by_producer_state[st] += 1
    detail_rows.append((s, nm, ";".join(sorted(ws))[:80], st))
for k, n in buckets.most_common():
    print(f"  {n:4d}  ({n/max(1,len(miss)):5.1%})  {k}")
print("\n  producer states among the dead-producer cases:")
for k, n in by_producer_state.most_common():
    print(f"     {n:4d}  {k}")

dead_prod = buckets["written ONLY by scripts that are themselves DEAD"]
none_prod = buckets["written by NO script in this corpus"]
live_prod = (buckets["written by a script that RUNS (build order)"]
             + buckets["written by several, at least one of which RUNS"])
tot = max(1, len(miss))

print("\n" + "=" * 100)
if not (POS and NEG):
    V = "**UNVERIFIED — the write-scan's own two-sided control did not pass, so no share is admissible.**"
else:
    lead = max([(dead_prod, "CASCADE"), (none_prod, "INDEPENDENT"), (live_prod, "BUILD-ORDER")])
    V = (f"**{lead[1]}** leads: dead-producer {dead_prod} ({dead_prod/tot:.0%}) · "
         f"no-producer {none_prod} ({none_prod/tot:.0%}) · live-producer {live_prod} "
         f"({live_prod/tot:.0%}), over {len(miss)} `MISSING-INPUT` scripts.\n"
         f"  ⇒ the answer to *are these all the same error* is **"
         + ("yes — one broken name upstream, and the rest is propagation**"
            if lead[1] == "CASCADE" else
            "no — the missing inputs were never produced here at all**"
            if lead[1] == "INDEPENDENT" else
            "no — the producers are alive and were simply never re-run; this is a build order**")
         + ".")
print(V)
print("\n⚠ **Scope**: this is a STATIC write-scan, so a file written through a variable filename is "
      "invisible to it and lands in `written by NO script` — the flattering direction for "
      "INDEPENDENT and the conservative one for CASCADE. The share is a lower bound on cascade.")

json.dump(dict(writers_known=len(writers), controls=dict(positive=POS, positive_file=known,
                                                         negative=NEG),
               n_missing=len(miss), buckets=dict(buckets),
               producer_states=dict(by_producer_state),
               rows=[dict(script=a, missing=b, writers=c, producer_state=d) for a, b, c, d in detail_rows],
               verdict=V),
          open(OUT / "is_it_one_cause_or_many.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  artifact -> {OUT/'is_it_one_cause_or_many.json'}")
