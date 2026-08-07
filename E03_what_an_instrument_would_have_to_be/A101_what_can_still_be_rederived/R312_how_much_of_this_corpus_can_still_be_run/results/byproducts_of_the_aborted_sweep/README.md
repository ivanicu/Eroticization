# Byproducts of the ABORTED sweep of `#874` — not evidence of anything

`#874` set out to run all 835 round scripts to measure how much of this corpus can still be
re-derived. **It wedged after ~2.5 hours with no output and was killed.** These 179 CSVs are what the
scripts it managed to run wrote into their own `results/` directories before that.

They are moved here rather than left in place (`L81`: never `rm`) for the same reason as `#871` and
`#872`: **left where they were written, today's output would read as those rounds' original
evidence.**

⚠ **They are NOT a partial result.** The sweep produced no classification, so nothing here supports
any statement about re-derivability. The diagnosis of why it wedged is in ledger entry `#874`.
