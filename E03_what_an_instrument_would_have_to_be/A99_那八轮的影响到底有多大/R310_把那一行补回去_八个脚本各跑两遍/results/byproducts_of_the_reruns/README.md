# Byproducts of R310's re-runs — NOT the original rounds' artifacts

`#871` re-ran five of `#869`'s eight scripts (as-committed, with the missing line restored, and with
one real item dropped). Those scripts write CSVs into **their own** `results/` directories, and those
CSVs had never been committed.

So after the re-runs, six untracked CSVs were sitting inside round directories whose ledger entries
date from far earlier. **Left there they would read as those rounds' original artifacts, which they
are not — they were produced today.** That is the "a rebuild launders staleness into currency"
failure, so they were moved here instead of being left in place or deleted (`L81`: never `rm`).

The filename encodes the path they were written to, with `/` replaced by `__`.

⚠ **They are the output of the AS-COMMITTED run** (the first of the three runs per script), so they
are faithful to the committed code — but they carry today's date, not the round's, and nothing in
this project should cite them as the original evidence.
