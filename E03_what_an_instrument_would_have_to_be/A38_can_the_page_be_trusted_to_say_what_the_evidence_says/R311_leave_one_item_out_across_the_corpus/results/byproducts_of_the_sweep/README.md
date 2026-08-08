# Byproducts of R311's corpus sweep — NOT the original rounds' artifacts

`#872` ran 35 scripts seven times each (as-committed twice, a no-op injection, and five single-item
drops). Those scripts write CSVs into **their own** `results/` directories, and those CSVs had never
been committed. After the sweep 35 untracked files were sitting inside round directories whose
ledger entries date from far earlier.

**Left in place, today's output would read as those rounds' original evidence. It is not.**
That is the "a rebuild launders staleness into currency" failure, so they were moved here rather
than left or deleted (`L81`: never `rm`).

The filename encodes the path they were written to, `/` replaced by `__`.

## The part worth keeping

This is the SECOND time (after `#871`). Last time the lesson was recorded as *"re-running an old
round is a write, not a read — check `git status` after"*. **I did check, and it happened again**,
because checking is a habit and the write is automatic. **A rule that depends on remembering to look
is one distraction from failing**; the durable form is for the sweep itself to relocate its
byproducts, which is what `#872`③ registers.
