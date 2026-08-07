r"""`#875`② — a registry key that no longer exists is a script that no longer runs. Check it.

**WHY THIS EXISTS, measured rather than imagined.** `lib/rounds.py` is a registry: it exists so that
a path is written once and referred to by key, instead of being copied into every caller. It solved
that. Then it created the failure one level up: **nothing checks that a caller's KEY still exists.**

Twice in three days a restructure rebuilt the key vocabulary and migrated no callers:
  · `5807d87` flattened `src/00..51.py` into one directory per round;
  · `4819b9b` flattened again and renamed `PATHS`→`ROUNDS`, `round_path`→`path`, turning 52
    numbered filenames into 635 stems with **zero keys in common**.
`#875` measured the result: **79 of 836 round scripts dead at one name, the only broken name in the
whole shared library**, and it went unnoticed for a day because **nothing in this project re-runs an
old script**. A full corpus sweep costs ~70 minutes. This gate costs one pass over the sources.

**WHAT IT CHECKS.** For every literal `round_path('k')` / `path('k')` in `E01`/`E02`/`E03`/`lib`/
`tools`, the key `k` must be present in the registry that accessor reads, and the path it maps to
must exist on disk.

⚠ **A grep is a measuring instrument** (`realstat`), so this one carries its own two-sided control
and **exits 2 on an empty population** — a gate that reports success having examined nothing is the
`empty population passes` failure, and it is why the control runs first.
⚠ **Scope, stated rather than implied**: only *literal* keys are visible. A caller that builds its
key at runtime is **UNSEEN, which is not CLEARED** — `#875` measured that as 30 of 80 importers, so
this gate covers the majority and says so rather than claiming coverage it does not have.
"""
import ast
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib import rounds  # noqa: E402

# ⚠⚠ **PARSED, NOT MATCHED — and the first version of this gate proved why in one run.**
# Written with a regex it reported **2 broken references, both inside its own source**, because its
# docstring shows `round_path('k')` as an illustration. That is the FOURTH instance in one session
# of a text scan counting a MENTION as a USE (a defect-token gate flagged the ledger for quoting the
# string it was retracting; a backtick pair spanning two code spans invented a token; an import scan
# counted a docstring quoting the broken import). An `ast.Call` cannot be written in a docstring, so
# parsing removes the whole class rather than excluding this file and hoping.
ACCESSORS = {"round_path": "LEGACY_PATHS", "path": "ROUNDS"}


def literal_keys(src):
    """(accessor, key) for every call `accessor("literal")` — a CALL, never a string that looks like one."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return []
    out = []
    for n in ast.walk(tree):
        if (isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                and n.func.id in ACCESSORS and len(n.args) == 1
                and isinstance(n.args[0], ast.Constant) and isinstance(n.args[0].value, str)):
            out.append((n.func.id, n.args[0].value))
    return out
SOURCES = [f for d in ("E01_sexual_as_a_value_not_a_category",
                       "E02_condemnation_is_not_rarity",
                       "E03_what_an_instrument_would_have_to_be", "lib", "tools")
           for f in sorted((ROOT / d).rglob("*.py"))
           if "__pycache__" not in str(f) and "/_archive/" not in f"/{f}"]


def scan():
    out = []
    for f in SOURCES:
        for acc, k in sorted(set(literal_keys(f.read_text(errors="replace")))):
            table = ACCESSORS[acc]
            reg = getattr(rounds, table, None)
            if reg is None:
                out.append((str(f.relative_to(ROOT)), acc, k,
                            f"the registry `{table}` this accessor reads does not exist"))
            elif k not in reg:
                out.append((str(f.relative_to(ROOT)), acc, k,
                            f"key absent from `{table}` ({len(reg)} keys)"))
            elif not (ROOT / reg[k]).exists():
                out.append((str(f.relative_to(ROOT)), acc, k,
                            f"key maps to `{reg[k]}`, which does not exist on disk"))
    return out


def controls():
    """Two-sided, and it must be able to fail: a known-good key passes, an invented one is caught."""
    good = next(iter(rounds.LEGACY_PATHS)), "LEGACY_PATHS"
    pos = good[0] in rounds.LEGACY_PATHS and (ROOT / rounds.LEGACY_PATHS[good[0]]).exists()
    neg = "__no_such_registry_key__" not in rounds.LEGACY_PATHS
    return pos, neg, good[0]


if __name__ == "__main__":
    pos, neg, gk = controls()
    n_lit = sum(len(set(literal_keys(f.read_text(errors="replace")))) for f in SOURCES)
    print("=== registry-key gate (`#875`②) ===")
    print(f"  sources scanned: {len(SOURCES)} · literal keys seen: **{n_lit}**")
    print(f"  control positive (a real key resolves and its file exists): {gk!r} -> "
          f"**{'PASS' if pos else 'FAIL'}**")
    print(f"  control negative (an invented key is absent): **{'PASS' if neg else 'FAIL'}**")
    if not (pos and neg):
        print("  => **UNVERIFIED** — the gate's own control did not pass; it has no standing to rule")
        sys.exit(2)
    if n_lit == 0:
        print("  => **EXIT 2** — an empty population must never be counted as a pass")
        sys.exit(2)
    bad = scan()
    print(f"\n  registries: LEGACY_PATHS {len(rounds.LEGACY_PATHS)} · ROUNDS {len(rounds.ROUNDS)}")
    print(f"  **broken literal references: {len(bad)}**")
    for f, acc, k, why in bad[:40]:
        print(f"     {acc}({k!r}) — {why}\n        {f}")
    print("\n  ⚠ **Scope**: only LITERAL keys are visible — a caller that builds its key at "
          "runtime is **UNSEEN, not CLEARED** (`#875` measured that as 30 of 80 importers). "
          "This gate parses rather than matches, so a key quoted in prose is correctly invisible.")
    sys.exit(1 if bad else 0)
