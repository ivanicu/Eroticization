#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plant_baseline_gate — a plant's g=0 arm must sit on the null the round judges against.

Built at `#922` to pay `#920`(2), which named five instances of one family — *I describe what I
meant rather than measure what I did* — and observed that **the gate caught every one and my own
resolutions caught none**. `#921` then carried two repairs forward BY HAND, which is the evidence
that they belong in a tool rather than in my memory of the last round.

WHAT IT CHECKS, and it is the mechanically checkable core of that family:
  a positive control plants an effect of size g and sweeps. At g=0 there is no plant, so the
  statistic must land on the SAME baseline the round uses to judge the observed value. If it does
  not, the plant arm and the null arm were drawn from **two different worlds**, and the control
  could only ever have failed (or, worse, passed for the wrong reason).
  ⇒ measured witness, `#920` before its repair: sweep at g=0 was **+0.6894** while the null it was
  judged against was **-0.0015 +/- 0.0549** — a gap of ~12.6 spreads, because the plant went into
  the OBSERVED data while the baseline came from the PERMUTED world. Same shape as `#905`.

⚠ P6 PROXY LEDGER
  PROPERTY    the plant's zero arm and the null are the same world
  PROXY       |sweep(g=0) - null_center| / null_spread
  IMPLICATION one direction only: **a large ratio ⇒ they really are different worlds** (reliable).
              The converse does NOT hold: a small ratio is not proof of a shared world — a plant
              into the observed data can sit near a null by coincidence.
  WITNESS     `#920` pre-repair (above), reproduced as a unit test below.
  SAFE SIDE   report "different worlds" only. **Never certify that a control is sound.**

⚠ AND THE SECOND-ORDER HONESTY THIS TOOL OWES ITSELF (`#913`'s lesson, turned inward): the round
  artifacts have **no common schema for "the null this round judged against"** — five different
  shapes across 19 files. An artifact whose null cannot be located is scored **UNREADABLE, never
  PASS**, and the unreadable count is printed beside the failures. A gate that silently passes what
  it cannot read is the "empty population passes" failure wearing a tool's clothes.

Exit codes: 0 clean · 1 at least one artifact fails · 2 empty population or a control failed.
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
K = 2.0                      # gap tolerated, in units of the null's own spread
CUTOFF = 922                 # from this entry on, an unlocatable null BLOCKS; earlier ones are
                             # named but not blocked — the project's own precedent (`readme_gate`
                             # grandfathers pre-cutoff entries; `#605`: never freeze a worse count
                             # into the baseline, but never retro-block either).

# ⚠ WHAT A ROUND MUST PERSIST FOR THIS GATE TO SEE IT AT ALL. Measured at `#922`: of 19 artifacts
#   carrying a plant sweep, **9 record no locatable null**, so a post-hoc gate is blind on half the
#   corpus. That is the finding, and it corrects `#920`(2): the mechanism that caught the five
#   instances was the IN-ROUND `lib/gates.py`, which sees live variables. A retrospective tool
#   cannot substitute for it — it can only enforce that the round wrote down what it judged against.
REQUIRED = "a `[[g, value], ...]` sweep AND a scalar `null_median` + `null_sd` (or `null_p95`)"


def _find_sweep(a):
    """Return (key, [(g, value), ...]) for the first plant sweep, or None."""
    for k, v in a.items():
        if "sweep" not in k.lower() or not isinstance(v, list) or not v:
            continue
        pairs = []
        for item in v:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                try:
                    pairs.append((float(item[0]), float(item[1])))
                except (TypeError, ValueError):
                    pairs = []
                    break
        if pairs:
            return k, sorted(pairs)
    return None


def _find_null(a):
    """Return (center, spread, where) for the null this round judged against, or None.

    ⚠ Five schemas are in use across the corpus. They are enumerated rather than guessed, and
    anything not matching returns None so the artifact is scored UNREADABLE."""
    def _num(x):
        # ⚠ a "null" field is not always a NUMBER: measured, `null_p95` is a dict in part of the
        #   corpus (per-cell nulls). A gate that assumes the shape crashes on the first exception;
        #   a gate that coerces it invents a value. Neither — it is UNREADABLE.
        return float(x) if isinstance(x, (int, float)) and not isinstance(x, bool) else None

    def _pair(d, prefix=""):
        c = _num(d.get("null_median", d.get("null_modal")))
        s = _num(d.get("null_sd"))
        if c is not None and s is not None:
            return c, s, f"{prefix}null_median/null_sd"
        p95 = _num(d.get("null_p95"))
        if p95 is not None:
            return 0.0, abs(p95) / 1.96, f"{prefix}null_p95"
        return None

    got = _pair(a)
    if got:
        return got
    for key in ("null", "nulls"):
        d = a.get(key)
        if isinstance(d, dict):
            c = _num(d.get("id_median", d.get("median")))
            s = _num(d.get("id_sd", d.get("sd")))
            if c is not None and s is not None:
                return c, s, f"{key}.median/sd"
    return None


def audit(paths):
    fails, unreadable, ok = [], [], []
    for p in paths:
        try:
            a = json.load(open(p))
        except Exception:                                   # noqa: BLE001
            unreadable.append((p, "unparseable"))
            continue
        if not isinstance(a, dict):
            unreadable.append((p, "not an object"))
            continue
        sw = _find_sweep(a)
        if not sw:
            continue                                        # no plant sweep: out of scope
        nl = _find_null(a)
        if not nl:
            unreadable.append((p, f"sweep `{sw[0]}` present but NO LOCATABLE NULL"))
            continue
        key, pairs = sw
        g0, v0 = pairs[0]
        center, spread, where = nl
        if spread <= 0:
            unreadable.append((p, "null spread is zero or negative"))
            continue
        gap = abs(v0 - center) / spread
        row = dict(path=str(p), entry=a.get("entry"), sweep_key=key, g0=g0, value_at_g0=v0,
                   null_center=center, null_spread=spread, null_where=where, gap_in_spreads=gap)
        (fails if (g0 == 0 and gap > K) else ok).append(row)
    return fails, unreadable, ok


def _unit_tests():
    """Positive: the measured `#920` pre-repair numbers must be CAUGHT. Negative: its repaired
    numbers must NOT be. A gate that has never fired is silence, not a clean corpus."""
    import tempfile
    bad = dict(entry=9001, positive_sweep=[[0.0, 0.6894], [0.6, 0.7528]],
               null_median=-0.0015, null_sd=0.0549)
    good = dict(entry=9002, positive_sweep=[[0.0, 0.0200], [0.6, 0.1728]],
                null_median=-0.0015, null_sd=0.0549)
    blind = dict(entry=9003, positive_sweep=[[0.0, 0.5], [0.6, 0.9]])
    out = []
    with tempfile.TemporaryDirectory() as td:
        for name, obj in (("bad", bad), ("good", good), ("blind", blind)):
            q = pathlib.Path(td) / f"{name}.json"
            q.write_text(json.dumps(obj))
            out.append(q)
        f, u, o = audit(out)
    return (len(f) == 1 and f[0]["entry"] == 9001,
            len(o) == 1 and o[0]["entry"] == 9002,
            len(u) == 1 and u[0][1].endswith("NO LOCATABLE NULL"))


def main():
    print("=== plant-baseline gate (`#922`, pays `#920`(2)) ===")
    pc, nc, bc = _unit_tests()
    print(f"  control positive (`#920` pre-repair, gap ~12.6 spreads, must be CAUGHT): "
          f"{'PASS' if pc else 'FAIL'}")
    print(f"  control negative (`#920` post-repair, must NOT be caught):               "
          f"{'PASS' if nc else 'FAIL'}")
    print(f"  control blind    (a sweep with no locatable null -> UNREADABLE, not PASS): "
          f"{'PASS' if bc else 'FAIL'}")
    if not (pc and nc and bc):
        print("  ⛔ a control failed — the gate's verdict on the corpus is INADMISSIBLE")
        return 2

    paths = sorted(ROOT.glob("E0*/A*/R*/results/*.json"))
    if not paths:
        print("  ⛔ EMPTY POPULATION — exit 2, never 0")
        return 2
    fails, unreadable, ok = audit(paths)
    scanned = len(fails) + len(unreadable) + len(ok)
    print(f"\n  artifacts on disk {len(paths)} · carrying a plant sweep {scanned}")
    print(f"  checked {len(ok) + len(fails)} · ⚠ UNREADABLE {len(unreadable)} · FAIL {len(fails)}")
    for r in sorted(fails, key=lambda r: -r["gap_in_spreads"]):
        print(f"    ⛔ #{r['entry']}  g=0 at {r['value_at_g0']:+.4f} vs null "
              f"{r['null_center']:+.4f} +/- {r['null_spread']:.4f}  = {r['gap_in_spreads']:.1f} spreads"
              f"  [{r['null_where']}]")
    blocking_unreadable = []
    for p, why in unreadable:
        try:
            e = json.load(open(p)).get("entry") or 0
        except Exception:                                   # noqa: BLE001
            e = 0
        tag = "⛔ BLOCKS" if e >= CUTOFF else "⚠ named, pre-cutoff"
        if e >= CUTOFF:
            blocking_unreadable.append((p, e))
        print(f"    {tag} {pathlib.Path(p).parent.parent.name}: {why}")
    if blocking_unreadable:
        print(f"\n  ⛔ {len(blocking_unreadable)} artifact(s) at or after #{CUTOFF} persist a sweep "
              f"with no locatable null. Required: {REQUIRED}")
    print(f"\n  ⚠ Scope: only rounds that PERSIST a sweep and a locatable null are visible. "
          f"{len(unreadable)} are UNREADABLE and are NOT cleared. And a small gap is not proof the "
          f"plant and the null share a world — this gate reports difference only (P6 safe side).")
    print(f"  ⚠ From #{CUTOFF} on, an unlocatable null BLOCKS. Required: {REQUIRED}")
    return 1 if (fails or blocking_unreadable) else 0


if __name__ == "__main__":
    sys.exit(main())
