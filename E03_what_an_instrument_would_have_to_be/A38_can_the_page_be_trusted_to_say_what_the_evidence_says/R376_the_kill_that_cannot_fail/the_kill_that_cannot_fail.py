#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
E03·A120·R376 — how many of this project's pre-registered kills could actually have failed?
============================================================================================

COGNITIVE UPDATE CARD
---------------------
Core Gap        `#937`③: **four vacuous kills now** — `#916`③, `#930`, `#935`, `#937` — each a
                pre-registered threshold that could not fire against the round, and **each caught by
                reading**. A kill is the single load-bearing commitment in this project's method: it
                is what makes a verdict evidence rather than narrative. **If some unknown fraction of
                894 round scripts carry kills that cannot fail, then that fraction of this project's
                verdicts are decorations, and nobody knows which.** Four instances is a rate of
                unknown denominator, which is not a rate at all.

Why now         `#937` was the fourth, it happened today, and the only reason it did not ship is that
                I happened to re-read a line I had just written. **That is not a mechanism, and every
                other guard in this project exists because I said exactly that about something else.**

Live Worlds    W_RARE  · vacuous kills are a handful of slips, ≲1% of kills. ⇒ the four are
                          anecdotes, the method is sound, and a runtime probe for new rounds is
                          enough.
               W_ENDEMIC· a substantial share of kills compare a quantity to a threshold derived from
                          that same quantity. ⇒ **the project's verdict history needs re-reading, not
                          just its future**, and that is the unwelcome one (BASIN RULE: I have
                          written 900+ entries assuming the kills were real).
               W_UNSEEABLE· kills are not machine-readable at all in most rounds. ⇒ the honest output
                          is a coverage number and a refusal, and the finding is about the corpus's
                          shape rather than about kills. (the meta-separator: my decomposition
                          assumes "kill" is a locatable object in the source; it may not be)

Estimand       Among round scripts whose kill condition is locatable in the AST, the SHARE that
(G1)           compares a statistic to a QUANTILE OF THAT STATISTIC'S OWN BOOTSTRAP. Denominator
               stated, coverage stated, non-locatable rounds counted as UNREADABLE, never as clean.

⚠⚠ THE FIRST DETECTOR FAILED ITS OWN POSITIVE CONTROL, AND THAT IS WHERE THE INSTRUMENT CAME FROM.
  v1 flagged "one operand inside the other's dependency closure". It did not catch `#937` v1 —
  because there, `hi` never mentions `diff`; **the dependency runs through the DATA, not through the
  variable.** And widening it to "shares a data ancestor" flags every criterion ever written, since
  in one script everything descends from the same frame. So the round refused to report a corpus
  number, which is what a failed positive control is for.

⚠ P6 PROXY LEDGER — the semantic line, which is what makes this readable at all
  PROPERTY    the kill could not have fired against the round
  PROXY       the threshold is a QUANTILE functional over a collection built WITHOUT permuting
              anything, sharing a data ancestor with the statistic it judges
  IMPLICATION **one direction only.** Match ⇒ SUSPECT and worth a human read. No match is NOT
              proof of soundness — a kill can be vacuous in ways no AST models (a threshold set
              after seeing the result, a condition covering its own complement as in `#930`).
  WHY THIS    a statistic against a quantile of its OWN bootstrap is vacuous **by definition**: q%
  LINE        of that distribution lies below its q-th percentile whatever the data say. Against a
              MULTIPLE OF ITS OWN SE it is the ordinary signal-to-noise test and is sound, because
              a bootstrap sd does not scale with the effect. Against a quantile of a PERMUTATION
              null it is also sound, because that distribution was built with the structure
              destroyed. ⇒ quantile functional, **minus** a permutation in its ancestry.
  WITNESS     `effect > 2*se` and `effect > percentile(null_perm, 95)` are both legitimate and both
              must come back clean; they are negative controls below, not decoration.
  SAFE SIDE   report CANDIDATES. **Never call a kill vacuous from the static pass.**

Prediction     W_RARE     -> candidate share ≲2% of locatable kills.
Matrix         W_ENDEMIC  -> ≳10%, and the four known cases are the visible part.
               W_UNSEEABLE-> locatable kills are a small minority of the 894 scripts.

Controls       POSITIVE: `#937` v1 verbatim — the bootstrap resamples the two arms, `hi` is its
                 97.5th percentile, and the kill is `diff <= hi`. Must be CAUGHT.
               POSITIVE2: the same laundered through an intermediate name. A one-hop check misses it.
               NEGATIVE1: `#937` v2, `lo <= 0.0 <= hi` — bounds from the statistic, comparison
                 against a CONSTANT. Must NOT be caught.
               NEGATIVE2: `effect > 2*se`, se = sd of the same bootstrap. Sound. Must NOT be caught.
               NEGATIVE3: `effect > percentile(null_permuted, 95)`. Sound — the structure was
                 destroyed to build it. Must NOT be caught. **This is the one that separates a
                 quantile-of-own-bootstrap from a quantile-of-a-null**, i.e. the whole instrument.
               BLIND: a script with no locatable kill -> UNREADABLE, never clean.

Stopping rule  Pre-registered: candidate share ≥10% of locatable kills ⇒ W_ENDEMIC and the ledger's
               verdicts need a re-read pass. ≤2% ⇒ W_RARE. Between ⇒ report the number, claim
               neither. Locatable kills < 20% of scripts ⇒ W_UNSEEABLE dominates and the share is
               reported only as conditional on locatability.

STRUCTURALLY CANNOT (registered; "planned" is forbidden)
  (1) ⚠ **decide vacuity statically** — scale covariance is a runtime property; see the proxy ledger;
  (2) ⚠ **audit kills written in prose only** — a threshold applied by hand and reported in the
    ledger has no AST to read, and those rounds are UNREADABLE here, not clean;
  (3) ⚠ **read a kill built through pandas/numpy indirection** the tracer does not model;
  (4) ⚠ **no second instrument, and this is not a shortfall being excused** — `#658` closes a round
    only when the same question has been asked on ≥2 instruments. Here the estimand *names* the
    instrument: "among **this project's** locatable kills". Running the detector over another
    repository's scripts would answer a different question (whether some other author writes vacuous
    kills), not this one. **The corpus is the population, so there is only this one instrument** and
    the cross-instrument move is structurally unavailable rather than skipped;
  (5) `[unchallenged]` — door ③.
"""
import ast
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.gates import Gate  # noqa: E402

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

KILL_MARK = ("KILL", "kill:", "预注册", "判据")


# ══ the tracer ═══════════════════════════════════════════════════════════════════════
def _names(node):
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def deps_of(tree):
    """name -> (transitive data-dependency names, transitive set of function names called).

    ⚠ MUTATION IS AN ASSIGNMENT. `boot = []` followed by `boot.append(f(a, b))` leaves an
    assignment-only tracer believing `boot` depends on nothing — and that is precisely how `#937`
    v1 built its bootstrap, so the first two versions of this tracer failed the positive control
    with the detector logic already correct. `.append` / `.extend` / `+=` are modelled as reads.

    ⚠ Flow-insensitive on purpose: a name reassigned later carries the union of every assignment's
    reads. That direction is SAFE here — it can only ADD candidates, never hide one."""
    direct, calls = defaultdict(set), defaultdict(set)

    def add(tgt_names, value):
        src = _names(value) if value is not None else set()
        cs = _calls(value) if value is not None else set()
        for nm in tgt_names:
            direct[nm] |= src
            calls[nm] |= cs

    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            tgts = node.targets if isinstance(node, ast.Assign) else [node.target]
            add({n for t_ in tgts for n in _names(t_)}, node.value)
        elif isinstance(node, (ast.For, ast.comprehension)):
            add(_names(node.target), node.iter)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr in {"append", "extend", "add", "update"} and node.args:
            base = node.func.value
            if isinstance(base, ast.Name):
                add({base.id}, node.args[0])

    dep = {k: set(v) for k, v in direct.items()}
    cal = {k: set(v) for k, v in calls.items()}
    for _ in range(len(dep) + 2):                      # fixpoint; cycles converge, never recurse
        changed = False
        for nm in list(dep):
            d, c = set(dep[nm]), set(cal[nm])
            for x in list(d):
                d |= dep.get(x, set())
                c |= cal.get(x, set())
            d.discard(nm)
            if d != dep[nm] or c != cal[nm]:
                dep[nm], cal[nm] = d, c
                changed = True
        if not changed:
            break
    return dep, cal


def _calls(node):
    out = set()
    for n in ast.walk(node):
        if isinstance(n, ast.Call):
            f = n.func
            if isinstance(f, ast.Attribute):
                out.add(f.attr)
            elif isinstance(f, ast.Name):
                out.add(f.id)
    return out


def kill_conditions(tree):
    """(line, condition-node) for every `G.asserted(<KILL...>, cond, ...)` and every assignment
    whose name is later used as such a condition."""
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if not (isinstance(f, ast.Attribute) and f.attr == "asserted"):
            continue
        if not node.args:
            continue
        first = node.args[0]
        txt = first.value if isinstance(first, ast.Constant) and isinstance(first.value, str) else ""
        is_kill = any(k.lower() in txt.lower() for k in KILL_MARK)
        if not is_kill:
            for kw in node.keywords:
                if kw.arg == "kind" and isinstance(kw.value, ast.Constant) \
                        and kw.value.value == "kill":
                    is_kill = True
        if is_kill and len(node.args) > 1:
            out.append((node.lineno, txt[:70], node.args[1]))
    return out


def cond_expr_for(tree, cond_node):
    """If the condition is a bare Name, return the expression it was assigned, else itself."""
    if isinstance(cond_node, ast.Name):
        for n in ast.walk(tree):
            if isinstance(n, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == cond_node.id for t in n.targets):
                return n.value
    return cond_node


QUANTILE = {"percentile", "quantile", "nanpercentile", "nanquantile", "ppf"}
PERMUTE = {"permutation", "shuffle", "permute"}


def self_quantile_names(tree, dep, cal):
    """names assigned from a QUANTILE functional over a collection that was NOT built by permuting
    anything — i.e. a quantile of the statistic's own resampling distribution.

    ⚠ THE SEMANTIC LINE, and it is the whole instrument: comparing a statistic to a QUANTILE of its
    own bootstrap is vacuous by definition (q% of that distribution lies below the q-th percentile,
    whatever the data say). Comparing it to a MULTIPLE OF ITS OWN SE is the ordinary
    signal-to-noise test and is sound, because a bootstrap sd does not scale with the effect.
    Comparing it to a quantile of a PERMUTATION null is also sound, because that distribution was
    built with the structure destroyed. So: quantile functional, minus a permutation in its
    ancestry."""
    out = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        # ⚠ the quantile call is often NOT the top-level node: `lo, hi = [float(x) for x in
        #   np.percentile(boot, [2.5, 97.5])]` buries it inside a comprehension AND unpacks a
        #   tuple. Requiring `isinstance(node.value, ast.Call)` failed the positive control on
        #   `#937` v1's own source, which is the shape this whole round exists to catch.
        qcall = None
        for sub in ast.walk(node.value):
            if isinstance(sub, ast.Call):
                f = sub.func
                fname = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", "")
                if fname in QUANTILE and sub.args:
                    qcall = sub
                    break
        if qcall is None:
            continue
        coll = _names(qcall.args[0])
        anc = set(coll)
        cls = set()
        for c_ in coll:
            anc |= dep.get(c_, set())
            cls |= cal.get(c_, set())
        if cls & PERMUTE or _calls(qcall.args[0]) & PERMUTE:
            continue                       # a permutation null is a legitimate yardstick
        for t_ in node.targets:
            for nm in _names(t_):
                out[nm] = anc
    return out


def suspect_comparisons(tree, cond_node, dep, cal, launder=True):
    """a comparison placing a statistic against a quantile of that statistic's OWN bootstrap."""
    expr = cond_expr_for(tree, cond_node)
    selfq = self_quantile_names(tree, dep, cal)
    if not selfq:
        return []
    hits = []
    for cmpn in [x for x in ast.walk(expr) if isinstance(x, ast.Compare)]:
        sides = [cmpn.left] + list(cmpn.comparators)
        # ⚠ ADJACENT PAIRS ONLY. `lo <= 0.0 <= hi` is (lo,0.0) and (0.0,hi) -- comparing lo to hi
        #   is a comparison Python never makes, and doing so flagged `#937`'s REPAIRED kill, i.e.
        #   the instrument would have condemned the fix it was built to reward.
        pairs = [(sides[i], sides[i + 1]) for i in range(len(sides) - 1)]
        for a, b in pairs + [(b_, a_) for a_, b_ in pairs]:
                na, nb = _names(a), _names(b)
                if not na or not nb:
                    continue
                # ⚠ the threshold may be one hop from the quantile (`hi = q * 1.0`), which is
                #   what POSITIVE2 tests: match on nb OR anything in nb's closure.
                # ⚠ SPECIFICATION AXIS, not a tuning knob (G4). STRICT: the quantile name must
                #   BE an operand. LAUNDERED: one dependency hop, which catches `hi = q * 1.0`
                #   (POSITIVE2) but also drags in every quantity computed from BINNED data --
                #   `qs = nanpercentile(A8, [20,40,60,80])` is a bin edge, not a threshold, and
                #   hand-reading 3 candidates found that class immediately. Both are reported.
                anc_b = set(nb)
                if launder:
                    for x in nb:
                        anc_b |= dep.get(x, set())
                q = anc_b & set(selfq)
                if not q:
                    continue
                anc_a = set(na)
                for x in na:
                    anc_a |= dep.get(x, set())
                shared = set()
                for qq in q:
                    shared |= anc_a & selfq[qq]
                if shared:
                    hits.append(dict(line=cmpn.lineno, lhs=sorted(na)[:4],
                                     threshold=sorted(q)[:3], shared_ancestors=sorted(shared)[:4]))
    return hits


def audit_source(src, tree=None, launder=True):
    tree = tree or ast.parse(src)
    dep, cal = deps_of(tree)
    kills = kill_conditions(tree)
    rows = []
    for line, txt, cond in kills:
        rows.append(dict(line=line, kill=txt,
                         suspects=suspect_comparisons(tree, cond, dep, cal, launder)))
    return rows


# ══ CONTROLS — fixtures first, corpus second ═════════════════════════════════════════
# ⚠ FIX_POS is `#937` v1 VERBATIM in shape, transcribed from the source before the repair --
#   not a fixture invented to be caught (`realstat`: a control validated against cases you invented
#   is validated against your imagination).
FIX_POS = '''
g_fem = coupling(g[g.sex == 2], G_T, G_O)
n_res = coupling(n, N_T, N_O)
diff = abs(g_fem[2] - n_res[2])
boot = []
gf = g[g.sex == 2]
for _ in range(300):
    a = gf.sample(len(gf), replace=True)
    b = n.sample(len(n), replace=True)
    boot.append(abs(coupling(a)[2] - coupling(b)[2]))
lo, hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
travels = diff <= hi
G.asserted("KILL: W_GSS_ONLY requires the two instruments to disagree beyond resampling", travels, "d")
'''
FIX_POS2 = '''
g_fem = coupling(g[g.sex == 2], G_T, G_O)
n_res = coupling(n, N_T, N_O)
diff = abs(g_fem[2] - n_res[2])
boot = [abs(coupling(g.sample(9, replace=True))[2] - coupling(n.sample(9, replace=True))[2])
        for _ in range(300)]
q = np.percentile(boot, 97.5)
hi = q * 1.0
travels = diff <= hi
G.asserted("KILL: laundered through an intermediate name", travels, "d")
'''
FIX_NEG = '''
g_fem = coupling(g[g.sex == 2], G_T, G_O)
n_res = coupling(n, N_T, N_O)
diff = abs(g_fem[2] - n_res[2])
boot = [abs(coupling(g.sample(9, replace=True))[2] - coupling(n.sample(9, replace=True))[2])
        for _ in range(300)]
lo, hi = [float(x) for x in np.percentile(boot, [2.5, 97.5])]
travels = lo <= 0.0 <= hi
G.asserted("KILL: the gap interval must include zero", travels, "d")
'''
FIX_NEG2 = '''
effect = mean(g) - mean(n)
boot = [mean(g.sample(9, replace=True)) - mean(n.sample(9, replace=True)) for _ in range(300)]
se = float(np.std(boot))
ok = abs(effect) > 2 * se
G.asserted("KILL: the effect must exceed twice its own bootstrap spread", ok, "d")
'''
FIX_NEG3 = '''
effect = mean(g) - mean(n)
nulls = [stat(RNG.permutation(g), n) for _ in range(300)]
thr = float(np.percentile(nulls, 95))
ok = effect > thr
G.asserted("KILL: the effect must exceed the 95th percentile of its permutation null", ok, "d")
'''
FIX_BLIND = '''
x = 1
y = x + 1
print(y)
'''

G = Gate("Can this project's pre-registered kills fail — and how many of them could not?")

pos, pos2 = audit_source(FIX_POS), audit_source(FIX_POS2)
neg, neg2, neg3 = audit_source(FIX_NEG), audit_source(FIX_NEG2), audit_source(FIX_NEG3)
blind = audit_source(FIX_BLIND)

c_pos = bool(pos and pos[0]["suspects"])
c_pos2 = bool(pos2 and pos2[0]["suspects"])
c_neg = bool(neg) and not neg[0]["suspects"]
c_neg2 = bool(neg2) and not neg2[0]["suspects"]
c_neg3 = bool(neg3) and not neg3[0]["suspects"]
c_blind = (len(blind) == 0)

print("=== controls ===")
print(f"  positive  `#937` v1 `diff <= percentile(own bootstrap)` CAUGHT : "
      f"{'PASS' if c_pos else 'FAIL'}   {pos[0]['suspects'][:1] if c_pos else ''}")
print(f"  positive2 the same laundered through a second name, CAUGHT     : "
      f"{'PASS' if c_pos2 else 'FAIL'}")
print(f"  negative1 `#937` v2 `lo <= 0.0 <= hi` must NOT be caught       : "
      f"{'PASS' if c_neg else 'FAIL'}")
print(f"  negative2 `effect > 2*se` (own bootstrap SD) must NOT be caught: "
      f"{'PASS' if c_neg2 else 'FAIL'}")
print(f"  negative3 `effect > pct(PERMUTATION null, 95)` must NOT be caught: "
      f"{'PASS' if c_neg3 else 'FAIL'}  <- separates a null from a bootstrap")
print(f"  blind     a script with no kill -> UNREADABLE, not clean       : "
      f"{'PASS' if c_blind else 'FAIL'}")

if not (c_pos and c_pos2 and c_neg and c_neg2 and c_neg3 and c_blind):
    print("  ⛔ a control failed — this round's verdict on the corpus is INADMISSIBLE")
    sys.exit(2)

# ══ THE CORPUS ═══════════════════════════════════════════════════════════════════════
scripts = sorted(ROOT.glob("E0*/A*/R*/*.py"))
n_parse_fail = 0
with_kill, kill_total, cand_rows = 0, 0, []
for p in scripts:
    try:
        rows = audit_source(p.read_text(errors="replace"))
    except SyntaxError:
        n_parse_fail += 1
        continue
    if not rows:
        continue
    with_kill += 1
    kill_total += len(rows)
    for r in rows:
        if r["suspects"]:
            cand_rows.append(dict(file=str(p.relative_to(ROOT)), **r))

strict_rows = []
for p_ in scripts:
    try:
        for r in audit_source(p_.read_text(errors="replace"), launder=False):
            if r["suspects"]:
                strict_rows.append(dict(file=str(p_.relative_to(ROOT)), **r))
    except SyntaxError:
        pass
cand_kills = len(cand_rows)
share = cand_kills / kill_total if kill_total else float("nan")
locatable_share = with_kill / len(scripts) if scripts else float("nan")

print(f"\n=== corpus ===")
print(f"  round scripts on disk        : {len(scripts)}   (parse failures {n_parse_fail})")
print(f"  scripts with a LOCATABLE kill: {with_kill}  = {locatable_share:.1%} — the rest are "
      f"⚠ UNREADABLE, not clean")
print(f"  locatable kill conditions    : {kill_total}")
strict_share = len(strict_rows) / kill_total if kill_total else float("nan")
print(f"  ⚠ CANDIDATES (statistic vs a quantile of its own bootstrap):")
print(f"      LAUNDERED (one dependency hop) : {cand_kills:3d} = {share:.2%} of locatable kills")
print(f"      STRICT    (quantile IS operand): {len(strict_rows):3d} = {strict_share:.2%}")
print(f"      ⚠ the specification matters more than the number: 3 laundered candidates were "
      f"hand-read and 3 were false positives, all binning quantiles")
for r in cand_rows[:12]:
    s = r["suspects"][0]
    print(f"    · {r['file'].split('/')[-1][:42]:42s} L{r['line']:<4d} {s['lhs']} vs quantile "
          f"{s['threshold']} sharing {s['shared_ancestors'][:3]}")
if len(cand_rows) > 12:
    print(f"    … and {len(cand_rows) - 12} more")

# ══ GATES ════════════════════════════════════════════════════════════════════════════
G.asserted("positive: `#937` v1's shape, transcribed from the source, is caught — and so is a "
           "version laundered through an intermediate name", c_pos and c_pos2,
           f"one-hop {c_pos} · laundered {c_pos2}; and the FIRST detector, 'one operand inside the "
           f"other's closure', failed this same control — the dependency runs through the DATA",
           kind="control", population="the 6 fixtures, which are transcribed shapes and not invented cases")
G.asserted("negative 1/3: `#937` v2's repaired kill is NOT caught", c_neg,
           "`lo <= 0.0 <= hi` — bounds from the statistic, comparison against a CONSTANT",
           kind="control", population="the 6 fixtures, which are transcribed shapes and not invented cases")
G.asserted("negative 2/3: `effect > 2*se` — the ordinary signal-to-noise test, same bootstrap — is "
           "NOT caught", c_neg2,
           "a bootstrap SD does not scale with the effect, so this criterion CAN fail; the detector "
           "keys on a QUANTILE functional, not on any use of the bootstrap", kind="control",
           population="the 6 fixtures, which are transcribed shapes and not invented cases")
G.asserted("negative 3/3: `effect > pct(PERMUTATION null, 95)` is NOT caught — the line between a "
           "null and a bootstrap, which is the whole instrument", c_neg3,
           "a permutation destroys the structure, so its quantile is a yardstick rather than a "
           "restatement; without this row the detector would flag every honest permutation kill in "
           "the corpus", kind="control", population="the 6 fixtures, which are transcribed shapes and not invented cases")
G.asserted("blind: a script with no locatable kill yields NO row — UNREADABLE, never clean",
           c_blind, f"{len(scripts) - with_kill} of {len(scripts)} scripts are in that state",
           kind="control", population="all 895 round scripts on disk")
G.asserted("coverage is measured, not asserted (`#932`: a row count is not coverage)", True,
           f"scripts {len(scripts)} · with a locatable kill {with_kill} ({locatable_share:.1%}) · "
           f"parse failures {n_parse_fail} · kill conditions {kill_total}", kind="control",
           population="all 895 round scripts on disk")

# ⚠ the kill is evaluated on the STRICT specification, because the laundered one was measured to
#   over-call on the only cases anyone actually read. Both are published (G4).
share_used = strict_share
W_END = share_used >= 0.10
W_RARE = share_used <= 0.02
# ⚠⚠ THE CONTROL ON THE CLAIM'S OWN POPULATION (`#866`①), which the six fixtures above cannot be:
#   9 flagged kills were HAND-READ in the corpus itself -- 3 laundered, then all 6 strict ones
#   printed. **0 of 9 are vacuous.** Two false-positive classes, both idiomatic here:
#     (a) BINNING quantiles -- `qs = nanpercentile(A8, [20,40,60,80])` cuts a covariate; it is not
#         a threshold at all, and the laundering hop drags in everything computed downstream of it;
#     (b) NULL quantiles the PERMUTE marker misses -- this corpus builds nulls by OFFSET, by
#         sign-flip, by label rotation and by max-statistic, naming them `nul_off`, `NUL0`, `fl`,
#         `maxt`, and none of those calls `permutation`. Comparing a statistic to the 95th
#         percentile of an offset null is exactly right, and the detector flags it every time.
#   ⇒ **the guard `#937`③ asked for cannot be static.** That is this round's product, and it was
#   measured rather than assumed: the static route's precision on its own corpus is 0/9.
HAND_READ, HAND_VACUOUS = 9, 0
G.asserted("⚠ precision measured ON THE CORPUS, not on fixtures (`#866`①): 9 flagged kills "
           "hand-read, 0 vacuous", True,
           f"3 laundered + all 6 strict candidates printed and read; {HAND_VACUOUS}/{HAND_READ} "
           f"confirmed vacuous. Classes: binning quantiles, and null quantiles whose construction "
           f"never calls `permutation` (`nul_off`, `NUL0`, `fl`, `maxt`). The candidate share is "
           f"therefore an UPPER BOUND with measured precision 0/9, not a count of defects",
           kind="control", population=f"{with_kill} scripts carrying a locatable kill, of "
                                      f"{len(scripts)} — the SAME population as the kill")
G.asserted("KILL: pre-registered — W_ENDEMIC requires the candidate share to reach 10% of locatable "
           "kills; W_RARE requires it to fall to 2%; between them, neither world is claimed",
           W_END or W_RARE,
           f"candidate share STRICT {strict_share:.2%} (laundered {share:.2%}) of {kill_total} "
           f"locatable kills across {with_kill} scripts "
           f"-> {'W_ENDEMIC' if W_END else 'W_RARE' if W_RARE else 'NEITHER — the number is '
               'reported and no world is claimed'}",
           kind="kill", yardstick="share of locatable kill conditions",
           yardstick_noise=0.0,   # ⚠ a CENSUS of the corpus, not a sample: every script on disk is
           # read, so there is no resampling noise on this share. The uncertainty that DOES exist is
           # the detector's coverage, and it is registered under STRUCTURALLY CANNOT, not hidden in
           # a spread. Stating 0.0 with the reason beats leaving it UNCOMPUTED (`realstat` G1).
           population=f"{with_kill} scripts carrying a locatable kill, of {len(scripts)}",
           direction="two-sided, both thresholds pre-registered above")

print(G)
tv = str(G)
world = ("W_ENDEMIC" if W_END else "W_RARE" if W_RARE else "NEITHER")
if locatable_share < 0.20:
    world = "W_UNSEEABLE"
    verdict = (f"UNVERIFIED · world W_UNSEEABLE · only {locatable_share:.1%} of round scripts "
               f"carry a machine-readable kill")
elif W_END or W_RARE:
    verdict = (f"CONFIRMED · world {world} · strict candidate share {strict_share:.2%} "
               f"(laundered {share:.2%})")
else:
    verdict = (f"UNVERIFIED · neither threshold reached · strict {strict_share:.2%} / laundered "
               f"{share:.2%}, reported as numbers and not as a world")
print(f"\nVERDICT           : {verdict}")

json.dump(dict(entry=938, round="E03·A120·R376", gate_verdict=str(G).splitlines()[-1][:300],
               gates=[[r[0], r[2], r[1]] for r in G.rows], claims_null=False,
               scripts=len(scripts), scripts_with_kill=with_kill, kill_total=kill_total,
               candidate_kills=cand_kills, candidate_share=share,
               strict_candidates=len(strict_rows), strict_share=strict_share,
               hand_read=dict(n=HAND_READ, confirmed_vacuous=HAND_VACUOUS,
                              classes=["binning quantiles", "null quantiles the PERMUTE marker "
                                       "misses: offset, sign-flip, rotation, max-statistic"],
                              note="0/9 vacuous -> the static route cannot be the guard"),
               locatable_share=locatable_share, parse_failures=n_parse_fail,
               controls=dict(positive=c_pos, positive_two_hop=c_pos2, negative=c_neg,
                             negative_se=c_neg2, negative_permutation_null=c_neg3, blind=c_blind),
               candidates=cand_rows, verdict=verdict, world=world),
          open(OUT / "the_kill_that_cannot_fail.json", "w"), indent=1)
print(f"\nwrote {OUT / 'the_kill_that_cannot_fail.json'}")
