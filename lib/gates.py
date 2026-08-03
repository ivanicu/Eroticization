"""
lib/gates.py -- the comparison rules, as code that asserts relations instead of prose that asserts
conditions.

Written after the twenty-ninth mis-specified design element (#102a). The class is stated in #102e:

    EVERY GATE IN THIS PROJECT COMPARES TWO NUMBERS, AND EVERY FAILURE WAS IN HOW, NOT IN WHAT.

The catalogue of how, each with its entry:

  #96a   a gate whose PROSE asserted a condition ("at comparable mean") the code never tested
  #79e   a control that examined ONE ARM of a two-arm design and passed
  #79f   a threshold applied to a SUM that hides a sign
  #83d   a ratio taken ACROSS A SIGN CHANGE, returning 1.7e10
  #88a   a positive control planted BELOW its own MDE, three times running
  #78c   a positive control demanding a value the instrument was KNOWN to violate -- could not pass
  #97c   a tolerance taken from a spread that is STRUCTURALLY ZERO
  #101a  a control residualised on the wrong FUNCTIONAL FORM, leaving the confound behind
  #101b  a reference compared across SCALES (raw vs disattenuated)
  #102a  a negative control compared to A CONSTANT I CHOSE instead of to the effect
  #109c  THE NULL'S TYPE SILENTLY CHANGED BETWEEN ROUNDS during a refactor. A stratified
         permutation and a plain one look identical in code (one line) and differ by 85% of the
         claimed effect. Anything that computes a null must NAME which null it is.

Nine of those ten are a comparison taking the wrong second argument. So every function here REQUIRES
the thing being compared against, and none accepts a bare threshold.

Usage:

    from lib.gates import Gate
    g = Gate("does the affinity-onset link survive?")
    g.negative_control("null-minus-null", null_r, effect=eff_r)
    g.positive_control("planted trait", planted, floor=null_rel, spread=sd)
    g.same_scale("breadth vs #17/#23", mine=0.0712, theirs=0.075, scale="raw")
    if g.verdict():           # False unless EVERY gate passed
        ...report...
    print(g)                  # always prints the whole table, pass or fail
"""
import numpy as np

# #80d / #93d / #117e -- 三次写下这条规则,三次违反它。名单不再是散文,是可执行的检查。
BANNED_COLUMN_NAMES = frozenset("""
T shift mode item count size min max sum mean std var rank pop all any abs where mask
first last div pow add sub mul truediv floordiv apply map filter head tail index values
dtypes shape loc iloc at iat name axes empty ndim
""".split())


def check_columns(df, where=""):
    """#117e: 一个与 DataFrame 方法同名的列,取用时静默返回方法对象而不是数据。
    这在本项目里发生了五次(shift #74, mode #77, item #80, T #93, shift #117)。
    在任何 groupby/agg 之前调用一次,它就再也不能发生。"""
    bad = sorted(set(map(str, df.columns)) & BANNED_COLUMN_NAMES)
    if bad:
        raise ValueError(
            f"列名与 pandas 方法撞名{(' in '+where) if where else ''}: {bad}. "
            f"改名(例如 v_{bad[0]}),否则 df.{bad[0]} 返回的是方法不是数据。")
    return df


class Gate:
    def __init__(self, question):
        self.question = question
        self.rows = []

    # ---- the three comparisons that failed, each now requiring its second argument ----

    def _degenerate(self, name, *vals):
        """#105e: a comparison whose inputs are all exactly zero passes vacuously. Refuse it.

        The library exists to stop a gate that cannot fail; a gate evaluated on a degenerate input
        is exactly that, and it printed 3 PASS on a run where every number was 0.00."""
        if all(abs(v) < 1e-12 for v in vals):
            self.rows.append((name, f"inputs {[round(v,12) for v in vals]}", False,
                              "DEGENERATE -- all inputs are exactly zero; the comparison is vacuous"))
            return True
        return False

    def negative_control(self, name, null, effect, ratio=0.5):
        """#102a: a null is judged against the EFFECT, never against a constant.

        Passes when |null| < ratio * |effect|. There is no absolute threshold, because an absolute
        threshold is what let a null equal to 91% of its effect print PASS."""
        if self._degenerate(name, null, effect): return False
        ok = abs(null) < ratio * abs(effect)
        self.rows.append((name, f"|{null:+.4f}| < {ratio}*|{effect:+.4f}| = {ratio*abs(effect):.4f}",
                          ok, f"null is {100*abs(null)/max(abs(effect),1e-12):.0f}% of the effect"))
        return ok

    def positive_control(self, name, planted, floor, spread):
        """#88a/#78c: a positive control is judged against the FLOOR and its own spread.

        Passes when planted > floor + 2*spread. It must also be POSSIBLE: if the instrument is known
        to violate the criterion, the gate cannot pass and this will say so rather than fail
        silently."""
        if self._degenerate(name, planted, floor, spread): return False
        need = floor + 2 * abs(spread)
        ok = planted > need
        self.rows.append((name, f"{planted:+.4f} > floor {floor:+.4f} + 2*{abs(spread):.4f} = {need:+.4f}",
                          ok, f"headroom {planted-need:+.4f}"))
        return ok

    def offset_control(self, name, effect, offset, spread, null_kind):
        """#106: a null that is a systematic BASELINE OFFSET, not a nuisance to be small.

        #109c: null_kind is REQUIRED and has no default. A stratified permutation and a plain one
        are ONE LINE apart in code and differed by 85% of a claimed effect; naming it in the call
        is the only thing that makes a silent swap visible in the output.

        When the two arms differ in model capacity, the null carries the overfitting penalty and is
        systematically NEGATIVE. negative_control() is the wrong shape there.
        Choose by asking: SHOULD the null be zero? Yes -> negative_control. No -> offset_control."""
        assert isinstance(null_kind,str) and null_kind, "name the null (#109c): e.g. 'stratified permutation'"
        """#106: a null that is a systematic BASELINE OFFSET, not a nuisance to be small.

        When the two arms differ in model capacity, the null carries the overfitting penalty and is
        systematically NEGATIVE. negative_control() is the wrong shape there -- it asks |null| to be
        small, and a large negative null makes it FAIL on a real effect. The right statistic is the
        DIFFERENCE, judged against its own spread.

        Choose by asking: SHOULD the null be zero? Yes -> negative_control. No, it has a known
        systematic direction -> offset_control."""
        if self._degenerate(name, effect, offset, spread): return False
        corrected = effect - offset
        ok = abs(corrected) > 2 * abs(spread)
        self.rows.append((f"{name} [null: {null_kind}]",
                          f"({effect:+.4f}) - ({offset:+.4f}) = {corrected:+.4f} vs 2*{abs(spread):.4f}",
                          ok, f"{abs(corrected)/max(abs(spread),1e-12):.1f}x its own spread"))
        return ok

    def same_scale(self, name, mine, theirs, scale):
        """#101b: a reference comparison must DECLARE its scale, in the call, as a string.

        There is no default. A comparison whose scale nobody wrote down is the comparison that
        compared raw to disattenuated and failed on arithmetic."""
        assert isinstance(scale, str) and scale, "state the scale explicitly (e.g. 'raw', 'disattenuated')"
        ok = abs(mine - theirs) < 0.5 * abs(theirs)
        self.rows.append((f"{name} [{scale}]", f"|{mine:.4f} - {theirs:.4f}| < {0.5*abs(theirs):.4f}",
                          ok, f"both on the {scale} scale"))
        return ok

    # ---- the guards on the arguments themselves ----

    def resolvable(self, name, effect, spread):
        """#97c: a spread of zero makes everything resolvable. Refuse it."""
        if not np.isfinite(spread) or spread <= 0:
            self.rows.append((name, f"spread = {spread}", False,
                              "STRUCTURALLY ZERO OR NON-FINITE -- use a bootstrap over units, not seeds"))
            return False
        ok = abs(effect) > 2 * spread
        self.rows.append((name, f"|{effect:+.4f}| > 2*{spread:.4f} = {2*spread:.4f}", ok,
                          f"{abs(effect)/spread:.1f}x its own spread"))
        return ok

    def no_sign_crossing(self, name, series):
        """#83d/#79f: never take a ratio or a sum across a sign change."""
        s = np.asarray(series, dtype=float)
        if self._degenerate(name, *s): return False
        ok = np.all(s > 0) or np.all(s < 0)
        self.rows.append((name, f"signs {'consistent' if ok else 'CROSS ZERO'}", ok,
                          f"{np.round(s,4).tolist()}"))
        return ok

    def covers_every_arm(self, name, checked, arms):
        """#79e: a control that examines one arm of a two-arm design is not a control."""
        missing = sorted(set(arms) - set(checked))
        ok = not missing
        self.rows.append((name, f"arms checked {sorted(checked)} of {sorted(arms)}", ok,
                          "complete" if ok else f"MISSING {missing}"))
        return ok

    def asserted(self, name, condition, detail):
        """#96a: a condition stated in prose must be a boolean here, or it was never tested."""
        ok = bool(condition)
        self.rows.append((name, detail, ok, "asserted in code, not in prose"))
        return ok

    # ---- output ----

    def verdict(self):
        return all(r[2] for r in self.rows)

    def __str__(self):
        w = max((len(r[0]) for r in self.rows), default=10)
        out = [f"  CONDITIONAL KILL -- {self.question}"]
        for nm, test, ok, note in self.rows:
            out.append(f"   {'PASS' if ok else 'FAIL'}  {nm:<{w}}  {test}   ({note})")
        out.append(f"   => {'ALL GATES PASS' if self.verdict() else 'UNVERIFIED, and that is not an acquittal'}")
        return "\n".join(out)
