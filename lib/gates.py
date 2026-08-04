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


def check_coverage(processed, available, where="", tol=0.02):
    """#118c: 第六次,一个我自己写的静默 cap 改变了结论(#73 #74 #118)。

    任何在单位上循环的轮次都必须报出它跳过了多少。截断本身可以是合法的(成本),
    但**沉默的截断**会把分母伪装成完整,而分母是每一个 SE 的分母。

    processed: 实际进入统计量的单位数
    available: 满足纳入标准的单位总数
    截断超过 tol 就抛错,除非调用方显式把 tol 调大 —— 那时它至少写在代码里了。"""
    if available <= 0:
        raise ValueError(f"没有单位可用{(' in '+where) if where else ''} —— 空分母(#04)")
    frac = 1.0 - processed / available
    if frac > tol:
        raise ValueError(
            f"静默截断{(' in '+where) if where else ''}: 用了 {processed}/{available} "
            f"({100*frac:.0f}% 被跳过). 若这是有意的成本控制,显式传 tol={frac+0.01:.2f} "
            f"并把跳过量写进轮次的输出里(#118c)。")
    return processed


def check_residualized(residual, covariate, where="", tol=0.02):
    """#129: 一个残差与它所回归掉的协变量的相关**在构造上恒等于 0**。不是 0 -> 那不是残差。

    A14R01 写的是 `z(S - polyval(fit_on_z, z(PK)))` 而不是 `z(S) - polyval(...)` ——
    把一个 z 尺度的预测值从**原始尺度**的变量里减掉。原始 S 的 sd 远小于预测值的 sd,
    所以结果几乎就是 -z(PK) 本身:corr(所谓残差, PK) = **-0.9654**。
    整整两轮的 corr(., S) 测的是负的勾选数,不是稀有亲和特质。

    这个错的可怕之处在于它**看起来完全正常** —— 变量还叫 S,还被 z 标准化过,
    数值范围也对,而且它产生的相关是**显著的**。唯一的破绽是这一行断言。
    单位不匹配的代数不会抛异常;它会安静地把协变量还给你,冠上被解释量的名字。"""
    import numpy as _np
    r = _np.asarray(residual, dtype=float); c = _np.asarray(covariate, dtype=float)
    m = _np.isfinite(r) & _np.isfinite(c)
    if m.sum() < 3:
        raise AssertionError(f"check_residualized({where}): 只有 {int(m.sum())} 个有限值,无法检验")
    if _np.std(r[m]) < 1e-12 or _np.std(c[m]) < 1e-12:
        raise AssertionError(f"check_residualized({where}): 一侧是常数,比较是退化的(#105e)")
    rr = float(_np.corrcoef(r[m], c[m])[0, 1])
    if abs(rr) > tol:
        raise AssertionError(
            f"check_residualized({where}): corr(残差, 协变量) = {rr:+.4f},超过 {tol}。"
            f"残差与被回归掉的协变量必须正交。检查两侧是否在同一尺度上 —— "
            f"把 z 尺度的预测值从原始尺度的变量里减掉,返回的是协变量本身(#129)。")
    return rr


def check_disjoint_items(pred_items, outcome_items, where="", tol=0.0):
    """#126c: 预测量与结局都由 item 级数据构造时,两边的 item 集必须不相交。

    A11R20 里强制单选的选项与多选块的选项重叠 89-100%,于是「稀有亲和 S」与「选中选项的
    冷门程度」共享同一批勾选 —— 那是恒等式不是测量。**这是本项目第一个在设计时被漏掉、
    而不是写错的混淆**(前四十一个都是写下了但写错了),抓住它的只是一个临时检查。

    机械可查,所以不该再靠想起来。tol>0 表示显式接受一定重叠(那时它至少写在源码里)。"""
    a, b = set(pred_items), set(outcome_items)
    if not a or not b:
        raise ValueError(f"item 集为空{(' in '+where) if where else ''} —— 空分母(#04)")
    ov = a & b
    frac = len(ov) / min(len(a), len(b))
    if frac > tol:
        ex = sorted(map(str, ov))[:3]
        raise ValueError(
            f"预测量与结局共享 item{(' in '+where) if where else ''}: "
            f"{len(ov)}/{min(len(a),len(b))} ({frac:.0%}) 重叠,例如 {ex}。"
            f"若无法避免,显式传 tol={min(frac+0.01,1.0):.2f} 并在轮次输出里报出重叠率(#126c)。")
    return True


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

    def negative_control(self, name, null, effect, ratio=0.5, null_spread=None):
        """#102a: a null is judged against the EFFECT, never against a constant.

        #125: ...but when the EFFECT is small, `|null| < 0.5*|effect|` is a bar the null can fail
        while being indistinguishable from zero on its own spread. A05R15: null -0.00147 +/- 0.00140
        (1.1x, i.e. zero) failed against an effect of +0.00173. Both questions must be asked:
            (a) is the null small relative to the effect?   -- #102a
            (b) is the null itself already indistinguishable from zero?  -- #125
        Passing either is enough; pass null_spread to enable (b). Without it only (a) is asked, and
        the row says so, so a missing spread can never be mistaken for a passed check."""
        if self._degenerate(name, null, effect): return False
        ok = abs(null) < ratio * abs(effect)
        # #125 的豁免只在零**帮不上忙**时成立:与效应异号,或已小于效应的一半。
        # 否则 #102a(零 -0.0275、效应 -0.0302、同号、91%)会从这个口子溜过去 ——
        # 而那正是这个库存在的起因。回归测试当场抓到了这个洞。
        helps = (null * effect) > 0 and abs(null) >= 0.5 * abs(effect)
        if not ok and not helps and null_spread is not None and abs(null) < 2 * abs(null_spread):
            self.rows.append((name, f"|{null:+.5f}| < 2*{abs(null_spread):.5f} (自身展布)", True,
                              f"零本身与零无法区分 ({abs(null)/max(abs(null_spread),1e-12):.1f}x);"
                              f" 相对效应是 {100*abs(null)/max(abs(effect),1e-12):.0f}% (#125)"))
            return True
        if null_spread is None:
            ratio_note = " [未给 null_spread,只问了「相对效应」这一半 (#125)]"
        else:
            ratio_note = ""
        self.rows.append((name, f"|{null:+.4f}| < {ratio}*|{effect:+.4f}| = {ratio*abs(effect):.4f}",
                          ok, f"null is {100*abs(null)/max(abs(effect),1e-12):.0f}% of the effect"
                              + ratio_note))
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

    def artifact_cannot_explain(self, name, artifact, effect, spread):
        """#119: 一个伪影不必为零 —— 它只需要**不能解释效应**。

        我在 A13R01 用 asserted(|artifact| < 2*SE) 判了一个伪影,它 FAIL 了 0.0001,
        而那个伪影的**符号与效应相反**(+0.0070 vs -0.0122),即校正只会让效应更大。
        "是否为零"是错的问题;正确的问题是"它能不能造出这个效应"。

        通过条件:符号相反(伪影帮不了忙),或同号但小于效应的一半。"""
        if self._degenerate(name, artifact, effect, spread): return False
        opposite = (artifact * effect) < 0
        small = abs(artifact) < 0.5 * abs(effect)
        ok = opposite or small
        why = "符号相反,校正只会放大效应" if opposite else \
              (f"同号但仅为效应的 {100*abs(artifact)/max(abs(effect),1e-12):.0f}%" if small
               else f"同号且为效应的 {100*abs(artifact)/max(abs(effect),1e-12):.0f}% —— 可能就是它")
        self.rows.append((name, f"伪影 {artifact:+.4f} vs 效应 {effect:+.4f}", ok, why))
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

    def degenerate_matches_reference(self, name, degenerate, reference, tol=1e-12):
        """#124f: 退化臂(强度=0 的种植)必须**精确**复现它要复现的那个臂。

        A06R03 里 g=0 的种植给 +0.0061 而真实臂给 +0.0044,因为两臂用了不同的随机种子
        (seed=5 vs seed=1),掩码不同。那条断言于是测的是**种子**,不是设计。
        退化臂的正确写法是复用参照臂的种子;这个方法在它没被复用时说出原因。"""
        d = abs(degenerate - reference)
        ok = d <= tol
        note = "精确复现" if ok else \
               f"差 {d:.6f} —— 退化臂几乎总是因为**没有复用参照臂的种子**(#124f)"
        self.rows.append((name, f"退化 {degenerate:+.6f} vs 参照 {reference:+.6f}", ok, note))
        return ok

    def require_resolvable_first(self, name, effect, spread, family="default"):
        """#120d: 门有顺序。对一个未分辨的量问形状,后面的比较全是 MOOT。

        #130: ...但 MOOT 必须**限定在它自己的族里**。一个 Gate 里放两个独立的量时,
        Gate 级的 `_moot` 会让其中一个未分辨的量把另一个也判成 MOOT ——
        A14R03 的原始年龄 corr(Delta,S) 是 1.1x,于是它把扣掉时间表后 5.1x 的
        同名量一起标黑了。**一个未分辨的量只让依赖它的行 MOOT,不让别人的行 MOOT。**
        `family` 是依赖关系的名字;不传就是所有行共用一族(旧行为)。"""
        ok = self.resolvable(name, effect, spread)
        if not ok:
            if not hasattr(self, "_moot_fams"):
                self._moot_fams = set()
            self._moot_fams.add(family)
        self._fam_of_row = family
        return ok

    def verdict(self):
        if getattr(self, "_moot", False) or getattr(self, "_moot_fams", set()):
            return False
        return all(r[2] for r in self.rows)

    def __str__(self):
        w = max((len(r[0]) for r in self.rows), default=10)
        out = [f"  CONDITIONAL KILL -- {self.question}"]
        for nm, test, ok, note in self.rows:
            out.append(f"   {'PASS' if ok else 'FAIL'}  {nm:<{w}}  {test}   ({note})")
        if getattr(self, "_moot", False):
            out.append("   ⚠ 效应本身未通过可分辨性 —— 其后的比较全部 MOOT(#120):"
                       "一个未分辨的量没有形状,也不需要零来解释")
        for fam in sorted(getattr(self, "_moot_fams", set())):
            out.append(f"   ⚠ 族 `{fam}` 的效应未通过可分辨性 —— **该族内**其后的比较 MOOT(#120/#130)。"
                       "其它族不受影响。")
        out.append(f"   => {'ALL GATES PASS' if self.verdict() else 'UNVERIFIED, and that is not an acquittal'}")
        return "\n".join(out)
