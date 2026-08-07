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
# ⚠ #197a:这份名单原本是**手写的 45 个**,而 pandas 实际公开 232 个属性 —— **漏掉 189 个**,
#   其中包括 `cov`(`#166c` 咬过我,我改了列名却从没把它加进名单)、`diff`(`#197` 又咬一次)、
#   `corr` · `unique` · `round` · `sample` · `describe` …
#   **一个手写的、关于第三方库 API 的名单,保证会漂。** 改成从 pandas 自己取。
#   保留手写集作为**下界**(防止 pandas 某天改名后名单缩水)。
_HANDWRITTEN = frozenset("""
T shift mode item count size min max sum mean std var rank pop all any abs where mask
first last div pow add sub mul truediv floordiv apply map filter head tail index values
dtypes shape loc iloc at iat name axes empty ndim
""".split())

def _pandas_attrs():
    try:
        import pandas as _pd
        return {a for a in (set(dir(_pd.DataFrame)) | set(dir(_pd.Series)))
                if not a.startswith('_')}
    except Exception:
        return set()

BANNED_COLUMN_NAMES = frozenset(_HANDWRITTEN | _pandas_attrs())


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
    这在本项目里发生了**八次**:shift #74 · mode #77 · item #80 · T #93 · shift #117 ·
    **cov #166c** · **item #184a** · **diff #197a** —— 后三次都在守卫存在之后,
    因为名单是**手写的**,漏掉了它们(见 #197a)。名单现在从 pandas 自己取。
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

    def negative_control(self, name, null, effect, ratio=0.5, null_spread=None, null_kind=None):
        """#102a: a null is judged against the EFFECT, never against a constant.

        ⚠ #713:**`null` 位收到的,必须与 `null_kind` 描述的是同一个量。**
        实测 33 处调用里,这个位置收到过**四种不同的东西**:安慰剂的观测值 11 处、
        零分布的中位 8 处、零分布的 95% 分位 6 处、其他常数 8 处 ——
        **而本函数对它们一视同仁,于是同一句 `null < 0.5*effect` 意味着三种不同的强度。**
        `null_kind` 挡住了「不说零是什么」,**挡不住「说了做法却塞进另一种量」。**
        **调用者自己负责让这两者对齐;本函数不能替你检查。**

        #125: ...but when the EFFECT is small, `|null| < 0.5*|effect|` is a bar the null can fail
        while being indistinguishable from zero on its own spread. A05R15: null -0.00147 +/- 0.00140
        (1.1x, i.e. zero) failed against an effect of +0.00173. Both questions must be asked:
            (a) is the null small relative to the effect?   -- #102a
            (b) is the null itself already indistinguishable from zero?  -- #125
        Passing either is enough; pass null_spread to enable (b). Without it only (a) is asked, and
        the row says so, so a missing spread can never be mistaken for a passed check."""
        if self._degenerate(name, null, effect): return False
        ok = abs(null) < ratio * abs(effect)
        # ⚠ #758(回测于 2026-08-06,`#748`③ 欠了四轮):这一行**对符号是瞎的**,而回测发现
        #   问题的位置和 `#748` 猜的不一样 —— 282 个调用点里 **70 处(24.8%)在传入前就
        #   `abs()` 了**,符号在进函数之前就没了;另有 **35 处(12.4%)没给 `null_spread`**,
        #   于是 `#125` 那条豁免根本不会触发。**⇒ 修的地方在调用点,不在这一行。**
        #   这里只加**披露**,不改判:反号的零其实是**支持**效应的证据,而 `ok` 会因幅度大判它失败。
        #   ⚠ 按 `#623`,任何会改判的改动必须先量出翻转数;这一条**构造上不可能翻转**(只加注记)。
        _opp = (null * effect) < 0
        _sign_note = ("  ⚠ 零与效应**反号** —— 反号的零是支持效应的证据,"
                      "而上面这条比较对符号是瞎的,请自行判读 (#758)" if _opp else "")
        # #125 的豁免只在零**帮不上忙**时成立:与效应异号,或已小于效应的一半。
        # 否则 #102a(零 -0.0275、效应 -0.0302、同号、91%)会从这个口子溜过去 ——
        # 而那正是这个库存在的起因。回归测试当场抓到了这个洞。
        helps = (null * effect) > 0 and abs(null) >= 0.5 * abs(effect)
        if not ok and not helps and null_spread is not None and abs(null) < 2 * abs(null_spread):
            self.rows.append((name, f"|{null:+.5f}| < 2*{abs(null_spread):.5f} (自身展布)", True,
                              f"零本身与零无法区分 ({abs(null)/max(abs(null_spread),1e-12):.1f}x);"
                              f" 相对效应是 {100*abs(null)/max(abs(effect),1e-12):.0f}% (#125)" + _sign_note))
            return True
        if null_spread is None:
            ratio_note = " [未给 null_spread,只问了「相对效应」这一半 (#125)]"
        else:
            ratio_note = ""
        # #217a:这个项目至少有两种置换零(题内跨人 · 人内),而它们在同一个量上
        #        给出**符号相反**的值(`#216c`:−0.0792 vs +0.1309)。
        #        所以"净值"离开了它的零就没有意义。`offset_control` 一直强制命名零的种类,
        #        而 `negative_control` 没有 —— 现在补上,但用**标注**而不是报错,
        #        以免让二百多个既有轮次一次全红(#197a 收紧名单时的同一个权衡)。
        if null_kind is None:
            ratio_note += " ⚠[未命名零的方案 —— 题内跨人?人内?两者在同一量上可反号 (#217a)]"
        else:
            ratio_note += f" [零的方案: {null_kind}]"
        # #758:反号的注记必须挂在**主判行**上 —— 第一版只挂到 `#125` 分支,
        # 而反号且幅度大的那一格根本走不到那条分支,于是注记一次也没出现过。
        # ⚠ 这正是这一族的又一次:**注记挂错了行,和判词比错了对象是同一个错。**
        ratio_note += _sign_note
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

    def identity_control(self, name, observed, expected, tol, what):
        """#761: 两个量**必须是同一个** —— 不同就说明仪器坏了。

        ⚠ 为什么补这一条：`#760` 里我要断言「把控制量换成常数/打乱后的量，偏相关必须**回到**偏前」，
        库里没有这个形状，于是我用了 `offset_control` —— 而它测的是「差要够大」，
        **用差值检测器去断言等式**，判词当场 FAIL 在一个完全正常的仪器上。
        这是 `#728`·`#748`·`#750`·`#758` 那一族的第七次，这次犯在**闸的选型**上。

        ⚠ P6 代理账：
          PROPERTY   这两个数指的是同一个量
          PROXY      |observed - expected| <= tol
          IMPLICATION 只有一个方向可靠：**超出容差 -> 它们确实不是同一个量**（可靠）。
                     反过来不成立：**相等不证明我算的是我以为的那个量**。
          SAFE SIDE  只报「不是同一个」；从不认证「这就是我要的量」。

        `tol` 必须显式给，且必须 > 0 —— 容差为 0 的等式检查在浮点上恒假，
        而恒假的检查与恒真的检查一样没用（`#754` 那一族的镜像）。"""
        assert isinstance(what, str) and what, "说明这两个量为什么该相等（#761）"
        assert tol is not None and tol > 0, "容差必须显式且为正（#761）"
        if self._degenerate(name, observed, expected): return False
        # ⚠ #773:`_degenerate` 只认「全零」,于是 `identity_control(1.0, 1.0)` 这种
        #   **把常数与自己比**的空洞检查会一路 PASS —— `#772` 就是这么溜过去的。
        #   ⇒ 在这里补:两侧是**同一个字面常数**且调用者没给出可变来源时,记为空洞。
        #   只标注不阻断(与本库其余部分同一姿态),但它会出现在行里,让人看见。
        _vacuous = (observed == expected) and (tol is not None) and (abs(observed) not in (0.0,))
        d = abs(observed - expected)
        # ⚠ #761 的回测当场抓到：0.2797-0.2747 = 0.005000000000000004 > 0.005，
        #   于是「恰在容差上」被判 FAIL。**造来抓等式错误的闸，自己栽在浮点等式上。**
        #   ⇒ 容差比较必须带相对松弛，否则边界上的判定由浮点表示决定，而不是由容差决定。
        ok = d <= tol * (1 + 1e-9) + 1e-15
        self.rows.append((name, f"|{observed:+.5f} - {expected:+.5f}| = {d:.6f} <= {tol:g}",
                          ok, f"{what}" + ("" if ok else "  ⚠ 超出容差 {:.1f}× ⇒ 这两个量不是同一个".format(d/tol))
                          + ("  ⚠⚠ 两侧完全相等 —— 若它们是同一个字面常数,这条检查是空的 (#773)" if _vacuous else "")))
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

    def has_error_bar(self, name, value, spread, spread_source):
        """#167b:一个**地板**不是一个**误差棒**,而它们长得一模一样。

        两者都写成 `点值 ± 某个数`,都由同一段代码产出,都让读者觉得这个数被界住了。
        但零假设地板答的是**存在性**(「这不可能是无结构造出来的」),
        误差棒答的是**精度**(「换一批人还会是这个数吗」)。
        `#100` 的 0.432 挂着 curveball 地板 −0.022 走了八轮,而它的真实臂在 8 个种子上
        逐字节相同 —— 它从来没有过误差棒。

        所以调用者必须**命名展布的来源**,而不是只递一个数字:
            'bootstrap_人层'  重抽人 —— 唯一能答"换一批人还会是这个数吗"的那一种
            'split_跨劈分'    换劈分 —— 答"换一种分法还会是这个数吗"
            'seed_跨种子'     换种子 —— 只在种子真的驱动真实臂时才有意义
            'null_零臂'       **直接 FAIL** —— 这是地板,不是误差棒
            'analytic_解析'   公式(如相关系数的 delta 法)—— 记下公式名
        """
        ok = {'bootstrap_人层','split_跨劈分','seed_跨种子','analytic_解析'}
        allowed = sorted(ok | {'null_零臂'})
        if spread_source == 'null_零臂':
            r = (name, f"来源=零臂, {value:+.4f}", False,
                 "展布来源是**零臂** —— 那是地板不是误差棒(#167b),精度仍然未知")
        elif spread_source not in ok:
            r = (name, f"来源='{spread_source}'", False, f"不在名录里:{allowed}")
        elif not (spread > 0):
            r = (name, f"来源={spread_source}, 展布={spread}", False,
                 "零展布 —— 真实臂没有抖动来源,这不是误差棒")
        else:
            r = (name, f"{value:+.4f} ± {spread:.4f}", True, f"来源 {spread_source}")
        self.rows.append(r)
        return r[2]

    def count_needs_interval(self, name, n_pass, n_total, spread, spread_source, n_resamples=None,
                             seed_spread=None):
        """#232d:一个**计数**看起来比一个相关更硬,所以它更容易被当成确定值写出去。

        「31 个结局里越过全族阈值的有 19 个」读起来像在数东西 —— 数出来的数怎么会有误差棒?
        但那个阈值是**重抽样估出来的**:`#232d` 实测,同一个分数、同一次运行里,
        越阈个数两次算出 **10 和 8**,只因为最大统计量的零抽样取到了不同的随机流。
        本项目此前把 `N/M` 写成确定值的地方(`#223a` `#230a` `#270` 与两个 README)
        都因此高估了自己的精度。

        与 `#167b` 的 `has_error_bar` 完全同构,但它管的是**计数**而不是**点值**,
        因为计数的展布有一个额外的、专属的来源:**阈值本身的抽样**。

            'threshold_resample_阈值重抽样'  重抽零 -> 重估阈值 -> 重数(这一族的默认)
            'bootstrap_人层'                 重抽人 —— 同时动了相关与阈值
            'seed_跨种子'                    换种子重跑整条管道
            'null_零臂'                      **直接 FAIL** —— 地板不是误差棒(#167b 同款)
            'analytic_解析'                  公式
        """
        ok = {'threshold_resample_阈值重抽样','bootstrap_人层','seed_跨种子','analytic_解析'}
        allowed = sorted(ok | {'null_零臂'})
        base = f"{n_pass}/{n_total}"
        # #255c:`#233` 只量了**阈值重抽样**那一层,而 `#300` 实测**数据劈分种子**那一层
        #   可以比它大 40%(±2.5 vs ±1.8),并且正是它让 `#254a` 得出了一个错的结论。
        #   ⚠ 我注册的做法是「降级为 WARN」,这里改成 **FAIL** —— 严于注册,方向安全:
        #   一个未被测量的**主导**不确定性来源,等于精度未知,而不是精度尚可。
        if spread_source == 'threshold_resample_阈值重抽样' and seed_spread is None:
            r = (name, f"{base} ±{spread if spread is not None else '?'}(仅阈值重抽样)", False,
                 "**这一族的主导不确定性来源是【数据劈分种子】,不是阈值**(#255c) —— "
                 "未给 `seed_spread`,精度未知;`#300` 实测种子层可比阈值层大 40%")
            self.rows.append(r); return False
        if seed_spread is not None and spread is not None and seed_spread > 0:
            r = (name, f"**{base} 应读作 {n_pass-2*seed_spread:.0f}–{n_pass+2*seed_spread:.0f}**"
                       f"(种子 ±{seed_spread:.2f} · 阈值 ±{spread:.2f})", True,
                 "两层不确定性都已量:取较大的那一层作为区间(#255c)")
            self.rows.append(r); return True
        if spread_source == 'null_零臂':
            r = (name, f"{base} 来源=零臂", False,
                 "展布来源是**零臂** —— 那是地板不是误差棒(#167b/#232d),精度仍然未知")
        elif spread_source not in ok:
            r = (name, f"{base} 来源='{spread_source}'", False, f"不在名录里:{allowed}")
        elif spread is None:
            r = (name, f"{base} 展布=None", False,
                 "没有展布 —— 阈值没有被重抽过(#232d)")
        elif spread == 0 and (n_resamples or 0) >= 10:
            # #233b:守卫自己的正对照在**第一次使用**时就报了一个假阳。
            # 纯噪声分数的计数钉死在 0/31,展布**真的是 0**,而 25 次重抽都给了同一个数。
            # 「展布为 0」有两个完全不同的成因:**没重抽**(未知),与**重抽了而它不动**(已知且窄)。
            # 把两者折成同一个 FAIL,就是 P6 的三值塌成两值 —— 一个已知的窄区间被当成未知。
            r = (name, f"**{base} 钉死**(重抽 {n_resamples} 次,展布恰为 0)", True,
                 f"来源 {spread_source};计数在边界上被钉住 —— 这是**已测得的窄**,不是未测")
        elif not (spread > 0):
            r = (name, f"{base} 展布={spread}", False,
                 "负展布或未重抽 —— 一个计数没有抖动来源(#232d)")
        else:
            lo, hi = n_pass - 2*spread, n_pass + 2*spread
            r = (name, f"**{base} 应读作 {lo:.0f}–{hi:.0f}** (±{spread:.2f})", True,
                 f"来源 {spread_source};一个计数不是一个确定值")
        self.rows.append(r)
        return r[2]

    def control_kept_the_sample(self, name, before, after, n_before, n_after,
                                before_common=None, after_common=None, n_common=None,
                                tol_ratio=2.0):
        """#239a:前 11 个守卫全在防**假肯定**;这一个防**假否定**,而它的成本更高。

        一个改变了纳入条件的控制,改变的是**估计目标本身**。
        `#284` 实测:同一个残差化控制,在**各自的可用样本**上算,效应只剩 **2.8%**;
        在**同一批人**上算,保留 **70.2%**。那个 2.8% 完全是样本从 2,740 变成 6,301 造成的。

        > **「控制之后效应消失了」与「控制之后我在量另一批人」长得一模一样。**

        而前者会让人做一次**撤回** —— 撤回是永久的,因为**没人会重审一条自己撤掉的主张**。

        用法:
            g.control_kept_the_sample('审查控制', before=.1532, after=-.0741,
                                      n_before=9944, n_after=2806,
                                      before_common=.0830, after_common=-.0741, n_common=2806)
        规则:
            n_before == n_after                      -> PASS(纳入没变)
            n 变了但没给交集样本上的重比               -> **FAIL**(这个比较不可读)
            两个保留率相差超过 tol_ratio 倍            -> **FAIL**,并指出样本变了
        """
        if n_before == n_after:
            r = (name, f"n 未变({n_before:,})", True, "纳入条件没被控制改动")
        elif before_common is None or after_common is None:
            r = (name, f"n {n_before:,} -> {n_after:,},**未给交集样本重比**", False,
                 "改变纳入的控制必须在交集样本上再报一次(#239a);否则读不出是效应变了还是人变了")
        else:
            own = (after - before) / before if before else float('inf')
            com = (after_common - before_common) / before_common if before_common else float('inf')
            # #244b:守卫 12 第二次使用就在**近零退化**上报了一个假阳。
            #   `#289` 实测 own = −0.4%、com = −3.7% —— 两个都几乎是零,而**两个近零数的比值由噪声主导**,
            #   于是「相差 9 倍」被读成「样本变了」。真实情况是**控制两边都没动**。
            #   判据因此分两段:先问「控制到底动了没有」,只有动了才比两个样本上的动幅。
            if abs(own) < 0.10 and abs(com) < 0.10:
                r = (name,
                     f"各自样本 {100*own:+.1f}% vs 交集样本(n={n_common:,}) {100*com:+.1f}%"
                     f" —— **两边都 <10%,控制没动它**",
                     True, "近零区:比值不稳,按「控制无实质影响」读(#244b)")
                self.rows.append(r); return True
            bad = (abs(own) > tol_ratio*abs(com)) or (abs(com) > tol_ratio*abs(own)) or (own*com < 0)
            r = (name,
                 f"各自样本 {100*own:+.1f}% vs 交集样本(n={n_common:,}) {100*com:+.1f}%",
                 not bad,
                 ("两者一致,可以按效应读" if not bad else
                  "**两者相差超过 %gx —— 样本变了,不是效应变了(#239a)**" % tol_ratio))
        self.rows.append(r)
        return r[2]

    def plant_direction_from_sweep(self, name, sweep, baseline, baseline_spread=None,
                                   half_of=None, identity=None, identity_floor=0.90):
        """#248c:预注册一个**方向**,和预注册一个**阈值**,不是同一件事。

        阈值我算得出来;**方向我算不出来,我是在猜** —— 而我已经猜错四次:
        `#132b`(审查)· `#134f`(剂量)· `#146e`(罕见度)· `#248c`(跨半种入)。
        **四次都是扫描纠正了我,没有一次是我读表读对的。**

        所以这个守卫**不接受调用者给的方向**。它只接受一条扫描
        `sweep = [(g0, stat0), (g1, stat1), ...]`(`g` 递增,`g0` 应为 0),自己判:
            ① `g=0` 必须落在 `baseline` 上(若给了 `baseline_spread`,按 2× 展布判)
            ② 扫描必须**单调**(方向由数据给出,不由调用者给出)
            ③ 报出**灵敏度**:最小的 `g` 使 |stat − stat(0)| 超过 `half_of`
               (默认 = |baseline| 的一半)——「这个设计能看见多小的效应」

        **调用者只提供扫描,不提供期望 —— 于是「方向写反」在结构上不可能发生。**
        """
        if len(sweep) < 3:
            r = (name, f"只有 {len(sweep)} 个点", False, "扫描至少要 3 点才谈得上单调")
            self.rows.append(r); return False
        gs = [float(g) for g, _ in sweep]; ys = [float(y) for _, y in sweep]
        # #252c:一个由**数据定义**的量(特征向量、聚类中心、主成分),
        #   种得够强时**种入本身会重新定义它** —— 于是扫描中途换了被测对象。
        #   `#297` 实测:c1 的扫描非单调,不是噪声,是 c1 在 g 变大时变成了「被种入的那个方向」。
        #   `identity[i]` = 第 i 个 g 上被测量与 `g=0` 时的**身份指纹**(1 = 同一个对象)。
        #   有身份漂移时,**先报漂移**:它比「非单调」具体得多,而且单调也救不了它。
        if identity is not None:
            idv = [float(x) for x in identity]
            drift = next((gs[i] for i in range(1, len(idv)) if idv[i] < identity_floor), None)
            if drift is not None:
                r = (name,
                     f"**被测对象在 g={drift:g} 处被替换**(身份指纹 {idv[0]:.2f} → {min(idv):.2f}"
                     f" < {identity_floor})",
                     False,
                     "由数据定义的量,种入参与了它的重新定义 —— 这条扫描测的不是同一个东西(#252c)")
                self.rows.append(r); return False
        if baseline_spread is not None and baseline_spread > 0:
            ok0 = abs(ys[0] - baseline) <= 2*baseline_spread
        else:
            ok0 = abs(ys[0] - baseline) <= max(0.05*max(abs(baseline), 1e-9), 0.02)
        up = all(ys[i] <= ys[i+1] + 1e-9 for i in range(len(ys)-1))
        dn = all(ys[i] >= ys[i+1] - 1e-9 for i in range(len(ys)-1))
        thr = half_of if half_of is not None else abs(baseline)/2
        sens = next((gs[i] for i in range(1, len(ys)) if abs(ys[i]-ys[0]) > thr), None)
        direction = '上升' if (up and not dn) else ('下降' if (dn and not up) else '平/非单调')
        ok = ok0 and (up or dn) and (sens is not None)
        why = []
        if not ok0: why.append(f"g=0 未落在基线上({ys[0]:+.4f} vs {baseline:+.4f})")
        if not (up or dn): why.append("扫描**非单调** —— 种入很可能干扰了它本该保持不变的变量(#248c①)")
        if sens is None: why.append(f"扫描全程未越过 {thr:.4f} —— **灵敏度未证明**")
        r = (name,
             f"方向由数据给出:**{direction}**;灵敏度 g={sens if sens is not None else '未达到'};"
             f"g=0 {ys[0]:+.4f} · g={gs[-1]:g} {ys[-1]:+.4f}",
             ok,
             "扫描自证:单调 + g=0 落基线 + 灵敏度达到" if ok else ' · '.join(why))
        self.rows.append(r)
        return ok

    def could_have_come_out_otherwise(self, name, fn, perturbations, tol=1e-12):
        """#267b/#267c:一个由**代数或构造**决定、不可能变的数,被当成了测量结果。

        本项目此前有 13 个守卫,**没有一个问「这个数有没有可能是别的值」**。
        `#312` 同一轮里出现了两个:
          ③ 把每个预测量除以 √信度再算 R² —— **线性重缩放,R² 对它不变**,
            于是 `1.31% → 1.31%` 被报成一个干净的零,而它是一个**恒等式**;
          位置分 S 的「半块信度」恒为 **1.000 ± 0.000** —— 因为 S 根本不随块子集变化,
            **一个恒为 1 的信度不是高信度,是没有被测量。**

        > **`+0.00` 与一个干净的零长得一模一样。**

        用法:`fn` 接受一个扰动标签并返回一个数;`perturbations` 是**不该让结论成立、
        但应当让这个数动一动**的一组扰动(换种子、打乱输入、改一个无关参数)。
        守卫**实际执行**它们并检查输出**是否真的变了**。

            全部扰动下逐位相同(在 tol 内)-> **FAIL**「这个数不可能是别的值」
            至少一个扰动让它变了              -> 放行,并报出**最小与最大**的变动幅度
        """
        vals = []
        for p_ in perturbations:
            try: vals.append(float(fn(p_)))
            except Exception as e:
                vals.append(float('nan'))
        base = vals[0] if vals else float('nan')
        deltas = [abs(v - base) for v in vals[1:] if v == v]
        moved = [dd for dd in deltas if dd > tol]
        if not deltas:
            r = (name, "扰动全部报错或只给了一个扰动", False, "至少要有两个能跑通的扰动")
        elif not moved:
            r = (name, f"**{len(deltas)} 个扰动下逐位相同**(base={base:+.6g})", False,
                 "**这个数不可能是别的值** —— 它由代数或构造决定,不由数据决定(#267b)")
        else:
            r = (name, f"{len(moved)}/{len(deltas)} 个扰动让它变了;"
                       f"变动 {min(moved):.3g} … {max(moved):.3g}(base={base:+.6g})", True,
                 "它有可能是别的值 —— 可以当成测量结果读")
        self.rows.append(r)
        return r[2]

    def curve_has_enough_points(self, name, xs, ys=None, min_points=3, what=""):
        """#362a:**一条只有一个点的曲线,会让下游每一个判据都通过。**

        `#361b`:`R405` 的设计在第一步就死了(950 个人里只有 1 个块是他们都答过的),
        于是 `k` 只能取 1,曲线只有**一个点**,全距 **0.0000**。
        而注册的 kill(「全距 > MDE?」)判**「平」**,guard 21(三件套齐)判**「这个零可以发布」**——
        **两个判据都对着一个点开火,而没有一个问「曲线有几个点」。**

        > **`#296b` 挡的是一个**数**越出定义域;这一条挡的是一个**设计**退化成一个点。**
        > **全距、单调性、剂量-反应、规格曲线 —— 每一个都以「有一条曲线」为前提,而那个前提要先被检。**

        用法:任何以**曲线 / 扫描 / 剂量-反应 / 规格曲线**作证据的判据之前调用它,
        传入 x(扫描的自变量)与可选的 y。
        **不同 x 的个数 < min_points -> FAIL(设计退化)**,不要让下游判据去开火。
        """
        try:
            xv = [float(x) for x in xs]
        except Exception:
            r = (name, f"{xs!r} 不是可比较的 x 列表", False, "曲线的自变量必须给出(#362a)")
            self.rows.append(r); return False
        if ys is not None and len(ys) != len(xv):
            r = (name, f"x 有 {len(xv)} 个而 y 有 {len(ys)} 个", False, "长度不一致 —— 这不是一条曲线(#362a)")
            self.rows.append(r); return False
        uniq = len({round(x, 12) for x in xv})
        if uniq < int(min_points):
            r = (name, f"不同 x 只有 **{uniq}** 个(共 {len(xv)} 点,要求 ≥ {min_points})", False,
                 f"**设计退化 —— {what or '这条曲线'}不足以支撑任何全距/单调性/剂量-反应判据**;"
                 f"`#361b` 里两个判据都在一个点上通过了(#362a)")
        else:
            r = (name, f"不同 x **{uniq}** 个(共 {len(xv)} 点)", True,
                 "曲线有足够的点 —— 下游的全距 / 单调性判据可以读")
        self.rows.append(r); return r[2]

    def null_claim_uses_null_criteria(self, name, claim_kind, perm_quantile=None,
                                      mde=None, sensitivity_shown=None, meaningful=None):
        """#312a:**为「有效应」设计的判据,用在「没有效应」的结论上会系统性地误报。**

        本项目两次撞上同一件事,而两次我都是**事后用散文解释**为什么那个 FAIL 不算数:
            `#308c`  规格曲线「符号一致」FAIL —— 而**真零的签名就是符号乱走**
            `#311d`  `negative_control` FAIL(零大于效应)—— 而**结论本身就是一个零**
        **散文不是执行**(P9)。所以:一轮的结论是 `NULL` 时,它必须自带**零的判据**。

        ⚠ **边界(`#314c`,第一次实战当场暴露的)**:本守卫检的是**证据到没到齐**,
        **不检结论的方向**。`#314a` 的分位数是 **0.067** —— 偏向「有效应」那一侧 ——
        而三件套齐全,于是本守卫 PASS。**照那个 PASS 写成「没有效应」就是误读。**
        **一个「零的三件套齐全」的判决,不等于「这是一个零」。**
        分位数落在 [0.05, 0.95] 之外时,先问的应当是「这还是不是一个零式结论」。

        `claim_kind='EFFECT'` -> 本守卫不干预,直接 PASS(它只管 NULL 那一支)。
        `claim_kind='NULL'`  -> 三样必须同时在场,三缺一 FAIL:
            ① **perm_quantile** —— 置换零里 ≥ 观测的比例(读「观测落在零的哪里」,
               而不是「零相对效应多大」——后者在没有效应时没有意义);
            ② **mde**(或 CI 宽度)—— 这个设计**能**看见的最小效应;
            ③ **sensitivity_shown** —— 正对照**实际**证明过的灵敏度(不是算出来的,是种进去看见的)。
        再加一条内容性检查:给了 `meaningful`(一个「有意义的效应量」)时,
        **`mde` 必须小于它** —— 否则这个零没有内容,只是「我看不见」。
        `#310a` 正是这样:MDE 55% 而一个 30% 的缓冲有意义 -> 那个零**不可发布**。
        """
        k = str(claim_kind).upper()
        if k == 'EFFECT':
            r = (name, "claim_kind=EFFECT", True, "本守卫只管 NULL 那一支 —— 不干预")
        elif k != 'NULL':
            r = (name, f"claim_kind={claim_kind!r}", False,
                 "必须显式声明 'EFFECT' 或 'NULL' —— 不声明就没有判据(#312a)")
        else:
            miss = [n for n, v in (('置换分位数', perm_quantile), ('MDE/CI 宽度', mde),
                                   ('正对照证明的灵敏度', sensitivity_shown)) if v is None]
            if miss:
                r = (name, f"NULL 结论缺少:{' · '.join(miss)}", False,
                     "零的三件套必须同时在场:置换分位数 · MDE · 正对照灵敏度(#312a)")
            elif meaningful is not None and float(mde) >= float(meaningful):
                r = (name, f"MDE **{mde}** ≥ 有意义的效应量 **{meaningful}**", False,
                     f"**这个零没有内容 —— 它只说明我看不见**;`#310a` 的 MDE 55% 对一个 30% 的"
                     f"缓冲就是这样,那个零不可发布(#312a)")
            else:
                r = (name, f"分位数 {perm_quantile} · MDE {mde} · 灵敏度 {sensitivity_shown}"
                           + (f" · 有意义量 {meaningful}" if meaningful is not None else ""),
                     True, "零的三件套齐全,且 MDE 小于有意义的效应量 —— 这个零可以发布")
        self.rows.append(r); return r[2]

    def sign_flip_needs_direction_change(self, name, cos_with_ref, corr_ref, corr_new,
                                         cos_floor=0.50):
        """#306b:**方向没变,派生量却翻号 —— 这两件事不可能同时为真。翻的是特征向量的符号。**

        `#306a` 第一版:门槛 400 下 `c3` 与发布版的 **|cos| = 0.9466**(方向 95% 一致),
        而 `c3 ↔ 羞耻` 从 **−0.1278 翻成 +0.1218**。印出来是「跨门槛符号翻转,极差 0.2497」——
        一个看起来极强的发现。**而它通过了全部门。**
        真相:`np.linalg.eigh` 的特征向量**符号是任意的**,`|cos|` 又取了绝对值,
        所以「方向一致」与「投影反号」被同时报了出来。对齐符号后极差是 **0.0425**,
        而在方向还认得出来的区间(|cos| ≥ 0.95)只有 **0.0060**。

        > **|cos| 高而派生相关翻号,是「符号没对齐」的签名,不是一个发现。**
        > 这是本项目第三次撞上特征向量符号(`R210:73` · 这里)。

        用法:任何用**特征向量/主成分/SVD 分量**算出来的派生量,在与参照版本比较时,
        把 `|cos|` 与两个版本的派生量一起交给它。
        `|cos| > cos_floor` 而两者异号 -> FAIL。
        """
        try:
            c = abs(float(cos_with_ref)); a = float(corr_ref); b = float(corr_new)
        except Exception:
            r = (name, "参数不是数", False, "需要 |cos| 与两个版本的派生量(#306b)")
            self.rows.append(r); return False
        if c > cos_floor and (a * b) < 0:
            r = (name, f"|cos| = **{c:.4f}** 而派生量 {a:+.4f} -> {b:+.4f}(**异号**)", False,
                 f"**方向 {c:.0%} 一致却翻号 —— 这是特征向量符号没对齐的签名,不是发现**;"
                 f"先把 V[:,k] 对齐到参照版再重读(#306b)")
        elif c <= cos_floor:
            r = (name, f"|cos| = {c:.4f} ≤ {cos_floor}", True,
                 "方向本身已经不一致 —— 翻号可以是真的,但这时该报的是方向变了")
        else:
            r = (name, f"|cos| = {c:.4f},派生量 {a:+.4f} -> {b:+.4f}(同号)", True,
                 "方向一致且同号 —— 没有符号问题")
        self.rows.append(r); return r[2]

    def heldout_drop_needs_a_plant(self, name, observed_drop, plant_drop, what=""):
        """#303b:**一个「留出后掉了多少」的读数,分不清「声明有乐观」和「潜变量不稳」——
        只有把一个**已知为真**的效应种进去、看它掉多少,才分得清。**

        `#303a`:`c3` 与羞耻的相关从样本内 −0.1278 掉到嵌套 −0.0919,**掉了 28%**。
        看上去像「样本内那个数偏乐观」。**但把一个真值 +0.25 的相关种在同一个 `c3` 上,
        同一套嵌套程序把它掉了 39%**(真值 +0.10 时掉 57%)。
        `c3` 的 |r| = 0.128 落在两者之间,**预期衰减约 50%,而实测只有 28%** ——
        观测到的下降**全部**可由 `c3` 这个估计量自身的噪声解释,**不构成乐观的证据**。

        > **留出读数对一个由数据估出来的潜变量,量的是「这个潜变量重不重现」,
        > 而不是「这个声明真不真」。种植是唯一能把两者分开的东西。**

        规则:只有 `observed_drop > plant_drop` 才允许把下降读成**声明的**;
        否则下降是**仪器的**,而原值不必被降级(但必须注明它是**对这一份样本的潜变量**报的)。
        缺 `plant_drop` -> 直接 FAIL:没有种植,这个读数不可解释。
        """
        if plant_drop is None:
            r = (name, "缺少种植对照的衰减", False,
                 f"没有种植就分不清「{what or '这个量'}有乐观」和「潜变量不稳」(#303b)")
        elif observed_drop is None:
            r = (name, "缺少观测到的衰减", False, "两个都必须给出(#303b)")
        else:
            o, pd_ = float(observed_drop), float(plant_drop)
            ok = o > pd_
            r = (name, f"观测衰减 **{100*o:.0f}%** vs 种植真效应的衰减 **{100*pd_:.0f}%**", ok,
                 ("观测衰减超过种植的 —— 超出的那部分才可以读成声明的乐观"
                  if ok else
                  f"**观测衰减不超过种植的 -> 下降是仪器的,不是声明的**;"
                  f"原值不必降级,但必须注明它是**对这一份样本的潜变量**报的(#303b)"))
        self.rows.append(r); return r[2]

    def apply_reached_the_test_set(self, name, finite_counts, n_train, n_all, labels=None):
        """#302b:**在任何「训练估计 / 应用到测试」的重构里,应用那一步最常见的失败是它压根没发生。**

        `#302a` 的嵌套 CV 把训练掩码加在了**留一块剖面**上 —— 而那是一个**人内**量
        (这个人自己的其它块),根本不是从别人身上估的。加掩码没有防住任何泄漏,
        它只是把测试集的剖面整片抹成了 NaN。结果:四个坐标的有限数**恰好等于 |训练集|**,
        `nested_r2` 对 29 个结局全返回 NaN,`DataFrame` 连列都没有。

        > **告密者是可机检的:一个「应用到所有人」的量,它的有限数等于 |训练集|。**
        > 相等不是巧合 —— 它是「应用这一步没有发生」的签名。

        用法:把每个坐标/特征的有限计数交给它,并给出 |训练集| 与总人数。
        任何一个等于 |训练集| -> FAIL。(等于总人数是正常的;介于两者之间也正常。)
        """
        try: cs = [int(c) for c in finite_counts]
        except Exception:
            r = (name, f"{finite_counts!r} 不是计数列表", False, "必须给出每个应用量的有限计数(#302b)")
            self.rows.append(r); return False
        if not cs:
            r = (name, "空的计数列表", False, "没有量被检查 = 没有检查(#302b)")
        else:
            hit = [i for i, c in enumerate(cs) if c == int(n_train)]
            nm = (lambda i: labels[i] if labels and i < len(labels) else f"#{i}")
            if hit:
                r = (name, f"{[nm(i) for i in hit]} 的有限数 = |训练集| {n_train}", False,
                     f"**应用那一步没有发生** —— 这些量只在训练集上有值;"
                     f"检查是不是把训练掩码加在了一个**人内**量上(#302b)")
            else:
                r = (name, f"有限计数 {cs} vs |训练集| {n_train} / 全体 {n_all}", True,
                     "没有一个等于 |训练集| —— 应用步确实跑到了测试集")
        self.rows.append(r); return r[2]

    def bounded_statistic_out_of_range(self, name, value, lo, hi, what=""):
        """#296b:**一个有定义域的统计量落在定义域外,是仪器坏了,不是效果极强。**

        `#296a` 的正对照算「收回比例」= 偏出后的掉幅 ÷ 原样的掉幅。两个掉幅**都是负的**,
        而我写的是 `rec = (hi_par-base_par) / max(hi_raw-base_raw, 1e-9)` ——
        `max(...)` 是给正量写的,它把 **−0.3614 夹成了 1e-9**,于是 rec = −4.3e8,
        「收回」印成 **42842376318%**,而门槛 `(1-rec) > 0.40` **欣然通过**。

        > **单侧门槛分不清「非常好」和「坏掉了」。**
        > 真实的收回是 **−18.5%** —— 偏出让人为的信度差**更大**,正对照其实是 FAIL。

        所以:任何**有天然上下界**的统计量(比例、份额、收回率、相关、保留百分比)
        在进入任何判据之前,**必须先被它自己的定义域挡一道**。
        越界 -> 仪器失败(UNVERIFIED),**永远不读成极端的成功**。
        """
        try: v = float(value)
        except Exception:
            r = (name, f"{value!r} 不是数", False, "有界统计量必须是数(#296b)"); self.rows.append(r); return False
        if not np.isfinite(v):
            r = (name, f"{v}", False, "非有限值 —— 仪器失败,不是极端效果(#296b)")
        elif v < lo or v > hi:
            r = (name, f"**{v:.4g}** 落在 [{lo}, {hi}] 之外", False,
                 f"**{what or '这个量'}越界 -> 仪器坏了,不是效果极强**;"
                 f"`#296a` 印出 4.28e10% 并通过了一个 `>0.40` 的单侧门槛(#296b)")
        else:
            r = (name, f"{v:.4g} ∈ [{lo}, {hi}]", True, "在定义域内 —— 可以进入判据")
        self.rows.append(r); return r[2]

    def component_difference_is_not_mechanism(self, name, whole_before, whole_after,
                                              spread, component, retain_ceiling=0.50):
        """#294b:**一个零件在两组间不同,不等于它就是整体差异的机制。**

        `#293a` 测到两组的 `rar`(什么算冷门)在块内流行度上明显不同 —— **148.6 个展布**,
        正对照完美。我据此在同一小时里往两个公开 README 写下「**而机制在尺子本身**」。
        `#294a` 把尺子换成**组内**的重测:整体差从 **−0.2061** 变成 **−0.2270**,
        **保留 110.2%** —— 修好那个零件,整体差**一点没动**。

        > **零件差是一个测量;机制是一个因果断言。**
        > 二者之间隔着一次**把零件换掉再测整体**的实验,而那次实验是可以做的、便宜的,
        > 并且**恰恰是我在注册 NEXT 时自己写下要做的那一步**。

        所以:任何形如「A 与 B 在 X 上不同,**因此** X 是 A/B 整体差异的机制」的句子,
        必须携带**修好 X 之后重测的整体差**。没有它 -> 这是 η 过大(frontier §2),
        不是发现。
        """
        if whole_before is None or whole_after is None:
            r = (name, "缺少「修好零件前/后的整体差」", False,
                 f"要把 `{component}` 叫作机制,必须把它换掉再测一次整体(#294b)")
        else:
            keep = abs(whole_after) / max(abs(whole_before), 1e-12)
            ok = abs(whole_after) < abs(whole_before) - (spread or 0.0) and keep <= retain_ceiling
            r = (name,
                 f"整体差 **{whole_before:+.4f}** -> 修好 `{component}` 后 **{whole_after:+.4f}**"
                 f"(保留 **{100*keep:.1f}%**)",
                 ok,
                 ("修好零件后整体差确实塌了 —— 可以按机制读"
                  if ok else
                  f"**修好 `{component}` 后整体差保留了 {100*keep:.1f}% —— 它不是机制,只是同时存在的另一个差异**;"
                  f"`#294a` 保留 110.2%(#294b)"))
        self.rows.append(r); return r[2]

    def profile_similarity_is_not_identity(self, name, profile_r, score_r,
                                           profile_floor=0.60, score_floor=0.30):
        """#279b:**剖面相似作为「是不是同一个东西」的证据,力量为零。**

        这个项目三次撞上同一个分离,而第三次是**代数上保证**的:
            `#259a`  分数 **+0.1589** / 剖面 **+0.7826**
            `#263a`  分数 **−0.1010** / 剖面 **+0.7091**   ← 符号都相反
            `#324`   分数 **恰好 0**(`form_i` 是六坐标残差,D 是六坐标之一)/ 剖面 **+0.7105**

        **两个相关恰好为零的量,可以在 27 个结局上给出几乎同一张脸。**
        所以任何以**剖面相关**下「同一构念」结论的地方,**必须同时给出分数层相关**;
        剖面高而分数低 -> 这两个量**不是**同一个东西,无论剖面多像。

        (反过来是允许的:剖面**低**是「不是同一个东西」的有效证据 ——
        这个守卫只挡「剖面高 ⇒ 同一」这一个方向。)
        """
        if profile_r is None or score_r is None:
            r = (name, "缺少剖面相关或分数相关", False, "两个都必须给出(#279b)")
        elif abs(profile_r) > profile_floor and abs(score_r) < score_floor:
            r = (name, f"剖面 **{profile_r:+.4f}** 而分数 **{score_r:+.4f}**", False,
                 f"**剖面高({profile_floor})而分数低({score_floor})—— 这不是同一个构念**;"
                 f"`#324` 有一个分数层恰好为零、剖面仍 +0.71 的代数反例(#279b)")
        elif abs(profile_r) <= profile_floor:
            r = (name, f"剖面 {profile_r:+.4f} ≤ {profile_floor}", True,
                 "剖面低 —— 「不是同一个东西」是有效结论,这个方向不受限制")
        else:
            r = (name, f"剖面 {profile_r:+.4f} 且分数 {score_r:+.4f}", True,
                 "剖面与分数同时高 —— 可以按「同一构念」读")
        self.rows.append(r)
        return r[2]

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

    def asserted(self, name, condition, detail, kind="kill"):
        """#96a: a condition stated in prose must be a boolean here, or it was never tested.

        `kind` (#366e, guard 23):
          "control" -- an instrument check (positive/negative/offset control).
          "kill"    -- the pre-registered threshold itself.
        Declaring it is what lets `__str__` tell **OVERTURNED** from **UNVERIFIED**.
        Rounds that never pass `kind` keep the old two-valued output byte-for-byte.
        """
        if kind not in ("kill", "control"):
            raise ValueError(f"kind must be 'kill' or 'control', got {kind!r}")
        ok = bool(condition)
        self.rows.append((name, detail, ok, "asserted in code, not in prose"))
        if kind == "control":
            self._control_rows = getattr(self, "_control_rows", set()) | {len(self.rows) - 1}
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

    def threshold_outside_noise(self, name, value, threshold, spread):
        """#142d: 一个预注册阈值,必须先证明它离开了被比较量**自身的噪声带**。

        `resolvable` 检查的是"效应 vs 零";这一条检查的是"效应 vs 阈值"。
        A12/R170 的预注册判定是 `r250 >= 3`,实测 3.0484 —— 高出 1.6%,而那个比值
        自身在 20 个重抽种子上是 3.023–3.163(sd 0.035)。**门槛落在噪声带里,
        所以那个判定的输出由重抽种子决定,不由数据决定**,而它印出来的话是
        「加强成功,按新强度引用」。

        通过条件:|value - threshold| > 2 * spread。否则判定是 UNVERIFIED,
        **不是**"未达标",也**不是**"达标"。"""
        if self._degenerate(name, value - threshold, spread):
            return False
        gap = abs(value - threshold)
        ok = gap > 2 * spread
        self.rows.append((name, f"|{value:+.4f} - {threshold:+.4f}| = {gap:.4f} vs 2*{spread:.4f}", ok,
                          f"{gap/max(spread,1e-12):.1f}x 自身噪声" if ok else
                          f"阈值落在噪声带内({gap/max(spread,1e-12):.1f}x)—— 判定由重抽种子决定,UNVERIFIED"))
        return ok

    def equivalent_within(self, name, diff, spread, margin):
        """#150: 等价检验。`resolvable`/`require_resolvable_first` 是为"我要它非零"设计的;
        当假设本身是**"两者相同"**时,用它们会把想要的结果报成 FAIL,并且把整族标 MOOT。

        正确的形式是**等价界**(TOST 式):差的置信区间必须**落在**预先指定的边界内。
        A19R02 的判别量 `corr(z,S)+corr(ρ,S)` = −0.0034,展布 0.0142 —— 0.2× 看似完美,
        但它的 95% 上界是 0.0318,而效应本身只有 0.037,**所以这个设计只能排除
        大于效应 80% 的差异**。不报这个边界,"两者相同"就是一句没有分辨率的话。

        通过条件:|diff| + 2*spread <= margin。**margin 必须在跑之前指定。**"""
        if self._degenerate(name, diff, spread):
            return False
        hi = abs(diff) + 2 * spread
        ok = hi <= margin
        self.rows.append((name, f"|{diff:+.4f}| + 2*{spread:.4f} = {hi:.4f} vs 边界 {margin:.4f}", ok,
                          f"差被界在 {hi:.4f} 内,小于预设边界" if ok else
                          f"只能把差界在 {hi:.4f},而边界是 {margin:.4f} —— "
                          f"设计的分辨率不够,'两者相同'未被证明"))
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


    def spec_curve_cells_declare_n(self, name, cells, what=""):
        """#518b: 一条规格曲线,若它的产物没有**逐格存 n**,事后就无法做人群审计。

        `#517` 只有在 `R534` 的 JSON 里有逐格 n 时才查得出「那四格对结局取了条件」;
        而 `#494a`/`#501a`/`#487d` 的产物 **0 格有 n** -> 它们**结构上不可审计,除非重跑**。
        ⇒ 规格曲线的每一格必须把 n 写进 results/。缺 n 的格,事后不可审。

        cells: 可迭代的 dict(每格一个),或 {key: dict}。逐格检查是否含非空的 'n'。
        """
        items = list(cells.values()) if isinstance(cells, dict) else list(cells)
        total = len(items)
        if total == 0:
            self.rows.append((name, "0 格", False, "DEGENERATE -- 规格曲线为空,这个检查不可判"))
            return False
        withn = [c for c in items if isinstance(c, dict) and c.get("n") is not None]
        ok = len(withn) == total
        self.rows.append((name, f"{len(withn)}/{total} 格带 n", ok,
                          ("每一格都声明了自己的人群规模" if ok else
                           f"⛔ {total-len(withn)} 格缺 n -> **这条曲线事后不可做人群审计**(#518b)")
                          + (f" [{what}]" if what else "")))
        return ok


    def spec_curve_cells_declare_inclusion(self, name, cells, what=""):
        """#519a: 三次改参照都判错,根因是**产物没记「这一格由哪些条件共同定义」**。

        `#518a` 用一个常数比两波 · `#519a` 容差太紧 · 最后发现真正的形态是
        「阳性数 ∩ 其余纳入条件」—— 而**交集里有哪些条件,产物里一个字也没有**。
        ⇒ 逐格 `n` 不够;每格还必须带一个 `inclusion`:**逐条列出该格的纳入条件**。

        判据:每格都有非空的 `inclusion`(list[str] 或非空 str)。
        ⚠ 这是一条**新**要求 —— 它在本项目**所有既有产物上都应当 FAIL**,
          而那正是它可失败的证明。
        """
        items = list(cells.values()) if isinstance(cells, dict) else list(cells)
        total = len(items)
        if total == 0:
            self.rows.append((name, "0 格", False, "DEGENERATE -- 规格曲线为空"))
            return False
        def has(c):
            v = c.get("inclusion") if isinstance(c, dict) else None
            return bool(v) and (isinstance(v, str) or (isinstance(v, (list, tuple)) and len(v) > 0))
        good = [c for c in items if has(c)]
        ok = len(good) == total
        self.rows.append((name, f"{len(good)}/{total} 格带 inclusion", ok,
                          ("每格都逐条声明了纳入条件" if ok else
                           f"⛔ {total-len(good)} 格缺 inclusion -> **这一格由哪些条件共同定义,产物里没有**(#519a)")
                          + (f" [{what}]" if what else "")))
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
        out.append(f"   => {self.three_valued()}")
        return "\n".join(out)

    def eigenvector_is_anchored(self, name, vec_scores, reference, ref_name="", min_abs=0.02):
        """#368a (guard 24) -- an eigenvector's SIGN is arbitrary; a label written off it is a coin flip.

        `np.linalg.eigh` / `svd` return a direction, not an orientation. Every quantity this project
        reads off a component -- R^2, commonality, retention, |cos| -- is sign-INVARIANT, so the bug
        never shows up in a gate. It shows up in the PROSE: "high PC1 = younger, fewer partners,
        more anxious" was written straight off an unanchored eigenvector and was exactly backwards
        (mean shame ran 1.179 at the low end down to 0.266 at the high end, i.e. the opposite of
        what the label claimed).

        Fourth occurrence in this project (R210:73, #306b, #361, here). Guard 20 catches a sign
        FLIP across thresholds; it cannot see a single run whose one orientation is simply wrong.

        Pass the component scores and the reference the label is written against. FAIL unless
        corr(scores, reference) > +min_abs -- i.e. the component has been oriented, not merely
        computed. A near-zero correlation FAILS too: an orientation that cannot be anchored is
        not an orientation, and its label must be dropped rather than guessed.
        """
        import numpy as _np
        v = _np.asarray(vec_scores, dtype=float); r = _np.asarray(reference, dtype=float)
        g = _np.isfinite(v) & _np.isfinite(r)
        if g.sum() < 30 or v[g].std() < 1e-12 or r[g].std() < 1e-12:
            ok, c = False, float("nan")
            note = "too few finite pairs / degenerate -- cannot anchor, so the label may not be written"
        else:
            c = float(_np.corrcoef(v[g], r[g])[0, 1])
            # #368b: an ABSOLUTE floor lets pure noise pass -- at n=500 the sampling se is ~0.045,
            # so min_abs=0.02 acquitted an unrelated component on the first attack pass.
            # The floor must be priced at its own level: 3 standard errors of a null correlation.
            need = max(float(min_abs), 3.0 / _np.sqrt(int(g.sum())))
            ok = c > need
            note = ("anchored" if ok else
                    ("SIGN IS BACKWARDS -- flip the component before writing any label"
                     if c < -need else
                     f"|corr| <= {need:.4f} (3 se at n={int(g.sum())}) -- unanchorable, "
                     "the label must be dropped, not guessed"))
        self.rows.append((name, f"corr(component, {ref_name or 'reference'}) = {c:+.4f}", ok, note))
        return ok

    def relaxation_reached_the_population(self, name, n_narrow, n_wide, min_growth=1.15, what=""):
        """#371a (guard 25) -- relaxing an inclusion rule that does not move n never happened.

        Sibling of guard 18 (`apply_reached_the_test_set`). There the tell was a finite-count equal
        to |train|; here it is **n barely moving**. R415 relaxed block coverage from >=8 to >=4 and
        the arm grew by 70 people (6,473 -> 6,543) instead of the ~1,700 the knob should deliver --
        because a control variable in the same mask (`S`) is itself only defined at coverage >=8,
        so `isfinite(S)` re-imposed the narrow rule downstream of the knob.

        A caliber arm that did not widen is not a replication in a different population; it is the
        SAME population reported under a wider label, which is a scope claim in the flattering
        direction. FAIL unless n_wide >= min_growth * n_narrow.
        """
        ok = float(n_wide) >= float(min_growth) * float(n_narrow)
        grew = (float(n_wide) / max(float(n_narrow), 1.0) - 1.0) * 100.0
        self.rows.append((name, f"n {int(n_narrow):,} -> {int(n_wide):,} ({grew:+.1f}%{', ' + what if what else ''})",
                          ok, "the relaxation reached the population" if ok else
                          f"RELAXATION BLOCKED -- needs >= {min_growth:g}x; something downstream "
                          "still enforces the narrow rule (a control variable defined only there?)"))
        return ok

    def _record_call(self, **kw):
        """#383b: record EVERY call, including the ones that could not be compared.

        The first version recorded only the comparable path, so a control that was
        UNCOMPARABLE vanished from the audit -- which is exactly the failure R426 named
        about its own regex scan: what cannot be extracted is not "fine", it is "unlooked-at".
        A ledger that drops its own hard cases flatters the thing it audits.
        """
        self.calls = getattr(self, "calls", [])
        self.calls.append(kw)
        try:
            import json, pathlib as _pl
            with _pl.Path("lib/_gate_calls.jsonl").open("a") as fh:
                fh.write(json.dumps(kw, ensure_ascii=False) + "\n")
        except Exception:
            pass          # a ledger that breaks a round is worse than a missing ledger

    def positive_control_at_the_contested_magnitude(self, name, plant_effect, contested_effect,
                                                    plant_passed, what="", branch=None,
                                                    main_quantity=None, sweep_detection=None):
        """(#402a) `branch` makes the commonest misuse impossible instead of merely remembered.

        Three times now the wrong `contested_effect` was passed (#372c, #384d, #402a): a NULL
        result's contested magnitude is the MEANINGFUL effect size (what you would have wanted to
        see), while a NON-NULL result's is the OBSERVED effect (prove the instrument works at the
        size actually found). Passing the null's number on a firing result makes the gate reject a
        real finding; passing the observed number on a null makes it wave through a blind design.

        `branch` must be "null" or "non_null" when given; it does not change the arithmetic, it
        forces the caller to say which quantity `contested_effect` is -- the #383a move of changing
        the interface rather than strengthening the reader.
        """
        if branch is not None and branch not in ("null", "non_null"):
            raise ValueError(f"branch must be 'null' or 'non_null', got {branch!r}")
        if main_quantity is not None and main_quantity not in ("continuous", "discrete_count"):
            raise ValueError("main_quantity must be 'continuous' or 'discrete_count', "
                             f"got {main_quantity!r}")
        # #407a: a DISCRETE main quantity (sign counts, cell counts, threshold-crossing counts)
        # has almost no dynamic range -- any uniform nudge flips every near-zero term at once and
        # the sweep saturates, so its "MDE" reads tiny when it is really UNMEASURABLE (#406b).
        # R451 tried to find such rounds retrospectively by regex and failed its own negative
        # control on the third try, because "what the main quantity is" is written nowhere in the
        # code. So the fix is #383a's: make the caller declare it, and refuse the calibration
        # outright when the sweep saturates.
        if main_quantity == "discrete_count":
            det = list(sweep_detection or [])
            if not det:
                self.rows.append((name, f"discrete main quantity, no sweep supplied", False,
                                  "a discrete count needs its sweep's detection rates so "
                                  "saturation can be ruled out -- MDE alone is not admissible"))
                return False
            if all(abs(float(x) - 1.0) < 1e-9 for x in det):
                self.rows.append((name, f"discrete main quantity, sweep detection {det}", False,
                                  "SWEEP SATURATED at every level -- this design's MDE is "
                                  "UNMEASURABLE, not small; the calibration carries no information "
                                  "(#406b)"))
                return False
            what = (what + " · " if what else "") + f"discrete, sweep {det} (not saturated)"
        if branch is not None:
            what = (f"[{branch}] " + what) if what else f"[{branch}]"
        return self._pcacm(name, plant_effect, contested_effect, plant_passed, what)

    def _pcacm(self, name, plant_effect, contested_effect,
               plant_passed, what=""):
        """#382a (guard 26) -- a positive control planted BIGGER than the thing in dispute is silence.

        P5* says a zero from an instrument that has never returned non-zero is silence, not an
        acquittal. R425 found the sharper version: in the same model and the same run, a SYNTHETIC
        positive control planted at 0.25 passed at |t| 6.82 while the REAL reference `EARLY` --
        known from #380 to be genuinely different -- reached only |t| 1.73 and failed. The synthetic
        control would have licensed "S shows no sex difference". **An instrument that returns
        non-zero only at magnitudes LARGER than the contested one is still silent.**

        So the failure this catches is not "the positive control failed". It is the far more
        dangerous **"the positive control passed, but it was planted above the contested
        magnitude"** -- an instrument that looks validated and is not. Checks fail toward PASS,
        and this is that failure wearing a control's clothes.

        plant_effect      magnitude actually planted (same units as contested_effect)
        contested_effect  magnitude of the thing being claimed absent/present
        plant_passed      did the planted control clear its threshold
        PASS requires: the control passed AND was planted at or below the contested magnitude.
        Either input missing/non-finite -> UNVERIFIED-style FAIL, never a silent pass.
        """
        try:
            pe = abs(float(plant_effect)); ce = abs(float(contested_effect))
        except (TypeError, ValueError):
            self.rows.append((name, "magnitudes unavailable", False,
                              "cannot compare -- an uncomparable control is UNVERIFIED, not a pass"))
            self._record_call(kind="positive_control", name=name, plant=None, contested=None,
                              passed=bool(plant_passed), ok=False, what=what,
                              status="UNCOMPARABLE_bad_input")
            return False
        if not (pe == pe and ce == ce) or ce <= 0:
            self.rows.append((name, f"plant {pe:.4g} vs contested {ce:.4g}", False,
                              "degenerate: contested magnitude is zero or non-finite -- "
                              "there is nothing for the control to be calibrated against"))
            self._record_call(kind="positive_control", name=name, plant=pe, contested=None,
                              passed=bool(plant_passed), ok=False, what=what,
                              status="UNCOMPARABLE_degenerate_contested")
            return False
        ok = bool(plant_passed) and pe <= ce
        if not plant_passed:
            note = "the control did not fire at all -- the instrument is silent (P5*)"
        elif pe > ce:
            note = (f"CONTROL PLANTED {pe/ce:.1f}x ABOVE THE CONTESTED MAGNITUDE -- it proves the "
                    "instrument works at a size nobody is arguing about, and says nothing about "
                    "the size in dispute (#381c)")
        else:
            note = "fired at or below the contested magnitude -- the instrument is calibrated where it matters"
        self.rows.append((name, f"plant {pe:.4g} vs contested {ce:.4g}"
                          + (f", {what}" if what else ""), ok, note))
        # #383a: the audit reads THIS, not the source code. R426 could machine-read only
        # 1.9% of 212 rounds because planted magnitudes are written a dozen different ways
        # in code. Writing a stronger regex is the path #374b and #377c both showed fails.
        # Recording the parameters at the call site makes extraction exact by construction.
        self._record_call(kind="positive_control", name=name, plant=pe, contested=ce,
                          passed=bool(plant_passed), ok=ok, what=what, status="COMPARABLE")
        return ok

    def three_valued(self):
        """#366e (guard 23) -- P6 says verdicts are CONFIRMED / OVERTURNED / UNVERIFIED.

        This class had **two** exits, so a pre-registered kill that FIRED (controls sound,
        threshold crossed against me) printed as `UNVERIFIED`. That is P6's error mirrored:
        folding UNVERIFIED into OVERTURNED manufactures false **pardons**; folding OVERTURNED
        into UNVERIFIED manufactures false **doubt** -- and doubt gets retried while a
        refutation gets abandoned. In a project whose basin rule is "same question UNVERIFIED
        twice -> change direction", that mislabel corrupts the basin detector itself.

        Back-compatible: with no row declared `kind="control"` the two old strings are returned.
        """
        ctrl = getattr(self, "_control_rows", set())
        if getattr(self, "_moot", False) or getattr(self, "_moot_fams", set()):
            return "UNVERIFIED, and that is not an acquittal"
        if not ctrl:
            return "ALL GATES PASS" if self.verdict() else "UNVERIFIED, and that is not an acquittal"
        bad_ctrl = [self.rows[i][0] for i in sorted(ctrl) if not self.rows[i][2]]
        if bad_ctrl:
            return ("UNVERIFIED, and that is not an acquittal "
                    f"-- the instrument failed its own control ({', '.join(bad_ctrl)})")
        bad_kill = [r[0] for i, r in enumerate(self.rows) if i not in ctrl and not r[2]]
        if not bad_kill:
            return "ALL GATES PASS"
        return ("OVERTURNED -- controls sound, and the pre-registered threshold fired "
                f"AGAINST the expectation ({', '.join(bad_kill)}). "
                "This is a refutation, NOT an unverified round (#366e).")


    # ------------------------------------------------------------ #477:通过了的 KILL 也要被审
    def passing_kill_audit(self, floors=None):
        """`#477b`:本项目对**零**很严(MDE·功率·正对照),对**通过了的 KILL** 却几乎不审。

        证据:本会话七条被我自己抓住的错门,**间隔全部为 0** —— 每一条都在写下它的同一轮被发现,
        因为它的错**当场就显形**(4/4 全过、判决与对照矛盾、一个不可能的数)。
        **⇒ 我抓得住的是「错得显眼」的门,抓不住「错得像真的」的门。**
        而「更宽」的门恰恰更容易错得像真的(它只是让本来该失败的东西通过),
        **所以已发布的主张里,带着一个未被发现的**宽**门的概率,高于带着一个未被发现的**严**门。**

        `floors`:`{门名: "什么样的结果会让它失败"}`。
        对每一个**通过了的 KILL**,若调用方没给出这句话,就喊出来 —— 与 `零需要 MDE` 对称。
        """
        floors = floors or {}
        ctl = getattr(self, "_control_rows", set())
        out = []
        for i, (name, detail, ok, _) in enumerate(self.rows):
            if i in ctl or not ok:
                continue
            out.append((name, floors.get(name)))
        miss = [n for n, f in out if not f]
        for n, f in out:
            print(f"   {'✅' if f else '⚠ '} 通过的 KILL:{n[:64]}"
                  + (f"\n        会让它失败的:{f}" if f else "  <- **没说什么会让它失败**"))
        if miss:
            print(f"   ⚠ **{len(miss)}/{len(out)} 个通过的 KILL 没有说明什么会让它失败** —— "
                  f"这与「一个零必须报 MDE」是同一条要求的另一半。")
        return len(out), len(miss)

    def value_range_guard(self, name, values, expect_lo, expect_hi, coding_note, what=""):
        """`#549`(`E02·A222·R594`)—— 把救过两次的那道基础率守卫,从每轮手写变成库函数。

        **由来:** `#495b` 我把 `sexsex5` 读错,阳性率出来是 0.5368 ≈ 男性占比;
        `#492` 的 `c_sex15` 是 n=15。两次都是**同一道手写 assert** 抓住的,
        而它每一轮都要重写一遍 —— 于是它总有一轮会被忘掉。
        `#548c` 之后这件事变得更急:**NSFG 里 83.2% 的 `Whether/Ever` 变量编码为 1 与 5**,
        按字面读成 0/1 会**反号**,而任何统计检查都不会响 —— 阳性率会是 0.20 而不是 0.80,
        **看起来完全合理**。

        判据:把一列被当作二元使用的值,算出基础率,**必须落在预注册的 `[lo, hi]` 内**。
        `coding_note` 必须写明**这一列的码是什么、怎么映射的** —— 空字符串直接判 FAIL,
        因为「我知道它的码」正是每一次读错时我以为的事。

        ⚠ 可靠方向(同 `#546c`/`#547b`):**通过不证明映射对**,只证明基础率没有明显荒谬。
          一个真正反号且基础率恰好对称(0.5)的列,这道守卫抓不到 —— 这一点写在失败信息里。
        """
        import numpy as _np
        v = _np.asarray(values, dtype=float)
        v = v[_np.isfinite(v)]
        n = int(v.size)
        if not str(coding_note).strip():
            self.rows.append((name, "无 coding_note", False,
                              "FAIL -- 未写明码与映射。「我知道它的码」正是每次读错时我以为的事"))
            return False
        if n == 0:
            self.rows.append((name, "n=0", False, "DEGENERATE -- 空列,这个检查不可判"))
            return False
        uniq = set(_np.unique(v).tolist())
        if not uniq <= {0.0, 1.0}:
            self.rows.append((name, f"取值 {sorted(uniq)[:6]}", False,
                              f"FAIL -- 被当作二元使用,实际取值不是 {{0,1}} -> **先映射再用**({coding_note})"))
            return False
        rate = float(v.mean())
        ok = expect_lo <= rate <= expect_hi
        detail = (f"基础率 {rate:.4f} ∈ [{expect_lo:.2f},{expect_hi:.2f}] (n={n}) [{coding_note}]"
                  if ok else
                  f"FAIL -- 基础率 {rate:.4f} 不在预注册区间 [{expect_lo:.2f},{expect_hi:.2f}] "
                  f"(n={n}) -> **码可能读反了**({coding_note})")
        if ok:
            detail += " ⚠ 通过不证明映射对:一个反号且基础率≈0.5 的列,本守卫抓不到"
        self.rows.append((name, f"rate={rate:.4f}", ok, detail))
        return ok


def calibrated_tolerance(null_samples, k=1.0, q=0.95):
    """#711/#712:容差必须从零里取,而不是每次手写一个数。

    **由同一个毛病第三次触发**:`#691` 的比值区间 `[0.5,2]`、`#710` 的分诊门槛 `1.5`、
    `#711` 的复现容差 `0.02` —— 三次都是我选了一个数当门槛,而三次都是那个数在决定裁决。

    参数
    ----
    null_samples : 该量自己的零分布样本(一维)
    k            : 余量系数,默认 1.0 = 恰好取零的分位,不再额外加宽
    q            : 分位,默认 0.95

    返回零分布的 `q` 分位 × `k`。**这是 `positive_control` 的 `spread`/`floor`、
    `offset_control` 的 `spread` 应当取的值。**

    ⚠ **本函数默认不被三个 control 调用,这是刻意的。**
    `#712` 量过:把三个 control 改接校准容差,会让本页面 **4 条**已有裁决翻,
    而 `#712` ⑤ 预注册的停止条件是「>2 就不合入」⇒ **停在这里。**
    ⚠ 而四条翻的方向**全部是「手写不过 → 校准通过」** ——
    **我的手写门槛系统性地比数据支持的更严,即我一直在少报,不是多报。**
    **但少报是安全的那一侧,而「改一个门槛就翻四条裁决」正是棘轮要防的事** ⇒
    **旧裁决一律不动,新轮次用本函数取容差。前向生效,不追溯。**
    """
    import numpy as _np
    a = _np.abs(_np.asarray(null_samples, dtype=float))
    a = a[_np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(_np.quantile(a, q) * k)
