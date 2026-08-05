"""页面算术核对器(`#376`,由 `E01·A140·R420` 抽出为常驻仪器)。

凡由恒等式相连、却分别印出的数,必须互相对得上。
**一个孤立的数没有对错可言;两个由恒等式相连的数,只要都印出来,就会互相检查。**

⚠ **区间算术,不是标量容差** —— 印出的数是四舍五入的,标量容差会把**排版**报成**错误**
   (第一版正是这样给了 6 条假阳性,而那一轮的 IMPOSSIBLE 栏预言过它)。
⚠ **覆盖率必须和结论一起报**:「0 条不一致」不等于「整页算术干净」。
⚠ **用前必须跑 `controls()`**(`P5★` + `#374b`):
   一个从未开火过的核对器,它的每一个「一致」都是沉默,不是无罪。
"""
import re, itertools

NUM=r'[-−+]?\d+(?:\.\d+)?'
def f(x): return float(str(x).replace('−','-').replace('+',''))
# ⚠ **第一版用一个标量容差,而它给了 6 条假阳性 —— 本轮 IMPOSSIBLE 栏预言的正是这个。**
#   印出的数是**四舍五入**的:`± 0.0034` 的真值在 [0.00335, 0.00345],
#   于是 (0.2219−0.0884)/sd 的真值范围是 **38.7–39.9**,而页面印的 **+38.9 在里面**。
#   一个标量容差把**排版**报成了**错误**,那是核对器自己的假阳性,而它的方向是**制造假警报**。
# 修法:**区间算术** —— 每个印出的数按它自己的位数展成区间,恒等式在区间上求值,
#   只有当**印出值的区间与恒等式的区间不相交**时才判不一致。
def ival(s_):
    """把一个印出的数展成它的四舍五入区间。"""
    t=str(s_).replace('−','-').replace('+','')
    dp=len(t.split('.')[1]) if '.' in t else 0
    h=0.5*10**(-dp); v=float(t)
    return (v-h,v+h)
def apply_iv(fn,*ivs):
    """在各输入区间的角点上求值,取 min/max —— 对本轮这些单调形式是精确的。"""
    import itertools
    vals=[]
    for corner in itertools.product(*[(a,b) for a,b in ivs]):
        try: vals.append(fn(*corner))
        except ZeroDivisionError: pass
    return (min(vals),max(vals)) if vals else (float('nan'),float('nan'))
def disjoint(a,b):
    return a[1]<b[0] or b[1]<a[0]

def check_sb(txt):
    out=[]
    for m in re.finditer(r'reliability is only \*\*('+NUM+r')\*\*[^)]*?= \*\*('+NUM+r')\*\*',txt):
        out.append(('SB',m.group(0)[:110],ival(m.group(1)),
                    apply_iv(lambda r:2*r/(1+r),ival(m.group(2)))))
    for m in re.finditer(r'信度只有 \*\*('+NUM+r')\*\*[^)]*?= \*\*('+NUM+r')\*\*',txt):
        out.append(('SB',m.group(0)[:110],ival(m.group(1)),
                    apply_iv(lambda r:2*r/(1+r),ival(m.group(2)))))
    return out
def check_ci(txt):
    out=[]
    pat=(r'\*\*('+NUM+r')\*\*[^.\n]{0,60}?se ('+NUM+r')[^\[\n]{0,60}?'
         r'(?:95% CI|95% 区间)\s*\[\s*('+NUM+r')\s*,\s*('+NUM+r')\s*\]')
    for m in re.finditer(pat,txt):
        E,S=ival(m.group(1)),ival(m.group(2))
        out.append(('CI-lo',m.group(0)[:110],ival(m.group(3)),
                    apply_iv(lambda e,s_:e-1.96*s_,E,S)))
        out.append(('CI-hi',m.group(0)[:110],ival(m.group(4)),
                    apply_iv(lambda e,s_:e+1.96*s_,E,S)))
    return out
def check_sd(txt):
    out=[]
    pat=(r'\*\*('+NUM+r')\*\*[^\n]{0,80}?('+NUM+r')\s*±\s*('+NUM+r')[^\n]{0,60}?'
         r'\*\*('+NUM+r')\s*sd\*\*')
    for m in re.finditer(pat,txt):
        if f(m.group(3))<=0: continue
        out.append(('sd',m.group(0)[:110],ival(m.group(4)),
                    apply_iv(lambda v,mu,sd:(v-mu)/sd,ival(m.group(1)),ival(m.group(2)),ival(m.group(3)))))
    return out
def check_drop(txt):
    out=[]
    pat=(r'\*\*('+NUM+r')\*\* -> (?:后 )?\*\*('+NUM+r')\*\*[^\n]{0,40}?'
         r'(?:掉幅|drop[^\n]{0,12}?)\s*\*?\*?('+NUM+r')%')
    for m in re.finditer(pat,txt):
        if abs(f(m.group(1)))<1e-12: continue
        out.append(('drop',m.group(0)[:110],ival(m.group(3)),
                    apply_iv(lambda b0,b1:100*(1-abs(b1)/abs(b0)),ival(m.group(1)),ival(m.group(2)))))
    return out
ALL=lambda t: check_sb(t)+check_ci(t)+check_sd(t)+check_drop(t)


def controls():
    """返回 (正对照通过, 负对照通过, 正对照条数, 负对照条数)。用前必须跑。"""
    BAD="its Spearman–Brown reliability is only **0.380** (`animated`↔`written` = **+0.468**), so"
    GOOD="加 `EARLY` 前 **-0.02873** -> 后 **-0.01481** · 掉幅 **48.44%**"
    rb=ALL(BAD); rg=ALL(GOOD)
    return (len(rb)>0 and all(disjoint(g,p) for _,_,g,p in rb),
            len(rg)>0 and all(not disjoint(g,p) for _,_,g,p in rg), len(rb), len(rg))

def coverage(txt):
    """(核到的三元组数, 页面上加粗数字的总数) —— 分母是覆盖率的诚实下限。"""
    bold=len(re.findall(r"\*\*[-−+]?\d+(?:\.\d+)?%?\*\*",txt))
    return len(ALL(txt)), bold
