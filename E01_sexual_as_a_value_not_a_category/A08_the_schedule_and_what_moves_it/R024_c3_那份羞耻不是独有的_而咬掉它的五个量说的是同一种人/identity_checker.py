import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A140 R420 -- 页面算术核对器:凡由恒等式相连、却分别印出的数,必须互相对上

⚠ **这是 CLOSURE,不是 FRONTIER**(§0.2)。它不分离任何世界,它**保护已有结论**。
本轮不冒充发现;它保护的是:页面上每一个**导出量**与它印在同一句话里的**输入量**。

`#375c` 靠的不是新数据,是让页面上两个数**互相约束**:
「可靠性 0.380(r = +0.468)」—— SB(r)=2r/(1+r) 是恒等式,两个数不能都对。
**一个孤立的数没有对错可言;两个由恒等式相连的数,只要都印出来,就会互相检查。**

ESTIMAND        从两页里抽出四类**恒等式三元组**,逐条核:
                ① Spearman–Brown:`reliability a` ↔ `r`      a = 2r/(1+r)
                ② 95% CI:`e` · `se s` · `[lo, hi]`           lo = e−1.96s · hi = e+1.96s
                ③ sd 距离:`v` · `零 m ± sd` · `k sd`          k = (v−m)/sd
                ④ 掉幅:`b0 -> b1` · `p%`                      p = 100(1−|b1|/|b0|)
KILL(条件式)  仅当正/负对照都过 -> 判:**是否存在任何一条不一致**。
POSITIVE CTRL  **前提,不是装饰**(`P5★` + `#374b`):在**合成的已知错**片段上必须开火 ——
                用 `#375c` 修好前的原句(「0.380 … +0.468」)。
NEGATIVE CTRL  在**合成的已知对**片段上必须沉默 —— 用 `#372a` 的掉幅(数都在账上,恒等式成立)。
⚠ 抽取用正则   而**正则是子串匹配的近亲**(`#374b`)-> **每一条抽到的都打印原文**,让错配当场可见。
⚠ 容差         页面上的数是**四舍五入印出的**,所以容差必须按印出的位数算,不能用一个固定 epsilon ——
                否则这个核对器会把**排版**报成**错误**,那是它自己的假阳性。
IMPOSSIBLE      抽不到的形式什么也说不了;漏掉的方向是**保守的**(少报不一致)。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

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

# ---- 对照:前提,不是装饰 ----
BAD="its Spearman–Brown reliability is only **0.380** (`animated`↔`written` = **+0.468**), so"
GOOD="加 `EARLY` 前 **-0.02873** -> 后 **-0.01481** · 掉幅 **48.44%**"
rb=ALL(BAD); rg_=ALL(GOOD)
print("对照(`P5★`:一个从未开火过的核对器,它的每一个「一致」都是沉默):")
for kind,src,gi,pi in rb:
    print(f"   正对照 [{kind}] 印区间 [{gi[0]:.5g}, {gi[1]:.5g}] vs 恒等式区间 [{pi[0]:.5g}, {pi[1]:.5g}] -> "
          f"{'**开火**' if disjoint(gi,pi) else '沉默'}")
for kind,src,gi,pi in rg_:
    print(f"   负对照 [{kind}] 印区间 [{gi[0]:.5g}, {gi[1]:.5g}] vs 恒等式区间 [{pi[0]:.5g}, {pi[1]:.5g}] -> "
          f"{'开火' if disjoint(gi,pi) else '**沉默**'}")
CP=len(rb)>0 and all(disjoint(g,p) for _,_,g,p in rb)
CN=len(rg_)>0 and all(not disjoint(g,p) for _,_,g,p in rg_)
print(f"   -> 正对照 **{CP}** · 负对照 **{CN}**\n")

# ---- 跑真页面 ----
rows=[]
for page in ('README.md','README_zh.md'):
    txt=pathlib.Path(page).read_text()
    for kind,src,gi,pi in ALL(txt):
        ok_=not disjoint(gi,pi)
        rows.append(dict(v_page=page,v_kind=kind,v_plo=gi[0],v_phi=gi[1],
                         v_ilo=round(pi[0],6),v_ihi=round(pi[1],6),
                         v_ok=bool(ok_),v_src=src.replace('\n',' ')[:100]))
T=pd.DataFrame(rows); check_columns(T,'R420')
T.to_csv(pathlib.Path(__file__).parent/'results'/'identity_checks.csv',index=False)
print(f"抽到 **{len(T)}** 条恒等式三元组(⚠ 每一条都打印原文,让错配当场可见):")
for r in T.itertuples():
    mark='✅' if r.v_ok else '❌ **不一致**'
    print(f"   {mark} [{r.v_kind}] 印区间 [{r.v_plo:.5g}, {r.v_phi:.5g}] vs "
          f"恒等式区间 [{r.v_ilo:.5g}, {r.v_ihi:.5g}]")
    print(f"        「{r.v_src}」")
BADN=int((~T.v_ok).sum())
print(f"\n   **一致 {int(T.v_ok.sum())} · 不一致 {BADN}**")

g=Gate('页面上由恒等式相连的数,是否互相对得上')
g.asserted('★ 正对照:在 `#375c` 修好前的原句上必须开火',CP,f"抽到 {len(rb)} 条,全部开火",kind='control')
g.asserted('★ 负对照:在已知一致的句子上必须沉默',CN,f"抽到 {len(rg_)} 条,全部沉默",kind='control')
g.asserted('★ 抽取非空(抽不到就什么也没核)',len(T)>0,f"{len(T)} 条",kind='control')
if CP and CN and len(T)>0:
    g.asserted('★ 注册的 kill:页面上不存在任何一条恒等式不一致',BADN==0,
               f"不一致 {BADN} 条")
else:
    g.asserted('★ 注册的 kill(对照未过或抽取为空 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
