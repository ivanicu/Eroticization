"""E03·A19·R103 —— 那 11–19 个百分点,是性题特有的,还是所有道德题都有

**类型:FRONTIER**。`#660` 量到:同年同问卷,面访比网络高 11–19 个点说婚外性「总是错的」,
而死刑不动。**但那只有一道题。** 若机制是「在人前说不出口」,它应当**在性题上大、在非性道德题上小**。

⚠ **BASIN**:`#660` 是我这一段第一句关于人的话,**我喜欢它** ⇒ **本轮下注 W2/W3**,即**下注它该被降级**。
W1 性题特有 ⇒ 机制成立。 W2 泛道德题效应 ⇒ `#660` 的标题该改成「道德题」而不是「性」。
W3 天花板伪影 ⇒ 性题边际更极端,份额差本来就有更大余地(`#647`)。
**Meta-separator**:若差随**边际**走而不随**内容**走,**「性 vs 非性」这个分解本身就是错的**。

G1 ESTIMAND:每题每年 **面访份额 − 网络份额**(份额 = 最保守/最谴责那一档),
  **主量 = 性题的差中位 − 非性题的差中位**。**并同时报天花板归一版**(份额差 ÷ 边际允许的最大差)。
G2:**正对照** `xmarsex` 必须复现 +0.110 / +0.186。**安慰剂** `cappun` 必须 ≈0。
  **这个零该不该是零?** 该 —— 死刑与性无关,方式效应本就该趋零 ⇒ `negative_control`。
G3:19 题 × 2 年全报,含**堕胎**这一组单列(它算不算性题有争议,**不由我判**)。
G4:{份额差, 天花板归一} × {2022, 2024, 合并}。
KILL:if 正对照复现 and 安慰剂≈0:
  性−非性 的 bootstrap 区间**不含零 且 两个版本同号** -> W1
  含零 -> **判不了**,`#660` 降级为「一题的观察」
  两版不同号 -> **W3 天花板伪影**
else UNVERIFIED
IMPOSSIBLE(不写 planned):**作答方式非随机分配** ⇒ 非因果 · 一个国家两波 ·
  **跨仪器:换不了仪器,只此一具** —— 没有第二份调查在同一年内同时给出两种作答方式且带道德题组 · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
SEX   = ["premarsx","xmarsex","homosex","teensex","pornlaw","sexeduc"]
NONSEX= ["cappun","gunlaw","letdie1","suicide1","suicide4","spanking","polhitok","obey","fefam","natcrime","helppoor"]
ABORT = ["abany","abnomore"]
ALL=SEX+NONSEX+ABORT
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
                         usecols=["year","mode"]+ALL, encoding="latin1")
d=df[df.year.isin([2022,2024])]

def cell(c,y):
    s=d[(d.year==y)&d[c].between(1,9)]
    a=s[s["mode"]==1][c]; b=s[s["mode"]==4][c]
    if len(a)<200 or len(b)<200: return None
    v=1                                   # 最保守/最谴责那一档 = 码 1(GSS 全部如此编)
    pa,pb=(a==v).mean(),(b==v).mean()
    diff=pa-pb
    # 天花板:给定两边总人数与合并份额,份额差的最大可能值 = min(1,·) —— 用边际允许的最大差
    p=( (a==v).sum()+(b==v).sum() )/(len(a)+len(b))
    cap=min(p*(len(a)+len(b))/len(a), 1.0) - max(0.0, (p*(len(a)+len(b))-len(a))/len(b))
    se=np.sqrt(pa*(1-pa)/len(a)+pb*(1-pb)/len(b))
    return dict(item=c,year=int(y),n_face=int(len(a)),n_web=int(len(b)),face=float(pa),web=float(pb),
                diff=float(diff),ci=[float(diff-1.96*se),float(diff+1.96*se)],
                cap=float(cap),norm=float(diff/cap) if cap>1e-9 else np.nan)
rows=[r for c in ALL for y in (2022,2024) if (r:=cell(c,y))]
R=pd.DataFrame(rows)
grp=lambda c: "性" if c in SEX else ("堕胎" if c in ABORT else "非性")
R["组"]=R.item.map(grp)
print("=== G3:19 题 × 2 年全报(份额 = 最谴责那一档)===")
print(f"  {'题':10s}{'组':>4s}{'年':>6s}{'面访':>8s}{'网络':>8s}{'差':>9s}{'天花板':>8s}{'归一':>8s}")
for r in rows:
    print(f"  {r['item']:10s}{grp(r['item']):>4s}{r['year']:>6d}{r['face']:>8.4f}{r['web']:>8.4f}"
          f"{r['diff']:>+9.4f}{r['cap']:>8.3f}{r['norm']:>+8.4f}")

def med(g,col="diff"): return float(R[R.组==g][col].median())
main = med("性")-med("非性"); main_n = med("性","norm")-med("非性","norm")
print(f"\n=== 主量 ===\n  性 中位差 {med('性'):+.4f} · 非性 {med('非性'):+.4f} · **性−非性 = {main:+.4f}**")
print(f"  归一版:性 {med('性','norm'):+.4f} · 非性 {med('非性','norm'):+.4f} · **差 = {main_n:+.4f}**")
print(f"  堕胎(单列,不由我判归属):中位差 {med('堕胎'):+.4f} · 归一 {med('堕胎','norm'):+.4f}")

def boot(col,n=4000):
    rng=np.random.default_rng(SEEDS[0]); s=R[R.组=="性"][col].values; t=R[R.组=="非性"][col].values; out=[]
    for _ in range(n):
        out.append(np.median(rng.choice(s,len(s)))-np.median(rng.choice(t,len(t))))
    return np.quantile(out,[.025,.975])
lo,hi=boot("diff"); lo2,hi2=boot("norm")
print(f"  95% CI(份额差)= [{lo:+.4f}, {hi:+.4f}] -> {'含零' if lo<0<hi else '**不含零**'}")
print(f"  95% CI(归一) = [{lo2:+.4f}, {hi2:+.4f}] -> {'含零' if lo2<0<hi2 else '**不含零**'}")

xm={r['year']:r['diff'] for r in rows if r['item']=='xmarsex'}
cp=[abs(r['diff']) for r in rows if r['item']=='cappun']
print(f"\n=== 控制 ===\n  正对照 xmarsex {xm} · `#660` 记的是 +0.1099 / +0.1860")
print(f"  安慰剂 cappun |差| = {[round(x,4) for x in cp]}")
G=Gate("那 11–19 个百分点,是性题特有还是所有道德题都有")
p1=G.positive_control("xmarsex 必须复现 #660(容差 0.01)",
                      planted=float(0.01-max(abs(xm[2022]-0.1099),abs(xm[2024]-0.1860))),floor=0.0,spread=0.001)
p2=G.negative_control("安慰剂:死刑的方式效应该趋零",null=float(np.median(cp)),effect=abs(main),
                      null_spread=0.02,null_kind="与性无关的道德题,同一面访/网络对比")
same = (main>0)==(main_n>0)
if p1 and p2:
    if (lo<0<hi) or (lo2<0<hi2): verdict=f"**判不了 —— 区间含零。`#660` 降级为「一题的观察」**"
    elif not same: verdict=f"**W3 —— 两版不同号,天花板伪影**"
    elif main>0: verdict=f"**W1 —— 性题特有:性−非性 = {main:+.4f}(归一 {main_n:+.4f}),两版同号且都不含零**"
    else: verdict=f"**方向相反:非性题的方式效应更大**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(rows=rows,median_sex=med("性"),median_nonsex=med("非性"),main=main,
               median_sex_norm=med("性","norm"),median_nonsex_norm=med("非性","norm"),main_norm=main_n,
               ci=[float(lo),float(hi)],ci_norm=[float(lo2),float(hi2)],abortion=med("堕胎"),
               verdict=verdict,unchallenged=True),
          open(OUT/"sexual_or_all_moral.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'sexual_or_all_moral.json'}")

# ── ⚠ 轮内自查:码 1 在不同题上不是同一个意思 ─────────────────────────────
# `cappun` 1=赞成死刑(严厉)· `letdie1`/`suicide1`/`abany` 1=允许(宽容)· `premarsx` 1=总是错的(严厉)。
# **有符号的中位会让非性组内部符号抵消,从而把「性−非性」抬高** —— 这是一个偏向我的偏差。
# ⇒ **用 |差| 重算。撑不住,我这个结果就是符号抵消造出来的。**
print("\n=== 轮内自查:改用 |差|(符号不再抵消)===")
R["absdiff"]=R["diff"].abs(); R["absnorm"]=R["norm"].abs()
ma=med("性","absdiff")-med("非性","absdiff"); man=med("性","absnorm")-med("非性","absnorm")
print(f"  |差| :性 {med('性','absdiff'):+.4f} · 非性 {med('非性','absdiff'):+.4f} · **差 = {ma:+.4f}**")
print(f"  |归一|:性 {med('性','absnorm'):+.4f} · 非性 {med('非性','absnorm'):+.4f} · **差 = {man:+.4f}**")
la,ha=boot("absdiff"); lan,han=boot("absnorm")
print(f"  95% CI(|差|)  = [{la:+.4f}, {ha:+.4f}] -> {'含零' if la<0<ha else '**不含零**'}")
print(f"  95% CI(|归一|)= [{lan:+.4f}, {han:+.4f}] -> {'含零' if lan<0<han else '**不含零**'}")
survive = (la>0 or ha<0) and (lan>0 or han<0) and (ma>0)==(man>0)
print(f"  ⇒ **{'四个版本全部同号且不含零 ⇒ 不是符号抵消造的' if survive else '撑不住 ⇒ 结果是符号抵消造的,降级'}**")
d=json.load(open(OUT/"sexual_or_all_moral.json"))
d["abs_check"]=dict(main_abs=ma,main_absnorm=man,ci_abs=[float(la),float(ha)],
                    ci_absnorm=[float(lan),float(han)],survives=bool(survive),
                    note="码 1 在不同题上含义不同,有符号中位会让非性组抵消并抬高对比;|差| 版是防这个的")
json.dump(d,open(OUT/"sexual_or_all_moral.json","w"),indent=1,ensure_ascii=False)
