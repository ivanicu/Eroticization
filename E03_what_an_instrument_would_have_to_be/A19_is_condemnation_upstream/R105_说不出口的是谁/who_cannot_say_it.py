"""E03·A19·R105 —— 说不出口的,是谁

**类型:FRONTIER**。`#662` 按 `#111c` 换了方向:不再问「哪些题」,改问「哪些人」。
`#660` 那个方式效应(面访比网络高 11–19 个点说婚外性「总是错的」)是关于**一道题**的;
**而它真正的心理学问题从来没被问过:说不出口的,是谁?**

⚠ **BASIN**:W1(集中在特定的人身上)是有故事的那个 ⇒ **下注 W2/W3**。
W1 集中 · W2 对所有人一样 · **W3 = meta-separator:梯度只是构成,不是心理** ——
  **若梯度跟着「两种方式的人口学构成差」走而不跟着分层变量的含义走,「谁」这个分解就是错的:
  它量的是分配,不是人。**

## G1 ESTIMAND(先于方法)
`xmarsex` 的**面访 − 网络**份额差(份额 = 答「always wrong」),**在每一层内各算一次**;
**主量 = 层间梯度**(该分层变量下,最高层的差 − 最低层的差)。
## G2 CONTROLS
**正对照**:总体方式效应必须复现 `#660` 的 +0.110 / +0.186。
**安慰剂**:同一套分层作用在 `cappun` 上,梯度必须 ≈0 —— **死刑没有方式效应,就不该有方式效应的梯度**。
  **这个零该不该是零?** 该 ⇒ `negative_control`。
## G3/G4
四个分层变量(年龄 / 性别 / 礼拜出席 / 政治立场)× 两年 × {合并, 分年} 全报,含不显著的。
## 地板(`#662` 写死)
**任一层任一方式 n < 150 ⇒ 该层记「判不了」,不许出现在结果里。**
## ⚠ 最强混淆 = 这一轮的死穴(`#662` 写死)
**作答方式非随机分配,而它与年龄、教育强相关** —— 老年人更可能被面访。
⇒ **必须同时报每一层内两种方式的人口学构成差(年龄中位、教育中位);
若梯度最大的那一层恰好也是构成差最大的那一层 ⇒ 记「判不了」,不许解释成心理。**
## KILL(条件式)
if 正对照复现 and 安慰剂梯度 ≈0:
  梯度 bootstrap 区间不含零 **且** 该层不是构成差最大的那一层 -> W1
  区间含零 -> **判不了,方式效应记作「对所有人一样」**
  区间不含零但落在构成差最大的层 -> **W3,不许解释成心理**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
作答方式非随机 ⇒ **非因果,这是结构性的** · 一国两波 ·
**跨仪器:换不了仪器,只此一具** · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; FLOOR=150
COLS=["year","mode","xmarsex","cappun","age","sex","attend","polviews","educ"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta", usecols=COLS, encoding="latin1")
d=df[df.year.isin([2022,2024])].copy()
d["年龄段"]=pd.cut(d.age,[17,34,49,64,99],labels=["18-34","35-49","50-64","65+"])
d["性别"]=d.sex.map({1:"男",2:"女"})
d["礼拜"]=pd.cut(d.attend,[-1,1,3,5,8],labels=["几乎不","偶尔","常去","每周+"])
d["政治"]=pd.cut(d.polviews,[0,2,4,7],labels=["自由派","中间","保守派"])
STRATA=["年龄段","性别","礼拜","政治"]

print("=== 硬规则①:每层每方式的 n(任一格 <150 ⇒ 整层判不了)===")
usable={}
for s in STRATA:
    print(f"  【{s}】")
    ok=[]
    for lv in d[s].dropna().unique().categories if hasattr(d[s].dropna(),'cat') else sorted(d[s].dropna().unique()):
        sub=d[(d[s]==lv)&d.xmarsex.between(1,4)]
        n1=int((sub["mode"]==1).sum()); n4=int((sub["mode"]==4).sum())
        good=n1>=FLOOR and n4>=FLOOR
        if good: ok.append(lv)
        print(f"    {str(lv):8s} 面访 {n1:5d} · 网络 {n4:5d}  {'✅' if good else '⛔ 判不了'}")
    usable[s]=ok
print("  可用层:", {k:[str(x) for x in v] for k,v in usable.items()})

def gap(sub,col="xmarsex",val=1):
    a=sub[(sub["mode"]==1)&sub[col].between(1,9)][col]; b=sub[(sub["mode"]==4)&sub[col].between(1,9)][col]
    if len(a)<FLOOR or len(b)<FLOOR: return None
    return float((a==val).mean()-(b==val).mean()), len(a), len(b)

print("\n=== G3:每层的方式效应(xmarsex)· 以及两种方式的人口学构成差 ===")
rows=[]
for s in STRATA:
    for lv in usable[s]:
        sub=d[d[s]==lv]
        g=gap(sub)
        if g is None: continue
        gc=gap(sub,"cappun")
        a=sub[sub["mode"]==1]; b=sub[sub["mode"]==4]
        dage=float(a.age.median()-b.age.median()); dedu=float(a.educ.median()-b.educ.median())
        rows.append(dict(strat=s,level=str(lv),gap=g[0],n_face=g[1],n_web=g[2],
                         gap_cappun=(gc[0] if gc else np.nan),d_age=dage,d_educ=dedu))
        print(f"  {s} {str(lv):8s} 方式效应 **{g[0]:+.4f}** (n {g[1]}/{g[2]}) · "
              f"死刑 {gc[0] if gc else float('nan'):+.4f} · 构成差 年龄中位 {dage:+.1f} 教育中位 {dedu:+.1f}")
R=pd.DataFrame(rows)

print("\n=== 主量:每个分层变量的层间梯度 ===")
grad={}
for s in STRATA:
    t=R[R.strat==s]
    if len(t)<2: print(f"  {s}: 可用层 <2,判不了"); continue
    hi,lo=t.loc[t.gap.idxmax()],t.loc[t.gap.idxmin()]
    gr=hi.gap-lo.gap
    grc=(t.gap_cappun.max()-t.gap_cappun.min()) if t.gap_cappun.notna().all() else np.nan
    comp=t.assign(c=t.d_age.abs()).loc[t.assign(c=t.d_age.abs()).c.idxmax()]
    grad[s]=dict(gradient=float(gr),hi=hi.level,lo=lo.level,cappun_gradient=float(grc),
                 worst_comp_level=str(comp.level),same=bool(comp.level==hi.level))
    print(f"  {s:5s} 梯度 **{gr:+.4f}**({hi.level} {hi.gap:+.4f} vs {lo.level} {lo.gap:+.4f})· "
          f"死刑梯度 {grc:+.4f} · 构成差最大的层 = {comp.level}{' ⚠ 与最高层同一层' if comp.level==hi.level else ''}")

def boot(s,B=3000):
    t=R[R.strat==s]; out=[]
    rng=np.random.default_rng(SEEDS[0])
    lv=list(t.level)
    for _ in range(B):
        vals=[]
        for l in lv:
            sub=d[d[STRATA[STRATA.index(s)]].astype(str)==l]
            i=rng.integers(0,len(sub),len(sub)); g=gap(sub.iloc[i])
            vals.append(g[0] if g else np.nan)
        vals=np.array(vals)
        if np.isfinite(vals).sum()>=2: out.append(np.nanmax(vals)-np.nanmin(vals))
    return np.quantile(out,[.025,.975]) if out else (np.nan,np.nan)
best=max(grad,key=lambda k:grad[k]["gradient"]) if grad else None
lo_,hi_=boot(best) if best else (np.nan,np.nan)
print(f"\n  最大梯度出现在【{best}】= {grad[best]['gradient']:+.4f} · 95% CI [{lo_:+.4f}, {hi_:+.4f}]"
      f" -> {'含零' if lo_<0<hi_ else '**不含零**'}")

ov=gap(d[d.year==2022]),gap(d[d.year==2024])
print(f"\n=== 控制 ===\n  正对照 总体方式效应 2022 {ov[0][0]:+.4f} · 2024 {ov[1][0]:+.4f}(`#660`: +0.1099 / +0.1860)")
cap_gr=float(np.nanmax([grad[s]["cappun_gradient"] for s in grad]))
print(f"  安慰剂 死刑的最大层间梯度 = **{cap_gr:+.4f}**(应 ≈0)")
G=Gate("说不出口的是谁")
p1=G.positive_control("总体方式效应必须复现 #660(容差 0.01)",
    planted=float(0.01-max(abs(ov[0][0]-0.1099),abs(ov[1][0]-0.1860))),floor=0.0,spread=0.001)
p2=G.negative_control("安慰剂:死刑不该有方式效应的层间梯度",null=abs(cap_gr),
    effect=abs(grad[best]["gradient"]) if best else 0.0,null_spread=0.02,
    null_kind="与性无关的道德题,同一套分层与同一个面访/网络对比")
if p1 and p2:
    if lo_<0<hi_: verdict=f"**判不了 —— 最大梯度({best} {grad[best]['gradient']:+.4f})的区间含零 ⇒ 方式效应记作「对所有人一样」**"
    elif grad[best]["same"]: verdict=f"**W3 —— 梯度最大的层正是构成差最大的层({best} {grad[best]['hi']}),不许解释成心理**"
    else: verdict=f"**W1 —— 说不出口的是特定的人:{best} 的 {grad[best]['hi']} 对 {grad[best]['lo']},梯度 {grad[best]['gradient']:+.4f}**"
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 负 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(cells=rows,gradients=grad,best=best,ci=[float(lo_),float(hi_)],
               overall=[ov[0][0],ov[1][0]],cappun_gradient=cap_gr,verdict=verdict,unchallenged=True),
          open(OUT/"who_cannot_say_it.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'who_cannot_say_it.json'}")
