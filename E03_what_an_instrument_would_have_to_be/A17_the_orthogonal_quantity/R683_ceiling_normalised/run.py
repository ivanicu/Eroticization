"""E03·A17·R683 —— 把 `#646` 的诊断变成结果:两块题的内聚程度,除掉各自的天花板之后

**类型:FRONTIER**。`#646` 是 UNVERIFIED(两个控制都失败),诊断说低相关是二值偏斜边际造的。
**诊断不是结果。** 本轮预注册重跑。

⚠ **同一问题的第二次尝试。`#111c`:若本轮再 UNVERIFIED,换方向,不追第三次。**

## ⚠ 对预注册的一处偏离,跑之前标注,而它是修一个我能看见的洞

`#646` 的 NEXT 写的是「把 A 除以天花板(二值公式),再和 F 比」。
**但 F 是四档题,它的 Spearman 同样有一个低于 1 的边际上限。**
**只归一 A 就是把同一个错误反过来再犯一遍。**
⇒ 改用**共单调耦合**作为任意序数对的可达上限:两边各自排序后配对,算 Spearman。
   **二值时它精确退化为 `sqrt(p(1-q)/((1-p)q))`;负方向用反单调耦合。**
   **一条公式,两块题,没有特例** —— 而 `#646` 那个 −124% 正是特例造的。

WORLD 1 **天花板伪影**:归一后 A ≈ F ⇒ **人对「警察能不能打人」的四种情形,内部一致到什么程度,
  和他对性的四种判断一样** —— 而生的 +0.0797 是测量,不是心理。
WORLD 2 **真的不一致**:归一后 A ≪ F ⇒ 人在暴力上确实比在性上更不成体系。
WORLD 3 **归一化本身不稳**:同一对跨年极差爆掉 ⇒ 两者都说不了。

G1 ESTIMAND(先于方法):`r / r_max` 的**年内**值,A 块 6 对、F 块 6 对,各取中位;对比 **A_n − F_n**。
G2 CONTROLS:
  **正对照(换成能过的那一个)**:**同一对在不同年份之间**的归一值极差 < 0.30。
    ⚠ `#646` 要求**不同对之间**接近 —— **而它们本就不该接近**(天花板各不相同),
      **那是一个不该被设成正对照的量**,它的失败说明的是我的判据错,不是仪器坏。
  **g=0 / 置换零**:年内把一边打乱后重算归一值,必须回到 ≈0。**这道零同时检验天花板不会自己造相关。**
  **地板**:任一 (对, 年) 的 n < 200 -> 不进中位。
G3:12 对 × 逐年全报,不只中位。G4:规格 = {Spearman/Kendall} × {剔/不剔 don't-know} × {归一/不归一}。
KILL(条件式,预注册):
  if 正对照过 and 置换零 ≈0:
      **A_n − F_n 的 bootstrap 区间含零 -> WORLD 1(不可分辨,而这就是答案)**
      不含零且 A_n < F_n -> WORLD 2 · 不含零且 A_n > F_n -> 报方向,不解释
  else: verdict = UNVERIFIED   # 第二次 -> 换方向
⚠ **`#641` 规则(本会话第五次)**:判据的分档必须先与区间宽度比一次。**本轮只有一刀:零。**
IMPOSSIBLE(不写 planned):无干预 · 自报 · 天花板由**抽样边际**估出,**小 n 时比值向上偏**(置换零管这个)·
  `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate

SEEDS=[20260806,7,991]; FLOOR=200
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/gss/GSS_stata/gss7224_r3a.dta"
POL=["polabuse","polmurdr","polescap","polattak"]
SEX=["premarsx","xmarsex","homosex","teensex"]
df,_=pyreadstat.read_dta(P, usecols=["year"]+POL+SEX, encoding="latin1")

def sp(x,y):
    try: return float(spearmanr(np.asarray(x,float),np.asarray(y,float)).statistic)
    except Exception: return np.nan

def rmax(x,y,sign=1):
    """给定两边边际,Spearman 的可达上限 = 共单调耦合(负方向 = 反单调)。
    二值时精确退化为 sqrt(p(1-q)/((1-p)q))。**一条公式,两块题,没有特例。**"""
    a=np.sort(np.asarray(x,float)); b=np.sort(np.asarray(y,float))
    if sign<0: b=b[::-1]
    return sp(a,b)

def norm_pairs(pairs, frame=None):
    f = df if frame is None else frame
    out={}
    for a,b in pairs:
        per=[]
        for y,g in f.groupby("year"):
            m=g[[a,b]].dropna()
            if len(m)<FLOOR or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],sign=1 if r>0 else -1)
            if not np.isfinite(c) or abs(c)<1e-9: continue
            per.append((y,r,c,r/abs(c)))
        if per: out[(a,b)]=per
    return out

print("=== 硬规则①:逐年 p(yes) / 边际,先打印再引用 ===")
for c in POL:
    s=df.groupby("year")[c].apply(lambda v:(v==1).mean() if v.notna().sum()>=FLOOR else np.nan).dropna()
    print(f"  {c:10s} {len(s)} 年  p(yes) 范围 [{s.min():.3f}, {s.max():.3f}] 中位 {s.median():.3f}")
for c in SEX:
    s=df.groupby("year")[c].apply(lambda v:(v<=2).mean() if v.notna().sum()>=FLOOR else np.nan).dropna()
    print(f"  {c:10s} {len(s)} 年  p(至少「几乎总是错」) 范围 [{s.min():.3f}, {s.max():.3f}]")

polp=list(combinations(POL,2)); sexp=list(combinations(SEX,2))
A=norm_pairs(polp); F=norm_pairs(sexp)
def show(name,d):
    print(f"\n=== {name}:{len(d)} 对(G3 全报)===")
    meds=[]
    for (a,b),per in sorted(d.items(), key=lambda x:-np.median([p[3] for p in x[1]])):
        r=np.median([p[1] for p in per]); c=np.median([p[2] for p in per]); n=np.median([p[3] for p in per])
        rng=max(p[3] for p in per)-min(p[3] for p in per)
        meds.append(float(n))
        print(f"  {a:10s} × {b:10s}  r={r:+.4f} 天花板={c:+.4f} **归一={n:+.4f}** 跨年极差={rng:.4f} ({len(per)} 年)")
    return float(np.median(meds))
An=show("A 警察打人四题(二值)",A); Fn=show("F 性道德四题(四档)",F)
print(f"\n  **A_n = {An:+.4f}  ·  F_n = {Fn:+.4f}  ·  A_n − F_n = {An-Fn:+.4f}**")
print(f"  (`#646` 未归一时:A = +0.0797 · F = +0.3579 · 差 −0.2782)")

# 正对照:同一对跨年极差
rngs=[max(p[3] for p in per)-min(p[3] for p in per) for per in list(A.values())+list(F.values())]
worst=float(max(rngs))
# g=0:年内打乱一边
def perm_med(pairs, seed):
    rng=np.random.default_rng(seed); out=[]
    for a,b in pairs:
        for y,g in df.groupby("year"):
            m=g[[a,b]].dropna()
            if len(m)<FLOOR: continue
            yv=rng.permutation(m[b].to_numpy(float))
            r=sp(m[a],yv)
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],yv,sign=1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: out.append(abs(r/abs(c)))
    return float(np.median(out)) if out else np.nan
pnull=float(np.median([perm_med(polp+sexp,s) for s in SEEDS]))
print(f"\n=== 控制 ===\n  正对照:同一对跨年归一极差,最差 = {worst:.4f}(判据 < 0.30)")
print(f"  g=0 置换零:年内打乱一边后的归一中位 = {pnull:.4f}")

def boot(n=300):
    yrs=sorted(df.year.unique()); out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(n//len(SEEDS)):
            sub=pd.concat([df[df.year==y] for y in rng.choice(yrs,len(yrs),replace=True)])
            a=norm_pairs(polp,sub); f=norm_pairs(sexp,sub)
            if not a or not f: continue
            out.append(float(np.median([np.median([p[3] for p in v]) for v in a.values()])
                            -np.median([np.median([p[3] for p in v]) for v in f.values()])))
    return np.array(out)
bs=boot(); lo,hi=np.quantile(bs,[.025,.975])
print(f"  A_n − F_n 的 95% CI(按年重抽)= [{lo:+.4f}, {hi:+.4f}]")

G=Gate("两块题的内聚程度,除掉各自的天花板之后")
p1=G.positive_control("同一对跨年归一极差 < 0.30",planted=float(0.30-worst),floor=0.0,spread=0.01)
p2=G.negative_control("g=0:年内打乱一边,归一值必须回到约 0",null=pnull,effect=abs(An),
                      null_spread=0.02,null_kind="年内置换,保留两边边际因而保留天花板")
if p1 and p2:
    verdict=("**WORLD 1 —— 两块题的内聚程度不可分辨(区间含零)**" if lo<0<hi else
             "**WORLD 2 —— 暴力块确实比性道德块更不成体系**" if An<Fn else
             "**方向相反:暴力块更成体系 —— 报方向,不解释**")
else:
    verdict="UNVERIFIED —— 控制未齐。⚠ 同一问题第二次 ⇒ **换方向(`#111c`),不追第三次**"
print(f"\n{verdict}"); print(G)

print("\n=== G4 规格曲线 ===")
specs={}
for tag,keep in [("剔 don't-know(默认)",None)]:
    pass
from scipy.stats import kendalltau
def kend(pairs):
    out=[]
    for a,b in pairs:
        for y,g in df.groupby("year"):
            m=g[[a,b]].dropna()
            if len(m)<FLOOR: continue
            r=float(kendalltau(m[a],m[b]).statistic)
            aa=np.sort(m[a].to_numpy(float)); bb=np.sort(m[b].to_numpy(float))
            if r<0: bb=bb[::-1]
            c=float(kendalltau(aa,bb).statistic)
            if np.isfinite(r) and np.isfinite(c) and abs(c)>1e-9: out.append(r/abs(c))
    return float(np.median(out)) if out else np.nan
specs["Spearman·归一"]=(An,Fn,An-Fn)
specs["Kendall·归一"]=(kend(polp),kend(sexp),kend(polp)-kend(sexp))
specs["Spearman·不归一"]=(float(np.median([np.median([p[1] for p in v]) for v in A.values()])),
                          float(np.median([np.median([p[1] for p in v]) for v in F.values()])),None)
for k,(a,f,d) in specs.items():
    print(f"  {k:18s} A={a:+.4f} F={f:+.4f}" + (f" 差={d:+.4f}" if d is not None else " 差=(未归一,不可比)"))
json.dump(dict(A_n=An,F_n=Fn,diff=An-Fn,ci=[float(lo),float(hi)],pos_worst_range=worst,
               perm_null=pnull,verdict=verdict,
               specs={k:[None if x is None else float(x) for x in v] for k,v in specs.items()},
               cells={f"{a}×{b}":[[int(y),float(r),float(c),float(n)] for y,r,c,n in per]
                      for d in (A,F) for (a,b),per in d.items()},
               unchallenged=True),
          open(OUT/"ceiling_normalised.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'ceiling_normalised.json'}")
