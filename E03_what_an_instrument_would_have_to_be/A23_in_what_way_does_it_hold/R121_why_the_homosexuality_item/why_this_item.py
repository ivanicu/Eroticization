"""E03·A23·R121 —— 为什么是同性恋这一题

**类型:FRONTIER。** `#678` 的 NEXT。一个三十六年恒定的教育效应,应该在别的题上留下同样的指纹。
⚠ **而 `#672`/`#673` 已经证明「我用眼睛看出的共同点」不是证据** ⇒ 共同点**先**预注册成可判定的清单:
**用 `#672`(`R114`)已经枚举好的 GSS 自有题组归属**,它先于本轮存在,不可能被我挑着拟合。
**而这份清单里就有决定性的一刀**:容忍三组各含 `spkhomo`/`colhomo`/`libhomo` ——
**同一题组、同一格式、同一年份,只有对象不同。**

## 硬规则①(已跑,在挑任何题之前)
性四题+educ 基底 **n = 15,000**。
⚠ **`政府该管`(4 题)与 `机构信任`(6 题)与基底的联合 n = 0** —— GSS 分票设计,从没问过同一批人
⇒ **这十题结构性地进不了网格**,不是被我剔掉的。(`#672` 用的是各题组**自己**的梯度,不受此影响。)
入网格 **31 题**,联合 n **6,864–14,634**:容忍三组 15 题 · 支出 5 · 自杀 4 · 堕胎 7。

## G1 ESTIMAND
每题 X:`C_X(e)` = X 与性四题各自的**天花板归一相关的中位**,核加权于 educ;**Δ_X = 两端各三格点之差**。
**预注册主量(一个检验):** `Δ_contrast = median(Δ | 三个 *homo 题) − median(Δ | 容忍组其余 12 题)`。
**同一题组、同一格式、同一年份,只有对象不同 ⇒ 这是这份数据能给的最干净的对照。**
## G2 CONTROLS
**正对照**:`homosex` 自己与性四题其余三题的耦合必须复现 `#677` 量级(该题在块内,>0.20)。
**零**:打乱 educ(`negative_control`;**这个零该不该是零?** 该 —— 若对象不重要,两组的 Δ 应无差别)。
## G3 多重性:31 题的 Δ **全部登**,BH(q=0.05, C=31)覆盖整个网格,**不许挑显著的报**。
## KILL(条件式)
if 正对照复现:
  `Δ_contrast` 超零 -> **机制是「关于同性恋的态度被绑在一起」** ·
  含零 -> **不是关于对象的,记「测不出机制」——这也是一个结果**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
容忍题是**二值**,天花板受边际卡死(`#647`:GSS 二值最高约 0.36)⇒ 已除天花板但残余偏斜无法完全消掉;
**跨仪器:MFQ 无容忍题组** ⇒ 结构性拿不到第二具仪器的对象对照;因果:横断面无干预。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]
B={"容忍·言论":["spkath","spkrac","spkcom","spkmil","spkhomo"],
   "容忍·任教":["colath","colrac","colcom","colmil","colhomo"],
   "容忍·藏书":["libath","librac","libcom","libmil","libhomo"],
   "支出":["natheal","nateduc","natcrime","natenvir","natrace"],
   "自杀":["suicide1","suicide2","suicide3","suicide4"],
   "堕胎":["abdefect","abnomore","abhlth","abpoor","abrape","absingle","abany"]}
ALL=[c for v in B.values() for c in v]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ"]+SEX+ALL, apply_value_formats=False, encoding="latin1")
grid=np.arange(8,19.01,1.0); BW=2.5
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def norm_at(x,y,W):
    r=wc(x,y,W)
    if not np.isfinite(r) or abs(r)<1e-12: return np.nan
    idx=W>np.quantile(W,0.5)
    xs=np.sort(x[idx]); ys=np.sort(y[idx]); ys=ys if r>0 else ys[::-1]
    c=np.corrcoef(xs,ys)[0,1]
    return r/abs(c) if np.isfinite(c) and abs(c)>1e-9 else np.nan
def item_delta(X,S,e):
    out=[]
    for g in grid:
        W=np.exp(-0.5*((e-g)/BW)**2)
        if W.sum()<150: out.append(np.nan); continue
        v=[norm_at(X,s,W) for s in S]; v=[u for u in v if np.isfinite(u)]
        out.append(float(np.median(v)) if v else np.nan)
    y=np.array(out); m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
CACHE={}
def prep(c):
    if c in CACHE: return CACHE[c]
    j=df.dropna(subset=SEX+["educ",c])
    S=[j[s].rank().to_numpy(float) for s in SEX if s!=c]
    CACHE[c]=(j[c].rank().to_numpy(float),S,j["educ"].to_numpy(float),len(j))
    return CACHE[c]
print(f"{'题':12s} {'题组':10s} {'n':>7s} {'Δ':>9s} {'p':>8s}  对象同性恋")
rows=[]
for bn,items in B.items():
    for c in items:
        X,S,e,n=prep(c); d=item_delta(X,S,e)
        rng=np.random.default_rng(20260806)
        nul=np.array([abs(item_delta(X,S,rng.permutation(e))) for _ in range(120)])
        p=float(np.nanmean(nul>=abs(d)))
        rows.append(dict(item=c,batt=bn,n=n,delta=d,p=p,homo=c.endswith("homo")))
        print(f"{c:12s} {bn:10s} {n:>7,} {d:>+9.4f} {p:>8.4f}  {'✅' if c.endswith('homo') else '—'}")
ps=sorted(r["p"] for r in rows); C=len(rows)
surv=[r["item"] for r in rows if r["p"]<=0.05*(ps.index(r["p"])+1)/C]
print(f"\n多重性 BH(q=0.05, C={C}):**存活 {len(surv)} 题** -> {surv}")
tol=[r for r in rows if r["batt"].startswith("容忍")]
dh=float(np.median([r["delta"] for r in tol if r["homo"]]))
do=float(np.median([r["delta"] for r in tol if not r["homo"]]))
print(f"\n预注册主量:容忍组内 *homo 三题 Δ中位 **{dh:+.4f}** − 其余 12 题 **{do:+.4f}** = **{dh-do:+.4f}**")
rng=np.random.default_rng(991); cn=[]
for _ in range(200):
    e2={c:rng.permutation(prep(c)[2]) for c in [r["item"] for r in tol]}
    a=[item_delta(*prep(r["item"])[:2],e2[r["item"]]) for r in tol if r["homo"]]
    b=[item_delta(*prep(r["item"])[:2],e2[r["item"]]) for r in tol if not r["homo"]]
    cn.append(abs(float(np.median(a))-float(np.median(b))))
cn=np.array(cn); q=float(np.nanquantile(cn,.95)); pc=float(np.nanmean(cn>=abs(dh-do)))
print(f"  打乱 educ 的零:95% 分位 **{q:.4f}** · p = **{pc:.4f}**  {'✅ 超零' if pc<0.05 else '⛔ 落在地板上'}")
Xh,Sh,eh,_=prep("homosex") if "homosex" in CACHE else (None,None,None,None)
jj=df.dropna(subset=SEX+["educ"])
Xo=jj["homosex"].rank().to_numpy(float); So=[jj[s].rank().to_numpy(float) for s in SEX if s!="homosex"]
W1=np.ones(len(jj)); pos=float(np.median([norm_at(Xo,s,W1) for s in So]))
G=Gate("为什么是同性恋这一题")
p1=G.positive_control("homosex 与性四题其余三题的耦合复现 #677 量级(>0.20)",planted=pos,floor=0.20,spread=0.01)
p2=G.negative_control("打乱 educ 后,*homo 与非 *homo 的 Δ 差应消失",
                      null=q,effect=abs(dh-do),null_spread=0.02,
                      null_kind="同一题组内打乱 educ —— 若对象不重要,两组 Δ 应无差别")
verdict=(("**机制是对象:容忍组内只有对象是同性恋的三题,与性四题的耦合随教育上升,"
          f"比同组其余 12 题高 {dh-do:+.4f}(p={pc:.4f})**") if (p1 and p2 and pc<0.05)
         else (f"**不是关于对象的:*homo 与其余的 Δ 差 {dh-do:+.4f},p = {pc:.4f} 落在零里 ⇒ 测不出机制,这也是一个结果**"
               if p1 else "UNVERIFIED —— 正对照失败"))
print(f"\n{verdict}"); print(G)
json.dump(dict(rows=rows,bh_survivors=surv,d_homo=dh,d_other=do,contrast=dh-do,q95=q,p=pc,
               verdict=verdict,unchallenged=True),open(OUT/"why_this_item.json","w"),indent=1,ensure_ascii=False)
