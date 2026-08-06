"""E03·A33·R174 —— 同一组题,换一个陪衬域,它的比值会挪多少

**类型:FRONTIER。这是 `#730`③ —— 而它的正结果会给这一页 28 行每一行加一个限定。**

**心理学的那一句(元层,但它决定这一页二十多句话怎么读):
我说「性道德的紧密度是随机题组的 4.56 倍」。那个「随机题组」是从**我选的池**里抽的。
换一个陪衬域,同一批人、同一组性题,这句话会变成几倍?**

## 缺口
`#730` 拆出「GSS 4.56× 对 NSFG 0.96×」的缺口 **84% 在分母**,而分母是**同池随机块** ——
**而池是我选的**:GSS 我一直用「性四 + 警察四」。**若换个陪衬域分母就变,那 4.56 这个数
就不是关于性道德的,而是关于「我恰好拿警察题来做背景」的。**

## G1 ESTIMAND
**分子固定**:性四题的最弱一环(天花板归一 · 最优符号)——**它不依赖陪衬域,所以整个展布都在分母上。**
**分母**:`池 = 性四 ∪ 陪衬域` 的同 k 全枚举零(排除任何整块落在同一真域内的块,与 k 无关 —— `#730` 的修法)。
**估计量 = 比值在五个陪衬域上的展布 max/min。**

## W1 / W2
| | max/min | 读法 |
|---|---|---|
| **W1 池无关** | **≤ 1.5** | 4.56 是关于性道德的,页上不必改 |
| **W2 池决定它** | **> 2** | **页上 28 行每一个比值都要带「相对哪个池」,而这是一个页级限定** |

⚠ **W2 的正结果我不高兴** —— 它给整页加限定。**而这正是本轮要设计成能出 W2 的理由。**

## G2 CONTROLS
**④ 正对照**:警察池必须复现 `#718`/`#730` 的 **零 0.0911 · 比值 4.56×**(容差 0.005)。
**零** = `negative_control`,**零的种类 = 同一批人、`性四 ∪ 陪衬域` 这个池、同样 k=4、
同样逐年取中位、同样最优符号,只打散「哪四题算一组」;全枚举,排除任何纯域块。**
**SHAM**:同一件事对**陪衬域自己**做一次 —— 陪衬域的比值也报出来,
**用来分开「性题特别」与「任何一个成套的域在任何池里都会高」。**
## G3:5 个陪衬域 × (性四 · 陪衬域自己) = **10 格全报**。G4:k=4 · 逐年中位 · 最优符号,并报贪心那一格。
## ⑤ 停止条件(跑之前写死)
- **警察池复现不到 0.005 ⇒ UNVERIFIED 并停。**
- **任一陪衬域与性四题的逐对联合 n < 2,000,或逐年 ≥100 的年数 < 5 ⇒ 那一格记「判不了」,不进展布。**
- **max/min > 2 ⇒ 判 W2**,页上加页级限定;**≤ 1.5 ⇒ 判 W1**;落在 1.5–2 ⇒ **记「判不了」,报区间。**
## IMPOSSIBLE(不写 planned)
**换不了仪器**:本轮问的就是 GSS 内部的池依赖性;NSFG 只有一个候选陪衬域(家庭七),
**做不出展布** ⇒ 本轮的结论**只对 GSS 成立**,而**页级限定是从它外推的,必须标明是外推。**
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from lib.blocks import pairmat, opt_batch, weakest_greedy
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]
COMP={"警察四":["polabuse","polmurdr","polescap","polattak"],
      "堕胎四":["abdefect","abnomore","abpoor","abrape"],
      "自杀四":["suicide1","suicide2","suicide3","suicide4"],
      "性别角色三":["fefam","fepresch","fepol"],
      "支出四":["natspac","natenvir","natheal","natcrime"]}
allc=sorted({c for v in COMP.values() for c in v})
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+SEX+allc,encoding="latin1")
print("=== 硬规则①:逐对联合 n(零真正用到的量),不是「八题同时非缺失」 ===")
ok={}
for nm,P in COMP.items():
    ns=[];yrs=[]
    for a in SEX:
        for b in P:
            m=g[["year",a,b]].dropna(); ns.append(len(m)); yrs.append(int((m.groupby("year").size()>=100).sum()))
    ok[nm]=(min(ns)>=2000 and min(yrs)>=5)
    print(f"  {nm:10s} 性×陪衬 {len(ns)} 对 · 联合 n [{min(ns):,}, {max(ns):,}] · 逐年≥100 的年数 [{min(yrs)}, {max(yrs)}]"
          f"  {'✅' if ok[nm] else '⚠ 判不了'}")
def ratio(pool,target):
    M=pairmat(g,pool,year="year",floor=100); I={c:i for i,c in enumerate(pool)}
    doms=[set(SEX)]+[set(v) for v in COMP.values()]
    def pure(c):
        st={pool[i] for i in c}
        return any(st<=d for d in doms)
    ix=[I[c] for c in target]
    obs=float(opt_batch(M,np.array([ix]))[0]); grd=weakest_greedy(M,ix)
    allb=[c for c in itertools.combinations(range(len(pool)),len(target)) if not pure(c)]
    v=opt_batch(M,np.array(allb)); v=v[np.isfinite(v)]
    q=float(np.quantile(v,0.95))
    return obs,grd,q,(obs/q if q>0 else np.nan),len(allb),int(v.size)
print("\n=== ④ 正对照:警察池必须复现 `#718`/`#730` 的 零 0.0911 · 4.56×(容差 0.005)===")
o,gr,q,r,nb,nv=ratio(SEX+COMP["警察四"],SEX)
d=max(abs(q-0.0911),abs(r-4.56)/100)
print(f"  最弱一环 {o:+.4f} · 零的 95% 分位 {q:.4f} · 比值 **{r:.2f}×** · 块 {nv}/{nb} · "
      f"差 {abs(q-0.0911):.4f} {'✅' if abs(q-0.0911)<=0.005 else '⛔ ⑤ 触发'}")
if abs(q-0.0911)>0.005:
    print("⛔ 停"); sys.exit(0)
print("\n=== G3 全格:5 个陪衬域 ×(性四 · 陪衬域自己)===")
print(f"{'陪衬域':12s}{'性四最弱':>10s}{'零95%':>9s}{'性四比值':>10s}{'陪衬自己':>10s}{'陪衬比值':>10s}{'块数':>7s}")
res={}
for nm,P in COMP.items():
    if not ok[nm]: print(f"{nm:12s}   ⚠ 判不了(地板)"); res[nm]=dict(undecidable=True); continue
    pool=SEX+P
    o,gr,q,r,nb,nv=ratio(pool,SEX)
    o2,gr2,q2,r2,_,_=ratio(pool,P)
    res[nm]=dict(obs=o,greedy=gr,null=q,ratio=r,comp_obs=o2,comp_ratio=r2,blocks=nv,undecidable=False)
    print(f"{nm:12s}{o:>+10.4f}{q:>9.4f}{r:>10.2f}{o2:>+10.4f}{r2:>+10.2f}{nv:>7d}")
rs=[v["ratio"] for v in res.values() if not v.get("undecidable")]
obs=[v["obs"] for v in res.values() if not v.get("undecidable")]
print(f"\n=== 判据 ===")
print(f"  性四题的**实测最弱一环**在五个池里:"+" ".join(f"{x:.4f}" for x in obs)+
      f"   ⇒ 展布 {max(obs)-min(obs):.4f} —— **分子几乎不动(它本来就不依赖陪衬域)**")
print(f"  性四题的**比值**:"+" ".join(f"{x:.2f}" for x in rs)+
      f"   ⇒ **max/min = {max(rs)/min(rs):.2f}**")
G=Gate("换一个陪衬域比值会挪多少")
p1=G.positive_control("警察池必须复现 #718/#730 的零 0.0911(容差 0.005)",
    planted=float(0.005-abs(q-0.0911)) if False else float(0.005-abs(res.get('警察四',{}).get('null',0.0911)-0.0911)),
    floor=0.0,spread=0.0002)
p2=G.negative_control("同池同 k 的随机题组应低于性四题",
    null=float(np.median([v["null"] for v in res.values() if not v.get("undecidable")])),
    effect=float(np.median(obs)),null_spread=0.005,
    null_kind="同一批人、`性四 ∪ 陪衬域` 这个池、同样 k=4、同样逐年取中位、同样最优符号,只打散哪四题算一组;全枚举,排除任何纯域块")
sp=max(rs)/min(rs)
if not p1: v="**UNVERIFIED:正对照没过**"
elif sp>2: v=(f"**W2:同一批人、同一组性题,只换陪衬域,比值从 {min(rs):.2f}× 走到 {max(rs):.2f}×(max/min = {sp:.2f})"
              f"—— 而实测的最弱一环几乎不动 ⇒ 比值是「相对哪个池」的,页上每一行都要带这个限定**")
elif sp<=1.5: v=f"**W1:比值 max/min = {sp:.2f} ≤ 1.5 ⇒ 池无关,页上不必改**"
else: v=f"**判不了:max/min = {sp:.2f} 落在 1.5–2 之间 ⇒ 报区间 [{min(rs):.2f}, {max(rs):.2f}]**"
print(f"\n{v}"); print(G)
json.dump(dict(cells=res,ratios=rs,obs=obs,spread=sp,verdict=v,unchallenged=True),
  open(OUT/"pool.json","w"),indent=1,ensure_ascii=False)
