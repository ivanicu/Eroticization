"""E03·A30·R157 —— 「性是一块」的零:任取同样多的题,最弱一环能到多少

**类型:FRONTIER。而这是 `#714`⑤ 点名的「正结果我会不高兴」的那一步。**

**心理学的那一句(本轮要判的):同一份问卷里随便挑几道题,它们本来就会挤成一团吗?
如果会,那「性道德是一件事」就不是关于性的发现,只是关于问卷的。**

## 为什么是这一条
`#653` 的表有两列(生 / 归一),**两列都不是零 —— 它的「团/链」判决从来没有零**。
而它的三条「团」里有两条只高出那个**选定**的 0.30 阈:**NSFG 0.346(+0.046)· MFQ·PURITY 0.335(+0.035)**。
`#653` 自己第②条是**这一页第一条跨仪器复制**,页上招牌行。**若零清得掉它们,这一行当场要削。**

## G1 ESTIMAND
块内 **最弱一环(min 归一对相关)**,与 `#653` 逐字同一个量、同一条代码路径(`align`+`pairmap`,逐年取中位)。

## G2 CONTROLS —— 零 = 保住仪器,毁掉「谁跟谁同域」的指派
在**同一具仪器的题池**里随机抽同样大小的块(排除真块),算同一个最弱一环。
**`negative_control`,零的种类 = 同一份问卷、同样题数、同样的对齐与逐年取中位,只把域指派打散。**
⚠ **这个零对我不利,是故意的:** 随机块也走 `align()`,等于白送它一个自由参数去自选翻向 ⇒ **零被抬高**。
**④ 正对照**:必须复现 `#653` 的四个归一值 —— SCCS·以身作则 0.895 · GSS 0.416 · NSFG 0.346 · MFQ·PURITY 0.335。

## ⑤ 停止条件(跑之前写死,不许跑完再找理由)
- **正对照复现不到 0.001** ⇒ 记「旧值不可复现」并停,本轮不判任何团。
- **逐块:零的 95% 分位 ≥ 该块实测最弱一环 ⇒ 该块的「团」撤回**(是撤回,不是降级)。
- **若 GSS 与 NSFG 两块同时被清掉 ⇒ `#653`② 那条跨仪器复制整条倒**,页上要删。

## 最强混淆(先写下来)
混合块会把**反向题**配到一起 ⇒ 最弱一环变成大负数 ⇒ 零被人为压低 ⇒ 实测显得特别。
**控制在同一迭代内:`align()` 对真块与随机块一视同仁地翻向**,所以这条混淆被同一步吃掉;
**并且另报一版「不对齐」的零作 G4 第二格**,两版都印。

## IMPOSSIBLE(不写 planned)
MFQ 非概率样本 · GSS 逐年再取中位 · 题池只有本页用过的题(不是全问卷)⇒ **零是条件在这个池上的**。
`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate, calibrated_tolerance
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)

# ── 与 #653 逐字同一条路径 ──
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,sign=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if sign<0: y=y[::-1]
    return sp(x,y)
def align(F,items):
    Z=F[items].dropna(); m=Z.mean(axis=1)
    A=Z.copy()
    for i in items:
        if sp(Z[i],m)<0: A[i]=-A[i]
    return A
def weakest(F,items,year=None,floor=150,do_align=True):
    A=align(F,items) if do_align else F[items].dropna()
    Fa=A.join(F[[year]]) if year else A
    groups=[(None,Fa)] if year is None else list(Fa.groupby(year))
    vals=[]
    for a,b in combinations(items,2):
        per=[]
        for _,g in groups:
            m=g[[a,b]].dropna()
            if len(m)<floor or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
        if not per: return np.nan
        vals.append(float(np.median(per)))
    return float(np.min(vals))

print("=== 硬规则①:题池、n、年份 —— 先印后引 ===")
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
Xm=mfq[list(ITEM)].dropna(); print(f"  MFQ  题池 {len(ITEM)} 题 / 5 域 · n = {len(Xm)} · 无年份(单次施测)")
POL=["polabuse","polmurdr","polescap","polattak"]; SEXG=["premarsx","xmarsex","homosex","teensex"]
gss,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+POL+SEXG,encoding="latin1")
yg=sorted(gss.dropna(subset=SEXG).year.unique())
print(f"  GSS  题池 {len(POL+SEXG)} 题(警察 4 + 性 4)· 性四题非缺失 n = {len(gss.dropna(subset=SEXG)):,} · "
      f"年份 {int(yg[0])}–{int(yg[-1])}({len(yg)} 个调查年)")
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"'); LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEXN+FAMN if n in LAY}; buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
Xn=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
print(f"  NSFG 题池 {len(SEXN+FAMN)} 题(性 3 + 家庭 7)· n = {len(Xn.dropna()):,} · 2011–2013 一轮")
S=pathlib.Path("data/external/dplace/repo/datasets/SCCS")
Dd=pd.read_csv(S/"data.csv"); W=Dd.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
SCC=[f"SCCS{429+i}" for i in range(4)]
print(f"  SCCS 以身作则四对象 {SCC} · 联合 n = {len(W[SCC].dropna())}(正对照用)")

# ── ④ 正对照:复现 #653 的四个归一值 ──
LED={"SCCS·以身作则":(W,SCC,None,30,0.895),"GSS·性道德四题":(gss,SEXG,"year",150,0.416),
     "NSFG·性三题":(Xn,SEXN,None,150,0.346),"MFQ·PURITY":(Xm,[k for k,v in ITEM.items() if v=="PURITY"],None,150,0.335)}
print("\n=== ④ 正对照:与 #653 账本值逐块对照(容差 0.001)===")
rep={}
for k,(F,it,yc,fl,led) in LED.items():
    v=weakest(F,it,year=yc,floor=fl); rep[k]=v
    print(f"  {k:16s} 实测 {v:+.4f}  账本 {led:+.3f}  差 {abs(v-led):.4f}  {'✅' if abs(v-led)<=0.001 else '⛔'}")
maxd=max(abs(rep[k]-LED[k][4]) for k in LED)
if maxd>0.001:
    print(f"\n⛔ ⑤ 触发:最大差 {maxd:.4f} > 0.001 ⇒ **旧值不可复现,本轮不判任何团**")
    json.dump(dict(stop="旧值不可复现",max_diff=maxd,reproduced=rep),open(OUT/"null_for_653.json","w"),indent=1,ensure_ascii=False); sys.exit(0)

# ── 零:同池随机块 ──
rng=np.random.default_rng(20260806)
POOLS={"GSS·性道德四题":(gss,POL+SEXG,SEXG,"year",150,4,None),
       "NSFG·性三题":(Xn,SEXN+FAMN,SEXN,None,150,3,None),
       "MFQ·PURITY":(Xm,list(ITEM),[k for k,v in ITEM.items() if v=="PURITY"],None,150,6,1500)}
print("\n=== 零:同一具仪器的题池里随机抽同样多的题(排除真块)===")
res={}
for name,(F,pool,true_b,yc,fl,ksz,nsamp) in POOLS.items():
    allb=[c for c in combinations(pool,ksz) if set(c)!=set(true_b)]
    if nsamp and len(allb)>nsamp:
        idx=rng.choice(len(allb),nsamp,replace=False); blocks=[allb[i] for i in idx]; how=f"随机抽 {nsamp} / 共 {len(allb):,}"
    else:
        blocks=allb; how=f"**全枚举 {len(allb)}**"
    for tag,do_al in (("对齐",True),("不对齐",False)):
        nul=np.array([weakest(F,list(b),year=yc,floor=fl,do_align=do_al) for b in blocks])
        nul=nul[np.isfinite(nul)]
        q=float(np.quantile(nul,0.95)); obs=rep[name]
        res[f"{name}·{tag}"]=dict(obs=obs,q95=q,median=float(np.median(nul)),n_blocks=int(nul.size),how=how,
                                  cleared=bool(q>=obs))
        print(f"  {name:16s} {tag:4s} {how:22s} 零 95% 分位 **{q:+.4f}** · 中位 {np.median(nul):+.4f} · "
              f"实测 {obs:+.4f} ⇒ {'⛔ **被零清掉,撤回「团」**' if q>=obs else f'✅ 站得住({obs/q:.2f} 倍)' if q>0 else '✅ 站得住'}")

G=Gate("「性是一块」是不是任取几题都这样")
p1=G.positive_control("必须复现 #653 的四个归一值(最大差 <0.001)",planted=float(0.001-maxd),floor=0.0,spread=0.00005)
key=[k for k in res if k.endswith("·对齐")]
worst=max(res[k]["q95"]/res[k]["obs"] for k in key)
p2=G.negative_control("同池随机块的最弱一环应低于真块",
    null=float(np.mean([res[k]["q95"] for k in key])),effect=float(np.mean([res[k]["obs"] for k in key])),
    null_spread=0.005,null_kind="同一份问卷、同样题数、同样对齐与逐年取中位,只把域指派打散(对齐版)")
cleared=[k for k in key if res[k]["cleared"]]
if not p1: v="**判不了:旧值不可复现**"
elif cleared: v=f"**撤回:{'、'.join(c.split('·')[0] for c in cleared)} 的「团」落在自己的零里**"
else: v=f"**三具仪器上「性是一块」都在零之上(最紧的一块是零的 {1/worst:.2f} 倍)**"
print(f"\n{v}"); print(G)
json.dump(dict(reproduced=rep,ledger={k:LED[k][4] for k in LED},max_diff=maxd,null=res,
               verdict=v,cleared=cleared,unchallenged=True),open(OUT/"null_for_653.json","w"),indent=1,ensure_ascii=False)
