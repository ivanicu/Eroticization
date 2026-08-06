"""E03·A30·R158 —— 「五个道德域里只有纯洁是团」是发现,还是那个 0.30 阈造出来的

**类型:FRONTIER。这是 `#715`② 点名的、方向对我不利的那一步 —— 它的正结果会削掉一条挂着五条声明的行。**

**心理学的那一句:五个道德域里,真的只有「纯洁」整个连成一件事吗?
还是说五个域都比随机凑的六道题更连成一片,只是纯洁最紧 —— 而「只有纯洁」是那把尺子造出来的?**

## 为什么现在做
`#715` 量出:**最弱一环随块大小系统性变化**(同池 k=3 零 0.3148 → k=6 零 0.1029),
而 `#653` 拿**一个固定的 0.30 阈**判了 k=3/4/6/7 的块。
⚠ **MFQ 五域都是 k=6,而 k=6 在 30 题池上的零只有 0.0162 —— 四条「链」是 0.062–0.170,全都高于它。**

## W1 / W2(预测矩阵)
| | 纯洁 vs 零 | 其余四域 vs 零 | 页上那行 |
|---|---|---|---|
| **W1 纯洁独一份** | 远高于零 | **落在零里** | 站得住 |
| **W2 是那个阈造的** | 高于零 | **也都高于零,只是低一些** | **「只有」要撤回,改成排序** |

**W2 的正结果是我不高兴的那个** —— `#653`④ 明写「这一页有五条声明立在纯洁是一个域上」。

## G1 ESTIMAND
每块的**最弱一环(归一)**÷ **它自己那一格的零**,而零是 **同题池 + 同块大小** 的随机块。
## G2 CONTROLS
**④ 正对照**:必须复现 `#653` 全表 14 块的归一最弱一环(容差 0.001)。
**零** = `negative_control`,**零的种类 = 同一具仪器的同一个题池、同样的块大小、同样的对齐,只把域指派打散**。
⚠ **① 先问池是否同一**(`#715`① 的要求):池不同则零不可互比,四具仪器各算各的。
## G3/G4:14 块全报,含不支持结论的;比值与阈两种判法并列。
## ⑤ 停止条件(跑之前写死)
- 正对照复现不到 0.001 ⇒ 记「旧值不可复现」并停。
- **若 MFQ 五域全部 > 各自的零 ⇒ 「只有纯洁是团」撤回**,改成「五域都高于零,纯洁最紧」+ 比值排序。
- **若其中任何一域 ≤ 零 ⇒ 「只有纯洁」在那一域上得到支持**,如实报。
## IMPOSSIBLE(不写 planned)
池只到本页用过的题 · MFQ 非概率 · SCCS 同一编码团队(池内 20 个变量全出自 Barry 1977)。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float)); return sp(x,y[::-1] if s<0 else y)
def weakest(F,items,year=None,floor=150):
    A=F[items].dropna(); m=A.mean(axis=1); A=A.copy()
    for i in items:
        if sp(A[i],m)<0: A[i]=-A[i]
    Fa=A.join(F[[year]]) if year else A
    groups=[(None,Fa)] if year is None else list(Fa.groupby(year))
    vals=[]
    for a,b in combinations(items,2):
        per=[]
        for _,g in groups:
            mm=g[[a,b]].dropna()
            if len(mm)<floor or mm[a].nunique()<2 or mm[b].nunique()<2: continue
            r=sp(mm[a],mm[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(mm[a],mm[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
        if not per: return np.nan
        vals.append(float(np.median(per)))
    return float(np.min(vals))

# ── 数据 ──
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
Xm=mfq[list(ITEM)].dropna()
POL=["polabuse","polmurdr","polescap","polattak"]; SEXG=["premarsx","xmarsex","homosex","teensex"]
gss,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+POL+SEXG,encoding="latin1")
NS=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"'); LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEXN+FAMN}; buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
Xn=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
S=pathlib.Path("data/external/dplace/repo/datasets/SCCS")
W=pd.read_csv(S/"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
SC={nm:[f"SCCS{b+i}" for i in range(4)] for nm,b in [("体罚",453),("以身作则",429),("讲课",437),("放任",465),("疼爱",469)]}
SCPOOL=[v for L in SC.values() for v in L]

# ── ①(#715 的要求):池是否同一 ──
print("=== ① 题池核对:池不同则零不可互比 ===")
POOL={"MFQ":(Xm,list(ITEM),None,150),"GSS":(gss,POL+SEXG,"year",150),"NSFG":(Xn,SEXN+FAMN,None,150),"SCCS":(W,SCPOOL,None,30)}
for k,(F,p,yc,fl) in POOL.items(): print(f"  {k:5s} 池 {len(p):2d} 题 · 逐年={'是' if yc else '否'} · 每对地板 {fl}")
BLK={}
for d in sorted(set(ITEM.values())): BLK[f"MFQ·{d}"]=("MFQ",[k for k,v in ITEM.items() if v==d])
BLK["GSS·性道德四题"]=("GSS",SEXG); BLK["GSS·警察四题"]=("GSS",POL)
BLK["NSFG·性三题"]=("NSFG",SEXN); BLK["NSFG·家庭七题"]=("NSFG",FAMN)
for nm,it in SC.items(): BLK[f"SCCS·{nm}四对象"]=("SCCS",it)
LED={"MFQ·PURITY":0.335,"MFQ·AUTHORITY":0.170,"MFQ·INGROUP":0.149,"MFQ·HARM":0.139,"MFQ·FAIRNESS":0.062,
 "GSS·性道德四题":0.416,"GSS·警察四题":-0.149,"NSFG·性三题":0.346,"NSFG·家庭七题":-0.220,
 "SCCS·以身作则四对象":0.895,"SCCS·体罚四对象":0.822,"SCCS·讲课四对象":0.786,"SCCS·疼爱四对象":0.785,"SCCS·放任四对象":0.651}
print("\n=== ④ 正对照:复现 #653 全表 14 块(容差 0.001)===")
obs={}
for nm,(pk,it) in BLK.items():
    F,_,yc,fl=POOL[pk]; obs[nm]=weakest(F,it,year=yc,floor=fl)
    print(f"  {nm:20s} 实测 {obs[nm]:+.4f}  账本 {LED[nm]:+.3f}  差 {abs(obs[nm]-LED[nm]):.4f}  {'✅' if abs(obs[nm]-LED[nm])<=0.001 else '⛔'}")
maxd=max(abs(obs[n]-LED[n]) for n in LED)
if maxd>0.001:
    print(f"\n⛔ ⑤ 触发:最大差 {maxd:.4f} > 0.001 ⇒ **旧值不可复现,停**")
    json.dump(dict(stop="旧值不可复现",max_diff=maxd,obs=obs),open(OUT/"size_matched.json","w"),indent=1,ensure_ascii=False); sys.exit(0)

# ── 零:同池 + 同 k,排除全部真块 ──
rng=np.random.default_rng(20260806); NUL={}; SAMP={("MFQ",6):3000,("SCCS",4):2500}
true={pk:[frozenset(it) for n,(p,it) in BLK.items() if p==pk] for pk in POOL}
print("\n=== 零:同池 + 同块大小,排除全部真块 ===")
for pk,(F,pool,yc,fl) in POOL.items():
    for ksz in sorted({len(it) for n,(p,it) in BLK.items() if p==pk}):
        allb=[c for c in combinations(pool,ksz) if frozenset(c) not in true[pk]]
        ns=SAMP.get((pk,ksz))
        if ns and len(allb)>ns:
            idx=rng.choice(len(allb),ns,replace=False); blocks=[allb[i] for i in idx]; how=f"随机 {ns:,}/{len(allb):,}"
        else: blocks=allb; how=f"**全枚举 {len(allb):,}**"
        v=np.array([weakest(F,list(b),year=yc,floor=fl) for b in blocks]); v=v[np.isfinite(v)]
        NUL[(pk,ksz)]=(float(np.quantile(v,0.95)),float(np.median(v)),int(v.size),how)
        print(f"  {pk:5s} k={ksz}  {how:22s} 零 95% 分位 **{NUL[(pk,ksz)][0]:+.4f}** · 中位 {NUL[(pk,ksz)][1]:+.4f}")

print("\n=== G3/G4:14 块全报 —— 固定阈判法 vs 同池同 k 的零判法 ===")
print(f"{'块':22s}{'最弱一环':>9s}{'它自己的零':>11s}{'比值':>8s}   {'0.30 阈':>8s}  {'对照零':>8s}")
ROWS={}
for nm,(pk,it) in BLK.items():
    q=NUL[(pk,len(it))][0]; o=obs[nm]; rat=o/q if q>0 else float("inf")
    old="团" if o>=0.30 else "链"; new="高于零" if o>q else "落在零里"
    ROWS[nm]=dict(obs=o,null_q95=q,ratio=rat,old=old,new=new,pool=pk,k=len(it))
    flag=" ⚠判决不同" if (old=="团")!=(new=="高于零") else ""
    print(f"{nm:22s}{o:>+9.4f}{q:>+11.4f}{rat:>8.2f}   {old:>8s}  {new:>8s}{flag}")

mfqd=[n for n in ROWS if n.startswith("MFQ·")]
above=[n for n in mfqd if ROWS[n]["obs"]>ROWS[n]["null_q95"]]
G=Gate("「只有纯洁是团」是发现还是那个阈造的")
p1=G.positive_control("必须复现 #653 全表 14 块(最大差 <0.001)",planted=float(0.001-maxd),floor=0.0,spread=0.00005)
p2=G.negative_control("同池同 k 的随机块应低于真块(以 MFQ 五域均值计)",
    null=float(np.mean([ROWS[n]["null_q95"] for n in mfqd])),effect=float(np.mean([ROWS[n]["obs"] for n in mfqd])),
    null_spread=0.005,null_kind="同一份问卷的同一个 30 题池、同样 k=6、同样对齐,只把域指派打散")
if not p1: v="**判不了:旧值不可复现**"
elif len(above)==5:
    order=" > ".join(f"{n.split('·')[1]} {ROWS[n]['ratio']:.1f}×" for n in sorted(mfqd,key=lambda x:-ROWS[x]['ratio']))
    v=f"**撤回「只有纯洁是团」:五域全部高于自己的零,差的是程度不是种类 —— {order}**"
else: v=f"**「只有纯洁」在 {5-len(above)} 个域上得到支持(它们落在零里)**"
print(f"\n{v}"); print(G)
json.dump(dict(rows=ROWS,null={f"{k[0]}·k{k[1]}":dict(q95=v0[0],median=v0[1],n=v0[2],how=v0[3]) for k,v0 in NUL.items()},
               ledger=LED,max_diff=maxd,verdict=v,unchallenged=True),open(OUT/"size_matched.json","w"),indent=1,ensure_ascii=False)
