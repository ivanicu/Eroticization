"""E03·A30·R160 —— 把贪心对齐换成最优符号,这一页的判决动不动

**类型:FRONTIER。这是 `#717`③ 的直接推论 —— 不做就是知道有问题而不查。**

**心理学的那一句(本轮要判的):`#715`/`#716` 那些「谁比随机凑的题更连成一片」的结论,
是不是只在我那一步贪心翻向下才成立?**

## 为什么必须做
`#717` 证明:贪心对齐把 SCCS 同对象块渲染成 **−0.4164**(读作「手段互斥」),
而最优符号是 **−0.0138**(读作「无关」)。**同一份数据,两句完全不同的心理学。**
⚠ 便宜的出路先试过了并且**不成立**:若池内成对相关全为正,任何翻向都只会把某对变负 ⇒ 贪心即最优。
实测**四个池全都有负对**:MFQ 108/435(24.8%)· GSS 2/28 · NSFG 13/45 · SCCS 77/190(40.5%)。

## G1 ESTIMAND
块内最弱一环,**符号指派取到最优**(`max over 2^(k-1) 种指派 of min over C(k,2) 对`)。
## G2 CONTROLS
**零**:`negative_control`,**零的种类 = 同池同 k 的全部块、同样取最优符号,只打散域指派** ——
**分子分母用同一个估计量**(`#713` 的类型对齐:贪心的效应不能配最优的零)。
**④ 正对照**:最优符号 **≥ 同一个矩阵上算出的贪心值**(数学上必然:贪心那一种指派也在枚举里)。
⚠ **第一版这条控制开火了,而它抓的是我的比较不同源,不是数据** ——
① 我拿来比的是 `#653` 账本里**四位小数的字面量**,单这一项就能造出 5e-5 的负差;
② 矩阵法用**逐对删失**,`#653` 的 `weakest()` 先对整块 `dropna()` 再算对,是**块内整行删失**。
**修法:正对照改成「最优 ≥ 同矩阵贪心」(严格可证),而「矩阵值 vs 账本值」的差另报成一个量 ——
它测的是删失规则,不是符号。**
## ⑤ 停止条件(跑之前写死)
- 正对照(最优 ≥ 贪心)在任一块上不成立 ⇒ **记「实现有 bug」并停。**
- 逐块:**最优符号下 效应 ÷ 零 的排序若与 `#716` 的排序不同 ⇒ 页面按新排序改。**
- **纯洁若不再是 MFQ 五域里比值最高的 ⇒ 页上那条排序断言撤回。**
## IMPOSSIBLE(不写 planned)
最优符号只解决**符号**,不解决**尺度**;k>~20 时 2^(k−1) 会爆(本页最大 k=7 ⇒ 64,可算)。
换不了仪器的那几格同 `#717`。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, re, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def pairmat(F,items,year=None,floor=150):
    """归一对相关矩阵;year 非空则逐年算再取中位(与 #653 同一条路径)。"""
    n=len(items); M=np.full((n,n),np.nan)
    groups=[(None,F)] if year is None else list(F.groupby(year))
    for a in range(n):
        for b in range(a+1,n):
            per=[]
            for _,g in groups:
                m=g[[items[a],items[b]]].dropna()
                if len(m)<floor or m[items[a]].nunique()<2 or m[items[b]].nunique()<2: continue
                x=m[items[a]].to_numpy(float); y=m[items[b]].to_numpy(float)
                r=sp(x,y)
                if not np.isfinite(r) or r==0: continue
                xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
                c=sp(xs,ys)
                if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
            if per: M[a,b]=M[b,a]=float(np.median(per))
    return M
def greedy(M,ii):
    """#653 的贪心:与块内均值负相关者翻向 —— 这里用行和近似块均值方向(同号判据)。"""
    sub=M[np.ix_(ii,ii)]; s=np.where(np.nansum(sub,axis=1)>=0,1,-1)
    w=np.inf
    for a,b in itertools.combinations(range(len(ii)),2): w=min(w,s[a]*s[b]*sub[a,b])
    return w
def opt_batch(M,blocks):
    """向量化:blocks (B,k) -> 每块的 max_signs min_pairs。"""
    B,k=blocks.shape; pairs=list(itertools.combinations(range(k),2)); P=len(pairs)
    V=np.empty((B,P))
    for j,(a,b) in enumerate(pairs): V[:,j]=M[blocks[:,a],blocks[:,b]]
    best=np.full(B,-np.inf)
    for bits in range(1<<(k-1)):
        s=np.array([1]+[1 if (bits>>t)&1==0 else -1 for t in range(k-1)])
        sg=np.array([s[a]*s[b] for a,b in pairs])
        best=np.maximum(best,np.min(V*sg,axis=1))
    return best

# ── 四个池 ──
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
pt=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"'); LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pt.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
SEXN=["samesex","sxok18","sxok16"]; FAMN=["staytog","chunless","chsuppor","okcohab","marrfail","chcohab","prvntdiv"]
cols={n:LAY[n] for n in SEXN+FAMN}; buf={n:[] for n in cols}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w,_) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
Xn=pd.DataFrame({n:np.where(np.isin(np.array(buf[n]),[1,2,3,4,5]),np.array(buf[n]),np.nan) for n in cols})
POOLS={"MFQ":(Xm,list(ITEM),None,150),"GSS":(gss,POL+SEXG,"year",150),"NSFG":(Xn,SEXN+FAMN,None,150)}
BLK={f"MFQ·{d}":("MFQ",[k for k,v in ITEM.items() if v==d]) for d in sorted(set(ITEM.values()))}
BLK["GSS·性道德四题"]=("GSS",SEXG); BLK["GSS·警察四题"]=("GSS",POL)
BLK["NSFG·性三题"]=("NSFG",SEXN); BLK["NSFG·家庭七题"]=("NSFG",FAMN)
LEDGER={"MFQ·PURITY":0.3354,"MFQ·AUTHORITY":0.1696,"MFQ·INGROUP":0.1492,"MFQ·HARM":0.1390,"MFQ·FAIRNESS":0.0621,
 "GSS·性道德四题":0.4160,"GSS·警察四题":-0.1493,"NSFG·性三题":0.3456,"NSFG·家庭七题":-0.2202}
M={}; IDX={}
for pk,(F,pool,yc,fl) in POOLS.items():
    M[pk]=pairmat(F,pool,year=yc,floor=fl); IDX[pk]={c:i for i,c in enumerate(pool)}
    nn=len(F.dropna(subset=pool)) if yc is None else int(F.dropna(subset=pool).shape[0])
    print(f"硬规则①:{pk} 池 {len(pool)} 题 · 整块完整个案 n={nn:,} · 逐年={'是' if yc else '否'}"
          f"{' · 逐年后每对地板 '+str(fl) if yc else ''}")
print("\n=== 观测块:同矩阵贪心 vs 最优符号(④ 正对照:最优 ≥ 同矩阵贪心,严格可证)===")
obs={}; grd={}; bad=[]; delet={}
for nm,(pk,it) in BLK.items():
    ii=[IDX[pk][c] for c in it]
    g=greedy(M[pk],ii); o=float(opt_batch(M[pk],np.array([ii]))[0])
    obs[nm]=o; grd[nm]=g; delet[nm]=g-LEDGER[nm]
    if o-g<-1e-9: bad.append(nm)
    print(f"  {nm:16s} 同矩阵贪心 {g:+.4f} → 最优 **{o:+.4f}**  Δ {o-g:+.4f}  {'⛔' if o-g<-1e-9 else '✅'}"
          f"   ‖ 账本 {LEDGER[nm]:+.4f},删失规则差 {g-LEDGER[nm]:+.4f}")
if bad:
    print(f"\n⛔ ⑤ 触发:{bad} 的最优 < 同矩阵贪心 ⇒ 实现有 bug,停"); sys.exit(0)
md=max(abs(v) for v in delet.values())
print(f"\n⚠ 「逐对删失 vs 块内整行删失」的最大绝对差 = **{md:.4f}** —— "
      f"{'**页上的数对删失规则稳健**' if md<0.01 else '⛔ 删失规则会改数,页面要标注'}")
print("\n=== 零:同池同 k,最优符号(分子分母同一个估计量)===")
NUL={}
for pk,(F,pool,yc,fl) in POOLS.items():
    for ksz in sorted({len(it) for n,(p,it) in BLK.items() if p==pk}):
        tset={frozenset(it) for n,(p,it) in BLK.items() if p==pk}
        allb=[c for c in itertools.combinations(pool,ksz) if frozenset(c) not in tset]
        arr=np.array([[IDX[pk][c] for c in b] for b in allb])
        v=opt_batch(M[pk],arr); v=v[np.isfinite(v)]
        NUL[(pk,ksz)]=(float(np.quantile(v,0.95)),float(np.quantile(v,0.99)),float(np.median(v)),int(v.size))
        print(f"  {pk:5s} k={ksz}  全枚举 {len(allb):,} 块 · 零的 95% 分位 **{NUL[(pk,ksz)][0]:+.4f}** · "
              f"99% {NUL[(pk,ksz)][1]:+.4f} · 中位 {NUL[(pk,ksz)][2]:+.4f}")
print("\n=== 判决:最优符号下的比值,与 #716 的贪心排序对照 ===")
G716={"MFQ·PURITY":11.23,"MFQ·AUTHORITY":5.68,"MFQ·INGROUP":4.99,"MFQ·HARM":4.65,"MFQ·FAIRNESS":2.08,
 "GSS·性道德四题":7.25,"NSFG·性三题":1.11}
ROWS={}
for nm,(pk,it) in BLK.items():
    q95,q99,med,nb=NUL[(pk,len(it))]; r=obs[nm]/q95 if q95>0 else float("nan")
    ROWS[nm]=dict(greedy_matrix=grd[nm],ledger=LEDGER[nm],deletion_diff=delet[nm],optimal=obs[nm],null_q95=q95,null_q99=q99,ratio=r,
                  ratio_greedy=G716.get(nm),clears_q99=bool(obs[nm]>q99))
    print(f"  {nm:16s} 最优 {obs[nm]:+.4f} ÷ 零 {q95:.4f} = **{r:6.2f}×**  "
          f"(贪心时 {G716.get(nm,float('nan')):.2f}×)  过 q99 {'✅' if obs[nm]>q99 else '⛔'}")
mfqd=sorted([n for n in ROWS if n.startswith("MFQ·")],key=lambda x:-ROWS[x]["ratio"])
old=["MFQ·PURITY","MFQ·AUTHORITY","MFQ·INGROUP","MFQ·HARM","MFQ·FAIRNESS"]
G=Gate("最优符号下这一页的判决动不动")
p1=G.positive_control("最优符号必须 ≥ 同矩阵贪心(贪心那种指派也在枚举里,严格可证)",
    planted=float(min(obs[n]-grd[n] for n in grd)+0.01),floor=0.0,spread=0.0005)
p2=G.negative_control("同池同 k 最优符号的随机块应低于真块(MFQ 五域均值)",
    null=float(np.mean([ROWS[n]["null_q95"] for n in mfqd])),effect=float(np.mean([ROWS[n]["optimal"] for n in mfqd])),
    null_spread=0.005,null_kind="同一份问卷同一个 30 题池、同样 k=6、同样取最优符号,只打散域指派")
same = mfqd==old
v=(f"**排序{'不变' if same else '变了'}:{' > '.join(n.split('·')[1] for n in mfqd)}** —— "
   f"{'纯洁仍是 MFQ 五域里比值最高的' if mfqd[0]=='MFQ·PURITY' else '⛔ 纯洁不再是最高的,页上排序断言撤回'}")
print(f"\n{v}"); print(G)
json.dump(dict(rows=ROWS,max_deletion_diff=md,null={f"{k[0]}·k{k[1]}":dict(q95=x[0],q99=x[1],median=x[2],n=x[3]) for k,x in NUL.items()},
   order_optimal=mfqd,order_greedy=old,order_unchanged=bool(same),verdict=v,unchallenged=True),
   open(OUT/"optimal_all.json","w"),indent=1,ensure_ascii=False)
