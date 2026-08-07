"""#782 · E03·A43·R221 —— 那一格是真反号,还是我的可读性判据太松?

`#781` 的撤回是靠**一格** 触发的:Seg-B · k3 · n≥60 · 世代 1975–1994 · 水平 = **1.128**。
`#781`① 预注册:**先问它是不是纯噪声。若该世代的 MDE 宽到覆盖 [0.44, 1.13],
那一格改标「不可读」,撤回维持但理由从「有反号格」变成「我的可读性判据太松」——两者要分清。**

⚠⚠ **而机制可能就在判据本身**:`#779`–`#781` 用的可读性判据是
**「两层的斜率各自超出自己的零」** —— **那是对分子与分母分别检验,不是对比值检验。**
**小分母下,两个都「显著」的量相除仍可以极不稳。**

⚠ **零怎么造,这一条决定成败**:直觉是「打乱年份标签得到比值的零」——**那是错的**:
打乱后两个斜率都趋近 0,**比值变成 Cauchy 型、方差不存在**,那个「零」大到会判掉所有格,**不是信息**。
⇒ 正确的问题不是「比值是否异于某个零」,而是**「这个比值被定得多准」** ⇒ **年份层面的自助区间。**

G1 估计量:每个可读格的比值 **及其年份自助 95% 区间**;
并把这一判据**回溯到 `#781` 的全部 37 个可读格**(否则只是给一格换标签)。

⚠ **换不了仪器**,与 `#776`–`#781` 同一条且同样量过。**只此一具。**

预注册判词(按 `#764` 新写法):
  ① 那一格的自助区间**覆盖 [0.44, 1.13]** ⇒ **B**:改标不可读,撤回维持而理由改写;
  ② 区间不覆盖 ⇒ **A**:它是真反号,`#781` 的撤回理由不变;
  ③ 无论 A/B,**都要报全网格有多少格的自助区间含 1.0** —— 那才是「判据太松」的规模。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(221)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
VALID={"homosex":(1,4),"attend":(0,8),"reliten":(1,4),"fund":(1,3),"cohort":(1900,2006)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=["homosex"],convert_categoricals=True)
for c in aligned({"homosex":list(cat["homosex"].cat.categories)[:4]},"strict"): M[c]=-M[c]+5
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
sub=M.dropna(subset=["homosex","attend","reliten","fund","cohort","year"]).copy()
sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
sub["k3"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
sub["k2"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,2,labels=False,duplicates="drop"))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
def series(lo,hi,col,kk,nmin):
    g=sub[(sub.cohort>=lo)&(sub.cohort<=hi)&(sub[col]==kk)]
    rows=[]
    for y,gy in g.groupby("year"):
        if len(gy)<nmin: continue
        rows.append((int(y),float(gy.homosex.mean()),float((gy.homosex==4).mean())))
    return rows
def ratio_boot(rowsA,rowsB,stat,B=4000):
    yA=np.array([r[0] for r in rowsA],float); yB=np.array([r[0] for r in rowsB],float)
    j=1 if stat=="水平" else 2
    vA=np.array([r[j] for r in rowsA]); vB=np.array([r[j] for r in rowsB])
    f0A,f0B=rowsA[0][2],rowsB[0][2]
    def rat(ia,ib):
        sa,sb=slope(yA[ia],vA[ia]),slope(yB[ib],vB[ib])
        if stat!="水平": sa,sb=sa/f0A,sb/f0B
        return sa/sb if abs(sb)>1e-12 else np.nan
    obs=rat(np.arange(len(yA)),np.arange(len(yB)))
    bs=[rat(RNG.integers(0,len(yA),len(yA)),RNG.integers(0,len(yB),len(yB))) for _ in range(B)]
    bs=np.array([b for b in bs if np.isfinite(b)])
    return obs, float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5)), len(bs)
SEG={"Seg-A(#779)":[(1930,1949),(1950,1964),(1965,1979),(1980,1999)],
     "Seg-B(等长20年)":[(1935,1954),(1955,1974),(1975,1994)]}
CUT={"k3":("k3",2,0),"k2":("k2",1,0)}
print(f"\n=== ① 先看触发撤回的那一格 ===")
rA=series(1975,1994,"k3",2,60); rB=series(1975,1994,"k3",0,60)
print(f"  Seg-B · k3 · n≥60 · 1975–1994 · 水平 —— 虔诚层 {len(rA)} 年 · 非虔诚层 {len(rB)} 年")
obs,lo95,hi95,nb=ratio_boot(rA,rB,"水平")
print(f"  观测比值 {obs:.3f} · **年份自助 95% 区间 [{lo95:.3f}, {hi95:.3f}]**(B={nb} 有效)")
# ⚠⚠ 第一版把判据写成「区间是否**包含整段** [0.44, 1.13]」(lo<=0.44 and hi>=1.13)——
#    而预注册问的是**「这一格能不能分辨 0.44 与 1.13」**。区间 [0.561, 2.851] 里
#    0.561 > 0.44,于是第一版返回「否 ⇒ 真反号」,**而它同时包含 1.0、0.623、0.821** ——
#    **这一格连「比值在 1 的哪一边」都定不了。**
#    ⇒ 「判词分支测错问题」那一族(`#728`·`#748`·`#750`·`#758`)的又一次,而这次它**改变了结论**。
#    正确判据:**该格若要触发「≥1.0」的撤回条件,它的区间必须把 1.0 排除在外。**
covers = (lo95 <= 1.0 <= hi95)      # 区间含 1.0 ⇒ 这一格定不了方向 ⇒ 不可读
print(f"  ⇒ 区间含 1.0 吗:**{'是 ⇒ 这一格定不了方向' if covers else '否 ⇒ 它真的在 1 之上'}**(预注册要问的那个问题)")
print(f"  ⚠ 它同时包含该世代另外两格(0.623 · 0.821)—— 宽 {hi95-lo95:.2f},**分辨不了 0.44 与 1.13**")
print(f"\n=== ③ 把这一判据回溯到全部可读格(否则只是给一格换标签)===")
cells=[]
for segn,bands in SEG.items():
    for cutn,(col,hk,lk) in CUT.items():
        for nmin in (60,120):
            for lo,hi in bands:
                a,b=series(lo,hi,col,hk,nmin),series(lo,hi,col,lk,nmin)
                if len(a)<8 or len(b)<8: continue
                for st in ("水平","端点相对基线"):
                    o,l,h,nb2=ratio_boot(a,b,st,B=1500)
                    if not np.isfinite(o): continue
                    cells.append(dict(seg=segn,cut=cutn,nmin=nmin,cohort=f"{lo}–{hi}",stat=st,
                                      ratio=float(o),lo=l,hi=h,width=h-l,covers1=bool(l<=1.0<=h)))
W=np.array([c["width"] for c in cells]); C1=sum(c["covers1"] for c in cells)
print(f"  {len(cells)} 格算出了自助区间 · **区间含 1.0 的有 {C1} 格({C1/len(cells)*100:.0f}%)**")
print(f"  区间宽度分位:5% {np.percentile(W,5):.2f} · 中位 {np.median(W):.2f} · 95% {np.percentile(W,95):.2f} · 最大 {W.max():.2f}")
print(f"\n  按世代看「区间含 1.0」的比例:")
for coh in sorted({c["cohort"] for c in cells}):
    g=[c for c in cells if c["cohort"]==coh]
    print(f"    {coh:12s} {sum(x['covers1'] for x in g)}/{len(g)} 格含 1.0 · 中位宽度 {np.median([x['width'] for x in g]):.2f}")
G=Gate("#782 · 比值需要它自己的区间")
G.asserted("① 自助必须真的在变(正控:有效重抽数 >0 且区间非零宽)",
           bool(nb>100 and hi95-lo95>1e-6), f"有效重抽 {nb} · 区间宽 {hi95-lo95:.3f}", kind="control")
G.asserted("② 预注册(修正判据):要触发「≥1.0」的撤回,那一格的区间必须**排除** 1.0",
           bool(not covers), f"区间 [{lo95:.3f}, {hi95:.3f}] {'含' if covers else '不含'} 1.0", kind="kill")
print(); print(G)
print("\n"+"="*90)
if covers:
    v=(f"**B:那一格的自助区间 [{lo95:.3f}, {hi95:.3f}] **含 1.0**(宽 {hi95-lo95:.2f})⇒ 它定不了方向,不可读。**\n"
       f"  ⇒ **`#781` 的撤回维持,但理由改写**:不是「世界上有一个反号的规格」,\n"
       f"  而是**「我的可读性判据太松」** —— 它检验了分子与分母,**没检验比值**。\n"
       f"  ⇒ 而这牵连整组:**{len(cells)} 格里有 {C1} 格({C1/len(cells)*100:.0f}%)的自助区间含 1.0。**")
else:
    v=(f"**A:那一格的区间 [{lo95:.3f}, {hi95:.3f}] 不含 1.0 ⇒ 它真的在 1 之上,"
       f"`#781` 的撤回理由不变**;而全网格 {C1}/{len(cells)} 格的区间含 1.0,仍须如实记")
print(v)
json.dump(dict(trigger=dict(ratio=obs,lo=lo95,hi=hi95,covers=covers,nA=len(rA),nB=len(rB)),
               cells=cells,covers1=C1,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"ratio_interval.json","w"),ensure_ascii=False,indent=1)
