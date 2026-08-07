"""#773 · E03·A41·R212 —— 唯一没到顶的那一题,它的方向站得住吗?

`#772` 判世界 C:`xmarsex`/`teensex`/`homosex` 在最虔诚层饱和到 80–88%,读数不可读。
**唯一没到顶的是 `premarsx`(端点 46.4%)**,而它是唯一方向相反的一题
(原始 ρ +0.1759 → +0.1267 → **+0.1970**)。⇒ 本轮**只在它上面**重做。

⚠⚠ **最强混淆写在跑之前,而它是压低别人那个机制的镜像**:
`premarsx` 的 sd **反向上升 41%**(0.888 → 1.102 → 1.252)。
**上升的方差会抬高相关** —— 所以「最虔诚层更强」可能只是方差变大。
⇒ 必须用**不吃 sd(Y)** 的量,而且**它们已经互相矛盾**:
`#772` 里 OLS 斜率末端更高(+0.1495→+0.1863),**边际匹配后的 ρ 末端更低**(+0.2647→+0.1539)。
**两个不吃方差的量给出相反方向,这本身就是本轮要面对的东西。**

G1 估计量:三个统计量 × 两种切法(三分位 · 五分位)的**方向**,各带**匹配同一相关矩阵的三元高斯零**。
⚠ **删失守门**:每层端点占比必须 <60%,否则该切法整格不可读(`#772` 的预注册规则,继续生效)。
⚠ **MDE 先算**:高斯零上 spread 的 95% 分位;观测低于它 ⇒ 判不了,不许跑完再解释。

预注册判词(按 `#764` 新写法):
  三个统计量在两种切法上**方向一致且超零** ⇒ 交互为真;
  **方向不一致** ⇒ **框架就是发现**(realstat §2.5),报出整张网格,不选边;
  全部未超零 ⇒ 判不了。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(212)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
VALID={"obey":(1,5),"premarsx":(1,4),"attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey","premarsx"],convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}
for c in aligned({"premarsx":cats["premarsx"]},"strict")|aligned({"obey":cats["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
D=M.dropna().copy(); D["REL"]=z(D[["attend","reliten","fund"]]).mean(axis=1)
print(f"\nn={len(D)}")
def sp(a,b): return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)),method="spearman"))
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
def matched_rho(g_list, pooled, col):
    out=[]
    for g in g_list:
        take=[]
        for cval,share in pooled.items():
            pool=g[g[col]==cval]; need=int(round(share*len(g)))
            if len(pool) and need: take.append(pool.iloc[RNG.integers(0,len(pool),need)])
        gg=pd.concat(take); out.append(sp(gg.obey,gg[col]))
    return out
GRID={}
for q,name in ((3,"三分位"),(5,"五分位")):
    D["k"]=pd.qcut(D.REL,q,labels=False,duplicates="drop")
    gs=[D[D.k==k] for k in sorted(D.k.unique())]
    tops=[float((g.premarsx==D.premarsx.max()).mean()) for g in gs]
    sds=[float(g.premarsx.std(ddof=1)) for g in gs]
    raw=[sp(g.obey,g.premarsx) for g in gs]
    sl=[slope(g.obey,g.premarsx) for g in gs]
    pooled=D.premarsx.value_counts(normalize=True)
    mm=np.median(np.array([matched_rho(gs,pooled,"premarsx") for _ in range(200)]),axis=0).tolist()
    # 高斯零:同 n、同相关矩阵
    rxy,rxz,ryz=sp(D.obey,D.premarsx),sp(D.obey,D.REL),sp(D.premarsx,D.REL)
    R=np.array([[1,rxy,rxz],[rxy,1,ryz],[rxz,ryz,1]]); w,V=np.linalg.eigh(R)
    L=V@np.diag(np.sqrt(np.clip(w,1e-9,None)))
    nul=[]
    for _ in range(400):
        G3=RNG.normal(size=(len(D),3))@L.T
        zt=pd.qcut(pd.Series(G3[:,2]),q,labels=False,duplicates="drop").to_numpy()
        rr=[sp(G3[zt==k,0],G3[zt==k,1]) for k in range(q)]
        nul.append(max(rr)-min(rr))
    q95=float(np.quantile(nul,.95))
    GRID[name]=dict(n=[len(g) for g in gs],tops=tops,sds=sds,raw=raw,slope=sl,matched=mm,
                    null_q95=q95,spread_raw=max(raw)-min(raw),spread_mm=max(mm)-min(mm),
                    spread_sl=max(sl)-min(sl))
    print(f"\n=== {name} · 端点占比 {[f'{t*100:.1f}%' for t in tops]} · sd {[round(v,3) for v in sds]} ===")
    print(f"  原始 ρ      {[round(v,4) for v in raw]}  spread {max(raw)-min(raw):.4f} · 高斯零95% {q95:.4f} ⇒ {(max(raw)-min(raw))/q95:.2f}×")
    print(f"  OLS 斜率    {[round(v,4) for v in sl]}  spread {max(sl)-min(sl):.4f}")
    print(f"  边际匹配 ρ  {[round(v,4) for v in mm]}  spread {max(mm)-min(mm):.4f} · ⇒ {(max(mm)-min(mm))/q95:.2f}×")
    print(f"  方向(首→末):原始 {'升' if raw[-1]>raw[0] else '降'} · 斜率 {'升' if sl[-1]>sl[0] else '降'} · 匹配 {'升' if mm[-1]>mm[0] else '降'}")

G=Gate("#773 · 唯一没到顶的那一题")
mt=max(max(GRID[n]["tops"]) for n in GRID)
G.asserted("① 删失守门:每层端点占比须 <60%(否则整格不可读)", bool(mt<0.60),
           f"最高端点占比 {mt*100:.1f}%(阈值 60%)", kind="control")
G.offset_control("② 原始 ρ 的 spread 须高于高斯零 95% 分位",
                 effect=GRID["三分位"]["spread_raw"], offset=GRID["三分位"]["null_q95"],
                 spread=GRID["三分位"]["null_q95"]*0.1,
                 null_kind="同 n、同相关矩阵的三元高斯世界,按同样分位切层并算同一个 spread")
print(); print(G)
dirs={n:(GRID[n]["raw"][-1]>GRID[n]["raw"][0], GRID[n]["slope"][-1]>GRID[n]["slope"][0],
         GRID[n]["matched"][-1]>GRID[n]["matched"][0]) for n in GRID}
allsame=len({v for t in dirs.values() for v in t})==1
print("\n"+"="*72)
if mt>=0.60: v=f"**删失守门未过(最高 {mt*100:.1f}%)⇒ 整轮不可读**"
elif allsame: v=f"**方向一致:六格(3 统计量 × 2 切法)全部 {'上升' if dirs['三分位'][0] else '下降'} ⇒ 交互为真**"
else:
    v=("**方向不一致 —— 框架就是发现(realstat §2.5)**:\n"
       +"\n".join(f"  {n}:原始 {'升' if a else '降'} · 斜率 {'升' if b else '降'} · 边际匹配 {'升' if c else '降'}"
                  for n,(a,b,c) in dirs.items())
       +"\n  ⇒ **三个都不吃同一样东西:原始 ρ 吃 sd(Y);斜率不吃 sd(Y) 但吃删失;边际匹配把 Y 的分布拉平。**\n"
        "  **它们分歧的地方,就是「这个交互到底定义在什么上」还没定下来。**")
print(v)
json.dump(dict(n=len(D),grid=GRID,dirs={k:list(v) for k,v in dirs.items()},verdict=v,
               gate_ok=all(r[2] for r in G.rows)),open(OUT/"premarsx.json","w"),ensure_ascii=False,indent=1)
