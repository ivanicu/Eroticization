"""#771 · E03·A41·R210 —— 交互是不是也只是二阶矩?

`#770` 证明「解释份额」是推导:它是 (r_XY, r_XZ, r_YZ) 的确定函数,
而一个**只匹配那两条相关的高斯**给出同一个份额(差 0.89pp),尽管偏度相差 0.38。
`#770`② 预注册的下一问:**那么这条线上什么才是非推导量?**
候选:**交互** —— 宗教是不是改变 `obey↔性态度` 的**斜率**,而不只是平移。

⚠⚠ **而这个问题自带一个漂亮的零,这是本轮设计的全部**:
**匹配全部二阶矩的高斯合成世界**。二元/三元高斯里,条件相关**不随条件变量取值改变** ——
所以若真实数据的分层相关差异**超过**高斯合成的同一统计量,那部分**不是二阶矩的函数**,
**即:它是一个份额永远看不见的东西。** 若它没超过,**交互也是推导**,这条线上就真的只剩相关。

G1 估计量:`spread = max_k ρ_k(obey, 态度) − min_k ρ_k`,k 为宗教三分位。
**零 = 同一 n、同一相关矩阵的高斯合成世界上的同一统计量。**

⚠ **最强混淆写在跑之前:分层会造出全距截断** —— 某一层若态度方差小,相关会被机械压低。
   而**高斯合成零同样带这个截断**(它按同一相关矩阵生成、同样分层),
   ⇒ **这个零把「分层本身造成的差异」减掉了,而那正是我需要的。**
   ⚠ 但它减不掉**方差随层变化的真实幅度差异** —— 所以另报**每层的 sd**,如实列出。

⚠ **MDE 先算(`#746`①)**:在高斯零上量出 spread 的 95% 分位;
   若观测 spread 落在它之下,**本设计判不了,不许跑完再解释。**
预注册判词(按 `#764` 新写法:只比已测量的量,各带自己的零):
  观测 spread 与**高斯零的 95% 分位**并排报;超出 ⇒ 交互不是二阶矩的函数(非推导量);
  未超出 ⇒ **交互也是推导**,这条线上只剩相关本身。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(210)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"obey":(1,5),"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
cat=pd.read_stata(gp,columns=["obey"]+SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
for c in aligned({c:cats[c] for c in SEX},"strict")|aligned({c:cats[c] for c in ["obey"]},"important"): M[c]=-M[c]
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
def sp(a,b): return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)),method="spearman"))

def strata_spread(x,y,zv,q=3):
    cut=pd.qcut(pd.Series(zv),q,labels=False,duplicates="drop")
    rs=[]; sds=[]
    for k in sorted(pd.Series(cut).dropna().unique()):
        m=(cut==k).to_numpy()
        if m.sum()<200: continue
        rs.append(sp(np.asarray(x)[m],np.asarray(y)[m])); sds.append(float(np.std(np.asarray(y)[m])))
    return (max(rs)-min(rs), rs, sds)

print(f"\n=== 观测:宗教三分位内的 ρ(obey, 态度) ===")
print(f"  {'题':10s}{'n':>8s}{'低':>9s}{'中':>9s}{'高':>9s}{'spread':>9s}   每层 sd")
obs={}
for s in SEX:
    sub=M[["obey","attend","reliten","fund",s]].dropna().copy()
    sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
    sp_,rs,sds=strata_spread(sub.obey.to_numpy(),sub[s].to_numpy(),sub.REL.to_numpy())
    obs[s]=dict(n=len(sub),rs=rs,spread=sp_,sds=sds,
                rxy=sp(sub.obey.to_numpy(),sub[s].to_numpy()),
                rxz=sp(sub.obey.to_numpy(),sub.REL.to_numpy()),
                ryz=sp(sub[s].to_numpy(),sub.REL.to_numpy()))
    print(f"  {s:10s}{len(sub):8d}"+"".join(f"{r:+9.4f}" for r in rs)+f"{sp_:9.4f}   {[round(v,3) for v in sds]}")

print(f"\n=== 零 = 匹配同一相关矩阵的高斯合成世界(同 n、同分层做法)===")
B=400
print(f"  {'题':10s}{'观测 spread':>12s}{'零中位':>9s}{'零 95%':>9s}{'倍数':>8s}")
res={}
for s in SEX:
    o=obs[s]; n=o["n"]
    R=np.array([[1,o["rxy"],o["rxz"]],[o["rxy"],1,o["ryz"]],[o["rxz"],o["ryz"],1]])
    w,V=np.linalg.eigh(R); w=np.clip(w,1e-9,None); L=V@np.diag(np.sqrt(w))
    nul=[]
    for _ in range(B):
        G3=RNG.normal(size=(n,3))@L.T
        nul.append(strata_spread(G3[:,0],G3[:,1],G3[:,2])[0])
    q95=float(np.quantile(nul,.95)); med=float(np.median(nul))
    res[s]=dict(spread=o["spread"],null_med=med,null_q95=q95,ratio=o["spread"]/q95)
    print(f"  {s:10s}{o['spread']:12.4f}{med:9.4f}{q95:9.4f}{o['spread']/q95:7.2f}×")

G=Gate("#771 · 交互是不是也只是二阶矩")
G.identity_control("① 高斯零的相关矩阵须复现目标(合成世界搭对了吗)",
                   observed=float(np.round(sp(*(lambda A: (A[:,0],A[:,1]))(
                       RNG.normal(size=(obs['premarsx']['n'],3))@(lambda R:(lambda w,V:V@np.diag(np.sqrt(np.clip(w,1e-9,None))))(*np.linalg.eigh(R)))(
                       np.array([[1,obs['premarsx']['rxy'],obs['premarsx']['rxz']],
                                 [obs['premarsx']['rxy'],1,obs['premarsx']['ryz']],
                                 [obs['premarsx']['rxz'],obs['premarsx']['ryz'],1]])).T)),2)),
                   expected=round(obs['premarsx']['rxy'],2), tol=0.02,
                   what="合成世界的 r_XY 必须等于目标 r_XY,否则零不是同一个世界")
mx=max(res[s]["ratio"] for s in SEX)
G.offset_control("② 观测 spread 须显著高于高斯零,才算交互不是二阶矩的函数",
                 effect=max(res[s]["spread"] for s in SEX),
                 offset=max(res[s]["null_q95"] for s in SEX),
                 spread=float(np.std([res[s]["null_q95"] for s in SEX])),
                 null_kind="同 n、同相关矩阵的三元高斯世界,按同样的三分位分层并算同一个 spread —— 它把「分层本身造成的差异」减掉")
print(); print(G)
print("\n"+"="*72)
above=[s for s in SEX if res[s]["ratio"]>1.0]
if len(above)==len(SEX):
    v=(f"**四题全部超出高斯零的 95% 分位(倍数 {min(res[s]['ratio'] for s in SEX):.2f}×–{mx:.2f}×)"
       f" ⇒ 交互不是二阶矩的函数 —— 这条线上存在份额看不见的东西**")
elif not above:
    v=(f"**四题全部落在高斯零之下(最大 {mx:.2f}×)⇒ 交互也是推导,这条线上只剩相关本身**")
else:
    v=(f"**{len(above)}/4 题超出(倍数 {min(res[s]['ratio'] for s in SEX):.2f}×–{mx:.2f}×)⇒ 不是普遍的,整张网格已全列**")
print(v)
json.dump(dict(obs={s:obs[s] for s in SEX},res=res,verdict=v,
               gate_ok=all(r[2] for r in G.rows)),open(OUT/"interaction.json","w"),ensure_ascii=False,indent=1)
