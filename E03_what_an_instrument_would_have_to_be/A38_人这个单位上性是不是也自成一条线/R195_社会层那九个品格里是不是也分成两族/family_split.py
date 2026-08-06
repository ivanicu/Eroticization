"""#753 · E03·A38·R195 —— 社会层那九个品格里,是不是也分成两族?

`#752` 在**人**层量到:性独立于伤害-公平,却是权威-内群体那一族的一部分。
`#750` 在**社会**层偏掉了九个育儿品格 —— **而那九个里只有「服从」属于权威族。**
⇒ 所以 `#750` 的 +0.5504 要**重算**,不是重述:把九个按族拆开,看代价是不是也分成两档。

⚠⚠ **最强的混淆写在跑之前:我已经看过 `#750` 的九格网格**(服从 −0.094、自我克制 −0.076 最大,
攻击性 +0.025 反向)。**所以任何我自己划的族都是被污染的。** 两条外部防御:
  ① **文献二分** —— Barry/Child/Bacon 1957–59 的 pressure toward COMPLIANCE(责任 · 服从)
     vs ASSERTION(成就 · 自立)。**它比我看那张网格早六十年,不可能被它拟合。**
     ⚠ 它只覆盖九个里的四个,这是它的代价,如实报。
  ② **盲因子拆分** —— 对九个非性品格做 PCA,**全程不看 `SCCS165`,也不看性克制**,
     按 PC1 载荷符号分两组。覆盖全部九个,而且与结果变量无关。

**预注册判词(阈值写在跑之前):**
  W-A 族的故事在社会层也成立:文献 compliance 组的代价 ≥ assertion 组的 2 倍
  W-B ⚠ **对上的是「服从」本身,不是一个族**:服从单独的代价 > 整个 compliance 组的代价
      —— 这一支会削弱跨层的漂亮故事,所以它写在前面
  W-C 两条都不满足 ⇒ 判不了
每个规格报**它自己的**零(社会间打乱 `SCCS165` 后重算偏相关)与 MDE。
"""
import pandas as pd, numpy as np, re, json, pathlib, sys
from scipy.stats import spearmanr
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
RNG=np.random.default_rng(195)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
B=ROOT/"data/external/dplace/repo/datasets/SCCS/"
V=pd.read_csv(B/"variables.csv",low_memory=False)
W=pd.read_csv(B/"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
P=V[V.source.astype(str).str.contains("barry1976traits",na=False)]
fam={}
for _,r in P.iterrows():
    m=re.match(r'(.+?):\s*(Early|Late)\s+(Boy|Girl)s?$',str(r.title))
    if m: fam.setdefault(m.group(1),[]).append(r.id)
FAM={k:sorted(v) for k,v in fam.items() if len(v)==4}
TG="SCCS165"; SR="Sexual Restraint"; NON=sorted([k for k in FAM if k!=SR])
comp=lambda cols:(W[cols].apply(lambda s:(s-s.mean())/s.std())).mean(axis=1)
C=pd.DataFrame({k:comp(v) for k,v in FAM.items()}); C[TG]=W[TG]
sp=lambda a,b: float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,ctrl):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    Cc=np.asarray(ctrl,float); Cc=Cc.reshape(-1,1) if Cc.ndim==1 else Cc
    rc=np.column_stack([r(Cc[:,j]) for j in range(Cc.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
print("九个非性品格:",NON)

# --- 防御②:盲因子拆分 —— 只用九个非性品格,不看 SCCS165,不看性克制 ---
Z=C[NON].dropna()
Zc=(Z-Z.mean())/Z.std(ddof=1)
u,s,vt=np.linalg.svd(Zc.to_numpy(),full_matrices=False)
load=pd.Series(vt[0],index=NON)
print(f"\n=== 防御②:盲因子拆分(n={len(Z)},只看九个非性品格)===")
print(f"  PC1 解释方差 {s[0]**2/np.sum(s**2):.3f}")
for k,v in load.sort_values().items(): print(f"    {k:18s} 载荷 {v:+.3f}")
gA=[k for k in NON if load[k]>=0]; gB=[k for k in NON if load[k]<0]
print(f"  组A(载荷≥0)= {gA}\n  组B(载荷<0)= {gB}")
print(f"  ⚠ 服从落在 {'组A' if 'Obedience' in gA else '组B'}")

def run(name,ctrl_cols,NPERM=2000):
    df=C[[SR,TG]+ctrl_cols].dropna()
    if len(df)<40: return None
    ct=df[ctrl_cols].to_numpy()
    obs=prho(df[SR].to_numpy(),df[TG].to_numpy(),ct)
    nul=[abs(prho(df[SR].to_numpy(),RNG.permutation(df[TG].to_numpy()),ct)) for _ in range(NPERM)]
    return dict(name=name,n=len(df),partial=obs,mde=float(np.quantile(nul,.95)),cols=ctrl_cols)

BIV=sp(*[C[[SR,TG]].dropna()[c] for c in (SR,TG)])
print(f"\n=== 双变量基线 ρ(性克制, {TG}) = {BIV:+.4f} · n={len(C[[SR,TG]].dropna())} ===")
LIT_C=[c for c in ["Responsibility","Obedience"] if c in NON]
LIT_A=[c for c in ["Achievement","Self-reliance"] if c in NON]
print(f"  ⚠ 文献二分覆盖 {len(LIT_C)+len(LIT_A)}/9 个 —— compliance {LIT_C} · assertion {LIT_A}")
SPECS={"文献 compliance(责任+服从)":LIT_C,"文献 assertion(成就+自立)":LIT_A,
       "只服从":["Obedience"],"只责任":["Responsibility"],
       f"盲因子 组A":gA,f"盲因子 组B":gB,"全部九个":NON}
res={}
print(f"\n{'规格':30s} {'n':>4s} {'偏相关':>9s} {'MDE':>7s} {'相对双变量掉了':>12s}")
for nm,cols in SPECS.items():
    r=run(nm,cols)
    if r is None: print(f"  {nm:28s}  n 太小,跳过"); continue
    res[nm]=r; drop=BIV-r["partial"]
    print(f"  {nm:28s} {r['n']:4d} {r['partial']:+9.4f} {r['mde']:7.3f} {drop:+12.4f}")

dC=BIV-res["文献 compliance(责任+服从)"]["partial"] if "文献 compliance(责任+服从)" in res else np.nan
dA=BIV-res["文献 assertion(成就+自立)"]["partial"] if "文献 assertion(成就+自立)" in res else np.nan
dO=BIV-res["只服从"]["partial"]
print("\n"+"="*66)
if dO>dC: v=(f"**W-B:服从单独的代价 {dO:+.4f} > 整个 compliance 组的 {dC:+.4f} ⇒ **对上的是「服从」本身,不是一个族** —— "
             f"跨层的族故事在社会层**不**成立,成立的是「服从」这一个品格**")
elif dA!=0 and dC>=2*abs(dA): v=(f"**W-A:compliance 组代价 {dC:+.4f} ≥ assertion 组 {dA:+.4f} 的两倍 ⇒ 族的故事在社会层也成立**")
else: v=f"**W-C:compliance {dC:+.4f} · assertion {dA:+.4f} · 只服从 {dO:+.4f} —— 落在两条预注册之间,判不了**"
print(v)
json.dump(dict(bivariate=BIV,loadings=load.to_dict(),groupA=gA,groupB=gB,
               specs={k:v for k,v in res.items()},verdict=v),
          open(OUT/"family_split.json","w"),ensure_ascii=False,indent=1)
