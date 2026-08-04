import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A54 R268 -- 起始年龄矩阵里,有几个人层维度

`#222c`:题目侧卡在 31 维的分辨率上。**但人的那一侧从没被问过。**
`Δ` 用的是每个人的**斜率** `rho_i`;而 ALS 同时给出每个人在**最佳方向**上的载荷 `u_i` ——
**`u_i` 从没被拿去对过任何外部变量。**

ESTIMAND        把 `u_i` 当人层量,跑 `#244` 的 20 道 Likert 面板 + `#251` 的 11 个非情色字段;
                并报 `corr(u_i, rho_i)`。
KILL            **若 `u_i` 越阈值的题与 `rho_i` 的那 7 道基本不重叠(交集 ≤2)->
                起始年龄矩阵里有两个人层维度,而本项目只用过一个。**
NEGATIVE CTRL   每个结局在分析样本内打乱(`#184b`)。
POSITIVE CTRL   合成一个由 `u_i` 造的结局 -> 必须强测到。
IMPOSSIBLE      `u_i` 与 `rho_i` 都来自同一个残差矩阵,**不可能真正独立** ——
                本轮判的是它们**指向的外部结局**重不重叠,不是它们本身正交不正交。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
KEEP0=(np.isfinite(V0).sum(1)>=8)
_,RHO=betas(V0)

def demean_np(Aa,iters=200,tol=1e-10):
    D=np.where(np.isfinite(Aa),Aa,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def als_ux(Aa,seed=0,iters=300):
    D=demean_np(Aa); W=np.isfinite(D)&KEEP0[:,None]; Z=np.where(W,D,0.0)
    rng=np.random.default_rng(seed); x=rng.standard_normal(Mc)
    for _ in range(iters):
        Xc=W*x[None,:]; den=(Xc*Xc).sum(1); u=np.where(den>1e-12,(Z*Xc).sum(1)/np.maximum(den,1e-12),0.0)
        Uc=W*u[:,None]; den2=(Uc*Uc).sum(0); x=np.where(den2>1e-12,(Z*Uc).sum(0)/np.maximum(den2,1e-12),0.0)
        n=np.linalg.norm(x)
        if n>0: x=x/n
    Xc=W*x[None,:]; den=(Xc*Xc).sum(1)
    u=np.where(den>1e-12,(Z*Xc).sum(1)/np.maximum(den,1e-12),np.nan)
    u[~(KEEP0&(den>1e-12))]=np.nan
    if np.corrcoef(x,rar0)[0,1]<0: u,x=-u,-x          # 符号锚到稀有度正向
    return u,x
U,xb=als_ux(V0)
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
mm=np.isfinite(U)&np.isfinite(RHO)&KEEP0
print(f"u_i 有效 {int(np.isfinite(U).sum()):,};corr(u_i, rho_i) = **{np.corrcoef(U[mm],RHO[mm])[0,1]:+.4f}**")
print(f"corr(ALS 的 x, 稀有度) = {np.corrcoef(xb,rar0)[0,1]:+.4f}")

lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EXTRA={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EXTRA.items()]
bi=np.flatnonzero(mm)
def cr(y,x,ii):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    return float(np.corrcoef(y[jj],x[jj])[0,1]), len(jj)
rng=np.random.default_rng(20260803)
rows=[]; nU=[]; nR=[]
for nm,y in OUT:
    ru,n=cr(y,U,bi); rr_,_=cr(y,RHO,bi)
    for store,x in ((nU,U),(nR,RHO)):
        ps=[]
        for _ in range(30):
            yp=y.copy(); yp[bi]=rng.permutation(y[bi]); v,_=cr(yp,x,bi)
            if np.isfinite(v): ps.append(abs(v))
        if len(ps)>=15: store.append(ps)
    rows.append(dict(q=nm[:56],n=n,r_u=ru,r_rho=rr_))
T=pd.DataFrame(rows); check_columns(T,'R268'); check_coverage(len(T),len(OUT),'R268 面板',tol=0.05)
th=lambda nl:(lambda L: float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nl]),axis=0),0.95)))(min(len(x) for x in nl))
tu,tr=th(nU),th(nR)
T['u_pass']=T.r_u.abs()>tu; T['rho_pass']=T.r_rho.abs()>tr
T.to_csv(pathlib.Path(__file__).parent/'results'/'panel.csv',index=False)
print(f"\n全族阈值:u_i {tu:.4f} · rho_i {tr:.4f}")
print(f"越阈值:u_i **{int(T.u_pass.sum())}/{len(T)}** · rho_i {int(T.rho_pass.sum())}/{len(T)}"
      f" · **交集 {int((T.u_pass&T.rho_pass).sum())}**")
print(f"\nu_i 最强的 6 个:")
for _,r in T.sort_values('r_u',key=abs,ascending=False).head(6).iterrows():
    print(f"  {r.r_u:>+8.4f}{'★' if r.u_pass else ' '}  (rho {r.r_rho:>+7.4f}{'★' if r.rho_pass else ' '})  {r.q[:50]}")

synth=np.full(len(U),np.nan); synth[bi]=U[bi]+rng.standard_normal(len(bi))*np.nanstd(U[bi])*2
r_pl,_=cr(synth,U,bi)
g=Gate('起始年龄矩阵里有几个人层维度')
g.asserted('正对照:由 u_i 造的合成结局必须强测到',abs(r_pl)>0.3,f"{r_pl:+.4f}")
g.asserted('⚠ u_i 与 rho_i 同源,不可能真正独立 —— 判的是外部结局重不重叠',True,
           f"corr(u_i, rho_i) = {np.corrcoef(U[mm],RHO[mm])[0,1]:+.4f}")
g.asserted('可判前提:rho_i 侧复现 `#199` 的 7/20 量级',abs(int(T[:len(lik)].rho_pass.sum())-7)<=2,
           f"Likert 段 rho_i 越阈值 {int(T[:len(lik)].rho_pass.sum())}/{len(lik)}")
g.asserted('注册的 kill:u_i 与 rho_i 越阈值的题交集 ≤2 -> 两个人层维度',
           int((T.u_pass&T.rho_pass).sum())<=2,
           f"交集 {int((T.u_pass&T.rho_pass).sum())};u_i {int(T.u_pass.sum())} · rho_i {int(T.rho_pass.sum())}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
