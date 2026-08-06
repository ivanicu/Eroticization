import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A47 R252 -- 「第三个维度」有多少是年龄

`#206b`:`rho_i × age = +0.1563` 是整张非情色面板最大的一格,比 rho_i 做的任何别的事大 4 倍。
**而 `#199` 那条「第三个维度」(rho_i 越阈值 7/20)的整条面板从没控过年龄。**

ESTIMAND        把 `#244` 的 20 道 Likert 面板**原样重跑**,**唯一改动是把 age 加进控制项**;
                并报 `rho_i` 在**每个年龄段内部**的越阈值数。
KILL            **若控制年龄后越阈值题数从 7 掉到 ≤2 -> 「第三个维度」有一大半是年龄,
                `#199`/`#200`/`#245` 里凡以 rho_i 为预测子的结论都要降级。**
NEGATIVE CTRL   每题在分析样本内打乱(`#184b`)。
POSITIVE CTRL   把 age 自己当结局塞进面板 -> 控制 age 后它必须塌到零(证明控制生效)。
段内对照        5 个年龄段各自跑一遍 —— 段内没有年龄变异,所以段内仍显著的题
                是**不可能由年龄解释**的那些。
IMPOSSIBLE      年龄只有 5 段(14–17 … 29–32),**段内仍有 3–4 年的残余变异**;
                所以"控制年龄"是不完全的,段内对照是对它的补充而不是替代。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
_,RHO=betas(V)
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=d['age'].map(AGE).values.astype(float)
m0=np.isfinite(RHO)&KEEP&np.isfinite(age); bi=np.flatnonzero(m0)
print(f"n = {len(bi):,};corr(rho_i, age) = {np.corrcoef(RHO[bi],age[bi])[0,1]:+.4f}(`#206b` 报 +0.1563)")

def cr(y,x,ii,ctrls=()):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    for c in ctrls: m2=np.isfinite(c[jj]); jj=jj[m2]
    X=np.c_[np.ones(len(jj)),*[c[jj] for c in ctrls]] if ctrls else np.ones((len(jj),1))
    ry=y[jj]-X@np.linalg.lstsq(X,y[jj],rcond=None)[0]
    rx=x[jj]-X@np.linalg.lstsq(X,x[jj],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1]), len(jj)

rng=np.random.default_rng(20260803)
OUT=[(c,d[c].values.astype(float)) for c in lik]+[('【正对照】age 自己',age)]
rows=[]; nul0=[]; nul1=[]
for nm,y in OUT:
    r0,n=cr(y,RHO,bi); r1,_=cr(y,RHO,bi,(age,))
    for store,ct in ((nul0,()),(nul1,(age,))):
        ps=[]
        for _ in range(40):
            yp=y.copy(); yp[bi]=rng.permutation(y[bi]); v,_=cr(yp,RHO,bi,ct)
            if np.isfinite(v): ps.append(abs(v))
        if len(ps)>=20 and '正对照' not in nm: store.append(ps)
    rows.append(dict(q=nm[:58],n=n,r_raw=r0,r_ctrl_age=r1))
T=pd.DataFrame(rows); check_columns(T,'R252'); check_coverage(len(T),len(OUT),'R252 面板',tol=0.0)
th=lambda nl:(lambda L: float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nl]),axis=0),0.95)))(min(len(x) for x in nl))
t0,t1=th(nul0),th(nul1)
T=T.sort_values('r_raw',key=abs,ascending=False)
T.to_csv(pathlib.Path(__file__).parent/'results'/'age_controlled.csv',index=False)
print(f"\n全族阈值:未控 {t0:.4f} · 控年龄后 {t1:.4f}\n")
print(f"{'未控':>9}{'控年龄后':>11}  题")
for _,r in T.head(10).iterrows():
    print(f"{r.r_raw:>+9.4f}{'★' if abs(r.r_raw)>t0 else ' '}{r.r_ctrl_age:>+10.4f}"
          f"{'★' if abs(r.r_ctrl_age)>t1 else ' '}  {r.q[:56]}")
real=T[~T.q.str.contains('正对照')]
n0=int((real.r_raw.abs()>t0).sum()); n1=int((real.r_ctrl_age.abs()>t1).sum())
print(f"\n越阈值:未控 **{n0}/{len(real)}** -> 控年龄后 **{n1}/{len(real)}**")

# 段内对照
print("\n---- 段内(每段内部没有年龄变异)----")
seg=[]
for lab,v in AGE.items():
    ii=bi[age[bi]==v]
    if len(ii)<800: continue
    k=0
    for c in lik:
        y=d[c].values.astype(float); r,_=cr(y,RHO,ii)
        if abs(r)>t0: k+=1
    seg.append(dict(band=lab,n=len(ii),n_pass=k))
    print(f"  {lab}  n={len(ii):>5,}  越(未控)阈值 {k}/20")
S=pd.DataFrame(seg); S.to_csv(pathlib.Path(__file__).parent/'results'/'by_band.csv',index=False)

pc=T[T.q.str.contains('正对照')].iloc[0]
g=Gate('第三个维度有多少是年龄')
g.asserted('正对照:age 自己在控制 age 后必须塌到零',abs(pc.r_ctrl_age)<0.02,
           f"未控 {pc.r_raw:+.4f} -> 控后 {pc.r_ctrl_age:+.6f}")
g.asserted('可判前提:未控时复现 `#199` 的 7/20',abs(n0-7)<=1,f"{n0}/20")
g.asserted('注册的 kill:控年龄后越阈值 ≤2 -> 第三个维度有一大半是年龄',n1<=2,
           f"{n0} -> {n1}")
g.asserted('⚠ 年龄只有 5 段,段内仍有 3–4 年残余 -> 控制不完全,段内对照是补充不是替代',True,
           f"段内越阈值:{[f'{r.band}:{r.n_pass}' for r in S.itertuples()]}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
