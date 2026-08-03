import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
ITER 10. Iter 5: pornhabit reorganises the grammar (block-matched congruence deficit .0871,
about sex-sized). Iter 9: pornhabit is the only thing touching theta (rho .17). Never modelled
together. Is consumption linked to HOW MANY (theta) or to WHICH ONES (coordinates)?

Separator: match groups on THETA as well as block count, and re-measure congruence.
  consumption -> more of everything : deficit collapses once theta is matched
  consumption -> different things   : deficit survives theta matching
POSITIVE CONTROL for the matching itself: run SEX through the identical theta-matching. If sex
also collapses, the procedure is over-matching and proves nothing.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(832040)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
Rt=df[rate].apply(pd.to_numeric,errors='coerce')
theta=((Rt>0).sum(1)/Rt.notna().sum(1).clip(lower=1)).reindex(pool)
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool)))
    if len(ppl)<500: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
male=df['biomale'].reindex(pool)
SPL={'PORNHABIT hi/lo':(pool[(ph>ph.median()).values], pool[(ph<ph.median()).values]),
     'SEX male/female [pos. control]':(pool[(male==1).values], pool[(male==0).values])}
tq=pd.qcut(theta,8,labels=False,duplicates='drop')
def matched(g1,g2,also_theta):
    i1=pd.Index(g1); i2=pd.Index(g2)
    s1=nblk.reindex(i1).astype(int); s2=nblk.reindex(i2).astype(int)
    if also_theta:
        k1=pd.Series(list(zip(s1,tq.reindex(i1).fillna(-1))),index=i1)
        k2=pd.Series(list(zip(s2,tq.reindex(i2).fillna(-1))),index=i2)
    else:
        k1,k2=pd.Series(list(s1),index=i1),pd.Series(list(s2),index=i2)
    a=[];b=[]
    for v in set(k1)|set(k2):
        x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
        if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
    return np.array(a),np.array(b)
for name,(g1,g2) in SPL.items():
    print(f"\n=== {name} ===")
    for lab,alsoT in [('block-count matched only',False),('block-count AND theta matched',True)]:
        a,b=matched(g1,g2,alsoT)
        if len(a)<500: print(f"   {lab}: too few after matching ({len(a)})"); continue
        c=cong(loadings(a),loadings(b))
        ceil=[cong(loadings(p[:len(a)]),loadings(p[len(a):len(a)+len(b)])) for p in [rng.permutation(pool) for _ in range(6)]]
        d=float(np.nanmean(ceil))-c
        gapT=abs(theta.reindex(a).mean()-theta.reindex(b).mean())
        print(f"   {lab:30s} n={len(a):,}/{len(b):,}  deficit={d:.4f} +/- {np.nanstd(ceil):.4f}"
              f"  |theta gap|={gapT:.4f}")
