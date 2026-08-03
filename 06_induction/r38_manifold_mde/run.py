import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Both misfit and extremity came out null at rho~0.02. The positive control spans only .943->.978
for a full 1-sd of injected noise, so the measure is compressed and the null may be low power.
MDE: inject misfit ONLY into the 'totally different' group -- the exact alternative hypothesis --
and read off the rho that would have been observed. Same for extremity.
Also coarsen the block-count matching into strata to recover sample size.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(514229)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
pm={p:i for i,p in enumerate(pool)}; cols=[]
for q in allq:
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
    mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
    cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
X=np.hstack(cols); X=X-X.mean(0)
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
b,*_=lstsq(D,X,rcond=None); X=X-D@b
K=5; half=rng.permutation(len(pool)); f,e=half[:len(half)//2],half[len(half)//2:]
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
gp=df[IND].map(ORD).reindex(pool[e]).values
nb=nblk.reindex(pool[e]).values
strata=pd.cut(nb,[0,4,6,8,12,40],labels=False)
ok=~np.isnan(gp)
keep=[]
for v in np.unique(strata[ok]):
    idx=[np.flatnonzero(ok&(strata==v)&(gp==g_)) for g_ in range(4)]
    m=min(len(i) for i in idx)
    if m: keep+=[x for i in idx for x in rng.choice(i,m,replace=False)]
keep=np.array(keep)
print(f"stratum-matched n = {len(keep):,}  (exact-match was 3,352)")
print("  mean blocks per group:", [round(float(nb[keep[gp[keep]==g_]].mean()),2) for g_ in range(4)])

def measure(Xm):
    U,S,Vt=svd(Xm[f]-Xm[f].mean(0),full_matrices=False); V=Vt[:K]
    E=Xm[e]-Xm[f].mean(0); on=E@V.T
    return np.linalg.norm(E-on@V,axis=1)/np.linalg.norm(E,axis=1), np.linalg.norm(on,axis=1)
r0,x0=measure(X)
print(f"\nOBSERVED  misfit rho = {stats.spearmanr(gp[keep],r0[keep]).statistic:+.4f}"
      f"   extremity rho = {stats.spearmanr(gp[keep],x0[keep]).statistic:+.4f}")
print("\nMDE -- inject the alternative into the 'totally different' group only:")
tgt=e[gp==3]
for lvl in [0.1,0.2,0.35,0.5,0.75]:
    Xn=X.copy(); Xn[tgt]=Xn[tgt]+rng.normal(0,lvl*X.std(),size=Xn[tgt].shape)
    rr,_=measure(Xn)
    print(f"   added misfit {lvl:.2f}x sd  -> misfit rho = {stats.spearmanr(gp[keep],rr[keep]).statistic:+.4f}")
print("\nMDE -- inject EXTREMITY (scale their on-manifold position) into that group only:")
U,S,Vt=svd(X[f]-X[f].mean(0),full_matrices=False); V=Vt[:K]
for mult in [1.1,1.25,1.5,2.0]:
    Xn=X.copy()
    proj=(Xn[tgt]-X[f].mean(0))@V.T
    Xn[tgt]=Xn[tgt]+ (mult-1)*(proj@V)
    _,xx=measure(Xn)
    print(f"   extremity x{mult:.2f}      -> extremity rho = {stats.spearmanr(gp[keep],xx[keep]).statistic:+.4f}")
