"""
Q1 measured a person's set as 1.1% MORE concentrated than a size-matched base-rate set. But the
coordinate loadings were fitted on the same people, which guarantees real sets align with them
better than random sets do. Refit the loadings OUT OF SAMPLE and use many null draws per person.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import svd
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(121393)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float).where(R.notna())
H=H[H.notna().sum(1)>=40]; Hf=H.fillna(H.mean()); V=Hf.values
n=len(V); idx=rng.permutation(n); f,e=idx[:n//2],idx[n//2:]
base=V[f].mean(0)
Z=V[f]-V[f].mean(0)
L=svd(Z,full_matrices=False)[2][:6].T
L=L/(np.linalg.norm(L,axis=1,keepdims=True)+1e-9)
def PR(mask):
    i=np.flatnonzero(mask)
    if len(i)<3: return np.nan
    Y=L[i]-L[i].mean(0); ev=np.linalg.eigvalsh(Y.T@Y)
    return (ev.sum()**2)/max((ev**2).sum(),1e-12)
p=base/base.sum()
obs=[];nul=[];bn=[]
for r_ in V[e]:
    m=r_>0.5; k_=int(m.sum())
    if k_<3: continue
    o=PR(m)
    ns=[]
    for _ in range(10):
        j=rng.choice(len(base),size=min(k_,len(base)),replace=False,p=p)
        mm=np.zeros(len(base),bool); mm[j]=True; ns.append(PR(mm))
    obs.append(o); nul.append(np.nanmean(ns)); bn.append(k_)
obs=np.array(obs); nul=np.array(nul); bn=np.array(bn)
ok=~(np.isnan(obs)|np.isnan(nul))
d=obs[ok]-nul[ok]
print(f"out-of-sample loadings, 10 null draws per person, n={ok.sum():,}")
print(f"  observed participation ratio : {obs[ok].mean():.3f}")
print(f"  size-matched null            : {nul[ok].mean():.3f}")
print(f"  difference                   : {d.mean():+.4f}  ({100*d.mean()/nul[ok].mean():+.2f}%)")
print(f"  paired t                     : {stats.ttest_rel(obs[ok],nul[ok]).statistic:+.1f}")
print(f"  in-sample version gave       : -0.0520  (-1.13%)")
print(f"\n  corr(breadth, own-minus-null concentration) = {stats.spearmanr(bn[ok],d).statistic:+.3f}")
lo,hi=bn[ok]<np.median(bn[ok]), bn[ok]>=np.median(bn[ok])
print(f"  narrow people (n<{int(np.median(bn[ok]))}): diff {d[lo].mean():+.4f}   broad people: diff {d[hi].mean():+.4f}")
print("\n  VERDICT:", "breadth sets are essentially base-rate random in coordinate diversity"
      if abs(d.mean())/nul[ok].mean()<0.03 else "breadth has real coordinate structure")
