import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
If breadth is a THRESHOLD parameter rather than a structural one -- the readout's gain moves and
everything above threshold gets endorsed -- then sets should be approximately NESTED: a broad
person's set should largely CONTAIN a narrow person's. If breadth instead means different tastes
at larger volume, containment should be no better than size-matched chance.
Null: base-rate-random sets of the identical sizes, which already produce substantial containment
because popular categories are shared. The question is the EXCESS over that.
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(196418)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float).where(R.notna())
H=H[H.notna().sum(1)>=40]; V=(H.fillna(0).values>0.5)
base=V.mean(0); p=base/base.sum(); K=V.shape[1]
n=len(V); sz=V.sum(1)
print(f"people {n:,}  categories {K}  median set size {int(np.median(sz))}")
def containment(A,Bm):
    if Bm.sum()==0: return np.nan
    return (A&Bm).sum()/Bm.sum()
pairs=[]
for _ in range(40000):
    i,j=rng.integers(0,n,2)
    if sz[i]<=sz[j] or sz[j]<5: continue
    pairs.append((i,j))
pairs=pairs[:20000]
obs=np.array([containment(V[i],V[j]) for i,j in pairs])
def rnd(k_):
    m=np.zeros(K,bool); m[rng.choice(K,size=min(k_,K),replace=False,p=p)]=True; return m
nul=np.array([containment(rnd(int(sz[i])),rnd(int(sz[j]))) for i,j in pairs])
print(f"\npairs (broader, narrower), n={len(pairs):,}")
print(f"  observed containment of the narrow set by the broad set : {obs.mean():.4f}")
print(f"  size-matched base-rate null                             : {nul.mean():.4f}")
print(f"  excess                                                  : {obs.mean()-nul.mean():+.4f}")
print(f"  paired t                                                : {stats.ttest_rel(obs,nul).statistic:+.1f}")
print(f"  perfect nesting would be 1.000; observed reaches "
      f"{(obs.mean()-nul.mean())/(1-nul.mean()):.1%} of the gap from chance to perfect")
r=sz[[i for i,_ in pairs]]/np.maximum(sz[[j for _,j in pairs]],1)
q=pd.qcut(r,4,labels=False,duplicates='drop')
print("\n  by size ratio (broad/narrow):")
for k_ in np.unique(q):
    m=q==k_
    print(f"    ratio quartile {k_+1}: observed {obs[m].mean():.3f}  null {nul[m].mean():.3f}  excess {obs[m].mean()-nul[m].mean():+.4f}")
