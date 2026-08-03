import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
ITER 8. Breadth drove the iter-5 congruence artifact, is what pornhabit split on, and is what
the induction report tracks. It has never been modelled. Is it one thing?

Q1  Given N endorsed categories, is a person's set MORE or LESS coordinate-diverse than a
    size-matched random set drawn from the same base rates?
      spreading along coordinates -> concentrated, LESS diverse than random
      occupying more coordinates  -> dispersed,    MORE diverse than random
      saying yes to everything    -> indistinguishable from base-rate random
Q2  Is breadth one trait? split-half over categories.
Q3  ATTACK: how much of breadth is acquiescence -- a tendency to endorse regardless of content?
    Measured on the survey's OWN mixed-valence Likert items, which are not erotic-category items.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import svd
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(75025)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float); H=H.where(R.notna())
have=H.notna().sum(1)>=40
H=H[have]; print(f"people with >=40 categories answered: {len(H):,}  categories: {H.shape[1]}")
Hf=H.fillna(H.mean())
base=Hf.mean(0).values
Z=Hf.values-base
U,S,Vt=svd(Z-Z.mean(0),full_matrices=False)
K=6; L=Vt[:K].T                       # category x K coordinate loadings
L=L/ (np.linalg.norm(L,axis=1,keepdims=True)+1e-9)
breadth=Hf.values.sum(1)

def dispersion(setmask):
    """participation ratio of the endorsed set's spread in coordinate space"""
    out=np.empty(len(setmask))
    for i,m in enumerate(setmask):
        idx=np.flatnonzero(m)
        if len(idx)<3: out[i]=np.nan; continue
        Y=L[idx]-L[idx].mean(0)
        ev=np.linalg.eigvalsh(Y.T@Y)
        out[i]=(ev.sum()**2)/max((ev**2).sum(),1e-12)
    return out
obs=dispersion(Hf.values>0.5)
# size-matched base-rate null: same N, drawn with probability proportional to base rate
nul=np.empty(len(obs))
p=base/base.sum()
for i,n_ in enumerate(breadth.astype(int)):
    if n_<3: nul[i]=np.nan; continue
    idx=rng.choice(len(base),size=min(n_,len(base)),replace=False,p=p)
    m=np.zeros(len(base),bool); m[idx]=True
    Y=L[m]-L[m].mean(0); ev=np.linalg.eigvalsh(Y.T@Y)
    nul[i]=(ev.sum()**2)/max((ev**2).sum(),1e-12)
ok=~(np.isnan(obs)|np.isnan(nul))
print(f"\nQ1 coordinate diversity of a person's own set vs a size-matched base-rate set")
print(f"   observed participation ratio : {np.nanmean(obs[ok]):.3f}")
print(f"   size-matched null            : {np.nanmean(nul[ok]):.3f}")
print(f"   difference                   : {np.nanmean(obs[ok]-nul[ok]):+.4f}  "
      f"(paired t = {stats.ttest_rel(obs[ok],nul[ok]).statistic:+.1f})")
d=obs-nul
print(f"   corr(breadth, own-minus-null diversity) = {stats.spearmanr(breadth[ok],d[ok]).statistic:+.3f}")
print("   -> " + ("sets are MORE diverse than chance: breadth occupies more coordinates"
      if np.nanmean(d[ok])>0 else "sets are LESS diverse than chance: breadth goes further along fewer"))

cats=np.arange(H.shape[1]); rng.shuffle(cats); a,b=cats[:len(cats)//2],cats[len(cats)//2:]
ba,bb=Hf.values[:,a].sum(1),Hf.values[:,b].sum(1)
r=np.corrcoef(ba,bb)[0,1]
print(f"\nQ2 is breadth one trait? split-half over categories r={r:+.3f}  Spearman-Brown={2*r/(1+r):+.3f}")

lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
Lk=df[lik].apply(pd.to_numeric,errors='coerce')[have]
acq=Lk.mean(axis=1)                     # signed mean over mixed-valence items = agree-bias
ext=Lk.abs().mean(axis=1)               # extreme responding
print(f"\nQ3 acquiescence attack ({len(lik)} non-category Likert items)")
for nm,v in [('agree-bias',acq),('extreme-responding',ext)]:
    m=v.notna()
    print(f"   corr(breadth, {nm:20s}) = {stats.spearmanr(breadth[m.values],v[m]).statistic:+.3f}")
m=acq.notna()&ext.notna()
X=np.c_[np.ones(m.sum()),acq[m].values,ext[m].values]
bb_,*_=np.linalg.lstsq(X,breadth[m.values],rcond=None)[:1]+(None,None,None)
res=breadth[m.values]-X@bb_
print(f"   share of breadth variance explained by response style = "
      f"{1-res.var()/breadth[m.values].var():.3f}")
np.save('data/derived/breadth.npy',breadth)
