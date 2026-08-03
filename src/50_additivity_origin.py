"""
The self/other contrast correlates +.75 between two male-origin substances and -.24 between a
male- and a female-origin one. Two readings:
  (A) genuine non-additivity -- 'receiving' has no substance-independent weight, Ivan's A_i term
  (B) a MISSING FEATURE -- partner gender. Someone attracted to men wants to receive male fluids.
      Under a basis that crosses role with source gender, additivity could be restored.
Separator: partial the orientation proxies out of every difference vector and recompute.
  (B) -> the male/female split collapses
  (A) -> it persists
This also settles a deeper point: additivity is BASIS-DEPENDENT, so "is phi additive" is only
answerable relative to a stated feature set.
"""
import pandas as pd, numpy as np, warnings, itertools
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
exec(open('src/49_additivity.py').read().split('print(f"\\ndifference vectors built')[0])
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
ORI=[c for c in df.columns if 'opposite gender to me' in c or 'gender identity' in c.lower()]
Q=df[['biomale']+ORI].copy()
for c in Q.columns:
    if Q[c].dtype==object: Q[c]=Q[c].astype('category').cat.codes.replace(-1,np.nan)
Q=Q.apply(pd.to_numeric,errors='coerce'); Q=Q.fillna(Q.median())
MALE_SRC={'precum','ejaculate'}; FEM_SRC={'squirt','breastmilk'}; NEU={'saliva','urine','sweat'}
Q=Q.replace([np.inf,-np.inf],np.nan).fillna(0.0)
def resid(v):
    v=v[np.isfinite(v.values)]
    idx=v.index.intersection(Q.index)
    v=v.reindex(idx)
    X=np.c_[np.ones(len(idx)),Q.loc[idx].values]
    good=np.isfinite(X).all(1)&np.isfinite(v.values)
    idx=idx[good]; X=X[good]; y=v.values[good]
    if len(idx)<50 or not np.isfinite(X).all(): return v
    b,*_=lstsq(X,y,rcond=None); return pd.Series(y-X@b,index=idx)
Dr={k:resid(v) for k,v in D.items()}
def pairs(Dx):
    same_within=[];same_across=[];diff=[]
    for k1,k2 in itertools.combinations(Dx,2):
        s1,a1=k1; s2,a2=k2
        if s1==s2: continue
        a,b=Dx[k1],Dx[k2]; c=a.index.intersection(b.index)
        if len(c)<250: continue
        x,y=a.reindex(c).values,b.reindex(c).values
        if x.std()==0 or y.std()==0: continue
        r=np.corrcoef(x,y)[0,1]
        if a1!=a2: diff.append(r); continue
        gendered = (s1 in MALE_SRC or s1 in FEM_SRC) and (s2 in MALE_SRC or s2 in FEM_SRC)
        if gendered:
            (same_within if ((s1 in MALE_SRC)==(s2 in MALE_SRC)) else same_across).append(r)
    return np.array(same_within),np.array(same_across),np.array(diff)
for label,Dx in [('RAW',D),('ORIENTATION PARTIALLED OUT',Dr)]:
    w,a,d=pairs(Dx)
    print(f"\n=== {label} ===")
    print(f"  same contrast, SAME-gender source pair   n={len(w):2d}  mean r = {w.mean():+.3f}")
    print(f"  same contrast, CROSS-gender source pair  n={len(a):2d}  mean r = {a.mean():+.3f}")
    print(f"  different contrast                       n={len(d):2d}  mean r = {d.mean():+.3f}")
    print(f"  same-gender minus cross-gender           = {w.mean()-a.mean():+.3f}"
          f"   Welch t={stats.ttest_ind(w,a,equal_var=False).statistic:+.2f}"
          f"  p={stats.ttest_ind(w,a,equal_var=False).pvalue:.3g}")
print("\n  (B) missing-feature predicts the gap collapses once orientation is removed")
print("  (A) non-additivity  predicts it persists")
