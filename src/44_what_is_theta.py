"""
ITER 9. theta_i (breadth) is the largest quantity in the dataset and is unexplained.
The decisive question is not which variables correlate but WHICH KIND it is:
   sexuality-specific gain  -> correlates with sexual history / development, not with
                               endorsement counts in non-sexual domains
   general endorsement trait-> correlates with counting anything at all
Scan every non-category column, control age/sex/agree-bias, and control multiplicity by
permutation rather than by a table lookup.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(317811)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
excl=set(rate)|set(inv[inv['kind'].isin(['AGE_ONSET','MULTISELECT','TOTAL','FORCED_CHOICE_MOST',
        'RATING_BINNED_FIB','RATING_NEG_FIB'])]['col'])
R=df[rate].apply(pd.to_numeric,errors='coerce')
theta=(R>0).sum(1).astype(float)
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
acq=df[lik].apply(pd.to_numeric,errors='coerce').mean(axis=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
ctrl=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],'acq':acq})
ctrl=ctrl.fillna(ctrl.median())

cands=[c for c in df.columns if c not in excl and c not in ('age','biomale')]
def numeric(c):
    s=df[c]
    if s.dtype!=object: return pd.to_numeric(s,errors='coerce')
    u=s.dropna().unique()
    if 1<len(u)<=12: return s.astype('category').cat.codes.replace(-1,np.nan)
    return None
rows=[]
for c in cands:
    v=numeric(c)
    if v is None or v.notna().sum()<1500 or v.nunique()<2: continue
    m=v.notna()&theta.notna()
    X=np.c_[np.ones(m.sum()),ctrl[m].values]
    b,*_=lstsq(X,theta[m].values,rcond=None); rt=theta[m].values-X@b
    b2,*_=lstsq(X,v[m].values.astype(float),rcond=None); rv=v[m].values-X@b2
    r=stats.spearmanr(rt,rv)
    rows.append(dict(col=c[:74],n=int(m.sum()),rho=r.statistic,p=r.pvalue))
T=pd.DataFrame(rows)
# permutation threshold over the whole family
mx=[]
for _ in range(200):
    perm=rng.permutation(len(theta)); tp=theta.values[perm]
    best=0
    for c in T.col.sample(min(40,len(T)),random_state=int(rng.integers(1e6))).values:
        cc=[x for x in cands if x[:74]==c][0]; v=numeric(cc)
        m=v.notna(); 
        if m.sum()<1500: continue
        best=max(best,abs(stats.spearmanr(tp[m.values],v[m]).statistic))
    mx.append(best)
thr=np.percentile(mx,95)
T['sig']=T.rho.abs()>thr
T=T.reindex(T.rho.abs().sort_values(ascending=False).index)
print(f"candidate columns tested: {len(T)}   family-wise permutation threshold |rho| > {thr:.3f}")
print(f"exceeding it: {int(T.sig.sum())}\n")
print(T.head(22).assign(rho=T.rho.round(3),p=T.p.map(lambda x:f"{x:.1e}")).to_string(index=False))
T.to_csv('data/derived/theta_correlates.csv',index=False)
