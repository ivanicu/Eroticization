import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The pornhabit coordinate deficit survives theta-matching at 78%. But men consume more and sex
reorganises the grammar at .084, so the residual could be sex leaking through group composition
even though demographics are partialled at the item level.
Match on sex AND block count AND theta simultaneously. Positive control stays: an unrelated
split matched the same way must not collapse, or the triple matching is just destroying signal.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(1346269)
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
neu=pd.to_numeric(df['neuroticismvariable'],errors='coerce').reindex(pool)
tq=pd.qcut(theta,6,labels=False,duplicates='drop')
def key(i,use):
    parts=[nblk.reindex(i).astype(int)]
    if 'theta' in use: parts.append(tq.reindex(i).fillna(-1).astype(int))
    if 'sex'  in use: parts.append(male.reindex(i).fillna(-1).astype(int))
    return pd.Series(list(zip(*[p.values for p in parts])),index=i)
def run(g1,g2,use,label):
    i1,i2=pd.Index(g1),pd.Index(g2); k1,k2=key(i1,use),key(i2,use)
    a=[];b=[]
    for v in set(k1)|set(k2):
        x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
        if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
    a,b=np.array(a),np.array(b)
    if len(a)<500: print(f"   {label}: too few ({len(a)})"); return
    c=cong(loadings(a),loadings(b))
    ceil=[cong(loadings(p[:len(a)]),loadings(p[len(a):len(a)+len(b)])) for p in [rng.permutation(pool) for _ in range(6)]]
    d=float(np.nanmean(ceil))-c
    print(f"   {label:38s} n={len(a):,}  deficit={d:.4f} +/- {np.nanstd(ceil):.4f}"
          f"  |male gap|={abs(male.reindex(a).mean()-male.reindex(b).mean()):.3f}"
          f"  |theta gap|={abs(theta.reindex(a).mean()-theta.reindex(b).mean()):.4f}")
    return d
G=(pool[(ph>ph.median()).values], pool[(ph<ph.median()).values])
print("=== PORNHABIT hi/lo, escalating match ===")
d1=run(*G,['blk'],'blocks only')
d2=run(*G,['blk','theta'],'+ theta')
d3=run(*G,['blk','theta','sex'],'+ theta + sex')
print("\n=== POSITIVE CONTROL: NEUROTICISM hi/lo through the identical triple match ===")
N=(pool[(neu>neu.median()).values], pool[(neu<neu.median()).values])
run(*N,['blk','theta','sex'],'+ theta + sex')
print("\n=== POSITIVE CONTROL 2: SEX through blocks+theta (must not collapse) ===")
S=(pool[(male==1).values], pool[(male==0).values])
run(*S,['blk','theta'],'+ theta')
if d1 and d3: print(f"\n  pornhabit deficit surviving the full match: {100*d3/d1:.0f}%")
