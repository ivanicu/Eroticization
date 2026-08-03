import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""Block-count-matched deficits for the two modality splits, the ones the A-vs-B claim rests on."""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(112358)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool_all=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool_all).fillna(0)
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool_all)))
    if len(ppl)<600: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.0)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
ANI=[c for c in df.columns if 'type of erotic content you prefer tends to be more' in c][0]
S={'MODALITY written/visual':(set(df.index[df[MOD].isin(['Mostly written','Entirely written'])]),
                              set(df.index[df[MOD].isin(['Mostly visual','Entirely visual'])])),
   'ANIMATION drawn/live':(set(df.index[df[ANI].isin(['Mostly animated/drawn','Entirely animated/drawn'])]),
                           set(df.index[df[ANI].isin(['Mostly live action vid/photos','Entirely live action vid/photos'])]))}
prev={'MODALITY written/visual':0.0360,'ANIMATION drawn/live':0.0123}
for name,(g1,g2) in S.items():
    g1=np.array(sorted(set(g1)&set(pool_all))); g2=np.array(sorted(set(g2)&set(pool_all)))
    s1,s2=nblk.reindex(g1),nblk.reindex(g2); k1=[];k2=[]
    for v in sorted(set(s1.unique())|set(s2.unique())):
        a=g1[s1.values==v]; b=g2[s2.values==v]; m=min(len(a),len(b))
        if m: k1+=list(rng.choice(a,m,replace=False)); k2+=list(rng.choice(b,m,replace=False))
    k1,k2=np.array(k1),np.array(k2)
    c=cong(loadings(k1),loadings(k2))
    ceil=[cong(loadings(p[:len(k1)]),loadings(p[len(k1):len(k1)+len(k2)])) for p in [rng.permutation(pool_all) for _ in range(6)]]
    d=float(np.nanmean(ceil))-c; sd=float(np.nanstd(ceil))
    print(f"{name}: n {len(k1):,} vs {len(k2):,}  mean blocks {nblk.reindex(k1).mean():.1f}/{nblk.reindex(k2).mean():.1f}")
    print(f"   matched deficit {d:.4f} +/- {sd:.4f} ({d/max(sd,1e-9):.1f} sd)   unmatched {prev[name]:.4f}"
          f"   survives {100*d/prev[name]:.0f}%\n")
