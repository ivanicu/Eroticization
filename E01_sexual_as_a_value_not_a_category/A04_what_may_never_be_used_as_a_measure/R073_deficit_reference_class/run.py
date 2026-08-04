import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
A deficit of 0.034 is only 'small' against a reference class. Locate modality among OTHER
person-variable splits of comparable kind, all measured the same way against their own matched
random-split ceiling. If modality sits at the bottom of the pack it is invariant; if it sits in
the middle, 0.034 is just what any real person variable costs.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(65537)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
pool_all=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
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
def cong(A_,B_):
    if A_ is None or B_ is None: return np.nan
    Qa,_=qr(A_.T,mode='reduced'); Qb,_=qr(B_.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
def median_split(col):
    s=pd.to_numeric(df[col],errors='coerce'); m=s.median()
    return set(df.index[s>m]), set(df.index[s<m])
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
ANI=[c for c in df.columns if 'type of erotic content you prefer tends to be more' in c][0]
REL=[c for c in df.columns if 'preferred relationship style' in c][0]
SPL={
 'MODALITY written/visual':(set(df.index[df[MOD].isin(['Mostly written','Entirely written'])]),
                            set(df.index[df[MOD].isin(['Mostly visual','Entirely visual'])])),
 'ANIMATION drawn/live'  :(set(df.index[df[ANI].isin(['Mostly animated/drawn','Entirely animated/drawn'])]),
                           set(df.index[df[ANI].isin(['Mostly live action vid/photos','Entirely live action vid/photos'])])),
 'SEX male/female'       :(set(df.index[df.biomale==1]), set(df.index[df.biomale==0])),
 'OPENNESS hi/lo'        :median_split('opennessvariable'),
 'NEUROTICISM hi/lo'     :median_split('neuroticismvariable'),
 'EXTROVERSION hi/lo'    :median_split('extroversionvariable'),
 'POWERLESSNESS hi/lo'   :median_split('powerlessnessvariable'),
 'PORNHABIT hi/lo'       :median_split('pornhabit'),
 'MONOGAMY yes/no'       :(set(df.index[df[REL]=='Monogamous']), set(df.index[df[REL]=='Not monogamous'])),
}
rows=[]
for name,(g1,g2) in SPL.items():
    n1,n2=len(g1&set(pool_all)),len(g2&set(pool_all))
    if min(n1,n2)<600: continue
    c=cong(loadings(g1),loadings(g2))
    ceil=[cong(loadings(p[:n1]),loadings(p[n1:n1+n2])) for p in [rng.permutation(pool_all) for _ in range(6)]]
    rows.append(dict(split=name,n1=n1,n2=n2,cong=round(c,3),ceiling=round(float(np.nanmean(ceil)),3),
                     deficit=round(float(np.nanmean(ceil))-c,4), sd=round(float(np.nanstd(ceil)),4)))
T=pd.DataFrame(rows).sort_values('deficit')
T['z']=(T.deficit/T.sd).round(1)
print(T.to_string(index=False))
T.to_csv('data/derived/deficit_reference.csv',index=False)
mod=T[T.split.str.startswith('MODALITY')].deficit.iloc[0]; sex=T[T.split.str.startswith('SEX')].deficit.iloc[0]
print(f"\n  modality deficit / sex deficit = {mod/sex:.2f}")
print(f"  modality rank among {len(T)} splits (1 = most invariant) : {int((T.deficit<=mod).sum())}")
print(f"  median deficit across all splits = {T.deficit.median():.4f}")
