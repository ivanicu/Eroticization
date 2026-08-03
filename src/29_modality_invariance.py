"""
ITER 5. Ivan's phase-4 cross-modal test, already run, sitting in the data.
If "sexual" is a VISUAL CONTENT category, the grammar should differ between people who consume
erotic content as pictures and as text -- and between live-action and drawn, where drawn content
contains no real bodies at all. If it is an abstract valuation on ordinary semantics, the
coordinates should be modality-invariant.

People differ between groups, so compare in the space they SHARE: the items. Fit item loadings
within each group, orthonormalise, and take the singular values of V1 V2^T = subspace congruence.
Ceiling = random splits of the same sizes. POSITIVE CONTROL = split by sex, a difference the
method must be able to see, otherwise a modality null is just blindness.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(8191)
exec(open('src/16_dimensionality.py').read().split("allq=list(B)")[0])
allq=list(B)
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
ANI=[c for c in df.columns if 'type of erotic content you prefer tends to be more' in c][0]
print("modality  :",df[MOD].value_counts().to_dict())
print("animation :",df[ANI].value_counts().to_dict())

def groups(col,a,b):
    s=df[col]
    return set(df.index[s.isin(a)]), set(df.index[s.isin(b)])
G_written = groups(MOD,['Mostly written','Entirely written'],['Mostly visual','Entirely visual'])
G_drawn   = groups(ANI,['Mostly animated/drawn','Entirely animated/drawn'],
                       ['Mostly live action vid/photos','Entirely live action vid/photos'])
G_sex     = (set(df.index[df.biomale==1]), set(df.index[df.biomale==0]))

pool_all=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
ITEMS=[(q,j) for q in allq for j in range(B[q]['R'].shape[1])]
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool_all)))
    if len(ppl)<600: return None
    pm={p:i for i,p in enumerate(ppl)}
    cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0) if np.isfinite(np.nanmean(Z,axis=0)).all() else 0.0,Z)
        cols.append(np.nan_to_num(Z))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    U,S,Vt=svd(Z,full_matrices=False)
    return Vt[:K]
def congruence(A_,B_):
    if A_ is None or B_ is None: return np.nan
    Qa,_=qr(A_.T,mode='reduced'); Qb,_=qr(B_.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))

res={}
for name,(g1,g2) in [('MODALITY written vs visual',G_written),
                     ('ANIMATION drawn vs live',G_drawn),
                     ('SEX male vs female  [pos. control]',G_sex)]:
    c=congruence(loadings(g1),loadings(g2))
    n1,n2=len(g1&set(pool_all)),len(g2&set(pool_all))
    # matched random-split ceiling at the SAME group sizes
    ceil=[]
    for _ in range(8):
        p=rng.permutation(pool_all); ceil.append(congruence(loadings(p[:n1]),loadings(p[n1:n1+n2])))
    res[name]=(c,np.nanmean(ceil),np.nanstd(ceil),n1,n2)
    print(f"\n{name}\n   n = {n1:,} vs {n2:,}"
          f"\n   subspace congruence           = {c:.3f}"
          f"\n   matched random-split ceiling  = {np.nanmean(ceil):.3f} +/- {np.nanstd(ceil):.3f}"
          f"\n   deficit vs ceiling            = {np.nanmean(ceil)-c:+.3f}"
          f"  ({(np.nanmean(ceil)-c)/max(np.nanstd(ceil),1e-9):+.1f} sd)")
