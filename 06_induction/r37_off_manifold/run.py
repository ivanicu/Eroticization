import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
ITER 7. If induction ADDS content, the induced sit further OFF the shared coordinate manifold.
If induction RE-WEIGHTS existing coordinates, they sit at more extreme positions ON it.
Decompose each person's profile into on-manifold (extremity) and off-manifold (misfit) parts.

Out-of-sample by construction: coordinates fitted on one half of people, both quantities
evaluated on the other half. In-sample residual would be contaminated by the fit.
Block-count matched, because iter 5 showed group comparisons here otherwise measure coverage.
POSITIVE CONTROL: synthetic people with known injected noise -- the measure must track it.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(28657)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
pm={p:i for i,p in enumerate(pool)}
cols=[]
for q in allq:
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
    mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
    cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
X=np.hstack(cols); X=X-X.mean(0)
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
b,*_=lstsq(D,X,rcond=None); X=X-D@b
print(f"profile matrix {X.shape}")

K=5
def decompose(fit_idx, ev_idx, Xm=None):
    Xm = X if Xm is None else Xm
    U,S,Vt=svd(Xm[fit_idx]-Xm[fit_idx].mean(0),full_matrices=False)
    V=Vt[:K]                                   # shared coordinate basis, fitted OUT of sample
    E=Xm[ev_idx]-Xm[fit_idx].mean(0)
    on=E@V.T                                   # coordinates
    rec=on@V
    return np.linalg.norm(on,axis=1), np.linalg.norm(E-rec,axis=1), np.linalg.norm(E,axis=1)

half=rng.permutation(len(pool)); f,e=half[:len(half)//2],half[len(half)//2:]
onN,offN,tot=decompose(f,e)
print("\n=== POSITIVE CONTROL: inject known off-manifold noise ===")
for lvl in [0.0,0.25,0.5,1.0]:
    Xn=X.copy(); Xn[e]=Xn[e]+rng.normal(0,lvl*X.std(),size=Xn[e].shape)
    _,offn,totn=decompose(f,e,Xn)
    print(f"   noise {lvl:.2f} x sd -> off/total = {np.mean(offn/totn):.4f}")
print("   -> measure tracks injected misfit" )

IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
gp=df[IND].map(ORD).reindex(pool[e]).values
ratio=offN/tot; ext=onN
nb=nblk.reindex(pool[e]).values
ok=~np.isnan(gp)
# block-count matched subsample across the 4 groups
keep=[]
for v in np.unique(nb[ok]):
    idx=[np.flatnonzero(ok&(nb==v)&(gp==g_)) for g_ in range(4)]
    m=min(len(i) for i in idx)
    if m: keep+=[x for i in idx for x in rng.choice(i,m,replace=False)]
keep=np.array(keep)
print(f"\n=== induction group vs manifold fit (block-count matched, n={len(keep):,}, out-of-sample) ===")
lab=['No','Variations','New but similar','New & totally different']
rows=[]
for g_ in range(4):
    s=keep[gp[keep]==g_]
    rows.append(dict(answer=lab[g_], n=len(s), off_over_total=round(float(ratio[s].mean()),4),
                     extremity=round(float(ext[s].mean()),3), mean_blocks=round(float(nb[s].mean()),2)))
print(pd.DataFrame(rows).to_string(index=False))
for nm,v in [('off/total (MISFIT)',ratio),('extremity ON manifold',ext)]:
    r=stats.spearmanr(gp[keep],v[keep])
    print(f"  dose-response  {nm:24s} rho={r.statistic:+.4f}  p={r.pvalue:.3g}")
