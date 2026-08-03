import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Attack each GCCA coordinate with its OWN matched control: find every block contributing to its
top-25 loadings, DELETE THEM ALL, refit, and ask whether the coordinate is still there.
Drop-one-block was too weak -- coord 4 reads as self/other but its loadings are the 7-block
fluid family, which shares one option template. If the coordinate is a template artifact it
cannot survive deleting the template.
Also: is coord 4 just the POWER composite from iteration 0 re-derived?
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(577)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index)); pm={p:i for i,p in enumerate(pool)}
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
def Zof(q):
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
    Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0),Z); Z=Z-Z.mean(0)
    b,*_=lstsq(D,Z,rcond=None); return Z-D@b
CACHE={q:Zof(q) for q in allq}
def gcca(blocks,K=6):
    U,s,_=svd(np.hstack([qr(CACHE[q],mode='reduced')[0] for q in blocks]),full_matrices=False)
    return U[:,:K], s[:K]**2/len(blocks)
G,ev=gcca(allq)

def loading_blocks(k,top=25):
    c=[];b=[]
    for q in allq:
        Z=CACHE[q]
        for j in range(Z.shape[1]): c.append(np.corrcoef(G[:,k],Z[:,j])[0,1]); b.append(q)
    c=np.array(c);b=np.array(b); t=np.argsort(-np.abs(c))[:top]
    return pd.Series(b[t]).value_counts()

FLUID=[q for q in allq if 'For ' in B[q]['col'][:4] or B[q]['col'].startswith('For ')]
print("fluid-family blocks:",FLUID,[B[q]['col'][:28] for q in FLUID])
print()
rows=[]
for k in range(1,5):
    lb=loading_blocks(k); contrib=list(lb.index)
    remain=[q for q in allq if q not in contrib]
    if len(remain)<6: remain=[q for q in allq if q not in contrib[:3]]
    G2,_=gcca(remain)
    stab=max(abs(np.corrcoef(G[:,k],G2[:,j])[0,1]) for j in range(6))
    # also: drop the whole fluid family regardless
    G3,_=gcca([q for q in allq if q not in FLUID])
    stabF=max(abs(np.corrcoef(G[:,k],G3[:,j])[0,1]) for j in range(6))
    rows.append(dict(coord=k+1, agree=round(ev[k],3), n_contrib_blocks=len(contrib),
                     stab_drop_contributors=round(stab,2), stab_drop_fluidfamily=round(stabF,2)))
    print(f"coord {k+1}: top-25 loadings come from {len(contrib)} blocks {contrib[:6]}")
T=pd.DataFrame(rows); print("\n"+T.to_string(index=False))

A=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in A.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=pd.concat([z(pd.to_numeric(A[c],errors='coerce'))*sg[c] for c in pc],axis=1).mean(axis=1).reindex(pool).values
m=~np.isnan(POWER)
print("\n=== is any coordinate just the iteration-0 POWER composite re-derived? ===")
for k in range(6):
    print(f"  coord {k+1}: |r(GCCA, POWER)| = {abs(np.corrcoef(G[m,k],POWER[m])[0,1]):.3f}")
T.to_csv('data/derived/coord_attack.csv',index=False)
