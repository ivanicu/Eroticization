import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
SVD on the concatenation is variance-weighted, so the largest block (qi=1, n=14,901) owns every
factor -- drop-top-block stability was 0.25-0.70. A SHARED coordinate is not the direction with
the most total variance, it is the direction best represented in EVERY block.
That is MAXVAR generalized CCA: eigenvectors of sum_b P_b, P_b = projector onto block b's column
space. Block size enters only through the projector's rank, not through variance.
Same attack as before: drop the biggest block and see if the coordinate survives.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore'); rng=np.random.default_rng(161803)
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])

allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index))
print(f"pool = people in >=8 of {len(allq)} blocks : {len(pool):,}")
pm={p:i for i,p in enumerate(pool)}
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]

def whitened(q):
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm])
    src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
    Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0),Z); Z=Z-Z.mean(0)
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b          # demographics out
    Qm,_=qr(Z, mode='reduced')                    # orthonormal basis => projector Q Q^T
    return Qm

def gcca(blocks,K=6):
    S=np.zeros((len(pool),len(pool))) if False else None
    # avoid n x n: accumulate via stacking bases and SVD of [Q1 Q2 ... Qb]
    Qs=np.hstack([whitened(q) for q in blocks])
    U,s,_=svd(Qs,full_matrices=False)
    return U[:,:K], s[:K]**2/len(blocks)          # eigenvalue/ nblocks = mean cos^2 agreement

G,ev=gcca(allq,6)
print("mean squared agreement per coordinate (1.0 = present in every block):")
print("  ", np.round(ev,3))

# --- attack: drop the biggest block ---
G2,_=gcca([q for q in allq if q!=1],6)
print("\ndrop-block-1 stability of each GCCA coordinate:")
stab=[max(abs(np.corrcoef(G[:,k],G2[:,j])[0,1]) for j in range(6)) for k in range(6)]
print("  ", np.round(stab,2))

# --- null: shuffle persons independently within each block, redo ---
def gcca_null(blocks,K=6):
    Qs=[]
    for q in blocks:
        Qm=whitened(q); Qs.append(Qm[rng.permutation(len(pool))])
    U,s,_=svd(np.hstack(Qs),full_matrices=False)
    return s[:K]**2/len(blocks)
nl=np.mean([gcca_null(allq) for _ in range(5)],axis=0)
print("\nnull (persons shuffled independently per block):", np.round(nl,3))
print("ratio observed/null                            :", np.round(ev/nl,2))

# --- name them ---
print("\n=== GCCA COORDINATES ===")
for k in range(4):
    corrs=[]; names=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm])
        src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
        Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0),Z)
        for j,o in enumerate(B[q]['opt']):
            corrs.append(np.corrcoef(G[:,k],Z[:,j])[0,1]); names.append(f"{o}")
    corrs=np.array(corrs); names=np.array(names); o_=np.argsort(corrs)
    print(f"\n--- coord {k+1}  (agreement {ev[k]:.3f}, stability to dropping block 1: {stab[k]:.2f}) ---")
    print("   +  "+" | ".join(n[:31] for n in names[o_[::-1]][:9]))
    print("   -  "+" | ".join(n[:31] for n in names[o_][:9]))
np.save('data/derived/gcca_G.npy',G)
