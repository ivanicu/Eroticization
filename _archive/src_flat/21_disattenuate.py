"""
"Three orthogonal axes" is exactly what measurement error produces. POWER is 3 indicators,
GAZE is 2. Disattenuate: r_true = r_obs / sqrt(rel_a * rel_b). If the axes are really one
folk axis measured badly, disattenuated r goes to ~1. Reliabilities measured, not assumed.
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore'); rng=np.random.default_rng(89)
exec(open('src/16_dimensionality.py').read().split("allq=list(B)")[0])
allq=list(B)
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts()
pool=np.array(sorted(nblk[nblk>=8].index)); pm={p:i for i,p in enumerate(pool)}
G=np.load('data/derived/gcca_G.npy')
A=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)

def alpha(M):
    M=M.dropna(); k=M.shape[1]
    return k/(k-1)*(1-M.var(axis=0,ddof=1).sum()/M.sum(axis=1).var(ddof=1))

pc=[c for c in A.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
Pm=pd.concat([z(pd.to_numeric(A[c],errors='coerce'))*sg[c] for c in pc],axis=1)
ex=[c for c in A.columns if 'exhibition' in c][0]; vo=[c for c in A.columns if 'voyeur' in c][0]
Gm=pd.concat([z(pd.to_numeric(A[ex],errors='coerce')),-z(pd.to_numeric(A[vo],errors='coerce'))],axis=1)
rel_P=float(alpha(Pm)); rel_G=float(alpha(Gm))
print(f"POWER  alpha ({Pm.shape[1]} indicators) = {rel_P:.3f}")
print(f"GAZE   alpha ({Gm.shape[1]} indicators) = {rel_G:.3f}")

# coord4 reliability: split blocks in half, refit GCCA on each half, correlate coord4 scores
from numpy.linalg import lstsq,svd,qr
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
def Zof(q):
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
    Z=np.where(np.isnan(Z),np.nanmean(Z,axis=0),Z); Z=Z-Z.mean(0)
    b,*_=lstsq(D,Z,rcond=None); return Z-D@b
C_={q:Zof(q) for q in allq}
def gc(bl,K=6):
    U,s,_=svd(np.hstack([qr(C_[q],mode='reduced')[0] for q in bl]),full_matrices=False); return U[:,:K]
rs=[]
for _ in range(8):
    p=rng.permutation(allq); h1,h2=list(p[:len(p)//2]),list(p[len(p)//2:])
    G1,G2=gc(h1),gc(h2)
    r=max(abs(np.corrcoef(G1[:,3],G2[:,j])[0,1]) for j in range(6))
    rs.append(2*r/(1+r))
rel_C=float(np.median(rs)); print(f"COORD4 split-half (blocks), Spearman-Brown = {rel_C:.3f}")

X=pd.DataFrame({'P':pd.Series(np.nan,index=range(len(pool))),'G':np.nan,'C':G[:,3]})
Ps=Pm.mean(axis=1).reindex(pool).values; Gs=Gm.mean(axis=1).reindex(pool).values
obs={('P','G'):np.corrcoef(*[v[~np.isnan(Ps)&~np.isnan(Gs)] for v in (Ps,Gs)])[0,1],
     ('P','C'):np.corrcoef(*[v[~np.isnan(Ps)] for v in (Ps,G[:,3])])[0,1],
     ('G','C'):np.corrcoef(*[v[~np.isnan(Gs)] for v in (Gs,G[:,3])])[0,1]}
rel={'P':rel_P,'G':rel_G,'C':rel_C}
print("\n           observed   disattenuated")
for (a,b),r in obs.items():
    d=r/np.sqrt(max(rel[a],1e-6)*max(rel[b],1e-6))
    print(f"  {a}-{b}  :   {r:+.3f}      {d:+.3f}")
mx=max(abs(r/np.sqrt(max(rel[a],1e-6)*max(rel[b],1e-6))) for (a,b),r in obs.items())
print(f"\n  largest disattenuated |r| = {mx:.3f}")
print("  VERDICT:", "still three axes -- orthogonality is not measurement error" if mx<0.4
      else "collapses -- the axes may be one construct measured badly")
