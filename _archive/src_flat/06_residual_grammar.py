"""
Q: after removing (a) item base rate and (b) each person's overall endorsement count,
   does ANY shared person-level structure remain -- and does it transfer across blocks?

No LLM features are used anywhere here, so LLM/corpus leakage cannot produce the result.
Nulls: degree-preserving (Chung-Lu style) rewiring that fixes BOTH person counts and item
counts, so anything the null reproduces is marginal structure, not grammar.
"""
import pandas as pd, numpy as np, json
from numpy.linalg import svd
rng = np.random.default_rng(20260802)

qm = pd.read_csv('data/derived/multiselect_questions.csv')
lg = pd.read_parquet('data/derived/endorsements_long.parquet')
df = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
keep = qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)]

def block_matrix(qi):
    sub = lg[lg.qi==qi]
    people = sorted(sub.person.unique()); opts = sorted(sub.option.unique())
    pi={p:i for i,p in enumerate(people)}; oi={o:i for i,o in enumerate(opts)}
    M=np.zeros((len(people),len(opts)),dtype=np.float64)
    M[sub.person.map(pi).values, sub.option.map(oi).values]=1.0
    return M, np.array(people), np.array(opts)

def curveball(M, swaps_per_edge=6):
    """Degree-preserving randomisation: keeps every person's count AND every item's count.
    Anything surviving this null is NOT explainable by 'this person says yes a lot' or
    'this item is popular'."""
    A = M.copy().astype(bool); n,k = A.shape
    n_sw = int(swaps_per_edge*A.sum())
    for _ in range(n_sw):
        i,j = rng.integers(0,n,2)
        if i==j: continue
        a, b = A[i], A[j]
        di = np.flatnonzero(a & ~b); dj = np.flatnonzero(b & ~a)
        m = min(len(di),len(dj))
        if m==0: continue
        t = rng.integers(1,m+1)
        ci = rng.choice(di,t,replace=False); cj = rng.choice(dj,t,replace=False)
        A[i,ci]=False; A[i,cj]=True; A[j,cj]=False; A[j,ci]=True
    return A.astype(np.float64)

def residualise(M):
    """remove item base rate AND person endorsement rate (the two nuisances)"""
    R = M - M.mean(0, keepdims=True)
    R = R - R.mean(1, keepdims=True)
    return R

def spectrum(R, ncomp=6):
    s = svd(R, compute_uv=False)
    tot = (R**2).sum()
    return (s[:ncomp]**2)/tot

rows=[]; scores={}
for _,q in keep.iterrows():
    M,people,opts = block_matrix(q.qi)
    if M.shape[0] < 1000: continue
    R = residualise(M)
    obs = spectrum(R)
    nulls = np.array([spectrum(residualise(curveball(M))) for _ in range(12)])
    z = (obs - nulls.mean(0)) / (nulls.std(0)+1e-12)
    U,S,Vt = svd(R, full_matrices=False)
    scores[q.qi] = dict(people=people, f=U[:,:3]*S[:3], load=Vt[:3], opts=opts)
    rows.append(dict(qi=q.qi, col=q.col[:52], n=M.shape[0], k=M.shape[1],
        pc1=round(obs[0],4), pc1_null=round(nulls[:,0].mean(),4), pc1_z=round(z[0],1),
        pc2=round(obs[1],4), pc2_z=round(z[1],1), pc3_z=round(z[2],1)))
res = pd.DataFrame(rows).sort_values('pc1_z', ascending=False)
print("=== Is there residual structure beyond BOTH marginals? (z vs degree-preserving null, 12 draws) ===")
print(res.to_string(index=False))
print("\nblocks with PC1 z>5 :", int((res.pc1_z>5).sum()), "of", len(res))
print("blocks with PC2 z>5 :", int((res.pc2_z>5).sum()), "of", len(res))
np.save('data/derived/block_scores.npy', np.array([1]))
import pickle; pickle.dump(scores, open('data/derived/block_scores.pkl','wb'))
res.to_csv('data/derived/residual_spectrum.csv', index=False)
