"""
ITER 1. The transferable subspace is real but unnamed. Name it from the data, not from an LLM
(the corpus that would label these options also contains the kink taxonomy = leakage).

Method: one joint low-rank factorization across ALL blocks at once, on OBSERVED cells only
(the survey is a gated tree, so "not asked" != "disliked"; imputing 0 would conflate them).
Within-block double-centering first, so factors are profile SHAPE, not propensity.
Held-out cells give an honest score. Extreme loadings NAME each factor.
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore'); rng=np.random.default_rng(1729)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)&(qm.mean_picks>1.5)]

# ---- build observed-cell table: every option of every block the person ENTERED ----
recs=[]; item_names=[]; iid=0; item_block={}
person_blocks={}
for _,q in keep.iterrows():
    sub=lg[lg.qi==q.qi]; ppl=np.array(sorted(sub.person.unique())); opt=np.array(sorted(sub.option.unique()))
    if len(ppl)<1000: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)      # shape only
    cols=np.arange(iid, iid+len(opt))
    for j,o in enumerate(opt): item_names.append(f"{q.qi}::{o}"); item_block[iid+j]=q.qi
    iid+=len(opt)
    P=np.repeat(ppl,len(opt)); I=np.tile(cols,len(ppl))
    recs.append(np.c_[P, I, R.ravel()])
    for p in ppl: person_blocks.setdefault(p,set()).add(q.qi)
D=np.vstack(recs)
persons=np.unique(D[:,0]); pmap={p:i for i,p in enumerate(persons)}
D[:,0]=np.vectorize(pmap.get)(D[:,0])
n_p, n_i = len(persons), iid
print(f"observed cells {len(D):,}   persons {n_p:,}   items {n_i}   blocks {len(keep)}")
print(f"density {len(D)/(n_p*n_i):.3f}")

def soft_impute(rows,cols,vals,K,iters=25):
    X=np.zeros((n_p,n_i)); mask=np.zeros((n_p,n_i),bool)
    X[rows,cols]=vals; mask[rows,cols]=True
    Z=X.copy()
    for _ in range(iters):
        U,S,Vt=np.linalg.svd(Z,full_matrices=False)
        S2=np.maximum(S-S[K]*0.9,0); S2[K:]=0
        L=(U*S2)@Vt
        Z=np.where(mask,X,L)
    U,S,Vt=np.linalg.svd(Z,full_matrices=False)
    return U[:,:K]*S[:K], Vt[:K], (U[:,:K]*S[:K])@Vt[:K]

rows,cols,vals=D[:,0].astype(int),D[:,1].astype(int),D[:,2]
# ---- held-out cell score ----
te=rng.random(len(rows))<0.2
F,L,Rec=soft_impute(rows[~te],cols[~te],vals[~te],K=8)
pred=Rec[rows[te],cols[te]]; truth=vals[te]
r2=1-((truth-pred)**2).sum()/(truth**2).sum()
print(f"\nheld-out cell R2 (rank 8, trained on 80% of cells) : {r2:+.4f}")

# ---- NULL that preserves the gating tree exactly ----
nb=np.array([len(person_blocks[p]) for p in persons])
strata=pd.qcut(nb,10,labels=False,duplicates='drop')
perm=np.arange(n_p)
for s in np.unique(strata):
    idx=np.flatnonzero(strata==s); perm[idx]=rng.permutation(idx)
rows_n=perm[rows]
te_n=rng.random(len(rows))<0.2
Fn,Ln,Recn=soft_impute(rows_n[~te_n],cols[~te_n],vals[~te_n],K=8)
r2n=1-((vals[te_n]-Recn[rows_n[te_n],cols[te_n]])**2).sum()/(vals[te_n]**2).sum()
print(f"gating-preserving permutation null              : {r2n:+.4f}   <- floor")
print(f"                                          ratio : {r2/max(r2n,1e-6):.1f}x")

# ---- name the factors ----
F,L,_=soft_impute(rows,cols,vals,K=8)
names=np.array(item_names); blk=np.array([item_block[i] for i in range(n_i)])
share=(np.abs(L)>0).sum(1)
print("\n=== FACTOR COORDINATES (top loadings, and how many BLOCKS each factor spans) ===")
for k in range(6):
    l=L[k]; o=np.argsort(l)
    nb_span=len(set(blk[np.argsort(-np.abs(l))[:25]]))
    print(f"\n--- factor {k+1}   spans {nb_span} distinct blocks in its top-25 loadings ---")
    print("  +  "+" | ".join(n.split('::')[1][:34] for n in names[o[::-1]][:7]))
    print("  -  "+" | ".join(n.split('::')[1][:34] for n in names[o][:7]))
np.save('data/derived/joint_F.npy',F); np.save('data/derived/joint_L.npy',L)
pd.Series(item_names).to_csv('data/derived/joint_items.csv',index=False,header=['item'])
