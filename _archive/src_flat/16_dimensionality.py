"""
ITER 2. How MANY shared coordinates are there? Naming an unshared coordinate is a story.
Test: split the 32 BLOCKS into two disjoint halves, fit person factors on each half
independently, then canonically correlate the two person-score sets on common people.
A coordinate that is domain-general must be recoverable from EITHER half of the domains.
The number of canonical correlations above the block-permutation floor IS the dimensionality.
Sex/age/personality/orientation projected out of the person space FIRST, so coordinate 1
is not just sex. Singleton options floor-filtered (they entered the last run as noise columns).
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq, svd
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore'); rng=np.random.default_rng(2718)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]

MINN=20
B={}; dropped=0
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); good=set(vc[vc>=MINN].index); dropped+=int((vc<MINN).sum())
    s=s[s.option.isin(good)]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values, s.option.map(oi).values]=1
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    B[q.qi]=dict(ppl=ppl,R=R,opt=opt,col=q.col)
print(f"blocks {len(B)}   singleton/rare options dropped (n<{MINN}): {dropped}")

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
ORI=[c for c in df.columns if 'opposite gender to me' in c or 'gender identity' in c.lower()]
cc=['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable',
    'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']+ORI
X=df[[c for c in cc if c in df.columns]].copy()
for c in X.columns:
    if X[c].dtype==object: X[c]=X[c].astype('category').cat.codes.replace(-1,np.nan)
X=X.apply(pd.to_numeric,errors='coerce'); X=X.fillna(X.median()); COV=((X-X.mean())/(X.std()+1e-9)).fillna(0.)

def factors(blocks, K=8, ppl_pool=None):
    ppl=ppl_pool if ppl_pool is not None else np.unique(np.concatenate([B[q]['ppl'] for q in blocks]))
    pm={p:i for i,p in enumerate(ppl)}
    cols=[]
    for q in blocks:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm])
        src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan); Z[idx]=B[q]['R'][src]
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]        # project out demographics FIRST
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    U,S,Vt=svd(Z,full_matrices=False)
    return ppl,(U[:,:K]*S[:K]),Vt[:K],np.hstack([[q]*B[q]['R'].shape[1] for q in blocks]),\
           np.hstack([B[q]['opt'] for q in blocks])

allq=list(B); K=8
pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
obs=[]; nul=[]
for rep in range(12):
    perm=rng.permutation(allq); h1,h2=list(perm[:len(perm)//2]),list(perm[len(perm)//2:])
    p1,F1,_,_,_=factors(h1,K,pool); p2,F2,_,_,_=factors(h2,K,pool)
    idx=rng.permutation(len(pool)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    c=CCA(n_components=K,max_iter=900).fit(F1[tr],F2[tr]); a,b_=c.transform(F1[te],F2[te])
    obs.append([abs(np.corrcoef(a[:,j],b_[:,j])[0,1]) for j in range(K)])
    sh=rng.permutation(len(pool))
    c2=CCA(n_components=K,max_iter=900).fit(F1[tr],F2[sh][tr]); a2,b2=c2.transform(F1[te],F2[sh][te])
    nul.append([abs(np.corrcoef(a2[:,j],b2[:,j])[0,1]) for j in range(K)])
O=np.array(obs).mean(0); N=np.array(nul).mean(0)
print("\n=== how many coordinates survive a BLOCK split-half? (12 random splits, held-out) ===")
print("  dim :  observed   floor   ratio")
for j in range(K): print(f"   {j+1}  :   {O[j]:.3f}    {N[j]:.3f}   {O[j]/max(N[j],1e-9):5.1f}x")
ndim=int((O>3*N).sum())
print(f"\n  coordinates above 3x floor : {ndim}")
print(f"  coordinates above .20      : {int((O>0.20).sum())}")
np.save('data/derived/dim_obs.npy',O); np.save('data/derived/dim_null.npy',N)
