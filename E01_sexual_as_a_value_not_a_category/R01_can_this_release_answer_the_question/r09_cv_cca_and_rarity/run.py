import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Three separate claims, measured separately:
 (0) PROPENSITY  : does endorsement RATE transfer across domains?  (Ivan's model 0, alpha_i)
 (1) GRAMMAR     : does the residual SUBSPACE transfer?  cross-validated CCA (in-sample CCA is
                   guaranteed positive, so only held-out canonical r is admissible)
 (2) RARITY/n_i  : does a person's taste for RARE options transfer, after removing how many
                   options they picked?  This is the transgression parameter.
"""
import pandas as pd, numpy as np, itertools, warnings
from numpy.linalg import svd, lstsq
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore')
rng=np.random.default_rng(11)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)&(qm.mean_picks>1.5)]

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
cc=['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable',
    'consciensiousnessvariable','agreeablenessvariable']
COV=df[[c for c in cc if c in df.columns]].apply(pd.to_numeric,errors='coerce')
COV=COV.fillna(COV.median()); COV=((COV-COV.mean())/(COV.std()+1e-9)).fillna(0.0)
def dres(y,idx):
    X=np.c_[np.ones(len(idx)),COV.loc[idx].values]
    b,*_=lstsq(X,y,rcond=None); return y-X@b

B={}
for _,q in keep.iterrows():
    sub=lg[lg.qi==q.qi]; ppl=np.array(sorted(sub.person.unique())); opt=np.array(sorted(sub.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    if len(ppl)<1000: continue
    rate=M.mean(0)                                   # item base rate
    cnt = M.sum(1)                                   # how many the person picked
    rar = (M@np.log(1/np.clip(rate,1e-3,1)))/np.clip(cnt,1,None)   # mean rarity of picks
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    B[q.qi]=dict(ppl=ppl,M=M,R=R,cnt=cnt,rar=rar,col=q.col,k=M.shape[1])
print("blocks:",len(B))

rows=[]
for a,b in itertools.combinations(B,2):
    A,Bb=B[a],B[b]
    common=np.intersect1d(A['ppl'],Bb['ppl'])
    if len(common)<600: continue
    ia=np.searchsorted(A['ppl'],common); ib=np.searchsorted(Bb['ppl'],common)
    # (0) propensity: endorsement RATE (fraction of options picked)
    pa=A['cnt'][ia]/A['k']; pb=Bb['cnt'][ib]/Bb['k']
    r_prop=np.corrcoef(dres(pa,common),dres(pb,common))[0,1]
    # (2) rarity taste, with count removed
    ra=dres(A['rar'][ia]-np.polyval(np.polyfit(pa,A['rar'][ia],2),pa), common)
    rb=dres(Bb['rar'][ib]-np.polyval(np.polyfit(pb,Bb['rar'][ib],2),pb), common)
    r_rar=np.corrcoef(ra,rb)[0,1]
    # (1) grammar: cross-validated CCA on residual subspaces
    Xa=A['R'][ia]; Xb=Bb['R'][ib]
    idx=rng.permutation(len(common)); tr,te=idx[:len(idx)//2], idx[len(idx)//2:]
    nc=min(3,Xa.shape[1]-1,Xb.shape[1]-1)
    try:
        c=CCA(n_components=nc,max_iter=800).fit(Xa[tr],Xb[tr])
        ua,ub=c.transform(Xa[te],Xb[te])
        cv=[abs(np.corrcoef(ua[:,j],ub[:,j])[0,1]) for j in range(nc)]
    except Exception: cv=[np.nan]
    # permutation null for CV-CCA
    nulls=[]
    for _ in range(20):
        pmt=rng.permutation(len(common))
        try:
            c2=CCA(n_components=nc,max_iter=300).fit(Xa[tr],Xb[pmt][tr])
            u2,v2=c2.transform(Xa[te],Xb[pmt][te])
            nulls.append(max(abs(np.corrcoef(u2[:,j],v2[:,j])[0,1]) for j in range(nc)))
        except Exception: pass
    rows.append(dict(a=a,b=b,n=len(common),r_prop=round(r_prop,3),r_rar=round(r_rar,3),
                     cca1=round(cv[0],3), cca_max=round(float(np.nanmax(cv)),3),
                     cca_null=round(float(np.mean(nulls)) if nulls else np.nan,3)))
T=pd.DataFrame(rows); T.to_csv('data/derived/cca_rarity.csv',index=False)
q=lambda s: (round(float(s.median()),3), round(float(s.quantile(.9)),3))
print(f"\npairs: {len(T)}   median n={int(T.n.median())}")
print(f"(0) PROPENSITY transfer  r : median {q(T.r_prop)[0]:+.3f}   p90 {q(T.r_prop)[1]:+.3f}   |r|>.2 in {int((T.r_prop.abs()>.2).sum())}/{len(T)}")
print(f"(2) RARITY-taste transfer r: median {q(T.r_rar)[0]:+.3f}   p90 {q(T.r_rar)[1]:+.3f}   |r|>.2 in {int((T.r_rar.abs()>.2).sum())}/{len(T)}")
print(f"(1) GRAMMAR  held-out CCA  : median {q(T.cca_max)[0]:.3f}   p90 {q(T.cca_max)[1]:.3f}")
print(f"    permutation null CCA   : median {round(float(T.cca_null.median()),3)}  <- the floor")
print(f"    pairs above null       : {int((T.cca_max>T.cca_null).sum())}/{len(T)}")
print(f"    pairs 2x above null    : {int((T.cca_max>2*T.cca_null).sum())}/{len(T)}")
