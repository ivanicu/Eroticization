"""
A: CONFOUND CONTROL -- does the cross-domain subspace transfer survive removing
   sex / age / personality / orientation-proxies from the item matrices themselves?
B: NAMED DIMENSION -- the survey asks the SAME agent-vs-patient contrast in 6+ separate
   domains. If a domain-general receptive/agentive coordinate exists, these cohere.
   No LLM, no invented ontology: the contrast is the survey's own.
"""
import pandas as pd, numpy as np, itertools, warnings
from numpy.linalg import lstsq
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore'); rng=np.random.default_rng(23)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')

AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
ORI=[c for c in df.columns if 'opposite gender to me' in c or 'gender identity' in c.lower()]
cc=['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable',
    'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']+ORI
X=df[[c for c in cc if c in df.columns]].copy()
for c in X.columns:
    if X[c].dtype==object: X[c]=X[c].astype('category').cat.codes.replace(-1,np.nan)
X=X.apply(pd.to_numeric,errors='coerce'); X=X.fillna(X.median()); X=((X-X.mean())/(X.std()+1e-9)).fillna(0.)
print("covariates used:", list(X.columns))

def strip(Mat, idx):
    D=np.c_[np.ones(len(idx)), X.loc[idx].values]
    b,*_=lstsq(D,Mat,rcond=None); return Mat-D@b

keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)&(qm.mean_picks>1.5)]
B={}
for _,q in keep.iterrows():
    sub=lg[lg.qi==q.qi]; ppl=np.array(sorted(sub.person.unique())); opt=np.array(sorted(sub.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    if len(ppl)<1000: continue
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    B[q.qi]=dict(ppl=ppl,R=R,col=q.col)

def cvcca(Xa,Xb,nc):
    idx=rng.permutation(len(Xa)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    c=CCA(n_components=nc,max_iter=800).fit(Xa[tr],Xb[tr]); ua,ub=c.transform(Xa[te],Xb[te])
    return max(abs(np.corrcoef(ua[:,j],ub[:,j])[0,1]) for j in range(nc))

rows=[]
for a,b in itertools.combinations(B,2):
    A_,B_=B[a],B[b]; common=np.intersect1d(A_['ppl'],B_['ppl'])
    if len(common)<600: continue
    ia=np.searchsorted(A_['ppl'],common); ib=np.searchsorted(B_['ppl'],common)
    Ra,Rb=A_['R'][ia],B_['R'][ib]
    nc=min(3,Ra.shape[1]-1,Rb.shape[1]-1)
    try:
        raw=cvcca(Ra,Rb,nc); adj=cvcca(strip(Ra,common),strip(Rb,common),nc)
        nul=np.mean([cvcca(Ra,Rb[rng.permutation(len(common))],nc) for _ in range(8)])
    except Exception: continue
    rows.append(dict(n=len(common),cca_raw=raw,cca_adj=adj,cca_null=nul))
T=pd.DataFrame(rows)
print(f"\n=== A. does cross-domain grammar survive demographics? ({len(T)} pairs) ===")
print(f"  held-out CCA, raw                 : {T.cca_raw.median():.3f}")
print(f"  held-out CCA, demographics removed: {T.cca_adj.median():.3f}")
print(f"  permutation floor                 : {T.cca_null.median():.3f}")
print(f"  ratio adj/floor                   : {T.cca_adj.median()/T.cca_null.median():.2f}x")
print(f"  share of raw transfer that is NOT demographic: {T.cca_adj.median()/T.cca_raw.median():.0%}")

# ---------------- B. the agent/patient axis ----------------
print("\n=== B. is 'am I the agent or the patient?' one domain-general coordinate? ===")
PAIRS=[('exhibitionself','exhibitionother'),('voyeurself','voyeurother'),
       ('worshipped','worshipping'),('receivepain','givepain')]
BEG=[c for c in df.columns if 'eagerly beg' in c]
if len(BEG)==2: PAIRS.append((BEG[0],BEG[1]))
ME=[c for c in df.columns if set(map(str,df[c].dropna().unique()))=={'Me','Someone else'}]
cols={}
for s,o in PAIRS:
    v=pd.to_numeric(df[s],errors='coerce')-pd.to_numeric(df[o],errors='coerce')
    cols[f"{s}-{o}"[:38]]=v
for c in ME:
    cols['ME:'+c[:34]]=(df[c]=='Me').astype(float).where(df[c].notna())
A=pd.DataFrame(cols)
print(f"  indicators: {A.shape[1]}   pairwise-complete correlation matrix:")
C=A.corr(min_periods=200)
print(C.round(2).to_string())
ev=np.linalg.eigvalsh(C.fillna(0).values)[::-1]
print(f"\n  eigenvalues: {np.round(ev,2)}")
print(f"  PC1 explains {ev[0]/len(ev):.1%} of {len(ev)} indicators  (chance = {1/len(ev):.1%})")
off=C.values[np.triu_indices_from(C.values,1)]; off=off[~np.isnan(off)]
print(f"  mean off-diagonal r = {off.mean():+.3f}   share positive = {(off>0).mean():.0%}   n_pairs={len(off)}")
A.to_csv('data/derived/agent_patient.csv',index=False)
