import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
THE test: does a person's residual position in block A predict their position in block B?
Shared grammar => yes. Block-local idiosyncratic preference => no.
Item redundancy cannot produce this (blocks share no items).
Controls: person endorsement rate already removed by double-centering; demographics
partialled out explicitly; person-identity permutation null; positive control on a
KNOWN dimension (sex) so a null here cannot be instrument blindness.
"""
import pandas as pd, numpy as np, pickle, itertools
from numpy.linalg import svd, lstsq
rng = np.random.default_rng(7)
df = pd.read_csv('data/raw/BKSPublic.csv', low_memory=False)
qm = pd.read_csv('data/derived/multiselect_questions.csv')
lg = pd.read_parquet('data/derived/endorsements_long.parquet')
keep = qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1000)&(qm.mean_picks>1.5)]

def block(qi):
    sub=lg[lg.qi==qi]; ppl=sorted(sub.person.unique()); opt=sorted(sub.option.unique())
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    U,S,Vt=svd(R,full_matrices=False)
    return np.array(ppl), U[:,:2]*S[:2], Vt[:2], np.array(opt), M

B={}
for _,q in keep.iterrows():
    ppl,f,load,opt,M = block(q.qi)
    if len(ppl)>=1000: B[q.qi]=dict(ppl=ppl,f=f,load=load,opt=opt,M=M,col=q.col)
print("blocks:", len(B))

# --- covariates for partialling ---
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
df['_age']=df['age'].map(AGEMAP)
cov_cols=['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable',
          'consciensiousnessvariable','agreeablenessvariable']
cov_cols=[c for c in cov_cols if c in df.columns]
COV=df[cov_cols].apply(pd.to_numeric, errors='coerce')
COV=COV.fillna(COV.median()); COV=(COV-COV.mean())/(COV.std()+1e-9); COV=COV.fillna(0.0)
def covmat(idx):
    X=COV.loc[idx].values
    return np.c_[np.ones(len(X)), X]
def resid(y,X):
    b,*_ = lstsq(X,y,rcond=None); return y-X@b

rows=[]
for a,b in itertools.combinations(B,2):
    sa,sb=B[a],B[b]
    common=np.intersect1d(sa['ppl'],sb['ppl'])
    if len(common)<400: continue
    ia=np.searchsorted(sa['ppl'],common); ib=np.searchsorted(sb['ppl'],common)
    fa=sa['f'][ia,0]; fb=sb['f'][ib,0]
    r_raw=np.corrcoef(fa,fb)[0,1]
    X=covmat(common)
    r_adj=np.corrcoef(resid(fa,X),resid(fb,X))[0,1]
    perm=[np.corrcoef(fa,rng.permutation(fb))[0,1] for _ in range(200)]
    rows.append(dict(a=a,b=b,n=len(common),r_raw=round(r_raw,3),r_adj=round(r_adj,3),
                     null_sd=round(float(np.std(perm)),3),
                     z=round(float((abs(r_adj))/ (np.std(perm)+1e-12)),1)))
T=pd.DataFrame(rows)
T.to_csv('data/derived/transfer.csv',index=False)
print("\nblock pairs tested :", len(T), " (median n =", int(T.n.median()), ")")
print("|r_adj| > 0.10     :", int((T.r_adj.abs()>0.10).sum()))
print("|r_adj| > 0.20     :", int((T.r_adj.abs()>0.20).sum()))
print("median |r_raw|     :", round(float(T.r_raw.abs().median()),3))
print("median |r_adj|     :", round(float(T.r_adj.abs().median()),3),
      "   <- after partialling out sex/age/personality")
print("permutation null sd:", round(float(T.null_sd.median()),4))
print("\n--- strongest 14 transfers (PC1 of A  vs  PC1 of B) ---")
for _,r in T.reindex(T.r_adj.abs().sort_values(ascending=False).index).head(14).iterrows():
    print(f"  r_adj={r.r_adj:+.3f} (raw {r.r_raw:+.3f}) n={int(r.n):5d}  {B[r.a]['col'][:40]:42s} <-> {B[r.b]['col'][:40]}")

# ---------- POSITIVE CONTROL: can this pipeline see a KNOWN person dimension? ----------
print("\n=== POSITIVE CONTROL: does block PC1 recover sex (biomale)? ===")
pc=[]
for qi,s in B.items():
    y=pd.to_numeric(df.loc[s['ppl'],'biomale'],errors='coerce').values
    m=~np.isnan(y)
    if m.sum()>500: pc.append(abs(np.corrcoef(s['f'][m,0],y[m])[0,1]))
print(f"  |r(PC1, biomale)| across {len(pc)} blocks: median {np.median(pc):.3f}  max {np.max(pc):.3f}")
print("  -> instrument is NOT blind to person-level dimensions" if np.max(pc)>0.15 else "  -> WARNING: instrument may be blind")
