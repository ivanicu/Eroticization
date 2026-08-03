"""
(a) Language contamination in block 1 (the biggest block, n=15,250): are the French-variant
    endorsers a DISJOINT subsample seeing a different instrument?
(b) The right estimand for "shared coordinates": LEAVE-ONE-BLOCK-OUT. Fit person factors on
    all other blocks, predict the held-out block's within-person profile. Cell-level R2 was
    saturated by item covariance and could not see the person effect; this can.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(4242)
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
qm=pd.read_csv('data/derived/multiselect_questions.csv')

print("=== (a) block 1 language contamination ===")
b1=lg[lg.qi==1]
FR=['Creampies (éjaculer dans un vagin)','Destination vaginale','Fessée','Rimmissement']
EN=['Creampies (ejaculating into a vagina)','Spanking','Rimming']
fr_p=set(b1[b1.option.isin(FR)].person); en_p=set(b1[b1.option.isin(EN)].person)
allp=set(b1.person)
print(f"  respondents in block 1        : {len(allp):,}")
print(f"  endorsed >=1 FRENCH variant   : {len(fr_p):,}")
print(f"  endorsed >=1 ENGLISH twin     : {len(en_p):,}")
print(f"  overlap                       : {len(fr_p&en_p):,}  ({len(fr_p&en_p)/max(len(fr_p),1):.1%} of French endorsers)")
# how many options did each subsample pick on average -> different instrument?
cnt=b1.groupby('person').size()
print(f"  mean picks, French endorsers  : {cnt.reindex(list(fr_p)).mean():.2f}")
print(f"  mean picks, others            : {cnt.reindex(list(allp-fr_p)).mean():.2f}")

print("\n=== (b) LEAVE-ONE-BLOCK-OUT: do person factors from OTHER blocks predict a new block? ===")
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
B={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values, s.option.map(oi).values]=1
    R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
    B[q.qi]=dict(ppl=ppl,R=R,col=q.col,prop=M.mean(1))

def person_features(target, K=6):
    """rank-K factors from all blocks EXCEPT target, over the union of people"""
    others=[q for q in B if q!=target]
    ppl=np.unique(np.concatenate([B[q]['ppl'] for q in others]))
    pm={p:i for i,p in enumerate(ppl)}
    cols=[]; 
    for q in others:
        idx=np.array([pm[p] for p in B[q]['ppl']])
        X=np.full((len(ppl),B[q]['R'].shape[1]),np.nan); X[idx]=B[q]['R']
        cols.append(X)
    X=np.hstack(cols)
    col_mu=np.nanmean(X,axis=0); X=np.where(np.isnan(X),col_mu,X)
    X=X-X.mean(0,keepdims=True)
    U,S,Vt=np.linalg.svd(X,full_matrices=False)
    return ppl,(U[:,:K]*S[:K])

rows=[]
for t in B:
    ppl_o,Fo=person_features(t)
    pmap={p:i for i,p in enumerate(ppl_o)}
    tgt=B[t]; common=np.array([p for p in tgt['ppl'] if p in pmap])
    if len(common)<800: continue
    ia=np.searchsorted(tgt['ppl'],common); Y=tgt['R'][ia]
    Fx=Fo[[pmap[p] for p in common]]
    prop=tgt['prop'][ia]
    idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
    def r2(X):
        b,*_=lstsq(X[tr],Y[tr],rcond=None); P=X[te]@b
        return 1-((Y[te]-P)**2).sum()/(Y[te]**2).sum()
    base=r2(np.c_[np.ones(len(common)),prop])
    full=r2(np.c_[np.ones(len(common)),prop,Fx])
    nul=np.mean([r2(np.c_[np.ones(len(common)),prop,Fx[rng.permutation(len(common))]]) for _ in range(12)])
    rows.append(dict(qi=t,col=tgt['col'][:44],n=len(common),base=round(base,4),
                     gain=round(full-base,4), null_gain=round(nul-base,4)))
T=pd.DataFrame(rows).sort_values('gain',ascending=False)
T.to_csv('data/derived/lobo.csv',index=False)
print(T.to_string(index=False))
print(f"\n  blocks: {len(T)}")
print(f"  median GAIN from other-block person factors : {T.gain.median():+.4f}")
print(f"  median PERMUTED-factor gain (floor)         : {T.null_gain.median():+.4f}")
print(f"  blocks with gain > 0                        : {int((T.gain>0).sum())}/{len(T)}")
print(f"  blocks with gain > 3x |floor|               : {int((T.gain>3*T.null_gain.abs()).sum())}/{len(T)}")
