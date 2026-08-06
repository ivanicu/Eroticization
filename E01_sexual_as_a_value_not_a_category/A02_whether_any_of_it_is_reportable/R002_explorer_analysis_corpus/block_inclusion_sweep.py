"""
E01 A07 R03 -- sweep the block-inclusion threshold that 20 rounds inherited and none chose.

#32 mapped the exec graph: 20 of 57 rounds inherit `16_dimensionality`'s block filter --
n_respondents >= 1200, n_options >= 10, mean_picks > 1.5 -- and that filter decides WHICH BLOCKS
EXIST AT ALL. It sits underneath the modality deficit, the consumption deficit, the coverage work
and the theta-vs-coordinates result simultaneously. A07 R01 swept K on these same quantities; the
block filter is the other unswept choice, and it is more fundamental because K reshapes an analysis
while this reshapes the corpus.

ESTIMAND        the congruence deficit for each split, as a function of which blocks are admitted.
IDENTIFICATION  identified at each cell; the random-split ceiling is recomputed inside every cell so
                a corpus that changes size cannot pass as an effect that changes size.
SCOPE           block thresholds this release can support: n>=600 to n>=2000, options>=8 to >=12.
WORLDS          A  the deficits are properties of the domain: stable in magnitude and ordering
                B  they are properties of my corpus cut: they move with the threshold
KILL (CONDITIONAL -- evaluated ONLY if the gate passes)
      gate: sex remains the largest deficit in EVERY cell AND the known-null placebo stays under
            0.03 in EVERY cell
      then: modality and consumption both stay within +/-50% of their published values -> STABLE
            either moves more than 2x, or their ordering flips in >1/3 of cells -> CORPUS-DEPENDENT
            otherwise -> UNVERIFIED
POSITIVE CTRL   sex, the largest documented effect in this dataset.
NEGATIVE CTRL   row-parity placebo, recomputed per cell.
NOISE FLOOR     random-split ceiling sd inside each cell.
MULTIPLICITY    5 respondent thresholds x 2 option thresholds x 4 splits x 2 seeds, all reported.
SEEDS           2 (each cell refits the whole corpus; 3 was over budget and this is stated rather
                than hidden -- see IMPOSSIBLE).
IMPOSSIBLE      3 seeds x the full grid would be ~240 corpus rebuilds; run at 2 and the seed spread
                is reported so the reader can price it.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
ORI=[c for c in df.columns if 'opposite gender to me' in c or 'gender identity' in c.lower()]
X=df[[c for c in ['biomale','_age','opennessvariable','neuroticismvariable','extroversionvariable',
    'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']+ORI if c in df.columns]].copy()
for c in X.columns:
    if X[c].dtype==object: X[c]=X[c].astype('category').cat.codes.replace(-1,np.nan)
X=X.apply(pd.to_numeric,errors='coerce'); X=X.fillna(X.median()); COV=((X-X.mean())/(X.std()+1e-9)).fillna(0.)
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
def corpus(minresp,minopt,MINN=20):
    keep=qm[(~qm.single_pick)&(qm.n_options>=minopt)&(qm.n_respondents>=minresp)&(qm.mean_picks>1.5)]
    B={}
    for _,q in keep.iterrows():
        s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=MINN].index)]
        ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
        if len(ppl)<minresp or len(opt)<8: continue
        pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
        M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
        R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        B[q.qi]=dict(ppl=ppl,R=R)
    return B
def loadings(B,allq,pool,people,K=5):
    ppl=np.array(sorted(set(people)&set(pool)))
    if len(ppl)<400: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
rows=[]
for minresp,minopt,seed in itertools.product([600,800,1000,1200,2000],[8,10],[7,17]):
    B=corpus(minresp,minopt); allq=list(B)
    if len(allq)<8: continue
    pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
    rng=np.random.default_rng(seed)
    male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
    mod=df[MOD].reindex(pool)
    SPL={'sex':(pool[(male==1).values],pool[(male==0).values]),
         'consumption':(pool[(ph>ph.median()).values],pool[(ph<ph.median()).values]),
         'modality':(pool[mod.isin(['Mostly written','Entirely written']).values],
                     pool[mod.isin(['Mostly visual','Entirely visual']).values]),
         'placebo_rowparity':(pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1])}
    # FIX (#33): #11 and #25 established that ANY group comparison on this release must be
    # block-count matched or it measures survey coverage. I built this sweep and did not carry
    # that forward, so the first run compared unmatched deficits against matched published values.
    nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
    def match(g1,g2):
        i1,i2=pd.Index(sorted(set(g1)&set(pool))),pd.Index(sorted(set(g2)&set(pool)))
        k1,k2=nblk.reindex(i1).astype(int),nblk.reindex(i2).astype(int)
        a=[];b=[]
        for v in set(k1)|set(k2):
            x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
            if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
        return np.array(a),np.array(b)
    for name,(g1,g2) in SPL.items():
        g1,g2=match(g1,g2)
        n1,n2=len(g1),len(g2)
        if min(n1,n2)<400: continue
        c=cong(loadings(B,allq,pool,g1),loadings(B,allq,pool,g2))
        ceil=[cong(loadings(B,allq,pool,p[:n1]),loadings(B,allq,pool,p[n1:n1+n2])) for p in [rng.permutation(pool) for _ in range(4)]]
        rows.append(dict(min_resp=minresp,min_opt=minopt,seed=seed,n_blocks=len(allq),
                         split=name,deficit=float(np.nanmean(ceil))-c))
G=pd.DataFrame(rows); G.to_csv(OUT/'block_sweep.csv',index=False)
print("=== blocks admitted at each threshold ===")
print(G.groupby(['min_resp','min_opt']).n_blocks.first().unstack().to_string())
print("\n=== congruence deficit by block-inclusion threshold (median over seeds) ===")
print(G.pivot_table(index=['min_resp','min_opt'],columns='split',values='deficit').round(4).to_string())
P=G.pivot_table(index=['min_resp','min_opt'],columns='split',values='deficit')
gate_pos=bool((P['sex']>=P[['consumption','modality','placebo_rowparity']].max(axis=1)).all())
gate_neg=bool((P['placebo_rowparity'].abs()<0.03).all())
pubm,pubc=0.0546,0.0398
print(f"\n  seed spread (max |diff| between seeds): {G.groupby(['min_resp','min_opt','split']).deficit.apply(lambda s: s.max()-s.min()).max():.4f}")
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  sex largest in every cell            : {'PASS' if gate_pos else 'FAIL'}")
print(f"  placebo under 0.03 in every cell     : {'PASS' if gate_neg else 'FAIL'} (max {P['placebo_rowparity'].abs().max():.4f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED, threshold not evaluated")
else:
    rm=P['modality']/pubm; rc=P['consumption']/pubc
    flips=int((P['consumption']>P['modality']).sum())
    print(f"  modality  range {P['modality'].min():+.4f}..{P['modality'].max():+.4f}  ({rm.min():.2f}x..{rm.max():.2f}x published)")
    print(f"  consumption range {P['consumption'].min():+.4f}..{P['consumption'].max():+.4f}  ({rc.min():.2f}x..{rc.max():.2f}x published)")
    print(f"  cells where consumption > modality: {flips}/{len(P)}")
    if max(rm.max(),rc.max())>2 or rm.min()<0.5 or rc.min()<0.5 or flips>len(P)/3:
        print("  -> CORPUS-DEPENDENT : these deficits move with which blocks are admitted")
    else: print("  -> STABLE : within +/-50% of published across the whole grid")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
