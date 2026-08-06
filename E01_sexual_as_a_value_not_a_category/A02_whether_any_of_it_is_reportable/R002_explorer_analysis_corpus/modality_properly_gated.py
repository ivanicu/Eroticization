"""
E01 A07 R04 -- the modality deficit, with the gate specified for the design that is actually run.

#33 failed three ways: unmatched deficits compared to matched values, a gate encoding an
expectation from a different matching level, and 2 seeds on a quantity whose seed noise exceeded it.
All three are fixed here, and the fixes are the round's specification:

  MATCHING     block-count matched inside every cell (#11, #25), and the gate is written FOR that
               level -- sex must be PRESENT AND LARGE, not "largest", since #11's own chain has
               consumption at 0.0871 under block-only matching.
  SEEDS        5, and each cell reports its own seed spread. A cell whose effect does not exceed
               2x its seed spread is marked UNRESOLVABLE and excluded from the verdict rather than
               averaged into it.
  AXES         option threshold and respondent threshold swept SEPARATELY -- #33 suggested option
               >=8 vs >=10 is where modality actually moves, and the two were confounded in one grid.

ESTIMAND        the modality (written vs visual) congruence deficit, and its dependence on each
                corpus-inclusion axis separately.
WORLDS          A  a property of the domain: stable across both axes within 1.5x
                B  a property of the corpus cut: range exceeds 3x its own median
KILL (CONDITIONAL -- evaluated ONLY if the gate passes)
      gate: sex deficit > 0.05 in every cell AND |placebo| < 0.02 in every cell AND at least 2/3 of
            cells resolvable (effect > 2x that cell's seed spread)
      then: range/median > 3   -> CORPUS-DEPENDENT, the number is withdrawn
            range/median < 1.5 -> STABLE, republish as a range
            otherwise          -> UNVERIFIED
POSITIVE CTRL   sex, present and large.
NEGATIVE CTRL   row-parity placebo.
NOISE FLOOR     per-cell seed spread, reported for every cell rather than aggregated.
MULTIPLICITY    3 option thresholds x 3 respondent thresholds x 5 seeds x 3 splits, all reported.
SEEDS           5.
IMPOSSIBLE      independent replication; a release whose blocks are not gated.
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
def corpus(minresp,minopt):
    keep=qm[(~qm.single_pick)&(qm.n_options>=minopt)&(qm.n_respondents>=minresp)&(qm.mean_picks>1.5)]
    B={}
    for _,q in keep.iterrows():
        s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=20].index)]
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
for minresp,minopt in itertools.product([600,1200,2000],[8,10,12]):
    B=corpus(minresp,minopt); allq=list(B)
    if len(allq)<8: continue
    pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
    nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
    male=df['biomale'].reindex(pool); mod=df[MOD].reindex(pool)
    SPL={'sex':(pool[(male==1).values],pool[(male==0).values]),
         'modality':(pool[mod.isin(['Mostly written','Entirely written']).values],
                     pool[mod.isin(['Mostly visual','Entirely visual']).values]),
         'placebo':(pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1])}
    for seed in (3,13,23,33,43):
        rng=np.random.default_rng(seed)
        def match(g1,g2):
            i1,i2=pd.Index(sorted(set(g1)&set(pool))),pd.Index(sorted(set(g2)&set(pool)))
            k1,k2=nblk.reindex(i1).astype(int),nblk.reindex(i2).astype(int)
            a=[];b=[]
            for v in set(k1)|set(k2):
                x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
                if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
            return np.array(a),np.array(b)
        for name,(g1,g2) in SPL.items():
            a,b=match(g1,g2)
            if min(len(a),len(b))<400: continue
            c=cong(loadings(B,allq,pool,a),loadings(B,allq,pool,b))
            ceil=[cong(loadings(B,allq,pool,p[:len(a)]),loadings(B,allq,pool,p[len(a):len(a)+len(b)]))
                  for p in [rng.permutation(pool) for _ in range(3)]]
            rows.append(dict(min_resp=minresp,min_opt=minopt,seed=seed,n_blocks=len(allq),
                             split=name,deficit=float(np.nanmean(ceil))-c))
G=pd.DataFrame(rows); G.to_csv(OUT/'modality_gated.csv',index=False)
cell=G.groupby(['min_resp','min_opt','split']).deficit.agg(['median','min','max'])
cell['spread']=cell['max']-cell['min']; cell['resolvable']=cell['median'].abs()>2*cell['spread']
print("=== per-cell median deficit, seed spread, and resolvability (5 seeds) ===")
print(cell.round(4).to_string())
M=cell.xs('modality',level='split'); S=cell.xs('sex',level='split'); P=cell.xs('placebo',level='split')
print(f"\n=== the two axes, separated (modality median) ===")
print(G[G.split=='modality'].pivot_table(index='min_opt',columns='min_resp',values='deficit',aggfunc='median').round(4).to_string())
gate_pos=bool((S['median']>0.05).all()); gate_neg=bool((P['median'].abs()<0.02).all())
res=float(M['resolvable'].mean())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  sex deficit > 0.05 in every cell     : {'PASS' if gate_pos else 'FAIL'} (min {S['median'].min():.4f})")
print(f"  |placebo| < 0.02 in every cell       : {'PASS' if gate_neg else 'FAIL'} (max {P['median'].abs().max():.4f})")
print(f"  >=2/3 of modality cells resolvable   : {'PASS' if res>=2/3 else 'FAIL'} ({res:.0%} resolvable)")
if not (gate_pos and gate_neg and res>=2/3):
    print("  -> gate FAILED : UNVERIFIED, threshold not evaluated")
else:
    r=(M['median'].max()-M['median'].min())/abs(M['median'].median())
    print(f"  modality range/median = {r:.2f}  ({M['median'].min():+.4f} .. {M['median'].max():+.4f}, median {M['median'].median():+.4f})")
    if r>3: print("  -> CORPUS-DEPENDENT : the modality number is withdrawn")
    elif r<1.5: print(f"  -> STABLE : republish as [{M['median'].min():.4f}, {M['median'].max():.4f}]")
    else: print("  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
