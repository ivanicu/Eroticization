"""
E01 A04 R12 -- does coverage matching CORRECT anything, or does it just subsample?

A04's decision has two halves. The first (self-reported induction is unusable as an outcome) rests
on timing nulls that already carry planted-effect controls. The second -- "group comparisons must be
coverage-matched" -- rests on the +0.815 coverage law that #20 WITHDREW as one influential point,
and #45 has since shown matching can shrink an effect below its own resolvability. The rule is
therefore standing on a retracted premise plus a known cost, and has never been tested directly.

The direct test: matching does two things at once. It removes coverage imbalance AND it discards
people, which adds noise and can move an estimate on its own. Separate them with a PLACEBO MATCH --
match on a random variable, discarding the same number of people while correcting nothing.

  if matching CORRECTS  -> the deficit change tracks the split's coverage gap, and the placebo match
                           changes little
  if it just SUBSAMPLES -> the placebo match changes the deficit as much as the real match does

ESTIMAND        change in congruence deficit from unmatched to matched, per split, against (a) the
                split's coverage gap and (b) a placebo match of equal sample cost.
IDENTIFICATION  identified; coverage gap is observed and the placebo is constructible.
WORLDS          A  matching corrects coverage: change ∝ coverage gap, placebo change ~0
                B  matching subsamples: placebo change comparable to real change
KILL (CONDITIONAL) gate: the split with the largest coverage gap (pornhabit, 2.435) must show the
                   largest real-match change; otherwise the premise fails and no verdict follows.
                   then: placebo change < 1/3 of real change -> matching CORRECTS, A04's rule holds
                         placebo >= real                     -> it SUBSAMPLES, the rule is withdrawn
POSITIVE CTRL   pornhabit, coverage gap 2.435, where #11 measured a 62% correction.
NEGATIVE CTRL   modality, coverage gap 0.069 -- matching should change it least.
PLACEBO         match on a random person-level variable with the same number of strata, so the
                sample loss is comparable and the correction is nil.
SEEDS           5.
MULTIPLICITY    4 splits x 3 match types x 5 seeds, all reported.
IMPOSSIBLE      a split with a large coverage gap and no other difference -- coverage gap is not
                manipulable in an observational release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd, qr
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
neu=pd.to_numeric(df['neuroticismvariable'],errors='coerce').reindex(pool)
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
mod=df[MOD].reindex(pool)
def loadings(people,K=5):
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
SPL={'pornhabit [POS CTRL]':(pool[(ph>ph.median()).values],pool[(ph<ph.median()).values]),
     'sex':(pool[(male==1).values],pool[(male==0).values]),
     'neuroticism':(pool[(neu>neu.median()).values],pool[(neu<neu.median()).values]),
     'modality [NEG CTRL]':(pool[mod.isin(['Mostly written','Entirely written']).values],
                            pool[mod.isin(['Mostly visual','Entirely visual']).values])}
GAP={k:abs(nblk.reindex(pd.Index(sorted(set(a)))).mean()-nblk.reindex(pd.Index(sorted(set(b)))).mean())
     for k,(a,b) in SPL.items()}
def deficit(g1,g2,mode,rng):
    i1,i2=pd.Index(sorted(set(g1))),pd.Index(sorted(set(g2)))
    if mode=='none': a,b=i1.values,i2.values
    else:
        if mode=='coverage': k1,k2=nblk.reindex(i1).astype(int),nblk.reindex(i2).astype(int)
        else:
            rv=pd.Series(rng.integers(0,int(nblk.nunique()),len(pool)),index=pool)
            k1,k2=rv.reindex(i1),rv.reindex(i2)
        a=[];b=[]
        for v in set(k1)|set(k2):
            x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
            if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
        a,b=np.array(a),np.array(b)
    if min(len(a),len(b))<400: return np.nan,0
    c=cong(loadings(a),loadings(b))
    ceil=[cong(loadings(p[:len(a)]),loadings(p[len(a):len(a)+len(b)])) for p in [rng.permutation(pool) for _ in range(3)]]
    return float(np.nanmean(ceil))-c, len(a)
rows=[]
for name,(g1,g2) in SPL.items():
    for mode,seed in itertools.product(['none','coverage','placebo'],[1,2,3,4,5]):
        d,n=deficit(g1,g2,mode,np.random.default_rng(seed))
        if np.isfinite(d): rows.append(dict(split=name,mode=mode,seed=seed,deficit=d,n=n,gap=GAP[name]))
G=pd.DataFrame(rows); G.to_csv(OUT/'matching.csv',index=False)
P=G.groupby(['split','mode']).agg(deficit=('deficit','median'),n=('n','median')).reset_index()
W=P.pivot(index='split',columns='mode',values='deficit')
N=P.pivot(index='split',columns='mode',values='n')
W['coverage_gap']=[GAP[i] for i in W.index]
W['real_change']=(W['coverage']-W['none']).abs()
W['placebo_change']=(W['placebo']-W['none']).abs()
print("=== deficit by matching mode ===")
print(W.round(4).to_string())
print("\n=== people retained ===")
print(N.round(0).to_string())
pos=W.loc['pornhabit [POS CTRL]']
gate=bool(pos['real_change']==W['real_change'].max())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  largest coverage gap shows largest real-match change : {'PASS' if gate else 'FAIL'}")
if not gate: print("  -> gate FAILED : UNVERIFIED")
else:
    ratio=W['placebo_change'].median()/max(W['real_change'].median(),1e-9)
    print(f"  median placebo change / median real change = {ratio:.2f}")
    if ratio<1/3: print("  -> MATCHING CORRECTS : A04's rule holds, and the correction is not subsampling")
    elif ratio>=1: print("  -> IT SUBSAMPLES : A04's matching rule is withdrawn")
    else: print(f"  -> UNVERIFIED : placebo does {100*ratio:.0f}% of what the real match does")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
