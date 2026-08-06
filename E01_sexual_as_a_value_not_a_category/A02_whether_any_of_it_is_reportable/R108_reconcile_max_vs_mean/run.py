"""
E01 A08 R02 -- reconcile #35: is 0.200 vs 0.149 a split-count effect or a different statistic?

#35 flagged the cross-domain CCA re-running at 0.1491 against a published 0.200, 3.3 seed spreads
apart, and I recorded it as an open discrepancy rather than guessing. Reading the code:

  A01 R09  cca_max = np.nanmax(cv)                  -> MAX over components
  A02 R01  cvcca(...) returns max(abs(corr))        -> MAX over components
  A08 R01  np.mean([abs(corr) for j in range(5)])   -> MEAN over components

So the two numbers are different statistics of the same fit. This round confirms that by computing
BOTH at every split count, which separates the two candidate causes rather than assuming one.

It also corrects a factual error in #27, where I wrote "this round reports the maximum canonical
correlation; the published 0.273/0.200 are the mean across components." Both published values are
the MAX. The note had it exactly backwards and has been sitting in the ledger since.

ESTIMAND        the held-out cross-domain canonical correlation under both aggregations, as a
                function of the number of half-splits averaged.
WORLDS          A  the gap is the statistic: max ~0.20 and mean ~0.15 at every split count
                B  the gap is averaging depth: both converge as splits increase
KILL (CONDITIONAL) gate: the two statistics must differ at EVERY split count (else they are not
                   distinguishable here and the comparison says nothing)
                   then: max stable within 15% across split counts AND mean stable within 15%
                         -> STATISTIC, #35 resolved and #27 corrected
                         either drifts more than 30% with split count -> AVERAGING DEPTH
POSITIVE CTRL   the max must be >= the mean in every cell, by construction. A cell violating that
                is a bug and invalidates the round.
NEGATIVE CTRL   a person-permuted fit, where both statistics must collapse toward the floor.
SEEDS           4 per split count.
MULTIPLICITY    6 split counts x 2 statistics x 4 seeds, all reported.
IMPOSSIBLE      nothing here -- this is a pure reconciliation and every quantity is computable.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, svd
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
def F(bl,K=5):
    pm={q:i for i,q in enumerate(pool)}; cols=[]
    for q in bl:
        idx=np.array([pm[x] for x in B[q]['ppl'] if x in pm]); src=np.array([i for i,x in enumerate(B[q]['ppl']) if x in pm])
        Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    U,S,_=svd(Z,full_matrices=False); return U[:,:K]*S[:K]
def one_split(rng,permute=False,K=5):
    p=rng.permutation(allq); h1,h2=list(p[:len(p)//2]),list(p[len(p)//2:])
    F1,F2=F(h1,K),F(h2,K)
    if permute: F2=F2[rng.permutation(len(pool))]
    idx=rng.permutation(len(pool)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    c=CCA(n_components=K,max_iter=600).fit(F1[tr],F2[tr]); a,b_=c.transform(F1[te],F2[te])
    cv=[abs(np.corrcoef(a[:,j],b_[:,j])[0,1]) for j in range(K)]
    return float(np.max(cv)),float(np.mean(cv))
rows=[]
for nsplit in [4,6,8,12,16,24]:
    for seed in (1,2,3,4):
        rng=np.random.default_rng(seed*1000+nsplit)
        vals=[one_split(rng) for _ in range(nsplit)]
        nul=[one_split(rng,permute=True) for _ in range(max(2,nsplit//3))]
        rows.append(dict(n_splits=nsplit,seed=seed,
                         cca_max=float(np.mean([v[0] for v in vals])),
                         cca_mean=float(np.mean([v[1] for v in vals])),
                         null_max=float(np.mean([v[0] for v in nul])),
                         null_mean=float(np.mean([v[1] for v in nul]))))
G=pd.DataFrame(rows); G.to_csv(OUT/'reconcile.csv',index=False)
P=G.groupby('n_splits')[['cca_max','cca_mean','null_max','null_mean']].median().round(4)
print("=== both statistics, by number of half-splits averaged ===")
print(P.to_string())
print(f"\n  published values: MAX-based 0.200 (A02 R01, demographics removed) · A08 audit used MEAN")
mx,mn=P['cca_max'],P['cca_mean']
gate_sep=bool((mx>mn).all()); gate_null=bool((P['null_max']<0.12).all())
dmx=(mx.max()-mx.min())/mx.median(); dmn=(mn.max()-mn.min())/mn.median()
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  max >= mean in every cell (construction) : {'PASS' if gate_sep else 'FAIL -- bug'}")
print(f"  permuted null collapses (max < 0.12)     : {'PASS' if gate_null else 'FAIL'} (max {P['null_max'].max():.3f})")
if not (gate_sep and gate_null): print("  -> gate FAILED : UNVERIFIED")
else:
    print(f"  drift with split count: max {100*dmx:+.0f}%   mean {100*dmn:+.0f}%")
    if dmx<0.15 and dmn<0.15:
        print(f"  -> THE GAP IS THE STATISTIC. max settles at {mx.median():.3f}, mean at {mn.median():.3f};")
        print("     split count changes neither. #35 resolved, #27's units note corrected.")
    elif max(dmx,dmn)>0.30: print("  -> AVERAGING DEPTH matters after all")
    else: print("  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
