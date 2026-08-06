"""
E01 A01 R17 -- the central claim under a deliberately different framing.

Forty-eight ledger entries, all self-inflicted, and ADVERSARY_FORECAST has never been scored
against an outside challenger. This session forbids dispatching one, so those rows stay
[unchallenged] -- not "clean" (frontier §5). Self-review is VOID, not weak: a reviewer sampled from
the weights that produced the material can only attack what was already anticipated.

The available substitute is not another reviewer but another FRAMING. realstat §2.5: three
independent designs test the FRAMING, where a re-implementation only tests the CODE, and most
retractions in this project were framing errors that correct code would have preserved.

  published framing : held-out CANONICAL CORRELATION between two blocks' residual matrices,
                      max over 3 components, median across 306 pairs. 0.269 raw / 0.198 adjusted.
  this framing      : held-out PREDICTION. Fit block A's residuals -> block B's residuals on 70% of
                      the shared people, score R2 on the held-out 30%. Different estimand family
                      (prediction vs correlation), different aggregation, different failure modes.

Agreement on sign and rough magnitude means the framing is robust. Disagreement means the framing
IS the finding, and must not be averaged away.

ESTIMAND        held-out cross-block predictive R2, per pair, median across pairs.
IDENTIFICATION  identified; the prediction is out-of-sample by construction.
WORLDS          A  transfer is real: prediction R2 clearly above its permuted-person null
                B  the CCA result was a property of canonical correlation: prediction finds nothing
KILL (CONDITIONAL) gate: the graded tiers from #38 must reproduce under THIS framing too --
                   precum<->ejaculate above within-family above all-other. If the ordering breaks,
                   the new framing is not measuring the same thing and no comparison follows.
                   then: median R2 > 3x its null -> the claim survives a framing change
                         R2 at the null           -> framing-dependent, and that is the finding
POSITIVE CTRL   the #38 tier ordering, reproduced under prediction rather than correlation.
NEGATIVE CTRL   permuted persons in block B before fitting.
SEEDS           3.
MULTIPLICITY    306 pairs x 3 seeds, tiers reported separately.
IMPOSSIBLE      genuine independence -- I designed both framings. This is labelled a framing check,
                NOT a replication, and the ledger rows stay [unchallenged].
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); FLUID={7,8,9,11,83,6,10}
def strip(M,idx):
    D=np.c_[np.ones(len(idx)),COV.loc[idx].values]
    b,*_=lstsq(D,M,rcond=None); return M-D@b
def predict_r2(Xa,Xb,rng,permute=False):
    n=len(Xa); idx=rng.permutation(n); tr,te=idx[:int(.7*n)],idx[int(.7*n):]
    Y=Xb[rng.permutation(n)] if permute else Xb
    A=np.c_[np.ones(len(tr)),Xa[tr]]
    b,*_=lstsq(A,Y[tr],rcond=None)
    P=np.c_[np.ones(len(te)),Xa[te]]@b
    return 1-((Y[te]-P)**2).sum()/max((Y[te]**2).sum(),1e-12)
def tier(a,b):
    if {a,b}=={7,83}: return 'tier1 precum<->ejaculate'
    if a in FLUID and b in FLUID: return 'tier2 within fluid family'
    if (a in FLUID)!=(b in FLUID): return 'tier3 fluid<->non-fluid'
    return 'tier4 all other pairs'
rows=[]
for seed in (1,2,3):
    rng=np.random.default_rng(seed)
    for a,b in itertools.combinations(allq,2):
        A_,B_=B[a],B[b]; common=np.intersect1d(A_['ppl'],B_['ppl'])
        if len(common)<600: continue
        ia=np.searchsorted(A_['ppl'],common); ib=np.searchsorted(B_['ppl'],common)
        Ra,Rb=A_['R'][ia],B_['R'][ib]
        rows.append(dict(seed=seed,a=a,b=b,tier=tier(a,b),
                         raw=predict_r2(Ra,Rb,rng),
                         adj=predict_r2(strip(Ra,common),strip(Rb,common),rng),
                         null=predict_r2(Ra,Rb,rng,permute=True)))
G=pd.DataFrame(rows); G.to_csv(OUT/'prediction_framing.csv',index=False)
print(f"pairs per seed: {G[G.seed==1].shape[0]}")
print("\n=== cross-block held-out predictive R2, by tier (the #38 ordering, under PREDICTION) ===")
T=G.groupby('tier')[['raw','adj','null']].median().round(4)
print(T.to_string())
med=G[['raw','adj','null']].median()
print(f"\n=== overall medians ===")
print(f"  raw prediction R2      {med['raw']:+.4f}   (CCA framing: 0.269)")
print(f"  demographics-removed   {med['adj']:+.4f}   (CCA framing: 0.198)")
print(f"  permuted-person null   {med['null']:+.4f}")
order=['tier1 precum<->ejaculate','tier2 within fluid family','tier4 all other pairs']
have=[t for t in order if t in T.index]; vals=[float(T.loc[t,'raw']) for t in have]
gate=all(vals[i]>vals[i+1] for i in range(len(vals)-1))
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  #38's tier ordering reproduces under prediction : {'PASS' if gate else 'FAIL'} ({[round(v,3) for v in vals]})")
if not gate: print("  -> gate FAILED : the new framing is not measuring the same thing, UNVERIFIED")
else:
    ratio=med['adj']/max(abs(med['null']),1e-9)
    print(f"  adjusted / |null| = {ratio:.1f}")
    if med['adj']>3*abs(med['null']): print("  -> SURVIVES A FRAMING CHANGE : prediction and correlation agree that transfer is real")
    else: print("  -> FRAMING-DEPENDENT : and that is the finding, not something to average away")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
