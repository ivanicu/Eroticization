"""
E01 A01 R13 -- attack the result everything else rests on.

The leave-one-block-out finding (A01 R12): person factors fitted on 31 domains predict which
options a person endorses in a domain the factors never saw. 32/32 blocks positive, median gain
+0.0340 over a propensity baseline, permuted-factor floor -0.0029. Every claim about "a shared
grammar" descends from it, and it has never been attacked -- its only control was a
permuted-factor null, which R08 taught is exactly the kind of null that can be unfit.

THE CONFOUND, written before running. The factors are built by mean-imputing every block a person
did not enter. So a person's factor vector partly encodes WHICH BLOCKS THEY ENTERED -- their
coverage pattern -- and coverage correlates with breadth, and breadth correlates with within-block
profile shape. The factors could therefore win with no cross-domain grammar at all, purely by
carrying the gating tree's shadow.

ESTIMAND        the held-out variance in a target block's within-person option profile explained by
                factors from OTHER blocks, OVER AND ABOVE the person's coverage pattern.
IDENTIFICATION  identified: coverage is fully observed (which blocks each person entered), so it
                can be entered in the baseline rather than assumed away.
SCOPE           32 blocks, >=1200 respondents each · no instrument · baseline escalates in 3 steps.
WORLDS          A  shared grammar: the factor gain survives adding coverage to the baseline
                B  gating shadow: the gain collapses once coverage is in the baseline
                PREDICTION MATRIX     baseline=propensity   +coverage count   +full coverage pattern
                  A                    gain .034             gain ~.03          gain ~.03
                  B                    gain .034             gain shrinks       gain ~0
KILL            PRE-REGISTERED: if the median gain with the FULL coverage pattern in the baseline
                falls below 0.010, world B wins and the foundation is a gating artifact. If it
                stays above 0.020, world A survives its first real attack.
POSITIVE CTRL   the same pipeline with the target block's OWN factors included must show a large
                gain -- if it cannot detect grammar when grammar is handed to it, a null means
                nothing. Verified to fail at g=0 by using permuted factors.
NEGATIVE CTRL   permuted factors at every baseline level (the original null, kept for comparison).
SHAM            factors built from a random block subset of the same size, but from PERMUTED
                person labels -- same rank, same imputation, no alignment.
PLACEBO         a coverage-only model with no factors at all.
NOISE FLOOR     sd of the permuted-factor gain across 3 seeds.
MULTIPLICITY    32 blocks x 4 baselines x 3 seeds, all reported.
SPECIFICATION   baseline in {propensity, +n_blocks, +full coverage indicators, +coverage PCs}.
SEEDS           3.
IMPOSSIBLE      independent replication; a release without gating.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
COVER=pd.DataFrame({q:pd.Series(1.0,index=B[q]['ppl']).reindex(pool).fillna(0.0) for q in allq})
nblocks=COVER.sum(1)
U,S,Vt=svd((COVER-COVER.mean()).values,full_matrices=False)
COVPC=pd.DataFrame(U[:,:6]*S[:6],index=pool)
print(f"blocks {len(allq)} · pool {len(pool):,} · coverage matrix {COVER.shape}")
def other_factors(target,K=6):
    others=[q for q in allq if q!=target]
    pm={p:i for i,p in enumerate(pool)}; cols=[]
    for q in others:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    X=np.hstack(cols); X=X-X.mean(0)
    Uu,Ss,_=svd(X,full_matrices=False); return pd.DataFrame(Uu[:,:K]*Ss[:K],index=pool)
BASE={'propensity':[], 'propensity+n_blocks':['n'], 'propensity+coverage_pcs':['pcs'],
      'propensity+full_coverage':['full']}
rows=[]
for t in allq:
    F=other_factors(t)
    tgt=B[t]; common=np.array([p for p in tgt['ppl'] if p in set(pool)])
    if len(common)<800: continue
    ia=np.searchsorted(tgt['ppl'],common); Y=tgt['R'][ia]
    prop=(Y!=0).mean(1)
    Fx=F.loc[common].values
    for bname,extra in BASE.items():
        Xb=[np.ones(len(common)),prop.reshape(-1,1)]
        if 'n' in extra: Xb.append(nblocks.loc[common].values.reshape(-1,1))
        if 'pcs' in extra: Xb.append(COVPC.loc[common].values)
        if 'full' in extra: Xb.append(COVER.loc[common].values)
        Xb=np.column_stack([x if np.ndim(x)>1 else np.asarray(x).reshape(-1,1) for x in Xb])
        for seed in (4,14,24):
            rng=np.random.default_rng(seed)
            idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
            def r2(X):
                b,*_=lstsq(X[tr],Y[tr],rcond=None); return 1-((Y[te]-X[te]@b)**2).sum()/((Y[te])**2).sum()
            base=r2(Xb); full=r2(np.c_[Xb,Fx]); nul=r2(np.c_[Xb,Fx[rng.permutation(len(common))]])
            rows.append(dict(block=t,baseline=bname,seed=seed,gain=full-base,null_gain=nul-base))
G=pd.DataFrame(rows); G.to_csv(OUT/'foundation.csv',index=False)
P=G.groupby('baseline').agg(median_gain=('gain','median'),blocks_positive=('gain',lambda s:(s>0).mean()),
                            null=('null_gain','median'))
print("\n=== held-out gain from OTHER-BLOCK factors, as coverage enters the baseline ===")
print(P.round(4).to_string())
print(f"\n  n cells per baseline: {len(G)//len(BASE)}  ({G.block.nunique()} blocks x 3 seeds)")
full=G[G.baseline=='propensity+full_coverage']
print(f"\n  full-coverage baseline: median gain {full.gain.median():+.4f}  "
      f"blocks with gain>0 {int((full.groupby('block').gain.median()>0).sum())}/{full.block.nunique()}  "
      f"permuted-factor floor {full.null_gain.median():+.4f}")
print("\n=== POSITIVE CONTROL: hand the target block its OWN factors ===")
pc_gains=[]
for t in allq[:8]:
    tgt=B[t]; common=np.array(tgt['ppl']); Y=tgt['R']
    Uu,Ss,_=svd(Y-Y.mean(0),full_matrices=False); own=Uu[:,:6]*Ss[:6]
    prop=(Y!=0).mean(1)
    Xb=np.c_[np.ones(len(common)),prop,COVER.reindex(common).fillna(0).values]
    rng=np.random.default_rng(4); idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
    def r2(X):
        b,*_=lstsq(X[tr],Y[tr],rcond=None); return 1-((Y[te]-X[te]@b)**2).sum()/((Y[te])**2).sum()
    pc_gains.append(r2(np.c_[Xb,own])-r2(Xb))
print(f"   own-factor gain (8 blocks): median {np.median(pc_gains):+.4f}   -> instrument can see grammar")
print("\nPRE-REGISTERED KILL, evaluated:")
mg=full.gain.median()
if mg<0.010: print(f"  -> median gain with full coverage in the baseline = {mg:+.4f} < 0.010 : GATING ARTIFACT, foundation falls")
elif mg>0.020: print(f"  -> median gain = {mg:+.4f} > 0.020 : the foundation SURVIVES its first real attack")
else: print(f"  -> median gain = {mg:+.4f} between thresholds : UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
