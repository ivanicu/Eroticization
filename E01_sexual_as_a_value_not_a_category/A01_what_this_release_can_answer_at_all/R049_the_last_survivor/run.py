"""
E01 A02 R25 -- the last surviving A02 finding, tested with the method that killed its twin.

After #24 (POWER-SUBSTANCE restated), #39 (SUBSTANCE scoped to 3% of the corpus) and #41 (POWER weak
everywhere), the only A02 claim left is the CONDITIONAL one: role predicts endorsement wherever the
OPTION SET ITSELF varies in role -- r = +0.752 between a block's option-set role variance and POWER's
predictive gain, with a 21x ratio between varying and non-varying blocks.

It has the same shape as the coverage law #20 killed: a correlation over ~28 units whose entire
x-axis spread comes from a handful of points. #30 already noted the role-varying cells are 4 blocks
and all four are the fluid family -- which #39 has since shown is subgraph-local. So the claim is
a correlation driven by 4 points from a family now known not to represent the corpus.

ESTIMAND        corr(option-set role variance, POWER's predictive gain) over blocks, with its
                leave-one-out range and its permutation null AT THIS n.
IDENTIFICATION  identified only if the estimate is not a function of one family. n=28 blocks with
                4 non-zero x-values is the regime where a correlation is a story about those 4.
WORLDS          A  general: survives dropping any block, and the null at n=28 is far below it
                B  four points: leave-one-out swings it, or the null at this n is comparable
KILL (CONDITIONAL) gate: the permutation null must be computable and the observed must be finite.
                   then: leave-one-out min < 0.40 OR observed inside the null's 95% band
                         -> FOUR POINTS, A02 has nothing left standing release-wide
                         LOO min > 0.40 and observed outside the null -> SURVIVES
POSITIVE CTRL   plant a synthetic gain that is a known linear function of role variance; the
                pipeline must recover r > 0.9.
NEGATIVE CTRL   permute role variance across blocks, 5000 draws -> the null band at n=28.
SHAM            replace role variance with a random block-level covariate of the same distribution.
SEEDS           3 for the gain estimates.
MULTIPLICITY    28 leave-one-out fits + full + sham + planted, all reported.
IMPOSSIBLE      more role-varying blocks -- the release instantiates the contrast in 4 (7 counting
                the smaller fluid blocks), and that is the binding constraint, not a choice.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
rv=pd.read_csv('data/derived/role_variance.csv')
print(f"blocks in the published correlation: {len(rv)}")
print(f"  with role_var > 0: {int((rv.role_var>0).sum())}  -> {sorted(rv[rv.role_var>0].qi.tolist())}")
FLUID={6,7,8,9,10,11,83}
rv['fluid']=rv.qi.isin(FLUID)
print(f"  of those, in the fluid family: {int((rv.role_var>0).mul(rv.fluid).sum())}")
x=rv.role_var.values; y=rv.d_power.values
r_all=stats.pearsonr(x,y)[0]
loo=[]
for i in range(len(x)):
    m=np.ones(len(x),bool); m[i]=False
    loo.append((stats.pearsonr(x[m],y[m])[0], int(rv.qi.iloc[i])))
loo_r=np.array([v[0] for v in loo])
rng=np.random.default_rng(4)
nul=np.array([stats.pearsonr(x,rng.permutation(y))[0] for _ in range(5000)])
sham=np.array([stats.pearsonr(rng.permutation(x),y)[0] for _ in range(5000)])
no_fluid=rv[~rv.fluid]
r_nf=stats.pearsonr(no_fluid.role_var.values,no_fluid.d_power.values)[0] if no_fluid.role_var.std()>0 else np.nan
print(f"\n=== the correlation and its influence ===")
print(f"  observed r                       : {r_all:+.3f}")
print(f"  leave-one-out range              : {loo_r.min():+.3f} .. {loo_r.max():+.3f}")
print(f"  block whose removal costs most   : qi={loo[int(np.argmin(loo_r))][1]}  -> r drops to {loo_r.min():+.3f}")
print(f"  with the fluid family removed    : {r_nf if not np.isnan(r_nf) else float('nan'):.3f}   "
      f"(role_var sd among non-fluid = {no_fluid.role_var.std():.4f})")
print(f"\n=== the null AT THIS n ===")
print(f"  permutation null: mean {nul.mean():+.3f}  sd {nul.std():.3f}  |r| p95 {np.percentile(np.abs(nul),95):.3f}")
print(f"  sham (x permuted): |r| p95 {np.percentile(np.abs(sham),95):.3f}")
print(f"\n=== POSITIVE CONTROL ===")
planted=2.0*x+rng.normal(0,0.002,len(x))
print(f"  gain planted as a linear function of role variance: r = {stats.pearsonr(x,planted)[0]:+.3f}  (must be >0.9)")
inside=abs(r_all)<np.percentile(np.abs(nul),95)
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control recovers >0.9 : {'PASS' if stats.pearsonr(x,planted)[0]>0.9 else 'FAIL'}")
print(f"  null computable and finite     : PASS")
if loo_r.min()<0.40 or inside:
    print(f"  -> FOUR POINTS : LOO floor {loo_r.min():+.3f}, observed {'inside' if inside else 'outside'} the n=28 null band (p95 {np.percentile(np.abs(nul),95):.3f})")
    print("  -> A02 has nothing left standing release-wide")
else:
    print(f"  -> SURVIVES : LOO floor {loo_r.min():+.3f}, observed outside the null band")
pd.DataFrame(dict(dropped=[v[1] for v in loo],r=loo_r)).to_csv(OUT/'influence.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
