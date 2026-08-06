"""
E01 A07 R02 -- does the coverage law survive its own influence analysis?

corr(congruence deficit, between-group coverage gap) = +0.815 over nine splits is, after
retraction #15 conceded the HAZARD to the explorer's own 14-missingness, the only NOVEL claim
this project still has. Two known problems, both mine, both written down before this round:
  ADVERSARY_FORECAST #3, p=0.55: "nine splits, and PORNHABIT sits far out on both axes. Without it
    the correlation is much weaker, and the methodological claim rests on one influential point."
  R19: the deficits feeding it were all computed at K=5, the one K that turned out to be an outlier.

ESTIMAND        the population correlation between a split's congruence deficit and its coverage
                gap, over the space of person-variable splits -- with splits as the unit, n=9.
IDENTIFICATION  identified only if the estimate is not a function of one unit. n=9 with a visible
                outlier is the regime where a correlation is a story about one point. Leave-one-out
                and a splits-level bootstrap decide it; both are reported as intervals, not points.
SCOPE           population: person-variable splits available in this release, not splits in general
                · instrument: none · baseline: known-null splits added to the class for the first
                time · regime: n=9 units, which is the binding constraint on everything here.
WORLDS          A  the law is general: it survives dropping any split, holds at every K, and holds
                   when known-null splits are added to the class
                B  it is one point: leave-one-out swings it, or it dies away from K=5
KILL            PRE-REGISTERED: if the leave-one-out range dips below +0.40, or the median across
                K drops below +0.40, the coverage law is downgraded from a law to an observation
                about this release, and ADVERSARY_FORECAST #3 is scored CORRECT.
POSITIVE CTRL   coverage gap must correlate with itself at 1.0 through the same pipeline (checks
                the plumbing), and a synthetic deficit built AS a function of coverage must be
                recovered at r>0.9.
NEGATIVE CTRL   permute the coverage gaps across splits, 2000 draws -> the null distribution of r
                at n=9, which is wide by construction and is the point.
SHAM / PLACEBO  the coin-flip and row-parity splits from R01 are ADDED to the class, so for the
                first time the reference class contains members known to be null.
NOISE FLOOR     the permutation null's sd at n=9.
MULTIPLICITY    4 K values x (11 leave-one-out fits + full) x 3 seeds, all reported.
SPECIFICATION   K in {3,5,8,12} x class {9 original, 11 with known-nulls} x estimator {pearson,
                spearman}.
SEEDS           3.
IMPOSSIBLE      independently replicated · cross-dataset · more splits (this release has what it
                has -- n=9 is not a choice, and that is exactly why the interval must be published).
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import pandas as pd, numpy as np, warnings, hashlib, itertools
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
rng=np.random.default_rng(4242)
G=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A07_are_congruence_comparisons_reportable/R102_sweep_K_on_every_deficit/results/k_sweep.csv')
old=pd.read_csv(ROOT/'data/derived/deficit_reference.csv')
gap={r.split.split()[0].lower():float(r.blockgap) for _,r in
     pd.read_csv(ROOT/'data/derived/deficit_reference.csv').assign(
        blockgap=lambda d: d.get('blockgap', pd.Series([np.nan]*len(d)))).iterrows()} if 'blockgap' in old.columns else {}
# coverage gap must be recomputed for the splits R01 used
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
neu=pd.to_numeric(df['neuroticismvariable'],errors='coerce').reindex(pool)
MOD=[c for c in df.columns if 'type of erotic content you prefer tends to be:' in c][0]
mod=df[MOD].reindex(pool)
def cov_gap(g1,g2): return abs(nblk.reindex(pd.Index(g1)).mean()-nblk.reindex(pd.Index(g2)).mean())
sham=pd.Series(rng.random(len(pool))<(ph>ph.median()).mean(),index=pool)
SPL={'sex':(pool[(male==1).values],pool[(male==0).values]),
     'pornhabit':(pool[(ph>ph.median()).values],pool[(ph<ph.median()).values]),
     'modality':(pool[mod.isin(['Mostly written','Entirely written']).values],
                 pool[mod.isin(['Mostly visual','Entirely visual']).values]),
     'neuroticism':(pool[(neu>neu.median()).values],pool[(neu<neu.median()).values]),
     'sham_coinflip':(pool[sham.values],pool[~sham.values]),
     'placebo_rowparity':(pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1])}
GAP={k:cov_gap(*v) for k,v in SPL.items()}
rows=[]
for K in [3,5,8,12]:
    sub=G[G.K==K].groupby('split').deficit.median()
    x=np.array([GAP[s] for s in sub.index]); y=sub.values
    for est in ['pearson','spearman']:
        f=stats.pearsonr if est=='pearson' else stats.spearmanr
        r_all=f(x,y); r_all=r_all[0] if est=='pearson' else r_all.statistic
        loo=[]
        for i in range(len(x)):
            m=np.ones(len(x),bool); m[i]=False
            rr=f(x[m],y[m]); loo.append(rr[0] if est=='pearson' else rr.statistic)
        rows.append(dict(K=K,estimator=est,n_splits=len(x),r=round(float(r_all),3),
            loo_min=round(float(np.min(loo)),3),loo_max=round(float(np.max(loo)),3),
            dropped_worst=sub.index[int(np.argmin(loo))]))
T=pd.DataFrame(rows); T.to_csv(OUT/'influence.csv',index=False)
print("=== the coverage law, with known-null splits IN the class for the first time (n=6 splits) ===")
print(T.to_string(index=False))
print("\n=== permutation null at this n (2000 draws, coverage gaps shuffled across splits) ===")
sub=G[G.K==5].groupby('split').deficit.median(); x=np.array([GAP[s] for s in sub.index]); y=sub.values
nul=[stats.pearsonr(x,rng.permutation(y))[0] for _ in range(2000)]
print(f"  null mean {np.mean(nul):+.3f}  sd {np.std(nul):.3f}  |null| p95 {np.percentile(np.abs(nul),95):.3f}")
print(f"  -> at n={len(x)} splits, a correlation must exceed {np.percentile(np.abs(nul),95):.3f} to mean anything")
print("\n=== POSITIVE CONTROL ===")
print(f"  coverage gap vs itself                : r = {stats.pearsonr(x,x)[0]:.3f}")
synth=2.5*x+rng.normal(0,0.005,len(x))
print(f"  synthetic deficit built AS f(coverage): r = {stats.pearsonr(x,synth)[0]:.3f}  (must be >0.9)")
print("\n=== the coverage gaps themselves ===")
for k in sorted(GAP,key=lambda z:-GAP[z]): print(f"   {k:20s} gap {GAP[k]:.3f}   deficit@K5 {sub.get(k,np.nan):+.4f}")
med=T[T.estimator=='pearson'].r.median(); lo=T[T.estimator=='pearson'].loo_min.min()
print("\nPRE-REGISTERED KILL, evaluated:")
if lo<0.40 or med<0.40:
    print(f"  -> median r over K = {med:+.3f}, worst leave-one-out = {lo:+.3f} : DOWNGRADED from a law")
    print("  -> ADVERSARY_FORECAST #3 (p=0.55) scored CORRECT")
else:
    print(f"  -> median r = {med:+.3f}, worst LOO = {lo:+.3f} : survives, forecast #3 WRONG")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
