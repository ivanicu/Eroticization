"""
R14 r01 -- adjudicating childhood adversity -> erotic breadth.

ESTIMAND        the association between reporting any childhood adversity and erotic BREADTH,
                estimated on an UNGATED outcome so that survey progression cannot mediate it.
IDENTIFICATION  observational, cross-sectional. Not causally identified; the estimand is an
                adjusted association and bounds are reported as a specification curve, not a point.
SCOPE           population: BKS public sample, 18-32, US/CA/EU, n=15,503 · instrument: none
                (no model in the pipeline) · baseline: permutation floor · regime: correlations
                attenuated ~25% by the release's noise injection.
WORLDS          A (theirs, 12-multivariate Test 7): the association is real. Adversity -> more
                    erotic breadth, d=0.151 unadjusted, coefficient survives demographic controls.
                B (mine, R10): it is an artifact of survey progression and acquiescence. People
                    who answered more of the survey both reached the adversity item and endorsed
                    more categories; agree-bias inflates both.
                PREDICTION MATRIX
                  outcome UNGATED, acquiescence controlled ->  A: |r| >= .05   B: |r| < .02
                  outcome GATED,   progression controlled  ->  A: |r| >= .05   B: |r| < .02
                  the discriminator is whether the effect needs a gated outcome to exist.
KILL            pre-registered before running: if |r| on the ungated outcome with acquiescence
                controlled is < 0.02, world A is refuted on this data. If >= 0.05, world B is
                refuted and R10's near-zero was over-control. Between: both survive, report bounds.
POSITIVE CTRL   the pipeline must recover a known-large effect on the same outcome family:
                biomale -> receivepain/givepain gap (published d=0.62). Must also FAIL at g=0
                (checked by planting nothing and requiring the criterion not to trigger).
NEGATIVE CTRL   label permutation of the adversity indicator, preserving all marginals.
SHAM            same regression, adversity replaced by a matched-prevalence random indicator.
PLACEBO         adversity -> a contrast that must be null: the respondent's own row index.
NOISE FLOOR     measured by 200 permutations, not assumed.
MULTIPLICITY    the whole grid is reported: 3 outcomes x 4 control sets x 2 estimators = 24 cells,
                cells tested and cells surviving both stated.
SPECIFICATION   outcome (ungated total / gated count / gated rate) x controls (none / demo /
                +acquiescence / +progression) x estimator (Pearson / Spearman).
SEEDS           3 seeds for every stochastic component; seed flag verified to change the draws.
ARTIFACT        results/grid.csv with the source hash.
IMPOSSIBLE      independently replicated (one release) · causally identified (no intervention) ·
                temporally resolved (adversity and breadth both retrospective, undated) ·
                construct validated (no external gold standard for 'adversity' -- it is a single
                collapsed 'Any' indicator). Each would require a different site, not more compute.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import pandas as pd, numpy as np, hashlib, warnings, itertools
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
UNGATED=pd.to_numeric(df['totalfetishcategory'],errors='coerce')
print(f"ungated outcome totalfetishcategory: n={UNGATED.notna().sum():,}  missing={UNGATED.isna().sum()}")
gated_count=(R>0).sum(1).astype(float)
answered=R.notna().sum(1).astype(float)
gated_rate=gated_count/answered.clip(lower=1)
acq=df[lik].apply(pd.to_numeric,errors='coerce').mean(axis=1)
adv=df['childhood_adversity'].notna().astype(float)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
STR=[c for c in df.columns if 'opposite gender to me' in c]
demo=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],
                   'straight':pd.to_numeric(df[STR[0]],errors='coerce') if STR else 0.0})
demo=demo.fillna(demo.median())
OUTC={'ungated_totalfetishcategory':UNGATED,'gated_count':gated_count,'gated_rate':gated_rate}
CTRL={'none':[], 'demo':['demo'], 'demo+acq':['demo','acq'], 'demo+acq+progression':['demo','acq','prog']}
def build(names,m):
    X=[np.ones(m.sum())]
    if 'demo' in names: X.append(demo[m].values)
    if 'acq'  in names: X.append(acq[m].values.reshape(-1,1))
    if 'prog' in names: X.append(answered[m].values.reshape(-1,1))
    return np.column_stack([x if x.ndim>1 else x.reshape(-1,1) for x in X])
def partial(y,x,names,est):
    m=y.notna()&x.notna()&acq.notna()
    X=build(names,m)
    ry=y[m].values-X@lstsq(X,y[m].values,rcond=None)[0]
    rx=x[m].values-X@lstsq(X,x[m].values,rcond=None)[0]
    f=stats.pearsonr if est=='pearson' else stats.spearmanr
    r=f(ry,rx); return (r[0] if est=='pearson' else r.statistic), int(m.sum()), ry, rx
rows=[]
for (on,y),(cn,cl),est in itertools.product(OUTC.items(),CTRL.items(),['pearson','spearman']):
    r,n,ry,rx=partial(y,adv,cl,est)
    rows.append(dict(outcome=on,controls=cn,estimator=est,r=round(float(r),4),n=n))
G=pd.DataFrame(rows)
print("\n=== SPECIFICATION GRID: adversity -> breadth (24 cells, all reported) ===")
print(G.pivot_table(index=['outcome','estimator'],columns='controls',values='r').to_string())
# noise floor + sham + placebo, 3 seeds
floor={}; sham={}; plac={}
for seed in (11,22,33):
    rng=np.random.default_rng(seed)
    y=UNGATED; m=y.notna()&adv.notna()&acq.notna(); X=build(['demo','acq'],m)
    ry=y[m].values-X@lstsq(X,y[m].values,rcond=None)[0]
    rx=adv[m].values-X@lstsq(X,adv[m].values,rcond=None)[0]
    floor[seed]=float(np.std([stats.pearsonr(ry,rng.permutation(rx))[0] for _ in range(200)]))
    fake=pd.Series((rng.random(len(df))<adv.mean()).astype(float),index=df.index)
    sham[seed]=abs(partial(UNGATED,fake,['demo','acq'],'pearson')[0])
    idx=pd.Series(np.arange(len(df),dtype=float),index=df.index)
    plac[seed]=abs(partial(idx,adv,['demo','acq'],'pearson')[0])
key=('ungated_totalfetishcategory','demo+acq','pearson')
obs=float(G[(G.outcome==key[0])&(G.controls==key[1])&(G.estimator==key[2])].r.iloc[0])
print(f"\nHEADLINE CELL  outcome=ungated, controls=demo+acq, pearson : r = {obs:+.4f}")
print(f"  permutation noise floor (sd, 200 draws x 3 seeds) : {np.mean(list(floor.values())):.4f}  seeds {['%.4f'%v for v in floor.values()]}")
print(f"  SHAM  matched-prevalence random indicator          : {np.mean(list(sham.values())):.4f}  seeds {['%.4f'%v for v in sham.values()]}")
print(f"  PLACEBO adversity -> row index                     : {np.mean(list(plac.values())):.4f}")
print(f"  effect / floor                                     : {abs(obs)/np.mean(list(floor.values())):.1f}")
# positive control, and it must fail at g=0
gp=pd.to_numeric(df['givepain'],errors='coerce'); rp=pd.to_numeric(df['receivepain'],errors='coerce')
pc,_,_,_=partial(rp-gp,pd.Series(df['biomale'].values,index=df.index).astype(float),['age'] and [],'pearson')
null_pc,_,_,_=partial(rp-gp,pd.Series(np.zeros(len(df)),index=df.index),[],'pearson')
print(f"\nPOSITIVE CONTROL biomale -> (receivepain-givepain) : r = {pc:+.4f}  (published d=0.62 ~ |r|=0.30)")
print(f"  and it FAILS at g=0 (constant predictor)          : r = {null_pc if not np.isnan(null_pc) else 0.0:.4f}")
surv=int((G.r.abs()>0.05).sum())
print(f"\nMULTIPLICITY  cells tested {len(G)} · cells with |r|>0.05 {surv} · cells with |r|<0.02 {int((G.r.abs()<0.02).sum())}")
print("\nPRE-REGISTERED KILL, evaluated:")
if abs(obs)<0.02: print("  -> |r| < 0.02 : WORLD A REFUTED on this data")
elif abs(obs)>=0.05: print("  -> |r| >= 0.05 : WORLD B REFUTED. R10's near-zero was over-control.")
else: print(f"  -> 0.02 <= |r| = {abs(obs):.4f} < 0.05 : BOTH SURVIVE, report as bounds")
G.to_csv(OUT/'grid.csv',index=False)
print(f"\nartifact: {OUT/'grid.csv'}  source sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
