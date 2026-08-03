"""
E01 A05 R10 -- was "theta has no external correlate" a reliability ceiling?

ADVERSARY_FORECAST #6, p=0.35: "breadth's 0.557 reliability makes theta's 'no external correlate'
partly a ceiling, and a better-measured theta would show the adversity correlations at 0.15 rather
than under 0.09." This became the most load-bearing unscored forecast when #17 showed one of those
correlates (childhood adversity) is real at r=+0.059 on an ungated outcome.

Classical test theory gives a sharp, falsifiable prediction rather than a hand-wave:
    r_observed = r_true * sqrt(reliability)
So build theta at a LADDER of reliabilities by compositing 1..K breadth indicators, measure the
reliability at each rung, and check whether r_observed / sqrt(rel) is FLAT.

ESTIMAND        r_true, the disattenuated association between each non-sexual variable and latent
                erotic breadth -- estimated by extrapolation along a measured reliability ladder
                rather than by dividing once by a single alpha.
IDENTIFICATION  identified if r_obs/sqrt(rel) is constant across rungs; if it drifts, the
                measurement model is wrong and only bounds are reportable.
SCOPE           BKS public n=15,503 · no instrument · baseline: known-null variables carried
                alongside · regime: reliability 0.3-0.9 achievable from this release's indicators.
WORLDS          A  CEILING: r_obs rises as sqrt(rel); r_obs/sqrt(rel) flat; extrapolates to ~0.15
                B  REALLY SMALL: r_obs flat as rel rises, so r_obs/sqrt(rel) DECLINES
                PREDICTION MATRIX          rung 1 (low rel)   rung K (high rel)
                  A  r_obs                  small              larger, by sqrt(rel) ratio
                  A  r_obs/sqrt(rel)        r_true             r_true (same)
                  B  r_obs                  small              same
                  B  r_obs/sqrt(rel)        higher             lower
KILL            PRE-REGISTERED: if r_obs/sqrt(rel) for childhood adversity declines by more than
                30% from the lowest to the highest reliability rung, world A is refuted and
                forecast #6 is scored WRONG. If flat within 30% AND the extrapolated r_true
                exceeds 0.12, #6 is scored CORRECT.
POSITIVE CTRL   a variable with a KNOWN large association with breadth must show the sqrt(rel)
                scaling clearly -- pornhabit, r about 0.22. If the ladder cannot reproduce the
                scaling law where the effect is large, it cannot test it where the effect is small.
NEGATIVE CTRL   two known-null variables (row-parity indicator, coin flip) carried on every rung;
                their r must stay ~0 at every reliability.
SHAM            a composite built from the same number of indicators drawn at random from
                NON-breadth columns -- its reliability ladder must not produce the scaling.
NOISE FLOOR     split-half sd of each rung's reliability estimate, 3 seeds.
MULTIPLICITY    6 rungs x 8 target variables x 3 seeds, all reported.
SPECIFICATION   rung size 1..6 indicators x seed x estimator {pearson, spearman}.
SEEDS           3.
IMPOSSIBLE      independent replication; a gold-standard breadth measure (none exists).
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
tot=[c for c in inv[inv['kind']=='TOTAL']['col'] if c in df.columns]
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
Rt=df[rate].apply(pd.to_numeric,errors='coerce')
def z(s): return (s-s.mean())/(s.std()+1e-9)
IND={ 'count_gated':(Rt>0).sum(1).astype(float),
      'rate_gated':((Rt>0).sum(1)/Rt.notna().sum(1).clip(lower=1)),
      'totalfetishcategory':pd.to_numeric(df['totalfetishcategory'],errors='coerce'),
      'answered_depth':Rt.notna().sum(1).astype(float),
      'total_cols_sum':df[tot].apply(pd.to_numeric,errors='coerce').sum(axis=1,min_count=1),
      'mean_intensity':Rt[Rt>0].mean(axis=1)}
IND={k:z(v.fillna(v.median())) for k,v in IND.items()}
names=list(IND)
acq=df[lik].apply(pd.to_numeric,errors='coerce').mean(axis=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
ctrl=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],'acq':acq})
ctrl=ctrl.fillna(ctrl.median())
rng0=np.random.default_rng(99)
TGT={'childhood_adversity':df['childhood_adversity'].notna().astype(float),
     'mental_illness':df['TotalMentalIllness'].notna().astype(float),
     'pornhabit_POSCTRL':pd.to_numeric(df['pornhabit'],errors='coerce'),
     'spanked':df[[c for c in df.columns if 'were you spanked' in c][0]].map({'Never':0,'Sometimes':1,'Often':2}),
     'assault':df[[c for c in df.columns if 'victim of sexual assault' in c][0]].map({'No':0,'Yes':1}),
     'openness':pd.to_numeric(df['opennessvariable'],errors='coerce'),
     'NULL_rowparity':pd.Series((np.arange(len(df))%2).astype(float),index=df.index),
     'NULL_coinflip':pd.Series((rng0.random(len(df))<0.5).astype(float),index=df.index)}
def rel_of(cols,seed):
    """split-half reliability of the composite over its own indicators, Spearman-Brown"""
    rng=np.random.default_rng(seed)
    if len(cols)<2: 
        # single indicator: reliability estimated by its correlation with the full pool mean
        full=pd.concat([IND[c] for c in names],axis=1).mean(axis=1)
        r=abs(np.corrcoef(IND[cols[0]],full)[0,1]); return float(r)
    p=rng.permutation(cols); h=len(p)//2
    a=pd.concat([IND[c] for c in p[:h]],axis=1).mean(axis=1)
    b=pd.concat([IND[c] for c in p[h:2*h]],axis=1).mean(axis=1)
    r=np.corrcoef(a,b)[0,1]; return float(2*r/(1+r))
def partial_r(theta,y):
    m=y.notna()&theta.notna()
    X=np.c_[np.ones(m.sum()),ctrl[m].values]
    rt=theta[m].values-X@lstsq(X,theta[m].values,rcond=None)[0]
    ry=y[m].values-X@lstsq(X,y[m].values,rcond=None)[0]
    return float(stats.pearsonr(rt,ry)[0])
rows=[]
for ksz,seed in itertools.product(range(1,7),[3,13,23]):
    rng=np.random.default_rng(seed*100+ksz)
    cols=list(rng.choice(names,ksz,replace=False))
    theta=pd.concat([IND[c] for c in cols],axis=1).mean(axis=1)
    rel=rel_of(cols,seed)
    for tname,y in TGT.items():
        rows.append(dict(k=ksz,seed=seed,rel=round(rel,3),target=tname,
                         r=round(partial_r(theta,y),4)))
G=pd.DataFrame(rows); G.to_csv(OUT/'ladder.csv',index=False)
S=G.groupby(['target','k']).agg(rel=('rel','median'),r=('r','median')).reset_index()
S['r_over_sqrt_rel']=S.r/np.sqrt(S.rel.clip(lower=0.05))
print("=== reliability ladder: does r_obs track sqrt(reliability)? ===")
print(S.pivot_table(index='k',columns='target',values='r').round(4).to_string())
print("\n=== measured reliability at each rung ===")
print(S.groupby('k').rel.median().round(3).to_string())
print("\n=== r_obs / sqrt(rel)  -- FLAT means a ceiling, DECLINING means the effect is really small ===")
print(S.pivot_table(index='k',columns='target',values='r_over_sqrt_rel').round(4).to_string())
print("\n=== POSITIVE CONTROL and NEGATIVE CONTROLS ===")
for t in ['pornhabit_POSCTRL','NULL_rowparity','NULL_coinflip']:
    s=S[S.target==t].sort_values('k')
    print(f"  {t:20s} r by rung: {[round(v,3) for v in s.r]}")
def drift(t):
    s=S[S.target==t].sort_values('rel'); 
    lo,hi=s.r_over_sqrt_rel.iloc[0],s.r_over_sqrt_rel.iloc[-1]
    return lo,hi,(hi-lo)/abs(lo) if lo!=0 else np.nan
print("\nPRE-REGISTERED KILL, evaluated on childhood_adversity:")
lo,hi,d=drift('childhood_adversity')
extrap=S[S.target=='childhood_adversity'].r_over_sqrt_rel.median()
print(f"  r/sqrt(rel) at lowest rel {lo:+.4f} -> at highest rel {hi:+.4f}   drift {100*d:+.0f}%")
print(f"  extrapolated r_true (median of r/sqrt(rel)) = {extrap:+.4f}")
if d<-0.30: print("  -> declines >30% : WORLD A REFUTED. Forecast #6 scored WRONG; the effect is really small.")
elif abs(d)<=0.30 and extrap>0.12: print("  -> flat and extrapolates above 0.12 : forecast #6 CORRECT, a ceiling was operating.")
else: print(f"  -> flat within 30% but extrapolates to {extrap:.3f} <= 0.12 : #6 WRONG on magnitude, ceiling small.")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
