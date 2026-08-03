"""
E01 A01 R14 -- how much of "73% of the transfer is not demographic" is an erotic item removed from
an erotic item matrix?

The cross-domain CCA drops from 0.273 raw to 0.200 adjusted, and I reported that 73% of the
transfer survives demographic adjustment. But two of the eight covariates are erotic items:
    "I find it erotic when two people of the opposite gender to me sexually interact with each other"
    "You're more sexually attracted to people whose gender identity ____ match their genitals"
#26 found the acquiescence index survived audit because its erotic content CANCELLED across
opposite-signed items. These two have no cancelling twin, so the same audit should come out
differently, and the prediction is recorded here before running.

ESTIMAND        the cross-domain held-out CCA under a LADDER of covariate sets, and specifically
                the share of the total adjustment attributable to the two erotic covariates.
IDENTIFICATION  identified -- the covariate sets are nested and the decomposition is a difference
                of measured quantities, not an inference.
SCOPE           32 blocks · pairs with >=600 common people · K swept, because A07 R01 established
                this machinery is monotone in K and a single K is not reportable.
WORLDS          A  the adjustment is demographic: sex+age do the work, the erotic pair adds little
                B  the adjustment is partly self-inflicted: the erotic pair carries a large share,
                   and "73% not demographic" conflates 'not demographic' with 'not erotic'
KILL            PRE-REGISTERED: if the two erotic covariates account for more than 40% of the total
                raw->fully-adjusted drop, world B wins and the figure is republished as a
                decomposition rather than a single percentage. Under 20% the original stands.
POSITIVE CTRL   the ladder must be monotone -- each added covariate set can only reduce the CCA,
                since these are nested projections. A non-monotone rung means a bug, not a finding.
NEGATIVE CTRL   a covariate set of the same size built from random numeric columns; it must remove
                almost nothing.
SHAM            the two erotic covariates replaced by two non-erotic Likert items.
NOISE FLOOR     spread across 3 seeds at each rung.
MULTIPLICITY    5 covariate sets x 3 K x 3 seeds, all reported.
SEEDS           3.
IMPOSSIBLE      a release whose orientation is measured non-erotically -- orientation IS erotic
                preference, so this confound may be irreducible rather than merely unmeasured here.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; df['_age']=df['age'].map(AGEMAP)
ERO=[c for c in df.columns if 'opposite gender to me' in c or 'gender identity' in c.lower()]
PERS=['opennessvariable','neuroticismvariable','extroversionvariable','consciensiousnessvariable',
      'agreeablenessvariable','powerlessnessvariable']
rng0=np.random.default_rng(5)
numcols=[c for c in df.columns if df[c].dtype!=object and c not in ERO+PERS+['biomale','_age']]
RANDC=list(rng0.choice(numcols,8,replace=False))
lik=[c for c in pd.read_csv('data/derived/inventory.csv').query("kind=='LIKERT_PM3'")['col'] if c in df.columns]
NONERO=[c for c in lik if 'narcissist' in c.lower() or 'harrass' in c.lower()][:2]
SETS={'none':[], 'sex_age':['biomale','_age'], 'sex_age_pers':['biomale','_age']+PERS,
      'PUBLISHED (+2 erotic)':['biomale','_age']+PERS+ERO,
      'SHAM (+2 non-erotic)':['biomale','_age']+PERS+NONERO,
      'NEGCTRL (+8 random)':['biomale','_age']+PERS+RANDC}
def mk(cols):
    if not cols: return None
    X=df[cols].copy()
    for c in X.columns:
        if X[c].dtype==object: X[c]=X[c].astype('category').cat.codes.replace(-1,np.nan)
    X=X.apply(pd.to_numeric,errors='coerce'); X=X.fillna(X.median())
    return ((X-X.mean())/(X.std()+1e-9)).fillna(0.)
COVS={k:mk(v) for k,v in SETS.items()}
def factors(blocks,K,pool,COV):
    pm={p:i for i,p in enumerate(pool)}; cols=[]
    for q in blocks:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    if COV is not None:
        D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
        b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    U,S,Vt=svd(Z,full_matrices=False); return U[:,:K]*S[:K]
pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
rows=[]
for name,COV in COVS.items():
    for K,seed in itertools.product([3,5,8],[6,16,26]):
        rng=np.random.default_rng(seed)
        p=rng.permutation(allq); h1,h2=list(p[:len(p)//2]),list(p[len(p)//2:])
        F1,F2=factors(h1,K,pool,COV),factors(h2,K,pool,COV)
        idx=rng.permutation(len(pool)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
        c=CCA(n_components=K,max_iter=700).fit(F1[tr],F2[tr]); a,b_=c.transform(F1[te],F2[te])
        rows.append(dict(covset=name,K=K,seed=seed,cca=max(abs(np.corrcoef(a[:,j],b_[:,j])[0,1]) for j in range(K))))
G=pd.DataFrame(rows); G.to_csv(OUT/'covariate_ladder.csv',index=False)
P=G.pivot_table(index='covset',columns='K',values='cca',aggfunc='median')
order=['none','sex_age','sex_age_pers','PUBLISHED (+2 erotic)','SHAM (+2 non-erotic)','NEGCTRL (+8 random)']
print("=== held-out cross-domain CCA by covariate set and K ===")
print(P.reindex(order).round(4).to_string())
med=P.median(axis=1)
raw=med['none']; pers=med['sex_age_pers']; pub=med['PUBLISHED (+2 erotic)']
tot=raw-pub; ero=pers-pub
print(f"\n  raw (no covariates)            {raw:.4f}")
print(f"  after sex+age                  {med['sex_age']:.4f}")
print(f"  after +personality             {pers:.4f}")
print(f"  after +2 EROTIC (published)    {pub:.4f}")
print(f"  SHAM +2 non-erotic instead     {med['SHAM (+2 non-erotic)']:.4f}")
print(f"  NEGCTRL +8 random instead      {med['NEGCTRL (+8 random)']:.4f}")
print(f"\n  total adjustment raw->published : {tot:+.4f}")
print(f"  share carried by the 2 erotic items : {100*ero/tot:.0f}%")
mono=all(med[order[i]]>=med[order[i+1]]-1e-6 for i in range(3))
print(f"  POSITIVE CONTROL monotone across nested rungs: {'YES' if mono else 'NO -- bug, not finding'}")
print("\nPRE-REGISTERED KILL, evaluated:")
sh=100*ero/tot
if sh>40: print(f"  -> the 2 erotic covariates carry {sh:.0f}% (>40%) : republish as a decomposition, not a single percentage")
elif sh<20: print(f"  -> they carry {sh:.0f}% (<20%) : the original 73% figure stands")
else: print(f"  -> they carry {sh:.0f}% : between thresholds, UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
