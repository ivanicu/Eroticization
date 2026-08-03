"""
The scan's top correlates were nearly all erotic items, i.e. theta correlating with itself.
Isolate the genuinely NON-SEXUAL side: personality, life history, childhood, mental health.
Two of those (TotalMentalIllness, childhood_adversity) were collapsed by the anonymisation to
a single value 'Any', so presence lives in MISSINGNESS -- and my nunique>=2 filter dropped them.
Recover them as binaries, but first check they are not just survey progression (iter 5's lesson).
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(514229)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
theta=(R>0).sum(1).astype(float)
answered=R.notna().sum(1).astype(float)          # survey progression
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
acq=df[lik].apply(pd.to_numeric,errors='coerce').mean(axis=1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
def code(c,order=None):
    s=df[c]
    if order: return s.map({v:i for i,v in enumerate(order)})
    return pd.to_numeric(s,errors='coerce')
NS={
 'openness':code('opennessvariable'),'conscientiousness':code('consciensiousnessvariable'),
 'extroversion':code('extroversionvariable'),'neuroticism':code('neuroticismvariable'),
 'agreeableness':code('agreeablenessvariable'),'powerlessness':code('powerlessnessvariable'),
 'has any mental illness':df['TotalMentalIllness'].notna().astype(float),
 'any childhood adversity':df['childhood_adversity'].notna().astype(float),
 'adult sexual assault victim':code([c for c in df.columns if 'victim of sexual assault' in c][0],['No','Yes']),
 'upbringing liberated':code([c for c in df.columns if 'sexually liberated' in c][0],['Repressed','Neutral','Liberated']),
 'spanked as a child':code([c for c in df.columns if 'were you spanked' in c][0],['Never','Sometimes','Often']),
 'childhood gender tolerance':code('childhood_gender_tolerance',['Intolerant','Medium','Tolerant']),
 'lifetime partner count':code('sexcount',['0','1-2','3-7','8-20','21+']),
 'monogamous':code([c for c in df.columns if 'preferred relationship style' in c][0],['Not monogamous','Monogamous']),
 'pornhabit':code('pornhabit'),
}
print("=== FIRST: are the missingness-coded variables just survey progression? ===")
for k in ['has any mental illness','any childhood adversity']:
    print(f"   corr({k:24s}, categories ANSWERED) = {stats.spearmanr(NS[k],answered).statistic:+.3f}")
print("   (a large value here means the variable is progression, not history)\n")

ctrl=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],'acq':acq,'answered':answered})
ctrl=ctrl.fillna(ctrl.median())
rows=[]
for k,v in NS.items():
    m=v.notna()&theta.notna()
    if m.sum()<1500: continue
    X=np.c_[np.ones(m.sum()),ctrl[m].values]
    b,*_=lstsq(X,theta[m].values,rcond=None); rt=theta[m].values-X@b
    b2,*_=lstsq(X,v[m].values.astype(float),rcond=None); rv=v[m].values-X@b2
    r=stats.spearmanr(rt,rv)
    rows.append(dict(variable=k,n=int(m.sum()),rho=round(r.statistic,4),p=f"{r.pvalue:.1e}",
                     var_pct=round(100*r.statistic**2,2)))
T=pd.DataFrame(rows).reindex(pd.DataFrame(rows).rho.abs().sort_values(ascending=False).index)
print("=== theta vs NON-SEXUAL variables (age, sex, agree-bias, survey progression controlled) ===")
print(T.to_string(index=False))
print(f"\n  largest |rho| among non-sexual variables : {T.rho.abs().max():.3f}"
      f"  = {100*T.rho.abs().max()**2:.1f}% of variance")
print(f"  total variance in theta explained by ALL {len(T)} together:")
m=theta.notna()
for k,v in NS.items(): m&=v.notna()
X=np.c_[np.ones(m.sum()),ctrl[m].values]+0.0
for k,v in NS.items(): X=np.c_[X,v[m].values.astype(float)]
b,*_=lstsq(X,theta[m].values,rcond=None); res=theta[m].values-X@b
Xc=np.c_[np.ones(m.sum()),ctrl[m].values]
bc,*_=lstsq(Xc,theta[m].values,rcond=None); resc=theta[m].values-Xc@bc
print(f"     R2 over and above age/sex/acq/progression = {1-res.var()/resc.var():.3f}   (n={m.sum():,})")
T.to_csv('data/derived/theta_nonsexual.csv',index=False)
