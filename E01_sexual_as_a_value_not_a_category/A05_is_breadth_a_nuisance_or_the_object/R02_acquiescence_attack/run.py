import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Iter 7 concluded: the induction self-report tracks breadth (rho .25) and nothing structural.
But breadth is 9.4% acquiescence measured on NON-erotic Likert items, and "yes, porn induced
new fetishes in me" is itself an agreeing answer. If agree-bias drives both, iter 7's headline
shrinks and part of it was response style.
Control it on items that have nothing to do with erotic categories, and report what survives.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
breadth=(R>0).sum(1).astype(float)
Lk=df[lik].apply(pd.to_numeric,errors='coerce')
acq=Lk.mean(axis=1); extr=Lk.abs().mean(axis=1)
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
g=df[IND].map(ORD)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
cov=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],
                  'porn':pd.to_numeric(df['pornhabit'],errors='coerce')})
cov=cov.fillna(cov.median())
print("=== does 'yes I was induced' itself ride on agree-bias? ===")
m=g.notna()&acq.notna()
print(f"  corr(induction answer, agree-bias)      = {stats.spearmanr(g[m],acq[m]).statistic:+.4f}  n={m.sum():,}")
print(f"  corr(induction answer, extreme-response)= {stats.spearmanr(g[m],extr[m]).statistic:+.4f}")
print(f"  corr(breadth, agree-bias)               = {stats.spearmanr(breadth[m],acq[m]).statistic:+.4f}")

def rho_after(controls,label):
    m=g.notna()&breadth.notna()
    for c in controls: m&=c.notna()
    X=np.c_[np.ones(m.sum()),cov[m].values]+0.0
    for c in controls: X=np.c_[X,c[m].values]
    b,*_=lstsq(X,breadth[m].values,rcond=None); rb=breadth[m].values-X@b
    b2,*_=lstsq(X,g[m].values.astype(float),rcond=None); rg=g[m].values-X@b2
    r=stats.spearmanr(rg,rb)
    print(f"  {label:44s} rho={r.statistic:+.4f}  p={r.pvalue:.2e}  n={m.sum():,}")
    return r.statistic
print("\n=== iter 7's headline, under increasing control ===")
r0=rho_after([],                 "breadth ~ induction | age,sex,porn        ")
r1=rho_after([acq],              "  + agree-bias                            ")
r2=rho_after([acq,extr],         "  + agree-bias + extreme-responding       ")
print(f"\n  survives response-style control: {100*r2/r0:.0f}% of the original rho")
print(f"  {'CONFIRMED, response style is a minority of it' if r2/r0>0.7 else 'DOWNGRADED -- response style carries much of it'}")

print("\n=== is agree-bias itself just 'has more of everything'? (reverse check) ===")
m=acq.notna()&breadth.notna()
X=np.c_[np.ones(m.sum()),breadth[m].values]
b,*_=lstsq(X,acq[m].values,rcond=None); ra=acq[m].values-X@b
print(f"  corr(induction, agree-bias | breadth)   = {stats.spearmanr(g[m].astype(float),ra,nan_policy='omit').statistic:+.4f}")
print("  -> if this stays large, agree-bias is a separate channel into the induction answer")
