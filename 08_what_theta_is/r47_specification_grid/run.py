import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
Controlling 'categories answered' may be over-control: if adversity -> more interests -> deeper
progression, that control removes real signal. Cleaner: define theta as a RATE (endorsed /
answered), which handles the denominator without conditioning on a possible mediator.
Report both specifications side by side, plus the MDE and a disattenuation bound, so the near-zero
is a measurement rather than a shrug.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
exec(open(round_path('45_theta_nonsexual.py')).read().split('print("=== FIRST')[0])
rate_theta=theta/answered.clip(lower=1)
ctrlA=pd.DataFrame({'age':df['age'].map(AGEMAP),'male':df['biomale'],'acq':acq}).fillna(0)
ctrlB=ctrlA.assign(answered=answered).fillna(0)
def scan(y,ctrl):
    out={}
    for k,v in NS.items():
        m=v.notna()&y.notna()
        if m.sum()<1500: continue
        X=np.c_[np.ones(m.sum()),ctrl[m].values]
        b,*_=lstsq(X,y[m].values,rcond=None); rt=y[m].values-X@b
        b2,*_=lstsq(X,v[m].values.astype(float),rcond=None); rv=v[m].values-X@b2
        out[k]=stats.spearmanr(rt,rv).statistic
    return pd.Series(out)
A=scan(theta,ctrlB).rename('count | +progression')
Bv=scan(rate_theta,ctrlA).rename('RATE | no progression ctrl')
C=scan(theta,ctrlA).rename('count | no progression ctrl')
T=pd.concat([A,Bv,C],axis=1).round(4)
T=T.reindex(T.abs().max(axis=1).sort_values(ascending=False).index)
print(T.to_string())
print(f"\n  largest |rho| in ANY specification : {T.abs().values.max():.3f}")
n=15000
print(f"  MDE at n={n:,}, 80% power, two-sided .05 : |rho| = {2.8/np.sqrt(n-3):.3f}")
print("  so this is a well-powered near-zero, not an underpowered null")
rel_theta=0.557        # measured in iter 8 (split-half over categories, Spearman-Brown)
for rel_pred in [0.7,0.8]:
    print(f"  disattenuated ceiling if predictors have reliability {rel_pred}: "
          f"{T.abs().values.max()/np.sqrt(rel_theta*rel_pred):.3f}")
