import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The same-gender-source set is only 7 pairs and the male ones are precum<->ejaculate, which are
nearly the same substance. "Same source gender" may just mean "same thing twice". Break it out:
squirt<->breastmilk are both female-origin but are NOT the same substance. If they also transfer,
source class is doing the work. If only precum<->ejaculate transfers, substance identity is.
Also add the neutral substances (saliva/urine/sweat), which have no source gender at all and were
excluded from the gendered comparison entirely.
"""
import pandas as pd, numpy as np, warnings, itertools
from numpy.linalg import lstsq
warnings.filterwarnings('ignore')
exec(open(round_path('50_additivity_origin.py')).read().split('for label,Dx in')[0])
def rr(Dx,k1,k2):
    a,b=Dx[k1],Dx[k2]; c=a.index.intersection(b.index)
    if len(c)<250: return np.nan,0
    x,y=a.reindex(c).values,b.reindex(c).values
    if x.std()==0 or y.std()==0: return np.nan,0
    return np.corrcoef(x,y)[0,1],len(c)
GROUPS={'MALE-origin, near-identical (precum/ejaculate)':[('precum','ejaculate')],
        'FEMALE-origin, distinct (squirt/breastmilk)'   :[('squirt','breastmilk')],
        'NEUTRAL pairs (saliva/urine/sweat)'            :[('saliva','urine'),('saliva','sweat'),('urine','sweat')],
        'CROSS source-gender'                           :[('precum','squirt'),('precum','breastmilk'),
                                                          ('ejaculate','squirt'),('ejaculate','breastmilk')],
        'NEUTRAL vs GENDERED'                           :[('saliva','precum'),('urine','ejaculate'),
                                                          ('sweat','squirt'),('saliva','squirt')]}
for label,Dx in [('RAW',D),('ORIENTATION PARTIALLED',Dr)]:
    print(f"\n=== {label}: same feature-contrast, by substance-pair type ===")
    for g,ps in GROUPS.items():
        vals=[]
        for s1,s2 in ps:
            for act in set(a for _,a in Dx):
                if (s1,act) in Dx and (s2,act) in Dx:
                    r,n=rr(Dx,(s1,act),(s2,act))
                    if not np.isnan(r): vals.append(r)
        if vals: print(f"   {g:48s} n={len(vals):2d}  mean r = {np.mean(vals):+.3f}")
print("\n  substance-identity account predicts: only precum/ejaculate high, everything else ~0")
print("  source-class account predicts:        squirt/breastmilk high too")
print("  pure additivity predicts:             ALL of them high, including cross and neutral")
