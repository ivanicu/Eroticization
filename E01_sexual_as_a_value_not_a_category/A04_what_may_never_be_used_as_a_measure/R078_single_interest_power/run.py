import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
The concentration null assumed induction touches ~35% of a person's interests. The realistic
case is ONE new fetish. Does the concentration measure (max gap minus own mean gap) still see a
single-interest induction? If yes, the observed null covers that case too. If no, the null is
blind to the very scenario the question is about and must be reported as UNVERIFIED.
"""
import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(6765)
exec(open(round_path('33_induction_timing.py')).read().split('print("\\n=== mean interest-onset')[0])
EDGES=[0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,25.5,99]; MIDS=[2,5.5,7.5,9.5,11.5,13.5,15.5,17.5,22,28]
def binify(x): return np.array(MIDS)[np.clip(np.digitize(x,EDGES)-1,0,len(MIDS)-1)]
Ov=O.values; mask=~np.isnan(Ov); gv=g.values; okv=ok.values
def sim(shift, n_affected):
    S=Ov.copy()
    for i in np.flatnonzero(okv&(gv==3)):
        j=np.flatnonzero(mask[i])
        sel=rng.choice(j,min(n_affected,len(j)),replace=False); S[i,sel]=S[i,sel]+shift
    S=np.where(mask,binify(S),np.nan)
    gp=pd.DataFrame(S).sub(porn,axis=0)
    c=(gp.max(axis=1)-gp.mean(axis=1))
    m=(gp.mean(axis=1))
    return (stats.spearmanr(gv[okv],c[okv],nan_policy='omit').statistic,
            stats.spearmanr(gv[okv],m[okv],nan_policy='omit').statistic)
print("planted induction in the 'New & totally different' group only")
print("  n_interests  shift(yr)   concentration_rho   mean_gap_rho")
for na in [1,2,3]:
    for sh in [1.0,2.0,3.0]:
        c,m=np.mean([sim(sh,na) for _ in range(5)],axis=0)
        print(f"      {na}          {sh:.1f}          {c:+.4f}            {m:+.4f}")
print(f"\n  OBSERVED                     +0.0151 (p=0.15)     +0.0614")
print("\n  power check: smallest planted case above vs observed concentration")
c1,_=np.mean([sim(1.0,1) for _ in range(8)],axis=0)
print(f"    one interest, 1 year -> concentration rho {c1:+.4f}  = {c1/0.0151:.1f}x the observed value")
print(f"    -> the concentration null {'DOES cover single-interest induction' if c1>0.05 else 'is BLIND to it; report UNVERIFIED'}")
