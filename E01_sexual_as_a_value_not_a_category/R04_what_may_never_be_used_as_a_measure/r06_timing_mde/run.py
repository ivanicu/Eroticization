import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
rho=0.03 is only meaningful with an MDE. Both onset and porn-onset are binned in ~2-year buckets
and the gap is a difference of two binned variables, so attenuation could be severe.
Simulate: plant a REAL induction of X years in the 'totally different' group only, push it
through the same binning and the same missingness, and read off the rho we would have observed.
That converts a small rho into a bound in YEARS, which is the thing worth knowing.
"""
import pandas as pd, numpy as np, warnings
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(9871)
exec(open(round_path('33_induction_timing.py')).read().split('print("\\n=== mean interest-onset')[0])
EDGES=[0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,25.5,99]; MIDS=[2,5.5,7.5,9.5,11.5,13.5,15.5,17.5,22,28]
def binify(x): return np.array(MIDS)[np.clip(np.digitize(x,EDGES)-1,0,len(MIDS)-1)]
Ov=O.values; mask=~np.isnan(Ov); pv=porn.values
gv=g.values; okv=ok.values
obs_rho=stats.spearmanr(gv[okv], (pd.DataFrame(np.where(mask,Ov,np.nan)).sub(porn,axis=0).mean(axis=1))[okv],
                        nan_policy='omit').statistic
print(f"observed rho (raw mean gap, unadjusted) = {obs_rho:+.4f}")

def simulate(shift_years, frac_affected=0.35):
    """plant: in the 'totally different' group, a fraction of that person's interests actually
    arrive `shift_years` later than they really did."""
    S=Ov.copy()
    tgt=np.flatnonzero(okv & (gv==3))
    for i in tgt:
        j=np.flatnonzero(mask[i]); k_=max(1,int(round(frac_affected*len(j))))
        sel=rng.choice(j,k_,replace=False); S[i,sel]=S[i,sel]+shift_years
    S=np.where(mask,binify(S),np.nan)
    gap_=pd.DataFrame(S).sub(porn,axis=0).mean(axis=1)
    return stats.spearmanr(gv[okv],gap_[okv],nan_policy='omit').statistic
print("\nplanted induction -> rho we would have seen (mean of 5 draws):")
for s in [0.5,1.0,1.5,2.0,3.0,4.0]:
    r=np.mean([simulate(s) for _ in range(5)])
    print(f"   shift {s:.1f} yr in 35% of a 'totally different' person's interests -> rho = {r:+.4f}"
          f"{'   <-- observed sits here' if abs(r-obs_rho)<0.012 else ''}")
print("\nalso: concentration under a planted REAL induction (the discriminator that came out null)")
def conc_rho(shift_years,frac=0.35):
    S=Ov.copy()
    for i in np.flatnonzero(okv&(gv==3)):
        j=np.flatnonzero(mask[i]); k_=max(1,int(round(frac*len(j))))
        sel=rng.choice(j,k_,replace=False); S[i,sel]=S[i,sel]+shift_years
    S=np.where(mask,binify(S),np.nan)
    gp=pd.DataFrame(S).sub(porn,axis=0)
    c=gp.max(axis=1)-gp.mean(axis=1)
    return stats.spearmanr(gv[okv],c[okv],nan_policy='omit').statistic
for s in [1.0,2.0,3.0]:
    print(f"   shift {s:.1f} yr -> concentration rho = {np.mean([conc_rho(s) for _ in range(5)]):+.4f}")
print(f"   OBSERVED concentration rho = +0.0151 (p=0.15)")
