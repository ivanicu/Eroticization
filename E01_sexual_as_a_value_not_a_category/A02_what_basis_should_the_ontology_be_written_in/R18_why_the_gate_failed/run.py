"""
E01 A02 R18 -- why the discriminant gate failed twice, and what would unfreeze it.

R16 and R17 both had their negative control return ~0.42-0.44 where it should return ~0, so both
pre-registered kills fired on an instrument that could not be trusted. Whitening the battery did
not help, which rules out the intercorrelation explanation I gave in R17. A freeze without a
diagnosis is abandonment with better manners, so: sample the noise control properly instead of once.

HYPOTHESIS  a profile is a vector of only 14 correlations, so the correlation BETWEEN two profiles
            has sampling sd about 1/sqrt(14-3) = 0.30 under the null. The single noise draw used in
            R16/R17 was not evidence of a biased instrument -- it was one draw from a very wide
            null, and I read it as a bias.
ESTIMAND    the null distribution of profile-profile correlation at this battery length.
KILL        PRE-REGISTERED: if the noise null has sd > 0.20, the design is declared underpowered by
            construction, the line is FROZEN, and the unfreeze condition is a battery long enough
            to bring that sd below 0.10. If sd < 0.20 the instrument is biased instead and the
            diagnosis is wrong.
CONTROLS    300 independent noise axes, each 4 series, through the identical pipeline.
IMPOSSIBLE  a longer battery -- this release has 14 usable non-role external variables.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
src=open(ROOT/'E01_sexual_as_a_value_not_a_category/A02_what_basis_should_the_ontology_be_written_in/R16_discriminant_validity/run.py').read()
body=src.split('def profile(')[0].split('"""')[2]
body='\n'.join(l for l in body.split('\n') if 'chdir' not in l and 'parents[3]' not in l and 'sys.path' not in l)
exec(body)
B=pd.DataFrame({k:pd.to_numeric(v,errors='coerce') for k,v in BATT.items()}); B=B.fillna(B.median())
Bz=(B-B.mean())/(B.std()+1e-9)
def prof(items):
    a=pd.concat(items,axis=1).mean(axis=1); m=a.notna()
    return np.array([float(stats.pearsonr(a[m].values,Bz.loc[m,k].values)[0]) for k in Bz.columns])
def pc(p1,p2):
    m=~(np.isnan(p1)|np.isnan(p2)); return float(np.corrcoef(p1[m],p2[m])[0,1])
Ppow=prof(AXES['POWER'])
L=len(Bz.columns)
print(f"profile length (battery size) = {L}")
d=np.array([pc(Ppow,prof([pd.Series(np.random.default_rng(s).normal(size=len(df)),index=df.index) for _ in range(4)])) for s in range(300)])
print(f"\nnoise-axis profile correlation with POWER, 300 INDEPENDENT noise draws")
print(f"  mean {d.mean():+.3f}   sd {d.std():.3f}   95% band [{np.percentile(d,2.5):+.3f}, {np.percentile(d,97.5):+.3f}]")
print(f"  the single draw R16/R17 relied on was +0.44 = the {100*(d<0.44).mean():.0f}th percentile")
print(f"  theoretical null sd at L={L}: {1/np.sqrt(L-3):.3f}")
print("\n  observed between-axis values against the properly sampled null:")
for k,v in {'POWER-GAZE':0.639,'POWER-SUBSTANCE':0.934,'GAZE-SUBSTANCE':0.825}.items():
    p=float((np.abs(d)>=abs(v)).mean())
    print(f"    {k:18s} {v:+.3f}   p={p:.3f}   {'ABOVE the noise band' if abs(v)>np.percentile(np.abs(d),95) else 'inside the noise band'}")
pd.Series(d,name='noise_profile_r').to_csv(OUT/'noise_null.csv',index=False)
print("\nPRE-REGISTERED KILL, evaluated:")
if d.std()>0.20:
    print(f"  -> noise null sd = {d.std():.3f} > 0.20 : UNDERPOWERED BY CONSTRUCTION, line FROZEN")
    print(f"  -> unfreeze condition: a battery of ~{int(1/0.10**2)+3} external variables (have {L})")
else:
    print(f"  -> noise null sd = {d.std():.3f} : the instrument is biased, not merely noisy; diagnosis wrong")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
