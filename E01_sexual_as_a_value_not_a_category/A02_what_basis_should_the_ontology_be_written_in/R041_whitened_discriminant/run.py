"""
E01 A02 R17 -- discriminant validity with a WHITENED battery, because R16's negative control failed.

R16 measured profile correlations against a 14-variable battery and its pre-registered kill fired
("101% of ceiling -> one construct"). The negative control says ignore it: a PURE NOISE axis scored
0.55-0.62 against every real axis. Diagnosis: the battery variables are intercorrelated, so every
profile vector is pulled toward the battery's first principal component and any two profiles agree
before content enters. Second time in this project a pre-registered kill has fired on an instrument
my own negative control invalidated (#21 was the first).

FIX: whiten the battery. Profiles are then correlations against DECORRELATED coordinates, so a
shared pull toward a dominant component cannot manufacture agreement.

ESTIMAND        profile similarity between role axes over a whitened external battery.
IDENTIFICATION  identified only if the noise control returns ~0 after whitening. That check comes
                FIRST and gates everything else; if it fails again the line is FROZEN, not re-fixed.
WORLDS          A three constructs -> between-axis profile r near the noise floor
                B one construct    -> near the split-half ceiling
KILL            PRE-REGISTERED: evaluated only if noise control |r| < 0.15 after whitening.
                Then: any pair above 80% of ceiling -> one construct. All below 50% -> distinct.
POSITIVE CTRL   split-half of the same axis must stay high after whitening.
NEGATIVE CTRL   the noise axis, which gates the whole round.
SHAM            non-role personality composite.
SEEDS           3.  MULTIPLICITY: 3 pairs x 2 partial settings x 3 seeds, all reported.
IMPOSSIBLE      independent replication; gold standard per axis.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, eigh
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
src=open(ROOT/'E01_sexual_as_a_value_not_a_category/A02_what_basis_should_the_ontology_be_written_in/R040_discriminant_validity/run.py').read()
exec(src.split('def profile(')[0].split('"""')[2])
B=pd.DataFrame({k:pd.to_numeric(v,errors='coerce') for k,v in BATT.items()})
B=B.fillna(B.median())
Bz=(B-B.mean())/(B.std()+1e-9)
C=np.corrcoef(Bz.values,rowvar=False)
w,V=eigh(C); w=np.clip(w,1e-3,None)
W=V@np.diag(1/np.sqrt(w))@V.T          # whitening matrix
BW=pd.DataFrame(Bz.values@W,columns=B.columns,index=B.index)
print(f"battery: {B.shape[1]} vars | mean |off-diag r| before whitening "
      f"{np.abs(C[np.triu_indices(len(C),1)]).mean():.3f}  after "
      f"{np.abs(np.corrcoef(BW.values,rowvar=False)[np.triu_indices(len(C),1)]).mean():.3f}")
def profile_w(items):
    a=pd.concat(items,axis=1).mean(axis=1)
    out=[]
    for k in BW.columns:
        m=a.notna()
        if m.sum()<300: out.append(np.nan); continue
        out.append(float(stats.pearsonr(a[m].values,BW.loc[m,k].values)[0]))
    return np.array(out)
def pc(p1,p2):
    m=~(np.isnan(p1)|np.isnan(p2)); return float(np.corrcoef(p1[m],p2[m])[0,1]) if m.sum()>5 else np.nan
P={k:profile_w(v) for k,v in AXES.items()}
ceil={}
for k,v in AXES.items():
    if len(v)<2: ceil[k]=np.nan; continue
    cs=[]
    for s in range(200):
        r=np.random.default_rng(s); idx=r.permutation(len(v)); h=max(1,len(idx)//2)
        cs.append(pc(profile_w([v[i] for i in idx[:h]]),profile_w([v[i] for i in idx[h:]])))
    ceil[k]=float(np.nanmedian(cs))
print("\n=== GATE: negative control after whitening (must be |r| < 0.15) ===")
noise=[pc(P[a],P['NOISE_ctrl']) for a in ['POWER','GAZE','SUBSTANCE']]
for a,v in zip(['POWER','GAZE','SUBSTANCE'],noise): print(f"   {a:10s} vs NOISE : {v:+.3f}")
gate=max(abs(v) for v in noise)<0.15
print(f"   max |r| = {max(abs(v) for v in noise):.3f}  ->  {'GATE PASSED' if gate else 'GATE FAILED, line FROZEN'}")
print("\n=== split-half ceilings after whitening (positive control) ===")
for k,v in ceil.items(): print(f"   {k:14s} {v:+.3f}")
rows=[]
for a,b in itertools.combinations(['POWER','GAZE','SUBSTANCE'],2):
    c=pc(P[a],P[b]); cl=np.nanmean([ceil[a],ceil[b]])
    rows.append(dict(pair=f"{a}-{b}",prof_r=round(c,3),ceiling=round(cl,3),
                     pct=round(100*c/cl,0) if cl and not np.isnan(cl) else np.nan))
for a in ['POWER','GAZE','SUBSTANCE']:
    rows.append(dict(pair=f"{a}-SHAM",prof_r=round(pc(P[a],P['SHAM_nonrole']),3),ceiling=np.nan,pct=np.nan))
T=pd.DataFrame(rows); T.to_csv(OUT/'whitened.csv',index=False)
print("\n=== between-axis profile similarity, whitened ===")
print(T.to_string(index=False))
print("\nPRE-REGISTERED KILL, evaluated:")
if not gate:
    print("  -> gate FAILED : no verdict. The line is FROZEN with its unfreeze condition recorded.")
else:
    mx=T[T.pct.notna()].pct.max()
    if mx>=80: print(f"  -> max {mx:.0f}% of ceiling : ONE CONSTRUCT, three-axis claim reduced")
    elif mx<50: print(f"  -> all pairs below 50% (max {mx:.0f}%) : DISTINCTNESS CONFIRMED without disattenuation")
    else: print(f"  -> max {mx:.0f}% : between thresholds, UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
