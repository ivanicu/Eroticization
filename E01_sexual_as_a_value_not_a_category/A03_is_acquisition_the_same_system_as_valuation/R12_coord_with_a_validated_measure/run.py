"""
E01 A03 R12 -- was COORD's null an absence, or a broken measure?

A03's decision -- model acquisition and valuation as TWO SYSTEMS -- rests on #9: with preference
similarity controlled, TEMPO (population arrival-time distance) predicts within-person onset
similarity at t=+4.59 while COORD (GCCA coordinate-loading similarity) gives t=-0.46. I concluded
COORD was absent rather than hidden, because the two predictors correlate only -0.028.

Since then the coordinates COORD was built from have been damaged: #39 scoped SUBSTANCE to 3% of
the corpus, #41 found POWER uniformly weak, #18/#19 withdrew the coordinate COUNT as K-dependent.
A null from a measure that has since failed its own validation is silence, not absence -- P6, and
the reason UNVERIFIED exists.

The person factors are the one construct here that PASSED a generalisation check (#40: they predict
blocks the filter discarded at 46% of the included effect, 2.8x a placebo). So rebuild COORD from
them and re-run the comparison.

ESTIMAND        the partial contribution of coordinate similarity to within-person onset similarity,
                measured with a coordinate construct that has passed a generalisation test.
IDENTIFICATION  identified; both predictors are computable for the same category pairs.
WORLDS          A  two systems: COORD stays null even with a validated coordinate measure
                B  measurement failure: COORD becomes significant, and A03's decision flips
KILL (CONDITIONAL) gate: TEMPO must reproduce its published t>4 (the positive control from #9) AND
                   the permuted-COORD null must be ~0. Otherwise UNVERIFIED.
                   then: |t_COORD| > 2.5 -> A03's decision FLIPS, acquisition tracks coordinates
                         |t_COORD| < 1.5 -> two systems CONFIRMED with a validated measure
                         otherwise       -> UNVERIFIED
POSITIVE CTRL   TEMPO, whose coefficient is published at t=+4.59.
NEGATIVE CTRL   permuted COORD values across category pairs.
SHAM            a random category-level vector standing in for the coordinate profile.
SEEDS           3.
MULTIPLICITY    2 coordinate constructions (GCCA as published, person-factor) x 3 seeds, both
                reported side by side.
IMPOSSIBLE      a coordinate construct validated OUTSIDE this release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('24_attack_rsa.py')).read().split('print("=== is onset a proxy')[0])
Ostr=Ores.copy()
for c in Ores.columns:
    m=Ores[c].notna()&Pres[c].notna()
    X=np.c_[np.ones(m.sum()),Pres.loc[m,c].values]; b,*_=lstsq(X,Ores.loc[m,c].values,rcond=None)
    Ostr.loc[m,c]=Ores.loc[m,c].values-X@b
CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan)
for i in range(k):
    for j in range(i+1,k):
        m=Ostr.iloc[:,i].notna()&Ostr.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
        if m.sum()<150: continue
        CO[i,j]=CO[j,i]=np.corrcoef(Ostr.iloc[:,i][m],Ostr.iloc[:,j][m])[0,1]
        CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
# --- build the PERSON-FACTOR coordinate profile per category (the construct that passed #40) ---
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
pm={p:i for i,p in enumerate(pool)}; cols=[]
for q in allq:
    idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
    Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
    if len(idx): Z[idx]=B[q]['R'][src]
    mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
    cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
Zall=np.hstack(cols); Zall=Zall-Zall.mean(0)
D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
b,*_=lstsq(D,Zall,rcond=None); Zall=Zall-D@b
U,S,_=svd(Zall,full_matrices=False); FAC=pd.DataFrame(U[:,:6]*S[:6],index=pool)
catcols=[v for _,v in pairs]
prof=np.full((k,6),np.nan)
for ci,cat in enumerate(catcols):
    y=pd.to_numeric(df[cat],errors='coerce').reindex(pool)
    m=y.notna()
    if m.sum()<400: continue
    prof[ci]=[np.corrcoef(FAC.loc[m,j],y[m])[0,1] for j in range(6)]
have=~np.isnan(prof).any(1)
print(f"categories with a person-factor profile: {have.sum()} of {k}")
meanons=np.array([O[c].mean() for c in catcols])
iu=np.triu_indices(k,1)
ok=(~np.isnan(CO[iu]))&(~np.isnan(CP[iu]))&have[iu[0]]&have[iu[1]]
def cos(a,b): return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12))
COORD_PF=np.array([cos(prof[i],prof[j]) for i,j in zip(iu[0][ok],iu[1][ok])])
TEMPO=-np.abs(meanons[iu[0][ok]]-meanons[iu[1][ok]])
y=CO[iu][ok]; pref=CP[iu][ok]
def fit(extra,labels):
    X=np.c_[np.ones(len(y)),pref,*extra]
    b,*_=lstsq(X,y,rcond=None); r=y-X@b
    se=np.sqrt((r**2).sum()/(len(y)-X.shape[1])*np.diag(np.linalg.pinv(X.T@X)))
    return {l:(b[2+i],b[2+i]/se[2+i]) for i,l in enumerate(labels)}
print(f"\npairs usable: {ok.sum()}")
res=fit([COORD_PF,TEMPO],['COORD_personfactor','TEMPO'])
rng=np.random.default_rng(5)
nulls=[fit([rng.permutation(COORD_PF),TEMPO],['C','T'])['C'][1] for _ in range(200)]
sham=[fit([rng.normal(size=len(COORD_PF)),TEMPO],['C','T'])['C'][1] for _ in range(200)]
print("\n=== COORD rebuilt from the person factors (the construct that passed #40) ===")
for lab,(coef,t) in res.items(): print(f"   {lab:22s} b={coef:+.4f}  t={t:+.2f}")
print(f"\n   permuted-COORD null t: mean {np.mean(nulls):+.2f}  sd {np.std(nulls):.2f}  |t| p95 {np.percentile(np.abs(nulls),95):.2f}")
print(f"   sham (gaussian) t    : |t| p95 {np.percentile(np.abs(sham),95):.2f}")
tT=res['TEMPO'][1]; tC=res['COORD_personfactor'][1]
gate_pos=tT>4; gate_neg=abs(np.mean(nulls))<1.0
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  TEMPO reproduces t>4 (published +4.59) : {'PASS' if gate_pos else 'FAIL'} (t={tT:+.2f})")
print(f"  permuted-COORD null centred near 0     : {'PASS' if gate_neg else 'FAIL'} ({np.mean(nulls):+.2f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
elif abs(tC)>2.5: print(f"  -> A03 FLIPS : COORD t={tC:+.2f} with a validated coordinate measure")
elif abs(tC)<1.5: print(f"  -> TWO SYSTEMS CONFIRMED : COORD t={tC:+.2f} even with a validated measure")
else: print(f"  -> UNVERIFIED : COORD t={tC:+.2f}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
