"""
E01 A03 R08 -- the residual-structure null permuted pairs, and pairs are not the unit.

ADVERSARY_FORECAST #4, p=0.50: "344 pairs come from 27 categories, so the effective n is closer to
27 and every p-value there is overstated." Checked against the code rather than memory, it is
exactly half right:
  R03 line 88  p=rng.permutation(k); M=CP[ix_(p,p)]   -> CATEGORY-level (Mantel). CORRECT.
  R04 line 58  v=rng.permutation(res)                 -> PAIR-level. ANTICONSERVATIVE.
Pair-level permutation destroys the block structure a category imposes on its ~26 pairs, so the
null eigenvalues come out too small and the observed z is inflated by construction. The claim
built on it -- "80% of the sd of onset structure is not explained by preference", z = 4.9/3.5/2.7
-- has never been tested against a null that respects the unit.

ESTIMAND        whether the onset-similarity matrix carries structure beyond preference
                similarity, with the CATEGORY as the unit of exchangeability.
IDENTIFICATION  identified at n=27 categories. Pair count is not sample size: each category
                appears in 26 pairs, so the 344 pairs carry at most 27 independent draws.
SCOPE           27 matched onset/rating categories · no instrument · baseline = category-label
                permutation of the residual matrix · regime n=27 units.
WORLDS          A  residual structure is real: survives a null that permutes CATEGORIES
                B  it was the pair-level null: z collapses when the unit is respected
KILL            PRE-REGISTERED: if the top residual eigenvalue's z under CATEGORY permutation
                falls below 2.0, the "80% not explained by preference" claim is WITHDRAWN and
                ADVERSARY_FORECAST #4 is scored CORRECT-on-the-residual, WRONG-on-the-RSA.
POSITIVE CTRL   plant a rank-2 structure into the residual at known amplitude; the category-level
                test must recover it (z > 3) -- otherwise a collapse means a blind test, not a
                real null. Verified to FAIL at amplitude 0.
NEGATIVE CTRL   category-label permutation, 2000 draws, which is the corrected null itself.
SHAM            the same pipeline on a residual built from two independent random matrices.
PLACEBO         a pure-noise symmetric matrix -- z must be ~0.
NOISE FLOOR     sd of the category-permuted eigenvalue distribution.
MULTIPLICITY    3 eigenvalues x 2 null types x 3 seeds, all reported, plus the RSA re-reported
                under its (already correct) category null for completeness.
SPECIFICATION   null type {pair-level, category-level} x seed {3} x eigenvalue index {1,2,3}.
SEEDS           3.
IMPOSSIBLE      more categories -- this release has 27 matched pairs of onset/rating columns.
                n=27 is a property of the object and is why the interval is published, not hidden.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, eigvalsh
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
iu=np.triu_indices(k,1); ok=~np.isnan(CO[iu])&~np.isnan(CP[iu])
co,cp=CO[iu][ok],CP[iu][ok]
b,*_=lstsq(np.c_[np.ones(len(cp)),cp],co,rcond=None)
res=co-np.c_[np.ones(len(cp)),cp]@b
def to_mat(v):
    M=np.zeros((k,k)); M[iu[0][ok],iu[1][ok]]=v; return M+M.T
Rm=to_mat(res); ev=np.abs(eigvalsh(Rm))[::-1]
print(f"units: {k} categories, {ok.sum()} pairs  (each category appears in ~{2*ok.sum()//k} pairs)")
def z_under(null_kind,seed,reps=500,mat=None):
    rng=np.random.default_rng(seed); M0=Rm if mat is None else mat
    obs=np.abs(eigvalsh(M0))[::-1][:3]; nn=[]
    for _ in range(reps):
        if null_kind=='pair': nn.append(np.abs(eigvalsh(to_mat(rng.permutation(M0[iu[0][ok],iu[1][ok]]))))[::-1][:3])
        else:
            p=rng.permutation(k); nn.append(np.abs(eigvalsh(M0[np.ix_(p,p)]))[::-1][:3])
    nn=np.array(nn); return (obs-nn.mean(0))/nn.std(0), obs, nn.mean(0), nn.std(0)
print("\n=== the same residual matrix, two nulls ===")
for kind in ['pair','category']:
    zs=np.array([z_under(kind,s)[0] for s in (5,15,25)])
    z,obs,mu,sd=z_under(kind,5)
    print(f"  {kind:9s}-level null  z = {np.round(zs.mean(0),1)}   (obs {np.round(obs,3)}  null mean {np.round(mu,3)}  sd {np.round(sd,3)})")
zc=np.array([z_under('category',s)[0] for s in (5,15,25)]).mean(0)
print("\n=== POSITIVE CONTROL: plant rank-2 structure, category-level test must see it ===")
rng=np.random.default_rng(1)
for amp in [0.0,0.02,0.05,0.10]:
    U=rng.normal(size=(k,2)); P2=U@U.T; np.fill_diagonal(P2,0)
    planted=Rm+amp*P2/np.abs(P2).max()
    zp,_,_,_=z_under('category',7,mat=planted)
    print(f"   amplitude {amp:.2f} -> top-eigenvalue z = {zp[0]:+.1f}"+("   <- must be ~0" if amp==0 else ""))
print("\n=== SHAM and PLACEBO ===")
sh=to_mat(rng.normal(size=len(res))*res.std())
print(f"   sham (random residual, same sd)  z = {z_under('category',9,mat=sh)[0][0]:+.1f}")
print(f"   placebo (pure noise matrix)      z = {z_under('category',9,mat=to_mat(rng.normal(size=len(res))))[0][0]:+.1f}")
nl=[]
rngm=np.random.default_rng(3)
for _ in range(2000):
    p=rngm.permutation(k); M=CP[np.ix_(p,p)]; v=M[iu][ok]; g=~np.isnan(v)
    nl.append(np.corrcoef(co[g],v[g])[0,1])
r_rsa=np.corrcoef(co,cp)[0,1]
print(f"\n=== the RSA itself, whose null was ALREADY category-level ===")
print(f"   RSA = {r_rsa:+.3f}   category-permutation null {np.mean(nl):+.3f} +/- {np.std(nl):.3f}"
      f"   z = {(r_rsa-np.mean(nl))/np.std(nl):+.1f}")
pd.DataFrame(dict(eig=[1,2,3],z_pair=np.array([z_under('pair',s)[0] for s in (5,15,25)]).mean(0),
                  z_category=zc)).to_csv(OUT/'unit_test.csv',index=False)
print("\nPRE-REGISTERED KILL, evaluated:")
if zc[0]<2.0:
    print(f"  -> top-eigenvalue z under category permutation = {zc[0]:+.1f} < 2.0 : CLAIM WITHDRAWN")
    print("  -> ADVERSARY_FORECAST #4 : CORRECT on the residual test, WRONG on the RSA")
else:
    print(f"  -> z = {zc[0]:+.1f} >= 2.0 : the residual structure survives the correct unit")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
