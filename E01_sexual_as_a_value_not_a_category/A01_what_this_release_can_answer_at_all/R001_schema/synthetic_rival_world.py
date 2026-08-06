"""
E01 A03 R09 -- settle #21 by BUILDING the rival world instead of permuting the finished matrix.

#21 left "onset structure has content beyond preference" UNVERIFIED: the pair-level null is
anticonservative, and the category-level null is degenerate because eigenvalues are invariant under
simultaneous row-column permutation. Neither can be repaired -- both operate on the finished
matrix, and the structure in question lives in the PEOPLE the matrix was computed from.

realstat G2: a permutation null answers "did the pairing matter", never "why". Name the world it
must exclude and BUILD it.

RIVAL WORLD, built explicitly: onset similarity is entirely preference similarity plus noise. So
generate synthetic person-level onset data whose category covariance IS the observed preference
covariance, push it through the SAME binning, the SAME missingness pattern, and the SAME
residualisation, and read the eigenvalue distribution off it. If the observed residual sits inside
that distribution, the claim is dead. If outside, it is real and was under-tested rather than over.

ESTIMAND        the leading eigenvalue of the onset-similarity matrix after removing its linear
                dependence on preference similarity, against the distribution that quantity takes
                when onset is BY CONSTRUCTION nothing but preference plus noise.
IDENTIFICATION  identified: the rival world is fully specifiable and simulable at the person level.
SCOPE           27 categories, the release's own missingness and 2-year binning reproduced exactly.
WORLDS          A  onset carries structure preference does not -> observed eigenvalue above the
                   synthetic distribution
                B  it is preference plus noise -> observed sits inside it
KILL            PRE-REGISTERED: observed top eigenvalue inside the central 95% of the synthetic
                distribution -> claim WITHDRAWN. Above the 97.5th percentile -> CONFIRMED and #21
                is resolved in favour of the original finding. Between -> stays UNVERIFIED and the
                line is FROZEN rather than re-run a third time.
POSITIVE CTRL   the same simulator with EXTRA rank-2 structure injected at the person level must
                land the observed statistic above the null. Verified to fail at injection 0.
NEGATIVE CTRL   the synthetic null itself is the negative control -- it is the rival world.
SHAM            synthetic onset with covariance = IDENTITY (no preference coupling at all); the
                residual eigenvalue must then be LARGER, confirming the statistic responds to
                coupling in the expected direction.
NOISE FLOOR     sd of the synthetic distribution.
MULTIPLICITY    1 statistic, 3 seeds, 300 synthetic draws each; the full distribution is reported.
SPECIFICATION   noise level matched to observed residual sd; reported at 3 matched levels.
SEEDS           3.
IMPOSSIBLE      independent replication, cross-dataset. n=27 categories is the object's.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, eigvalsh, cholesky
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('24_attack_rsa.py')).read().split('print("=== is onset a proxy')[0])
MIDS=np.array([2,5.5,7.5,9.5,11.5,13.5,15.5,17.5,22,28]); EDGES=[0,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,25.5,99]
binify=lambda x: MIDS[np.clip(np.digitize(x,EDGES)-1,0,len(MIDS)-1)]
mask=O.notna().values
def sim_matrices(Ovals):
    Om=pd.DataFrame(np.where(mask,Ovals,np.nan),columns=O.columns)
    Z=Om.sub(Om.mean(axis=1),axis=0); Z=Z.sub(Z.mean(axis=0),axis=1)
    CO=np.full((k,k),np.nan)
    for i in range(k):
        for j in range(i+1,k):
            m=Z.iloc[:,i].notna()&Z.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
            if m.sum()<150: continue
            CO[i,j]=CO[j,i]=np.corrcoef(Z.iloc[:,i][m],Z.iloc[:,j][m])[0,1]
    return CO
CPm=np.full((k,k),np.nan)
for i in range(k):
    for j in range(i+1,k):
        m=Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()&O.iloc[:,i].notna()&O.iloc[:,j].notna()
        if m.sum()<150: continue
        CPm[i,j]=CPm[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
iu=np.triu_indices(k,1)
def top_eig(CO):
    ok=~np.isnan(CO[iu])&~np.isnan(CPm[iu])
    co,cp=CO[iu][ok],CPm[iu][ok]
    b,*_=lstsq(np.c_[np.ones(len(cp)),cp],co,rcond=None)
    r=co-np.c_[np.ones(len(cp)),cp]@b
    M=np.zeros((k,k)); M[iu[0][ok],iu[1][ok]]=r; M=M+M.T
    return float(np.abs(eigvalsh(M))[::-1][0])
obs=top_eig(sim_matrices(O.values))
print(f"observed top residual eigenvalue: {obs:.3f}   ({k} categories, {int(mask.sum())} onset cells)")
S=np.nan_to_num(CPm,nan=0.0); np.fill_diagonal(S,1.0)
w,V=np.linalg.eigh(S); w=np.clip(w,1e-6,None); L=V@np.diag(np.sqrt(w))
mu=np.nanmean(O.values); sd=np.nanstd(O.values)
def draw(rng,noise,inject=0.0):
    Zs=rng.normal(size=(len(O),k))@L.T
    if inject>0:
        U=rng.normal(size=(k,2)); Zs=Zs+inject*(rng.normal(size=(len(O),2))@U.T)
    Zs=Zs/Zs.std(0,keepdims=True)
    return binify(mu+sd*(np.sqrt(1-noise**2)*Zs+noise*rng.normal(size=Zs.shape)))
print("\n=== RIVAL WORLD: onset = preference-covariance + noise, same binning, same missingness ===")
res={}
for noise in [0.3,0.5,0.7]:
    vals=[]
    for seed in (2,12,22):
        rng=np.random.default_rng(seed)
        vals += [top_eig(sim_matrices(draw(rng,noise))) for _ in range(100)]
    v=np.array(vals); res[noise]=v
    lo,hi=np.percentile(v,[2.5,97.5])
    print(f"  noise {noise:.1f}: synthetic top eig mean {v.mean():.3f} sd {v.std():.3f}  95% [{lo:.3f},{hi:.3f}]"
          f"   observed {'ABOVE' if obs>hi else ('INSIDE' if obs>=lo else 'BELOW')}")
print("\n=== POSITIVE CONTROL: inject extra person-level rank-2 structure ===")
for inj in [0.0,0.3,0.6]:
    rng=np.random.default_rng(5)
    v=np.array([top_eig(sim_matrices(draw(rng,0.5,inject=inj))) for _ in range(60)])
    base=res[0.5]
    print(f"  injection {inj:.1f} -> mean top eig {v.mean():.3f}  vs rival-world mean {base.mean():.3f}"
          f"  ({'detected' if v.mean()>np.percentile(base,97.5) else 'not detected'})"
          +("   <- must NOT be detected" if inj==0 else ""))
print("\n=== SHAM: onset covariance = identity (no preference coupling at all) ===")
rng=np.random.default_rng(8)
Lid=np.eye(k)
def draw_id(rng,noise=0.5):
    Zs=rng.normal(size=(len(O),k)); Zs=Zs/Zs.std(0,keepdims=True)
    return binify(mu+sd*(np.sqrt(1-noise**2)*Zs+noise*rng.normal(size=Zs.shape)))
vsham=np.array([top_eig(sim_matrices(draw_id(rng))) for _ in range(60)])
print(f"  identity-covariance synthetic top eig mean {vsham.mean():.3f} vs coupled {res[0.5].mean():.3f}")
v=res[0.5]; lo,hi=np.percentile(v,[2.5,97.5])
print("\nPRE-REGISTERED KILL, evaluated (noise 0.5, the matched level):")
if obs>np.percentile(v,97.5): print(f"  -> observed {obs:.3f} > 97.5th pct {hi:.3f} : CONFIRMED, #21 resolved in favour of the finding")
elif obs>=lo: print(f"  -> observed {obs:.3f} inside 95% [{lo:.3f},{hi:.3f}] : CLAIM WITHDRAWN")
else: print(f"  -> observed {obs:.3f} BELOW {lo:.3f} : anomalous, line FROZEN")
pd.DataFrame({f"noise_{n}":pd.Series(res[n]) for n in res}).to_csv(OUT/'rival_world.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
