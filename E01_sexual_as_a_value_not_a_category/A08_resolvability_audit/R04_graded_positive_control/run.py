"""
E01 A08 R04 -- a graded positive control for the pairwise-block CCA, replacing the degenerate one.

#37 reproduced the headline (raw 0.2686 +/- 0.0137, adjusted 0.1980 +/- 0.0169, null 0.0509) and
found both RESOLVABLE -- but its positive control was a block split against itself, which returned
exactly 1.0000 because row-centring forces sum(half1) = -sum(half2). A control fixed by geometry
tests nothing, and unity is the signature of a constraint rather than a strong effect.

A binary control ("must be high") cannot distinguish detection from constraint anyway. This one is
GRADED: pairs ordered by independently known similarity, and the pipeline must reproduce the
ordering, not merely clear a bar.

  tier 1  precum <-> ejaculate      same option template, same source class, role scores r=0.729
                                    in A02 R21 -- the most similar pair in the release
  tier 2  other within-fluid-family same template, different substance
  tier 3  fluid <-> non-fluid       different template
  tier 4  all other pairs           the population the headline is a median over
  floor   person-permuted           the null

ESTIMAND        the held-out max canonical correlation per tier, and whether the tier ordering holds.
WORLDS          A  the pipeline detects alignment: tier1 > tier2 > tier4 > floor, all sub-unity
                B  it reads constraint or noise: ordering breaks, or a tier hits unity
KILL (CONDITIONAL) gate: no tier may reach 0.99 (that would signal a constraint, as in #37) AND the
                   floor must sit below tier 4.
                   then: tiers strictly ordered 1>2>4>floor -> the control is VALID and #37's
                         resolvability verdict is re-gated on it
                         ordering broken -> the pipeline cannot grade similarity; #37 stays on two
                         of three controls and says so
POSITIVE CTRL   this round IS the control being built; its own validity is judged by the ordering.
NEGATIVE CTRL   person-permutation within each tier.
SEEDS           4.
MULTIPLICITY    4 tiers x 4 seeds, every pair's value retained.
IMPOSSIBLE      a pair with known ground-truth similarity from outside the release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
FLUID={7,8,9,11,83,6,10}
def cvcca(Xa,Xb,nc,rng):
    idx=rng.permutation(len(Xa)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    try:
        c=CCA(n_components=nc,max_iter=800).fit(Xa[tr],Xb[tr]); ua,ub=c.transform(Xa[te],Xb[te])
        return max(abs(np.corrcoef(ua[:,j],ub[:,j])[0,1]) for j in range(nc))
    except Exception: return np.nan
def tier_of(a,b):
    if {a,b}=={7,83}: return 'tier1 precum<->ejaculate'
    if a in FLUID and b in FLUID: return 'tier2 within fluid family'
    if (a in FLUID) != (b in FLUID): return 'tier3 fluid<->non-fluid'
    return 'tier4 all other pairs'
rows=[]
for seed in (1,2,3,4):
    rng=np.random.default_rng(seed)
    for a,b in itertools.combinations(allq,2):
        A_,B_=B[a],B[b]; common=np.intersect1d(A_['ppl'],B_['ppl'])
        if len(common)<600: continue
        ia=np.searchsorted(A_['ppl'],common); ib=np.searchsorted(B_['ppl'],common)
        Ra,Rb=A_['R'][ia],B_['R'][ib]
        nc=min(3,Ra.shape[1]-1,Rb.shape[1]-1)
        rows.append(dict(seed=seed,a=a,b=b,tier=tier_of(a,b),
                         cca=cvcca(Ra,Rb,nc,rng),
                         null=cvcca(Ra,Rb[rng.permutation(len(common))],nc,rng)))
G=pd.DataFrame(rows); G.to_csv(OUT/'graded_control.csv',index=False)
T=G.groupby('tier').agg(n_pairs=('cca','size'),median=('cca','median'),
                        p90=('cca',lambda s: s.quantile(.9)),max=('cca','max'),
                        null=('null','median')).round(4)
print("=== graded positive control: does the pipeline reproduce known similarity ordering? ===")
print(T.to_string())
order=['tier1 precum<->ejaculate','tier2 within fluid family','tier4 all other pairs']
have=[t for t in order if t in T.index]
vals=[float(T.loc[t,'median']) for t in have]
floor=float(G['null'].median())
print(f"\n  floor (person-permuted, all pairs): {floor:.4f}")
print(f"  tier medians in order: {[round(v,4) for v in vals]}")
gate_unity=bool(T['max'].max()<0.99)
gate_floor=bool(floor<float(T.loc['tier4 all other pairs','median']))
ordered=all(vals[i]>vals[i+1] for i in range(len(vals)-1))
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  no tier reaches 0.99 (would signal a constraint) : {'PASS' if gate_unity else 'FAIL'} (max {T['max'].max():.4f})")
print(f"  floor below tier4                                : {'PASS' if gate_floor else 'FAIL'}")
if not (gate_unity and gate_floor): print("  -> gate FAILED : UNVERIFIED")
elif ordered:
    print(f"  -> tiers strictly ordered {' > '.join(str(round(v,3)) for v in vals)} > floor {floor:.3f}")
    print("  -> CONTROL VALID. #37's resolvability verdict is now gated on a working positive control.")
else:
    print(f"  -> ordering BROKEN {vals} : the pipeline cannot grade similarity; #37 stays on two of three controls")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
