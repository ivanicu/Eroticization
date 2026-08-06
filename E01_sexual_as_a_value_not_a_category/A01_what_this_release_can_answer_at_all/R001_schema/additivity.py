import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
ITER 11. Ivan's model B is v_i = w_i^T phi(s). That REQUIRES additivity: adding feature f to a
scene shifts erotic value by a consistent amount, whatever the scene. His A_i interaction term is
the admission that it might not.
Direct test: the fluid blocks instantiate ONE option template across 7 substances, so the same
feature contrast (self vs other) appears seven times with the substance held out.
  additive phi      -> the self-minus-other difference vector points the same way for saliva as
                       for urine as for ejaculate. Parallel across substances.
  interactive phi   -> 'receiving' means something different depending on what is received.
Null: the same computation for MISMATCHED acts, which shares the substance and the people but
not the feature.
"""
import pandas as pd, numpy as np, re, warnings, itertools
from scipy import stats
warnings.filterwarnings('ignore'); rng=np.random.default_rng(2178309)
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
o=pd.read_csv('data/derived/options.csv')
FLUID={7:'precum',8:'saliva',9:'squirt',11:'urine',83:'ejaculate',6:'breastmilk',10:'sweat'}
# the template: four act-pairs, each with a self pole and an other pole
ACTS={'consume':(r'consuming it myself', r'others consuming it'),
      'produce':(r'^(making|ejaculating|squirting) .*myself|myself$', r'^others (making|ejaculating|squirting)'),
      'play'   :(r'playing with it myself', r'others playing with it'),
      'orifice':(r'into my orifices', r"into others' orifices")}
D={}
for qi,name in FLUID.items():
    sub=lg[lg.qi==qi]
    if not len(sub): continue
    ppl=np.array(sorted(sub.person.unique())); opts=np.array(sorted(sub.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={x:i for i,x in enumerate(opts)}
    M=np.zeros((len(ppl),len(opts))); M[sub.person.map(pi).values, sub.option.map(oi).values]=1
    low=pd.Series(opts).str.lower()
    for act,(rs,ro) in ACTS.items():
        a=np.flatnonzero(low.str.contains(rs,regex=True).values)
        b=np.flatnonzero(low.str.contains(ro,regex=True).values)
        if len(a)==0 or len(b)==0: continue
        D[(name,act)]=pd.Series(M[:,a].mean(1)-M[:,b].mean(1), index=ppl)
print(f"difference vectors built: {len(D)}  ({len(set(k[0] for k in D))} substances x {len(set(k[1] for k in D))} acts)")
def r_between(k1,k2,minn=250):
    a,b=D[k1],D[k2]; c=a.index.intersection(b.index)
    if len(c)<minn: return np.nan,0
    x,y=a.reindex(c).values,b.reindex(c).values
    if x.std()==0 or y.std()==0: return np.nan,0
    return np.corrcoef(x,y)[0,1], len(c)
same_act=[]; diff_act=[]
for k1,k2 in itertools.combinations(D,2):
    s1,a1=k1; s2,a2=k2
    if s1==s2: continue                       # different substances only
    r,n=r_between(k1,k2)
    if np.isnan(r): continue
    (same_act if a1==a2 else diff_act).append((r,n,k1,k2))
sa=np.array([x[0] for x in same_act]); da=np.array([x[0] for x in diff_act])
print(f"\nSAME feature-contrast, different substance : n_pairs={len(sa):3d}  mean r={sa.mean():+.3f}  median={np.median(sa):+.3f}")
print(f"DIFFERENT contrast, different substance    : n_pairs={len(da):3d}  mean r={da.mean():+.3f}  median={np.median(da):+.3f}")
print(f"  difference = {sa.mean()-da.mean():+.3f}   Welch t = {stats.ttest_ind(sa,da,equal_var=False).statistic:+.2f}"
      f"  p={stats.ttest_ind(sa,da,equal_var=False).pvalue:.2e}")
print(f"\n  additivity would put SAME-contrast pairs near the reliability ceiling and")
print(f"  DIFFERENT-contrast pairs near zero. Observed same={sa.mean():.3f} diff={da.mean():+.3f}")
print("\n  strongest same-contrast pairs:")
for r,n,k1,k2 in sorted(same_act,key=lambda x:-x[0])[:6]:
    print(f"    r={r:+.3f} n={n:5d}  {k1[1]:8s} {k1[0]:11s} <-> {k2[0]}")
print("  weakest same-contrast pairs:")
for r,n,k1,k2 in sorted(same_act,key=lambda x:x[0])[:4]:
    print(f"    r={r:+.3f} n={n:5d}  {k1[1]:8s} {k1[0]:11s} <-> {k2[0]}")
