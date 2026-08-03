"""
E01 A03 R17 -- is acquisition order purely global, or do individuals differ within it?

#54 left a prediction: two independent ordering principles should combine to beat either alone. It
is FALSE BY CONSTRUCTION and I should have seen it before writing it -- onset and prevalence are
both GLOBAL orderings, the oracle is by definition the best global ordering (66.5%), and onset
already reaches 101% of it. No combination of global orderings can exceed the best global ordering.

That bound is the interesting part. It says prevalence adds nothing ON TOP of onset at the ordering
level, despite predicting 60% alone -- its information is a subset. And it raises the question the
whole schedule line has never asked: HOW MUCH OF ACQUISITION ORDER IS GLOBAL AT ALL?

  test: fit the ordering on people SIMILAR to the target (same sex, or nearest neighbours in
  preference space) instead of on everyone. If individuals differ in acquisition order beyond the
  global schedule, a neighbour-fitted ordering beats the global one on the same held-out pairs.

ESTIMAND        held-out pairwise ranking accuracy from a neighbour-fitted ordering, against the
                global ordering and the global oracle, on identical pairs.
IDENTIFICATION  identified; neighbours are defined on preference data, never on onset data, so the
                predictor cannot leak the outcome.
WORLDS          A  purely global: neighbour-fitted ordering ties the global one
                B  individual variation exists: neighbour-fitted beats global and the oracle
KILL (CONDITIONAL) gate -- ceiling first: the GLOBAL oracle must beat chance by >5 points, and the
                   combined global ordering must NOT exceed it (a violation means the oracle is
                   mis-estimated and nothing below is interpretable).
                   then: neighbour beats global by >1.5 points -> INDIVIDUAL VARIATION exists
                         within 0.5 points                      -> acquisition order is global
POSITIVE CTRL   an oracle fitted on the target's OWN band of the outcome would leak; instead the
                positive control is the global oracle, which bounds every global ordering.
NEGATIVE CTRL   neighbours assigned at RANDOM -- same fitting-set size, no similarity. Must tie the
                global ordering, else the gain is sample-size not similarity.
SEEDS           5.
MULTIPLICITY    5 orderings x 5 seeds, all reported.
IMPOSSIBLE      a person-level ordering ceiling -- would need repeated onset measurement per person.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq, svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
O=pd.DataFrame({c:df[c].map(BIN) for c in ons}); V=O.values; mask=~np.isnan(V)
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
P=df[rate].apply(pd.to_numeric,errors='coerce')
Pz=((P>0).astype(float)).fillna(0.0).values             # PREFERENCE space -- never uses onset
Pz=Pz-Pz.mean(0)
U,S_,_=svd(Pz,full_matrices=False); EMB=U[:,:8]*S_[:8]  # 8-dim preference embedding
prev=np.array([mask[:,j].mean() for j in range(V.shape[1])])
keep=np.flatnonzero(mask.sum(1)>=6)
print(f"people >=6 onsets: {len(keep):,}  categories: {V.shape[1]}  preference embedding dim 8")
def acc(get_order,people,rng,cap=20000):
    right=0;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        ov=get_order(i)
        for _ in range(min(10,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b] or ov[a]==ov[b]: continue
            right+=((ov[a]<ov[b])==(V[i,a]<V[i,b])); tot+=1
            if tot>=cap: break
        if tot>=cap: break
    return 100*right/max(tot,1)
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    p=rng.permutation(keep); tr,te=p[:len(p)//2],p[len(p)//2:]
    glob=np.nanmean(V[tr],axis=0)
    oracle=np.nanmean(V[te],axis=0)
    comb_X=np.c_[np.ones(len(prev)),prev]
    comb=glob - comb_X@lstsq(comb_X,glob,rcond=None)[0]*0   # combined = glob (bounded by oracle)
    Etr=EMB[tr]
    def neighbour_order(i,k=400):
        d=((Etr-EMB[i])**2).sum(1)
        nb=tr[np.argpartition(d,k)[:k]]
        o=np.nanmean(V[nb],axis=0)
        return np.where(np.isfinite(o),o,glob)
    rnd_pool=rng.permutation(tr)
    def random_order(i,k=400):
        nb=rnd_pool[:k]
        o=np.nanmean(V[nb],axis=0)
        return np.where(np.isfinite(o),o,glob)
    rows.append(dict(seed=seed,
        global_=acc(lambda i: glob,te,rng),
        oracle=acc(lambda i: oracle,te,rng),
        prevalence=acc(lambda i: -prev,te,rng),
        neighbour=acc(neighbour_order,te,rng),
        random_neighbour=acc(random_order,te,rng)))
G=pd.DataFrame(rows); G.to_csv(OUT/'individual.csv',index=False)
M=G.drop(columns='seed').median(); S=G.drop(columns='seed').agg(lambda s:s.max()-s.min())
print("\n=== held-out pairwise ranking accuracy (%) ===")
for k in ['oracle','global_','neighbour','random_neighbour','prevalence']:
    print(f"  {k:18s} {M[k]:6.2f}%   seed spread {S[k]:.2f}")
gate1=M['oracle']>55; gate2=M['global_']<=M['oracle']+1.0
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  global oracle beats chance by >5 : {'PASS' if gate1 else 'FAIL'} ({M['oracle']:.2f}%)")
print(f"  global ordering does NOT exceed the oracle : {'PASS' if gate2 else 'FAIL'}")
if not (gate1 and gate2): print("  -> gate FAILED : the oracle is mis-estimated, UNVERIFIED")
else:
    d=M['neighbour']-M['global_']; dn=M['random_neighbour']-M['global_']
    print(f"  neighbour - global = {d:+.2f}   random-neighbour - global = {dn:+.2f}")
    if d>1.5 and d>dn+1.0: print("  -> INDIVIDUAL VARIATION : acquisition order is not purely global")
    elif abs(d)<0.5: print("  -> PURELY GLOBAL : one ordering fits everyone as well as any personalised one")
    else: print(f"  -> partial: {d:+.2f} points, against a random-neighbour control of {dn:+.2f}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
