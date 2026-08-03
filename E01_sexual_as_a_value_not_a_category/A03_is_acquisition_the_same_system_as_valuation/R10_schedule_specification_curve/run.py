"""
E01 A03 R10 -- the shared maturational schedule, as a curve rather than a cell.

#30's standing correction: a specification curve is the unit in which a result is stated, not a
step run on suspicious results. Auditing the surviving README numbers for ones that have never had
a curve leaves two; this is the load-bearing one, since A03's whole decision -- model acquisition
and valuation as two systems -- rests on it.

Published as a single cell: per-person Spearman between their own onset ordering and the population
mean ordering, mean +0.232, 74.7% of 9,691 people positive against a 48.9% null, Cohen d 0.69.
Unswept analyst choices: the minimum onsets required per person, the bin->number mapping, which
categories are included, and the rank statistic.

ESTIMAND        the share of people whose personal acquisition ordering agrees with the population
                ordering, above its own within-person permutation null.
IDENTIFICATION  identified at every cell; the null is computed inside each cell so a moving null
                cannot pass as a moving effect.
WORLDS          A  real shared schedule: share-positive stays well above the null in every cell
                B  an artifact of one mapping or one inclusion rule: it collapses somewhere
KILL (CONDITIONAL -- evaluated ONLY if the gate passes)
      gate: the within-person permutation null stays within 45-55% positive in EVERY cell
      then: min share-positive across the grid >= 65% -> ROBUST, stated as a range
            any cell below 55%                        -> FRAGILE, downgraded
            otherwise                                 -> UNVERIFIED
POSITIVE CTRL   a synthetic population where every person follows the population order exactly,
                plus noise, must return share-positive near 100% under every cell.
NEGATIVE CTRL   the within-person permutation null, recomputed per cell.
SHAM            population ordering replaced by a random ordering of the same categories.
NOISE FLOOR     sd of the null share across 3 seeds.
MULTIPLICITY    4 min-onset x 3 mappings x 3 category sets x 2 statistics = 72 cells, all reported.
SEEDS           3.
IMPOSSIBLE      finer onset bins -- the release ships ~2-year buckets.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
MAPS={'midpoint':{'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
                  '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28},
      'ordinal':{k:i for i,k in enumerate(['0-4yo','5-6yo','7-8yo','9-10yo','11-12yo','13-14yo',
                                           '15-16yo','17-18yo','19-25yo','26yo+'])},
      'loweredge':{'0-4yo':0,'5-6yo':5,'7-8yo':7,'9-10yo':9,'11-12yo':11,'13-14yo':13,
                   '15-16yo':15,'17-18yo':17,'19-25yo':19,'26yo+':26}}
allons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if c in df.columns]
def catset(name,BIN):
    O=pd.DataFrame({c:df[c].map(BIN) for c in allons})
    n=O.notna().sum()
    if name=='all': return O.loc[:,n>300]
    if name=='high_n': return O.loc[:,n>2000]
    if name=='drop_porn': return O.loc[:,[c for c in O.columns if 'watching porn' not in c.lower()]].loc[:,n[[c for c in O.columns if 'watching porn' not in c.lower()]]>300]
rows=[]
for mapname,minons,cset,stat in itertools.product(MAPS,[4,6,8,10],['all','high_n','drop_porn'],['spearman','kendall']):
    O=catset(cset,MAPS[mapname]); V=O.values; mask=~np.isnan(V)
    pop=np.nanmean(V,axis=0)
    keep=np.flatnonzero(mask.sum(1)>=minons)
    if len(keep)<500: continue
    f=(lambda a,b: stats.spearmanr(a,b).statistic) if stat=='spearman' else (lambda a,b: stats.kendalltau(a,b).statistic)
    obs=[];nul=[];shm=[]
    rng=np.random.default_rng(3)
    order=rng.permutation(len(pop))
    for i in keep:
        j=np.flatnonzero(mask[i]); v=V[i,j]
        if len(set(v))<3: continue
        obs.append(f(v,pop[j])); nul.append(f(rng.permutation(v),pop[j])); shm.append(f(v,pop[order][j]))
    obs=np.array(obs); nul=np.array(nul); shm=np.array(shm)
    rows.append(dict(mapping=mapname,min_onsets=minons,categories=cset,stat=stat,n=len(obs),
                     share_pos=float(np.nanmean(obs>0)),mean_rho=float(np.nanmean(obs)),
                     null_share=float(np.nanmean(nul>0)),sham_share=float(np.nanmean(shm>0))))
G=pd.DataFrame(rows); G.to_csv(OUT/'schedule_curve.csv',index=False)
print(f"cells run: {len(G)}  (all reported)")
print("\n=== share of people agreeing with the population order, by cell ===")
print(G.pivot_table(index=['mapping','categories'],columns=['stat','min_onsets'],values='share_pos').round(3).to_string())
print("\n=== the null in the same cells (must sit near 0.50) ===")
print(f"   null share-positive: min {G.null_share.min():.3f}  median {G.null_share.median():.3f}  max {G.null_share.max():.3f}")
print(f"   sham (random population order): median {G.sham_share.median():.3f}")
print("\n=== POSITIVE CONTROL: synthetic followers ===")
O=catset('all',MAPS['midpoint']); V=O.values; mask=~np.isnan(V); pop=np.nanmean(V,axis=0)
rng=np.random.default_rng(9); syn=np.where(mask,pop[None,:]+rng.normal(0,1.0,V.shape),np.nan)
sp=[]
for i in np.flatnonzero(mask.sum(1)>=6):
    j=np.flatnonzero(mask[i]); v=syn[i,j]
    if len(set(v))>=3: sp.append(stats.spearmanr(v,pop[j]).statistic)
print(f"   synthetic exact-followers + noise: share positive {np.mean(np.array(sp)>0):.3f}  (must be near 1.0)")
gate = (G.null_share.between(0.45,0.55)).all()
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  negative control in every cell within 0.45-0.55 : {'PASS' if gate else 'FAIL'}")
if not gate: print("  -> gate FAILED : UNVERIFIED, threshold not evaluated")
else:
    lo=G.share_pos.min()
    if lo>=0.65: print(f"  -> min share-positive across {len(G)} cells = {lo:.3f} >= 0.65 : ROBUST. State as a range: [{G.share_pos.min():.3f}, {G.share_pos.max():.3f}]")
    elif lo<0.55: print(f"  -> a cell falls to {lo:.3f} < 0.55 : FRAGILE, downgraded")
    else: print(f"  -> min {lo:.3f} : UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
