"""
E01 A03 R15 -- is the maturational schedule a developmental fact, or just rarity?

#52 made the schedule the strongest surviving claim: a population ordering predicts held-out
people's pairwise acquisition order at 66.71%, which is 101% of what the best possible single global
ordering achieves. It has never been asked its most obvious rival.

RIVAL: rare interests are REPORTED later. Someone notices an uncommon interest later, or is less
certain when it began, or acquires it after the common ones simply because it is uncommon. If so the
"maturational schedule" is a PREVALENCE ordering wearing a developmental name, and prevalence is
directly observable in this release.

  discriminating test: run the ranking task with PREVALENCE ordering as a rival predictor, on the
  identical held-out pairs. If prevalence does as well, the schedule is rarity. If the onset
  ordering wins clearly, it is not.

ESTIMAND        held-out pairwise ranking accuracy from three orderings -- onset, prevalence, and
                onset-with-prevalence-partialled -- on the same pairs.
IDENTIFICATION  identified; prevalence is observed per category, no proxy needed.
WORLDS          A  developmental: onset ordering >> prevalence ordering, and the residual onset
                   ordering (prevalence removed) still beats chance clearly
                B  rarity: prevalence matches the onset ordering, and the residual collapses
KILL (CONDITIONAL) gate -- ceiling first (#50): the oracle must beat chance by >5 points.
                   then: residual-onset accuracy within 3 points of full onset -> DEVELOPMENTAL
                         prevalence within 3 points of onset AND residual collapsing to chance
                            -> RARITY, the strongest claim in the project is an artifact
                         otherwise -> partial, report the decomposition
POSITIVE CTRL   the in-sample oracle ordering.
NEGATIVE CTRL   a random ordering.
SHAM            a category-level nuisance with no developmental meaning -- the number of options in
                the category's block -- which must not predict acquisition order.
SEEDS           5.
MULTIPLICITY    5 orderings x 5 seeds, tie share reported.
IMPOSSIBLE      prospective onset measurement, which is the only thing that fully separates
                "acquired later" from "noticed later". Stated, not planned.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
br=pd.read_csv('data/derived/branching.csv'); qm=pd.read_csv('data/derived/multiselect_questions.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
O=pd.DataFrame({c:df[c].map(BIN) for c in ons}); V=O.values; mask=~np.isnan(V)
# PREVALENCE per category: how many people answered its onset item at all (= passed its gate)
prev=np.array([mask[:,j].mean() for j in range(V.shape[1])])
# SHAM: option count of the block this category gates into, if any
gate={str(r.gate).strip('"'):int(r.qi) for _,r in br.iterrows()}
nopt=[]
for c in ons:
    hit=[q for g,q in gate.items() if g and g.split('(')[0].strip().lower()[:6] in c.lower()]
    nopt.append(float(qm[qm.qi.isin(hit)].n_options.mean()) if hit else np.nan)
nopt=np.array(nopt); nopt=np.where(np.isfinite(nopt),nopt,np.nanmean(nopt))
keep=np.flatnonzero(mask.sum(1)>=6)
print(f"people >=6 onsets: {len(keep):,}  categories: {V.shape[1]}")
print(f"corr(prevalence, population mean onset) = {stats.pearsonr(prev,np.nanmean(V,axis=0))[0]:+.3f}")
def acc(order_vals,people,rng,cap=40000):
    right=0;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(12,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b] or order_vals[a]==order_vals[b]: continue
            right+=((order_vals[a]<order_vals[b])==(V[i,a]<V[i,b])); tot+=1
            if tot>=cap: break
        if tot>=cap: break
    return 100*right/max(tot,1)
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    idx=rng.permutation(keep); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    onset_tr=np.nanmean(V[tr],axis=0)
    X=np.c_[np.ones(len(prev)),prev]
    resid=onset_tr-X@lstsq(X,onset_tr,rcond=None)[0]      # onset ordering, prevalence removed
    rows.append(dict(seed=seed,
        onset=acc(onset_tr,te,rng),
        prevalence=acc(-prev,te,rng),                      # rarer -> later, so order by -prevalence
        onset_minus_prev=acc(resid,te,rng),
        oracle=acc(np.nanmean(V[te],axis=0),te,rng),
        sham_noptions=acc(nopt,te,rng),
        random=acc(rng.permutation(onset_tr),te,rng)))
G=pd.DataFrame(rows); G.to_csv(OUT/'schedule_or_rarity.csv',index=False)
M=G.drop(columns='seed').median(); S=G.drop(columns='seed').agg(lambda s:s.max()-s.min())
print("\n=== held-out pairwise ranking accuracy by ordering (%) ===")
for k in ['oracle','onset','onset_minus_prev','prevalence','sham_noptions','random']:
    print(f"  {k:18s} {M[k]:6.2f}%   seed spread {S[k]:.2f}")
gate_ok=M['oracle']>55
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  oracle beats chance by >5 points : {'PASS' if gate_ok else 'FAIL'} ({M['oracle']:.2f}%)")
if not gate_ok: print("  -> gate FAILED : UNVERIFIED")
else:
    dev = M['onset_minus_prev'] >= M['onset']-3
    rar = (M['prevalence']>=M['onset']-3) and (M['onset_minus_prev']<=52)
    if rar: print("  -> RARITY : the strongest claim in the project is a prevalence ordering")
    elif dev: print(f"  -> DEVELOPMENTAL : removing prevalence costs only {M['onset']-M['onset_minus_prev']:.2f} points")
    else: print(f"  -> PARTIAL : onset {M['onset']:.2f}%, prevalence {M['prevalence']:.2f}%, residual {M['onset_minus_prev']:.2f}%")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
