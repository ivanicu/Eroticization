"""
E01 A03 R14 -- the maturational schedule as a ranking task, completing the #51 sort.

#51 sorted every surviving headline by predictable variance except one. The schedule's statistic is
a within-person RANK agreement (share of people whose own onset ordering agrees with the
population's, [0.747, 0.860] across 72 cells), and forcing that through an R2 pipeline would be a
category error. Its proper predictive analogue is a RANKING task, and this round builds it.

  task     for a held-out person and a pair of their categories, predict which one they acquired
           first, using only the POPULATION ordering.
  chance   50% -- the pair is unordered a priori
  ceiling  the best a SINGLE GLOBAL ORDERING can do, estimated by fitting the ordering on the
           held-out people themselves (in-sample oracle). No global ordering can beat it, so the
           fraction (observed-50)/(oracle-50) is the ranking analogue of "% of ceiling".

Ties matter and are handled explicitly: onsets are 2-year bins, so many within-person pairs are
tied. Tied pairs are excluded from accuracy and their share is reported, because scoring them as
half-credit would hide how much of the task the binning has already decided.

ESTIMAND        held-out pairwise ranking accuracy from the population ordering, and its position
                between chance and the single-ordering oracle.
WORLDS          A  the schedule carries orderable information: accuracy well above 50 and a
                   substantial fraction of the oracle
                B  it is agreement without predictive content: accuracy at chance on held-out people
KILL (CONDITIONAL) gate -- CEILING FIRST (#50): the oracle must exceed chance by >5 points, else no
                   global ordering can do this task and nothing is sortable.
                   then: (obs-50)/(oracle-50) > 25% -> CARRIES, completing the sort
                         obs <= chance + 1           -> agreement without predictive content
POSITIVE CTRL   the in-sample oracle ordering.
NEGATIVE CTRL   a random ordering of the same categories.
NOISE FLOOR     across-seed spread, 5 seeds.
MULTIPLICITY    5 seeds x 3 orderings, tie share reported.
IMPOSSIBLE      test-retest onsets, which would give a true individual ceiling. The oracle used here
                bounds a GLOBAL ordering only, and is stated as such.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
O=pd.DataFrame({c:df[c].map(BIN) for c in ons})
V=O.values; mask=~np.isnan(V)
keep=np.flatnonzero(mask.sum(1)>=6)
print(f"people with >=6 onsets: {len(keep):,}   categories: {V.shape[1]}")
def accuracy(order_vals, people, rng, max_pairs=60000):
    right=0; total=0; tied=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(12,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b]: tied+=1; continue
            if order_vals[a]==order_vals[b]: continue
            pred_a_first = order_vals[a]<order_vals[b]
            true_a_first = V[i,a]<V[i,b]
            right += (pred_a_first==true_a_first); total+=1
            if total>=max_pairs: break
        if total>=max_pairs: break
    return (100*right/max(total,1), total, tied/max(total+tied,1))
rows=[]
for seed in (1,2,3,4,5):
    rng=np.random.default_rng(seed)
    idx=rng.permutation(keep); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    pop_tr=np.nanmean(V[tr],axis=0)                   # ordering fitted on TRAINING people
    oracle=np.nanmean(V[te],axis=0)                   # in-sample on the held-out half: the bound
    rnd=rng.permutation(pop_tr)
    a1,n1,t1=accuracy(pop_tr,te,rng); a2,_,_=accuracy(oracle,te,rng); a3,_,_=accuracy(rnd,te,rng)
    rows.append(dict(seed=seed,population=a1,oracle=a2,random=a3,pairs=n1,tie_share=t1))
G=pd.DataFrame(rows); G.to_csv(OUT/'ranking.csv',index=False)
M=G[['population','oracle','random']].median(); S=G[['population','oracle','random']].agg(lambda s:s.max()-s.min())
print("\n=== held-out pairwise ranking accuracy (%) ===")
for k in M.index: print(f"  {k:12s} {M[k]:6.2f}%   seed spread {S[k]:.2f}")
print(f"  tied pairs excluded: {100*G.tie_share.median():.1f}% of all within-person pairs")
print(f"  pairs scored per seed: {int(G.pairs.median()):,}")
gate=M['oracle']>55
frac=100*(M['population']-50)/max(M['oracle']-50,1e-9)
print(f"\nCONDITIONAL KILL -- gate first (ceiling before threshold, per #50)")
print(f"  oracle ordering beats chance by >5 points : {'PASS' if gate else 'FAIL'} ({M['oracle']:.2f}%)")
if not gate: print("  -> gate FAILED : no global ordering can do this task, UNSORTABLE")
else:
    print(f"  population ordering reaches {frac:.0f}% of the oracle's margin over chance")
    if frac>25: print("  -> CARRIES ORDERABLE INFORMATION : the sort in #51 is complete")
    elif M['population']<=51: print("  -> AGREEMENT WITHOUT PREDICTIVE CONTENT")
    else: print(f"  -> partial, {frac:.0f}% of the oracle margin")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
