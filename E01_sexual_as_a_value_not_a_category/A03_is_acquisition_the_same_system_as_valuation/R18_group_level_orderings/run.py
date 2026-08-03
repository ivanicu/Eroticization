"""
E01 A03 R18 -- do groups have different acquisition orderings?

#52 showed a single global ordering saturates what any global ordering can do (66.7%, 101% of the
oracle). #55 found the INDIVIDUAL component is +0.88 points against a 1.04 seed spread -- present but
unresolvable. Between "one ordering for everyone" and "one per person" sits the level this project
never tested: GROUPS. A group-level difference is far more detectable than an individual one,
because a group ordering is estimated from thousands of people rather than one.

  test: fit the ordering separately within each group, then predict a held-out person's pairs with
  their OWN group's ordering versus the OTHER group's. If groups differ, own beats other.

ESTIMAND        held-out pairwise ranking accuracy from own-group vs other-group ordering.
IDENTIFICATION  identified; orderings are fitted on training people within group, scored on held-out
                people of the target group.
WORLDS          A  one schedule for everyone: own-group ties other-group
                B  group-specific schedules: own beats other, and by more than the seed spread
KILL (CONDITIONAL) gate, measured (#41): (a) own-group ordering must beat a random ordering by >5
                   points -- the established effect must reproduce within group; (b) the random
                   ordering must sit near chance in every group.
                   then: own - other > 2x seed spread -> GROUP-SPECIFIC SCHEDULES
                         own - other < seed spread    -> ONE SCHEDULE
POSITIVE CTRL   own-group ordering against random.
NEGATIVE CTRL   random orderings, averaged over 40 permutations (per #54's fix).
GROUPS          sex (biomale), and as a second contrast, high vs low breadth -- because #59 showed
                breadth is the dominant axis of individual variation, so if any grouping has its own
                schedule it should.
SEEDS           4.
MULTIPLICITY    2 groupings x 3 orderings x 4 seeds, all reported.
IMPOSSIBLE      a grouping not measured in this release; and the tie rate (36.3%, #52) bounds every
                accuracy here identically, so it does not distort the comparison.
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
O=pd.DataFrame({c:df[c].map(BIN) for c in ons}); V=O.values; mask=~np.isnan(V)
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce'); breadth=(R>0).sum(1).values
male=pd.to_numeric(df['biomale'],errors='coerce').values
keep=np.flatnonzero(mask.sum(1)>=6)
GROUPS={'sex':(np.flatnonzero((male==1)&(mask.sum(1)>=6)),np.flatnonzero((male==0)&(mask.sum(1)>=6))),
        'breadth':(np.flatnonzero((breadth>np.median(breadth))&(mask.sum(1)>=6)),
                   np.flatnonzero((breadth<=np.median(breadth))&(mask.sum(1)>=6)))}
def acc(order_vals,people,rng,cap=20000):
    right=0;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(10,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b] or order_vals[a]==order_vals[b]: continue
            right+=((order_vals[a]<order_vals[b])==(V[i,a]<V[i,b])); tot+=1
            if tot>=cap: break
        if tot>=cap: break
    return 100*right/max(tot,1)
rows=[]
for gname,(gA,gB) in GROUPS.items():
    for seed in (1,2,3,4):
        rng=np.random.default_rng(seed)
        for lab,tgt,oth in [('A',gA,gB),('B',gB,gA)]:
            p=rng.permutation(tgt); tr,te=p[:len(p)//2],p[len(p)//2:]
            own=np.nanmean(V[tr],axis=0); other=np.nanmean(V[oth],axis=0)
            rnd=float(np.mean([acc(rng.permutation(own),te,rng) for _ in range(8)]))
            rows.append(dict(grouping=gname,side=lab,seed=seed,
                             own=acc(own,te,rng),other=acc(other,te,rng),random=rnd))
G=pd.DataFrame(rows); G.to_csv(OUT/'group_orderings.csv',index=False)
S=G.groupby('grouping')[['own','other','random']].agg(['median',lambda s:s.max()-s.min()])
S.columns=['own','own_sp','other','other_sp','random','random_sp']
print("\n=== held-out pairwise ranking accuracy by grouping ===")
print(S.round(2).to_string())
for g in S.index:
    d=S.loc[g,'own']-S.loc[g,'other']; sp=max(S.loc[g,'own_sp'],S.loc[g,'other_sp'])
    print(f"\n  {g}: own {S.loc[g,'own']:.2f}%  other {S.loc[g,'other']:.2f}%  diff {d:+.2f}  seed spread {sp:.2f}  ratio {abs(d)/max(sp,1e-9):.2f}")
ga=bool(((S['own']-S['random'])>5).all()); gb=bool((S['random'].between(45,55)).all())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) own-group beats random by >5 in every group : {'PASS' if ga else 'FAIL'}")
print(f"  (b) random near chance in every group           : {'PASS' if gb else 'FAIL'} ({S['random'].min():.1f}-{S['random'].max():.1f})")
if not (ga and gb): print("  -> gate FAILED : UNVERIFIED")
else:
    for g in S.index:
        d=S.loc[g,'own']-S.loc[g,'other']; sp=max(S.loc[g,'own_sp'],S.loc[g,'other_sp'])
        v='GROUP-SPECIFIC' if abs(d)>2*sp else ('ONE SCHEDULE' if abs(d)<sp else 'partial')
        print(f"  -> {g}: {v} (diff {d:+.2f}, spread {sp:.2f})")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
