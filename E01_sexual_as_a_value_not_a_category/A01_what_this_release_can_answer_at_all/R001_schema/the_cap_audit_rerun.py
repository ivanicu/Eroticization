import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A03 R21 -- RE-PRICING EVERY ROUND THAT SHARED THE CAP.

#73 found that #55's null was produced by `cap=20000` with `if tot>=cap: break`, where the break
exits the PERSON loop. A project-wide audit found the same helper copied into FIVE rounds:

  R14 schedule as a ranking task   -> the 66.5% headline itself
  R15 schedule or rarity
  R16 rarity or censoring
  R17 is the order individual      -> #55, already overturned by #73
  R18 group-level orderings        -> #62

Every one of them therefore ran on roughly 2,000 of 12,459 eligible people. The obvious reading is
"random truncation, so location is unbiased and only precision suffers" -- but that is an ASSUMPTION
about my own code, and assumptions about my own code are what produced this entry. It is tested here
rather than asserted.

ESTIMAND        for each capped comparison: the point estimate and the seed spread, computed at the
                cap and with the cap lifted, on identical seeds and identical orderings.
IDENTIFICATION  identified: the only thing that changes between arms is how many people are scored.
SCOPE           12,459 people with >=6 recorded onsets. n_eff is PEOPLE.
WORLDS          precision-only  location moves less than the uncapped seed spread; spreads shrink
                biased          location moves MORE than the uncapped spread -> every capped number
                                in five rounds is wrong, not merely imprecise
KILL            threshold-free: |location shift| compared against the UNCAPPED seed spread, per
                comparison, all comparisons published.
POSITIVE CTRL   the shrinkage itself must appear -- if uncapping does not reduce the seed spread,
                the cap was not the binding constraint and the diagnosis in #73 is wrong.
NEGATIVE CTRL   a random ordering, scored under both arms: must be at chance in both.
NOISE FLOOR     6 seeds per arm.
MULTIPLICITY    2 groupings x 2 sides x 3 orderings x 2 arms x 6 seeds, published whole.
IMPOSSIBLE      recovering the exact original RNG draw sequence -- the capped runs consumed the
                generator differently. Arms are matched on SEED, not on draw sequence, so a small
                residual difference is expected and is bounded by the reported spreads.
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'; OUT.mkdir(exist_ok=True)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values; mask=~np.isnan(V)
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
breadth=((df[rate].apply(pd.to_numeric,errors='coerce')>0).sum(1)).values
keep=np.flatnonzero(mask.sum(1)>=6)
GROUPS={'sex':(np.flatnonzero((male==1)&(mask.sum(1)>=6)),np.flatnonzero((male==0)&(mask.sum(1)>=6))),
        'breadth':(np.flatnonzero((breadth>np.median(breadth))&(mask.sum(1)>=6)),
                   np.flatnonzero((breadth<=np.median(breadth))&(mask.sum(1)>=6)))}
print(f"eligible people {len(keep):,}   sex groups {len(GROUPS['sex'][0]):,}/{len(GROUPS['sex'][1]):,}",flush=True)

def acc(order_vals,people,rng,cap):
    """VERBATIM the capped helper from R14-R18, with cap as a parameter. cap=None lifts it."""
    right=0;tot=0
    for i in people:
        j=np.flatnonzero(mask[i])
        if len(j)<2: continue
        for _ in range(min(10,len(j))):
            a,b=rng.choice(j,2,replace=False)
            if V[i,a]==V[i,b] or order_vals[a]==order_vals[b]: continue
            right+=((order_vals[a]<order_vals[b])==(V[i,a]<V[i,b])); tot+=1
            if cap and tot>=cap: break
        if cap and tot>=cap: break
    return 100*right/max(tot,1),tot

rows=[]
for arm,cap in [('capped',20000),('uncapped',None)]:
    for gname,(gA,gB) in GROUPS.items():
        for seed in range(1,7):
            rng=np.random.default_rng(seed)
            for lab,tgt,oth in [('A',gA,gB),('B',gB,gA)]:
                p=rng.permutation(tgt); tr,te=p[:len(p)//2],p[len(p)//2:]
                own=np.nanmean(V[tr],axis=0); other=np.nanmean(V[oth],axis=0)
                glob=np.nanmean(V[rng.permutation(keep)[:len(keep)//2]],axis=0)
                o,n=acc(own,te,rng,cap); ot,_=acc(other,te,rng,cap)
                gl,_=acc(glob,te,rng,cap); rd,_=acc(rng.permutation(own),te,rng,cap)
                rows.append(dict(arm=arm,grouping=gname,side=lab,seed=seed,own=o,other=ot,
                                 glob=gl,random=rd,n_pairs=n,n_people=len(te)))
    print(f"  {arm} done",flush=True)
D=pd.DataFrame(rows); D.to_csv(OUT/'grid.csv',index=False)

print("\n=== DID THE CAP EVEN BIND? pairs each round's design generates, vs the 20,000 ceiling ===")
for nm,pool in [('R14/R15/R16/R17  full pool',keep),
                ('R18  sex group A',GROUPS['sex'][0]),('R18  sex group B',GROUPS['sex'][1]),
                ('R18  breadth hi',GROUPS['breadth'][0]),('R18  breadth lo',GROUPS['breadth'][1])]:
    te=len(pool)//2
    est=sum(min(10,int(mask[i].sum())) for i in pool[:te])
    print(f"  {nm:24s} held-out people {te:6,d}   pairs ~{est:7,d}   "
          f"cap binds: {'YES  (%.0f%% of people dropped)'%(100*(1-20000/max(est,1))) if est>20000 else 'no'}")
print("\n=== HOW MUCH OF THE POPULATION EACH ARM ACTUALLY SCORED ===")
print(D.groupby('arm')[['n_pairs','n_people']].mean().round(0).to_string())

print("\n=== NEGATIVE CONTROL: a random ordering must sit at chance in BOTH arms ===")
print(D.groupby('arm')['random'].agg(['mean','std']).round(3).to_string())

print("\n=== LOCATION AND SPREAD, capped vs uncapped, matched seeds ===")
out=[]
for gname in GROUPS:
    for col in ['own','other','glob']:
        c=D[(D.arm=='capped')&(D.grouping==gname)][col]
        u=D[(D.arm=='uncapped')&(D.grouping==gname)][col]
        out.append(dict(grouping=gname,quantity=col,cap_mean=c.mean(),unc_mean=u.mean(),
                        loc_shift=u.mean()-c.mean(),cap_sd=c.std(),unc_sd=u.std(),
                        sd_ratio=c.std()/max(u.std(),1e-9)))
T=pd.DataFrame(out); T['biased']=T.loc_shift.abs()>T.unc_sd
print(T.round(3).to_string(index=False))

print("\n=== #62 RE-PRICED: own-group vs other-group ordering ===")
res=[]
for gname in GROUPS:
    for arm in ['capped','uncapped']:
        d=D[(D.arm==arm)&(D.grouping==gname)]
        diff=(d.own-d.other)
        res.append(dict(grouping=gname,arm=arm,diff=diff.mean(),spread=diff.std(),
                        ratio=abs(diff.mean())/max(diff.std(),1e-9),
                        own=d.own.mean(),other=d.other.mean()))
P=pd.DataFrame(res); P['resolvable']=P.ratio>2
print(P.round(3).to_string(index=False))

print("\n  CONDITIONAL KILL -- gates first")
g_rand=abs(D.random.mean()-50)<3
g_shrink=(T.sd_ratio.median()>1.5)
print(f"   (a) random ordering at chance in both arms : {'PASS' if g_rand else 'FAIL'} "
      f"({D.random.mean():.2f})")
print(f"   (b) uncapping shrinks the spread           : {'PASS' if g_shrink else 'FAIL -- the cap was not binding'}"
      f" (median sd ratio {T.sd_ratio.median():.2f}x)")
if not(g_rand and g_shrink):
    print("   -> UNVERIFIED.")
else:
    nb=int((~T.biased).sum())
    print(f"\n   location shifts smaller than the uncapped spread: {nb}/{len(T)} quantities")
    print(f"   median |shift| {T.loc_shift.abs().median():.3f}   median uncapped sd {T.unc_sd.median():.3f}")
    if nb==len(T):
        print("\n   -> PRECISION-ONLY, and now tested rather than assumed. Every point estimate in")
        print("      R14-R18 stands; every NULL and every resolvability verdict in them was computed")
        print(f"      against a spread {T.sd_ratio.median():.1f}x too large and must be re-read.")
    else:
        print(f"\n   -> BIASED in {len(T)-nb} quantities: the cap moved the ANSWER, not just its")
        print("      precision, and those numbers are wrong rather than imprecise.")
    for _,r in P[P.arm=='uncapped'].iterrows():
        was=P[(P.grouping==r.grouping)&(P.arm=='capped')].iloc[0]
        print(f"     {r.grouping:9s} own-other {r['diff']:+.2f} ratio {r.ratio:.2f} "
              f"({'RESOLVABLE' if r.resolvable else 'below 2x spread'})   "
              f"[capped: {was['diff']:+.2f} ratio {was.ratio:.2f}]")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
