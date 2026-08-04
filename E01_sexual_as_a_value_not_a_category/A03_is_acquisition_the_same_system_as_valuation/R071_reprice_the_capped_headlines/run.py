import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A03 R22 -- RE-PRICING THE THREE ROUNDS WHERE THE CAP ACTUALLY BOUND.

#74 measured which designs the cap reached: the four FULL-POOL rounds, not R18. R17 was re-run in
#73. This round re-prices the remaining three, and R14 among them carries the project's strongest
surviving claim.

Random truncation leaves LOCATION unbiased and inflates SPREAD, so every downstream resolvability
verdict was computed against a denominator ~1.7x too large -- i.e. CONSERVATIVE. The point of
re-pricing is therefore not to defend those verdicts but to find the ones that were called
UNRESOLVABLE only because the cap made the ruler coarse. A verdict that flips here was never a fact
about the world.

ESTIMAND        for R14 (population / oracle / random ordering), R15 (onset / prevalence /
                onset-minus-prevalence) and R16 (onset vs prevalence within onset bands): the point
                estimate and seed spread, capped vs uncapped, matched seeds.
IDENTIFICATION  identified; the only difference between arms is how many people are scored.
SCOPE           12,459 people with >=6 onsets. n_eff is PEOPLE. Bands in R16 are onset strata.
WORLDS          conservative-only  spreads shrink, no verdict flips -> the cap cost precision only
                verdict-flipping   a comparison previously below 2x spread crosses it -> that
                                   verdict was an artefact of the ruler, not a finding
KILL            threshold-free: each comparison's |effect| / seed spread reported in both arms, and
                a flip is declared only when the uncapped ratio crosses 2 while the capped did not.
POSITIVE CTRL   the spread must actually shrink here (unlike R18 in #74, where it could not). If it
                does not, the cap was not binding and #74's measurement is wrong.
NEGATIVE CTRL   a random ordering: at chance in both arms, in every stratum.
NOISE FLOOR     8 seeds per arm.
MULTIPLICITY    3 rounds x their comparisons x 2 arms x 8 seeds, published whole.
IMPOSSIBLE      unchanged -- no per-person onset ceiling exists in this release.
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
prev=mask.mean(0)
keep=np.flatnonzero(mask.sum(1)>=6)
print(f"eligible {len(keep):,}  categories {V.shape[1]}",flush=True)

def acc(order_vals,people,rng,cap):
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
    for seed in range(1,9):
        rng=np.random.default_rng(seed)
        p=rng.permutation(keep); tr,te=p[:len(p)//2],p[len(p)//2:]
        pop=np.nanmean(V[tr],axis=0); orc=np.nanmean(V[te],axis=0)
        X=np.c_[np.ones(len(prev)),prev]
        resid=pop-X@np.linalg.lstsq(X,pop,rcond=None)[0]      # onset with prevalence projected out
        r14p,n=acc(pop,te,rng,cap); r14o,_=acc(orc,te,rng,cap); r14r,_=acc(rng.permutation(pop),te,rng,cap)
        r15pr,_=acc(-prev,te,rng,cap); r15rs,_=acc(resid,te,rng,cap)
        rows.append(dict(arm=arm,seed=seed,band='all',n_people=len(te),pairs=n,
                         population=r14p,oracle=r14o,random=r14r,prevalence=r15pr,onset_minus_prev=r15rs))
        # R16: within onset bands (censoring strata) -- people binned by their own mean onset
        pm=np.nanmean(np.where(mask,V,np.nan),axis=1)
        for bi,(lo,hi) in enumerate([(0,11),(11,14),(14,18),(18,99)]):
            idx=te[(pm[te]>=lo)&(pm[te]<hi)]
            if len(idx)<200: continue
            a_on,nb=acc(pop,idx,rng,cap); a_pr,_=acc(-prev,idx,rng,cap)
            a_rd,_=acc(rng.permutation(pop),idx,rng,cap)
            rows.append(dict(arm=arm,seed=seed,band=f"onset {lo}-{hi}",n_people=len(idx),pairs=nb,
                             population=a_on,oracle=np.nan,random=a_rd,prevalence=a_pr,
                             onset_minus_prev=np.nan))
    print(f"  {arm} done",flush=True)
D=pd.DataFrame(rows); D.to_csv(OUT/'grid.csv',index=False)

A=D[D.band=='all']
print("\n=== DID THE SPREAD SHRINK? (positive control for #74's measurement) ===")
print(A.groupby('arm')[['pairs','n_people']].mean().round(0).to_string())
sp=A.groupby('arm')[['population','oracle','random','prevalence','onset_minus_prev']].std()
print(sp.round(3).to_string())
ratio=(sp.loc['capped']/sp.loc['uncapped']).median()
print(f"  median sd ratio capped/uncapped = {ratio:.2f}x  -> "
      f"{'PASS, the cap was binding here' if ratio>1.3 else 'FAIL -- #74 mis-measured'}")

print("\n=== NEGATIVE CONTROL: random ordering, both arms, every stratum ===")
print(D.groupby('arm')['random'].agg(['mean','std']).round(2).to_string())

print("\n=== R14 / R15 HEADLINES, capped vs uncapped ===")
res=[]
for col in ['population','oracle','prevalence','onset_minus_prev']:
    for arm in ['capped','uncapped']:
        d=A[A.arm==arm]
        res.append(dict(quantity=col,arm=arm,mean=d[col].mean(),sd=d[col].std(),
                        vs_chance=d[col].mean()-50,ratio=(d[col].mean()-50)/max(d[col].std(),1e-9)))
H=pd.DataFrame(res); H['resolvable']=H.ratio>2
print(H.round(3).to_string(index=False))

print("\n=== R16: onset vs prevalence WITHIN onset bands (the censoring test) ===")
res=[]
for b in sorted(D.band.unique()):
    if b=='all': continue
    for arm in ['capped','uncapped']:
        d=D[(D.band==b)&(D.arm==arm)]
        if not len(d): continue
        gap=d.population-d.prevalence
        res.append(dict(band=b,arm=arm,n=d.n_people.mean(),onset=d.population.mean(),
                        prevalence=d.prevalence.mean(),gap=gap.mean(),sd=gap.std(),
                        ratio=abs(gap.mean())/max(gap.std(),1e-9)))
B=pd.DataFrame(res); B['resolvable']=B.ratio>2
print(B.round(3).to_string(index=False))

print("\n  CONDITIONAL KILL -- gate first")
g1=ratio>1.3; g2=abs(D.random.mean()-50)<3
print(f"   (a) spread shrinks when uncapped : {'PASS' if g1 else 'FAIL'} ({ratio:.2f}x)")
print(f"   (b) random ordering at chance     : {'PASS' if g2 else 'FAIL'} ({D.random.mean():.2f})")
if not(g1 and g2): print("   -> UNVERIFIED.")
else:
    flips=[]
    for df_,key in [(H,'quantity'),(B,'band')]:
        for k in df_[key].unique():
            c=df_[(df_[key]==k)&(df_.arm=='capped')]; u=df_[(df_[key]==k)&(df_.arm=='uncapped')]
            if len(c) and len(u) and (not c.resolvable.iloc[0]) and u.resolvable.iloc[0]:
                flips.append((k,float(c.ratio.iloc[0]),float(u.ratio.iloc[0])))
    print(f"\n   verdicts that FLIP from unresolvable to resolvable when uncapped: {len(flips)}")
    for k,a,b_ in flips: print(f"     {k}: ratio {a:.2f} -> {b_:.2f}")
    if not flips:
        print("   -> CONSERVATIVE-ONLY. The cap cost precision and cost NO verdict. Every")
        print("      unresolvable call in R14/R15/R16 was unresolvable for a reason other than the cap.")
    else:
        print("   -> VERDICTS FLIPPED. Those calls were artefacts of a coarse ruler, not findings.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
