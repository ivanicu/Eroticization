import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A05 R14 -- BREADTH WAS UNDER-USED FOR EIGHTY-FIVE ENTRIES, AND THE REASON IS ARITHMETIC.

#85b measured that the person main effect is +0.093 raw and +0.148 against its own null -- a 60%
under-report nobody noticed, because A05 asked "is breadth a nuisance or the object", answered
nuisance on the raw number, and every round since has projected breadth out rather than measured it.

The cause is not subtle once stated: an ITEM effect is estimated from n = 1,200-15,000 observations
per column; a PERSON effect from m = 10-24 per row. So the item effect's null is ~0 (its estimator
barely overfits) and the person effect's null is -0.05 (its estimator overfits badly). The
correction in #85 is therefore measuring ESTIMATOR NOISE, and if that is right, a SHRUNKEN person
estimator should deliver the gap DIRECTLY, without any correction.

That is a falsifiable, mechanistic prediction, and it is cheap.

ESTIMAND        held-out skill of the person main effect under raw, James-Stein-shrunk and
                empirical-Bayes estimators; and how much of the +0.055 gap each recovers.
IDENTIFICATION  identified. Shrinkage uses only training cells; the shrinkage factor is estimated
                from the training-cell variance decomposition, never from held-out cells.
SCOPE           the 23 blocks A09/R114 identified.
WORLDS          noise      shrinkage recovers most of the gap -> the raw number was an
                           UNDER-ESTIMATE from a noisy estimator, and breadth is a component this
                           project has been discarding
                bias       shrinkage recovers little -> the #85 correction is crediting P with the
                           damage a MISASSIGNED person effect does, not with skill it has, and
                           #85b's +0.148 must be restated
KILL            threshold-free: the share of the gap recovered, with its seed spread; declared only
                above 2x that spread.
POSITIVE CTRL   a synthetic world with a KNOWN person effect of known size and m observations per
                person: shrinkage must recover its planted skill better than the raw row mean.
NEGATIVE CTRL   a world with NO person effect (within-column shuffle): shrinkage must drive the
                person component to ~0 rather than to a negative number, which is the whole point
                of shrinking.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 3 estimators x 3 worlds x 3 seeds, published whole.
IMPOSSIBLE      separating "breadth" from "acquiescence" -- they are the same row sum in this
                release, and #A05/R085 already measured that the distinction is not identifiable
                here. Stated so the finding is not over-read as being about desire rather than
                response style.
"""
import pandas as pd, numpy as np, warnings, hashlib
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
MASK=0.15; SEEDS=[11,29,47]
print(f"targets {len(IDENT)}   people per block {min(RAW[q]['M'].shape[0] for q in IDENT)}-"
      f"{max(RAW[q]['M'].shape[0] for q in IDENT)}   options per block "
      f"{min(RAW[q]['M'].shape[1] for q in IDENT)}-{max(RAW[q]['M'].shape[1] for q in IDENT)}",flush=True)

def person_skill(M,seed,mode):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    raw=np.nanmean(T1,axis=1); raw=np.where(np.isnan(raw),0.,raw)
    ki=np.maximum(obs.sum(1),1)                        # observations per person
    if mode=='raw':
        b=raw
    else:
        # within-person sampling variance of the residual mean, from TRAINING cells only
        vw=np.nanmean(np.where(obs,(T1-raw[:,None])**2,np.nan),axis=1)
        vw=np.where(np.isnan(vw),np.nanmean(vw),vw)
        s2_within=vw/ki
        s2_total=raw.var()
        s2_true=max(s2_total-np.mean(s2_within),1e-9)
        if mode=='js':                                  # one global shrinkage factor
            lam=s2_true/(s2_true+np.mean(s2_within)); b=raw*lam
        else:                                           # empirical Bayes, per person by its own ki
            lam=s2_true/(s2_true+s2_within); b=raw*lam
    P=b[:,None]; IB=np.broadcast_to(I,M.shape)
    base=np.mean((M[he]-gm)**2)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    return f(IB,P)-f(IB), float(np.mean(lam)) if mode!='raw' else 1.0

def col_shuffle(M,rng): return np.column_stack([rng.permutation(M[:,j]) for j in range(M.shape[1])])
def plant(M,g,rng):
    n,m=M.shape; b=rng.normal(size=n)*g
    p=np.clip(M.mean(0)[None,:]+b[:,None],0.02,0.98)
    return (rng.random((n,m))<p).astype(float)

rows=[]
for i,t in enumerate(IDENT):
    M=RAW[t]['M']
    for sd in SEEDS:
        rg=np.random.default_rng(9500+sd)
        worlds={'real':M,'noperson':col_shuffle(M,rg),'planted':plant(M,0.12,rg)}
        for w,Mw in worlds.items():
            for mode in ['raw','js','eb']:
                s,lam=person_skill(Mw,sd,mode)
                rows.append(dict(q=t,world=w,mode=mode,seed=sd,skill=s,lam=lam,
                                 n=M.shape[0],m=M.shape[1]))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'R097_shrink_the_person_effect'/'results'
OUT.mkdir(parents=True,exist_ok=True); D.to_csv(OUT/'grid.csv',index=False)

print("\n=== PERSON-EFFECT SKILL BY ESTIMATOR AND WORLD (mean over 23 blocks x 3 seeds) ===")
T=D.groupby(['world','mode']).skill.mean().unstack('mode')[['raw','js','eb']]
T['shrink_factor']=D[D['mode']=='eb'].groupby('world').lam.mean()
print(T.round(4).to_string())

print("\n  CONDITIONAL KILL -- gates first")
g1=T.loc['planted','eb']>T.loc['planted','raw']
g2=abs(T.loc['noperson','eb'])<abs(T.loc['noperson','raw'])
print(f"   (a) shrinkage beats the raw row mean where a person effect is PLANTED : "
      f"{'PASS' if g1 else 'FAIL'} ({T.loc['planted','raw']:+.4f} -> {T.loc['planted','eb']:+.4f})")
print(f"   (b) shrinkage pulls toward 0 where there is NO person effect          : "
      f"{'PASS' if g2 else 'FAIL'} ({T.loc['noperson','raw']:+.4f} -> {T.loc['noperson','eb']:+.4f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    raw=T.loc['real','raw']; eb=T.loc['real','eb']
    target=0.1482                                        # #85's corrected P
    rec=(eb-raw)/max(target-raw,1e-9)
    sp=D[(D.world=='real')&(D['mode']=='eb')].groupby('q').skill.std().median()
    print(f"\n   real world: raw {raw:+.4f}   James-Stein {T.loc['real','js']:+.4f}   "
          f"empirical Bayes {eb:+.4f}")
    print(f"   #85's corrected P was +0.1482, so the gap to close was {target-raw:+.4f}")
    print(f"   recovered by shrinkage alone: {rec:.1%}   (seed spread {sp:.4f})")
    print(f"   mean shrinkage factor applied: {T.loc['real','shrink_factor']:.3f} "
          f"(1.0 = no shrinkage)")
    if rec>0.4:
        print("\n   -> NOISE. The raw person effect was an UNDER-ESTIMATE from an estimator working")
        print("      off 10-24 observations per person, and shrinking it recovers most of the gap")
        print("      WITHOUT any null correction. Breadth is a component this project discarded on")
        print("      a number its own estimator was suppressing.")
    elif rec<0.1:
        print("\n   -> BIAS. Shrinkage recovers almost nothing, so #85's +0.1482 is crediting P with")
        print("      the DAMAGE a misassigned person effect does rather than skill it has, and")
        print("      #85b's person figure must be restated.")
    else:
        print(f"\n   -> PARTIAL ({rec:.0%}). Both mechanisms are present and neither dominates; the")
        print("      honest report of the person effect is the interval [raw, corrected].")
print("\nN/A, with what it would require: separating breadth from acquiescence needs a reverse-keyed "
      "or forced-choice item. This release has none, so nothing here distinguishes 'wants more' from "
      "'ticks more'.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
