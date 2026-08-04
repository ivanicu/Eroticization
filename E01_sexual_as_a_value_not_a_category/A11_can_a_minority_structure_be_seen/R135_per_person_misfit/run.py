import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R01 -- THE FIRST METHOD IN THIS PROJECT THAT COULD SEE A MINORITY.

#91 measured a capability boundary: a person x option structure carried by 5% of people is
undetectable at +/-50 percentage points, because every method used in eleven arcs is a
VARIANCE-EXPLAINED method and those weight by prevalence. A fetish -- a minority with intense
specific attachments -- is exactly that object.

The fix is not a better factorisation. It is a statistic that weights every PERSON equally and needs
no per-person parameters, because a person contributes only ~18 observations and nothing can be fit
from them. A per-person GOODNESS-OF-FIT test does both:

  p_hat_ij   = clip(grand mean + item effect + SHRUNK person effect), from TRAINING cells only
  z_ij       = (M_ij - p_hat_ij) / sqrt(p_hat(1-p_hat))          on HELD-OUT cells only
  T_i        = mean of z^2 over every held-out cell that person has, pooled across ALL their blocks

Under a population model with no personal structure, T_i concentrates around 1. A person carrying
specific structure has T_i inflated. The signal is then in the UPPER TAIL of the T distribution, and
a minority shows up there while contributing almost nothing to any variance-explained measure.

ESTIMAND        the upper quantiles of the per-person misfit distribution T, relative to a
                parametric null generated from the fitted population model itself.
IDENTIFICATION  identified; no per-person parameter is estimated, so there is nothing to overfit.
SCOPE           people appearing in >=4 of the 23 identified blocks, so T has enough held-out cells.
WORLDS          minority-present  real upper tail exceeds the null's -> a structure exists that
                                  every previous round was blind to
                none              real tail matches the null -> the blind spot is real but this
                                  release has nothing hiding in it
KILL            threshold-free: the real-vs-null gap at each quantile against its own seed spread.
POSITIVE CTRL   THE EXACT WORLD #91 SHOWED WAS INVISIBLE -- 5% of people carrying rank-5 structure at
                scale 0.50. If this statistic cannot see that either, the method is not the fix and
                the boundary stands.
NEGATIVE CTRL   the parametric null itself, regenerated from p_hat: T must concentrate at 1.
NOISE FLOOR     3 masks x 3 null draws.
MULTIPLICITY    quantiles {50,75,90,95,99} x 3 worlds x 3 seeds, published whole.
IMPOSSIBLE      attributing an inflated T_i to a CAUSE -- misfit is misfit, and a careless responder
                and an intense specific interest both produce it. Stated up front.
"""
import pandas as pd, numpy as np, warnings, hashlib
warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
MASK=0.15; SEEDS=[11,29,47]; QS=[50,75,90,95,99]
print(f"targets {len(IDENT)}",flush=True)
def shrunk_rows(T1,obs):
    rmn=np.nanmean(T1,axis=1); rmn=np.where(np.isnan(rmn),0.,rmn)
    k=np.maximum(obs.sum(1),1)
    vw=np.nanmean(np.where(obs,(T1-rmn[:,None])**2,np.nan),axis=1)
    vw=np.where(np.isnan(vw),np.nanmean(vw),vw); s2w=vw/k
    s2t=max(np.var(rmn)-np.mean(s2w),1e-9)
    return rmn*(s2t/(s2t+s2w))
def per_person(M,seed,regen=False):
    """returns (person_index_in_block, sum z^2, count) on held-out cells."""
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    P=shrunk_rows(T-gm-I,obs)[:,None]
    ph=np.clip(gm+I+P,0.02,0.98)
    Y=M
    if regen: Y=(np.random.default_rng(seed+31).random(M.shape)<ph).astype(float)
    z2=(Y-ph)**2/(ph*(1-ph))
    return np.where(he,z2,0.).sum(1), he.sum(1)
def plant_sparse(M,dens,sc,rng):
    n,m=M.shape; F=rng.normal(size=(n,5)); L=rng.normal(size=(5,m))*sc
    dev=F@L; on=rng.random(n)<dens
    dev=dev*on[:,None]/np.sqrt(max(dens,1e-9))
    p=np.clip(M.mean(0)[None,:]+dev,0.02,0.98)
    return (rng.random((n,m))<p).astype(float)
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
rows=[]
for sd in SEEDS:
    acc={w:(np.zeros(len(ALLP)),np.zeros(len(ALLP))) for w in ['real','null','sparse5','sparse5null']}
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        Msp=plant_sparse(M,0.05,0.50,np.random.default_rng(6600+sd))
        for w,(Mw,regen) in {'real':(M,False),'null':(M,True),
                             'sparse5':(Msp,False),'sparse5null':(Msp,True)}.items():
            s,c=per_person(Mw,sd,regen); a,b=acc[w]; a[idx]+=s; b[idx]+=c
    ok=acc['real'][1]>=40
    for w in acc:
        a,b=acc[w]; Tv=a[ok]/np.maximum(b[ok],1)
        for qv in QS: rows.append(dict(seed=sd,world=w,q=qv,T=float(np.percentile(Tv,qv))))
        rows.append(dict(seed=sd,world=w,q=-1,T=float(Tv.mean())))
        rows.append(dict(seed=sd,world=w,q=-2,T=float((Tv>2.0).mean())))
    print(f"  seed {sd} done  people with >=40 held-out cells: {int(ok.sum()):,}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
P=D.groupby(['world','q']).T.agg(['mean','std']).unstack('world')['mean']
S=D.groupby(['world','q']).T.std().unstack('world')
lab={-2:'share T>2',-1:'mean T',50:'p50',75:'p75',90:'p90',95:'p95',99:'p99'}
print("\n=== PER-PERSON MISFIT T, BY QUANTILE ===")
out=P.rename(index=lab)[['real','null','sparse5','sparse5null']]
print(out.round(4).to_string())
print("\n  CONDITIONAL KILL -- gates first")
g1=abs(P.loc[50,'null']-1)<0.15
g2=(P.loc[99,'sparse5']-P.loc[99,'sparse5null'])>2*S.loc[99,'sparse5']
print(f"   (a) parametric null concentrates at 1 : {'PASS' if g1 else 'FAIL'} (median {P.loc[50,'null']:.3f})")
print(f"   (b) SEES the world #91 showed was invisible (5% carriers, +/-50pp): "
      f"{'PASS' if g2 else 'FAIL'} (p99 {P.loc[99,'sparse5']:.3f} vs its null {P.loc[99,'sparse5null']:.3f})")
if not(g1 and g2): print("   -> UNVERIFIED: the method is not the fix, and #91's boundary stands.")
else:
    print(f"\n   THE METHOD WORKS ON THE OBJECT THE OLD ONE COULD NOT SEE.")
    print(f"   real vs its own null, by quantile:")
    for qv in QS:
        d=P.loc[qv,'real']-P.loc[qv,'null']
        print(f"     p{qv:<3d} real {P.loc[qv,'real']:.4f}  null {P.loc[qv,'null']:.4f}  "
              f"diff {d:+.4f}  ({'RESOLVABLE' if abs(d)>2*S.loc[qv,'real'] else 'within 2x spread'})")
    d99=P.loc[99,'real']-P.loc[99,'null']
    if d99>2*S.loc[99,'real']:
        print(f"\n   -> A MINORITY STRUCTURE IS PRESENT. The upper tail of per-person misfit exceeds")
        print(f"      what the population model can generate, by {d99:+.4f} at p99.")
    else:
        print(f"\n   -> NO MINORITY STRUCTURE DETECTED, and this time the null is admissible: the")
        print(f"      instrument demonstrably sees a 5%-carrier structure (gate b). #91's blind spot")
        print(f"      is real as a property of the OLD methods, and this release has nothing hiding")
        print(f"      in that particular place.")
print("\nN/A: an inflated T_i cannot be attributed to a cause. A careless responder and an intense "
      "specific interest produce the same misfit.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
