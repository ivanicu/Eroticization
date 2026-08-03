import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A10 R18 -- TESTING #90's OWN SCOPE LIMITATION.

#90 converted the interaction from a skill number (+0.036) into a magnitude (+/-23.7 pp) by inverting
a DENSE Gaussian low-rank plant, and named the assumption in its own scope paragraph: a SPARSE
structure -- a few people with strong specific tastes rather than everyone deviating a little --
would produce the same skill at a different sd.

That is a falsifiable statement about my own headline, and it costs one plant ladder.

  dense   every person deviates (f = 1.0), the #90 specification
  sparse  only a fraction f of people deviate at all; the rest sit exactly at the base rate

Each shape is calibrated to the SAME observed skill, and the question is whether they imply the same
realized per-cell probability sd.

ESTIMAND        implied realized (post-clip) per-cell sd at matched skill, as a function of the
                fraction of people carrying the structure.
IDENTIFICATION  identified per shape; the point is precisely that the shape is NOT identified by
                skill alone, and this measures how much that costs.
SCOPE           the 23 blocks A09/R03 identified, soft-thresholded estimator (#88).
WORLDS          shape-robust  implied sd varies little with f -> #90a's +/-23.7 pp stands unqualified
                shape-bound   implied sd varies strongly -> #90a is CONDITIONAL on density and must
                              be reported as a range indexed by an unmeasurable shape parameter
KILL            threshold-free: the spread of implied sd across f, relative to its mean, IS the answer.
POSITIVE CTRL   the dense arm (f = 1.0) must reproduce #90's +/-23.7 pp, or this pipeline is not
                measuring the same thing #90 measured.
NEGATIVE CTRL   scale 0 at every f must give skill ~0 and implied sd 0.
NOISE FLOOR     2 masks.
MULTIPLICITY    4 densities x 6 scales x 23 blocks x 2 seeds, published whole.
IMPOSSIBLE      measuring the real structure's density directly -- that is the whole difficulty; if
                density were observable the inversion would not be needed.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd
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
    RAW[q.qi]=M
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
MASK=0.15; SEEDS=[11,29]; LAMS=[0.5,1.,2.,4.,8.]
DENS=[0.05,0.15,0.30,1.0]; SCALES=[0.,0.08,0.12,0.20,0.30,0.50]
REAL_SKILL=0.0362
print(f"targets {len(IDENT)}",flush=True)
def wskill(M,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    val=(rng.random(M.shape)<0.2)&obs; fit=obs&~val
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I; rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    R=T1-P; Rf=np.where(np.isnan(R),0.,R)
    base=np.mean((M[he]-gm)**2); IB=np.broadcast_to(I,M.shape)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    b0=f(IB,P); best=None
    for lam in LAMS:
        F=np.where(fit,Rf,0.)
        for _ in range(20):
            U,S,V=svd(F,full_matrices=False); F=np.where(fit,Rf,(U*np.maximum(S-lam,0.))@V)
        U,S,V=svd(F,full_matrices=False); W=(U*np.maximum(S-lam,0.))@V
        e=np.mean((Rf[val]-W[val])**2)
        if best is None or e<best[0]: best=(e,f(IB,P,W)-b0)
    return best[1]
rows=[]
for i,t in enumerate(IDENT):
    M=RAW[t]; n,m=M.shape; base=M.mean(0)[None,:]
    for sd in SEEDS:
        for dens in DENS:
            for sc in SCALES:
                rg=np.random.default_rng(5500+sd)
                F=rg.normal(size=(n,5)); L=rg.normal(size=(5,m))*sc
                dev=F@L
                if dens<1.0:
                    on=rg.random(n)<dens
                    dev=dev*on[:,None]/np.sqrt(max(dens,1e-9))   # same total variance at same sc
                p=np.clip(base+dev,0.02,0.98); realized=p-base
                Mw=(rg.random((n,m))<p).astype(float)
                rows.append(dict(q=t,dens=dens,scale=sc,seed=sd,skill=wskill(Mw,sd),
                                 pp=100*realized.std(),
                                 pp_carriers=100*(realized[on].std() if dens<1.0 else realized.std())))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby(['dens','scale']).agg(skill=('skill','mean'),pp=('pp','mean'),ppc=('pp_carriers','mean'))
print("\n=== SKILL AND REALIZED MAGNITUDE BY DENSITY AND SCALE ===")
print(G.round(3).to_string())
print(f"\n=== INVERSION AT THE REAL SKILL ({REAL_SKILL:+.4f}) ===")
fam=[]
for d_ in DENS:
    s=G.loc[d_].sort_values('skill')
    if s.skill.min()<=REAL_SKILL<=s.skill.max():
        fam.append(dict(density=d_,implied_pp=float(np.interp(REAL_SKILL,s.skill.values,s.pp.values)),
                        implied_pp_carriers=float(np.interp(REAL_SKILL,s.skill.values,s.ppc.values))))
F=pd.DataFrame(fam); print(F.round(2).to_string(index=False) if len(F) else "  outside the ladder")
print("\n  CONDITIONAL KILL -- gates first")
g1=abs(G.loc[(1.0,0.0),'skill'])<0.02
g2=len(F)>=3
dense=F[F.density==1.0].implied_pp.iloc[0] if (len(F) and (F.density==1.0).any()) else np.nan
g3=abs(dense-23.7)<6 if not np.isnan(dense) else False
print(f"   (a) scale 0 gives skill ~0            : {'PASS' if g1 else 'FAIL'} ({G.loc[(1.0,0.0),'skill']:+.4f})")
print(f"   (b) real skill inside >=3 densities   : {'PASS' if g2 else 'FAIL'}")
print(f"   (c) dense arm reproduces #90's 23.7pp : {'PASS' if g3 else 'FAIL'} ({dense:.1f} pp)")
if not(g1 and g2 and g3): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    cv=F.implied_pp.std()/F.implied_pp.mean()
    print(f"\n   POPULATION-WIDE implied sd across densities {F.implied_pp.round(1).tolist()} pp  CV {cv:.1%}")
    print(f"   AMONG CARRIERS ONLY                        {F.implied_pp_carriers.round(1).tolist()} pp")
    if cv<0.25:
        print(f"\n   -> SHAPE-ROBUST. #90a's +/-{F.implied_pp.mean():.1f} pp holds whether 5% of people carry")
        print("      the structure or all of them, because the POPULATION sd is what skill pins down.")
        print(f"      What shape DOES change is the per-carrier magnitude: "
              f"{F.implied_pp_carriers.min():.0f}-{F.implied_pp_carriers.max():.0f} pp.")
    else:
        print(f"\n   -> SHAPE-BOUND ({cv:.0%} spread). #90a is CONDITIONAL on density and must be")
        print("      reported as a range indexed by a shape parameter this release cannot measure.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
