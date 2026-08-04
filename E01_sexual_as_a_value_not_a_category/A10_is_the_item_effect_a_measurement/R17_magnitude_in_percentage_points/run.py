import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A10 R17 -- WHAT THE THREE COMPONENTS ARE, IN PERCENTAGE POINTS.

Every number in this project is held-out R2 or Shapley skill -- quantities whose meaning depends on
the estimator that produced them, which is exactly what #86, #87, #88 and #89 spent five entries
discovering. #89d left the inversion open: a planted rank-5 world at loading scale 0.30 yields skill
+0.166 and the real data yields +0.035, so what magnitude does the real interaction correspond to?

Inverting the plant converts skill into a PROBABILITY-SCALE effect: how many percentage points does
a person's endorsement probability for a given option move away from that option's base rate?

The inversion is one-to-many -- skill is one number and (rank, loading scale) is two -- so the answer
is a FAMILY of worlds, not a point. But if the family lies along constant s*sqrt(r), the implied
per-cell probability sd is identified even though rank and scale separately are not. That is the
result this round is testing for, and it is stated before the run rather than after.

ESTIMAND        the (rank, scale) family of synthetic worlds whose soft-thresholded W skill matches
                the real data's, and the implied per-cell probability sd along it; plus the same
                quantity read directly for the item and person components.
IDENTIFICATION  the family is identified; individual (rank, scale) pairs are NOT, and the round
                reports the invariant along the family rather than a point on it.
SCOPE           the 23 blocks A09/R03 identified. Soft thresholding throughout (#88's estimator).
WORLDS          invariant  s*sqrt(r) is constant along the matching family -> the magnitude is
                           identified in percentage points
                not        different (r,s) with the same skill imply different magnitudes -> only
                           the skill is reportable and this conversion is impossible here
KILL            threshold-free: the spread of implied per-cell sd along the matching family, relative
                to its own mean, IS the verdict.
POSITIVE CTRL   the implied sd recovered for a planted world must match the sd actually planted.
NEGATIVE CTRL   the fixed-margin world must invert to ~0 percentage points.
NOISE FLOOR     2 masks x 2 draws per cell.
MULTIPLICITY    3 ranks x 5 scales x 23 blocks x 2 seeds, published whole.
IMPOSSIBLE      a per-PERSON magnitude -- this is a population sd, and #A05/R02 measured that breadth
                and acquiescence are the same row sum here, so nothing distinguishes "wants" from
                "ticks" at the individual level.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd
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
    RAW[q.qi]=M
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
MASK=0.15; SEEDS=[11,29]; LAMS=[0.5,1.,2.,4.,8.]
RANKS=[2,5,10]; SCALES=[0.05,0.08,0.12,0.20,0.30]
print(f"targets {len(IDENT)}",flush=True)
def curveball(M,rng,per_row=5.):
    A=[set(np.flatnonzero(r).tolist()) for r in M]; n=len(A)
    for _ in range(int(per_row*n)):
        i,j=int(rng.integers(n)),int(rng.integers(n))
        if i==j: continue
        ai,aj=A[i],A[j]; inter=ai&aj
        di=list(ai-inter); dj=list(aj-inter); L=di+dj
        if not L: continue
        rng.shuffle(L); k=len(di)
        A[i]=inter|set(L[:k]); A[j]=inter|set(L[k:])
    out=np.zeros_like(M)
    for i,s in enumerate(A): out[i,list(s)]=1.
    return out
def plant(M,r,sc,rng):
    n,m=M.shape; F=rng.normal(size=(n,r)); L=rng.normal(size=(r,m))*sc
    dev=F@L
    base=M.mean(0)[None,:]
    p=np.clip(base+dev,0.02,0.98)
    # ⚠ #157:第一版返回 `dev.std()` —— **clip 之前**的种植扰动 sd。
    #   `#90c` 已经把它判为「第十七个 mis-specified statistic」,并把校正后的数
    #   (交互 ±30.8 -> ±23.7 pp,族 CV 22.4% -> 15.9%)写进了账本 —— **但代码从没改**。
    #   于是这一轮至今仍在打印被自己账本判为错的那一版,而且**静默地**(exit 0)。
    return (rng.random((n,m))<p).astype(float), float((p-base).std())
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
    M=RAW[t]
    for sd in SEEDS:
        rows.append(dict(q=t,world='real',rank=np.nan,scale=np.nan,seed=sd,
                         skill=wskill(M,sd),planted_sd=np.nan,
                         item_pp=100*M.mean(0).std(),person_pp=100*M.mean(1).std(),
                         # #157:同一循环里记下二项噪声,#90d 的校正才有可用的分母
                         person_noise_pp=100*float(np.sqrt(np.mean(
                             M.mean(1)*(1-M.mean(1))/M.shape[1])))))
        rg=np.random.default_rng(3300+sd)
        rows.append(dict(q=t,world='margin',rank=np.nan,scale=np.nan,seed=sd,
                         skill=wskill(curveball(M,rg),sd),planted_sd=0.,
                         item_pp=np.nan,person_pp=np.nan,person_noise_pp=np.nan))
        for r in RANKS:
            for sc in SCALES:
                Mw,psd=plant(M,r,sc,np.random.default_rng(4400+sd))
                rows.append(dict(q=t,world='plant',rank=r,scale=sc,seed=sd,
                                 skill=wskill(Mw,sd),planted_sd=psd,
                                 item_pp=np.nan,person_pp=np.nan,person_noise_pp=np.nan))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
real=D[D.world=='real'].skill.mean(); marg=D[D.world=='margin'].skill.mean()
print(f"\n=== real W skill {real:+.4f}   fixed-margin {marg:+.4f} ===")
P=D[D.world=='plant'].groupby(['rank','scale']).agg(skill=('skill','mean'),
                                                    planted_pp=('planted_sd',lambda s:100*s.mean()))
print("\n=== PLANT LADDER: skill and the per-cell probability sd actually planted ===")
print(P.round(4).to_string())
print("\n=== INVERSION: which planted magnitude reproduces the real skill? ===")
fam=[]
for r in RANKS:
    sub=P.loc[r].sort_values('skill')
    if sub.skill.min()<=real<=sub.skill.max():
        pp=float(np.interp(real,sub.skill.values,sub.planted_pp.values))
        sc=float(np.interp(real,sub.skill.values,np.array(sub.index,dtype=float)))
        fam.append(dict(rank=r,scale=sc,implied_pp=pp))
F=pd.DataFrame(fam)
print(F.round(3).to_string(index=False) if len(F) else "  real skill outside the ladder's range")
print("\n  CONDITIONAL KILL -- gates first")
g1=abs(marg)<0.03
g2=len(F)>=2
print(f"   (a) fixed-margin world inverts to ~0 : {'PASS' if g1 else 'FAIL'} ({marg:+.4f})")
print(f"   (b) real skill inside the ladder for >=2 ranks : {'PASS' if g2 else 'FAIL'}")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    cv=F.implied_pp.std()/F.implied_pp.mean()
    print(f"\n   implied per-cell probability sd along the family: "
          f"{F.implied_pp.round(2).tolist()} pp   CV {cv:.1%}")
    ip=D[D.world=='real'].item_pp.mean(); pp_raw=D[D.world=='real'].person_pp.mean()
    # ⚠ #157:`#90d` 的校正也只活在账本散文里 —— 观测到的人层展布里有一部分是二项噪声,
    #   校正量 = sqrt(观测^2 - 噪声^2)。代码一直打印**未校正**的那个。
    noise=float(D[D.world=='real'].person_noise_pp.mean())
    pp=float(np.sqrt(max(pp_raw**2-noise**2,0.)))
    print(f"\n   人层展布的二项噪声校正(#90d):观测 {pp_raw:.1f} pp,噪声 {noise:.1f} pp "
          f"-> 校正后 {pp:.1f} pp")
    print(f"\n   THE THREE COMPONENTS IN PERCENTAGE POINTS")
    print(f"     option base rates differ by      +/- {ip:.1f} pp  (sd of column means)")
    print(f"     people differ in overall rate by +/- {pp:.1f} pp  (sd of row means)")
    print(f"     person x option interaction      +/- {F.implied_pp.mean():.1f} pp  (inverted plant)")
    if cv<0.25:
        print(f"\n   -> IDENTIFIED. The implied magnitude varies only {cv:.0%} across ranks 2-10, so the")
        print(f"      percentage-point statement holds even though rank and scale separately do not.")
    else:
        print(f"\n   -> NOT IDENTIFIED: {cv:.0%} spread across the family. Only the skill is reportable.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
