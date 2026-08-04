import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A10 R16 -- DIMENSIONALITY, RE-ASKED ON THE ESTIMATOR THAT SURVIVED #88.

#72 ("no cliff in the cross-block spectrum") and #83 ("knee at 5 within blocks") were both measured
with HARD rank truncation, whose floor #88 measured at -0.09 on a structureless world -- an artifact,
not a floor. Soft singular-value thresholding has a floor of -0.02 and beats hard everywhere the
comparison is fair (#87a, settled in #88).

Under soft thresholding, rank is not chosen -- it EMERGES: the effective rank is the number of
singular values that survive the threshold, and the threshold is tuned on training cells. So the
question "how many dimensions" becomes a measurement rather than a sweep, and it can be calibrated
the same way #83 calibrated the knee: against worlds whose true rank is known.

ESTIMAND        effective rank = #{singular values > lambda} at the training-tuned lambda, per block;
                and the shape of the skill-vs-lambda curve, whose flatness is itself the answer if
                no optimum is identified.
IDENTIFICATION  identified relative to the controls. A bare effective rank means nothing; an
                effective rank that reads 2 on a rank-2 world and 5 on a rank-5 world means something.
SCOPE           the 23 blocks A09/R114 identified. Within-block structure only -- the cross-block
                analogue needs the wide basis and #84 showed that comparison is unidentifiable here.
WORLDS          low-rank    real effective rank lands near the controls' -> #83's knee at 5 survives
                            the method change
                high-rank   real effective rank far above -> #83 was an artifact of hard truncation
                unidentified the skill-vs-lambda curve is flat -> no rank is estimable and both #72
                            and #83 were reading structure into an estimator's tuning parameter
KILL            threshold-free: the real effective rank is reported beside the two controls' and the
                whole lambda curve is published, including its flatness.
POSITIVE CTRL   known within-block rank 2 and 5; recovered effective rank must be monotone in true
                rank and near it.
NEGATIVE CTRL   fixed-margin curveball: effective rank must collapse toward 0.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 4 worlds x 8 lambdas x 3 seeds, published whole.
IMPOSSIBLE      the cross-block analogue at matched budget -- #84 measured that the instrument goes
                blind there. Not attempted, and that is why this round is about W alone.
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
MASK=0.15; SEEDS=[11,29,47]; LAMS=[0.5,1.,2.,3.,4.,6.,8.,12.]
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
def plant(M,r,rng,sc=0.30):
    n,m=M.shape; F=rng.normal(size=(n,r)); L=rng.normal(size=(r,m))*sc
    p=np.clip(M.mean(0)[None,:]+F@L,0.02,0.98)
    return (rng.random((n,m))<p).astype(float)
def run(M,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    val=(rng.random(M.shape)<0.2)&obs; fit=obs&~val
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I; rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    R=T1-P; Rf=np.where(np.isnan(R),0.,R)
    base=np.mean((M[he]-gm)**2); IB=np.broadcast_to(I,M.shape)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    b0=f(IB,P); out=[]
    for lam in LAMS:
        F=np.where(fit,Rf,0.)
        for _ in range(20):
            U,S,V=svd(F,full_matrices=False); F=np.where(fit,Rf,(U*np.maximum(S-lam,0.))@V)
        U,S,V=svd(F,full_matrices=False)
        W=(U*np.maximum(S-lam,0.))@V; er=int((S>lam).sum())
        out.append(dict(lam=lam,val=np.mean((Rf[val]-W[val])**2),skill=f(IB,P,W)-b0,eff_rank=er))
    return out
rows=[]
for i,t in enumerate(IDENT):
    M=RAW[t]
    for sd in SEEDS:
        rg=np.random.default_rng(2200+sd)
        for w,Mw in {'real':M,'r2':plant(M,2,rg),'r5':plant(M,5,rg),'margin':curveball(M,rg)}.items():
            for r in run(Mw,sd): rows.append(dict(q=t,world=w,seed=sd,**r))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print("\n=== SKILL AND EFFECTIVE RANK ACROSS THE LAMBDA CURVE ===")
print(D.groupby(['world','lam'])[['skill','eff_rank']].mean().unstack('world').round(3).to_string())
print("\n=== AT THE TRAINING-TUNED LAMBDA (validation MSE minimised, never held-out) ===")
best=D.loc[D.groupby(['world','q','seed']).val.idxmin()]
B=best.groupby('world').agg(lam=('lam','mean'),eff_rank=('eff_rank','mean'),
                            rank_sd=('eff_rank','std'),skill=('skill','mean'))
print(B.round(3).to_string())
print("\n  CONDITIONAL KILL -- gates first")
g1=B.loc['r5','eff_rank']>B.loc['r2','eff_rank']
g2=B.loc['margin','eff_rank']<B.loc['r2','eff_rank']
print(f"   (a) effective rank monotone in TRUE rank : {'PASS' if g1 else 'FAIL'} "
      f"(r2 {B.loc['r2','eff_rank']:.2f}, r5 {B.loc['r5','eff_rank']:.2f})")
print(f"   (b) collapses on the structureless world : {'PASS' if g2 else 'FAIL'} "
      f"({B.loc['margin','eff_rank']:.2f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    er=B.loc['real','eff_rank']
    print(f"\n   real effective rank {er:.2f} (sd {B.loc['real','rank_sd']:.2f})   "
          f"controls: r2 -> {B.loc['r2','eff_rank']:.2f}, r5 -> {B.loc['r5','eff_rank']:.2f}, "
          f"no-structure -> {B.loc['margin','eff_rank']:.2f}")
    if er<=B.loc['r5','eff_rank']+1:
        print(f"\n   -> #83's KNEE AT 5 SURVIVES THE METHOD CHANGE. The within-block structure reads")
        print(f"      as effective rank {er:.1f} on an estimator whose floor is not an artifact.")
    else:
        print(f"\n   -> #83 WAS AN ARTIFACT of hard truncation. Effective rank {er:.1f} is well above")
        print(f"      what a true-rank-5 world produces here.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
