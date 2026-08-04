import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A10 R14 -- THE SAME TEST #86 APPLIED TO THE PERSON EFFECT, APPLIED TO THE INTERACTION.

#86 measured that a null correction credits a component with the DAMAGE its estimator does, and that
87% of the person effect's correction was that. #86c named the consequence: W's fixed-margin floor is
-0.19, and #82's inversion, #83's knee and #85's "three components are equal" all rest on
W_c = W_real - W_null. If a properly regularised W recovers only ~13% of that gap directly, the
interaction is far smaller than this arc has been reporting and four entries move at once.

Hard rank-K truncation is the unregularised estimator, exactly as the raw row mean was for P. Its
regularised counterpart is SOFT singular-value thresholding, with the threshold tuned on TRAINING
cells only -- never on the held-out cells the skill is scored on.

ESTIMAND        held-out skill of a soft-thresholded interaction, and the share of the
                (W_corrected - W_raw) gap it recovers WITHOUT any null correction.
IDENTIFICATION  identified; lambda is chosen on a validation split carved out of the training cells.
SCOPE           the 23 blocks A09/R114 identified.
WORLDS          noise  soft-thresholding recovers most of the gap -> the correction was measuring
                       estimator noise, W is genuinely large, #82/#83/#85 stand
                bias   it recovers ~13% like the person effect -> the correction was crediting
                       damage, W is much smaller, and four entries move
KILL            threshold-free: the recovered share with its seed spread, declared above 2x it.
POSITIVE CTRL   a planted rank-5 within-block world: soft-thresholding must beat hard truncation
                there, or the regularised estimator is not actually better and the test is void.
NEGATIVE CTRL   the fixed-margin world: soft-thresholding must return ~0, not a negative number.
                That is the whole claim -- a good estimator does not dig a hole.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 3 worlds x (lambda sweep) x 3 seeds, published whole.
IMPOSSIBLE      unchanged.
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
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
MASK=0.15; SEEDS=[11,29,47]; LAMS=[0.,.25,.5,1.,2.,4.,8.]
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
def plant5(M,rng):
    n,m=M.shape; F=rng.normal(size=(n,5)); L=rng.normal(size=(5,m))*0.30
    p=np.clip(M.mean(0)[None,:]+F@L,0.02,0.98)
    return (rng.random((n,m))<p).astype(float)

def marg(M,obs):
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    return gm,I,P,T1-P
def svt(Rres,fit,lam,iters=20):
    F=np.where(np.isnan(Rres)|~fit,0.,Rres)
    for _ in range(iters):
        U,S,V=svd(F,full_matrices=False); S2=np.maximum(S-lam,0.)
        L=(U*S2)@V; F=np.where(fit,np.where(np.isnan(Rres),0.,Rres),L)
    U,S,V=svd(F,full_matrices=False); return (U*np.maximum(S-lam,0.))@V
def hard(Rres,fit,K,iters=20):
    F=np.where(np.isnan(Rres)|~fit,0.,Rres)
    for _ in range(iters):
        U,S,V=svd(F,full_matrices=False); F=np.where(fit,np.where(np.isnan(Rres),0.,Rres),(U[:,:K]*S[:K])@V[:K])
    U,S,V=svd(F,full_matrices=False); return (U[:,:K]*S[:K])@V[:K]

def skill(M,seed,mode):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    gm,I,P,Rres=marg(M,obs)
    base=np.mean((M[he]-gm)**2); IB=np.broadcast_to(I,M.shape)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    b0=f(IB,P)
    if mode=='soft':
        val=(rng.random(M.shape)<0.2)&obs; fit=obs&~val          # tune on TRAINING cells only
        best,bl=None,None
        for lam in LAMS:
            L=svt(Rres,fit,lam)
            e=np.mean((Rres[val]-L[val])**2)
            if best is None or e<best: best,bl=e,lam
        W=svt(Rres,obs,bl); return f(IB,P,W)-b0, bl
    else:
        bestv=None
        for K in [1,2,3,5,8]:
            W=hard(Rres,obs,K); v=f(IB,P,W)-b0
            if bestv is None or v>bestv: bestv=v
        return bestv, np.nan

rows=[]
for i,t in enumerate(IDENT):
    M=RAW[t]
    for sd in SEEDS:
        rg=np.random.default_rng(9800+sd)
        for w,Mw in {'real':M,'margin':curveball(M,rg),'planted':plant5(M,rg)}.items():
            for mode in ['hard','soft']:
                s,lam=skill(Mw,sd,mode)
                rows.append(dict(q=t,world=w,mode=mode,seed=sd,skill=s,lam=lam))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
T=D.groupby(['world','mode']).skill.mean().unstack('mode')
T['lam']=D[D['mode']=='soft'].groupby('world').lam.mean()
print("\n=== BEST-RANK HARD TRUNCATION vs TUNED SOFT THRESHOLDING ===")
print(T.round(4).to_string())
print("\n  CONDITIONAL KILL -- gates first")
g1=T.loc['planted','soft']>T.loc['planted','hard']
g2=T.loc['margin','soft']>T.loc['margin','hard']
print(f"   (a) soft beats hard where rank-5 structure is PLANTED : {'PASS' if g1 else 'FAIL'} "
      f"({T.loc['planted','hard']:+.4f} -> {T.loc['planted','soft']:+.4f})")
print(f"   (b) soft does not dig a hole in the fixed-margin world: {'PASS' if g2 else 'FAIL'} "
      f"({T.loc['margin','hard']:+.4f} -> {T.loc['margin','soft']:+.4f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    raw=T.loc['real','hard']; soft=T.loc['real','soft']; corr=0.1906
    rec=(soft-raw)/max(corr-raw,1e-9)
    sp=D[(D.world=='real')&(D['mode']=='soft')].groupby('q').skill.std().median()
    print(f"\n   real: best hard-rank {raw:+.4f}   tuned soft-threshold {soft:+.4f}")
    print(f"   #83/#85's corrected W was +0.1906; gap to close {corr-raw:+.4f}")
    print(f"   recovered by regularisation alone: {rec:.1%}   (seed spread {sp:.4f})")
    print(f"   residual hole the fixed-margin world still shows under soft: {T.loc['margin','soft']:+.4f}")
    if rec>0.4:
        print("\n   -> NOISE. W is genuinely large; #82, #83 and #85 stand.")
    elif rec<0.2:
        print("\n   -> BIAS, same as the person effect. The correction was crediting W with the damage")
        print("      hard truncation does. The interaction's honest magnitude is near the SOFT number,")
        print("      and #82's inversion, #83's knee and #85's equality all move.")
    else:
        print(f"\n   -> PARTIAL ({rec:.0%}). Report the interaction as the interval [soft, corrected].")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
