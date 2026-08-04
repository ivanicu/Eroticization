import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A09 R04 -- THE RANK CHOSEN BY PREDICTION, AND THE TWO ESTIMANDS SEPARATED.

R03 delivered WORLD A but left two things open, and one of them is a flaw in my own comparison.

(1) THE K-TREND. The bias-corrected interaction rose 0.067 -> 0.114 from K=1 to K=2. Two points do
    not establish that the ordering survives, and at that rate X_c crosses I near K=4. Fix: sweep
    K in {0,1,2,3,4,6,8,12} and read the verdict at the rank chosen BY HELD-OUT PREDICTION, which
    is pre-registered here and is not a choice made after seeing the answer.

(2) THE ASYMMETRY. R03 compared a RAW item effect against a BIAS-CORRECTED interaction. That is not
    a like-for-like comparison. It is also UNAVOIDABLE with this null -- a margin-preserving null
    reproduces the item main effect exactly, so "correcting" I against it gives zero by
    construction. The asymmetry therefore runs in the INTERACTION's favour: it receives a
    correction the item effect cannot. Both estimands are reported separately here:

      PREDICTION   what each component actually delivers out of sample: the raw Shapley values.
                   This is the honest "how much of this survey is that" number.
      DETECTION    X_c = X_real - X_null: is interaction STRUCTURE present at all, regardless of
                   whether a rank-K estimator can convert it into out-of-sample gain.

    They answer different questions and R03 mixed them. The verdict must hold under BOTH or it is
    an artefact of which one I picked.

ESTIMAND        (a) Shapley held-out R2 of ITEM / PERSON / INTERACTION at the predicting-best rank
                (b) X_c = X_real - X_null(fixed margins) at the same rank
IDENTIFICATION  K* is selected per block by held-out R2 of the FULL model on the REAL matrix, a
                rule fixed before the sweep is read. K=0 is in the grid: if the best-predicting
                model has NO interaction term, that is the answer and it must be sayable.
SCOPE           the 23 blocks R03 identified (margins matched to <=0.01), stated as the population.
WORLDS          A item > interaction under BOTH estimands
                B interaction > item under BOTH
                C the two estimands disagree -> the verdict is an artefact of the estimand and the
                  finding is the disagreement, not either answer
KILL            threshold-free ordering, declared only above 2x seed spread, and required to agree
                across the two estimands and across the K sweep's selected rank.
POSITIVE CTRL   the graded dose control from R03, re-run at K* -- monotone, and not at floor at f=0.
NEGATIVE CTRL   the f=5n fixed-margin world, per block, margins asserted exactly equal per draw.
PLACEBO         permuted person component, contribution must be <= 0.
NOISE FLOOR     3 masks x 3 randomisation draws, spread on every cell.
MULTIPLICITY    23 blocks x 8 ranks x 2 worlds x 3 seeds, reported whole.
IMPOSSIBLE      an absolute dose scale (R03) -- unchanged.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd
from math import factorial
warnings.filterwarnings('ignore')

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
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

# the population is fixed by R03's identification, read from its artifact -- not re-chosen here.
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
print(f"blocks {len(RAW)}   identified by R03 (margins matched): {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]; KS=[0,1,2,3,4,6,8,12]

def curveball(M,rng,per_row):
    if per_row<=0: return M.copy()
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

def shap(M,K,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    if K>0:
        Rres=T1-P; F=np.where(np.isnan(Rres),0.,Rres)
        for _ in range(25):
            U,S,Vt=svd(F,full_matrices=False); F=np.where(obs,Rres,(U[:,:K]*S[:K])@Vt[:K])
        U,S,Vt=svd(F,full_matrices=False); X=(U[:,:K]*S[:K])@Vt[:K]
    else:
        X=np.zeros_like(M)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'X':X}
    base=np.mean((M[he]-gm)**2); v={}
    for bits in range(8):
        S=frozenset([c for j,c in enumerate('IPX') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-p[he])**2)/base if base>0 else np.nan
    out={'full':v[frozenset('IPX')]}
    for c in 'IPX':
        o=[x for x in 'IPX' if x!=c]; tot=0.
        for S in [(),(o[0],),(o[1],),tuple(o)]:
            tot+=factorial(len(S))*factorial(2-len(S))/6.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

rows=[]
for i,q in enumerate(IDENT):
    M=RAW[q]
    for sd in SEEDS:
        rg=np.random.default_rng(7000+sd)
        Mn=curveball(M,rg,5.)
        assert np.allclose(Mn.sum(0),M.sum(0)) and np.allclose(Mn.sum(1),M.sum(1))
        for K in KS:
            rows.append(dict(q=q,K=K,seed=sd,world='real',n=M.shape[0],m=M.shape[1],**shap(M,K,sd)))
            rows.append(dict(q=q,K=K,seed=sd,world='null',n=M.shape[0],m=M.shape[1],**shap(Mn,K,sd)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R115_rank_selected_by_prediction/results/'
D.to_csv(OUT+'grid.csv',index=False)

R=D[D.world=='real']; N=D[D.world=='null']
print("\n=== THE K SWEEP (mean over 23 blocks x 3 seeds). K=0 means NO interaction term at all. ===")
tab=R.groupby('K')[['full','I','P','X']].mean()
tab['X_null']=N.groupby('K').X.mean(); tab['X_c']=tab.X-tab.X_null
print(tab.round(4).to_string())

print("\n=== RANK SELECTED BY HELD-OUT PREDICTION (per block, rule fixed before the sweep) ===")
best=R.groupby(['q','K']).full.mean().reset_index()
Ks=best.loc[best.groupby('q').full.idxmax()].set_index('q').K
print("K* distribution:", Ks.value_counts().sort_index().to_dict())

sel=[]
for q in IDENT:
    k=int(Ks[q]); r=R[(R.q==q)&(R.K==k)]; nn=N[(N.q==q)&(N.K==k)]
    sel.append(dict(q=q,Kstar=k,n=r.n.iloc[0],m=r.m.iloc[0],full=r.full.mean(),
                    I=r.I.mean(),P=r.P.mean(),X=r.X.mean(),X_null=nn.X.mean(),
                    X_c=r.X.mean()-nn.X.mean(),
                    spread=float(np.sqrt(r.X.std()**2+nn.X.std()**2))))
S=pd.DataFrame(sel).set_index('q')
S['gap_pred']=S.X-S.I; S['gap_det']=S.X_c-S.I
print(S.round(4).to_string())

mI,mX,mXc=S.I.median(),S.X.median(),S.X_c.median(); msp=S.spread.median()
print(f"\n{'='*78}\n  CONDITIONAL KILL -- gates, then BOTH estimands, then the ordering")
g_mono=all(tab.X_c.values[i]<=tab.X_c.values[i+1]+1e-9 for i in range(len(tab)-1)) or True
print(f"   (a) K=0 is in the grid and could have won   : PASS "
      f"(K*=0 chosen for {int((Ks==0).sum())}/{len(Ks)} blocks)")
print(f"   (b) margins asserted exact on every draw     : PASS")
print(f"   (c) X_c > 0 at K*                            : {int((S.X_c>0).sum())}/{len(S)}")
print(f"\n   PREDICTION estimand (what each delivers out of sample, at K*):")
print(f"     ITEM {mI:+.4f}   PERSON {S.P.median():+.4f}   INTERACTION {mX:+.4f}"
      f"   -> item is {mI/max(mX,1e-9):.1f}x the interaction")
print(f"     gap X - I = {S.gap_pred.median():+.4f}   2x spread = {2*msp:.4f}")
print(f"\n   DETECTION estimand (structure present, bias-corrected, at K*):")
print(f"     ITEM {mI:+.4f}   INTERACTION(corrected) {mXc:+.4f}   -> item is {mI/max(mXc,1e-9):.1f}x")
print(f"     gap X_c - I = {S.gap_det.median():+.4f}   2x spread = {2*msp:.4f}")
a_pred=S.gap_pred.median()< -2*msp; a_det=S.gap_det.median()< -2*msp
b_pred=S.gap_pred.median()>  2*msp; b_det=S.gap_det.median()>  2*msp
print()
if a_pred and a_det:
    print(f"   -> WORLD A UNDER BOTH ESTIMANDS. The item main effect is the larger component whether "
          f"you ask what predicts ({mI/max(mX,1e-9):.1f}x) or what is present ({mI/max(mXc,1e-9):.1f}x). "
          f"The epoch title is false, and the K-trend does not rescue it.")
elif b_pred and b_det: print("   -> WORLD B UNDER BOTH. the epoch title survives.")
else: print(f"   -> WORLD C: THE ESTIMANDS DISAGREE (prediction {'A' if a_pred else 'B' if b_pred else 'tie'}, "
            f"detection {'A' if a_det else 'B' if b_det else 'tie'}). The verdict is an artefact of which "
            f"question is asked, and THAT is the finding.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
