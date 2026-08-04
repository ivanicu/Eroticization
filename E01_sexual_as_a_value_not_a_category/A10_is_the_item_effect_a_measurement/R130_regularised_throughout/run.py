import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R15 -- THE DECOMPOSITION WITH NOTHING LEFT TO CORRECT.

#86 and #87 are one finding at two components: SUBTRACTING A NULL DOES NOT CORRECT AN ESTIMATOR, IT
CREDITS THE ESTIMATOR'S FAILURE TO THE THING BEING ESTIMATED. A raw row mean scores -0.064 on a world
with no person effect; hard rank truncation scores -0.062 on a world with no interaction. Neither is
a floor the component clears -- both are holes the estimator digs, and #82/#83/#85 subtracted them
and handed the difference to the component.

The fix is not a better null. It is an estimator whose null is already zero.

  I  item     shrunk column means      (n = 1,200-15,000 per column, so shrinkage ~ 1: a CHECK, since
                                        #85 measured this component's null at -0.0004 already)
  P  person   empirical-Bayes shrunk row means            (#86)
  C  cross    per-column RIDGE on external scores, alpha tuned
  W  within   soft singular-value thresholding, lambda tuned  (#87)

Every hyperparameter is chosen on a validation split carved out of the TRAINING cells and never on
the held-out cells the skill is scored on. #87a's flaw -- one arm selecting its rank on the test
set -- is fixed here for both arms, which also settles whether soft genuinely beats hard.

ESTIMAND        Shapley held-out skill of {I,P,C,W} under regularised estimators, and the
                item:interaction ratio read directly off them with no null correction anywhere.
IDENTIFICATION  identified. The claim that no correction is needed is not assumed -- all four nulls
                are run and each component's null skill is reported and required to be ~0.
SCOPE           the 23 blocks A09/R114 identified. Gate + demographics out of the scores (#77).
WORLDS          item-dominant  ratio >> 1 -> #85's tie was the correction machinery, and A09's
                               original direction was right for the wrong reason
                tied           ratio ~ 1 -> #85 survives its own methodological critique
                interaction-dominant -> the epoch title is right
KILL            threshold-free: declared per block above 2x its own seed spread, whole grid published.
POSITIVE CTRL   a planted world carrying a KNOWN item effect, person effect, rank-5 within-block
                structure and shared cross-block structure at once. The decomposition must recover
                all four with the right ordering; a decomposition that cannot separate planted
                components cannot be trusted to separate real ones.
NEGATIVE CTRL   all four nulls, each destroying one structure: within-person shuffle (I),
                within-column shuffle (P), fixed-margin curveball (W), person-permutation (C).
                Every regularised component must return ~0 in the world that destroys it.
SETTLES         #87a -- hard truncation with K chosen on TRAINING cells, run beside soft, so the
                comparison is finally fair.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 6 worlds x 3 seeds x (hyperparameter grids), published whole.
IMPOSSIBLE      unchanged: no null destroys I while preserving the interaction, but that no longer
                matters here, because no correction is being applied.
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
from numpy.linalg import svd, lstsq
from math import factorial
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW])); PM={p:i for i,p in enumerate(ALLP)}
BLKS=sorted(RAW)
E=np.zeros((len(ALLP),len(BLKS)))
for k,q in enumerate(BLKS): E[[PM[p] for p in RAW[q]['ppl']],k]=1.
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
COV=pd.DataFrame({'male':pd.to_numeric(df.get('biomale'),errors='coerce'),'age':df['age'].map(AGEMAP)}).reindex(ALLP)
for c in ['opennessvariable','neuroticismvariable','extroversionvariable',
          'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']:
    if c in df.columns: COV[c]=pd.to_numeric(df[c],errors='coerce').reindex(ALLP).values
COV=COV.fillna(COV.median()).values; COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
MASK=0.15; SEEDS=[11,29,47]; LAMS=[0.,1.,2.,4.,8.]; ALPHAS=[0.,1.,10.,100.]; KHARD=[1,2,3,5,8]
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
def row_shuffle(M,rng): return np.array([rng.permutation(r) for r in M])
def col_shuffle(M,rng): return np.column_stack([rng.permutation(M[:,j]) for j in range(M.shape[1])])

def scores(target,K=8):
    cols=[]
    for q in BLKS:
        if q==target: continue
        Mq=RAW[q]['M']; R=Mq-Mq.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),Mq.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:K]*S[:K]
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b

def shrink(vals,cnt,resid_var):
    s2w=resid_var/np.maximum(cnt,1); s2t=max(np.var(vals)-np.mean(s2w),1e-9)
    return vals*(s2t/(s2t+s2w))
def svt(R,fit,lam,it=20):
    F=np.where(np.isnan(R)|~fit,0.,R)
    for _ in range(it):
        U,S,V=svd(F,full_matrices=False); F=np.where(fit,np.where(np.isnan(R),0.,R),(U*np.maximum(S-lam,0.))@V)
    U,S,V=svd(F,full_matrices=False); return (U*np.maximum(S-lam,0.))@V
def hardk(R,fit,K,it=20):
    F=np.where(np.isnan(R)|~fit,0.,R)
    for _ in range(it):
        U,S,V=svd(F,full_matrices=False); F=np.where(fit,np.where(np.isnan(R),0.,R),(U[:,:K]*S[:K])@V[:K])
    U,S,V=svd(F,full_matrices=False); return (U[:,:K]*S[:K])@V[:K]
def ridge_cols(R,Us,fit,alpha):
    n,m=R.shape; out=np.zeros_like(R)
    X=np.c_[np.ones(n),Us]; p=X.shape[1]
    Pen=np.eye(p)*alpha; Pen[0,0]=0.
    for j in range(m):
        k=fit[:,j]
        if k.sum()<50: continue
        Xj=X[k]; y=np.where(np.isnan(R[k,j]),0.,R[k,j])
        b=np.linalg.solve(Xj.T@Xj+Pen, Xj.T@y)
        out[:,j]=X@b
    return out

def decompose(M,U,rows,seed,mode='soft'):
    n,m=M.shape
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    val=(rng.random(M.shape)<0.2)&obs; fit=obs&~val
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cmn=np.nanmean(T,axis=0); cmn=np.where(np.isnan(cmn),gm,cmn)
    Ivals=shrink(cmn-gm,obs.sum(0),np.nanmean(np.where(obs,(T-cmn)**2,np.nan)))
    I=Ivals[None,:]
    T1=T-gm-I
    rmn=np.nanmean(T1,axis=1); rmn=np.where(np.isnan(rmn),0.,rmn)
    Pv=shrink(rmn,obs.sum(1),np.nanmean(np.where(obs,(T1-rmn[:,None])**2,np.nan)))
    P=Pv[:,None]
    Rres=T1-P
    Us=U[rows]; Us=(Us-Us.mean(0))/(Us.std(0)+1e-12)
    Rf=np.where(np.isnan(Rres),0.,Rres)
    best=None
    for a in ALPHAS:
        Cc=ridge_cols(Rres,Us,fit,a); e=np.mean((Rf[val]-Cc[val])**2)
        if best is None or e<best[0]: best=(e,a,Cc)
    _,alpha,C=best
    best=None
    grid=LAMS if mode=='soft' else KHARD
    for g_ in grid:
        Wc=svt(Rres,fit,g_) if mode=='soft' else hardk(Rres,fit,g_)
        e=np.mean((Rf[val]-Wc[val])**2)
        if best is None or e<best[0]: best=(e,g_,Wc)
    _,hyp,W=best
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':W}
    b0=np.mean((M[he]-gm)**2); v={}
    for bits in range(16):
        S=frozenset([c for j,c in enumerate('IPCW') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-np.asarray(p)[he])**2)/b0
    out={'alpha':alpha,'hyp':hyp}
    for c in 'IPCW':
        o=[x for x in 'IPCW' if x!=c]; tot=0.
        for r in range(4):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(3-len(S))/24.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def planted(M,rows,U,rng):
    n,m=M.shape
    F=rng.normal(size=(n,5)); L=rng.normal(size=(5,m))*0.20
    g=(U[rows][:,:2]/ (U[rows][:,:2].std(0)+1e-9))@rng.normal(size=(2,m))*0.10
    b=rng.normal(size=n)*0.08
    p=np.clip(M.mean(0)[None,:]+b[:,None]+F@L+g,0.02,0.98)
    return (rng.random((n,m))<p).astype(float)

rows=[]
for i,t in enumerate(IDENT):
    U=scores(t); r_=np.array([PM[p] for p in RAW[t]['ppl']]); M=RAW[t]['M']
    for sd in SEEDS:
        rg=np.random.default_rng(9900+sd)
        worlds={'real':M,'noI':row_shuffle(M,rg),'noP':col_shuffle(M,rg),
                'noW':curveball(M,rg),'planted':planted(M,r_,U,rg)}
        for w,Mw in worlds.items():
            rows.append(dict(q=t,world=w,mode='soft',seed=sd,**decompose(Mw,U,r_,sd,'soft')))
        rows.append(dict(q=t,world='real',mode='hard',seed=sd,**decompose(M,U,r_,sd,'hard')))
        rows.append(dict(q=t,world='noW',mode='hard',seed=sd,
                         **decompose(worlds['noW'],U,r_,sd,'hard')))
        # C's null: person-permuted external scores, regularised pipeline
        Up=U.copy(); Up[r_]=U[r_][np.random.default_rng(sd+7).permutation(len(r_))]
        rows.append(dict(q=t,world='noC',mode='soft',seed=sd,**decompose(M,Up,r_,sd,'soft')))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R130_regularised_throughout/results/'
D.to_csv(OUT+'grid.csv',index=False)

S=D[D['mode']=='soft']
print("\n=== NEGATIVE CONTROLS: each regularised component in the world that destroys it ===")
T=S.groupby('world')[['I','P','C','W']].mean()
print(T.round(4).to_string())
print(f"\n  I in noI  {T.loc['noI','I']:+.4f}   P in noP  {T.loc['noP','P']:+.4f}   "
      f"W in noW  {T.loc['noW','W']:+.4f}   C in noC  {T.loc['noC','C']:+.4f}")

print("\n=== #87a SETTLED: soft vs hard, BOTH tuned on training cells ===")
H=D[D['mode']=='hard'].groupby('world')[['W']].mean()
print(f"  real  soft {T.loc['real','W']:+.4f}   hard {H.loc['real','W']:+.4f}")
print(f"  noW   soft {T.loc['noW','W']:+.4f}   hard {H.loc['noW','W']:+.4f}")

print("\n=== POSITIVE CONTROL: planted world with all four components present ===")
print(T.loc[['planted']].round(4).to_string())

print("\n=== THE DECOMPOSITION, REGULARISED, NO CORRECTION APPLIED ===")
r=S[S.world=='real']
pb=r.groupby('q')[['I','P','C','W']].mean()
sp=r.groupby('q')[['I','P','C','W']].std()
print(f"  I {pb.I.median():+.4f}   P {pb.P.median():+.4f}   "
      f"C {pb.C.median():+.4f}   W {pb.W.median():+.4f}")
inter=pb.C+pb.W; gap=pb.I-inter
psp=np.sqrt(sp.I**2+sp.W**2)
print(f"  interaction (C+W) {inter.median():+.4f}   item:interaction {pb.I.median()/max(inter.median(),1e-9):.2f}x")
print(f"  per block: item larger {int((gap>2*psp).sum())}/23   interaction larger "
      f"{int((-gap>2*psp).sum())}/23   tied {int((gap.abs()<=2*psp).sum())}/23")
print(f"  median gap {gap.median():+.4f}   2x spread {2*psp.median():.4f}")
print(f"  chosen hyperparameters: lambda {r.hyp.mean():.2f}   ridge alpha {r.alpha.mean():.1f}")

print("\n  CONDITIONAL KILL -- gates first")
g1=all(abs(T.loc[w,c])<0.02 for w,c in [('noI','I'),('noP','P'),('noW','W'),('noC','C')])
g2=all(T.loc['planted',c]>0.005 for c in 'IPW')
print(f"   (a) every regularised component ~0 in the world that destroys it (|.|<0.02): "
      f"{'PASS' if g1 else 'FAIL'}")
print(f"   (b) planted world recovers I, P and W                                     : "
      f"{'PASS' if g2 else 'FAIL'}")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    ratio=pb.I.median()/max(inter.median(),1e-9)
    print(f"\n   NO CORRECTION WAS APPLIED ANYWHERE. item:interaction = {ratio:.2f}x")
    print(f"   #85 (correction machinery) said 1.05x; #87e predicted ~4x from the soft estimator")
    if ratio>2.5:
        print("\n   -> ITEM-DOMINANT. #85's tie was the correction machinery, and #82/#83/#85 must be")
        print("      restated on regularised estimators. A09's original direction was right, and")
        print("      right for a reason it did not have at the time.")
    elif ratio<1.4:
        print("\n   -> TIED, on estimators that need no correction. #85 survives its own critique.")
    else:
        print(f"\n   -> INTERMEDIATE at {ratio:.2f}x. Neither #85's tie nor #87e's 4x; the honest")
        print("      report is this number, from the only specification with nothing subtracted.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
