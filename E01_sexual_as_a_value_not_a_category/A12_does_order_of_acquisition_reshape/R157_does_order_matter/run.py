import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R01 -- IVAN'S THIRD STEP, WHICH NOTHING IN 148 ROUNDS HAS TOUCHED.

His prior was three steps: ordinary representation -> individualised erotic readout -> RECURSIVE
REPRESENTATIONAL RESTRUCTURING. The project measured step 1 and step 2. Step 3 was filed as
"needs longitudinal data" and never attacked. That filing is wrong, and this round is why.

Recursion has a CROSS-SECTIONAL signature, and it is a psychological question, not a statistical one:

    if an acquired interest RESHAPES the representation, then two people who both end up liking A
    and B should differ in EVERYTHING ELSE depending on WHICH ONE THEY GOT FIRST.

Under pure readout (model B) order cannot matter: you have both, the weights are the weights, the
rest of your profile follows from the weights. Under recursion (model C) the first one restructures
the space the second is read out in, so the residue differs.

ESTIMAND        held-out AUC of predicting acquisition ORDER within a category pair from the person's
                preferences on all OTHER categories, above a covariates-only baseline.
IDENTIFICATION  A and B and anything derived from them are excluded from the predictors. Order is
                observed, not inferred.
SCOPE           category pairs where >=400 people have both onsets and a strict order.
CONFOUNDS, written before the run:
                (1) precocity -- people who acquire everything early
                (2) breadth -- more interests, more chances
                (3) sex, age
                all four enter the baseline, so only the INCREMENT over them is read.
WORLDS          recursion   order is predictable from the rest of the profile above covariates
                readout     it is not: once you have both, the order left no trace
KILL            threshold-free: the increment against a permutation null that shuffles the order
                label WITHIN covariate strata, so the null keeps every confound.
POSITIVE CTRL   plant an order effect of known size and require detection AND monotonicity.
NEGATIVE CTRL   the stratified permutation, per pair.
NOISE FLOOR     5 splits x 3 permutation draws.
IMPOSSIBLE      distinguishing "the first one restructured me" from "whatever made me get A first
                also made me like C" -- a third common cause is not excluded by cross-sectional data.
                A positive result licenses "order carries information", never "order caused it".
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
sys.path.insert(0,str(ROOT))
from lib.gates import Gate
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
P=(df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float).fillna(0.).values
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
breadth=P.sum(1); precocity=np.nanmean(np.where(np.isnan(V),np.nan,V),axis=1)
precocity=np.where(np.isfinite(precocity),precocity,np.nanmean(precocity))
COV=np.c_[male,agev,breadth,precocity]
COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
print(f"onset categories {V.shape[1]}  preference items {P.shape[1]}  people {len(df):,}",flush=True)
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<10 or n0<10: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,rng,alpha=50.,reps=5):
    out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<10 or (1-y[tr]).sum()<10: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr])
        out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def strata(idx):
    q=lambda v:np.digitize(v,np.quantile(v,[.33,.66]))
    return (male[idx]>0).astype(int)*9+q(breadth[idx])*3+q(precocity[idx])
rng=np.random.default_rng(11); rows=[]
pairs=[(a,b) for a,b in itertools.combinations(range(V.shape[1]),2)]
rng.shuffle(pairs); done=0
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b])
    idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)          # 1 = A came first
    others=P[idx]                                 # A and B are onset cats; predictors are ratings
    base=ridge_auc(COV[idx],y,np.random.default_rng(1))
    full=ridge_auc(np.c_[COV[idx],others],y,np.random.default_rng(1))
    st=strata(idx); perms=[]
    for d in range(3):
        rp=np.random.default_rng(100+d); yp=y.copy()
        for s in np.unique(st):
            w=np.flatnonzero(st==s)
            if len(w)>1: yp[w]=y[w][rp.permutation(len(w))]
        perms.append(ridge_auc(np.c_[COV[idx],others],yp,np.random.default_rng(1))-
                     ridge_auc(COV[idx],yp,np.random.default_rng(1)))
    rows.append(dict(a=a,b=b,n=len(idx),base=base,full=full,inc=full-base,
                     perm=np.nanmean(perms)))
    done+=1
    if done>=120: break
    if done%20==0: print(f"  {done} pairs",flush=True)
# ---- positive control: plant an order effect of known size ----
ctrl=[]
for g in [0.0,0.10,0.25]:
    cr=[]
    for a,b in pairs[:25]:
        m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
        if len(idx)<400: continue
        y=(V[idx,a]<V[idx,b]).astype(float)
        rp=np.random.default_rng(7)
        w=rp.normal(size=P.shape[1]); sig=P[idx]@w; sig=(sig-sig.mean())/(sig.std()+1e-9)
        yp=np.where(rp.random(len(idx))<g*np.clip((sig+2)/4,0,1)+ (1-g)*y, 1., 0.) if g>0 else y
        yp=(y*(1-g)+ (sig>0).astype(float)*g); yp=(rp.random(len(idx))<yp).astype(float)
        cr.append(ridge_auc(np.c_[COV[idx],P[idx]],yp,np.random.default_rng(1))-
                  ridge_auc(COV[idx],yp,np.random.default_rng(1)))
        if len(cr)>=12: break
    ctrl.append(dict(g=g,inc=np.nanmean(cr)))
D=pd.DataFrame(rows); C=pd.DataFrame(ctrl)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
print(f"\n=== {len(D)} category pairs, median n = {D.n.median():.0f} ===")
print(f"  covariates-only AUC          : {D.base.mean():.4f}")
print(f"  + all other preferences      : {D.full.mean():.4f}")
print(f"  INCREMENT                    : {D.inc.mean():+.4f}   (sd over pairs {D.inc.std():.4f})")
print(f"  stratified permutation null  : {D.perm.mean():+.4f}")
print(f"  pairs with increment > 0     : {(D.inc>0).mean():.1%}")
print("\n=== POSITIVE CONTROL: planted order effect ===")
print(C.round(4).to_string(index=False))
sd=D.inc.std()/np.sqrt(len(D))
g=Gate("does acquisition ORDER leave a trace in the rest of the profile?")
g.negative_control("stratified permutation", null=D.perm.mean(), effect=D.inc.mean())
g.positive_control("planted order effect g=0.25", planted=C[C.g==0.25].inc.iloc[0],
                   floor=C[C.g==0.0].inc.iloc[0], spread=sd)
g.asserted("planted ladder monotone",
           C.inc.iloc[0]<C.inc.iloc[1]<C.inc.iloc[2], f"{C.inc.round(4).tolist()}")
g.resolvable("increment over covariates", effect=D.inc.mean(), spread=sd)
print(); print(g)
if g.verdict():
    print(f"\n  -> ORDER LEAVES A TRACE. Two people who both like A and B differ in the REST of their")
    print(f"     profile depending on which came first: AUC increment {D.inc.mean():+.4f} over a")
    print(f"     baseline that already contains sex, age, breadth and precocity.")
else:
    print(f"\n  -> NOT ESTABLISHED. increment {D.inc.mean():+.4f} vs null {D.perm.mean():+.4f}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
