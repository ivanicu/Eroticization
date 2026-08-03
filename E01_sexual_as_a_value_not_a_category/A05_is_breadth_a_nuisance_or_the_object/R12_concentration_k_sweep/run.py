"""
E01 A05 R12 -- is "breadth has no shape" a property of breadth, or of K=6?

A05's decision -- model the scalar gain as the object rather than control it away -- rests on
breadth being QUANTITY WITHOUT SHAPE: a person's endorsed set is only 0.88% more concentrated in
coordinate space than a size-matched base-rate set (out-of-sample loadings, #A05 R03).

That was computed at K=6 coordinate loadings and never swept. #18/#19 established this project's
CCA machinery is monotone in K everywhere, and #34 established that small effects here are often
smaller than their own seed noise. A 0.88% effect has both problems and has had neither check.

ESTIMAND        the participation-ratio gap between a person's own endorsed set and a size-matched
                base-rate set, as a function of K, with its seed spread.
IDENTIFICATION  identified at each K; loadings fitted out of sample so the gap is not guaranteed.
WORLDS          A  no shape: the gap stays near zero at every K -- breadth is a scalar gain
                B  a K artifact: |gap| grows with K, and "no shape" was a statement about K=6
KILL (CONDITIONAL) gate: the positive control must fire (a synthetically CONCENTRATED set must show
                   a large negative gap at every K) AND the base-rate null must be unbiased.
                   then: |gap| < 3% at every K and resolvable-or-null -> NO SHAPE, A05 holds
                         |gap| > 5% at any K -> SHAPE, and A05's decision needs restating
POSITIVE CTRL   synthetic people whose sets are drawn from a single coordinate direction -- they
                MUST show a large concentration gap, else the measure cannot see shape at all.
NEGATIVE CTRL   synthetic people drawn exactly from base rates -- gap must be ~0.
NOISE FLOOR     across-seed spread at each K.
MULTIPLICITY    5 K values x 4 seeds x 3 populations (real, concentrated, base-rate), all reported.
SEEDS           4.
IMPOSSIBLE      a ground-truth coordinate space -- the loadings are estimated, which is why the
                positive and negative controls are synthetic populations rather than known answers.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import svd
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce')
H=(R>0).astype(float).where(R.notna()); H=H[H.notna().sum(1)>=40]
V=(H.fillna(H.mean()).values>0.5)
print(f"people {len(V):,}  categories {V.shape[1]}")
def run(Vmat,K,seed):
    rng=np.random.default_rng(seed)
    n=len(Vmat); idx=rng.permutation(n); f,e=idx[:n//2],idx[n//2:]
    base=Vmat[f].mean(0); Z=Vmat[f].astype(float)-base
    L=svd(Z,full_matrices=False)[2][:K].T
    L=L/(np.linalg.norm(L,axis=1,keepdims=True)+1e-9)
    p=base/max(base.sum(),1e-9); Kc=Vmat.shape[1]
    def PR(mask):
        i=np.flatnonzero(mask)
        if len(i)<3: return np.nan
        Y=L[i]-L[i].mean(0); ev=np.linalg.eigvalsh(Y.T@Y)
        return (ev.sum()**2)/max((ev**2).sum(),1e-12)
    obs=[];nul=[]
    for r_ in Vmat[e]:
        m=r_>0.5; k_=int(m.sum())
        if k_<3: continue
        o=PR(m)
        ns=[]
        for _ in range(6):
            j=rng.choice(Kc,size=min(k_,Kc),replace=False,p=p)
            mm=np.zeros(Kc,bool); mm[j]=True; ns.append(PR(mm))
        if np.isfinite(o): obs.append(o); nul.append(np.nanmean(ns))
    obs=np.array(obs); nul=np.array(nul); ok=np.isfinite(obs)&np.isfinite(nul)
    return float(np.mean(obs[ok]-nul[ok])), float(np.mean(nul[ok]))
# synthetic populations for the two controls
rngp=np.random.default_rng(3); Kc=V.shape[1]; base=V.mean(0); p=base/base.sum()
sizes=V.sum(1)
CONC=np.zeros_like(V)
Ldir=rngp.normal(size=(Kc,2))
for i,k_ in enumerate(sizes):
    w=np.abs(Ldir@rngp.normal(size=2)); w=w/w.sum()
    j=rngp.choice(Kc,size=min(int(k_),Kc),replace=False,p=w)
    CONC[i,j]=True
BASE=np.zeros_like(V)
for i,k_ in enumerate(sizes):
    j=rngp.choice(Kc,size=min(int(k_),Kc),replace=False,p=p); BASE[i,j]=True
rows=[]
for K,seed in itertools.product([2,4,6,10,16],[1,2,3,4]):
    for name,M in [('real',V),('CONCENTRATED [pos ctrl]',CONC),('base-rate [neg ctrl]',BASE)]:
        g,nl=run(M,K,seed)
        rows.append(dict(pop=name,K=K,seed=seed,gap=g,null_PR=nl,pct=100*g/max(nl,1e-9)))
G=pd.DataFrame(rows); G.to_csv(OUT/'concentration_k.csv',index=False)
P=G.groupby(['pop','K']).pct.agg(['median','min','max'])
P['spread']=P['max']-P['min']
print("\n=== concentration gap as % of the size-matched null, by K ===")
print(P.round(2).to_string())
real=P.xs('real',level='pop'); pos=P.xs('CONCENTRATED [pos ctrl]',level='pop'); neg=P.xs('base-rate [neg ctrl]',level='pop')
gate_pos=bool((pos['median']<-5).all()); gate_neg=bool((neg['median'].abs()<2).all())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  concentrated population shows a big negative gap at every K : {'PASS' if gate_pos else 'FAIL'} (max {pos['median'].max():.2f}%)")
print(f"  base-rate population near zero at every K                   : {'PASS' if gate_neg else 'FAIL'} (max |{neg['median'].abs().max():.2f}|%)")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
else:
    mx=real['median'].abs().max()
    print(f"  real population: {real['median'].min():.2f}% .. {real['median'].max():.2f}%  (published -0.88% at K=6)")
    if mx<3: print(f"  -> NO SHAPE at every K (max |{mx:.2f}|%) : A05's decision holds")
    elif mx>5: print(f"  -> SHAPE appears at some K (max |{mx:.2f}|%) : A05 needs restating")
    else: print(f"  -> UNVERIFIED : max |{mx:.2f}|%")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
