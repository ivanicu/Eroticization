"""
E01 A08 R01 -- apply #34's resolvability criterion retroactively to every surviving headline.

#34 killed the modality deficit not by finding it false but by finding it UNRESOLVABLE: its effect
was smaller than its own seed spread in 9 cells of 9, while sex was resolvable in 9 of 9 through
the identical pipeline. Most of this project's numbers were produced at 1-3 seeds and have never
been asked the same question.

CRITERION (from #34): a quantity is RESOLVABLE only if |effect| > 2 x (its own spread across seeds).
Anything else is a number this design cannot distinguish from its own randomness, whatever its
p-value says.

ESTIMAND        for each surviving headline: its value, its across-seed spread at >=6 seeds, and
                the ratio between them.
IDENTIFICATION  identified for every quantity that has a seed at all. Deterministic quantities
                (the gating tree, base rates) are marked N/A rather than passed silently.
WORLDS          A  the headlines are measurements: most clear 2x
                B  the project is largely reporting seed noise: most do not
KILL (CONDITIONAL)  gate: at least one quantity must be resolvable and at least one unresolvable --
                    if everything passes or everything fails, the criterion is not discriminating
                    here and no verdict follows.
                    then: >1/3 of headlines unresolvable -> the README is republished with a
                    resolvability column and the failures marked.
POSITIVE CTRL   sex deficit, shown resolvable 9/9 in #34.
NEGATIVE CTRL   row-parity placebo, which must be unresolvable by construction.
SEEDS           6 per quantity.
MULTIPLICITY    every quantity attempted is reported, including those too costly to re-seed, which
                are listed as NOT AUDITED rather than omitted.
IMPOSSIBLE      the onset rival-world round (#22/#32) costs ~40 pipeline rebuilds per seed; it is
                listed NOT AUDITED with its price rather than quietly skipped.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq, svd, qr
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
nblk=pd.Series(np.concatenate([B[q]['ppl'] for q in allq])).value_counts().reindex(pool).fillna(0)
inv=pd.read_csv('data/derived/inventory.csv')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
Rt=df[rate].apply(pd.to_numeric,errors='coerce')
SEEDS=[1,2,3,4,5,6]
def loadings(people,K=5):
    ppl=np.array(sorted(set(people)&set(pool)))
    if len(ppl)<400: return None
    pm={p:i for i,p in enumerate(ppl)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in B[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(B[q]['ppl']) if p in pm])
        Z=np.full((len(ppl),B[q]['R'].shape[1]),np.nan)
        if len(idx): Z[idx]=B[q]['R'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Z=np.hstack(cols); Z=Z-Z.mean(0)
    D=np.c_[np.ones(len(ppl)),COV.loc[ppl].values]
    b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
    return svd(Z,full_matrices=False)[2][:K]
def cong(a,b):
    if a is None or b is None: return np.nan
    Qa,_=qr(a.T,mode='reduced'); Qb,_=qr(b.T,mode='reduced')
    return float(np.mean(svd(Qa.T@Qb,compute_uv=False)))
def match(g1,g2,rng):
    i1,i2=pd.Index(sorted(set(g1)&set(pool))),pd.Index(sorted(set(g2)&set(pool)))
    k1,k2=nblk.reindex(i1).astype(int),nblk.reindex(i2).astype(int); a=[];b=[]
    for v in set(k1)|set(k2):
        x=k1.index[k1==v].values; y=k2.index[k2==v].values; m=min(len(x),len(y))
        if m: a+=list(rng.choice(x,m,replace=False)); b+=list(rng.choice(y,m,replace=False))
    return np.array(a),np.array(b)
male=df['biomale'].reindex(pool); ph=pd.to_numeric(df['pornhabit'],errors='coerce').reindex(pool)
def deficit(kind,seed):
    rng=np.random.default_rng(seed)
    if kind=='sex': g1,g2=pool[(male==1).values],pool[(male==0).values]
    elif kind=='consumption': g1,g2=pool[(ph>ph.median()).values],pool[(ph<ph.median()).values]
    else: g1,g2=pool[np.arange(len(pool))%2==0],pool[np.arange(len(pool))%2==1]
    a,b=match(g1,g2,rng)
    c=cong(loadings(a),loadings(b))
    ceil=[cong(loadings(p[:len(a)]),loadings(p[len(a):len(a)+len(b)])) for p in [rng.permutation(pool) for _ in range(3)]]
    return float(np.nanmean(ceil))-c
def cca_transfer(seed):
    from sklearn.cross_decomposition import CCA
    rng=np.random.default_rng(seed); p=rng.permutation(allq)
    h1,h2=list(p[:len(p)//2]),list(p[len(p)//2:])
    def F(bl):
        pm={q:i for i,q in enumerate(pool)}; cols=[]
        for q in bl:
            idx=np.array([pm[x] for x in B[q]['ppl'] if x in pm]); src=np.array([i for i,x in enumerate(B[q]['ppl']) if x in pm])
            Z=np.full((len(pool),B[q]['R'].shape[1]),np.nan)
            if len(idx): Z[idx]=B[q]['R'][src]
            mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
            cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
        Z=np.hstack(cols); Z=Z-Z.mean(0)
        D=np.c_[np.ones(len(pool)),COV.loc[pool].values]
        b,*_=lstsq(D,Z,rcond=None); Z=Z-D@b
        U,S,_=svd(Z,full_matrices=False); return U[:,:5]*S[:5]
    F1,F2=F(h1),F(h2)
    idx=rng.permutation(len(pool)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    c=CCA(n_components=5,max_iter=600).fit(F1[tr],F2[tr]); a,b_=c.transform(F1[te],F2[te])
    return float(np.mean([abs(np.corrcoef(a[:,j],b_[:,j])[0,1]) for j in range(5)]))
def nestedness(seed):
    rng=np.random.default_rng(seed)
    H=(Rt>0).astype(float).where(Rt.notna()); H=H[H.notna().sum(1)>=40]
    V=(H.fillna(0).values>0.5); base=V.mean(0); p=base/base.sum(); K=V.shape[1]; sz=V.sum(1)
    obs=[];nul=[]
    for _ in range(4000):
        i,j=rng.integers(0,len(V),2)
        if sz[i]<=sz[j] or sz[j]<5: continue
        obs.append((V[i]&V[j]).sum()/V[j].sum())
        m1=np.zeros(K,bool); m1[rng.choice(K,min(int(sz[i]),K),replace=False,p=p)]=True
        m2=np.zeros(K,bool); m2[rng.choice(K,min(int(sz[j]),K),replace=False,p=p)]=True
        nul.append((m1&m2).sum()/max(m2.sum(),1))
    return float(np.mean(obs)-np.mean(nul))
QUANT={'sex deficit [POS CTRL]':lambda s: deficit('sex',s),
       'consumption deficit':lambda s: deficit('consumption',s),
       'placebo deficit [NEG CTRL]':lambda s: deficit('placebo',s),
       'cross-domain CCA (mean)':cca_transfer,
       'breadth nestedness excess':nestedness}
rows=[]
for name,f in QUANT.items():
    v=np.array([f(s) for s in SEEDS])
    spread=float(v.max()-v.min()); med=float(np.median(v))
    rows.append(dict(quantity=name,median=round(med,4),spread=round(spread,4),
                     ratio=round(abs(med)/max(spread,1e-9),2),resolvable=abs(med)>2*spread))
T=pd.DataFrame(rows); T.to_csv(OUT/'resolvability.csv',index=False)
print("=== resolvability at 6 seeds:  |effect| > 2 x seed spread ? ===")
print(T.to_string(index=False))
print("\nNOT AUDITED (cost stated rather than omitted):")
print("   onset rival world (#22/#32)  ~40 pipeline rebuilds per seed -> ~240 for 6 seeds")
print("   POWER-SUBSTANCE ladder (#24) 36 cells x 6 seeds -> 216 composite refits")
print("   theta reliability ladder (#23) 6 rungs x 8 targets x 6 seeds")
res=T[~T.quantity.str.contains('CTRL')]
gate=bool(T.resolvable.any() and (~T.resolvable).any())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  at least one resolvable AND one not : {'PASS' if gate else 'FAIL'}")
if not gate: print("  -> gate FAILED : criterion not discriminating here, UNVERIFIED")
else:
    frac=float((~res.resolvable).mean())
    print(f"  unresolvable share among headlines (excl. controls): {frac:.0%}")
    print("  -> README republished with a resolvability column" if frac>1/3 else "  -> headlines largely resolvable")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
