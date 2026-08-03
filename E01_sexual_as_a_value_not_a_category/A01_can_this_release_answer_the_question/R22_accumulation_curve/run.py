"""
E01 A01 R22 -- how does the cross-block signal accumulate with the number of source domains?

#59 established that pairwise block->block prediction is ~0 (#49) while factors from 31 other blocks
give +0.0409 -- the largest structural increment in the project. The signal ACCUMULATES. The curve
of that accumulation has never been measured, and it is the one number a phase-1 collection would
need to size itself: how many domains before the shared structure is recoverable, and does it
saturate.

ESTIMAND        held-out R2 increment over a marginals base from cross-block factors fitted on n
                source blocks, as a function of n.
IDENTIFICATION  identified; source blocks are sampled at random and never include the target.
WORLDS          A  saturating: the curve flattens, and a phase-1 design can be sized
                B  linear in log n: more domains always help, and no finite design suffices
                C  no accumulation: flat from n=1, contradicting #59 and #49 both
KILL (CONDITIONAL) gate, all measured (#41): (a) n=31 must reproduce #59's +0.0409 within its seed
                   spread -- the established endpoint; (b) the permuted-label control must be ~0 at
                   EVERY n, since a control that only holds at one n is not a control.
                   then: increment at n=8 >= 80% of n=31 -> SATURATES BY 8
                         increment still rising >20% from n=16 to n=31 -> NOT SATURATED
                         otherwise -> report the curve and where it flattens
POSITIVE CTRL   the n=31 endpoint from #59.
NEGATIVE CTRL   person labels permuted inside the cross-block factor, at every n.
SEEDS           4, with source blocks resampled per seed.
MULTIPLICITY    6 values of n x 2 arms x 4 seeds, all reported.
IMPOSSIBLE      n beyond 31 -- this release has 32 usable blocks, so the curve cannot be shown to
                saturate beyond its own right edge. Stated rather than extrapolated.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import svd, lstsq
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
M={}
for q in allq:
    s=lg[lg.qi==q]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=20].index)]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    X=np.zeros((len(ppl),len(opt))); X[s.person.map(pi).values,s.option.map(oi).values]=1
    if X.shape[0]>=800 and X.shape[1]>=8: M[q]=dict(X=X,ppl=ppl)
pool=np.unique(np.concatenate([v['ppl'] for v in M.values()]))
pmap={p:i for i,p in enumerate(pool)}
print(f"blocks {len(M)}   pool {len(pool):,}   max sources per target {len(M)-1}")
def factors(srcs,K=8):
    cols=[]
    for q in srcs:
        v=M[q]; X=v['X']; idx=np.array([pmap[p] for p in v['ppl']])
        Z=np.full((len(pool),X.shape[1]),np.nan); Z[idx]=X
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Zc=np.hstack(cols); Zc=Zc-Zc.mean(0)
    kk=min(K,Zc.shape[1]-1)
    U,S_,_=svd(Zc,full_matrices=False); return U[:,:kk]*S_[:kk]
def run(nsrc,seed,permute=False):
    rng=np.random.default_rng(seed*97+nsrc)
    D={'im':[],'pm':[],'cross':[]}; Y=[]
    for q,v in M.items():
        others=[x for x in M if x!=q]
        srcs=list(rng.choice(others,min(nsrc,len(others)),replace=False))
        F=factors(srcs)
        Fi=F[np.array([pmap[p] for p in v['ppl']])]
        X=v['X']; n,k=X.shape
        if permute: Fi=Fi[rng.permutation(n)]
        mask=rng.random((n,k))<0.30
        Xtr=np.where(mask,np.nan,X)
        im=np.nanmean(Xtr,axis=0); pm=np.nanmean(Xtr,axis=1)
        im=np.where(np.isfinite(im),im,np.nanmean(X)); pm=np.where(np.isfinite(pm),pm,np.nanmean(X))
        filled=np.where(np.isnan(Xtr),im[None,:],Xtr)
        A=np.c_[np.ones(n),Fi]; b,*_=lstsq(A,filled,rcond=None); crec=A@b
        ii,jj=np.where(mask); Y.append(X[ii,jj])
        D['im'].append(im[jj]); D['pm'].append(pm[ii]); D['cross'].append(crec[ii,jj])
    y=np.concatenate(Y); D={k:np.concatenate(v) for k,v in D.items()}
    h=rng.random(len(y))<0.5
    def sc(names):
        Xm=np.c_[np.ones(len(y)),np.column_stack([D[n] for n in names])]
        b,*_=lstsq(Xm[h],y[h],rcond=None); p=Xm[~h]@b
        return 1-((y[~h]-p)**2).sum()/((y[~h]-y[~h].mean())**2).sum()
    return sc(['im','pm','cross'])-sc(['im','pm'])
rows=[]
for nsrc in [1,2,4,8,16,31]:
    for seed in (1,2,3,4):
        rows.append(dict(n_sources=nsrc,seed=seed,increment=run(nsrc,seed),
                         null=run(nsrc,seed,permute=True)))
G=pd.DataFrame(rows); G.to_csv(OUT/'accumulation.csv',index=False)
S=G.groupby('n_sources')[['increment','null']].agg(['median',lambda s:s.max()-s.min()])
S.columns=['increment','inc_spread','null','null_spread']
print("\n=== cross-block increment by number of source domains ===")
print(S.round(4).to_string())
end=float(S.loc[31,'increment']); e8=float(S.loc[8,'increment']); e16=float(S.loc[16,'increment'])
ga=abs(end-0.0409)<0.010; gb=bool((S['null'].abs()<0.004).all())
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) n=31 reproduces #59's +0.0409 : {'PASS' if ga else 'FAIL'} ({end:+.4f})")
print(f"  (b) permuted-label null ~0 at EVERY n : {'PASS' if gb else 'FAIL'} (max |{S['null'].abs().max():.4f}|)")
if not (ga and gb): print("  -> gate FAILED : UNVERIFIED")
else:
    print(f"\n  n=1 {float(S.loc[1,'increment']):+.4f}  ->  n=8 {e8:+.4f} ({100*e8/end:.0f}% of n=31)  ->  n=31 {end:+.4f}")
    rise=(end-e16)/max(e16,1e-9)
    if e8>=0.8*end: print(f"  -> SATURATES BY 8 SOURCES : a phase-1 design needs about 8 domains, not 31")
    elif rise>0.20: print(f"  -> NOT SATURATED : still rising {100*rise:.0f}% from n=16 to n=31, so 32 blocks is not enough to see the ceiling")
    else: print(f"  -> flattening: n=16 to n=31 adds {100*rise:.0f}%")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
