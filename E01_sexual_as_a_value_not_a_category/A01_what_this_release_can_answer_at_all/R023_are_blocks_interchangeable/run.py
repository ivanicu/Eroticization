"""
E01 A01 R23 -- are source blocks interchangeable, as the sqrt(n) mechanism requires?

#60 fitted increment = 0.00723 * sqrt(n_sources) at CV 6.4%, and read sqrt(n) as the signature of
ONE SHARED LATENT measured with independent per-block noise. That reading makes a testable
prediction inside this release: if every block is an independent noisy estimate of the same thing,
blocks are INTERCHANGEABLE -- the increment depends on HOW MANY sources, not WHICH.

  interchangeable  -> variance across different SUBSETS at fixed n equals variance across SEEDS at
                      a fixed subset. All the spread is estimation noise.
  not              -> subset variance exceeds seed variance, and some block property predicts which
                      subsets do better. The independence assumption behind sqrt(n) fails.

ESTIMAND        the variance of the increment across source subsets at fixed n, decomposed against
                the variance across seeds at fixed subset.
IDENTIFICATION  identified; both variances are directly measurable by resampling one and holding
                the other.
WORLDS          A  interchangeable: variance ratio near 1
                B  structured: ratio >> 1, and a block property predicts subset quality
KILL (CONDITIONAL) gate, all measured (#41): (a) the mean increment at n=8 reproduces #60's +0.0225
                   within its spread; (b) permuted-label null ~0 for every subset.
                   then: subset/seed variance ratio < 1.5 -> INTERCHANGEABLE, mechanism holds
                         ratio > 3                        -> NOT, and sqrt(n) needs another reading
                         between                          -> partial, report the excess
POSITIVE CTRL   #60's n=8 value.
NEGATIVE CTRL   permuted person labels inside the cross-block factor.
SECOND TEST     regress subset increment on subset composition -- mean block size, mean option
                count, fluid-family share -- to see whether any property predicts subset quality.
SEEDS           3 per subset, 14 subsets.
MULTIPLICITY    14 subsets x 3 seeds x 2 arms, all reported.
IMPOSSIBLE      subsets larger than 31 or disjoint respondent pools; the gating tree shares people
                across every block.
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
FLUID={6,7,8,9,10,11,83}
M={}
for q in allq:
    s=lg[lg.qi==q]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=20].index)]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    X=np.zeros((len(ppl),len(opt))); X[s.person.map(pi).values,s.option.map(oi).values]=1
    if X.shape[0]>=800 and X.shape[1]>=8: M[q]=dict(X=X,ppl=ppl)
pool=np.unique(np.concatenate([v['ppl'] for v in M.values()]))
pmap={p:i for i,p in enumerate(pool)}
NS=8
def factors(srcs,K=8):
    cols=[]
    for q in srcs:
        v=M[q]; X=v['X']; idx=np.array([pmap[p] for p in v['ppl']])
        Z=np.full((len(pool),X.shape[1]),np.nan); Z[idx]=X
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Zc=np.hstack(cols); Zc=Zc-Zc.mean(0)
    U,S_,_=svd(Zc,full_matrices=False); return U[:,:K]*S_[:K]
def evaluate(subset,seed,permute=False):
    rng=np.random.default_rng(seed)
    D={'im':[],'pm':[],'cross':[]}; Y=[]
    for q,v in M.items():
        srcs=[s for s in subset if s!=q]
        if len(srcs)<NS-1: continue
        F=factors(srcs)[np.array([pmap[p] for p in v['ppl']])]
        X=v['X']; n,k=X.shape
        if permute: F=F[rng.permutation(n)]
        mask=rng.random((n,k))<0.30
        Xtr=np.where(mask,np.nan,X)
        im=np.nanmean(Xtr,axis=0); pm=np.nanmean(Xtr,axis=1)
        im=np.where(np.isfinite(im),im,np.nanmean(X)); pm=np.where(np.isfinite(pm),pm,np.nanmean(X))
        filled=np.where(np.isnan(Xtr),im[None,:],Xtr)
        A=np.c_[np.ones(n),F]; b,*_=lstsq(A,filled,rcond=None); crec=A@b
        ii,jj=np.where(mask); Y.append(X[ii,jj])
        D['im'].append(im[jj]); D['pm'].append(pm[ii]); D['cross'].append(crec[ii,jj])
    y=np.concatenate(Y); D={k:np.concatenate(v) for k,v in D.items()}
    h=rng.random(len(y))<0.5
    def sc(names):
        Xm=np.c_[np.ones(len(y)),np.column_stack([D[n] for n in names])]
        b,*_=lstsq(Xm[h],y[h],rcond=None); p=Xm[~h]@b
        return 1-((y[~h]-p)**2).sum()/((y[~h]-y[~h].mean())**2).sum()
    return sc(['im','pm','cross'])-sc(['im','pm'])
rng0=np.random.default_rng(11)
subsets=[list(rng0.choice(allq,NS,replace=False)) for _ in range(14)]
rows=[]
for si,sub in enumerate(subsets):
    for seed in (1,2,3):
        rows.append(dict(subset=si,seed=seed,inc=evaluate(sub,seed),
                         null=evaluate(sub,seed,permute=True),
                         mean_n=np.mean([M[q]['X'].shape[0] for q in sub]),
                         mean_k=np.mean([M[q]['X'].shape[1] for q in sub]),
                         fluid_share=np.mean([q in FLUID for q in sub])))
G=pd.DataFrame(rows); G.to_csv(OUT/'interchangeable.csv',index=False)
per=G.groupby('subset').inc.agg(['median','std'])
across=float(per['median'].std()); within=float(per['std'].mean())
print(f"\n=== {len(subsets)} random subsets of n={NS} sources, 3 seeds each ===")
print(f"  mean increment across subsets       : {per['median'].mean():+.4f}")
print(f"  SD ACROSS subsets (which blocks)    : {across:.5f}")
print(f"  SD WITHIN subset across seeds       : {within:.5f}")
print(f"  variance ratio (across/within)^2    : {(across/max(within,1e-9))**2:.2f}")
print(f"  permuted-label null                 : {G.null.median():+.5f}  (max |{G.null.abs().max():.5f}|)")
comp=G.groupby('subset').agg(inc=('inc','median'),mean_n=('mean_n','first'),
                             mean_k=('mean_k','first'),fluid=('fluid_share','first'))
print("\n=== does any subset property predict its increment? ===")
for c in ['mean_n','mean_k','fluid']:
    print(f"  corr(subset {c:8s}, increment) = {comp['inc'].corr(comp[c]):+.3f}")
ga=abs(per['median'].mean()-0.0225)<0.010; gb=abs(G.null).max()<0.004
ratio=(across/max(within,1e-9))**2
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) n=8 mean reproduces #60's +0.0225 : {'PASS' if ga else 'FAIL'} ({per['median'].mean():+.4f})")
print(f"  (b) permuted null ~0 for every subset : {'PASS' if gb else 'FAIL'}")
if not (ga and gb): print("  -> gate FAILED : UNVERIFIED")
elif ratio<1.5: print("  -> INTERCHANGEABLE : which blocks you use does not matter, only how many. The sqrt(n) mechanism holds.")
elif ratio>3: print(f"  -> NOT INTERCHANGEABLE : subset variance is {ratio:.1f}x seed variance; sqrt(n) needs another reading")
else: print(f"  -> partial : subset variance is {ratio:.1f}x seed variance")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
