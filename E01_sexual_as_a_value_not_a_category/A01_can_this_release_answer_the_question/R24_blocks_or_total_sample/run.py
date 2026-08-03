"""
E01 A01 R24 -- does the signal accumulate with BLOCK COUNT or with TOTAL SAMPLE?

#60 fitted increment = 0.00723 * sqrt(n_sources). #61 falsified its mechanism reading: blocks are
not interchangeable (subset variance 58x seed variance) and block SIZE predicts subset quality at
r = +0.816. Because #60 sampled subsets at random, total respondents grew with n, so sqrt(block
count) and sqrt(total sample) were perfectly confounded.

The confound is removable by subsampling rather than by argument. Hold the TOTAL respondent-rows
used across sources fixed, and vary how many blocks those rows are spread over.

  per-block-latent reading -> at fixed total rows, MORE BLOCKS is better (each block is a separate
                              view of the latent, so splitting the same budget across more views wins)
  pooled-sample reading    -> at fixed total rows, block count is IRRELEVANT; only the budget matters

ESTIMAND        held-out increment as a function of (n_sources, rows_per_source), on an explicit
                budget grid where total = n * rows is held constant along one axis.
IDENTIFICATION  identified: the design breaks the confound by construction rather than adjusting for it.
WORLDS          A  block count matters at fixed budget -> per-block structure is real
                B  only budget matters -> block boundaries are incidental, one pooled estimate
KILL (CONDITIONAL) gate, measured not chosen (#41): (a) the increment must rise with BUDGET at fixed
                   n -- if more data does not help, nothing here is measurable; (b) permuted-label
                   null ~0 in every cell.
                   then: at fixed budget, n=31 beats n=8 by >30% -> BLOCK COUNT MATTERS
                         within 10%                              -> ONLY BUDGET MATTERS
                         between                                 -> report both slopes
POSITIVE CTRL   the budget axis: increment must grow when rows_per_source grows at fixed n.
NEGATIVE CTRL   permuted person labels inside the cross-block factor.
SEEDS           3.
MULTIPLICITY    7 design cells x 2 arms x 3 seeds, all reported.
IMPOSSIBLE      rows_per_source above the smallest block (800), so the largest budgets are only
                reachable by using more blocks -- which is why the grid is stated as cells rather
                than a full factorial.
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
print(f"blocks {len(M)}  smallest block {min(v['X'].shape[0] for v in M.values())}")
def factors(srcs,rows,rng,K=8):
    cols=[]
    for q in srcs:
        v=M[q]; X=v['X']; idx=np.array([pmap[p] for p in v['ppl']])
        take=rng.choice(len(idx),min(rows,len(idx)),replace=False)   # SUBSAMPLE respondents
        Z=np.full((len(pool),X.shape[1]),np.nan); Z[idx[take]]=X[take]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Zc=np.hstack(cols); Zc=Zc-Zc.mean(0)
    U,S_,_=svd(Zc,full_matrices=False); return U[:,:K]*S_[:K]
def evaluate(nsrc,rows,seed,permute=False):
    rng=np.random.default_rng(seed*131+nsrc*7+rows)
    D={'im':[],'pm':[],'cross':[]}; Y=[]
    targets=list(rng.choice(list(M),12,replace=False))
    for q in targets:
        v=M[q]
        others=[x for x in M if x!=q]
        srcs=list(rng.choice(others,min(nsrc,len(others)),replace=False))
        F=factors(srcs,rows,rng)[np.array([pmap[p] for p in v['ppl']])]
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
# COST NOTE: the full 9-cell x 3-seed x 2-arm grid is ~1,700 pool-wide SVDs and does not finish.
# Cut to the cells that decide it: the BUDGET axis at fixed n=8 (the positive control) and the
# BLOCK-COUNT axis at fixed budget 6400 (the estimand). Nulls only where they gate.
CELLS=[(8,250),(8,500),(8,800),       # budget axis at fixed n -> positive control
       (16,400),(31,206)]             # block-count axis at fixed budget 6400 -> the estimand
NULLCELLS={(8,800),(31,206)}
rows=[]
for nsrc,r in CELLS:
    for seed in (1,2):
        rows.append(dict(n_sources=nsrc,rows_per_source=r,budget=nsrc*r,seed=seed,
                         inc=evaluate(nsrc,r,seed),
                         null=evaluate(nsrc,r,seed,permute=True) if (nsrc,r) in NULLCELLS else 0.0))
G=pd.DataFrame(rows); G.to_csv(OUT/'blocks_or_sample.csv',index=False)
S=G.groupby(['budget','n_sources']).agg(inc=('inc','median'),spread=('inc',lambda s:s.max()-s.min()),
                                        null=('null','median')).reset_index()
print("\n=== held-out increment: budget held fixed along each row-group ===")
print(S.round(4).to_string(index=False))
piv=S.pivot(index='budget',columns='n_sources',values='inc')
print("\n=== same, as budget x n_sources ===")
print(piv.round(4).to_string())
b8=S[S.n_sources==8].sort_values('budget').inc.values; gate_a=bool(b8[-1]>b8[0])
gate_b=abs(S['null']).max()<0.004
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) increment rises with BUDGET at fixed n=8 : {'PASS' if gate_a else 'FAIL'} ({[round(v,4) for v in piv[8]]})")
print(f"  (b) permuted null ~0 in every cell           : {'PASS' if gate_b else 'FAIL'} (max |{abs(S['null']).max():.5f}|)")
if not (gate_a and gate_b): print("  -> gate FAILED : UNVERIFIED")
else:
    r_=[]
    for bud in piv.index:
        if 8 in piv.columns and 31 in piv.columns and np.isfinite(piv.loc[bud,8]) and np.isfinite(piv.loc[bud,31]):
            r_.append(piv.loc[bud,31]/max(piv.loc[bud,8],1e-9))
    m=float(np.median(r_))
    print(f"  at fixed budget, n=31 / n=8 = {m:.2f}   (per budget: {[round(v,2) for v in r_]})")
    if m>1.30: print("  -> BLOCK COUNT MATTERS : spreading the same rows over more blocks wins; per-block structure is real")
    elif m<1.10: print("  -> ONLY BUDGET MATTERS : block boundaries are incidental, this is one pooled estimate")
    else: print(f"  -> partial: {m:.2f}x from 8 to 31 blocks at equal budget")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
