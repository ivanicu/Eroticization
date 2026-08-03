"""
E01 A01 R21 -- within-block structure and cross-domain transfer are different objects, and this
project has been quoting them as one.

#58 closed on an apparent contradiction: this within-block structure predicts at ratio 12
(+0.0293), while #49 showed the CROSS-DOMAIN direction carries no predictable variance (pairwise
block->block R2 -0.0022). Reading the code rather than remembering it resolves half of it
immediately -- R19/R20 fit the SVD INSIDE each block (`Z=filled-filled.mean(0)` within the loop),
so their "person margin" is a block-internal reconstruction, never cross-block transfer. The two
numbers were never measuring the same thing and I quoted them in one breath.

What is still open is the decomposition: of the within-block predictable structure, how much does a
factor built ONLY FROM OTHER BLOCKS supply, and how much is irreducibly block-internal?

ESTIMAND        out-of-sample R2 increments over a marginals base from (a) cross-block factors
                fitted on OTHER blocks only, (b) within-block low-rank reconstruction, (c) both.
IDENTIFICATION  identified; the cross-block factor never sees the target block, by construction.
WORLDS          A  one structure: cross-block factors supply most of the within-block increment
                B  two objects: within-block dominates and cross-block adds little, which is what
                   #49's thin direction predicts
KILL (CONDITIONAL) gate, all measured not chosen (#41): (a) within-block increment > 0 -- the
                   established quantity reproduces; (b) random control near zero; (c) combined not
                   below its components.
                   then: cross-block increment >= 50% of within-block -> ONE STRUCTURE
                         cross-block <= 20% of within-block           -> TWO OBJECTS
                         between                                      -> report the fraction
POSITIVE CTRL   the within-block reconstruction, ratio 12 in #58.
NEGATIVE CTRL   random-neighbour column, and cross-block factors with person labels permuted.
SEEDS           5.
MULTIPLICITY    5 models x 5 seeds, all reported.
IMPOSSIBLE      a block whose people are disjoint from every other block -- the gating tree makes
                overlap unavoidable, so "cross-block" always shares respondents. Stated, not fixed.
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
print(f"blocks {len(M)}  pool {len(pool):,}")
def cross_factors(exclude,K=8):
    cols=[]
    for q,v in M.items():
        if q==exclude: continue
        X=v['X']; idx=np.array([pmap[p] for p in v['ppl']])
        Z=np.full((len(pool),X.shape[1]),np.nan); Z[idx]=X
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Zc=np.hstack(cols); Zc=Zc-Zc.mean(0)
    U,S_,_=svd(Zc,full_matrices=False); return U[:,:K]*S_[:K]
def run(seed):
    rng=np.random.default_rng(seed)
    D={'im':[],'pm':[],'within':[],'cross':[],'crossperm':[],'rand':[]}; Y=[]
    for q,v in M.items():
        X=v['X']; n,k=X.shape
        F=cross_factors(q)[np.array([pmap[p] for p in v['ppl']])]
        mask=rng.random((n,k))<0.30
        Xtr=np.where(mask,np.nan,X)
        im=np.nanmean(Xtr,axis=0); pm=np.nanmean(Xtr,axis=1)
        im=np.where(np.isfinite(im),im,np.nanmean(X)); pm=np.where(np.isfinite(pm),pm,np.nanmean(X))
        filled=np.where(np.isnan(Xtr),im[None,:],Xtr)
        Z=filled-filled.mean(0); U,S_,Vt=svd(Z,full_matrices=False)
        rec=(U[:,:8]*S_[:8])@Vt[:8]+filled.mean(0)
        A=np.c_[np.ones(n),F]; b,*_=lstsq(A[~mask.any(1)] if False else A,filled,rcond=None)
        crec=A@b
        Fp=F[rng.permutation(n)]; Ap=np.c_[np.ones(n),Fp]
        bp,*_=lstsq(Ap,filled,rcond=None); cperm=Ap@bp
        C=np.corrcoef(filled.T); C=np.nan_to_num(C); np.fill_diagonal(C,-9)
        rnd=rng.integers(0,k,size=(k,5))
        ii,jj=np.where(mask); Y.append(X[ii,jj])
        D['im'].append(im[jj]); D['pm'].append(pm[ii])
        D['within'].append(rec[ii,jj]); D['cross'].append(crec[ii,jj])
        D['crossperm'].append(cperm[ii,jj]); D['rand'].append(filled[ii[:,None],rnd[jj]].mean(1))
    y=np.concatenate(Y); D={k:np.concatenate(v) for k,v in D.items()}
    h=rng.random(len(y))<0.5
    def sc(names):
        Xm=np.c_[np.ones(len(y)),np.column_stack([D[n] for n in names])]
        b,*_=lstsq(Xm[h],y[h],rcond=None); p=Xm[~h]@b
        return 1-((y[~h]-p)**2).sum()/((y[~h]-y[~h].mean())**2).sum()
    B0=['im','pm']
    return dict(base=sc(B0),within=sc(B0+['within']),cross=sc(B0+['cross']),
                both=sc(B0+['within','cross']),crossperm=sc(B0+['crossperm']),
                rand=sc(B0+['rand']),cells=len(y))
rows=[run(s) for s in (1,2,3,4,5)]
G=pd.DataFrame(rows); G.to_csv(OUT/'within_vs_cross.csv',index=False)
Md=G.median(); Sd=G.agg(lambda s:s.max()-s.min())
inc={k:Md[k]-Md['base'] for k in ['within','cross','both','crossperm','rand']}
print("\n=== increments over the marginals base ===")
for k in ['within','cross','both','crossperm','rand']:
    print(f"  base+{k:10s} {Md[k]:+.4f}   increment {inc[k]:+.4f}   spread {Sd[k]:.4f}")
print(f"  held-out cells {int(Md['cells']):,}")
ga=inc['within']>0; gb=abs(inc['rand'])<0.003 and abs(inc['crossperm'])<max(0.3*inc['within'],0.003)
gc=inc['both']>=max(inc['within'],inc['cross'])-Sd['both']
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) within-block increment > 0        : {'PASS' if ga else 'FAIL'} ({inc['within']:+.4f})")
print(f"  (b) both null controls near zero      : {'PASS' if gb else 'FAIL'} (rand {inc['rand']:+.4f}, cross-perm {inc['crossperm']:+.4f})")
print(f"  (c) combined not below its components : {'PASS' if gc else 'FAIL'} ({inc['both']:+.4f})")
if not (ga and gb and gc): print("  -> gate FAILED : UNVERIFIED")
else:
    frac=inc['cross']/inc['within']
    print(f"\n  cross-block increment is {100*frac:.0f}% of the within-block increment")
    if frac>=0.5: print("  -> ONE STRUCTURE")
    elif frac<=0.2: print("  -> TWO OBJECTS : within-block structure is not cross-domain transfer, and this project quoted them as one")
    else: print(f"  -> partial: {100*frac:.0f}%")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
