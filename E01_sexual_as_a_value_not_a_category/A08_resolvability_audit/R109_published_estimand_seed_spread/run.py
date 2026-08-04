"""
E01 A08 R03 -- resolvability of the estimand that actually carries the headline.

#36 found that A08 R01's audit measured a quantity nobody published: half-split factor CCA, where
the headline is PAIRWISE-BLOCK CCA -- 321 block pairs, nc=min(3,...), max over components, median
across pairs, reported raw (0.273) and demographics-removed (0.200). That estimand has never had
its seed spread measured, so #34's resolvability criterion has never been applied to it.

Reproduce the published computation exactly, at 6 seeds, and apply |effect| > 2 x spread.

ESTIMAND        the median across block pairs of the held-out max canonical correlation, raw and
                demographics-removed -- i.e. the two numbers in README's headline table.
IDENTIFICATION  identified; this is a re-run of a published computation with the seed varied.
WORLDS          A  the headline is a measurement: ratio well above 2
                B  it is seed noise: ratio at or below 2, as modality turned out to be
KILL (CONDITIONAL) gate: the positive control must exceed the observed value AND the permuted null
                   must fall below a third of it; otherwise UNVERIFIED.
                   then: ratio > 2 -> RESOLVABLE, the headline stands with a stated spread
                         ratio <= 2 -> UNRESOLVABLE, the headline is withdrawn like modality's
POSITIVE CTRL   a block split against ITSELF -- two halves of one block's own options. The same
                pipeline must return a high canonical correlation there, since it is the same people
                measured twice. If it does not, the pipeline cannot detect alignment at all.
NEGATIVE CTRL   person-permuted pairs: one block's rows shuffled before the fit.
NOISE FLOOR     the permuted value, recomputed at every seed.
MULTIPLICITY    2 adjustments x 6 seeds x ~320 pairs each, all summarised; per-seed values printed.
SEEDS           6.
IMPOSSIBLE      independent replication; a second release.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq
from sklearn.cross_decomposition import CCA
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
def strip(Mat, idx):
    D=np.c_[np.ones(len(idx)), COV.loc[idx].values]
    b,*_=lstsq(D,Mat,rcond=None); return Mat-D@b
def cvcca(Xa,Xb,nc,rng):
    idx=rng.permutation(len(Xa)); tr,te=idx[:len(idx)//2],idx[len(idx)//2:]
    try:
        c=CCA(n_components=nc,max_iter=800).fit(Xa[tr],Xb[tr]); ua,ub=c.transform(Xa[te],Xb[te])
        return max(abs(np.corrcoef(ua[:,j],ub[:,j])[0,1]) for j in range(nc))
    except Exception: return np.nan
def sweep(seed):
    rng=np.random.default_rng(seed)
    raw=[];adj=[];nul=[]
    for a,b in itertools.combinations(allq,2):
        A_,B_=B[a],B[b]; common=np.intersect1d(A_['ppl'],B_['ppl'])
        if len(common)<600: continue
        ia=np.searchsorted(A_['ppl'],common); ib=np.searchsorted(B_['ppl'],common)
        Ra,Rb=A_['R'][ia],B_['R'][ib]
        nc=min(3,Ra.shape[1]-1,Rb.shape[1]-1)
        raw.append(cvcca(Ra,Rb,nc,rng))
        adj.append(cvcca(strip(Ra,common),strip(Rb,common),nc,rng))
        nul.append(cvcca(Ra,Rb[rng.permutation(len(common))],nc,rng))
    return (float(np.nanmedian(raw)),float(np.nanmedian(adj)),float(np.nanmedian(nul)),len(raw))
rows=[]
for s in (1,2,3,4,5,6):
    r,a,n,k=sweep(s); rows.append(dict(seed=s,raw=r,adjusted=a,permuted_null=n,pairs=k))
G=pd.DataFrame(rows); G.to_csv(OUT/'published_estimand.csv',index=False)
print("=== the published estimand, one row per seed (321-pair medians) ===")
print(G.round(4).to_string(index=False))
def stat(col):
    v=G[col].values; return float(np.median(v)),float(v.max()-v.min())
print("\n=== resolvability ===")
res={}
for col,pub in [('raw',0.273),('adjusted',0.200),('permuted_null',None)]:
    med,spr=stat(col); ratio=abs(med)/max(spr,1e-9)
    res[col]=(med,spr,ratio)
    tag='' if pub is None else f'   published {pub}'
    print(f"  {col:14s} median {med:.4f}  seed spread {spr:.4f}  ratio {ratio:5.1f}  "
          f"{'RESOLVABLE' if abs(med)>2*spr else 'UNRESOLVABLE'}{tag}")
print("\n=== POSITIVE CONTROL: a block split against itself ===")
rng=np.random.default_rng(99); pcs=[]
for q in allq[:10]:
    R=B[q]['R']; k=R.shape[1]
    if k<8: continue
    p=rng.permutation(k); h=k//2
    pcs.append(cvcca(R[:,p[:h]],R[:,p[h:2*h]],min(3,h-1),rng))
print(f"   same block, two halves of its own options: median {np.nanmedian(pcs):.4f}  (n={len(pcs)} blocks)")
gate_pos=float(np.nanmedian(pcs))>res['raw'][0]
gate_neg=res['permuted_null'][0]<res['raw'][0]/3
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control exceeds observed : {'PASS' if gate_pos else 'FAIL'} ({np.nanmedian(pcs):.3f} vs {res['raw'][0]:.3f})")
print(f"  permuted null below a third       : {'PASS' if gate_neg else 'FAIL'} ({res['permuted_null'][0]:.3f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
else:
    ok=res['adjusted'][2]>2
    print(f"  -> adjusted ratio {res['adjusted'][2]:.1f} : {'RESOLVABLE, headline stands with spread ±%.4f'%res['adjusted'][1] if ok else 'UNRESOLVABLE, withdraw like modality'}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
