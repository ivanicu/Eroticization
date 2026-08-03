"""
E01 A01 R19 -- the item margin, third attempt, built the way the first two should have been.

#56: attempt 1 compared a bare item-neighbour average against a baseline carrying both marginals
(not information-matched, and it printed "EMPTY"); attempt 2 nested residuals without fitting their
weights (the positive control fell below baseline and the gate refused). Neither measured anything.

This design fixes both faults by construction:
  - CELL-LEVEL masking, not person-level. Every feature is computed from UNMASKED cells only.
  - All four models NESTED on the same base and JOINTLY FITTED on training cells, scored on the
    identical held-out cells. Each is base + one candidate, so the increment is what is compared.

ESTIMAND        the out-of-sample R2 increment over a marginals-only base, from (a) item-neighbour
                structure and (b) person-factor structure, on identical held-out cells.
IDENTIFICATION  identified; leakage is impossible because features never see a masked cell.
WORLDS          A  item margin carries structure the person side missed: item increment > person
                B  redundant: the two increments are comparable
                C  empty: item increment at the random-neighbour control
KILL (CONDITIONAL) gate -- both controls, per #33: the PERSON increment must be positive (the
                   established margin must reproduce) AND the random-neighbour increment must be
                   near zero. Either failing -> UNVERIFIED, no comparison.
                   then: item > 1.5x person -> NEW STRUCTURE
                         item within 25% of person -> REDUNDANT
                         item <= random-neighbour  -> EMPTY
POSITIVE CTRL   the person-factor increment, the margin measured exhaustively for 50 rounds.
NEGATIVE CTRL   item neighbours chosen at random, same count.
BASE            item marginal + person rate, both from unmasked cells.
SEEDS           4.
MULTIPLICITY    4 models x 4 seeds, every block contributing, all reported.
IMPOSSIBLE      cross-block item neighbours for gated pairs -- most option pairs are never observed
                together because the tree routes people apart. Neighbours are within-block, stated.
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
    if X.shape[0]>=800 and X.shape[1]>=8: M[q]=X
print(f"blocks usable {len(M)}   cells {sum(v.size for v in M.values()):,}")
def run(seed):
    rng=np.random.default_rng(seed)
    F={'base':[],'item':[],'person':[],'rand':[]}; Y=[]
    for q,X in M.items():
        n,k=X.shape
        mask=rng.random((n,k))<0.30                    # held-out CELLS
        Xtr=np.where(mask,np.nan,X)
        im=np.nanmean(Xtr,axis=0); pm=np.nanmean(Xtr,axis=1)
        im=np.where(np.isfinite(im),im,np.nanmean(X)); pm=np.where(np.isfinite(pm),pm,np.nanmean(X))
        filled=np.where(np.isnan(Xtr),im[None,:],Xtr)
        C=np.corrcoef(filled.T); C=np.nan_to_num(C); np.fill_diagonal(C,-9)
        nbr=np.argsort(-C,axis=1)[:,:5]
        rnd=rng.integers(0,k,size=(k,5))
        Z=filled-filled.mean(0); U,S_,Vt=svd(Z,full_matrices=False)
        rec=(U[:,:8]*S_[:8])@Vt[:8]+filled.mean(0)
        ii,jj=np.where(mask)
        Y.append(X[ii,jj])
        F['base'].append(np.c_[im[jj],pm[ii]])
        F['item'].append(np.c_[im[jj],pm[ii],filled[ii][:,None][:,0,:][np.arange(len(ii))[:,None],nbr[jj]].mean(1)])
        F['rand'].append(np.c_[im[jj],pm[ii],filled[ii][:,None][:,0,:][np.arange(len(ii))[:,None],rnd[jj]].mean(1)])
        F['person'].append(np.c_[im[jj],pm[ii],rec[ii,jj]])
    y=np.concatenate(Y)
    out={}
    for kname,parts in F.items():
        Xm=np.vstack(parts); Xm=np.c_[np.ones(len(Xm)),Xm]
        h=rng.random(len(y))<0.5                        # fit the weights on half the held-out cells,
        b,*_=lstsq(Xm[h],y[h],rcond=None)               # score on the other half -- weights never see it
        p=Xm[~h]@b
        out[kname]=1-((y[~h]-p)**2).sum()/((y[~h]-y[~h].mean())**2).sum()
    out['cells']=len(y); return out
rows=[run(s) for s in (1,2,3,4)]
G=pd.DataFrame(rows); G.to_csv(OUT/'item_margin_v3.csv',index=False)
Md=G.median(); Sd=G.agg(lambda s:s.max()-s.min())
print("\n=== out-of-sample R2, all models jointly fitted on the same base ===")
for k in ['base','person','item','rand']:
    print(f"  {k:8s} {Md[k]:+.4f}   seed spread {Sd[k]:.4f}   increment over base {Md[k]-Md['base']:+.4f}")
print(f"  held-out cells: {int(Md['cells']):,}")
ip=Md['person']-Md['base']; ii_=Md['item']-Md['base']; ir=Md['rand']-Md['base']
gate_pos=ip>0; gate_neg=abs(ir)<max(0.3*abs(ip),0.002)
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  person increment positive (established margin reproduces) : {'PASS' if gate_pos else 'FAIL'} ({ip:+.4f})")
print(f"  random-neighbour increment near zero                      : {'PASS' if gate_neg else 'FAIL'} ({ir:+.4f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
elif ii_>1.5*ip: print(f"  -> NEW STRUCTURE : item increment {ii_:+.4f} is {ii_/ip:.1f}x the person increment")
elif abs(ii_-ip)<0.25*ip: print(f"  -> REDUNDANT : item {ii_:+.4f} vs person {ip:+.4f}")
elif ii_<=ir: print("  -> EMPTY : the item margin adds nothing over random neighbours")
else: print(f"  -> partial : item increment is {100*ii_/ip:.0f}% of the person increment")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
