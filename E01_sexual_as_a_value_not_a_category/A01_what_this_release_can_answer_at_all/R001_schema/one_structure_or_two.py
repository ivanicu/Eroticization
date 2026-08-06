"""
E01 A01 R20 -- are the person margin and the item margin one structure or two?

#57 measured the item margin for the first time in fifty-seven rounds: +0.0206 held-out R2 over a
marginals base, against the person factors' +0.0289 and a random-neighbour control's +0.0006. Each
was tested SEPARATELY against the base, so "71%" says item-alone recovers 71% of person-alone. It
does not say whether the item margin adds anything ON TOP of the person one.

Two margins can each be large and be the same structure seen from two sides. Fitting them jointly
decides it, and both outcomes are possible -- unlike #54's prediction, which was bounded by an
oracle and could not have been wrong.

  combined ~ max(item, person)   -> ONE structure, and this project has been describing it twice
  combined ~ item + person       -> TWO structures, and half of it has never been reported

ESTIMAND        the out-of-sample R2 increment of base+item+person over base, compared to each
                margin's own increment.
IDENTIFICATION  identified; nested models, common base, weights fitted on cells the scoring never
                sees.
KILL (CONDITIONAL) gate -- three conditions, all measured, none chosen (per #41):
                   (a) person increment > 0, the established margin reproduces
                   (b) random-neighbour increment near zero
                   (c) combined >= max(item, person) - seed spread; a combined model BELOW its own
                       components means the joint fit is broken, not that they cancel
                   then: combined <= max + 25% of min -> ONE STRUCTURE
                         combined >= max + 75% of min -> TWO STRUCTURES
                         between                      -> partial, report the overlap fraction
POSITIVE CTRL   the person margin.
NEGATIVE CTRL   random neighbours, and base+person+random (must not exceed base+person).
SEEDS           5.
MULTIPLICITY    6 models x 5 seeds, all reported.
IMPOSSIBLE      separating "same structure" from "two structures that happen to align" -- this
                measures shared explained variance, not shared mechanism.
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
def run(seed):
    rng=np.random.default_rng(seed)
    cols={'im':[],'pm':[],'item':[],'person':[],'rand':[]}; Y=[]
    for q,X in M.items():
        n,k=X.shape
        mask=rng.random((n,k))<0.30
        Xtr=np.where(mask,np.nan,X)
        im=np.nanmean(Xtr,axis=0); pm=np.nanmean(Xtr,axis=1)
        im=np.where(np.isfinite(im),im,np.nanmean(X)); pm=np.where(np.isfinite(pm),pm,np.nanmean(X))
        filled=np.where(np.isnan(Xtr),im[None,:],Xtr)
        C=np.corrcoef(filled.T); C=np.nan_to_num(C); np.fill_diagonal(C,-9)
        nbr=np.argsort(-C,axis=1)[:,:5]; rnd=rng.integers(0,k,size=(k,5))
        Z=filled-filled.mean(0); U,S_,Vt=svd(Z,full_matrices=False)
        rec=(U[:,:8]*S_[:8])@Vt[:8]+filled.mean(0)
        ii,jj=np.where(mask); Y.append(X[ii,jj])
        cols['im'].append(im[jj]); cols['pm'].append(pm[ii])
        cols['item'].append(filled[ii[:,None],nbr[jj]].mean(1))
        cols['rand'].append(filled[ii[:,None],rnd[jj]].mean(1))
        cols['person'].append(rec[ii,jj])
    y=np.concatenate(Y); D={k:np.concatenate(v) for k,v in cols.items()}
    h=rng.random(len(y))<0.5
    def score(names):
        Xm=np.c_[np.ones(len(y)),np.column_stack([D[n] for n in names])]
        b,*_=lstsq(Xm[h],y[h],rcond=None); p=Xm[~h]@b
        return 1-((y[~h]-p)**2).sum()/((y[~h]-y[~h].mean())**2).sum()
    B0=['im','pm']
    return dict(base=score(B0),item=score(B0+['item']),person=score(B0+['person']),
                both=score(B0+['item','person']),rand=score(B0+['rand']),
                person_rand=score(B0+['person','rand']),cells=len(y))
rows=[run(s) for s in (1,2,3,4,5)]
G=pd.DataFrame(rows); G.to_csv(OUT/'joint.csv',index=False)
Md=G.median(); Sd=G.agg(lambda s:s.max()-s.min())
inc={k:Md[k]-Md['base'] for k in ['item','person','both','rand','person_rand']}
print("\n=== out-of-sample R2 increments over the common base ===")
for k in ['person','item','both','rand','person_rand']:
    print(f"  base+{k:12s} {Md[k]:+.4f}   increment {inc[k]:+.4f}   seed spread {Sd[k]:.4f}")
print(f"  held-out cells {int(Md['cells']):,}")
mx=max(inc['item'],inc['person']); mn=min(inc['item'],inc['person'])
ga=inc['person']>0; gb=abs(inc['rand'])<max(0.3*abs(inc['person']),0.002)
gc=inc['both']>=mx-Sd['both']
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  (a) person increment > 0                : {'PASS' if ga else 'FAIL'} ({inc['person']:+.4f})")
print(f"  (b) random-neighbour near zero          : {'PASS' if gb else 'FAIL'} ({inc['rand']:+.4f})")
print(f"  (c) combined not below its components   : {'PASS' if gc else 'FAIL'} ({inc['both']:+.4f} vs max {mx:+.4f})")
if not (ga and gb and gc): print("  -> gate FAILED : UNVERIFIED")
else:
    extra=inc['both']-mx
    print(f"\n  combined {inc['both']:+.4f} = max {mx:+.4f} + {extra:+.4f}, where the smaller margin alone is {mn:+.4f}")
    print(f"  overlap: the smaller margin contributes {100*extra/mn:.0f}% of itself on top of the larger")
    if extra<=0.25*mn: print("  -> ONE STRUCTURE : this project has been describing it twice")
    elif extra>=0.75*mn: print("  -> TWO STRUCTURES : half of it has never been reported")
    else: print(f"  -> PARTIAL OVERLAP : {100*(1-extra/mn):.0f}% shared")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
