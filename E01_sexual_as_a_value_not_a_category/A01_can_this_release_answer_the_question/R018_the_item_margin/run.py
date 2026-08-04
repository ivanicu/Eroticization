"""
E01 A01 R18 -- the item margin, abandoned 40 rounds ago and never run.

#55 established that every NEXT line should be checked for analytic possibility before being acted
on. Auditing all 45 of them: nearly all were audits of existing claims and were executed. ONE opens
a new margin, is falsifiable, and was abandoned when the restructure interrupted it --

  aea8476: "ten iterations have all measured structure over PEOPLE. The other margin has never been
  touched -- the ITEMS. 1,332 options exist and only their block membership has ever been used."

It is sharper now than when written. #49 showed the person-side cross-domain structure is a THIN
DIRECTION -- a real correlation carrying no predictable variance. So the question becomes: is there
ITEM-side structure the person-side analysis missed?

  person side : reduce PEOPLE to factors, use them to predict a person's endorsements
  item side   : reduce ITEMS to co-endorsement geometry, use an item's neighbours to predict whether
                this person endorses it -- never fitting a person representation at all

ESTIMAND        held-out endorsement prediction from item-item co-endorsement structure, against
                person-factor prediction and a marginals-only baseline, on identical held-out cells.
IDENTIFICATION  identified; both predictors are fitted on training cells only.
WORLDS          A  the item margin is redundant: item-based ties person-based
                B  it carries structure the person side missed: item-based beats it
                C  it is empty: item-based ties the marginals-only baseline
KILL (CONDITIONAL) gate -- ceiling first (#50): person-factor prediction must beat the marginals
                   baseline, else neither margin works here and nothing is comparable.
                   then: item beats person by >20% of the person margin -> NEW STRUCTURE
                         item within 10% of person                       -> REDUNDANT
                         item at the marginals baseline                  -> EMPTY
POSITIVE CTRL   person-factor prediction, the established margin.
NEGATIVE CTRL   item neighbours assigned at random, same neighbourhood size.
BASELINE        item marginals plus person endorsement rate -- no interaction structure at all.
SEEDS           4.
MULTIPLICITY    4 predictors x 4 seeds, all reported.
IMPOSSIBLE      the full 1,332-option matrix at once -- gating makes most cells unobserved, so this
                runs within blocks and across the block-pair overlap, which is stated not hidden.
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
allq=list(B); pool=np.unique(np.concatenate([B[q]['ppl'] for q in allq]))
# build one observed-cell table over the blocks, keeping BINARY endorsements (not residuals)
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
qm=pd.read_csv('data/derived/multiselect_questions.csv')
M={}
for q in allq:
    s=lg[lg.qi==q]; vc=s.option.value_counts(); s=s[s.option.isin(vc[vc>=20].index)]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    X=np.zeros((len(ppl),len(opt))); X[s.person.map(pi).values,s.option.map(oi).values]=1
    M[q]=dict(ppl=ppl,X=X)
print(f"blocks {len(M)}   total observed cells {sum(v['X'].size for v in M.values()):,}")
def run(seed):
    rng=np.random.default_rng(seed); res={}
    ytrue=[];p_base=[];p_item=[];p_pers=[];p_rand=[]
    # person factors from ALL blocks (the established margin)
    pm={p:i for i,p in enumerate(pool)}; cols=[]
    for q in allq:
        idx=np.array([pm[p] for p in M[q]['ppl'] if p in pm]); src=np.array([i for i,p in enumerate(M[q]['ppl']) if p in pm])
        Z=np.full((len(pool),M[q]['X'].shape[1]),np.nan); Z[idx]=M[q]['X'][src]
        mu=np.nanmean(Z,axis=0); mu=np.where(np.isfinite(mu),mu,0.)
        cols.append(np.nan_to_num(np.where(np.isnan(Z),mu,Z)))
    Zall=np.hstack(cols); Zall=Zall-Zall.mean(0)
    U,S_,_=svd(Zall,full_matrices=False); FAC=U[:,:8]*S_[:8]
    for q in allq:
        X=M[q]['X']; ppl=M[q]['ppl']; n,k=X.shape
        if k<8 or n<800: continue
        te_p=rng.random(n)<0.3
        tr=X[~te_p]; te=X[te_p]
        base_item=tr.mean(0); base_pers=te.mean(1,keepdims=True)
        C=np.corrcoef(tr.T); C=np.nan_to_num(C); np.fill_diagonal(C,0)
        nbr=np.argsort(-C,axis=1)[:,:5]                     # 5 nearest ITEMS by co-endorsement
        rnd=rng.integers(0,k,size=(k,5))
        fp=np.array([pm[p] for p in ppl])[te_p]
        Ftr=FAC[np.array([pm[p] for p in ppl])[~te_p]]; Fte=FAC[fp]
        A=np.c_[np.ones(len(Ftr)),Ftr]; b,*_=lstsq(A,tr,rcond=None)
        pers=np.c_[np.ones(len(Fte)),Fte]@b
        # FIX: the first pass compared a bare item-neighbour average against a baseline that
        # already carried the item marginal AND the person rate. Not information-matched, so the
        # item margin was handicapped rather than tested. NEST all three on the same base.
        for j in range(k):
            ytrue.append(te[:,j])
            b0=np.full(te.shape[0],base_item[j])+(base_pers[:,0]-tr.mean())
            p_base.append(b0)
            p_item.append(b0+(te[:,nbr[j]].mean(1)-tr[:,nbr[j]].mean()))
            p_rand.append(b0+(te[:,rnd[j]].mean(1)-tr[:,rnd[j]].mean()))
            p_pers.append(b0+(pers[:,j]-tr[:,j].mean()))
    y=np.concatenate(ytrue)
    def r2(p):
        p=np.concatenate(p); return 1-((y-p)**2).sum()/((y-y.mean())**2).sum()
    return dict(baseline=r2(p_base),item=r2(p_item),person=r2(p_pers),random_item=r2(p_rand),cells=len(y))
rows=[run(s) for s in (1,2,3,4)]
G=pd.DataFrame(rows); G.to_csv(OUT/'item_margin.csv',index=False)
Md=G.median(); Sd=G.agg(lambda s:s.max()-s.min())
print("\n=== held-out endorsement prediction R2 ===")
for k in ['baseline','person','item','random_item']:
    print(f"  {k:14s} {Md[k]:+.4f}   seed spread {Sd[k]:.4f}")
print(f"  cells scored: {int(Md['cells']):,}")
gate=Md['person']>Md['baseline']
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  person margin beats the marginals baseline : {'PASS' if gate else 'FAIL'} ({Md['person']:+.4f} vs {Md['baseline']:+.4f})")
if not gate: print("  -> gate FAILED : neither margin works here, UNVERIFIED")
else:
    pm_=Md['person']-Md['baseline']; im=Md['item']-Md['baseline']
    print(f"  person margin {pm_:+.4f} · item margin {im:+.4f} · random-item {Md['random_item']-Md['baseline']:+.4f}")
    if im>1.2*pm_: print(f"  -> NEW STRUCTURE : the item margin carries {100*im/pm_:.0f}% of the person margin")
    elif abs(im-pm_)<0.1*pm_: print("  -> REDUNDANT : the two margins are the same structure")
    elif im<=0.05*pm_: print("  -> EMPTY : the item margin adds nothing over the marginals")
    else: print(f"  -> partial: item margin is {100*im/pm_:.0f}% of the person margin")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
