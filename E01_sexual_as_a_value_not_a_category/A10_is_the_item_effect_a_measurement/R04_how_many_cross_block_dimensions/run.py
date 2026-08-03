import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R04 -- HOW MANY CROSS-BLOCK DIMENSIONS, MEASURED BY PREDICTION ON BLOCKS NEVER SEEN.

R03 left C still rising at Kc=8. #18 and #49 both tried to settle the coordinate count from WITHIN
the shared space -- congruence between halves, canonical correlations -- and both failed, one on
K-dependence and one because the direction was thin. This estimator is different in kind: the person
scores come from 31 blocks, the loadings from the target's training cells, and the score is held-out
prediction on the target's masked cells. A dimension only counts if it PREDICTS somewhere it was
never fit.

A monotone rising curve is not by itself an answer, because more dimensions means more loadings and
a rising curve could be the estimator flattering itself. So the curve is only readable against a
world whose TRUE rank is known.

ESTIMAND        the largest Kc at which C_corrected still gains more than its own seed spread; and
                the same quantity measured on synthetic worlds of KNOWN rank 2 and 5.
IDENTIFICATION  identified only relative to the control: the answer is "the real curve behaves like
                a world of true rank r", never a bare number off the real curve.
SCOPE           23 identified blocks (real); 8 target blocks for each control world (cost).
WORLDS          low-dim   the curve saturates by Kc~5, like the r=5 control
                high-dim  the curve keeps rising past the controls' saturation
                artifact  the controls ALSO rise without saturating -> the estimator cannot count
                          dimensions at all and no number is reportable
KILL            the r=2 and r=5 controls must SATURATE near their true rank. If they do not, the
                real curve is uninterpretable and this round reports UNVERIFIED -- which is the
                outcome that would kill the whole question, and it is the reason the controls are
                here rather than in a follow-up.
POSITIVE CTRL   the r=5 world: saturation must occur later than in the r=2 world (dose-response
                in the true rank, not just "it saturates").
NEGATIVE CTRL   person-permutation at every Kc, on the real data.
NOISE FLOOR     2 masks; per-Kc spread reported.
MULTIPLICITY    9 ranks x 23 blocks x 2 seeds real, plus 9 x 8 x 2 for each of 2 control worlds.
IMPOSSIBLE      naming the dimensions -- #18/#39/#49 measured that naming fails here. Counting and
                naming are different tasks and only the first is attempted.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd, lstsq
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW])); PM={p:i for i,p in enumerate(ALLP)}
MASK=0.15; SEEDS=[11,29]; KCS=[1,2,3,5,8,12,16,24,32]
print(f"targets {len(IDENT)}  people {len(ALLP)}  ranks {KCS}",flush=True)

def build_scores(BLK,target,K):
    cols=[]
    for q in BLK:
        if q==target: continue
        M=BLK[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in BLK[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,Vt=svd(Z,full_matrices=False)
    return U[:,:K]*S[:K]

def cross_r2(BLK,target,U_all,Kc,seed,permute=False):
    M=BLK[target]['M']; n,m=M.shape; rows=[PM[p] for p in BLK[target]['ppl']]
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    U=U_all[rows][:,:Kc]
    if permute: U=U[np.random.default_rng(seed+555).permutation(n)]
    U=(U-U.mean(0))/(U.std(0)+1e-12)
    C=np.zeros_like(M)
    for j in range(m):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),U[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(n),U]@b
    base=np.mean((M[he]-gm)**2)
    f=lambda *parts: 1.-np.mean((M[he]-np.clip(gm+sum(parts),0,1)[he])**2)/base
    IB=np.broadcast_to(I,M.shape)
    return f(IB,P,C)-f(IB,P)            # the increment C adds on top of both marginals

def synth_world(r,seed):
    """all blocks regenerated from a shared person space of TRUE rank r, matched shapes+prevalences."""
    rng=np.random.default_rng(3000+seed)
    F=rng.normal(size=(len(ALLP),r))
    out={}
    for q in RAW:
        M=RAW[q]['M']; n,m=M.shape; rows=[PM[p] for p in RAW[q]['ppl']]
        L=rng.normal(size=(r,m))*0.30
        p=np.clip(M.mean(0)[None,:]+F[rows]@L,0.02,0.98)
        out[q]=dict(M=(rng.random((n,m))<p).astype(float),ppl=RAW[q]['ppl'])
    return out

rows=[]
def sweep(BLK,targets,tag,arms=('real','perm')):
    for i,t in enumerate(targets):
        U=build_scores(BLK,t,max(KCS))
        for Kc in KCS:
            for sd in SEEDS:
                for a in arms:
                    rows.append(dict(world=tag,q=t,Kc=Kc,seed=sd,arm=a,
                                     C=cross_r2(BLK,t,U,Kc,sd,permute=(a=='perm'))))
        print(f"  [{tag}] {i+1}/{len(targets)}",flush=True)

CTRL_T=IDENT[:8]
sweep(RAW,IDENT,'real')
for r in [2,5]:
    sweep(synth_world(r,1),CTRL_T,f'r{r}',arms=('real',))
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R04_how_many_cross_block_dimensions/results/'
D.to_csv(OUT+'grid.csv',index=False)

def curve(w,arm='real',ts=None):
    d=D[(D.world==w)&(D.arm==arm)]
    if ts is not None: d=d[d.q.isin(ts)]
    return d.groupby('Kc').C.mean(), d.groupby(['Kc','q']).C.mean().groupby('Kc').std()

print("\n=== THE CURVES ===")
cr,_=curve('real'); cp,_=curve('real','perm')
c2,_=curve('r2'); c5,_=curve('r5')
crc=cr-cp
T=pd.DataFrame({'real':cr,'perm':cp,'real_corrected':crc,'ctrl_r2':c2,'ctrl_r5':c5})
T['gain']=T.real_corrected.diff()
sp=D[(D.world=='real')&(D.arm=='real')].groupby(['Kc','q']).C.std().groupby('Kc').mean()
T['seed_spread']=sp; T['gain_over_spread']=(T.gain/sp).round(2)
print(T.round(4).to_string())

def sat(c,tol=0.05):
    """first rank whose remaining gain is below tol of the total gain achieved."""
    tot=c.max()-c.iloc[0]
    for k in c.index:
        if (c.max()-c[k])<=tol*max(tot,1e-12): return k
    return c.index[-1]
s2,s5,sr=sat(c2),sat(c5),sat(crc)
print(f"\n=== CONTROL CALIBRATION (saturation rank at 95% of total gain) ===")
print(f"   true rank 2 world -> saturates at Kc={s2}")
print(f"   true rank 5 world -> saturates at Kc={s5}")
print(f"   dose-response in the true rank: {'PASS' if s5>s2 else 'FAIL -- the estimator cannot count'}")
print(f"   controls saturate before the sweep ends: {'PASS' if max(s2,s5)<max(KCS) else 'FAIL'}")
print(f"\n=== THE REAL DATA ===")
print(f"   saturates at Kc={sr}   corrected C at saturation {crc[sr]:+.4f}   at Kc={max(KCS)} {crc.iloc[-1]:+.4f}")
last=T.gain_over_spread.iloc[-1]
print(f"   final increment / seed spread = {last}")
print("\n  CONDITIONAL KILL")
if not(s5>s2 and max(s2,s5)<max(KCS)):
    print("   -> UNVERIFIED: the controls do not calibrate the estimator, so the real curve carries")
    print("      no rank information. No number is reportable and #18/#49 remain unsettled.")
else:
    print(f"   -> the estimator CAN count (r=2 saturates at {s2}, r=5 at {s5}).")
    if sr>s5:
        print(f"      The real cross-block space saturates LATER than a true-rank-5 world (Kc={sr} vs {s5}):")
        print(f"      the eroticization readout carries MORE than 5 domain-general dimensions.")
    else:
        print(f"      The real space saturates at Kc={sr}, i.e. between the r=2 and r=5 controls.")
    print(f"      Reported as a comparison to known-rank worlds, never as a bare count.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
