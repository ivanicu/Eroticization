import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R12 -- SETTLING #83c's OWN CONFOUND: WAS "GENERAL IS HIGH-DIMENSIONAL" JUST A BIGGER BUDGET?

#83 found the within-block spectrum knees at 5 (calibrated against known-rank-2 and rank-5 worlds),
while #72 found the cross-block spectrum has no knee through rank 32. #83c called that the sharpest
structural statement in the project and then named the reason it might be false:

  the within-block estimator draws its person scores from the target's own m = 10-24 columns.
  the cross-block estimator draws them from ~500 columns pooled over 31 other blocks.

A knee at 5 out of <=24 and no knee at 32 out of ~500 are not like-for-like. If the difference is the
BUDGET rather than the structure, #83c is withdrawn.

The match: build the cross-block scores from a RANDOM SUBSET of exactly m other-block columns -- the
same number the within-block estimator gets -- so both spectra are searched in spaces of equal size.
Three independent subsets per block, because one subset is one draw and #83c is not worth deciding
on a lucky one.

ESTIMAND        C corrected against its person-permutation null, swept over rank, with the score
                basis restricted to m randomly chosen other-block columns; and its knee.
IDENTIFICATION  identified. The only thing that changes from #72's design is the number of columns
                the basis is estimated from; ranks, projections and evaluation are unchanged.
SCOPE           the 23 blocks A09/R114 identified. Gate + demographics projected out (#77).
WORLDS          structural  budget-matched C still has no knee -> #83c stands, and the two levels
                            are genuinely different kinds of object
                budgetary   budget-matched C knees near 5 like W -> #83c is WITHDRAWN and the
                            difference between the levels was the size of the search space
KILL            threshold-free: the knee is read by the same per-dimension rule as #83, against the
                same known-rank calibration, and all three column subsets are published separately
                so a subset-driven answer is visible rather than averaged away.
POSITIVE CTRL   a shared rank-5 world at the SAME restricted budget: the estimator must still knee
                at 5 when only m columns are available, or the restriction itself destroys the
                instrument and no conclusion about the real curve is licensed.
NEGATIVE CTRL   person-permutation at every rank and every subset.
NOISE FLOOR     2 masks x 3 column subsets.
MULTIPLICITY    23 blocks x 3 subsets x 7 ranks x 2 arms x 2 seeds, published whole.
IMPOSSIBLE      matching the NUMBER OF PEOPLE as well -- both estimators already use the same people,
                so the budgets differ only in columns. Stated so the match is not over-claimed.
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd, lstsq
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
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
BLKS=sorted(RAW)
E=np.zeros((len(ALLP),len(BLKS)))
for k,q in enumerate(BLKS): E[[PM[p] for p in RAW[q]['ppl']],k]=1.
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
COV=pd.DataFrame({'male':pd.to_numeric(df.get('biomale'),errors='coerce'),'age':df['age'].map(AGEMAP)}).reindex(ALLP)
for c in ['opennessvariable','neuroticismvariable','extroversionvariable',
          'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']:
    if c in df.columns: COV[c]=pd.to_numeric(df[c],errors='coerce').reindex(ALLP).values
COV=COV.fillna(COV.median()).values; COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
MASK=0.15; SEEDS=[11,29]; SUBSETS=[0,1,2]; KCS=[1,2,3,5,8,12,16]
print(f"targets {len(IDENT)}   option counts {min(RAW[q]['M'].shape[1] for q in IDENT)}-"
      f"{max(RAW[q]['M'].shape[1] for q in IDENT)}",flush=True)

def pooled(target,BLK):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=BLK[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in BLK[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0)
    return np.where(np.isnan(Z),mu,Z)

def scores_budget(Zp,target,m,sub,K):
    rng=np.random.default_rng(3300+sub)
    pick=rng.choice(Zp.shape[1],size=min(m,Zp.shape[1]),replace=False)
    Z=Zp[:,pick]; Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:K]*S[:K]
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b

def c_skill(M,U,rows,Kc,seed,permute=False):
    n,m=M.shape
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    Us=U[rows][:,:Kc]
    if permute: Us=Us[np.random.default_rng(seed+7).permutation(n)]
    Us=(Us-Us.mean(0))/(Us.std(0)+1e-12)
    C=np.zeros_like(M)
    for j in range(m):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),Us[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(n),Us]@b
    base=np.mean((M[he]-gm)**2)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    IB=np.broadcast_to(I,M.shape)
    return f(IB,P,C)-f(IB,P)

def synth_shared5(seed):
    rng=np.random.default_rng(7700+seed); F=rng.normal(size=(len(ALLP),5)); out={}
    for q in BLKS:
        M=RAW[q]['M']; n,m=M.shape; rows=[PM[p] for p in RAW[q]['ppl']]
        L=rng.normal(size=(5,m))*0.30
        p=np.clip(M.mean(0)[None,:]+F[rows]@L,0.02,0.98)
        out[q]=dict(M=(rng.random((n,m))<p).astype(float),ppl=RAW[q]['ppl'])
    return out

rows=[]
def sweep(BLK,targets,tag):
    for i,t in enumerate(targets):
        M=BLK[t]['M']; m=M.shape[1]; r_=np.array([PM[p] for p in BLK[t]['ppl']])
        Zp=pooled(t,BLK)
        for sub in SUBSETS:
            U=scores_budget(Zp,t,m,sub,max(KCS))
            for Kc in [k for k in KCS if k<=m]:
                for sd in SEEDS:
                    rows.append(dict(world=tag,q=t,sub=sub,Kc=Kc,seed=sd,arm='real',
                                     C=c_skill(M,U,r_,Kc,sd),m=m))
                    rows.append(dict(world=tag,q=t,sub=sub,Kc=Kc,seed=sd,arm='perm',
                                     C=c_skill(M,U,r_,Kc,sd,True),m=m))
        print(f"  [{tag}] {i+1}/{len(targets)}",flush=True)
sweep(RAW,IDENT,'real')
sweep(synth_shared5(1),IDENT[:8],'shared5')
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R127_matched_dimensional_budget/results/'
D.to_csv(OUT+'grid.csv',index=False)

def corr(w,sub=None):
    d=D[D.world==w]
    if sub is not None: d=d[d['sub']==sub]
    return d[d.arm=='real'].groupby('Kc').C.mean()-d[d.arm=='perm'].groupby('Kc').C.mean()
def perdim(c):
    k=np.array(c.index,dtype=float); return np.diff(c.values)/np.diff(k)
def knee(c,frac=0.10):
    k=np.array(c.index); inc=perdim(c)
    if inc.max()<=0: return int(k[0])
    for i,g_ in enumerate(inc):
        if g_<frac*inc.max(): return int(k[i])
    return int(k[-1])

print("\n=== BUDGET-MATCHED CROSS-BLOCK CURVE (basis from m random other-block columns) ===")
T=pd.DataFrame({f"subset {s}":corr('real',s) for s in SUBSETS})
T['mean']=T.mean(axis=1); T['ctrl_shared5']=corr('shared5')
print(T.round(4).to_string())
print("\n=== GAIN PER DIMENSION, and the knee, per subset (published separately) ===")
for s in SUBSETS:
    c=corr('real',s); print(f"  subset {s}: knee={knee(c):2d}   {np.round(perdim(c),5).tolist()}")
cm=corr('real'); c5=corr('shared5')
print(f"  MEAN     : knee={knee(cm):2d}   {np.round(perdim(cm),5).tolist()}")
print(f"  ctrl r=5 : knee={knee(c5):2d}   {np.round(perdim(c5),5).tolist()}")

print("\n  CONDITIONAL KILL -- gate first")
g=knee(c5)>=4 and knee(c5)<=8
print(f"   (a) the restriction preserves the instrument (shared rank-5 world still knees near 5): "
      f"{'PASS' if g else 'FAIL -- restricting the budget destroyed the estimator'} (knee {knee(c5)})")
if not g:
    print("   -> UNVERIFIED: no conclusion about the real curve is licensed at this budget.")
else:
    kk=[knee(corr('real',s)) for s in SUBSETS]
    print(f"\n   budget-matched cross-block knee, per subset: {kk}   mean-curve knee: {knee(cm)}")
    print(f"   #83 within-block knee: 5   |   #72 full-budget cross-block: no knee through 32")
    if knee(cm)>=max(KCS) or all(k>=8 for k in kk):
        print("\n   -> #83c STANDS. Even at an equal dimensional budget the general spectrum does not")
        print("      knee where the specific one does. The two levels are different kinds of object,")
        print("      and the difference is not the size of the search space.")
    elif abs(knee(cm)-5)<=1:
        print("\n   -> #83c WITHDRAWN. At equal budget the general spectrum knees where the specific")
        print("      one does. The difference #83c called structural was the search-space size.")
    else:
        print(f"\n   -> PARTIAL: knee at {knee(cm)} against the specific side's 5. The levels differ,")
        print("      but by less than #83c claimed, and the claim must be restated as a comparison")
        print("      of knees rather than presence-versus-absence of one.")
print("\nN/A, with what it would require: matching the NUMBER OF PEOPLE as well is unnecessary here -- "
      "both estimators already use the same people, so the budgets differ only in columns.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
