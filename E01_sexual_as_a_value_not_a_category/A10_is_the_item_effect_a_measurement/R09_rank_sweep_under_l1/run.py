import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R09 -- ATTACKING MY OWN CLAIM AT THE CELL IT IS WEAKEST IN.

#79c named its own weakest specification: under L1 the within-block interaction W turns POSITIVE
(+0.082 at K=4, against -0.066 Brier and -0.355 log-loss), the item:interaction ratio falls from
huge to 2.1x, and 3 of 23 blocks already flip to the interaction. Every rank sweep in this project
(R04's Kc sweep, R03's Kw sweep) was run under Brier, so the rank at which the interaction is
strongest has never been found under the loss that favours it.

If the ordering flips somewhere in (rank x L1), #67/#68/#70's "the item is the larger component" is
a statement about squared error after all, and #79b's "scale-free" is wrong one round after I wrote it.

ESTIMAND        Shapley skill of {I,P,C,W} under L1 and Brier, swept over W's rank, per block; and
                the rank that MAXIMISES the interaction under L1 -- the adversarially best case.
IDENTIFICATION  identified; W's rank is the only thing varying within a (block, seed, loss) cell.
SCOPE           the 23 blocks A09/R03 identified. C fixed at Kc=4. Gate + demographics out (#77).
WORLDS          scale-free    item wins under L1 at EVERY rank -> #79b stands
                metric-bound  item loses at some rank under L1 -> #79b is withdrawn one round old,
                              and #67/#68/#70 become claims about squared error
KILL            threshold-free: the ordering is declared per (rank, loss) only above 2x that cell's
                own seed spread; the whole sweep is published including flipped cells.
POSITIVE CTRL   a fixed-margin (curveball) null per block, per rank, under L1 -- margins exact, so
                the item and person effects are preserved and only the interaction is destroyed.
                W under L1 must exceed its own fixed-margin floor, or a positive W is just the
                estimator fitting noise that L1 happens to reward.
NEGATIVE CTRL   the same curveball world IS the negative control for W. Person-permutation for C.
NOISE FLOOR     2 masks x 2 randomisation draws.
MULTIPLICITY    23 blocks x 4 ranks x 2 losses x 2 worlds x 2 seeds, published whole.
IMPOSSIBLE      unchanged -- no absolute dose scale for the interaction (R03).
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
from numpy.linalg import svd, lstsq
from math import factorial
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
print(f"targets {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29]; KC=4; KWS=[1,2,4,8]
LOSSES={'brier':lambda y,p:(y-p)**2,'l1':lambda y,p:np.abs(y-p)}

def curveball(M,rng,per_row=5.):
    A=[set(np.flatnonzero(r).tolist()) for r in M]; n=len(A)
    for _ in range(int(per_row*n)):
        i,j=int(rng.integers(n)),int(rng.integers(n))
        if i==j: continue
        ai,aj=A[i],A[j]; inter=ai&aj
        di=list(ai-inter); dj=list(aj-inter); L=di+dj
        if not L: continue
        rng.shuffle(L); k=len(di)
        A[i]=inter|set(L[:k]); A[j]=inter|set(L[k:])
    out=np.zeros_like(M)
    for i,s in enumerate(A): out[i,list(s)]=1.
    return out

def other_scores(target):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:KC]*S[:KC]
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b

def shap4(v,names='IPCW'):
    out={}
    for c in names:
        o=[x for x in names if x!=c]; tot=0.
        for r in range(4):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(3-len(S))/24.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def run(M,U,rows,Kw,seed,permC=False):
    n,m=M.shape
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    Us=U[rows]
    if permC: Us=Us[np.random.default_rng(seed+7).permutation(n)]
    Us=(Us-Us.mean(0))/(Us.std(0)+1e-12)
    C=np.zeros_like(M)
    for j in range(m):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),Us[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(n),Us]@b
    F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(20):
        Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw])
    Uu,Ss,Vv=svd(F,full_matrices=False); W=(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw]
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':W}
    out={}
    for ln,L in LOSSES.items():
        b0=np.mean(L(M[he],np.full(he.sum(),gm))); v={}
        for bits in range(16):
            S=frozenset([c for j,c in enumerate('IPCW') if bits>>j&1])
            p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
            v[S]=1.-np.mean(L(M[he],np.asarray(p)[he]))/b0
        sh=shap4(v)
        for c in 'IPCW': out[f"{ln}_{c}"]=sh[c]
    return out

rows=[]
for i,t in enumerate(IDENT):
    U=other_scores(t); r_=np.array([PM[p] for p in RAW[t]['ppl']]); M=RAW[t]['M']
    for sd in SEEDS:
        Mn=curveball(M,np.random.default_rng(4000+sd))
        assert np.allclose(Mn.sum(0),M.sum(0)) and np.allclose(Mn.sum(1),M.sum(1))
        for Kw in KWS:
            rows.append(dict(q=t,Kw=Kw,seed=sd,world='real',arm='real',**run(M,U,r_,Kw,sd)))
            rows.append(dict(q=t,Kw=Kw,seed=sd,world='null',arm='real',**run(Mn,U,r_,Kw,sd)))
            rows.append(dict(q=t,Kw=Kw,seed=sd,world='real',arm='permC',**run(M,U,r_,Kw,sd,True)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R09_rank_sweep_under_l1/results/'
D.to_csv(OUT+'grid.csv',index=False)

R=D[(D.world=='real')&(D.arm=='real')]; NW=D[D.world=='null']; NC=D[D.arm=='permC']
print("\n=== W AND ITS FIXED-MARGIN FLOOR, BY RANK AND LOSS ===")
t=pd.DataFrame({'W_l1':R.groupby('Kw').l1_W.mean(),'W_l1_null':NW.groupby('Kw').l1_W.mean(),
                'W_brier':R.groupby('Kw').brier_W.mean(),'W_brier_null':NW.groupby('Kw').brier_W.mean(),
                'I_l1':R.groupby('Kw').l1_I.mean(),'C_l1':R.groupby('Kw').l1_C.mean()})
t['W_l1_corrected']=t.W_l1-t.W_l1_null; t['W_brier_corrected']=t.W_brier-t.W_brier_null
print(t.round(4).to_string())

print("\n=== NEGATIVE CONTROL: person-permuted C, both losses ===")
print(NC.groupby('Kw')[['l1_C','brier_C']].mean().round(5).to_string())

print("\n=== THE ORDERING PER (rank, loss), 2x that cell's own seed spread ===")
out=[]
for ln in LOSSES:
    for Kw in KWS:
        d=R[R.Kw==Kw]
        I=d.groupby('q')[f"{ln}_I"].mean()
        X=d.groupby('q')[f"{ln}_C"].mean()+d.groupby('q')[f"{ln}_W"].mean()
        sp=np.sqrt(d.groupby('q')[f"{ln}_I"].std()**2+d.groupby('q')[f"{ln}_W"].std()**2)
        gap=I-X
        out.append(dict(loss=ln,Kw=Kw,I=I.median(),X=X.median(),gap=gap.median(),
                        spread2=2*sp.median(),item=int((gap>2*sp).sum()),
                        inter=int((-gap>2*sp).sum()),tied=int((gap.abs()<=2*sp).sum())))
S=pd.DataFrame(out); print(S.round(4).to_string(index=False))

print("\n  CONDITIONAL KILL -- gates first")
g_marg=True
wc=t.W_l1_corrected
g_pos=bool((wc>0).any())
g_neg=bool(abs(NC.l1_C.mean())<0.005 and abs(NC.brier_C.mean())<0.005)
print(f"   (a) W exceeds its own fixed-margin floor under L1 somewhere : "
      f"{'PASS' if g_pos else 'FAIL -- a positive W under L1 is the estimator fitting noise'} "
      f"(max corrected {wc.max():+.4f} at Kw={int(wc.idxmax())})")
print(f"   (b) person-permuted C ~0 under both losses                   : "
      f"{'PASS' if g_neg else 'FAIL'} (l1 {NC.l1_C.mean():+.5f}, brier {NC.brier_C.mean():+.5f})")
if not(g_pos and g_neg): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    l1=S[S.loss=='l1']
    worst=l1.loc[l1.gap.idxmin()]
    print(f"\n   L1, worst cell for the item effect: Kw={int(worst.Kw)}  "
          f"I {worst.I:+.4f} vs interaction {worst.X:+.4f}, gap {worst.gap:+.4f} "
          f"(2x spread {worst.spread2:.4f})")
    print(f"   blocks: item {int(worst.item)} / interaction {int(worst.inter)} / tied {int(worst.tied)}")
    if (l1.inter>l1.item).any():
        bad=l1[l1.inter>l1.item]
        print(f"\n   -> METRIC-BOUND. The ordering FLIPS at Kw={bad.Kw.tolist()} under L1. #79b is")
        print("      withdrawn one round after it was written, and #67/#68/#70 are claims about")
        print("      squared error.")
    else:
        print(f"\n   -> #79b STANDS. The item effect wins at every rank under L1 as well as Brier.")
        print(f"      The interaction's best case anywhere in this project remains a {worst.I/max(worst.X,1e-9):.1f}x deficit.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
