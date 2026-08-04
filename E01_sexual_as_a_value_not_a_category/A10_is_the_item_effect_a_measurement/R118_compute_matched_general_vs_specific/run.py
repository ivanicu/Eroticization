import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R03 -- DOMAIN-GENERAL vs DOMAIN-SPECIFIC, PAID FOR IN PARAMETERS.

R02 measured C (cross-block, person-side) = +0.0172 against W (within-block) = +0.0078 at Kw=1, and
refused the ordering because C used 4 factors and W used 1. That refusal was right and the fix is
not a rank match -- it is a PARAMETER COUNT, and the two are not the same thing here:

  C at rank k  fits k x m loadings on the target block. The PERSON SCORES ARE FREE: they come
               entirely from the other 31 blocks and no target cell is used to estimate them.
  W at rank k  fits k x (n + m) parameters on the target block -- both scores and loadings.

With n in the thousands and m ~ 20, W is one to three ORDERS OF MAGNITUDE more expensive per unit
of rank. So a rank match would have handed W a large hidden advantage, and reporting it as "matched"
would have been an unavailability claim in the flattering direction.

ESTIMAND        C and W as Shapley held-out R2, each reported against its own degrees of freedom
                fit to the TARGET block's training cells, swept over both ranks.
IDENTIFICATION  identified. df is exact, not estimated: df_C = Kc*m, df_W = Kw*(n+m).
SCOPE           the 23 blocks A09/R114 identified.
WORLDS          general-wins   C exceeds W at every df where both are measured
                specific-wins  W exceeds C at matched df
                mixed          the ordering flips somewhere in the grid -> report where
KILL            threshold-free: the ordering is declared per (Kc,Kw) cell only where the gap
                exceeds 2x the cell's own seed spread, and the WHOLE grid is published including
                the cells that disagree. No cell is reported alone.
POSITIVE CTRL   W must be positive somewhere in the sweep -- it is the object of 105 prior rounds.
                If W is nowhere positive the instrument, not the world, is the finding.
NEGATIVE CTRL   person-permutation on the cross-block scores, run at every Kc.
PLACEBO         inherited (R02): I and P unmoved by the permutation.
NOISE FLOOR     2 masks per cell; spread reported per cell.
MULTIPLICITY    23 blocks x 4 Kc x 4 Kw x 2 seeds, published whole.
IMPOSSIBLE      a cross-block ITEM channel (per-option content annotations) -- unchanged from R02.
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools
from numpy.linalg import svd, lstsq
from math import factorial
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
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW])); PM={p:i for i,p in enumerate(ALLP)}
print(f"identified targets {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29]; KCS=[1,2,4,8]; KWS=[0,1,2,4]

def other_scores(target,K=8):
    cols=[]
    for q in RAW:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,Vt=svd(Z,full_matrices=False)
    return U[:,:K]*S[:K]

def decompose(target,U_all,Kc,Kw,seed,permute=False):
    M=RAW[target]['M']; n,m=M.shape; rows=[PM[p] for p in RAW[target]['ppl']]
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
    if Kw>0:
        F=np.where(np.isnan(Rres),0.,Rres)
        for _ in range(15):
            Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw])
        Uu,Ss,Vv=svd(F,full_matrices=False); W=(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw]
    else: W=np.zeros_like(M)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':W}
    base=np.mean((M[he]-gm)**2); names='IPCW'; v={}
    for bits in range(16):
        S=frozenset([c for j,c in enumerate(names) if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-p[he])**2)/base
    out={'full':v[frozenset(names)],'df_C':Kc*m,'df_W':Kw*(n+m)}
    for c in names:
        o=[x for x in names if x!=c]; tot=0.
        for r in range(4):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(3-len(S))/24.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

rows=[]
for i,t in enumerate(IDENT):
    U_all=other_scores(t)
    n,m=RAW[t]['M'].shape
    for Kc in KCS:
        for sd in SEEDS:
            rows.append(dict(q=t,Kc=Kc,Kw=0,seed=sd,arm='perm',n=n,m=m,
                             **decompose(t,U_all,Kc,0,sd,permute=True)))
        for Kw in KWS:
            for sd in SEEDS:
                rows.append(dict(q=t,Kc=Kc,Kw=Kw,seed=sd,arm='real',n=n,m=m,
                                 **decompose(t,U_all,Kc,Kw,sd)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R118_compute_matched_general_vs_specific/results/'
D.to_csv(OUT+'grid.csv',index=False)

R=D[D.arm=='real']; N=D[D.arm=='perm']
print("\n=== THE PARAMETER PRICE, per block median ===")
pr=R.groupby(['Kc','Kw'])[['df_C','df_W']].median()
pr['W_costs_x']=(pr.df_W/pr.df_C).round(1)
print(pr.astype(int,errors='ignore').to_string())

print("\n=== C and W across the whole grid (mean over 23 blocks x 2 seeds) ===")
print(R.groupby(['Kc','Kw'])[['C','W','full']].mean().round(4).unstack('Kw').to_string())

print("\n=== C corrected against the person-permutation null, by Kc ===")
cc=R[R.Kw==0].groupby('Kc').C.mean()-N.groupby('Kc').C.mean()
print(pd.DataFrame({'C_real':R[R.Kw==0].groupby('Kc').C.mean(),
                    'C_perm':N.groupby('Kc').C.mean(),'C_corrected':cc}).round(4).to_string())

print("\n=== THE ORDERING, per cell, declared only above 2x the cell's own seed spread ===")
res=[]
for Kc in KCS:
    for Kw in KWS:
        if Kw==0: continue
        d=R[(R.Kc==Kc)&(R.Kw==Kw)]
        pb=d.groupby('q')[['C','W','df_C','df_W']].mean()
        sp=np.sqrt(d.groupby('q').C.std()**2+d.groupby('q').W.std()**2)
        gap=pb.C-pb.W; ok=gap.abs()>2*sp
        res.append(dict(Kc=Kc,Kw=Kw,C=pb.C.median(),W=pb.W.median(),gap=gap.median(),
                        spread2=2*sp.median(),
                        n_general=int(((gap>0)&ok).sum()),n_specific=int(((gap<0)&ok).sum()),
                        n_tied=int((~ok).sum()),
                        df_ratio=float((pb.df_W/pb.df_C).median())))
G=pd.DataFrame(res)
print(G.round(4).to_string(index=False))

wpos=R[R.Kw>0].groupby(['Kc','Kw']).W.mean().max()
print(f"\n  CONDITIONAL KILL -- gate first")
print(f"   (a) W positive somewhere in the sweep : {'PASS' if wpos>0 else 'FAIL'} (max mean W {wpos:+.4f})")
print(f"   (b) C above the permutation null      : {'PASS' if (cc>0).all() else 'FAIL'} "
      f"(corrected C at every Kc: {cc.round(4).tolist()})")
tot=G.n_general.sum()+G.n_specific.sum()+G.n_tied.sum()
print(f"\n   ACROSS THE WHOLE GRID ({len(G)} cells x 23 blocks = {tot} comparisons):")
print(f"     domain-GENERAL larger : {G.n_general.sum()}  ({G.n_general.sum()/tot:.1%})")
print(f"     domain-SPECIFIC larger: {G.n_specific.sum()}  ({G.n_specific.sum()/tot:.1%})")
print(f"     not distinguishable   : {G.n_tied.sum()}  ({G.n_tied.sum()/tot:.1%})")
print(f"   and W pays {G.df_ratio.median():.0f}x more parameters than C to do it.")
if G.n_general.sum()>G.n_specific.sum()*2:
    print("\n   -> DOMAIN-GENERAL WINS, and it wins while spending far fewer target-block parameters.")
    print("      The person-side structure that transfers across content domains is larger than the")
    print("      structure specific to a domain -- so the readout is not assembled per-domain.")
elif G.n_specific.sum()>G.n_general.sum()*2:
    print("\n   -> DOMAIN-SPECIFIC WINS. the readout is assembled per-domain and R02's ordering was rank.")
else:
    print("\n   -> MIXED: the ordering flips inside the grid. Both exist; neither dominates, and the")
    print("      cell that is reported alone is the finding, not the effect.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
