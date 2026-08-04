import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R02 -- THE DECOMPOSITION AT THE LEVEL WHERE THE MODELS ACTUALLY DIFFER.

#69 showed the item main effect is a DERIVATION from prevalence dispersion, and that within a block
"content" and "base rate" are THE SAME NUMBER -- so A09's contrast cannot adjudicate model A vs
model B. It measured a real magnitude and answered a different question.

The level where they differ is ACROSS blocks. Base rates are block-local: an option in block 1 and
an option in block 2 are different options, so there is no item-level quantity that can transfer.
A person-side readout weight CAN transfer. So:

  if a person's deviation from prevalence in one block predicts their deviation in another,
  that shared part CANNOT be a base rate, and it is the only thing in this release that
  distinguishes "value assigned by a person" from "property of the content".

This round measures all four on ONE scale, on the SAME held-out cells:

  I  item prevalence          (block-local, a derivation -- #69)
  P  person breadth           (block-local scalar)
  C  CROSS-BLOCK interaction  person factors fit ONLY on the OTHER 31 blocks, never on the target
  W  WITHIN-BLOCK interaction low-rank structure in the target's own residual

ESTIMAND        Shapley held-out R2 of {I,P,C,W} over all 24 orderings, per block.
IDENTIFICATION  C is identified because the person scores come entirely from other blocks and the
                loadings from the target's TRAINING cells -- no held-out cell touches either.
SCOPE           the 23 blocks A09/R114 identified. Population/instrument/regime as A09.
WORLDS          A content property dominates -> C ~ 0: nothing person-side survives a block change
                B individualised readout     -> C > 0 and resolvable
                C' mixed                     -> C > 0 but small next to W
KILL            C is real iff it exceeds the person-permutation null by > 2x seed spread (#34).
                THEN, and only then, the ordering C vs W is reported with its own spread.
                No threshold is chosen; every comparison is against a measured floor.
POSITIVE CTRL   the WITHIN-block component W is the positive control for the machinery: it is the
                quantity 105 rounds have measured and it must come back positive.
NEGATIVE CTRL   person-permutation on the cross-block scores -- destroys the person correspondence
                between the other blocks and the target, preserves everything else exactly.
                World it excludes: "C is an artefact of both sides sharing a rank-K basis".
PLACEBO         I and P recomputed under the permutation must be UNCHANGED (they do not use C).
NOISE FLOOR     3 masks x 3 permutation draws.
MULTIPLICITY    23 blocks x 2 K_within x 3 seeds x 16 subset models, reported whole.
IMPOSSIBLE      a cross-block ITEM property -- it would require item-level annotations (content
                codes for each option). This release carries option text only, and #28/#44 measured
                that string proxies here lose to their own shams. N/A, with that requirement.
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
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW]))
PM={p:i for i,p in enumerate(ALLP)}
print(f"blocks {len(RAW)}   identified targets {len(IDENT)}   people {len(ALLP)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]; KC=4

def other_scores(target,K=KC):
    """person factors from every block EXCEPT the target. no target cell is ever seen."""
    cols=[]
    for q in RAW:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan)
        Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,Vt=svd(Z,full_matrices=False)
    return U[:,:K]*S[:K]

def decompose(target,U_all,Kw,seed,permute=False):
    M=RAW[target]['M']; rows=[PM[p] for p in RAW[target]['ppl']]
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    U=U_all[rows]
    if permute: U=U[np.random.default_rng(seed+555).permutation(len(rows))]
    U=(U-U.mean(0))/(U.std(0)+1e-12)
    # loadings fit on TRAINING cells only, one column at a time (each column has its own mask)
    C=np.zeros_like(M)
    for j in range(M.shape[1]):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),U[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(len(U)),U]@b
    if Kw>0:
        F=np.where(np.isnan(Rres),0.,Rres)
        for _ in range(25):
            Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw])
        Uu,Ss,Vv=svd(F,full_matrices=False); W=(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw]
    else: W=np.zeros_like(M)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':W}
    base=np.mean((M[he]-gm)**2); names='IPCW'; v={}
    for bits in range(16):
        S=frozenset([c for j,c in enumerate(names) if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-p[he])**2)/base
    out={'full':v[frozenset(names)]}
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
    for Kw in [0,1]:
        for sd in SEEDS:
            rows.append(dict(q=t,Kw=Kw,seed=sd,arm='real',n=RAW[t]['M'].shape[0],
                             m=RAW[t]['M'].shape[1],**decompose(t,U_all,Kw,sd)))
            rows.append(dict(q=t,Kw=Kw,seed=sd,arm='perm',n=RAW[t]['M'].shape[0],
                             m=RAW[t]['M'].shape[1],**decompose(t,U_all,Kw,sd,permute=True)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R117_the_four_way_decomposition/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== FOUR-WAY SHAPLEY DECOMPOSITION (mean over 23 blocks x 3 seeds) ===")
print(D.groupby(['Kw','arm'])[['full','I','P','C','W']].mean().round(4).to_string())

print("\n=== PLACEBO: I and P must be unchanged by the person permutation (they do not use C) ===")
for Kw in [0,1]:
    d=D[D.Kw==Kw]
    dI=abs(d[d.arm=='real'].groupby('q').I.mean()-d[d.arm=='perm'].groupby('q').I.mean()).median()
    dP=abs(d[d.arm=='real'].groupby('q').P.mean()-d[d.arm=='perm'].groupby('q').P.mean()).median()
    print(f"  Kw={Kw}  median |dI| {dI:.5f}   median |dP| {dP:.5f}   -> {'PASS' if max(dI,dP)<0.005 else 'FAIL'}")

for Kw in [0,1]:
    d=D[D.Kw==Kw]
    r=d[d.arm=='real'].groupby('q')[['I','P','C','W','full']].mean()
    pm=d[d.arm=='perm'].groupby('q').C.mean()
    sp=np.sqrt(d[d.arm=='real'].groupby('q').C.std()**2+d[d.arm=='perm'].groupby('q').C.std()**2)
    r['C_perm']=pm; r['C_c']=r.C-pm; r['spread']=sp; r['real']=r.C_c>2*sp
    print(f"\n{'='*78}\nK_within={Kw}")
    print(r.round(4).sort_values('C_c').to_string())
    print(f"\n  CONDITIONAL KILL -- gate first")
    print(f"   (a) W > 0 (the machinery reproduces what 105 rounds measured): "
          f"{'PASS' if (Kw==0 or r.W.median()>0) else 'FAIL'}  median W {r.W.median():+.4f}")
    print(f"   (b) C above the person-permutation null by >2x spread: "
          f"{int(r.real.sum())}/{len(r)} blocks   median C_c {r.C_c.median():+.4f}"
          f"  median 2x spread {2*r.spread.median():.4f}")
    if r.real.sum()<len(r)*0.5:
        print("   -> C UNVERIFIED at this K_within."); continue
    print(f"\n   ON ONE SCALE, SAME HELD-OUT CELLS:")
    print(f"     I (block-local, a DERIVATION)   {r.I.median():+.4f}")
    print(f"     P (block-local breadth)         {r.P.median():+.4f}")
    print(f"     C (CROSS-BLOCK, person-side)    {r.C.median():+.4f}   corrected {r.C_c.median():+.4f}")
    print(f"     W (within-block interaction)    {r.W.median():+.4f}")
    print(f"\n   At the CROSS-BLOCK level -- the only level where models A and B differ -- the item")
    print(f"   contribution is 0 BY CONSTRUCTION (no option is shared between blocks), and the")
    print(f"   measured person-side transfer is {r.C_c.median():+.4f}, resolvable in "
          f"{int(r.real.sum())}/{len(r)} blocks.")
print("\nN/A, with what it would require: a cross-block ITEM property needs content annotations per "
      "option. This release carries option text only, and #28/#44 measured that string proxies here "
      "lose to their own shams.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
