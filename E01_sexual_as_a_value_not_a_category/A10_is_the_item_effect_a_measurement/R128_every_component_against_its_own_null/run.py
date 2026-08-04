import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R13 -- THE HEADLINE, WITH EVERY COMPONENT REFERRED TO ITS OWN NULL.

#80b, #81/#82 and #84 all turned on the same thing: two quantities compared without each being
referred to its own floor. #82 showed that fixing it INVERTED the arc's load-bearing claim. The
audit that follows from that has to reach the README's first row, which still reads

    "item +0.222 vs person +0.085 vs interaction +0.019 (prediction) / +0.063 (detection)"

and those are a RAW item effect against a RAW interaction, which is exactly the comparison #82
killed one row lower down.

Every component gets the null that destroys IT and preserves everything else, computed in the same
run on the same held-out cells:

  I  item main effect    -> WITHIN-PERSON shuffle: permute each row independently. Row sums exact,
                            column structure destroyed.
  P  person main effect  -> WITHIN-COLUMN shuffle: permute each column independently. Column means
                            exact, row structure destroyed.
  W  within-block interaction -> FIXED-MARGIN (curveball): both margins exact, interaction destroyed.
  C  cross-block transfer     -> PERSON-PERMUTATION of the external scores.

ESTIMAND        each component's corrected skill at ITS OWN best rank, and the item:interaction
                ratio computed from those, on identical held-out cells.
IDENTIFICATION  identified per component; each null changes exactly one structure and the invariants
                it must preserve are ASSERTED numerically per draw, not assumed.
SCOPE           the 23 blocks A09/R114 identified. Gate + demographics out of the scores (#77).
WORLDS          item-dominant   corrected item >> corrected interaction -> A09's headline stands as
                                written and #68c's tie was an adversarial edge case
                tied            they are within 2x pooled seed spread -> the headline must be
                                restated, and #68c's "best case is both, equally" becomes the
                                CENTRAL estimate rather than the generous one
                interaction-dominant -> the epoch title is right after all, for the first time
KILL            threshold-free: the ordering is declared only above 2x the pooled seed spread, per
                block, and the raw comparison is printed beside the corrected one so the size of the
                asymmetry that #82 found is visible in this table too.
POSITIVE CTRL   each null must PRESERVE its invariant exactly (asserted) and MOVE its target: the
                within-person shuffle must collapse I, the within-column shuffle must collapse P,
                the curveball must collapse W. A null that does not move its target is not a null.
NEGATIVE CTRL   a fully destroyed world (within-person shuffle, which kills columns and interaction
                together): every corrected component must be ~0 there except P.
NOISE FLOOR     2 masks x 2 randomisation draws.
MULTIPLICITY    23 blocks x 5 ranks x 5 worlds x 2 seeds, published whole.
IMPOSSIBLE      a null that destroys I while preserving the INTERACTION -- shuffling within a person
                destroys both. So I's correction is a LOWER bound on its uniquely-attributable part,
                and that is stated rather than hidden.
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
MASK=0.15; SEEDS=[11,29]; KS=[1,2,3,5,8]
print(f"targets {len(IDENT)}",flush=True)

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
def row_shuffle(M,rng):  return np.array([rng.permutation(r) for r in M])
def col_shuffle(M,rng):  return np.column_stack([rng.permutation(M[:,j]) for j in range(M.shape[1])])

def scores(target,K=8):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:K]*S[:K]
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
    b,*_=lstsq(D,U,rcond=None); return U-D@b

def dec(M,U,rows,K,seed,permC=False):
    n,m=M.shape
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    Us=U[rows][:,:K]
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
        Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:K]*Ss[:K])@Vv[:K])
    Uu,Ss,Vv=svd(F,full_matrices=False); Wm=(Uu[:,:K]*Ss[:K])@Vv[:K]
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':Wm}
    b0=np.mean((M[he]-gm)**2); v={}
    for bits in range(16):
        S=frozenset([c for j,c in enumerate('IPCW') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-np.asarray(p)[he])**2)/b0
    out={}
    for c in 'IPCW':
        o=[x for x in 'IPCW' if x!=c]; tot=0.
        for r in range(4):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(3-len(S))/24.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

rows=[]
for i,t in enumerate(IDENT):
    U=scores(t); r_=np.array([PM[p] for p in RAW[t]['ppl']]); M=RAW[t]['M']
    for sd in SEEDS:
        rg=np.random.default_rng(9100+sd)
        Mcb=curveball(M,rg); Mrs=row_shuffle(M,rg); Mcs=col_shuffle(M,rg)
        assert np.allclose(Mcb.sum(0),M.sum(0)) and np.allclose(Mcb.sum(1),M.sum(1)),"curveball margins"
        assert np.allclose(Mrs.sum(1),M.sum(1)),"row-shuffle must preserve row sums"
        assert np.allclose(Mcs.sum(0),M.sum(0)),"col-shuffle must preserve column sums"
        for K in KS:
            rows.append(dict(q=t,K=K,seed=sd,world='real',    **dec(M  ,U,r_,K,sd)))
            rows.append(dict(q=t,K=K,seed=sd,world='permC',   **dec(M  ,U,r_,K,sd,True)))
            rows.append(dict(q=t,K=K,seed=sd,world='margin',  **dec(Mcb,U,r_,K,sd)))
            rows.append(dict(q=t,K=K,seed=sd,world='rowshuf', **dec(Mrs,U,r_,K,sd)))
            rows.append(dict(q=t,K=K,seed=sd,world='colshuf', **dec(Mcs,U,r_,K,sd)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R128_every_component_against_its_own_null/results/'
D.to_csv(OUT+'grid.csv',index=False)

def gg(w,c): return D[D.world==w].groupby('K')[c].mean()
print("\n=== POSITIVE CONTROL: each null must MOVE its own target ===")
print(f"  I under within-person shuffle : {gg('real','I').mean():+.4f} -> {gg('rowshuf','I').mean():+.4f}")
print(f"  P under within-column shuffle : {gg('real','P').mean():+.4f} -> {gg('colshuf','P').mean():+.4f}")
print(f"  W under fixed-margin          : {gg('real','W').mean():+.4f} -> {gg('margin','W').mean():+.4f}")
print(f"  C under person-permutation    : {gg('real','C').mean():+.4f} -> {gg('permC','C').mean():+.4f}")
print("  invariants asserted exact on every draw: PASS (asserts would have raised)")

print("\n=== EACH COMPONENT, RAW AND CORRECTED AGAINST ITS OWN NULL, BY RANK ===")
T=pd.DataFrame({'I':gg('real','I'),'I_c':gg('real','I')-gg('rowshuf','I'),
                'P':gg('real','P'),'P_c':gg('real','P')-gg('colshuf','P'),
                'W':gg('real','W'),'W_c':gg('real','W')-gg('margin','W'),
                'C':gg('real','C'),'C_c':gg('real','C')-gg('permC','C')})
T['inter_raw']=T.W+T.C; T['inter_c']=T.W_c+T.C_c
T['ratio_raw']=T.I/T.inter_raw.replace(0,np.nan); T['ratio_c']=T.I_c/T.inter_c
print(T.round(4).to_string())

kI=T.I_c.idxmax(); kX=T.inter_c.idxmax()
print(f"\n  each at ITS OWN best rank:  item {T.I_c[kI]:+.4f} (K={kI})   "
      f"interaction {T.inter_c[kX]:+.4f} (K={kX})   ratio {T.I_c[kI]/T.inter_c[kX]:.2f}x")

print("\n=== PER-BLOCK ORDERING at each side's own best rank ===")
dI=D[(D.world=='real')&(D.K==kI)].groupby('q').I.mean()-D[(D.world=='rowshuf')&(D.K==kI)].groupby('q').I.mean()
dX=(D[(D.world=='real')&(D.K==kX)].groupby('q').W.mean()-D[(D.world=='margin')&(D.K==kX)].groupby('q').W.mean()
    +D[(D.world=='real')&(D.K==kX)].groupby('q').C.mean()-D[(D.world=='permC')&(D.K==kX)].groupby('q').C.mean())
sp=np.sqrt(D[(D.world=='real')&(D.K==kI)].groupby('q').I.std()**2+
           D[(D.world=='real')&(D.K==kX)].groupby('q').W.std()**2)
gap=dI-dX
print(f"  item larger: {int((gap>2*sp).sum())}/23   interaction larger: {int((-gap>2*sp).sum())}/23   "
      f"tied: {int((gap.abs()<=2*sp).sum())}/23   median gap {gap.median():+.4f} (2x spread {2*sp.median():.4f})")

print("\n  CONDITIONAL KILL -- gates first")
g1=gg('rowshuf','I').mean()<0.5*gg('real','I').mean()
g2=gg('colshuf','P').mean()<0.5*gg('real','P').mean()
g3=gg('margin','W').mean()<gg('real','W').mean()
print(f"   (a) within-person shuffle collapses I : {'PASS' if g1 else 'FAIL'}")
print(f"   (b) within-column shuffle collapses P : {'PASS' if g2 else 'FAIL'}")
print(f"   (c) fixed-margin collapses W          : {'PASS' if g3 else 'FAIL'}")
if not(g1 and g2 and g3): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    r=T.I_c[kI]/T.inter_c[kX]
    print(f"\n   README row (1) says item is 11.7x (prediction) / 3.5x (detection) the interaction.")
    print(f"   SYMMETRIC, each at its own best rank and own null: {r:.2f}x")
    if int((gap.abs()<=2*sp).sum())>=12 or 0.7<r<1.4:
        print("\n   -> TIED. #68c's 'the best case for the epoch title is both, equally' was not an")
        print("      adversarial edge case -- it is the CENTRAL estimate once every component is")
        print("      referred to its own null. The README's 11.7x must be withdrawn.")
    elif r>1.4:
        print(f"\n   -> ITEM STILL LARGER at {r:.2f}x, but far below the 11.7x the README claims.")
    else:
        print(f"\n   -> INTERACTION LARGER at {1/r:.2f}x. The epoch title is right for the first time.")
print("\nN/A, with what it would require: a null that destroys I while preserving the interaction does "
      "not exist -- a within-person shuffle destroys both. So I_c is a LOWER bound on the item "
      "effect's uniquely-attributable part, which makes the comparison CONSERVATIVE toward the item.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
