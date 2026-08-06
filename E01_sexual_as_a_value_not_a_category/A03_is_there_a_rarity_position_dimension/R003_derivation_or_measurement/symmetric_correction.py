import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R10 -- THE SYMMETRIC RUN, DESIGNED TO OVERTURN #71 RATHER THAN DEFEND IT.

#71 concluded "domain-general beats domain-specific, 186/276 comparisons" from `gap = pb.C - pb.W`
with NEITHER side corrected (#81). #80 then measured that the two sides have wildly different nulls:

    C  cross-block   own null = person-permutation      about -0.002   -> nearly unbiased
    W  within-block  own null = fixed-margin curveball  -0.070..-0.287 -> severely under-reported

So #71's comparison handicapped the DOMAIN-SPECIFIC side by up to 0.29 of held-out skill. If the
symmetric comparison inverts, the load-bearing claim of arc A10 is wrong and the README row that
cites it must be withdrawn.

This round is built so that outcome CAN happen: the estimand is the corrected gap, the grid is the
same Kc x Kw grid #71 used, and both nulls are computed per block, per rank, per seed, in the same
run as the quantity they correct -- never carried over from another round.

ESTIMAND        (C - C_null) - (W - W_null), per block, per (Kc,Kw), on identical held-out cells.
IDENTIFICATION  identified. Each null is the same estimator applied to a world that destroys exactly
                one structure: person-permutation destroys the person<->block correspondence that C
                needs; fixed-margin randomisation destroys the interaction W needs while preserving
                both margins exactly (asserted per draw).
SCOPE           the 23 blocks A09/R114 identified. Gate + demographics projected out of the person
                scores throughout (#77), so C enters at its corrected magnitude, not its inflated one.
WORLDS          general-wins   corrected gap > 0 -> #71 survives its own correction
                specific-wins  corrected gap < 0 -> #71 INVERTS and is withdrawn
                mixed          the sign flips inside the grid -> neither dominates and #71 must be
                               restated as "both exist", which is a different claim from the one made
KILL            threshold-free: the corrected gap is declared per cell only above 2x that cell's own
                pooled seed spread; every cell is published including the ones that disagree, and the
                RAW gap is printed beside the corrected one so the size of the handicap is visible.
POSITIVE CTRL   each null must MOVE its own quantity: |C - C_null| and |W - W_null| both nonzero,
                and the fixed-margin draw must preserve margins exactly. A null that changes nothing
                is not a null.
NEGATIVE CTRL   a doubly-destroyed world -- fixed-margin randomised AND person-permuted -- in which
                both corrected quantities must fall to zero. This is the control that catches a
                correction which manufactures signal by subtracting a too-negative floor.
NOISE FLOOR     2 masks x 2 randomisation draws per cell.
MULTIPLICITY    23 blocks x 3 Kc x 3 Kw x 2 seeds x 4 worlds, published whole.
IMPOSSIBLE      a null that destroys W while preserving C exactly, or vice versa -- the two share the
                same residual matrix. Reported N/A: the corrections are therefore each valid for
                their own quantity and are NOT claimed to be independent of one another.
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
               'R114_fixed_margin_null/results/grid.csv')
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

MASK=0.15; SEEDS=[11,29]; KCS=[1,4,8]; KWS=[1,2,4]

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

def shap4(v):
    out={}
    for c in 'IPCW':
        o=[x for x in 'IPCW' if x!=c]; tot=0.
        for r in range(4):
            for S in itertools.combinations(o,r):
                tot+=factorial(len(S))*factorial(3-len(S))/24.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

def dec(M,U,rows,Kc,Kw,seed,permute=False):
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
    F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(20):
        Uu,Ss,Vv=svd(F,full_matrices=False); F=np.where(obs,Rres,(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw])
    Uu,Ss,Vv=svd(F,full_matrices=False); Wm=(Uu[:,:Kw]*Ss[:Kw])@Vv[:Kw]
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'C':C,'W':Wm}
    b0=np.mean((M[he]-gm)**2); v={}
    for bits in range(16):
        S=frozenset([c for j,c in enumerate('IPCW') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-np.asarray(p)[he])**2)/b0
    return shap4(v)

rows=[]
for i,t in enumerate(IDENT):
    U=scores(t); r_=np.array([PM[p] for p in RAW[t]['ppl']]); M=RAW[t]['M']
    for sd in SEEDS:
        Mn=curveball(M,np.random.default_rng(6000+sd))
        assert np.allclose(Mn.sum(0),M.sum(0)) and np.allclose(Mn.sum(1),M.sum(1)),"margins broken"
        for Kc in KCS:
            for Kw in KWS:
                rows.append(dict(q=t,Kc=Kc,Kw=Kw,seed=sd,world='real',   **dec(M ,U,r_,Kc,Kw,sd)))
                rows.append(dict(q=t,Kc=Kc,Kw=Kw,seed=sd,world='permC',  **dec(M ,U,r_,Kc,Kw,sd,True)))
                rows.append(dict(q=t,Kc=Kc,Kw=Kw,seed=sd,world='margin', **dec(Mn,U,r_,Kc,Kw,sd)))
                rows.append(dict(q=t,Kc=Kc,Kw=Kw,seed=sd,world='both',   **dec(Mn,U,r_,Kc,Kw,sd,True)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R125_symmetric_correction/results/'
D.to_csv(OUT+'grid.csv',index=False)

def g(w,col): return D[D.world==w].groupby(['Kc','Kw'])[col].mean()
print("\n=== POSITIVE CONTROL: each null must MOVE its own quantity ===")
print(f"  |C - C_permC| mean {abs(g('real','C')-g('permC','C')).mean():.5f}"
      f"   |W - W_margin| mean {abs(g('real','W')-g('margin','W')).mean():.5f}")
print("  margins asserted exact on every curveball draw: PASS (assert would have raised)")

print("\n=== NEGATIVE CONTROL: doubly-destroyed world, both corrected quantities must vanish ===")
dc=g('both','C')-g('permC','C'); dw=g('both','W')-g('margin','W')
print(f"  corrected C in the doubly-destroyed world: mean {dc.mean():+.5f}  max |{dc.abs().max():.5f}|")
print(f"  corrected W in the doubly-destroyed world: mean {dw.mean():+.5f}  max |{dw.abs().max():.5f}|")

print("\n=== RAW vs SYMMETRICALLY CORRECTED, over the Kc x Kw grid ===")
T=pd.DataFrame({'C_raw':g('real','C'),'C_null':g('permC','C'),
                'W_raw':g('real','W'),'W_null':g('margin','W')})
T['C_c']=T.C_raw-T.C_null; T['W_c']=T.W_raw-T.W_null
T['gap_raw']=T.C_raw-T.W_raw; T['gap_corrected']=T.C_c-T.W_c
print(T.round(4).to_string())

print("\n=== PER-CELL VERDICT, 2x that cell's own pooled seed spread ===")
out=[]
for (Kc,Kw) in T.index:
    d=lambda w,c: D[(D.world==w)&(D.Kc==Kc)&(D.Kw==Kw)].groupby('q')[c]
    cc=d('real','C').mean()-d('permC','C').mean()
    ww=d('real','W').mean()-d('margin','W').mean()
    sp=np.sqrt(d('real','C').std()**2+d('real','W').std()**2+
               d('permC','C').std()**2+d('margin','W').std()**2)
    gap=cc-ww
    out.append(dict(Kc=Kc,Kw=Kw,C_c=cc.median(),W_c=ww.median(),gap=gap.median(),
                    spread2=2*sp.median(),general=int((gap>2*sp).sum()),
                    specific=int((-gap>2*sp).sum()),tied=int((gap.abs()<=2*sp).sum())))
S=pd.DataFrame(out); print(S.round(4).to_string(index=False))

tot=int(S.general.sum()+S.specific.sum()+S.tied.sum())
print("\n  CONDITIONAL KILL -- gates first")
g1=abs(g('real','C')-g('permC','C')).mean()>0.001 and abs(g('real','W')-g('margin','W')).mean()>0.001
g2=(dc.abs().max()<0.01) and (dw.abs().max()<0.02)
print(f"   (a) both nulls move their own quantity        : {'PASS' if g1 else 'FAIL'}")
print(f"   (b) doubly-destroyed world returns ~0 for both: {'PASS' if g2 else 'FAIL'}")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"\n   #71 reported, UNCORRECTED : general 186 / specific 32 / tied 58  (67.4% general)")
    print(f"   SYMMETRICALLY CORRECTED   : general {int(S.general.sum())} / specific "
          f"{int(S.specific.sum())} / tied {int(S.tied.sum())}  "
          f"({100*S.general.sum()/tot:.1f}% general) over {tot} comparisons")
    print(f"   median raw gap {T.gap_raw.median():+.4f}  ->  median corrected gap {T.gap_corrected.median():+.4f}")
    print(f"   the handicap #71 gave the specific side: {T.W_null.median():+.4f} of held-out skill")
    if S.specific.sum()>S.general.sum():
        print("\n   -> #71 INVERTS. Domain-SPECIFIC structure is larger once both sides are referred")
        print("      to their own nulls. The claim is withdrawn and the README row with it.")
    elif S.general.sum()>S.specific.sum()*2:
        print("\n   -> #71 SURVIVES its own correction. The conclusion was right and the method was")
        print("      lucky: it handicapped the side that lost anyway.")
    else:
        print("\n   -> #71 must be RESTATED as 'both exist, neither dominates'. The 186/276 margin was")
        print("      bought by the handicap, and the corrected grid does not support a ranking.")
print("\nN/A, with what it would require: a null destroying W while preserving C exactly (or the "
      "reverse) does not exist -- both live in the same residual matrix. The two corrections are "
      "each valid for their own quantity and are NOT claimed to be mutually independent.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
