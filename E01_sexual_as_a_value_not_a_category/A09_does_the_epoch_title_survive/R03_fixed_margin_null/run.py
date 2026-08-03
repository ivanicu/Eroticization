import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A09 R03 -- THE FLOOR, BUILT FROM A WORLD THAT IS ACTUALLY LIKE THE DATA.

R02's identification gate refused all 32 blocks and it was right to. Its synthetic null drew
Bernoulli cells from an additive probability model, and that world is far WEAKER than the data on
both main effects: item -0.157 below real (median), person effect only 0.3x real. A floor built
from a world unlike the data is not a floor.

The correct null for a binary matrix is FIXED-MARGIN RANDOMISATION (curveball). It preserves every
row sum and every column sum EXACTLY, so both main effects are matched BY CONSTRUCTION rather than
in expectation, and it destroys the person x item interaction and nothing else.

That construction also buys a graded positive control that R02 could not have: run only a FRACTION
of the mixing trades. f=0 is the real matrix (all interaction), f=5n is fully mixed (none), and
margins are exact at every f. The dose is "fraction of the real structure retained" -- a relative
scale, which is the right one, because it calibrates against the structure actually present rather
than against a magnitude I invented.

ESTIMAND        X_c = X_real - X_null(fixed margins, same block), on the held-out R2 scale,
                directly comparable to the ITEM main effect I_real measured on the same cells.
IDENTIFICATION  identified by construction: the null has identical margins, so I and P are matched.
                VERIFIED not assumed -- |I_real - I_null| and |P_real - P_null| are reported.
SCOPE           as R01/R02. K in {1,2}: R01 measured that overfit grows monotonically with K.
WORLDS          A content/base-rate dominates -> I - X_c > 2x spread
                B individualised value        -> X_c - I > 2x spread
                C neither                     -> |X_c - I| <= 2x spread
KILL            THRESHOLD-FREE: the verdict is the ordering of two quantities on one scale,
                declared only when the gap exceeds 2x its own seed spread (#34).
POSITIVE CTRL   graded, margin-preserving: X must fall monotonically as trades increase
                (more real structure destroyed -> less interaction recovered). It must NOT
                already be at the floor at f=0, or the instrument is blind.
NEGATIVE CTRL   the f=5n cell IS the negative control, per block, with exact margins.
PLACEBO         permuted person component must contribute <= 0 (direction corrected in R02).
NOISE FLOOR     3 masks x 3 randomisation draws.
MULTIPLICITY    32 blocks x 2 K x 5 doses x 3 seeds, reported whole.
IMPOSSIBLE      an ABSOLUTE dose scale -- that needs a margin-preserving generator with a planted
                interaction of known magnitude, which fixed-margin randomisation cannot provide.
                Reported as N/A with what it would require, not as "planned".
"""
import pandas as pd, numpy as np, warnings, hashlib
from numpy.linalg import svd
from math import factorial
warnings.filterwarnings('ignore')

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
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
    RAW[q.qi]=M
print(f"blocks {len(RAW)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]; DOSES=[0.,0.25,0.5,1.,5.]

def curveball(M,rng,per_row):
    """fixed-margin randomisation. every row sum and column sum is preserved EXACTLY."""
    if per_row<=0: return M.copy()
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

def components(M,obs,K):
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P; F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(25):
        U,S,Vt=svd(F,full_matrices=False); F=np.where(obs,Rres,(U[:,:K]*S[:K])@Vt[:K])
    U,S,Vt=svd(F,full_matrices=False)
    return gm,I,P,(U[:,:K]*S[:K])@Vt[:K]

def shap(M,K,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    gm,I,P,X=components(M,obs,K)
    comp={'I':np.broadcast_to(I,M.shape),'P':np.broadcast_to(P,M.shape),'X':X}
    base=np.mean((M[he]-gm)**2); v={}
    for bits in range(8):
        S=frozenset([c for j,c in enumerate('IPX') if bits>>j&1])
        p=np.clip(gm+sum(comp[c] for c in S),0.,1.) if S else np.full(M.shape,gm)
        v[S]=1.-np.mean((M[he]-p[he])**2)/base if base>0 else np.nan
    out={}
    for c in 'IPX':
        o=[x for x in 'IPX' if x!=c]; tot=0.
        for S in [(),(o[0],),(o[1],),tuple(o)]:
            tot+=factorial(len(S))*factorial(2-len(S))/6.*(v[frozenset(S+(c,))]-v[frozenset(S)])
        out[c]=tot
    return out

rows=[]
for qi,(q,M) in enumerate(RAW.items()):
    for sd in SEEDS:
        rgen=np.random.default_rng(7000+sd)
        for f in DOSES:
            Mf=curveball(M,rgen,f) if f>0 else M
            assert np.allclose(Mf.sum(0),M.sum(0)) and np.allclose(Mf.sum(1),M.sum(1)),"margins broken"
            for K in [1,2]:
                rows.append(dict(q=q,K=K,seed=sd,f=f,n=M.shape[0],m=M.shape[1],**shap(Mf,K,sd)))
    print(f"  block {qi+1}/{len(RAW)}",flush=True)
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== MARGIN MATCH (verified, not assumed): main effects at f=0 vs f=5 ===")
for K in [1,2]:
    d=D[D.K==K].groupby(['q','f'])[['I','P']].mean().unstack('f')
    print(f"  K={K}  |I_real - I_null| median {abs(d[('I',0.)]-d[('I',5.)]).median():.5f}"
          f"   |P_real - P_null| median {abs(d[('P',0.)]-d[('P',5.)]).median():.5f}")

print("\n=== GRADED POSITIVE CONTROL: interaction recovered vs fraction of structure destroyed ===")
print(D.groupby(['K','f'])[['I','P','X']].mean().round(4).to_string())

for K in [1,2]:
    d=D[D.K==K]
    mono=d.groupby('f').X.mean().values
    g_mono=all(mono[i]>=mono[i+1]-1e-9 for i in range(len(mono)-1))
    g_notfloor=(mono[0]-mono[-1])>2*d[d.f==0].groupby('q').X.std().median()
    piv=d.groupby(['q','f'])[['I','X']].mean().unstack('f')
    res=pd.DataFrame({'n':d.groupby('q').n.first(),'m':d.groupby('q').m.first(),
                      'I':piv[('I',0.)],'X_real':piv[('X',0.)],'X_null':piv[('X',5.)]})
    res['X_c']=res.X_real-res.X_null; res['gap']=res.X_c-res.I
    res['dI']=(piv[('I',0.)]-piv[('I',5.)]).abs()
    res['spread']=np.sqrt(d[d.f==0].groupby('q').X.std()**2+d[d.f==5.].groupby('q').X.std()**2)
    res['ok']=res.dI<=0.01
    print(f"\n{'='*78}\nK={K}")
    print(res.sort_values('X_c').round(4).to_string())
    print(f"\n  CONDITIONAL KILL -- gates first")
    print(f"   (a) graded control monotone in f          : {'PASS' if g_mono else 'FAIL'}  {mono.round(4)}")
    print(f"   (b) not already at floor at f=0           : {'PASS' if g_notfloor else 'FAIL'}")
    print(f"   (c) margins matched (|dI|<=0.01, per block): {int(res.ok.sum())}/{len(res)}")
    v=res[res.ok]
    if not(g_mono and g_notfloor) or len(v)<len(res)//2:
        print("   -> UNVERIFIED, and that is not an acquittal."); continue
    mI,mX,mg,ms=v.I.median(),v.X_c.median(),v.gap.median(),v.spread.median()
    print(f"\n   identified blocks {len(v)}/{len(res)}   medians:")
    print(f"     ITEM main effect        I   = {mI:+.4f}")
    print(f"     interaction, corrected  X_c = {mX:+.4f}     ratio X_c/I = {mX/mI:.3f}")
    print(f"     gap X_c - I = {mg:+.4f}    2x seed spread = {2*ms:.4f}")
    print(f"     blocks with X_c>0: {int((v.X_c>0).sum())}/{len(v)}   with X_c>I: {int((v.gap>0).sum())}/{len(v)}")
    if abs(mg)<=2*ms:
        print("   -> WORLD C: item and interaction NOT DISTINGUISHABLE in size.")
    elif mg>0:
        print(f"   -> WORLD B: interaction exceeds the item main effect by {mg:+.4f}. epoch title SURVIVES.")
    else:
        print(f"   -> WORLD A: the ITEM MAIN EFFECT exceeds the interaction by {-mg:+.4f} "
              f"({mI/max(mX,1e-9):.1f}x). 'a value, NOT a category' IS FALSE AS STATED -- and 105 "
              f"rounds were run on the residual left after deleting the larger component.")
print("\nN/A, with what it would require: an ABSOLUTE dose scale needs a margin-preserving generator "
      "with a planted interaction of known magnitude. Fixed-margin randomisation cannot provide one.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
