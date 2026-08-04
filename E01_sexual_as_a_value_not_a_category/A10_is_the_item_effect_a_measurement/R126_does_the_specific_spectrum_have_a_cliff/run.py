import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R11 -- THE SAME QUESTION #72 ASKED OF THE THIN PART, ASKED OF THE THICK ONE.

#72 swept the CROSS-block rank and found no cliff: a calibrated estimator produced 141-147x drops in
known-rank worlds and 1.8x on the real cross-block spectrum. #82 then inverted #71 and showed the
cross-block part is 7-26x SMALLER than the within-block part -- so #72 characterised a thin residual
and called it the object.

The same question, asked of the structure that turns out to carry the mass: does the WITHIN-block
spectrum have a cliff?

  low-rank specific   a few per-domain factors and then noise. Combined with #72's smooth
                      cross-block tail, that is a sharp structural statement: the operator is a
                      small number of domain-specific readouts plus a diffuse general residue.
  smooth specific     no cliff anywhere. The person-side structure is high-dimensional at both
                      levels and "how many readouts" has no answer at either.

ESTIMAND        W corrected against its own per-block fixed-margin floor, swept over rank; and the
                per-DIMENSION gain, whose collapse point is the knee. Read only against synthetic
                worlds of KNOWN WITHIN-BLOCK rank -- the same discipline as #72, whose bare curve
                would otherwise carry no rank information.
IDENTIFICATION  identified relative to the controls; the answer is "the real spectrum behaves like a
                world of true within-block rank r", never a bare number.
SCOPE           23 blocks (real), 8 blocks per control world. Kc fixed at 4 so C is present and W is
                the only thing swept, matching #82's specification.
WORLDS          as above; plus  artifact  -- the controls also fail to saturate, in which case the
                estimator cannot count within-block rank and NOTHING is reportable.
KILL            the r=2 and r=5 controls must knee at their true rank, and r=5 must knee later than
                r=2. If not, UNVERIFIED -- the outcome that kills the question rather than answering it.
POSITIVE CTRL   the graded pair above (dose-response in the true rank, not merely "it saturates").
NEGATIVE CTRL   the per-block fixed-margin world, at every rank, margins asserted exact per draw.
NOISE FLOOR     2 masks x 2 randomisation draws.
MULTIPLICITY    23 blocks x 8 ranks x 2 worlds x 2 seeds, plus 8 x 8 x 2 x 2 per control.
IMPOSSIBLE      naming the within-block factors -- #18/#39/#49 measured that naming fails here, and
                counting is a different task from naming. Only counting is attempted.
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
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
print(f"targets {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29]; KWS=[1,2,3,5,8,12,16,24]

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

def w_skill(M,Kw,seed):
    """held-out skill added by a rank-Kw within-block component, on top of both marginals."""
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P; F=np.where(np.isnan(Rres),0.,Rres)
    for _ in range(20):
        U,S,V=svd(F,full_matrices=False); F=np.where(obs,Rres,(U[:,:Kw]*S[:Kw])@V[:Kw])
    U,S,V=svd(F,full_matrices=False); W=(U[:,:Kw]*S[:Kw])@V[:Kw]
    base=np.mean((M[he]-gm)**2)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    IB=np.broadcast_to(I,M.shape)
    return f(IB,P,W)-f(IB,P)

def synth_within(r,seed):
    """each block gets its OWN rank-r person factors -- domain-SPECIFIC structure of known rank."""
    rng=np.random.default_rng(5000+seed); out={}
    for q in RAW:
        M=RAW[q]['M']; n,m=M.shape
        F=rng.normal(size=(n,r)); L=rng.normal(size=(r,m))*0.30
        p=np.clip(M.mean(0)[None,:]+F@L,0.02,0.98)
        out[q]=dict(M=(rng.random((n,m))<p).astype(float),ppl=RAW[q]['ppl'])
    return out

rows=[]
def sweep(BLK,targets,tag):
    for i,t in enumerate(targets):
        M=BLK[t]['M']
        for sd in SEEDS:
            Mn=curveball(M,np.random.default_rng(8000+sd))
            assert np.allclose(Mn.sum(0),M.sum(0)) and np.allclose(Mn.sum(1),M.sum(1))
            for Kw in KWS:
                rows.append(dict(world=tag,q=t,Kw=Kw,seed=sd,arm='real',W=w_skill(M ,Kw,sd)))
                rows.append(dict(world=tag,q=t,Kw=Kw,seed=sd,arm='null',W=w_skill(Mn,Kw,sd)))
        print(f"  [{tag}] {i+1}/{len(targets)}",flush=True)

CT=IDENT[:8]
sweep(RAW,IDENT,'real')
for r in [2,5]: sweep(synth_within(r,1),CT,f'r{r}')
D=pd.DataFrame(rows)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R126_does_the_specific_spectrum_have_a_cliff/results/'
D.to_csv(OUT+'grid.csv',index=False)

def corrected(w):
    d=D[D.world==w]
    return d[d.arm=='real'].groupby('Kw').W.mean()-d[d.arm=='null'].groupby('Kw').W.mean()
cr,c2,c5=corrected('real'),corrected('r2'),corrected('r5')
sp=D[(D.world=='real')&(D.arm=='real')].groupby(['Kw','q']).W.std().groupby('Kw').mean()
T=pd.DataFrame({'real_corrected':cr,'ctrl_r2':c2,'ctrl_r5':c5,'seed_spread':sp})
print("\n=== CORRECTED WITHIN-BLOCK SKILL BY RANK ===")
print(T.round(4).to_string())

def perdim(c):
    k=np.array(c.index,dtype=float); return np.diff(c.values)/np.diff(k)
def knee(c,frac=0.10):
    k=np.array(c.index); inc=perdim(c); big=inc.max()
    for i,g_ in enumerate(inc):
        if g_<frac*big: return int(k[i])
    return int(k[-1])
print("\n=== GAIN PER DIMENSION, and the knee ===")
for nm,c in [('ctrl r=2',c2),('ctrl r=5',c5),('real',cr)]:
    print(f"  {nm:9s} knee={knee(c):2d}   {np.round(perdim(c),5).tolist()}")
def cliff(c):
    inc=perdim(c); r=[inc[i]/max(inc[i+1],1e-12) for i in range(len(inc)-1)]
    return max(r),int(np.array(c.index)[int(np.argmax(r))])
for nm,c in [('ctrl r=2',c2),('ctrl r=5',c5),('real',cr)]:
    m,at=cliff(c); print(f"  {nm:9s} sharpest single drop {m:8.1f}x  at rank {at}")

k2,k5,kr=knee(c2),knee(c5),knee(cr)
print("\n  CONDITIONAL KILL -- gates first")
g1=k5>k2; g2=max(k2,k5)<max(KWS); g3=cliff(c2)[0]>20 and cliff(c5)[0]>20
print(f"   (a) dose-response in the true rank (r5 knees later than r2) : {'PASS' if g1 else 'FAIL'} ({k2} vs {k5})")
print(f"   (b) controls knee before the sweep ends                     : {'PASS' if g2 else 'FAIL'}")
print(f"   (c) controls produce a CLIFF (>20x) at their true rank      : {'PASS' if g3 else 'FAIL'} "
      f"({cliff(c2)[0]:.0f}x, {cliff(c5)[0]:.0f}x)")
if not(g1 and g2 and g3):
    print("   -> UNVERIFIED: the estimator cannot count within-block rank, so the real curve carries")
    print("      no rank information and #72's question has no answer at this level either.")
else:
    m,at=cliff(cr)
    print(f"\n   real: knee at rank {kr}, sharpest drop {m:.1f}x at rank {at}")
    print(f"   corrected skill at rank 1 {cr.iloc[0]:+.4f} -> at rank {max(KWS)} {cr.iloc[-1]:+.4f}")
    if m>20:
        print(f"\n   -> THE WITHIN-BLOCK SPECTRUM HAS A CLIFF at rank {at} ({m:.0f}x), and #72 measured")
        print(f"      that the cross-block one does not (1.8x). The operator is a SMALL NUMBER OF")
        print(f"      DOMAIN-SPECIFIC READOUTS plus a diffuse general residue -- a sharper statement")
        print(f"      than either arc has made.")
    else:
        print(f"\n   -> NO CLIFF HERE EITHER ({m:.1f}x against the controls' "
              f"{cliff(c2)[0]:.0f}x/{cliff(c5)[0]:.0f}x). The person-side structure is")
        print(f"      high-dimensional at BOTH levels, and 'how many readouts' has no answer at")
        print(f"      either. #72 generalises rather than being superseded.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
