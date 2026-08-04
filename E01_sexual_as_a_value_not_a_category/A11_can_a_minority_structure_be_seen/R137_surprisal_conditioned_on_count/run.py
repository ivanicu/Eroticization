import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R03 -- THE FIRST STATISTIC WHOSE ESTIMAND AND NULL ARE COMPATIBLE.

#93b: the per-person misfit T is confounded with the row sum -- a person with specific structure
picks a different NUMBER of options, and the only admissible null (fixed-margin) preserves exactly
that. So T could not work.

The statistic that survives that argument conditions on the count and asks about the CHOICE:

    S_i = mean over this person's picked options of  -log(base rate of that option)

S is invariant to HOW MANY options were picked and sensitive to WHICH. Fixed-margin randomisation
preserves each person's count and each option's base rate while destroying which person picked what
-- so it is not merely an admissible null for S, it is the exactly-matched one.

And the plant is finally the right shape too. #91's Gaussian low-rank plant was a mathematical
convenience; a fetish is an intense interest in something UNCOMMON. So the positive control plants a
minority whose extra picks land on the RAREST options in each block.

ESTIMAND        upper quantiles of the per-person mean surprisal S, real vs fixed-margin null.
IDENTIFICATION  identified: curveball preserves the conditioning variable (count) and the reference
                marginal (base rates) exactly, asserted per draw, and destroys only the assignment.
SCOPE           people with >=15 picks pooled across the 23 identified blocks.
WORLDS          rare-minority   real upper tail exceeds the null -> a minority concentrating on
                                uncommon options exists, which is the thing #91 said was invisible
                none            it does not, and the null is admissible because the control fires
KILL            threshold-free, per quantile, above 2x that quantile's own seed spread.
POSITIVE CTRL   GRADED: 5% of people given extra picks on the 3 rarest options of each block, at
                three strengths. Must be monotone and must NOT fire at strength 0.
NEGATIVE CTRL   the fixed-margin null itself, per draw, margins asserted.
NOISE FLOOR     6 seeds.
MULTIPLICITY    quantiles x 3 plant strengths x 6 seeds, published whole.
IMPOSSIBLE      attributing high surprisal to a cause -- unchanged; rare-option preference and
                careless responding both raise S.
"""
import pandas as pd, numpy as np, warnings, hashlib
warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,ppl=ppl)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
QS=[50,75,90,95,99]
print(f"targets {len(IDENT)}  people {len(ALLP)}",flush=True)
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
def surprisal(M):
    """mean -log(base rate) over each person's picked options. Base rates from the SAME matrix, so
    real and null are each scored against their own margins -- which curveball makes identical."""
    p=np.clip(M.mean(0),1e-4,1.); s=-np.log(p)
    tot=M@s; k=M.sum(1)
    return tot,k
def plant_rare(M,dens,extra,rng):
    """a minority whose EXTRA picks land on the rarest options -- what a fetish looks like."""
    Mw=M.copy(); n,m=M.shape
    order=np.argsort(M.mean(0))[:3]              # the 3 rarest options in this block
    on=np.flatnonzero(rng.random(n)<dens)
    for i in on:
        for j in order:
            if rng.random()<extra: Mw[i,j]=1.
    return Mw
rows=[]
for sd in range(1,7):
    acc={}
    arms=['real','cb']+[f'p{int(100*e)}' for e in [0.,0.3,0.6]]+[f'p{int(100*e)}cb' for e in [0.,0.3,0.6]]
    for w in arms: acc[w]=[np.zeros(len(ALLP)),np.zeros(len(ALLP))]
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        rg=np.random.default_rng(1500+sd)
        pack={'real':M,'cb':curveball(M,rg)}
        for e in [0.,0.3,0.6]:
            Mp=plant_rare(M,0.05,e,np.random.default_rng(1600+sd))
            pack[f'p{int(100*e)}']=Mp
            pack[f'p{int(100*e)}cb']=curveball(Mp,np.random.default_rng(1601+sd))
        assert np.allclose(pack['cb'].sum(0),M.sum(0)) and np.allclose(pack['cb'].sum(1),M.sum(1))
        for w,Mw in pack.items():
            tot,k=surprisal(Mw); acc[w][0][idx]+=tot; acc[w][1][idx]+=k
    ok=acc['real'][1]>=15
    for w in acc:
        S=acc[w][0][ok]/np.maximum(acc[w][1][ok],1)
        for qv in QS: rows.append(dict(seed=sd,world=w,qq=qv,S=float(np.percentile(S,qv))))
        rows.append(dict(seed=sd,world=w,qq=-1,S=float(S.mean())))
    print(f"  seed {sd}  n={int(ok.sum()):,}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def cmp(a,b,qv):
    r=D[(D.world==a)&(D.qq==qv)]['S']; n=D[(D.world==b)&(D.qq==qv)]['S']
    return r.mean(),n.mean(),r.mean()-n.mean(),2*np.sqrt(r.std()**2+n.std()**2)
lab={-1:'mean',50:'p50',75:'p75',90:'p90',95:'p95',99:'p99'}
print("\n=== GRADED POSITIVE CONTROL: 5% of people given extra picks on the 3 rarest options ===")
for e in [0,30,60]:
    r,n,d,s=cmp(f'p{e}',f'p{e}cb',95)
    r9,n9,d9,s9=cmp(f'p{e}',f'p{e}cb',99)
    print(f"  strength {e/100:.1f}:  p95 diff {d:+.4f} (2xsp {s:.4f}) {'FIRES' if d>s else '-'}   "
          f"p99 diff {d9:+.4f} (2xsp {s9:.4f}) {'FIRES' if d9>s9 else '-'}")
print("\n=== REAL vs FIXED-MARGIN NULL (6 seeds) ===")
for qv in [-1,50,75,90,95,99]:
    r,n,d,s=cmp('real','cb',qv)
    print(f"  {lab[qv]:5s} real {r:.4f}  null {n:.4f}  diff {d:+.4f}  2xspread {s:.4f}  "
          f"{'RESOLVABLE' if abs(d)>s else 'no'}")
_,_,d0,s0=cmp('p0','p0cb',95); _,_,d6,s6=cmp('p60','p60cb',95)
print("\n  CONDITIONAL KILL -- gates first")
g1=d0<=s0; g2=d6>s6
print(f"   (a) control does NOT fire at strength 0 : {'PASS' if g1 else 'FAIL'} ({d0:+.4f} vs {s0:.4f})")
print(f"   (b) control FIRES at strength 0.6       : {'PASS' if g2 else 'FAIL'} ({d6:+.4f} vs {s6:.4f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    hits=[qv for qv in [50,75,90,95,99] if (lambda t:(t[2]>t[3]))(cmp('real','cb',qv))]
    print(f"\n   quantiles where the real tail exceeds the null: {hits if hits else 'NONE'}")
    if hits:
        print("   -> A MINORITY CONCENTRATING ON RARE OPTIONS IS PRESENT, measured by a statistic")
        print("      whose null preserves exactly what it conditions on.")
    else:
        print("   -> NO SUCH MINORITY, and this is the first admissible null in the arc: the")
        print("      instrument fires on a planted rare-option minority and the null preserves both")
        print("      the pick count and the base rates by construction.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
