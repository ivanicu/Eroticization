import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R05 -- THE CONTROL AT A MAGNITUDE THAT WAS CALCULATED INSTEAD OF GUESSED.

#94c: three positive controls in a row failed because I chose the plant strength by intuition. #94's
closing paragraph derives it instead: p95 has a 2x seed spread of 0.0017 and the real effect is
+0.0503, so an informative plant must move p95 by >= 0.0017 and a comparable one by ~0.05. With ~35
pooled picks per carrier that is a surprisal shift of ~1.75 nats per carrier -- two swaps from a
median-prevalence option to a much rarer one IN EVERY BLOCK THE CARRIER ENTERS, not two swaps in one.

The round reports the ACHIEVED per-carrier surprisal shift beside the intended one, so a
below-MDE plant can never again be mistaken for a null result.

ESTIMAND        does the p95 of the surprisal distribution rise when a 5% minority is given a
                rare-option preference of a magnitude computed to be detectable?
IDENTIFICATION  identified; the achieved shift is measured, not assumed.
WORLDS          works    the plant fires -> #94a is licensed and the real signal is a finding
                broken   a plant at a magnitude computed to be detectable does NOT fire -> the
                         statistic is broken in a way three guessed plants could not have revealed
KILL            threshold-free: the plant fires iff its p95 exceeds the no-plant arm's by more than
                2x the pooled seed spread, and the ACHIEVED per-carrier shift is published beside it.
POSITIVE CTRL   the graded ladder IS it; 0 swaps must be a genuine no-op on a separate matrix.
NEGATIVE CTRL   fixed-margin null per arm, row sums asserted.
NOISE FLOOR     6 seeds.
IMPOSSIBLE      unchanged.
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
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in IDENT])); PM={p:i for i,p in enumerate(ALLP)}
QS=[50,75,90,95,99]; DENS=0.05
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
def plant_every_block(M,carriers,nsw,rng):
    """IN THIS BLOCK: each carrier swaps nsw picks from above-median-prevalence to the 3 rarest."""
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in carriers:
        done=0
        for _ in range(6*nsw):
            if done>=nsw: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; done+=1
    return Mw
rows=[]; shifts=[]
for sd in range(1,7):
    rgc=np.random.default_rng(3300+sd)
    carrier=rgc.random(len(ALLP))<DENS                     # the SAME people in every block
    arms=['real','cb']+[f'n{k}' for k in [0,1,2]]+[f'n{k}cb' for k in [0,1,2]]
    acc={w:[np.zeros(len(ALLP)),np.zeros(len(ALLP))] for w in arms}
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        loc=np.flatnonzero(carrier[idx])
        rg=np.random.default_rng(3400+sd)
        pack={'real':M,'cb':curveball(M,rg)}
        for k in [0,1,2]:
            Mp=plant_every_block(M,loc,k,np.random.default_rng(3500+sd))
            pack[f'n{k}']=Mp; pack[f'n{k}cb']=curveball(Mp,np.random.default_rng(3501+sd))
        assert np.allclose(pack['cb'].sum(1),M.sum(1))
        for w,Mw in pack.items():
            acc[w][0][idx]+=Mw@ref; acc[w][1][idx]+=Mw.sum(1)
    ok=acc['real'][1]>=15
    base=acc['n0'][0][ok]/np.maximum(acc['n0'][1][ok],1)
    for k in [1,2]:
        v=acc[f'n{k}'][0][ok]/np.maximum(acc[f'n{k}'][1][ok],1)
        c=carrier[ok]
        shifts.append(dict(seed=sd,nsw=k,achieved=float((v-base)[c].mean()),
                           achieved_nats=float(((acc[f'n{k}'][0][ok]-acc['n0'][0][ok])[c]).mean())))
    for w in acc:
        S=acc[w][0][ok]/np.maximum(acc[w][1][ok],1)
        for qv in QS: rows.append(dict(seed=sd,world=w,qq=qv,S=float(np.percentile(S,qv))))
    print(f"  seed {sd}  n={int(ok.sum()):,}  carriers={int(carrier[ok].sum()):,}",flush=True)
D=pd.DataFrame(rows); SH=pd.DataFrame(shifts)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def cmp(a,b,qv):
    r=D[(D.world==a)&(D.qq==qv)]['S']; n=D[(D.world==b)&(D.qq==qv)]['S']
    return r.mean(),n.mean(),r.mean()-n.mean(),2*np.sqrt(r.std()**2+n.std()**2)
print("\n=== ACHIEVED PER-CARRIER SHIFT (target from #94: ~1.75 nats, ~0.05 in mean surprisal) ===")
print(SH.groupby('nsw')[['achieved','achieved_nats']].mean().round(4).to_string())
print("\n=== PLANT vs NO-PLANT, at p95 and p99 ===")
for k in [0,1,2]:
    r,n,d,s=cmp(f'n{k}',f'n{k}cb',95); r9,n9,d9,s9=cmp(f'n{k}',f'n{k}cb',99)
    print(f"  {k} swaps/block: p95 real-vs-null {d:+.4f} (2xsp {s:.4f})   p99 {d9:+.4f} (2xsp {s9:.4f})")
_,_,d0,s0=cmp('n0','n0cb',95); _,_,d1,s1=cmp('n1','n1cb',95); _,_,d2,s2=cmp('n2','n2cb',95)
_,_,dr,sr=cmp('real','cb',95)
print("\n  CONDITIONAL KILL -- gates first")
ach=SH[SH.nsw==2].achieved.mean()
g1=abs(d0-dr)<max(s0,sr)
g2=ach>0.0017
print(f"   (a) 0-swap arm is a genuine no-op        : {'PASS' if g1 else 'FAIL'} ({d0:+.4f} vs real {dr:+.4f})")
print(f"   (b) the plant ACHIEVED a detectable shift: {'PASS' if g2 else 'FAIL'} "
      f"({ach:+.4f} in mean surprisal vs the 0.0017 needed)")
if not(g1 and g2): print("   -> UNVERIFIED: the plant is still below the resolution, and #94a stays provisional.")
else:
    fires=(d2-d0)>2*max(s0,s2)
    print(f"\n   plant effect on p95: {d2-d0:+.4f}   2x spread {2*max(s0,s2):.4f}  -> "
          f"{'FIRES' if fires else 'does NOT fire'}")
    if fires:
        print("\n   -> #94a IS LICENSED. The statistic responds to a planted rare-option minority at a")
        print("      magnitude computed to be detectable, and the real data shows the same signature")
        print(f"      at p95 (+{dr:.4f}, {dr/(sr/2):.0f}x its own spread).")
    else:
        print("\n   -> THE STATISTIC IS BROKEN. A plant that achieved a detectable per-carrier shift")
        print("      does not move the population quantile, so the real p95 excess cannot be read as")
        print("      a minority either. Three guessed plants could not have shown this.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
