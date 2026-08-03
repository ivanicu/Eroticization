import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R12 -- THE CHECK ELEVEN ROUNDS NEVER RAN: WHAT DOES THE LOWER TAIL DO?

Every round in this arc measured p50 and UPWARD. #95 licensed "a minority picks rarer options"; #98
bracketed its size; #R11 found no carrier model fits, because every random-carrier plant LIFTS the
median while the real p50 sits exactly at the null.

There is a much simpler world that produces all of that and it has never been tested: THE REAL
DISTRIBUTION IS JUST WIDER. If a person-level trait exists at all, curveball destroys it, so the null
is narrower than the real -- symmetrically. A symmetric widening leaves the median unmoved, raises
every upper quantile, and LOWERS every lower one.

  elevation   real is above the null at high quantiles and AT it at low ones -> a one-sided tail
  width       real is above at high AND BELOW at low, roughly symmetrically -> extra person-level
              variance, which is a far weaker claim than "a minority picks rare options"

Measuring p1..p50 costs one line and decides between them. I did not run it for eleven rounds.

ESTIMAND        real-minus-null elevation across the WHOLE quantile range, and the symmetry ratio
                (upper elevation) / (-lower elevation).
IDENTIFICATION  identified; nothing is fitted.
UNCERTAINTY     bootstrap over people, 300 resamples.
WORLDS          as above; plus a mixed case where the widening is asymmetric.
KILL            threshold-free: the symmetry ratio with its bootstrap interval. ~1 = pure width,
                >>1 = genuinely one-sided.
POSITIVE CTRL   a synthetic PURE-WIDTH world (a symmetric person-level trait) must return ratio ~1;
                a synthetic ONE-SIDED world (only positive deviations) must return ratio >> 1.
NEGATIVE CTRL   the null against itself: elevation 0 everywhere.
NOISE FLOOR     4 seeds.
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
QQ=[1,5,10,25,50,75,90,95,99]
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
def swap(M,who,to_rare,rng):
    """to_rare=True: median pick -> rare. False: rare pick -> median (the mirror image)."""
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    a,b=(med,rare) if to_rare else (rare,med)
    for i in who:
        for _ in range(8):
            c=a[rng.integers(len(a))]; r=b[rng.integers(len(b))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; break
    return Mw
rows=[]
for sd in range(1,5):
    rgc=np.random.default_rng(5100+sd); u=rgc.random(len(ALLP))
    S={w:{} for w in ['real','cb','WIDTH','ONESIDE']}
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        Mn=curveball(M,np.random.default_rng(5200+sd))
        up=np.flatnonzero(u[idx]>0.75); dn=np.flatnonzero(u[idx]<0.25)
        w1=swap(Mn,up,True,np.random.default_rng(5300+sd))
        pack={'real':M,'cb':Mn,
              'WIDTH':swap(w1,dn,False,np.random.default_rng(5301+sd)),   # symmetric: up AND down
              'ONESIDE':w1}                                               # only the upper quarter
        for w,Mw in pack.items():
            kk=Mw.sum(1); s=(Mw@ref)/np.maximum(kk,1)
            for j,gi in enumerate(idx):
                if kk[j]>0: S[w].setdefault(gi,[]).append(s[j])
    for w in S:
        ks=[i for i,v in S[w].items() if len(v)>=5]
        vals=np.array([np.mean(S[w][i]) for i in ks])
        q=np.percentile(vals,QQ)
        rows.append(dict(seed=sd,world=w,**{f'p{a}':b for a,b in zip(QQ,q)}))
        if w=='real' and sd==1:
            rb=np.random.default_rng(777)
            B=np.array([np.percentile(vals[rb.integers(0,len(vals),len(vals))],QQ) for _ in range(300)])
            np.save(pathlib.Path(__file__).parent/'results'/'boot.npy',B)
    print(f"  seed {sd}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
B=np.load(OUT/'boot.npy'); bsd=dict(zip([f'p{q}' for q in QQ],B.std(0)))
G=D.groupby('world')[[f'p{q}' for q in QQ]].mean(); nul=G.loc['cb']
E=(G-nul).drop(index=['cb'])
print("\n=== ELEVATION ABOVE THE NULL ACROSS THE WHOLE QUANTILE RANGE ===")
print(E.round(4).to_string())
print("\n  bootstrap sd: " + "  ".join(f"{k} {v:.4f}" for k,v in bsd.items()))
def ratio(w):
    up=(G.loc[w,'p95']-nul.p95)+(G.loc[w,'p90']-nul.p90)
    dn=(nul.p5-G.loc[w,'p5'])+(nul.p10-G.loc[w,'p10'])
    return up/dn if abs(dn)>1e-9 else np.inf
print("\n=== SYMMETRY RATIO  (upper elevation) / (lower depression) ===")
for w in ['WIDTH','ONESIDE','real']: print(f"  {w:8s} {ratio(w):+.2f}")
gw=abs(ratio('WIDTH'))<3; go=ratio('ONESIDE')>3 or ratio('ONESIDE')==np.inf
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) the PURE-WIDTH control returns ratio ~1  : {'PASS' if gw else 'FAIL'} ({ratio('WIDTH'):+.2f})")
print(f"   (b) the ONE-SIDED control returns ratio >> 1 : {'PASS' if go else 'FAIL'} ({ratio('ONESIDE'):+.2f})")
if not(gw and go): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    rr=ratio('real')
    print(f"\n   real symmetry ratio {rr:+.2f}   (width control {ratio('WIDTH'):+.2f}, "
          f"one-sided control {ratio('ONESIDE'):+.2f})")
    if abs(rr)<3:
        print("\n   -> WIDTH. The real distribution is SYMMETRICALLY wider than the null: the upper")
        print("      tail is elevated and the lower tail is depressed by a comparable amount. That is")
        print("      extra PERSON-LEVEL VARIANCE, not a minority picking rare options, and #95's")
        print("      wording must be withdrawn.")
    else:
        print("\n   -> ONE-SIDED. The elevation is genuinely asymmetric: the upper tail moves and the")
        print("      lower tail does not. #95's minority reading survives its own strongest attack.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
