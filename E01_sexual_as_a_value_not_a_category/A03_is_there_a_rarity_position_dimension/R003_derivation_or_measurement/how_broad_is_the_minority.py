import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A11 R11 -- HOW BROAD IS IT? THE TWO CONSTRAINTS PIN THE FRACTION.

#98 excluded both a gradient across everyone (it lifts the median, and the real median sits exactly
at the null) and a sharp 3% minority (it overshoots the held-out tail shape 4x). What survives is a
BROAD, MILD minority, and #98's two facts pin its size before a single run:

    p50 sits at the null       ->  the median person is NOT a carrier  ->  f < 0.5
    p90 is elevated by +0.034  ->  the p90 person IS a carrier         ->  f > 0.10

So the answer is inside [0.10, 0.50], and inside that window fraction and intensity trade off along
exactly the axis the two constraints separate. This round sweeps that window on both axes and reads
the pair off, then tests the fit on the quantiles it did NOT match.

ESTIMAND        the (carrier fraction f, per-carrier intensity c) that reproduces the real p50 and
                p90; and the held-out error at p95 and p99.
IDENTIFICATION  identified by the two constraints being independent -- f moves WHERE the elevation
                starts, c moves HOW BIG it is. Checked by requiring the grid to separate on both.
SCOPE           people entering >=5 of the 23 identified blocks. n ~ 12,000.
UNCERTAINTY     bootstrap over people (#97c), 200 resamples.
WORLDS          narrow    best f near 0.10 -> a small minority, mild
                broad     best f near 0.40 -> a large minority, very mild
                nofit     no (f,c) reproduces both anchors -> the one-swap model is wrong in shape,
                          not merely in size, and that is the finding
KILL            threshold-free: the held-out error at p95 and p99 in units of the bootstrap sd.
POSITIVE CTRL   f=0 must reproduce the null exactly; the grid must separate on p50 across f (the
                constraint that caps f) and on p90 across c.
NEGATIVE CTRL   the fixed-margin null, shared origin for every arm.
NOISE FLOOR     4 seeds.
MULTIPLICITY    5 fractions x 3 intensities x 4 seeds, published whole.
IMPOSSIBLE      distinguishing "a real subgroup" from "a smooth trait truncated at zero" -- both
                produce a population where the median is unmoved and a tail is elevated. This round
                measures the SIZE of the elevated part, not whether it is a natural kind.
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

def one_swap(M,who,rng):
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in who:
        for _ in range(8):
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; break
    return Mw

FRACS=[0.0,0.10,0.15,0.22,0.32,0.45]; COVER=[0.25,0.5,1.0]
rows=[]
for sd in range(1,5):
    rgc=np.random.default_rng(3100+sd)
    u=rgc.random(len(ALLP))                                  # each person's rank on the latent
    blockdraw={c:rgc.random((len(ALLP),len(IDENT)))<c for c in COVER}
    arms=['real','cb']+[f'f{f}c{c}' for f in FRACS for c in COVER]
    S={w:{} for w in arms}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        Mn=curveball(M,np.random.default_rng(3200+sd))
        pack={'real':M,'cb':Mn}
        for f in FRACS:
            car=u[idx]<f
            for c in COVER:
                who=np.flatnonzero(car & blockdraw[c][idx,bi]) if f>0 else np.array([],dtype=int)
                pack[f'f{f}c{c}']=one_swap(Mn,who,np.random.default_rng(3300+sd)) if f>0 else Mn
        for w,Mw in pack.items():
            kk=Mw.sum(1); s=(Mw@ref)/np.maximum(kk,1)
            for j,gi in enumerate(idx):
                if kk[j]>0: S[w].setdefault(gi,[]).append(s[j])
    for w in arms:
        ks=[i for i,v in S[w].items() if len(v)>=5]
        vals=np.array([np.mean(S[w][i]) for i in ks])
        q=np.percentile(vals,[50,75,90,95,99])
        rows.append(dict(seed=sd,world=w,n=len(vals),p50=q[0],p75=q[1],p90=q[2],p95=q[3],p99=q[4]))
        if w=='real' and sd==1:
            rb=np.random.default_rng(4242)
            B=np.array([np.percentile(vals[rb.integers(0,len(vals),len(vals))],[50,75,90,95,99])
                        for _ in range(200)])
            np.save(pathlib.Path(__file__).parent/'results'/'boot.npy',B)
    print(f"  seed {sd}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
B=np.load(OUT/'boot.npy'); bsd=dict(zip(['p50','p75','p90','p95','p99'],B.std(0)))
G=D.groupby('world')[['p50','p75','p90','p95','p99']].mean()
r=G.loc['real']; nul=G.loc['cb']
print("\n=== BOOTSTRAP SD OF THE REAL QUANTILES (200 resamples over people) ===")
print("  " + "   ".join(f"{k} {v:.4f}" for k,v in bsd.items()))
print("\n=== THE GRID, as ELEVATION ABOVE THE NULL ===")
E=(G-nul).drop(index=['cb'])
E.loc['REAL']=r-nul
print(E.round(4).to_string())
cand=[w for w in G.index if w.startswith('f') and not w.startswith('f0.0')]
def anchor(w): return abs(G.loc[w,'p50']-r.p50)/bsd['p50']+abs(G.loc[w,'p90']-r.p90)/bsd['p90']
best=min(cand,key=anchor)
f_hat=float(best.split('c')[0][1:]); c_hat=float(best.split('c')[1])
held=abs(G.loc[best,'p95']-r.p95)/bsd['p95']+abs(G.loc[best,'p99']-r.p99)/bsd['p99']
print(f"\n  best anchor match: {best}  -> fraction {f_hat:.0%} of people, elevated in "
      f"{c_hat:.0%} of their blocks")
print(f"     anchors (p50,p90) off by {anchor(best):.1f} bootstrap sd; "
      f"HELD-OUT (p95,p99) off by {held:.1f}")
zero=abs(G.loc['f0.0c1.0','p50']-nul.p50)<1e-9
sep_f=(G.loc['f0.45c1.0','p50']-G.loc['f0.1c1.0','p50'])/bsd['p50']
sep_c=(G.loc['f0.22c1.0','p90']-G.loc['f0.22c0.25','p90'])/bsd['p90']
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) f=0 reproduces the null exactly     : {'PASS' if zero else 'FAIL'}")
print(f"   (b) the grid separates on f (via p50)   : {'PASS' if abs(sep_f)>3 else 'FAIL'} ({sep_f:.1f} sd)")
print(f"   (c) the grid separates on c (via p90)   : {'PASS' if abs(sep_c)>3 else 'FAIL'} ({sep_c:.1f} sd)")
if not(zero and abs(sep_f)>3 and abs(sep_c)>3):
    print("   -> UNVERIFIED, and that is not an acquittal.")
elif anchor(best)>6:
    print(f"\n   -> NO FIT: the best (f,c) still misses the anchors by {anchor(best):.1f} sd. The")
    print("      one-swap carrier model is wrong in SHAPE, not merely in size, and that is the finding.")
else:
    print(f"\n   -> THE ELEVATED GROUP IS {f_hat:.0%} OF PEOPLE, each elevated in {c_hat:.0%} of their")
    print(f"      domains by a single median->rare substitution.")
    print(f"      Held-out error {held:.1f} sd -> {'the model also predicts the tail' if held<6 else 'but the tail is NOT predicted, so the size is conditional on the shape being right'}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
