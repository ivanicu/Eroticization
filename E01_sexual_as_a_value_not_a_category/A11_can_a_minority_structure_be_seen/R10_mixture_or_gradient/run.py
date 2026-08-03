import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A11 R10 -- IS IT A MINORITY AT ALL, OR A GRADIENT ACROSS EVERYONE?

#95 licensed "a minority concentrating on rare options". #97d says the word MINORITY may be the wrong
one: the real high-surprisal group is elevated far more EVENLY across its blocks than any
discrete-carrier plant produces at the same mean, and that is what selecting the tail of a CONTINUUM
looks like. A distinct subgroup and a gradient are different objects with different mechanisms, and
nothing so far separates them.

  MIXTURE   a fraction f of people carry rare-option affinity; the rest carry none
  GRADIENT  everyone carries some, drawn from a continuous latent trait; nobody is special

Both can hit the same top-5% statistics. They differ in the SHAPE OF THE WHOLE DISTRIBUTION, so the
design matches them on TWO quantiles and tests them on the OTHERS -- an out-of-sample shape test
rather than another tail number.

ESTIMAND        the shape profile (p_q - p50)/(p75 - p25) at q in {90,95,99}, for each family tuned
                to reproduce the real p50 and p95, compared against the real profile.
IDENTIFICATION  identified iff the two families, once matched on p50 and p95, still DIFFER on the
                held-out quantiles; checked numerically, never asserted (#96a).
SCOPE           people entering >=5 of the 23 identified blocks.
UNCERTAINTY     #97c: the real arm's seed spread is structurally ZERO because the real matrix does not
                vary with the seed. The real point's uncertainty here is a BOOTSTRAP OVER PEOPLE.
WORLDS          mixture   real profile matches the mixture family -> a distinct subgroup exists and
                          "minority" in #95 is the right word
                gradient  real profile matches the gradient family -> there is no subgroup, only a
                          continuum, and #95 must be reworded as a tail of a gradient
                neither   real profile matches neither -> both models are wrong and the shape is the
                          finding
KILL            threshold-free: the real profile's distance to each family, in units of the real
                bootstrap spread; declared only when one distance is smaller by more than that spread.
POSITIVE CTRL   the two families must SEPARATE on the held-out quantiles after matching. If they do
                not, the design cannot tell them apart and no reading is licensed.
NEGATIVE CTRL   both families at zero intensity must reproduce the fixed-margin null.
NOISE FLOOR     5 seeds; 200 bootstrap resamples for the real point.
MULTIPLICITY    2 families x parameter grids x 5 seeds, published whole.
IMPOSSIBLE      distinguishing a gradient from a mixture with MANY components -- a mixture of enough
                components IS a gradient. The test separates "one distinct carrier group" from
                "smooth", not every conceivable latent structure.
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

def swaps_per_person(M,nsw_by_row,rng):
    """nsw_by_row[i] swaps for row i: a median-prevalence pick becomes a rare one."""
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in np.flatnonzero(nsw_by_row>0):
        d=0; want=int(nsw_by_row[i])
        for _ in range(8*max(want,1)):
            if d>=want: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; d+=1
    return Mw

def profile(S):
    q=np.percentile(S,[25,50,75,90,95,99]); iqr=max(q[2]-q[0],1e-9)
    return dict(p50=q[1],p90=q[3],p95=q[4],p99=q[5],
                s90=(q[3]-q[1])/iqr,s95=(q[4]-q[1])/iqr,s99=(q[5]-q[1])/iqr)

MIX=[(0.03,4),(0.05,3),(0.05,4),(0.08,3),(0.10,2),(0.15,2)]   # (fraction, swaps per elevated block)
GRAD=[0.4,0.7,1.0,1.5,2.2]                                     # scale of a continuous latent
rows=[]
for sd in range(1,6):
    rgc=np.random.default_rng(7700+sd)
    lat=rgc.exponential(1.0,size=len(ALLP))                    # the continuous trait, same all round
    arms=['real','cb']+[f'M{i}' for i in range(len(MIX))]+[f'G{i}' for i in range(len(GRAD))]
    S={w:{} for w in arms}
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        Mn=curveball(M,np.random.default_rng(7800+sd))
        pack={'real':M,'cb':Mn}
        for i,(f,k) in enumerate(MIX):
            car=(rgc.random(len(idx))<f) if False else (np.random.default_rng(7900+sd).random(len(ALLP))[idx]<f)
            nsw=np.where(car,k,0).astype(float)
            pack[f'M{i}']=swaps_per_person(Mn,nsw,np.random.default_rng(8000+sd))
        for i,g in enumerate(GRAD):
            nsw=np.floor(lat[idx]*g).astype(float)
            pack[f'G{i}']=swaps_per_person(Mn,nsw,np.random.default_rng(8000+sd))
        for w,Mw in pack.items():
            kk=Mw.sum(1); s=(Mw@ref)/np.maximum(kk,1)
            for j,gi in enumerate(idx):
                if kk[j]>0: S[w].setdefault(gi,[]).append(s[j])
    for w in arms:
        ks=[i for i,v in S[w].items() if len(v)>=5]
        vals=np.array([np.mean(S[w][i]) for i in ks])
        rows.append(dict(seed=sd,world=w,n=len(vals),**profile(vals)))
        if w=='real' and sd==1:
            bs=[]
            rb=np.random.default_rng(9999)
            for _ in range(200):
                r=profile(vals[rb.integers(0,len(vals),len(vals))])
                bs.append([r['p50'],r['p95'],r['s90'],r['s99']])
            B=np.array(bs)
            np.save(pathlib.Path(__file__).parent/'results'/'boot.npy',B)
    print(f"  seed {sd}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
B=np.load(OUT/'boot.npy'); bsd={'p50':B[:,0].std(),'p95':B[:,1].std(),'s90':B[:,2].std(),'s99':B[:,3].std()}
G=D.groupby('world')[['p50','p90','p95','p99','s90','s95','s99']].mean()
print("\n=== BOOTSTRAP SPREAD OF THE REAL POINT (200 resamples over PEOPLE, #97c) ===")
print("  " + "   ".join(f"{k} {v:.4f}" for k,v in bsd.items()))
print("\n=== EVERY ARM: matched-on quantiles (p50,p95) and held-out SHAPE (s90,s99) ===")
print(G.round(4).to_string())
r=G.loc['real']
def dist_match(w): return abs(G.loc[w,'p50']-r.p50)/bsd['p50']+abs(G.loc[w,'p95']-r.p95)/bsd['p95']
mixw=[f'M{i}' for i in range(len(MIX))]; gradw=[f'G{i}' for i in range(len(GRAD))]
bm=min(mixw,key=dist_match); bg=min(gradw,key=dist_match)
print(f"\n  best MIXTURE match : {bm} {MIX[int(bm[1:])]}   p50 {G.loc[bm,'p50']:.4f} p95 {G.loc[bm,'p95']:.4f}"
      f"  (match distance {dist_match(bm):.2f} bootstrap sd)")
print(f"  best GRADIENT match: {bg} scale {GRAD[int(bg[1:])]}   p50 {G.loc[bg,'p50']:.4f} p95 {G.loc[bg,'p95']:.4f}"
      f"  (match distance {dist_match(bg):.2f} bootstrap sd)")
sep=abs(G.loc[bm,'s99']-G.loc[bg,'s99'])/bsd['s99']
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) both families matched to within 3 bootstrap sd on (p50,p95): "
      f"{'PASS' if max(dist_match(bm),dist_match(bg))<3 else 'FAIL'}")
print(f"   (b) the two families SEPARATE on the held-out s99      : "
      f"{'PASS' if sep>2 else 'FAIL'} (gap {sep:.2f} bootstrap sd)")
if not(max(dist_match(bm),dist_match(bg))<3 and sep>2):
    print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    dm=abs(r.s99-G.loc[bm,'s99'])/bsd['s99']; dg=abs(r.s99-G.loc[bg,'s99'])/bsd['s99']
    print(f"\n   held-out s99: real {r.s99:.4f}   mixture {G.loc[bm,'s99']:.4f} ({dm:.2f} sd away)"
          f"   gradient {G.loc[bg,'s99']:.4f} ({dg:.2f} sd away)")
    if abs(dm-dg)<1:
        print("   -> INDISTINGUISHABLE on the held-out shape; both models fit and the fork stands open.")
    elif dg<dm:
        print("\n   -> GRADIENT. There is no distinct carrier group; rare-option affinity is a")
        print("      CONTINUUM across everyone, and #95's 'minority' must be reworded as the tail")
        print("      of a gradient.")
    else:
        print("\n   -> MIXTURE. A distinct subgroup carries the structure and #95's wording stands.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
