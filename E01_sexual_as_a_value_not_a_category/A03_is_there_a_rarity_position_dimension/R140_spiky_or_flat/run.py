import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R06 -- WHICH OF THE TWO THINGS #95 COULD NOT TELL APART IS IT?

#95 licensed a minority whose picks concentrate on rare options, and named the limit: surprisal
cannot distinguish an idiosyncratic erotic attachment from a careless responder. But those two make
DIFFERENT predictions about SHAPE ACROSS BLOCKS, and that is testable here:

  careless    picks poorly everywhere -> high mean S, LOW between-block variance. FLAT.
  attachment  intense about a few specific things -> high S in a FEW blocks, ordinary elsewhere. SPIKY.

So the separator is the within-person spread of block-level surprisal, conditioned on the mean.

ESTIMAND        for the high-S minority: the between-block sd of their block-level S, relative to
                what each of two calibrated plants produces at the SAME mean S.
IDENTIFICATION  identified relative to the two plants; a bare sd means nothing, an sd read against a
                careless world and an attachment world at matched mean means something.
SCOPE           people entering >=5 of the 23 identified blocks, so a between-block sd exists.
WORLDS          spiky  real high-S people look like the ATTACHMENT plant -> the minority is
                       content-specific, which is the eroticization operator's own signature
                flat   they look like the CARELESS plant -> #95's signal is response quality and
                       must be reported as such
                neither -> the two plants do not separate and the question is unanswerable here
KILL            threshold-free: real sd placed between the two plants' sds, each with its own spread.
POSITIVE CTRL   the two plants must SEPARATE from each other at matched mean S. If they do not, no
                reading of the real data is licensed -- this is the gate that matters.
NEGATIVE CTRL   fixed-margin null: its high-S people are high by chance, so their sd is the chance
                reference.
NOISE FLOOR     6 seeds.
IMPOSSIBLE      naming WHAT the attachment is about -- that needs option semantics, and #28/#44
                measured that string proxies here lose to their own shams.
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
def swap_to_rare(M,who,nsw,rng):
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in who:
        d=0
        for _ in range(6*max(nsw,1)):
            if d>=nsw: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; d+=1
    return Mw
rows=[]
for sd in range(1,7):
    rgc=np.random.default_rng(4400+sd)
    carrier=rgc.random(len(ALLP))<0.05
    # ATTACHMENT: each carrier is intense in exactly 2 randomly chosen blocks (many swaps there)
    home={i:set(rgc.choice(len(IDENT),2,replace=False).tolist()) for i in np.flatnonzero(carrier)}
    arms=['real','cb','flat','spiky']
    S={w:{} for w in arms}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.))
        loc=np.flatnonzero(carrier[idx]); rg=np.random.default_rng(4500+sd)
        flat_loc=loc                                             # careless: every block, 1 swap
        spiky_loc=np.array([j for j in loc if bi in home.get(idx[j],set())],dtype=int)
        pack={'real':M,'cb':curveball(M,rg),
              'flat':swap_to_rare(M,flat_loc,1,np.random.default_rng(4600+sd)),
              'spiky':swap_to_rare(M,spiky_loc,6,np.random.default_rng(4600+sd))}
        for w,Mw in pack.items():
            k=Mw.sum(1); s=(Mw@ref)/np.maximum(k,1)
            for j,gi in enumerate(idx):
                if k[j]>0: S[w].setdefault(gi,[]).append(s[j])
    for w in arms:
        keep_=[i for i,v in S[w].items() if len(v)>=5]
        mean=np.array([np.mean(S[w][i]) for i in keep_])
        sdv =np.array([np.std(S[w][i]) for i in keep_])
        car =np.array([carrier[i] for i in keep_])
        hi=mean>=np.percentile(mean,95)
        rows.append(dict(seed=sd,world=w,n=len(keep_),
                         mean_hi=float(mean[hi].mean()),sd_hi=float(sdv[hi].mean()),
                         mean_all=float(mean.mean()),sd_all=float(sdv.mean()),
                         mean_car=float(mean[car].mean()) if car.any() else np.nan,
                         sd_car=float(sdv[car].mean()) if car.any() else np.nan))
    print(f"  seed {sd}  n={rows[-1]['n']:,}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby('world')[['mean_hi','sd_hi','mean_car','sd_car','mean_all','sd_all']].agg(['mean','std'])
print("\n=== BETWEEN-BLOCK SD OF SURPRISAL, top-5% by mean S ===")
print(D.groupby('world')[['mean_hi','sd_hi','mean_car','sd_car']].mean().round(4).to_string())
def g(w,c): return D[D.world==w][c]
print("\n  CONDITIONAL KILL -- gates first")
fm,fs=g('flat','mean_car').mean(),g('flat','sd_car').mean()
sm,ss=g('spiky','mean_car').mean(),g('spiky','sd_car').mean()
sep=abs(ss-fs); spread=2*np.sqrt(g('flat','sd_car').std()**2+g('spiky','sd_car').std()**2)
print(f"   carriers under FLAT  plant: mean S {fm:.4f}  between-block sd {fs:.4f}")
print(f"   carriers under SPIKY plant: mean S {sm:.4f}  between-block sd {ss:.4f}")
print(f"   (a) the two plants SEPARATE on sd at comparable mean: "
      f"{'PASS' if sep>spread else 'FAIL'} (gap {sep:.4f}, 2x spread {spread:.4f})")
if sep<=spread:
    print("   -> UNVERIFIED: the plants do not separate, so no reading of the real data is licensed.")
else:
    rh,rs=g('real','mean_hi').mean(),g('real','sd_hi').mean()
    nh,ns=g('cb','mean_hi').mean(),g('cb','sd_hi').mean()
    print(f"\n   real top-5%: mean S {rh:.4f}  sd {rs:.4f}")
    print(f"   null top-5%: mean S {nh:.4f}  sd {ns:.4f}   (high by chance -- the reference)")
    pos=(rs-ns)/max(ss-fs,1e-9)
    print(f"   real sd sits {rs-ns:+.4f} above the chance reference; the spiky-minus-flat gap is "
          f"{ss-fs:+.4f}")
    if rs>ns+0.5*spread:
        print("\n   -> SPIKY. The real high-surprisal minority is concentrated in a FEW blocks, which")
        print("      is the attachment signature and not the carelessness one.")
    elif rs<ns-0.5*spread:
        print("\n   -> FLATTER THAN CHANCE. Consistent with a uniform response-quality effect.")
    else:
        print("\n   -> INDISTINGUISHABLE from chance on this axis; #95's limit stands unresolved.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
