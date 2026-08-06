import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R07 -- READING SD AS A FUNCTION OF MEAN, SO THE MEAN CONFOUND CANNOT DECIDE THE ANSWER.

#96 built the right separator -- careless is FLAT across blocks, an attachment is SPIKY -- and could
not read it, because the two plants differed in mean surprisal by 36% and sd scales with mean. Raw sd
said spiky; sd/mean said flat. Choosing between two defensible normalisations is not a measurement.

The fix is not to hit a matched mean by hand. It is to SWEEP both arms and compare sd AS A FUNCTION
OF MEAN, then place the real point on that plane. If the flat and spiky curves are separated in sd at
every mean they share, the real point falls on one of them and the mean confound is gone by
construction.

  flat arm   each carrier gets 1 swap in a FRACTION f of the blocks they enter    f in {0.25,0.5,1.0}
  spiky arm  each carrier gets k swaps in exactly 2 blocks                        k in {2,4,6,10}

Every arm is summarised on its TOP 5% BY MEAN S -- the same selection applied to the real data --
because selecting on mean induces its own correlation with sd, and comparing a selected group to an
unselected one is how #96's confound would come back in a different costume.

ESTIMAND        the (mean, sd) curve of block-level surprisal for each plant family, and where the
                real top-5% falls relative to them.
IDENTIFICATION  identified iff the curves separate in sd over the mean range containing the real
                point; that condition is CHECKED numerically, not asserted in a label (#96a).
SCOPE           people entering >=5 of the 23 identified blocks.
WORLDS          attachment  real point lies on the spiky curve
                careless    real point lies on the flat curve
                neither     curves do not separate where the real point is -> unanswerable here
KILL            threshold-free: the real point's distance to each curve, in units of its own spread.
POSITIVE CTRL   the curves must separate somewhere; and each family must be monotone in its own
                intensity parameter, or the sweep is not a sweep.
NEGATIVE CTRL   fixed-margin null, whose top-5% is high by chance -- the origin of the plane.
NOISE FLOOR     5 seeds.
IMPOSSIBLE      naming what the attachment is about -- unchanged (#28/#44).
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
def swaps(M,who,nsw,rng):
    Mw=M.copy(); p=M.mean(0); o=np.argsort(p); rare=o[:3]; med=o[len(o)//2:]
    for i in who:
        d=0
        for _ in range(8*max(nsw,1)):
            if d>=nsw: break
            c=med[rng.integers(len(med))]; r=rare[rng.integers(len(rare))]
            if Mw[i,c]==1 and Mw[i,r]==0: Mw[i,c]=0.; Mw[i,r]=1.; d+=1
    return Mw
FLAT=[0.25,0.5,1.0]; SPIKY=[2,4,6,10]
rows=[]
for sd in range(1,6):
    rgc=np.random.default_rng(5500+sd)
    carrier=rgc.random(len(ALLP))<0.05
    ci=np.flatnonzero(carrier)
    home={i:set(rgc.choice(len(IDENT),2,replace=False).tolist()) for i in ci}
    fl_on={f:{i:set(np.flatnonzero(rgc.random(len(IDENT))<f).tolist()) for i in ci} for f in FLAT}
    arms=['real','cb']+[f'F{f}' for f in FLAT]+[f'S{k}' for k in SPIKY]
    S={w:{} for w in arms}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.)); rg=np.random.default_rng(5600+sd)
        pack={'real':M,'cb':curveball(M,rg)}
        loc=np.flatnonzero(carrier[idx])
        for f in FLAT:
            w=np.array([j for j in loc if bi in fl_on[f][idx[j]]],dtype=int)
            pack[f'F{f}']=swaps(M,w,1,np.random.default_rng(5700+sd))
        for k in SPIKY:
            w=np.array([j for j in loc if bi in home[idx[j]]],dtype=int)
            pack[f'S{k}']=swaps(M,w,k,np.random.default_rng(5700+sd))
        for w,Mw in pack.items():
            kk=Mw.sum(1); s=(Mw@ref)/np.maximum(kk,1)
            for j,gi in enumerate(idx):
                if kk[j]>0: S[w].setdefault(gi,[]).append(s[j])
    for w in arms:
        ks=[i for i,v in S[w].items() if len(v)>=5]
        mean=np.array([np.mean(S[w][i]) for i in ks]); sdv=np.array([np.std(S[w][i]) for i in ks])
        hi=mean>=np.percentile(mean,95)
        rows.append(dict(seed=sd,world=w,mean=float(mean[hi].mean()),sd=float(sdv[hi].mean())))
    print(f"  seed {sd}",flush=True)
D=pd.DataFrame(rows); OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
G=D.groupby('world').agg(mean=('mean','mean'),sd=('sd','mean'),
                         mean_e=('mean','std'),sd_e=('sd','std'))
print("\n=== (mean, sd) OF THE TOP-5% BY MEAN SURPRISAL, EVERY ARM ===")
print(G.reindex(['cb','real']+[f'F{f}' for f in FLAT]+[f'S{k}' for k in SPIKY]).round(4).to_string())
fm=np.array([G.loc[f'F{f}','mean'] for f in FLAT]); fs=np.array([G.loc[f'F{f}','sd'] for f in FLAT])
sm=np.array([G.loc[f'S{k}','mean'] for k in SPIKY]); ss=np.array([G.loc[f'S{k}','sd'] for k in SPIKY])
rm_,rs_=G.loc['real','mean'],G.loc['real','sd']; re_=G.loc['real','sd_e']
lo,hi=max(fm.min(),sm.min()),min(fm.max(),sm.max())
print(f"\n  overlapping mean range of the two families: [{lo:.4f}, {hi:.4f}]   real mean {rm_:.4f}")
mono_f=all(fs[i]<=fs[i+1]+1e-9 for i in range(len(fs)-1)) or all(fs[i]>=fs[i+1]-1e-9 for i in range(len(fs)-1))
mono_s=all(ss[i]<=ss[i+1]+1e-9 for i in range(len(ss)-1)) or all(ss[i]>=ss[i+1]-1e-9 for i in range(len(ss)-1))
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) each family monotone in its own intensity : flat {'PASS' if mono_f else 'FAIL'}, "
      f"spiky {'PASS' if mono_s else 'FAIL'}")
inrange=lo<=rm_<=hi
print(f"   (b) the real mean lies in the overlapping range: {'PASS' if inrange else 'FAIL'}")
if inrange:
    fp=np.interp(rm_,fm[np.argsort(fm)],fs[np.argsort(fm)])
    sp=np.interp(rm_,sm[np.argsort(sm)],ss[np.argsort(sm)])
    gap=abs(sp-fp)
    print(f"   (c) the curves SEPARATE at the real mean      : "
          f"{'PASS' if gap>2*re_ else 'FAIL'} (flat {fp:.4f} vs spiky {sp:.4f}, gap {gap:.4f}, "
          f"2x real spread {2*re_:.4f})")
    if mono_f and mono_s and gap>2*re_:
        df_,ds_=abs(rs_-fp),abs(rs_-sp)
        print(f"\n   real sd {rs_:.4f}   distance to flat {df_:.4f}   to spiky {ds_:.4f}")
        print("   ->","ATTACHMENT-LIKE (spiky)" if ds_<df_ else "CARELESS-LIKE (flat)",
              f"-- and by {abs(df_-ds_)/max(2*re_,1e-9):.1f}x the real spread")
    else:
        print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"   -> UNVERIFIED: the real mean {rm_:.4f} is outside where the families overlap, so no")
    print("      interpolated comparison is licensed. Extend the ladders until it is inside.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
