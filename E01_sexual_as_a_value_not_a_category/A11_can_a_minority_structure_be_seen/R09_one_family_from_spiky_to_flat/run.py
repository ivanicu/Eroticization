import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R09 -- ONE FAMILY WITH A SINGLE PARAMETER, BECAUSE "SPIKY vs FLAT" WAS TWO LADDERS THAT
COULD NOT MEET.

R08 rebased both plant families on the fixed-margin null and gate (b) still refused -- but for a
reason worth more than the gate: THE SPIKY FAMILY SATURATES IN MEAN. Intensity inside 2 blocks tops
out at mean S = 0.976 (S2 0.9659 -> S10 0.9759) because there are only 3 rare options per block and a
carrier who already holds all of them cannot become more surprising. The real top-5% sits at 1.0277,
ABOVE that structural ceiling. Twenty-third mis-specified design element: I parameterised the spiky
arm on the axis that saturates.

The axis that does not saturate is BREADTH -- how many of a carrier's blocks are elevated. And breadth
turns two ladders that could never meet into ONE family with a single parameter:

    b = number of elevated blocks per carrier,  from 1 (maximally spiky) to all (flat)

Every point on that family has the same per-block intensity; only the concentration changes. So sd
can finally be read at matched mean, because the family sweeps mean continuously.

ESTIMAND        the (mean, sd) curve of the top-5% as a function of BREADTH b, and where the real
                point falls on it -- read as an implied b.
IDENTIFICATION  identified iff the curve is monotone in b and the real mean lies inside its range;
                both CHECKED numerically (#96a's lesson: never assert a condition in a label).
SCOPE           people entering >=5 of the 23 identified blocks.
WORLDS          narrow   real point implies small b -> a minority intense about FEW domains
                broad    real point implies large b -> elevated across most domains, which is the
                         response-quality reading
KILL            threshold-free: the implied b, with the spread of the real point carried through.
POSITIVE CTRL   the curve must be monotone in b, and b=all must reproduce R08's flat arm.
NEGATIVE CTRL   b=0 must reproduce the fixed-margin null exactly.
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
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/R03_fixed_margin_null/results/grid.csv')
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
BS=[0,1,2,3,5,8,13,23]; INTENSITY=3
rows=[]
for sd in range(1,6):
    rgc=np.random.default_rng(6600+sd)
    carrier=rgc.random(len(ALLP))<0.05; ci=np.flatnonzero(carrier)
    perm={i:rgc.permutation(len(IDENT)) for i in ci}       # each carrier's block priority order
    arms=['real','cb']+[f'b{b}' for b in BS]
    S={w:{} for w in arms}
    for bi,t in enumerate(IDENT):
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        ref=-np.log(np.clip(M.mean(0),1e-4,1.)); rg=np.random.default_rng(6700+sd)
        Mn=curveball(M,rg)                                  # shared origin
        pack={'real':M,'cb':Mn}
        loc=np.flatnonzero(carrier[idx])
        for b in BS:
            w=np.array([j for j in loc if bi in perm[idx[j]][:b]],dtype=int) if b>0 else np.array([],dtype=int)
            pack[f'b{b}']=swaps(Mn,w,INTENSITY,np.random.default_rng(6800+sd)) if b>0 else Mn
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
print("\n=== ONE FAMILY: (mean, sd) OF THE TOP-5% vs BREADTH b (blocks elevated per carrier) ===")
print(G.reindex(['cb','b0','b1','b2','b3','b5','b8','b13','b23','real']).round(4).to_string())
bm=np.array([G.loc[f'b{b}','mean'] for b in BS]); bs=np.array([G.loc[f'b{b}','sd'] for b in BS])
rm_,rs_,re_=G.loc['real','mean'],G.loc['real','sd'],G.loc['real','sd_e']
mono=all(bm[i]<=bm[i+1]+1e-9 for i in range(len(bm)-1))
zero=abs(G.loc['b0','mean']-G.loc['cb','mean'])<1e-9
print("\n  CONDITIONAL KILL -- gates first")
print(f"   (a) b=0 reproduces the null exactly : {'PASS' if zero else 'FAIL'}")
print(f"   (b) mean is monotone in b           : {'PASS' if mono else 'FAIL'}  {bm.round(4).tolist()}")
inr=bm.min()<=rm_<=bm.max()
print(f"   (c) the real mean is inside the range: {'PASS' if inr else 'FAIL'} "
      f"(real {rm_:.4f}, range [{bm.min():.4f}, {bm.max():.4f}])")
if not(zero and mono and inr):
    print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    b_hat=float(np.interp(rm_,bm,np.array(BS,dtype=float)))
    sd_at=float(np.interp(rm_,bm,bs))
    print(f"\n   implied breadth b_hat = {b_hat:.2f} of {len(IDENT)} blocks "
          f"({100*b_hat/len(IDENT):.0f}% of a carrier's domains)")
    print(f"   sd predicted by the family at that mean: {sd_at:.4f}   real sd {rs_:.4f}  "
          f"(2x real spread {2*re_:.4f})")
    fit = abs(rs_-sd_at)<=2*re_
    print(f"   the real point LIES ON the family: {'YES' if fit else 'NO -- real sd is off the curve'}")
    if not fit:
        print("   -> the real minority is not a uniform-intensity structure at ANY breadth; its sd")
        print("      does not match what this family produces at its own mean.")
    elif b_hat<=3:
        print("\n   -> NARROW. The observed level is reproduced by a minority elevated in only")
        print(f"      ~{b_hat:.1f} of {len(IDENT)} domains -- concentrated, which is the attachment shape.")
    elif b_hat>=0.6*len(IDENT):
        print("\n   -> BROAD. The observed level requires elevation across most domains, which is the")
        print("      response-quality reading rather than a specific attachment.")
    else:
        print(f"\n   -> INTERMEDIATE at b_hat={b_hat:.1f}: neither a narrow attachment nor a uniform")
        print("      elevation; the minority is elevated across a substantial minority of domains.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
