import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R02 -- THE ONE ADMISSIBLE COMPARISON, WITH ITS TWO KNOWN FAULTS FIXED.

#92 left exactly one admissible quantile (p99, the only one where the positive control fires) and two
fixable faults: a parametric null that is mis-specified at the median because respondents pick a
roughly stable NUMBER of options per block, and three seeds giving a 2x spread (0.96) wider than the
effect it had to judge (0.60).

Both fixes are mechanical. The null becomes FIXED-MARGIN (curveball) -- row sums exact, so the
response-format under-dispersion is reproduced rather than mistaken for signal -- and the seed count
goes to 8.

ESTIMAND        per-person misfit T at p99, real vs a row-sum-preserving null.
IDENTIFICATION  identified; curveball preserves each person's pick count and each option's base rate
                EXACTLY (asserted per draw) and destroys only person x option structure.
SCOPE           people with >=40 held-out cells across the 23 identified blocks.
WORLDS          minority-present  real p99 exceeds the curveball null's beyond 2x spread
                none              it does not, and now the null is admissible because #92d showed
                                  the instrument fires at p99 on a 5%-carrier world
KILL            threshold-free, and ONLY p99 is interpreted -- the quantiles where the instrument is
                blind are printed but explicitly not read.
POSITIVE CTRL   the 5%-carrier world again, against its own curveball null, at 4 seeds.
NEGATIVE CTRL   the curveball null IS it, and #92c's median mismatch is the check: real and null p50
                must now AGREE, or the new null is mis-specified too.
NOISE FLOOR     8 seeds real/null, 4 seeds control.
MULTIPLICITY    quantiles printed whole; one interpreted, declared in advance.
IMPOSSIBLE      attributing an inflated T to a cause -- unchanged from #92.
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
MASK=0.15; QS=[50,75,90,95,99]
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
def shrunk_rows(T1,obs):
    rmn=np.nanmean(T1,axis=1); rmn=np.where(np.isnan(rmn),0.,rmn)
    k=np.maximum(obs.sum(1),1)
    vw=np.nanmean(np.where(obs,(T1-rmn[:,None])**2,np.nan),axis=1)
    vw=np.where(np.isnan(vw),np.nanmean(vw),vw); s2w=vw/k
    s2t=max(np.var(rmn)-np.mean(s2w),1e-9)
    return rmn*(s2t/(s2t+s2w))
def per_person(M,seed):
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    P=shrunk_rows(T-gm-I,obs)[:,None]
    ph=np.clip(gm+I+P,0.02,0.98)
    z2=(M-ph)**2/(ph*(1-ph))
    return np.where(he,z2,0.).sum(1), he.sum(1)
def plant_sparse(M,dens,sc,rng):
    n,m=M.shape; F=rng.normal(size=(n,5)); L=rng.normal(size=(5,m))*sc
    dev=(F@L)*(rng.random(n)<dens)[:,None]/np.sqrt(dens)
    return (rng.random((n,m))<np.clip(M.mean(0)[None,:]+dev,0.02,0.98)).astype(float)
rows=[]
for sd in range(1,9):
    acc={w:[np.zeros(len(ALLP)),np.zeros(len(ALLP))] for w in ['real','cb']}
    ctrl=sd<=4
    if ctrl: acc.update({w:[np.zeros(len(ALLP)),np.zeros(len(ALLP))] for w in ['sp5','sp5cb']})
    for t in IDENT:
        M=RAW[t]['M']; idx=np.array([PM[p] for p in RAW[t]['ppl']])
        rg=np.random.default_rng(7700+sd)
        Mcb=curveball(M,rg)
        assert np.allclose(Mcb.sum(0),M.sum(0)) and np.allclose(Mcb.sum(1),M.sum(1)),"margins"
        pack={'real':M,'cb':Mcb}
        if ctrl:
            Msp=plant_sparse(M,0.05,0.50,np.random.default_rng(8800+sd))
            pack.update({'sp5':Msp,'sp5cb':curveball(Msp,np.random.default_rng(8801+sd))})
        for w,Mw in pack.items():
            s,c=per_person(Mw,sd); acc[w][0][idx]+=s; acc[w][1][idx]+=c
    ok=acc['real'][1]>=40
    for w in acc:
        Tv=acc[w][0][ok]/np.maximum(acc[w][1][ok],1)
        for qv in QS: rows.append(dict(seed=sd,world=w,qq=qv,T=float(np.percentile(Tv,qv))))
        rows.append(dict(seed=sd,world=w,qq=-2,T=float((Tv>2.0).mean())))
    print(f"  seed {sd}  n={int(ok.sum()):,}",flush=True)
D=pd.DataFrame(rows)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
def cmp(a,b,qv):
    r=D[(D.world==a)&(D.qq==qv)].T; n=D[(D.world==b)&(D.qq==qv)].T
    return r.mean(),n.mean(),r.mean()-n.mean(),2*np.sqrt(r.std()**2+n.std()**2)
lab={-2:'share T>2',50:'p50',75:'p75',90:'p90',95:'p95',99:'p99'}
print("\n=== REAL vs FIXED-MARGIN NULL (8 seeds) -- only p99 is interpreted ===")
for qv in [-2,50,75,90,95,99]:
    r,n,d,s=cmp('real','cb',qv)
    print(f"  {lab[qv]:10s} real {r:.4f}  null {n:.4f}  diff {d:+.4f}  2xspread {s:.4f}"
          f"  {'RESOLVABLE' if abs(d)>s else 'no'}{'   <- INTERPRETED' if qv==99 else ''}")
print("\n=== POSITIVE CONTROL: 5% carriers vs its own fixed-margin null (4 seeds) ===")
for qv in [90,95,99]:
    r,n,d,s=cmp('sp5','sp5cb',qv)
    print(f"  {lab[qv]:10s} {r:.4f} vs {n:.4f}  diff {d:+.4f}  2xspread {s:.4f}  "
          f"{'SEES IT' if d>s else 'blind'}")
r50,n50,d50,s50=cmp('real','cb',50)
r99,n99,d99,s99=cmp('real','cb',99)
cr,cn,cd,cs=cmp('sp5','sp5cb',99)
print("\n  CONDITIONAL KILL -- gates first")
g1=abs(d50)<=s50
g2=cd>cs
print(f"   (a) the new null reproduces the median (#92c's fault fixed): {'PASS' if g1 else 'FAIL'} "
      f"(diff {d50:+.4f}, 2x spread {s50:.4f})")
print(f"   (b) positive control still fires at p99                    : {'PASS' if g2 else 'FAIL'} "
      f"(diff {cd:+.4f})")
if not(g1 and g2): print("   -> UNVERIFIED, and that is not an acquittal.")
else:
    print(f"\n   p99: real {r99:.4f}  null {n99:.4f}  diff {d99:+.4f}  2x spread {s99:.4f}")
    print(f"   (#92 had diff +0.5958 against a 2x spread of 0.9558)")
    if d99>s99:
        print("\n   -> A MINORITY STRUCTURE IS PRESENT, at the one quantile where this instrument is")
        print("      validated and against a null that reproduces the response format.")
    else:
        print("\n   -> NO MINORITY STRUCTURE at p99, and this null IS admissible: the instrument fires")
        print("      on a 5%-carrier world and the null reproduces the median. The strongest negative")
        print("      result this project has produced, because for once every precondition holds.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
