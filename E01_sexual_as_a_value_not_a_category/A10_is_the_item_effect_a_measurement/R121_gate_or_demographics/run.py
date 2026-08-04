"""
E01 A10 R06 -- WHICH OF THE TWO SENTENCES DOES #76b LICENSE?

#76 measured that projecting entry pattern AND demographics out of the cross-block person scores
shrinks C by 60.5%. It bundled them, so it cannot say whether the transfer is confounded with SURVEY
STRUCTURE (who answers which block) or with ORDINARY DEMOGRAPHICS (sex, age, personality). Those
license very different sentences and only one of them is a threat to the finding: the A02-era loader
always projected demographics out, so demographic variance was never part of what this project
claimed C was.

Four arms, identical cells, identical masks, one variable changed at a time:
  none  the raw scores                         (the #70/#72 specification)
  gate  entry indicators for the other 31 blocks only
  demo  sex, age, five personality scales, powerlessness only
  both  the #76 specification

ESTIMAND        C under each projection, and the DECOMPOSITION of the 60.5% shrinkage into a gate
                share and a demographic share, with the overlap named rather than assigned.
IDENTIFICATION  identified. Each projection is exact and its residual orthogonality is verified.
                The overlap between gate and demographics is NOT identified as a split and is
                reported as an interval, not apportioned.
SCOPE           the 23 blocks A09/R114 identified, Kc=8.
WORLDS          structural  gate alone reproduces most of the shrinkage -> C is partly an artefact
                            of the gated survey tree and #70/#72 must carry that caveat
                demographic demo alone reproduces most of it -> C was measuring sex/age/personality,
                            which this project had already decided to remove everywhere else
                shared      neither alone reproduces it -> the two are collinear here and the
                            honest report is the interval
KILL            threshold-free: shrinkage attributed to each factor is bounded below by its own
                marginal effect and above by the joint, and both bounds are published.
POSITIVE CTRL   each projection must remove variance and leave an orthogonal residual (verified).
NEGATIVE CTRL   the free permutation, retained in every arm.
NOISE FLOOR     3 masks.
MULTIPLICITY    23 blocks x 4 arms x 3 seeds, published whole.
IMPOSSIBLE      splitting the shared variance between gate and demographics -- they are collinear by
                construction, since who passes a gate depends on who they are. Reported as an
                interval with that reason, never as a point.
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
               'R03_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
ALLP=np.unique(np.concatenate([RAW[q]['ppl'] for q in RAW])); PM={p:i for i,p in enumerate(ALLP)}
BLKS=sorted(RAW)
E=np.zeros((len(ALLP),len(BLKS)))
for k,q in enumerate(BLKS): E[[PM[p] for p in RAW[q]['ppl']],k]=1.
NENT=E.sum(1)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
COV=pd.DataFrame({'male':pd.to_numeric(df.get('biomale'),errors='coerce'),
                  'age':df['age'].map(AGEMAP)}).reindex(ALLP)
for c in ['opennessvariable','neuroticismvariable','extroversionvariable',
          'consciensiousnessvariable','agreeablenessvariable','powerlessnessvariable']:
    if c in df.columns: COV[c]=pd.to_numeric(df[c],errors='coerce').reindex(ALLP).values
COV=COV.fillna(COV.median()).values
COV=(COV-COV.mean(0))/(COV.std(0)+1e-9)
print(f"targets {len(IDENT)}  people {len(ALLP)}  blocks {len(BLKS)}  "
      f"entry counts {NENT.min():.0f}-{NENT.max():.0f} (median {np.median(NENT):.0f})",flush=True)

MASK=0.15; SEEDS=[11,29,47]; KC=8

# ---- strata for the permutation ladder ----
cnt_stratum=NENT.astype(int)
pat=np.array([''.join(map(str,row.astype(int))) for row in E])
vc=pd.Series(pat).value_counts(); big=set(vc[vc>=20].index)
pat_stratum=np.array([p if p in big else '__unique__' for p in pat])
n_excl=int((pat_stratum=='__unique__').sum())
print(f"pattern strata >=20 people: {len(big)}   people with no usable stratum: {n_excl:,} "
      f"({100*n_excl/len(ALLP):.1f}%)",flush=True)

def strat_perm(idx,strata,rng):
    out=np.arange(len(idx))
    s=strata[idx]
    for v in np.unique(s):
        w=np.flatnonzero(s==v)
        if len(w)>1: out[w]=w[rng.permutation(len(w))]
    return out

def scores(target,mode):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:KC]*S[:KC]
    if mode=='none': return U,0.,np.nan
    oth=[k for k,q in enumerate(BLKS) if q!=target]
    parts=[np.ones((len(ALLP),1))]
    if mode in ('gate','both'): parts.append(E[:,oth])
    if mode in ('demo','both'): parts.append(COV)
    Dm=np.hstack(parts)
    b,*_=lstsq(Dm,U,rcond=None); Uf=U-Dm@b
    removed=1-Uf.var(0).sum()/max(U.var(0).sum(),1e-12)
    orth=float(np.abs((Dm[:,1:]-Dm[:,1:].mean(0)).T@(Uf-Uf.mean(0))).max()/len(ALLP))
    return Uf,removed,orth

def cval(target,U_all,seed,perm=None):
    M=RAW[target]['M']; n,m=M.shape; rows=np.array([PM[p] for p in RAW[target]['ppl']])
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm); I=(cm-gm)[None,:]
    T1=T-gm-I
    rm=np.nanmean(T1,axis=1); rm=np.where(np.isnan(rm),0.,rm); P=rm[:,None]
    Rres=T1-P
    U=U_all[rows]
    if perm is not None: U=U[perm]
    U=(U-U.mean(0))/(U.std(0)+1e-12)
    C=np.zeros_like(M)
    for j in range(m):
        k=obs[:,j]
        if k.sum()<50: continue
        b,*_=lstsq(np.c_[np.ones(k.sum()),U[k]],Rres[k,j],rcond=None)
        C[:,j]=np.c_[np.ones(n),U]@b
    base=np.mean((M[he]-gm)**2)
    f=lambda *p: 1.-np.mean((M[he]-np.clip(gm+sum(p),0,1)[he])**2)/base
    IB=np.broadcast_to(I,M.shape)
    return f(IB,P,C)-f(IB,P)

rows=[]; diag=[]
for i,t in enumerate(IDENT):
    for mode in ['none','gate','demo','both']:
        U,rem,orth=scores(t,mode)
        diag.append(dict(q=t,mode=mode,var_removed=rem,max_orth=orth))
        idx=np.array([PM[p] for p in RAW[t]['ppl']])
        for sd in SEEDS:
            rng=np.random.default_rng(2000+sd)
            rows.append(dict(q=t,mode=mode,arm='real',seed=sd,C=cval(t,U,sd,None)))
            rows.append(dict(q=t,mode=mode,arm='perm',seed=sd,
                             C=cval(t,U,sd,rng.permutation(len(idx)))))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows); G=pd.DataFrame(diag)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R121_gate_or_demographics/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== POSITIVE CONTROL: every projection exact and non-vacuous ===")
print(G.groupby('mode').agg(var_removed=('var_removed','median'),max_orth=('max_orth','max')).round(6).to_string())
print(f"  all residuals orthogonal: {'PASS' if G.max_orth.max(skipna=True)<1e-8 else 'FAIL'}")

print("\n=== C BY PROJECTION (mean over 23 blocks x 3 seeds) ===")
piv=D.groupby(['mode','arm']).C.mean().unstack('arm')
base=float(piv.loc['none','real'])
piv['shrinkage']=1-piv.real/base
print(piv.reindex(['none','gate','demo','both']).round(5).to_string())

print("\n=== PER-BLOCK survival against the free permutation ===")
out=[]
for m in ['none','gate','demo','both']:
    d=D[D.mode==m]; r=d[d.arm=='real'].groupby('q').C; nn=d[d.arm=='perm'].groupby('q').C
    gap=r.mean()-nn.mean(); sp=np.sqrt(r.std()**2+nn.std()**2)
    out.append(dict(mode=m,C=r.mean().median(),beats=int((gap>2*sp).sum()),n=len(gap)))
S=pd.DataFrame(out); print(S.to_string(index=False))

g_=float(piv.loc['gate','shrinkage']); d_=float(piv.loc['demo','shrinkage']); b_=float(piv.loc['both','shrinkage'])
print("\n  CONDITIONAL KILL -- gates first")
ok=G.max_orth.max(skipna=True)<1e-8 and G[G['mode']!='none'].var_removed.min()>0.01
print(f"   (a) projections exact and non-vacuous : {'PASS' if ok else 'FAIL'}")
if not ok: print("   -> UNVERIFIED.")
else:
    print(f"\n   shrinkage of C:  gate alone {g_:.1%}   demographics alone {d_:.1%}   both {b_:.1%}")
    print(f"   overlap (collinear, NOT apportioned): {g_+d_-b_:+.1%}")
    print(f"   gate's share is bounded by [{max(b_-d_,0):.1%}, {min(g_,b_):.1%}] of C")
    print(f"   demographics' share is bounded by [{max(b_-g_,0):.1%}, {min(d_,b_):.1%}] of C")
    if g_>d_*1.5:
        print("\n   -> STRUCTURAL. The gated survey tree drives most of the shrinkage, so #70/#72")
        print("      must carry that caveat: part of what looked like a shared readout is shared")
        print("      exposure to the same entry conditions.")
    elif d_>g_*1.5:
        print("\n   -> DEMOGRAPHIC. C was partly sex/age/personality, which this project removes")
        print("      everywhere else. Not a threat to the finding, but the A09/A10 scores should")
        print("      have carried the same projection the A02 loader always did.")
    else:
        print("\n   -> SHARED and collinear: who passes a gate depends on who they are, so no split")
        print("      is identified. The interval above is the honest report.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
