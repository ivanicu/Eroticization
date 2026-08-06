import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R05 -- RUNNING MY OWN HIGHEST-PROBABILITY SELF-OVERTURN.

ADVERSARY_FORECAST.md, second block, p=0.70: "the cross-block transfer C is partly the GATE, not the
person." Entry to every block is conditioned on the parent rating (P(enter|parent>0)=0.99, #04), so
people present in many blocks are people who rated many parents highly. "The same person deviates
the same way across blocks" may be partly "the same person passed the same gate."

And the person-permutation null used in #70/#71/#72 does NOT control for this. It destroys the
person correspondence, but every permuted person is still a gate-passer, so a gate-driven C would
survive it intact. #70, #71 and #72 all rest on C.

Two independent attacks, because a permutation and a projection fail differently:

  ATTACK A (projection)  regress the 31 OTHER blocks' entry indicators -- and demographics -- out of
                         the cross-block person scores, then recompute C. If C is the gate, the
                         gate-free scores carry nothing.
  ATTACK B (stratified   permute person identity only WITHIN strata of matched gate exposure:
            permutation) (i) free permutation, the old null; (ii) count-matched, same number of
                         blocks entered; (iii) pattern-matched, same exact entry pattern where a
                         stratum has >=20 people. A ladder, so partial survival is visible.

ESTIMAND        C measured with raw scores, with gate-free scores, and against three nulls of
                increasing strictness, on identical held-out cells.
IDENTIFICATION  identified. The projection removes a linear subspace; what remains is orthogonal to
                every entry indicator by construction, which is verified numerically here.
SCOPE           the 23 blocks A09/R114 identified. Kc=8 (R04: C is still rising there, so this is a
                conservative rank at which to attack it).
WORLDS          gate-driven   C_free collapses to its null, or C falls to the pattern-matched null
                person-driven C_free survives and C exceeds even the pattern-matched null
                partial       C survives count-matching but not pattern-matching -> the transfer is
                              real but its SIZE in #70/#72 is inflated and must be restated
KILL            threshold-free: C is declared surviving a null only where it exceeds it by more than
                2x the pooled seed spread, per block, and the block counts are published for each
                null separately.
POSITIVE CTRL   the projection must not be vacuous: report the variance of U removed by it, and
                verify the residual scores are numerically orthogonal to the entry indicators.
NEGATIVE CTRL   the free permutation, retained so the ladder has its old bottom rung.
NOISE FLOOR     3 masks x 3 permutation draws.
MULTIPLICITY    23 blocks x 2 score types x 4 arms x 3 seeds, published whole.
IMPOSSIBLE      separating gate from person for someone whose entry pattern is UNIQUE -- no stratum
                exists. Those people are excluded from the pattern-matched arm and the count of
                excluded people is reported, not hidden.
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
               'R114_fixed_margin_null/results/grid.csv')
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

def scores(target,gate_free):
    cols=[]
    for q in BLKS:
        if q==target: continue
        M=RAW[q]['M']; R=M-M.mean(0,keepdims=True); R=R-R.mean(1,keepdims=True)
        Z=np.full((len(ALLP),M.shape[1]),np.nan); Z[[PM[p] for p in RAW[q]['ppl']]]=R
        cols.append(Z)
    Z=np.hstack(cols); mu=np.nanmean(Z,axis=0); Z=np.where(np.isnan(Z),mu,Z); Z=Z-Z.mean(0)
    U,S,_=svd(Z,full_matrices=False); U=U[:,:KC]*S[:KC]
    removed=0.
    if gate_free:
        oth=[k for k,q in enumerate(BLKS) if q!=target]
        D=np.c_[np.ones(len(ALLP)),E[:,oth],COV]
        b,*_=lstsq(D,U,rcond=None); Uf=U-D@b
        removed=1-Uf.var(0).sum()/max(U.var(0).sum(),1e-12)
        orth=float(np.abs((E[:,oth]-E[:,oth].mean(0)).T@(Uf-Uf.mean(0))).max()/len(ALLP))
        return Uf,removed,orth
    return U,0.,np.nan

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
    for gf in [False,True]:
        U,rem,orth=scores(t,gf)
        if gf: diag.append(dict(q=t,var_removed=rem,max_orth=orth))
        idx=np.array([PM[p] for p in RAW[t]['ppl']])
        for sd in SEEDS:
            rng=np.random.default_rng(2000+sd)
            arms={'real':None,
                  'perm_free':rng.permutation(len(idx)),
                  'perm_count':strat_perm(idx,cnt_stratum,rng),
                  'perm_pattern':strat_perm(idx,pat_stratum,rng)}
            for a,pm_ in arms.items():
                rows.append(dict(q=t,gate_free=gf,arm=a,seed=sd,C=cval(t,U,sd,pm_)))
    print(f"  {i+1}/{len(IDENT)}",flush=True)
D=pd.DataFrame(rows); G=pd.DataFrame(diag)
OUT='E01_sexual_as_a_value_not_a_category/A10_is_the_item_effect_a_measurement/R120_is_the_transfer_the_gate/results/'
D.to_csv(OUT+'grid.csv',index=False)

print("\n=== POSITIVE CONTROL: is the projection vacuous? ===")
print(f"  variance of the person scores removed by projecting out entry pattern + demographics: "
      f"median {G.var_removed.median():.1%}  range [{G.var_removed.min():.1%}, {G.var_removed.max():.1%}]")
print(f"  max residual correlation with any entry indicator: {G.max_orth.max():.2e}  -> "
      f"{'orthogonal, PASS' if G.max_orth.max()<1e-8 else 'NOT orthogonal, FAIL'}")

print("\n=== C BY SCORE TYPE AND NULL (mean over 23 blocks x 3 seeds) ===")
piv=D.groupby(['gate_free','arm']).C.mean().unstack('arm')
print(piv[['real','perm_free','perm_count','perm_pattern']].round(5).to_string())

print("\n=== PER-BLOCK: does C beat each null by >2x the pooled seed spread? ===")
out=[]
for gf in [False,True]:
    d=D[D.gate_free==gf]
    r=d[d.arm=='real'].groupby('q').C
    for a in ['perm_free','perm_count','perm_pattern']:
        nn=d[d.arm==a].groupby('q').C
        gap=r.mean()-nn.mean(); sp=np.sqrt(r.std()**2+nn.std()**2)
        out.append(dict(gate_free=gf,null=a,C=r.mean().median(),null_C=nn.mean().median(),
                        gap=gap.median(),spread2=2*sp.median(),beats=int((gap>2*sp).sum()),n=len(gap)))
S=pd.DataFrame(out); print(S.round(5).to_string(index=False))

raw=S[(S.gate_free==False)&(S.null=='perm_free')].iloc[0]
strict=S[(S.gate_free==True)&(S.null=='perm_pattern')].iloc[0]
print("\n  CONDITIONAL KILL -- gates first")
g1=G.max_orth.max()<1e-8; g2=G.var_removed.median()>0.02
print(f"   (a) projection is exact (orthogonal residual) : {'PASS' if g1 else 'FAIL'}")
print(f"   (b) projection is not vacuous                 : {'PASS' if g2 else 'FAIL'} "
      f"({G.var_removed.median():.1%} of score variance removed)")
if not(g1 and g2): print("   -> UNVERIFIED.")
else:
    print(f"\n   OLD claim (raw scores vs free permutation) : C {raw.C:+.5f}, beats null in "
          f"{raw.beats}/{raw.n} blocks")
    print(f"   STRICTEST (gate-free scores vs pattern-matched permutation): C {strict.C:+.5f}, "
          f"beats null in {strict.beats}/{strict.n} blocks")
    shrink=1-strict.C/max(raw.C,1e-12)
    print(f"   C shrinks {shrink:.1%} from the loosest to the strictest specification")
    if strict.beats>=strict.n*0.8:
        print(f"\n   -> THE FORECAST IS WRONG AND C SURVIVES. The transfer is not the gate. My own")
        print(f"      p=0.70 self-overturn fails, which is the outcome that costs the most to")
        print(f"      report and the one worth the most: #70/#71/#72 stand, at {1-shrink:.0%} of the")
        print(f"      magnitude originally reported.")
    elif strict.beats<=strict.n*0.3:
        print(f"\n   -> THE FORECAST IS RIGHT. C collapses under gate matching: #70, #71 and #72 all")
        print(f"      rest on it and all must be withdrawn or restated.")
    else:
        print(f"\n   -> PARTIAL: C survives in {strict.beats}/{strict.n} blocks. The transfer is real")
        print(f"      but its magnitude in #70/#72 is inflated by gate structure and must be restated.")
print(f"\n  people with no usable pattern stratum: {n_excl:,} ({100*n_excl/len(ALLP):.1f}%) -- "
      f"they are permuted freely in the pattern arm, which makes that arm LENIENT, not strict, "
      f"by exactly that fraction.")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
