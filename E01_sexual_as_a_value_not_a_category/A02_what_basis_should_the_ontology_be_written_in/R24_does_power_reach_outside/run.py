"""
E01 A02 R24 -- does POWER survive the check that scoped SUBSTANCE?

#39 scoped SUBSTANCE to seven blocks. #40 cleared the person factors. POWER is now the only role
axis in this project claimed to hold release-wide, and A02's decision -- what basis the ontology
should be written in -- rests on it. It has never had the same check.

POWER is built from three indicator families that come from SPECIFIC blocks: receivepain/givepain
(the sadomasochism gate), the "eagerly beg" items, and the worship items. Blocks gated on those
categories CONTRIBUTED to its construction; the rest did not.

ESTIMAND        held-out variance in a block's within-person option profile explained by POWER, in
                blocks that contributed to POWER vs blocks that did not.
IDENTIFICATION  identified -- contribution is determined by which gate column routes into a block,
                which is observable (A01 R04's branching map).
WORLDS          A  POWER is release-wide: it predicts non-contributing blocks above the placebo
                B  it is local like SUBSTANCE: no gain outside its own blocks
KILL (CONDITIONAL) gate: POWER must predict CONTRIBUTING blocks, and the permuted null must be
                   STABLE across blocks (the df-cost correction from #40 -- not required to be ~0).
                   then: non-contributing effect > 2x placebo effect -> RELEASE-WIDE, A02's basis
                         decision keeps a coordinate
                         effect <= placebo -> LOCAL, and A02 has NO surviving release-wide role axis
POSITIVE CTRL   contributing blocks.
NEGATIVE CTRL   permuted POWER, per block.
PLACEBO         the person's endorsement rate in the contributing blocks -- same source, no role
                content. This is the control that decided #39 and #40.
SEEDS           4.
MULTIPLICITY    every block reported by contribution status.
IMPOSSIBLE      a POWER measure built outside any block -- the indicators are block-gated by design.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
from lib.rounds import round_path
exec(open(round_path('16_dimensionality.py')).read().split("allq=list(B)")[0])
allq=list(B)
Ax=pd.read_csv('data/derived/agent_patient.csv')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in Ax.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=pd.concat([z(pd.to_numeric(Ax[c],errors='coerce'))*sg[c] for c in pc],axis=1).mean(axis=1)
br=pd.read_csv('data/derived/branching.csv')
gate={int(r.qi):str(r.gate) for _,r in br.iterrows()}
CONTRIB_GATES=('receivepain','givepain','sadomasochism','humiliation','nonconsent','bondage','worship','powerdynamic','obedience')
contrib={q: any(g in gate.get(q,'') for g in CONTRIB_GATES) for q in allq}
print(f"blocks contributing to POWER's construction: {sum(contrib.values())} of {len(allq)}")
RATE=pd.Series(np.nan,index=POWER.index)
inc=[q for q in allq if contrib[q]]
acc={}
for q in inc:
    s=pd.Series((B[q]['R']!=0).mean(1),index=B[q]['ppl'])
    for p,v in s.items(): acc.setdefault(p,[]).append(v)
RATE=pd.Series({p:np.mean(v) for p,v in acc.items()})
rows=[]
for t in allq:
    blk=B[t]; common=np.array([p for p in blk['ppl'] if p in POWER.index and np.isfinite(POWER.get(p,np.nan))
                               and p in RATE.index])
    if len(common)<500: continue
    ia=np.searchsorted(blk['ppl'],common); Y=blk['R'][ia]
    prop=(Y!=0).mean(1)
    Xb=np.c_[np.ones(len(common)),prop,COV.loc[common].values]
    pw=POWER.reindex(common).values.reshape(-1,1); rt=RATE.reindex(common).values.reshape(-1,1)
    for seed in (1,2,3,4):
        rng=np.random.default_rng(seed)
        idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
        def r2(X):
            b,*_=lstsq(X[tr],Y[tr],rcond=None); return 1-((Y[te]-X[te]@b)**2).sum()/((Y[te])**2).sum()
        base=r2(Xb)
        rows.append(dict(block=t,contributes=contrib[t],seed=seed,
                         gain=r2(np.c_[Xb,pw])-base, placebo=r2(np.c_[Xb,rt])-base,
                         null=r2(np.c_[Xb,pw[rng.permutation(len(common))]])-base))
G=pd.DataFrame(rows); G.to_csv(OUT/'power_reach.csv',index=False)
S=G.groupby('contributes')[['gain','placebo','null']].median().round(5)
S.index=['NON-contributing blocks','contributing blocks']
print("\n=== held-out gain from POWER ===")
print(S.to_string())
print(f"  blocks: {G[~G.contributes].block.nunique()} non-contributing, {G[G.contributes].block.nunique()} contributing")
co=S.loc['contributing blocks']; no=S.loc['NON-contributing blocks']
eff_c=co['gain']-co['null']; eff_n=no['gain']-no['null']; eff_p=no['placebo']-no['null']
null_sd=float(G.groupby('block').null.median().std())
gate_pos=eff_c>0.002; gate_neg=null_sd<0.02
print(f"\nCONDITIONAL KILL -- gate first (effect = gain - null, per #40)")
print(f"  POWER predicts contributing blocks (effect>0.002) : {'PASS' if gate_pos else 'FAIL'} ({eff_c:.5f})")
print(f"  permuted null stable across blocks (sd<0.02)      : {'PASS' if gate_neg else 'FAIL'} (sd {null_sd:.5f})")
print(f"\n  contributing effect {eff_c:.5f} | non-contributing {eff_n:.5f} | placebo there {eff_p:.5f}")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
elif eff_n>2*abs(eff_p) and eff_n>0.001:
    print(f"  -> RELEASE-WIDE : POWER reaches blocks that did not build it, at {100*eff_n/max(eff_c,1e-9):.0f}% of its home effect")
elif eff_n<=abs(eff_p):
    print("  -> LOCAL like SUBSTANCE : A02 has NO surviving release-wide role axis")
else: print("  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
