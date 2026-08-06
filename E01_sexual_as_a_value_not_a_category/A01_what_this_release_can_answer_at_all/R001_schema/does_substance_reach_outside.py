"""
E01 A02 R23 -- does the SUBSTANCE coordinate reach outside the family it was built in?

#38 found the seven fluid blocks are partially ISOLATED: fluid<->non-fluid pairs transfer at 0.203
against 0.272 for ordinary pairs. Most of this project's role findings live in that family --
additivity (R11-R13, R20), source-gender transfer (R21-R22), the SUBSTANCE axis itself (R09, R19).
If the family is a subgraph, those are subgraph conclusions and their generality has never been
tested.

The role contrast only exists as explicitly paired options INSIDE the fluid blocks, which is why it
was built there. So generality cannot be tested by rebuilding the contrast elsewhere -- it has to be
tested by asking whether the coordinate PREDICTS elsewhere.

ESTIMAND        held-out variance in a NON-fluid block's within-person option profile explained by
                the fluid-derived SUBSTANCE score, over demographics + propensity.
IDENTIFICATION  identified for the 25 non-fluid blocks; people must appear in both, which the
                overlap allows.
WORLDS          A  SUBSTANCE is a general person coordinate: it predicts non-fluid profiles above
                   baseline and above the fluid-propensity placebo
                B  it is subgraph-local: no gain outside the family, or no better than the placebo
KILL (CONDITIONAL) gate: SUBSTANCE must predict FLUID-block profiles (where it is defined) AND the
                   permuted-SUBSTANCE null must be ~0. Otherwise UNVERIFIED.
                   then: median non-fluid gain > 2x the placebo gain -> REACHES OUTSIDE
                         gain <= placebo -> SUBGRAPH-LOCAL, and A02's role results are scoped to
                         the fluid family in README
POSITIVE CTRL   SUBSTANCE predicting fluid blocks -- the coordinate is defined there.
NEGATIVE CTRL   permuted SUBSTANCE scores.
PLACEBO         the person's overall fluid-block ENDORSEMENT RATE -- same source blocks, same
                people, no role content. Anything SUBSTANCE achieves that this also achieves is
                not about role.
SEEDS           4.
MULTIPLICITY    25 non-fluid + 7 fluid blocks x 3 predictors x 4 seeds, all reported.
IMPOSSIBLE      rebuilding the self/other contrast outside the fluid family -- no other block
                instantiates it in its option set. That is a property of the release.
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
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
FLUIDQ=[7,8,9,11,83,6,10]
ACTS={'consume':(r'consuming it myself',r'others consuming it'),
      'produce':(r'^(making|ejaculating|squirting).*myself|myself$',r'^others (making|ejaculating|squirting)'),
      'play':(r'playing with it myself',r'others playing with it'),
      'orifice':(r'into my orifices',r"into others' orifices")}
sub_parts=[]; rate_parts=[]
for qi in FLUIDQ:
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower(); acc=[]
    for act,(rs,ro) in ACTS.items():
        a=np.flatnonzero(lo.str.contains(rs,regex=True).values); b=np.flatnonzero(lo.str.contains(ro,regex=True).values)
        if len(a) and len(b): acc.append(M[:,a].mean(1)-M[:,b].mean(1))
    if acc: sub_parts.append(pd.Series(np.mean(acc,axis=0),index=ppl))
    rate_parts.append(pd.Series(M.mean(1),index=ppl))     # PLACEBO: fluid endorsement rate, no role
SUB=pd.concat(sub_parts,axis=1).mean(axis=1)
RATE=pd.concat(rate_parts,axis=1).mean(axis=1)
print(f"SUBSTANCE defined for {SUB.notna().sum():,} people from {len(sub_parts)} fluid blocks")
rows=[]
for t in allq:
    tgt=B[t]; ppl=tgt['ppl']
    common=np.array([p for p in ppl if p in SUB.index and np.isfinite(SUB.get(p,np.nan))])
    if len(common)<500: continue
    ia=np.searchsorted(ppl,common); Y=tgt['R'][ia]
    prop=(Y!=0).mean(1)
    Xb=np.c_[np.ones(len(common)),prop,COV.loc[common].values]
    s=SUB.reindex(common).values.reshape(-1,1); r=RATE.reindex(common).values.reshape(-1,1)
    for seed in (1,2,3,4):
        rng=np.random.default_rng(seed)
        idx=rng.permutation(len(common)); tr,te=idx[:int(.7*len(idx))],idx[int(.7*len(idx)):]
        def r2(X):
            b,*_=lstsq(X[tr],Y[tr],rcond=None); return 1-((Y[te]-X[te]@b)**2).sum()/((Y[te])**2).sum()
        base=r2(Xb)
        rows.append(dict(block=t,fluid=t in FLUIDQ,seed=seed,
                         gain_substance=r2(np.c_[Xb,s])-base,
                         gain_placebo=r2(np.c_[Xb,r])-base,
                         gain_null=r2(np.c_[Xb,s[rng.permutation(len(common))]])-base))
G=pd.DataFrame(rows); G.to_csv(OUT/'reach.csv',index=False)
S=G.groupby('fluid')[['gain_substance','gain_placebo','gain_null']].median().round(5)
S.index=['non-fluid blocks','fluid blocks']
print("\n=== held-out gain from the fluid-derived SUBSTANCE score ===")
print(S.to_string())
print(f"\n  blocks: {G[~G.fluid].block.nunique()} non-fluid, {G[G.fluid].block.nunique()} fluid")
fl=S.loc['fluid blocks']; nf=S.loc['non-fluid blocks']
gate_pos=fl['gain_substance']>3*abs(fl['gain_null']); gate_neg=abs(nf['gain_null'])<0.002
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  SUBSTANCE predicts fluid blocks (>3x its null) : {'PASS' if gate_pos else 'FAIL'} ({fl['gain_substance']:.5f} vs null {fl['gain_null']:.5f})")
print(f"  permuted null ~0 outside                       : {'PASS' if gate_neg else 'FAIL'} ({nf['gain_null']:.5f})")
if not (gate_pos and gate_neg): print("  -> gate FAILED : UNVERIFIED")
else:
    ratio=nf['gain_substance']/max(abs(nf['gain_placebo']),1e-9)
    print(f"  non-fluid: substance {nf['gain_substance']:.5f} vs placebo {nf['gain_placebo']:.5f}  ratio {ratio:.2f}")
    if nf['gain_substance']>2*abs(nf['gain_placebo']) and nf['gain_substance']>0.002:
        print("  -> REACHES OUTSIDE : SUBSTANCE is a general person coordinate")
    elif nf['gain_substance']<=abs(nf['gain_placebo']):
        print("  -> SUBGRAPH-LOCAL : A02's role results are scoped to the fluid family")
    else: print("  -> UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
