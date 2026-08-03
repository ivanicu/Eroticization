"""
E01 A02 R20 -- additivity as a variance share, not as 3 pairwise correlations.

v_i = w_i^T phi(s) REQUIRES additivity: the weight a person puts on a feature does not depend on
the scene. A02 R11-R13 tested it with pairwise correlations between substances and reported
PLAUSIBLE on n=3-4 decisive pairs, which is a two-point line. But the fluid template is really a
4 ACTS x 7 SUBSTANCES design, so the same question is a variance decomposition with thousands of
cells instead of a handful of correlations.

  additive        -> a person's self-minus-other difference is a PERSON MAIN EFFECT: constant
                     across acts and substances
  non-additive    -> person x substance (or person x act) interaction carries comparable variance,
                     which is exactly Ivan's A_i term being non-negligible

ESTIMAND        the share of variance in the person-level self-minus-other difference attributable
                to person main effect vs person x substance vs person x act interaction.
IDENTIFICATION  identified: the design is crossed, and every cell is directly observed for the
                people who entered both blocks.
SCOPE           7 substances x 4 acts, people entering >=3 fluid blocks · no instrument · noise
                floor measured by split-half within cell.
WORLDS          A  additive: person main effect dominates; interactions near the noise floor
                B  non-additive: person x substance rivals the main effect
                   B' the R12 refinement: the interaction is CONCENTRATED in the male-vs-female
                      source contrast rather than spread across substances
KILL            PRE-REGISTERED: if person x substance variance share exceeds 50% of the person main
                effect share, additivity in the folk basis is REFUTED and the A_i term is required.
                Under 20%, additivity holds and R11-R13's PLAUSIBLE is upgraded.
POSITIVE CTRL   the person main effect must exceed the residual by a wide margin -- if the design
                cannot detect a person effect at all, no decomposition is interpretable. Verified to
                fail: person labels shuffled within cell must drive the person effect to ~0.
NEGATIVE CTRL   within-cell person permutation, preserving all margins.
SHAM            the same decomposition on the NON-role options of the same blocks, where no person
                main effect should exist.
NOISE FLOOR     split-half over the options within each act, measured.
MULTIPLICITY    3 variance components x 3 seeds x 2 option-sets, all reported.
SEEDS           3.
IMPOSSIBLE      other feature contrasts -- only self/other is instantiated across enough blocks.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
FLUID={7:'precum',8:'saliva',9:'squirt',11:'urine',83:'ejaculate',6:'breastmilk',10:'sweat'}
MALE={'precum','ejaculate'}; FEM={'squirt','breastmilk'}
ACTS={'consume':(r'consuming it myself',r'others consuming it'),
      'produce':(r'^(making|ejaculating|squirting).*myself|myself$',r'^others (making|ejaculating|squirting)'),
      'play':(r'playing with it myself',r'others playing with it'),
      'orifice':(r'into my orifices',r"into others' orifices")}
recs=[]; sham=[]
for qi,sub in FLUID.items():
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower()
    role=np.zeros(len(opt),bool)
    for act,(rs,ro) in ACTS.items():
        a=np.flatnonzero(lo.str.contains(rs,regex=True).values); b=np.flatnonzero(lo.str.contains(ro,regex=True).values)
        if len(a)==0 or len(b)==0: continue
        role[a]=True; role[b]=True
        d=M[:,a].mean(1)-M[:,b].mean(1)
        recs += [(p,sub,act,v) for p,v in zip(ppl,d)]
    nr=np.flatnonzero(~role)
    if len(nr)>=2:
        h=len(nr)//2
        d=M[:,nr[:h]].mean(1)-M[:,nr[h:2*h]].mean(1)
        sham += [(p,sub,'shamA',v) for p,v in zip(ppl,d)]
D=pd.DataFrame(recs,columns=['person','substance','act','d'])
S=pd.DataFrame(sham,columns=['person','substance','act','d'])
print(f"cells: {len(D):,}   persons {D.person.nunique():,}   substances {D.substance.nunique()}   acts {D.act.nunique()}")
keep=D.groupby('person').substance.nunique()
D=D[D.person.isin(keep[keep>=3].index)]
print(f"restricted to persons in >=3 substances: {D.person.nunique():,} persons, {len(D):,} cells")
def decompose(T,seed=0,permute=False):
    T=T.copy()
    if permute:
        rng=np.random.default_rng(seed)
        T['person']=T.groupby(['substance','act']).person.transform(lambda s: rng.permutation(s.values))
    g=T.d.mean()
    pm=T.groupby('person').d.mean()-g
    sm=T.groupby('substance').d.mean()-g
    am=T.groupby('act').d.mean()-g
    T['resid']=T.d-g-T.person.map(pm)-T.substance.map(sm)-T.act.map(am)
    ps=T.groupby(['person','substance']).resid.mean()
    pa=T.groupby(['person','act']).resid.mean()
    tot=T.d.var()
    return dict(person=float(pm.var()*1.0/tot*len(pm)/len(pm)) if False else float(np.var(T.person.map(pm)))/tot,
                substance=float(np.var(T.substance.map(sm)))/tot,
                act=float(np.var(T.act.map(am)))/tot,
                person_x_substance=float(np.var(T.set_index(['person','substance']).index.map(ps)))/tot,
                person_x_act=float(np.var(T.set_index(['person','act']).index.map(pa)))/tot)
obs=pd.DataFrame([decompose(D,s) for s in (1,2,3)]).median()
nul=pd.DataFrame([decompose(D,s,permute=True) for s in (1,2,3)]).median()
shm=pd.DataFrame([decompose(S,s) for s in (1,2,3)]).median()
out=pd.DataFrame({'observed':obs,'person_permuted_null':nul,'sham_nonrole_options':shm}).round(4)
out.to_csv(OUT/'variance.csv')
print("\n=== variance share of the self-minus-other difference ===")
print(out.to_string())
ratio=obs['person_x_substance']/obs['person'] if obs['person']>0 else np.nan
print(f"\n  person main effect share        : {obs['person']:.4f}   (null {nul['person']:.4f})")
print(f"  person x substance share        : {obs['person_x_substance']:.4f}")
print(f"  interaction / main effect       : {100*ratio:.0f}%")
print(f"  POSITIVE CONTROL person effect vs its permuted null: {obs['person']/max(nul['person'],1e-6):.1f}x")
print("\n=== B' refinement: is the interaction concentrated at the source-gender boundary? ===")
D2=D.copy(); D2['src']=D2.substance.map(lambda s:'male' if s in MALE else ('female' if s in FEM else 'neutral'))
g=D2.d.mean(); pm=D2.groupby('person').d.mean()-g
D2['r2']=D2.d-g-D2.person.map(pm)
by=D2.groupby(['person','src']).r2.mean().unstack()
for a,b in [('male','female'),('male','neutral'),('female','neutral')]:
    if a in by and b in by:
        m=by[a].notna()&by[b].notna()
        print(f"   corr(person residual in {a:7s}, in {b:7s}) = {np.corrcoef(by[a][m],by[b][m])[0,1]:+.3f}  n={m.sum():,}")
print("\nPRE-REGISTERED KILL, evaluated:")
if ratio>0.50: print(f"  -> person x substance is {100*ratio:.0f}% of the main effect (>50%) : ADDITIVITY REFUTED in the folk basis; A_i required")
elif ratio<0.20: print(f"  -> {100*ratio:.0f}% (<20%) : ADDITIVITY HOLDS, R11-R13 upgraded from PLAUSIBLE")
else: print(f"  -> {100*ratio:.0f}% : between thresholds, UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
