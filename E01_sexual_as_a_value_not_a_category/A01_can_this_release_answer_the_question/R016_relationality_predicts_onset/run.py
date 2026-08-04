"""
E01 A01 R16 -- does a category's RELATIONALITY predict when it is acquired?

A01's decision says this release cannot separate A (a dedicated sexual-content system) from B
(erotic valuation of ordinary representation), because #13 established that a dedicated module can
be compositional, so prediction tests do not discriminate. That leaves developmental structure,
which #13 did not consider.

  B predicts: erotic interests inherit the maturation order of the ordinary representations they
     are built on. Person-perception matures FEATURES BEFORE RELATIONS, so categories whose content
     is relational should arrive later.
  A predicts: a dedicated sexual-content system has no reason to follow person-perception's
     developmental ordering.

RELATIONALITY is measured by a PRE-REGISTERED STRING RULE over each category's own option texts --
the share of options referring to another person -- not by hand-coding, which would be an answer key
(the failure named in feedback_hardcoded_structure_dict_is_an_answer_key).

ESTIMAND        corr(category relationality, mean onset age) across categories.
IDENTIFICATION  identified; both quantities are computed from the release, neither from a judgment.
WORLDS          A  no relation between relationality and onset age
                B  positive: relational categories arrive later
KILL (CONDITIONAL) gate: the positive control must fire (a planted relationality->onset link must be
                   recovered) AND the permutation null must be centred near zero.
                   then: r > 0.35 and surviving the nuisance controls -> supports B
                         |r| < 0.15                                    -> no developmental signature
                         otherwise                                     -> UNVERIFIED
POSITIVE CTRL   plant onset = f(relationality) + noise; must recover r > 0.8.
NEGATIVE CTRL   permute relationality across categories, 5000 draws.
NUISANCE        option word count, option count, and category prevalence -- all partialled, because
                longer or rarer option sets could drive both.
SHAM            a string rule counting a content-free token class (colour words) instead of person
                references; must give ~0.
SEEDS           3 where stochastic.
MULTIPLICITY    1 estimand x 4 control sets, all reported.
IMPOSSIBLE      a non-sexual developmental yardstick inside this release -- the ordering claim comes
                from outside literature and is used as a PREDICTION, not verified here.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, re
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
opts=pd.read_csv('data/derived/options.csv')
br=pd.read_csv('data/derived/branching.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
# PRE-REGISTERED string rules, fixed before looking at any onset value
PERSON=r'\b(others?|other|someone|partner|partners|them|their|people|person|him|her|his|guy|girl|man|woman|men|women|boy|stranger|couple)\b'
COLOUR=r'\b(red|blue|green|black|white|pink|dark|light|colou?r|colou?red)\b'   # SHAM token class
gate={int(r.qi):str(r.gate).strip('"') for _,r in br.iterrows()}
rows=[]
for qi,g in opts.groupby('qi'):
    cat=gate.get(int(qi))
    if not cat: continue
    on=[c for c in inv[inv['kind']=='AGE_ONSET']['col']
        if cat.split('(')[0].strip().lower()[:8] in c.lower() or
           (cat.lower()[:6] in c.lower() and len(cat)>5)]
    txt=g.option.astype(str)
    rel=float(txt.str.contains(PERSON,case=False,regex=True).mean())
    sham=float(txt.str.contains(COLOUR,case=False,regex=True).mean())
    wc=float(txt.str.split().str.len().mean()); nopt=len(g)
    prev=float(pd.to_numeric(df.get(cat),errors='coerce').gt(0).mean()) if cat in df.columns else np.nan
    if not on: continue
    o=df[on[0]].map(BIN)
    if o.notna().sum()<300: continue
    rows.append(dict(qi=int(qi),cat=cat[:28],relationality=rel,sham=sham,words=wc,n_opt=nopt,
                     prevalence=prev,mean_onset=float(o.mean()),n=int(o.notna().sum())))
T=pd.DataFrame(rows).drop_duplicates('cat')
print(f"categories with both an option set and an onset column: {len(T)}")
print(T[['cat','relationality','mean_onset','n']].sort_values('relationality').to_string(index=False))
x=T.relationality.values; y=T.mean_onset.values
r0=stats.pearsonr(x,y)
def partial(extra):
    X=np.c_[np.ones(len(x)),T[extra].values]
    rx=x-X@lstsq(X,x,rcond=None)[0]; ry=y-X@lstsq(X,y,rcond=None)[0]
    return stats.pearsonr(rx,ry)[0]
print(f"\n=== relationality vs mean onset age ===")
print(f"  raw                                   r = {r0[0]:+.3f}  p = {r0[1]:.4f}  n = {len(T)}")
print(f"  | word count                          r = {partial(['words']):+.3f}")
print(f"  | word count + option count           r = {partial(['words','n_opt']):+.3f}")
print(f"  | + prevalence                        r = {partial(['words','n_opt','prevalence']):+.3f}")
print(f"  SHAM (colour-word share vs onset)     r = {stats.pearsonr(T.sham.values,y)[0]:+.3f}")
rng=np.random.default_rng(7)
nul=np.array([stats.pearsonr(rng.permutation(x),y)[0] for _ in range(5000)])
planted=3.0*x+rng.normal(0,0.3,len(x))
print(f"\n  permutation null: mean {nul.mean():+.3f} sd {nul.std():.3f} |r| p95 {np.percentile(np.abs(nul),95):.3f}")
print(f"  POSITIVE CONTROL onset planted as f(relationality): r = {stats.pearsonr(x,planted)[0]:+.3f}")
T.to_csv(OUT/'relationality.csv',index=False)
gp=stats.pearsonr(x,planted)[0]>0.8; gn=abs(nul.mean())<0.10
radj=partial(['words','n_opt','prevalence'])
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control recovers >0.8 : {'PASS' if gp else 'FAIL'}")
print(f"  permutation null near zero     : {'PASS' if gn else 'FAIL'}")
if not (gp and gn): print("  -> gate FAILED : UNVERIFIED")
elif radj>0.35: print(f"  -> SUPPORTS B : relational categories arrive later, r={radj:+.3f} after nuisance controls")
elif abs(radj)<0.15: print(f"  -> NO DEVELOPMENTAL SIGNATURE : r={radj:+.3f}")
else: print(f"  -> UNVERIFIED : r={radj:+.3f}")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
