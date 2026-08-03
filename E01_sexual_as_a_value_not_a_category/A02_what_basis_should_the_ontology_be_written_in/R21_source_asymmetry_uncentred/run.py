"""
E01 A02 R21 -- is the male/female source opposition real, or induced by centring?

R20 reported person residuals across source classes at male-female -0.706, male-neutral -0.486,
female-neutral -0.263, and I flagged it as "not explained by the arithmetic". That flag was itself
unchecked: I removed a person mean across 3 classes, which forces negative correlations among the
residuals, and the exact induced value depends on each person's cell counts -- so eyeballing it
against a nominal -0.5 proves nothing.

Two fixes in one round.
  (1) UNCENTRED estimate: correlate the raw self-minus-other score computed WITHIN male-source
      substances against the same score computed WITHIN female-source substances. No person mean
      removed, so no mechanical negativity to argue about.
  (2) The centred version kept, but compared to a null that preserves the centring geometry exactly
      -- persons permuted within (substance, act) cells, which destroys real person structure while
      leaving every cell count and the whole centring operation untouched.

ESTIMAND        the association between a person's role position for male-source and female-source
                substances, estimated uncentred; and the centred residual correlation against a
                geometry-preserving null.
IDENTIFICATION  identified. The uncentred quantity has no induced component by construction.
WORLDS          A  real opposition: uncentred correlation clearly negative, and the centred value
                   is more negative than the geometry-preserving null
                B  arithmetic: uncentred correlation near zero, centred value matches its null
KILL (CONDITIONAL, per the rule added to P16 after three overrides -- the threshold is evaluated
      ONLY IF the controls pass, otherwise UNVERIFIED)
      gate: positive control fires AND negative control returns null
      then: uncentred |r| > 0.15 and centred below its null's 5th pct -> REAL opposition
            uncentred |r| < 0.05 and centred inside its null                -> ARITHMETIC
            otherwise                                                        -> UNVERIFIED
POSITIVE CTRL   within-source-class correlation (male substance vs the OTHER male substance) must be
                clearly positive -- the same person measured twice on the same class. If that is not
                positive the score is not measuring anything and no verdict follows.
NEGATIVE CTRL   persons permuted within (substance, act): all correlations must go to ~0 for the
                uncentred estimate, and to the induced geometry for the centred one.
SHAM            the same computation on the non-role options of the same blocks.
NOISE FLOOR     sd across 200 permutations.
MULTIPLICITY    3 class pairs x 2 estimators (centred/uncentred) x 3 seeds, all reported.
SEEDS           3.
IMPOSSIBLE      more substances per source class -- 2 male, 2 female, 3 neutral is what exists.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
FLUID={7:'precum',8:'saliva',9:'squirt',11:'urine',83:'ejaculate',6:'breastmilk',10:'sweat'}
SRC={'precum':'male','ejaculate':'male','squirt':'female','breastmilk':'female',
     'saliva':'neutral','urine':'neutral','sweat':'neutral'}
ACTS={'consume':(r'consuming it myself',r'others consuming it'),
      'produce':(r'^(making|ejaculating|squirting).*myself|myself$',r'^others (making|ejaculating|squirting)'),
      'play':(r'playing with it myself',r'others playing with it'),
      'orifice':(r'into my orifices',r"into others' orifices")}
rec=[]; shm=[]
for qi,sub in FLUID.items():
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower(); role=np.zeros(len(opt),bool)
    for act,(rs,ro) in ACTS.items():
        a=np.flatnonzero(lo.str.contains(rs,regex=True).values); b=np.flatnonzero(lo.str.contains(ro,regex=True).values)
        if len(a)==0 or len(b)==0: continue
        role[a]=True; role[b]=True
        rec += [(p,sub,act,v) for p,v in zip(ppl,M[:,a].mean(1)-M[:,b].mean(1))]
    nr=np.flatnonzero(~role); h=len(nr)//2
    if h>=1: shm += [(p,sub,'sham',v) for p,v in zip(ppl,M[:,nr[:h]].mean(1)-M[:,nr[h:2*h]].mean(1))]
D=pd.DataFrame(rec,columns=['person','substance','act','d']); D['src']=D.substance.map(SRC)
S=pd.DataFrame(shm,columns=['person','substance','act','d']); S['src']=S.substance.map(SRC)
print(f"cells {len(D):,}  persons {D.person.nunique():,}  |  sham cells {len(S):,}")
def uncentred(T):
    w=T.groupby(['person','src']).d.mean().unstack()
    out={}
    for a,b in itertools.combinations(['male','female','neutral'],2):
        if a in w and b in w:
            m=w[a].notna()&w[b].notna()
            out[f"{a}-{b}"]=(float(np.corrcoef(w[a][m],w[b][m])[0,1]),int(m.sum()))
    return out
def centred(T):
    T=T.copy(); g=T.d.mean(); pm=T.groupby('person').d.mean()-g
    T['r']=T.d-g-T.person.map(pm)
    w=T.groupby(['person','src']).r.mean().unstack()
    out={}
    for a,b in itertools.combinations(['male','female','neutral'],2):
        if a in w and b in w:
            m=w[a].notna()&w[b].notna()
            out[f"{a}-{b}"]=(float(np.corrcoef(w[a][m],w[b][m])[0,1]),int(m.sum()))
    return out
def permute(T,seed):
    T=T.copy(); rng=np.random.default_rng(seed)
    T['person']=T.groupby(['substance','act']).person.transform(lambda s: rng.permutation(s.values))
    return T
UO,CO=uncentred(D),centred(D)
print("\n=== POSITIVE CONTROL: same source class, two different substances ===")
w=D.groupby(['person','substance']).d.mean().unstack()
for a,b in [('precum','ejaculate'),('squirt','breastmilk'),('saliva','urine')]:
    if a in w and b in w:
        m=w[a].notna()&w[b].notna()
        print(f"   {a:11s} vs {b:11s} r = {np.corrcoef(w[a][m],w[b][m])[0,1]:+.3f}  n={int(m.sum()):,}")
UN=[uncentred(permute(D,s)) for s in range(200)]
CN=[centred(permute(D,s)) for s in range(200)]
print("\n=== between source classes ===")
print(f"{'pair':16s} {'UNCENTRED':>10s} {'null':>16s} | {'CENTRED':>9s} {'null (geometry-preserving)':>28s}")
rows=[]
for k in UO:
    un=np.array([d[k][0] for d in UN if k in d]); cn=np.array([d[k][0] for d in CN if k in d])
    print(f"{k:16s} {UO[k][0]:+10.3f} {un.mean():+8.3f}±{un.std():.3f} | {CO[k][0]:+9.3f} {cn.mean():+18.3f}±{cn.std():.3f}")
    rows.append(dict(pair=k,uncentred=UO[k][0],unc_null=un.mean(),unc_sd=un.std(),
                     centred=CO[k][0],cen_null=cn.mean(),cen_sd=cn.std(),n=UO[k][1]))
T=pd.DataFrame(rows); T.to_csv(OUT/'source_asymmetry.csv',index=False)
US=uncentred(S)
print(f"\n=== SHAM (non-role options) uncentred: " + "  ".join(f"{k} {v[0]:+.3f}" for k,v in US.items()))
pos=np.corrcoef(*[w[c][w['precum'].notna()&w['ejaculate'].notna()] for c in ('precum','ejaculate')])[0,1]
mf=T[T.pair=='male-female'].iloc[0]
gate_pos = pos>0.20
gate_neg = abs(mf.unc_null)<0.05
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control (precum vs ejaculate r={pos:+.3f} > 0.20) : {'PASS' if gate_pos else 'FAIL'}")
print(f"  negative control (uncentred null |{mf.unc_null:+.3f}| < 0.05) : {'PASS' if gate_neg else 'FAIL'}")
if not (gate_pos and gate_neg):
    print("  -> gate FAILED : verdict UNVERIFIED, threshold NOT evaluated")
else:
    beyond = mf.centred < mf.cen_null-1.64*mf.cen_sd
    if abs(mf.uncentred)>0.15 and beyond: print(f"  -> uncentred {mf.uncentred:+.3f}, centred beyond its null : REAL OPPOSITION")
    elif abs(mf.uncentred)<0.05 and not beyond: print(f"  -> uncentred {mf.uncentred:+.3f}, centred at its null ({mf.cen_null:+.3f}) : ARITHMETIC. R20's flag withdrawn.")
    else: print(f"  -> uncentred {mf.uncentred:+.3f}, centred {mf.centred:+.3f} vs null {mf.cen_null:+.3f} : UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
