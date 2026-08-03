"""
E01 A02 R22 -- is the male/female non-transfer a role fact or a BODY FEASIBILITY artifact?

R21 established the sharpest structural fact in the project: the self/other role feature transfers
male<->neutral at +0.479 and male<->female at +0.063, below even a meaningless sham (+0.266). I
named orientation as the rival. There is a blunter one I had not noticed:

    "ejaculating myself" requires a penis.  "squirting myself" requires a vulva.

So on GENDERED substances the self-pole is partly gated by the respondent's own body, in opposite
directions for the two classes -- which would suppress the male-female correlation with no role
content involved at all. Neutral substances (saliva, urine, sweat) have no such gate, which is
exactly why they bridge to both sides.

TWO ORTHOGONAL FIXES, either of which removes the artifact if it is one:
  (1) stratify by biomale -- within a sex, feasibility is constant
  (2) drop the 'produce' act entirely and rebuild from consume/play/orifice, all of which are
      feasible for any body

ESTIMAND        the male-female transfer of the role feature, purged of body feasibility.
WORLDS          A  role fact: +0.063 stays low within sex AND without the produce act
                B  feasibility artifact: it rises materially under either fix
                B' orientation: it rises under an orientation stratification but not under (1) or (2)
KILL (CONDITIONAL -- threshold evaluated ONLY if the gate passes)
      gate: positive control fires (within-class transfer stays high under the same fix) AND
            negative control null stays ~0
      then: male-female transfer > 0.25 under either fix -> ARTIFACT, R21's fact withdrawn
            stays < 0.12 under both                       -> ROLE FACT, R21 upheld
            otherwise                                     -> UNVERIFIED
POSITIVE CTRL   male-neutral transfer must stay high under each fix; if a fix destroys the
                signal everywhere, it is destroying data rather than removing a confound.
NEGATIVE CTRL   person-permutation within (substance, act) under each fix.
SHAM            the non-role option contrast under each fix.
SEEDS           3.  MULTIPLICITY: 3 pairs x 4 fixes x 3 seeds, all reported.
IMPOSSIBLE      a body-neutral gendered substance -- none exists.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
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
male=df['biomale']
def transfer(T):
    w=T.groupby(['person','src']).d.mean().unstack()
    out={}
    for a,b in itertools.combinations(['male','female','neutral'],2):
        if a in w and b in w:
            m=w[a].notna()&w[b].notna()
            out[f"{a}-{b}"]=(float(np.corrcoef(w[a][m],w[b][m])[0,1]) if m.sum()>150 else np.nan,int(m.sum()))
    return out
def perm(T,seed):
    T=T.copy(); rng=np.random.default_rng(seed)
    T['person']=T.groupby(['substance','act']).person.transform(lambda s: rng.permutation(s.values)); return T
FIXES={'ALL (as R21)':D,
       'men only':D[D.person.map(male)==1],
       'women only':D[D.person.map(male)==0],
       'no produce act':D[D.act!='produce'],
       'men only + no produce':D[(D.person.map(male)==1)&(D.act!='produce')]}
print(f"{'fix':24s} {'male-female':>12s} {'male-neutral':>13s} {'female-neutral':>15s} {'n(m-f)':>8s} {'null(m-f)':>10s}")
rows=[]
for name,T in FIXES.items():
    if len(T)<2000: continue
    o=transfer(T); nl=np.nanmean([transfer(perm(T,s)).get('male-female',(np.nan,0))[0] for s in range(30)])
    mf=o.get('male-female',(np.nan,0)); mn=o.get('male-neutral',(np.nan,0)); fn=o.get('female-neutral',(np.nan,0))
    print(f"{name:24s} {mf[0]:+12.3f} {mn[0]:+13.3f} {fn[0]:+15.3f} {mf[1]:8,d} {nl:+10.3f}")
    rows.append(dict(fix=name,male_female=mf[0],male_neutral=mn[0],female_neutral=fn[0],n=mf[1],null=nl))
R=pd.DataFrame(rows); R.to_csv(OUT/'feasibility.csv',index=False)
print("\n=== SHAM (non-role options), same fixes ===")
for name in ['ALL (as R21)','men only','no produce act']:
    T=S if name=='ALL (as R21)' else (S[S.person.map(male)==1] if name=='men only' else S)
    o=transfer(T); print(f"   {name:22s} male-female {o.get('male-female',(np.nan,0))[0]:+.3f}")
base=R[R.fix=='ALL (as R21)'].iloc[0]
gate_pos = R.male_neutral.min()>0.20
gate_neg = abs(R.null).max()<0.06
print(f"\nCONDITIONAL KILL -- gate first")
print(f"  positive control (male-neutral stays >0.20 under every fix; min {R.male_neutral.min():+.3f}) : {'PASS' if gate_pos else 'FAIL'}")
print(f"  negative control (|null| < 0.06 under every fix; max {abs(R.null).max():.3f}) : {'PASS' if gate_neg else 'FAIL'}")
if not (gate_pos and gate_neg):
    print("  -> gate FAILED : UNVERIFIED, threshold not evaluated")
else:
    fixed=R[R.fix!='ALL (as R21)'].male_female
    if fixed.max()>0.25: print(f"  -> male-female rises to {fixed.max():+.3f} under a fix : FEASIBILITY ARTIFACT, R21's fact withdrawn")
    elif fixed.max()<0.12: print(f"  -> stays below 0.12 under every fix (max {fixed.max():+.3f}) : ROLE FACT, R21 upheld")
    else: print(f"  -> max under fixes {fixed.max():+.3f} : UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
