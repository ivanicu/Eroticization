"""
E01 A02 R19 -- settle POWER vs SUBSTANCE distinctness on a reliability ladder.

The three-axis claim rests on disattenuated correlations of 0.233-0.362, and ADVERSARY_FORECAST #5
correctly called that fragile because GAZE has alpha 0.163. R16-R18 tried to bypass the
disattenuation entirely and froze: the discriminant design is underpowered by construction.

So stop trying to settle all three at once. GAZE genuinely has only 4 indicators in this release
(searched the option level: the other hits were regex false positives -- "Wetlook", "Video game
characters"). But POWER-SUBSTANCE does not involve GAZE, POWER has 3 indicators and SUBSTANCE has
7, and the reliability-ladder method already worked on theta in A05 R10. Settle the pair that can
be settled and scope the claim to it.

ESTIMAND        r_true(POWER, SUBSTANCE), estimated by extrapolation along a MEASURED reliability
                ladder rather than by dividing once by two alphas.
IDENTIFICATION  identified if r_obs / sqrt(rel_P * rel_S) is flat across rungs. If it drifts, the
                measurement model is wrong and only bounds are reportable.
SCOPE           POWER 3 indicators, SUBSTANCE 7 (one per fluid block) · n varies by rung overlap.
WORLDS          A  distinct axes: r_obs/sqrt(rel_P*rel_S) flat and extrapolates well below 1
                B  one construct measured twice: extrapolates toward 1.0
KILL            PRE-REGISTERED: extrapolated r_true > 0.70 -> one construct, the three-axis claim
                loses its POWER-SUBSTANCE leg. r_true < 0.45 with drift under 30% -> DISTINCT
                confirmed for this pair. Between -> UNVERIFIED.
POSITIVE CTRL   THE DECISIVE ONE: run the identical pipeline on SUBSTANCE against ITSELF, split
                into two disjoint halves of its own 7 indicators. The method MUST extrapolate to
                ~1.0 there. A disattenuation procedure that cannot recover 1.0 for a measure
                against itself cannot be trusted to report 0.23 for two different ones.
                Also verified to fail: two independent noise composites must extrapolate to ~0.
NEGATIVE CTRL   noise-vs-noise and noise-vs-POWER at every rung.
SHAM            POWER against a non-role composite of the same indicator count.
NOISE FLOOR     sd of r_obs across 3 seeds at each rung.
MULTIPLICITY    rungs x 3 seeds x 5 pairings, all reported.
SEEDS           3.
IMPOSSIBLE      GAZE cannot be laddered here -- 4 indicators. That leg stays on its single
                unstable disattenuation and is reported as such, not folded into this verdict.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
Ax=pd.read_csv('data/derived/agent_patient.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in Ax.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=[z(pd.to_numeric(Ax[c],errors='coerce'))*sg[c] for c in pc]
SUB=[]
for qi in [7,8,9,11,83,6,10]:
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower()
    a=np.flatnonzero(lo.str.contains(r'(myself|\bmy\b)',regex=True).values)
    b=np.flatnonzero(lo.str.contains(r'(others|other )',regex=True).values)
    if len(a) and len(b): SUB.append(pd.Series(M[:,a].mean(1)-M[:,b].mean(1),index=ppl).reindex(df.index))
rng0=np.random.default_rng(7)
NOISE=[pd.Series(rng0.normal(size=len(df)),index=df.index) for _ in range(7)]
NOISE2=[pd.Series(rng0.normal(size=len(df)),index=df.index) for _ in range(7)]
SHAM=[z(pd.to_numeric(df[c],errors='coerce')) for c in
      ['opennessvariable','neuroticismvariable','extroversionvariable','consciensiousnessvariable',
       'agreeablenessvariable','powerlessnessvariable','pornhabit']]
print(f"indicators -- POWER {len(POWER)}  SUBSTANCE {len(SUB)}  NOISE {len(NOISE)}  SHAM {len(SHAM)}")
def comp(items): return pd.concat(items,axis=1).mean(axis=1)
def rel(items,seed):
    if len(items)<2: return np.nan
    r=np.random.default_rng(seed); idx=r.permutation(len(items)); h=len(idx)//2
    a=comp([items[i] for i in idx[:h]]); b=comp([items[i] for i in idx[h:2*h]])
    m=a.notna()&b.notna(); rr=np.corrcoef(a[m],b[m])[0,1]
    return float(2*rr/(1+rr)) if rr>-1 else np.nan
def robs(A_,B_):
    a,b=comp(A_),comp(B_); m=a.notna()&b.notna()
    return float(stats.pearsonr(a[m],b[m])[0]) if m.sum()>300 else np.nan
rows=[]
PAIRS={'POWER-SUBSTANCE':(POWER,SUB),'SUBSTANCE-SUBSTANCE_posctrl':(SUB[:3],SUB[3:]),
       'NOISE-NOISE2_negctrl':(NOISE,NOISE2),'POWER-NOISE_negctrl':(POWER,NOISE),
       'POWER-SHAM':(POWER,SHAM)}
for name,(X,Y) in PAIRS.items():
    for kx,ky,seed in itertools.product(range(2,len(X)+1),range(2,len(Y)+1),[5,15,25]):
        if kx>len(X) or ky>len(Y): continue
        r=np.random.default_rng(seed*31+kx*7+ky)
        xs=[X[i] for i in r.choice(len(X),kx,replace=False)]
        ys=[Y[i] for i in r.choice(len(Y),ky,replace=False)]
        rx,ry=rel(xs,seed),rel(ys,seed)
        if not (np.isfinite(rx) and np.isfinite(ry)) or rx<=0.02 or ry<=0.02: continue
        ro=robs(xs,ys)
        if not np.isfinite(ro): continue
        rows.append(dict(pair=name,kx=kx,ky=ky,seed=seed,rel_prod=round(rx*ry,4),
                         r_obs=round(ro,4),r_true=round(ro/np.sqrt(rx*ry),4)))
T=pd.DataFrame(rows); T.to_csv(OUT/'ladder.csv',index=False)
print("\n=== r_true = r_obs / sqrt(rel_x * rel_y), by pairing and reliability tercile ===")
T['tercile']=T.groupby('pair').rel_prod.transform(lambda s: pd.qcut(s,3,labels=['low','mid','high'],duplicates='drop'))
print(T.pivot_table(index='pair',columns='tercile',values='r_true',aggfunc='median',observed=False).round(3).to_string())
print("\n=== r_obs (should RISE with reliability if the model holds) ===")
print(T.pivot_table(index='pair',columns='tercile',values='r_obs',aggfunc='median',observed=False).round(3).to_string())
print("\n=== summary ===")
for name in PAIRS:
    s=T[T.pair==name]
    if not len(s): print(f"  {name:30s} no admissible rungs"); continue
    print(f"  {name:30s} r_true median {s.r_true.median():+.3f}  IQR [{s.r_true.quantile(.25):+.3f},{s.r_true.quantile(.75):+.3f}]  cells {len(s)}")
ps=T[T.pair=='POWER-SUBSTANCE']; pcl=T[T.pair=='SUBSTANCE-SUBSTANCE_posctrl']
lo=ps.groupby('tercile',observed=False).r_true.median()
drift=(lo.get('high',np.nan)-lo.get('low',np.nan))/abs(lo.get('low',np.nan)) if len(lo)>1 else np.nan
print(f"\nPOSITIVE CONTROL SUBSTANCE vs itself: r_true = {pcl.r_true.median():+.3f}  (must be ~1.0)")
print(f"NEGATIVE CONTROL noise vs noise      : r_true = {T[T.pair=='NOISE-NOISE2_negctrl'].r_true.median():+.3f}  (must be ~0)")
print("\nPRE-REGISTERED KILL, evaluated:")
if not (0.75<pcl.r_true.median()<1.35):
    print(f"  -> positive control returned {pcl.r_true.median():+.3f}, not ~1.0 : METHOD UNFIT, no verdict")
else:
    rt=ps.r_true.median()
    if rt>0.70: print(f"  -> r_true = {rt:+.3f} > 0.70 : ONE CONSTRUCT, POWER-SUBSTANCE leg lost")
    elif rt<0.45 and abs(drift)<0.30: print(f"  -> r_true = {rt:+.3f} < 0.45, drift {100*drift:+.0f}% : DISTINCT CONFIRMED for this pair")
    else: print(f"  -> r_true = {rt:+.3f}, drift {100*drift:+.0f}% : UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
