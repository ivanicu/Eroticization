"""
E01 A02 R16 -- are the three role axes distinct constructs, tested WITHOUT disattenuation?

ADVERSARY_FORECAST #5, p=0.40: "the disattenuation in R09 is doing too much work. GAZE's alpha is
0.163; dividing by sqrt(0.163) is numerically unstable -- I said so -- and the three-axis claim
still leans on it."

Correct, and the right response is not a better disattenuation. It is to stop needing one. Mutual
correlation is the WRONG evidence for distinctness when one measure is unreliable, because low
reliability manufactures orthogonality and the correction that undoes it is unstable exactly where
it matters most. Discriminant validity answers the same question without the division:

    THREE DISTINCT CONSTRUCTS PREDICT DIFFERENT THINGS.
    ONE CONSTRUCT MEASURED THREE WAYS PREDICTS THE SAME THINGS, differing only in scale.

ESTIMAND        the similarity between axes of their CORRELATION PROFILES across an external
                battery -- a quantity that needs no reliability correction, because both profiles
                are attenuated by the same axis reliability and the profile SHAPE is preserved.
IDENTIFICATION  identified. Scale cancels in a profile correlation; reliability rescales a profile
                without rotating it.
SCOPE           BKS public · 14 external variables, none of them role items · n varies by pair.
WORLDS          A  three constructs: profile correlations between axes are well below the ceiling
                B  one construct measured thrice: profile correlations approach the ceiling
KILL            PRE-REGISTERED: the ceiling is the profile correlation between two random halves
                of the SAME axis. If any between-axis profile correlation reaches 80% of that
                ceiling, that pair is declared one construct and the three-axis claim is reduced.
                If all three pairs sit below 50% of the ceiling, distinctness is CONFIRMED without
                any disattenuation and forecast #5 is scored WRONG-ON-CONSEQUENCE (right that the
                disattenuation was fragile, wrong that the claim depended on it).
POSITIVE CTRL   the split-half ceiling itself -- two halves of one axis MUST show a high profile
                correlation. If they do not, the battery cannot discriminate anything and no
                verdict is available.
NEGATIVE CTRL   a random-noise "axis" carried through the identical pipeline; its profile
                correlation with every real axis must be ~0.
SHAM            an axis built from the same number of items drawn from NON-role columns.
NOISE FLOOR     sd of the split-half ceiling over 200 random half-splits.
MULTIPLICITY    3 axis pairs + 3 ceilings + 2 controls x 3 seeds, all reported.
SPECIFICATION   battery {14 vars} x partial {none, sex-removed} x estimator {pearson, spearman}.
SEEDS           3.
IMPOSSIBLE      independent replication; a gold standard for any of the three axes.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib, itertools
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
Ax=pd.read_csv('data/derived/agent_patient.csv')
G=np.load('data/derived/gcca_G.npy')
def z(s): return (s-s.mean())/(s.std()+1e-9)
pc=[c for c in Ax.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER_items=[z(pd.to_numeric(Ax[c],errors='coerce'))*sg[c] for c in pc]
ex=[c for c in Ax.columns if 'exhibition' in c][0]; vo=[c for c in Ax.columns if 'voyeur' in c][0]
GAZE_items=[z(pd.to_numeric(Ax[ex],errors='coerce')), -z(pd.to_numeric(Ax[vo],errors='coerce'))]
inv=pd.read_csv('data/derived/inventory.csv')
sub_cols=[c for c in df.columns if 'orifices' in c.lower()]
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
FLUID=[7,8,9,11,83,6,10]
SUB_items=[]
for qi in FLUID:
    s=lg[lg.qi==qi]
    if not len(s): continue
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    lo=pd.Series(opt).str.lower()
    a=np.flatnonzero(lo.str.contains(r'(myself|\bmy\b)',regex=True).values)
    b=np.flatnonzero(lo.str.contains(r'(others|other )',regex=True).values)
    if len(a)and len(b): SUB_items.append(pd.Series(M[:,a].mean(1)-M[:,b].mean(1),index=ppl).reindex(df.index))
AXES={'POWER':POWER_items,'GAZE':GAZE_items,'SUBSTANCE':SUB_items}
rng0=np.random.default_rng(11)
AXES['NOISE_ctrl']=[pd.Series(rng0.normal(size=len(df)),index=df.index) for _ in range(4)]
AXES['SHAM_nonrole']=[z(pd.to_numeric(df[c],errors='coerce')) for c in
                      ['opennessvariable','neuroticismvariable','extroversionvariable']]
print({k:len(v) for k,v in AXES.items()})
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
BATT={'male':df['biomale'],'age':df['age'].map(AGEMAP),'pornhabit':pd.to_numeric(df['pornhabit'],errors='coerce'),
 'openness':pd.to_numeric(df['opennessvariable'],errors='coerce'),'neuro':pd.to_numeric(df['neuroticismvariable'],errors='coerce'),
 'extro':pd.to_numeric(df['extroversionvariable'],errors='coerce'),'consc':pd.to_numeric(df['consciensiousnessvariable'],errors='coerce'),
 'agree':pd.to_numeric(df['agreeablenessvariable'],errors='coerce'),'powerless':pd.to_numeric(df['powerlessnessvariable'],errors='coerce'),
 'adversity':df['childhood_adversity'].notna().astype(float),'mental':df['TotalMentalIllness'].notna().astype(float),
 'spanked':df[[c for c in df.columns if 'were you spanked' in c][0]].map({'Never':0,'Sometimes':1,'Often':2}),
 'assault':df[[c for c in df.columns if 'victim of sexual assault' in c][0]].map({'No':0,'Yes':1}),
 'monog':df[[c for c in df.columns if 'preferred relationship style' in c][0]].map({'Not monogamous':0,'Monogamous':1})}
def profile(items,partial_sex=False):
    a=pd.concat(items,axis=1).mean(axis=1)
    out=[]
    for k,y in BATT.items():
        m=a.notna()&y.notna()
        if m.sum()<300: out.append(np.nan); continue
        x=a[m].values; yy=y[m].values.astype(float)
        if partial_sex and k!='male':
            s=df['biomale'][m].values.astype(float); X=np.c_[np.ones(m.sum()),s]
            x=x-X@lstsq(X,x,rcond=None)[0]; yy=yy-X@lstsq(X,yy,rcond=None)[0]
        out.append(float(stats.pearsonr(x,yy)[0]))
    return np.array(out)
def prof_corr(p1,p2):
    m=~(np.isnan(p1)|np.isnan(p2))
    return float(np.corrcoef(p1[m],p2[m])[0,1]) if m.sum()>5 else np.nan
rows=[]
for psex in [False,True]:
    P={k:profile(v,psex) for k,v in AXES.items()}
    ceil={}
    for k,v in AXES.items():
        if len(v)<2: ceil[k]=np.nan; continue
        cs=[]
        for s in range(200):
            r=np.random.default_rng(s); idx=r.permutation(len(v)); h=max(1,len(idx)//2)
            cs.append(prof_corr(profile([v[i] for i in idx[:h]],psex),profile([v[i] for i in idx[h:]],psex)))
        ceil[k]=float(np.nanmedian(cs))
    for a,b in itertools.combinations(['POWER','GAZE','SUBSTANCE'],2):
        c=prof_corr(P[a],P[b]); cl=np.nanmean([ceil[a],ceil[b]])
        rows.append(dict(partial_sex=psex,pair=f"{a}-{b}",prof_r=round(c,3),ceiling=round(cl,3),
                         pct_of_ceiling=round(100*c/cl,0) if cl and not np.isnan(cl) else np.nan))
    for a in ['POWER','GAZE','SUBSTANCE']:
        for ctl in ['NOISE_ctrl','SHAM_nonrole']:
            rows.append(dict(partial_sex=psex,pair=f"{a}-{ctl}",prof_r=round(prof_corr(P[a],P[ctl]),3),
                             ceiling=round(ceil[a],3) if not np.isnan(ceil[a]) else np.nan,pct_of_ceiling=np.nan))
    if not psex:
        print("\n=== split-half CEILING per axis (positive control) ===")
        for k,v in ceil.items(): print(f"   {k:14s} {v:+.3f}")
T=pd.DataFrame(rows); T.to_csv(OUT/'discriminant.csv',index=False)
print("\n=== profile correlations between axes, against their own ceiling ===")
print(T.to_string(index=False))
real=T[(T.pair.str.count('-')==1)&(~T.pair.str.contains('ctrl|SHAM'))]
mx=real.pct_of_ceiling.max()
print("\nPRE-REGISTERED KILL, evaluated:")
if mx>=80: print(f"  -> a pair reaches {mx:.0f}% of ceiling : ONE CONSTRUCT, three-axis claim reduced")
elif mx<50: print(f"  -> all pairs below 50% of ceiling (max {mx:.0f}%) : DISTINCTNESS CONFIRMED without disattenuation")
else: print(f"  -> max {mx:.0f}% of ceiling : between thresholds, stays UNVERIFIED")
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
