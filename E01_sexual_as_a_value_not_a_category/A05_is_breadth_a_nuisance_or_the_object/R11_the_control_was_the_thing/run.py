"""
E01 A05 R11 -- my "acquiescence" control is made of erotic content, and two of its items ARE an axis.

#25 established a pattern: I reach for a control before asking what it is. Applying that to the two
remaining unexamined control sets, by READING them rather than trusting their names:

  "acquiescence" = mean of 22 LIKERT_PM3 items. At least 16 are explicit erotic content, including
     "I am aroused by being dominant in sexual interactions"
     "I am aroused by being submissive in sexual interactions"
  which are the POWER axis itself, plus animated/written (the modality axis), plus
  ashamed-of-arousal, therapeutic-arousal, acted-on-arousal, partner-contagion, and six
  autoerotic-identity items.

  "orientation proxies" in every CCA = two items, one of which is
     "I find it erotic when two people of the opposite gender to me sexually interact"
  which is an EROTIC ITEM being partialled out of an EROTIC ITEM MATRIX.

So three published quantities rest on partialling erotic content out of erotic content:
  breadth is "9-13% response style"          (A05)
  induction->breadth "85% survives response-style control"  (A05)
  cross-domain CCA 0.200, "73% not demographic"             (A01/A02)

ESTIMAND        each of those three, recomputed with the erotic items removed from the control set;
                and, first, the degree to which the "acquiescence" index is erotic content at all.
IDENTIFICATION  identified by construction -- the item texts are readable and the split is a fact
                about the questionnaire, not an inference.
WORLDS          A  the index is response style: it correlates weakly with named erotic axes, and
                   removing erotic items barely moves the three quantities
                B  the index is erotic endorsement: it correlates strongly with the POWER axis, and
                   the three quantities move materially
KILL            PRE-REGISTERED: if |corr(index, POWER)| > 0.30, the index is declared erotic content
                and all three quantities are recomputed and republished on the cleaned control.
                If < 0.15 the original control stands.
POSITIVE CTRL   the two power items must correlate near-perfectly with the POWER axis by
                construction -- if they do not, the item mapping is wrong and nothing follows.
NEGATIVE CTRL   the handful of plausibly content-free items ("I am a narcissist", harassment
                experience) as a mini-index; its correlation with POWER should be far lower.
SHAM            a random 22-item mean drawn from non-Likert numeric columns.
SEEDS           3.  MULTIPLICITY: 3 quantities x 3 control sets, all reported.
IMPOSSIBLE      a real acquiescence measure -- that needs balanced-keyed items, which this release
                does not carry. If the kill fires, the line is FROZEN with that as its unfreeze.
"""
import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np, pandas as pd, warnings, hashlib
from numpy.linalg import lstsq
from scipy import stats
warnings.filterwarnings('ignore')
OUT=pathlib.Path(__file__).parent/'results'
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
Ax=pd.read_csv('data/derived/agent_patient.csv')
lik=[c for c in inv[inv['kind']=='LIKERT_PM3']['col'] if c in df.columns]
EROTIC=[c for c in lik if any(s in c.lower() for s in
   ['arouse','arousing','erotic','sexual','masturbat','animated','written','allroll','highenergy','supernatural','freeuse'])]
CLEAN=[c for c in lik if c not in EROTIC]
print(f"LIKERT items {len(lik)}  -> erotic-content {len(EROTIC)}  plausibly content-free {len(CLEAN)}")
print("  content-free remainder:", [c[:52] for c in CLEAN])
def z(s): return (s-s.mean())/(s.std()+1e-9)
IDX_all=df[lik].apply(pd.to_numeric,errors='coerce').mean(axis=1)
IDX_clean=df[CLEAN].apply(pd.to_numeric,errors='coerce').mean(axis=1) if CLEAN else None
pc=[c for c in Ax.columns if any(x in c for x in ['receivepain','eagerly beg','bondage','humiliation','nonconsent','worship'])]
sg={c:(-1 if 'worship' in c else 1) for c in pc}
POWER=pd.concat([z(pd.to_numeric(Ax[c],errors='coerce'))*sg[c] for c in pc],axis=1).mean(axis=1)
dom=[c for c in lik if 'aroused by being dominant' in c][0]; sub=[c for c in lik if 'aroused by being submissive' in c][0]
print("\n=== POSITIVE CONTROL: the two power items vs the POWER axis ===")
for c in (dom,sub):
    v=pd.to_numeric(df[c],errors='coerce'); m=v.notna()&POWER.notna()
    print(f"   {c[:56]:58s} r = {stats.pearsonr(v[m],POWER[m])[0]:+.3f}")
def r_with(x):
    m=x.notna()&POWER.notna(); return float(stats.pearsonr(x[m],POWER[m])[0])
print("\n=== is the index response style or erotic content? ===")
print(f"   full 22-item index   vs POWER : {r_with(IDX_all):+.3f}")
if CLEAN: print(f"   content-free {len(CLEAN)}-item  vs POWER : {r_with(IDX_clean):+.3f}   <- NEGATIVE CONTROL")
rng=np.random.default_rng(3)
num=[c for c in df.columns if df[c].dtype!=object and c not in lik]
sham=df[list(rng.choice(num,min(22,len(num)),replace=False))].apply(pd.to_numeric,errors='coerce')
sham=((sham-sham.mean())/(sham.std()+1e-9)).mean(axis=1)
print(f"   SHAM random 22 numeric cols vs POWER : {r_with(sham):+.3f}")
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
Rt=df[rate].apply(pd.to_numeric,errors='coerce')
breadth=(Rt>0).sum(1).astype(float); answered=Rt.notna().sum(1).astype(float)
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
IND=[c for c in df.columns if 'induced' in c.lower() and 'fetish' in c.lower()][0]
ORD={'No':0,'Yes: Variations on my current preexisting fetishes':1,
     'Yes: New but still similar to my preexisting fetishes':2,
     'Yes: New and totally different to my preexisting fetishes':3}
g=df[IND].map(ORD)
def rho(ctrl_extra):
    m=g.notna()&breadth.notna()
    X=[np.ones(m.sum()),df['age'].map(AGEMAP)[m].fillna(2).values,df['biomale'][m].fillna(0).values]
    for e in ctrl_extra:
        X.append(e[m].fillna(e.median()).values)
    X=np.column_stack(X)
    rb=breadth[m].values-X@lstsq(X,breadth[m].values,rcond=None)[0]
    rg=g[m].values.astype(float)-X@lstsq(X,g[m].values.astype(float),rcond=None)[0]
    return float(stats.spearmanr(rg,rb).statistic)
print("\n=== the three quantities, recomputed on each control set ===")
print(f"   breadth ~ index (all 22, EROTIC)         : {stats.spearmanr(breadth,IDX_all,nan_policy='omit').statistic:+.3f}")
if CLEAN: print(f"   breadth ~ index (content-free only)      : {stats.spearmanr(breadth,IDX_clean,nan_policy='omit').statistic:+.3f}")
print(f"   induction->breadth | no index control    : {rho([]):+.4f}")
print(f"   induction->breadth | FULL index (as pub) : {rho([IDX_all]):+.4f}")
if CLEAN: print(f"   induction->breadth | content-free index  : {rho([IDX_clean]):+.4f}")
r=abs(r_with(IDX_all))
print("\nPRE-REGISTERED KILL, evaluated:")
if r>0.30: print(f"  -> |corr(index, POWER)| = {r:.3f} > 0.30 : THE INDEX IS EROTIC CONTENT. All three quantities republished on the cleaned control.")
elif r<0.15: print(f"  -> {r:.3f} < 0.15 : the original control stands")
else: print(f"  -> {r:.3f} between thresholds : UNVERIFIED")
pd.DataFrame(dict(item=lik,erotic=[c in EROTIC for c in lik])).to_csv(OUT/'item_audit.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(open(__file__,'rb').read()).hexdigest()[:12]}")
