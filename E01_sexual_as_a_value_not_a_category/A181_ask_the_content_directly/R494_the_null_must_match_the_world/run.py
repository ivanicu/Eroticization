import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #449c's sham control failed -- an arbitrary six-block split cleared the bar in 7 of 20
   tries -- so the null result was inadmissible. The diagnosis: the bar was built from the
   WRONG NULL. Permuting people asks "is this quantity related to the outcome at all"; the
   sham varies the BLOCK SPLIT and asks "would any six blocks do". Those exclude different
   worlds, and the second is the one that matters here.

Worlds
  A  the fluid split beats arbitrary splits -> the content DOES predict independently, and
     #449's unresolved becomes a result.
  B  it does not -> then at n = 2,490 any six blocks would do, and the fluid split has no
     special claim in this arm. That is a real answer, not a failure.

Null : max |b| over the four outcomes, across >= 500 RANDOM six-block splits -- the same
   family as the sham control that broke #449. Bar = its 95th percentile.
⚠ The CRITERION does not change with the bar: it stays "max |b| over the four outcomes clears
   the 95th percentile". Only the null changes. Changing both would be moving the goalposts.
CONTROL : the person-permutation bar from #449 is recomputed here too, so the two nulls are
   printed side by side and the difference is visible rather than asserted.
CONTROL2: the random-split null must be centred well below the observed -- otherwise it has no
   power and world B would be unreadable.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R494 match the null to the world")
_R489=(ROOT/'E01_sexual_as_a_value_not_a_category/A179_two_ways_of_being_involved'
            /'R489_does_the_fluid_pattern_replicate/run.py').read_text()
exec(_R489.split('"""',2)[2].split('FLU=re.compile')[0])
FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)   # ⚠ 一字未改
fl=np.array([bool(FLU.search(n)) for n in NAMES]); nf=int(fl.sum())
RATE=(np.where(np.isfinite(A),A,0)+np.where(np.isfinite(B),B,0))/2
RATE=np.where(np.isfinite(A)|np.isfinite(B),RATE,np.nan)
_RATE=RATE.copy(); _NB=int(NB)                       # #449a:splice 之前私有化
def content(mk):
    return np.nanmean(_RATE[mk],axis=0)-np.nanmean(_RATE[~mk],axis=0)
FLUIDPREF=content(fl)

_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
rawcsv=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(rawcsv['age'].dropna().astype(str).unique()))}
AGE=rawcsv['age'].astype(str).map(AGEmap).values.astype(float)
BASE = M & np.isfinite(AGE)
MM = BASE & np.isfinite(FLUIDPREF)
print(f"n(体液口径)= **{int(MM.sum()):,}**")

def coef(v,y,idx):
    X=np.column_stack([np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
                       z(ncat,idx), z(AGE,idx)])
    return float(np.linalg.lstsq(X,z(np.asarray(y,dtype=float),idx),rcond=None)[0][1])
def maxb(v):
    g=BASE&np.isfinite(v)
    if g.sum()<500: return np.nan
    return max(abs(coef(v,y,g)) for y in OUT.values())

obs_each={nm: coef(FLUIDPREF,y,MM) for nm,y in OUT.items()}
obs=max(abs(x) for x in obs_each.values())

# ---- 零 A(#449 用的):打乱人
NP_=400
nulP=np.array([max(abs(coef(perm_in(FLUIDPREF,MM,seed=17000+i),y,MM)) for y in OUT.values())
               for i in range(NP_)])
# ---- 零 B(本轮):随机六块划分 —— 与假对照同族
rg=np.random.default_rng(61); nulS=[]
for i in range(500):
    mk=np.zeros(_NB,bool); mk[rg.choice(_NB,nf,replace=False)]=True
    v=maxb(content(mk))
    if v==v: nulS.append(v)
nulS=np.array(nulS)
barP=float(np.percentile(nulP,95)); barS=float(np.percentile(nulS,95))
T=pd.DataFrame([dict(q='观测 max|b|',v=obs),
                dict(q='零A 打乱人 · 95 分位',v=barP),
                dict(q='零A 均值',v=float(nulP.mean())),
                dict(q='零B 随机六块划分 · 95 分位',v=barS),
                dict(q='零B 均值',v=float(nulS.mean())),
                dict(q='零B 中位',v=float(np.median(nulS)))])
show(T, HERE/'results/two_nulls.csv', n=6, label="两个零")
print("\n四个结局各自的系数:")
for nm,b in obs_each.items(): print(f"   {nm:<8} **{b:+.4f}**")

GATE.asserted("CONTROL the two nulls are visibly different",
              abs(barS-barP)>1e-4, f"person-permutation bar {barP:.5f} vs block-split bar {barS:.5f}",
              kind="control")
GATE.asserted("CONTROL2 the block-split null has power (its mean is below the observation)",
              float(nulS.mean())<obs, f"null mean {nulS.mean():.5f} vs observed {obs:.5f}", kind="control")
beats = obs>barS
GATE.asserted("KILL the fluid split beats arbitrary six-block splits", beats,
              f"observed {obs:.5f} vs block-split 95th pct {barS:.5f}")
verdict = "FLUID_SPLIT_IS_SPECIAL" if beats else "ANY_SIX_WOULD_DO"
print(f"\n观测 **{obs:.5f}** · 零A 阈 **{barP:.5f}** · **零B 阈 {barS:.5f}** -> 判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=int(MM.sum()),obs=obs,bar_person=barP,bar_split=barS,
               null_split_mean=float(nulS.mean()),n_split_null=len(nulS),
               coefs=obs_each), open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
