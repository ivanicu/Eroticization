import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #451d(1) named recall bias as the only serious threat to "what arrived early feels more
   fixed". #451's NEXT said to reuse the page's existing recall-bias control rather than
   invent one, and to READ what it actually is first.

⚠ TWO THINGS FOUND BY READING, BEFORE ANY COMPUTATION:
  (1) the page cites `#289` for "recall bias controlled". **`#289` is about measurement
      invariance by sex.** The project's actual recall-bias work is **`#114`** (via R166):
      people report their more highly rated categories as starting earlier. The citation on
      the page is wrong and is corrected by this round.
  (2) `#114`'s bias is defined as: reported onset MINUS the population schedule for that
      category, **minus that person's own overall precocity**. So it is a WITHIN-PERSON tilt
      across categories, and the person-level mean is removed BY CONSTRUCTION. `EARLY` -- the
      focal in #451 -- IS that person-level mean. **So the bias as measured cannot, by its own
      definition, be what produced #451.**

Worlds
  A  by construction and in fact orthogonal : a person's recall tilt does not track their mean
     onset -> #451's warning can be narrowed to name what remains untested.
  B  they do track each other : then the construction argument is wrong somewhere and #451
     needs the warning it currently carries, or stronger.

⚠ THE RULE I WROTE FOR MYSELF IN #451's NEXT: if the existing control cannot be applied in
   this arm, **write it down and stop -- do not swap in an approximation and pass it off as
   the real thing.** So the rating-intensity covariate below is run and reported AS A
   STAND-IN, explicitly labelled, never as `#114`'s control.
CONTROL : the stand-in must actually be related to onset, or controlling it is vacuous.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R496 can recall bias produce it")
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
_O=np.array(O,dtype=float).copy()          # 人 × 类别 的起始年龄
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
BEL=np.asarray(OUT['能不能改'],dtype=float)

# ---- ① 把「按构造正交」量出来,而不是只论证它
# 人群时间表 = 每个类别的中位起始年龄;偏离 = 报告 − 时间表 − 该人整体早熟度
sched=np.nanmedian(_O,axis=0)
dev=_O-sched[None,:]
prec=np.nanmean(dev,axis=1)                       # 该人的整体早熟度
tilt_input=dev-prec[:,None]                       # #114 的量所在的空间
resid_spread=np.nanstd(tilt_input,axis=1)
r_orth=float(np.corrcoef(*[x[np.isfinite(prec)&np.isfinite(_EARLY)]
                           for x in (prec,_EARLY)])[0,1])
r_tilt_early=float(np.corrcoef(*[x[np.isfinite(resid_spread)&np.isfinite(_EARLY)]
                                 for x in (resid_spread,_EARLY)])[0,1])
print(f"① 整体早熟度 与 `EARLY` 相关 **{r_orth:+.4f}**(应当很高 —— 它们几乎是同一个量)")
print(f"   `#114` 的偏离**扣掉早熟度之后**,其人内展布 与 `EARLY` 相关 **{r_tilt_early:+.4f}**")
print(f"   ⇒ **`#114` 测的那一层,与 `EARLY` 所在的那一层,是分开的。**")

# ---- ② 替代品(明确标注),两臂并报
RAT=[c for c in pd.read_csv('data/derived/inventory.csv').query("kind=='RATING_0_5'").col
     if c in raw.columns]
RI=np.nanmean(np.column_stack([pd.to_numeric(raw[c],errors='coerce').values for c in RAT]),axis=1)
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY) & np.isfinite(RI)
print(f"\n② 替代品 = **整体评分强度**(`RATING_0_5` 67 列的人均,方向已由 `#424b` 定为一个方向)")
print(f"   ⚠ **这不是 `#114` 的控制,是替代品,如实标注**(`#451` NEXT 的自订规则)")
print(f"   n = **{int(MM.sum()):,}** · 评分强度 与 `EARLY` 相关 **"
      f"{np.corrcoef(RI[MM],_EARLY[MM])[0,1]:+.4f}**")
GATE.asserted("CONTROL the stand-in is actually related to onset",
              abs(np.corrcoef(RI[MM],_EARLY[MM])[0,1])>0.05,
              f"r(rating intensity, EARLY) = {np.corrcoef(RI[MM],_EARLY[MM])[0,1]:+.4f}", kind="control")

def coef(v,y,idx,extra=None):
    cols=[np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx), z(C3,idx), z(ncat,idx), z(AGE,idx)]
    if extra is not None: cols.append(z(extra,idx))
    return float(np.linalg.lstsq(np.column_stack(cols),z(np.asarray(y,dtype=float),idx),rcond=None)[0][1])

b0=coef(_EARLY,BEL,MM); b1=coef(_EARLY,BEL,MM,extra=RI)
NP_=400
nl0=np.array([abs(coef(perm_in(_EARLY,MM,seed=19000+i),BEL,MM)) for i in range(NP_)])
nl1=np.array([abs(coef(perm_in(_EARLY,MM,seed=19000+i),BEL,MM,extra=RI)) for i in range(NP_)])
t0=float(np.percentile(nl0,95)); t1=float(np.percentile(nl1,95))
T=pd.DataFrame([dict(arm='不控评分强度',b=b0,thr=t0,sig=abs(b0)>t0),
                dict(arm='**控评分强度(替代品)**',b=b1,thr=t1,sig=abs(b1)>t1)])
show(T, HERE/'results/two_arms.csv', n=4, label="两臂")
print(f"   保留比 = **{abs(b1)/max(abs(b0),1e-12):.3f}**")

GATE.asserted("KILL the effect survives the stand-in control", bool(abs(b1)>t1),
              f"b {b0:+.4f} -> {b1:+.4f}, threshold {t1:.5f}")
verdict = "SURVIVES" if abs(b1)>t1 else "COLLAPSES"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,b_before=b0,b_after=b1,thr_before=t0,thr_after=t1,
               r_prec_early=r_orth,r_tilt_early=r_tilt_early,n=int(MM.sum()),
               citation_fix="page cited #289 for recall bias; the real work is #114 (R166)"),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
