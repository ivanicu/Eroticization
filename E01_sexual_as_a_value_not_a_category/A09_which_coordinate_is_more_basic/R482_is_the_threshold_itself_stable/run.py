import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #437b left "how much acted on" at z = -2.25 against a family-wise threshold of 2.37 -- a
   5% gap, with the threshold estimated from 400 permutations. Is the gap real, or is the
   THRESHOLD itself noisier than the gap?

Worlds
  A  the threshold is stable : it settles in a narrow band and z stays under it -> "not
     separable from shame at this design" becomes a resolved statement, not an open one.
  B  the threshold is noisy  : it moves materially with the number of permutations, or its
     own resampling interval spans the gap -> then EVERY conclusion on this page that uses
     the same family-wise construction is stated more precisely than it can be, which is a
     sentence about the whole page, not about one cell.

Design: 4,000 permutations (10x), threshold recomputed at 400/1000/2000/4000, PLUS the
threshold's own bootstrap interval (resample the 4,000 permutation draws). The observed z is
recomputed once and does not depend on the null.
KILL   : world A requires (i) the threshold's own 95% interval to EXCLUDE the observed |z|,
         and (ii) the threshold at 4,000 to sit within 5% of the threshold at 2,000.
CONTROL: the observed coefficient must be identical to #437's (nothing but the null changed).
FRONTIER: world B forces a page-wide qualification.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show

GATE=Gate("R482 threshold stability")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(SHAME)
OUTS={k:v for k,v in OUT.items() if k!='羞耻'}

def fit(y, focal):
    X=np.column_stack([np.ones(int(MM.sum())), z(focal,MM), z(A,MM), z(Bv,MM),
                       z(ncat,MM), z(AGE,MM), z(SHAME,MM)])
    return float(np.linalg.lstsq(X,z(np.asarray(y,dtype=float),MM),rcond=None)[0][1])

NP_=4000
raw_null={nm: np.array([fit(y, perm_in(C3,MM,seed=9900+i)) for i in range(NP_)])
          for nm,y in OUTS.items()}
obs={nm: fit(y,C3) for nm,y in OUTS.items()}
sd ={nm: float(raw_null[nm].std()) for nm in OUTS}
zobs={nm: obs[nm]/sd[nm] for nm in OUTS}
NZ=np.column_stack([raw_null[nm]/sd[nm] for nm in OUTS])
mx=np.abs(NZ).max(1)

rows=[dict(nperm=k, thr=float(np.percentile(mx[:k],95))) for k in (400,1000,2000,4000)]
S=pd.DataFrame(rows)
rng=np.random.default_rng(23)
bs=np.array([float(np.percentile(rng.choice(mx,len(mx),replace=True),95)) for _ in range(400)])
tlo,thi=np.percentile(bs,[2.5,97.5])
S['thr_lo']=tlo; S['thr_hi']=thi
show(S, HERE/'results/threshold_sweep.csv', n=6, label="阈随置换次数")

zi=zobs['实践了多少']
print(f"\n观测 z(实践了多少,控羞耻)= **{zi:+.4f}**  (`#437` 报 −2.25;系数 {obs['实践了多少']:+.4f})")
print(f"阈的自助 95% 区间 = **[{tlo:.3f}, {thi:.3f}]** · 4,000 次的点估 = **{S.thr.iloc[-1]:.3f}**")
drift=abs(S.thr.iloc[-1]-S.thr.iloc[-2])/S.thr.iloc[-2]
print(f"2,000 -> 4,000 的漂移 = **{drift:.1%}**")

GATE.asserted("CONTROL the observed coefficient is unchanged from #437",
              abs(obs['实践了多少']-(-0.027281))<5e-4,
              f"b = {obs['实践了多少']:+.5f} vs #437 -0.02728", kind="control")
excl = not (tlo <= abs(zi) <= thi)
GATE.asserted("KILL-i the threshold's own interval excludes the observed |z|",
              excl, f"|z| = {abs(zi):.3f} vs threshold interval [{tlo:.3f}, {thi:.3f}]")
GATE.asserted("KILL-ii the threshold is converged (2000 -> 4000 within 5%)",
              drift<0.05, f"drift = {drift:.1%}")
verdict = "RESOLVED" if (excl and drift<0.05) else "THRESHOLD_IS_THE_LIMIT"
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,z_acted=zi,b_acted=obs['实践了多少'],
               thr4000=float(S.thr.iloc[-1]),thr_lo=tlo,thr_hi=thi,drift=drift,
               z_all={k:float(v) for k,v in zobs.items()}),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
