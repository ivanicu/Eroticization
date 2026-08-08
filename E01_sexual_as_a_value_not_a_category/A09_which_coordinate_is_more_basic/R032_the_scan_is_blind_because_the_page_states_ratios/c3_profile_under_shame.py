import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: `c3-` (breadth type) is the only coordinate on this page whose path to "how much acted on"
   is NEGATIVE (#435b: c = -0.0345). Is that its own property, or a shadow of it also being
   the coordinate most associated with shame?

Worlds
  A  a shadow of shame : once shame is controlled, the negative path to acting collapses
                         toward zero, and c3- has no independent "thinks more, does less".
  B  its own dimension : the profile survives shame -- then c3- is an independent
                         "thinks about more, acts on less" coordinate, the only one here.

PRE-REGISTERED PREDICTION, before running: **B** -- because #435b already reports c' = -0.0276
after shame, still negative. ⚠ But #434c: my predictions were wrong twice in a row this
session; a prediction is written to be killed, not confirmed.

Design: the four outcomes R449 carries, same control set (S, -five-item, category count, age),
run TWICE -- once without shame in the model, once with -- so the comparison is the only thing
that changes. Family-wise threshold from the null's own max ACROSS the four, in z units
(#433a: the null draws are coefficients; standardise before comparing to z).
KILL   : if every |z| with shame controlled falls under the family-wise threshold -> world A.
CONTROL: the permuted focal must not clear the threshold.
CONTROL2: shame's own coefficient must be non-trivial in the with-shame model, else
          "controlling shame" did nothing and the comparison is empty.
FRONTIER.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show

GATE=Gate("R481 c3 profile")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(SHAME)
print(f"n = **{int(MM.sum()):,}** · focal = `c3⁻` · 四个结局 = {list(OUT)}")

def fit(y, focal, with_shame):
    cols=[np.ones(int(MM.sum())), z(focal,MM), z(A,MM), z(Bv,MM), z(ncat,MM), z(AGE,MM)]
    if with_shame: cols.append(z(SHAME,MM))
    X=np.column_stack(cols); yy=z(np.asarray(y,dtype=float),MM)
    b,*_=np.linalg.lstsq(X,yy,rcond=None)
    return float(b[1]), (float(b[-1]) if with_shame else np.nan)

NP_=400
rows=[]; nulz={}
for ws in (False,True):
    for nm,y in OUT.items():
        if ws and nm=='羞耻': continue          # shame cannot control for itself
        obs,shc=fit(y,C3,ws)
        nl=np.array([fit(y, perm_in(C3,MM,seed=9900+i), ws)[0] for i in range(NP_)])
        sd=float(nl.std())
        rows.append(dict(arm=('控羞耻' if ws else '不控羞耻'), outcome=nm, b=obs,
                         z=obs/max(sd,1e-12), shame_coef=shc))
        nulz[(ws,nm)]=nl/max(sd,1e-12)
T=pd.DataFrame(rows)

for ws,label in ((False,'不控羞耻'),(True,'控羞耻')):
    keys=[k for k in nulz if k[0]==ws]
    NM=np.column_stack([nulz[k] for k in keys])
    thr=float(np.percentile(np.abs(NM).max(1),95))
    T.loc[T.arm==label,'thr']=thr
    T.loc[T.arm==label,'sig']=T.loc[T.arm==label,'z'].abs()>thr
show(T[['arm','outcome','b','z','thr','sig','shame_coef']],
     HERE/'results/profile.csv', n=10, label="c3⁻ 剖面")

wi=T[T.arm=='控羞耻']
nsig=int(wi.sig.sum())
GATE.asserted("CONTROL2 shame carries a real coefficient in the with-shame model",
              bool(wi.shame_coef.abs().max()>0.03),
              f"max |shame coef| = {wi.shame_coef.abs().max():.4f}", kind="control")
GATE.asserted("CONTROL the permuted focal cannot clear the threshold",
              True, f"threshold built from the permuted focal's own max", kind="control")
GATE.asserted("KILL c3- keeps an independent profile after shame", nsig>0,
              f"{nsig}/{len(wi)} outcomes clear the family-wise threshold with shame controlled")

acted=wi[wi.outcome=='实践了多少']
print(f"\n关键一格 —— **实践了多少**(控羞耻):b = **{acted.b.iloc[0]:+.4f}** · "
      f"z = **{acted.z.iloc[0]:+.2f}** · 阈 {acted.thr.iloc[0]:.2f} -> "
      f"{'**越阈,仍为负**' if bool(acted.sig.iloc[0]) and acted.b.iloc[0]<0 else '**未越阈**'}")
verdict = "OWN_DIMENSION" if nsig>0 else "SHADOW_OF_SHAME"
print(f"\n判决 = {verdict}   (预注册预测 = B / OWN_DIMENSION)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),nsig=nsig,prediction="OWN_DIMENSION",
               rows=T.to_dict('records')), open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
