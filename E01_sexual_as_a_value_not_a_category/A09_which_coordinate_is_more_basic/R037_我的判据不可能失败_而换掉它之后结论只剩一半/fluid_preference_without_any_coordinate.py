import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #448a established the fluid content does not depend on the coordinate carrying it. So stop
   using the coordinate. What does a purely content-defined quantity -- how much a person
   picks in the six fluid blocks relative to the other twenty-six -- predict?

This is the first round in this session that uses NONE of the eigenvector machinery, so it is
an INDEPENDENT INSTRUMENT on the same claim. If #446b's "the content is the sturdier of the
two" is true, this arm should be CLEANER, not messier.

Worlds
  A  it reproduces the coordinate's profile -> the coordinate was only ever a noisy version of
     this content quantity, and the "two halves" material on the page can be simplified.
  B  the profile differs -> content and coordinate carry different things and both stay.

⚠ Direction is fixed BY CONSTRUCTION (a difference of two pick rates: higher = relatively more
   drawn to the fluid blocks), so NO ANCHORING IS NEEDED. That is the concrete benefit of
   leaving the eigenvectors behind, and it is worth writing down: five of this project's
   published errors came from an eigenvector's arbitrary sign.
⚠ The fluid regex is untouched, character for character.
CONTROL : the new quantity must correlate with the coordinate-based fluid direction -- if it
          does not, they are not the same construct and world A is not even testable.
CONTROL2: a SHAM content quantity -- six blocks chosen at random instead of the fluid six --
          must NOT clear the bar, or the design would flag any arbitrary block split.
MULTIPLICITY: four outcomes -> the null of the maximum (#440b, standing practice).
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show
from lib.nulls import perm_in

GATE=Gate("R493 content without a coordinate")
_R489=(ROOT/'E01_sexual_as_a_value_not_a_category/A179_two_ways_of_being_involved'
            /'R489_does_the_fluid_pattern_replicate/run.py').read_text()
exec(_R489.split('"""',2)[2].split('FLU=re.compile')[0])
FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)   # ⚠ 一字未改
fl=np.array([bool(FLU.search(n)) for n in NAMES]); nf=int(fl.sum())
print(f"六个体液块(纯题面):{[i for i in range(NB) if fl[i]]}")

# 纯内容量:六块的勾选率 − 其余 26 块的勾选率。方向由构造固定 -> **不需要锚定**。
RATE=(np.where(np.isfinite(A),A,0)+np.where(np.isfinite(B),B,0))/2
RATE=np.where(np.isfinite(A)|np.isfinite(B),RATE,np.nan)
# ⚠ #449a: third namespace collision this session (#427e, #443a). Stop diagnosing them one at
# a time -- PRIVATISE everything the later code needs BEFORE any splice runs. The rule is
# structural: a cross-round exec imports a namespace, so nothing that must survive it may be
# left under a name the other round also uses.
_RATE=RATE.copy(); _NB=int(NB)
def content(mask_blocks):
    a=np.nanmean(_RATE[mask_blocks],axis=0); b=np.nanmean(_RATE[~mask_blocks],axis=0)
    return a-b
FLUIDPREF=content(fl)
print(f"`FLUIDPREF` n={int(np.isfinite(FLUIDPREF).sum()):,} · "
      f"值域 [{np.nanmin(FLUIDPREF):+.3f}, {np.nanmax(FLUIDPREF):+.3f}] · "
      f"中位 {np.nanmedian(FLUIDPREF):+.4f} · **方向由构造固定,无需锚定**")

_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
rawcsv=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(rawcsv['age'].dropna().astype(str).unique()))}
AGE=rawcsv['age'].astype(str).map(AGEmap).values.astype(float)
MM = M & np.isfinite(AGE) & np.isfinite(FLUIDPREF)
gg=MM&np.isfinite(C3)
r_c3=float(np.corrcoef(FLUIDPREF[gg],C3[gg])[0,1])
print(f"n = **{int(MM.sum()):,}** · `FLUIDPREF` 与 `c3⁻` 相关 **{r_c3:+.4f}**")
GATE.asserted("CONTROL the content quantity and the coordinate are related at all",
              abs(r_c3)>0.10, f"r = {r_c3:+.4f}", kind="control")

def coef(v,y,idx):
    X=np.column_stack([np.ones(idx.sum()), z(v,idx), z(A,idx), z(Bv,idx),
                       z(ncat,idx), z(AGE,idx)])
    return float(np.linalg.lstsq(X,z(np.asarray(y,dtype=float),idx),rcond=None)[0][1])

rows=[dict(outcome=nm, b=coef(FLUIDPREF,y,MM)) for nm,y in OUT.items()]
NP_=400
nul=np.zeros((NP_,len(OUT)))
for i in range(NP_):
    pv=perm_in(FLUIDPREF,MM,seed=16000+i)
    for j,(nm,y) in enumerate(OUT.items()): nul[i,j]=abs(coef(pv,y,MM))
thr=float(np.percentile(nul.max(1),95))
T=pd.DataFrame(rows); T['sig']=T.b.abs()>thr
show(T, HERE/'results/fluidpref.csv', n=6, label="纯内容量")
print(f"   **多重性阈(四个里最大 |b| 的零分布 95 分位)= {thr:.5f}**")

rg=np.random.default_rng(29); sham=[]
for s in range(20):
    mk=np.zeros(_NB,bool); mk[rg.choice(_NB,nf,replace=False)]=True
    v=content(mk); g=MM&np.isfinite(v)
    sham.append(max(abs(coef(v,y,g)) for y in OUT.values()))
sham=np.array(sham)
print(f"   **CONTROL2 假内容量(随机六块 ×20)最大 |b| 中位 {np.median(sham):.5f} · "
      f"越阈 {int((sham>thr).sum())}/20**")
GATE.asserted("CONTROL2 a sham six-block split does not clear the bar",
              int((sham>thr).sum())<=2, f"{int((sham>thr).sum())}/20 sham splits cleared", kind="control")

nsig=int(T.sig.sum())
GATE.asserted("KILL the content quantity predicts something", nsig>0,
              f"{nsig}/4 clear: {list(T[T.sig].outcome)}")
verdict = "CONTENT_PREDICTS" if nsig>0 else "NOTHING"
print(f"\n判决 = {verdict}({nsig}/4)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),thr=thr,r_with_c3=r_c3,nsig=nsig,
               sham_clear=int((sham>thr).sum()),rows=T.to_dict('records')),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
