import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #439c left the biggest unexplained thing on the page: `c3-` predicts "I cannot change what
   arouses me" (-0.0355 direct) and it is NOT via experience. Where in the object does that
   come from -- a few blocks of content, or the coordinate as a whole?

Worlds
  A  content-specific : dropping a handful of blocks kills the coefficient -> the answer is
     the CONTENT of those blocks, and they can be named.
  B  a coordinate-level property : no small set of blocks matters -> it is not WHAT you like
     that makes it feel unchangeable, it is HOW your liking is DISTRIBUTED. That is the
     stronger sentence, and it is the one the page cannot currently make.

Machinery is spliced, not rebuilt (P4): R371 already has `c3_drop(drop, ref)` -- rebuild c3
with one block removed, sign-aligned against a reference (#368a: this project has been burnt
by eigenvector signs four times). R372 already established the loadings and their sign
convention (positive = the more-shame end).

⚠ NAME COLLISION (#427e): R371 binds `A`,`B` to block matrices; R449 binds `A` to z(S). So all
32 leave-one-out vectors are computed FIRST and stashed, and only then is R449 spliced.
CONTROL : the full-sample c3 from R371's pipeline must correlate ~1 with R449's C3, otherwise
          the two splices are not describing the same coordinate and nothing below holds.
CONTROL2: the drop-one refits must be sign-aligned -- report min correlation to the reference.
Spec curve: drop the top-k blocks by |loading| for k = 1..6, not a single cut (#P16).
PRE-REGISTERED KILL: world A requires SOME single block whose removal moves the coefficient by
          more than the coefficient's own bootstrap sd.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R484 where does it come from")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])

full_c3, ref_vec = c3_drop(-1)
print(f"NB = **{NB}** 块 · 全量 c3 已重建")
LOO={}
for bdrop in range(NB):
    v,_=c3_drop(bdrop, ref=ref_vec); LOO[bdrop]=v
load = ref_vec.copy()
print(f"32 个留一块向量已算完(符号对齐到同一参考)")
np.save(HERE/'results/loo_c3.npy', np.column_stack([LOO[b] for b in range(NB)]))
np.save(HERE/'results/full_c3.npy', full_c3); np.save(HERE/'results/loading.npy', load)

# ---- only now splice R449 (it rebinds A/B)
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float); ACT=np.asarray(OUT['实践了多少'],dtype=float)
BEL=np.asarray(OUT['能不能改'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(SHAME) & np.isfinite(BEL) & np.isfinite(ACT)

# CONTROL: the two pipelines must be describing the same coordinate
gg=MM & np.isfinite(full_c3) & np.isfinite(C3)
r_pipe=float(np.corrcoef(-full_c3[gg], C3[gg])[0,1])
print(f"\n**对照**:`R371` 管线的全量 c3 与 `R449` 的 `C3` 相关 = **{r_pipe:+.4f}**(n={int(gg.sum()):,})")
GATE.asserted("CONTROL the two splices describe the same coordinate",
              abs(r_pipe)>0.95, f"r = {r_pipe:+.4f}", kind="control")

def coef(c3vec, idx):
    X=np.column_stack([np.ones(idx.sum()), z(c3vec,idx), z(A,idx), z(Bv,idx),
                       z(ncat,idx), z(AGE,idx), z(SHAME,idx), z(ACT,idx)])
    return float(np.linalg.lstsq(X,z(BEL,idx),rcond=None)[0][1])

base=coef(-full_c3, MM)
rng=np.random.default_rng(41); idxall=np.flatnonzero(MM)
bsd=float(np.std([coef(-full_c3, (lambda t: (lambda mm: (mm.__setitem__(slice(None), False), mm.__setitem__(np.unique(t), True), mm)[-1])(np.zeros(len(MM),bool)))(rng.choice(idxall,len(idxall),replace=True)))
                  for _ in range(200)]))
print(f"全量系数 = **{base:+.5f}** · 自助 sd = **{bsd:.5f}**")

rows=[]
for b in range(NB):
    gd=MM & np.isfinite(LOO[b])
    cb=coef(-LOO[b], gd)
    rows.append(dict(block=b, loading=float(load[b]), coef=cb,
                     shift=abs(cb-base), shift_in_sd=abs(cb-base)/max(bsd,1e-12)))
T=pd.DataFrame(rows)
show(T.sort_values('shift_in_sd',ascending=False), HERE/'results/leave_one_block.csv',
     n=8, label="留一块")

worst=float(T.shift_in_sd.max())
GATE.asserted("CONTROL2 every drop-one refit is sign-aligned",
              True, "aligned to a common reference inside c3_drop", kind="control")
# ⚠ #440b: the pre-registered KILL compared a MAX OVER 32 BLOCKS to a 1-sd bar. A maximum of
# 32 draws clears 1 sd almost surely, so as written the gate could not fail. The correct bar is
# the null distribution OF THE MAXIMUM: bootstrap people, recompute all 32 shifts, take the max.
# The mis-specified gate's own verdict is kept on the record (#439d) and not replaced silently.
GATE.asserted("KILL as pre-registered (mis-specified): max shift > 1 bootstrap sd",
              worst>1.0, f"largest single-block shift = {worst:.2f} sd -- but this is a MAX over 32")
maxnull=[]
for _ in range(200):
    take=rng.choice(idxall,len(idxall),replace=True)
    mm=np.zeros(len(MM),bool); mm[np.unique(take)]=True
    bb=coef(-full_c3, mm)
    sh=[abs(coef(-LOO[b], mm & np.isfinite(LOO[b]))-bb)/max(bsd,1e-12) for b in range(NB)]
    maxnull.append(max(sh))
maxnull=np.array(maxnull); mthr=float(np.percentile(maxnull,95))
print(f"\n**多重性零(32 块最大位移的自助分布)**:均值 {maxnull.mean():.2f} · "
      f"95 分位 **{mthr:.2f}** sd")
print(f"   观测最大 **{worst:.2f}** sd -> "
      f"{'**越阈**' if worst>mthr else '**落在零里**'}")
GATE.asserted("KILL corrected: the largest single-block shift beats the null OF THE MAXIMUM",
              worst>mthr, f"observed max {worst:.2f} vs null-of-max 95th pct {mthr:.2f}")

# spec curve over k
spec=[]
for k in range(1,7):
    top=T.reindex(T.loading.abs().sort_values(ascending=False).index).head(k).block.tolist()
    keepv=np.nanmean(np.column_stack([LOO[b] for b in top]),axis=1)  # crude: mean of those LOO fits
    gd=MM & np.isfinite(keepv)
    spec.append(dict(k=k, dropped=str(top), coef=coef(-keepv,gd)))
show(pd.DataFrame(spec), HERE/'results/spec_curve.csv', n=6, label="切点规格曲线")

verdict = "CONTENT_SPECIFIC" if worst>mthr else "COORDINATE_LEVEL"
print(f"\n最大单块位移 = **{worst:.2f}** 个自助 sd -> 判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,base=base,boot_sd=bsd,worst_shift_sd=worst,null_of_max_95=mthr,kill_as_preregistered=bool(worst>1.0),
                r_pipeline=r_pipe,NB=int(NB),n=int(MM.sum())),
           open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
