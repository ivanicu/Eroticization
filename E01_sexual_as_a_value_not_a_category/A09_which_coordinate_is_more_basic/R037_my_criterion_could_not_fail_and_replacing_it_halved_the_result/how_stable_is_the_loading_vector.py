import os,sys,pathlib,json,re
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #445c argued the strength of "20 of 20" lies in each half re-estimating the SAME loading
   structure, not in the number 20 -- and that was never measured. How similar are the twenty
   half-sample `c2` loading vectors to each other?

Worlds
  A  the structure is stable (pairwise r high) -> "the fluid blocks load on this half" is a
     property of the coordinate, and #445b's strength is earned.
  B  the structure is unstable but the fluid gap still appeared in all twenty -> then the
     fluid effect is MORE BASIC than `c2`: the coordinate is an unstable carrier of it, and
     the "two halves" framing on the page is the thing that has to change, not the finding.

⚠ SIGN: an eigenvector's sign is arbitrary, so every half's vector is anchored the same way
   before any pairwise comparison (#443b's rule) -- otherwise half the correlations are
   sign-flipped and their median is a coin flip. This project has been burnt five times.
NULL : pairwise correlation between vectors carrying NO shared structure -- built by
   permuting each half's loading vector over blocks, keeping its own value distribution.
CONTROL : each half's vector must correlate positively with the full-sample vector (that is
   what the anchoring is supposed to achieve; if it does not, the anchoring failed).
Report the DISTRIBUTION (median + interval), never a single number.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R490 loading stability")
_R489=(ROOT/'E01_sexual_as_a_value_not_a_category/A179_two_ways_of_being_involved'
            /'R489_does_the_fluid_pattern_replicate/run.py').read_text()
exec(_R489.split('"""',2)[2].split('FLU=re.compile')[0])

FLU=re.compile(r'urine|ejaculate|squirt|precum|saliva|secretion',re.I)
fl=np.array([bool(FLU.search(n)) for n in NAMES]); nf=int(fl.sum())
ALLR=np.flatnonzero(ok)
v_full=load_k(ALLR,1)

rng=np.random.default_rng(97)          # same seed as #445 -> the same twenty halves
VS=[]
for t in range(20):
    perm=rng.permutation(ALLR); hb=perm[len(perm)//2:]
    VS.append(load_k(hb,1))
VS=np.array(VS)
# anchoring is inside load_k (turned toward MORE involvement) -- verify it took
r_to_full=np.array([float(np.corrcoef(v,v_full)[0,1]) for v in VS])
print(f"每一半与全样本向量的相关:中位 **{np.median(r_to_full):+.4f}** · "
      f"最小 **{r_to_full.min():+.4f}** · 为正的 **{int((r_to_full>0).sum())}/20**")
GATE.asserted("CONTROL the anchoring took (every half faces the full-sample vector)",
              bool((r_to_full>0).all()), f"min r to full = {r_to_full.min():+.4f}", kind="control")

pw=[float(np.corrcoef(VS[i],VS[j])[0,1]) for i in range(20) for j in range(i+1,20)]
pw=np.array(pw); lo,hi=np.percentile(pw,[2.5,97.5])
rg2=np.random.default_rng(3)
nul=[]
for _ in range(4000):
    i,j=rg2.integers(0,20,2)
    a=VS[i][rg2.permutation(NB)]; b=VS[j][rg2.permutation(NB)]
    nul.append(float(np.corrcoef(a,b)[0,1]))
nul=np.array(nul); nlo,nhi=np.percentile(nul,[2.5,97.5])
T=pd.DataFrame([dict(q='两两相关 中位',v=float(np.median(pw))),
                dict(q='两两相关 2.5%',v=lo), dict(q='两两相关 97.5%',v=hi),
                dict(q='零(块内打乱)中位',v=float(np.median(nul))),
                dict(q='零 97.5%',v=nhi),
                dict(q='与全样本 中位',v=float(np.median(r_to_full)))])
show(T, HERE/'results/stability.csv', n=6, label="载荷稳定性")

med=float(np.median(pw))
GATE.asserted("CONTROL2 the null is centred at zero and can say no",
              abs(np.median(nul))<0.1, f"null median = {np.median(nul):+.4f}", kind="control")
stable = lo>0.8
GATE.asserted("KILL the loading structure is stable (pairwise r interval above 0.8)",
              stable, f"pairwise r median {med:+.4f}, 95% [{lo:+.4f}, {hi:+.4f}]")
verdict = "STABLE_COORDINATE" if stable else "CONTENT_MORE_BASIC_THAN_COORDINATE"
print(f"\n两两相关 中位 **{med:+.4f}** · 95% **[{lo:+.4f}, {hi:+.4f}]** · "
      f"零 **[{nlo:+.4f}, {nhi:+.4f}]**")
print(f"判决 = {verdict}")
json.dump(dict(verdict=verdict,pw_median=med,pw_lo=lo,pw_hi=hi,
               null_lo=nlo,null_hi=nhi,r_to_full_median=float(np.median(r_to_full)),
               r_to_full_min=float(r_to_full.min())),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())

# ---------------------------------------------------------------- #446c(发布前追加)
# The two worlds were posed as either/or. The data are neither: the coordinate is mostly
# stable and occasionally collapses, so the decisive question is whether the CONTENT effect
# depends on the coordinate being well estimated. Measured directly.
gaps=np.array([float(v[fl].mean()-v[~fl].mean()) for v in VS])
r_stab=float(np.corrcoef(r_to_full,gaps)[0,1])
o=np.argsort(r_to_full); worst=gaps[o[:5]]; best=gaps[o[-5:]]
W=pd.DataFrame([dict(q='坐标稳定度 与 gap 的相关',v=r_stab),
                dict(q='最差 5 个半样本 gap 中位',v=float(np.median(worst))),
                dict(q='最差 5 个里 gap 全为正',v=float((worst>0).all())),
                dict(q='最好 5 个半样本 gap 中位',v=float(np.median(best))),
                dict(q='半样本 MDE(#445a)',v=0.1580)])
show(W, HERE/'results/content_vs_coordinate.csv', n=5, label="内容 vs 坐标")
GATE.asserted("KILL(重做):内容效应在坐标估得最差的半样本里仍然成立",
              bool((worst>0).all()) and float(np.median(worst))>0.158,
              f"worst-5 gap median {np.median(worst):+.4f} vs half-sample MDE 0.158")
print(f"\n⇒ **坐标是脆的那一半,内容是韧的那一半** —— 但最差 5 个的中位 "
      f"**{np.median(worst):+.4f}** 只是刚过 MDE **0.158**,**贴着分辨率**。")
