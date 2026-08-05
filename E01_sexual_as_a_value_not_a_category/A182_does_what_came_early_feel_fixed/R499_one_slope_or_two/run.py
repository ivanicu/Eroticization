import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #454b found a clean gradient from onset stratum to feeling able to change. The page also
   carries "earlier onset -> more shame". Is that a gradient too, and if so, is it the SAME
   gradient -- one antecedent pushing both -- or two different shapes?

Worlds
  A  one antecedent : both outcomes are monotone across the four strata and their steepness is
     comparable. Then "early -> shame" and "early -> feels unmovable" are two readings of one
     thing, and the page should say so once rather than twice.
  B  two paths      : the shapes differ (one monotone and one not, or steepness differing by
     more than its own spread). Then they cannot be treated as the same fact, and the page's
     two arcs stay two.

Everything is held identical to #454 on purpose: the same four strata (the page's own cuts),
the same control set, the same monotonicity null. Only the outcome changes -- so a difference
in shape cannot be a difference in method.
⚠ Both outcomes are z-scored, so the two steepnesses are on one scale and can be compared
directly; the DIFFERENCE gets its own person-level bootstrap, because a difference of two
noisy slopes is noisier than either (#413b).
CONTROL : the changeability arm must reproduce #454's numbers exactly -- if it does not, the
          pipeline changed and the comparison is void.
CONTROL2: the shame arm's overall direction must match the page (earlier -> more shame).
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R499 one slope or two")
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_EARLY=np.array(EARLY,dtype=float).copy()
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY)
n=int(MM.sum())
q=np.quantile(_EARLY[MM],[0.25,0.5,0.75]); sub=np.digitize(_EARLY[MM],q)
X=np.column_stack([np.ones(n), z(A,MM), z(Bv,MM), z(C3,MM), z(ncat,MM), z(AGE,MM)])

def profile(y, idx=None, s=None):
    yy=z(np.asarray(y,dtype=float),MM); Xs=X
    if idx is not None: yy=yy[idx]; Xs=X[idx]
    r=yy-Xs@np.linalg.lstsq(Xs,yy,rcond=None)[0]
    ss=sub if s is None else s
    return np.array([r[ss==k].mean() for k in range(4)])

OUTS={'能不能改':OUT['能不能改'], '羞耻':OUT['羞耻']}
prof={k:profile(v) for k,v in OUTS.items()}
rho={k:float(spearmanr(np.arange(4),v).statistic) for k,v in prof.items()}
# 陡度 = 末层 − 首层(两个结局都已 z 化 -> 同尺度)
steep={k:float(v[3]-v[0]) for k,v in prof.items()}
rg=np.random.default_rng(83)
nul={k:np.array([float(spearmanr(np.arange(4),profile(OUTS[k],s=p)).statistic)
                 for p in (rg.permutation(sub) for _ in range(3000))]) for k in OUTS}
T=pd.DataFrame([dict(outcome=k, s1=prof[k][0], s2=prof[k][1], s3=prof[k][2], s4=prof[k][3],
                     rho=rho[k], null_p95=float(np.percentile(nul[k],95)),
                     steepness=steep[k]) for k in OUTS])
show(T, HERE/'results/two_profiles.csv', n=4, label="两条斜坡")

GATE.asserted("CONTROL the changeability arm reproduces #454 exactly",
              abs(prof['能不能改'][0]-(-0.059692))<1e-4 and abs(rho['能不能改']-1.0)<1e-9,
              f"s1 = {prof['能不能改'][0]:+.6f} vs #454 -0.059692; rho = {rho['能不能改']:+.3f}",
              kind="control")
GATE.asserted("CONTROL2 the shame arm points the way the page says (earlier -> more shame)",
              steep['羞耻']<0, f"shame steepness (late - early) = {steep['羞耻']:+.4f}", kind="control")

# 陡度差,及其自己的自助
d_obs=abs(steep['能不能改'])-abs(steep['羞耻'])
idxall=np.arange(n); rg2=np.random.default_rng(97); bs=[]
for _ in range(400):
    take=rg2.choice(idxall,n,replace=True)
    p1=profile(OUTS['能不能改'],idx=take,s=sub[take])
    p2=profile(OUTS['羞耻'],idx=take,s=sub[take])
    bs.append(abs(p1[3]-p1[0])-abs(p2[3]-p2[0]))
bs=np.array(bs); lo,hi=np.percentile(bs,[2.5,97.5])
print(f"\n陡度(末层 − 首层,同尺度):能不能改 **{steep['能不能改']:+.4f}** · "
      f"羞耻 **{steep['羞耻']:+.4f}**")
print(f"|陡度| 之差 = **{d_obs:+.4f}** · 自助 95% **[{lo:+.4f}, {hi:+.4f}]** "
      f"-> {'**不含 0:两条不一样陡**' if (lo>0)==(hi>0) else '**含 0:陡度分不开**'}")
# ⚠ #455b: written as `rho[k] >= 0.8`, which fails on shame's PERFECT but NEGATIVE monotonicity
# (-1.0) -- and "earlier -> more shame" is exactly what a negative rho means here. Sixth gate
# this session that tested something other than its own sentence (#433a #439d #440b #444a #451b).
# The word in the sentence was "monotone", which is |rho|.
# ⚠ #455c: `|rho| > two-sided p95` CANNOT FIRE at four strata -- |rho| = 1 arises in 2 of 4! = 8.3%
# of random orderings, so the two-sided 95th percentile IS 1.0. Both directions here were
# pre-specified by the page BEFORE this round (earlier -> more shame; earlier -> less able to
# change), so the correct test is ONE-SIDED in each outcome's own predicted direction.
PRED={'能不能改':+1, '羞耻':-1}          # 方向由页面既有结论事先给定
one_sided={k: float(np.mean(np.sign(PRED[k])*nul[k] >= np.sign(PRED[k])*rho[k])) for k in OUTS}
print(f"⚠ 双侧不可用(4 层下 |rho|=1 的随机概率 = 2/24 = 8.3%);"
      f"两个方向都由页面**事先**给定 -> 单侧 p = "
      f"{ {k: round(one_sided[k],4) for k in OUTS} }")
both_mono = all(abs(rho[k])>=0.8 and np.sign(rho[k])==PRED[k] and one_sided[k]<0.05 for k in OUTS)
mde_diff = float(hi-lo)
print(f"⚠ 陡度差的区间宽 **{mde_diff:.4f}**,而两条陡度本身只有 "
      f"**{abs(steep['能不能改']):.4f}** 与 **{abs(steep['羞耻']):.4f}** —— "
      f"**区间比被比较的量还宽 -> 「分不开」是低功率,不是「一样」**")
same_steep = not ((lo>0)==(hi>0))
GATE.asserted("KILL one antecedent (both monotone AND steepness indistinguishable)",
              both_mono and same_steep,
              f"rho {[round(rho[k],2) for k in OUTS]}; steepness diff CI [{lo:+.4f},{hi:+.4f}]")
verdict = "ONE_ANTECEDENT" if (both_mono and same_steep) else "TWO_PATHS"
print(f"\n单调性:{ {k:round(rho[k],3) for k in OUTS} } -> 判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,profiles={k:prof[k].tolist() for k in OUTS},
               rho=rho,steepness=steep,diff=d_obs,boot=[lo,hi],
               null_p95={k:float(np.percentile(np.abs(nul[k]),95)) for k in OUTS}),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
