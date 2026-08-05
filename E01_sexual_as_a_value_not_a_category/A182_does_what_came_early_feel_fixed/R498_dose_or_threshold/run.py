import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #451/#452 established that people whose interests arrived earlier feel less able to change
   them. That is one coefficient. Is it a DOSE -- more of it the earlier you go -- or a
   THRESHOLD, where only the earliest stratum differs?

Worlds
  A  dose      : the adjusted outcome rises monotonically across the four onset strata. Then
     "the earlier it arrived, the less movable it feels" is a gradient, which is a much
     stronger statement than a single slope.
  B  threshold : only the earliest stratum stands apart. Then the sentence changes: not
     "earlier is more fixed" but "before some point it is fixed", which is a different claim
     about people and one the page cannot currently make.

⚠ THE CUTS ARE NOT MINE: the page already carries a four-way onset split with stratum medians
   **12.5 / 14.7 / 16.7 / 19.9**. Inventing new cuts would be inventing a new hypothesis. This
   round reproduces that split and reports the stratum medians as a CONTROL that it did.
Method: the outcome is residualised on #451's control set (S, -(five-item), c3-, category
   count, current age) POOLED, then averaged within stratum -- so strata differ only in onset.
Monotonicity statistic: Spearman rho between stratum index and adjusted mean (4 points).
NULL : permute stratum labels across people, recompute -- names the world "onset stratum
   carries nothing", which is the world to exclude here (#450b).
CONTROL2: the stratum sizes must be roughly equal, else "quartile" is the wrong word.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R498 dose or threshold")
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
BEL=np.asarray(OUT['能不能改'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(_EARLY)
n=int(MM.sum()); print(f"n = **{n:,}**")

q=np.quantile(_EARLY[MM],[0.25,0.5,0.75])
strat=np.full(len(MM),-1); strat[MM]=np.digitize(_EARLY[MM],q)
meds=[float(np.median(_EARLY[MM&(strat==k)])) for k in range(4)]
sizes=[int((MM&(strat==k)).sum()) for k in range(4)]
print(f"四层中位起始年龄 = **{[round(x,1) for x in meds]}** "
      f"(页面既有:12.5 / 14.7 / 16.7 / 19.9)· 层大小 {sizes}")
PAGE=[12.5,14.7,16.7,19.9]
GATE.asserted("CONTROL the four strata reproduce the page's existing split",
              all(abs(a-b)<0.8 for a,b in zip(meds,PAGE)),
              f"medians {[round(x,1) for x in meds]} vs page {PAGE}", kind="control")
GATE.asserted("CONTROL2 the strata are quartiles (sizes within 5%)",
              (max(sizes)-min(sizes))/np.mean(sizes)<0.05, f"sizes {sizes}", kind="control")

# 结局对控制集残差化(合并),再按层取均值 -> 层与层只差起始年龄
X=np.column_stack([np.ones(n), z(A,MM), z(Bv,MM), z(C3,MM), z(ncat,MM), z(AGE,MM)])
y=z(BEL,MM)
res=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
sub=strat[MM]
adj=np.array([res[sub==k].mean() for k in range(4)])
se =np.array([res[sub==k].std()/np.sqrt((sub==k).sum()) for k in range(4)])
from scipy.stats import spearmanr
rho=float(spearmanr(np.arange(4),adj).statistic)
T=pd.DataFrame([dict(stratum=k+1, median_onset=round(meds[k],1), n=sizes[k],
                     adj_mean=adj[k], lo=adj[k]-1.96*se[k], hi=adj[k]+1.96*se[k])
                for k in range(4)])
show(T, HERE/'results/strata.csv', n=4, label="四层")
print(f"\n单调性(层序 × 调整后均值的 Spearman rho)= **{rho:+.3f}**")

rg=np.random.default_rng(71)
nul=np.array([spearmanr(np.arange(4),
              np.array([res[p==k].mean() for k in range(4)])).statistic
              for p in (rg.permutation(sub) for _ in range(4000))])
p_hi=float(np.mean(nul>=rho))
print(f"零(打乱层标签)rho 分布:中位 {np.median(nul):+.3f} · "
      f"95 分位 **{np.percentile(nul,95):+.3f}** · **观测的零上侧比例 = {p_hi:.4f}**")

# 阈值检验:只有最早一层不同?
d_first=float(adj[0]-adj[1:].mean()); d_rest=float(np.ptp(adj[1:]))
print(f"最早一层 vs 其余三层的差 = **{d_first:+.4f}** · 其余三层内部的极差 = **{d_rest:.4f}**")
mono = (rho>=0.8) and (p_hi<0.05)
GATE.asserted("KILL the relation is a dose (monotone across all four strata)", mono,
              f"rho = {rho:+.3f}, null upper-tail p = {p_hi:.4f}")
verdict = "DOSE" if mono else ("THRESHOLD" if abs(d_first)>d_rest else "NEITHER_RESOLVED")
print(f"\n判决 = {verdict}")
json.dump(dict(verdict=verdict,n=n,medians=meds,sizes=sizes,adj=adj.tolist(),
               rho=rho,p_hi=p_hi,d_first=d_first,d_rest=d_rest),
          open(HERE/'results/verdict.json','w'),indent=1)
print(GATE.verdict())
