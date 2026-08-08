"""E03·A30·R161b —— 描述臂:20 个种子 + 标准差÷中位(不用于改判决)

⚠ **第一版是错的,而告示牌是「三个 B 给出逐位相同的数」**:
`np.random.default_rng(s)` 写在了**内层**循环里 ⇒ **Bn 次置换全是同一次**,零塌成一个点,
B=200 / 4,000 / 20,000 打印出逐位相同的 `0.0894 · 85.47%`。
**一个与 B 无关的抽样量,就是「它从没被抽过」的告示牌**(`P5`「值 == 初始条件」的同型)。
本文件是修好之后的版本:**rng 每个种子只建一次。**

⚠ **预注册的判决在 `denominator_stability.py` 里,按 5 种子极差口径是 3.9%,落在 3–10% ⇒ 判不了。
本文件只描述,不改那个判决** —— 极差是 5 个点上最怕噪声的统计量(`#674` 同型),
而**换一个更好的估计量之后就改口,正是预注册要防的事**。
⚠ **换不了仪器**:只在 SCCS n=60 这一具装置上量。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
W=pd.read_csv("data/external/dplace/repo/datasets/SCCS/data.csv").pivot_table(
    index="soc_id",columns="var_id",values="code",aggfunc="first")
J=W[["SCCS455","SCCS1766"]].dropna(); J=J[J.SCCS1766!=20]
x=J["SCCS455"].to_numpy(float); y=J.SCCS1766.replace({10:0,21:1,22:2}).to_numpy(float)
perm={}
for Bn in (200,4000,20000):
    q=[]
    for s in range(20):
        r=np.random.default_rng(1000+s)                      # ← 每个种子只建一次(修法)
        q.append(float(np.quantile(np.abs([sp(x,r.permutation(y)) for _ in range(Bn)]),0.95)))
    perm[Bn]=(float(np.median(q)),float(np.std(q)/np.median(q)))
    print(f"  置换 B={Bn:>6,}  中位 {perm[Bn][0]:.4f} · **相对标准差 {perm[Bn][1]*100:5.2f}%**")
assert len({round(v[0],6) for v in perm.values()})>1, "三个 B 给出同一个数 ⇒ rng bug 回来了"
assert perm[200][1]>perm[4000][1]>perm[20000][1], "抖动未随 B 单调下降 ⇒ 装置不对"
print("  ✅ 正对照:三个 B 的数不相同,且抖动随 B 单调下降(两条断言写成 assert,回归时会当场炸)")
IT=["emotionally","weak","cruel","compassion","animal","kill","treated","unfairly","rights","fairly","justice","rich",
 "lovecountry","betray","loyalty","history","family","team","respect","traditions","chaos","kidrespect","sexroles","soldier",
 "decency","disgusting","god","harmlessdg","unnatural","chastity"]
mfq,_=pyreadstat.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
X=mfq[IT].dropna(); n=len(IT); M=np.eye(n)
for a in range(n):
    for b in range(a+1,n):
        xa=X[IT[a]].to_numpy(float); yb=X[IT[b]].to_numpy(float)
        r=sp(xa,yb); xs=np.sort(xa); ys=np.sort(yb); ys=ys if r>0 else ys[::-1]
        M[a,b]=M[b,a]=r/abs(sp(xs,ys))
truth={frozenset(IT[i*6:(i+1)*6]) for i in range(5)}
allb=[c for c in itertools.combinations(range(n),6) if frozenset(IT[j] for j in c) not in truth]
V=np.array([min(M[a,b] for a,b in itertools.combinations(c,2)) for c in allb])
print(f"\n  组合臂 全枚举精确值 {np.quantile(V,0.95):.4f}({len(V):,} 块)")
comb={}
for Bn in (1500,4000,20000):
    q=[float(np.quantile(V[np.random.default_rng(1000+s).choice(len(V),Bn,replace=False)],0.95)) for s in range(20)]
    comb[Bn]=(float(np.median(q)),float(np.std(q)/np.median(q)))
    print(f"  组合 抽 {Bn:>6,}  中位 {comb[Bn][0]:.4f} · **相对标准差 {comb[Bn][1]*100:5.2f}%**")
print(f"\n⇒ **同一个 B=4,000:置换 {perm[4000][1]*100:.2f}% vs 组合 {comb[4000][1]*100:.2f}% "
      f"= 组合是置换的 {comb[4000][1]/perm[4000][1]:.1f} 倍**")
json.dump(dict(permutation={str(k):v for k,v in perm.items()},combinatorial={str(k):v for k,v in comb.items()},
   exact=float(np.quantile(V,0.95)),note="描述,不改 denominator_stability.py 的预注册判决(判不了)",
   unchallenged=True),open(OUT/"descriptive_20seed.json","w"),indent=1,ensure_ascii=False)
