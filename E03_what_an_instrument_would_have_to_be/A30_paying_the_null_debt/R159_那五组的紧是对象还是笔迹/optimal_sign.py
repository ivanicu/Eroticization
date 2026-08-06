"""E03·A30·R159b —— 最优符号指派:把「贪心对齐失败」与「数据里真有负关系」分开

**为什么必须有这一步。** `align()` 是**贪心**的:每一项按它与块均值的相关翻向。
k=4 的同对象块在贪心下中位 **−0.4164**,读起来像「**手段之间互相排斥**」;
而 k=4 只有 **2^(k−1) = 8** 种本质不同的符号指派,**穷举到最优是 −0.0138** ——
**那个负号是贪心步造的,不是数据里的。「互斥」与「无关」是两句完全不同的心理学。**

⚠ **换不了仪器**:同 `sccs_exact_null.py`。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
S=pathlib.Path("data/external/dplace/repo/datasets/SCCS")
W=pd.read_csv(S/"data.csv").pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
TECH={"以身作则":429,"讲课":437,"体罚":453,"放任":465,"疼爱":469}
GRP={k:[f"SCCS{b+i}" for i in range(4)] for k,b in TECH.items()}
POOL=[v for L in GRP.values() for v in L]; idx={c:i for i,c in enumerate(POOL)}
n=len(POOL); M=np.full((n,n),np.nan)
for a in range(n):
    for b in range(a+1,n):
        m=W[[POOL[a],POOL[b]]].dropna()
        if len(m)<30: continue
        x=m[POOL[a]].to_numpy(float); y=m[POOL[b]].to_numpy(float)
        r=sp(x,y); xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        c=sp(xs,ys)
        if abs(c)>1e-9: M[a,b]=M[b,a]=r/abs(c)
def best(cols):
    ii=[idx[c] for c in cols]; bv=-2.0
    for bits in range(1<<(len(ii)-1)):
        s=[1]+[1 if (bits>>t)&1==0 else -1 for t in range(len(ii)-1)]
        bv=max(bv,min(s[a]*s[b]*M[ii[a],ii[b]] for a,b in itertools.combinations(range(len(ii)),2)))
    return bv
tru={k:best(v) for k,v in GRP.items()}
print("真块(同一手段 × 四对象),最优符号:"+" · ".join(f"{k} {v:+.4f}" for k,v in tru.items()))
st=[best(list(c)) for i in range(4) for c in itertools.combinations([GRP[k][i] for k in TECH],4)]
print(f"同对象块(20 块):中位 **{np.median(st):+.4f}** · 范围 [{min(st):+.4f}, {max(st):+.4f}]"
      f"  ← 贪心下是 −0.4164")
truset={frozenset(v) for v in GRP.values()}
allb=[c for c in itertools.combinations(POOL,4) if frozenset(c) not in truset]
v=np.array([best(list(c)) for c in allb]); v=v[np.isfinite(v)]
q95=float(np.quantile(v,0.95))
print(f"零(全枚举 {v.size:,} 块):**零的 95% 分位 {q95:+.4f}** · 中位 {np.median(v):+.4f} · 最大 {v.max():+.4f}")
print(f"\n真块中位 ÷ 零 = **{np.median(list(tru.values()))/q95:.2f}×** · "
      f"同对象块中位 ÷ 零 = **{np.median(st)/q95:.2f}×** ⇒ "
      f"{'**任何符号指派都救不回来 —— 不是贪心的锅**' if np.median(st)<q95 else '⚠ 最优符号救回来了'}")
json.dump(dict(true_optimal=tru,same_target_median=float(np.median(st)),
   same_target_range=[float(min(st)),float(max(st))],null_q95=q95,n_blocks=int(v.size),
   greedy_same_target_median=-0.4164,
   note="零的种类 = 同一编码团队 Barry 1977 的 20 变量池、同样 k=4、最优符号,只打散分组",
   unchallenged=True),open(OUT/"optimal_sign.json","w"),indent=1,ensure_ascii=False)
