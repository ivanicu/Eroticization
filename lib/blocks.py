"""块统计的共享实现 —— 由 `#722` 从 `#653`/`#715`–`#721` 的重复代码提出来。

**为什么存在(三个证据,都在两轮之内):**
- `#718`:贪心对齐把 SCCS 同对象块渲染成 −0.4164,最优符号是 −0.0138 ——
  **「互斥」与「无关」是两句完全不同的心理学。**
- `#720`:一具全新仪器上**我自己又忘了用**最优符号,三个面因此是负的。
- `#721`:**「逐组中位的中位」被拿去比「合并中位」,正对照因此蒙对。**

⇒ 本模块把三件事定死:**符号取最优、口径必须显式命名、天花板归一是默认。**
每个函数的名字里就带着口径,**调用处读不出口径的函数一律不提供**。
"""
import itertools
import numpy as np
from scipy.stats import spearmanr

__all__=["spearman","ceiling","pairmat","weakest_greedy","weakest_optimal",
         "opt_batch","pooled_median","median_of_group_medians"]

def spearman(a,b):
    return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)

def ceiling(a,b):
    """同调耦合下这一对能达到的最大 |ρ| —— 边际决定的天花板。"""
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    r=spearman(a,b)
    return abs(spearman(x,y if r>0 else y[::-1]))

def pairmat(F,items,year=None,floor=150,normalise=True):
    """归一(或生)对相关矩阵。

    ⚠ **删失规则:逐对删失**(每对用它自己的完整个案)。
    `#653` 的 `weakest()` 用的是**块内整行删失** —— `#718` 实测两者在那 9 块上最大差 **0.0044**,
    但**这是那一批数据上的量,不是保证**。需要整行删失时,调用方先 `F[items].dropna()` 再传进来。
    """
    n=len(items); M=np.full((n,n),np.nan)
    groups=[(None,F)] if year is None else list(F.groupby(year))
    for a in range(n):
        for b in range(a+1,n):
            per=[]
            for _,g in groups:
                m=g[[items[a],items[b]]].dropna()
                if len(m)<floor or m[items[a]].nunique()<2 or m[items[b]].nunique()<2: continue
                r=spearman(m[items[a]],m[items[b]])
                if not np.isfinite(r) or r==0: continue
                if not normalise: per.append(r); continue
                c=ceiling(m[items[a]],m[items[b]])
                if np.isfinite(c) and c>1e-9: per.append(r/c)
            if per: M[a,b]=M[b,a]=float(np.median(per))
    np.fill_diagonal(M,1.0)
    return M

def weakest_greedy(M,ix):
    """**贪心**:按行和定翻向。⚠ `#718` 证明它会把「无关」渲染成「负关系」。保留只为复现旧数。"""
    sub=M[np.ix_(ix,ix)]; s=np.where(np.nansum(sub,axis=1)>=0,1,-1)
    return float(min(s[a]*s[b]*sub[a,b] for a,b in itertools.combinations(range(len(ix)),2)))

def _sign_table(k):
    pr=list(itertools.combinations(range(k),2))
    S=[[1]+[1 if (b>>t)&1==0 else -1 for t in range(k-1)] for b in range(1<<(k-1))]
    return pr,np.array([[s[a]*s[b] for a,b in pr] for s in S])

def opt_batch(M,blocks):
    """向量化最优符号:blocks (B,k) -> 每块 `max over 2^(k-1) signs of min over C(k,2) pairs`。

    ⚠ 2^(k-1) 指数增长:k=7 → 64 可算,k≈20 会爆。**调用方负责 k 的上界。**
    """
    blocks=np.asarray(blocks); B,k=blocks.shape
    pr,SG=_sign_table(k)
    V=np.stack([M[blocks[:,a],blocks[:,b]] for a,b in pr],axis=1)
    return np.max(np.min(V[:,None,:]*SG[None,:,:],axis=2),axis=1)

def weakest_optimal(M,ix):
    """单块的最优符号最弱一环。**严格 ≥ `weakest_greedy`**(贪心那种指派也在枚举里)。"""
    return float(opt_batch(M,np.asarray([list(ix)]))[0])

def pooled_median(M,ix):
    """**把这一组的全部对倒进一个池再取中位。** `#542` 用的是这一种。"""
    return float(np.median([M[a,b] for a,b in itertools.combinations(ix,2)]))

def median_of_group_medians(M,groups):
    """**先每组取中位,再对组取中位。** `#720` 用的是这一种,而它拿去比了 `pooled_median`。

    ⚠ **两者不是同一个数**(`#721` 实测:0.5645 对 0.5856)。名字里带着口径就是为了不再比错。
    """
    return float(np.median([pooled_median(M,g) for g in groups]))
