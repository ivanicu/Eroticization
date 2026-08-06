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
import numpy as _np
from scipy.stats import spearmanr

__all__=["spearman","ceiling","pairmat","weakest_greedy","weakest_optimal",
         "opt_batch","pooled_median","aligned_pooled_median","median_of_group_medians"]

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

def aligned_pooled_median(M,ix,optimal=False):
    """**先定翻向,再取块内全部对的中位。**

    ⚠ **为什么必须有这个函数**:`pooled_median` 不做符号对齐,而 NSFG 的 `okcohab`
    题干本就反向("should **not** live together unless married")⇒ 不对齐时块内中位是 **−0.0332**,
    对齐后是 **+0.4257**。`#650` 用的是对齐版,而 `#734` 第一次重跑时漏了这一步,
    **正对照当场开火** —— 这是 `#718`/`#720` 之后**第三次**漏掉同一步。
    ⇒ **凡是要和 `#650`/`#653` 系的数比较的中位,都必须走这里,不能走 `pooled_median`。**

    `optimal=False` 走贪心(行和定号,与 `#653`/`#650` 同路);`optimal=True` 穷举 2^(k-1)
    取**使中位最大**的那一组符号(注意:与 `weakest_optimal` 的目标函数不同,那个最大化的是 min)。
    """
    ix=list(ix); k=len(ix); sub=M[_np.ix_(ix,ix)]
    def med(s):
        return float(_np.median([s[a]*s[b]*sub[a,b] for a,b in itertools.combinations(range(k),2)]))
    if not optimal:
        s=_np.where(_np.nansum(sub,axis=1)>=0,1,-1)
        return med(s)
    best=-2.0
    for bits in range(1<<(k-1)):
        s=_np.array([1]+[1 if (bits>>t)&1==0 else -1 for t in range(k-1)])
        best=max(best,med(s))
    return best

def median_of_group_medians(M,groups):
    """**先每组取中位,再对组取中位。** `#720` 用的是这一种,而它拿去比了 `pooled_median`。

    ⚠ **两者不是同一个数**(`#721` 实测:0.5645 对 0.5856)。名字里带着口径就是为了不再比错。
    """
    return float(np.median([pooled_median(M,g) for g in groups]))


# ============================================================================
# #759 · 符号方向 —— 方向未知时 raise,不默认
# ============================================================================
# ⚠ 动机:符号这一族已经六次(`#718`·`#720`·`#734`·`#751`·`#756`·`#758`)。
#   而 `#755` 刚证明「为语义缺陷造仪器」会失败,所以**先分类,只造能造的那半**:
#     机械的(可机读):`#718`/`#720`/`#734`(合成时忘了反向题)· `#751`(两套电池方向相反)
#                      —— **方向就写在值标签里**;
#     语义的(不可机读):`#756`(恒等式管的是哪个统计量)—— 只能靠约定,不在此处。
#
# ⚠ P6 代理账:
#   PROPERTY   这一列的高值代表的是我以为的那一端
#   PROXY      该列**最后一个**值标签里出现的词
#   IMPLICATION 只有一个方向可靠:**词表匹配不上 -> 我确实不知道方向**(可靠,于是 raise)。
#              匹配上**不**证明我用对了 —— 只证明标签里有那个词。**从不认证方向正确。**
#   SAFE SIDE  不确定就 **raise**,绝不猜一个默认值。

_POLES = {
    # 高值端的词 -> 这一端叫什么
    "permissive": ("not wrong at all", "not wrong", "not at all wrong", "never wrong"),
    "strict":     ("seriously wrong", "always wrong", "very wrong"),
    "unimportant":("least important", "not at all important", "not important"),
    "important":  ("most important", "very important", "extremely important"),
    "agree":      ("strongly agree", "agree strongly"),
    "disagree":   ("strongly disagree", "disagree strongly"),
}

def label_pole(categories):
    """给定按码序排列的值标签,回答**高值那一端叫什么**。认不出就 raise。"""
    if not categories: raise ValueError("没有值标签 —— 方向不可知,不许默认")
    hi = str(categories[-1]).strip().lower()
    for pole, words in _POLES.items():
        if any(w in hi for w in words): return pole
    raise ValueError(f"最高码的标签 {categories[-1]!r} 不在词表里 —— **方向未知,去看值标签,不要猜**")

def aligned(categories_by_col, high_means):
    """检查每一列的高值端是否都是 `high_means`;返回需要取负的列名集合。

    `categories_by_col`: {列名: 按码序排列的值标签}
    `high_means`: 期望的高值端(`_POLES` 的键之一)

    ⚠ 它**不改数据**,只回答「哪几列要翻」—— 翻不翻由调用者写出来,
      因为 `#734` 的教训是**我忘了翻**,不是我翻错了方向。
    ⚠ 任何一列认不出方向就 raise,**整个调用失败**,而不是跳过那一列。
    """
    if high_means not in _POLES:
        raise ValueError(f"high_means={high_means!r} 不是已知的极点:{sorted(_POLES)}")
    flip = set()
    for col, cats in categories_by_col.items():
        pole = label_pole(cats)          # 认不出就在这里 raise
        if pole != high_means: flip.add(col)
    return flip


# ============================================================================
# #766 · 缺失码 —— 笼统过滤会静默删掉最大的一档
# ============================================================================
# ⚠ 动机(`#765`,代价是一个已发表的数):`.where(x>0)` 是我加的笼统缺失过滤,
#   而 GSS 的 `.dta` **本来就把 DK/NA 设成缺失了** —— 于是它删掉的是真数据:
#   `attend` 码 0 = "never",**n=14,883,最大的一档**,整档消失,而那正是宗教分析的参照组。
#
# ⚠ 与 `#755` 那个被否掉的 lint 的差别,这一条是它能成立的原因:
#   **真值来源是可机读的**(值标签),而不是我的意图。`#755` 要判「判词比错了对象」——语义,造不出来。
#
# ⚠ P6 代理账:
#   PROPERTY   我保留的码集合 == 这一列真正有效的码集合
#   PROXY      被排除的码里,有没有**带标签且样本量可观**的
#   IMPLICATION 只有一个方向可靠:**排除了一个带标签的大档 -> 我确实删了真数据**(可靠)。
#              反过来不成立:**没排除大档不证明我的码集合对** —— 小档同样可能是真数据。
#   SAFE SIDE  只报「你删掉了带标签的档」;**从不认证「这个范围是对的」。**

def labelled_codes(dta_path, col, encoding="latin1"):
    """读出一列**带标签的码**与各自的 n。真值来源是值标签,不是我的记忆。"""
    import pandas as pd
    num = pd.read_stata(dta_path, columns=[col], convert_categoricals=False)[col]
    lab = pd.read_stata(dta_path, columns=[col], convert_categoricals=True)[col]
    j = pd.DataFrame({"code": pd.to_numeric(num, errors="coerce"), "label": lab.astype(str)}).dropna()
    g = j.groupby(["code", "label"]).size().reset_index(name="n").sort_values("code")
    return [(float(r.code), str(r.label), int(r.n)) for r in g.itertuples()]

def check_kept_codes(dta_path, col, keep, min_share=0.01, encoding="latin1"):
    """`keep` 是一个判定函数或 (lo,hi);报出**被排除的带标签档**及其占比。

    返回 (dropped, total) —— `dropped` 是 [(码, 标签, n, 占比)],按 n 降序。
    ⚠ 不 raise、不改数据:**它只把「你正在删什么」摆到面前**,决定权在调用者。
      理由与 `#759` 的 `aligned()` 相同 —— `#734`/`#765` 的教训都是**我没看见**,不是我看错。
    ⚠ `min_share` 只影响**打印时的醒目程度**,不影响返回值;默认 1%。
    """
    rows = labelled_codes(dta_path, col, encoding)
    total = sum(n for _, _, n in rows)
    pred = keep if callable(keep) else (lambda c, lo=keep[0], hi=keep[1]: lo <= c <= hi)
    dropped = [(c, l, n, n / total) for c, l, n in rows if not pred(c)]
    return sorted(dropped, key=lambda t: -t[2]), total
