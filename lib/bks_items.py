"""BKS 的 likert 题集合 —— 一条**鸭子类型发现规则加一个例外**,被 `import`,不被拷贝。

⚠⚠ **这个模块存在的理由,是 `#869` 的审计实测出来的一次「已经发生了一半的灭绝」:**

同一段规则在 `E01` 里出现了 **38 份逐字相同的拷贝**:
```
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']              # ⚠ 剔除:它是协变量
```
**而第二行只出现在 32 份里** —— **8 份只有发现规则,没有那个例外。**

**而那个例外不是可有可无的**:在真实的 `data/raw/BKSPublic.csv` 上,
鸭子类型规则选出 **20 列**,而 **`biomale` 就在里面**(取值 {0,1} ⊂ {−3..3},非缺 15,503 > 10,000)。
⇒ **那 8 份把一个性别协变量当成了 likert 题,占它们题集的 5%(1/20)。**

**这正是 `colcom` 的形状,只是它已经死了一半:**
一条规则活在被逐轮拷贝的两行里,拷贝到第 33 份时有人漏了第二行,而**没有任何东西会报错**。
⇒ 搬到这里:**被 `import`,不被拷贝**。

⚠ **本模块不回改那 8 轮的结果** —— 逐轮重跑与重估影响是另一件工作,
**已登记为 `#869`①,而它的影响至今 UNASSESSED,不许写成「影响很小」。**
"""
import pandas as pd

# 鸭子类型发现规则的三个参数,与 38 份拷贝逐字相同(`#869` 用 `uniq -c` 查过:38/38 一致)
LIKERT_VALUES = {-3., -2., -1., 0., 1., 2., 3.}
MIN_NONNULL = 10000
# 通过鸭子类型但**不是** likert 题的列。每一条都要写清**为什么**,否则它就是下一个 `biomale`。
NOT_ITEMS = {
    "biomale": "性别协变量:取值 {0,1} ⊂ {−3..3} 且非缺 15,503,**碰巧通过**鸭子类型",
}


def likert_columns(d, values=None, min_nonnull=None, exclude=None):
    """返回 BKS 的 likert 题列名。

    ⚠ **例外是默认生效的。** 想看「不剔除会怎样」必须**显式**传 `exclude=set()`,
    因为那是 `#869` 查出的那 8 轮实际做的事,而**默认值不该是那个**。
    """
    # ⚠ **空总体不许静默通过**(`realstat`:a gate that reports success having examined nothing）。
    # 攻击向量 ⑦ 抓到的:空 DataFrame 原本安静地返回 `[]`,而一个返回空表的题集
    # 会让下游每一个统计量变成 nan 或 0,**而没有任何一行会说出原因**。
    if d is None or getattr(d, "shape", (0, 0))[1] == 0:
        raise ValueError("likert_columns 收到一个没有列的框 —— 空总体不许当作「没有 likert 题」")
    vals = LIKERT_VALUES if values is None else values
    mn = MIN_NONNULL if min_nonnull is None else min_nonnull
    ex = NOT_ITEMS.keys() if exclude is None else exclude
    lik = [c for c in d.columns
           if d[c].dtype != object
           and set(pd.Series(d[c]).dropna().unique()) <= vals
           and d[c].notna().sum() > mn]
    return [c for c in lik if c not in ex]


def audit(d):
    """把「鸭子类型选中了谁、其中哪些被剔除、为什么」印出来 —— 让例外**看得见**。"""
    raw = likert_columns(d, exclude=set())
    kept = likert_columns(d)
    dropped = [c for c in raw if c not in kept]
    return dict(n_raw=len(raw), n_kept=len(kept), dropped=dropped,
                reasons={c: NOT_ITEMS.get(c, "⚠ 未说明理由 —— 不许无理由剔除") for c in dropped})
