"""零的构造(`#394`,由 `#393b` 要求)。

**三次同族之后建的基础设施**(`#385c` · `#391c` · `#437`):
每一次我都在写新代码时**从零重写零的构造**,而每一次重写都是一次重新犯错的机会。
三次里两次的表象是**多个门同时失败**,也就是**最劝我停下的那种表象**。

⚠ **用前必须跑 `controls()`**(`P5★` + `#374b`):
一个从未开火过的零构造器,它给出的每一个「未越阈」都是沉默,不是无罪。
"""
import numpy as np


def perm_in(v, mask, seed):
    """掩码内打乱,**保住缺失格局**。

    ⚠ 这是三次栽跟头的那一个。手写版通常是 `v[rng.permutation(len(v))]` ——
    它打乱**整个含 NaN 的数组**,于是掩码内出现 NaN,阈变成 NaN,
    **后面每一个门都恒为 False**,读起来像四个问题,其实是一个。
    """
    v = np.asarray(v, dtype=float)
    out = v.copy()
    j = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(v))
    if j.size:
        out[j] = v[j][np.random.default_rng(seed).permutation(j.size)]
    return out


def sham_control(x, target_corr, mask, seed):
    """合成一个与 `x` 相关度 = `target_corr`、**其余独立**的控制列。

    用途:`offset_control` 的地板 —— 「加任何与 x 相关的控制都会吃掉一点」,
    所以掉幅的零不是零(`#372a` · `#437` 都在用)。
    """
    m = np.asarray(mask, dtype=bool)
    x = np.asarray(x, dtype=float)
    n = int(m.sum())
    rg = np.random.default_rng(seed)
    zx = (x[m] - x[m].mean()) / max(x[m].std(), 1e-12)
    a = float(target_corr)
    out = np.full(x.shape, np.nan)
    out[m] = a * zx + np.sqrt(max(1.0 - a * a, 1e-9)) * rg.standard_normal(n)
    return out


def row_perm(cols, mask, seed):
    """**整行**打乱:把同一批人的整行一起搬走,**保住列与列之间的相关与维数**。

    用途:多列联合进模型时的零(`#390` 在用)。逐列独立打乱会**破坏列间相关**,
    于是零比真实情况更宽松 —— 那是一个**偏向通过**的零。
    """
    m = np.asarray(mask, dtype=bool)
    idx = np.flatnonzero(m)
    p = np.random.default_rng(seed).permutation(idx)
    out = []
    for c in cols:
        c = np.asarray(c, dtype=float)
        o = c.copy()
        o[idx] = c[p]
        out.append(o)
    return out


def controls():
    """三个构造器各自的对照。返回 dict,全部为 True 才可用。"""
    rg = np.random.default_rng(0)
    n = 4000
    m = np.zeros(n, bool); m[:3000] = True
    v = rg.standard_normal(n); v[np.arange(0, n, 7)] = np.nan     # 故意掺 NaN

    # perm_in:含 NaN 的输入必须仍给出**有限**的统计量(三次都栽在这一点)
    pv = perm_in(v, m, 1)
    finite_ok = bool(np.isfinite(pv[m & np.isfinite(v)]).all())
    pattern_ok = bool((np.isfinite(pv) == np.isfinite(v)).all())   # 缺失格局不变
    shuffled_ok = bool(not np.allclose(pv[m & np.isfinite(v)], v[m & np.isfinite(v)]))

    # sham_control:实现的相关度必须等于指定的
    x = rg.standard_normal(n)
    s = sham_control(x, 0.4, m, 2)
    got = float(np.corrcoef(x[m], s[m])[0, 1])
    sham_ok = abs(got - 0.4) < 0.05

    # row_perm:列间相关必须被保住
    a = rg.standard_normal(n); b = 0.6 * a + 0.8 * rg.standard_normal(n)
    ra = float(np.corrcoef(a[m], b[m])[0, 1])
    a2, b2 = row_perm([a, b], m, 3)
    rb = float(np.corrcoef(a2[m], b2[m])[0, 1])
    row_ok = abs(ra - rb) < 0.02
    # 负对照:逐列独立打乱**必须**破坏它(证明这个检验有分辨力)
    a3 = perm_in(a, m, 4); b3 = perm_in(b, m, 5)
    rc = float(np.corrcoef(a3[m], b3[m])[0, 1])
    row_neg_ok = abs(rc) < 0.1

    return dict(perm_in_finite=finite_ok, perm_in_pattern=pattern_ok,
                perm_in_shuffled=shuffled_ok, sham_corr=sham_ok,
                row_perm_keeps_corr=row_ok, row_perm_control_discriminates=row_neg_ok,
                _detail=dict(sham_got=round(got, 4), row_before=round(ra, 4),
                             row_after=round(rb, 4), row_indep=round(rc, 4)))
