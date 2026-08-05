"""#432:把 `#398d` 从「靠记住」搬进「靠接口」。

`#398d` 的**措辞**是「从 `results/` 读数」,但 `#431c`/`#432a` 表明那是**目的的代理**:
458 轮里只有 **51** 轮(11%)真的从 `results/` 读回,而其余 89% 也没有违反它 ——
它们在内存里算完就下结论,**根本没经过终端**。

**真正的失败是:让一段**可能被截断**的终端输出成为证据。**
而截断只在输出**长**的时候发生,所以可机械化的那一半是:**让轮次的输出永远短**。

`show()` 是唯一被允许把一张表送上终端的入口:
- 最多打印 `n` 行,**并且在同一次调用里把整张表写进产物**;
- 打印的最后一行**永远**说明「完整的在哪里」,所以屏幕上的东西**自称是节选**;
- 被截掉的行数是**打印出来的**,不是静默的。

⚠ 它管不到的那一半:我在会话里临时跑的 shell。那一半只能靠 `#398d` 本身,
**而这正是「不可机械检查」的那一类还剩下的部分** —— 写在这里,不假装已经解决。
"""
import pathlib as _p

def show(df, path, n=12, sort=None, ascending=False, label=""):
    """打印至多 n 行,同时把整张表写到 `path`。返回写出的路径。"""
    path=_p.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    d = df.sort_values(sort, ascending=ascending) if sort else df
    head = d.head(n)
    print(head.to_string(index=False))
    hidden = len(d)-len(head)
    tag = f"{label} " if label else ""
    if hidden > 0:
        print(f"   … {tag}另有 **{hidden}** 行未显示 —— **完整的 {len(d)} 行在 `{path}`**")
    else:
        print(f"   ({tag}全部 {len(d)} 行已显示 · 产物 `{path}`)")
    return path

def controls():
    """四个自检:截断被报出 · 全显时不谎称截断 · 产物真的写了 · 行数与源一致。"""
    import pandas as pd, io, contextlib, tempfile
    T=pd.DataFrame(dict(a=range(30)))
    out=[]
    with tempfile.TemporaryDirectory() as td:
        for k,frame in (("trunc",T),("full",T.head(5))):
            buf=io.StringIO()
            with contextlib.redirect_stdout(buf):
                p=show(frame, f"{td}/{k}.csv", n=12)
            out.append((buf.getvalue(), p, frame))
        trunc_reported = "另有 **18** 行未显示" in out[0][0]
        full_no_lie    = "未显示" not in out[1][0]
        written        = out[0][1].exists() and out[1][1].exists()
        same_len       = len(pd.read_csv(out[0][1]))==30 and len(pd.read_csv(out[1][1]))==5
    return trunc_reported, full_no_lie, written, same_len


# ---------------------------------------------------------------- #436:比值的守门
def share(num, den_boot, num_boot=None, name=""):
    """把一个「占比」变成**可以拒绝发布**的东西。

    `#435c`:`pornhabit` 的间接效应稳(越零阈、自身区间不含 0),而它的**占比**
    `a·b / c` 的自助区间是 `[−3.59, +2.16]` —— 因为分母 `c = 0.0122` 贴着零。
    **一个效应可以是稳的,而它的占比同时是不可估的,因为占比的分母是另一个量。**

    `#436`:扫描抓不到这件事 —— 页面上写的是**比值**,分母的区间从来不在句子里。
    所以修法不是审计,是**让分母的区间成为计算占比的前提**(`#383a`:改接口)。

    返回 `(ok, value, lo, hi, reason)`;`ok=False` 时 **`value` 是 None** ——
    调用方拿不到那个数,而不是拿到之后被提醒不要用。
    """
    import numpy as _np
    d=_np.asarray(den_boot,dtype=float); d=d[_np.isfinite(d)]
    if d.size < 20:
        return (False, None, _np.nan, _np.nan, f"{name}分母自助样本太少({d.size})")
    dlo,dhi=_np.percentile(d,[2.5,97.5])
    if dlo<=0<=dhi:
        return (False, None, dlo, dhi,
                f"{name}**分母的区间含零** [{dlo:+.4g}, {dhi:+.4g}] -> 占比不可估")
    if num_boot is not None:
        n_=_np.asarray(num_boot,dtype=float)
        k=min(len(n_),len(d)); r=n_[:k]/d[:k]
        lo,hi=_np.percentile(r,[2.5,97.5])
    else:
        lo,hi=_np.nan,_np.nan
    return (True, float(num)/float(_np.median(d)), lo, hi, "")

def share_controls():
    """三个对照:近零分母被拒 · 远离零的分母通过 · 拒绝时不返回数值。"""
    import numpy as _np
    rng=_np.random.default_rng(3)
    bad =rng.normal(0.0,0.02,400)          # 分母贴着零
    good=rng.normal(0.50,0.02,400)         # 分母远离零
    ok1,v1,*_ = share(0.05, bad,  name="bad ")
    ok2,v2,*_ = share(0.05, good, name="good ")
    return (ok1 is False), (ok2 is True), (v1 is None)
