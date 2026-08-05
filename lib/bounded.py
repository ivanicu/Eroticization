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


# ---------------------------------------------------------------- #456:中介的事前上界
def mediation_headroom(x, y1, y2, mask, min_move=0.05, name=""):
    """**在跑中介之前**算出「控制 `y2` 最多能把 `x -> y1` 的系数移动多少」。

    `#456b`:`R500` 的中介检验「通过」了,而它几乎是**被逼出来**的 ——
    两个结局彼此相关只有 −0.0022,所以控制其中任何一个,本来就拿不走另一个的任何东西。
    **一个通过得太轻松的检验,和一个不可能失败的检验,只差一点点。**
    而分辨两者只需要**三个相关系数**,不需要跑那一轮。

    在标准化下,OLS 的恒等式是 `c = c' + a·b`,而
        a = r(x, y2)
        b = (r(y1,y2) − r(y1,x)·r(x,y2)) / (1 − r(x,y2)²)
    所以 **间接项 a·b 完全由三个相关决定**,可以事前算出。
    `r(y1,y2) ≈ 0` ⇒ `a·b ≈ −r(x,y2)²·r(y1,x)/(1−r(x,y2)²)` ⇒ 二阶小量 ⇒ **无空间可拿**。

    返回 `(ok, indirect, share, why)`;`ok=False` 时 **`indirect` 是 None** ——
    调用方**拿不到那个数**,而不是拿到之后被提醒它没有意义(与 `share()` 同一形状)。
    `min_move` = 相对于总效应 `c` 的最小可辨移动比例。
    """
    import numpy as _np
    m=_np.asarray(mask,dtype=bool)
    def rr(u,v):
        g=m&_np.isfinite(u)&_np.isfinite(v)
        return float(_np.corrcoef(_np.asarray(u)[g],_np.asarray(v)[g])[0,1])
    r_x_y2=rr(x,y2); r_y1_y2=rr(y1,y2); r_y1_x=rr(y1,x)
    den=1.0-r_x_y2**2
    if abs(den)<1e-9:
        return (False, None, _np.nan, f"{name}x 与 y2 共线(r={r_x_y2:+.4f})")
    a=r_x_y2; b=(r_y1_y2-r_y1_x*r_x_y2)/den
    ind=a*b; c=r_y1_x
    share=ind/c if abs(c)>1e-9 else _np.nan
    if abs(c)<1e-9:
        return (False, None, _np.nan, f"{name}总效应贴零(r={c:+.4g})-> 占比无意义")
    if abs(share)<min_move:
        return (False, None, share,
                f"{name}**事前上界只有总效应的 {abs(share):.2%}**(< {min_move:.0%}) -> "
                f"这个中介检验**跑不出信息**:r(y1,y2)={r_y1_y2:+.4f}")
    return (True, float(ind), float(share), "")

def headroom_controls():
    """三个自检:正交结局被拒 · 真中介放行 · 拒绝时不返回数值。"""
    import numpy as _np
    rng=_np.random.default_rng(9); n=4000; m=_np.ones(n,dtype=bool)
    x=rng.normal(size=n)
    y2=0.3*x+rng.normal(size=n)                     # 与 x 相关
    y1_orth=rng.normal(size=n)+0.3*x                # 与 y2 几乎无关
    y1_med =0.5*y2+rng.normal(size=n)               # 真的经 y2
    ok1,v1,_,_=mediation_headroom(x,y1_orth,y2,m,name="orth ")
    ok2,v2,_,_=mediation_headroom(x,y1_med ,y2,m,name="med ")
    return (ok1 is False), (ok2 is True), (v1 is None)
