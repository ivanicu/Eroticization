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
