"""退化等式检查的 lint —— 「把差与零比」永远不会失败

⚠ 动机:这一族已经**四次**(`#769` 残差和 vs 0 · `#772` 常数 1.0 vs 1.0 ·
`#774` maxd vs 0.0,以及 `#770` 的第一版)。而 `#770` 自己写下过规则
——「等式检查要比两个值,不要比它们的差与零」—— **然后下一个脚本就破了它。**
⇒ 规则写进账本**不等于**规则进入下一个脚本;它需要一个**在写的时候拦住我**的东西。

⚠ 与 `#755` 那个被否掉的 lint 的差别,这是它能成立的原因:
**PROPERTY 是机械的** —— 「两个实参在语法上是同一个表达式」「其一是字面 0 而另一是一个差」
都在 AST 上直接可判。`#755` 要判「判词比错了对象」是语义的,造不出来。

⚠⚠ **边界(`#776` 当场量到的,写在这里而不是只写在账本里)**:**它只认语法形状,认不出值上的退化。**
`identity_control(observed=float(vmax(1.0)), expected=0.0)` —— `vmax(1.0)` 恰好返回 0,
而这具 lint 在跑之前扫过那个脚本、返回 **0 命中**;抓住它的是库运行时的 `_degenerate`。
⇒ **两层互补:这具 lint 在写的时候拦语法上的那部分,库在跑的时候拦全零那部分。**
⇒ 而正确的做法不是加深溯源(`#755` 已证明那条路是封的),是**改写法:直接把两个非零的值交给闸。**

⚠ P6 代理账:
  PROPERTY   这条 identity_control 有可能失败
  PROXY      两个实参的 AST 形状:同一表达式 / 一侧字面 0 且另一侧是减法或 abs(减法)
  IMPLICATION 只有一个方向可靠:**命中 -> 它确实不可能失败**(可靠)。
             反过来不成立:**没命中不证明这条检查有意义。**
  SAFE SIDE  只报「这条不可能失败」;从不认证「这条检查是好的」。
"""
import ast, sys, pathlib

def _norm(node):
    try: return ast.unparse(node)
    except Exception: return None

def _is_zero(node):
    return isinstance(node,ast.Constant) and isinstance(node.value,(int,float)) and node.value==0

def _is_difference(node):
    if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Sub): return True
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in ("abs","float"):
        return any(_is_difference(a) for a in node.args)
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Attribute) and node.func.attr in ("sum","max","min","mean","item"):
        # ⚠ 回测抓到的真 bug:`(A@b - y).sum()` 里那个差在 `func.value`,不在 `args`。
        return any(_is_difference(a) for a in node.args) or _is_difference(node.func.value)
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in ("max","min","sum","abs","float"):
        return any(_is_difference(a) for a in node.args)
    return False

# ── `#784` 追加的第二条规则:`X / X` 与 `X - X` ──────────────────────────────────
# ⚠ 为什么现在才加:`#784` 第一版把正控写成 `eq.append(d0/d0)` —— 一个恒等于 1 的量,
#   而**这具 lint 扫过那个脚本返回 0 命中**。原因是它只看 `identity_control` 的实参,
#   而这次的退化发生在一个**普通表达式**里,之后才喂给 `G.asserted`。
# ⚠ 它与既有规则同样机械(`unparse(左) == unparse(右)`),所以属于这具 lint 能造的那一半;
#   `#755` 证明过语义那一半造不出来。
# ⚠⚠ **回测把 `X - X` 那一半砍掉了,而这是回测该做的事(`#623`)。**
#   第一版同时收 Div 与 Sub,在 761 个历史脚本上命中 **2 处,人工逐个看过:两处都是有意的 g=0 臂**
#   (`R003` 种植阶梯的第 0 级 `-(C.slope.iloc[0]-C.slope.iloc[0])-1e-9`;
#    `R060` 的「同一分布与自身之差,必为 0」,连 `null_kind` 都写着)。⇒ **Sub 的精度 0/2。**
#   而**一个 g=0 臂本来就是「差为零」,那是本项目应该写的东西** —— 报它等于罚正确做法。
#   `X / X` 没有这个歧义:**没有哪种 g=0 臂是「比为一」。** 历史语料上 Div 命中 **0 处**
#   ⇒ 它对既有代码零误报,而对触发它的那个真缺陷(`d0/d0`)命中。
#   ⚠ 但 0/0 不是精度:**Div 这条规则在真实代码上从未返回过非零**,按 `P5` 的 ★ 条,
#   它目前只有合成夹具的正对照,**下一次它真的开火时必须人工复核,不许直接当缺陷。**
# ⚠ P6 代理账(与文件头同一张表,只补这一行):
#   PROXY      二元运算两侧的 AST 反解字符串相同,且运算是除或减
#   IMPLICATION 命中 -> 该表达式恒为 1 或 0(可靠)。**不命中仍不证明任何检查是好的。**
#   SAFE SIDE  只在两侧都不是字面常数时报(`2-2` 是有意为之的写法,不算缺陷)
def self_ops(tree):
    out=[]
    for n in ast.walk(tree):
        if not isinstance(n,ast.BinOp) or not isinstance(n.op,ast.Div): continue
        if isinstance(n.left,ast.Constant) or isinstance(n.right,ast.Constant): continue
        l,r=_norm(n.left),_norm(n.right)
        if l is None or l!=r: continue
        out.append((n.lineno, "`X / X` 恒为 1", f"{l} vs {r}"))
    return out


def audit(path):
    try: tree=ast.parse(path.read_text())
    except SyntaxError: return []
    out=[]
    for n in ast.walk(tree):
        if not (isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
                and n.func.attr=="identity_control"): continue
        kw={k.arg:k.value for k in n.keywords if k.arg}
        o,e=kw.get("observed"),kw.get("expected")
        if o is None or e is None:
            pos=[a for a in n.args if not isinstance(a,ast.Constant) or not isinstance(a.value,str)]
            if len(pos)>=2: o,e=pos[0],pos[1]
        if o is None or e is None: continue
        so,se=_norm(o),_norm(e)
        why=None
        if so is not None and so==se: why=f"两个实参是同一个表达式 `{so}`"
        elif (_is_zero(e) and _is_difference(o)) or (_is_zero(o) and _is_difference(e)):
            why=f"一侧是字面 0,另一侧是一个差 `{so if _is_zero(e) else se}`"
        elif isinstance(o,ast.Constant) and isinstance(e,ast.Constant) and o.value==e.value:
            why=f"两侧是同一个字面常数 `{o.value}`"
        if why: out.append((n.lineno,why,f"{so} vs {se}"))
    out.extend(self_ops(tree))
    return sorted(set(out))

if __name__=="__main__":
    root=pathlib.Path(__file__).resolve().parents[1]
    paths=[pathlib.Path(a) for a in sys.argv[1:]] or sorted(root.glob("E0*/A*_*/R*_*/*.py"))
    tot=0; files=0
    for p in paths:
        r=audit(p)
        if r:
            files+=1; tot+=len(r)
            print(f"\n{p.relative_to(root) if root in p.parents else p}")
            for ln,why,expr in r: print(f"  line {ln}: ⚠ **不可能失败** —— {why}\n      {expr[:90]}")
    print(f"\n扫了 {len(paths)} 个脚本,{files} 个文件、{tot} 处不可能失败的等式检查。"
          "**只报「不可能失败」,从不认证「这条检查是好的」。**")
