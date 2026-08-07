"""退化等式检查的 lint —— 「把差与零比」永远不会失败

⚠ 动机:这一族已经**四次**(`#769` 残差和 vs 0 · `#772` 常数 1.0 vs 1.0 ·
`#774` maxd vs 0.0,以及 `#770` 的第一版)。而 `#770` 自己写下过规则
——「等式检查要比两个值,不要比它们的差与零」—— **然后下一个脚本就破了它。**
⇒ 规则写进账本**不等于**规则进入下一个脚本;它需要一个**在写的时候拦住我**的东西。

⚠ 与 `#755` 那个被否掉的 lint 的差别,这是它能成立的原因:
**PROPERTY 是机械的** —— 「两个实参在语法上是同一个表达式」「其一是字面 0 而另一是一个差」
都在 AST 上直接可判。`#755` 要判「判词比错了对象」是语义的,造不出来。

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
    return out

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
