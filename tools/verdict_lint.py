"""判词体检 —— ⛔ **回测不合格:不是闸,永远不要把它接进 `readme_gate.py`。**

⚠⚠ **`#623` 回测结果(2026-08-06,#755):recall 0/4 · precision 0/1。**
  催生它的四处判词缺陷(`#728` 比计数不比集合 · `#748` 分支必然触发 · `#750` 比错对象 ·
  `#754` 阈值留白),**它一处也没抓到**;而它在 730 个脚本上唯一开的那一枪(R080)
  经人工核对是**假的** —— 那处的外层 `else:` 是覆盖了的。
  ⇒ **按 `#733` 降级 `internal_consistency` 的同一先例,此文件保留为诊断脚本,不得成为规则。**

⚠ **它换来的是一条能力边界,而边界才是这轮的产出:**
  「判词比错了对象」**不可由代码溯源机检** —— 它是语义的,不是名字上的。
  第一版把变量传递展开地追源,连格式化字符串都算「同源」,命中 48/730 全是噪声;
  把展开压到 0 层,`#750` 那个真缺陷也就看不见了。**两端都不成立 ⇒ 这条路是封的。**
  唯一剩下的补救是**约定**:分支里每个被比较的量,人工写出它的**单位**。

⚠⚠ **而最有内容的一条:我为了抓某一族缺陷而建的仪器,自己犯了那一族的两个。**
  ① 它的 PROXY 是「赋值给叫 `v` 的东西」,而 `person_level.py` 里的 `v` 是一列 pandas 数据
     —— **仪器的单位 ≠ 主张的单位**,正是它要抓的那个错。
  ② 它把同一条 if/elif 链的每一节当成新链,同一处报三遍。
  两个都是我自己读输出读出来的,不是它报出来的。

原始动机(保留):判词是一句话,而一句话不会因为写错而报错。

⚠ 动机:同一族已经四次。`#728`(比的是计数不是集合)· `#748`(分支在两个单调序列上必然触发)·
`#750`(比的是「性 vs 性」而我想比「性 vs 非性」)· `#754`(W-1 与 W-2 之间留了一条谁也够不着的缝)。
**数据错了会崩,判词错了只会印出来。**

⚠ P6 代理账:
  PROPERTY   这条判词比的是我想比的东西,且它的分支覆盖了全部可能结果
  PROXY      ① 判词分支里每个比较的两侧,各自能追到哪个变量(**同源 = 可能在自比**)
             ② 同一个左值上的阈值链,中间有没有没被任何分支命名的区间
  IMPLICATION 只有一个方向可靠:**两侧同源 -> 确实值得人看一眼**(可靠)。
             反过来不成立:**不同源不证明比对了**。**只报可疑,从不认证判词正确。**
  SAFE SIDE  只打印,不判定,不阻断。它是一张必读清单,不是裁决。

用法:python3 tools/verdict_lint.py [路径…]   (默认扫全部 E0*/A*/R* 下的 .py)
"""
import ast, sys, pathlib, re

def _names(node):
    out=set()
    for n in ast.walk(node):
        if isinstance(n,ast.Name): out.add(n.id)
        elif isinstance(n,ast.Attribute): out.add(n.attr)
        elif isinstance(n,ast.Constant) and isinstance(n.value,str): out.add(n.value)
    return out

def _origin(name, assigns):
    """把一个变量名追到它最初来自的名字集合(一层展开就够,再深会引入噪声)。"""
    seen=set()
    def go(n,d=0):
        if d>3 or n in seen: return {n}
        seen.add(n)
        if n not in assigns: return {n}
        out=set()
        for src in assigns[n]:
            out|=go(src,d+1)
        return out or {n}
    return go(name)

def audit(path):
    try: tree=ast.parse(path.read_text())
    except SyntaxError as e: return [f"  ⚠ 解析失败:{e}"]
    assigns={}
    for n in ast.walk(tree):
        if isinstance(n,ast.Assign):
            for t in n.targets:
                if isinstance(t,ast.Name): assigns.setdefault(t.id,set()).update(_names(n.value)-{t.id})
    rows=[]
    # ⚠ 第一版对同一条 if/elif 链的每一节都当成新链,于是同一处被报三遍。
    #   先记下所有「作为别人 orelse 出现」的节点,它们不是链的开头。
    inner={id(c.orelse[0]) for c in ast.walk(tree)
           if isinstance(c,ast.If) and len(c.orelse)==1 and isinstance(c.orelse[0],ast.If)}
    for n in ast.walk(tree):
        if not isinstance(n,ast.If) or id(n) in inner: continue
        # 只看「给判词赋值」的分支链
        def assigns_verdict(body):
            # ⚠ 第一版只看变量名叫不叫 v —— 而 `person_level.py` 里的 `v` 是一列 pandas 数据。
            #   **PROPERTY 是「判词分支」,PROXY 是「赋值给叫 v 的东西」,两个单位不同。**
            #   这正是这具 lint 自己要抓的那一族,发生在它自己身上。收紧:**赋的值必须是字符串**。
            def is_str(x):
                if isinstance(x,ast.Constant): return isinstance(x.value,str)
                if isinstance(x,ast.JoinedStr): return True
                if isinstance(x,ast.BinOp): return is_str(x.left) or is_str(x.right)
                if isinstance(x,ast.IfExp): return is_str(x.body) or is_str(x.orelse)
                return False
            return any(isinstance(s,ast.Assign) and is_str(s.value)
                       and any(isinstance(t,ast.Name) and t.id in ("v","verdict","VERDICT") for t in s.targets)
                       for s in body)
        chain=[]; cur=n
        while True:
            chain.append(cur)
            if len(cur.orelse)==1 and isinstance(cur.orelse[0],ast.If): cur=cur.orelse[0]
            else: break
        if not any(assigns_verdict(c.body) for c in chain): continue
        rows.append(f"  判词链 @ line {n.lineno} · {len(chain)} 个分支" + ("" if chain[-1].orelse else "  ⚠ **没有 else —— 有结果落不进任何分支**"))
        thresholds={}
        for c in chain:
            for cmpn in [x for x in ast.walk(c.test) if isinstance(x,ast.Compare)]:
                L=_names(cmpn.left); R=set().union(*[_names(x) for x in cmpn.comparators]) if cmpn.comparators else set()
                oL=set().union(*[_origin(x,assigns) for x in L]) if L else set()
                oR=set().union(*[_origin(x,assigns) for x in R]) if R else set()
                src=ast.unparse(cmpn)[:74]
                # ⚠⚠ **同源检查已按 `#623` 回测后砍掉,不是调参。**
                #   第一版把变量传递展开地追到源头,结果格式化字符串也算「同源」:
                #   `dO > dC` 被标出 30 个共享名,而 dO 与 dC 是两个真正不同的量。
                #   命中 48/730,人工可核的几条全是假的 ⇒ **精度不合格。**
                #   而把展开压到 0 层,`#750` 那个真缺陷(`others` 装的是同一个网格的切片)
                #   **也就看不见了** —— 因为它是语义的,不是名字上的。
                #   ⇒ **记下能力边界:「判词比错了对象」不可由代码溯源机检。**
                #      它只能靠约定:分支里每个被比较的量,人工写出它的**单位**。
                # 阈值链:同一左值上的数值比较
                nums=[x.value for x in cmpn.comparators if isinstance(x,ast.Constant) and isinstance(x.value,(int,float))]
                if nums and L:
                    key=tuple(sorted(L))
                    thresholds.setdefault(key,[]).append((ast.unparse(cmpn.ops[0]).strip() if hasattr(cmpn.ops[0],'__doc__') else type(cmpn.ops[0]).__name__, nums[0], src))
        for key,ts in thresholds.items():
            if len(ts)<2: continue
            lo=[t for t in ts if t[0] in ("LtE","Lt")]; hi=[t for t in ts if t[0] in ("GtE","Gt")]
            if lo and hi:
                a=max(x[1] for x in lo); b=min(x[1] for x in hi)
                if b-a>1e-9:
                    rows.append(f"    ⚠ 阈值缝 {sorted(key)}:≤{a} 与 ≥{b} 之间的 ({a}, {b}) 没有任何分支命名")
                    rows.append(f"       -> 这是 `#754` 那一族:留白等于把裁量权还给跑完之后的我")
    return rows

def main(argv):
    root=pathlib.Path(__file__).resolve().parents[1]
    paths=[pathlib.Path(a) for a in argv[1:]] or sorted(root.glob("E0*/A*_*/R*_*/*.py"))
    tot=0
    for p in paths:
        r=audit(p)
        flags=[x for x in r if "⚠" in x]
        if flags:
            tot+=1; print(f"\n{p.relative_to(root) if root in p.parents else p}")
            for line in r: print(line)
    print(f"\n扫了 {len(paths)} 个脚本,{tot} 个有可疑判词。**只报可疑,从不认证判词正确。**")
    return 0
if __name__=="__main__": sys.exit(main(sys.argv))
