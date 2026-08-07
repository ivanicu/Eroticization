"""还有几条守则只靠拷贝活着 —— 把项目里的「除了 X」数出来,并说出它们住在哪里。

⚠⚠ **这个工具存在的理由是一次实测的灭绝:**
`#680` 查码本查出 `colcom` 极性反了,写进账本,并在代码里写下 `c != "colcom"`;
`R123` · `R127` · `R151` 等约 180 轮**都带着那个例外,一次没错** ——
**而 `#866` 从头重写题目构造时,它静静地消失了。**
**那条守则活在一个被逐轮拷贝的 lambda 里。靠拷贝续命的规则,离灭绝永远只差一次重写。**
⇒ `#868` 把它搬进 `lib/gss_polarity.py`。**本工具问的是:还有多少条这样的。**

⚠⚠⚠ **`realstat`:一个 grep 是一具测量仪器,而它必须有正控。**
本模块因此**先在答案已知的地方跑一遍**(`colcom` 例外,已知存在于若干轮),
再在**不该有的地方**跑一遍(`colhomo` 从来不是例外,必须返回 0)。
**两条都过,零才是测量;否则那是沉默,不是无罪。**

⚠ **P6 代理账 —— 仪器的单位与主张的单位不相等,必须写下来:**
  PROPERTY    一条**规则**:把某一项/某一码/某一年从一条本来统一的处理里单拎出来,
              且**住在轮次脚本里而不是共享模块里**(所以下一次重写就断)
  PROXY       一**行** Python 命中下面的形状之一
  IMPLICATION 只有一个方向可靠:**命中 ⇒ 那一行确实在做特例**(可靠)。
              反过来不成立:**没命中不证明没有例外** —— 例外可以跨行写、可以藏在数据字典里、
              可以用变量名间接表达。**所以本工具只报下界。**
  WITNESS     `lib/gss_polarity.py` 里的 `ITEM_EXCEPTIONS` 是一条**规则**,
              但它**不以行内 `!=` 的形状出现** ⇒ 形状表看不见它 ⇒ 证明 PROXY ⊊ PROPERTY
  SAFE SIDE   只报「至少有这些」,**从不报「例外一共就这么多」**
"""
import pathlib
import re
import sys
import tokenize

ROOT = pathlib.Path(__file__).resolve().parents[1]
SHARED = ("lib/", "tools/")          # 住在这里 = 被 import,不被拷贝
SCAN = ("E01_sexual_as_a_value_not_a_category", "E02_condemnation_is_not_rarity",
        "E03_what_an_instrument_would_have_to_be", "lib", "tools")

# 形状表。每一条都**窄**,而且**每一条都必须自带正控**:
# 一个命中一切的模式不是保守,是没被检验过(`realstat`:一个 grep 是仪器)。
# ⚠ 第一版把 `category_pick`(裸 `== <整数>`)也算了进来,它在 **169 轮**里命中 `== 1` ——
#   那不是「例外规则」,那是每一个整数比较。**它没有正控,而它一个人就淹掉了整张表。**
#   ⇒ 现在:**没有正控的形状不计入,只列出并标 UNCONTROLLED。**
# `positive` = 一个**答案已知存在**的 token;`negative` = 一个**必须为 0** 的 token。
SHAPES = {
    "item_exception":  (re.compile(r'!=\s*[\'"]([A-Za-z_][A-Za-z0-9_]{2,})[\'"]'),
                        "colcom", "colhomo"),
    "value_whitelist": (re.compile(r'isin\(\s*\[([^\]]{1,60})\]\s*\)'),
                        "4, 5", "42, 43"),
    "year_guard":      (re.compile(r'year\s*[<>]=?\s*(\d{4})'), None, None),
    "exclusion_list":  (re.compile(r'not\s+in\s*[\(\[]([^\)\]]{1,60})[\)\]]'), None, None),
    # ⛔ `category_pick` 已移除 —— 无正控且命中一切,见上。
}


def _files():
    for d in SCAN:
        p = ROOT / d
        if not p.exists():
            continue
        yield from sorted(p.rglob("*.py"))


def scan():
    """返回 [(rel_path, shape, token, lineno, is_shared)]。"""
    hits = []
    for f in _files():
        rel = str(f.relative_to(ROOT))
        if "/_archive/" in f"/{rel}" or "/__pycache__/" in f"/{rel}":
            continue
        shared = rel.startswith(SHARED)
        try:
            lines = f.read_text(errors="replace").split("\n")
        except OSError:
            continue
        # ⚠⚠ **字符串字面量(含 docstring)里的不算,而这一条是实测出来的,不是想出来的。**
        # 本模块第一版只剥注释,不剥字符串 ⇒ `lib/bks_items.py` 的 docstring 里**引用**了那行
        # 原始代码 `lik=[c for c in lik if c!='biomale']`,于是审计把 `biomale` 判成
        # 「已经有家了」,at_risk 从 19 掉到 18 —— **而那个家是一句引文,不是一个实现。**
        # ⚠ **它错在讨好的方向**:光是在 docstring 里提一句,就能让待办数字变小。
        # ⇒ 用 `tokenize` 跳过 STRING 与 COMMENT,只看**会被执行的**那部分。
        try:
            import io as _io
            src = "\n".join(lines)
            code_lines = {}
            for tok in tokenize.generate_tokens(_io.StringIO(src).readline):
                if tok.type in (tokenize.STRING, tokenize.COMMENT):
                    continue
                code_lines.setdefault(tok.start[0], []).append(tok.line)
        except (tokenize.TokenError, IndentationError, SyntaxError):
            code_lines = {i: [ln] for i, ln in enumerate(lines, 1)}   # 解析不了就退回整行,如实偏保守
        for i, ln in enumerate(lines, 1):
            if i not in code_lines:          # 该行只有字符串/注释 ⇒ 不会被执行
                continue
            s = ln.split("#", 1)[0]
            for shape, (rx, _pos, _neg) in SHAPES.items():
                for m in rx.finditer(s):
                    hits.append((rel, shape, m.group(1).strip(), i, shared))
    return hits


def controls(hits):
    """⚠ grep 是仪器 ⇒ **每一个形状**各自的正控在答案已知处跑,负控在不该有的地方跑。

    返回每形状的 `(正控命中文件数, 负控命中文件数, 是否受控)`。
    **不受控的形状不计入总数** —— 它的计数只是列出来,标 `UNCONTROLLED`。
    """
    out = {}
    for shape, (_rx, pos, neg) in SHAPES.items():
        if pos is None:
            out[shape] = dict(controlled=False, positive_n=None, negative_n=None)
            continue
        p = {h[0] for h in hits if h[1] == shape and h[2] == pos}
        n = {h[0] for h in hits if h[1] == shape and h[2] == neg}
        out[shape] = dict(controlled=bool(p) and not n, positive_n=len(p), negative_n=len(n),
                          positive_token=pos, negative_token=neg)
    return out


def _homed_tokens():
    """一个 token **有没有家**:它是否在 `lib/` 的任何地方出现过 —— **字符串字面量也算**。

    ⚠⚠ **两侧必须用两条不同的规则,而这是实测逼出来的:**
    · 「它是不是轮次脚本里的行内特例」⇒ **只看会被执行的代码**(字符串/注释剥掉),
      否则 docstring 里引一句原始代码就能让待办数字变小 —— **讨好方向的假阴性**。
    · 「它有没有家」⇒ **要连字符串字面量一起看**,因为**一个登记表的家就长成字典的键**:
      `lib/gss_polarity.py` 的 `ITEM_EXCEPTIONS = {"colcom": +1}`、
      `lib/bks_items.py` 的 `NOT_ITEMS = {"biomale": ...}` —— 它们**是**家,却不是 `!=` 的形状。
    **只用一条规则时,两个方向各错一次:先是 18(把引文当家),后是 20(看不见字典键)。**
    ⚠ 这正是本模块 `P6` 代理账里写下的那个 WITNESS —— **它是跑之前就预测到的,不是事后解释。**
    """
    toks = set()
    for d in SHARED:
        p = ROOT / d
        if not p.exists():
            continue
        for f in sorted(p.rglob("*.py")):
            if "__pycache__" in str(f):
                continue
            txt = f.read_text(errors="replace")
            for rx, _pos, _neg in SHAPES.values():
                for m in rx.finditer(txt):
                    toks.add(m.group(1).strip())
            for m in re.finditer(r"""^\s*['"]([A-Za-z_][A-Za-z0-9_]{2,})['"]\s*:""", txt, re.M):
                toks.add(m.group(1))          # 登记表的键 = 家
    return toks


def at_risk(hits):
    """**靠拷贝活着** = 同一个 token 在 ≥2 个轮次脚本的**可执行代码**里,且**在 `lib/` 里没有家**。

    这正是 `colcom` 在 `#866` 之前的形状:180 轮拷贝、零个 import。
    """
    from collections import defaultdict
    homed = _homed_tokens()
    per = defaultdict(set)
    for rel, shape, tok, _ln, shared in hits:
        if not shared:
            per[(shape, tok)].add(rel)
    out = []
    for (shape, tok), files in per.items():
        if len(files) >= 2 and tok not in homed:
            out.append(dict(shape=shape, token=tok, n_rounds=len(files), files=sorted(files)))
    return sorted(out, key=lambda d: -d["n_rounds"])


if __name__ == "__main__":
    h = scan()
    c = controls(h)
    ctrl = [s for s, v in c.items() if v["controlled"]]
    print(f"扫了 {len(list(_files()))} 个 .py · 命中 {len(h)} 行")
    for s, v in c.items():
        if v["controlled"]:
            print(f"  {s:16s} **受控**(正控 `{v['positive_token']}` 命中 {v['positive_n']} 文件 · "
                  f"负控 `{v['negative_token']}` = {v['negative_n']})")
        else:
            n = len({x[0] for x in h if x[1] == s})
            print(f"  {s:16s} **UNCONTROLLED —— 不计入**(出现在 {n} 个文件,但没有已知答案可校)")
    r = [d for d in at_risk(h) if d["shape"] in ctrl]
    print(f"\n**靠拷贝活着的 token(只数受控形状):{len(r)}**")
    for d in r[:15]:
        print(f"  {d['shape']:16s} {d['token']:24s} {d['n_rounds']:3d} 个轮次脚本 · 0 个共享模块")
    ok = bool(ctrl)
    sys.exit(0 if ok else 2)
