#!/usr/bin/env python3
"""#840① —— 把「绝不手抄上一轮的数字」从决心变成工具。

**为什么需要工具而不是决心**:`#836`① 教过同一条道理 ——
**教训不落到工具上就是没落地**。而 `#840` 当场证了它:
`R279` 第一版把 `#838` 的四位小数手抄进脚本(`D838 = {1990: (+0.1137, +0.6266), …}`),
于是**恒等式控制在 2.3e−5 上失败** —— 失败的不是那个恒等式,是**转抄的精度**。

**检测什么(P6 代理账):**
  PROPERTY    一个数是从别的轮次的产物里**手抄**过来的
  PROXY       脚本里出现一个**小数位 ≥4** 的字面量,而某个 `results/*.json` 里
              存在一个值 `v`,使得该字面量恰好是 `v` 的四舍五入(或与 `v` 精确相等)
  IMPLICATION 只有一个方向可靠:**命中 ⇒ 这个数确实与某份产物对得上**(可靠)。
              反过来**不成立**:没命中不证明没手抄 —— 产物可能没存那个量,
              或作者抄的是打印输出而不是 JSON。**只报命中,从不报「本文件干净」。**
  SAFE SIDE   报「疑似转抄 + 它对应的产物路径」,由作者改成从产物读;**从不自动改写代码。**

⚠ **小数位 ≥4 这条门槛是有理由的,不是随手定的**:预注册阈值、容差、q 值几乎总是
  **人选的圆数**(0.05 · 0.30 · 1.5 · 120),而**从产物抄回来的量几乎总是四位以上**。
  ⇒ 这条门槛把「作者的选择」与「机器算出来的量」分开,而这正是要区分的两类。

⚠ **默认棘轮而非零容忍**(`L81` + `readme_gate` 的先例):
  历史轮次里已有的命中不该让今天的提交失败;基线存在 `tools/no_transcribed_baseline.json`。

⚠⚠⚠ **第一版用「数值比对」做代理,而它被自己的零杀掉了 —— 记在这里,因为这是本工具的主要发现:**
  第一版:脚本里小数位 ≥4 的字面量,若与某个 `results/*.json` 里的值相等或是它的四舍五入 ⇒ 判为手抄。
  正控过了(种一个抄来的值,命中),**于是我差点就发了它。**
  ⇒ 但 `realstat` 有一条专门的:**「正控只问『这仪器看得见吗』,从不问『它看见的是不是我要主张的东西』。」**
  ⇒ 于是量了它的**误报率**:拿 4000 个**随机**字面量比对真实产物(11,376 个有限值,90.6% 落在 |v|<1):
     **4 位小数 → 误报 32.0% · 5 位 → 4.2% · 6 位 → 0.4% · 8 位 → 0.0%**
  **而第一版报了 1057 处 —— 其中一大半是巧合。** 证据不用统计也看得见:
  `E01/R001` 的字面量被判为抄自 `E03/R195` 与 `R265`,**而那些轮次在它之后好几百轮才存在。**
  ⇒ **更要命的是把门槛提到 6 位也救不了它**:6 位命中的是**全精度复制**,
     而全精度复制在数值上无害;`#840` 真正出事的那个是 `0.6266 ← 0.626612`,**一个四位四舍五入** ——
     **恰好落在误报率 32% 的那一档。**
  ⇒ **结论:数值比对这个代理,结构上抓不到它被造出来要抓的那种错。不是调参数能救的。**

**⇒ 第二版换代理,而新代理不看数值,看出处(P6 代理账):**
  PROPERTY    一个数是从别的轮次的产物里**手抄**过来的
  PROXY       脚本里出现对**另一轮编号**的引用(`#NNN` / `DNNN` / `RNNN`),
              **同一处附近有小数位 ≥4 的字面量**,而该脚本**从不 `json.load` 任何路径含该轮编号的产物**
  IMPLICATION 只有一个方向可靠:**命中 ⇒ 这脚本确实在提别的轮次、且确实没去读它的产物**(可靠)。
              反过来**不成立**:作者可以不提编号就手抄。**只报命中,从不报「干净」。**
  SAFE SIDE   报「提到了 `#NNN` 却没读它的产物 + 附近有高精度字面量」,由作者改成从产物读。
  ⚠ 这个代理**不依赖数值巧合**,所以它的误报来源是「引用了但确实不需要读产物」——
    **那是一个人可以一眼判掉的误报,不是一个需要统计才能发现的误报。**
"""

import ast, io, json, pathlib, re, sys, tokenize, argparse

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASE = ROOT/"tools"/"no_transcribed_baseline.json"
LIT = re.compile(r"(?<![\w.])[-+]?(\d+)\.(\d{4,})(?![\w.])")
REF = re.compile(r"[#DR](\d{3,4})\b")
LOADS = re.compile(r"json\.load\s*\(\s*open\s*\(([^)]*)\)")

def scan(paths=None):
    """第二版:不看数值,看出处。见文件头 —— 为什么数值比对被它自己的零杀掉。"""
    hits = []
    files = paths if paths else sorted(ROOT.rglob("E0*/**/*.py"))
    for f in files:
        f = pathlib.Path(f)
        if "_archive" in f.parts: continue
        try: txt = f.read_text(encoding="utf-8")
        except Exception: continue
        # ⚠⚠ **只看可执行代码里的字面量,不看模块 docstring 里的。** 理由是 `#840` 那个 bug 的形状:
        #    出事的是**一个进了控制期望值的数**,不是一段引用旧结果的散文。
        #    第一版不分,于是 890 处里绝大多数是 docstring 在引用上一轮 —— **那是引用,不是消费。**
        #    ⇒ 用 `ast` 把模块 docstring 的字符区间切掉;解析失败就整个文件跳过(**宁可漏,不可冤**)。
        try:
            _tree = ast.parse(txt)
        except SyntaxError:
            continue
        # ⚠⚠ **第一版这里只跳 `body[0]` 的 docstring,而实测它经常不在 body[0]:**
        #    `R003/amplify_with_forced_choice.py` 前三行是 `import`/`os.chdir`,docstring 在第 5 行 ——
        #    于是 `body[0]` 是 `Import`,跳过从没生效,而那段散文里的 `#150b` 引用被报了三次。
        #    **一个只覆盖「标准写法」的跳过规则,在非标准写法上静默失效** —— 而失效方向是**多报**。
        #    ⇒ 改成跳过**模块里所有**「字符串常量表达式语句」(任何位置的 docstring),
        #      并用 `tokenize` 额外跳掉**注释**——两者都是**引用**,不是**消费**。
        _skip = []
        for _n in ast.walk(_tree):
            if (isinstance(_n, ast.Expr) and isinstance(getattr(_n, "value", None), ast.Constant)
                    and isinstance(_n.value.value, str)):
                _ls = txt.split("\n")
                _skip.append((sum(len(x)+1 for x in _ls[:_n.lineno-1]),
                              sum(len(x)+1 for x in _ls[:_n.end_lineno])))
        try:
            _ls2 = txt.split("\n")
            for _tok in tokenize.generate_tokens(io.StringIO(txt).readline):
                if _tok.type == tokenize.COMMENT:
                    _b = sum(len(x)+1 for x in _ls2[:_tok.start[0]-1])
                    _skip.append((_b+_tok.start[1], _b+_tok.end[1]))
        except Exception:
            continue                      # 分词失败 ⇒ 整个文件跳过(宁可漏,不可冤)
        loaded = " ".join(LOADS.findall(txt))
        m0 = re.match(r"R(\d+)_", f.parent.name)
        own = m0.group(1) if m0 else None
        for m in LIT.finditer(txt):
            if any(a <= m.start() < b for a, b in _skip): continue
            # ⚠ **补零的圆数是作者的选择,不是机器的输出**:`-0.2000` 其实是 `-0.2`,
            #   `0.0500` 其实是 `0.05`。本仪器要抓的是**从产物抄回来的量**,而那种量
            #   几乎不会在第 3 位以后全是零。⇒ 去掉尾零后有效小数位 <3 的一律跳过。
            if len(m.group(2).rstrip("0")) < 3: continue
            line = txt.count("\n", 0, m.start())+1
            lo = max(0, m.start()-260); hi = min(len(txt), m.end()+260)
            refs = {r for r in REF.findall(txt[lo:hi]) if r != own}
            unread = sorted(r for r in refs if r not in loaded)
            if unread:
                hits.append(dict(file=str(f.relative_to(ROOT)), line=line, literal=m.group(0),
                                 exact=False, artifact_value=None,
                                 sources=[f"引用 #{unread[0]} 却没读它的产物"]))
    return hits

def controls():
    """★ P5:先证它会开火,再证它不会对合规写法开火。"""
    import tempfile
    d = pathlib.Path(tempfile.mkdtemp(dir=str(ROOT), prefix="_ctl_"))
    (d/"R900_ctl").mkdir()
    (d/"R900_ctl"/"a.py").write_text(
        "D838 = {1990: (+0.1137, +0.6266)}   # 抄自 #838\n", encoding="utf-8")
    pc = len(scan([d/"R900_ctl"/"a.py"])) >= 1
    (d/"R900_ctl"/"b.py").write_text(
        'g = json.load(open("E03/R277_x/results/who_838.json"))   # #838\n'
        "X = g['a'] + 0.1137\n", encoding="utf-8")
    nc1 = len(scan([d/"R900_ctl"/"b.py"])) == 0
    (d/"R900_ctl"/"c.py").write_text("Q = 0.05\nTOL = 0.30\nB = 6000\n", encoding="utf-8")
    nc2 = len(scan([d/"R900_ctl"/"c.py"])) == 0
    for x in (d/"R900_ctl").iterdir(): x.unlink()
    (d/"R900_ctl").rmdir(); d.rmdir()
    return pc, (nc1 and nc2), f"负控α(引用且已读产物){nc1} · 负控β(只有人选圆数){nc2}"

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--precommit", action="store_true")
    ap.add_argument("--rebaseline", action="store_true")
    a = ap.parse_args()
    pc, nc, note = controls()
    print(f"仪器控制:正控(抄一个产物值必须命中)**{pc}** · 负控(只有人选圆数不许命中)**{nc}** —— {note}")
    if not (pc and nc):
        print("⛔ **控制没过 ⇒ 本次扫描不可采**(★ P5:仪器没证明会开火,它的零是沉默)"); sys.exit(2)
    hits = scan()
    print(f"疑似从别轮产物**手抄**的字面量:**{len(hits)}** 处")
    for h in hits[:12]:
        print(f"  {h['file']}:{h['line']}  `{h['literal']}`  ⇐ {h['sources'][0]}")
    if len(hits) > 12: print(f"  …另 {len(hits)-12} 处")
    print("⚠ **只报命中,从不报「本文件干净」** —— 没命中不证明没手抄(产物可能没存那个量)。")
    if a.rebaseline:
        BASE.write_text(json.dumps({"count": len(hits)}, indent=1)); print(f"基线写入 {len(hits)}"); sys.exit(0)
    if a.precommit:
        old = json.loads(BASE.read_text())["count"] if BASE.exists() else None
        if old is None:
            BASE.write_text(json.dumps({"count": len(hits)}, indent=1))
            print(f"首次建立基线 = {len(hits)}"); sys.exit(0)
        if len(hits) > old:
            print(f"🔒 PRE-COMMIT BLOCK:手抄字面量 {old} -> {len(hits)}(+{len(hits)-old})"); sys.exit(1)
        if len(hits) < old:
            BASE.write_text(json.dumps({"count": len(hits)}, indent=1))
            print(f"棘轮收紧:{old} -> {len(hits)}")
    sys.exit(0)
