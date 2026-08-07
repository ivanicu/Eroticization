"""#850 · E03·A81·R289 —— 把 `#807`/`#828` 的判词也重算:方法早在手,只是没跑

**还 `#843`①。** `#843` 把 `#798` 的判词在**当时页面 × 今天页面 × 旧切法 × 新切法**四格上
重算过(全为 0),并**明确登记**剩下两条**「不是做不到,是我只做了最可能翻的那一条」**。
⇒ **本轮把那句登记兑现 —— 而这正是 `#843` 自己立的规矩:
一条「改不了」要写成可检验的形式,否则它就变成下一堵没查过的墙。**

G1 估计量:`#807` 与 `#828` 各自的判词,在 **{当时页面, 今天页面} × {旧切法, 新切法}** 四格上的值。
   `#807`:两条更正(**顶对够不着** · **一般化缩回 `homosex`**)有没有落到页面上 ——
          判据是「提到该主张、却没有对应标记」的单位数。
   `#828`:同一形状,目标是 `sexeduc/1990s` 与 `teensex/2010s` 的**未校正标记**。

三个世界:
   A **四格全一致** ⇒ 两条判词都不受切法影响,`#833` 的缺陷没碰到它们。
   B **当时旧 ≠ 当时新** ⇒ **那一条当时就算错了,判词要撤。**
   C **当时一致而今天不同** ⇒ **判词没错,页面退化了** —— 要修的是页面。

预注册判词(条件式):
  if 正控开火(**在当时页面植入一处「提到主张但无标记」的单位,两种切法都必须数到它**)
     and 负控开火(**给同一处补上标记后,两种切法都必须不再数它**):
      任一条 当时旧 != 当时新 -> B
      当时一致而今天 > 0      -> C
      四格全 0/一致           -> A
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**谓词与页面必须同代** —— `#807`/`#828` 的标记正则里含它们自己的编号,
  而「当时的页面」正是刚写完它们那一版。⇒ 控制:**谓词也从那次提交的脚本里读**,并印出两代是否相同。
⚠ 本轮对象是页面自己 ⇒ 换不了仪器。⚠ 标注 **Production**,不产生新的关于人的判断。
"""
import ast, json, pathlib, re, subprocess, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.page_units import units
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
E3 = "E03_what_an_instrument_would_have_to_be"
JOBS = {
 807: dict(script=f"{E3}/A57_两条更正要走到页面上/R246_不是没跑是跑不了而页面还没这么说/two_corrections_must_land.py",
           pairs=[("TOP", "MARK_TOP"), ("GEN", "MARK_GEN")]),
 828: dict(script=f"{E3}/A72_两条更正落地/R267_未校正的标记与从没重数的欠账/land_two_corrections.py",
           pairs=[("TARGET", "MARKED")]),
}
def git(*a): return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout
def pats(src):
    out = {}
    for n in ast.parse(src).body:
        if (isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)
                and isinstance(n.value, ast.Call) and getattr(n.value.func, "attr", None) == "compile"):
            try: out[n.targets[0].id] = eval(compile(ast.Expression(n.value), "<p>", "eval"), {"re": re})
            except Exception: pass
    return out

print("=== ⓪ 从对象读:每条是哪一次提交引入的,以及谓词两代是否相同 ===")
ST = {}
for n, J in JOBS.items():
    c = git("log", "--format=%H", "-S", f"Entry {n} ·", "--", "RETRACTIONS.md").strip().split("\n")[-1]
    P_then, P_now = pats(git("show", f"{c}:{J['script']}")), pats((ROOT/J["script"]).read_text(encoding="utf-8"))
    same = {k: (P_now[k].pattern == P_then.get(k, re.compile("")).pattern)
            for pr in J["pairs"] for k in pr}
    ST[n] = dict(commit=c, then=P_then, now=P_now, same=same)
    print(f"  `#{n}` ← {c[:12]} · {git('log','-1','--format=%ad','--date=short',c).strip()} · "
          f"谓词两代相同 {same}")

def naked(txt, P, pairs):
    """提到某主张(claim 正则)却没有对应标记(mark 正则)的单位数。"""
    tot, det = 0, []
    for kind, owner, body in units(txt, fixed=True) + [(k, o, b) for k, o, b in units(txt, fixed=False)][:0]:
        pass
    return tot, det
def count(txt, P, pairs, fixed):
    tot, det = 0, []
    for kind, owner, body in units(txt, fixed=fixed):
        for cl, mk in pairs:
            if cl in P and P[cl].search(body) and not (mk in P and P[mk].search(body)):
                tot += 1; det.append((kind, owner, cl)); break
    return tot, det

print(f"\n=== ① 四格重算(当时/今天 × 旧切法/新切法)===")
G4 = {}
for n, J in JOBS.items():
    pages = {"当时": {f: git("show", f"{ST[n]['commit']}:{f}") for f in ("README_zh.md", "README.md")},
             "今天": {f: (ROOT/f).read_text(encoding="utf-8") for f in ("README_zh.md", "README.md")}}
    P = {"当时": ST[n]["then"], "今天": ST[n]["now"]}
    row = []
    for era in ("当时", "今天"):
        for cut, fx in (("旧", False), ("新", True)):
            tot = sum(count(pages[era][f], P[era], J["pairs"], fx)[0] for f in pages[era])
            G4[(n, era, cut)] = tot; row.append(f"{era}×{cut} **{tot}**")
    print(f"  `#{n}`:" + " · ".join(row))
    ST[n]["pages_then_len"] = {f: len(pages["当时"][f]) for f in pages["当时"]}

print(f"\n=== ② 控制(在**当时的页面**上做,与判词同代)===")
n0 = 828; J0 = JOBS[n0]; P0 = ST[n0]["then"]
base_txt = git("show", f"{ST[n0]['commit']}:README_zh.md")
cl0 = J0["pairs"][0][0]
# ⚠⚠ **第一版取的是命中单位的前 60 个字符,而那 60 个字符里未必包含命中的那一段** ——
#    于是植入进去的文本根本不触发 `TARGET`,正控 +0。**夹具没有编码它要检验的性质,
#    而这是本会话第三次同一形状**(`#846` 的 `⟨无作用域⟩`、`#848` 的正控夹具)。
#    ⇒ 改成**把真正命中的那一段原样取出来**(`.group(0)`),而不是取前缀。
_m = next((P0[cl0].search(b) for k, o, b in units(base_txt, fixed=True) if P0[cl0].search(b)), None)
_seed = _m.group(0) if _m else "sexeduc/1990s"
print(f"  ⚠ 正控夹具用的是**真正命中的那一段**(`.group(0)`),不是前缀:`{_seed[:40]}`")
INJ = "\n植入对照:" + _seed + " `[#700「植入」]`\n"
pc = {c: count(base_txt+INJ, P0, J0["pairs"], fx)[0] - count(base_txt, P0, J0["pairs"], fx)[0]
      for c, fx in (("旧", False), ("新", True))}
INJ2 = INJ.replace("`[#700「植入」]`", "而它在 `#826` 校正后不存活 `[#701「植入带标记」]`")
assert P0[cl0].search(INJ), "夹具必须真的触发 claim 正则,否则正控测的不是它该测的东西"
nc = {c: count(base_txt+INJ2, P0, J0["pairs"], fx)[0] - count(base_txt, P0, J0["pairs"], fx)[0]
      for c, fx in (("旧", False), ("新", True))}
print(f"  正控:在当时页面植入**一处提到主张但无标记**的单位 ⇒ 旧 +{pc['旧']} · 新 +{pc['新']}(都该 +1)")
print(f"  负控:同一处**补上标记** ⇒ 旧 +{nc['旧']} · 新 +{nc['新']}(都该 **+0**)—— "
      f"⚠ **「这个零该不该是零?」该**:带了标记按定义就不算裸着")

flip = [n for n in JOBS if G4[(n, "当时", "旧")] != G4[(n, "当时", "新")]]
now_bad = [n for n in JOBS if G4[(n, "今天", "新")] > 0]
G = Gate("#850 · 把 `#807`/`#828` 的判词也重算")
G.asserted("① 前提(跑前写下的最强混淆):标记正则里含它们自己的编号,而当时的页面正是刚写完它们那一版 ⇒ "
           "**谓词必须与页面同代**,谓词也从那次提交的脚本里读",
           bool(all(ST[n]["then"] for n in JOBS)),
           " · ".join(f"#{n} 两代相同 {ST[n]['same']}" for n in JOBS), kind="control")
G.asserted("② 正控:在**当时页面**植入一处「提到主张但无标记」的单位,两种切法都必须 +1",
           bool(pc["旧"] == 1 and pc["新"] == 1), f"旧 +{pc['旧']} · 新 +{pc['新']}", kind="control")
G.asserted("③ 负控:同一处**补上标记**后,两种切法都必须 **+0**(⚠ **这个零该是零**)",
           bool(nc["旧"] == 0 and nc["新"] == 0), f"旧 +{nc['旧']} · 新 +{nc['新']}", kind="control")
G.asserted("④ kill(预注册):「这两条当时没算错」要成立,需**当时页面上旧切法与新切法给出同一个数**",
           bool(not flip), f"两版不一致的条目 {flip or '无'}", kind="kill",
           yardstick="同一页面同一谓词下,两种切法各自的「裸着」计数", yardstick_noise=0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif flip:
    VERD = f"**B 当时就算错了:{flip} ⇒ 那些判词要撤。**"
elif now_bad:
    # ⚠⚠ **第一版这里写「今天有退化 ⇒ 要修的是页面」,而那超出了它自己的输出。**
    #    `#807` 的 `TOP` 正则是 `sexeduc|racmar|顶对|top pair` —— **它匹配的是题名本身,
    #    不是那条主张**;而 `#846`/`#849` 整轮都在正当地谈 `sexeduc`。
    #    ⇒ **43 不是「43 个退化」,是一个未分辨的存量**(`#841`/`#848` 的同一条:
    #      多报的工具会被无视;存量不是缺陷计数)。⇒ 报归属分布,不下「要修」的判词。
    own = {}
    for f in ("README_zh.md", "README.md"):
        _, det = count((ROOT/f).read_text(encoding="utf-8"), ST[807]["now"], JOBS[807]["pairs"], True)
        for k, o, cl in det: own[o] = own.get(o, 0) + 1
    late = sum(v for o, v in own.items() if o >= 846)
    VERD = (f"**A/C 之间,而本轮只有资格说前半句:`#807` 当时没算错(四格里当时两格都是 "
            f"{G4[(807,'当时','旧')]}),今天的计数升到 {G4[(807,'今天','新')]}。**\n"
            f"  ⚠⚠ **但 {G4[(807,'今天','新')]} 不是「{G4[(807,'今天','新')]} 处退化」**:"
            f"`#807` 的 `TOP` 正则是 `sexeduc|racmar|顶对|top pair` —— **它匹配的是题名本身,"
            f"不是那条主张**,而 `#846`/`#849` 整轮都在正当地谈 `sexeduc`。\n"
            f"  实测归属:**`#846` 及之后的条目占 {late} 处**,其余 "
            f"{G4[(807,'今天','新')]-late} 处分布在更早的条目上。\n"
            f"  ⇒ **登记为「未分辨的存量」,不是缺陷计数**(`#841`/`#848` 的同一条)。\n"
            f"  ⇒ **一句关于方法的话:一个用题名做代理的检查,在题名后来成了研究对象的时候,\n"
            f"  会把正当的讨论数成退化 —— 而这不是它坏了,是它的代理从一开始就只在\n"
            f"  「没人再谈这道题」的世界里成立。**")
else:
    VERD = (f"**A 两条判词都不受切法影响,而这现在是量出来的。** "
            f"`#807` 四格 {[G4[(807,e,c)] for e in ('当时','今天') for c in ('旧','新')]} · "
            f"`#828` 四格 {[G4[(828,e,c)] for e in ('当时','今天') for c in ('旧','新')]}。\n"
            f"  ⇒ **`#843` 登记的那句「不是做不到,是我只做了最可能翻的那一条」现在兑现了 ——\n"
            f"  而兑现它花的,和 `#843` 那次一样,是一条 `git show`。**")
print(VERD)
json.dump(dict(grid={f"{n}|{e}|{c}": v for (n, e, c), v in G4.items()},
               commits={n: ST[n]["commit"] for n in JOBS},
               same_generation={n: ST[n]["same"] for n in JOBS},
               pos_control=pc, neg_control=nc, flipped=flip, degraded_today=now_bad,
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), action="Production"),
          open(OUT/"other_two_verdicts.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'other_two_verdicts.json'}")
