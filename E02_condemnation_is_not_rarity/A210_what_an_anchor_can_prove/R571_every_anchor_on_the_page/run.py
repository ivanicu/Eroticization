"""E02·A210·R571 — 页面上每一个锚:有多少个根本无法被检查?

`#526` 的 NEXT,**留在 `A210` 内**,同样一轮打包多个操作:
识别性(逐锚打印所引记号)+ 主检验(可检查率 × 碰撞率)+ 三道对照 + 逐格规格曲线 + 页面兑现。

行动类型:**FRONTIER**。`#524d` 留下两个洞,`#526` 只给第一个洞(语义)一个数;
本轮量**第二个洞**:**一个不含可核验记号的断言,它的锚是不可检查的** —— 到底有几个?

G1 ESTIMAND(先于方法),两个量,**都先于方法命名**:
  ① **可检查率** = 页面上带有至少一个可核验记号(数字或引号内短语)的锚 / 全部锚
  ② **零碰撞率** = 在可检查的锚里,其记号在账本中**只出现在被引条目**的锚 / 可检查的锚
  记号的定义(预注册,先于看结果):数字 = 正则 `−?\d+(?:\.\d+)?%?`,取该锚**同段落**内、
  长度 ≥ 3 字符的;短语 = 该锚 `[...]` 括号内成对引号里的串。

WORLDS:
  W-CHECKED   页面基本可检查:可检查率 ≥ 0.8 且零碰撞率 ≥ 0.8
  W-DECORATIVE 锚大体是装饰:可检查率 < 0.5 ⇒ **页面上大多数「出处」并不能被任何机器核对**
  W-WEAK      可检查但撞得厉害:可检查率高、零碰撞率 < 0.5
⚠ **W-CHECKED 的阳性结局是我欢迎的,所以它不是那个该被设计的步。**
  按 BASIN RULE,本轮真正下注的是 **W-DECORATIVE** —— 它意味着我上一轮刚在页面上
  夸口的那套「每条断言都带出处」大部分**不可核对**。预先写明:**我不希望它为真。**

CONTROLS(G2):
  正对照:把一个可检查锚的记号改掉一个字符 -> 必须被抓(g=0 时必须不报)
  安慰剂:一个字面不存在的记号 -> 碰撞数必须恰为 0(「这个零该不该是零?」该 ⇒ negative_control)
  ⚠ 第三道,针对本轮自己的仪器:**记号抽取器必须能在一个人造的、显然可检查的锚上抽到记号**
    —— 否则「可检查率低」测的是我的正则,不是页面。
KILL(条件式,预注册):
  if 篡改被抓 and g=0 不报 and 安慰剂为 0 and 抽取器在人造锚上抽到记号:
      按上面三个世界判
  else: UNVERIFIED
IMPOSSIBLE:可检查 ≠ 引对(语义仍不可自动验,`#526e`)· 单页面 ⇒ 无跨站点 · [unchallenged]
"""
import os, sys, pathlib, json, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
LED = (ROOT / "RETRACTIONS.md").read_text()
mk = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', LED, re.M)]
BODY = {n: LED[s:(mk[i+1][1] if i+1 < len(mk) else len(LED))] for i, (n, s) in enumerate(mk)}
AUDIT = {526}   # 本线自己的审计条目,不计入被测总体(见 `#527b`)
NUM = re.compile(r'[−-]?\d+(?:\.\d+)?%?')
PHR = re.compile(r'[「"“]([^「」"“”]{4,})[」"”]')

def tokens(para, bracket):
    """预注册的记号抽取:括号内引号短语优先,其次同段落内 ≥3 字符的数字。"""
    ph = [p for p in PHR.findall(bracket)]
    nu = [t for t in NUM.findall(para) if len(t.strip("−-")) >= 3]
    return ph, nu

def scan(page):
    t = pathlib.Path(page).read_text()
    out = []
    for m in re.finditer(r'`\[([^\]]*#\d+[^\]]*)\]`', t):
        br = m.group(1)
        # ⚠ 修:原版用空行切段,而项目符号列表整块无空行 -> 每个锚都拿到列表首个数 `−0.29`。
        # 记号必须来自**锚自己那一行**,否则测的是列表,不是锚(硬规则 1:在标题处下结论)。
        ls = t.rfind("\n", 0, m.start()) + 1; le = t.find("\n", m.end())
        para = t[ls:(le if le > 0 else len(t))]
        for e in re.findall(r'#(\d+)', br):
            ph, nu = tokens(para, br)
            out.append(dict(page=page, entry=int(e), bracket=br[:60], phrases=ph, numbers=nu[:6]))
    return out

rows = scan("README.md") + scan("README_zh.md")
print(f"页面锚引用总数 = {len(rows)}  (README {len(scan('README.md'))} · README_zh {len(scan('README_zh.md'))})")
print("\n=== 规则①:逐锚打印它带的可核验记号 ===")
for r in rows:
    e, B = r["entry"], BODY.get(r["entry"], "")
    toks = [t for t in r["phrases"]] or [t for t in r["numbers"]]
    r["kind"] = "phrase" if r["phrases"] else ("number" if r["numbers"] else "none")
    r["token"] = toks[0] if toks else None
    r["checkable"] = bool(toks)
    if r["checkable"]:
        # 记号必须先在被引条目里 —— 否则锚本身就是错的(`#526d`)
        r["in_entry"] = sum(1 for t in toks if t in B)
        # ⚠ 修:`#526` 的表格逐字引用了全部四个短语,于是**审计的记录本身**制造了碰撞。
        # 一个仪器把自己的记录算进被测总体,它测的就不再是页面。排除审计条目,并把它单独报出来。
        best = min((sum(1 for n in BODY if n != e and n not in AUDIT and t in BODY[n]), t) for t in toks)
        r["n_collide"], r["token"] = best[0], best[1]
        r["rate"] = round(best[0] / (len(BODY) - 1), 5)
    else:
        r["in_entry"], r["n_collide"], r["rate"] = 0, None, None
    r["n"] = len(BODY) - 1
    r["inclusion"] = [f"账本 {len(BODY)} 条除 #{e}", "预注册记号规则:引号短语优先,其次 ≥3 字符数字"]
    tag = "不可检查 ⛔" if not r["checkable"] else (f"碰撞 {r['n_collide']:3d}/{r['n']} = {r['rate']:.4f} "
                                                f"{'✅' if r['n_collide'] == 0 else '⛔'}")
    print(f"  {r['page'][:9]:9s} #{e:3d} {r['kind']:6s} {str(r['token'])[:26]:28s} {tag}")

chk = [r for r in rows if r["checkable"]]
zero = [r for r in chk if r["n_collide"] == 0]
CR, ZR = len(chk) / len(rows), (len(zero) / len(chk) if chk else 0.0)
print(f"\n  **可检查率 = {len(chk)}/{len(rows)} = {CR:.4f}**   "
      f"**零碰撞率 = {len(zero)}/{len(chk)} = {ZR:.4f}**")

# ---- 对照
G = Gate("页面上每一个锚:有多少个根本无法被检查?")
def nfail(led):
    b = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
    bd = {n: led[s:(b[i+1][1] if i+1 < len(b) else len(led))] for i, (n, s) in enumerate(b)}
    return sum(1 for r in chk if r["token"] not in bd.get(r["entry"], ""))
clean = nfail(LED)
victim = next(r for r in chk if r["kind"] == "phrase")
tamp = nfail(LED.replace(victim["token"], victim["token"][:-1] + "x", 1))
print(f"\n=== 对照 ===\n  g=0 未篡改 -> 失败 {clean}(必须 0) · 篡改 {victim['token'][:14]!r} -> 失败 {tamp}(必须 ≥1)")
G.positive_control("篡改一个记号必须被抓", planted=float(tamp), floor=0.5, spread=1e-9)
G.negative_control("g=0:未篡改必须一个都不报", null=float(clean), effect=float(max(tamp, 1)),
                   null_spread=1e-9, null_kind="不做任何篡改的同一份账本")
PL = "绝不存在的记号-qz93"
G.negative_control("安慰剂:不存在的记号,碰撞必须恰为 0",
                   null=float(sum(1 for n in BODY if PL in BODY[n])), effect=1.0,
                   null_spread=1e-9, null_kind="字面不存在的串")
# 第三道:抽取器自检 —— 人造的、显然可检查的锚必须抽到记号
fake_ph, fake_nu = tokens("这一格是 −0.1234 的人造段落 `[#999「人造短语记号」]`", "#999「人造短语记号」")
G.positive_control("抽取器自检:人造锚上必须抽到记号",
                   planted=float(len(fake_ph) + len(fake_nu)), floor=0.5, spread=1e-9)
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['page'][:2]}#{r['entry']}": r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {f"{r['page'][:2]}#{r['entry']}": r for r in rows})

print("\n" + "=" * 72)
ok = tamp >= 1 and clean == 0 and (len(fake_ph) + len(fake_nu)) > 0
if ok:
    if CR < 0.5:
        world, verdict = "W-DECORATIVE", f"**可检查率 {CR:.4f} < 0.5 -> 页面上大多数出处不可被任何机器核对**"
    elif CR >= 0.8 and ZR >= 0.8:
        world, verdict = "W-CHECKED", f"可检查率 {CR:.4f}、零碰撞率 {ZR:.4f} -> **页面基本可检查**"
    else:
        world, verdict = "W-WEAK", f"可检查率 {CR:.4f} 但零碰撞率 {ZR:.4f} -> **可检查而不可分辨**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:可检查 ≠ 引对。一个记号在被引条目里、又不与别处碰撞的锚,"
          "仍可能指向一条讲别的事的条目 —— **语义不可自动验(`#526e`),本轮没有改变这一点。**")
else:
    world, verdict = "UNVERIFIED", f"控制未齐 clean={clean} tamp={tamp} extractor={len(fake_ph)+len(fake_nu)}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(n_anchor_refs=len(rows), checkable_rate=CR, zero_collision_rate=ZR,
               n_entries=len(BODY), rows=rows, world=world, verdict=verdict,
               controls=dict(clean=clean, tampered=tamp, extractor_selftest=len(fake_ph) + len(fake_nu)),
               estimand=["可检查率", "零碰撞率"], instrument="纯子串匹配 + 预注册记号抽取正则",
               impossible=["可检查 ≠ 引对(语义不可自动验)", "单页面无跨站点", "无干预非因果"],
               unchallenged=True), open(OUT / "every_anchor.json", "w"), indent=1)
print(f"\nwrote {OUT/'every_anchor.json'}")
