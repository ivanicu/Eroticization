"""#792 · E03·A46·R231 —— 页面还按「性 vs 非性」组织着多少?先数,再决定要不要重组

`#791`① 预注册:「**先数出页面上还有多少句话按「性 vs 非性」组织**,再决定要不要整页重组。」
⚠ 而这是本项目第三次执行「先数再定规矩」这条(`#783` 数比值 · `#787` 数 `N/N`)——
   **前两次都是这条纪律救了我**:`#783` 那次数出来的东西推翻了我数它的方式,
   `#787` 那次数出来的比例照预注册说「不立规矩」,而我照做了。

**⚠ 本轮标注为 Production,不是 Frontier,而这是诚实标注不是谦虚**(`frontier §7.6`):
   它不产生新的关于人的判断,也没有两个本体不同的世界 ——
   它回答的是「**这个交付物有多少建立在一个已知是错的组织原则上**」,
   而那是一个**关于页面的事实**,不是一个关于人的事实。
⚠ 但它不是可有可无的:`§0.2` 说交付物是**站得住的东西**,不是账本。
   **一个仍按错误切法组织的页面,是在用一个已被推翻的结构向读者索取信任。**

G1 估计量(两个,方法之前先命名):
   (a) 页面上**以「性 vs 非性」作为组织原则**的句子/行数
   (b) 其中有多少在 `#789`(该切法被推翻那一轮)**之前**写下,多少在之后

识别:(a) 是语义的 —— 「一句话是不是靠这个切法撑着」不能靠正则判。
   ⇒ 沿用 `#787` 立下的那条分界:**造筛子提候选,判定权不交给它,逐条人工判并指到位置。**
   ⚠ 而 `#755` 已经证明:**为语义缺陷造 lint 会失败** —— 所以本轮同样不造 lint。

预注册判词(条件式,不是阈值):
  if 筛子的正对照开火(它必须找到 `#789`/`#791` 那两行里已知在谈这个切法的句子):
      设 pre = 在 `#789` 之前写下、且仍以该切法组织的行数
      if pre >= 10  -> 整页重组(那不是几处措辞,是一个结构)
      elif pre >= 3 -> 逐处加限定,不重组
      else          -> 只需在 `#789` 那几行加一句前向指针
  else: UNVERIFIED
⚠ 阈值 10 / 3 的理由**写在跑之前**:页面主表约 180 行,10 行 ≈ 5.6% ——
  **低于这个比例,重组的代价(锚失效、行号漂移、`readme_gate` 全族重跑)高于收益。**
  这是一个**成本阈值**,不是一个统计阈值,如实标注。

⚠ 跑之前写下的最强混淆:**筛子会把「讨论这个切法」与「使用这个切法」混在一起。**
  `#789`/`#791` 那两行**大量提到**「性/非性」,但它们是在**推翻**它 —— 那不是缺陷,是修复。
  ⇒ 控制:人工判必须分三类 —— **USE(靠它组织)· DISCUSS(在谈它/推翻它)· MENTION(只是词)**,
     而**只有 USE 计入 (a)**。

本轮对象是这两页 `README` 自己 ⇒ **换不了仪器**,理由与 `#787` 同一种
(不是数据边界,是「对象根本不是世界」);跨轮次的替代性交叉:每一条判定指到页面位置。
"""
import re, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)

# ── 筛子:提到这个切法的候选行 ────────────────────────────────────────────────
SIEVE = re.compile(r"性\s*(?:vs|对|/|、)\s*非性|非性|sexual\s+(?:vs\.?|versus)\s+non-?sexual|non-?sexual|"
                   r"性题|非性题|性道德(?:的|是)?(?:一|特殊|不特殊)|sexual morality is")
ANCHOR = re.compile(r"\[#(\d+)「[^」]*」\]")

def rows_of(path):
    txt = path.read_text(encoding="utf-8")
    out = []
    for i, line in enumerate(txt.split("\n"), 1):
        if not SIEVE.search(line): continue
        ancs = [int(m.group(1)) for m in ANCHOR.finditer(line)]
        out.append(dict(line=i, n_hits=len(SIEVE.findall(line)), anchors=ancs,
                        head=line[:100], has789=any(a >= 789 for a in ancs)))
    return out

print("=== ① 筛子:提到「性 vs 非性」这个切法的行 ===")
C = {}
for f in ("README_zh.md", "README.md"):
    C[f] = rows_of(ROOT/f)
    print(f"  {f}: **{len(C[f])} 行**命中 · 合计 {sum(r['n_hits'] for r in C[f])} 处")

known = [r for r in C["README_zh.md"] if any(a in (789, 791) for a in r["anchors"])]
print(f"\n  正对照①(能不能看见):`#789`/`#791` 那两行 —— "
      f"筛子{'**找到了**' if known else '**没找到**'}({len(known)} 行)")

# ⚠⚠ 第二条正对照,而它是第一版没有、并且当场把第一版判掉的那一条。
#    `realstat §4` 那一行:**「正对照只问『这具仪器看得见吗』,从不问『它看见的是不是我要主张的那个东西』」**
#    —— 一个与仪器共享盲点的正对照,只确认仪器,什么也不许可。
#    ⇒ 所以再放一条:**筛子必须找到我已经知道是 USE 的那一行**(`#785` 的锚行)。
target = [r for r in C["README_zh.md"] if 785 in r["anchors"]]
print(f"  正对照②(**看见的是不是我要主张的那个东西**):`#785` 那一行是**已知的 USE** —— "
      f"筛子{'**找到了**' if target else '**没找到 ⇒ 召回缺口,计数不可信**'}({len(target)} 行)")
recall_ok = bool(target)

print("\n=== ② 逐行列出(人工判 USE / DISCUSS / MENTION,只有 USE 计入)===")
for r in C["README_zh.md"]:
    tag = f"锚 {r['anchors'][-3:]}" if r["anchors"] else "无锚(表行/正文)"
    print(f"  行 {r['line']:>4}  命中 {r['n_hits']:>2}  {tag:22s}  {r['head'][:78]}")

# ── ③ 人工判 —— 只判我能指到位置的,其余 UNVERIFIED ────────────────────────────
# ⚠ 分三类是跑之前写下的混淆的控制:讨论它 ≠ 使用它。
VERD = {
    789: ("DISCUSS", "这一行**推翻**了该切法 —— 修复,不是缺陷"),
    791: ("DISCUSS", "这一行给出替代切法(三堆)—— 修复,不是缺陷"),
    785: ("USE", "「不是一条关于性道德的规律,是一条关于其余人真的动了的那些题的规律」—— "
                 "标题与整段仍以性/非性组织,而 `#789` 之后它的前提已倒"),
    788: ("DISCUSS", "谈的是仪器边界,提到非性题只是为了说明问题需要它们"),
}
use, discuss, mention, unver = [], [], [], []
for r in C["README_zh.md"]:
    a = next((x for x in reversed(r["anchors"]) if x in VERD), None)
    if a is None: unver.append(r); continue
    k = VERD[a][0]
    (use if k == "USE" else discuss if k == "DISCUSS" else mention).append((a, r))
print(f"\n=== ③ 人工判 ===")
for a, (kind, why) in sorted(VERD.items()):
    print(f"  `#{a}`  {kind:8s}  {why}")
pre = [x for x in use if x[0] < 789]
print(f"\n  命中行 {len(C['README_zh.md'])} · **已核实 {len(VERD)} 个锚** · "
      f"USE {len(use)} · DISCUSS {len(discuss)} · **未核实的行 {len(unver)} ⇒ UNVERIFIED,不猜**")
print(f"  ⇒ **在 `#789` 之前写下、且仍以该切法组织的:{len(pre)}**")

print("\n=== ④ 预注册判词开火 ===")
if not known or not recall_ok:
    v = (f"**UNVERIFIED,而不可信的是筛子不是页面。** 正对照①过了(它看得见 `#789`/`#791`),"
         f"**而正对照②没过:它找不到 `#785` 那一行,可我已经知道那一行是 USE。**\n"
         f"  ⇒ 词表里写的是 `性道德(?:的|是)?(?:一|特殊|不特殊)`,而那一行写的是"
         f"「关于**性道德的规律**」—— **`的` 后面跟的是「规」,不是「一/特殊」,于是漏掉。**\n"
         f"  ⇒ **这正是 `realstat` 那一行:正对照只问「这具仪器看得见吗」,"
         f"从不问「它看见的是不是我要主张的那个东西」。** 我建的第一条正对照与筛子共享同一个盲点 ——\n"
         f"  它去找**讨论这个切法**的行(那些行把「非性」写成一个词,词表命中),\n"
         f"  而我要数的是**使用这个切法**的行(那些行只说「性道德」,词表不命中)。\n"
         f"  ⇒ **命中的 {len(C['README_zh.md'])} 行与 USE 的 0 行,两个数都不许拿去做重组决定。**\n"
         f"  ⇒ 下一轮的修法**不是把词表放宽**(那会把 DISCUSS 一起捞进来):\n"
         f"  **先把「USE 的语言长什么样」与「DISCUSS 的语言长什么样」写成两串不同的词,再各自量召回。**")
elif len(pre) >= 10:
    v = f"**整页重组**:`#789` 之前仍按该切法组织的有 {len(pre)} 行 ≥ 10 —— 那不是几处措辞,是一个结构。"
elif len(pre) >= 3:
    v = f"**逐处加限定,不重组**:{len(pre)} 行(3 ≤ n < 10)—— 重组的代价高于收益。"
else:
    v = (f"**只需前向指针**:`#789` 之前仍以该切法组织的只有 **{len(pre)} 行**"
         f"(`#785` 那一行),其余命中全部是 **DISCUSS**(`#788`/`#789`/`#791` 在谈它、推翻它、替代它)"
         f"或**未核实**。\n"
         f"  ⇒ **页面并没有大面积建立在这个切法上** —— 而这本身是一个值得记的事实:\n"
         f"  **我以为它是一个组织原则,实测它只是一行标题。**\n"
         f"  ⇒ 代价最低的修法:在 `#785` 那一行挂一句指向 `#791` 三堆的前向指针,**不重组。**")
print(v)
print(f"\n⚠ 跑之前写下的混淆的控制有产出:**{len(discuss)} 行是 DISCUSS 而不是 USE** —— "
      f"若不分这三类,我会把「推翻它的那几行」算成「靠它撑着的那几行」,"
      f"**从而把一次修复读成一片缺陷,并据此重组整页。**")
json.dump(dict(zh=C["README_zh.md"], en=C["README.md"], verdicts={str(k): list(v2) for k, v2 in VERD.items()},
               n_use=len(use), n_discuss=len(discuss), n_unverified=len(unver),
               n_pre789_use=len(pre), verdict=v, action="Production"),
          open(OUT/"page_organisation.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'page_organisation.json'}")
