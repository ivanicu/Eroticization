"""#787 · E03·A44·R226 —— 页面上有几个 `N/N`,它的分母是**活下来的格**而不是**尝试过的格**?

`#786` 逮到 `#785` 表里的一行:**「`premarsx` 可读 2 格 · 全在 1.0 之下 2/2」** ——
读起来像一致,**而实测四个规格里只有 1 个可读**,那 2 格是同一规格的两种连接函数。
⇒ **「2/2」的分母是活下来的格,不是尝试过的格。**
`#786`② 预注册:**这是一个可机械检出的形状 —— 先数,再决定要不要一条 lint。**

⚠⚠ **而「先数」这件事本身有一条已经付过学费的规矩(`#783`)**:
   抽取器是一具仪器,必须**同时报 recall 与 precision**,并且**手标名单写在抽取器跑之前**。
   ⚠ 还有一条更硬的(`#755`):**为语义缺陷造 lint 会失败。**
   「分母是不是尝试过的格」是**语义的** —— 所以本轮**不造 lint**,
   只造一具**候选筛子**(`N/N` 这个语法形状),然后**逐个人工判**,并把命中率报出来。
   **筛子只提出候选,判定权不交给它** —— 这正是 `#775` 那条 lint 与 `#755` 那条被否掉的 lint 的分界。

G1 估计量(两个,方法之前先命名):
   (a) 页面上 `N/N` 形状的**总数**
   (b) 其中**分母是活下来的格**的个数,以及**分母是尝试过的格**的个数(= 筛子的假阳性)

识别:(a) 机械可判;(b) **只能人工判**,而人工判的是我自己写的页面 ⇒ **它是我自己的仪器**。
   ⇒ 所以本轮报的不是「真相」,是**一张逐条可复核的表**,每一条都指到页面上的位置。

预注册判词(条件式):
  if 筛子的正对照开火(它必须找到 `#785` 那个已知的 `2/2`):
      if 分母是活下来的格 的比例 >= 0.5 -> 这是一个族,值得一条**写作规矩**(不是 lint)
      else                              -> 那是孤例,`#786` 修掉它就够了,不立规矩
  else: UNVERIFIED(筛子找不到已知的那一个 ⇒ 它是瞎的)

⚠⚠ **本轮换不了仪器,而它的理由与 `R223` 那次是两种不同的东西 —— 这一点必须分清,
   否则同一句豁免语会掩盖两种完全不同的处境:**
   · `R223` 的豁免:**对象是世界**(美国人的性道德随年代的变化),而承载它的第二具仪器
     在本机六份外部数据里不存在 —— 那是一条**关于数据的、可失效的**声明,
     由 `R223/instrument_search.py` **跑出来**,加一份多波调查它立刻失效。
   · **本轮的豁免:对象根本不是世界,是这两页 `README` 自己。**
     一个「页面上有几处 N/N」的问题**没有第二具仪器可换**,因为它问的不是自然,是我写的字。
   ⇒ 而本轮能做、也做了的**替代性交叉**是:每一条判定都必须落到**另一个轮次已持久化的产物**上
     (`R221`/`R224`/`R225` 的 json),**不许只靠读页面**。判定表里每一行都带着它的产物出处。
   ⚠ **这不是跨仪器,是跨轮次,而它弱得多** —— 那三个产物是同一具 GSS、同一个我造的。如实标注。

⚠ 跑之前写下的最强混淆:**`N/N` 里绝大多数是无辜的** —— 「8/8」「4/4」在本项目里通常是
   「尝试了 8 格,8 格都过」,那是**正确写法**。⇒ 筛子的 precision 必然低,
   **而低 precision 正是本轮要量的东西**,不是缺陷。
"""
import re, json, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)

# ── 筛子:`N/N`(两数相等),排除日期、版本号、分数形式的比例 ──────────────────────
SIEVE = re.compile(r"(?<![\d./])(\d{1,3})\s*/\s*(\d{1,3})(?![\d./])")

def candidates(path):
    txt = path.read_text(encoding="utf-8")
    out = []
    for m in SIEVE.finditer(txt):
        a, b = int(m.group(1)), int(m.group(2))
        if a != b or a == 0: continue
        s = max(0, m.start()-160); e = min(len(txt), m.end()+60)
        out.append(dict(value=f"{a}/{b}", pos=m.start(), ctx=txt[s:e].replace("\n", " ")))
    return out

print("=== ① 筛子:页面上所有 `N/N`(两数相等)===")
cand = {}
for f in ("README_zh.md", "README.md"):
    cand[f] = candidates(ROOT/f)
    print(f"  {f}: **{len(cand[f])}** 处")

# ── 正对照:筛子必须找到 `#785` 那个已知的 `2/2` ────────────────────────────────
known = [c for c in cand["README_zh.md"] if c["value"] == "2/2" and "premarsx" in c["ctx"]]
print(f"\n  正对照:`#785` 那个已知的 `premarsx … 2/2` —— 筛子{'**找到了**' if known else '**没找到 ⇒ 它是瞎的**'}"
      f"({len(known)} 处)")

print("\n=== ② 逐条列出,人工判(筛子只提候选,判定权不交给它 —— `#755`)===")
for i, c in enumerate(cand["README_zh.md"], 1):
    print(f"\n  [{i}] `{c['value']}` @ {c['pos']}")
    print(f"      …{c['ctx'][-190:]}")

json.dump(dict(zh=cand["README_zh.md"], en=cand["README.md"],
               n_zh=len(cand["README_zh.md"]), n_en=len(cand["README.md"]),
               positive_control_found=len(known)),
          open(OUT/"nn_candidates.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  候选表 → {OUT/'nn_candidates.json'}")
print("\n⚠ 到此为止只**提出候选**;判定在下面的 ③,而它是**人工**的,逐条指到产物。"
      "两段分开写,是为了让「机器筛的」与「我判的」在输出与产物里都分得开。")

# ─────────────────────────────────────────────────────────────────────────────
# ③ 人工判(只判我能从产物里核实的那些;其余如实标 UNVERIFIED,不猜)
# ─────────────────────────────────────────────────────────────────────────────
VERDICTS = {
    "premarsx 2/2":   ("SURVIVOR", "R224 产物:premarsx 四个规格里只有 1 个可读,"
                                   "那 2 格是同一规格的两种连接 ⇒ 分母是活下来的格。方向:**显得一致**(讨好)"),
    "1975–94 4/4":    ("SURVIVOR", "R221 产物:该世代 8 格设计里只有 4 格算得出来 ⇒ 分母是算得出的格。"
                                   "方向:**显得全军覆没**(不讨好)"),
    "1980–99 2/2":    ("SURVIVOR", "R221 产物:8 格设计里只有 2 格算得出来。方向:**不讨好**"),
    "homosex 8/8":    ("ATTEMPTED", "R224 产物:2 分层 × 2 统计量 × 2 连接 = 8 格全部尝试且全部可读"),
    "homosex 4/4 ×3": ("ATTEMPTED", "R225 产物:三个窗口各 4 格全部尝试且全部可读"),
    "1965–79 6/8":    ("ATTEMPTED", "分母写的就是 8 —— **这一处写法是对的**,可作为正确写法的样例"),
}
n_surv = sum(1 for v in VERDICTS.values() if v[0] == "SURVIVOR")
n_att = sum(1 for v in VERDICTS.values() if v[0] == "ATTEMPTED")
print("\n=== ③ 人工判(逐条可复核,每条指到产物)===")
for k, (kind, why) in VERDICTS.items():
    print(f"  {kind:9s}  {k:16s}  {why}")
n_cand = len(cand["README_zh.md"])
print(f"\n  候选 {n_cand} · **已核实 {len(VERDICTS)}** · 其中分母是活下来的格 **{n_surv}** · "
      f"分母是尝试过的格 {n_att} · **其余 {n_cand-len(VERDICTS)} 条本轮没核实 ⇒ UNVERIFIED,不猜**")
print(f"  ⚠ **两个 SURVIVOR 指向相反方向**:`premarsx` 那个让结论显得一致(讨好),"
      f"世代那两个让不可读显得全军覆没(不讨好)⇒ **这不是一个有偏的毛病,是一个写法的毛病。**")
print(f"  ⚠ 而这一条是从 `realstat` 那张表抄来的纪律:**「它们全都在讨好我」本身是一个叙事主张,"
      f"要被计数推翻** —— 这里就被推翻了。")

# ── 预注册判词 ────────────────────────────────────────────────────────────────
print("\n=== ④ 预注册判词开火 ===")
rate = n_surv/max(1, n_cand)
fires_family = rate >= 0.5
print(f"  预注册写的是:「分母是活下来的格 的比例 >= 0.5 ⇒ 立一条写作规矩」。"
      f"实测 {n_surv}/{n_cand} = **{rate:.3f} < 0.5 ⇒ 按预注册:孤例,不立规矩。**")
print("  ⚠⚠ **而我现在要说的是:那个判据我写错了,并且我不推翻它。**")
print("     它问的是**筛子的精度**,而决定该不该立规矩的是**真实实例的个数**"
     f"({n_surv} 处,来自 2 个不同的轮 —— R221 两处、R224 一处,方向还相反)—— 这两件事我在跑之前把它们混成了一个数。")
print("     ⇒ 这正是「判词分支测错问题」那一族(`#728`·`#748`·`#750`·`#758`·`#782`)的又一次,"
     "**而这次它发生在判据本身,不是在分支上。**")
print("     ⇒ **本轮照预注册执行:不立规矩。** 已知的两处照常修(修缺陷不是立规矩)。"
     "要不要一条写作规矩,留给一个把判据写对的轮次。")

json.dump(dict(zh=cand["README_zh.md"], en=cand["README.md"], n_zh=n_cand,
               verdicts={k: list(v) for k, v in VERDICTS.items()},
               n_survivor=n_surv, n_attempted=n_att, n_unverified=n_cand-len(VERDICTS),
               preregistered_rate=rate, fires_family=fires_family,
               criterion_was_misspecified=True),
          open(OUT/"nn_candidates.json", "w"), ensure_ascii=False, indent=1)
