#!/usr/bin/env python3
"""#847① —— 把「作用域」从我记得变成机械检查。

**为什么**:`#846` 把「不信教的人走得更快」收窄到 `homosex` 一道题(八题实测,
世俗层排零 **1/8**、虔诚层 **1/8**)。而 `#839`/`#840` 也**只跑过 `homosex`**,
页面上却没说 —— `#847` 花了**三次**才把注记贴全,三次的错法是同一种:
**按我记得的句子去找,而不是按结构把单位枚举出来。**
⇒ **`#836`①/`#840`① 的同一句:教训不落到工具上就是没落地。**

**检测什么(P6 代理账):**
  PROPERTY    页面上一处引用了 `homosex`-only 的结论,却没说它只到那一道题
  PROXY       一个**单位**(`lib/page_units.units()` 定的表行/叙事段)引用了
              `#838`/`#839`/`#840`/`#846`,而该单位里**不含任何作用域词**
  IMPLICATION 只有一个方向可靠:**命中 ⇒ 该单位确实引用了它们且确实没写作用域**(可靠)。
              反过来**不成立**:写了 `homosex` 三个字不等于作用域说对了。
              **只报命中,从不报「本单位的作用域正确」。**
  SAFE SIDE   报「引用了但没写作用域」,由作者补;**从不自动改写页面。**

⚠⚠ **一个已知且量出来的误报类别,写在这里而不是藏起来:**
  `#840` 这个编号在页面上有**两个所指** —— ① 它关于那条缝的**结论**(只到 `homosex`),
  ② 它那次**手抄数字的事故**(`0.6266 ← 0.626612`,`#841`/`#842`/`#843`/`#844` 都在引它)。
  第一版试图用「欠账写成 `#NNN①`、结论写成裸 `#NNN`」把两者分开,**而实测分不开**:
  那 4 条(×2 版 = 8 处)在引事故时**也用裸 `#840`**。
  ⇒ **本仪器无法区分这两个所指,基线因此定在 8。**
  ⚠ **那 8 处不是「8 个错」** —— 它们是**已知误报**;
  棘轮只拦**新增**的,而新增一处只需在句子里写一个 `homosex` 或「作用域」就能消掉。
  **(`#841` 的教训:多报的工具最终会被无视,那等于没有工具 —— 所以误报要被命名,不是被容忍。)**

⚠ **总体由结构定,不由词表定** —— 这正是 `#847` 三次失败的原因,也是 `#842` 修 `units()` 的目的。
⚠ **棘轮而非零容忍**:历史单位里的引用不该让今天的提交失败(`L81`)。
"""
import json, pathlib, re, sys, argparse
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from lib.page_units import units

ROOT = pathlib.Path(__file__).resolve().parents[1]
BASE = ROOT/"tools"/"scope_baseline.json"
SCOPED = (838, 839, 840, 846)                       # 结论只到 `homosex` 的那几条
# ⚠⚠ **第一版用裸串 `"#840" in body`,于是 8 处全是同一种误报:**
#    那些单位引的是 **`#840`①**(那笔「绝不手抄上一轮的数字」的**欠账 id**),
#    **不是 `#840` 关于那条缝的结论。** 一个编号在这个项目里有两个所指
#    (`#170b` 那类「同一个串两个所指」)—— 而它们的写法恰好不同:
#    **欠账一律写成 `#NNN①`,结论写成裸 `#NNN`。**
#    ⇒ 用「后面不跟圈码」把两者分开,这是从对象读出来的区别,不是我发明的约定。
REFS = re.compile(r"#(?:" + "|".join(str(n) for n in SCOPED) + r")(?![①②③④⑤⑥⑦⑧⑨\d])")
SCOPE_WORDS = ("homosex", "这一道题", "这道题", "one item", "single item",
               "作用域", "scope", "1/8")

def scan(pages=("README_zh.md", "README.md")):
    hits = []
    for f in pages:
        p = ROOT/f
        if not p.exists():
            continue  # Skip missing pages (e.g., archived README_zh.md)
        t = p.read_text(encoding="utf-8")
        for kind, owner, body in units(t):
            if owner in SCOPED: continue            # 它自己那一条不必自证
            m = REFS.findall(body) if hasattr(REFS, "findall") else []
            if not REFS.search(body): continue
            if any(w in body for w in SCOPE_WORDS): continue
            hits.append(dict(file=f, kind=kind, owner=owner,
                             refs=sorted(set(REFS.findall(body))), head=body[:90].replace("\n", " ")))
    return hits

def controls():
    """★P5:先证它会开火,再证它不会对合规写法开火。"""
    # ⚠ 第一版这里写成 `⟨无作用域⟩`,而那个串**包含 `作用域` 这个词** ——
    #   于是正控的样本被判为「已写作用域」,**正控自己把自己解除了武装**。
    #   `feedback_check_encodes_instance_not_property` 的同一形状:夹具决定了结论。
    bad = "| 甲引用了 #839 的结论 `[#900「x」]`(Entry 900) ⟨没有说范围⟩|\n"
    good = "| 乙引用了 #839,而它只在 homosex 这一道题上成立 `[#901「y」]`(Entry 901) ⟨⟩|\n"
    none = "| 丙什么都没引用 `[#902「z」]`(Entry 902) ⟨⟩|\n"
    debt = "| 丁只引了欠账 #840① 那条工具 `[#903「w」]`(Entry 903) ⟨⟩|\n"
    def n(txt):
        c = 0
        for kind, owner, body in units(txt):
            if owner in SCOPED: continue
            if not REFS.search(body): continue
            if any(w in body for w in SCOPE_WORDS): continue
            c += 1
        return c
    return dict(pos=n(bad) >= 1, neg_scoped=n(good) == 0, neg_unrelated=n(none) == 0,
                neg_debt_id=n(debt) == 0,
                ok=bool(n(bad) >= 1 and n(good) == 0 and n(none) == 0 and n(debt) == 0))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--precommit", action="store_true")
    ap.add_argument("--rebaseline", action="store_true"); a = ap.parse_args()
    c = controls()
    print(f"仪器控制:正控(引用而无作用域必须命中)**{c['pos']}** · "
          f"负控α(引用且写了作用域)**{c['neg_scoped']}** · 负控β(什么都没引用)**{c['neg_unrelated']}** · "
          f"**负控γ(只引欠账 id `#840①`)**{c['neg_debt_id']}**")
    if not c["ok"]:
        print("⛔ **控制没过 ⇒ 本次扫描不可采**(★P5:仪器没证明会开火,它的零是沉默)"); sys.exit(2)
    h = scan()
    print(f"引用了 `#838`/`#839`/`#840`/`#846` 却**没写作用域**的单位:**{len(h)}** 个")
    for x in h[:10]: print(f"  {x['file']} [{x['kind']} #{x['owner']}] {x['refs']} …{x['head']}…")
    if len(h) > 10: print(f"  …另 {len(h)-10} 个")
    print("⚠ **只报命中,从不报「本单位的作用域正确」** —— 写了 `homosex` 不等于作用域说对了。")
    if a.rebaseline:
        BASE.write_text(json.dumps({"count": len(h)}, indent=1)); print("基线已写"); sys.exit(0)
    if a.precommit:
        if not BASE.exists():
            BASE.write_text(json.dumps({"count": len(h)}, indent=1))
            print(f"首次建立基线 = {len(h)}"); sys.exit(0)
        old = json.loads(BASE.read_text())["count"]
        if len(h) > old:
            print(f"🔒 PRE-COMMIT BLOCK:无作用域的引用 {old} -> {len(h)}(+{len(h)-old})"); sys.exit(1)
        if len(h) < old:
            BASE.write_text(json.dumps({"count": len(h)}, indent=1)); print(f"棘轮收紧:{old} -> {len(h)}")
    sys.exit(0)
