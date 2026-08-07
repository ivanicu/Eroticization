#!/usr/bin/env python3
"""#831① —— `debts_gate` 保证 `DEBTS.tsv` **自洽**,从不保证它**被更新过**。

**缺口**(`#831` 记下,一直 `OPEN`):账本里一条 `NEXT` 写下 `#844①`,
**而表里可以根本没有这一行,`debts_gate` 一个字都不会说** ——
它只检查表内部的一致性(状态合法、`SETTLED` 有 `settled_in`、无重复 id)。
⇒ **一张永远自洽的空表,和一张自洽且完整的表,在那道闸下长得一模一样。**

**检测什么(P6 代理账):**
  PROPERTY    账本里被承诺过的欠账,表里都有一行
  PROXY       `RETRACTIONS.md` 里形如 `#NNN①` 的**欠账 id**,与 `DEBTS.tsv` 第一列做差集
  IMPLICATION 只有一个方向可靠:**差集非空 ⇒ 那个 id 确实在账本里出现过而表里没有**(可靠)。
              反过来**不成立**:表里有行不证明那笔欠账被认真对待过。
              **只报缺失,从不报「表是完整的」。**
  SAFE SIDE   报「账本提过而表里没有」,由作者补行;**从不自动写表。**

⚠ **方向是单向的,这一点必须写死**:表里有而账本没提**不算缺失**
  (欠账可以先建行、后写账;也可能来自更早的、已归档的叙述)。
⚠⚠ **零容忍,但只在表的作用域之内 —— 而这个分界是从对象推出来的,不是我挑的:**
  `DEBTS.tsv` **建于 `#830`**,并只导入了 `#800` 那份「26 笔 `UNVERIFIED`」的清单。
  实测:账本里被承诺过的欠账共 **172** 笔,表里 **42** 行,**差 130 笔** ——
  而其中绝大多数在 `#830` **之前**,表从来就没打算覆盖它们(`#800` 已如实登记它们
  「状态真的不知道」,`#829` 又证过从散文反推不可行)。
  ⇒ **`#830` 及之后被承诺的欠账:零容忍**(表存在就是为了它们);
    **`#830` 之前的:登记为已知历史存量,只报数不阻断**(`L81`)。
  ⚠ **那 130 不是「130 个没还的债」** —— 是**未分辨的历史存量**
  (`#841` 的 357、`#848` 的 8、`#850` 的 39,同一条)。

⚠ **本工具的第一版把这个数报成 3。** 正则写成「圈码紧跟数字」,
  而账本的写法是 `` `#844`① `` —— **中间隔着一个反引号**。
  **正控证明了它会开火,却没证明它看得见页面真正用的那种写法**
  (`realstat`:「正控只问『这仪器看得见吗』,从不问『它看见的是不是我要主张的东西』」)。
  **本会话第四次同一形状。** ⇒ 紧式/宽式两个计数现在都印出来,永远并列。
"""
TABLE_BORN = 830
import pathlib, re, sys, argparse

ROOT = pathlib.Path(__file__).resolve().parents[1]
# ⚠⚠ **第一版写成 `#(\d{3,4})([①-⑨])`,要求圈码**紧跟**数字 —— 而账本里的写法是
#    `` `#844`① ``:**数字和圈码之间隔着一个反引号。** 于是它只找到 **3** 笔,
#    而本会话一轮就写过六笔。**正控证明了它会开火,却没证明它看得见页面真正用的那种写法** ——
#    `realstat` 点名的那条:「正控只问『这仪器看得见吗』,从不问『它看见的是不是我要主张的东西』」。
#    **本会话第四次同一形状。** ⇒ 允许中间有反引号/空格,并把两种写法的计数都印出来。
ID = re.compile(r"#(\d{3,4})`?\s*([①②③④⑤⑥⑦⑧⑨])")
ID_TIGHT = re.compile(r"#(\d{3,4})([①②③④⑤⑥⑦⑧⑨])")

def promised():
    """账本里出现过的欠账 id(带圈码的才算 —— 裸 `#NNN` 是结论引用,不是欠账)。"""
    t = (ROOT/"RETRACTIONS.md").read_text(encoding="utf-8")
    # ⚠⚠ **`#999①` 被数成了「被承诺的欠账」,而它是更早一轮的**负控夹具**
    #    (原文:「`#999`① 查无此项 ⇒ 匹配器不乱开火」)。⇒ 又一次:**代理把控制当成了对象。**
    #    结构性的排除法,不靠词表:**账本里没有 `## Entry NNN` 的编号,就不是一笔承诺。**
    real = {m.group(1) for m in re.finditer(r"^## Entr(?:y|ies) (\d+)", t, re.M)}
    loose = {f"#{a}{b}" for a, b in ID.findall(t) if a in real}
    tight = {f"#{a}{b}" for a, b in ID_TIGHT.findall(t) if a in real}
    if len(loose) > len(tight):
        print(f"  ⚠ 写法对照:紧式(圈码紧跟数字)**{len(tight)}** 笔 · "
              f"宽式(允许中间有反引号)**{len(loose)}** 笔 —— "
              f"**账本主用的是宽式,紧式会漏掉 {len(loose)-len(tight)} 笔。**")
    return loose

def tabled():
    rows = (ROOT/"DEBTS.tsv").read_text(encoding="utf-8").strip().split("\n")
    return {r.split("\t")[0].strip() for r in rows if r.strip()}

def missing():
    """⚠ 排掉表头行 `debt`,以及不是真编号的东西。"""
    return sorted(x for x in (promised() - tabled()) if x.startswith("#"))
def split_by_scope(m):
    """表建于 `#830` ⇒ 之后的阻断,之前的只登记。分界从对象推出,不是我挑的。"""
    inn, out = [], []
    for x in m:
        n = int(x[1:-1])
        (inn if n >= TABLE_BORN else out).append(x)
    return inn, out

def controls():
    """★P5:先证它会开火,再证它不会对合规状态开火。"""
    p, tb = promised(), tabled()
    fake = "#9991"                                  # 一定不在账本里
    pos = bool(({fake+"①"} | p) - tb)               # 造一个「承诺了但没入表」⇒ 必须命中
    neg = not bool(p - (tb | p))                    # 全部承诺都入表的世界 ⇒ 不许命中
    one_way = not bool((tb - p) & missing_set(p, tb))  # 表里有而账本没提 ⇒ 不许算缺失
    return dict(pos=pos, neg=neg, one_way=one_way, ok=bool(pos and neg and one_way))

def missing_set(p, tb): return p - tb

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--precommit", action="store_true")
    a = ap.parse_args()
    c = controls()
    print(f"仪器控制:正控(造一个承诺但未入表的 id ⇒ 必须命中)**{c['pos']}** · "
          f"负控(全部承诺都入表 ⇒ 不许命中)**{c['neg']}** · "
          f"**单向性(表里有而账本没提 ⇒ 不算缺失)**{c['one_way']}**")
    if not c["ok"]:
        print("⛔ **控制没过 ⇒ 本次扫描不可采**(★P5)"); sys.exit(2)
    m = missing()
    extra = sorted(tabled() - promised())
    print(f"账本承诺过的欠账 **{len(promised())}** 笔 · 表里 **{len(tabled())}** 行")
    print(f"**账本提过而表里没有的:{len(m)}** ⇒ {m or '无'}")
    print(f"(表里有而账本没提的 {len(extra)} 笔 —— ⚠ **不算缺失,单向**:"
          f"欠账可以先建行后写账){' ⇒ ' + str(extra[:6]) if extra else ''}")
    print("⚠ **只报缺失,从不报「表是完整的」** —— 有行不证明那笔欠账被认真对待过。")
    inn, out = split_by_scope(m)
    print(f"  ⇒ **表的作用域内(`#{TABLE_BORN}` 及之后)缺行:{len(inn)}** ⇒ {inn or '无'}")
    print(f"  ⇒ 作用域外(`#{TABLE_BORN}` 之前)缺行:**{len(out)}** —— "
          f"⚠ **不是「{len(out)} 个没还的债」,是未分辨的历史存量**"
          f"(`#800` 已登记它们状态真的不知道,`#829` 证过从散文反推不可行)")
    if a.precommit and inn:
        print(f"🔒 PRE-COMMIT BLOCK:{len(inn)} 笔**表建立之后**被承诺的欠账没有行。"
              f"**零容忍,不是棘轮**:表存在就是为了它们,而修它只需一行 `printf`。")
        sys.exit(1)
    sys.exit(0)
