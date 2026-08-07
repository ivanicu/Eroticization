"""#833① —— `units()` 把表行切了两次,而两次的边界不一样。

**缺陷**(`#833` 记下,`#841` 之后才还):
  旧版把页面切成两类单位:
    `row`  —— 每个带 `(Entry N)` 的表行,取**整行**;
    `para` —— 相邻两个行内锚 `[#NNN「…」]` 之间的文本,取到**锚为止**。
  **而表行自己也带行内锚** ⇒ 表行的文本**同时进了 `row` 和 `para`**,
  **且两次的右边界不同**:`row` 到行尾(含锚之后的 `⟨…⟩` 注记块),`para` 只到锚。
  ⇒ 两个后果,方向相反,都要说:
    ① **总体规模虚高**(同一段文本被数两次)—— 自 `#798` 起的每个「总体 N」都受影响;
    ② **锚后追加的注记被算到下一条账上** ——
       ⚠⚠ **而这一条比 `#833` 当初写的更糟,是本轮的控制逼出来的更正:**
       `#833` 写的是「锚后追加的注记对 `para` 不可见」。**实测不是不可见,是归错了主。**
       旧版 `para` 的左边界是**上一个锚的结尾**,所以 `#N` 锚之后的文本落进了 **`#N+1` 的单位**。
       实测(见 `controls()`):`⟨注记:三堆⟩` 写在 `#100` 的锚后,旧版把它交给了 **`#101`**。
       ⇒ 后果不是「少数了」,而是**「数到别人头上」** —— 对 `#798` 这种
       **按 owner 把单位分成「撤回侧 / 使用侧」**的轮次,一次错归就可能把单位换边。
       **一个「不可见」的缺陷只让你漏;一个「错归」的缺陷会让你反向定罪。**
       而本项目的更正**恰恰**追加在锚之后(`#836`/`#837` 的就地收窄就是这么写的)。

**修法**:先把表行整行从文本里摘出去,`para` 只在**剩下的文本**上找锚。
  ⇒ `row` 与 `para` 的文本**不再重叠**,且 `row` 保留整行(含锚后注记)。

⚠ **三份拷贝**(`R237`/`R246`/`R267`)实测**语义相同**(格式不同,逻辑一致)——
  所以三轮用的是同一个总体定义,缺陷是**统一的**,不是三种。
  ⇒ 修在一处(本文件),三轮的重算都指向它。`feedback_fix_lands_on_one_path`:
  **命名不变量,枚举调用者** —— 调用者就是这三个,已枚举。
"""
import re

ANC = re.compile(r"\[#(\d+)「[^」]*」\]")
ROW = re.compile(r"^\|.*\(Entry\s+(\d+)")

def units(t, fixed=True):
    """把页面切成互不重叠的判定单位。`fixed=False` 复现旧版,供逐轮对照。"""
    out = []
    keep = []
    for line in t.split("\n"):
        m = ROW.match(line)
        if m:
            out.append(("row", int(m.group(1)), line))
            if fixed: continue          # ⇐ 修好的地方:表行不再流进 para
        keep.append(line)
    src = "\n".join(keep) if fixed else t
    a = list(ANC.finditer(src))
    for i, m in enumerate(a):
        out.append(("para", int(m.group(1)), src[(a[i-1].end() if i else 0):m.end()]))
    return out

def controls():
    """★P5:先证这把尺子会开火,再证它不会对合规写法开火。"""
    page = ("| x `[#100「a」]`(Entry 100) ⟨注记:三堆⟩|\n"
            "\n正文一段 `[#101「b」]`\n")
    old, new = units(page, fixed=False), units(page, fixed=True)
    # 正控:锚后注记「三堆」在新版只出现在 row 里;**在旧版被归给了下一个锚**
    old_owner = {n for k, n, b in old if k == "para" and "三堆" in b}
    new_para_sees = any(k == "para" and "三堆" in b for k, _, b in new)
    row_sees = any(k == "row" and "三堆" in b for k, _, b in new)
    misattributed = old_owner == {101}          # ⇐ 写在 #100 锚后,旧版却记到 #101 名下
    pc = row_sees and not new_para_sees and misattributed
    # 正控 β:旧版的 para 必须**含有表行文本**(这就是重复计数),新版必须不含
    old_dup = any(k == "para" and "(Entry 100)" in b for k, _, b in old)
    new_dup = any(k == "para" and "(Entry 100)" in b for k, _, b in new)
    # 负控:没有表行的页面上,两版必须给出**相同**的 para
    plain = "甲 `[#200「p」]`\n乙 `[#201「q」]`\n"
    nc = ([b for k, _, b in units(plain, False) if k == "para"]
          == [b for k, _, b in units(plain, True) if k == "para"])
    return dict(anchor_note_misattributed_by_old=pc, old_misattributed_to=sorted(old_owner),
                old_para_duplicated_row=old_dup,
                new_para_duplicated_row=new_dup, no_table_identical=nc,
                ok=bool(pc and old_dup and not new_dup and nc))

if __name__ == "__main__":
    c = controls()
    for k, v in c.items(): print(f"  {k:28s} {v}")
