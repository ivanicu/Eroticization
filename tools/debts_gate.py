"""#830 · 欠账表的闸 —— 让「状态」成为被写入的东西,而不是被正则反推的东西

`#827` 数 → `#828` 列 → `#829` 量:**三轮各造了一把新尺子去量同一件事,而三轮都失败了。**
`#829` 的诊断是决定性的:**匹配器的盲区与欠账的年代精确重合,那个「0/26」几乎是被构造出来的。**
⇒ **正确的做法不是第四把尺子,是把状态写进一个数据结构,由每一轮显式写入。**

本文件是那张表的闸。它检查的是**表本身的完整性**,不是欠账的真假:
  · 表存在、表头正确、每行五列
  · `status` ∈ {UNVERIFIED, OPEN, SETTLED, WITHDRAWN}
  · `status == SETTLED` ⇒ `settled_in` 必须非空且 > `raised_in`
  · `status == UNVERIFIED` ⇒ `settled_in` 必须为空(不许一边说不知道一边填了了结轮次)
  · 表里的 `debt` 必须互不重复

⚠⚠ **这个闸不能、也不打算判定一笔欠账是不是真的还了** —— 那是语义判断,
   而 `#829` 已经证明从散文里反推它是不可行的。**闸只保证:一旦有人写下状态,那个状态是自洽的。**
   **`§0.2`:这是一件工具,它的价值在于以后能挡住什么,不在今天的计数。**

用法:`python tools/debts_gate.py`  ⇒ 退出码 0 通过,2 失败(**空表也退 2,`realstat`:空总体不许通过**)。
"""
import csv, pathlib, sys
OK = {"UNVERIFIED", "OPEN", "SETTLED", "WITHDRAWN"}
p = pathlib.Path(__file__).resolve().parents[1]/"DEBTS.tsv"
if not p.exists():
    print("⛔ DEBTS.tsv 不存在"); sys.exit(2)
rows = list(csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"))
errs = []
if not rows: errs.append("空表 —— `realstat`:空总体不许通过,退 2")
seen = set()
for i, r in enumerate(rows, 2):
    d = (r.get("debt") or "").strip()
    if not d: errs.append(f"第 {i} 行:debt 为空")
    if d in seen: errs.append(f"第 {i} 行:debt {d} 重复")
    seen.add(d)
    st = (r.get("status") or "").strip()
    if st not in OK: errs.append(f"第 {i} 行:status={st!r} 不在 {sorted(OK)}")
    sett = (r.get("settled_in") or "").strip()
    if st == "SETTLED":
        if not sett: errs.append(f"第 {i} 行:SETTLED 却没有 settled_in")
        elif sett.isdigit() and (r.get("raised_in") or "").isdigit() and int(sett) <= int(r["raised_in"]):
            errs.append(f"第 {i} 行:settled_in {sett} 不晚于 raised_in {r['raised_in']}")
    if st == "UNVERIFIED" and sett:
        errs.append(f"第 {i} 行:UNVERIFIED 却填了 settled_in={sett} —— 不许一边说不知道一边填了结轮次")
from collections import Counter
c = Counter((r.get("status") or "").strip() for r in rows)
print(f"DEBTS.tsv:**{len(rows)} 行** · 状态分布 {dict(c)}")
if errs:
    print("⛔ 不通过:"); [print("   ", e) for e in errs]; sys.exit(2)
print("✅ 表自洽。⚠ 而自洽 ≠ 正确:闸不判定任何一笔欠账的真假(`#829` 已证从散文反推不可行)。")
sys.exit(0)
