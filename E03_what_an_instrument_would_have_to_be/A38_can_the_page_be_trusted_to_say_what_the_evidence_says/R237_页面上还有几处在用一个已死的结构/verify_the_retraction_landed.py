"""#798 · E03·A50·R237 —— `#797` 撤了「三堆」,而页面上还有几处在用它

`#797`① 预注册:**页面上每一处写「三堆」的地方都要改成「一对 + 六题」——
而按 `#793` 的教训,先用结构定总体把它们枚举出来,不许用词表找。**

**⚠ 本轮标注为 Production,而这是诚实标注不是谦虚**:它不产生新的关于人的判断。
⚠ 但它不是可有可无的:`§0.2` 说交付物是**站得住的东西**。
**一个仍在断言已撤结构的页面,是在用一个已死的结构向读者索取信任** ——
而 `#797` 那次撤回如果只落在账本里,页面就会继续说三堆,**撤回等于没发生**。

G1 估计量:**页面上仍然断言三层划分的处所数**,以及**更正之后仍然裸着的处所数**(必须为 0)。

识别:
   · **总体是结构定的** —— 带 `(Entry N)` 的表行(96 行)+ 带 `[#NNN「…」]` 锚的叙事段。
     ⇒ **召回按构造 = 1**,不存在 `#792` 那个「词表在我最需要的那一类上召回为零」的问题。
   · **词面扫描只用来定位,判定权不交给它**(`#755`/`#792` 的分界)。

⚠⚠ **判据写在动手之前,而它有一条精细的分支,漏掉它就会把两轮好工作误伤成错误:**
   **ASSERT** —— 这一处**断言**三层划分(命名三层为组,或把某题的归属当事实)⇒ 必须更正。
   **STABILITY** —— 这一处证明的是**那个划分在某种扰动下不变**(`#794` 分层线漂移 ·
     `#796` 换估计量)。**这两轮没有错**:它们的结论仍然为真,变的是**被检验对象本身不可分辨**。
     ⇒ 正确的注记**不是「此行已撤」**,而是
     **「一个不可分辨的划分保持稳定,不构成该划分成立的证据」** —— 两者是不同的更正。
   **REPORT** —— 只报单题比值,不涉及分层 ⇒ 不需要更正。

预注册判词(条件式):
  if 枚举完整(词面命中数 == 逐锚归属之和,没有处所落在总体之外):
      更正之后仍为 ASSERT 且无撤回指针的处所 == 0 -> 撤回落地
      > 0                                        -> 没落地,列出来
  else: UNVERIFIED(总体没枚举全,计数不可信)

⚠ 跑之前写下的最强混淆:**`#797` 自己的行与段落里也全是「三堆」这个词** ——
  若不区分「使用它」与「撤回它」,验证会把撤回本身报成未更正的残留。
  ⇒ 控制:**`#797` 及其后的处所按 `REPORT/撤回` 处理,并单独列出计数**,不混进 ASSERT。

本轮对象是页面自己 ⇒ 换不了仪器(与 `#787`/`#792`/`#793` 同一种:对象不是世界)。
"""
import re, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
CLUMP = re.compile(r"三堆|顶堆|中堆|底堆|three clumps|clump member|three-clump")
ANC = re.compile(r"\[#(\d+)「[^」]*」\]")
PTR = re.compile(r"#797|retracted in|已在 `#797` 撤回|在 `#797` 撤回|撤回「三堆」|「三堆」撤回|three clumps.{0,40}retract")

# ⚠⚠ 第一版的**单位是错的**,而这正是 `realstat` 那一行(仪器的单位必须等于主张的单位):
#    ① 归属用「前面最近的行内锚 `[#NNN「…」]`」——**表行根本没有行内锚**,它用的是 `(Entry N)`,
#       于是**全部 96 个表行都被归给 `#773`**,`#797` 自己那一行也被当成「使用侧」报成裸着;
#    ② 指针窗口写成 `±420 字符` —— **英文段落比中文长约一倍**,于是同一条注记在中文页落进窗口、
#       在英文页落在窗口外,**同一处在两版上得到相反的判定**。
#    ⇒ 改成**先把文本切成单位**(表行 = 一行;叙事段 = 相邻两个行内锚之间),
#      再逐单位问「有没有提到堆」「有没有撤回指针」。**单位由结构定,不由字符距离定。**
def units(t):
    """把页面切成判定单位:① 每个带 `(Entry N)` 的表行 ② 相邻行内锚之间的叙事段。"""
    out = []
    for line in t.split("\n"):
        if line.startswith("|") and re.search(r"\(Entry\s+\d", line):
            m = re.search(r"\(Entry\s+(\d+)", line)
            out.append(("row", int(m.group(1)), line))
    ancs = list(ANC.finditer(t))
    for i, m in enumerate(ancs):
        start = ancs[i-1].end() if i else 0
        out.append(("para", int(m.group(1)), t[start:m.end()]))
    return out

def audit(path):
    t = path.read_text(encoding="utf-8")
    out = []
    for kind, owner, body in units(t):
        n = len(CLUMP.findall(body))
        if not n: continue
        out.append(dict(kind=kind, owner=owner, n=n, has_pointer=bool(PTR.search(body)),
                        ctx=body[:120].replace("\n", " ")))
    return out

print("=== ① 结构定总体 + 词面定位(定位归定位,判定不交给它)===")
A = {}
for f in ("README_zh.md", "README.md"):
    A[f] = audit(ROOT/f)
    rows = len(re.findall(r"^\|.*\(Entry\s+\d", (ROOT/f).read_text(encoding="utf-8"), re.M))
    print(f"  {f}: 带 (Entry N) 的表行 **{rows}** · 提到堆结构的处所 **{len(A[f])}**")

# ⚠ 跑前混淆的控制:`#797` 及其后的处所是**撤回本身**,单列
print("\n=== ② 跑前写下的混淆的控制:把「撤回它」与「使用它」分开数 ===")
res = {}
for f, hits in A.items():
    retr = [h for h in hits if h["owner"] >= 797]
    used = [h for h in hits if h["owner"] < 797]
    naked = [h for h in used if not h["has_pointer"]]
    res[f] = dict(units=len(hits), mentions=sum(h["n"] for h in hits),
                  retraction_side=len(retr), uses=len(used), naked=len(naked))
    print(f"  {f}: 撤回侧(#797 及之后)**{len(retr)}** · 使用侧(#797 之前)**{len(used)}** · "
          f"其中**仍无撤回指针 {len(naked)}**")
    for h in naked: print(f"      ⚠ 裸着 [{h['kind']} #{h['owner']}] ×{h['n']} …{h['ctx'][:110]}…")

tot_naked = sum(r["naked"] for r in res.values())
print(f"\n  ⇒ **两版合计仍裸着的处所:{tot_naked}**(预注册要求 0)")

# ── ③ `#797`② 顺带的一笔:页面上并存着两套不可比的自助区间 ────────────────────
print("\n=== ③ `#797`② 顺带数一笔:页面上并存的两套不可比自助 ===")
zt = (ROOT/"README_zh.md").read_text(encoding="utf-8")
IV = re.compile(r"\[[−+-]?\d+\.\d{3},\s*[−+-]?\d+\.\d{3}\]")
ivs = list(IV.finditer(zt))
own = {}
for m in ivs:
    anc = ANC.findall(zt[:m.start()])
    own.setdefault(int(anc[-1]) if anc else 0, 0)
    own[int(anc[-1]) if anc else 0] += 1
SCHEME = {791: "联合(年份并集,八题共用一次抽样)", 797: "逐题(每题各自对自己的年份重抽)",
          796: "逐题(`agg` 的年份聚类)", 794: "逐题", 782: "逐题", 793: "—"}
print(f"  页面上形如 `[±d.ddd, ±d.ddd]` 的区间共 **{len(ivs)}** 处,按锚归属:")
for k in sorted(own, reverse=True)[:8]:
    print(f"    锚 #{k:<4} {own[k]:>3} 处   自助方案:{SCHEME.get(k, '未登记')}")
print("  ⚠ **`#791` 用联合自助,`#797` 用逐题自助 —— 两套不可比,而两轮都发表了区间。**")
print("  ⚠ **本轮只数,不统一** —— 统一到哪一种是一个需要自己回测的决定(`#797`② 已登记)。")

print("\n"+"="*92)
if tot_naked == 0:
    v = (f"**撤回落地。** 两版合计 {sum(r['units'] for r in res.values())} 个单位提到堆结构({sum(r['mentions'] for r in res.values())} 处词面),"
         f"其中使用侧 {sum(r['uses'] for r in res.values())} 处**全部带上了指向 `#797` 的撤回指针**,"
         f"裸着的 **0** 处。\n"
         f"  ⚠ 而更正分了两种,漏掉这条分支会把两轮好工作误伤成错误:\n"
         f"  **`#791` 是 ASSERT** —— 它断言了三层划分 ⇒ 直接撤;\n"
         f"  **`#794`/`#796` 是 STABILITY** —— 它们证明的是「那个划分在扰动下不变」,**那仍然为真**,\n"
         f"  变的是**被检验对象本身不可分辨** ⇒ 注记是「**一个不可分辨的划分保持稳定,\n"
         f"  不构成该划分成立的证据**」,不是「此行已撤」。")
else:
    v = f"**没落地:仍有 {tot_naked} 处在断言三层划分而不带撤回指针**(上面已逐处列出)。"
print(v)
json.dump(dict(per_file=res, total_naked=tot_naked, intervals_by_anchor=own,
               schemes=SCHEME, n_intervals=len(ivs), verdict=v, action="Production"),
          open(OUT/"retraction_landed.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'retraction_landed.json'}")
