"""#807 · E03·A57·R246 —— 「不是没跑,是跑不了」;而页面还没这么说

`#806` 的 ① 与 ②,两条更正,**都必须走到页面上,而不是停在账本里。**
⚠ **本轮标注 Production —— 诚实标注,它不产生新的关于人的判断。**
但它不是可有可无的:**一条只写在账本里的更正,对读页面的人等于没发生。**

**两个总体,一个仪器:**
**① 顶对**:`#806` 发现 `sexeduc`(42%/31%)与 `racmar`(72%/71%)的两条定标臂都超过 30% 可达幅度
   ⇒ **这套共同位移对照对它们零个可用格,而且是结构性的。**
   `#801` 当时写的是「那套对照**从来没跑到**这一对上」——**语气是欠账;而真相是够不着。**
   ⚠ **两者对读者的意思完全相反**:欠账意味着「以后会补」,够不着意味着「这条路没有」。
   ⇒ 页面上每一处**断言顶对**的地方,都缺这个标注。
**② `#805` 的一般化**:页面上写了「即使两群人心里改变的幅度一模一样,他们在问卷上的答案也会越离越远」,
   而 `#806` 实测 `spanking`(+0.008~+0.036)与 `helpblk`(−0.011)上**共同位移几乎解释不了任何一部分**
   ⇒ **那句话只在 `homosex` 上成立,页面把它写成了一般规律。**

G1 估计量:**两个总体各自「仍然裸着」的处所数**(必须为 0)。

识别(沿用 `#793`→`#798` 已经付过代价的方法):
   · **总体是结构定的** —— 带 `(Entry N)` 的表行 + 相邻行内锚 `[#NNN「…」]` 之间的叙事段。
     **召回按构造 = 1**,不存在 `#792` 那个「词表在最需要的那一类上召回为零」的问题。
   · **词面扫描只用来定位,判定权不交给它**(`#755`/`#792` 的分界)。

⚠⚠ 判据的分支写在动手之前(`#798` 的教训:漏掉分支会把好工作误伤成错误):
   **ASSERT** —— 这一处**拿顶对当结论在用**(报它的 `r`、它的水平、它的差距)⇒ 需要标注。
   **LIMIT** —— 这一处本来就在说这一对的局限(`#801`/`#806` 自己)⇒ **不需要再标,它就是那条标注。**
   **MENTION** —— 只是提到题名(如列举八题)⇒ 不需要标注。

预注册判词(条件式):
  if 枚举完整(每个处所都落在结构定的总体里,没有落在总体之外的)
     and 正控开火(**一个我确知在断言顶对的处所,必须被判成 ASSERT**)
     and 负控开火(**`#806` 自己那一行必须被判成 LIMIT,不能被算成待修**):
      两个总体裸着的处所都为 0 -> 更正落地
      > 0                     -> 没落地,逐处列出
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`#806` 与 `#801` 自己的行里全是「顶对」「sexeduc」「racmar」这些词** ——
  若不区分「使用它」与「限制它」,验证会把更正本身报成待修的残留(**`#798` 那一轮踩的正是这个**)。
  ⇒ 控制:**`#801` 及其后的处所按 LIMIT 处理并单独计数**,不混进 ASSERT。

⚠ 本轮对象是页面自己 ⇒ 换不了仪器(与 `#787`/`#792`/`#793`/`#798` 同一种:对象不是世界)。
"""
import re, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
ANC = re.compile(r"\[#(\d+)「[^」]*」\]")
TOP = re.compile(r"sexeduc|racmar|顶对|top pair")
GEN = re.compile(r"即使两群人心里改变的幅度一模一样|even if the two groups had changed by exactly the same amount")
# 已落地的标注长什么样(两版各自的说法)
MARK_TOP = re.compile(r"够不着|结构性不可能|跑不了|out of this design's reach|structurally impossible|cannot run")
MARK_GEN = re.compile(r"缩回? *`?homosex|narrowed (back )?to `?homosex|只在 `?homosex|a property of `?homosex")

def units(t):
    out = []
    for line in t.split("\n"):
        if line.startswith("|") and re.search(r"\(Entry\s+\d", line):
            out.append(("row", int(re.search(r"\(Entry\s+(\d+)", line).group(1)), line))
    ancs = list(ANC.finditer(t))
    for i, m in enumerate(ancs):
        out.append(("para", int(m.group(1)), t[(ancs[i-1].end() if i else 0):m.end()]))
    return out

print("=== ① 结构定总体(带 `(Entry N)` 的表行 + 相邻行内锚之间的叙事段)===")
U = {}
for f in ("README_zh.md", "README.md"):
    t = (ROOT/f).read_text(encoding="utf-8")
    U[f] = units(t)
    rows = sum(1 for k, _, _ in U[f] if k == "row"); paras = sum(1 for k, _, _ in U[f] if k == "para")
    print(f"  {f}: 表行 **{rows}** · 叙事段 **{paras}** ⇒ 总体 **{len(U[f])}** 个单位")

def classify(owner, body, pat, mark):
    if not pat.search(body): return None
    if owner >= 801: return "LIMIT"          # `#801` 及其后 = 限制侧(跑前写下的混淆的控制)
    return "MARKED" if mark.search(body) else "ASSERT"

print("\n=== ② 两个总体,逐处判(⚠ 跑前混淆的控制:`#801` 及其后按 LIMIT,单独计数)===")
RES = {}
for f, us in U.items():
    RES[f] = {}
    for tag, pat, mark in (("顶对", TOP, MARK_TOP), ("`#805` 的一般化", GEN, MARK_GEN)):
        hits = [(k, o, classify(o, b, pat, mark), b) for k, o, b in us]
        hits = [h for h in hits if h[2]]
        naked = [h for h in hits if h[2] == "ASSERT"]
        RES[f][tag] = dict(total=len(hits), limit=sum(1 for h in hits if h[2] == "LIMIT"),
                           marked=sum(1 for h in hits if h[2] == "MARKED"), naked=len(naked),
                           naked_owners=sorted({h[1] for h in naked}))
        r = RES[f][tag]
        print(f"  {f} / {tag}: 命中 {r['total']} · 限制侧 {r['limit']} · 已标注 {r['marked']} · "
              f"**仍裸着 {r['naked']}**" + (f" ⇒ 条目 {r['naked_owners']}" if r['naked'] else ""))

tot = {tag: sum(RES[f][tag]["naked"] for f in RES) for tag in ("顶对", "`#805` 的一般化")}
print(f"\n  ⇒ **两版合计仍裸着:顶对 {tot['顶对']} 处 · `#805` 的一般化 {tot['`#805` 的一般化']} 处**(预注册要求各为 0)")

# ── 控制 ──────────────────────────────────────────────────────────────────
print("\n=== ③ 控制 ===")
zh = (ROOT/"README_zh.md").read_text(encoding="utf-8")
pc_body = next((b for k, o, b in U["README_zh.md"] if o == 801 and TOP.search(b)), "")
pc_syn = "本页把 `sexeduc` 的比值 2.192 与 `racmar` 的 1.841 当作结论在用。"
pc_ok = classify(700, pc_syn, TOP, MARK_TOP) == "ASSERT"
nc_ok = classify(806, "顶对这个设计够不着,结构性不可能", TOP, MARK_TOP) == "LIMIT"
print(f"  正控:一个**确知在断言顶对**的合成处所(条目 700,无标注)⇒ 判为 "
      f"**{classify(700, pc_syn, TOP, MARK_TOP)}**(该是 ASSERT)")
print(f"  负控:`#806` 自己那种限制侧文字 ⇒ 判为 "
      f"**{classify(806, '顶对这个设计够不着,结构性不可能', TOP, MARK_TOP)}**(该是 LIMIT,不能算待修)")
mk_ok = classify(700, pc_syn + " ⚠ 而这个设计对它够不着。", TOP, MARK_TOP) == "MARKED"
print(f"  正控②:同一处**加上标注之后**必须变成 **MARKED** ⇒ "
      f"**{classify(700, pc_syn + ' ⚠ 而这个设计对它够不着。', TOP, MARK_TOP)}**(该是 MARKED)"
      f" —— ⚠ **这一条才让「标注」这个动作可被检验:没有它,判据无法区分「改过」与「没改」**")

G = Gate("#807 · 「不是没跑,是跑不了」;而页面还没这么说")
G.asserted("① 正控:一个确知在断言顶对、且无标注的处所必须判为 ASSERT", pc_ok,
           "合成处所(条目 700)判为 ASSERT", kind="control")
G.asserted("② 负控:`#801` 及其后的限制侧文字必须判为 LIMIT,不得算成待修"
           "(⚠ 跑前写下的混淆的控制 —— `#798` 正是踩在这里)", nc_ok,
           "`#806` 式限制文字判为 LIMIT", kind="control")
G.asserted("③ 正控②:同一处加上标注之后必须变成 MARKED"
           "(否则判据分不清「改过」与「没改」,整轮的计数就没有意义)", mk_ok,
           "加标注后判为 MARKED", kind="control")
G.asserted("④ 前提:总体由结构定(表行 + 锚间叙事段),不由词表定 —— 召回按构造 = 1",
           bool(all(len(u) > 0 for u in U.values())),
           f"两版总体 {[len(u) for u in U.values()]} 个单位", kind="control")
G.asserted("⑤ kill(预注册):两个总体裸着的处所都必须为 0",
           bool(tot["顶对"] == 0 and tot["`#805` 的一般化"] == 0),
           f"顶对 {tot['顶对']} · 一般化 {tot['`#805` 的一般化']}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*94)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 计数不可信。**"
elif tot["顶对"] == 0 and tot["`#805` 的一般化"] == 0:
    V = "**两条更正都已落地:两个总体裸着的处所各为 0。**"
else:
    V = (f"**没落地。仍需标注:顶对 {tot['顶对']} 处 · `#805` 的一般化 {tot['`#805` 的一般化']} 处。**\n"
         f"  ⚠ **而「顶对」那一类的意思要改的是语气,不只是补一句**:\n"
         f"  `#801` 写的是「从来没跑到这一对上」(=欠账,读者会以为以后会补),\n"
         f"  `#806` 实测的是「跑不了」(=这条路没有)。**两者对读者的意思相反。**")
print(V)
json.dump(dict(per_file=RES, naked_total=tot, action="Production",
               pos_control=pc_ok, neg_control=nc_ok, marked_control=mk_ok,
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"two_corrections.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'two_corrections.json'}")
