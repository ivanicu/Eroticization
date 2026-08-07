"""#828 · E03·A72·R267 —— 把两条更正落到页面上;而第二条的仪器单位与主张单位不相等

`#827` 的 NEXT,两件:
① `#826` 校正之后 **`sexeduc`/1990s 与 `teensex`/2010s 都不存活**,
   **而它们此刻还写在页面的网格里当作「偏离」。**
② **重数欠账,一笔一笔,并写成一个可以被 `grep` 检验的清单**,
   而不是又一个抄在收尾行的数字(`#827` 测出那个数字被抄了 25 遍从没重数)。

⚠ **本轮标注 Production —— 诚实标注,不是谦虚**:它没有 ≥2 个本体世界,
  所以**不填 Cognitive Update Card**(`frontier §7` 规则 6)。它产生的是**落地**与**清单**,不是新判断。
  ⚠ 但它仍然报数,所以 `realstat` 照常适用,控制照常要跑。

**总体 ①(结构定,`#793`→`#798`→`#807` 已付过学费的方法)**:
   带 `(Entry N)` 的表行 + 相邻行内锚之间的叙事段。**召回按构造 = 1。**
   判定三分支:**ASSERT**(把那两格当作「偏离」在用,且没有 `#826` 的更正注记)⇒ 必须改;
   **MARKED**(已带 `#826` 更正)⇒ 不需要改;**LIMIT**(`#826` 及其后,即更正本身)⇒ 不算待修。
   ⚠ **最后这一分支是 `#798` 踩过的坑**:一次更正必然把被更正的措辞写进页面,
   **不分侧的话,验证会把更正本身报成残留。**

**总体 ②:`#800` 那条规范清单的 26 笔欠账,逐笔查后文有没有声称执行过。**

⚠⚠ **而这里必须先写下一件会决定结论措辞的事(`realstat`「搜索是仪器且没有正对照」那一行的补丁):**
   **仪器的单位 = 「账本里出现一句声称 `#NNN`① 被执行」;**
   **主张的单位 = 「那笔欠账真的被了结」。**
   **这两个单位不相等 —— 一句「本轮是 `#NNN`① 的执行」不等于那件事做完了。**
   ⇒ **所以本轮只能产出「声称已执行」的计数,不能产出「已了结」的计数;**
   **真实了结数登记为 UNVERIFIED,而不是四舍五入成前者。**

预注册判词(条件式):
  if 正控开火(**`#802`① 是已知被还的**(`#808` 明写「`#802`① 已还」)⇒ 匹配器必须找到它)
     and 负控开火(**一个不存在的欠账号 `#999`① 必须查无此项** ⇒ 匹配器不乱开火)
     and 正控②开火(**同一处加上更正注记之后,判定必须从 ASSERT 变成 MARKED** ——
        `#807` 立的第三条:**没有它,「裸着 0」与「一个字没改」在判据眼里一样**):
      总体① 裸着的处所 == 0 -> 落地
      并逐笔输出总体② 的清单(声称已执行 / 未见声称),**写成可 grep 的行**
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`#826`/`#827` 自己的行里全是 `sexeduc`/1990s 与 `teensex`/2010s** ——
  不分侧就会把更正本身报成待修(`#798` 原样的坑)。⇒ 控制:**拥有者条目号 ≥ 826 一律算 LIMIT。**

⚠ 本轮**换不了仪器**:对象是页面与账本自己。⚠ 总判由 `Gate.admissible()` 决定。
"""
import re, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
ANC = re.compile(r"\[#(\d+)「[^」]*」\]")
# 那两格的写法:`sexeduc`/1990s · `teensex`/2010s(页面上以多种排版出现)
TARGET = re.compile(r"`?sexeduc`?\s*[/／]\s*1990s|`?teensex`?\s*[/／]\s*2010s|"
                    r"`sexeduc` 1990s|`teensex` 2010s")
MARKED = re.compile(r"#826|校正后(都)?不存活|不存活|did not survive|BH|多重性校正")

def units(t):
    out = []
    for line in t.split("\n"):
        if line.startswith("|") and re.search(r"\(Entry\s+\d", line):
            out.append(("row", int(re.search(r"\(Entry\s+(\d+)", line).group(1)), line))
    a = list(ANC.finditer(t))
    for i, m in enumerate(a):
        out.append(("para", int(m.group(1)), t[(a[i-1].end() if i else 0):m.end()]))
    return out

def classify(owner, body):
    if not TARGET.search(body): return None
    if owner >= 826: return "LIMIT"          # 跑前混淆的控制:更正本身不算待修
    return "MARKED" if MARKED.search(body) else "ASSERT"

print("=== ① 总体一:页面上还把那两格当作「偏离」的处所(结构定总体)===")
POP = {}
for f in ("README_zh.md", "README.md"):
    t = (ROOT/f).read_text(encoding="utf-8"); us = units(t)
    hits = [(k, o, classify(o, b), b) for k, o, b in us]
    hits = [h for h in hits if h[2]]
    naked = [h for h in hits if h[2] == "ASSERT"]
    POP[f] = dict(units=len(us), hits=len(hits), naked=len(naked),
                  limit=sum(1 for h in hits if h[2] == "LIMIT"),
                  marked=sum(1 for h in hits if h[2] == "MARKED"),
                  owners=sorted({h[1] for h in naked}))
    r = POP[f]
    print(f"  {f}: 总体 {r['units']} 单位 · 命中那两格 {r['hits']} · 限制侧 {r['limit']} · "
          f"已标注 {r['marked']} · **仍裸着 {r['naked']}**" + (f" ⇒ 条目 {r['owners']}" if r['naked'] else ""))
    for h in naked: print(f"      ⚠ 裸着 [{h[0]} #{h[1]}] …{h[3][:100].replace(chr(10),' ')}…")
tot_naked = sum(POP[f]["naked"] for f in POP)
print(f"  ⇒ **两版合计仍裸着 {tot_naked} 处**(预注册要求 0)")

print("\n=== ② 总体二:`#800` 的 26 笔欠账,逐笔查后文有没有**声称**执行过 ===")
led = (ROOT/"RETRACTIONS.md").read_text(encoding="utf-8")
marks = [(int(m.group(1)), m.start()) for m in re.finditer(r'^## Entry (\d+)', led, re.M)]
bodies = {n: led[s:(marks[i+1][1] if i+1 < len(marks) else len(led))] for i, (n, s) in enumerate(marks)}
best = max(bodies.items(), key=lambda kv: len(re.findall(r'`#(\d+)`([①②③])', kv[1])))
DEBTS = sorted(set(re.findall(r'`#(\d+)`([①②③])', best[1])), key=lambda x: (-int(x[0]), x[1]))
print(f"  规范清单来自条目 **#{best[0]}**,去重后 **{len(DEBTS)}** 笔")

CLAIM = lambda a, b: re.compile(rf'`#{a}`{b}[^。\n]{{0,50}}?(?:的执行|的直接执行|已还|还上|本轮就是|一起做|在此撤回|作废)')
def claimed(a, b, exclude_owner=None):
    for n, body in bodies.items():
        if exclude_owner and n <= int(a): continue
        m = CLAIM(a, b).search(body)
        if m and n > int(a): return n, m.group(0)[:44]
    return None, None
LIST = []
for a, b in DEBTS:
    n, txt = claimed(a, b, exclude_owner=True)
    LIST.append(dict(debt=f"#{a}{b}", claimed_by=n, evidence=txt))
n_claim = sum(1 for x in LIST if x["claimed_by"])
print(f"  **可 grep 的清单**(每行一笔;`DEBT` 前缀便于 `grep '^  DEBT'`):")
for x in LIST:
    st = f"声称已执行 @#{x['claimed_by']}" if x["claimed_by"] else "**未见声称**"
    print(f"  DEBT {x['debt']:<8s} {st}" + (f"  …{x['evidence']}…" if x["evidence"] else ""))
print(f"  ⇒ **{len(DEBTS)} 笔中,{n_claim} 笔在后文被声称执行过,{len(DEBTS)-n_claim} 笔未见声称**")
print(f"  ⚠⚠ **而这不是「已了结 {n_claim} 笔」**:仪器的单位是「账本里出现一句声称执行」,"
      f"主张的单位是「那笔欠账真的被了结」——**两个单位不相等**,")
print(f"     **所以真实了结数登记为 UNVERIFIED,不许四舍五入成 {n_claim}。**")

print("\n=== ③ 控制 ===")
pc1_n, pc1_t = claimed("802", "①", exclude_owner=True)
pc1 = pc1_n is not None
print(f"  正控①:`#802`① 是**已知被还的**(`#808` 明写「`#802`① 已还」)⇒ 匹配器找到:"
      f"**{pc1}**" + (f" @#{pc1_n} …{pc1_t}…" if pc1 else ""))
nc_n, _ = claimed("999", "①", exclude_owner=True)
nc = nc_n is None
print(f"  负控:一个**不存在**的欠账号 `#999`① ⇒ 查无此项:**{nc}**(匹配器不乱开火)")
syn_bare = "本页把 `sexeduc`/1990s 当作一个偏离在用。"
syn_marked = syn_bare + " ⚠ 而 `#826` 校正后它并不存活。"
c_bare, c_marked = classify(700, syn_bare), classify(700, syn_marked)
pc2 = (c_bare == "ASSERT" and c_marked == "MARKED")
print(f"  正控②(`#807` 立的第三条):合成处所加注记前 **{c_bare}** → 加注记后 **{c_marked}** ⇒ "
      f"判定确实会变:**{pc2}**")
print(f"     ⚠ **没有它,「裸着 0」与「一个字没改」在判据眼里是同一件事。**")
nc2 = classify(826, "`#826` 校正后 `sexeduc`/1990s 不存活") == "LIMIT"
print(f"  负控②(跑前写下的混淆):`#826` 及其后的更正文字必须判 LIMIT,不算待修:**{nc2}**")

G = Gate("#828 · 把两条更正落到页面上")
G.asserted("① 正控:`#802`① 是已知被还的(`#808` 明写)⇒ 匹配器必须找到它"
           "(`realstat`:一个 grep 是测量仪器,必须在答案已知处跑过)",
           bool(pc1), f"找到 @#{pc1_n}:{pc1_t}", kind="control")
G.asserted("② 负控:不存在的欠账号 `#999`① 必须查无此项 —— 匹配器不乱开火",
           bool(nc), "`#999`① 未命中", kind="control")
G.asserted("③ 正控②(`#807` 的第三条,最容易漏):同一处**加上更正注记之后判定必须从 ASSERT 变 MARKED**"
           " —— **没有它,「裸着 0」与「一个字没改」在判据眼里一样**",
           bool(pc2), f"{c_bare} → {c_marked}", kind="control")
G.asserted("④ 负控②(跑前写下的混淆):`#826` 及其后的更正文字判 LIMIT,不算待修"
           "(`#798` 原样的坑:一次更正必然把被更正的措辞写进页面)",
           bool(nc2), "`#826` 式文字 ⇒ LIMIT", kind="control")
G.asserted("⑤ 前提(`realstat` 单位相等):**仪器单位「一句声称执行」≠ 主张单位「已了结」** ⇒ "
           "本轮只产出前者的计数,**后者登记 UNVERIFIED**",
           True, f"声称已执行 {n_claim}/{len(DEBTS)} · 真实了结数 = UNVERIFIED", kind="control")
G.asserted("⑥ kill(预注册):两条更正落地要成立,需总体① 裸着的处所 == 0",
           bool(tot_naked == 0), f"裸着 {tot_naked} 处", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 计数不可信。**"
elif tot_naked == 0:
    V = (f"**两条更正都已落地。** 总体① 裸着 0 处;总体② 产出 {len(DEBTS)} 行可 grep 的清单,"
         f"其中 **{n_claim} 笔被声称执行过、{len(DEBTS)-n_claim} 笔未见声称**。\n"
         f"  ⚠⚠ **而「已了结」的真实笔数仍是 UNVERIFIED —— 仪器只能看见「声称」。**")
else:
    V = (f"**没落地:总体① 仍有 {tot_naked} 处把那两格当作偏离在用**(上面已逐处列出)。\n"
         f"  ⇒ **而这正是 `#826` 的结论没有走到页面上的证据。**")
print(V)
json.dump(dict(action="Production", pop1=POP, total_naked=tot_naked,
               debts_source_entry=best[0], debts=LIST, n_debts=len(DEBTS), n_claimed=n_claim,
               n_discharged="UNVERIFIED — instrument unit != claim unit",
               pos_control=bool(pc1), neg_control=bool(nc), pos_control2=bool(pc2), neg_control2=bool(nc2),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"two_corrections_landed.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'two_corrections_landed.json'}")
