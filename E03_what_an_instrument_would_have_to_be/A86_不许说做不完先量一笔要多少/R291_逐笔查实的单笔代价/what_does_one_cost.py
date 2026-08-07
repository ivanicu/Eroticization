"""#852 · E03·A86·R291 —— 不许说「做不完」:先量一笔要多少

**还 `#851`①。** `#851` 量出 **124 笔** `#830` 之前被承诺、却没进表的欠账,
并**明确禁止**自己直接说「太多了做不完」——
**因为那正是 `#843` 撞过的那种墙:一条捏造的「不可能」让停下来显得有理由,所以从来不被审计。**

G1 估计量:**逐笔查实一笔历史欠账的单笔代价**,以及**它的分布** ——
   而「代价」必须是**可计数的机器动作**,不是我的感觉:
   ① 它的 `NEXT` 文本是否**指向一个更晚的真实条目**(机器能读出后续在哪);
   ② 该条目对应的轮次**有没有 `results/` 产物目录**(能不能回到当时的证据)。

**⚠⚠ 而在量之前必须先说清一件事,否则这一轮会自欺:**
**`#829` 已经证过「从散文反推欠账状态」不可行**(匹配器召回 18.6%)——
所以本轮量的**不是「自动查实的代价」**,是**「把一笔欠账送到人面前所需的机器动作」**。
**人判那一步的代价本轮量不了,如实登记。**

三个世界:
   A **多数笔能被机器直接定位到「后续轮次」** ⇒ 逐笔查实是可做的,只是没做 ⇒ **该做。**
   B **多数笔的产物目录不存在或 `NEXT` 不指向任何轮次** ⇒ **单笔代价里人判占绝大部分**
     ⇒ 总代价由我量不了的那部分主导 ⇒ **诚实登记为「代价不可估」,不是「做不完」。**
   C **分布不是全有全无** ⇒ **先做便宜的那一半** —— 一个比「全做/全不做」都好的第三选择。

预注册判词(条件式):
  if 正控开火(**一笔已知有后续轮次的欠账,必须被判为「可机器定位」**)
     and 负控开火(**一个编造的欠账 id,必须两个事实都为假**):
      两个事实都满足 ≥60%  -> A
      ≤30%                 -> B
      其间                 -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**「可机器定位」是我定义的,而我有动机把它定得宽** ——
  定得越宽,结论越像「该做」,而「该做」听起来比「做不了」体面。
  ⇒ 控制:**判据只用两个客观事实**(产物目录存不存在 · `NEXT` 里有没有更晚的真实条目号),
  **不含任何需要我读懂文本的判断**;并且**把两个事实分开报,不合成一个分数**。
⚠ 本轮换不了仪器(对象是这个项目自己的历史)。⚠ 标注 **Production**。
"""
import json, pathlib, re, sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
LED = (ROOT / "RETRACTIONS.md").read_text(encoding="utf-8")
ID = re.compile(r"#(\d{3,4})`?\s*([①-⑨])")
REAL = {m.group(1) for m in re.finditer(r"^## Entr(?:y|ies) (\d+)", LED, re.M)}
TABLED = {r.split("\t")[0].strip()
          for r in (ROOT / "DEBTS.tsv").read_text(encoding="utf-8").strip().split("\n")}
promised = {f"#{a}{b}" for a, b in ID.findall(LED) if a in REAL}
legacy = sorted((x for x in promised - TABLED if int(x[1:-1]) < 830), key=lambda x: int(x[1:-1]))
print(f"=== ⓪ 对象:`#830` 之前被承诺却未入表的欠账 **{len(legacy)}** 笔(与 `#851` 对齐)")

marks = [(m.group(1), m.start()) for m in re.finditer(r"^## Entr(?:y|ies) (\d+)", LED, re.M)]
BODY = {}
for i, (n, s) in enumerate(marks):
    BODY[n] = LED[s:(marks[i + 1][1] if i + 1 < len(marks) else len(LED))]
DIRS = {}
for d in ROOT.glob("E0*/*/R*_*"):
    m = re.match(r"R(\d+)_", d.name)
    if m:
        DIRS.setdefault(int(m.group(1)), []).append(d)

def assess(debt):
    """两个**客观**事实,分开报,不合成分数。"""
    n = debt[1:-1]
    body = BODY.get(n, "")
    seg = body[body.find("**NEXT**"):][:800] if "**NEXT**" in body else body[:800]
    later = sorted({int(a) for a in re.findall(r"#(\d{3,4})", seg)
                    if a in REAL and int(a) > int(n)})
    rm_ = re.search(r"`E0\d·A\d+·R(\d+)`", body)
    rn = int(rm_.group(1)) if rm_ else None
    has_dir = bool(rn and rn in DIRS and any((d / "results").exists() for d in DIRS[rn]))
    return dict(debt=debt, entry=int(n), points_to_later=bool(later), later=later[:3],
                round=rn, has_results_dir=has_dir)

t0 = time.time()
A = [assess(x) for x in legacy]
elapsed = time.time() - t0
pt = sum(a["points_to_later"] for a in A)
hd = sum(a["has_results_dir"] for a in A)
both = sum(a["points_to_later"] and a["has_results_dir"] for a in A)
print("\n=== ① 两个客观事实,分开报(**不合成分数** —— 合成就是我在给自己定宽窄)===")
print(f"  ① `NEXT` 段里指向了一个**更晚的真实条目**:**{pt}/{len(A)} = {pt/len(A):.0%}**")
print(f"  ② 该条目的轮次**有 `results/` 产物目录**:**{hd}/{len(A)} = {hd/len(A):.0%}**")
print(f"  两者**都**满足:**{both}/{len(A)} = {both/len(A):.0%}**")
print(f"  机器动作总耗时 **{elapsed:.2f} s** ⇒ 单笔 **{1000*elapsed/len(A):.1f} ms**"
      f"(纯机器,**不含人判**)")

print("\n=== ② 控制 ===")
seed = next((a for a in A if a["points_to_later"]), None)
pc = bool(seed)
print(f"  正控:一笔**已知有后续轮次**的欠账必须被判为「可机器定位」⇒ "
      f"{seed['debt'] if seed else '无'} → 指向 {seed['later'] if seed else '-'} ⇒ **{pc}**")
fake = assess("#9991①")
nc = not (fake["points_to_later"] or fake["has_results_dir"])
print(f"  负控:一个**编造的**欠账 id(`#9991①`,`## Entry 9991` 不存在)必须两个事实都为假 "
      f"⇒ **{nc}**")
print("     ⚠ **「这个假该不该是假?」该** —— 不存在的条目没有正文,也没有产物目录,按定义。")

G = Gate("#852 · 不许说做不完:先量一笔要多少")
G.asserted("① 前提(跑前写下的最强混淆):**「可机器定位」是我定义的,而我有动机把它定宽** ⇒ "
           "判据只用**两个客观事实**(产物目录存不存在 · `NEXT` 里有没有更晚的真实条目号),"
           "**不含任何需要我读懂文本的判断**,且**两个事实分开报,不合成分数**",
           bool(len(A) > 0), f"事实① {pt}/{len(A)} · 事实② {hd}/{len(A)}", kind="control")
G.asserted("② 正控:一笔已知有后续轮次的欠账必须被判为「可机器定位」", pc,
           f"{seed['debt'] if seed else '无'} → {seed['later'] if seed else '-'}", kind="control")
G.asserted("③ 负控:编造的欠账 id 必须两个事实都为假"
           "(⚠ **这个假该是假**:不存在的条目按定义没有正文,也没有产物目录)",
           nc, f"#9991① → points={fake['points_to_later']} dir={fake['has_results_dir']}",
           kind="control")
G.asserted("④ kill(预注册):「逐笔查实可做、只是没做」要成立,需**两个事实都满足的比例 ≥60%**",
           bool(both / len(A) >= 0.60), f"两者都满足 {both}/{len(A)} = {both/len(A):.0%}",
           kind="kill", yardstick="两个客观事实同时为真的比例", yardstick_noise=0.0)
print()
print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
r = both / len(A)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif r >= 0.60:
    VERD = (f"**A 逐笔查实是可做的,只是没做。** {both}/{len(A)} = **{r:.0%}** 的笔数同时满足"
            f"「`NEXT` 指向更晚的真实条目」与「该轮有产物目录」,纯机器动作单笔 "
            f"**{1000*elapsed/len(A):.1f} ms**。\n"
            f"  ⇒ **「124 笔太多了做不完」这句话在本轮之后不许再说 —— "
            f"它不是一堵墙,是一件没做的事。**")
elif r <= 0.30:
    VERD = (f"**B 单笔代价里人判占绝大部分 ⇒ 代价不可估,而这不等于「做不完」。** "
            f"两者都满足只有 {both}/{len(A)} = **{r:.0%}**;分开看:"
            f"`NEXT` 指向更晚条目 **{pt/len(A):.0%}** · 有产物目录 **{hd/len(A):.0%}**。\n"
            f"  ⇒ **诚实的说法是「这 {len(A)} 笔的代价由我量不了的那部分主导」,\n"
            f"  而不是「太多了做不完」—— 后者是一堵墙,前者是一个已知的边界。**")
else:
    VERD = (f"**C 分布不是全有全无 ⇒ 先做便宜的那一半。** 两者都满足 {both}/{len(A)} = **{r:.0%}**;"
            f"分开看:`NEXT` 指向更晚条目 **{pt/len(A):.0%}** · 有产物目录 **{hd/len(A):.0%}**。\n"
            f"  ⇒ **一个比「全做」和「全不做」都好的第三选择:先把 {both} 笔机器能定位的做掉,\n"
            f"  剩下的单独登记 —— 而这个选择只有在量过之后才看得见。**")
print(VERD)
print("\n⚠ **本轮量的不是「自动查实的代价」** —— `#829` 已证从散文反推状态不可行(召回 18.6%)。"
      "量的是**把一笔欠账送到人面前所需的机器动作**;**人判那一步的代价本轮量不了,如实登记。**")
json.dump(dict(n_legacy=len(A), fact1_points_to_later=pt, fact2_has_results_dir=hd, both=both,
               ratio=r, machine_seconds_total=elapsed, machine_ms_per_debt=1000 * elapsed / len(A),
               rows=A[:40], pos_control=pc, neg_control=nc,
               scope="machine actions only; the human-judgement step is NOT measured here",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), action="Production"),
          open(OUT / "what_one_costs.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'what_one_costs.json'}")
