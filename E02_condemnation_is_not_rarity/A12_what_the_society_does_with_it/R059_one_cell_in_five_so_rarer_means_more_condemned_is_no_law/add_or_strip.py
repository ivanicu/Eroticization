"""E02·A12·R661 —— 那三十次措辞更正,是在「补」还是在「删」?

`#624` 的 NEXT。**行动类型:CLOSURE**(如实标注)。

⚠ **§3 梯度检查改了做法,而这一条写在前面:**
   `#624` 的 NEXT 写的是「**由我逐条归类**」——**而那本身就是一台判断仪器**,
   它没有正对照、不可复现、而且我是被测对象的作者。
   ⇒ 改成**机械口径**:用**已有的** `_QUAL`(`readme_ledger_audit.py`)数每次 diff 里
     **增加/删除的限定语行**,`net = 加 − 删`。
   ⇒ 真正需要语义判断的那一类(**换掉一个范围词**)**机械上分不开**,
     **单独标 `UNCOMPUTED`,并说明它需要什么** —— 不许标「planned」。

G1 ESTIMAND(先于方法):对每一次措辞更正提交,在 `README.md` 与 `README_zh.md` 的 diff 上:
  `q_add` = 匹配 `_QUAL` 的 `+` 行数 · `q_del` = 匹配 `_QUAL` 的 `−` 行数 · `net = q_add − q_del`。
  归类(**先写死**):`net>0` -> **(b) 补限定** · `net<0` -> **(c) 删限定** ·
  `net==0 且 q_del>0` -> **(a/d) 改写或换范围词(机械上分不开)** · `net==0 且 q_del==0` -> **不涉及限定语**。
KILL(条件式,预注册于 `#624`,但因 (d) 不可计算而如实调整,调整写在这里):
  原判据是 `(c)+(d) >= 5` / `<= 2`。**(d) 机械上不可计算** ⇒
  **改为对 (c) 单独评判,并把 (a/d) 桶的大小作为 (d) 的上界一起报**:
  `(c) >= 5` -> 「措辞缺陷不可检测」为假,值得再造工具;
  `(c) <= 2` **且** `(a/d) 桶 <= 2` -> **如实记「这个项目的措辞更正几乎全是补,不是删」**;
  否则 -> 报区间,不下判决。
CONTROLS:
  正对照:合成一段 diff,**删掉一条带限定语的行** -> 必须判为 **(c)**。
  **g=0**:合成一段**不含任何限定语行**的 diff -> 必须判为 **不涉及限定语**。
  安慰剂:合成一段**只增加**限定语行的 diff -> 必须判为 **(b)**,不得判为 (c)。
G3:30 条逐条发布(`q_add`/`q_del`/`net`/类别),**让下一个人能推翻我**。
G4:是否把中文页一起计入 {只英文 / 两版合计} 两档。
IMPOSSIBLE(不写 planned):**`_QUAL` 是一个词表,它认得的「限定语」不是「限定语」本身** ·
  **(d) 换范围词** 需要**逐句语义比对**,机械口径给不了 —— 它需要的是**一个能判断两句话是不是同一主张、
  而范围词被换掉了的判定器**,本仓库没有 · 条目→提交映射 `git log -S` 已知 6 次失配 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, subprocess, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r'(?<![\w.])\d+\.\d{2,4}(?![\w])')
QUAL = A._QUAL
sh = lambda *a: subprocess.run(a, capture_output=True, text=True).stdout


def classify(diff):
    add = [l[1:] for l in diff.splitlines() if l.startswith("+") and not l.startswith("+++")]
    dele = [l[1:] for l in diff.splitlines() if l.startswith("-") and not l.startswith("---")]
    qa = sum(1 for l in add if QUAL.search(l)); qd = sum(1 for l in dele if QUAL.search(l))
    net = qa - qd
    if net > 0: k = "(b) 补限定"
    elif net < 0: k = "(c) 删限定"
    elif qd > 0: k = "(a/d) 改写或换范围词(机械上分不开)"
    else: k = "不涉及限定语"
    return qa, qd, net, k


# ── 控制先跑 ─────────────────────────────────────────────────
G = Gate("那三十次措辞更正,是在补还是在删?")
POS = "--- a\n+++ b\n-⚠ 这个 +0.432 only 在同一份问卷内部成立,cannot 外推。\n 别的句子。"
G0 = "--- a\n+++ b\n-一句普通陈述句。\n+另一句普通陈述句。"
PLA = "--- a\n+++ b\n+⚠ 新增:这个 +0.432 only 在同一份问卷内部成立。\n 别的句子。"
cp, cg, cl = classify(POS), classify(G0), classify(PLA)
print(f"=== 控制 ===")
print(f"  正对照(删一条限定语行)-> {cp[3]}  (q_add {cp[0]} · q_del {cp[1]} · net {cp[2]})")
print(f"  g=0(无任何限定语行)  -> {cg[3]}  (q_add {cg[0]} · q_del {cg[1]})")
print(f"  安慰剂(只增加限定语行)-> {cl[3]}  (q_add {cl[0]} · q_del {cl[1]})")
pos_ok = G.positive_control("正对照:删一条限定语行必须判为 (c)",
                            planted=float(cp[3].startswith("(c)")), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:无限定语行的 diff 不得判为 (c)",
                            null=float(cg[3].startswith("(c)")), effect=1.0,
                            null_spread=0.4, null_kind="不含任何限定语行的 diff")

# ── 主测 ────────────────────────────────────────────────────
led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
PAT = re.compile(r'两份 README 已改|两版 README|页面已改|页面.{0,6}更正|已更正.{0,10}页面|正确写法')
CAND = sorted({ent[j] for j, l in enumerate(led) if PAT.search(l) and ent[j] > 0})
rows = []
for N in CAND:
    log = sh("git", "log", "--format=%H", "-S", f"## Entry {N} ", "--", "RETRACTIONS.md").split()
    if not log: rows.append(dict(entry=N, cls="找不到提交")); continue
    c = log[-1]
    d_en = sh("git", "show", c, "--", "README.md")
    rem = {n for l in d_en.splitlines() if l.startswith("-") and not l.startswith("---") for n in NUM.findall(l)}
    add = {n for l in d_en.splitlines() if l.startswith("+") and not l.startswith("+++") for n in NUM.findall(l)}
    if rem - add: rows.append(dict(entry=N, cls="数字更正(不在本轮总体)")); continue
    both = d_en + sh("git", "show", c, "--", "README_zh.md")
    qa, qd, net, k = classify(both)
    qa1, qd1, net1, k1 = classify(d_en)
    rows.append(dict(entry=N, q_add=qa, q_del=qd, net=net, cls=k, cls_en_only=k1))
T = pd.DataFrame(rows)
W = T[T.cls.isin(["(b) 补限定", "(c) 删限定", "(a/d) 改写或换范围词(机械上分不开)", "不涉及限定语"])]
print(f"\n=== G3:{len(T)} 条候选,其中措辞更正 {len(W)} 条 —— 逐条发布 ===")
print(T.cls.value_counts().to_string())
print("\n  逐条(让下一个人能推翻我):")
for r in W.itertuples():
    print(f"    #{int(r.entry):4d} q_add {int(r.q_add):3d} · q_del {int(r.q_del):3d} · net {int(r.net):+4d} · {r.cls}")
n_c = int((W.cls == "(c) 删限定").sum()); n_ad = int((W.cls.str.startswith("(a/d)")).sum())
n_b = int((W.cls == "(b) 补限定").sum()); n_none = int((W.cls == "不涉及限定语").sum())
print(f"\n**(b) 补 {n_b} · (c) 删 **{n_c}** · (a/d) 分不开 **{n_ad}**(= (d) 的上界)· 不涉及 {n_none}**")

if pos_ok and pla_ok:
    if n_c >= 5: verdict = f"**(c) = {n_c} ≥ 5 ⇒ 「措辞缺陷不可检测」为假,值得再造工具**"
    elif n_c <= 2 and n_ad <= 2: verdict = (f"**(c) = {n_c} ≤ 2 且 (a/d) = {n_ad} ≤ 2 ⇒ "
                                            f"这个项目的措辞更正几乎全是补,不是删;`#623d` 的缺口是理论上的,不是实测到的**")
    else: verdict = f"报区间:(c) = {n_c} · (a/d) 上界 {n_ad} ⇒ **不下判决**"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4:只英文 vs 两版合计 ===")
for col, nm in (("cls_en_only", "只英文"), ("cls", "两版合计")):
    vc = W[col].value_counts().to_dict()
    print(f"  {nm}: (b) {vc.get('(b) 补限定',0)} · (c) **{vc.get('(c) 删限定',0)}** · "
          f"(a/d) {vc.get('(a/d) 改写或换范围词(机械上分不开)',0)} · 不涉及 {vc.get('不涉及限定语',0)}")
json.dump(dict(table=T.to_dict("records"), n_b=n_b, n_c=n_c, n_ad=n_ad, n_none=n_none,
               controls=dict(pos=cp[3], g0=cg[3], placebo=cl[3]), verdict=verdict,
               d_is="UNCOMPUTED —— 机械上分不开,需要一个能判断「同一主张而范围词被换掉」的判定器",
               unchallenged=True),
          open(OUT/"add_or_strip.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'add_or_strip.json'}")
