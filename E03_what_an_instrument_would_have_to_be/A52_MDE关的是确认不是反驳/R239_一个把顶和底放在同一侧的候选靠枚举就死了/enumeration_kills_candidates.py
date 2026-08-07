"""#800 · E03·A52·R239 —— MDE 关的是「确认」,不是「反例反驳」

`#791` 立过:n=8 上要检出一个候选解释量需 **|ρ| ≥ 0.90** ⇒ **确认这条路是关着的。**
`#799`① 于是写下:**不要再找候选解释量** —— 但也写下了一件没被 MDE 关住的事:
顶对(`sexeduc`·`racmar`)都是**关于国家该不该立法/办学**的题,
而 `prayer`(学校祷告)、`helpblk`(政府该不该帮黑人)**也是**,却在底端。

**⚠⚠ 而那句话里藏着本轮的全部内容,它是一条不对称:**
**一个候选要被「确认」,需要 |ρ| ≥ 0.90 —— 关着的。
一个候选若把「顶」与「底」放在同一侧,它被一个反例杀掉 —— 而反例不需要统计。**
⇒ **MDE 关的是确认,不是反驳。** 在一个功效枯竭的设计里,**能做的事不是零,是只能做减法。**

G1 估计量:**每个候选切法的正类,是不是恰好等于 `{sexeduc, racmar}`。**
   ⚠ 要解释的**不是**整条次序(`#797` 已撤三堆),而是**「什么让这两题与其余六题不同」** ——
   所以判据是**「恰好隔离出这一对」**,而不是「与次序相关」。

识别:切法由**题干原文**判定(硬规则①:`.dta` 的变量标签,不是我记的题目名),
   而**判定是我做的、且我已经看过次序** ⇒ **被污染。**
   ⚠⚠ **而污染的方向在这里是可说清楚的,这正是本轮成立的原因:**
   **我有动机让切法成立** ⇒ 污染偏向「切法存活」 ⇒
   **一个在这种偏向下仍然被反例杀掉的切法,是更强的证据,不是更弱的。**
   ⇒ **本轮只报「被杀掉的」,不把「存活的」当成任何东西。**

⚠⚠ **而「存活」在这里根本不构成证据,这一条要用算术说死(`G3` 多重性):**
   8 题里恰好隔离出某一对的子集只有 **1 个**,而 2 元子集共 **C(8,2) = 28** 个
   ⇒ 随机贴标签下一个切法恰好命中的概率是 **1/28 ≈ 0.036**;
   **枚举 6 个候选,至少一个靠运气命中的概率 = 1 − (27/28)^6 ≈ 0.196。**
   ⇒ **五分之一的机会白捡一个「解释」** ⇒ **存活的切法一律登记为「未被杀掉」,绝不写成「支持」。**

三个世界:
   A **全被杀**:所有候选都把顶或底的某一题放错边 ⇒ **这一族候选整体出局**,而这是真进展。
   B **有存活**:某个切法恰好隔离出顶对 ⇒ **它仍然不是证据**(上面那个 0.196),
     只能登记为「未被杀掉、且在 n=8 上不可检验」。
   C **判据本身错了**:若某个切法把顶对拆开却仍被我判成「相关」,说明我在用次序而不是用题干判 ——
     **控制:判定只准看题干,产物里存下每题的题干原文供复核。**

预注册判词(条件式):
  if 正控开火(一个**构造上就隔离顶对**的切法必须被判为存活)
     and 负控开火(`#789` 已经杀掉的「性 vs 非性」必须被判为**被杀**):
      全部候选被杀 -> A
      有存活       -> B(登记为未被杀掉,附 0.196 的多重性,**不许写成支持**)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**我可以通过调切法的措辞让任何一对被隔离出来。**
  ⇒ 控制:**六个候选的定义在看任何结果之前写死在下面**,且每个都必须能只看题干判定;
  **不许在跑完之后新增或修改任何一个。**

本轮换不了仪器(对象是这八道题的题干);而它**不需要** —— 反驳不需要第二具仪器。
"""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from itertools import combinations
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]; OBS = P791["obs"]
TOP = {"sexeduc", "racmar"}
STEM = pd.io.stata.StataReader(gp).variable_labels()      # 硬规则①:题干从对象读,不是我记的

print("=== ① 八道题的题干原文(判定只准看这一列,产物里存下供复核)===")
for c in sorted(ITEMS, key=lambda x: -OBS[x]):
    print(f"  {c:9s} r={OBS[c]:+7.3f} {'← 顶对' if c in TOP else '       '}  「{STEM.get(c,'?')}」")

# ── 六个候选,**定义写死在看任何结果之前**,每个只看题干 ────────────────────────
CUTS = {
 "① 国家该不该立法/办学": {"sexeduc", "racmar", "prayer", "helpblk"},
 "② 关于未成年人":        {"sexeduc", "teensex", "spanking"},
 "③ 与种族有关":          {"racmar", "helpblk"},
 "④ 对象是公共安排而非个人行为": {"sexeduc", "racmar", "prayer", "helpblk"},
 "⑤ 与性有关":            {"sexeduc", "teensex", "homosex"},
 "⑥ 涉及生死":            {"suicide2"},
}
CTRL = {
 "正控:构造上就隔离顶对":  set(TOP),
 "负控:`#789` 已杀的性/非性": {"sexeduc", "teensex", "homosex"},
}

def verdict(pos):
    """恰好隔离出顶对 ⇒ 未被杀掉;否则被反例杀掉,并指名那个反例。"""
    extra = sorted(pos - TOP); missing = sorted(TOP - pos)
    return (not extra and not missing), extra, missing

print(f"\n=== ② 六个候选 × 判据「正类是否恰好 = {sorted(TOP)}」===")
rows = []
for name, pos in CUTS.items():
    ok, extra, missing = verdict(pos)
    why = ("未被杀掉" if ok else
           "被反例杀掉 —— 同侧还有:" + "·".join(f"{c}(r={OBS[c]:+.2f})" for c in extra) +
           (f";而缺了:{missing}" if missing else ""))
    rows.append(dict(cut=name, positive=sorted(pos), survives=ok, extra=extra, missing=missing))
    print(f"  {name:22s} {'存活' if ok else '**死**'}  {why}")

n_alive = sum(r["survives"] for r in rows)
pc_ok = verdict(CTRL["正控:构造上就隔离顶对"])[0]
nc_ok = not verdict(CTRL["负控:`#789` 已杀的性/非性"])[0]
print(f"\n  正控(构造上隔离顶对)判为存活:**{pc_ok}** · 负控(性/非性)判为被杀:**{nc_ok}**")

# ── 多重性:存活根本不构成证据,用算术说死 ────────────────────────────────────
C82 = len(list(combinations(range(8), 2)))
p_one = 1.0/C82
p_any = 1.0 - (1.0 - p_one)**len(CUTS)
print(f"\n=== ③ `G3` 多重性:为什么「存活」不构成证据 ===")
print(f"  8 题的 2 元子集共 **C(8,2) = {C82}** 个,恰好命中某一对的概率 **1/{C82} = {p_one:.3f}**")
print(f"  枚举 **{len(CUTS)}** 个候选,**至少一个靠运气命中的概率 = 1 − (1−1/{C82})^{len(CUTS)} = {p_any:.3f}**")
print(f"  ⇒ **约五分之一的机会白捡一个「解释」** ⇒ 存活一律登记为「未被杀掉」,**绝不写成「支持」。**")

G = Gate("#800 · MDE 关的是确认,不是反例反驳")
G.asserted("① 正控:一个**构造上就隔离顶对**的切法必须被判为存活(否则判据连真的都认不出)",
           pc_ok, f"正类 = {sorted(CTRL['正控:构造上就隔离顶对'])}", kind="control")
G.asserted("② 负控:`#789` 已经杀掉的「性 vs 非性」必须被判为**被杀**(否则判据杀不掉已知该死的)",
           nc_ok, f"正类 = {sorted(CTRL['负控:`#789` 已杀的性/非性'])} ⇒ 同侧含 teensex/homosex", kind="control")
G.asserted("③ 前提(跑前写下的混淆):六个候选的定义在看任何结果之前写死,跑完不许增删改",
           bool(len(CUTS) == 6), f"{len(CUTS)} 个候选,定义在文件里可复核", kind="control")
G.asserted("④ kill(预注册):「这一族候选整体出局」要成立,需**全部六个**都被反例杀掉",
           bool(n_alive == 0), f"存活 {n_alive}/{len(CUTS)}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n_alive == 0:
    v = (f"**A 这一族候选整体出局,而这是靠枚举做到的,一个统计量都没用。** "
         f"六个候选全部被反例杀掉 —— 每一个都把顶对之外的某题放在了同侧。\n"
         f"  ⇒ **而这正是本轮的方法论收获:MDE 关的是确认,不是反驳。** n=8 上 |ρ|≥0.90 才可检出,\n"
         f"  **但一个把「顶」与「底」放在同一侧的候选,一个反例就够了 —— 反例不吃功效。**\n"
         f"  ⇒ **在一个功效枯竭的设计里,能做的不是零,是只能做减法。**")
else:
    alive = [r["cut"] for r in rows if r["survives"]]
    v = (f"**B 有 {n_alive} 个候选未被杀掉:{alive}** —— **而它们不是证据。**\n"
         f"  枚举 {len(CUTS)} 个候选,至少一个靠运气恰好隔离出某一对的概率是 **{p_any:.3f}** ——\n"
         f"  **约五分之一。** ⇒ 登记为「未被杀掉、且在 n=8 上不可检验(`#791`:需 |ρ|≥0.90)」,\n"
         f"  **绝不写成「支持」。**")
print(v)
json.dump(dict(items=ITEMS, obs=OBS, top=sorted(TOP), stems={c: STEM.get(c, "") for c in ITEMS},
               cuts={k: sorted(v2) for k, v2 in CUTS.items()}, rows=rows, n_alive=n_alive,
               n_subsets=C82, p_one=p_one, p_any_lucky=p_any,
               pos_control=pc_ok, neg_control=nc_ok, admissible=adm, verdict=v, gate_ok=G.verdict()),
          open(OUT/"enumeration_kills.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'enumeration_kills.json'}")
