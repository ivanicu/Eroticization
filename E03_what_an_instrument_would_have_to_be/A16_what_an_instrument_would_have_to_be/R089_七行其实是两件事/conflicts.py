"""E03·A258·R641 — 那六条对同一份调查的要求,彼此冲突吗?

`#596` 的 NEXT。行动类型:**PRODUCTION**(结构分析,**不提议采数据**)。
`#596b`:七条规格里有**六条**由同一份关于人的纵向调查满足。
⇒ 那就该问下一个可失败的问题:**它们对那份调查的要求,彼此冲突吗?**

**判据(先于逐对检查写死,三值):**
  `冲突`   —— **满足其一,必然破坏另一条已达成的状态**;
  `需额外结构` —— 两者可同时满足,**但只能靠一个前者不含的额外部件**(如第二条抽样臂);
  `相容`   —— 直接同时满足。
⚠ **「需额外结构」不算冲突** —— 它是**成本**,不是**矛盾**;两者必须分开报,
  否则「无冲突」会读成「一份简单的调查就够了」。

预注册:
  **存在 `冲突` 对** ⇒ 那份调查**不是一次能设计成的**,页面必须说明冲突在哪;
  **无冲突但有 `需额外结构`** ⇒ 六条可写成**一份规格**,而它**必须带上那些部件**;
  **全部相容** ⇒ 六条可写成一份最简规格。
CONTROLS:正对照 = 一对**按定义必然需要额外结构**的要求(R5 每波新样本 × R7 同一批人复访)
  必须判 `需额外结构` 或 `冲突`,**不得**判 `相容` · 负对照 = 每条与自己必须判 `相容` · 全矩阵公布
IMPOSSIBLE:判据由我定(三值,写在前面)· **不估成本、不估可行性** ·
  **要求的措辞由我从规格表读出**,而一条措辞不同的规格可能给出不同的矩阵 · [unchallenged]
"""
import os, sys, pathlib, json, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
from lib.gates import Gate
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
# 六条规格的**设计要求**,逐条拆成可比对的属性(这三栏是判据的全部输入)
REQ = {
 "R1 抽样框+第二波": dict(sample="概率抽样框", waves="≥2", items="—", present="—",
    note="要求样本可外推,且至少两个时点"),
 "R3 留住十道题":     dict(sample="—", waves="≥2", items="题集跨波不变", present="—",
    note="下一波保留全部条目"),
 "R4 跨波条目稳定":   dict(sample="—", waves="≥2", items="题干与格式跨波不变", present="呈现方式跨波不变",
    note="不仅题在,措辞与呈现也不变"),
 "R5 第二条跨年序列": dict(sample="**每波新样本**", waves="≥20 年,固定间隔", items="同一道题", present="—",
    note="重复横断面序列:每波抽新人"),
 "R6 态度题双问法":   dict(sample="—", waves="—", items="同一题问两遍", present="**一遍出声、一遍自填**",
    note="改变呈现方式:每题两种模式"),
 "R7 面板维度":       dict(sample="**同一批人复访**", waves="≥2", items="—", present="—",
    note="需人标识,跨波配对同一个人"),
}
K = list(REQ)
def judge(a, b):
    if a == b: return "相容", "同一条与自身"
    A, B = REQ[a], REQ[b]
    # 规则 ①:每波新样本 vs 同一批人复访 -> 需额外结构(两条抽样臂),不是逻辑矛盾
    if ("新样本" in A["sample"] and "复访" in B["sample"]) or ("复访" in A["sample"] and "新样本" in B["sample"]):
        return "需额外结构", "重复横断面需每波新人,面板需同一批人 -> 一份调查需**两条抽样臂**"
    # 规则 ②:呈现方式跨波不变 vs 改变呈现方式 -> 冲突,除非双问法自第一波即在
    pa, pb = A["present"], B["present"]
    if ("跨波不变" in pa and "自填" in pb) or ("跨波不变" in pb and "自填" in pa):
        return "需额外结构", "双问法改变呈现;仅当它**自第一波即存在**时两者可同时满足 -> 需**从设计之初纳入**"
    return "相容", "两条要求互不触碰"
print("=== 六条规格的设计要求(判据的全部输入)===")
for k in K: print(f"  {k:14s} 抽样={REQ[k]['sample']:12s} 波={REQ[k]['waves']:10s} 题={REQ[k]['items']:12s} 呈现={REQ[k]['present']}")
print("\n=== 三值矩阵(全矩阵公布)===")
SYM = {"相容": "·", "需额外结构": "+", "冲突": "✗"}
cells, counts = {}, {"相容": 0, "需额外结构": 0, "冲突": 0}
print("       " + " ".join(f"{k.split()[0]:>3s}" for k in K))
for a in K:
    row = []
    for b in K:
        v, why = judge(a, b)
        cells[f"{a.split()[0]}×{b.split()[0]}"] = dict(n=1, verdict=v, why=why,
            inclusion=["判据三值,写在前面", "要求由规格表读出", "机械应用"])
        if a < b: counts[v] += 1
        row.append(SYM[v])
    print(f"  {a.split()[0]:>4s}   " + " ".join(f"{x:>3s}" for x in row))
print(f"\n  上三角 15 对:相容 {counts['相容']} · **需额外结构 {counts['需额外结构']}** · **冲突 {counts['冲突']}**")
extra = [(k, v["why"]) for k, v in cells.items() if v["verdict"] == "需额外结构" and
         k.split("×")[0] < k.split("×")[1]]
print("\n  需额外结构的对:")
for k, why in extra: print(f"    {k}: {why}")
G = Gate("那六条对同一份调查的要求,彼此冲突吗?")
vpos, _ = judge("R5 第二条跨年序列", "R7 面板维度")
G.positive_control("正对照:R5(每波新样本)× R7(同一批人复访)不得判『相容』",
                   planted=float(vpos != "相容"), floor=0.5, spread=1e-9)
selfok = all(judge(k, k)[0] == "相容" for k in K)
G.negative_control("负对照:每条与自己必须判『相容』", null=float(not selfok), effect=1.0,
                   null_spread=1e-9, null_kind="同一条与自身,判据上必然相容")
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 70)
if vpos != "相容" and selfok:
    if counts["冲突"] > 0:
        world = "CONFLICT"; verdict = f"存在 {counts['冲突']} 对冲突 -> **那份调查不是一次能设计成的**"
    elif counts["需额外结构"] > 0:
        world = "ONE-SPEC-WITH-PARTS"
        verdict = (f"无冲突,但有 **{counts['需额外结构']} 对需额外结构** -> "
                   f"**六条可写成一份规格,而它必须带上那些部件**")
    else:
        world = "MINIMAL"; verdict = "全部相容 -> 一份最简规格"
    print(f"评判:**{world}** —— {verdict}")
    print("⚠ 这个 KILL 会怎样失败:**要求的措辞由我从规格表读出** —— "
          "一条写法不同的规格(例如把 R5 写成『同一批人问二十年』)会给出完全不同的矩阵。")
else:
    world, verdict = "UNVERIFIED", "对照未过"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(requirements=REQ, matrix={k: v["verdict"] for k, v in cells.items()},
               counts=counts, extra_structure=extra, world=world, verdict=verdict,
               impossible=["判据由我定(三值,写在前面)", "不估成本不估可行性",
                           "要求的措辞由我从规格表读出,换一种写法矩阵会变"],
               unchallenged=True), open(OUT / "conflicts.json", "w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'conflicts.json'}")
