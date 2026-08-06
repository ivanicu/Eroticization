"""E02·A217·R587 — 「性」在更细的尺度上还是一整块吗?

`#541` 的 NEXT(541d③ 的转向)。行动类型:**FRONTIER**。
`#536`/`#539` 说「性是一个更紧的领域」。本轮把尺度调细一格:**性内部有没有子面?**
如果在一份**自选、非概率、非政府**的样本里它清楚地分成若干块,
那么「性是一个整块」就是**尺度依赖的**,不是一个事实。

⚠ **按 `#541c` 自己的教训,归类必须从题目文本推,不能从我记得的量表结构推。**
   码本**没有**写子面结构,**但 100 道题的文本全在**。
   可验证的机械规则:**题 i 与题 i+20 语义配对**(Q1/Q21 焦虑 · Q2/Q22 胜任 · Q3/Q23 觉察 ·
   Q6/Q26 全神贯注 · Q17/Q37 恐惧 · Q20/Q100 主责)⇒ **20 个子面 × 5 题**。
   **本轮把 20 组的题目文本全部打印进产物,作为归类的证据** —— 归类来自文本,检验来自数据,**不循环**。

G1 ESTIMAND(先于方法):`ρ_in` = 子面**内**两两 |ρ| 中位(20 组 × C(5,2)=10 对 = 200 对);
   `ρ_out` = 子面**间**两两 |ρ| 中位(4950−200 = 4750 对)。**主量 = `ρ_in − ρ_out`。**

⚠ 硬规则 2:**自选样本的选择效应会同时抬高所有相关**(愿意做在线性量表并同意存档的人不是随机的)。
   ⇒ **只比 `ρ_in − ρ_out`,绝不与 NSFG 的绝对值比。**

WORLDS:
  W-FACETED  `ρ_in` 明显 > `ρ_out` ⇒ **性内部有清楚的子面 ⇒「一整块」是尺度依赖的**
  W-ONE-BLOCK 两者接近 ⇒ **在这具仪器上性就是一块,子面结构不显**
  W-NOISE    随机归类也给出同样的差 ⇒ 差是**组大小的产物**,不是结构
⚠ BASIN:`W-FACETED` 会让「尺度依赖」这个漂亮说法成立,**不是**下注方向。本轮下注 `W-ONE-BLOCK`。

CONTROLS(G2):
  **关键零:随机重新归类** —— 把 100 题随机打成 20×5,重算 `ρ_in − ρ_out`,300 次。
     「这个零该不该是零?」**该** ⇒ `negative_control`,零的方案 = 题目→子面标签置换;
  正对照 文本上最像重复的那一对(Q3 觉察感受 × Q23 觉察动机)必须处在高位;
  ⚠ 反向计分题(题干带 (R))**不翻转** —— 本轮全部取 |ρ|,方向不进入任何量。
KILL(条件式):if 正对照 > `ρ_out` and 随机归类的差 ≈ 0:
     实际差 > 随机归类差的 q95 -> W-FACETED;落在其内 -> W-ONE-BLOCK else UNVERIFIED
IMPOSSIBLE:自选样本 ⇒ 无人群外推 · 单一量表 ⇒ 子面定义由该量表作者决定 ·
   横断面非因果 · 归类规则由我从文本读出,**未经第二人复核** ⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, re, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
BASE = ROOT / "data/external/openpsych/MSSCQ/MSSCQ"
txt = (BASE / "codebook.txt").read_text(errors="replace")
ITEMS = {}
for m in re.finditer(r'^(\d{1,3})\.\s+(.+)$', txt, re.M):
    k = int(m.group(1))
    if 1 <= k <= 100 and k not in ITEMS: ITEMS[k] = m.group(2).strip()
print(f"=== 硬规则 1:从码本读出 {len(ITEMS)}/100 道题的文本 ===")
assert len(ITEMS) == 100, f"只读到 {len(ITEMS)} 道题 —— 解析不全,不得继续"
D = pd.read_csv(BASE / "data.csv", sep="\t" if "\t" in open(BASE / "data.csv").readline() else ",")
QC = [c for c in D.columns if re.fullmatch(r"Q\d{1,3}", c)]
print(f"  data.csv:{len(D)} 行 · {len(QC)} 个 Q 列")
X = D[QC].apply(pd.to_numeric, errors="coerce")
X = X.where(X.isin([1, 2, 3, 4, 5]))
print(f"  有效作答率(全 100 题非缺失的人)= {int(X.notna().all(1).sum())}/{len(X)}")
FACET = {q: (q - 1) % 20 for q in range(1, 101)}   # 题 i 与 i+20 同面 —— 由文本推出的机械规则
print("\n=== 归类证据:20 个子面,每面 5 道题的文本(归类来自文本,不来自记忆)===")
groups = {}
for f in range(20):
    qs = [q for q in range(1, 101) if FACET[q] == f]
    groups[f] = qs
    print(f"  面 {f:2d}: " + " | ".join(f"Q{q}:{ITEMS[q][:34]}" for q in qs[:3]))
V = {q: X[f"Q{q}"].values.astype(float) for q in range(1, 101) if f"Q{q}" in X}
def rho(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 500: return np.nan
    return abs(float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]))
print("\n  计算 4950 对相关…", flush=True)
R = {}
for a, b in itertools.combinations(sorted(V), 2):
    r = rho(V[a], V[b])
    if np.isfinite(r): R[(a, b)] = r
inn = [r for (a, b), r in R.items() if FACET[a] == FACET[b]]
out = [r for (a, b), r in R.items() if FACET[a] != FACET[b]]
IN, OUT_ = float(np.median(inn)), float(np.median(out))
print(f"  子面**内** {len(inn)} 对,中位 |ρ| = **{IN:.4f}**")
print(f"  子面**间** {len(out)} 对,中位 |ρ| = **{OUT_:.4f}**")
print(f"  **差 = {IN-OUT_:+.4f}**")
# 关键零:随机重新归类
null = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(100):
        perm = rng.permutation(100) + 1
        FP = {int(perm[i]): i % 20 for i in range(100)}
        i2 = [r for (a, b), r in R.items() if FP[a] == FP[b]]
        o2 = [r for (a, b), r in R.items() if FP[a] != FP[b]]
        if i2 and o2: null.append(np.median(i2) - np.median(o2))
null = np.array(null); Q95 = float(np.quantile(null, .95))
print(f"  随机重新归类(300 次)的差:中位 {np.median(null):+.5f} · q95 {Q95:+.5f}")
G = Gate("「性」在更细的尺度上还是一整块吗?(MSSCQ,自选网络样本)")
pc = R.get((3, 23), np.nan)
G.positive_control("正对照:Q3 觉察感受 × Q23 觉察动机(文本上最像重复)",
                   planted=pc, floor=OUT_, spread=1e-9)
G.negative_control("关键零:随机重新归类后的差必须 ≈ 0", null=float(np.median(null)),
                   effect=IN - OUT_, null_spread=float(np.std(null)), null_kind="题目→子面标签置换")
cells = {f"面{f}": dict(n=int(np.isfinite(V[groups[f][0]]).sum()), facet=f,
                        items=[f"Q{q}" for q in groups[f]],
                        median_in=float(np.median([R[(a, b)] for a, b in itertools.combinations(groups[f], 2)
                                                   if (a, b) in R])),
                        inclusion=[f"面 {f} 的 5 道题,C(5,2)=10 对", "每对 n>=500", "由 i≡i+20 文本规则归类"])
         for f in range(20)}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n=== 20 个子面各自的面内中位(全格公布)===")
for f in range(20):
    print(f"  面 {f:2d}: {cells[f'面{f}']['median_in']:.4f}  ({ITEMS[groups[f][0]][:44]})")
print("\n" + "=" * 76)
if pc > OUT_ and abs(np.median(null)) < 0.5 * abs(IN - OUT_):
    world = "W-FACETED" if (IN - OUT_) > Q95 else "W-ONE-BLOCK"
    verdict = (f"差 {IN-OUT_:+.4f} {'>' if world=='W-FACETED' else '<='} 随机归类 q95 {Q95:+.5f} -> "
               + ("**性内部有清楚的子面 ⇒「一整块」是尺度依赖的**" if world == "W-FACETED"
                  else "**在这具仪器上性就是一块**"))
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:子面的定义来自**该量表作者的设计**,"
          "而我只是从文本把它读了出来 —— 若作者的子面本身就是按「相关高」编出来的,"
          "那么面内高于面间是**设计的复现**,不是关于性心理的发现。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(median_in=IN, median_out=OUT_, diff=IN - OUT_, null_median=float(np.median(null)),
               null_q95=Q95, n_pairs_in=len(inn), n_pairs_out=len(out), n_resp=int(len(X)),
               world=world, verdict=verdict, facet_cells=cells, item_texts=ITEMS,
               grouping_rule="题 i 与 i+20 同面(由题目文本读出,非记忆);20 面 × 5 题",
               positive_control=dict(pair="Q3×Q23", rho=float(pc) if np.isfinite(pc) else None),
               seeds=SEEDS, instrument="Open Psychometrics MSSCQ,网络自选非概率样本",
               impossible=["自选样本无人群外推", "子面定义由量表作者决定", "横断面非因果",
                           "归类规则未经第二人复核"], unchallenged=True),
          open(OUT / "facets_inside_sex.json", "w"), indent=1)
print(f"\nwrote {OUT/'facets_inside_sex.json'}")
