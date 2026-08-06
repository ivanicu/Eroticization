"""E02·A249·R632 — 她的性道德,和她的其他道德,是两个东西吗?

`#587` 的 NEXT。行动类型:**FRONTIER**。**回到对象侧**(墙的清单已收口)。
`#534`/`#579` 只算了**性题内**与**非性题内**两个量,**从没算过跨两者**。
而那正是这条结论最后一个没被问过的大问题:
**「一个人有一套性道德」——它是她整套道德的一部分,还是另一件东西?**

G1 ESTIMAND(先于方法):在**每一个调查年内**(年份固定 ⇒ 无趋势污染),
   `性内` · `非性内` · **`跨两者`**(四道性题 × 五道非性道德题 = 20 对)三个 |ρ| 中位,再取年份中位。
   全部题目二值化(`#535` 的 S-A 严切点),使三个量格式一致。
预注册:
   **跨 ≈ 0(< 安慰剂 q95)** ⇒ **两个独立的东西**;
   **跨 ≈ 非性内** ⇒ **性只是道德的一部分**,「一套性道德」不是一个单独的东西;
   **介于两者** ⇒ **报比值**,并写成「相关但可分」。
CONTROLS:正对照 = `premarsx`×`teensex`(**两道最接近的性题**)必须最高 ·
   安慰剂 = 每道题 × `zodiac`(星座)· 逐年 n 全部打印
IMPOSSIBLE:非性道德题**全部原生二值**(死刑/大麻/安乐死/自杀/女性从政)⇒
   它们**不构成一个「领域」**,只是「其他道德」的一个便利集合(`#536a` 里 NSFG 的家庭题才是一个领域)·
   一国一仪器 · 横断面 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEX = ["homosex", "premarsx", "xmarsex", "teensex"]
NON = ["cappun", "grass", "letdie1", "suicide1", "fepol"]
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "zodiac"] + SEX + NON, convert_categoricals=False)
B = {}
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
for c in NON:
    v = g[c].where(g[c].isin([1, 2])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
YR = g.year.values.astype(int)
print("=== 硬规则 1:逐题 n 与二值阳性率 ===")
for c in SEX + NON:
    print(f"  {c:9s} n={int(np.isfinite(B[c]).sum()):6d} 阳性率={np.nanmean(B[c]):.4f}")
def med(pairs, ym):
    o = []
    for a, b in pairs:
        m = ym & np.isfinite(B[a]) & np.isfinite(B[b])
        if m.sum() < 200 or np.std(B[a][m]) == 0 or np.std(B[b][m]) == 0: continue
        o.append(abs(float(np.corrcoef(rankdata(B[a][m]), rankdata(B[b][m]))[0, 1])))
    return float(np.median(o)) if o else np.nan
P_IN = list(itertools.combinations(SEX, 2))
P_NON = list(itertools.combinations(NON, 2))
P_X = [(a, b) for a in SEX for b in NON]
print(f"\n配对数:性内 {len(P_IN)} · 非性内 {len(P_NON)} · **跨 {len(P_X)}**")
vals = {"性内": [], "非性内": [], "跨两者": []}
years = []
for y in sorted(set(YR)):
    ym = YR == y
    if ym.sum() < 500: continue
    a, b, c = med(P_IN, ym), med(P_NON, ym), med(P_X, ym)
    if not all(np.isfinite(x) for x in (a, b, c)): continue
    years.append(y); vals["性内"].append(a); vals["非性内"].append(b); vals["跨两者"].append(c)
M = {k: float(np.median(v)) for k, v in vals.items()}
print(f"\n=== {len(years)} 个可算年份({min(years)}–{max(years)})===")
for k in vals: print(f"  {k:6s} 中位 = **{M[k]:.4f}**")
G = Gate("她的性道德,和她的其他道德,是两个东西吗?")
pc = []
for y in years:
    ym = YR == y
    v = med([("premarsx", "teensex")], ym)
    if np.isfinite(v): pc.append(v)
G.positive_control("正对照:premarsx×teensex(两道最接近的性题)为上限",
                   planted=float(np.median(pc)), floor=M["性内"], spread=1e-9)
zs = []
for c in SEX + NON:
    m = np.isfinite(B[c]) & g.zodiac.between(1, 12).values
    if m.sum() > 1000: zs.append(abs(float(np.corrcoef(rankdata(B[c][m]), rankdata(g.zodiac.values[m]))[0, 1])))
Z95 = float(np.quantile(zs, .95))
G.negative_control("安慰剂:每道题 × 星座", null=float(np.median(zs)), effect=M["性内"],
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
cells = {k: dict(n=len(years), median=M[k], values=[round(x, 4) for x in vals[k]],
                 inclusion=[k, f"{len(years)} 个调查年", "全部题目二值化(S-A 严切点)", "每对 n>=200"])
         for k in vals}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print(f"\n  安慰剂(题×星座)中位 {np.median(zs):.4f} · q95 {Z95:.4f}")
print("\n" + "=" * 70)
if np.median(pc) > M["性内"] and np.median(zs) < 0.5 * M["性内"]:
    x, nn = M["跨两者"], M["非性内"]
    if x < Z95:
        world = "TWO-THINGS"; verdict = f"跨两者 {x:.4f} 落在安慰剂 q95 {Z95:.4f} 内 -> **两个独立的东西**"
    elif abs(x - nn) < 0.2 * nn:
        world = "ONE-MORALITY"; verdict = f"跨两者 {x:.4f} ≈ 非性内 {nn:.4f} -> **性只是道德的一部分**"
    else:
        world = "RELATED-BUT-SEPARABLE"
        verdict = (f"跨两者 {x:.4f} 介于安慰剂 {Z95:.4f} 与非性内 {nn:.4f} 之间"
                   f"(占非性内 {x/nn:.0%},占性内 {x/M['性内']:.0%})-> **相关但可分**")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:非性道德题**全部原生二值**且彼此不属一个领域"
          "(死刑/大麻/安乐死/自杀/女性从政)—— 它们是「其他道德」的一个**便利集合**,不是一个领域;"
          "所以「跨两者」比的是「性 vs 一堆杂题」,不是「性 vs 另一个领域」。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(medians=M, per_year={k: [round(x, 4) for x in v] for k, v in vals.items()},
               years=[int(y) for y in years], n_pairs=dict(sex=len(P_IN), non=len(P_NON), cross=len(P_X)),
               placebo_median=float(np.median(zs)), placebo_q95=Z95, world=world, verdict=verdict,
               seeds=SEEDS,
               impossible=["非性题全部原生二值且不构成一个领域,是便利集合", "一国一仪器", "横断面"],
               unchallenged=True), open(OUT / "cross_domain.json", "w"), indent=1)
print(f"\nwrote {OUT/'cross_domain.json'}")
