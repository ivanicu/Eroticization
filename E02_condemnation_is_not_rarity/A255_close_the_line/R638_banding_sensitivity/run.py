"""E02·A255·R638 — 那两条斜率,换一种分段还在吗?

`#593` 的 NEXT ③(`#593d` 欠着的那次)。行动类型:**CLOSURE**(补做敏感性,不分离新世界)。
`#593b` 的两条斜率(A 随年龄 +0.001857 · B 随队列 +0.001585)建立在**我划的**年龄段与队列段上,
而 `#575` 已证明:**划分会改变结论,而且能让符号翻转。**

预注册(先于计算):用 `lib/assignment_sensitivity`(符号优先于极差),
   **两条斜率各自的跨方案极差 ≥ 0.0005/年** ⇒ `#593b` 必须带上「分段敏感」;
   **< 0.0005** ⇒ 分段不承重;**任一方案下符号翻转** ⇒ 该斜率**必须撤回**。
三种站得住的分段(先于计算写死):
   S1 现行(年龄 4 段 / 队列 5 段)· S2 年龄细分(5 段)· S3 边界整体平移(年龄 4 段,断点各挪 4–5 岁)
IMPOSSIBLE:CLOSURE 不分离世界 · APC 共线未变(`#593a`)· 横断面 · 格 n<400 剔除 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
from lib.assignment_sensitivity import assignment_sensitivity
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEX = ["homosex", "premarsx", "xmarsex", "teensex"]
DOMS = {"D1": ["cappun", "grass", "letdie1", "suicide1", "fepol"],
        "D2": ["abany", "abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle"],
        "D3": ["letdie1", "suicide1"]}
ALL = SEX + sorted({x for v in DOMS.values() for x in v})
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "age", "cohort"] + ALL, convert_categoricals=False)
CO = g.cohort.where(g.cohort.between(1880, 2006)).values.astype(float)
AG = g.age.where(g.age.between(18, 89)).values.astype(float)
B = {}
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
for c in set(x for v in DOMS.values() for x in v):
    v = g[c].where(g[c].isin([1, 2])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
def med(pairs, m):
    o = []
    for a, b in pairs:
        k = m & np.isfinite(B[a]) & np.isfinite(B[b])
        if k.sum() < 150 or np.std(B[a][k]) == 0 or np.std(B[b][k]) == 0: continue
        o.append(abs(float(np.corrcoef(rankdata(B[a][k]), rankdata(B[b][k]))[0, 1])))
    return float(np.median(o)) if o else np.nan
def cell(m):
    v = [med([(s, d) for s in SEX for d in D], m) for D in DOMS.values()]
    v = [x for x in v if np.isfinite(x)]
    return float(np.median(v)) if v else np.nan
def slopes(AB, CB, minn=400):
    grid = {}
    for lo, hi in AB:
        grid[(lo, hi)] = [None if ((AG >= lo) & (AG <= hi) & (CO >= c0) & (CO <= c1)).sum() < minn
                          else cell((AG >= lo) & (AG <= hi) & (CO >= c0) & (CO <= c1)) for c0, c1 in CB]
    amid = np.array([(a + b) / 2 for a, b in AB], float)
    cmid = np.array([(a + b) / 2 for a, b in CB], float)
    A, Bs = [], []
    for j in range(len(CB)):
        col = [(amid[i], grid[AB[i]][j]) for i in range(len(AB)) if grid[AB[i]][j] is not None]
        if len(col) >= 3: A.append(np.polyfit([x for x, _ in col], [y for _, y in col], 1)[0])
    for i in range(len(AB)):
        row = [(cmid[j], grid[AB[i]][j]) for j in range(len(CB)) if grid[AB[i]][j] is not None]
        if len(row) >= 3: Bs.append(np.polyfit([x for x, _ in row], [y for _, y in row], 1)[0])
    return (float(np.median(A)) if A else np.nan, float(np.median(Bs)) if Bs else np.nan,
            sum(1 for v in grid.values() for x in v if x is not None))
S1 = ([(18, 34), (35, 49), (50, 64), (65, 89)],
      [(1900, 1929), (1930, 1944), (1945, 1959), (1960, 1974), (1975, 2006)])
S2 = ([(18, 29), (30, 39), (40, 49), (50, 64), (65, 89)],
      [(1900, 1929), (1930, 1944), (1945, 1959), (1960, 1974), (1975, 2006)])
S3 = ([(18, 30), (31, 45), (46, 60), (61, 89)],
      [(1900, 1934), (1935, 1949), (1950, 1964), (1965, 1979), (1980, 2006)])
SCH = {"S1 现行": S1, "S2 年龄细分": S2, "S3 边界平移": S3}
print("=== 三种分段(全格公布)===")
res = {}
for k, (AB, CB) in SCH.items():
    a, b, n = slopes(AB, CB)
    res[k] = dict(slope_age=a, slope_cohort=b, n_cells=n,
                  inclusion=[k, f"{len(AB)}×{len(CB)} 格", f"可算 {n} 格", "格 n>=400"])
    print(f"  {k:10s} 可算 {n:2d} 格 · A(随年龄)={a:+.6f} · B(随队列)={b:+.6f}")
rA = assignment_sensitivity(lambda m: res[m]["slope_age"], {k: k for k in SCH}, tol=0.0005)
rB = assignment_sensitivity(lambda m: res[m]["slope_cohort"], {k: k for k in SCH}, tol=0.0005)
print(f"\n  A 极差 = {rA['range']:.6f} · 判决 **{rA['verdict']}** · 符号 {rA['signs']}")
print(f"  B 极差 = {rB['range']:.6f} · 判决 **{rB['verdict']}** · 符号 {rB['signs']}")
G = Gate("那两条斜率,换一种分段还在吗?")
G.positive_control("三种分段都算得出(非退化)",
                   planted=float(sum(np.isfinite(res[k]["slope_age"]) for k in SCH)), floor=2.5, spread=1e-9)
G.spec_curve_cells_declare_n("规格曲线逐格 n", {k: dict(n=v["n_cells"], **v) for k, v in res.items()})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", res)
print("\n" + "=" * 70)
flip = "方向翻转" in (rA["verdict"], rB["verdict"])
sens = rA["range"] >= 0.0005 or rB["range"] >= 0.0005
if flip:
    world, verdict = "RETRACT", f"某一方案下符号翻转(A {rA['signs']} · B {rB['signs']})-> **该斜率必须撤回**"
elif sens:
    world, verdict = "BANDING-SENSITIVE", f"极差 A {rA['range']:.6f} · B {rB['range']:.6f} ≥ 0.0005 -> **`#593b` 必须带上「分段敏感」**"
else:
    world, verdict = "ROBUST", f"极差 A {rA['range']:.6f} · B {rB['range']:.6f} < 0.0005 -> **分段不承重**"
print(f"CLOSURE 结论:**{world}** —— {verdict}")
print("⚠ 这是 CLOSURE:它没有分离任何新世界,只检验了一个我做过的划分是否承重;"
      "而 APC 共线(`#593a`)不受本轮影响。")
print(G)
json.dump(dict(schemes=res, range_age=rA["range"], range_cohort=rB["range"],
               verdict_age=rA["verdict"], verdict_cohort=rB["verdict"],
               world=world, decision=verdict, tol=0.0005, seeds=SEEDS,
               impossible=["CLOSURE 不分离世界", "APC 共线未变", "横断面", "格 n<400 剔除"],
               unchallenged=True), open(OUT / "banding.json", "w"), indent=1)
print(f"\nwrote {OUT/'banding.json'}")
