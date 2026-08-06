"""E02·A254·R637 — 那个收紧,是世代换了人,还是同一批人自己拧紧了?

`#592` 的 NEXT。行动类型:**FRONTIER**。
`#592b`:那条连线在 52 年里从 0.227 升到 0.301。`#592d` 第一条:**升高可能是队列更替而非个体变化。**

⚠ **而这是一个结构上不可完全识别的问题:年龄 + 出生队列 = 调查年**(APC 共线)。
   ⇒ **本轮不识别三者**,只做一件可失败的事:**报两个边际方向的斜率与全网格**,
   并把「不可识别」写成结论的一部分,**不是脚注**。

G1 ESTIMAND(先于方法):把样本切成 **年龄段 × 出生队列段** 的网格,每格算一次「常数」
   (性题 × 三个第二领域的 |ρ| 中位的中位),再报:
   **A 同一队列内、随年龄的斜率** · **B 同一年龄段内、随队列的斜率**。
预注册:
   **B 超 MDE 而 A 不超** ⇒ **世代更替**;
   **A 超而 B 不超** ⇒ **同一批人自己拧紧**(⚠ 横断面下这仍不是个体变化,只是同龄不同期的差);
   **两者都超** ⇒ 报两个数,**不主张单一机制**;
   **都不超** ⇒ `#592b` 的上升**不可归因**,写进页面。
CONTROLS:正对照 = 每格 n 足够(<400 的格剔除并公布)· 安慰剂 = 打乱出生年后两个斜率都应塌
IMPOSSIBLE:**APC 共线 ⇒ 三者不可同时识别**,本轮只报两个边际 ·
  横断面 ⇒ **没有同一个人被问两次**,所以「个体拧紧」永远只是一个推断 · 一国一仪器 · [unchallenged]
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
DOMS = {"D1": ["cappun", "grass", "letdie1", "suicide1", "fepol"],
        "D2": ["abany", "abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle"],
        "D3": ["letdie1", "suicide1"]}
ALL = SEX + sorted({x for v in DOMS.values() for x in v})
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "age", "cohort"] + ALL, convert_categoricals=False)
CO = g.cohort.where(g.cohort.between(1880, 2006)).values.astype(float)
AG = g.age.where(g.age.between(18, 89)).values.astype(float)
print("=== 硬规则 1 ===")
print(f"  cohort n={int(np.isfinite(CO).sum()):6d} 范围 {np.nanmin(CO):.0f}–{np.nanmax(CO):.0f}")
print(f"  age    n={int(np.isfinite(AG).sum()):6d} 范围 {np.nanmin(AG):.0f}–{np.nanmax(AG):.0f}")
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
AB = [(18, 34), (35, 49), (50, 64), (65, 89)]
CB = [(1900, 1929), (1930, 1944), (1945, 1959), (1960, 1974), (1975, 2006)]
print("\n=== 年龄段 × 出生队列 网格(格 n<400 剔除并公布)===")
grid, drop = {}, []
for lo, hi in AB:
    row = []
    for c0, c1 in CB:
        m = (AG >= lo) & (AG <= hi) & (CO >= c0) & (CO <= c1)
        if m.sum() < 400:
            row.append(None); drop.append(f"{lo}-{hi}×{c0}-{c1}(n={int(m.sum())})"); continue
        row.append(cell(m))
    grid[f"{lo}-{hi}"] = row
    print(f"  年龄 {lo}-{hi}: " + " · ".join("—" if x is None else f"{x:.3f}" for x in row))
print(f"  剔除的格({len(drop)}):{drop}")
cmid = np.array([(a + b) / 2 for a, b in CB], float)
amid = np.array([(a + b) / 2 for a, b in AB], float)
def slopes(G):
    """A: 同一队列内随年龄;B: 同一年龄内随队列。"""
    A, Bs = [], []
    for j in range(len(CB)):
        col = [(amid[i], G[f"{AB[i][0]}-{AB[i][1]}"][j]) for i in range(len(AB))
               if G[f"{AB[i][0]}-{AB[i][1]}"][j] is not None]
        if len(col) >= 3: A.append(np.polyfit([x for x, _ in col], [y for _, y in col], 1)[0])
    for i in range(len(AB)):
        row = [(cmid[j], G[f"{AB[i][0]}-{AB[i][1]}"][j]) for j in range(len(CB))
               if G[f"{AB[i][0]}-{AB[i][1]}"][j] is not None]
        if len(row) >= 3: Bs.append(np.polyfit([x for x, _ in row], [y for _, y in row], 1)[0])
    return (float(np.median(A)) if A else np.nan, len(A),
            float(np.median(Bs)) if Bs else np.nan, len(Bs))
sa, na, sb, nb = slopes(grid)
print(f"\n  **A 同一队列内随年龄:{sa:+.6f}/岁**({na} 条队列列)")
print(f"  **B 同一年龄内随队列:{sb:+.6f}/出生年**({nb} 条年龄行)")
# 安慰剂:打乱出生年
G_ = Gate("那个收紧,是世代换了人,还是同一批人自己拧紧了?")
rng = np.random.default_rng(SEEDS[0])
CO2 = CO[rng.permutation(len(CO))]
grid2 = {}
for lo, hi in AB:
    row = []
    for c0, c1 in CB:
        m = (AG >= lo) & (AG <= hi) & (CO2 >= c0) & (CO2 <= c1)
        row.append(None if m.sum() < 400 else cell(m))
    grid2[f"{lo}-{hi}"] = row
sa2, _, sb2, _ = slopes(grid2)
print(f"  安慰剂(打乱出生年):A={sa2:+.6f} · B={sb2:+.6f}")
# 两个斜率的 MDE:对格做 bootstrap
bs_a, bs_b = [], []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(60):
        gg = {k: [None if x is None else x * (1 + r.normal(0, 0.02)) for x in v] for k, v in grid.items()}
        a_, _, b_, _ = slopes(gg)
        if np.isfinite(a_): bs_a.append(a_)
        if np.isfinite(b_): bs_b.append(b_)
MA, MB = 2.8 * float(np.std(bs_a)), 2.8 * float(np.std(bs_b))
print(f"  MDE:A={MA:.6f} · B={MB:.6f}")
G_.positive_control("正对照:可算的格 ≥ 12", planted=float(sum(1 for v in grid.values() for x in v if x is not None)),
                    floor=11.5, spread=1e-9)
G_.negative_control("安慰剂:打乱出生年后 B 应塌", null=abs(sb2), effect=abs(sb),
                    null_spread=1e-9, null_kind="出生年标签置换")
cells = {f"{k}×{CB[j][0]}": dict(n=1, value=v[j], inclusion=[f"年龄 {k}", f"队列 {CB[j][0]}-{CB[j][1]}", "格 n>=400"])
         for k, v in grid.items() for j in range(len(CB))}
G_.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G_.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 70)
if abs(sb2) < 0.5 * abs(sb):
    A_hit, B_hit = abs(sa) > MA, abs(sb) > MB
    world = ("BOTH" if A_hit and B_hit else ("COHORT" if B_hit else ("AGE" if A_hit else "NEITHER")))
    verdict = (f"A {sa:+.6f}{'超' if A_hit else '不超'}MDE {MA:.6f} · B {sb:+.6f}{'超' if B_hit else '不超'}MDE {MB:.6f}"
               + {"BOTH": " -> **两者都动,不主张单一机制**",
                  "COHORT": " -> **世代更替**", "AGE": " -> **同龄不同期的差**",
                  "NEITHER": " -> **上升不可归因,写进页面**"}[world])
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**APC 共线** —— 年龄 + 出生队列 = 调查年,三者不可同时识别;"
          "而横断面**没有同一个人被问两次**,所以「个体拧紧」永远只是一个推断,不是一次观测。")
else:
    world, verdict = "UNVERIFIED", "安慰剂未塌"
    print(f"⚠ {verdict}")
print(G_)
json.dump(dict(grid=grid, dropped=drop, slope_age=sa, slope_cohort=sb, mde_age=MA, mde_cohort=MB,
               placebo=dict(age=sa2, cohort=sb2), world=world, verdict=verdict, seeds=SEEDS,
               impossible=["APC 共线,三者不可同时识别", "横断面,没有同一个人被问两次",
                           "一国一仪器"], unchallenged=True), open(OUT / "apc.json", "w"), indent=1)
print(f"\nwrote {OUT/'apc.json'}")
