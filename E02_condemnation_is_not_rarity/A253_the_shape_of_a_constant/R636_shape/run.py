"""E02·A253·R636 — 那个常数有形状吗?

`#591` 的 NEXT。行动类型:**FRONTIER**。
`#591b` 留下一个说不出来源的连接(扣掉政治与作答习惯后仍有 0.15–0.21)——
**它现在只有一个数,没有形状。** 两个可测的形状:

  ① **它是不是被少数几道题扛着的?** 逐题剔除(leave-one-out),报最大跌幅。
     预注册:某一题剔除让常数跌 **>30%** ⇒ **它不是一个「面对面」的连接,是一两道题的事**。
  ② **它随年代变吗?** 逐年画常数,报斜率与 MDE。
     预注册:斜率 **超 MDE** ⇒ 它是一个**历史现象**,不是结构常数。

G1 ESTIMAND:`常数 = median over {D1,D2,D3} of (性题 × D 题 的 |ρ| 中位)`,逐年算再取年份中位。
CONTROLS:正对照 = **剔除一道 D 题**(第二领域内部)对常数的影响应**小于**剔除一道性题 ——
  若不然,「性侧扛着」这个说法就不成立 · 安慰剂 = 每题 × 星座 · 逐格 n 打印
IMPOSSIBLE:一国一仪器 · 横断面 · 三个第二领域仍是我挑的 ·
  **逐年常数的年数不等**(各题年份覆盖不同)⇒ 斜率对年份加权敏感 · [unchallenged]
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
                  columns=["year", "zodiac"] + ALL, convert_categoricals=False)
B = {}
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
for c in set(x for v in DOMS.values() for x in v):
    v = g[c].where(g[c].isin([1, 2])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
YR = g.year.values.astype(int)
YEARS = [y for y in sorted(set(YR)) if (YR == y).sum() >= 500]
print(f"=== 硬规则 1:{len(YEARS)} 个可算年份 {YEARS[0]}–{YEARS[-1]} ===")
def med(pairs, ym):
    o = []
    for a, b in pairs:
        m = ym & np.isfinite(B[a]) & np.isfinite(B[b])
        if m.sum() < 200 or np.std(B[a][m]) == 0 or np.std(B[b][m]) == 0: continue
        o.append(abs(float(np.corrcoef(rankdata(B[a][m]), rankdata(B[b][m]))[0, 1])))
    return float(np.median(o)) if o else np.nan
def constant(sex=SEX, doms=DOMS, years=None):
    per = []
    for y in (years or YEARS):
        ym = YR == y
        vals = [med([(s, d) for s in sex for d in D], ym) for D in doms.values()]
        vals = [x for x in vals if np.isfinite(x)]
        if vals: per.append(float(np.median(vals)))
    return (float(np.median(per)) if per else np.nan), per
C0, per0 = constant()
print(f"\n常数(全题)= **{C0:.4f}**({len(per0)} 年)")
print("\n=== ① 逐题剔除 ===")
loo = {}
for q in SEX:
    c, _ = constant(sex=[x for x in SEX if x != q])
    loo[f"性:{q}"] = dict(c=c, drop=(C0 - c) / C0, side="sex",
                          inclusion=[f"剔除 {q}", f"{len(per0)} 年", "常数 = 三领域中位"])
    print(f"  剔除性题 {q:9s} -> {c:.4f}  跌幅 {(C0-c)/C0:+.1%}")
for D, items in DOMS.items():
    for q in items:
        d2 = {k: [x for x in v if x != q] or v for k, v in DOMS.items()}
        if any(len(v) == 0 for v in d2.values()): continue
        c, _ = constant(doms=d2)
        loo[f"{D}:{q}"] = dict(c=c, drop=(C0 - c) / C0, side="dom",
                               inclusion=[f"从全部领域剔除 {q}", f"{len(per0)} 年", "常数 = 三领域中位"])
sexdrop = max(abs(v["drop"]) for k, v in loo.items() if v["side"] == "sex")
domdrop = max(abs(v["drop"]) for k, v in loo.items() if v["side"] == "dom")
print(f"  领域题剔除的最大跌幅 = {domdrop:.1%}(逐题见产物)")
print(f"\n  **性题最大跌幅 {sexdrop:.1%} · 领域题最大跌幅 {domdrop:.1%}**(预注册门槛 30%)")
print("\n=== ② 逐年剖面 ===")
yy = np.array([y for y in YEARS][:len(per0)], dtype=float)
pp = np.array(per0)
slope = float(np.polyfit(yy, pp, 1)[0])
bs = []
for sd in SEEDS:
    r = np.random.default_rng(sd)
    for _ in range(2000):
        i = r.integers(0, len(yy), len(yy))
        if np.ptp(yy[i]) > 0: bs.append(np.polyfit(yy[i], pp[i], 1)[0])
MDE = 2.8 * float(np.std(bs))
print(f"  斜率 = {slope:+.6f}/年 · MDE = {MDE:.6f} · {'**超 MDE**' if abs(slope) > MDE else '看不见'}")
print(f"  首末:{int(yy[0])}={pp[0]:.4f} … {int(yy[-1])}={pp[-1]:.4f}")
G = Gate("那个常数有形状吗?")
G.positive_control("正对照:剔除一道性题的影响应大于剔除一道领域题",
                   planted=float(sexdrop), floor=float(domdrop), spread=1e-9)
zs = []
for c in ALL:
    m = np.isfinite(B[c]) & g.zodiac.between(1, 12).values
    if m.sum() > 1000: zs.append(abs(float(np.corrcoef(rankdata(B[c][m]), rankdata(g.zodiac.values[m]))[0, 1])))
G.negative_control("安慰剂:每题 × 星座", null=float(np.median(zs)), effect=C0,
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {k: dict(n=len(per0), **v) for k, v in loo.items()})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", loo)
print("\n" + "=" * 70)
if sexdrop > domdrop and np.median(zs) < 0.5 * C0:
    carried = sexdrop > 0.30
    hist = abs(slope) > MDE
    world = ("CARRIED-BY-FEW" if carried else "DISTRIBUTED") + ("+HISTORICAL" if hist else "+STABLE")
    verdict = (f"最大单题跌幅 {sexdrop:.1%} {'>' if carried else '≤'} 30% · "
               f"斜率 {slope:+.6f} {'超' if hist else '不超'} MDE {MDE:.6f}")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:逐年常数的**年数不等**(各题年份覆盖不同),"
          "斜率对年份加权敏感;而剔除一道题同时改变了配对数与人群。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(constant=C0, loo=loo, max_drop_sex=float(sexdrop), max_drop_dom=float(domdrop),
               slope=slope, mde=float(MDE), per_year=[float(x) for x in pp],
               years=[int(x) for x in yy], world=world, verdict=verdict, seeds=SEEDS,
               impossible=["一国一仪器", "横断面", "三个第二领域仍是我挑的",
                           "逐年常数年数不等,斜率对年份加权敏感"],
               unchallenged=True), open(OUT / "shape.json", "w"), indent=1)
print(f"\nwrote {OUT/'shape.json'}")
