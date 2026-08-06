"""E02·A252·R635 — 那个常数(0.19–0.26)是从哪来的?

`#590` 的 NEXT。行动类型:**FRONTIER**。
`#590b`:换哪个第二领域,性态度与它的相关都落在 **0.19–0.26** —— **对内容惊人地不敏感**。
一个对内容不敏感的常数,通常有两个来源,而两个都可测:

  **候选 1「一般保守性」** —— 所有道德题都载在一条总体保守-自由轴上,
     于是任何两组的相关都被这条共同轴撑到同一水平。**可测:偏掉 `polviews` 后常数应塌。**
  **候选 2「共同方法方差 / 极端作答」** —— 同一份问卷同一批人,
     一个与内容无关的地板。**可测:偏掉该受访者的「极端作答倾向」后常数应塌。**

G1 ESTIMAND(先于方法):对每个第二领域 D,
   `跨(D)` = 性题 × D 题的 |ρ| 中位(逐年算,再取年份中位),在三种口径下:
   **RAW** · **偏掉 `polviews`** · **偏掉极端作答倾向**。
   极端作答倾向 = 该受访者在**全部纳入题**上取**极端选项**的比例(与内容无关的作答风格)。
预注册:
   某一口径下三个 `跨` **全部落到安慰剂 q95 以内** ⇒ **该候选成立**;
   **两个口径都不塌** ⇒ **常数另有来源**,写进页面「做不到什么」;
   **部分塌** ⇒ 报每个候选各解释掉多少。
CONTROLS:正对照 = 最紧配对(`abnomore×abpoor`)在每个口径下必须仍最高 ·
   安慰剂 = 每题 × 星座(逐口径重算)· 逐格 n 打印
IMPOSSIBLE:`polviews` 是**自评**,它本身可能是被这些态度决定的 ⇒
   **偏掉它可能偏掉了要解释的东西本身**(过度控制),本轮无法排除 · 一国一仪器 · 横断面 · [unchallenged]
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
D1 = ["cappun", "grass", "letdie1", "suicide1", "fepol"]
D2 = ["abany", "abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle"]
D3 = ["letdie1", "suicide1"]
ALL = SEX + sorted(set(D1 + D2 + D3))
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "zodiac", "polviews"] + ALL, convert_categoricals=False)
pv = g.polviews.where(g.polviews.between(1, 7)).values.astype(float)
print("=== 硬规则 1 ===")
print(f"  polviews n={int(np.isfinite(pv).sum()):6d} 取值 {sorted(set(pv[np.isfinite(pv)].astype(int)))} · "
      f"年 {g.year[np.isfinite(pv)].nunique()} 个 {int(g.year[np.isfinite(pv)].min())}-{int(g.year[np.isfinite(pv)].max())}")
B, EXTR = {}, []
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
    EXTR.append(np.where(v.isna(), np.nan, np.isin(v, [1, 4]).astype(float)))   # 四级题的两端
for c in set(D1 + D2 + D3):
    v = g[c].where(g[c].isin([1, 2])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
EXT = np.nanmean(np.vstack(EXTR), 0)          # 极端作答倾向:只用四级题(二值题无「极端」可言)
print(f"  极端作答倾向 n={int(np.isfinite(EXT).sum()):6d} 均值={np.nanmean(EXT):.4f}(仅用四级性题的两端)")
YR = g.year.values.astype(int)
def resid(x, ctrl, m):
    Y = rankdata(x[m]); A = np.column_stack([np.ones(m.sum()), rankdata(ctrl[m])])
    b, *_ = np.linalg.lstsq(A, Y, rcond=None)
    return Y - A @ b
def med(pairs, ym, ctrl=None):
    o = []
    for a, b in pairs:
        m = ym & np.isfinite(B[a]) & np.isfinite(B[b])
        if ctrl is not None: m = m & np.isfinite(ctrl)
        if m.sum() < 200: continue
        if ctrl is None:
            xa, xb = rankdata(B[a][m]), rankdata(B[b][m])
        else:
            xa, xb = resid(B[a], ctrl, m), resid(B[b], ctrl, m)
        if np.std(xa) == 0 or np.std(xb) == 0: continue
        o.append(abs(float(np.corrcoef(xa, xb)[0, 1])))
    return float(np.median(o)) if o else np.nan
DOMS = {"D1 便利集合": D1, "D2 堕胎题组": D2, "D3 生死题组": D3}
SPECS = {"RAW": None, "偏掉 polviews": pv, "偏掉极端作答": EXT}
res = {}
print("\n=== 三个第二领域 × 三种口径:『跨』的年份中位 ===")
for sp, ctrl in SPECS.items():
    row = {}
    for name, D in DOMS.items():
        vals = []
        for y in sorted(set(YR)):
            ym = YR == y
            if ym.sum() < 500: continue
            v = med([(s, d) for s in SEX for d in D], ym, ctrl)
            if np.isfinite(v): vals.append(v)
        row[name] = float(np.median(vals)) if vals else np.nan
    res[sp] = row
    print(f"  {sp:14s} " + " · ".join(f"{k.split()[0]}={v:.4f}" for k, v in row.items()))
G = Gate("那个常数是从哪来的?")
zs = {}
for sp, ctrl in SPECS.items():
    z = []
    for c in ALL:
        m = np.isfinite(B[c]) & g.zodiac.between(1, 12).values
        if ctrl is not None: m = m & np.isfinite(ctrl)
        if m.sum() < 1000: continue
        if ctrl is None: xa, xb = rankdata(B[c][m]), rankdata(g.zodiac.values[m])
        else: xa, xb = resid(B[c], ctrl, m), resid(g.zodiac.values.astype(float), ctrl, m)
        z.append(abs(float(np.corrcoef(xa, xb)[0, 1])))
    zs[sp] = float(np.quantile(z, .95))
    print(f"  安慰剂 q95[{sp}] = {zs[sp]:.4f}")
pcs = {}
for sp, ctrl in SPECS.items():
    v = [med([("abnomore", "abpoor")], YR == y, ctrl) for y in sorted(set(YR)) if (YR == y).sum() >= 500]
    v = [x for x in v if np.isfinite(x)]
    pcs[sp] = float(np.median(v))
    G.positive_control(f"正对照[{sp}]:最紧配对仍最高", planted=pcs[sp],
                       floor=max(res[sp].values()), spread=1e-9)
G.negative_control("安慰剂:RAW 口径下每题 × 星座", null=zs["RAW"],
                   effect=float(np.median(list(res["RAW"].values()))), null_spread=1e-9,
                   null_kind="无关的个体层标签")
cells = {f"{sp}|{d}": dict(n=len(DOMS[d]), cross=res[sp][d],
                           inclusion=[sp, d, "逐年算再取年份中位", "每对 n>=200"])
         for sp in SPECS for d in DOMS}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 72)
collapsed = {sp: all(v < zs[sp] for v in res[sp].values()) for sp in SPECS if sp != "RAW"}
if all(pcs[sp] > max(res[sp].values()) for sp in SPECS):
    if any(collapsed.values()):
        who = [k for k, v in collapsed.items() if v]
        world = "EXPLAINED"; verdict = f"{who} 口径下三个『跨』全部落进安慰剂 -> **该候选成立**"
    else:
        world = "UNEXPLAINED"
        drop = {sp: 1 - np.median(list(res[sp].values())) / np.median(list(res["RAW"].values()))
                for sp in SPECS if sp != "RAW"}
        world_extra = " · ".join(f"{k} 解释掉 {v:.0%}" for k, v in drop.items())
        verdict = f"两个候选都没让它塌到安慰剂内({world_extra})-> **常数另有来源**"
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**`polviews` 是自评,它本身可能是被这些态度决定的** ——"
          "偏掉它可能偏掉了要解释的东西本身(过度控制),本轮无法排除。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(cross=res, placebo_q95=zs, positive=pcs, world=world, verdict=verdict, seeds=SEEDS,
               impossible=["polviews 是自评,偏掉它可能是过度控制", "一国一仪器", "横断面",
                           "极端作答倾向只用四级性题的两端"],
               unchallenged=True), open(OUT / "partial_out.json", "w"), indent=1)
print(f"\nwrote {OUT/'partial_out.json'}")
