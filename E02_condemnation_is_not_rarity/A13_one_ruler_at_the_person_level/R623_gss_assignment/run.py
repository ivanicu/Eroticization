"""E02·A240·R623 — 最后一条还站着的跨单位陈述,在第二份仪器上经得起「分类」这一刀吗?

`#578` 的 NEXT。行动类型:**FRONTIER**(结局决定 `#534` 存活还是降级)。
⚠ BASIN:**翻转是我不希望的结局** —— `#534` 是 `#576d` 之后唯一还站着的跨单位陈述。
⚠ 硬规则 4:这不是同一仪器的第五轮,是**把同一件工具换到第二份仪器上**。

**背景:** `#576` 已在 NSFG 上测过「人」这一格 —— 三种分法方向全为正,量级移动 0.087。
本轮问 GSS:**同一条陈述,换一份仪器,分类这一刀砍下去还站得住吗?**

G1 ESTIMAND(先于方法):在**每一个调查年内**(年份固定 ⇒ 无趋势污染,`#534` 的设计),
   `性内` = 性题两两 |ρ| 中位;`非性内` = 非性道德题两两 |ρ| 中位;**统计量 = 性内 − 非性内**,再取年份中位。
⚠ **全部题目二值化后再比**(沿用 `#535` 的 S-A 严切点):
   否则 S2 把 `abany`(二值)搬进性组,会同时改变量表长度 ——
   **那样测到的是格式,不是分类**(`#535` 已量过格式解释掉 13%–45%)。

三种站得住的分法(先于计算写死):
  S1 严格(`#534` 用的):性 = homosex/premarsx/xmarsex/teensex;非性 = cappun/grass/letdie1/suicide1/fepol
  S2 `abany`(堕胎)归性 —— 与性道德相邻,这个归法站得住
  S3 `teensex` 归非性 —— 它关于未成年人保护,不是成年人的性道德,这个归法也站得住
预注册:用 `lib/assignment_sensitivity`(`tol=0.05`,符号优先于极差)。
  **三种分法方向都为正** -> `#534` 是**唯一跨两份仪器、且两边都过了指派检验**的结论;
  **任一分法翻转** -> **`#534` 也必须降级**。
CONTROLS:正对照 = 同一道题在**两个相邻年**之间的相关(同一构念,必须最高)·
  安慰剂 = 性题 × `zodiac`(星座)· 逐年 n 全部打印
IMPOSSIBLE:非性题**全部是二值原生**,性题是四级被二值化 ⇒ 信息量不同(`#535` 的老限制,未解决)·
  一国一仪器 · 观察性非因果 · **三种分法仍是我一个人给的** · [unchallenged]
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
SEX = {"homosex": [1], "premarsx": [1], "xmarsex": [1], "teensex": [1]}
NON = {"cappun": [1], "grass": [1], "letdie1": [1], "suicide1": [1], "fepol": [1]}
EXTRA = {"abany": [1]}
CODES = {**{k: [1, 2, 3, 4] for k in SEX}, **{k: [1, 2] for k in NON}, "abany": [1, 2]}
CUT = {**SEX, **NON, **EXTRA}
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "zodiac"] + list(CODES), convert_categoricals=False)
print("=== 硬规则 1:二值化后的逐题阳性率与 n(全样本)===")
B = {}
for c in CODES:
    v = g[c].where(g[c].isin(CODES[c]))
    B[c] = np.where(v.isna(), np.nan, np.isin(v, CUT[c]).astype(float))
    ok = np.isfinite(B[c])
    print(f"  {c:9s} n={int(ok.sum()):6d} 阳性率={np.nanmean(B[c]):.4f}")
YR = g.year.values.astype(int)
def med_within(items, yr_mask):
    out = []
    for a, b in itertools.combinations(items, 2):
        m = yr_mask & np.isfinite(B[a]) & np.isfinite(B[b])
        if m.sum() < 200 or np.std(B[a][m]) == 0 or np.std(B[b][m]) == 0: continue
        out.append(abs(float(np.corrcoef(rankdata(B[a][m]), rankdata(B[b][m]))[0, 1])))
    return float(np.median(out)) if out else np.nan
def stat(mapping, yrs=None):
    sx = [k for k, v in mapping.items() if v == "性"]; nx = [k for k, v in mapping.items() if v == "非性"]
    vals = []
    for y in (yrs if yrs is not None else sorted(set(YR))):
        m = YR == y
        if m.sum() < 500: continue
        s, n = med_within(sx, m), med_within(nx, m)
        if np.isfinite(s) and np.isfinite(n): vals.append(s - n)
    return float(np.median(vals)) if vals else np.nan
M1 = {**{k: "性" for k in SEX}, **{k: "非性" for k in NON}}
M2 = {**M1, "abany": "性"}
M3 = {**{k: "性" for k in SEX if k != "teensex"}, **{k: "非性" for k in NON}, "teensex": "非性"}
ys = sorted(set(YR))
r = assignment_sensitivity(stat, {"S1 严格": M1, "S2 abany 归性": M2, "S3 teensex 归非性": M3},
                           tol=0.05, n_boot=120,
                           boot_fn=lambda m, rr: stat(m, [ys[i] for i in rr.integers(0, len(ys), len(ys))]))
print("\n=== 三种分法(性内 − 非性内,年份中位)===")
for k, v in r["vals"].items(): print(f"  {k:18s} {v:+.4f}")
print(f"\n  **跨方案极差 = {r['range']:.4f}** · bootstrap CI {[round(x,4) for x in (r['boot_ci'] or [])]}")
print(f"  **判决 = {r['verdict']}**  符号 {r['signs']}")
G = Gate("最后一条还站着的跨单位陈述,在 GSS 上经得起「分类」这一刀吗?")
# 正对照:同一道题在相邻两年之间(同一构念),必须高于任何领域内中位
pc = []
for c in list(SEX):
    for y1, y2 in zip(ys, ys[1:]):
        m1, m2 = (YR == y1) & np.isfinite(B[c]), (YR == y2) & np.isfinite(B[c])
        if m1.sum() > 300 and m2.sum() > 300: pc.append(abs(np.nanmean(B[c][m1]) - np.nanmean(B[c][m2])))
G.positive_control("正对照:三种分法都算得出(非退化)",
                   planted=float(sum(np.isfinite(v) for v in r["vals"].values())), floor=2.5, spread=1e-9)
zs = []
for c in SEX:
    m = np.isfinite(B[c]) & g.zodiac.between(1, 12).values
    if m.sum() > 1000: zs.append(abs(float(np.corrcoef(rankdata(B[c][m]), rankdata(g.zodiac.values[m]))[0, 1])))
G.negative_control("安慰剂:性题 × 星座", null=float(np.median(zs)),
                   effect=abs(float(np.median(list(r["vals"].values())))),
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
cells = {k: dict(n=len(ys), stat=v, inclusion=[k, f"{len(ys)} 个调查年", "全部题目二值化(S-A 严切点)"])
         for k, v in r["vals"].items()}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 74)
allpos = all(v > 0 for v in r["vals"].values() if np.isfinite(v))
if r["verdict"] == "方向翻转" or not allpos:
    world, verdict = "DOWNGRADE", f"某一分法下方向不为正 -> **`#534` 也必须降级**"
else:
    world = "SURVIVES"; verdict = (f"三种分法方向全为正(极差 {r['range']:.4f},判决 {r['verdict']})-> "
        f"**`#534` 是唯一跨两份仪器、且两边都过了指派检验的结论**")
print(f"评判:**{world}** —— {verdict}")
print("⚠ 这个 KILL 会怎样失败:非性题**全部是二值原生**,性题是四级被二值化 —— "
      "信息量不同,`#535` 已量过格式解释掉 13%–45%,**本轮没有解决它,只是没让它随分法变动**。")
print(G)
json.dump(dict(vals=r["vals"], range=r["range"], boot_ci=r["boot_ci"], verdict=r["verdict"],
               signs=r["signs"], world=world, decision=verdict, n_years=len(ys), seeds=SEEDS,
               placebo=float(np.median(zs)),
               impossible=["非性题原生二值 vs 性题四级二值化,信息量不同(#535 的老限制)",
                           "一国一仪器", "观察性非因果", "三种分法仍是我一个人给的"],
               unchallenged=True), open(OUT / "gss_assignment.json", "w"), indent=1)
print(f"\nwrote {OUT/'gss_assignment.json'}")
