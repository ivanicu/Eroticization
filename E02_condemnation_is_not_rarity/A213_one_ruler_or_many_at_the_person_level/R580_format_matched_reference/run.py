"""E02·A213·R580 — 把性题降成二值:那把尺还在吗?(`A213` 收口轮)

`#534` 的 NEXT。行动类型:**FRONTIER**(结局决定 `#534` 保留还是降级)。

**要堵的洞(`#534c`+`#534d`,我自己写的):** 参照题多为**二值**,性题是**四级**;
二值对的秩相关被机械压低,所以 `+0.2736` 有一部分是**量表长度差**,不是道德结构。

⚠ 硬规则 1 已先做,而它改变了修法:**GSS 里几乎不存在第二道「四级·非性·道德判断」题。**
   全库扫描,级数=4 且 n>3000 的只有 `spanking`(体罚)与 `pillok`(给 14–16 岁避孕药)——
   **而 `pillok` 本身是性题**,放进参照会把要测的东西扣掉。
   ⇒ **换参照这条路在这具仪器上不存在。** 改成把**性题降成二值**,向参照的格式靠。

G1 ESTIMAND:`差_二值 = median ρ(性题二值两两) − median ρ(非性二值参照两两)`,逐年,取中位。
   二值切点**预注册两套**:S-A 严格(仅 always wrong 记 1)· S-B 宽(always+almost always 记 1)。
   **两套都报,不挑。**

WORLDS:
  W-SURVIVES 降成二值后差仍 >> 参照散度 ⇒ **量表长度不是解释,`#534` 保留**
  W-FORMAT   差塌到参照散度内 ⇒ **那把尺大半是格式**,`#534` 降级
  W-PARTIAL  差明显缩小但仍在散度之上 ⇒ **报缩小的比例**,把 `#534` 的数从上界收成区间
⚠ BASIN:`W-SURVIVES` 保住我刚写上页面的第八条,**不是**本轮下注方向。本轮下注 `W-FORMAT`。

CONTROLS:正对照 = **堕胎题组本身就是二值**且已达 0.4371 —— 它证明**二值并不封顶**,
   这是本轮最便宜的那个反驳(gauge test);安慰剂 = 二值性题 × 星座;
   参照散度 = 非性参照的年间标准差。
KILL(条件式):if 二值正对照 > 参照中位 and 安慰剂 ≈ 0: 按三分判 else UNVERIFIED
IMPOSSIBLE:格式匹配只能靠**降级性题**,不能靠升级参照 ⇒ 信息被丢弃,差可能被低估 ·
   一国一仪器 · 观察性非因果 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEX = ["homosex", "premarsx", "xmarsex", "teensex"]
REF = {"cappun": [1, 2], "grass": [1, 2], "letdie1": [1, 2], "suicide1": [1, 2], "fepol": [1, 2]}
POS = {k: [1, 2] for k in ["abany", "abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle"]}
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "zodiac"] + SEX + list(REF) + list(POS), convert_categoricals=False)
CUTS = {"S-A 仅 always wrong": [1], "S-B always+almost always": [1, 2]}

def r2(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 200 or np.std(x[m]) == 0 or np.std(y[m]) == 0: return np.nan
    return float(np.corrcoef(rankdata(x[m]), rankdata(y[m]))[0, 1])

def med_pairs(d, spec):
    out = []
    for a, b in itertools.combinations(spec, 2):
        va = d[a].where(d[a].isin(spec[a])).values.astype(float)
        vb = d[b].where(d[b].isin(spec[b])).values.astype(float)
        r = r2(va, vb)
        if np.isfinite(r): out.append(abs(r))
    return (float(np.median(out)), len(out)) if out else (np.nan, 0)

rows = []
print("=== 逐年 × 两套二值切点。性题已降成二值,与参照同格式 ===")
for cname, cut in CUTS.items():
    for y, d in g.groupby("year"):
        dd = d.copy()
        for s in SEX:
            v = dd[s].where(dd[s].isin([1, 2, 3, 4]))
            dd[s + "_b"] = np.where(v.isna(), np.nan, np.isin(v, cut).astype(float))
        sx = {s + "_b": [0.0, 1.0] for s in SEX}
        ms, ks = med_pairs(dd, sx); mr, kr = med_pairs(dd, REF); mp, kp = med_pairs(dd, POS)
        if not (np.isfinite(ms) and np.isfinite(mr)): continue
        rows.append(dict(cut=cname, year=int(y), rho_sex_bin=ms, rho_ref=mr,
                         rho_pos=mp if np.isfinite(mp) else None, diff=ms - mr,
                         k_sex=ks, k_ref=kr, n=int(len(d)),
                         inclusion=[f"{int(y)} 年 n={len(d)}", cname, "性题已二值化,与参照同格式",
                                    "每对 n>=200"]))
for cname in CUTS:
    sub = [r for r in rows if r["cut"] == cname]
    SX = float(np.median([r["rho_sex_bin"] for r in sub])); RF = float(np.median([r["rho_ref"] for r in sub]))
    SP = float(np.std([r["rho_ref"] for r in sub])); DF = float(np.median([r["diff"] for r in sub]))
    print(f"  {cname:26s} {len(sub):2d} 年  性题(二值)={SX:.4f}  参照={RF:.4f}  "
          f"**差={DF:+.4f}**  参照散度={SP:.4f}")
POSM = float(np.median([r["rho_pos"] for r in rows if r["rho_pos"] is not None]))
A = [r for r in rows if r["cut"] == "S-A 仅 always wrong"]
B = [r for r in rows if r["cut"] == "S-B always+almost always"]
DA, DB = float(np.median([r["diff"] for r in A])), float(np.median([r["diff"] for r in B]))
SPREAD = float(np.std([r["rho_ref"] for r in rows]))
FOUR_LEVEL_DIFF = 0.2736      # `#534a`,四级时的差
print(f"\n  四级时的差(`#534a`)= {FOUR_LEVEL_DIFF:+.4f}  ->  二值后 S-A {DA:+.4f} · S-B {DB:+.4f}")
print(f"  保留比例:S-A {DA/FOUR_LEVEL_DIFF:.3f} · S-B {DB/FOUR_LEVEL_DIFF:.3f}   参照散度={SPREAD:.4f}")
print(f"  正对照(堕胎题组,**本身就是二值**)= {POSM:.4f} —— 二值并不封顶,这是最便宜的那个反驳")

G = Gate("把性题降成二值,那把尺还在吗?")
G.positive_control("正对照:二值堕胎题组(证明二值不封顶)", planted=POSM,
                   floor=float(np.median([r["rho_ref"] for r in rows])), spread=1e-9)
zs = []
for s in SEX:
    v = g[s].where(g[s].isin([1, 2, 3, 4]))
    b = np.where(v.isna(), np.nan, np.isin(v, [1]).astype(float))
    z = g.zodiac.where(g.zodiac.between(1, 12)).values.astype(float)
    r = r2(b, z)
    if np.isfinite(r): zs.append(abs(r))
G.negative_control("安慰剂:二值性题 × 星座", null=float(np.median(zs)), effect=POSM,
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['cut'][:3]}|{r['year']}": r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {f"{r['cut'][:3]}|{r['year']}": r for r in rows})
print("\n" + "=" * 76)
if POSM > np.median([r["rho_ref"] for r in rows]) and np.median(zs) < 0.5 * POSM:
    mn = min(DA, DB)
    if mn <= SPREAD:
        world = "W-FORMAT"; verdict = f"二值后差塌到参照散度 {SPREAD:.4f} 内(min {mn:+.4f}) -> **`#534` 降级**"
    elif mn > 0.7 * FOUR_LEVEL_DIFF:
        world = "W-SURVIVES"; verdict = (f"二值后差仍为 S-A {DA:+.4f} / S-B {DB:+.4f},保留 "
            f"{min(DA,DB)/FOUR_LEVEL_DIFF:.0%} 以上 -> **量表长度不是解释,`#534` 保留**")
    else:
        world = "W-PARTIAL"; verdict = (f"二值后差缩到 S-A {DA:+.4f} / S-B {DB:+.4f}"
            f"(保留 {DA/FOUR_LEVEL_DIFF:.0%} / {DB/FOUR_LEVEL_DIFF:.0%})但仍远高于散度 {SPREAD:.4f}"
            f" -> **格式解释掉一部分,不是全部;把 `#534` 的上界收成区间**")
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:二值化**丢掉了信息**,所以二值后的差可能**低估**真实余量;"
          "而参照题的内容(死刑/大麻/安乐死)彼此本就不属于同一道德领域,"
          "**「同领域」与「同格式」在本轮仍然混在一起**,只是格式那一半被处理了。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, diff_S_A=DA, diff_S_B=DB, four_level_diff=FOUR_LEVEL_DIFF,
               retained=[DA / FOUR_LEVEL_DIFF, DB / FOUR_LEVEL_DIFF], ref_spread=SPREAD,
               pos_median=POSM, world=world, verdict=verdict, placebo=[float(x) for x in zs],
               instrument="GSS 个体层,性题二值化以匹配参照格式",
               impossible=["格式匹配只能降级性题,不能升级参照(GSS 无第二道四级非性道德题)",
                           "二值化丢信息,差可能被低估", "同领域与同格式仍混在一起",
                           "一国一仪器", "观察性非因果"], unchallenged=True),
          open(OUT / "format_matched.json", "w"), indent=1)
print(f"\nwrote {OUT/'format_matched.json'}")
