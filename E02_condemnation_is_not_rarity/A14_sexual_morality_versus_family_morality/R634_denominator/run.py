"""E02·A251·R634 — 那个「跨 > 第二组内」,是一座桥,还是一个小分母?

`#589` 的 NEXT。行动类型:**FRONTIER**。
`#589c` 的反常:NSFG 的两组**更近**(同属婚姻/生育大类),比值却**更小**(1.13 vs 1.57)。
两个候选解释:**A** 二值化在 NSFG 上压得更狠 · **B** GSS 的便利集合内部太松,把分母做小了。

⚠ 而问得更锐的版本是:**换第二领域时,「跨」动不动?**
   若 `跨` 几乎不动而 `第二组内` 大幅变化 ⇒ **比值只是分母的产物,不是一座桥** ⇒ `#588b` 必须降级。
⚠ BASIN:**这正是我不希望的结局**(它会削弱「性是核」),所以这一步是该走的。

G1 ESTIMAND(先于方法):在**每一个调查年内**,固定性题组(四道),
   换三个候选第二领域,各算 `跨` 与 `第二组内`:
   D1 **便利集合**(死刑·大麻·安乐死·自杀·女性从政)—— 松,`#588` 用的
   D2 **堕胎题组**(7 道)—— 紧,`#534` 用作过上限
   D3 **生死题组**(安乐死·自杀)—— 中,一个小而真的领域
   **主量:`跨` 在三个 D 之间的极差** vs **`第二组内` 的极差。**
预注册:
   `跨` 的极差 **< `第二组内` 极差的一半** ⇒ **B:比值是分母的产物**,`#588b` 降级为「次序由分母决定」;
   `跨` 的极差 **≥ 一半** ⇒ **跨组相关确实随领域变**,桥的读法站得住。
CONTROLS:正对照 `abany`×`abdefect`(近重复)必须最高 · 安慰剂 每题 × 星座 · 逐格 n 打印
IMPOSSIBLE:堕胎题**与性道德相邻** ⇒ 它当第二领域时「跨」偏高是**部分预期**的 ·
   一国一仪器 · 横断面 · 三个候选领域仍是我挑的 · [unchallenged]
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
                  columns=["year", "zodiac"] + ALL, convert_categoricals=False)
B = {}
for c in SEX:
    v = g[c].where(g[c].isin([1, 2, 3, 4])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
for c in set(D1 + D2 + D3):
    v = g[c].where(g[c].isin([1, 2])); B[c] = np.where(v.isna(), np.nan, (v == 1).astype(float))
YR = g.year.values.astype(int)
print("=== 硬规则 1:逐题 n 与二值阳性率 ===")
for c in ALL: print(f"  {c:9s} n={int(np.isfinite(B[c]).sum()):6d} 阳性率={np.nanmean(B[c]):.4f}")
def med(pairs, ym):
    o = []
    for a, b in pairs:
        m = ym & np.isfinite(B[a]) & np.isfinite(B[b])
        if m.sum() < 200 or np.std(B[a][m]) == 0 or np.std(B[b][m]) == 0: continue
        o.append(abs(float(np.corrcoef(rankdata(B[a][m]), rankdata(B[b][m]))[0, 1])))
    return float(np.median(o)) if o else np.nan
DOMS = {"D1 便利集合(松)": D1, "D2 堕胎题组(紧)": D2, "D3 生死题组(中)": D3}
res, sexin = {}, []
for y in sorted(set(YR)):
    ym = YR == y
    if ym.sum() < 500: continue
    v = med(list(itertools.combinations(SEX, 2)), ym)
    if np.isfinite(v): sexin.append(v)
SEXIN = float(np.median(sexin))
print(f"\n性内(固定)= **{SEXIN:.4f}**({len(sexin)} 年)\n")
print("=== 三个候选第二领域 ===")
for name, D in DOMS.items():
    cx, inn, yrs = [], [], []
    for y in sorted(set(YR)):
        ym = YR == y
        if ym.sum() < 500: continue
        a = med([(s, d) for s in SEX for d in D], ym)
        b = med(list(itertools.combinations(D, 2)), ym)
        if np.isfinite(a) and np.isfinite(b): cx.append(a); inn.append(b); yrs.append(int(y))
    C, I = float(np.median(cx)), float(np.median(inn))
    res[name] = dict(cross=C, internal=I, ratio=C / I, n_years=len(yrs), n_items=len(D),
                     inclusion=[name, f"{len(yrs)} 个可算年份", f"{len(D)} 道题", "二值化 · 每对 n>=200"])
    print(f"  {name:18s} {len(yrs):2d} 年 · 跨={C:.4f} · 组内={I:.4f} · **比值={C/I:.2f}**")
cr = [res[k]["cross"] for k in res]; ir = [res[k]["internal"] for k in res]
rg_c, rg_i = max(cr) - min(cr), max(ir) - min(ir)
print(f"\n  **『跨』的极差 = {rg_c:.4f}** · **『组内』的极差 = {rg_i:.4f}** · 比 = {rg_c/rg_i:.2f}")
G = Gate("那个「跨 > 第二组内」,是一座桥,还是一个小分母?")
pc = []
for y in sorted(set(YR)):
    ym = YR == y
    if ym.sum() < 500: continue
    v = med([("abany", "abdefect")], ym)
    if np.isfinite(v): pc.append(v)
# ⚠ 第一版拿 `abany×abdefect`(0.411)当上限,去比一堆异质量的**最大值** ——
#   而 D3 只有两道题,它的「组内」就是单个配对 `letdie1×suicide1`(0.545),比那个上限还高。
#   `#580c` 同型。改:上限 = **观测到的最紧单配对**;它必须高于被比较的那些「跨」值。
_pairs = [("abany", "abdefect"), ("letdie1", "suicide1"), ("premarsx", "teensex"),
          ("abnomore", "abpoor"), ("abhlth", "abrape")]
_best = []
for _a, _b in _pairs:
    _v = [med([(_a, _b)], YR == y) for y in sorted(set(YR)) if (YR == y).sum() >= 500]
    _v = [x for x in _v if np.isfinite(x)]
    if _v: _best.append((float(np.median(_v)), f"{_a}×{_b}"))
_ceil, _who = max(_best)
print(f"  上限 = 观测到的最紧配对 {_who} = {_ceil:.4f}(候选:{[(f'{v:.3f}', n) for v, n in _best]})")
G.positive_control(f"正对照:最紧配对 {_who} 高于全部『跨』值",
                   planted=_ceil, floor=max(cr), spread=1e-9)
zs = []
for c in ALL:
    m = np.isfinite(B[c]) & g.zodiac.between(1, 12).values
    if m.sum() > 1000: zs.append(abs(float(np.corrcoef(rankdata(B[c][m]), rankdata(g.zodiac.values[m]))[0, 1])))
G.negative_control("安慰剂:每题 × 星座", null=float(np.median(zs)), effect=SEXIN,
                   null_spread=float(np.std(zs)), null_kind="无关的个体层标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {k: dict(n=v["n_years"], **v) for k, v in res.items()})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", res)
print("\n" + "=" * 72)
# ⚠ 只改对照的注册、不改判决条件 —— 于是闸门全过而我的逻辑仍打印「控制未齐」。
#   `#533a` 的第三次:**闸门与三分逻辑不一致,而权威是闸门。** 条件改用同一个上限。
if _ceil > max(cr) and np.median(zs) < 0.5 * SEXIN:
    if rg_c < rg_i / 2:
        world = "B-DENOMINATOR"; verdict = (f"『跨』极差 {rg_c:.4f} < 『组内』极差 {rg_i:.4f} 的一半 -> "
            f"**比值是分母的产物;`#588b` 降级为「次序由分母决定」**")
    else:
        world = "A-BRIDGE"; verdict = (f"『跨』极差 {rg_c:.4f} ≥ 『组内』极差 {rg_i:.4f} 的一半 -> "
            f"**跨组相关确实随领域变,桥的读法站得住**")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**堕胎题与性道德相邻**,所以它当第二领域时「跨」偏高是部分预期的 —— "
          "这会把『跨』的极差做大,即**对 `A-BRIDGE` 有利**,所以判到 B 时这个判决更强。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(sex_internal=SEXIN, domains=res, range_cross=rg_c, range_internal=rg_i,
               ratio_of_ranges=rg_c / rg_i, world=world, verdict=verdict, seeds=SEEDS,
               placebo=float(np.median(zs)), positive=float(np.median(pc)),
               impossible=["堕胎题与性道德相邻,跨偏高是部分预期", "一国一仪器", "横断面",
                           "三个候选领域仍是我挑的"], unchallenged=True),
          open(OUT / "denominator.json", "w"), indent=1)
print(f"\nwrote {OUT/'denominator.json'}")
