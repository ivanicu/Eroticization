"""E02·A246·R629 — 二值化到底丢了多少信息?

`#584` 的 NEXT。行动类型:**FRONTIER**。
`#535` 的条款:**「格式匹配只能靠降级性题,不能靠升级参照 ⇒ 二值化丢信息」**。
`#580` 已答了一半(两边一起二值化,那条差保留 93%/95%);**没答的一半是「丢了多少」**。

⚠ BASIN(连续第四轮在「我声明的限制不可靠」这个故事上):
   本轮下注的是**不受欢迎的一侧** —— 若二值化丢得多,`#580` 的「保留 93%」就要打折。

G1 ESTIMAND(先于方法):对每个量表(性 4 题 · 家庭 7 题),
   **分半信度(Spearman–Brown 校正),枚举全部可能的对半划分,取中位**;
   分别在**五级**与**二值**(S-A 严切点)下算。**主量 = 二值 / 五级 的比值。**
预注册:
   比值 **> 0.8** ⇒ 二值化丢的信息小,`#535` 这条降级为「已量且小」;
   **≤ 0.8** ⇒ 丢失真实,**`#580` 的「保留 93%」必须带上这个折扣**。
CONTROLS:正对照 = 五级下的信度必须高于二值下的(**若不然,判据不成立,整轮 UNVERIFIED**)·
   安慰剂 = 把量表随机拆成两半后**跨量表**混装(性 2 题 + 家庭 2 题),信度应明显更低 ·
   逐格 n 与划分数全部打印
IMPOSSIBLE:仅女性 · 单一波 · **分半信度衡量的是「同一量表内部一致」,不是「测到了什么」** ——
   一个二值化后仍然内部一致、却丢掉了区分度的量表,这个指标看不出来 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings, re
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
LAY = {}
for line in open(NS / "setup/2011_2013_FemRespSetup.dct", errors="replace"):
    m = pat.search(line)
    if m: LAY[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
SEX = ["samesex", "sxok18", "sxok16", "gayadopt"]
FAM = ["staytog", "chunless", "chsuppor", "okcohab", "marrfail", "chcohab", "prvntdiv"]
ALL = SEX + FAM
buf = {n: [] for n in ALL}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n in ALL:
        s, w, _ = LAY[n]; v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
X5 = {n: np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in ALL}
XB = {n: np.where(np.isnan(X5[n]), np.nan, (X5[n] == 1).astype(float)) for n in ALL}
print("=== 硬规则 1:逐题 n 与二值化后的阳性率 ===")
for n in ALL:
    print(f"  {n:10s} n={int(np.isfinite(X5[n]).sum()):5d}  二值阳性率={np.nanmean(XB[n]):.4f}")

def align(D, items):
    """⚠ 第一版直接把题相加,**没有对齐极性** —— 而 NSFG 的家庭题里
    `okcohab`(「**不该**同居」)· `marrfail` · `chunless` 是反向题干。
    结果:家庭量表二值化后信度反而近乎翻倍(0.245 -> 0.456),正对照当场失败。
    **`#557a` 抓到过同一条(方向必须先定),我在 YRBS 上做了,在 NSFG 上没做。**
    这里按第一主成分载荷符号对齐(与 `#583`/`R583` 同一做法)。"""
    m = np.all([np.isfinite(D[x]) for x in items], 0)
    M = np.column_stack([D[x][m] for x in items])
    M = (M - M.mean(0)) / (M.std(0) + 1e-12)
    _u, _s, vt = np.linalg.svd(M - M.mean(0), full_matrices=False)
    sg = np.sign(vt[0]); sg[sg == 0] = 1
    return {x: D[x] * sg[i] for i, x in enumerate(items)}, sg


def split_half(D, items):
    """枚举全部对半划分,Spearman–Brown 校正后取中位。"""
    out, k = [], len(items)
    for size in [k // 2]:
        for a in itertools.combinations(items, size):
            b = [x for x in items if x not in a]
            if not b: continue
            m = np.all([np.isfinite(D[x]) for x in items], 0)
            if m.sum() < 500: continue
            A = np.sum([D[x][m] for x in a], 0); Bv = np.sum([D[x][m] for x in b], 0)
            if np.std(A) == 0 or np.std(Bv) == 0: continue
            r = float(np.corrcoef(A, Bv)[0, 1])
            sb = 2 * r / (1 + r) if r > -1 else np.nan
            if np.isfinite(sb): out.append(sb)
    return (float(np.median(out)), len(out)) if out else (np.nan, 0)

res = {}
print("\n=== 分半信度(Spearman–Brown,枚举全部对半划分,取中位)===")
for name, items in (("性 4 题", SEX), ("家庭 7 题", FAM)):
    A5, sg5 = align(X5, items); AB, sgb = align(XB, items)
    print(f"    极性(五级)={sg5.astype(int)} · (二值)={sgb.astype(int)}")
    r5, k5 = split_half(A5, items)
    rb, kb = split_half(AB, items)
    res[name] = dict(n_items=len(items), n_splits=k5, rel5=r5, relb=rb,
                     ratio=rb / r5 if r5 else None,
                     inclusion=[f"{len(items)} 题,{k5} 种对半划分", "该量表全部题非缺失的人", "n>=500"])
    print(f"  {name:8s} {k5:2d} 种划分 · 五级={r5:.4f} · 二值={rb:.4f} · **比值={rb/r5:.3f}**")
ratios = [res[k]["ratio"] for k in res]
print(f"\n  **最小比值 = {min(ratios):.3f}**(预注册门槛 0.8)")
G = Gate("二值化到底丢了多少信息?")
G.positive_control("正对照:五级信度必须高于二值(否则判据不成立)",
                   planted=float(np.median([res[k]["rel5"] for k in res])),
                   floor=float(np.median([res[k]["relb"] for k in res])), spread=1e-9)
rng = np.random.default_rng(SEEDS[0])
mixed = [SEX[0], SEX[1], FAM[0], FAM[1]]
_Am, _sg = align(X5, mixed)
rm5, _ = split_half(_Am, mixed)
G.negative_control("安慰剂:跨量表混装(性 2 + 家庭 2)信度应明显更低",
                   null=float(rm5), effect=float(res["性 4 题"]["rel5"]),
                   null_spread=1e-9, null_kind="把两个不同量表的题混成一个假量表")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {k: dict(n=v["n_splits"], **v) for k, v in res.items()})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", res)
print(f"\n  安慰剂(混装 4 题)五级信度 = {rm5:.4f}  vs 性 4 题 {res['性 4 题']['rel5']:.4f}")
print("\n" + "=" * 70)
if np.median([res[k]["rel5"] for k in res]) > np.median([res[k]["relb"] for k in res]) and rm5 < res["性 4 题"]["rel5"]:
    if min(ratios) > 0.8:
        world = "A-SMALL"; verdict = f"最小比值 {min(ratios):.3f} > 0.8 -> **二值化丢的信息小;`#535` 这条降级为「已量且小」**"
    else:
        world = "B-DISCOUNT"; verdict = (f"最小比值 {min(ratios):.3f} ≤ 0.8 -> "
            f"**丢失真实:`#580` 的「保留 93%」必须带上这个折扣**")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**分半信度衡量的是「同一量表内部一致」,不是「测到了什么」** —— "
          "一个二值化后仍然内部一致、却丢掉了区分度的量表,这个指标看不出来。")
else:
    world, verdict = "UNVERIFIED", "控制未齐(五级未高于二值,或混装未更低)"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(scales=res, ratios=ratios, min_ratio=float(min(ratios)), placebo_mixed=float(rm5),
               world=world, verdict=verdict, seeds=SEEDS,
               impossible=["仅女性", "单一波", "分半信度衡量内部一致,不衡量区分度"],
               unchallenged=True), open(OUT / "reliability_ratio.json", "w"), indent=1)
print(f"\nwrote {OUT/'reliability_ratio.json'}")
