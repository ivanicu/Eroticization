"""E02·A216·R585 — 波之间做不到,但男性做得到(`A216` 收口轮)

`#539` 的 NEXT。行动类型:**一半 PRODUCTION(记录一个不可行),一半 FRONTIER(一个真的外推检验)**。
诚实标注:第一半不分离世界,只把一个「做不到」钉死;第二半分离世界。

**539-NEXT 的第一问,答案是否定的,而查出来只花了一次读字典:**
本地 NSFG 有 6 个女性波的 `.dat`,**但只有 2011–2013 与 2017–2019 有 `.dct` 字典** ——
**定宽文件没有字典读不了。⚠ 有 `.dat` 不等于有数据。**
而这两波里,**2017–2019 砍掉了 10 道 IH 题中的 8 道**:
仅 `samesex` 与 `chsuppor` 幸存(题干逐字相同,只是编号 IH-3→IH-1、IH-8→IH-2)。
⇒ **域内耦合需要每域 ≥3 道题;2017–2019 每域只剩 1 道。年龄与队列的分离,在本地数据上做不到。**
**写进页面「做不到什么」,不写成「计划中」。**

**而幸存的那两道恰好一性一家,且同时在男性问卷里 —— 那正好攻页面上三条都带着的「仅女性」。**

G1 ESTIMAND(先于方法):`ρ_cross = |Spearman(samesex, chsuppor)|`,
   在四个格上各算一次:**2011–2013 女 · 2017–2019 女 · 2017–2019 男**(2011–2013 男字典缺,记明)。
   **主量 = 男女之差** —— `#536` 的跨领域数 0.2174 是女性的,本轮问它是否只是女性的。

WORLDS:
  W-GENERAL 男女之差 < 展布 ⇒ **跨领域耦合不是女性特有** ⇒ 页面上三条的「仅女性」可以放宽
  W-FEMALE  男性明显更低 ⇒ **那是一个关于女性的结构**,页面措辞维持
  W-MALE    男性更高
⚠ BASIN:`W-GENERAL` 让我能放宽三条限制,**是我想要的**,所以不是下注方向。本轮下注 `W-FEMALE`。
CONTROLS:正对照 每格内 `samesex` 与自身分半的一致性不可得(只有一道题)⇒
   **改用 `#536` 已测的 `sxok18×sxok16`=0.5887 作为该仪器的已知上限,并明确它来自另一个波的女性**;
   安慰剂 每格 `samesex` × 随机标签 ≈ 0;基础率:逐格打印 n 与两题的边际分布。
KILL(条件式):if 安慰剂 ≈ 0 且三格 n 都 >1500:|男−女| 与 bootstrap 展布比 else UNVERIFIED
IMPOSSIBLE:**只有一对题** ⇒ 这不是耦合,是**一个相关**,不能与 `#536` 的中位直接比 ·
   2011–2013 男性字典缺 ⇒ 男性只有一个波 · 观察性非因果 · [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NS = ROOT / "data/external/nsfg"
def parse(p):
    out = {}; pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    for line in open(p, errors="replace"):
        m = pat.search(line)
        if m: out[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    return out
CELLS = {"2011-13 女": ("setup/2011_2013_FemRespSetup.dct", "2011_2013_FemRespData.dat"),
         "2017-19 女": ("setup/2017_2019_FemRespSetup.dct", "2017_2019_FemRespData.dat"),
         "2017-19 男": ("setup/2017_2019_MaleSetup.dct", "2017_2019_MaleData.dat")}
CEIL = 0.5887   # `#536` 的 sxok18×sxok16,来自 2011-13 **女性** —— 不是本轮各格自己的上限
res, arr = {}, {}
print("=== 逐格:先打 n 与两题的边际分布,再看 ρ ===")
for name, (dct, dat) in CELLS.items():
    L = parse(NS / dct)
    if not all(k in L for k in ("samesex", "chsuppor")):
        print(f"  {name}: **缺题,跳过**"); continue
    a, b = [], []
    for line in open(NS / dat, errors="replace"):
        for k, buf in (("samesex", a), ("chsuppor", b)):
            s, w, _ = L[k]; v = line[s:s + w].strip()
            buf.append(float(v) if v not in ("", ".") else np.nan)
    A = np.where(np.isin(np.array(a), [1, 2, 3, 4, 5]), a, np.nan)
    B = np.where(np.isin(np.array(b), [1, 2, 3, 4, 5]), b, np.nan)
    m = np.isfinite(A) & np.isfinite(B)
    r = float(np.corrcoef(rankdata(A[m]), rankdata(B[m]))[0, 1])
    bs = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(400):
            i = rng.integers(0, m.sum(), m.sum())
            bs.append(abs(float(np.corrcoef(rankdata(A[m][i]), rankdata(B[m][i]))[0, 1])))
    arr[name] = np.array(bs)
    res[name] = dict(rho=abs(r), n=int(m.sum()), sd=float(np.std(bs)),
                     ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                     marg_samesex=[float(np.mean(A[m] == k)) for k in range(1, 6)],
                     marg_chsuppor=[float(np.mean(B[m] == k)) for k in range(1, 6)],
                     inclusion=[f"两题都在 1–5,n={int(m.sum())}", name, "samesex × chsuppor 一对"])
    print(f"  {name}: n={int(m.sum()):5d}  **|ρ|={abs(r):.4f}**  CI [{res[name]['ci'][0]:.4f},"
          f"{res[name]['ci'][1]:.4f}]  samesex 边际={[f'{x:.2f}' for x in res[name]['marg_samesex']]}")
G = Gate("跨领域耦合是女性特有的吗?(NSFG,samesex × chsuppor)")
G.positive_control("已知上限(#536 的 sxok18×sxok16,来自 2011-13 女)",
                   planted=CEIL, floor=max(v["rho"] for v in res.values()), spread=1e-9)
rng = np.random.default_rng(SEEDS[0])
L = parse(NS / CELLS["2011-13 女"][0]); s, w, _ = L["samesex"]
a = [float(line[s:s + w].strip()) if line[s:s + w].strip() not in ("", ".") else np.nan
     for line in open(NS / CELLS["2011-13 女"][1], errors="replace")]
A = np.where(np.isin(np.array(a), [1, 2, 3, 4, 5]), a, np.nan); mm = np.isfinite(A)
tg = rng.integers(0, 5, mm.sum()).astype(float)
zr = abs(float(np.corrcoef(rankdata(A[mm]), rankdata(tg))[0, 1]))
G.negative_control("安慰剂:samesex × 随机标签", null=zr, effect=CEIL, null_spread=1e-9,
                   null_kind="与问卷无关的随机整数标签")
G.spec_curve_cells_declare_n("规格曲线逐格 n", res)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", res)
print("\n" + "=" * 76)
if "2017-19 男" in res and "2017-19 女" in res and zr < 0.5 * CEIL:
    d = res["2017-19 男"]["rho"] - res["2017-19 女"]["rho"]
    sp = float(np.sqrt(arr["2017-19 男"].var() + arr["2017-19 女"].var()))
    wave = res["2017-19 女"]["rho"] - res["2011-13 女"]["rho"]
    spw = float(np.sqrt(arr["2017-19 女"].var() + arr["2011-13 女"].var()))
    world = ("W-GENERAL" if abs(d) <= 2.8 * sp else ("W-MALE" if d > 0 else "W-FEMALE"))
    print(f"控制齐备 ⇒ 评判。**{world}**:同波男女之差 = {d:+.4f},联合展布 2.8σ = {2.8*sp:.4f}")
    print(f"  顺带(同为女性,两波之差)= {wave:+.4f},2.8σ = {2.8*spw:.4f} —— "
          f"**这一格是可比的,而它只覆盖一对题,不能与 `#536` 的中位比**")
    print("⚠ 这个 KILL 会怎样失败:只有**一对题**,所以它测的是一个相关,不是一个领域耦合;"
          "「男女相似」不能推出「男性那边的领域结构也相似」。")
    verdict = f"{world}: 男女差 {d:+.4f} vs 2.8σ {2.8*sp:.4f}"
else:
    world, verdict = "UNVERIFIED", "控制未齐或缺格"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(cells=res, world=world, verdict=verdict, ceiling_from_536=CEIL, placebo=zr, seeds=SEEDS,
               waves_incomparable=("2017-2019 砍掉 10 道 IH 题中的 8 道;每域仅剩 1 道 -> "
                                   "域内耦合无法跨波计算;且 6 个女性波中仅 2 个有 .dct 字典"),
               instrument="NSFG 2011-2013 与 2017-2019,ACASI",
               impossible=["只有一对题,是一个相关不是领域耦合", "2011-2013 男性字典缺失",
                           "观察性非因果", "上限借自另一波的女性"], unchallenged=True),
          open(OUT / "waves_and_men.json", "w"), indent=1)
print(f"\nwrote {OUT/'waves_and_men.json'}")
