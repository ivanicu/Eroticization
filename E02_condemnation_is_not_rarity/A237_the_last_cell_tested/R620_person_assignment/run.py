"""E02·A237·R620 — 「人」那一格的领域划分,也是我给的吗?

`#575` 的 NEXT ②。行动类型:**CLOSURE**(保护 `#536`,不分离新世界)。
`#575` 撤回了「年代·行为」那一格,理由是**领域划分由我一个人给,而换一种分法符号就翻**。
⇒ **「人」那一格(`#536`:性内 0.425 vs 家内 0.140)的划分同样是我给的** —— 必须同样验。
⚠ BASIN:**它翻转是我不希望的结局**,所以这一步正是该走的那一步。

三种站得住的分法(先于计算写死):
  S1 严格(`#536` 的臂 A:性 = samesex/sxok18/sxok16;家庭 = 七道;gayadopt 剔除)
  S2 `gayadopt` 归性(`#536` 的臂 B)
  S3 `chsuppor`(未婚女性生养孩子)归性 —— **非婚性行为的后果**,这个归法站得住
预注册:域内−跨域(性内、家内、跨)的**跨方案极差 ≥ 0.05** -> **`#536` 也必须降级**;
  **< 0.05** -> 它是四格里**唯一经过指派检验**的一格。
⚠ 先给极差本身做 bootstrap(`#574d`/`#575a` 的教训)。
IMPOSSIBLE:CLOSURE 不分离世界 · 仅女性 · 单一波 · 三种分法仍是我一个人给的 · [unchallenged]
"""
import os, sys, pathlib, json, re, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import rankdata
from lib.gates import Gate
NS = ROOT / "data/external/nsfg"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SEEDS = [20260805, 7, 991]
pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
LAY = {}
for line in open(NS / "setup/2011_2013_FemRespSetup.dct", errors="replace"):
    m = pat.search(line)
    if m: LAY[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
SEX = ["samesex", "sxok18", "sxok16"]
FAM = ["staytog", "chunless", "chsuppor", "okcohab", "marrfail", "chcohab", "prvntdiv"]
BND = ["gayadopt"]
ALL = SEX + FAM + BND
buf = {n: [] for n in ALL}
for line in open(NS / "2011_2013_FemRespData.dat", errors="replace"):
    for n in ALL:
        s, w, _ = LAY[n]; v = line[s:s + w].strip()
        buf[n].append(float(v) if v not in ("", ".") else np.nan)
X = {n: np.where(np.isin(np.array(buf[n]), [1, 2, 3, 4, 5]), buf[n], np.nan) for n in ALL}
print(f"=== 硬规则 1:逐题 n ===")
for n in ALL: print(f"  {n:10s} n={int(np.isfinite(X[n]).sum()):5d}  {LAY[n][2][:48]}")
def rho(a, b, idx=None):
    A_, B_ = (X[a], X[b]) if idx is None else (X[a][idx], X[b][idx])
    m = np.isfinite(A_) & np.isfinite(B_)
    if m.sum() < 200: return np.nan
    return abs(float(np.corrcoef(rankdata(A_[m]), rankdata(B_[m]))[0, 1]))
def gap(sx, fm, idx=None):
    i = [rho(a, b, idx) for a, b in itertools.combinations(sx, 2)]
    f = [rho(a, b, idx) for a, b in itertools.combinations(fm, 2)]
    i = [x for x in i if np.isfinite(x)]; f = [x for x in f if np.isfinite(x)]
    if not i or not f: return np.nan
    return float(np.median(i) - np.median(f))
SCH = {"S1 严格(gayadopt 剔除)": (SEX, FAM),
       "S2 gayadopt 归性": (SEX + BND, FAM),
       "S3 chsuppor 归性": (SEX + ["chsuppor"], [f for f in FAM if f != "chsuppor"])}
g = {}
print("\n=== 三种分法(性内 − 家内)===")
for k, (sx, fm) in SCH.items():
    g[k] = gap(sx, fm)
    si = [rho(a, b) for a, b in itertools.combinations(sx, 2)]
    fi = [rho(a, b) for a, b in itertools.combinations(fm, 2)]
    print(f"  {k:22s} 性内={np.median([x for x in si if np.isfinite(x)]):.4f} "
          f"家内={np.median([x for x in fi if np.isfinite(x)]):.4f}  **差={g[k]:+.4f}**")
rng_ = max(g.values()) - min(g.values())
print(f"\n  **跨方案极差 = {rng_:.4f}**(预注册门槛 0.05)")
N = len(X[ALL[0]]); bs = []
for sd in SEEDS:
    rr = np.random.default_rng(sd)
    for _ in range(150):
        idx = rr.integers(0, N, N)
        gg = [gap(sx, fm, idx) for sx, fm in SCH.values()]
        gg = [x for x in gg if np.isfinite(x)]
        if len(gg) == len(SCH): bs.append(max(gg) - min(gg))
bs = np.array(bs)
print(f"  极差的 bootstrap:中位 {np.median(bs):.4f} · 95% CI [{np.quantile(bs,.025):.4f},{np.quantile(bs,.975):.4f}]")
G = Gate("「人」那一格的领域划分,也是我给的吗?")
cells = {k: dict(n=N, gap=g[k], inclusion=[k, "受访者 bootstrap", "每对 n>=200"]) for k in SCH}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
G.positive_control("三种方案都算得出(非退化)",
                   planted=float(sum(np.isfinite(v) for v in g.values())), floor=2.5, spread=1e-9)
print("\n" + "=" * 72)
sens = rng_ >= 0.05
verdict = (f"极差 {rng_:.4f} {'≥' if sens else '<'} 0.05 -> "
           + ("**`#536` 也必须降级**" if sens else "**划分不承重;这是四格里唯一经过指派检验的一格**"))
print(f"CLOSURE 结论:{verdict}")
lo, hi = np.quantile(bs, [.025, .975])
print(f"⚠ 极差的 95% CI [{lo:.4f},{hi:.4f}] —— "
      f"{'跨过门槛,判决本身不稳' if lo < 0.05 < hi else '整段在门槛一侧,判决稳'}")
print(G)
json.dump(dict(gaps=g, range=rng_, threshold=0.05, sensitive=bool(sens),
               range_boot_median=float(np.median(bs)), range_boot_ci=[float(lo), float(hi)],
               n=int(N), verdict=verdict, action_type="CLOSURE",
               impossible=["CLOSURE 不分离世界", "仅女性", "单一波", "三种分法仍是我一个人给的"],
               unchallenged=True), open(OUT / "person_assignment.json", "w"), indent=1)
print(f"\nwrote {OUT/'person_assignment.json'}")
