"""E02·A235·R617 — 把那条「不可分辨」定价:基线从 22 对推到最大

`#572` 的 NEXT。行动类型:**FRONTIER**(结局决定 `#529`/`#530` 是升级还是降级)。

**要处理的东西:** `#529` 判「严厉不跨做法耦合」,`#530` 把它精确成
**「与仪器自身的耦合**不可分辨**」** —— 而那条基线只有 **k=22 对**,自身 CI 宽达 **[0.122, 0.334]**。
**一个宽到这个程度的零,判什么都容易判成『不可分辨』。**

⚠ 硬规则 1 已先做:`broude1976cross` 20 个变量里,**取值≥3 且 n≥30** 的有 **18** 个;
   两两配对中 **n≥25 的有 74 对**(`#530` 用了 22 对)。
⚠ **但不能直接用 74** —— 那会把**谴责量表之间的对**算进基线,而它正是被测目标(循环)。
   基线池 = **非目标**的序数变量,且沿用 `#529a` 对 `SCCS172`(名义类型学)的剔除。

G1 ESTIMAND(先于方法):
   **基线 = 同源、非目标、跨做法配对的 |ρ| 中位**;目标 = `#529` 的谴责跨做法非对角中位 **+0.1249**。
   **主量 = 基线的 CI 宽度**(定价)与 **|基线 − 目标| 是否仍小于基线自身 MDE**(判决)。

预注册(写在跑之前):
  **A** 基线 CI 宽度 **≤ 上次的一半(≤0.106)** 且判决仍是「不可分辨」
      -> `#529`/`#530` **从「不可分辨」升级为「已定价的不可分辨」**;
  **B** 判决翻转(|基线 − 目标| > 基线 MDE)-> **`#529`/`#530` 必须降级**;
  **C** CI 宽度没有明显收窄 -> **说明 k 不是瓶颈**,写进页面「做不到什么」。
CONTROLS:切分有效性(同做法对 > 跨做法对,`#530b` 的可失败检查)· 安慰剂(基线变量 × 社会纬度)·
  块 bootstrap(10°×10° 经纬格)给基线中位自身的展布
IMPOSSIBLE:一个编码团队(`#521`)· 无系统发生树 · 无干预 · 基线池的「做法」指派由我从标题读出
  (`#531` 已验它不承重,但那是在 22 对上验的)· [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
SEEDS = [20260805, 7, 991]
TARGET = 0.1249          # `#529` 的谴责跨做法非对角中位
OLD_CI = (0.122, 0.334)  # `#530` 的基线 CI
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False); S = pd.read_csv(SC / "societies.csv").set_index("id")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
G_ = S.reindex(W.index)
BLK = (np.floor(G_.Lat.values / 10).astype(str) + "_" + np.floor(G_.Long.values / 10).astype(str))
# 目标量表(#529 测的那六个谴责/限制量表)与被剔除者
TARGETV = {"SCCS165", "SCCS169", "SCCS173", "SCCS176", "SCCS159", "SCCS161"}
DROP = {"SCCS172", "SCCS177", "SCCS178"}     # 名义类型学 · 取值仅 2
PRACTICE = {"SCCS160": "婚内", "SCCS162": "婚内", "SCCS163": "着衣年龄", "SCCS164": "着衣年龄",
            "SCCS166": "婚前", "SCCS167": "婚前", "SCCS168": "婚前", "SCCS170": "婚外",
            "SCCS171": "婚外", "SCCS174": "强奸", "SCCS175": "男性主动"}
POOL = [v for v in PRACTICE if v in W.columns and v not in TARGETV and v not in DROP]
print(f"=== 基线池 {len(POOL)} 个(已排除 {len(TARGETV)} 个目标量表与 {len(DROP)} 个不合格)===")
for v in POOL: print(f"  {v} [{PRACTICE[v]}]")

def sp(a, b, mn=25):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < mn: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

V = {v: pd.to_numeric(W[v], errors="coerce").values.astype(float) for v in POOL}
same, cross, cells = [], [], []
for a, b in itertools.combinations(POOL, 2):
    r, n = sp(V[a], V[b])
    if not np.isfinite(r): continue
    sm = PRACTICE[a] == PRACTICE[b]
    (same if sm else cross).append(abs(r))
    cells.append(dict(pair=f"{a}×{b}", practices=f"{PRACTICE[a]}|{PRACTICE[b]}", same_practice=sm,
                      rho=abs(r), n=n, inclusion=[f"两列都非缺失 (n={n})", "非目标序数变量", "n>=25"]))
print(f"\n可用配对 {len(cells)}:同做法 {len(same)} · **跨做法 {len(cross)}**(`#530` 是 22)")
NEW = float(np.median(cross))
meds = []
ub = pd.unique(BLK)
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(300):
        idx = np.concatenate([np.where(BLK == x)[0] for x in rng.choice(ub, len(ub))])
        vals = []
        for a, b in itertools.combinations(POOL, 2):
            if PRACTICE[a] == PRACTICE[b]: continue
            r, _n = sp(V[a][idx], V[b][idx])
            if np.isfinite(r): vals.append(abs(r))
        if vals: meds.append(np.median(vals))
meds = np.array(meds); MDE = 2.8 * meds.std()
lo, hi = np.quantile(meds, [.025, .975]); width = float(hi - lo)
print(f"\n  **新基线 = {NEW:.4f}** · 自身 MDE = {MDE:.4f} · CI [{lo:.4f},{hi:.4f}] **宽 {width:.4f}**")
print(f"  旧基线(`#530`)= 0.1874 · CI [{OLD_CI[0]:.3f},{OLD_CI[1]:.3f}] 宽 {OLD_CI[1]-OLD_CI[0]:.4f}")
print(f"  目标(`#529` 谴责中位)= {TARGET:.4f} · |差| = {abs(NEW-TARGET):.4f}")
G = Gate("把那条「不可分辨」定价:基线从 22 对推到最大")
G.positive_control("切分有效:同做法对必须高于跨做法对",
                   planted=float(np.median(same)), floor=NEW, spread=1e-9)
lat = G_.Lat.values.astype(float)
plc = [abs(sp(V[v], lat)[0]) for v in POOL]; plc = [x for x in plc if np.isfinite(x)]
G.negative_control("安慰剂:基线变量 × 社会纬度", null=float(np.median(plc)),
                   effect=float(np.median(same)), null_spread=float(np.std(plc)),
                   null_kind="任意地理坐标(应无关)")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {c["pair"]: c for c in cells})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {c["pair"]: c for c in cells})
print("\n" + "=" * 74)
if np.median(same) > NEW and np.median(plc) < 0.5 * np.median(same):
    indist = abs(NEW - TARGET) < MDE
    halved = width <= (OLD_CI[1] - OLD_CI[0]) / 2
    if not indist:
        world = "B-DOWNGRADE"; verdict = f"|基线−目标| {abs(NEW-TARGET):.4f} > 基线 MDE {MDE:.4f} -> **`#529`/`#530` 必须降级**"
    elif halved:
        world = "A-PRICED"; verdict = (f"CI 宽 {width:.4f} ≤ 旧宽的一半 {(OLD_CI[1]-OLD_CI[0])/2:.4f},"
            f"且仍不可分辨 -> **升级为「已定价的不可分辨」**")
    else:
        world = "C-NOT-K"; verdict = (f"CI 宽 {width:.4f} 未达旧宽的一半 {(OLD_CI[1]-OLD_CI[0])/2:.4f} "
            f"-> **k 不是瓶颈**,写进页面「做不到什么」")
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print(f"  跨做法对 {len(cross)} 个(旧 22),k 增至 **{len(cross)/22:.1f}×**")
    print("⚠ 这个 KILL 会怎样失败:基线池的「做法」指派由我从标题读出;"
          "`#531` 验过它不承重,**但那是在 22 对上验的**,本轮未重验。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(pool=POOL, n_pairs=len(cells), k_cross=len(cross), k_same=len(same),
               baseline=NEW, mde=float(MDE), ci=[float(lo), float(hi)], ci_width=width,
               old_baseline=0.1874, old_ci=list(OLD_CI), target=TARGET,
               world=world, verdict=verdict, cells=cells, seeds=SEEDS,
               impossible=["一个编码团队", "无系统发生树", "无干预", "做法指派未在 74 对上重验"],
               unchallenged=True), open(OUT / "wider_baseline.json", "w"), indent=1)
print(f"\nwrote {OUT/'wider_baseline.json'}")
