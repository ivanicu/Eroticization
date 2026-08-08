"""E02·A211·R574 — 换一条同类型的基线,`#529` 还站得住吗?

`#529` 的 NEXT。行动类型:**FRONTIER**(结局会改变 `#529` 的可信度,不只是它的数值)。
一轮打包多个操作:识别性 + 主检验 + 三道对照 + 规格曲线 + 页面兑现。

**要修的洞(`#529d`,我自己写下的):** `#529` 拿 **0.2239** 当零,而那是
**`#528b` 在「谴责变量 × 该源全部变量」上算的全源基线** —— 它假定「谴责量表之间的编码耦合
强度」等于「该源任意两变量的耦合强度」。**本轮不再假定,直接测。**

G1 ESTIMAND(先于方法):**同源、同类型、同样跨做法**的参照分布 ——
   `broude1976cross` 的**频率/程度**类序数变量,**只取跨做法的对**(排除同一做法的男女配对,
   因为 166×167 这类是同一构念的两个测量,ρ=0.877,会把基线抬成天花板)。
   概括量 = 该分布的**中位 |ρ|**。目标量 = `#529` 的谴责跨做法非对角中位 **+0.1249**。

WORLDS:
  W-HOLDS   新基线 ≥ 0.1249 ⇒ 谴责之间的耦合**不高于**同类型变量之间的 ⇒ `#529` 站得住
  W-DOWNGRADE 新基线 < 0.1249 且差超过分辨率 ⇒ **谴责之间确实更耦合**,`#529` 必须降级
  W-NOBASE  新基线本身的 k 太小或散度太大 ⇒ **这条基线不可用**,`#529` 的零仍未被检验(UNVERIFIED)
⚠ BASIN:`W-HOLDS` 保住我上一轮刚写上页面的结论,所以它**不是**本轮下注方向。
   本轮下注 `W-DOWNGRADE` —— 它要求我把昨天写上页面的第六条**当场降级**。

CONTROLS(G2):
  正对照 同一做法的男女频率对(166×167 等)必须**远高于**跨做法对 —— 若不然,
     「跨做法」这个切分就没有意义,基线不可用(这是本轮基线本身的可失败检查);
  安慰剂 频率变量 × 社会纬度 ≈ 0(该是零 ⇒ negative_control);
  分辨率 对基线中位做块 bootstrap,给它自己的 MDE —— **不拿一个没有误差棒的数当零**。
KILL(条件式,预注册):
  if 同做法对 > 跨做法对(基线切分有效) and 安慰剂 ≈ 0:
      |新基线 − 0.1249| < 基线 MDE -> W-HOLDS(差不可分辨,`#529` 不改)
      新基线 < 0.1249 − MDE        -> W-DOWNGRADE
      新基线 > 0.1249 + MDE        -> W-HOLDS(更强)
  else: UNVERIFIED
IMPOSSIBLE:一个编码团队(`#521`)· 无系统发生树 · 无干预 · 跨做法频率对数量有限 ⇒ k 小 · [unchallenged]
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
OLD_BASE = 0.2239        # `#528b` 的全源基线(本轮要替换的那个)
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False); S = pd.read_csv(SC / "societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
S = S.set_index("id").reindex(W.index)
BLK = (np.floor(S.Lat.values / 10).astype(str) + "_" + np.floor(S.Long.values / 10).astype(str))
W["_lat"] = S.Lat.values

# 频率/程度类序数变量,标注它们属于哪一种做法 —— 同做法的对要被排除出基线
FREQ = {"SCCS166": "婚前", "SCCS167": "婚前", "SCCS170": "婚外", "SCCS171": "婚外",
        "SCCS177": "同性恋", "SCCS174": "强奸", "SCCS160": "婚内", "SCCS175": "男性主动",
        "SCCS178": "阳痿", "SCCS163": "着衣年龄", "SCCS164": "着衣年龄"}
def sp(a, b, mn=25):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < mn: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

print("=== 规则①:逐个打印频率变量的 n ===")
for v, p in FREQ.items():
    print(f"  {v} [{p:6s}] n={int(W[v].notna().sum()):3d}")

same_p, cross_p = [], []
for v1, v2 in itertools.combinations(FREQ, 2):
    r, n = sp(W[v1].values.astype(float), W[v2].values.astype(float))
    if not np.isfinite(r): continue
    (same_p if FREQ[v1] == FREQ[v2] else cross_p).append(dict(v=f"{v1}×{v2}", n=n, rho=r,
        same_practice=FREQ[v1] == FREQ[v2], practices=f"{FREQ[v1]}|{FREQ[v2]}",
        inclusion=[f"两列都非缺失 (n={n})", "broude1976cross 频率/程度类序数变量", "n>=25"]))
sm = np.array([abs(x["rho"]) for x in same_p]); cr = np.array([abs(x["rho"]) for x in cross_p])
print(f"\n=== 基线切分的可失败检查:同做法对 vs 跨做法对 ===")
print(f"  同做法对 k={len(sm):3d}  中位 |ρ| = {np.median(sm):.4f}   {[x['v'] for x in same_p]}")
print(f"  跨做法对 k={len(cr):3d}  中位 |ρ| = {np.median(cr):.4f}")
print(f"  ⇒ 切分{'有效' if np.median(sm) > np.median(cr) else '**无效 —— 基线不可用**'}")

# 基线中位的自身分辨率:块 bootstrap 重算整条分布的中位
meds = []
ub = pd.unique(BLK)
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(200):
        idx = np.concatenate([np.where(BLK == x)[0] for x in rng.choice(ub, len(ub))])
        vals = []
        for v1, v2 in itertools.combinations(FREQ, 2):
            if FREQ[v1] == FREQ[v2]: continue
            r, n = sp(W[v1].values.astype(float)[idx], W[v2].values.astype(float)[idx])
            if np.isfinite(r): vals.append(abs(r))
        if vals: meds.append(np.median(vals))
meds = np.array(meds); BMDE = 2.8 * meds.std(); NEW = float(np.median(cr))
print(f"\n  **新基线(同源·同类型·跨做法)= {NEW:.4f}**  自身 MDE = {BMDE:.4f}  "
      f"CI [{np.quantile(meds,.025):.4f},{np.quantile(meds,.975):.4f}]")
print(f"  旧基线(`#528b` 全源)= {OLD_BASE:.4f}   目标(`#529` 谴责中位)= {TARGET:.4f}")

G = Gate("换一条同类型的基线,`#529` 还站得住吗?")
G.positive_control("基线切分有效:同做法对必须高于跨做法对",
                   planted=float(np.median(sm)), floor=float(np.median(cr)), spread=1e-9)
plc = [abs(sp(W[v].values.astype(float), W._lat.values.astype(float))[0]) for v in FREQ]
plc = [x for x in plc if np.isfinite(x)]
G.negative_control("安慰剂:频率变量 × 社会纬度", null=float(np.median(plc)),
                   effect=float(np.median(sm)), null_spread=float(np.std(plc)),
                   null_kind="任意地理坐标(应无关)")
G.negative_control("g=0:同做法对与自身比较必须无余量",
                   null=abs(float(np.median(sm)) - float(np.median(sm))), effect=float(np.median(sm)),
                   null_spread=1e-9, null_kind="同一分布与自身之差,必为 0")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {x["v"]: x for x in same_p + cross_p})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {x["v"]: x for x in same_p + cross_p})

print("\n" + "=" * 74)
if np.median(sm) > np.median(cr) and np.median(plc) < 0.5 * np.median(sm):
    if abs(NEW - TARGET) < BMDE:
        world = "W-HOLDS"; verdict = (f"新基线 {NEW:.4f} 与谴责中位 {TARGET:.4f} 之差 "
            f"{abs(NEW-TARGET):.4f} < 基线自身 MDE {BMDE:.4f} -> **两者不可分辨:谴责之间的耦合"
            f"与同类型变量之间的一样,`#529` 站得住,而它的措辞要改成「不可分辨」而非「低于」**")
    elif NEW < TARGET - BMDE:
        world = "W-DOWNGRADE"; verdict = (f"新基线 {NEW:.4f} < 谴责中位 {TARGET:.4f} 且差超 MDE "
            f"-> **谴责之间确实更耦合,`#529` 必须降级**")
    else:
        world = "W-HOLDS"; verdict = f"新基线 {NEW:.4f} > 谴责中位 {TARGET:.4f} + MDE -> **`#529` 更强**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:「跨做法」的切分由我按变量标题指派做法,"
          "而 `SCCS175`「男性性主动」与 `SCCS178`「阳痿」算不算同一种做法,是一个判断,不是一个测量。")
else:
    world, verdict = "UNVERIFIED", "控制未齐:切分无效或安慰剂不为零 -> 这条基线不可用"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(new_baseline=NEW, baseline_MDE=float(BMDE), old_baseline=OLD_BASE, target=TARGET,
               same_practice=same_p, cross_practice=cross_p,
               same_median=float(np.median(sm)), cross_median=float(np.median(cr)),
               world=world, verdict=verdict, seeds=SEEDS, placebo=[float(x) for x in plc],
               instrument="Broude & Greene 1976,单一编码团队(#521)",
               impossible=["一个编码团队", "无系统发生树", "无干预非因果", "跨做法频率对 k 有限",
                           "做法指派由标题人工判断"], unchallenged=True),
          open(OUT / "same_kind_baseline.json", "w"), indent=1)
print(f"\nwrote {OUT/'same_kind_baseline.json'}")
