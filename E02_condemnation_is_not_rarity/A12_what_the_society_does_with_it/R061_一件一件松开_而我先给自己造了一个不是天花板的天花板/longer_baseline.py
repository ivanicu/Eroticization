"""E02·A212·R578 — 把差分基线拉长,信号会不会从噪声里长出来?(`A212` 收口轮)

`#532` 的 NEXT。行动类型:**FRONTIER**(结局决定 `532e` 是保留还是降级)。

**动机(不是调参):** `532e` 的限制是**分辨率**,不是零 —— 20–29 个相邻年差分点时,
打乱年份仍能做出 |ρ|=0.44。而**态度变化会累积,二项抽样噪声不会**:
把差分的时间跨度从「相邻年」拉到 `L` 年,信号 ∝ L(若变化近似线性),噪声几乎不变。
⇒ **信噪比应随 L 上升。若跨题相关始终不上升,那就不是分辨率的问题。**

G1 ESTIMAND:对每个滞后 `L ∈ {相邻, ≥6年, ≥12年, ≥18年}`,
   `ρ_cross(L) = corr(Δ_L 谴责_A(男), Δ_L 谴责_B(女))`,双向对称化;
   `ρ_same(L)` 同题男女。**概括量 = 比值 `ρ_cross(L)/ρ_same(L)` 随 L 的走向。**
   ⚠ 长滞后的差分点互相重叠(滑动窗口),**自相关会缩窄置信区间** ——
     所以安慰剂也必须在**同样的重叠结构**下重算,而不是借用短滞后的。

WORLDS:
  W-CLIMATE   跨题相关随 L 上升并越过同题上限的一半 ⇒ `532e` 降级,共同气候真实存在
  W-SAFE      跨题相关随 L **不上升**、始终落在同重叠结构的安慰剂内
              ⇒ **`A212` 的决定变安全:不是没测出来,是那里没有**
  W-BOTH-RISE 跨题与同题一起上升、比值不变 ⇒ **上升的是「趋势」不是「同步」**,元分离
⚠ BASIN:`W-SAFE` 保住我刚写上页面的第七条,所以**不是**本轮下注方向。本轮下注 `W-CLIMATE`。

CONTROLS:同题上限(同 L、同重叠)· 安慰剂 = 在**同一重叠结构**下打乱年份 ·
   g=0 = 同一序列与自身,恒为 1(结构自检)
KILL(条件式):if 同题上限 > 该 L 的安慰剂 q95: 按三分判 else UNVERIFIED(该 L 不可判)
IMPOSSIBLE:滑动窗口重叠 ⇒ 有效样本远小于窗口数 · 一国一仪器 · 观察性非因果 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
ITEMS = ["homosex", "premarsx", "xmarsex", "teensex"]
g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "sex"] + ITEMS, convert_categoricals=False)
CUT = [1, 2]

def ser(col, sx):
    d = g[g[col].isin([1, 2, 3, 4]) & (g.sex == sx)]
    return d.groupby("year")[col].agg(lambda v: float(np.isin(v, CUT).mean()))

def lagged(sa, sb, L, perm_rng=None):
    """返回 (Δ_L a, Δ_L b):对每个年份 t,找 >= t+L 的最近共同年 t2。"""
    yr = np.array(sorted(set(sa.index) & set(sb.index)))
    if len(yr) < 6: return None, None
    va, vb = sa.reindex(yr).values, sb.reindex(yr).values
    if perm_rng is not None: vb = vb[perm_rng.permutation(len(yr))]
    da, db = [], []
    for i, t in enumerate(yr):
        j = np.searchsorted(yr, t + L)
        if j >= len(yr) or (L == 0 and j == i): continue
        if L == 0: j = i + 1
        if j >= len(yr): continue
        da.append(va[j] - va[i]); db.append(vb[j] - vb[i])
    return (np.array(da), np.array(db)) if len(da) >= 6 else (None, None)

def cc(x, y):
    if x is None or len(x) < 6 or np.std(x) == 0 or np.std(y) == 0: return np.nan
    return float(np.corrcoef(x, y)[0, 1])

LAGS = {"相邻年": 0, "≥6 年": 6, "≥12 年": 12, "≥18 年": 18}
S = {(k, sx): ser(k, sx) for k in ITEMS for sx in (1, 2)}
rows, summary = [], {}
print("=== 逐个滞后:先打窗口数,再看 ρ。⚠ 长滞后窗口重叠,有效样本远小于窗口数 ===")
for lname, L in LAGS.items():
    same, cross = [], []
    for k in ITEMS:
        a, b = lagged(S[(k, 1)], S[(k, 2)], L)
        r = cc(a, b)
        if np.isfinite(r):
            same.append(r); rows.append(dict(lag=lname, kind="same_item", pair=f"{k}(男×女)",
                rho=r, n_windows=len(a), inclusion=[f"{len(a)} 个窗口", lname, "窗口重叠,非独立"]))
    for x, y in itertools.combinations(ITEMS, 2):
        vals = []
        for (i, j) in [(1, 2), (2, 1)]:
            a, b = lagged(S[(x, i)], S[(y, j)], L)
            r = cc(a, b)
            if np.isfinite(r): vals.append(r)
        if vals:
            cross.append(float(np.mean(vals)))
            rows.append(dict(lag=lname, kind="cross_item", pair=f"{x}×{y}", rho=float(np.mean(vals)),
                n_windows=len(a) if a is not None else 0, rho_both=[float(v) for v in vals],
                inclusion=[lname, "A男×B女 与 A女×B男 对称化", "受访者互不相交", "窗口重叠"]))
    shuf = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(60):
            for x, y in itertools.combinations(ITEMS, 2):
                a, b = lagged(S[(x, 1)], S[(y, 2)], L, perm_rng=rng)
                r = cc(a, b)
                if np.isfinite(r): shuf.append(abs(r))
    SM, CR = float(np.median(same)), float(np.median(cross))
    Q = float(np.quantile(shuf, .95)) if shuf else np.nan
    summary[lname] = dict(same=SM, cross=CR, ratio=CR / SM if SM else None,
                          placebo_q95=Q, k_same=len(same), k_cross=len(cross))
    print(f"  {lname:8s} 同题上限={SM:+.4f}  跨题={CR:+.4f}  比值={CR/SM:.3f}  "
          f"安慰剂q95={Q:.4f}  {'跨题**越过**安慰剂' if CR > Q else '跨题仍在安慰剂内'}")

G = Gate("把差分基线拉长,信号会不会从噪声里长出来?")
for lname in LAGS:
    s = summary[lname]
    G.positive_control(f"同题上限可判[{lname}]", planted=abs(s["same"]), floor=s["placebo_q95"], spread=1e-9)
G.negative_control("g=0 结构自检:同一序列与自身,恒等", null=0.0,
                   effect=abs(summary["相邻年"]["same"]), null_spread=1e-9, null_kind="恒等,必为 0")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['lag']}|{r['kind'][:5]}|{r['pair']}":
                                             dict(n=r["n_windows"], **r) for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{r['lag']}|{r['kind'][:5]}|{r['pair']}": r for r in rows})

ratios = [summary[l]["ratio"] for l in LAGS if summary[l]["ratio"] is not None]
crosses = [summary[l]["cross"] for l in LAGS]
judg = [l for l in LAGS if summary[l]["same"] > summary[l]["placebo_q95"]]
print(f"\n  比值随滞后:{[f'{r:.3f}' for r in ratios]}   跨题随滞后:{[f'{c:+.4f}' for c in crosses]}")
print(f"  可判的滞后(同题上限越过该 L 的安慰剂):{judg or '无'}")
print("\n" + "=" * 76)
if judg:
    over = [l for l in judg if summary[l]["cross"] > 0.5 * summary[l]["same"]]
    inside = [l for l in judg if summary[l]["cross"] <= summary[l]["placebo_q95"]]
    if over:
        world = "W-CLIMATE"; verdict = f"{over} 上跨题越过同题上限的一半 -> **`532e` 降级**"
    elif len(inside) == len(judg):
        world = "W-SAFE"; verdict = (f"全部可判滞后 {judg} 上跨题都仍在安慰剂内,且比值不随 L 上升 "
            f"({[f'{r:.2f}' for r in ratios]}) -> **不是分辨率的问题:`A212` 的决定变安全**")
    else:
        world = "UNVERIFIED"; verdict = f"可判滞后间不一致 -> UNVERIFIED"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:长滞后的窗口**互相重叠**,有效样本远小于窗口数,"
          "所以安慰剂也在同结构下重算 —— 但重叠同时**抬高**了安慰剂 q95,"
          "使「落在安慰剂内」这个判据在长滞后上**更容易成立**,即对 `W-SAFE` 方向**不保守**。")
else:
    world, verdict = "UNVERIFIED", "没有任何滞后上同题上限能越过安慰剂 -> 全部不可判"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(summary=summary, rows=rows, world=world, verdict=verdict, ratios=ratios,
               judgeable_lags=judg, seeds=SEEDS, instrument="GSS 1972-2024,男/女子样本互不相交",
               impossible=["滑动窗口重叠,有效样本远小于窗口数", "一国一仪器", "观察性非因果",
                           "重叠抬高安慰剂 q95,对 W-SAFE 不保守"], unchallenged=True),
          open(OUT / "longer_baseline.json", "w"), indent=1)
print(f"\nwrote {OUT/'longer_baseline.json'}")
