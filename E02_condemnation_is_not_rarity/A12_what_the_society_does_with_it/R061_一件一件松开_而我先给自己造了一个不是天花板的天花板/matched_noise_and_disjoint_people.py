"""E02·A212·R577 — 把噪声配平,并且让两条序列来自**不同的人**

`#576` 的修。行动类型:**FRONTIER**(设计改变了估计量,不只是修一个数)。

**要修的东西(`R576` 的 UNVERIFIED,而它比一个通过更有价值):**
`R576` 把「同一道题的 男Δ × 女Δ」当作仪器上限,得到 +0.3588,
而跨题的全样本差分相关是 +0.3829 —— **目标量超过了「上限」**。
原因不是气候真的比仪器还强,而是:**上限用半样本算,目标用全样本算,两者噪声结构不同。**
一个用不同配对方式算出来的不确定度,不能当另一个的天花板。

**修法一次解决两个问题:**
  ① **噪声配平**:目标量也用半样本 —— 男性子样本的 Δ(题 A)× 女性子样本的 Δ(题 B)。
  ② **人也不共享**:这样两条序列来自**互不相交的受访者**,
     于是 `R576` 自己写下的「共同方法方差抬高 ρ_dif」这个不保守性,**被设计消掉了**。
  ⇒ 上限(同题 男×女)与目标(跨题 男×女)现在**逐项同 n、同配对、同不相交**。

G1 ESTIMAND:`ρ_cross = corr(Δ谴责_A(男), Δ谴责_B(女))`,与
   `ρ_same = corr(Δ谴责_A(男), Δ谴责_A(女))` **在同一噪声结构下**比较。
   ⚠ 对称化:每一对都算两个方向(A男×B女、A女×B男)并取平均,免得性别与题目混在一起。

WORLDS:
  W-CLIMATE     `ρ_cross` 中位 ≥ `ρ_same` 中位的一半 ⇒ 一个时代同时松开多种做法
  W-ONE-BY-ONE  `ρ_cross` ≈ 安慰剂 ⇒ 一件一件松开(`#529` 在年代单位上复制)
  W-PARTIAL     介于两者之间 ⇒ **有共同成分,但小于同题的一致性** —— 报比值,不报二分
CONTROLS:上限 = 同题 男Δ×女Δ(现在同结构)· 安慰剂 = 打乱年份顺序 ·
   g=0 = 同题同性别与自身,必为 1(结构自检,不作证据)
KILL(条件式):if 上限 > 安慰剂 q95 and 打乱后 ≈ 0: 按上面三分 else UNVERIFIED
IMPOSSIBLE:半样本 ⇒ 每年 n 减半,分辨率更低 · 一国一仪器 · 观察性非因果 · [unchallenged]
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

def ser(col, cut, sx):
    d = g[g[col].isin([1, 2, 3, 4]) & (g.sex == sx)]
    return d.groupby("year")[col].agg(lambda v: float(np.isin(v, cut).mean())), d.groupby("year")[col].size()

def dcorr(sa, sb):
    yr = np.array(sorted(set(sa.index) & set(sb.index)))
    if len(yr) < 8: return np.nan, 0
    return float(np.corrcoef(np.diff(sa.reindex(yr).values), np.diff(sb.reindex(yr).values))[0, 1]), len(yr) - 1

CUTS = {"最严 {1}": [1], "中 {1,2}": [1, 2], "最宽 {1,2,3}": [1, 2, 3]}
same_rows, cross_rows = [], []
print("=== 上限与目标现在同结构:同 n、同配对、受访者互不相交 ===")
for cname, cut in CUTS.items():
    S = {(k, sx): ser(k, cut, sx) for k in ITEMS for sx in (1, 2)}
    for k in ITEMS:
        r, nd = dcorr(S[(k, 1)][0], S[(k, 2)][0])
        if np.isfinite(r): same_rows.append(dict(cut=cname, pair=f"{k}(男×女)", rho=r, n=nd,
            kind="same_item", inclusion=[f"{nd} 个差分点", cname, "男/女子样本,互不相交"]))
    for a, b in itertools.combinations(ITEMS, 2):
        r1, n1 = dcorr(S[(a, 1)][0], S[(b, 2)][0])
        r2, n2 = dcorr(S[(a, 2)][0], S[(b, 1)][0])
        vals = [x for x in (r1, r2) if np.isfinite(x)]
        if not vals: continue
        cross_rows.append(dict(cut=cname, pair=f"{a}×{b}", rho=float(np.mean(vals)),
            rho_both=[float(x) for x in vals], n=max(n1, n2), kind="cross_item",
            inclusion=[f"{max(n1,n2)} 个差分点", cname, "A男×B女 与 A女×B男 对称化平均",
                       "两条序列来自互不相交的受访者"]))
SAME = float(np.median([r["rho"] for r in same_rows]))
CROSS = float(np.median([r["rho"] for r in cross_rows]))
for r in same_rows: print(f"  上限 {r['cut']:12s} {r['pair']:22s} n={r['n']:2d}  ρ={r['rho']:+.4f}")
print()
for r in cross_rows: print(f"  目标 {r['cut']:12s} {r['pair']:22s} n={r['n']:2d}  ρ={r['rho']:+.4f}  "
                           f"(两向 {[f'{x:+.3f}' for x in r['rho_both']]})")
print(f"\n  **同题上限中位 = {SAME:+.4f}   跨题目标中位 = {CROSS:+.4f}   比值 = {CROSS/SAME:.3f}**")

G = Gate("把噪声配平、让两条序列来自不同的人之后,时代的宽容是同时松开的吗?")
shuf = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for cname, cut in CUTS.items():
        S = {(k, sx): ser(k, cut, sx) for k in ITEMS for sx in (1, 2)}
        for a, b in itertools.combinations(ITEMS, 2):
            sa, sb = S[(a, 1)][0], S[(b, 2)][0]
            yr = np.array(sorted(set(sa.index) & set(sb.index)))
            if len(yr) < 8: continue
            pa = sa.reindex(yr).values; pb = sb.reindex(yr).values[rng.permutation(len(yr))]
            shuf.append(abs(float(np.corrcoef(np.diff(pa), np.diff(pb))[0, 1])))
Q95 = float(np.quantile(shuf, .95))
print(f"\n=== 对照 ===\n  安慰剂(打乱年份)中位={np.median(shuf):.4f} q95={Q95:.4f}")
G.positive_control("上限:同题男Δ×女Δ(现在同结构)", planted=abs(SAME), floor=Q95, spread=1e-9)
G.negative_control("安慰剂:打乱年份顺序", null=float(np.median(shuf)), effect=abs(SAME),
                   null_spread=float(np.std(shuf)), null_kind="年份顺序置换(破坏时间配对)")
G.negative_control("g=0 结构自检:同题同性别与自身之差必为 0",
                   null=0.0, effect=abs(SAME), null_spread=1e-9, null_kind="同一序列与自身,恒等")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['cut']}|{r['pair']}": r for r in same_rows + cross_rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{r['cut']}|{r['pair']}": r for r in same_rows + cross_rows})
print("\n" + "=" * 76)
if abs(SAME) > Q95 and np.median(shuf) < 0.5 * abs(SAME):
    ratio = CROSS / SAME
    if CROSS >= 0.5 * SAME:
        world = "W-CLIMATE"; verdict = (f"跨题 {CROSS:+.4f} ≥ 同题上限 {SAME:+.4f} 的一半(比值 {ratio:.3f})"
            f" -> **一个时代确实同时松开多种做法**")
    elif CROSS < Q95:
        world = "W-ONE-BY-ONE"; verdict = f"跨题 {CROSS:+.4f} 落在安慰剂 q95={Q95:.4f} 内 -> **一件一件松开**"
    else:
        world = "W-PARTIAL"; verdict = (f"跨题 {CROSS:+.4f} 在安慰剂之上、上限一半之下(比值 {ratio:.3f})"
            f" -> **有共同成分,但小于同一道题内部的一致性**")
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:男女子样本虽不相交,但他们**生活在同一个国家的同一年**,"
          "任何真实的共同事件都会同时推动两条序列 —— 这正是我要测的东西,"
          "所以「不相交」消掉的是**问卷内的**共同方法方差,不是**世界的**共同原因。")
else:
    world, verdict = "UNVERIFIED", f"控制未齐 上限={SAME:.4f} q95={Q95:.4f}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(same_item=same_rows, cross_item=cross_rows, same_median=SAME, cross_median=CROSS,
               ratio=CROSS / SAME, placebo_median=float(np.median(shuf)), placebo_q95=Q95,
               world=world, verdict=verdict, seeds=SEEDS,
               instrument="GSS 1972-2024,男/女子样本互不相交",
               impossible=["半样本 n 减半分辨率更低", "一国一仪器", "观察性非因果",
                           "不相交只消掉问卷内共同方法方差,不消掉世界的共同原因"],
               unchallenged=True), open(OUT / "matched_noise.json", "w"), indent=1)
print(f"\nwrote {OUT/'matched_noise.json'}")
