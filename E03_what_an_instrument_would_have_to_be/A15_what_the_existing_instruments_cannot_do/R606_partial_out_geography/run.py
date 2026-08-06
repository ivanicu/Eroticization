"""E03·A229·R606 — 把地理偏掉之后,年代还剩下什么?

`#560` 的 NEXT(`#111c` 允许的第一次修)。行动类型:**FRONTIER**。
`#560b`:安慰剂(性变量 × 经度)= 0.0930,是正对照的 51% ⇒ 地域与年代分不开 ⇒ `UNVERIFIED`。
**修法(预注册于 `#560d`):把纬度与经度偏掉,在残差上重算全部三个量。**

G1 ESTIMAND:`d_partial(v) = |corr(resid(rank(code_v) ~ lat, long), resid(rank(year) ~ lat, long))|`。
判据(预注册,三个都必须过):
  ① 安慰剂 `code × 经度` 的偏相关**必须归零**(经度已被偏掉 ⇒ 结构上应为 0,这是**结构自检**);
  ② 正对照(发展类变量)在偏相关下**仍高于参照中位**;
  ③ 三者齐 -> 报性变量与参照之差,按符号判 `W-ERA` / `W-SAME` / `W-LOWER`。
⚠ ① 是**结构自检不是证据**:偏掉经度后它必然为 0,若不为 0 说明**残差算错了**。
   真正的安慰剂改为 **`code × 一个与年代无关的地理派生量`** —— 用**纬度的平方**(已被线性偏掉一次,
   其非线性成分仍在),它**不应**与年代残差相关。
IMPOSSIBLE:同 `#560e` · 偏掉线性地理**不等于**偏掉文化区域(区域是离散的,不是坐标的线性函数)·
  [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
D = pd.read_csv(SC / "data.csv", low_memory=False)
V = pd.read_csv(SC / "variables.csv"); S = pd.read_csv(SC / "societies.csv").set_index("id")
YR = D.dropna(subset=["year"]).groupby("soc_id")["year"].median()
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first").reindex(YR.index)
G_ = S.reindex(W.index)
lat, lon = G_.Lat.values.astype(float), G_.Long.values.astype(float)
yr = YR.values.astype(float)

def resid(y, X, m):
    """秩化后对 X 做最小二乘,返回残差(只在 m 上)。"""
    Y = rankdata(y[m]); A = np.column_stack([np.ones(m.sum())] + [rankdata(x[m]) for x in X])
    b, *_ = np.linalg.lstsq(A, Y, rcond=None)
    return Y - A @ b

def dpart(col, target, ctrl):
    v = pd.to_numeric(W[col], errors="coerce").values.astype(float)
    m = np.isfinite(v) & np.isfinite(target) & np.all([np.isfinite(c) for c in ctrl], 0)
    if m.sum() < 40 or len(np.unique(v[m])) < 3: return np.nan, int(m.sum())
    a = resid(v, ctrl, m); b = resid(target, ctrl, m)
    if np.std(a) == 0 or np.std(b) == 0: return np.nan, int(m.sum())
    return abs(float(np.corrcoef(a, b)[0, 1])), int(m.sum())

CTRL = [lat, lon]
BR = [v for v in V[V.source == "broude1976cross"].id if v in W.columns]
OTH = [v for v in W.columns if v not in BR]
bd = [(v, *dpart(v, yr, CTRL)) for v in BR]; bd = [(v, d, n) for v, d, n in bd if np.isfinite(d)]
od = [(v, *dpart(v, yr, CTRL)) for v in OTH]; od = [(v, d, n) for v, d, n in od if np.isfinite(d)]
BM = float(np.median([d for _, d, _ in bd])); OM = float(np.median([d for _, d, _ in od]))
SP = float(np.std([d for _, d, _ in od]))
print(f"=== 偏掉纬度与经度之后 ===")
print(f"  Broude 性变量 k={len(bd)}  偏相关中位 = **{BM:.4f}**")
print(f"  其余变量     k={len(od)}  中位 = **{OM:.4f}**  散度 = {SP:.4f}   **差 {BM-OM:+.4f}**")
G = Gate("把地理偏掉之后,年代还剩下什么?")
# ① 结构自检:偏掉经度后 code × 经度 必须为 0
sc = [dpart(v, lon, CTRL)[0] for v in BR]; sc = [x for x in sc if np.isfinite(x)]
print(f"\n=== 控制 ===\n  ① 结构自检 code × 经度(偏掉后)中位 = {np.median(sc):.6f}(必须 ≈0)")
# ② 正对照
tech = [v for v in V.id if any(k in str(V[V.id == v].title.iloc[0]).lower()
        for k in ("writing", "money", "urban", "population", "political integration")) and v in W.columns]
pc = [dpart(v, yr, CTRL)[0] for v in tech[:12]]; pc = [x for x in pc if np.isfinite(x)]
print(f"  ② 正对照(发展类 × 年代,偏相关)中位 = {np.median(pc):.4f}(须 > 参照 {OM:.4f})")
# ③ 真安慰剂:纬度的非线性成分
lat2 = lat ** 2
pl = [dpart(v, lat2, CTRL)[0] for v in BR]; pl = [x for x in pl if np.isfinite(x)]
print(f"  ③ 安慰剂(性变量 × 纬度²的非线性残差)中位 = {np.median(pl):.4f}(须 << 正对照)")
G.negative_control("① 结构自检:偏掉后 code×经度 必为 0", null=float(np.median(sc)),
                   effect=float(np.median(pc)), null_spread=1e-9, null_kind="已被偏掉的变量,结构上必为 0")
G.positive_control("② 正对照:发展类变量 × 年代(偏相关)", planted=float(np.median(pc)), floor=OM, spread=1e-9)
G.negative_control("③ 安慰剂:性变量 × 纬度² 残差", null=float(np.median(pl)),
                   effect=float(np.median(pc)), null_spread=float(np.std(pl)),
                   null_kind="地理的非线性成分,与年代无关")
cells = {v: dict(n=n, d=d, inclusion=[f"{v} · 年份 · 经纬度都非缺失 (n={n})", "秩化后对 lat/long 取残差"])
         for v, d, n in bd}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 74)
ok = np.median(sc) < 1e-6 and np.median(pc) > OM and np.median(pl) < 0.5 * np.median(pc)
if ok:
    world = ("W-SAME" if abs(BM - OM) <= SP else ("W-ERA" if BM > OM else "W-LOWER"))
    verdict = f"{world}: 性 {BM:.4f} vs 参照 {OM:.4f} ± {SP:.4f}(差 {BM-OM:+.4f})"
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:**偏掉线性地理 ≠ 偏掉文化区域** —— 区域是离散的,"
          "不是坐标的线性函数;一个沿区域边界跳变的效应,线性偏相关偏不掉。")
else:
    world, verdict = "UNVERIFIED", f"控制未齐 自检={np.median(sc):.2e} 正={np.median(pc):.4f} 安慰={np.median(pl):.4f}"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(broude_median=BM, other_median=OM, other_spread=SP, diff=BM - OM,
               k_broude=len(bd), k_other=len(od), world=world, verdict=verdict,
               ctrl=dict(structural=float(np.median(sc)), positive=float(np.median(pc)),
                         placebo=float(np.median(pl))),
               top=[[v, d, n] for v, d, n in sorted(bd, key=lambda x: -x[1])[:10]],
               impossible=["偏掉线性地理≠偏掉文化区域", "year 是焦点年份不是编码年份",
                           "真实社会变迁与观察者效应不可分"], unchallenged=True),
          open(OUT / "partial_geo.json", "w"), indent=1)
print(f"\nwrote {OUT/'partial_geo.json'}")
