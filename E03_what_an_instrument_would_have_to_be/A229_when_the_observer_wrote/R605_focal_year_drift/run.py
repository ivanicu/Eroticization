"""E03·A229·R605 — 那些码,有多少是在记录「观察者何时写的」?

`#559` 的 NEXT。行动类型:**FRONTIER**。
`#559` 的第二条规格是「**第二组编码者**」,而 SCCS 里没有编码者 ID(`references` 在 Broude 变量上**全空**)。
⇒ 换一个**同样要紧、且用现有数据可测**的量:**民族志的焦点年份 `year`**。
若编码受**观察者写作的年代**影响(早期民族志对性的描述方式不同),
则码会与 `year` 相关,**超出该社会实际做法所能解释的程度**。
**这是「仪器带来的方差」的一个下界**,而它不需要第二组编码者。

G1 ESTIMAND(先于方法):对每个变量 `v`,`d(v) = |Spearman(code_v, focal_year)|`,跨社会。
   **主量 = Broude 性变量的 `d` 中位,减去同数据集**非 Broude** 变量的 `d` 中位。**
   ⚠ 参照必须来自**同一个数据集、同一批社会、同样的焦点年份** —— 否则扣掉的不是仪器。

WORLDS:
  W-ERA      性变量的 `d` 显著高于参照 ⇒ **性的编码比别的变量更受写作年代影响**
             ⇒ `#529`/`#530` 的结论必须再带一层降级
  W-SAME     两者不可分辨 ⇒ 性编码**不比其他编码更受年代影响**(这**不等于**没有影响)
  W-LOWER    性变量反而更低
⚠ BASIN:`W-SAME` 省事,**不是**下注方向。本轮下注 `W-ERA`。
CONTROLS(G2):
  正对照 一个**必然与年代相关**的变量(如「有无电力/现代技术」类,先按标题找候选并打印)必须高;
  安慰剂 `code` 与**社会的经度**(与写作年代无关)必须 ≈0;
  参照分布 非 Broude 变量的 `d` 分布(中位与散度),**不是一个挑出来的零**(RULE-v3)。
KILL(条件式):if 正对照 > 参照中位 and 安慰剂 ≈ 0:
   |性中位 − 参照中位| > 参照散度 -> W-ERA / W-LOWER(按符号);否则 W-SAME
   else UNVERIFIED
IMPOSSIBLE:`year` 是**焦点年份**,不是编码年份 ⇒ 它混了「何时被观察」与「何时被编码」·
  无编码者 ID ⇒ **这不是编码者一致性,只是年代关联的下界** · [unchallenged]
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
print(f"=== 焦点年份 `year`:{len(YR)} 个社会,{YR.min():.0f}–{YR.max():.0f},中位 {YR.median():.0f} ===")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
W = W.reindex(YR.index)
yr = YR.values.astype(float)
lon = S.reindex(W.index).Long.values.astype(float)
def d_of(col, ref):
    v = pd.to_numeric(W[col], errors="coerce").values.astype(float)
    m = np.isfinite(v) & np.isfinite(ref)
    if m.sum() < 40 or len(np.unique(v[m])) < 3: return np.nan, int(m.sum())
    return abs(float(np.corrcoef(rankdata(v[m]), rankdata(ref[m]))[0, 1])), int(m.sum())
BR = [v for v in V[V.source == "broude1976cross"].id if v in W.columns]
OTH = [v for v in W.columns if v not in BR]
bd = [(v, *d_of(v, yr)) for v in BR]; bd = [(v, d, n) for v, d, n in bd if np.isfinite(d)]
od = [(v, *d_of(v, yr)) for v in OTH]; od = [(v, d, n) for v, d, n in od if np.isfinite(d)]
BM = float(np.median([d for _, d, _ in bd])); OM = float(np.median([d for _, d, _ in od]))
SP = float(np.std([d for _, d, _ in od]))
print(f"\n  Broude 性变量 k={len(bd)}  |ρ(code, 焦点年)| 中位 = **{BM:.4f}**")
print(f"  同数据集其他变量 k={len(od)}  中位 = **{OM:.4f}**  散度 = {SP:.4f}")
print(f"  **差 = {BM-OM:+.4f}**")
for v, d, n in sorted(bd, key=lambda x: -x[1])[:6]:
    print(f"    {v:9s} d={d:.4f} n={n:3d}  {V[V.id==v].title.iloc[0][:40]}")
G = Gate("那些码,有多少是在记录「观察者何时写的」?(SCCS 焦点年份)")
tech = [v for v in V.id if any(k in str(V[V.id == v].title.iloc[0]).lower()
        for k in ("writing", "money", "urban", "population", "political integration")) and v in W.columns]
pc = [d_of(v, yr)[0] for v in tech[:12]]
pc = [x for x in pc if np.isfinite(x)]
print(f"\n=== 对照 ===\n  正对照(必然与年代/发展相关的变量,k={len(pc)})中位 = {np.median(pc):.4f}")
G.positive_control("正对照:发展类变量与焦点年份", planted=float(np.median(pc)), floor=OM, spread=1e-9)
plc = [d_of(v, lon)[0] for v in BR]; plc = [x for x in plc if np.isfinite(x)]
G.negative_control("安慰剂:性变量 × 社会经度", null=float(np.median(plc)),
                   effect=float(np.median(pc)), null_spread=float(np.std(plc)),
                   null_kind="与写作年代无关的地理坐标")
cells = {v: dict(n=n, d=d, inclusion=[f"{v} 与焦点年份都非缺失 (n={n})", "至少 3 个不同码", "秩相关"])
         for v, d, n in bd}
G.spec_curve_cells_declare_n("规格曲线逐格 n", cells)
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", cells)
print("\n" + "=" * 74)
if np.median(pc) > OM and np.median(plc) < 0.5 * np.median(pc):
    if abs(BM - OM) <= SP:
        world = "W-SAME"; verdict = (f"性中位 {BM:.4f} 与参照中位 {OM:.4f} 之差 {abs(BM-OM):.4f} "
            f"≤ 参照散度 {SP:.4f} -> **性编码不比其他编码更受年代影响**(⚠ 不等于没有影响)")
    elif BM > OM:
        world = "W-ERA"; verdict = f"性中位高出参照 {BM-OM:+.4f} > 散度 {SP:.4f} -> **性的编码更受写作年代影响**"
    else:
        world = "W-LOWER"; verdict = f"性中位低于参照 {BM-OM:+.4f} -> **性编码反而更少受年代影响**"
    print(f"控制齐备 ⇒ 评判。**{world}**:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:`year` 是**焦点年份**(民族志描述的那个时点),"
          "**不是编码年份** —— 它混了「何时被观察」与「何时被编码」;"
          "而真正的社会变迁也会让码随年代变,**本轮无法把这两者分开**。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(broude_median=BM, other_median=OM, other_spread=SP, diff=BM - OM,
               k_broude=len(bd), k_other=len(od), world=world, verdict=verdict,
               top_broude=[[v, d, n] for v, d, n in sorted(bd, key=lambda x: -x[1])[:10]],
               positive_control=float(np.median(pc)), placebo=float(np.median(plc)),
               impossible=["year 是焦点年份不是编码年份", "无编码者 ID,这是年代关联的下界不是编码者一致性",
                           "真实社会变迁与观察者年代效应不可分"], unchallenged=True),
          open(OUT / "focal_year_drift.json", "w"), indent=1)
print(f"\nwrote {OUT/'focal_year_drift.json'}")
