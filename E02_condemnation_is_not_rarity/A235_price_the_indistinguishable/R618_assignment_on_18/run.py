"""E02·A235·R618 — 指派在 18 对上还不承重吗?(`A235` 收口)

`#573` 的 NEXT。行动类型:**CLOSURE**(诚实标注:保护 `#573c` 的结论,不分离新世界)。
`#573d` 的最后一条:基线池的「做法」指派由我从标题读出,`#531` 验过它不承重,
**但那是在 22 对上验的** —— 本轮在**本轮的 18 对**上重验。

预注册(`#573` NEXT ②,写在跑之前):
  基线的跨方案极差 **< 0.02** -> 指派仍不承重,`#573c` 的结论**不带这条保留**;
  **≥ 0.02** -> 结论必须带上「指派敏感」并写进页面。
方案(沿用 `#531` 的形状,适配本轮的池):
  P1 严格按标题 · P2 把「男性主动」并入「强奸」(同为侵犯性) · P3 剔除「着衣年龄」
IMPOSSIBLE:CLOSURE 不分离世界 · 一个编码团队 · 18 对基数小 ⇒ 极差本身也有抽样波动 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False)
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
BASE = {"SCCS160": "婚内", "SCCS162": "婚内", "SCCS163": "着衣年龄", "SCCS164": "着衣年龄",
        "SCCS166": "婚前", "SCCS167": "婚前", "SCCS168": "婚前", "SCCS170": "婚外",
        "SCCS171": "婚外", "SCCS174": "强奸", "SCCS175": "男性主动"}
P2 = dict(BASE); P2["SCCS175"] = "强奸"
P3 = {k: v for k, v in BASE.items() if v != "着衣年龄"}
SCHEMES = {"P1 严格按标题": BASE, "P2 主动并入强奸": P2, "P3 剔除着衣年龄": P3}
V = {v: pd.to_numeric(W[v], errors="coerce").values.astype(float) for v in BASE}
def sp(a, b, mn=25):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < mn: return np.nan
    return abs(float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]))
res, cells = {}, []
print("=== 三种指派方案(全格公布)===")
for name, MAP in SCHEMES.items():
    cross, same = [], []
    for a, b in itertools.combinations(MAP, 2):
        r = sp(V[a], V[b])
        if not np.isfinite(r): continue
        (same if MAP[a] == MAP[b] else cross).append(r)
        cells.append(dict(scheme=name, pair=f"{a}×{b}", same=MAP[a] == MAP[b], rho=r,
                          n=int((np.isfinite(V[a]) & np.isfinite(V[b])).sum()),
                          inclusion=[name, "n>=25", "非目标序数变量"]))
    res[name] = dict(baseline=float(np.median(cross)), k_cross=len(cross), k_same=len(same),
                     same_median=float(np.median(same)) if same else None)
    print(f"  {name:16s} 跨做法 {len(cross):2d} 对 · 基线 = **{np.median(cross):.4f}** · "
          f"同做法 {len(same)} 对 中位 {np.median(same) if same else float('nan'):.4f}")
vals = [res[k]["baseline"] for k in res]
rng_ = max(vals) - min(vals)
print(f"\n  **跨方案极差 = {rng_:.4f}**(预注册门槛 0.02)")
G = Gate("指派在 18 对上还不承重吗?")
for name in res:
    if res[name]["same_median"] is not None:
        G.positive_control(f"切分有效[{name[:2]}]", planted=res[name]["same_median"],
                           floor=res[name]["baseline"], spread=1e-9)
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{c['scheme'][:2]}|{c['pair']}": c for c in cells})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {f"{c['scheme'][:2]}|{c['pair']}": c for c in cells})
print("\n" + "=" * 72)
if rng_ < 0.02:
    world = "不承重"; verdict = f"极差 {rng_:.4f} < 0.02 -> **指派仍不承重,`#573c` 不带这条保留**"
else:
    world = "指派敏感"; verdict = f"极差 {rng_:.4f} >= 0.02 -> **结论必须带上「指派敏感」并写进页面**"
print(f"CLOSURE 结论:{world} —— {verdict}")
print("⚠ 这是 CLOSURE:它没有分离任何新世界,只检验了一个我做过的判断在新基数上是否仍然承重。")
print(G)
json.dump(dict(schemes=res, range=rng_, threshold=0.02, verdict=verdict, world=world,
               action_type="CLOSURE", cells=cells,
               impossible=["CLOSURE 不分离世界", "一个编码团队", "18 对基数小,极差本身有抽样波动"],
               unchallenged=True), open(OUT / "assignment_18.json", "w"), indent=1)
print(f"\nwrote {OUT/'assignment_18.json'}")
