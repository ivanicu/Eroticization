"""E02·A211·R575 — 那个指派是我做的判断,它撑得住吗?(`A211` 收口轮)

`#530` 的 NEXT。行动类型:**CLOSURE**(诚实标注:它保护 `#529`/`#530` 的结论,
**不分离任何新世界**)。按 §0 的三分法,这不是 Frontier,不假装是。
它要保护的具体东西:**「严厉附着在做法上,不附着在社会上」在页面上的那一条。**

**要检的洞(`#530d`,我自己写的):** 「跨做法」的指派是我**按变量标题人工判断**的。
若指派一动,基线就动,判决可能翻。**本轮把指派做成三种合理方案,各自重判。**

G1 ESTIMAND:每一种指派方案下,
   `新基线 = 跨做法频率对的中位 |ρ|`,以及判决 `|新基线 − 0.1249| < 基线自身 MDE ?`
   **概括量 = 三种方案的判决是否一致。**

方案(先于结果写死):
  P1 严格按标题(`#530` 用的那个)
  P2 把「男性主动」与「阳痿」合并为**性功能**(它们都不是一种「做法」,而是一种能力/倾向)
  P3 把「着衣年龄」整个剔除(它是身体暴露规范,不是性实践)
判决规则(预注册):**三种方案的世界标签一致 -> `A211` 的决定变安全,写进页面;
  任一翻转 -> 本条降级,并把翻转条件写进页面。**
CONTROLS:每种方案都要自己通过「同做法对 > 跨做法对」的切分检查;安慰剂同 `#530`。
IMPOSSIBLE:CLOSURE 不分离世界 · 一个编码团队 · k=20–22,基线 CI 宽 · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate

SEEDS = [20260805, 7, 991]; TARGET = 0.1249
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False); S = pd.read_csv(SC / "societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
S = S.set_index("id").reindex(W.index)
BLK = (np.floor(S.Lat.values / 10).astype(str) + "_" + np.floor(S.Long.values / 10).astype(str))

BASE = {"SCCS166": "婚前", "SCCS167": "婚前", "SCCS170": "婚外", "SCCS171": "婚外",
        "SCCS177": "同性恋", "SCCS174": "强奸", "SCCS160": "婚内", "SCCS175": "男性主动",
        "SCCS178": "阳痿", "SCCS163": "着衣年龄", "SCCS164": "着衣年龄"}
P2 = dict(BASE); P2["SCCS175"] = P2["SCCS178"] = "性功能"
P3 = {k: v for k, v in BASE.items() if v != "着衣年龄"}
SCHEMES = {"P1 严格按标题": BASE, "P2 主动+阳痿合并为性功能": P2, "P3 剔除着衣年龄": P3}

def sp(a, b, mn=25):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < mn: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

G = Gate("那个「做法」指派是我做的判断,它撑得住吗?")
res = {}
for sname, MAP in SCHEMES.items():
    sm, cr, cells = [], [], []
    for v1, v2 in itertools.combinations(MAP, 2):
        r, n = sp(W[v1].values.astype(float), W[v2].values.astype(float))
        if not np.isfinite(r): continue
        same = MAP[v1] == MAP[v2]
        (sm if same else cr).append(abs(r))
        cells.append(dict(v=f"{v1}×{v2}", n=n, rho=r, same_practice=same,
                          practices=f"{MAP[v1]}|{MAP[v2]}", scheme=sname,
                          inclusion=[f"两列都非缺失 (n={n})", sname, "n>=25"]))
    meds = []
    ub = pd.unique(BLK)
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(200):
            idx = np.concatenate([np.where(BLK == x)[0] for x in rng.choice(ub, len(ub))])
            vals = []
            for v1, v2 in itertools.combinations(MAP, 2):
                if MAP[v1] == MAP[v2]: continue
                r, n = sp(W[v1].values.astype(float)[idx], W[v2].values.astype(float)[idx])
                if np.isfinite(r): vals.append(abs(r))
            if vals: meds.append(np.median(vals))
    meds = np.array(meds); MDE = 2.8 * meds.std(); NEW = float(np.median(cr))
    split_ok = np.median(sm) > NEW
    world = ("W-HOLDS(不可分辨)" if abs(NEW - TARGET) < MDE else
             ("W-DOWNGRADE" if NEW < TARGET - MDE else "W-HOLDS(更强)")) if split_ok else "UNVERIFIED"
    res[sname] = dict(new_baseline=NEW, MDE=float(MDE), k_cross=len(cr), k_same=len(sm),
                      same_median=float(np.median(sm)), split_ok=bool(split_ok), world=world,
                      gap=float(abs(NEW - TARGET)), cells=cells)
    print(f"  {sname:24s} 跨做法 k={len(cr):3d} 基线={NEW:.4f} MDE={MDE:.4f}  "
          f"同做法中位={np.median(sm):.4f}  |差|={abs(NEW-TARGET):.4f}  -> **{world}**")
    G.positive_control(f"切分有效[{sname[:2]}]", planted=float(np.median(sm)), floor=NEW, spread=1e-9)

worlds = {v["world"] for v in res.values()}
print(f"\n  三种方案的世界标签:{sorted(worlds)}")
print(f"  基线跨方案极差 = {max(v['new_baseline'] for v in res.values()) - min(v['new_baseline'] for v in res.values()):.4f}")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{c['scheme'][:2]}|{c['v']}": c
                                             for v in res.values() for c in v["cells"]})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{c['scheme'][:2]}|{c['v']}": c for v in res.values() for c in v["cells"]})
print("\n" + "=" * 74)
consistent = len(worlds) == 1
if consistent:
    verdict = (f"三种指派方案给出**同一个世界** {list(worlds)[0]} -> "
               f"**`A211` 的决定变安全:指派不是结论的承重墙**")
else:
    verdict = f"方案之间**不一致** {sorted(worlds)} -> **本条降级,翻转条件写进页面**"
print(f"CLOSURE 结论:{verdict}")
print("⚠ 这是 CLOSURE,不是 Frontier:它没有分离任何新世界,只检验了一个我做过的判断是否承重。")
print(G)
json.dump(dict(schemes=res, consistent=bool(consistent), worlds=sorted(worlds), verdict=verdict,
               target=TARGET, action_type="CLOSURE", seeds=SEEDS,
               impossible=["CLOSURE 不分离世界", "一个编码团队", "k=20-22 基线 CI 宽"],
               unchallenged=True), open(OUT / "assignment_sensitivity.json", "w"), indent=1)
print(f"\nwrote {OUT/'assignment_sensitivity.json'}")
