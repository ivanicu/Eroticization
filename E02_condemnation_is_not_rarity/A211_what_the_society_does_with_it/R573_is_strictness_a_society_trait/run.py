"""E02·A211·R573 — 严厉是社会的性格,还是做法的性格?

`#528` 的 NEXT,留在 `A211`,一轮打包多个操作。行动类型:**FRONTIER**。
这是 `#519`(人层「不是一把尺,是两把」)在**社会层**的复制尝试 —— 硬规则 4。

⚠ 硬规则 1 在写代码之前就已经杀掉两个变量,记在这里因为它本身是结果:
   **`SCCS172`「妻共享」不是谴责量表** —— 码是「为任何理由 / 对特定群体 / 为丈夫的经济利益 / …」,
   是**名义类型学**,不是序数严厉度。它本来是三个有力对里的两个(n=81、80),**全部作废**。
   ⇒ 真正序数的谴责/限制量表只有:SCCS165 婚前(女)· SCCS169 婚外 · SCCS173 强奸 ·
     SCCS176 同性恋 · SCCS159 谈性 · SCCS161 性被认为危险。
   ⚠ `SCCS169` 的码是「单一/双重标准」类型学,只是**近似**序数(2 与 3 差在**种类**不在程度)——
     所以它进**规格曲线的一个轴**,而不是被默认接受。

G1 ESTIMAND(先于方法):
   **σ = 序数谴责量表之间的跨社会秩相关矩阵**;概括量 = **非对角中位 ρ**。
   判据不是「> 0」,而是 **「> 同源参照分布的中位」**(RULE-v3,`#528b` 已测同源中位 0.2239)——
   因为同一团队编的任意两个变量本来就相关 0.22,**低于它就什么都没说**。

WORLDS:
  W-TRAIT     非对角中位 ρ 显著高于同源基线 ⇒ **有些社会对什么都严** —— 严厉是社会的性格
  W-SPECIFIC  ≈ 同源基线 ⇒ 除编码者笔迹外无额外耦合 —— 严厉是**做法**的性格
  W-SPLIT     符号混杂 ⇒ **社会层也是两把尺**(`#519` 的跨单位复制)
⚠ BASIN:`W-TRAIT` 是最省事、最像故事的结局,所以它不是本轮下注方向。
   本轮下注 `W-SPECIFIC` —— 它意味着「严厉的社会」这个概念在这份数据里**不存在**。

CONTROLS(G2):
  正对照 SCCS166×SCCS167(同一构念不同性别)必须强正,且置换后不通过(g=0);
  安慰剂 每个谴责量表 × 社会**纬度**必须 ≈ 0(该是零 ⇒ negative_control);
  基线   同源参照分布中位(`#528b`),**不是一个挑出来的零**;
  n_eff  10°×10° 经纬块 bootstrap,3 seed。
KILL(条件式,预注册):
  if 正对照通过 and 置换后不通过 and 安慰剂 ≈ 0:
      非对角中位 > 2×同源中位 且全部同号 -> W-TRAIT
      符号混杂且至少两格超各自 MDE       -> W-SPLIT
      |非对角中位 − 同源中位| < MDE       -> W-SPECIFIC
      否则                               -> UNVERIFIED-by-power
  else: UNVERIFIED
IMPOSSIBLE:一个编码团队 ⇒ 无跨仪器复制(`#521`)· 无系统发生树 ⇒ 非 Galton 校正 ·
  无干预 ⇒ 非因果 · n=26–90,多数格 MDE ≈ 0.4 ⇒ **本轮大概率是 UNVERIFIED-by-power,已预先承认** ·
  未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
SAME_SOURCE_MEDIAN = 0.2239        # `#528b`,同源基线,不是挑出来的零
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False); S = pd.read_csv(SC / "societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
S = S.set_index("id").reindex(W.index)
W["_lat"] = S.Lat.values
BLK = (np.floor(S.Lat.values / 10).astype(str) + "_" + np.floor(S.Long.values / 10).astype(str))
W.loc[W.SCCS176 == 2, "SCCS176"] = np.nan       # 码 2 码本无描述

SCALES = {"婚前(女)": "SCCS165", "婚外": "SCCS169", "强奸": "SCCS173",
          "同性恋": "SCCS176", "谈性": "SCCS159", "性危险": "SCCS161"}
EXCLUDED = {"妻共享 SCCS172": "名义类型学(情形),不是序数严厉度 —— 硬规则 1 在写代码前杀掉"}

def sp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 25: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

def blkboot(a, b, k=300):
    ub = pd.unique(BLK); out = []
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        for _ in range(k):
            idx = np.concatenate([np.where(BLK == x)[0] for x in rng.choice(ub, len(ub))])
            v, _ = sp(a[idx], b[idx])
            if np.isfinite(v): out.append(v)
    return np.array(out)

# ---- 规格曲线的两个轴:SCCS169 是否保留 / 双重标准码是否剔除
SPECS = {"S1 全部量表,169 原码": dict(drop169=False, dropdouble=False),
         "S2 剔除 169 的双重标准码 2,3": dict(drop169=False, dropdouble=True),
         "S3 完全剔除 169": dict(drop169=True, dropdouble=False)}
print("=== 规则①:逐格先打 n。判据是「高于同源基线 0.2239」,不是「> 0」===")
allspec = {}
for sname, cfg in SPECS.items():
    Ws = W.copy()
    if cfg["dropdouble"]: Ws.loc[Ws.SCCS169.isin([2, 3]), "SCCS169"] = np.nan
    names = {k: v for k, v in SCALES.items() if not (cfg["drop169"] and v == "SCCS169")}
    cells = []
    for (n1, v1), (n2, v2) in itertools.combinations(names.items(), 2):
        a, b = Ws[v1].values.astype(float), Ws[v2].values.astype(float)
        r, n = sp(a, b)
        if not np.isfinite(r): continue
        bs = blkboot(a, b); MDE = 2.8 * bs.std() if len(bs) else np.nan
        cells.append(dict(pair=f"{n1}×{n2}", v=f"{v1}×{v2}", n=n, rho=r, MDE=float(MDE),
                          seen=bool(abs(r) > MDE),
                          inclusion=[f"两列都非缺失 (n={n})", sname,
                                     "SCCS176 码 2 已剔除", "块 bootstrap 10°格"]))
    med = float(np.median([c["rho"] for c in cells]))
    allspec[sname] = dict(cells=cells, offdiag_median=med, k=len(cells))
    print(f"\n  --- {sname} --- {len(cells)} 格,非对角中位 ρ = **{med:+.4f}** "
          f"(同源基线 {SAME_SOURCE_MEDIAN:.4f})")
    for c in sorted(cells, key=lambda x: -x["n"]):
        print(f"    {c['pair']:16s} n={c['n']:3d}  ρ={c['rho']:+.4f}  MDE={c['MDE']:.3f}  "
              f"{'超MDE' if c['seen'] else '看不见'}")

base = allspec["S1 全部量表,169 原码"]
meds = [v["offdiag_median"] for v in allspec.values()]
print(f"\n  规格曲线上的非对角中位:{[f'{m:+.4f}' for m in meds]}  "
      f"跨规格极差 = {max(meds)-min(meds):.4f}")
signs = [c["rho"] > 0 for c in base["cells"]]
print(f"  符号:{sum(signs)} 正 / {len(signs)-sum(signs)} 负")

# ---- 对照
G = Gate("严厉是社会的性格,还是做法的性格?(SCCS/Broude 序数谴责量表)")
pc, pcn = sp(W.SCCS166.values.astype(float), W.SCCS167.values.astype(float))
aa, ff = W.SCCS166.values.astype(float), W.SCCS167.values.astype(float)
m = np.isfinite(aa) & np.isfinite(ff); rng = np.random.default_rng(SEEDS[0])
perm = [abs(sp(aa[m][rng.permutation(m.sum())], ff[m])[0]) for _ in range(300)]
G.positive_control("正对照-v3:同一构念不同性别", planted=abs(pc),
                   floor=max(float(np.quantile(perm, .95)), SAME_SOURCE_MEDIAN), spread=1e-9)
G.negative_control("g=0:置换后必须不通过", null=float(np.median(perm)), effect=abs(pc),
                   null_spread=float(np.std(perm)), null_kind="社会层标签置换")
plc = [abs(sp(W[v].values.astype(float), W._lat.values.astype(float))[0]) for v in SCALES.values()]
plc = [x for x in plc if np.isfinite(x)]
G.negative_control("安慰剂:谴责量表 × 社会纬度", null=float(np.median(plc)), effect=abs(pc),
                   null_spread=float(np.std(plc)), null_kind="任意地理坐标(应无关)")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{s}|{c['pair']}": c
                                             for s, v in allspec.items() for c in v["cells"]})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件",
                                     {f"{s}|{c['pair']}": c for s, v in allspec.items() for c in v["cells"]})

print("\n" + "=" * 74)
if abs(pc) > max(np.quantile(perm, .95), SAME_SOURCE_MEDIAN) and np.median(plc) < 0.5 * abs(pc):
    seen = [c for c in base["cells"] if c["seen"]]
    med = base["offdiag_median"]
    mde_med = float(np.median([c["MDE"] for c in base["cells"]]))
    if med > 2 * SAME_SOURCE_MEDIAN and all(signs):
        world = "W-TRAIT"; verdict = f"非对角中位 {med:+.4f} > 2×基线 且全同号 -> **有些社会对什么都严**"
    elif len(seen) >= 2 and not all(c["rho"] > 0 for c in seen):
        world = "W-SPLIT"; verdict = f"{len(seen)} 格超 MDE 且符号混杂 -> **社会层也是两把尺**"
    elif abs(med - SAME_SOURCE_MEDIAN) < mde_med:
        world = "W-SPECIFIC"; verdict = (f"非对角中位 {med:+.4f} 与同源基线 {SAME_SOURCE_MEDIAN:.4f} "
                                         f"之差 {abs(med-SAME_SOURCE_MEDIAN):.4f} < 中位 MDE {mde_med:.4f} "
                                         f"-> **除编码者笔迹外无额外耦合:「严厉的社会」在这份数据里没有证据**")
    else:
        world = "UNVERIFIED"; verdict = f"只有 {len(seen)}/{len(base['cells'])} 格超 MDE -> UNVERIFIED-by-power"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:同源基线 0.2239 是**另一轮、另一组变量**上算的,"
          "把它当本轮的零,假定了两组变量的编码耦合强度相同 —— 这一点本轮没有检验。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(specs=allspec, spec_medians=meds, baseline_same_source=SAME_SOURCE_MEDIAN,
               excluded=EXCLUDED, world=world, verdict=verdict, seeds=SEEDS,
               positive_control=dict(rho=pc, n=pcn, perm_q95=float(np.quantile(perm, .95))),
               placebo=[float(x) for x in plc],
               instrument="Broude & Greene 1976,单一编码团队(#521)",
               impossible=["一个编码团队无跨仪器复制", "无系统发生树非 Galton 校正", "无干预非因果",
                           "n=26-90,多数格 MDE≈0.4"], unchallenged=True),
          open(OUT / "strictness_trait.json", "w"), indent=1)
print(f"\nwrote {OUT/'strictness_trait.json'}")
