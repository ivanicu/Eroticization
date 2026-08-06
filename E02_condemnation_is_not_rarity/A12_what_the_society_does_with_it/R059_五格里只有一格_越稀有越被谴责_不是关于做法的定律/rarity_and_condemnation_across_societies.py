"""E02·A211·R572 — 换一个单位:在 186 个社会之间,越稀有的做法是否越被谴责?

`#527` 的 NEXT。**方向已改:回到对象。** 一轮打包多个操作(Ivan 2026-08-05)。
行动类型:**FRONTIER**。这是 E01 那个 `corr(稀有, 羞耻) = 0.758` 的**跨仪器、跨单位复制**
(硬规则 4:跨仪器复制优先于同一仪器的第五轮)。单位从**人**换成**社会**。

⚠ 硬规则 1,已先做:20 个 `broude1976cross` 变量各有一行 186,但**实际编码数远少于 186**:
   SCCS165 婚前性态度(女)130 · SCCS167 婚前频率(女)109 · SCCS169 婚外性态度 109 ·
   SCCS171 婚外频率(女)53 · SCCS176 同性恋态度 **40** · SCCS177 同性恋频率 69。
   SCCS176 的码 `2` **在码本里没有描述**(`nan`)—— 本轮把它剔除并声明。

⚠ 硬规则 2,最强混淆,写在跑之前:**态度与频率是同一批人(Broude & Greene 1976)
   从同一份民族志里编出来的。** 一个民族志作者若写「此地此事罕见且遭鄙夷」,
   两个码就同时被那一句话决定。**任何 ρ 都可能是编码者的笔迹,不是社会的事实。**
   本轮**不试图去掉它**,而是**量它**:同源对(Broude×Broude)的 |ρ| 分布
   vs 跨源对(Broude 态度 × Murdock 社会结构)的 |ρ| 分布 —— `#482c` 的参照分布法。

G1 ESTIMAND(先于方法):对每一种做法 p,
   **ρ_p = Spearman(谴责度_p, 稀有度_p),跨社会,中位秩处理并列(`#482a`)。**
   方向已由码本钉死:态度码越大越谴责;频率码 1=Universal…4=Uncommon,**越大越稀有**。
   ⇒ **ρ > 0 就是「越稀有越被谴责」**,即 E01 那条被降级的关系在社会层面的形式。

WORLDS(本体不同):
  W-REPLICATE 三种做法 ρ 都 > 0 且超 MDE ⇒ 这条关系**跨单位、跨仪器存在**
  W-DIES      ρ ≈ 0 或为负 ⇒ **人层面的 0.758 在社会之间不出现**,它是关于人的,不是关于社会的
  W-CODER     同源 |ρ| 的分布显著高于跨源 ⇒ **我们测到的是那一个编码团队**(`#521`)
⚠ BASIN:`W-REPLICATE` 是我**想要**的结局,所以它不是本轮下注的方向。
   本轮真正下注 `W-CODER` —— 它一旦为真,**E02 至今为止所有 SCCS 结论都只能写成关于编码者的**。

CONTROLS(G2):
  正对照 SCCS166×SCCS167(婚前性频率:男 × 女,同一构念不同性别)必须强正,
     且在**置换后不通过**(g=0 必须失败);
  安慰剂 谴责度 × 社会**经度**(任意地理坐标)必须 ≈ 0 ——「这个零该不该是零?」**该**⇒ negative_control;
  参照分布 跨源对 |ρ| 的中位数(RULE-v3 用中位,不用 q95);
  n_eff Galton 问题:社会不独立。用 **10°×10° 经纬网格块**作聚类做 bootstrap,报块数。
KILL(条件式,预注册):
  if 正对照在真数据上通过 and 置换后不通过 and 安慰剂 ≈ 0:
      同源中位 |ρ| > 跨源中位 |ρ| 的 2 倍 -> W-CODER
      elif 三个 ρ 都 > 0 且都超各自 MDE     -> W-REPLICATE
      elif 任一 ρ < 0 且超 MDE              -> W-DIES
      else                                  -> UNVERIFIED-by-power
  else: UNVERIFIED
IMPOSSIBLE(结构上做不到):无干预 ⇒ 非因果 · 无系统发生树 ⇒ 只有空间聚类,不是真正的
  Galton 校正 · **一个编码团队 ⇒ 无跨仪器复制**(`#521`)· 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SC = ROOT / "data/external/dplace/repo/datasets/SCCS"
D = pd.read_csv(SC / "data.csv", low_memory=False)
V = pd.read_csv(SC / "variables.csv")
S = pd.read_csv(SC / "societies.csv")
W = D.pivot_table(index="soc_id", columns="var_id", values="code", aggfunc="first")
S = S.set_index("id").reindex(W.index)
W["_lat"], W["_lon"] = S.Lat.values, S.Long.values
W["_blk"] = (np.floor(W._lat / 10).astype("Int64").astype(str) + "_" +
             np.floor(W._lon / 10).astype("Int64").astype(str))
print(f"社会数 {len(W)} · 变量数 {W.shape[1]} · 10°网格块 {W._blk.nunique()}")

def sp(a, b):
    """中位秩 Spearman(`#482a`:argsort(argsort()) 给并列任意秩)。"""
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 20: return np.nan, int(m.sum())
    return float(np.corrcoef(rankdata(a[m]), rankdata(b[m]))[0, 1]), int(m.sum())

PAIRS = [("婚前性(女)", "SCCS165", "SCCS167"), ("婚前性(男)", "SCCS165", "SCCS166"),
         ("婚外性(女)", "SCCS169", "SCCS171"), ("婚外性(男)", "SCCS169", "SCCS170"),
         ("同性恋",     "SCCS176", "SCCS177")]
W.loc[W.SCCS176 == 2, "SCCS176"] = np.nan       # 码 2 无描述 -> 剔除并声明

print("\n=== 规则①:逐格打印 n,再看 ρ。谴责码大 = 更谴责;频率码大 = 更稀有 ⇒ ρ>0 即「越稀有越谴责」===")
rows = []
for nm, av, fv in PAIRS:
    a, f = W[av].values.astype(float), W[fv].values.astype(float)
    r, n = sp(a, f)
    bs = []
    blocks = W._blk.values
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        ub = pd.unique(blocks)
        for _ in range(400):
            pick = rng.choice(ub, len(ub))
            idx = np.concatenate([np.where(blocks == b)[0] for b in pick])
            v, _n = sp(a[idx], f[idx])
            if np.isfinite(v): bs.append(v)
    bs = np.array(bs); MDE = 2.8 * bs.std() if len(bs) else np.nan
    ci = np.quantile(bs, [.025, .975]) if len(bs) else [np.nan, np.nan]
    rows.append(dict(pair=nm, att=av, freq=fv, n=n, rho=r, MDE=float(MDE),
                     ci=[float(ci[0]), float(ci[1])], sd=float(bs.std()) if len(bs) else None,
                     n_blocks=int(W._blk.nunique()),
                     inclusion=[f"两列都非缺失的社会 (n={n})", "SCCS176 码 2(码本无描述)已剔除",
                                "块 bootstrap:10°×10° 经纬格"]))
    print(f"  {nm:12s} {av}×{fv}  n={n:3d}  **ρ={r:+.4f}**  MDE={MDE:.4f}  "
          f"CI [{ci[0]:+.4f},{ci[1]:+.4f}]  {'超MDE' if abs(r)>MDE else '**看不见**'}")

# ---- 共享仪器:同源 vs 跨源参照分布
BR = list(V[V.source == "broude1976cross"].id)
MU = list(V[V.source == "murdock1973factors"].id)[:60]
OTH = [v for v in V[~V.source.isin(["broude1976cross"])].id if v in W.columns][:120]
def dist(av, others):
    out = []
    for o in others:
        if o == av or o not in W.columns: continue
        r, n = sp(W[av].values.astype(float), pd.to_numeric(W[o], errors="coerce").values.astype(float))
        if np.isfinite(r) and n >= 30: out.append(abs(r))
    return np.array(out)
same = np.concatenate([dist(a, BR) for a in ["SCCS165", "SCCS169", "SCCS176"]])
cross = np.concatenate([dist(a, OTH) for a in ["SCCS165", "SCCS169", "SCCS176"]])
print(f"\n=== 共享仪器(最强混淆,写在跑之前)===")
print(f"  同源对 Broude×Broude  k={len(same):4d}  中位 |ρ| = {np.median(same):.4f}")
print(f"  跨源对 Broude×其他来源 k={len(cross):4d}  中位 |ρ| = {np.median(cross):.4f}")
ratio = float(np.median(same) / np.median(cross))
print(f"  **比值 = {ratio:.3f}**  (预注册:> 2.0 -> W-CODER)")

# ---- 对照
G = Gate("186 个社会之间,越稀有的做法是否越被谴责?(SCCS/Broude)")
pc, pcn = sp(W.SCCS166.values.astype(float), W.SCCS167.values.astype(float))
rng = np.random.default_rng(SEEDS[0])
aa, ff = W.SCCS166.values.astype(float), W.SCCS167.values.astype(float)
m = np.isfinite(aa) & np.isfinite(ff)
perm = [abs(sp(aa[m][rng.permutation(m.sum())], ff[m])[0]) for _ in range(300)]
print(f"\n=== 对照 ===\n  正对照 SCCS166×SCCS167(婚前频率 男×女)n={pcn} ρ={pc:+.4f} · "
      f"置换 q95={np.quantile(perm,.95):.4f} · 跨源参照中位={np.median(cross):.4f}")
G.positive_control("正对照-v3:同一构念不同性别", planted=abs(pc),
                   floor=max(float(np.quantile(perm, .95)), float(np.median(cross))), spread=1e-9)
G.negative_control("g=0:置换后必须不通过",
                   null=float(np.median(perm)), effect=abs(pc), null_spread=float(np.std(perm)),
                   null_kind="社会层标签置换")
plc = []
for av in ["SCCS165", "SCCS169", "SCCS176"]:
    r, n = sp(W[av].values.astype(float), W._lon.values.astype(float)); plc.append(abs(r))
G.negative_control("安慰剂:谴责度 × 社会经度", null=float(np.median(plc)), effect=abs(pc),
                   null_spread=float(np.std(plc)), null_kind="任意地理坐标(应无关)")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {r["pair"]: r for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {r["pair"]: r for r in rows})

print("\n" + "=" * 74)
if all([abs(pc) > max(np.quantile(perm, .95), np.median(cross)), np.median(plc) < 0.5 * abs(pc)]):
    seen = [r for r in rows if np.isfinite(r["rho"]) and abs(r["rho"]) > r["MDE"]]
    if ratio > 2.0:
        world = "W-CODER"; verdict = f"同源/跨源中位 |ρ| 比值 {ratio:.3f} > 2 -> **测到的是那一个编码团队**"
    elif seen and all(r["rho"] > 0 for r in seen) and len(seen) >= 3:
        world = "W-REPLICATE"; verdict = f"{len(seen)} 格超 MDE 且全为正 -> **跨单位跨仪器复制成功**"
    elif any(r["rho"] < 0 for r in seen):
        world = "W-DIES"; verdict = (f"{[r['pair'] for r in seen if r['rho']<0]} 为**负**且超 MDE -> "
                                     f"**人层面的关系在社会之间不出现,甚至反号**")
    else:
        world = "UNVERIFIED"; verdict = f"只有 {len(seen)}/5 格超 MDE -> **UNVERIFIED-by-power**"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:态度与频率同源同笔,一句民族志可以同时决定两个码;"
          "而块 bootstrap 只处理空间邻近,**不处理语言谱系** —— 真正的 Galton 校正这里做不到。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, world=world, verdict=verdict, same_median=float(np.median(same)),
               cross_median=float(np.median(cross)), ratio=ratio, k_same=len(same), k_cross=len(cross),
               positive_control=dict(rho=pc, n=pcn, perm_q95=float(np.quantile(perm, .95))),
               placebo=[float(x) for x in plc], seeds=SEEDS, n_blocks=int(W._blk.nunique()),
               instrument="Broude & Greene 1976,单一编码团队(#521)",
               impossible=["无干预非因果", "无系统发生树,块 bootstrap 不是 Galton 校正",
                           "一个编码团队,无跨仪器复制"], unchallenged=True),
          open(OUT / "sccs_rarity_condemnation.json", "w"), indent=1)
print(f"\nwrote {OUT/'sccs_rarity_condemnation.json'}")
