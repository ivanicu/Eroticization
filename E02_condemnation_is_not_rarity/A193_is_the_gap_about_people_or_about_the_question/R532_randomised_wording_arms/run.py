"""E02·A193·R532 — 同一群人随机劈半,一半说 32% 一半说 57.5%,而我上一轮把其中一半发上了页面

`#486` 的 NEXT 是用 `porn30` 问「差距是关于人还是关于题面」。
⚠ **规则 ① 先跑,然后它撞出一个比 NEXT 更好的东西**:
  `porn30` 与 `pornlaw` 只在 **2 年、n=620** 上重叠;`xmovie1` 只有 **1994 一年**;
  但 **`xmoviey` 在 2022/2024 与 `xmovie` 同时存在,却从不落在同一个人身上**
  -> GSS **把样本劈成了两种问法**,ballot 2 与 3 里两臂都有,人口学平衡
  (|Δ/se| ≤ 1.96,5 项比较里最大 1.96)-> **一次问卷内部的随机题面实验。**

⛔ **而这直接指向我上一轮发上页面的东西。** `#486b` 的五十一年零用的是 `xmovie` 这个**二值题**。
   若二值格式系统性压缩了报告,那个零可能是**格式的性质**,不是人的性质。

G1 ESTIMAND(先于方法):
  `Δ_lnOR = lnOR(谴责, 行为阳性 | 臂 = xmoviey) − lnOR(谴责, 行为阳性 | 臂 = xmovie)`,
  合并 2022+2024。`xmoviey` 是 10 级频率,**预注册主切点 = 「从不」vs「非从不」**
  (最贴近「过去一年看过」的语义);切点作为规格轴全格公布。

⚠ 题面(读过全文,不是猜的):
  `xmovie`  = "Seen x-rated movie in last year"      二值 yes/no
  `xmoviey` = "SEEN PORNOGRAPHY IN PAST YEAR"        10 级频率
⛔ **STRONGEST CONFOUND,写在跑之前:这个操纵是混淆的。**
  **对象**(「X 级电影」vs「色情内容」)与**作答格式**(二值 vs 频率)**同时变了**。
  任何 Δ **不可归因于其中任一个**。这是 GSS 设计的性质,在此不可拆开。

⚠ P6 代理账本 —— **本设计只在一个方向上可靠**:
  两臂都是**同一份问卷里的自报**,共享方法方差。
  **「两臂一致」几乎不授权任何东西;「两臂不一致」才是有信息的那个方向。**
  ⇒ 一个零在这里**不能**读成「差距是关于人的」,只能读成一个**大小的上界**。

WORLDS:
  W-FORMAT  差距有一部分是仪器造的 -> `|Δ_lnOR|` 大,或随切点系统性变化
  W-PERSON  差距是人的性质         -> `Δ_lnOR ≈ 0` 且跨切点稳定
  | World    | now | \|Δ\|>MDE | Δ≈0 |
  | W-FORMAT | 0.5 | 0.85      | 0.15 |
  | W-PERSON | 0.5 | 0.15      | 0.85 |
  ⚠ 但按上一段,只有 `|Δ|>MDE` 这一格真的移动信念。

CONTROLS:
  正对照 谴责×礼拜出席,**两臂各一次**,门槛 = 同问卷参照分布 q95(`#485a`)
  安慰剂 谴责×星座,两臂各一次 —— 必须为零
  平衡   已算(age/sex/degree/attend/polviews),随结果一起报
  精度   **人层 bootstrap**,臂内重抽

KILL(条件式,预注册):
  if 正对照在**两臂**都触发 and 安慰剂在**两臂**都为零:
      |Δ| > MDE 且 CI 不含 0 -> W-FORMAT:`#486b` 的零是**仪器绑定**的,页面必须重新界定范围
      |Δ| < MDE 且 CI 含 0 且 MDE < 0.5 -> **不授权 W-PERSON**(共享方法),只报一个上界
  else: UNVERIFIED

IMPOSSIBLE:对象与格式不可拆(GSS 设计)· 无干预 ⇒ 非因果 ·
  未派对抗 agent(会话约束)⇒ `[unchallenged]` · 只有 2 个年份 ⇒ 无时代分辨
"""
import os, sys, pathlib, json, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate, check_columns

SEEDS = [20260805, 7, 991]
MEANINGFUL = 0.5
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
REF529 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/"
                        "A192_does_condemnation_suppress_the_report_or_the_act/"
                        "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))["reference"]
T_POS = float(np.quantile([abs(x["r"]) for x in REF529], .95))

cols = ["year", "age", "sex", "degree", "attend", "polviews", "zodiac",
        "pornlaw", "xmovie", "xmoviey", "ballot"]
it = pd.read_stata(DTA, iterator=True); have = set(it.variable_labels())
d = pd.read_stata(DTA, columns=[c for c in cols if c in have], convert_categoricals=False)
check_columns(d, where="R532")
d = d[d.year.isin([2022, 2024])].copy()
d["arm"] = np.where(d.xmovie.notna(), "xmovie",
                    np.where(d.xmoviey.notna(), "xmoviey", None))
d = d[d.arm.notna() & d.pornlaw.notna()].copy()
d["condemn"] = (d.pornlaw == 1).astype(float)
print(f"n={len(d)}  两臂 {d.arm.value_counts().to_dict()}  年份 {sorted(d.year.unique().astype(int))}")

# ---------------------------------------------------------------- 平衡
print("\n=== 平衡(随机分配则应无差;5 项比较,期望最大 |z| ≈ 1.9)===")
bal = []
for v in ["age", "sex", "degree", "attend", "polviews"]:
    a = d[d.arm == "xmovie"][v].dropna(); b = d[d.arm == "xmoviey"][v].dropna()
    se = math.sqrt(a.var() / len(a) + b.var() / len(b))
    z = (a.mean() - b.mean()) / se
    bal.append(dict(var=v, mean_xmovie=float(a.mean()), mean_xmoviey=float(b.mean()), z=float(z)))
    print(f"  {v:9s} {a.mean():7.3f} vs {b.mean():7.3f}   z={z:+.2f}")
print(f"  max|z| = {max(abs(x['z']) for x in bal):.2f}")

# ---------------------------------------------------------------- 基率:纯随机对比
CUTS = [("never", 1), ("ge_lt_once_yr", 2), ("ge_once_twice_yr", 3),
        ("ge_several_yr", 4), ("ge_monthly", 5), ("ge_weekly", 7)]
# xmoviey 的码序需要先看:打印一次(规则①)
lab = pd.read_stata(DTA, columns=["xmoviey"], convert_categoricals=True)
raw = pd.read_stata(DTA, columns=["xmoviey"], convert_categoricals=False)
mapping = pd.DataFrame({"code": raw.xmoviey, "label": lab.xmoviey}).dropna().drop_duplicates()
mapping = mapping.sort_values("code")
print("\n=== xmoviey 码 → 标签(规则①:先打印)===")
for _, r in mapping.iterrows(): print(f"  {int(r['code'])} = {r['label']}")

A = d[d.arm == "xmovie"]
B = d[d.arm == "xmoviey"]
pa = (A.xmovie == 1).mean()
print(f"\n=== 随机劈半的基率对比(控制无关,这是分配本身给的)===")
print(f"  xmovie  「过去一年看过 X 级电影」= {pa:.4f}  (n={len(A)})")
for nm, c in CUTS[:2]:
    pb = (B.xmoviey > c - 1).mean() if nm != "never" else (B.xmoviey > mapping.code.min()).mean()
print(f"  xmoviey 「过去一年看过色情内容」非「从不」= {(B.xmoviey > mapping.code.min()).mean():.4f}  (n={len(B)})")
rng0 = np.random.default_rng(SEEDS[0])
bd = np.array([(B.xmoviey.sample(len(B), replace=True, random_state=int(rng0.integers(1e9)))
                > mapping.code.min()).mean()
               - (A.xmovie.sample(len(A), replace=True, random_state=int(rng0.integers(1e9))) == 1).mean()
               for _ in range(1000)])
print(f"  Δ基率 = {bd.mean():+.4f}  95% CI [{np.quantile(bd,.025):+.4f}, {np.quantile(bd,.975):+.4f}]")


def lnor(dd, pos):
    a, b = dd[dd.condemn == 1], dd[dd.condemn == 0]
    if len(a) < 30 or len(b) < 30: return np.nan
    p1, p0 = pos(a).mean(), pos(b).mean()
    if min(p1, p0) <= 0 or max(p1, p0) >= 1: return np.nan
    return math.log((p1 / (1 - p1)) / (p0 / (1 - p0)))


MINC = mapping.code.min()
posA = lambda x: (x.xmovie == 1)
def posB_at(cut): return lambda x: (x.xmoviey > cut)

lnA = lnor(A, posA)
lnB = lnor(B, posB_at(MINC))
main_delta = lnB - lnA
print(f"\n=== 主:两臂的 lnOR ===")
print(f"  xmovie  lnOR = {lnA:+.4f}")
print(f"  xmoviey lnOR = {lnB:+.4f}   (切点=非「从不」)")
print(f"  **Δ = {main_delta:+.4f}**")


def boot_delta(B_=800, seed=0, cut=MINC):
    rng = np.random.default_rng(seed); out = []
    for _ in range(B_):
        aa = A.iloc[rng.integers(0, len(A), len(A))]
        bb = B.iloc[rng.integers(0, len(B), len(B))]
        x, y = lnor(aa, posA), lnor(bb, posB_at(cut))
        if np.isfinite(x) and np.isfinite(y): out.append(y - x)
    return np.array(out)


bs = [boot_delta(800, s) for s in SEEDS]
ball = np.concatenate(bs)
se = float(ball.std()); MDE = 2.8 * se
lo_ci, hi_ci = np.quantile(ball, [.025, .975])
print(f"  95% CI [{lo_ci:+.4f}, {hi_ci:+.4f}]  se={se:.4f}  **MDE={MDE:.4f}**  "
      f"seed_spread={np.std([b.mean() for b in bs]):.5f}")

# ---------------------------------------------------------------- G4 切点规格曲线
print("\n=== G4 规格曲线:切点 + 子群(全格公布)===")
spec = []
codes = sorted(mapping.code.astype(int).tolist())
for i, cut in enumerate(codes[:-2]):
    v = lnor(B, posB_at(cut))
    if np.isfinite(v):
        spec.append(dict(kind="cut", cut=int(cut),
                         label=str(mapping[mapping.code == cut].label.iloc[0])[:28],
                         lnor_B=v, delta=v - lnA))
        print(f"  切点 > {int(cut)} ({spec[-1]['label']:28s}) lnOR_B={v:+.4f}  Δ={v-lnA:+.4f}")
for nm, m in [("age<50", d.age < 50), ("age>=50", d.age >= 50),
              ("HS+", d.degree >= 1), ("attend_lo", d.attend <= 2)]:
    aa, bb = d[(d.arm == "xmovie") & m], d[(d.arm == "xmoviey") & m]
    x, y = lnor(aa, posA), lnor(bb, posB_at(MINC))
    if np.isfinite(x) and np.isfinite(y):
        spec.append(dict(kind="subgroup", name=nm, lnor_A=x, lnor_B=y, delta=y - x))
        print(f"  {nm:10s} lnOR_A={x:+.4f} lnOR_B={y:+.4f}  Δ={y-x:+.4f}")
ds = [s["delta"] for s in spec]
sg = [np.sign(v) for v in ds]; dom = max(set(sg), key=sg.count)
print(f"\nspec_survival: {sg.count(dom)}/{len(sg)} = {sg.count(dom)/len(sg):.0%} 同号 ({dom:+.0f});"
      f" 超过 MDE 的格数 = {sum(1 for v in ds if abs(v) > MDE)}/{len(ds)}")

# ---------------------------------------------------------------- 控制,两臂各一次
G = Gate("差距是关于人的,还是关于题面的?(GSS 2022/2024 随机题面臂)")
oks = {}
for nm, arm in [("xmovie", A), ("xmoviey", B)]:
    s = arm.dropna(subset=["attend"])
    r = float(np.corrcoef(s.condemn, s.attend)[0, 1])
    oks[f"pos_{nm}"] = G.positive_control(f"正对照[{nm}]:谴责×礼拜出席(门槛=参照 q95)",
                                          planted=abs(r), floor=T_POS, spread=1e-9)
    z = arm.dropna(subset=["zodiac"])
    zr = float(np.corrcoef(z.condemn, z.zodiac)[0, 1])
    rngz = np.random.default_rng(SEEDS[0])
    zn = np.array([np.corrcoef(z.condemn.values[rngz.permutation(len(z))], z.zodiac)[0, 1]
                   for _ in range(300)])
    oks[f"plac_{nm}"] = G.negative_control(f"安慰剂[{nm}]:谴责×星座", null=zr, effect=r,
                                           null_spread=float(zn.std()), null_kind="个体层标签置换")
G.has_error_bar("Δ_lnOR", value=main_delta, spread=se, spread_source="bootstrap_人层")

all_pos = oks["pos_xmovie"] and oks["pos_xmoviey"]
all_plac = oks["plac_xmovie"] and oks["plac_xmoviey"]
ci_excl = not (lo_ci <= 0 <= hi_ci)
print("\n" + "=" * 70)
if all_pos and all_plac:
    if abs(main_delta) > MDE and ci_excl:
        verdict = (f"|Δ|={abs(main_delta):.4f} > MDE={MDE:.4f} 且 CI 不含 0 -> "
                   f"**W-FORMAT:`#486b` 的零是仪器绑定的,页面必须重新界定范围**")
    elif abs(main_delta) < MDE and not ci_excl and MDE < MEANINGFUL:
        verdict = (f"|Δ|={abs(main_delta):.4f} < MDE={MDE:.4f},CI 含 0 -> "
                   f"⚠ **不授权 W-PERSON**(两臂共享自报方法);只报一个上界 |Δ| < {MDE:.3f}")
    else:
        verdict = f"Δ={main_delta:+.4f},MDE={MDE:.4f} -> 未决"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:操纵是混淆的 —— 对象与格式同时变了,"
          "Δ 不可归因于任一个;而一个零仍然只是上界,因为两臂共享同一份问卷的自报方法。")
else:
    verdict = f"UNVERIFIED —— 控制未齐 {oks}"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(balance=bal, n_arms=d.arm.value_counts().to_dict(),
               base_xmovie=float(pa), base_xmoviey=float((B.xmoviey > MINC).mean()),
               base_delta=float(bd.mean()),
               base_delta_ci=[float(np.quantile(bd, .025)), float(np.quantile(bd, .975))],
               lnor_xmovie=lnA, lnor_xmoviey=lnB, delta=main_delta,
               ci=[float(lo_ci), float(hi_ci)], se=se, MDE=MDE, meaningful=MEANINGFUL,
               spec=spec, controls=oks, verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "randomised_wording_arms.json", "w"), indent=1)
print(f"\nwrote {OUT/'randomised_wording_arms.json'}")
