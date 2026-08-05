"""E02·A192·R529 — 谴责崩塌的五十一年里,谴责者和不谴责者报告的差距有没有跟着变

⚠ 本轮**不执行 `#483` 注册的 NEXT**,并说明为什么(frontier §3 盆地):
  `A190·R526` UNVERIFIED -> `A190·R527` UNVERIFIED -> `A191·R528` 评判 —— **三轮同一个问题**
  (共享编码者混淆是否足以解释 +0.7966),而三轮的产出**都是「我的仪器坏在哪」**。
  `#483` 的 NEXT 是**第四轮**同一个问题。
  更要命的是 `#483f` 我自己写下:**W-CODER 与 W-DOMAIN 在此同义** ——
  **我一直在分离一个退化的世界分解**,而 meta-separator 要求的正是发现这件事。
  ⇒ 换方向,并且换**仪器**(规则 ④)。`#483` 的 NEXT 记为**未执行**,不是失败。

SAMPLE:`#466` 里那个**从来没被测过**的世界。
  `#466` 写:本设计把 {W-A 稀有→谴责, W-C 谴责→稀有} 与 {W-B 无关} 分开,**分不开 A 和 C**。
  SCCS 结构上分不开:它没有「同一个人的谴责」与「同一个人的行为」。
  **GSS 有,而且有五十一年,期间社会谴责崩塌。**

G1 ESTIMAND(先于方法):**每一年内**,谴责者与非谴责者在**自报行为**上的风险差
  `RD_t = P(看过 | 谴责) − P(看过 | 不谴责)`;
  然后:`RD_t` 随该年**社会谴责率**下降如何变化(斜率)。

WORLDS:
  W-CONCEAL   谴责压低的是**报告**:谴责规范强时,谴责者更不敢承认 -> |RD| 被隐瞒抬高
              -> **社会谴责率下降时 |RD| 应当收缩**
  W-PREFERENCE 谴责反映**真实不想要** -> **|RD| 跨时代稳定**
  W-BOTH       部分收缩
  | World       | now | |RD| 收缩 | |RD| 稳定 |
  | W-CONCEAL   | 0.35| 0.85     | 0.05     |
  | W-PREFERENCE| 0.35| 0.05     | 0.85     |
  | W-BOTH      | 0.30| 0.40     | 0.20     |
  没有平行行。

⚠ STRONGEST CONFOUND,写在跑之前(三个):
  ① **成分**:谴责率下降时,仍在谴责的人变成更被选择的少数,|RD| 可因此变化而与隐瞒无关。
     控制:在 宗教×出席×学历×世代 分层内重算;并检验 |RD| 是否**只由谴责者占比**预测。
  ② **技术**:X 片的可得性在 51 年里剧变(录像带、互联网),基率会因技术移动。
     -> 报 |RD| 时**同时报基率**;RD 是组间差,对整体基率的敏感度低于水平本身。
  ③ **题面**:只用 `xmovie`(1973–2024 连续),不混 `xmovie1`/`xmoviey`。

CONTROLS:
  正对照   `pornlaw` × `attend`(礼拜出席)—— 教科书级相关,必须强;算 floor **和** ceiling(`#482d`)
  阴性     **不挑一个零,测一个参照分布**(`#482c`/`#483d` 的教训):
           `pornlaw` 与一批与性无关的 GSS 变量的关联分布
  安慰剂   `zodiac`(星座,由出生月导出)—— 这个零**必须**是零 ⇒ negative_control

KILL(条件式,预注册,写在跑之前):
  if 正对照触发 and 安慰剂为零:
      slope(|RD| vs 该年谴责率) 显著为**正**(谴责率高 -> |RD| 大)-> W-CONCEAL
      slope 落在自身零展布内                                     -> W-PREFERENCE
      slope 显著为负                                             -> 两个世界都没预测,MIXED
  else: UNVERIFIED

IMPOSSIBLE(附所需):
  causally identified     -> 需要对社会谴责的干预;这里只有观测的时间变异
  independently replicated-> **未派对抗 agent(会话约束)** ⇒ `[unchallenged]`
  个体内纵向              -> GSS 是重复横截面,不是面板;**「同一个人变了」在此不可测**,
                             只能测「同一世代的人在不同年份」
"""
import os, sys, pathlib, json, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate, check_columns

SEEDS = [20260805, 7, 991]
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)

CORE = ["year", "age", "cohort", "sex", "degree", "relig", "attend", "polviews", "region",
        "pornlaw", "xmovie", "wtssps"]
REF = ["tvhours", "happy", "health", "satfin", "sibs", "childs", "hompop",
       "zodiac", "res16", "mobile16", "wrkstat", "dwelown"]

it = pd.read_stata(DTA, iterator=True)
have = set(it.variable_labels())
REF = [c for c in REF if c in have]
cols = [c for c in CORE + REF if c in have]
df = pd.read_stata(DTA, columns=cols, convert_categoricals=False)
check_columns(df, where="R529")

# ---------------------------------------------------------------- 规则①:先验码,不假设
print("=== 规则①:码与 n,先打印再相信 ===")
lab = pd.read_stata(DTA, columns=["pornlaw", "xmovie", "zodiac"], convert_categoricals=True)
for c in ["pornlaw", "xmovie"]:
    vc = lab[c].value_counts()
    print(f"  {c}: " + " | ".join(f"{k}={v}" for k, v in vc.items()))
assert df["pornlaw"].dropna().isin([1, 2, 3]).all(), "pornlaw 码不是 {1,2,3}"
assert df["xmovie"].dropna().isin([1, 2]).all(), "xmovie 码不是 {1,2}"
d = df.dropna(subset=["pornlaw", "xmovie"]).copy()
d["condemn"] = (d.pornlaw == 1).astype(float)      # 「应对所有人非法」= 谴责
d["saw"] = (d.xmovie == 1).astype(float)           # 已由上面的标签核对
print(f"  分析样本 n={len(d)}  年份 {int(d.year.min())}-{int(d.year.max())}  "
      f"{d.year.nunique()} 个年份")
print(f"  谴责率 {d.groupby('year').condemn.mean().iloc[0]:.3f} -> "
      f"{d.groupby('year').condemn.mean().iloc[-1]:.3f};"
      f" 看过率 {d.groupby('year').saw.mean().iloc[0]:.3f} -> "
      f"{d.groupby('year').saw.mean().iloc[-1]:.3f}")


def rd_by_year(dd, wcol=None):
    """每年的风险差 RD = P(saw|condemn) - P(saw|~condemn),以及该年谴责率。"""
    out = []
    for y, g in dd.groupby("year"):
        a, b = g[g.condemn == 1], g[g.condemn == 0]
        if len(a) < 30 or len(b) < 30: continue
        if wcol:
            pa = np.average(a.saw, weights=a[wcol]); pb = np.average(b.saw, weights=b[wcol])
            share = np.average(g.condemn, weights=g[wcol])
        else:
            pa, pb, share = a.saw.mean(), b.saw.mean(), g.condemn.mean()
        out.append(dict(year=int(y), rd=float(pa - pb), share=float(share),
                        base=float(g.saw.mean()), n=len(g), na=len(a), nb=len(b)))
    return pd.DataFrame(out)


base = rd_by_year(d)
print("\n=== 每年 RD、谴责率、基率(全部 30 年公布)===")
print(base.assign(rd=base.rd.round(4), share=base.share.round(3), base=base.base.round(3))
      .to_string(index=False))

# ---------------------------------------------------------------- 主:斜率
def slope_of(bb):
    x = bb.share.values; y = np.abs(bb.rd.values)
    if len(x) < 8: return np.nan
    return float(np.polyfit(x, y, 1)[0])


main_slope = slope_of(base)
# 零:年内打乱 condemn 标签(破坏个体层配对,保留年份结构与谴责率)
def slope_null(dd, n=400, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(n):
        dd2 = dd.copy()
        dd2["condemn"] = dd2.groupby("year").condemn.transform(
            lambda s: s.values[rng.permutation(len(s))])
        bb = rd_by_year(dd2)
        s = slope_of(bb)
        if np.isfinite(s): out.append(s)
    return np.array(out)


nulls = [slope_null(d, 400, s) for s in SEEDS]
null_all = np.concatenate(nulls)
sd = float(null_all.std()); seed_spread = float(np.std([x.mean() for x in nulls]))
print(f"\n主:slope(|RD| ~ 该年谴责率) = {main_slope:+.4f}")
print(f"    零(年内打乱 condemn)mean={null_all.mean():+.4f} sd={sd:.4f}  "
      f"seed_spread={seed_spread:.5f}  |slope|/sd = {abs(main_slope)/max(sd,1e-12):.2f}x")
print(f"    单边 p = {(null_all >= main_slope).mean():.4f}")

# ---------------------------------------------------------------- G4 规格曲线
print("\n=== G4 规格曲线(全格公布,含反号格)===")
spec = []
for wname, wcol in [("unweighted", None), ("wtssps", "wtssps")]:
    for cname, mask in [("all", d.index),
                        ("HS+", d.index[d.degree >= 1]),
                        ("<HS", d.index[d.degree == 0]),
                        ("attend_hi", d.index[d.attend >= 5]),
                        ("attend_lo", d.index[d.attend <= 2]),
                        ("cohort<1950", d.index[d.cohort < 1950]),
                        ("cohort>=1950", d.index[d.cohort >= 1950])]:
        dd = d.loc[mask]
        if wcol: dd = dd.dropna(subset=[wcol])
        bb = rd_by_year(dd, wcol)
        s = slope_of(bb)
        spec.append(dict(weight=wname, stratum=cname, slope=s, n_years=len(bb), n=len(dd)))
        print(f"  {wname:11s} {cname:13s} slope={s:+.4f}  years={len(bb):2d} n={len(dd):6d}"
              if np.isfinite(s) else
              f"  {wname:11s} {cname:13s} slope=  n/a   years={len(bb):2d} n={len(dd):6d}")
fin = [s["slope"] for s in spec if np.isfinite(s["slope"])]
signs = [np.sign(v) for v in fin]; dom = max(set(signs), key=signs.count)
print(f"\nspec_survival: {signs.count(dom)}/{len(fin)} = {signs.count(dom)/len(fin):.0%} "
      f"同号 (sign={dom:+.0f})")

# ---------------------------------------------------------------- 混淆①:占比本身能否解释
print("\n=== 混淆①:|RD| 是不是只由谴责者占比机械地决定 ===")
# 若 RD 完全由占比驱动,则在占比相同的年份对上,|RD| 应无残差变异
r_share = np.corrcoef(base.share, np.abs(base.rd))[0, 1]
r_base = np.corrcoef(base.base, np.abs(base.rd))[0, 1]
r_year = np.corrcoef(base.year, np.abs(base.rd))[0, 1]
print(f"  corr(|RD|, 谴责率)={r_share:+.4f}   corr(|RD|, 看过基率)={r_base:+.4f}   "
      f"corr(|RD|, 年份)={r_year:+.4f}")
# 把基率回归掉之后,谴责率还剩多少
X = np.c_[np.ones(len(base)), base.base.values]
res = np.abs(base.rd.values) - X @ np.linalg.lstsq(X, np.abs(base.rd.values), rcond=None)[0]
r_share_resid = float(np.corrcoef(base.share, res)[0, 1])
print(f"  扣掉看过基率后 corr(|RD|残差, 谴责率) = {r_share_resid:+.4f}")

# ---------------------------------------------------------------- 控制
G = Gate("谴责压低的是报告还是行为?(GSS 1973-2024)")

# 正对照:pornlaw × attend
pc = d.dropna(subset=["attend"])
pc_r = float(np.corrcoef(pc.condemn, pc.attend)[0, 1])
rngp = np.random.default_rng(SEEDS[0])
pc_null = np.array([np.corrcoef(pc.condemn.values[rngp.permutation(len(pc))], pc.attend)[0, 1]
                    for _ in range(400)])
pc_floor = float(np.quantile(np.abs(pc_null), .95))
# ceiling:两列各自排序后配对(给定边际的上界)
pc_ceil = float(np.corrcoef(np.sort(pc.condemn.values), np.sort(pc.attend.values))[0, 1])
t = (pc_floor + pc_ceil) / 2
print(f"\n正对照 谴责 × 礼拜出席 r={pc_r:+.4f} n={len(pc)}  floor={pc_floor:.4f} "
      f"ceiling={pc_ceil:.4f} 门槛={t:.4f}")
pc_ok = G.positive_control("正对照:谴责 × 礼拜出席(门槛在 floor<t<ceiling 带内)",
                           planted=abs(pc_r), floor=t, spread=1e-9)

# 安慰剂:zodiac —— 这个零必须是零
if "zodiac" in d.columns:
    z = d.dropna(subset=["zodiac"])
    z_r = float(np.corrcoef(z.condemn, z.zodiac)[0, 1])
    rngz = np.random.default_rng(SEEDS[0])
    z_null = np.array([np.corrcoef(z.condemn.values[rngz.permutation(len(z))], z.zodiac)[0, 1]
                       for _ in range(400)])
    print(f"安慰剂 谴责 × 星座 r={z_r:+.4f} n={len(z)}  零 sd={z_null.std():.4f}")
    nc_ok = G.negative_control("安慰剂:谴责 × 星座(由出生月导出,必须为零)",
                               null=z_r, effect=pc_r, null_spread=float(z_null.std()),
                               null_kind="个体层标签置换(保留两侧边际)")
else:
    z_r, nc_ok = float("nan"), False
    print("⚠ zodiac 不在这份释放里 -> 安慰剂缺席")

# 参照分布(测量出来的同问卷基线,而不是挑一个零)
print("\n=== 参照分布:谴责 × 一批与性无关的 GSS 变量(测量,非挑选)===")
ref = []
for c in REF:
    if c not in d.columns: continue
    dd = d.dropna(subset=[c])
    if len(dd) < 2000 or dd[c].nunique() < 3: continue
    r = float(np.corrcoef(dd.condemn, dd[c])[0, 1])
    ref.append(dict(var=c, r=r, n=len(dd)))
ref.sort(key=lambda x: -abs(x["r"]))
if ref:
    ar = np.array([abs(x["r"]) for x in ref])
    print(f"  {len(ref)} 个参照变量  |r| 中位={np.median(ar):.4f}  max={ar.max():.4f}")
    for x in ref[:6]: print(f"    |{x['r']:+.4f}| n={x['n']:6d}  {x['var']}")

G.has_error_bar("斜率", value=main_slope, spread=sd, spread_source="null_零臂")
G.resolvable("斜率可分辨", effect=abs(main_slope), spread=sd)

print("\n" + "=" * 70)
if pc_ok and nc_ok:
    if main_slope > 2 * sd:
        verdict = "斜率显著为正(谴责率高的年份 |RD| 大)-> W-CONCEAL 得到支持"
    elif abs(main_slope) <= 2 * sd:
        verdict = "斜率落在自身零展布内 -> W-PREFERENCE:差距跨时代稳定"
    else:
        verdict = "斜率显著为负 -> 两个世界都没预测这个方向;MIXED/未知"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:GSS 是**重复横截面**,不是面板 ——"
          " 「同一个人变了」在此不可测,任何解释都只在**同世代不同年份**这一层成立。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} placebo={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(by_year=base.to_dict("records"), slope=main_slope, null_sd=sd,
               null_mean=float(null_all.mean()), seed_spread=seed_spread,
               p_one_sided=float((null_all >= main_slope).mean()),
               spec=[{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                      for k, v in s.items()} for s in spec],
               confound_share=dict(corr_rd_share=float(r_share), corr_rd_base=float(r_base),
                                   corr_rd_year=float(r_year), corr_resid_share=r_share_resid),
               positive=dict(r=pc_r, floor=pc_floor, ceiling=pc_ceil, threshold=t, ok=bool(pc_ok)),
               placebo=dict(r=z_r, ok=bool(nc_ok)), reference=ref,
               verdict=verdict, seeds=SEEDS, n=len(d), unchallenged=True),
          open(OUT / "gss_gap_over_51_years.json", "w"), indent=1)
print(f"\nwrote {OUT/'gss_gap_over_51_years.json'}")
