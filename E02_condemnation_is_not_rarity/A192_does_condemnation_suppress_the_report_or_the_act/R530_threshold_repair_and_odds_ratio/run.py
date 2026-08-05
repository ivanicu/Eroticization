"""E02·A192·R530 — 门槛规则修好之后,三轮的 UNVERIFIED 里有几轮是我自伤的;以及基率不变地重问

`#484` 的 NEXT,两件,按顺序,分别标注行动类型:

  ① CLOSURE  修正门槛规则并**回溯**它污染了什么。不是 Frontier —— 它保护既有结论的可读性,
             不开新世界。诚实标注(frontier §7 规则 6)。
  ② FRONTIER 用**基率不变**的统计量、**在固定年龄带内**重问 `#484c`。

⚠ frontier §3 要求设计一个**正面结果我会不愿意看到**的步骤。这里两个都是:
  ①的不愿意结果 = 门槛修复让 `R527`/`R529` 的正对照通过 ⇒ **三轮 UNVERIFIED 是我自伤的**;
  ②的不愿意结果 = W-PREFERENCE(差距跨时代稳定)⇒ **我从 `#466c` 起隐隐偏好的隐瞒故事是空的**,
     而那个「我从没测过的漏洞」原来没东西。**两个都让它出来。**

------------------------------------------------------------------ ① CLOSURE
新规则(预注册,写在跑之前):
  正对照门槛 t := **同一仪器测出的参照分布的 q95**,而不是 `(floor+ceiling)/2`。
  理由(`#484b`):中点是我 `R527` 发明的,realstat 只要求 `floor < t < ceiling`;
  二值 × 多级的数学上界与「一个真实社会效应能有多大」无关。
回溯对象:`R527` 的 `SCCS1776×SCCS1777` · `R528` 的 `SCCS166×SCCS167` · `R529` 的 `谴责×attend`。
⛔ **不得追认那几轮的结论**(`#482e`):门槛不可用 ⇒ 问题**回到未决**,
   需要一次**重新预注册**的运行,而不是把旧输出重读一遍。旧输出的重读只作**标注为事后、不具约束力**。

------------------------------------------------------------------ ② FRONTIER
G1 ESTIMAND(先于方法):在**固定年龄带**内,每一年
  `lnOR_t = ln[ (P(saw|condemn)/(1-P)) / (P(saw|~condemn)/(1-P)) ]`,
  然后 `slope( |lnOR_t| ~ 该年谴责率 )`。
  **优势比对基率位移不变,风险差不是** —— `#484d` 量到 `corr(|RD|,基率)=+0.79`。
  ⚠ 本轮**先验证这一点**:若 `corr(|lnOR|,基率)` 不比 `corr(|RD|,基率)` 小得多,换统计量没用。

⚠ STRONGEST CONFOUND,写在跑之前:**世代 × 年份 ≈ 年龄**。
  同一世代随时间变老,而谴责与观看都随年龄变。`#484c` 的世代分层**没有控制年龄**。
  ⇒ 控制:**固定年龄带**(同一人生阶段,不同年代),这是本轮与 `#484c` 的关键差别。

WORLDS:
  W-CONCEAL    谴责规范强时谴责者更不敢承认 -> |lnOR| 随该年谴责率**上升**  -> slope > 0
  W-PREFERENCE 谴责反映真实不想要           -> |lnOR| 跨时代**稳定**        -> slope ≈ 0
  | World       | now | slope>0 | slope≈0 |
  | W-CONCEAL   | 0.4 | 0.85    | 0.05    |
  | W-PREFERENCE| 0.4 | 0.05    | 0.85    |
  | W-BOTH      | 0.2 | 0.10    | 0.10    |

CONTROLS:
  正对照  谴责 × 礼拜出席,门槛 = **GSS 参照分布 q95**(新规则)
  安慰剂  谴责 × 星座 —— 必须为零
  精度    **人层 bootstrap**(`#484e`:零臂是地板,不是误差棒)
  ⚠ 零式声明纪律(guard 21):若判 W-PREFERENCE,必须同时给
    置换分位 + **MDE** + 已演示的灵敏度,否则不是零而是沉默。

KILL(条件式,预注册):
  if 正对照触发 and 安慰剂为零:
      slope > 0 且 bootstrap CI 不含 0 -> W-CONCEAL
      |slope| < MDE 且 CI 含 0         -> W-PREFERENCE(带 MDE 陈述)
      slope < 0 且 CI 不含 0           -> 两个世界都没预测;MIXED
  else: UNVERIFIED

IMPOSSIBLE:无干预 ⇒ 非因果 · 未派对抗 agent(会话约束)⇒ `[unchallenged]` ·
  GSS 是重复横截面 ⇒ 个体内变化结构上不可测
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
E02 = ROOT / "E02_condemnation_is_not_rarity"

# ================================================================== ① CLOSURE
print("=" * 72)
print("① CLOSURE — 门槛规则 (floor+ceiling)/2  ->  同仪器参照分布 q95;并回溯")
print("=" * 72)
r527 = json.load(open(E02 / "A190_is_the_coupling_the_coder_or_the_norm/"
                      "R527_repair_both_controls/results/repair_both_controls.json"))
r528 = json.load(open(E02 / "A191_is_the_headline_special_inside_its_own_source/"
                      "R528_broude_reference_distribution/results/broude_reference_distribution.json"))
r529 = json.load(open(E02 / "A192_does_condemnation_suppress_the_report_or_the_act/"
                      "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))

retro = []
# R527: 正对照 SCCS1776×SCCS1777,参照 = lang 同源分布
retro.append(dict(round="R527", pair="SCCS1776×SCCS1777 (lang)",
                  observed=abs(r527["positive"]["rho"]),
                  t_old=r527["positive"]["threshold_new"],
                  t_new=r527["reference"]["q95"],
                  ref_k=r527["reference"]["k"], ref_median=r527["reference"]["median"]))
# R528: 正对照 SCCS166×SCCS167,参照 = broude 主格分布
retro.append(dict(round="R528", pair="SCCS166×SCCS167 (broude)",
                  observed=abs(r528["positive"]["rho"]),
                  t_old=r528["positive"]["threshold"],
                  t_new=r528["main"]["q95"],
                  ref_k=r528["main"]["n_pairs"], ref_median=r528["main"]["median"]))
# R529: 正对照 condemn×attend,参照 = 12 个无关 GSS 变量
ref529 = np.array([abs(x["r"]) for x in r529["reference"]])
retro.append(dict(round="R529", pair="condemn×attend (GSS)",
                  observed=abs(r529["positive"]["r"]),
                  t_old=r529["positive"]["threshold"],
                  t_new=float(np.quantile(ref529, .95)),
                  ref_k=len(ref529), ref_median=float(np.median(ref529))))

print(f"{'round':6s} {'pair':26s} {'观测':>8s} {'旧门槛(中点)':>13s} {'新门槛(参照q95)':>16s} "
      f"{'旧':>5s} {'新':>5s}")
flipped = []
for x in retro:
    old = "PASS" if x["observed"] > x["t_old"] else "FAIL"
    new = "PASS" if x["observed"] > x["t_new"] else "FAIL"
    x["verdict_old"], x["verdict_new"] = old, new
    if old != new: flipped.append(x["round"])
    print(f"{x['round']:6s} {x['pair']:26s} {x['observed']:8.4f} {x['t_old']:13.4f} "
          f"{x['t_new']:16.4f} {old:>5s} {new:>5s}   (参照 k={x['ref_k']}, 中位={x['ref_median']:.4f})")
print(f"\n⇒ 门槛规则改变了 **{len(flipped)}** 轮的正对照判定:{flipped or '无'}")
print("⛔ 但这**不追认**那几轮的结论:预注册的 kill 挂在一个现已知不可用的门槛上,")
print("   ⇒ 那些问题**回到未决**,需要重新预注册的运行。下面只作**事后、不具约束力**的标注:")
if "R527" in flipped:
    print(f"   · R527 事后:范围匹配 ρ 在同源分布 q{100*r527['main']['percentile_in_same_source']:.0f}"
          f"(其 kill 门槛为 q90)—— **不具约束力**")

# ================================================================== ② FRONTIER
print("\n" + "=" * 72)
print("② FRONTIER — 基率不变的统计量,固定年龄带内")
print("=" * 72)
cols = ["year", "age", "cohort", "degree", "attend", "zodiac", "pornlaw", "xmovie", "wtssps"]
it = pd.read_stata(DTA, iterator=True); have = set(it.variable_labels())
cols = [c for c in cols if c in have]
df = pd.read_stata(DTA, columns=cols, convert_categoricals=False)
check_columns(df, where="R530")
d = df.dropna(subset=["pornlaw", "xmovie", "age"]).copy()
d["condemn"] = (d.pornlaw == 1).astype(float)
d["saw"] = (d.xmovie == 1).astype(float)
print(f"n={len(d)}  年份 {int(d.year.min())}-{int(d.year.max())} ({d.year.nunique()} 个)  "
      f"年龄 {int(d.age.min())}-{int(d.age.max())}")


def stats_by_year(dd, lo=None, hi=None, cond_col="condemn"):
    if lo is not None: dd = dd[(dd.age >= lo) & (dd.age <= hi)]
    out = []
    for y, g in dd.groupby("year"):
        a, b = g[g[cond_col] == 1], g[g[cond_col] == 0]
        if len(a) < 30 or len(b) < 30: continue
        pa, pb = a.saw.mean(), b.saw.mean()
        if min(pa, pb) <= 0 or max(pa, pb) >= 1: continue
        lnor = math.log((pa / (1 - pa)) / (pb / (1 - pb)))
        q = ((pa / (1 - pa)) - (pb / (1 - pb))) / ((pa / (1 - pa)) + (pb / (1 - pb)))  # Yule's Q
        out.append(dict(year=int(y), lnor=lnor, yq=q, rd=pa - pb,
                        share=g[cond_col].mean(), base=g.saw.mean(), n=len(g)))
    return pd.DataFrame(out)


AGE_BANDS = [("18-34", 18, 34), ("25-54", 25, 54), ("35-64", 35, 64), ("all", 18, 99)]
print("\n=== 先验证前提:优势比是不是真的比风险差更基率不变 ===")
for nm, lo, hi in AGE_BANDS:
    bb = stats_by_year(d, lo, hi)
    if len(bb) < 8: print(f"  {nm}: 年份不足 ({len(bb)}),跳过"); continue
    c_rd = np.corrcoef(bb.base, np.abs(bb.rd))[0, 1]
    c_or = np.corrcoef(bb.base, np.abs(bb.lnor))[0, 1]
    print(f"  {nm:6s} years={len(bb):2d}  corr(|RD|,基率)={c_rd:+.4f}   "
          f"corr(|lnOR|,基率)={c_or:+.4f}   {'✅ 更不变' if abs(c_or)<abs(c_rd) else '⛔ 没更不变'}")

MAIN = ("25-54", 25, 54)
bb = stats_by_year(d, MAIN[1], MAIN[2])
print(f"\n=== 主格 年龄 {MAIN[0]},{len(bb)} 年全公布 ===")
print(bb.assign(lnor=bb.lnor.round(4), yq=bb.yq.round(4), rd=bb.rd.round(4),
                share=bb.share.round(3), base=bb.base.round(3)).to_string(index=False))


def slope(bb, col="lnor"):
    x = bb.share.values; y = np.abs(bb[col].values)
    return float(np.polyfit(x, y, 1)[0]) if len(x) >= 8 else np.nan


main_slope = slope(bb, "lnor")
# 精度:人层 bootstrap(`#484e`)
def boot_slope(dd, lo, hi, col="lnor", B=400, seed=0):
    rng = np.random.default_rng(seed); out = []
    sub = dd[(dd.age >= lo) & (dd.age <= hi)]
    idx = np.arange(len(sub))
    for _ in range(B):
        s = sub.iloc[rng.choice(idx, len(idx), replace=True)]
        b = stats_by_year(s)
        v = slope(b, col)
        if np.isfinite(v): out.append(v)
    return np.array(out)


boots = [boot_slope(d, MAIN[1], MAIN[2], "lnor", 400, s) for s in SEEDS]
ball = np.concatenate(boots)
lo_ci, hi_ci = np.quantile(ball, [.025, .975])
seed_spread = float(np.std([b.mean() for b in boots]))
print(f"\n主 slope(|lnOR| ~ 谴责率) = {main_slope:+.4f}")
print(f"  人层 bootstrap 95% CI = [{lo_ci:+.4f}, {hi_ci:+.4f}]  sd={ball.std():.4f}  "
      f"seed_spread={seed_spread:.5f}")
# MDE:能被这个设计检出的最小斜率(2×bootstrap sd,单边 80% 功效近似)
MDE = 2.8 * float(ball.std())
print(f"  MDE (2.8×bootstrap sd) = {MDE:.4f}")

print("\n=== G4 规格曲线(全格公布,含反号格)===")
spec = []
for nm, lo, hi in AGE_BANDS:
    for col in ("lnor", "yq", "rd"):
        b2 = stats_by_year(d, lo, hi)
        if len(b2) < 8: continue
        s = slope(b2, col)
        spec.append(dict(band=nm, stat=col, slope=s, years=len(b2)))
for cdef, name in [((d.pornlaw == 1), "illegal_all"), ((d.pornlaw < 3), "illegal_any")]:
    d2 = d.copy(); d2["condemn"] = cdef.astype(float)
    b2 = stats_by_year(d2, 25, 54)
    if len(b2) >= 8:
        spec.append(dict(band="25-54", stat=f"lnor/{name}", slope=slope(b2, "lnor"), years=len(b2)))
for s in spec:
    print(f"  {s['band']:6s} {s['stat']:16s} slope={s['slope']:+.4f}  years={s['years']}")
fin = [s["slope"] for s in spec if np.isfinite(s["slope"])]
sg = [np.sign(v) for v in fin]; dom = max(set(sg), key=sg.count)
print(f"\nspec_survival: {sg.count(dom)}/{len(fin)} = {sg.count(dom)/len(fin):.0%} 同号 ({dom:+.0f})")

# ---------------------------------------------------------------- 控制
G = Gate("谴责压低的是报告还是行为?(GSS,固定年龄带,优势比)")
pc = d.dropna(subset=["attend"])
pc_r = float(np.corrcoef(pc.condemn, pc.attend)[0, 1])
t_new = float(np.quantile(ref529, .95))
print(f"\n正对照 谴责×礼拜出席 r={pc_r:+.4f}  新门槛(GSS 参照 q95)={t_new:.4f}")
pc_ok = G.positive_control("正对照:门槛 = 同问卷参照分布 q95(新规则)",
                           planted=abs(pc_r), floor=t_new, spread=1e-9)
z = d.dropna(subset=["zodiac"])
z_r = float(np.corrcoef(z.condemn, z.zodiac)[0, 1])
rngz = np.random.default_rng(SEEDS[0])
z_null = np.array([np.corrcoef(z.condemn.values[rngz.permutation(len(z))], z.zodiac)[0, 1]
                   for _ in range(300)])
nc_ok = G.negative_control("安慰剂:谴责×星座(必须为零)", null=z_r, effect=pc_r,
                           null_spread=float(z_null.std()), null_kind="个体层标签置换")
G.has_error_bar("主斜率", value=main_slope, spread=float(ball.std()), spread_source="bootstrap_人层")

ci_excl_0 = not (lo_ci <= 0 <= hi_ci)
if pc_ok and nc_ok:
    if main_slope > 0 and ci_excl_0:
        verdict = "slope>0 且 CI 不含 0 -> W-CONCEAL"
    elif abs(main_slope) < MDE and not ci_excl_0:
        verdict = f"|slope|={abs(main_slope):.4f} < MDE={MDE:.4f} 且 CI 含 0 -> W-PREFERENCE(带 MDE)"
        G.null_claim_uses_null_criteria("零式声明:差距跨时代稳定", claim_kind="NULL",
                                        perm_quantile=float((ball >= abs(main_slope)).mean()),
                                        mde=MDE, sensitivity_shown=True, meaningful=0.5)
    elif main_slope < 0 and ci_excl_0:
        verdict = "slope<0 且 CI 不含 0 -> 两个世界都没预测;MIXED"
    else:
        verdict = f"|slope|={abs(main_slope):.4f} 与 MDE={MDE:.4f} 不可区分且 CI 含 0 -> 未决"
    print(f"\n控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:年龄带固定了人生阶段,但**没有固定接触机会** ——"
          " X 片的可得性在 51 年里剧变,而可得性变化对谴责者/非谴责者可能不对称。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} placebo={nc_ok})"
    print(f"\n⚠ {verdict}")
print(G)

json.dump(dict(closure=dict(retro=retro, flipped=flipped,
                            rule_old="(floor+ceiling)/2", rule_new="same-instrument reference q95"),
               frontier=dict(by_year=bb.to_dict("records"), slope=main_slope,
                             ci=[float(lo_ci), float(hi_ci)], boot_sd=float(ball.std()),
                             seed_spread=seed_spread, MDE=MDE, age_band=MAIN[0],
                             spec=[{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                                    for k, v in s.items()} for s in spec]),
               positive=dict(r=pc_r, threshold=t_new, ok=bool(pc_ok)),
               placebo=dict(r=z_r, ok=bool(nc_ok)),
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "threshold_repair_and_odds_ratio.json", "w"), indent=1)
print(f"\nwrote {OUT/'threshold_repair_and_odds_ratio.json'}")
