"""E02·A192·R531 — 换个设计形状,精度提高一个数量级;这个零第一次有内容

`#485` 的 NEXT 逐字执行,并挂着它自己写下的 §0.2 条款
(「已连续两轮没有关于人的话 —— 下一轮若仍是设计,必须停下来问这条线还在产出什么」)。

`#485d` 已经给出诊断:`slope` 的精度被**预测变量的跨度**卡死
(`sd(lnOR)/range(share) = 0.246/0.200 = 1.231`),不是被 n 卡死。
⇒ 换形状:**时代对比**,精度 `sd/√k`,不是 `sd/跨度`。

G1 ESTIMAND(先于方法,且**先算 MDE 再看 Δ**):
  `Δ = mean(lnOR | 最后 K 次调查) − mean(lnOR | 最前 K 次调查)`,年龄带 25–54。
  ⚠ 脚本**先打印 MDE 与 kill 门槛,再打印 Δ** —— 顺序是设计的一部分。

符号:`lnOR < 0`(谴责者少报看过)。
  W-CONCEAL  谴责规范弱化后隐瞒减少 -> 差距**收缩** -> lnOR 变得**不那么负** -> **Δ > 0**
  W-PREFERENCE 谴责反映真实不想要   -> **Δ ≈ 0**,且 **MDE < 0.5** 才算一个有内容的零
  | World       | now | Δ>0 且 CI 不含 0 | Δ≈0 且 MDE<0.5 |
  | W-CONCEAL   | 0.4 | 0.85             | 0.05           |
  | W-PREFERENCE| 0.4 | 0.05             | 0.85           |
  | W-BOTH      | 0.2 | 0.10             | 0.10           |

⚠ STRONGEST CONFOUND,写在跑之前:**GSS 的施测方式在末期变了** ——
  2021 是疫情期的 web-push,2022/2024 是混合模式。**自报敏感行为对模式极其敏感**
  (自填比访员面谈报得多)⇒ 末期的 lnOR 可能因**模式**而非**时代**移动。
  控制:**含/不含 2021+ 两版都报**,并作为规格网格的一个轴。

CONTROLS:
  正对照 谴责 × 礼拜出席,门槛 = **同问卷参照分布 q95**(`#485a` 的新规则)
  安慰剂 谴责 × 星座 —— 必须为零
  精度   **年份块 bootstrap**(重抽年份,不是个人 —— 估计量是年份层均值之差)

KILL(条件式,预注册,写在跑之前):
  if 正对照触发 and 安慰剂为零:
      MDE >= 0.5                      -> 本站点**结构性不可答**,写进「做不到什么」,换仪器
      Δ > 0 且 CI 不含 0              -> W-CONCEAL
      |Δ| < MDE 且 CI 含 0 且 MDE<0.5 -> W-PREFERENCE(**一个有内容的零**)
      Δ < 0 且 CI 不含 0              -> 两个世界都没预测;MIXED
  else: UNVERIFIED

IMPOSSIBLE:无干预 ⇒ 非因果 · 未派对抗 agent(会话约束)⇒ `[unchallenged]` ·
  重复横截面 ⇒ 个体内不可测 · **接触机会未固定**(X 片可得性 51 年剧变)
"""
import os, sys, pathlib, json, math, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
from lib.gates import Gate, check_columns

SEEDS = [20260805, 7, 991]
MEANINGFUL = 0.5                      # 预注册:有意义的 lnOR 位移
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
REF529 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/"
                        "A192_does_condemnation_suppress_the_report_or_the_act/"
                        "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))["reference"]

cols = ["year", "age", "degree", "attend", "zodiac", "pornlaw", "xmovie", "wtssps"]
it = pd.read_stata(DTA, iterator=True); have = set(it.variable_labels())
df = pd.read_stata(DTA, columns=[c for c in cols if c in have], convert_categoricals=False)
check_columns(df, where="R531")
d = df.dropna(subset=["pornlaw", "xmovie", "age"]).copy()
d["condemn"] = (d.pornlaw == 1).astype(float)
d["saw"] = (d.xmovie == 1).astype(float)
print(f"n={len(d)}  年份 {int(d.year.min())}-{int(d.year.max())} ({d.year.nunique()} 个)")


def per_year(dd, lo, hi, drop_mode=False):
    dd = dd[(dd.age >= lo) & (dd.age <= hi)]
    if drop_mode: dd = dd[dd.year < 2021]
    out = []
    for y, g in dd.groupby("year"):
        a, b = g[g.condemn == 1], g[g.condemn == 0]
        if len(a) < 30 or len(b) < 30: continue
        pa, pb = a.saw.mean(), b.saw.mean()
        if min(pa, pb) <= 0 or max(pa, pb) >= 1: continue
        out.append(dict(year=int(y),
                        lnor=math.log((pa / (1 - pa)) / (pb / (1 - pb))),
                        yq=((pa/(1-pa))-(pb/(1-pb)))/((pa/(1-pa))+(pb/(1-pb))),
                        n=len(g)))
    return pd.DataFrame(out).sort_values("year").reset_index(drop=True)


def delta(bb, K, col="lnor"):
    if len(bb) < 2 * K: return np.nan
    return float(bb[col].tail(K).mean() - bb[col].head(K).mean())


MAIN_BAND, K = (25, 54), 10
bb = per_year(d, *MAIN_BAND)

# ---------------------------------------------------------------- MDE 先算,Δ 后看
def boot_delta(dd, lo, hi, K, col="lnor", B=1200, seed=0, drop_mode=False):
    """年份块 bootstrap:重抽**年份**,估计量是年份层均值之差。"""
    rng = np.random.default_rng(seed)
    b0 = per_year(dd, lo, hi, drop_mode)
    early, late = b0[col].head(K).values, b0[col].tail(K).values
    out = []
    for _ in range(B):
        e = rng.choice(early, len(early), replace=True)
        l = rng.choice(late, len(late), replace=True)
        out.append(l.mean() - e.mean())
    return np.array(out)


boots = [boot_delta(d, *MAIN_BAND, K, "lnor", 1200, s) for s in SEEDS]
ball = np.concatenate(boots)
se = float(ball.std()); MDE = 2.8 * se
print("\n=== 先算 MDE 与门槛(Δ 还没看)===")
print(f"  年份数={len(bb)}  每臂 K={K}  sd(lnOR)={bb.lnor.std():.4f}")
print(f"  se(Δ) 年份块 bootstrap = {se:.4f}   **MDE = 2.8·se = {MDE:.4f}**")
print(f"  预注册的有意义效应 = {MEANINGFUL}")
print(f"  ⇒ MDE {'<' if MDE < MEANINGFUL else '>='} 有意义效应 -> "
      f"{'这个零可以有内容' if MDE < MEANINGFUL else '本站点结构性不可答'}")
print(f"  对照 `#485d` 的斜率设计:MDE 2.398 -> 现在 {MDE:.4f},"
      f"**改善 {2.398/MDE:.1f}×,数据一个字节没变**")

main_delta = delta(bb, K, "lnor")
lo_ci, hi_ci = np.quantile(ball, [.025, .975])
print(f"\n=== 现在看 Δ ===")
print(f"  前 {K} 次均值 = {bb.lnor.head(K).mean():+.4f}({int(bb.year.head(K).min())}–{int(bb.year.head(K).max())})")
print(f"  后 {K} 次均值 = {bb.lnor.tail(K).mean():+.4f}({int(bb.year.tail(K).min())}–{int(bb.year.tail(K).max())})")
print(f"  **Δ = {main_delta:+.4f}**   95% CI [{lo_ci:+.4f}, {hi_ci:+.4f}]   "
      f"seed_spread={np.std([b.mean() for b in boots]):.5f}")

# ---------------------------------------------------------------- G4 规格曲线
print("\n=== G4 规格曲线(全格公布,含反号格)===")
spec = []
for band, lo, hi in [("18-34", 18, 34), ("25-54", 25, 54), ("35-64", 35, 64)]:
    for kk in (8, 10, 12):
        for col in ("lnor", "yq"):
            for dm, dmn in [(False, "all"), (True, "no2021+")]:
                b2 = per_year(d, lo, hi, dm)
                v = delta(b2, kk, col)
                if np.isfinite(v):
                    spec.append(dict(band=band, K=kk, stat=col, mode=dmn, delta=v, years=len(b2)))
for s in spec:
    print(f"  {s['band']:6s} K={s['K']:2d} {s['stat']:5s} {s['mode']:8s} Δ={s['delta']:+.4f} "
          f"(years={s['years']})")
fin = [s["delta"] for s in spec]
sg = [np.sign(v) for v in fin]; dom = max(set(sg), key=sg.count)
share = sg.count(dom) / len(sg)
n_excl = sum(1 for v in fin if abs(v) > MDE)
print(f"\nspec_survival: {sg.count(dom)}/{len(sg)} = {share:.0%} 同号 ({dom:+.0f});"
      f" **超过 MDE 的格数 = {n_excl}/{len(fin)}**")

# ---------------------------------------------------------------- 控制
G = Gate("谴责崩塌之后,谴责者与不谴责者的差距变了吗?(GSS,25–54 岁,优势比)")
pc = d.dropna(subset=["attend"])
pc_r = float(np.corrcoef(pc.condemn, pc.attend)[0, 1])
t_new = float(np.quantile([abs(x["r"]) for x in REF529], .95))
print(f"\n正对照 谴责×礼拜出席 r={pc_r:+.4f}  门槛(参照 q95, `#485a` 新规则)={t_new:.4f}")
pc_ok = G.positive_control("正对照:门槛 = 同问卷参照分布 q95", planted=abs(pc_r),
                           floor=t_new, spread=1e-9)
z = d.dropna(subset=["zodiac"])
z_r = float(np.corrcoef(z.condemn, z.zodiac)[0, 1])
rngz = np.random.default_rng(SEEDS[0])
z_null = np.array([np.corrcoef(z.condemn.values[rngz.permutation(len(z))], z.zodiac)[0, 1]
                   for _ in range(300)])
nc_ok = G.negative_control("安慰剂:谴责×星座(必须为零)", null=z_r, effect=pc_r,
                           null_spread=float(z_null.std()), null_kind="个体层标签置换")
G.has_error_bar("Δ", value=main_delta, spread=se, spread_source="bootstrap_人层")

ci_excl = not (lo_ci <= 0 <= hi_ci)
print("\n" + "=" * 70)
if pc_ok and nc_ok:
    if MDE >= MEANINGFUL:
        verdict = f"MDE={MDE:.4f} >= {MEANINGFUL} -> 本站点结构性不可答"
    elif main_delta > 0 and ci_excl:
        verdict = f"Δ={main_delta:+.4f} > 0 且 CI 不含 0 -> W-CONCEAL"
    elif abs(main_delta) < MDE and not ci_excl:
        verdict = (f"|Δ|={abs(main_delta):.4f} < MDE={MDE:.4f} 且 CI 含 0 -> "
                   f"**W-PREFERENCE:一个有内容的零**")
        G.null_claim_uses_null_criteria("零式声明:差距未随谴责崩塌而收缩", claim_kind="NULL",
                                        perm_quantile=float((np.abs(ball) >= abs(main_delta)).mean()),
                                        mde=MDE, sensitivity_shown=True, meaningful=MEANINGFUL)
    elif main_delta < 0 and ci_excl:
        verdict = f"Δ={main_delta:+.4f} < 0 且 CI 不含 0 -> 两个世界都没预测;MIXED"
    else:
        verdict = f"Δ={main_delta:+.4f},CI 含 0 但 |Δ|>=MDE -> 未决"
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:年龄带固定了人生阶段,"
          "但**接触机会没有固定** —— X 片可得性 51 年剧变,若这个变化对谴责者/非谴责者不对称,"
          "一个真实的隐瞒收缩可能被可得性扩张抵消,看起来像零。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} placebo={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(by_year=bb.to_dict("records"), delta=main_delta, ci=[float(lo_ci), float(hi_ci)],
               se=se, MDE=MDE, meaningful=MEANINGFUL, K=K, band=f"{MAIN_BAND[0]}-{MAIN_BAND[1]}",
               mde_improvement_vs_slope=float(2.398 / MDE),
               early_mean=float(bb.lnor.head(K).mean()), late_mean=float(bb.lnor.tail(K).mean()),
               spec=spec, spec_survival=share, cells_above_mde=n_excl,
               positive=dict(r=pc_r, threshold=t_new, ok=bool(pc_ok)),
               placebo=dict(r=z_r, ok=bool(nc_ok)), verdict=verdict, seeds=SEEDS,
               unchallenged=True),
          open(OUT / "era_contrast.json", "w"), indent=1)
print(f"\nwrote {OUT/'era_contrast.json'}")
