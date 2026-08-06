"""E02·A212·R576 — 一个时代的宽容,是同时松开所有做法,还是一件一件松开?

`#531` 的 NEXT。行动类型:**FRONTIER**。换仪器:SCCS(社会)→ GSS(年代),
检验 `#529`「严厉附着在做法上,不附着在社会上」在**另一个单位、另一份仪器**上是否成立(硬规则 4)。

⚠ 硬规则 1,已先做,四条序列的年份覆盖**并不相同**:
   homosex n=44,644 / 30 年 (1973–2024) · premarsx n=45,697 / 30 年 (1972–2024) ·
   xmarsex n=46,266 / 30 年 (1973–2024) · teensex n=33,901 / **22 年 (1986–2024)** ·
   pornlaw n=46,245 / 30 年。**premarsx 与 homosex 的年份表不同**(1972/75/78/83/86 vs 1973/76/77/80/…)。

⚠ 硬规则 2 + 最强混淆,写在跑之前:**四条序列都在长期自由化。**
   两条同向趋势的序列,水平相关**必然虚高**(Yule 的伪相关)。
   所以**水平相关不是证据**,本轮把它算出来只是为了展示那个陷阱;
   真正的量是**一阶差分**:在同一段时间里,它们是否**一起动**。

G1 ESTIMAND(先于方法):
   ① 水平:`ρ_lev = corr(谴责占比_A(t), 谴责占比_B(t))`(**陷阱量,不作证据**)
   ② **差分:`ρ_dif = corr(Δ谴责占比_A, Δ谴责占比_B)`,Δ 取相邻共同调查年之间。**
   ③ **噪声地板**:每年占比的二项抽样误差 `sqrt(p(1-p)/n_year)`,差分的地板 = `sqrt(2)×` 它;
      用它对 `ρ_dif` 做**衰减校正**,并报校正前后两个数。

WORLDS(本体不同):
  W-CLIMATE     `ρ_dif` 显著为正 ⇒ **一个时代的宽容同时松开多种做法** —— 与 `#529` 相反
  W-ONE-BY-ONE  `ρ_dif` ≈ 0 ⇒ **一件一件松开** —— `#529` 在新单位上复制
  W-TREND-ONLY  `ρ_lev` 高而 `ρ_dif` ≈ 0 ⇒ **「共同气候」的外观完全由共同趋势制造**
⚠ `W-TREND-ONLY` 是**元分离**:它说的不是哪个世界赢,而是**「水平/差分」这个区分本身
   才是问题所在** —— 我原来的世界分解(气候 vs 逐件)在水平上根本无法被区分。
⚠ BASIN:`W-ONE-BY-ONE` 会让 `#529` 漂亮地复制,所以它**不是**本轮下注方向。
   本轮下注 `W-CLIMATE` —— 它一旦为真,`#529` 就只是**关于社会这个单位的**,不是关于严厉的。

CONTROLS(G2):
  正对照 **同一道题**的男性子样本 Δ × 女性子样本 Δ —— 同一现实的两个独立测量,必须强正,
     且它给出**这具仪器在差分上能达到的上限**(不是 1.0);
  安慰剂 **打乱年份顺序**后重算 Δ 相关(该是零 ⇒ negative_control);
  噪声地板 二项 SE,measured 不是 assumed;
  规格曲线 谴责切点 {1} / {1,2} / {1,2,3} × 六对题 × 校正前后。
KILL(条件式,预注册):
  if 正对照通过 and 打乱年份后 ≈ 0:
      `ρ_dif` 中位 > 正对照的一半 -> W-CLIMATE
      `ρ_dif` 中位 ≈ 0 且 `ρ_lev` 中位 > 0.8 -> W-TREND-ONLY
      `ρ_dif` 中位 ≈ 0 且 `ρ_lev` 也不高 -> W-ONE-BY-ONE
      else UNVERIFIED-by-power
  else: UNVERIFIED
IMPOSSIBLE:一国一仪器 ⇒ 无跨国复制 · 观察性 ⇒ 非因果 · 共同年份最多 30 点 ⇒ 分辨率有限 ·
  同一份问卷同一批受访者 ⇒ **共同方法方差会抬高 ρ_dif**,即对 `W-CLIMATE` 方向**不保守** · [unchallenged]
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
ITEMS = {"homosex": [1, 2, 3, 4], "premarsx": [1, 2, 3, 4],
         "xmarsex": [1, 2, 3, 4], "teensex": [1, 2, 3, 4]}
g = pd.read_stata(DTA, columns=["year", "sex"] + list(ITEMS), convert_categoricals=False)

def series(col, cut, sub=None):
    d = g[g[col].isin(ITEMS[col])]
    if sub is not None: d = d[d.sex == sub]
    o = d.groupby("year")[col].agg(lambda v: float(np.isin(v, cut).mean()))
    n = d.groupby("year")[col].size()
    return o, n

CUTS = {"最严 {1}": [1], "中 {1,2}": [1, 2], "最宽 {1,2,3}": [1, 2, 3]}
rows = []
print("=== 规则①:每一格先打共同年数与每年最小 n,再看 ρ ===")
for cname, cut in CUTS.items():
    S = {k: series(k, cut) for k in ITEMS}
    for a, b in itertools.combinations(ITEMS, 2):
        ya = set(S[a][0].index) & set(S[b][0].index)
        yr = np.array(sorted(ya))
        if len(yr) < 8: continue
        pa, pb = S[a][0].reindex(yr).values, S[b][0].reindex(yr).values
        na, nb = S[a][1].reindex(yr).values, S[b][1].reindex(yr).values
        lev = float(np.corrcoef(pa, pb)[0, 1])
        da, db = np.diff(pa), np.diff(pb)
        dif = float(np.corrcoef(da, db)[0, 1])
        # 噪声地板:二项 SE,差分 = sqrt(SE_t^2 + SE_{t-1}^2)
        sea = np.sqrt(pa * (1 - pa) / na); seb = np.sqrt(pb * (1 - pb) / nb)
        fa = np.sqrt(sea[1:] ** 2 + sea[:-1] ** 2); fb = np.sqrt(seb[1:] ** 2 + seb[:-1] ** 2)
        rel_a = max(0.0, 1 - float(np.mean(fa ** 2) / np.var(da)))
        rel_b = max(0.0, 1 - float(np.mean(fb ** 2) / np.var(db)))
        dis = dif / np.sqrt(rel_a * rel_b) if rel_a * rel_b > 0 else np.nan
        rows.append(dict(cut=cname, pair=f"{a}×{b}", n_years=len(yr), n_diff=len(da),
                         min_n=int(min(na.min(), nb.min())), rho_level=lev, rho_diff=dif,
                         reliability=[round(rel_a, 4), round(rel_b, 4)],
                         rho_diff_disattenuated=float(dis) if np.isfinite(dis) else None,
                         inclusion=[f"两题共同调查年 {len(yr)} 个", cname,
                                    "Δ 取相邻共同年之差", "噪声地板 = 二项 SE"]))
        print(f"  {cname:12s} {a:9s}×{b:9s} 共同年={len(yr):2d} 最小n={int(min(na.min(),nb.min())):4d}  "
              f"水平ρ={lev:+.4f}  **差分ρ={dif:+.4f}**  信度=[{rel_a:.2f},{rel_b:.2f}]  "
              f"校正后={dis:+.4f}" if np.isfinite(dis) else "")

lev_med = float(np.median([r["rho_level"] for r in rows]))
dif_med = float(np.median([r["rho_diff"] for r in rows]))
print(f"\n  **水平 ρ 中位 = {lev_med:+.4f}(陷阱量)   差分 ρ 中位 = {dif_med:+.4f}(证据量)**")

# ---- 对照
G = Gate("一个时代的宽容,是同时松开所有做法,还是一件一件松开?(GSS 1972-2024)")
pcs = []
for k in ITEMS:
    m, f = series(k, [1, 2], 1)[0], series(k, [1, 2], 2)[0]
    yr = np.array(sorted(set(m.index) & set(f.index)))
    if len(yr) < 8: continue
    pcs.append(float(np.corrcoef(np.diff(m.reindex(yr).values), np.diff(f.reindex(yr).values))[0, 1]))
PC = float(np.median(pcs))
print(f"\n=== 对照 ===\n  正对照:同一道题的 男Δ × 女Δ,{len(pcs)} 题,中位 ρ = {PC:+.4f} "
      f"(这是本仪器在差分上的**上限**,不是 1.0)")
G.positive_control("正对照:同题男Δ×女Δ(仪器差分上限)", planted=abs(PC), floor=abs(dif_med), spread=1e-9)
shuf = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for cname, cut in CUTS.items():
        S = {k: series(k, cut) for k in ITEMS}
        for a, b in itertools.combinations(ITEMS, 2):
            yr = np.array(sorted(set(S[a][0].index) & set(S[b][0].index)))
            if len(yr) < 8: continue
            pa = S[a][0].reindex(yr).values; pb = S[b][0].reindex(yr).values[rng.permutation(len(yr))]
            shuf.append(abs(float(np.corrcoef(np.diff(pa), np.diff(pb))[0, 1])))
G.negative_control("安慰剂:打乱年份顺序后的 Δ 相关", null=float(np.median(shuf)), effect=abs(PC),
                   null_spread=float(np.std(shuf)), null_kind="年份顺序置换(破坏时间配对)")
G.spec_curve_cells_declare_n("规格曲线逐格 n", {f"{r['cut']}|{r['pair']}": dict(n=r["n_diff"], **r) for r in rows})
G.spec_curve_cells_declare_inclusion("规格曲线逐格纳入条件", {f"{r['cut']}|{r['pair']}": r for r in rows})

print("\n" + "=" * 76)
if abs(PC) > abs(dif_med) and np.median(shuf) < 0.5 * abs(PC):
    if dif_med > 0.5 * abs(PC):
        world = "W-CLIMATE"; verdict = f"差分中位 {dif_med:+.4f} > 正对照上限 {PC:+.4f} 的一半 -> **同时松开**"
    elif abs(dif_med) < np.quantile(shuf, .95) and lev_med > 0.8:
        world = "W-TREND-ONLY"; verdict = (f"水平中位 {lev_med:+.4f} 高而差分中位 {dif_med:+.4f} "
            f"落在打乱年份的 q95={np.quantile(shuf,.95):.4f} 以内 -> "
            f"**「共同气候」的外观完全由共同趋势制造**")
    elif abs(dif_med) < np.quantile(shuf, .95):
        world = "W-ONE-BY-ONE"; verdict = f"差分中位 {dif_med:+.4f} 落在安慰剂内,水平也不高 -> **一件一件松开**"
    else:
        world = "UNVERIFIED"; verdict = f"差分中位 {dif_med:+.4f} 介于安慰剂与上限之间 -> UNVERIFIED-by-power"
    print(f"控制齐备 ⇒ 评判。{world}:{verdict}")
    print("⚠ 这个 KILL 会怎样失败:四道题问的是**同一批受访者**在**同一份问卷**里,"
          "共同方法方差会**抬高** ρ_dif —— 也就是说它对 `W-CLIMATE` 方向**不保守**,"
          "若结论是「不同时松开」,那是在一个偏向反面的设计上得到的,更强;反之则更弱。")
else:
    world, verdict = "UNVERIFIED", "控制未齐"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, level_median=lev_med, diff_median=dif_med, world=world, verdict=verdict,
               positive_control=dict(per_item=pcs, median=PC), placebo_median=float(np.median(shuf)),
               placebo_q95=float(np.quantile(shuf, .95)), seeds=SEEDS,
               instrument="GSS 1972-2024,同一份问卷同一批受访者",
               impossible=["一国一仪器无跨国复制", "观察性非因果", "共同年份<=30 点",
                           "共同方法方差抬高 rho_diff,对 W-CLIMATE 不保守"], unchallenged=True),
          open(OUT / "levels_vs_differences.json", "w"), indent=1)
print(f"\nwrote {OUT/'levels_vs_differences.json'}")
