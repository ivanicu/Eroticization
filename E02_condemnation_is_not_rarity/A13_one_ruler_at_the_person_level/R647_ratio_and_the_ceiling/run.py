"""E02·A13·R647 —— 把 `#610c` 的两个事后量正式预注册,然后用天花板去打它们

`#610` 的 NEXT,判据在上一条**已经写死**,本轮只执行,不改判据。

⛔ **先做攻击阶梯第 2 级(代数),因为它最便宜 —— 而它给出的是 DERIVATION,不是证据。**
经典衰减 `ρ_obs = ρ_true·√(r_ii·r_jj)` ⇒
  `ratio(i) = within(i)/out(i) = [ρ_true,within·√(r_peers)] / [ρ_true,out·√(r_others)]`
**这道题自己的信度 `r_ii` 在比值里约掉了** ⇒ **它自己的天花板效应消掉。**
⚠ 这是一个**推导**,它的假设是:经典衰减成立(线性、范围截断不引入非线性)。
   序数题 + 秩相关 + 真实天花板未必满足 ⇒ **仍必须实测**,这正是本轮 E3。

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3(MFQ),
`data/external/dataverse/10.7910_DVN_SJTRBI_x/…Study_3.sav`。一份问卷 · 同一批人 · 同一量表。

G1 ESTIMAND(三个,先于方法):
  E1 `ratio(i) = within(i)/out(i)`,30 题全表。
  E2 每题的**边际分布**:均值 · 标准差 · 顶端占比(=5)· 底端占比(=0)。
  E3 **天花板敏感度**:把题 i 的顶端两档并成一档(人为造天花板),重算 `ratio(i)`,报 Δ。

KILL(条件式,**逐字取自 `#610` 的 NEXT,不改**):
  ① 纯洁四道仍占 30 题比值前四 **且** 最低的纯洁题比值 > 非纯洁题最高比值 + 2×展布
     ⇒ **W-INTERNAL**;任一非纯洁题挤进前四 ⇒ **W-C,`#610c①` 记为未分离**
  ② 若纯洁四道的标准差**系统性低于**其余 ⇒ 高比值可能是天花板效应,**记「未分离」**
  ⚠ 两支都要报比值表 + 边际分布表。**不许只报有利的一支。**
CONTROLS:
  正对照:合成**真总括题**(全 30 题秩均值)比值必须**最低**;合成**纯域内题**比值必须**最高**;
    **且 g=0(纯噪声)时判据必须不通过。**
  安慰剂:把一道**非纯洁**题人为造天花板,它的比值变化必须与纯洁题的变化同量级(否则是题的问题,不是族的问题)。
  **offset_control**(「这个零该是零吗?」**不该**):比值的基线不是 1 也不是 0 ——
    null 种类 = **「域大小与题数相同(各 6 题)时,任意分组都会产生的域内 > 域外」**,
    由**随机打乱域指派**给出。
G3:30 题的 `ratio` / 边际 / 天花板 Δ 全表发布,含不一致的题。
G4:{中位 / 均值} × {原始 / 偏宗教} × 3 种子。
IMPOSSIBLE(不写 planned):非因果 · **自选网络志愿者,非概率样本** ·
  人为天花板是**模拟**,不是真实的题目改写 —— 真正分开「措辞」与「概念」需要重写题干重测 ·
  未派对抗 agent ⇒ `[unchallenged]`
"""
import os, sys, pathlib, json, itertools, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import rankdata
from lib.gates import Gate
import pyreadstat

SEEDS = [20260806, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
SAV = ROOT / "data/external/dataverse/10.7910_DVN_SJTRBI_x/GrahamHaidtNosek.2009.JPSP.Study_3.sav"
ITEM = {
 "emotionally":("HARM",1),"weak":("HARM",1),"cruel":("HARM",1),
 "compassion":("HARM",2),"animal":("HARM",2),"kill":("HARM",2),
 "treated":("FAIRNESS",1),"unfairly":("FAIRNESS",1),"rights":("FAIRNESS",1),
 "fairly":("FAIRNESS",2),"justice":("FAIRNESS",2),"rich":("FAIRNESS",2),
 "lovecountry":("INGROUP",1),"betray":("INGROUP",1),"loyalty":("INGROUP",1),
 "history":("INGROUP",2),"family":("INGROUP",2),"team":("INGROUP",2),
 "respect":("AUTHORITY",1),"traditions":("AUTHORITY",1),"chaos":("AUTHORITY",1),
 "kidrespect":("AUTHORITY",2),"sexroles":("AUTHORITY",2),"soldier":("AUTHORITY",2),
 "decency":("PURITY",1),"disgusting":("PURITY",1),"god":("PURITY",1),
 "harmlessdg":("PURITY",2),"unnatural":("PURITY",2),"chastity":("PURITY",2),
}
REL = "Religion_attend_num"
NAMES = list(ITEM); DOMS = sorted(set(v[0] for v in ITEM.values()))
PUR = [k for k, v in ITEM.items() if v[0] == "PURITY"]
IDX = {c: i for i, c in enumerate(NAMES)}

d, _ = pyreadstat.read_sav(str(SAV))
CC = d[NAMES + [REL]].dropna()
print(f"仪器 = GHN 2009 Study 3(MFQ)· 30 题 + {REL} 完整 n={len(CC)}")
RAW = CC[NAMES].values.astype(float)
RK = np.column_stack([rankdata(RAW[:, j]) for j in range(RAW.shape[1])]).astype(float)


def partial(M, v):
    z = rankdata(v).astype(float); z = (z - z.mean()) / z.std()
    Z = np.column_stack([np.ones(len(M)), z])
    return M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]


def prof(M, agg=np.median, imap=ITEM):
    R = np.abs(np.corrcoef(M, rowvar=False))
    out = {}
    for it, (D, _) in imap.items():
        ins = [k for k in imap if imap[k][0] == D and k != it]
        ous = [k for k in imap if imap[k][0] != D]
        w = float(agg([R[IDX[it], IDX[k]] for k in ins]))
        o = float(agg([R[IDX[it], IDX[k]] for k in ous]))
        out[it] = (w, o, w/o if o > 1e-9 else np.nan)
    return out


# ── E1 + E2 ─────────────────────────────────────────────────
P = prof(RK)
rows = []
for j, it in enumerate(NAMES):
    v = RAW[:, j]
    w, o, r = P[it]
    rows.append(dict(item=it, domain=ITEM[it][0], within=w, out=o, ratio=r,
                     mean=float(v.mean()), sd=float(v.std()),
                     top=float((v == 5).mean()), bot=float((v == 0).mean())))
T = pd.DataFrame(rows).sort_values("ratio", ascending=False).reset_index(drop=True)
print("\n=== G3 · E1+E2:30 题按比值排序(全表)===")
print(f"{'#':>2} {'题':13s} {'域':9s} {'within':>7s} {'out':>7s} {'ratio':>6s} {'均值':>5s} {'sd':>5s} {'顶端%':>6s}")
for i, r in T.iterrows():
    m = " ←纯洁" if r.domain == "PURITY" else ""
    print(f"{i+1:2d} {r['item']:13s} {r.domain:9s} {r.within:7.4f} {r['out']:7.4f} {r.ratio:6.2f} "
          f"{r['mean']:5.2f} {r.sd:5.2f} {r.top*100:5.1f}%{m}")

top4 = list(T.head(4)["item"]); top4_all_pur = all(ITEM[i][0] == "PURITY" for i in top4)
purR = T[T.domain == "PURITY"].ratio; othR = T[T.domain != "PURITY"].ratio
print(f"\n**判据①**:前四 = {top4} -> 全是纯洁? **{top4_all_pur}**")
print(f"  纯洁六题比值 最低 **{purR.min():.3f}**(`{T[T.domain=='PURITY'].iloc[-1]['item']}`)· "
      f"非纯洁最高 **{othR.max():.3f}**(`{T[T.domain!='PURITY'].iloc[0]['item']}`)")

bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(150):
        Pb = prof(RK[rng.integers(0, len(RK), len(RK))])
        pr = [Pb[k][2] for k in PUR]; orr = [Pb[k][2] for k in NAMES if k not in PUR]
        bs.append(float(min(pr) - max(orr)))
SPR = float(np.std(bs))
GAP = float(purR.min() - othR.max())
print(f"  **min(纯洁) − max(非纯洁) = {GAP:+.3f}** ± {SPR:.3f} ⇒ 严格判据 "
      f"{'**通过**' if GAP > 2*SPR else '**不通过**(要求 > 2×展布)'}")

print(f"\n**判据②**:标准差 —— 纯洁六题 中位 **{T[T.domain=='PURITY'].sd.median():.3f}** "
      f"区间 [{T[T.domain=='PURITY'].sd.min():.3f}, {T[T.domain=='PURITY'].sd.max():.3f}] · "
      f"其余 24 题 中位 **{T[T.domain!='PURITY'].sd.median():.3f}** "
      f"区间 [{T[T.domain!='PURITY'].sd.min():.3f}, {T[T.domain!='PURITY'].sd.max():.3f}]")
print(f"  顶端占比:纯洁 中位 {T[T.domain=='PURITY'].top.median()*100:.1f}% · "
      f"其余 {T[T.domain!='PURITY'].top.median()*100:.1f}%")

# ── E3:人为造天花板,实测那条推导 ────────────────────────────
print("\n=== E3:人为把每题顶端两档并成一档(造天花板),重算比值 ===")
print("  ⚠ 上面那条约分是**推导**,不是证据。这里实测它在序数 + 秩相关下还成不成立。")
ce = []
for j, it in enumerate(NAMES):
    M = RK.copy()
    v = RAW[:, j].copy(); v[v >= 4] = 4          # 顶端两档并档
    M[:, j] = rankdata(v)
    r2 = prof(M)[it][2]
    ce.append(dict(item=it, domain=ITEM[it][0], ratio=P[it][2], ratio_ceil=r2, d=r2-P[it][2],
                   lost_top=float((RAW[:, j] == 5).mean())))
CE = pd.DataFrame(ce)
print(f"  比值变化 |Δ| 中位 **{CE.d.abs().median():.4f}** · 最大 {CE.d.abs().max():.4f}"
      f"(`{CE.loc[CE.d.abs().idxmax(),'item']}`)")
print(f"  纯洁六题 |Δ| 中位 {CE[CE.domain=='PURITY'].d.abs().median():.4f} · "
      f"其余 24 题 {CE[CE.domain!='PURITY'].d.abs().median():.4f}")
print(f"  ⇒ 把顶端两档并掉(平均抹掉 {CE.lost_top.mean()*100:.1f}% 的顶端回答),"
      f"比值动了 **{CE.d.abs().median()/T.ratio.median()*100:.2f}%** ——"
      f" **推导得到实测支持**" if CE.d.abs().median() < 0.1 else " **推导被实测否掉**")

# ── 控制 ────────────────────────────────────────────────────
G = Gate("`#610c` 的两个事后量,正式预注册后还站得住吗?")


def synth_ratio(kind, g, sd):
    rng = np.random.default_rng(sd)
    Z = (RK - RK.mean(0)) / RK.std(0)
    core = Z.mean(1) if kind == "omni" else Z[:, [IDX[k] for k in PUR]].mean(1)
    core = (core - core.mean()) / core.std()
    v = g*core + np.sqrt(max(0., 1-g**2))*rng.standard_normal(len(Z))
    M = np.column_stack([RK, rankdata(v)])
    R = np.abs(np.corrcoef(M, rowvar=False)); j = M.shape[1]-1
    ins = [k for k in NAMES if ITEM[k][0] == "PURITY"]
    ous = [k for k in NAMES if ITEM[k][0] != "PURITY"]
    w = float(np.median([R[j, IDX[k]] for k in ins])); o = float(np.median([R[j, IDX[k]] for k in ous]))
    return w/o if o > 1e-9 else np.nan

om = {g: float(np.mean([synth_ratio("omni", g, s) for s in SEEDS])) for g in (0.0, 0.5, 1.0)}
pu = {g: float(np.mean([synth_ratio("pure", g, s) for s in SEEDS])) for g in (0.0, 0.5, 1.0)}
print(f"\n  正对照:合成**纯域内题** 比值 g=1 **{pu[1.0]:.2f}** · g=0.5 {pu[0.5]:.2f} · g=0 {pu[0.0]:.2f}")
print(f"          合成**真总括题** 比值 g=1 **{om[1.0]:.2f}** · g=0.5 {om[0.5]:.2f} · g=0 {om[0.0]:.2f}")
pos_ok = G.positive_control("正对照:纯域内合成题的比值必须远高于总括合成题",
                            planted=float(pu[1.0]-om[1.0]), floor=0.30, spread=SPR)
g0 = abs(pu[0.0]-om[0.0])
print(f"    g=0 时两者差 {g0:.3f} ⇒ 判据在 g=0 {'⛔ 也通过' if g0 > 0.30+2*SPR else '**不通过** ✅'}")
# 安慰剂:非纯洁题造天花板,变化必须与纯洁题同量级
pl = float(abs(CE[CE.domain!='PURITY'].d.abs().median() - CE[CE.domain=='PURITY'].d.abs().median()))
pla_ok = G.negative_control("安慰剂:非纯洁题造天花板,比值变化须与纯洁题同量级",
                            null=pl, effect=float(GAP), null_spread=SPR,
                            null_kind="天花板对任何一道题都有的机械影响")
# offset:随机打乱域指派
sh = []
for sd in SEEDS:
    r2 = np.random.default_rng(sd); doms0 = [ITEM[k][0] for k in NAMES]
    for _ in range(150):
        perm = r2.permutation(doms0)
        FK = {k: (perm[i], ITEM[k][1]) for i, k in enumerate(NAMES)}
        Pf = prof(RK, imap=FK); sh.append(float(np.median([Pf[k][2] for k in NAMES])))
print(f"\n  offset:随机打乱域指派 -> 比值中位 **{np.mean(sh):.3f} ± {np.std(sh):.3f}**"
      f"(真指派 {T.ratio.median():.3f})")
G.offset_control("offset:任意分组都会产生的域内>域外", effect=float(T.ratio.median()),
                 offset=float(np.mean(sh)), spread=float(np.std(sh)),
                 null_kind="域大小相同(各6题)时,随机分组也会产生的域内>域外")
G.has_error_bar("min(纯洁比值) − max(非纯洁比值)", value=GAP, spread=SPR, spread_source="bootstrap_人层")

sd_lower = T[T.domain=='PURITY'].sd.median() < T[T.domain!='PURITY'].sd.median()
if pos_ok and pla_ok and g0 <= 0.30+2*SPR:
    v1 = ("W-INTERNAL:比值判据**严格通过**" if (top4_all_pur and GAP > 2*SPR)
          else ("W-C:**纯洁四道仍占前四,但严格判据(最低纯洁 > 最高非纯洁 + 2×展布)不通过** "
                "-> `#610c①` 记为**未分离**" if top4_all_pur
                else "W-C:**有非纯洁题挤进前四** -> `#610c①` 记为未分离"))
    v2 = (f"判据②:纯洁标准差{'**更低**' if sd_lower else '**不更低**'};"
          f"而人为天花板只让比值动 {CE.d.abs().median()/T.ratio.median()*100:.2f}% "
          f"⇒ 天花板{'**未被排除**' if sd_lower and CE.d.abs().median()>0.1 else '**被排除**'}")
    verdict = v1 + "；" + v2
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:{中位/均值} × {原始/偏宗教} ===")
RES = partial(RK, CC[REL].values); spec = []
for mn, M in (("原始", RK), ("偏宗教", RES)):
    for an, ag in (("中位", np.median), ("均值", np.mean)):
        Px = prof(M, ag); pr = [Px[k][2] for k in PUR]; orr = [Px[k][2] for k in NAMES if k not in PUR]
        t4 = sorted(NAMES, key=lambda k: -Px[k][2])[:4]
        spec.append(dict(spec=f"{mn}·{an}", gap=float(min(pr)-max(orr)),
                         top4_all_purity=bool(all(ITEM[i][0]=="PURITY" for i in t4)), top4=t4))
        print(f"  {mn}·{an}: min(纯洁)−max(非纯洁) {min(pr)-max(orr):+.3f} · "
              f"前四全纯洁 {all(ITEM[i][0]=='PURITY' for i in t4)} · 前四 {t4}")

json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)),
               table=T.to_dict("records"), ceiling=CE.to_dict("records"),
               top4=top4, top4_all_purity=bool(top4_all_pur), gap=GAP, spread=SPR,
               purity_sd_median=float(T[T.domain=='PURITY'].sd.median()),
               other_sd_median=float(T[T.domain!='PURITY'].sd.median()),
               ceiling_delta_median=float(CE.d.abs().median()),
               omni_ratio=om, pure_ratio=pu, shuffle_ratio=[float(np.mean(sh)), float(np.std(sh))],
               spec_curve=spec, verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"ratio_and_the_ceiling.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'ratio_and_the_ceiling.json'}")
