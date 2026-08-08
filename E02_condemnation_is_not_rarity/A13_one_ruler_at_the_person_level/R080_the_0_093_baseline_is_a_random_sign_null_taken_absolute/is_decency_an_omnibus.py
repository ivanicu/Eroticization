"""E02·A13·R646 —— `decency` 的载重,是纯洁内部的,还是一道措辞最泛的总括题?

`#609` 的 NEXT。**BASIN RULE**:这是连续第三轮确认「纯洁特殊」,所以本轮**就是去攻击它**,
而 **W-A 的成立是我不希望的**(它把 `#609c` 那一行降级为「未分离」)。

§3 GRADIENT CHECK(在踩上去之前):一道总括题会**同时**抬高本域一致性**和**跨域相关 ——
所以「删掉它领先下降」这件事的方向**不是代数强制的**。两个领先口径都要量:
  `lead_max = within(PURITY) − max(其余四域 within)`(`#609` 用的)
  `lead_cross = within(PURITY) − cross`(跨域基线口径)
**只量一个,就是把一个未定的符号当成证据。**

INSTRUMENT(硬规则②):GHN 2009 JPSP Study 3(MFQ),
`data/external/dataverse/10.7910_DVN_SJTRBI_x/…Study_3.sav`。一份问卷 · 同一批人 · 同一量表。

G1 ESTIMAND(先于方法):`out(i)` = 题 i 与**不属于 i 所在域**的题两两 |ρ| 的中位。
  比较三个量(判据要求三个都报,不许只报有利的一支):
  ① `out(decency)` ② 其余五道纯洁题的 `out` ③ 其余四域 24 题的 `out`(基线)。
  辅助量:`ratio(i) = within(i)/out(i)` —— **总括题的 ratio 低**。

WORLDS:A 总括题(`out(decency)` 显著高于其余五道)· B 载重是纯洁内部的 · C 之间 -> 报区间

CONTROLS:
  正对照:种一道**真总括题** = 全部 30 题秩的均值 + 噪声。它的 `out` 必须远高于任何真题;
    **且 g=0(纯噪声)时必须不通过。**
  安慰剂:种一道**纯域内题** = 仅与纯洁六题相关。它的 `out` 必须落在低端。
  **offset_control**(「这个零该是零吗?」**不该**):每道题都有**一般共性**,`out` 从来不是 0。
    null 种类 = **「同一份问卷上每道题都带的一般共性」**,由 30 题 `out` 的分布给出。
KILL(条件式,预注册,写在跑之前):
  if 正对照(g=1)触发 and 正对照在 g=0 不通过 and 安慰剂落在低端:
      `out(decency)` > 其余五道纯洁题的最大值 + 2×展布 -> **W-A:总括题,`#609c` 降级为「未分离」**
      `out(decency)` 落在其余五道的区间内                -> **W-B:载重是纯洁内部的**
      否则                                              -> **W-C:报区间**
  else: UNVERIFIED
G3:30 题的 `out`/`within`/`ratio` 全表发布,含不一致的题。
G4:{原始 / 偏宗教} × {中位 / 均值} × 3 种子,并同时报两个领先口径。
IMPOSSIBLE(不写 planned):非因果 · **自选网络志愿者,非概率样本** ·
  `out` 只量「这道题有多泛」,**不量它为什么泛**(措辞?概念?)—— 要分开需要改写题干重测,
  而那需要一次新的数据采集 · 未派对抗 agent ⇒ `[unchallenged]`
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

d, _ = pyreadstat.read_sav(str(SAV))
CC = d[NAMES + [REL]].dropna()
print(f"仪器 = GHN 2009 Study 3(MFQ)· 30 题 + {REL} 完整 n={len(CC)}")
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
IDX = {c: i for i, c in enumerate(NAMES)}


def partial(M, v):
    z = rankdata(v).astype(float); z = (z - z.mean()) / z.std()
    Z = np.column_stack([np.ones(len(M)), z])
    return M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]


def profile(M, agg=np.median):
    """每题的 within(本域其余题)· out(域外全部题)· ratio。"""
    R = np.abs(np.corrcoef(M, rowvar=False))
    rows = []
    for it, (D, P) in ITEM.items():
        ins = [k for k in NAMES if ITEM[k][0] == D and k != it]
        outs = [k for k in NAMES if ITEM[k][0] != D]
        w = float(agg([R[IDX[it], IDX[k]] for k in ins]))
        o = float(agg([R[IDX[it], IDX[k]] for k in outs]))
        rows.append(dict(item=it, domain=D, part=P, within=w, out=o, ratio=w/o if o > 1e-9 else np.nan))
    return pd.DataFrame(rows)


def leads(M, subset=None):
    R = np.abs(np.corrcoef(M, rowvar=False))
    w = {}
    for D in DOMS:
        its = [k for k, v in ITEM.items() if v[0] == D]
        if subset is not None and D == "PURITY": its = subset
        w[D] = float(np.median([R[IDX[a], IDX[b]] for a, b in itertools.combinations(its, 2)]))
    xs = []
    for a, b in itertools.combinations(NAMES, 2):
        if ITEM[a][0] == ITEM[b][0]: continue
        if subset is not None and (a in PUR and a not in subset or b in PUR and b not in subset): continue
        xs.append(R[IDX[a], IDX[b]])
    cross = float(np.median(xs))
    return (w["PURITY"] - max(w[D] for D in DOMS if D != "PURITY"), w["PURITY"] - cross, w, cross)


T = profile(RK)
print("\n=== G3 · 30 题全表(按 out 从高到低,总括题应排在最上)===")
S = T.sort_values("out", ascending=False)
for r in S.itertuples():
    mark = " ←" if r.item == "decency" else ("  *" if r.domain == "PURITY" else "")
    print(f"  {r.item:13s} {r.domain:9s} within {r.within:.4f}  out {r.out:.4f}  ratio {r.ratio:5.2f}{mark}")
dec = T[T.item == "decency"].iloc[0]
pur5 = T[(T.domain == "PURITY") & (T.item != "decency")]
oth = T[T.domain != "PURITY"]
print(f"\n**判据要求的三个数(都报)**")
print(f"  ① `decency` 的域外相关 out = **{dec.out:.4f}**   (ratio {dec.ratio:.2f})")
print(f"  ② 其余五道纯洁题 out 中位 **{pur5.out.median():.4f}** 区间 [{pur5.out.min():.4f}, {pur5.out.max():.4f}]")
print(f"  ③ 其余四域 24 题 out 中位 **{oth.out.median():.4f}** 区间 [{oth.out.min():.4f}, {oth.out.max():.4f}](基线)")
print(f"  `decency` 在 30 题里的 out 排名 = **第 {int((T.out > dec.out).sum())+1} 位**")

# 展布
bs = []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(150):
        Tb = profile(RK[rng.integers(0, len(RK), len(RK))])
        db = Tb[Tb.item == "decency"].out.iloc[0]
        p5 = Tb[(Tb.domain == "PURITY") & (Tb.item != "decency")].out
        bs.append(float(db - p5.max()))
SPR = float(np.std(bs))
GAP = float(dec.out - pur5.out.max())
print(f"\n  **out(decency) − max(其余五道) = {GAP:+.4f}** ± {SPR:.4f}(bootstrap 人层 × {len(SEEDS)} 种子)")

# ── §3 梯度检查兑现:两个领先口径都量 ─────────────────────────────
print("\n=== §3 GRADIENT CHECK 兑现:两个领先口径,原始 / 偏宗教 ===")
RES = partial(RK, CC[REL].values)
print(f"  {'删掉':12s} {'lead_max原始':>12s} {'lead_cross原始':>14s} {'lead_max偏后':>12s} {'lead_cross偏后':>14s}")
LL = []
for drop in [None] + PUR:
    sub = [k for k in PUR if k != drop] if drop else PUR
    a0, b0, _, _ = leads(RK, sub); a1, b1, _, _ = leads(RES, sub)
    LL.append(dict(dropped=drop or "—", lead_max_raw=a0, lead_cross_raw=b0, lead_max_par=a1, lead_cross_par=b1))
    print(f"  {(drop or '—(全六道)'):12s} {a0:12.4f} {b0:14.4f} {a1:12.4f} {b1:14.4f}")
L = pd.DataFrame(LL)
base = L[L.dropped == "—"].iloc[0]; dl = L[L.dropped == "decency"].iloc[0]
print(f"\n  删 `decency`:`lead_max` {base.lead_max_raw:+.4f}->{dl.lead_max_raw:+.4f} · "
      f"`lead_cross` {base.lead_cross_raw:+.4f}->{dl.lead_cross_raw:+.4f}")
print(f"  ⇒ 两个口径{'**同向**' if np.sign(dl.lead_max_raw-base.lead_max_raw)==np.sign(dl.lead_cross_raw-base.lead_cross_raw) else '**反向 —— 那个符号是口径造的,不是数据**'}")

# ── 控制 ────────────────────────────────────────────────────
G = Gate("`decency` 的载重,是纯洁内部的,还是一道措辞最泛的总括题?")


def synth_out(kind, g, sd):
    rng = np.random.default_rng(sd)
    Z = (RK - RK.mean(0)) / RK.std(0)
    core = Z.mean(1) if kind == "omni" else Z[:, [IDX[k] for k in PUR]].mean(1)
    core = (core - core.mean()) / core.std()
    v = g*core + np.sqrt(max(0., 1-g**2))*rng.standard_normal(len(Z))
    M = np.column_stack([RK, rankdata(v)])
    R = np.abs(np.corrcoef(M, rowvar=False))
    j = M.shape[1]-1
    outs = [k for k in NAMES if ITEM[k][0] != "PURITY"]     # 视合成题为纯洁域成员
    return float(np.median([R[j, IDX[k]] for k in outs]))


print("\n  正对照:种一道**真总括题**(全 30 题秩的均值 + 噪声),它的 out 必须远高于任何真题")
omni = {g: float(np.mean([synth_out("omni", g, sd) for sd in SEEDS])) for g in (0.0, 0.5, 1.0)}
for g, v in omni.items(): print(f"    g={g:.1f} -> out {v:.4f}   (真题最大 out = {T.out.max():.4f})")
pos_ok = G.positive_control("正对照:合成总括题(g=1)的 out 必须超过任何真题",
                            planted=float(omni[1.0] - T.out.max()), floor=0.02, spread=SPR)
g0_pass = (omni[0.0] - T.out.max()) > 0.02 + 2*SPR
print(f"    g=0 时 out {omni[0.0]:.4f} ⇒ 判据在 g=0 {'⛔ 也通过,不合格' if g0_pass else '**不通过** ✅'}")
pla = {g: float(np.mean([synth_out("pure", g, sd) for sd in SEEDS])) for g in (0.5, 1.0)}
print(f"\n  安慰剂:种一道**纯域内题**(只与纯洁六题相关)-> out g=1.0 {pla[1.0]:.4f} · g=0.5 {pla[0.5]:.4f}"
      f"  (须落在真题低端 {T.out.min():.4f})")
pla_ok = G.negative_control("安慰剂:纯域内合成题的 out 必须落在低端",
                            null=float(pla[1.0] - T.out.min()), effect=float(omni[1.0] - T.out.min()),
                            null_spread=SPR, null_kind="只与本域相关、与域外无结构的合成题")
G.offset_control("offset:同一份问卷上每道题都带的一般共性", effect=float(dec.out),
                 offset=float(T.out.median()), spread=SPR,
                 null_kind="每道题的一般共性;out 从来不是 0,基线由 30 题 out 的分布给出")
G.has_error_bar("out(decency) − max(其余五道纯洁题)", value=GAP, spread=SPR, spread_source="bootstrap_人层")

if pos_ok and pla_ok and not g0_pass:
    if GAP > 2*SPR: verdict = "W-A:**`decency` 是总括题 —— `#609c` 那行降级为「未分离」**"
    elif dec.out <= pur5.out.max() and dec.out >= pur5.out.min(): verdict = "W-B:**载重是纯洁内部的**"
    else: verdict = "W-C:介于其间 —— 报区间,不下二选一"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok} · g=0 不通过 {not g0_pass})"
    print(f"\n⚠ {verdict}")
print(G)

# G4 规格曲线
print("\n=== G4 规格曲线:{原始/偏宗教} × {中位/均值} ===")
spec = []
for mn, M in (("原始", RK), ("偏宗教", RES)):
    for an, ag in (("中位", np.median), ("均值", np.mean)):
        Tx = profile(M, ag); dx = Tx[Tx.item == "decency"].out.iloc[0]
        p5 = Tx[(Tx.domain == "PURITY") & (Tx.item != "decency")].out
        spec.append(dict(spec=f"{mn}·{an}", out_decency=float(dx), out_p5_max=float(p5.max()),
                         gap=float(dx - p5.max()), rank=int((Tx.out > dx).sum())+1))
        print(f"  {mn}·{an}: out(decency) {dx:.4f} · 其余五道最大 {p5.max():.4f} · "
              f"差 {dx-p5.max():+.4f} · 排名第 {int((Tx.out>dx).sum())+1}")
json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)),
               per_item=T.to_dict("records"), leads=L.to_dict("records"),
               out_decency=float(dec.out), out_purity5_median=float(pur5.out.median()),
               out_purity5_range=[float(pur5.out.min()), float(pur5.out.max())],
               out_other24_median=float(oth.out.median()),
               out_other24_range=[float(oth.out.min()), float(oth.out.max())],
               rank_of_decency=int((T.out > dec.out).sum())+1, gap=GAP, spread=SPR,
               omni_control=omni, placebo=pla, spec_curve=spec, verdict=verdict,
               seeds=SEEDS, unchallenged=True),
          open(OUT/"is_decency_an_omnibus.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'is_decency_an_omnibus.json'}")
