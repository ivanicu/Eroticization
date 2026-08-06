"""E02·A13·R645 —— 纯洁领先里「不是宗教」的那一半,是纯洁本身,还是一道提上帝的题在撑?

`#608` 的 NEXT,**并加宽**:原判据只说剔 `god` 一道 ——
**那正是「剔掉对我有利的那一道」的做法**。这里 **六道每一道都剔一次**,整条曲线发布(G4)。

INSTRUMENT(硬规则②):**Graham–Haidt–Nosek 2009 JPSP Study 3(MFQ)**,
`data/external/dataverse/10.7910_DVN_SJTRBI_x/…Study_3.sav`。一份问卷 · 同一批人 · 同一量表。
题目文本取自仪器发布方(moralfoundations.org self-scorable MFQ),不是变量名字面,也不是记忆。

G1 ESTIMAND(先于方法,三个,都在选统计量之前命名):
  E1 `retain(i)` = 题 i 与本域其余题的 |ρ| 中位,偏掉礼拜出席后 / 偏掉之前。**逐题**,五域全做。
  E2 `lead(S)` = `within(PURITY|S) − max(其余四域 within)`,S 遍历纯洁的 **6 个留一子集** × {原始, 偏宗教}。
  E3 纯洁在**每一个**留一变体下是否仍是五域第一(布尔)。

WORLDS:
  A `god` 的保留率显著低于其余五道 -> 领先由一道**宗教题**携带,剔掉它领先应大幅塌
  B 六道保留率相当                -> 那一半是**纯洁本身**,与题目提不提上帝无关
  C 剔掉任一道纯洁就掉出第一      -> **领先是单题伪影**
**A 与 C 都是我不希望的**(两者都削掉 `#608` 的一半)。

CONTROLS:
  安慰剂:改偏**年龄** -> 所有保留率必须 ≈ 100%(`#608d` 已测得纯洁领先保留 103.0%)。
  **offset_control**(「这个零该是零吗?」**不该**):偏掉任何与两端都相关的变量都有**机械衰减**。
    null 种类 = **「其余四个道德域在同一次偏相关下的保留率」** —— 它们同样被偏,同样衰减,
    但**没有一道题提上帝**,所以它们量的正是「与宗教内容无关的那部分衰减」。
  正对照:在纯洁域里**种一道合成的纯宗教题**(= 礼拜出席的秩 + 噪声),
    它的保留率必须**显著低于**任何真实题;且 **g=0(零权重)时必须不通过**。
KILL(条件式,预注册,写在跑之前):
  if 安慰剂(年龄)保留 ∈ [0.95,1.05] and 正对照种植题保留 < 真题最低 − 2×展布:
      `god` 保留 < 其余五道最低 − 2×展布            -> **W-A**
      六道保留全距 < 2×展布                          -> **W-B**
      任一留一变体下纯洁掉出第一                     -> **W-C**(与上面可同时成立,分别报)
  else: UNVERIFIED
G3:5 域 × 30 题的**逐题保留率全表**发布,含不一致的题。
G4:6 个留一子集 × {原始/偏宗教} × 3 种子 的整条曲线。
IMPOSSIBLE(不写 planned):非因果 · 无法定序 · **自选网络志愿者,非概率样本** ·
  六道纯洁题只有一道点到性行为 ⇒ 这仍是关于**纯洁这一族**的话 · 未派对抗 agent ⇒ `[unchallenged]`
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
DOMS = sorted(set(v[0] for v in ITEM.values()))
PUR = [k for k, v in ITEM.items() if v[0] == "PURITY"]

d, meta = pyreadstat.read_sav(str(SAV))
need = list(ITEM) + [REL, "Age"]
CC = d[need].dropna(subset=list(ITEM) + [REL])
print(f"仪器 = GHN 2009 Study 3(MFQ)· 全表 n={len(d)} · 30 题 + {REL} 完整 n={len(CC)}")
print(f"{REL} 取值 {sorted(CC[REL].unique())}  纯洁六题 = {PUR}")

NAMES = list(ITEM)
RK = np.column_stack([rankdata(CC[c].values) for c in NAMES]).astype(float)
IDX = {c: i for i, c in enumerate(NAMES)}


def partial(M, v):
    z = rankdata(v).astype(float); z = (z - z.mean()) / z.std()
    Z = np.column_stack([np.ones(len(M)), z])
    return M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]


def item_coh(M, item, peers):
    R = np.abs(np.corrcoef(M, rowvar=False))
    return float(np.median([R[IDX[item], IDX[p]] for p in peers if p != item]))


def within(M, subset=None):
    R = np.abs(np.corrcoef(M, rowvar=False))
    out = {}
    for D in DOMS:
        its = [k for k, v in ITEM.items() if v[0] == D]
        if subset is not None and D == "PURITY": its = subset
        out[D] = float(np.median([R[IDX[a], IDX[b]] for a, b in itertools.combinations(its, 2)]))
    return out


RES_REL = partial(RK, CC[REL].values)
RES_AGE = partial(RK, CC["Age"].values)

# ── E1:逐题保留率,五域全做(G3 整表)────────────────────────────
print("\n=== E1 · G3:逐题保留率(偏掉礼拜出席)—— 30 题全表 ===")
rows = []
for it, (D, P) in ITEM.items():
    peers = [k for k, v in ITEM.items() if v[0] == D]
    a = item_coh(RK, it, peers); b = item_coh(RES_REL, it, peers); c = item_coh(RES_AGE, it, peers)
    rows.append(dict(item=it, domain=D, part=P, raw=a, par_rel=b, retain=b/a, retain_age=c/a))
T = pd.DataFrame(rows)
for D in DOMS:
    s = T[T.domain == D].sort_values("retain")
    print(f"  {D:10s} 保留率 中位 {s.retain.median()*100:5.1f}%  " +
          " · ".join(f"{r.item}={r.retain*100:.1f}%" for r in s.itertuples()))
G_ret = T[T.domain != "PURITY"].retain
P_ret = T[T.domain == "PURITY"].retain
print(f"\n  **纯洁六题保留率 中位 {P_ret.median()*100:.1f}% · 全距 {(P_ret.max()-P_ret.min())*100:.1f} 个百分点**")
print(f"  **其余四域 24 题 中位 {G_ret.median()*100:.1f}%**(offset:与宗教内容无关的那部分衰减)")
god_r = float(T.loc[T.item == "god", "retain"].iloc[0])
oth_p = T[(T.domain == "PURITY") & (T.item != "god")].retain
print(f"  **`god` 保留 {god_r*100:.1f}% · 其余五道最低 {oth_p.min()*100:.1f}% 中位 {oth_p.median()*100:.1f}%**")

# 展布:bootstrap
bs_gap, bs_lead, bs_rng = [], [], []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(200):
        i = rng.integers(0, len(RK), len(RK))
        M = RK[i]; Mr = partial(M, CC[REL].values[i])
        rr = {}
        for it in PUR:
            a = item_coh(M, it, PUR); rr[it] = item_coh(Mr, it, PUR)/a
        v = np.array(list(rr.values()))
        bs_rng.append(float(v.max()-v.min()))
        bs_gap.append(float(np.min([rr[k] for k in PUR if k != "god"]) - rr["god"]))
        w = within(Mr); bs_lead.append(float(w["PURITY"] - max(w[D] for D in DOMS if D != "PURITY")))
SPR = float(np.std([x for x in bs_rng]))
print(f"  展布(bootstrap 人层 × {len(SEEDS)} 种子,{len(bs_rng)} 次):保留率全距 sd = {SPR:.4f}")

# ── E2/E3 · G4:六个留一子集 ─────────────────────────────────
print("\n=== E2/E3 · G4:纯洁六题**每一道都剔一次**,原始 / 偏宗教两条 ===")
loo = []
for drop in [None] + PUR:
    sub = [k for k in PUR if k != drop] if drop else PUR
    w0, w1 = within(RK, sub), within(RES_REL, sub)
    o0 = max(w0[D] for D in DOMS if D != "PURITY"); o1 = max(w1[D] for D in DOMS if D != "PURITY")
    rank1_0 = w0["PURITY"] == max(w0.values()); rank1_1 = w1["PURITY"] == max(w1.values())
    loo.append(dict(dropped=drop or "—(全六道)", k=len(sub), pur_raw=w0["PURITY"], lead_raw=w0["PURITY"]-o0,
                    pur_par=w1["PURITY"], lead_par=w1["PURITY"]-o1, rank1_raw=bool(rank1_0), rank1_par=bool(rank1_1)))
L = pd.DataFrame(loo)
print(L.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print(f"\n  **剔掉 `god`:领先 {L.loc[L.dropped=='god','lead_raw'].iloc[0]:+.4f}(原始)· "
      f"{L.loc[L.dropped=='god','lead_par'].iloc[0]:+.4f}(偏宗教后)**")
print(f"  **全六道:{L.loc[L.dropped=='—(全六道)','lead_raw'].iloc[0]:+.4f} · "
      f"{L.loc[L.dropped=='—(全六道)','lead_par'].iloc[0]:+.4f}**")
print(f"  **纯洁在 {int(L.rank1_raw.sum())}/{len(L)} 个变体下仍是原始第一 · "
      f"{int(L.rank1_par.sum())}/{len(L)} 个变体下偏宗教后仍是第一**")

# ── 控制 ────────────────────────────────────────────────────
G = Gate("纯洁领先里「不是宗教」的那一半,是纯洁本身,还是一道提上帝的题在撑?")
age_ok = G.negative_control("安慰剂:改偏年龄(所有题保留率须≈100%)",
                            null=float(abs(1-T.retain_age.median())), effect=float(abs(1-T.retain.median())),
                            null_spread=float(T.retain_age.std()), null_kind="与道德域结构无关的人口学变量")
print(f"\n  安慰剂:改偏年龄 -> 30 题保留率中位 {T.retain_age.median()*100:.1f}%"
      f"(纯洁 {T[T.domain=='PURITY'].retain_age.median()*100:.1f}%)")

# 正对照:在纯洁域里种一道合成的纯宗教题
print("\n  正对照:种一道**合成纯宗教题**进纯洁域,它的保留率必须显著低于任何真题")
plant_ret = {}
for gval in [0.0, 0.5, 1.0]:
    r2 = np.random.default_rng(SEEDS[0])
    z = rankdata(CC[REL].values).astype(float); z = (z-z.mean())/z.std()
    synth = gval*z + np.sqrt(max(0.,1-gval**2))*r2.standard_normal(len(CC))
    M2 = np.column_stack([RK, rankdata(synth)])
    IDX["_plant"] = M2.shape[1]-1
    a = item_coh(M2, "_plant", PUR + ["_plant"])
    M2r = partial(M2, CC[REL].values)
    b = item_coh(M2r, "_plant", PUR + ["_plant"])
    plant_ret[gval] = b/a if a > 1e-9 else np.nan
    print(f"    g={gval:.1f} -> 种植题保留 {plant_ret[gval]*100:6.1f}%  (原始一致 {a:.4f})")
real_min = float(T[T.domain == "PURITY"].retain.min())
pos_ok = G.positive_control("正对照:合成纯宗教题(g=1)的保留率必须远低于任何真题",
                            planted=float(real_min - plant_ret[1.0]), floor=0.05, spread=SPR)
print(f"    g=0 时保留 {plant_ret[0.0]*100:.1f}% ⇒ 判据在 g=0 **不通过** "
      f"({'✅' if (real_min - plant_ret[0.0]) < 0.05 + 2*SPR else '⛔ 正对照在 g=0 也通过,不合格'})")
G.offset_control("offset:其余四域在同一次偏相关下的保留率",
                 effect=float(1-P_ret.median()), offset=float(1-G_ret.median()), spread=SPR,
                 null_kind="偏掉任何与两端都相关的变量都有的机械衰减;这四域没有一道题提上帝")
G.has_error_bar("纯洁领先(偏宗教后,全六道)", value=float(L.loc[L.dropped=='—(全六道)','lead_par'].iloc[0]),
                spread=float(np.std(bs_lead)), spread_source="bootstrap_人层")

if age_ok and pos_ok:
    wA = god_r < oth_p.min() - 2*SPR
    wB = (P_ret.max()-P_ret.min()) < 2*SPR
    wC = not bool(L.rank1_par.all())
    verdict = (("W-A:领先由 `god` 那一道宗教题携带" if wA else "") +
               ("；" if wA and wB else "") + ("W-B:六道保留率相当 -> 那一半是**纯洁本身**" if wB else "") +
               ("；W-C:某个留一变体下纯洁掉出第一 -> 领先含单题成分" if wC else
                "；纯洁在**每一个**留一变体下都仍是第一 -> **领先不是任何单题的伪影**"))
    if not (wA or wB): verdict = "两支都不成立:`god` 低但未超 2×展布,且六道保留率全距也超了 —— 报区间,不下二选一" + \
        ("；纯洁在每个留一变体下都仍是第一" if not wC else "；某变体下掉出第一")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(安慰剂 {age_ok} · 正对照 {pos_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(instrument="GHN 2009 JPSP Study 3 (MFQ)", n=int(len(CC)),
               per_item=T.to_dict("records"), leave_one_out=L.to_dict("records"),
               god_retain=god_r, purity_other_min=float(oth_p.min()), purity_other_median=float(oth_p.median()),
               purity_retain_median=float(P_ret.median()), purity_retain_range=float(P_ret.max()-P_ret.min()),
               other_domains_retain_median=float(G_ret.median()), spread=SPR,
               lead_spread=float(np.std(bs_lead)), plant_retain=plant_ret,
               placebo_age_retain_median=float(T.retain_age.median()),
               rank1_raw=int(L.rank1_raw.sum()), rank1_par=int(L.rank1_par.sum()), variants=len(L),
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT/"which_item_carries_it.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'which_item_carries_it.json'}")
