"""E02·A13·R644 —— 该切的是「性 vs 非性」,还是「同一个道德域 vs 跨域」?

`#607` 的 NEXT 换仪器。**判据在跑之前就已经判掉了它自己的第一问**:
MFQ 只有 **1 道**明确的性道德题(`chastity`),不是 ≥3 ⇒ 「那把尺子是不是宗教」**在本地数据上问不了**。
而同一次读码打开了一个更该问的问题,**它的阳性结果同样是我不希望的**。

INSTRUMENT(硬规则②,先命名):**Graham–Haidt–Nosek 2009 JPSP Study 3**,
`data/external/dataverse/10.7910_DVN_SJTRBI_x/…Study_3.sav`,n=8,193(逐题 n≈7,400)。
**一份问卷、同一批人、无分票** —— 这正是 GSS 给不了的(`#607b`:30.7% 组合零重叠)。
题目文本取自仪器发布方(moralfoundations.org 的 self-scorable MFQ),**不是我的记忆**。

G1 ESTIMAND(先于方法):对每个道德域 D,
  `within(D) = median|ρ|(D 内题两两)` · `cross = median|ρ|(不同域的题两两)`。
  关心的量:**`within` 在五个域之间是否可分辨**,以及 **`within − cross`**。
  ⚠ 这不是「性题 vs 参照题」的重做 —— 换的是**分类轴本身**,这是 meta-separator。

WORLDS:
  A 性/纯洁特殊  -> `within(PURITY)` 明显高于另外四个域
  B **域才是刀口** -> 五个域的 `within` 不可分辨,而 `cross` 明显更低  ⇒ **我的世界分解是错的**
  C 格式主导     -> 域内/跨域差别小,而**同 Part / 跨 Part** 的差别大

⚠ MFQ 自带一个 GSS 结构上没有的东西:**Part 1 是「相关性」量表,Part 2 是「同意度」量表**
  ⇒ **格式的贡献可以直接量**,而 `#535` 记的正是「GSS 里不存在这样一组对照」。

CONTROLS:
  正对照:**同域且同 Part** 的题对必须是最高的一档(且 g=0 —— 随机指派域时**必须不通过**)。
  安慰剂:`astrology`(仪器自己的 catch 题)与所有真题的相关,必须≈0。
  **offset_control**(「这个零该是零吗?」**不该**):任何两道同一份问卷上的题都共享
    **默认同意倾向 + 一般道德化程度**。null 种类 = **「同一仪器上任意两题的共同响应风格」**,
    由 catch 题 `astrology` 的相关给出下界,由**跨域跨 Part** 给出上界。
KILL(条件式,预注册):控制齐备才评判 ——
  `within(PURITY) − max(其余四域 within) > 2×展布` -> W-A
  五域 `within` 全距 < 2×展布 且 `median(within) − cross > 2×展布` -> **W-B(我不希望的那个)**
  否则 -> W-C / 未定
G3:5 域 × {同 Part / 跨 Part} × {域内 / 跨域} 整张网格发布,含不一致的格。
G4:规格曲线 = 是否剔除答题不认真者(`MFQ_failed`)× 是否用秩 × 三个种子。
IMPOSSIBLE(不写 planned):
  · 单一时点、无干预 ⇒ **非因果** · 自选样本(yourmorals.org 网络志愿者)⇒ **非概率抽样**,
    分布外推做不到,需要一份带抽样框的 MFQ · 只有 1 道性题 ⇒ **「性域」本身在这具仪器上量不了**,
    需要一份含 ≥3 道同格式性行为道德题的问卷 · 未派对抗 agent ⇒ `[unchallenged]`
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

# 题 -> (域, Part)。**全部来自仪器发布方的题目文本**,不是变量名的字面猜测。
ITEM = {
 "emotionally":("HARM",1), "weak":("HARM",1), "cruel":("HARM",1),
 "compassion":("HARM",2), "animal":("HARM",2), "kill":("HARM",2),
 "treated":("FAIRNESS",1), "unfairly":("FAIRNESS",1), "rights":("FAIRNESS",1),
 "fairly":("FAIRNESS",2), "justice":("FAIRNESS",2), "rich":("FAIRNESS",2),
 "lovecountry":("INGROUP",1), "betray":("INGROUP",1), "loyalty":("INGROUP",1),
 "history":("INGROUP",2), "family":("INGROUP",2), "team":("INGROUP",2),
 "respect":("AUTHORITY",1), "traditions":("AUTHORITY",1), "chaos":("AUTHORITY",1),
 "kidrespect":("AUTHORITY",2), "sexroles":("AUTHORITY",2), "soldier":("AUTHORITY",2),
 "decency":("PURITY",1), "disgusting":("PURITY",1), "god":("PURITY",1),
 "harmlessdg":("PURITY",2), "unnatural":("PURITY",2), "chastity":("PURITY",2),
}
CATCH = "astrology"
# ⚠ Study 3 还有 8 道 MFQ41 的旧题(harm·dominate·yourgroup·duties·desires·victim·shutup·temples)——
#   **我手上没有它们的官方文本,所以全部排除**,而不是按变量名猜域。这是硬规则①。
EXCLUDED = ["harm","dominate","yourgroup","duties","desires","victim","shutup","temples"]

d, meta = pyreadstat.read_sav(str(SAV))
print(f"仪器 = GHN 2009 JPSP Study 3 · n={len(d)} · 变量 {len(d.columns)}")
print(f"纳入 {len(ITEM)} 道(五域各 6 道,同 Part 各 3 道)· 排除 {len(EXCLUDED)} 道无官方文本的旧题")
print("=== 硬规则①:逐题打印 n 与取值 ===")
for c in list(ITEM) + [CATCH]:
    s = d[c].dropna(); assert set(s.unique()) <= {0.,1.,2.,3.,4.,5.}, f"{c} 取值超范围"
    print(f"  {c:13s} {ITEM.get(c,('CATCH',0))[0]:9s} P{ITEM.get(c,('',0))[1]} n={len(s):5d} 均值 {s.mean():.2f}")


# ⚠ 第一版超时:每次 bootstrap 都重算 435 个两两秩相关。
#   改法:**一次建秩矩阵,bootstrap 只重算相关阵** —— 同一个估计量,只是不再重复排序。
ITEMS = list(ITEM) + [CATCH]
CC = d[ITEMS].dropna()
print(f"\n完整作答(30 题 + catch 全部非缺失)n = {len(CC)} / {len(d)}  "
      f"⇒ 用 listwise,**所有相关都在同一批人上**(GSS 结构上做不到的那一点)")
RK = np.column_stack([rankdata(CC[c].values) for c in ITEMS]).astype(float)
IDX = {c: i for i, c in enumerate(ITEMS)}
DOMS = sorted(set(v[0] for v in ITEM.values()))

def stats_from(M, imap=ITEM):
    R = np.abs(np.corrcoef(M, rowvar=False))
    W = {D: [] for D in sorted(set(v[0] for v in imap.values()))}
    Ws, Wc, X, Xs, Xc = {D: [] for D in W}, {D: [] for D in W}, [], [], []
    for a, b in itertools.combinations(imap, 2):
        r = R[IDX[a], IDX[b]]
        (Da, Pa), (Db, Pb) = imap[a], imap[b]
        sp = (Pa == Pb)
        if Da == Db: W[Da].append(r); (Ws if sp else Wc)[Da].append(r)
        else: X.append(r); (Xs if sp else Xc).append(r)
    return W, Ws, Wc, X, Xs, Xc

def catch_rhos(M):
    R = np.abs(np.corrcoef(M, rowvar=False))
    return [R[IDX[CATCH], IDX[c]] for c in ITEM]

W, Wsame, Wcross, X, Xsame, Xcross = stats_from(RK)
print(f"\n=== G3 整张网格(listwise n={len(CC)})===")
print(f"{'域':10s} {'域内对':>5s} {'within':>8s} {'同Part':>8s} {'跨Part':>8s}")
rows = []
for D in DOMS:
    a, b, c = np.median(W[D]), np.median(Wsame[D]), np.median(Wcross[D])
    rows.append(dict(domain=D, k=len(W[D]), within=float(a), same_part=float(b), cross_part=float(c)))
    print(f"{D:10s} {len(W[D]):5d} {a:8.4f} {b:8.4f} {c:8.4f}")
CROSS = float(np.median(X))
print(f"{'跨域':10s} {len(X):5d} {CROSS:8.4f} {np.median(Xsame):8.4f} {np.median(Xcross):8.4f}")

wv = np.array([r["within"] for r in rows]); names = [r["domain"] for r in rows]
PUR = float(wv[names.index("PURITY")]); OTH = float(np.max(np.delete(wv, names.index("PURITY"))))
bs_w, bs_x, bs_rng, bs_pd = [], [], [], []
for sd in SEEDS:
    rng = np.random.default_rng(sd)
    for _ in range(200):
        M = RK[rng.integers(0, len(RK), len(RK))]
        Wb, _, _, Xb, _, _ = stats_from(M)
        wb = np.array([np.median(Wb[D]) for D in DOMS])
        bs_w.append(float(np.median(wb))); bs_x.append(float(np.median(Xb)))
        bs_rng.append(float(wb.max()-wb.min()))
        bs_pd.append(float(wb[DOMS.index("PURITY")] - np.max(np.delete(wb, DOMS.index("PURITY")))))
SPREAD = float(np.std(bs_w))
print(f"\n展布(bootstrap 人层 × {len(SEEDS)} 种子,{len(bs_w)} 次):within 中位 sd = {SPREAD:.4f}")
print(f"**五域 within 全距 = {wv.max()-wv.min():.4f}** · 其 bootstrap 95%CI "
      f"[{np.quantile(bs_rng,.025):.4f}, {np.quantile(bs_rng,.975):.4f}]")
print(f"**PURITY {PUR:.4f} · 其余最高 {OTH:.4f} · 差 {PUR-OTH:+.4f}** "
      f"95%CI [{np.quantile(bs_pd,.025):+.4f}, {np.quantile(bs_pd,.975):+.4f}]")
print(f"**median(within) {np.median(wv):.4f} − cross {CROSS:.4f} = {np.median(wv)-CROSS:+.4f}**")

fs = float(np.median([r["same_part"] for r in rows])); fc = float(np.median([r["cross_part"] for r in rows]))
print(f"\n=== 格式的贡献(`#535` 说 GSS 里不存在这组对照)===")
print(f"  域内:同 Part {fs:.4f} vs 跨 Part {fc:.4f} -> **换一次答题量表拿走 {(fs-fc)/fs*100:.1f}%**")
print(f"  跨域:同 Part {np.median(Xsame):.4f} vs 跨 Part {np.median(Xcross):.4f}")

G = Gate("该切的是「性 vs 非性」,还是「同一个道德域 vs 跨域」?(MFQ,GHN 2009 Study 3)")
shuf = []
keys = list(ITEM); doms0 = [ITEM[k][0] for k in keys]
for sd in SEEDS:
    r2 = np.random.default_rng(sd)
    for _ in range(200):
        perm = r2.permutation(doms0)
        FAKE = {k: (perm[i], ITEM[k][1]) for i, k in enumerate(keys)}
        Wf, _, _, Xf, _, _ = stats_from(RK, FAKE)
        shuf.append(float(np.median([np.median(Wf[D]) for D in sorted(Wf)]) - np.median(Xf)))
SHUF = float(np.mean(shuf))
pos_ok = G.positive_control("正对照:真域指派的 within−cross 必须超过随机指派",
                            planted=float(np.median(wv)-CROSS), floor=SHUF, spread=float(np.std(shuf)))
print(f"\n  g=0 检验:**随机打乱域指派 -> within−cross = {SHUF:+.4f} ± {np.std(shuf):.4f}**"
      f"({len(shuf)} 次)⇒ 判据在 g=0 **不通过** ✅")
cat = catch_rhos(RK)
pla_ok = G.negative_control(f"安慰剂:仪器自带的 catch 题 `{CATCH}` × 全部真题",
                            null=float(np.median(cat)), effect=float(np.median(wv)-CROSS),
                            null_spread=float(np.std(cat)), null_kind="仪器自带的无意义题(响应风格的下界)")
G.offset_control("offset:同一仪器上任意两题的共同响应风格", effect=float(np.median(wv)),
                 offset=CROSS, spread=SPREAD,
                 null_kind="默认同意倾向 + 一般道德化程度;下界由 catch 题给出,上界由跨域跨 Part 给出")
G.has_error_bar("within 中位", value=float(np.median(wv)), spread=SPREAD, spread_source="bootstrap_人层")
print(f"  catch 题 `{CATCH}` 与真题的 |ρ| 中位 = {np.median(cat):.4f}(下界)· 跨域跨 Part = {np.median(Xcross):.4f}(上界)")

if pos_ok and pla_ok:
    if PUR - OTH > 2*SPREAD: verdict = "W-A:**纯洁/性域确实特殊**"
    elif (wv.max()-wv.min()) < 2*SPREAD and (np.median(wv)-CROSS) > 2*SPREAD:
        verdict = "W-B:**刀切错了 —— 该切的是「域内 vs 跨域」,不是「性 vs 非性」**"
    else: verdict = "W-C / 未定:五域可分辨,但纯洁不是独高 —— 分解需要重写而不是翻转"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

# G4 规格曲线
print("\n=== G4 规格曲线:剔不剔除答题不认真者 ===")
spec = []
ok = d.loc[CC.index, "MFQ_failed"].fillna(1).values == 0
for nm, sel in [("全样本", np.ones(len(CC), bool)), ("剔除 MFQ_failed", ok)]:
    if sel.sum() < 500: print(f"  {nm}: n={sel.sum()} 太小,跳过"); continue
    Wq, _, _, Xq, _, _ = stats_from(RK[sel])
    v = np.array([np.median(Wq[D]) for D in DOMS])
    spec.append(dict(spec=nm, n=int(sel.sum()), within=float(np.median(v)), cross=float(np.median(Xq)),
                     purity=float(v[DOMS.index("PURITY")]), rng=float(v.max()-v.min())))
    print(f"  {nm:16s} n={int(sel.sum()):5d} within {np.median(v):.4f} · cross {np.median(Xq):.4f} · "
          f"PURITY {v[DOMS.index('PURITY')]:.4f} · 五域全距 {v.max()-v.min():.4f}")


# ── 最强混淆,控制在同一轮内(realstat G2)────────────────────────────
# PURITY 的六道题里有一道是 "acted in a way that God would approve of" ->
# **纯洁域的高一致,会不会就是宗教?** 这正是 `#607` 在 GSS 上答不了的那个问题:
# 那里因为分票 + 子样本随偏相关移动而失败(安慰剂也抹掉 80%)。
# **这里没有分票,所有题同一批人、同一格式、listwise —— 那两个病都不存在。**
print("\n=== 最强混淆:纯洁域的高一致是不是宗教?(`#607` 在 GSS 上答不了的那一问)===")
REL = "Religion_attend_num"
rr = d.loc[CC.index, REL]
mrel = rr.notna().values
print(f"  {REL} n={int(mrel.sum())}/{len(CC)} 取值 {sorted(rr.dropna().unique())}")
RK2 = RK[mrel]
relr = rankdata(rr[mrel].values).astype(float); relr = (relr - relr.mean())/relr.std()
Z = np.column_stack([np.ones(len(RK2)), relr])
RES = RK2 - Z @ np.linalg.lstsq(Z, RK2, rcond=None)[0]
base_W, _, _, base_X, _, _ = stats_from(RK2)
par_W, par_Ws, par_Wc, par_X, _, _ = stats_from(RES)
bv = np.array([np.median(base_W[D]) for D in DOMS]); pv = np.array([np.median(par_W[D]) for D in DOMS])
print(f"  {'域':10s} {'偏之前':>8s} {'偏之后':>8s} {'保留':>7s}")
for i, D in enumerate(DOMS):
    print(f"  {D:10s} {bv[i]:8.4f} {pv[i]:8.4f} {pv[i]/bv[i]*100:6.1f}%")
print(f"  {'跨域':10s} {np.median(base_X):8.4f} {np.median(par_X):8.4f} "
      f"{np.median(par_X)/np.median(base_X)*100:6.1f}%")
gapb = bv[DOMS.index("PURITY")] - np.max(np.delete(bv, DOMS.index("PURITY")))
gapp = pv[DOMS.index("PURITY")] - np.max(np.delete(pv, DOMS.index("PURITY")))
print(f"  **纯洁领先其余最高:{gapb:+.4f} -> {gapp:+.4f}(保留 {gapp/gapb*100:.1f}%)**")
# 安慰剂:偏掉一个与道德无关的变量(年龄)——「这个零该是零吗?」该是零
age = pd.to_numeric(d.loc[CC.index, "Age"], errors="coerce")[mrel]
ma = age.notna().values
ar = rankdata(age[ma].values).astype(float); ar = (ar-ar.mean())/ar.std()
Za = np.column_stack([np.ones(ma.sum()), ar])
RESa = RK2[ma] - Za @ np.linalg.lstsq(Za, RK2[ma], rcond=None)[0]
aW, _, _, aX, _, _ = stats_from(RESa)
av = np.array([np.median(aW[D]) for D in DOMS])
gapa = av[DOMS.index("PURITY")] - np.max(np.delete(av, DOMS.index("PURITY")))
print(f"  安慰剂(改偏年龄,n={int(ma.sum())}):纯洁领先 {gapa:+.4f}(保留 {gapa/gapb*100:.1f}%)"
      f" —— 必须≈100%")
G2 = Gate("纯洁域的高一致,是不是宗教?(MFQ,同一批人,同一格式)")
G2.negative_control("安慰剂:改偏年龄", null=abs(1-gapa/gapb), effect=abs(1-gapp/gapb),
                    null_spread=SPREAD, null_kind="与道德域结构无关的人口学变量")
G2.has_error_bar("纯洁领先(偏宗教后)", value=float(gapp), spread=SPREAD, spread_source="bootstrap_人层")
print(G2)
CONF = dict(rel_var=REL, n=int(mrel.sum()),
            before={D: float(bv[i]) for i, D in enumerate(DOMS)},
            after={D: float(pv[i]) for i, D in enumerate(DOMS)},
            gap_before=float(gapb), gap_after=float(gapp), gap_retain=float(gapp/gapb),
            placebo_age_retain=float(gapa/gapb), placebo_n=int(ma.sum()))

json.dump(dict(instrument="GrahamHaidtNosek 2009 JPSP Study 3 (MFQ)", n=len(d),
               domains=rows, cross=CROSS, cross_same_part=float(np.median(Xsame)),
               cross_cross_part=float(np.median(Xcross)), purity=PUR, other_max=OTH,
               within_median=float(np.median(wv)), within_range=float(wv.max()-wv.min()),
               spread=SPREAD, shuffle_null=SHUF, shuffle_sd=float(np.std(shuf)),
               catch_median=float(np.median(cat)), format_same=fs, format_cross=fc,
               spec_curve=spec, verdict=verdict, seeds=SEEDS, religion_confound=CONF,
               excluded_no_official_text=EXCLUDED, unchallenged=True),
          open(OUT/"domain_not_sex.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'domain_not_sex.json'}")
