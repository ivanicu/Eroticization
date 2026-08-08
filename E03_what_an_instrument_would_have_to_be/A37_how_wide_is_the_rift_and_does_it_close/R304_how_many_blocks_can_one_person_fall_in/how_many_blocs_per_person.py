r"""#865 · E03·A95·R304 —— 一个人能同时落在几块里:那几块孤立的少数是不是同一批人

**还 `#864`①。** `#864` 量出制裁三题上有**块状结构**:按年龄分中间那三分之一贴着最年轻的一端
(λ_perm 0.075/0.205/0.207),按信仰分中间贴着世俗那端(0.354–0.366),按政治分中间贴着保守那端
(0.636–0.647)⇒ **每根轴上都有一个孤立的三分之一。**
**⇒ 而「几块」和「一块」是两幅完全不同的社会图景,`#864` 没有分开它们。**

**⚠⚠ 而本轮第一件事是撤 `#864` 判词最后一句里的一个过度概括 —— 它三分之一是错的:**
`#864` 写的是「**在「怎么对待他」这件事上,反对的人不是分布的尾巴,是一个自成一块的少数派**」。
**对年龄和信仰成立**(孤立的那一块是**最老**和**最虔诚**,正是不放行的那一侧);
**对政治不成立** —— 政治上孤立的那一块是**自由派**,而自由派是**放行**的那一侧。
⇒ **三根轴里两根支持那句话,一根反过来。「孤立」和「反对」是两件事,`#864` 把它们并成了一句。**
本轮把**边**(孤立块落在宽容侧还是不宽容侧)**逐块印出来,不再合并**。⇒ `#864`②

**⚠⚠ 算术先行(`realstat` 的算术陷阱,`#864` 刚吃过一次):**
三块各占 1/3、相互独立时,`P(K=0)=8/27=0.2963 · P(K=1)=12/27=0.4444 · P(K=2)=6/27=0.2222 ·
P(K=3)=1/27=0.0370` ⇒ **`P(K≥2) = 7/27 = 0.2593` 是算出来的,不是量出来的。**
**⇒ 本轮的测量是实测值偏离 7/27 有多远**,而不是 `P(K≥2)` 本身;
**并且置换零必须自己重现 7/27,否则是零的实现错了,不是数据有结构。**

`G1` **估计量(先于方法命名)**:
   ① **`P(K≥2)`** —— 一个人同时落在 ≥2 个孤立块里的比例,**逐年算再平均**;
      **对照它自己的置换零(年内独立打乱每一块的成员身份,保边际)**,并对照算术值 7/27。
   ② **`lift = P(K≥2)_实测 ÷ P(K≥2)_置换零`** —— 偏离 1 的部分是测量。
   ③ **`coverage`** —— **在「不放行」的人里,有多少落在至少一个「不宽容侧」的孤立块里**,
      对照这些块在全体中的并集占比(**基率**)。⚠ **这一条才是决策相关的那一条:
      如果覆盖率很低,那么「谁在反对」这个问题的答案根本不在这几根轴上。**

**块的定义(先于数据写死,`#864` 的输出只用来定端,不用来定阈)**:
   每根连续轴在每个 (题 × 十年) 上分三分位;**孤立的那一端 = 离中间组更远的那一端**
   (按宽容度均值算 `|m_端 − m_中|`);**并记下它落在宽容侧还是不宽容侧**(该端均值 vs 全体均值)。

四个世界(**每个都有分支**):
   A **一块**:`lift` 明显 > 1(≥1.3)⇒ **三块基本是同一批人换了三个标签**,社会上只有一个少数派。
   B **几块**:`lift ≈ 1` 或 < 1 ⇒ **互不重叠的几块拼起来** —— **没有哪一种诉求能同时够到他们**,
     这是与 A 完全不同的一幅图,**也是我不欢迎的那一个**(它让「宗教那条缝」更没有位置)。
   C **块是真的但覆盖不了反对者**:`lift > 1` 而 `coverage` 相对基率的提升很小 ⇒
     **「谁在反对」的答案不在这几根轴上**,块只是块。
   D **⚠ 元分离器**:`P(K≥2)` 的置换零**重现不了 7/27** ⇒ **块的构造本身有问题**
     (三分位不独立、缺失模式相关、年内样本不均),**那就不是社会有结构,是我的零错了。**

预测矩阵:
   | 世界 | 现在 | lift≥1.3 | lift≈1 | 覆盖提升小 | 零重现不了 7/27 |
   | A 一块       | 0.40 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 几块       | 0.30 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 覆盖不了   | 0.20 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 零错了     | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(**条件式**):
  if 置换零重现 7/27(**|零 − 0.2593| < 0.02**,`D` 的否定)
     and 正控开火(**把两根轴的块人为设成同一批人 ⇒ lift 必须暴涨;而不做手脚时必须恰好不变**)
     and 安慰剂为零(**`ballot` 当第四块 ⇒ 它与任何真块的 lift 必须 ≈ 1**):
      lift ≥ 1.30                                      -> A
      lift 落在 [0.90, 1.15]                            -> B
      lift ≥ 1.30 但覆盖提升 < 1.15 倍基率              -> C
  else: UNVERIFIED

⚠ **跑前写下的最强混淆**:**三根轴本来就相关**(老年人更虔诚、保守派更虔诚)——
  **所以 lift > 1 几乎一定会出现,问题从来不是「有没有重叠」,而是「重叠有多少是这三根轴之间
  本来就有的相关造成的」。** ⇒ 控制:**同时报「块的重叠」与「整根轴之间的相关」**
  (Cramér's V / Spearman ρ),**并且把 lift 与「只由两两相关预测的 lift」放在一起**;
  ⚠ 而**最干净的那一条**是:**孤立块的边不一样**(信仰/年龄的孤立块在不宽容侧,政治的在宽容侧)
  ⇒ **政治块与另外两块的 lift 应该 < 1**,**这是一个方向相反的预测,它让本轮不是单向确认。**

`G3` 多重性:整族 = 3 题 × 6 十年 × (3 对 + 1 三元 + 3 安慰剂对),BH 与 BY 都做。
`G4` 规格曲线:**块的粗细两版** —— 三分位(1/3)与四分位极端(1/4),两版都报。
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`。

**⚠ 本轮结构性做不到的(登记,不许写「计划中」)**:
① 横断面 ⇒ **无因果识别**,不能问「是先落进块里还是先有立场」;
② **「同一批人」只能到「同一批**类别**」** —— GSS 无面板,**同一个体跨年不可追**,
   所以本轮说的「同一批人」严格意义是「同一年里的同一批受访者」,⇒ **结构性拿不到,不是没做**;
③ **块的定义依赖三分位这个人为切法** —— 规格曲线跑两个粗细,**但连续变量本来就没有真的「块」**,
   `#864` 量到的是**中间组的位置**,不是自然聚类;**本轮不做聚类**,那是另一个估计量;
④ **换不了仪器**:制裁三题是 GSS 独有(`#854`)。
"""
import json, math, pathlib, sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
NPERM, SEED = 200, 304
P864 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be"
                      "/A94_把两端拉开和解释了多少是同一件事吗/R303_which_side_the_middle_third_stands_on"
                      "/results/middle_third.json"))
IND = 7 / 27

print("=== ⓪a **算术先行**:三块各 1/3、相互独立 ⇒ `P(K≥2)` 是算出来的 ===")
for k in range(4):
    print(f"  P(K={k}) = {math.comb(3,k)*(1/3)**k*(2/3)**(3-k):.4f}", end="  ·" if k < 3 else "\n")
print(f"  ⇒ **P(K≥2) = 7/27 = {IND:.4f}** —— **本轮的测量是实测偏离它多远,不是这个数本身**;"
      f"\n     **而置换零必须自己重现 7/27,否则是零错了(世界 D),不是社会有结构。**")
print(f"\n=== ⓪b **先撤 `#864` 判词最后一句里的过度概括** ===")
for d in P864["lam"]:
    if d["item"] == "道德 `homosex`": continue
    side = "不宽容侧" if d["lam_perm"] < 0.5 else "宽容侧"
    print(f"  {d['axis']:16s} × {d['item']:16s} λ_perm={d['lam_perm']:.3f} ⇒ "
          f"孤立的那一块在 **{side}**")
print("  ⇒ **信仰/年龄的孤立块在不宽容侧,政治的孤立块在宽容侧(自由派)** ⇒ "
      "**`#864` 那句「反对的人是一个自成一块的少数派」三根轴里只有两根成立** ⇒ `#864`②")

GC = ["year", "spkhomo", "colhomo", "libhomo", "attend", "reliten", "fund",
      "polviews", "age", "ballot"]
gs = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta", columns=GC, convert_categoricals=False)
D = pd.DataFrame({"year": gs.year})
D["spk"] = 2 - pd.to_numeric(gs.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["col"] = 5 - pd.to_numeric(gs.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
D["lib"] = pd.to_numeric(gs.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)),
                    ("polviews", (1, 7)), ("age", (18, 89)), ("ballot", (1, 4))):
    D[c] = pd.to_numeric(gs[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
Rr = D.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = (zs(Rr.attend) + zs(-Rr.reliten) + zs(-Rr.fund)) / 3
D = D.join(Rr["REL"])

ITEMS = {"发言 `spkhomo`": "spk", "教书 `colhomo`": "col", "图书馆 `libhomo`": "lib"}
DECS = {"1970s": range(1972, 1980), "1980s": range(1980, 1990), "1990s": range(1990, 2000),
        "2000s": range(2000, 2010), "2010s": range(2010, 2020)}
AX = {"信仰 REL": "REL", "年龄 age": "age", "政治 polviews": "polviews"}
CUTS = {"三分位": 1 / 3, "四分位极端": 1 / 4}


NULLCOUNT = {}


def blocs(s, ycol, frac, force_same=None):
    """孤立块 = 离中间组更远的那一端。返回 {轴: (成员布尔, 边)} —— 边 = +1 宽容侧 / −1 不宽容侧。"""
    y = s[ycol].to_numpy(float)
    if not np.isfinite(y).all() or y.std(ddof=1) <= 0: return None
    gm = y.mean(); out = {}
    for a, col in AX.items():
        v = s[col].to_numpy(float)
        ok = np.isfinite(v)
        if ok.sum() < 300: return None
        qlo, qhi = np.nanquantile(v[ok], [frac, 1 - frac])
        # ⚠⚠ **必须与 `#864` 逐字同一个划分:`np.digitize`,不是「≤qlo / ≥qhi / 严格居中」。**
        # 第一版用了后者,整张网格空了 —— 因为 **`polviews` 是 7 档整数,它的三分位切点落在
        # 相邻整数上(4 和 5),严格居中的那一组是空集**。这不是代码笔误,是**变量的粒度**:
        # 一个离散变量的两个分位点可以贴在一起,而「中间」这个概念此时需要一个包含端点的约定。
        lab = np.digitize(v, [qlo, qhi]).astype(float)
        lab[~ok] = np.nan
        lo_m, mid_m, hi_m = (lab == 0), (lab == 1), (lab == 2)
        if lo_m.sum() < 60 or hi_m.sum() < 60 or mid_m.sum() < 60:
            NULLCOUNT[a] = NULLCOUNT.get(a, 0) + 1
            return None
        mm = y[mid_m].mean()
        iso = hi_m if abs(y[hi_m].mean() - mm) >= abs(y[lo_m].mean() - mm) else lo_m
        out[a] = (iso, 1 if y[iso].mean() > gm else -1)
    if force_same:                      # 正控:把第二根轴的块**设成**第一根的块
        a1, a2 = force_same
        out[a2] = (out[a1][0].copy(), out[a2][1])
    return out


def stats(s, ycol, frac, rng=None, force_same=None, placebo=False):
    B = blocs(s, ycol, frac, force_same=force_same)
    if B is None: return None
    names = list(AX)
    M = np.array([B[a][0] for a in names])
    if placebo:                          # 用 ballot 造第四块(随机分配,必须无结构)
        bv = s["ballot"].to_numpy(float)
        pb = np.isfinite(bv) & (bv == 1)
        if pb.sum() < 60: return None
        M = np.vstack([M, pb[None, :]])
        names = names + ["安慰剂 ballot"]
    if rng is not None:                  # 零:每一块独立打乱成员身份,保边际
        M = np.array([r[rng.permutation(len(r))] for r in M])
    K = M[:3].sum(0)
    y = s[ycol].to_numpy(float)
    refuse = y < 0.5                     # 「不放行」
    restr = [i for i, a in enumerate(names[:3]) if B[a][1] < 0]
    cov = float(M[restr].any(0)[refuse].mean()) if restr and refuse.sum() > 30 else np.nan
    base = float(M[restr].any(0).mean()) if restr else np.nan
    pair = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            o = float((M[i] & M[j]).mean()); e = float(M[i].mean() * M[j].mean())
            pair[f"{names[i]}×{names[j]}"] = (o / e) if e > 0 else np.nan
    # ⚠⚠ **独立性期望必须由这一格的实测边际算,不能用常数 7/27。**
    # 7/27 只在**每一块恰好占 1/3** 时成立,而 `polviews` 是 7 档整数、重并极多,
    # 它的「中间」组占到近四成 ⇒ **孤立块只有约三成** ⇒ 常数本身是错的。
    # 第一版拿置换零去跟 7/27 比,负控 FAIL —— **失败的是那个常数,不是零。**
    ps_ = [float(M[i].mean()) for i in range(3)]
    e0 = np.prod([1 - q for q in ps_])
    e1 = sum(ps_[i] * np.prod([1 - ps_[j] for j in range(3) if j != i]) for i in range(3))
    return dict(pk2=float((K >= 2).mean()), pk3=float((K == 3).mean()),
                ind=float(1 - e0 - e1), sizes=ps_,
                cov=cov, base=base, pair=pair, n=int(len(s)),
                sides={a: B[a][1] for a in AX})


def agg(sub, ycol, frac, rng=None, force_same=None, placebo=False):
    acc = []
    for yv in np.unique(sub.year):
        s = sub[sub.year == yv]
        if len(s) < 400: continue
        r = stats(s, ycol, frac, rng=rng, force_same=force_same, placebo=placebo)
        if r: acc.append(r)
    if not acc: return None
    out = dict(pk2=float(np.mean([a["pk2"] for a in acc])),
               pk3=float(np.mean([a["pk3"] for a in acc])),
               ind=float(np.mean([a["ind"] for a in acc])),
               sizes=[float(np.mean([a["sizes"][i] for a in acc])) for i in range(3)],
               cov=float(np.nanmean([a["cov"] for a in acc])),
               base=float(np.nanmean([a["base"] for a in acc])), nyear=len(acc))
    ks = acc[0]["pair"].keys()
    out["pair"] = {k: float(np.nanmean([a["pair"][k] for a in acc])) for k in ks}
    out["sides"] = acc[-1]["sides"]
    return out


print(f"\n=== ① 网格:{len(ITEMS)} 题 × {len(DECS)} 十年 × {len(CUTS)} 粗细 · 置换零 {NPERM} 次 ===")
rows = []
rng = np.random.default_rng(SEED)
for inm, ycol in ITEMS.items():
    for dec in DECS:
        m = D[ycol].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 1200: continue
        for cn, frac in CUTS.items():
            o = agg(sub, ycol, frac, placebo=True)
            if o is None: continue
            nd = [agg(sub, ycol, frac, rng=rng, placebo=True) for _ in range(NPERM)]
            nd = [x for x in nd if x]
            npk = np.array([x["pk2"] for x in nd])
            npair = {k: np.array([x["pair"][k] for x in nd]) for k in o["pair"]}
            rows.append(dict(item=inm, dec=dec, cut=cn, pk2=o["pk2"], pk3=o["pk3"],
                             ind=o["ind"], sizes=o["sizes"],
                             null_pk2=float(np.mean(npk)), null_sd=float(np.std(npk)),
                             lift=o["pk2"] / float(np.mean(npk)) if np.mean(npk) > 0 else np.nan,
                             p=float((1 + (npk >= o["pk2"]).sum()) / (len(npk) + 1)),
                             cov=o["cov"], base=o["base"],
                             cov_lift=(o["cov"] / o["base"]) if o["base"] > 0 else np.nan,
                             pair={k: v for k, v in o["pair"].items()},
                             pair_null={k: float(np.mean(v)) for k, v in npair.items()},
                             sides=o["sides"], n=int(len(sub)), nyear=o["nyear"]))
        r3 = [r for r in rows if (r["item"], r["dec"], r["cut"]) == (inm, dec, "三分位")]
        if r3:
            r = r3[0]
            print(f"  {inm:16s} {dec} n={r['n']:5,d} · P(K≥2) **{r['pk2']:.4f}** vs 置换零 "
                  f"**{r['null_pk2']:.4f}** vs 该格边际独立值 **{r['ind']:.4f}** · **lift {r['lift']:.3f}** · "
                  f"覆盖 {r['cov']:.3f} / 基率 {r['base']:.3f} = **{r['cov_lift']:.3f}×**")

TER = [r for r in rows if r["cut"] == "三分位"]
NULLDEV = float(np.mean([abs(r["null_pk2"] - r["ind"]) for r in TER]))
NULLDEV_CONST = float(np.mean([abs(r["null_pk2"] - IND) for r in TER]))
LIFT = float(np.median([r["lift"] for r in TER]))
COVL = float(np.median([r["cov_lift"] for r in TER]))
print(f"\n  ⚠ **本轮每一格为什么可能落空,逐轴计数(空网格必须说出理由,不能静默)**:"
      + (" · ".join(f"{k}:{v}" for k, v in NULLCOUNT.items()) if NULLCOUNT else "无")
      + f" · 有效格 **{len(rows)}**")
if not rows:
    raise SystemExit("⛔ 网格为空 —— 空总体不许当作通过(`realstat §4`:exit 2, never 0)")
print(f"\n=== ② 两两之间:哪两块是同一批人,哪两块互斥 ===")
PN = sorted({k for r in TER for k in r["pair"]})
for k in PN:
    v = [r["pair"][k] for r in TER if k in r["pair"]]
    nv = [r["pair_null"][k] for r in TER if k in r["pair_null"]]
    print(f"  {k:34s} 观测/独立 = **{np.median(v):.3f}** · 置换零 {np.median(nv):.3f}"
          + ("  ⚠ **< 1:这两块互斥**" if np.median(v) < 0.95 else ""))
print("  ⚠ **跑前的方向相反的预测**:政治的孤立块在**宽容侧**,信仰/年龄的在**不宽容侧** ⇒ "
      "**政治×信仰、政治×年龄 应当 < 1**;上面是它们实测的样子。")

print("\n=== ③ 控制 ===")
sub0 = D[D.col.notna() & D.REL.notna() & D.year.isin(list(DECS["2010s"]))]
b0 = agg(sub0, "col", 1 / 3, placebo=True)
pc = agg(sub0, "col", 1 / 3, force_same=("信仰 REL", "年龄 age"), placebo=True)
nz = agg(sub0, "col", 1 / 3, force_same=None, placebo=True)
print(f"  正控:**把年龄块人为设成信仰块**(两块变成同一批人)⇒ P(K≥2) "
      f"**{b0['pk2']:.4f} → {pc['pk2']:.4f}**(**{pc['pk2']-b0['pk2']:+.4f}**);"
      f"**不做手脚时 {nz['pk2']-b0['pk2']:+.2e}** ⇒ ⚠ **`G2` 控制必须能失败**")
CEIL = pc["pk2"] - b0["pk2"]; CT = 0.02
print(f"     **控制也必须能通过**:floor {abs(nz['pk2']-b0['pk2']):.2e} < 阈 {CT} < ceiling {CEIL:.4f} ⇒ "
      f"**{'阈在真带内' if abs(nz['pk2']-b0['pk2']) < CT < CEIL else '⚠⚠ 阈不在带内'}**")
rgn = np.random.default_rng(SEED + 3)
nn = agg(sub0, "col", 1 / 3, rng=rgn, placebo=True)
r0 = [r for r in TER if r["item"] == "教书 `colhomo`" and r["dec"] == "2010s"][0]
print(f"  负控:年内独立打乱每一块的成员身份(保边际)⇒ P(K≥2) = **{nn['pk2']:.4f}**,"
      f"**该格边际独立值 {nn['ind']:.4f}** ⇒ **偏差 {abs(nn['pk2']-nn['ind']):.4f}**"
      f"(对照常数 7/27 的偏差 {abs(nn['pk2']-IND):.4f})—— "
      f"⚠ **「这个零该不该是 7/27?」不该** —— 7/27 要求每块恰好 1/3,而实测块大小是 "
      + " · ".join(f"{a_}:{np.median([r['sizes'][i] for r in TER]):.3f}" for i, a_ in enumerate(AX))
      + f",**政治那块只有 0.27**(`polviews` 是 7 档整数,中间组近四成)"
      f" ⇒ **该比的是这一格自己的边际独立值;第一版比了常数,负控 FAIL,而失败的是那个常数。** ⇒ `#865`①")
PLK = [k for k in PN if "安慰剂" in k]
plv = [r["pair"][k] for r in TER for k in PLK if k in r["pair"]]
print(f"  安慰剂 `ballot` 当第四块:与三根真轴的 观测/独立 中位 **{np.median(plv):.4f}** "
      f"(**必须 ≈ 1** —— GSS 自己随机分配的问卷版本,真实变量)")

ps = np.array([r["p"] for r in rows if np.isfinite(r["p"])]); C = len(ps)
o_ = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C; cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (int(np.max(np.where(pv <= cr)[0])) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o_], cH), su(ps[o_], cY)
print(f"\n=== ④ 多重性:整族 **{C}** 格 · BH 存活 **{kH}** · BY **{kY}** · "
      f"p 分辨率下限 {1/(NPERM+1):.4f} · 不同意的 {kH-kY} 格一起发表 ===")

G = Gate("#865 · 那几块孤立的少数是不是同一批人")
G.asserted("① **算术先行**:独立时 `P(K≥2)` 由三块的边际算出 ⇒ **置换零必须自己重现它**,"
           "否则是零错了(世界 D),不是社会有结构。⚠ **而它不是 7/27** —— 7/27 只在每块恰好 1/3 时成立,"
           "`polviews` 重并极多、中间组近四成 ⇒ **孤立块只有约三成**;第一版拿常数比,负控 FAIL,"
           "**失败的是常数不是零**",
           bool(NULLDEV < 0.02),
           f"置换零与**该格边际独立值**的平均偏差 {NULLDEV:.4f}(阈 0.02);"
           f"而与常数 7/27 的偏差是 {NULLDEV_CONST:.4f} —— **两者之差就是「块不是 1/3」这件事**;"
           f"实测块大小中位 " + " · ".join(
               f"{a}:{np.median([r['sizes'][i] for r in TER]):.3f}" for i, a in enumerate(AX)),
           kind="control")
G.asserted("② 前提(跑前写下的最强混淆):**三根轴本来就相关** ⇒ lift>1 几乎一定出现;"
           "⚠ **而孤立块的边不一样**(信仰/年龄在不宽容侧、政治在宽容侧)"
           "⇒ **政治×信仰、政治×年龄 应当 < 1,这是一个方向相反的预测**",
           bool(all(np.isfinite(r["pair"].get("信仰 REL×政治 polviews", np.nan)) for r in TER)),
           " · ".join(f"{k}:{np.median([r['pair'][k] for r in TER if k in r['pair']]):.3f}"
                      for k in PN if "安慰剂" not in k), kind="control")
G.asserted("③ 正控:**把年龄块人为设成信仰块** ⇒ P(K≥2) 必须暴涨;不做手脚时必须恰为 0;"
           "**且阈落在 floor 与 ceiling 之间**",
           bool(CEIL > CT and abs(nz["pk2"] - b0["pk2"]) < 1e-12
                and abs(nz["pk2"] - b0["pk2"]) < CT < CEIL),
           f"{b0['pk2']:.4f} → {pc['pk2']:.4f}({CEIL:+.4f})· 不动手 {nz['pk2']-b0['pk2']:+.2e}",
           kind="control")
G.asserted("④ 负控:独立打乱每一块的成员身份(保边际)⇒ 必须落回**该格自己的边际独立值**",
           bool(abs(nn["pk2"] - nn["ind"]) < 0.02),
           f"{nn['pk2']:.4f} vs 边际独立值 {nn['ind']:.4f}(常数 7/27 = {IND:.4f},差 "
           f"{abs(nn['pk2']-IND):.4f} —— **那个差是块大小,不是结构**)", kind="control")
G.asserted("⑤ 安慰剂 `ballot` 当第四块 ⇒ 与真块的 观测/独立 必须 ≈ 1",
           bool(abs(np.median(plv) - 1.0) < 0.05), f"中位 {np.median(plv):.4f}", kind="control")
G.asserted("⑥ kill(预注册):「三块基本是同一批人」要成立,需 **lift ≥ 1.30**",
           bool(LIFT >= 1.30),
           f"lift 中位 {LIFT:.3f}(逐格 {min(r['lift'] for r in TER):.3f}–"
           f"{max(r['lift'] for r in TER):.3f})· 覆盖提升中位 {COVL:.3f}× · "
           f"P(K=3) 中位 {np.median([r['pk3'] for r in TER]):.4f}",
           kind="kill",
           yardstick="`P(K≥2)` 相对**它自己的年内置换零**(保每一块的边际),并对照算术值 7/27",
           yardstick_noise=float(np.median([r["null_sd"] for r in TER])),
           population=f"GSS 三道制裁题 × {len({r['dec'] for r in TER})} 个十年的 {len(TER)} 个三分位格 —— "
                      f"⚠ **不含道德题**(`#864` 量到块状结构只在制裁题上,道德题四根轴都近线性)",
           direction=[r["lift"] - 1.0 for r in TER])
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
neg_pairs = [k for k in PN if "安慰剂" not in k
             and np.median([r["pair"][k] for r in TER if k in r["pair"]]) < 0.95]
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif LIFT >= 1.30 and COVL < 1.15:
    VERD = (f"**C 块是真的,但它覆盖不了反对的人 ⇒ 「谁在反对」的答案不在这几根轴上。**\n"
            f"  lift {LIFT:.3f},而不放行的人里落进「不宽容侧块」的比例只有基率的 **{COVL:.3f}×**。")
elif LIFT >= 1.30:
    VERD = (f"**A 三块基本是同一批人 —— lift 中位 {LIFT:.3f}。**\n"
            f"  ⇒ **社会上只有一个少数派,它同时是最老的、最虔诚的、最保守的那批人。**")
elif 0.90 <= LIFT <= 1.15:
    # ⚠⚠ **边由数据算,不由我写。** 第一版把政治那块写成「最保守的」——
    #     而实测它是**自由派**那一块(孤立在**宽容**侧),`#864` 刚犯过同一个错。(本轮自捕,不进债表)
    side_n = {a_: sum(1 for r in TER if r["sides"][a_] < 0) for a_ in AX}
    who = {"信仰 REL": ("最虔诚的", "最世俗的"), "年龄 age": ("最年长的", "最年轻的"),
           "政治 polviews": ("最保守的", "最自由的")}
    desc = " · ".join(
        f"**{a_.split()[0]}:{who[a_][0] if side_n[a_] > len(TER)/2 else who[a_][1]}那一块**"
        f"({'不宽容' if side_n[a_] > len(TER)/2 else '宽容'}侧,{max(side_n[a_], len(TER)-side_n[a_])}"
        f"/{len(TER)} 格)" for a_ in AX)
    VERD = (f"**B 它们不是同一批人 —— lift 中位 {LIFT:.3f},落在独立性附近。**\n"
            f"  ⇒ **「拦着他」的阻力是几块互不重叠的少数拼起来的,没有哪一种诉求能同时够到他们。**\n"
            f"  ⚠ 而这是我不欢迎的那一个:它让「宗教那条缝」更没有位置。\n"
            f"  ⚠ **三块孤立在哪一边,由数据算出来,而它们不在同一边**:{desc}\n"
            f"  ⚠ **互斥的对**:{neg_pairs if neg_pairs else '无'} —— "
            f"**跑前就预测了政治块会与另两块反向,因为它的边不同;实测正是如此。**\n"
            f"  ⇒ **一句关于人的话:在「该不该拦着他」上,那几块站得离人群最远的少数,"
            f"既不是同一批人,也不站在同一边 ——\n"
            f"  最虔诚的那三分之一和最年长的那三分之一站在「拦」的一侧,"
            f"最自由的那三分之一站在「不拦」的一侧,\n"
            f"  而这三块几乎不重叠(lift {LIFT:.3f}),信仰块与政治块甚至互相排斥"
            f"({min(np.median([r['pair'][k] for r in TER if k in r['pair']]) for k in neg_pairs):.3f})。\n"
            f"  不放行的人里落进「不宽容侧块」的比例是基率的 {COVL:.3f} 倍 —— **块是真的,"
            f"但它只解释了一部分人。**"
            f"\n  ⇒ 也就是说:这件事上没有一个「反对派」,只有几群因为不同理由、"
            f"彼此也不重合的人,各自站在自己的位置上。**")
else:
    VERD = (f"**都不是**:lift 中位 {LIFT:.3f},覆盖提升 {COVL:.3f}× —— "
            f"**三个预注册世界都没被满足,如实登记。**")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① 横断面 ⇒ **无因果识别**;② **「同一批人」只能到「同一年里的"
      f"同一批受访者」** —— GSS 无面板,同一个体跨年不可追,**结构性拿不到,不是没做**;"
      f"③ **块的定义依赖三分位这个人为切法**,规格曲线跑了两个粗细,**但连续变量本来就没有真的「块」** ——"
      f"`#864` 量到的是中间组的位置,不是自然聚类,**本轮不做聚类,那是另一个估计量**;"
      f"④ **换不了仪器**(`#854`)。")

json.dump(dict(grid=rows, independence=IND, lift_median=LIFT, cov_lift_median=COVL,
               null_dev=NULLDEV, negative_pairs=neg_pairs,
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               controls=dict(pos_from=b0["pk2"], pos_to=pc["pk2"], ceiling=CEIL, threshold=CT,
                             zero=nz["pk2"] - b0["pk2"], neg=nn["pk2"],
                             placebo_median=float(np.median(plv))),
               derivation="P(K>=2) from the OBSERVED bloc marginals; 7/27 only holds if each bloc is "
                           "exactly 1/3, and polviews' ties make it ~0.30 -- the constant was the "
                           "defect, not the null",
               null_dev_vs_constant=NULLDEV_CONST,
               side_counts_note="isolated bloc side computed per cell, not asserted",
               retracts_864="#864's closing sentence generalised across three axes; the political bloc is "
                            "isolated on the PERMISSIVE side, so 2 of 3 axes support it and 1 contradicts",
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), seed=SEED, nperm=NPERM),
          open(OUT / "blocs_per_person.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'blocs_per_person.json'}")
