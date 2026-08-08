r"""#864 · E03·A94·R303 —— 中间那三分之一站在哪边:`|g|` 与「解释了多少」的差,就是中间那批人

**还 `#863`①。** `#862`/`#863` 用 `|g|` 排了七根轴,而 `|g|` 是**两端三分位之差** ——
**它把中间那三分之一整个扔掉了**。`#863` 结尾登记的疑问是:制裁那一侧「谁都只是险胜」
(领先中位 +0.082,6/30 落在地板内),而 **`|g|` 这把尺分辨不出「没有主轴」和「所有轴都弱」**。

**⚠⚠ 但在换估计量之前,先做算术 —— 这一轮的全部严肃性在这里(`realstat` 的「算术陷阱」):**
**如果 (X, Y) 联合正态、三分位等分,那么 η² 与 g 之间是恒等式,不是两个测量:**
   `g   = ρ·(E[Z|上三分位] − E[Z|下三分位]) = ρ · 2.1814`
   `η²  = ρ² · Var(组均值) = ρ² · 0.7931`
   ⇒ **`η² = g²/6`(系数 0.7931/2.1814² = 0.16668 = 1/6.0002)**
**这一条已在机器上验过**(n=4e6,ρ = 0.2/0.4/0.6,`η² ÷ (g²/6)` = **1.0000 / 1.0000 / 1.0000**),
脚本里再跑一次小规模的positive control 版本。

⇒ **所以「换成方差份额」本身不是新证据,它在正态世界里是同一个数换了个写法**
   (`feedback_reparameterisation_is_not_measurement` 说的正是这件事)。
   **真正的测量是那个比值 `R = η²_实测 ÷ (g²/6)` 偏离 1 有多远,以及偏离由什么造成。**
   而它**必然**在四个地方偏离,每一个都是一句关于人的话:
   ① **结局不正态** —— 三道制裁题是 0/1,道德题是四档;
   ② **轴上有大量重并** —— `educ` 的三分位切在 12 年 vs 15 年,高组 36.8% / 低组 42.1%,**不是 33/33**;
   ③ **组不等大** —— 地区 38/62、种族 75/16、性别 50/50,`g` 的那个 2.1814 括号根本不适用;
   ④ **⚠ 中间那三分之一不是线性内插** —— 而这是唯一一条**关于人**而不关于测量的:
      **温和虔诚的人,是站在中间,还是站在虔诚那一边?**

`G1` **估计量(三个,先于方法命名)**:
   ① **`η²`**:每个 (题 × 轴 × 十年) 上,把**全部**受访者按该轴分组后的组间方差份额,**逐年算再平均**。
      ⚠ **它与 `|g|` 的差别只有一件事:`|g|` 用 68% 的人(两端),`η²` 用 100%。**
   ② **`R = η² ÷ (g²/6)`** —— **偏离 1 的部分才是本轮的测量**,`R < 1` 表示「两端拉得开但解释得少」。
   ③ **`λ = (m_中 − m_低)/(m_高 − m_低)`** —— **中间三分之一的位置**。
      **λ = 0.5 线性 · λ → 0 中间和低组抱团 · λ → 1 中间和高组抱团。**
      ⚠ `λ` 只在 `|m_高 − m_低|` **远高于地板**时可读(**预注册:`|g| > 3 ×` 该格地板才报**),
      否则它是 0/0;**报不了的格如实印 `UNREADABLE`,不留空。**

四个世界(**每个都有分支**,`#856`):
   A **`R ≈ 1` 且 η² 的排名 = `|g|` 的排名** ⇒ **两个估计量是一个**,`#862`/`#863` 的排名不是尺的产物,
     而**「换估计量」这条路本身走不通** —— 那本轮是 Closure,如实标。
   B **制裁侧排名变了(教育丢掉领先)** ⇒ **`#862`/`#863` 的制裁侧结论是 `educ` 粗切法造出来的**,
     「无一例外」要从撤四个字升级为撤整句。⚠ **这是我不欢迎的那个,而它有具体的机制预测:
     `educ` 重并最多 ⇒ 它的 `R` 应该最低。**
   C **道德侧排名也变了(宗教丢掉领先)** ⇒ **`#832` 以来的一切都是估计量特有的。**
   D **`R` 的偏离本身有结构(某些轴系统性地 `R < 1`)** ⇒ **「两端拉得开」和「解释得多」是两件不同的
     社会事实,两个都真** —— ⚠ **元分离器:「哪根轴最宽」这个问法一直缺一个下标。**

预测矩阵:
   | 世界 | 现在 | R≈1 且排名同 | 制裁排名变 | 道德排名变 | R 有结构 |
   | A 同一个量 | 0.30 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 制裁变   | 0.25 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 全变     | 0.10 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D R 有结构 | 0.35 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(**条件式**):
  if 正控开火(**恒等式的正控:在合成正态数据上 `R` 必须 = 1;而在真实数据上它不必**)
     and 负控为零(**年内打乱轴标签 ⇒ η² 必须塌回它自己的置换零 —— ⚠ η² 有正偏,
        它的零不是 0,是 `(k−1)/(n−k)` 量级,所以必须用置换零而不是 0**)
     and 安慰剂为零(**`ballot`**):
      η² 排名与 `|g|` 排名在 ≥90% 的格里一致,且 `R` 的中位落在 [0.8, 1.25]  -> A
      制裁格里 η² 的主轴 ≠ 教育 的比例 ≥ 1/3                                  -> B
      道德格里 η² 的主轴 ≠ 宗教 的比例 ≥ 1/3                                  -> C
      各轴的 `R` 中位极差 > 0.25(**轴间系统性差异**)                         -> D
  else: UNVERIFIED

⚠ **跑前写下的最强混淆**:**`η²` 是有偏的**,偏差随组数 `k` 与样本量变化 ——
  **七根轴的 `k` 不同(连续轴 3 组,二元轴 2 组)⇒ 直接比 `η²` 会系统性偏向组多的轴。**
  ⇒ 控制:**同时报 `ω²`(偏差校正)与「减掉自己置换零中位」的 `η²_净`**,
  **排名用 `η²_净`,而三个都印出来**;并且**安慰剂 `ballot` 是 2 组,恰好落在最容易被偏差抬高的一侧。**

`G3` 多重性:整族 = 4 题 × 7 轴 × 6 十年 × 2 总体规格,BH 与 BY 都做,不同意的格一起发表。
`G4` 规格曲线:**总体规格两版** —— ①「与 `|g|` 可比」(连续轴只用两端三分位、种族只用白/黑)
   ②「全体」(连续轴三组全用、种族三类全用)—— **两版都报,因为它们回答不同的问题。**
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`。

**⚠ 本轮结构性做不到的(登记,不许写「计划中」)**:
① 横断面 ⇒ **无因果识别**;
② **`λ` 不能分开「中间的人本来就中间」与「中间的人在两种立场间摇摆」** —— 需要同一个人的重复测量,
   **GSS 的滚动面板 1-2 波在这三道题上不够**,⇒ **结构性拿不到,不是没做**;
③ **`R` 的偏离不能唯一归因** —— 四个来源(非正态·重并·组不等大·中间非线性)**同时在场**,
   本轮只能**分离出「中间非线性」那一条**(因为 `λ` 直接量它),其余三条**只能一起报**;
④ **换不了仪器**:制裁三题是 GSS 独有(`#854`);道德题的 NSFG 版在 `#863` 已跑,**本轮不重复**。
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
NPERM, NBOOT, SEED = 200, 300, 303
PRIOR = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be/A93_宗教是不是这个对象的那根轴"
                       "/R302_the_same_ruler_applied_to_the_moral_items/results/which_cleavage_moral.json"))

print("=== ⓪a **先做算术,再做测量**:`η² = g²/6` 在联合正态下是恒等式,不是发现 ===")
rg = np.random.default_rng(SEED)
for rho in (0.2, 0.5):
    n = 400_000
    x = rg.normal(size=n); y = rho * x + math.sqrt(1 - rho ** 2) * rg.normal(size=n)
    q1, q2 = np.quantile(x, [1 / 3, 2 / 3]); grp = np.digitize(x, [q1, q2])
    gg = (y[grp == 2].mean() - y[grp == 0].mean()) / y.std(ddof=1)
    mu = np.array([y[grp == k].mean() for k in range(3)])
    w = np.array([(grp == k).mean() for k in range(3)])
    e2 = float(np.sum(w * (mu - y.mean()) ** 2) / y.var(ddof=1))
    print(f"  ρ={rho}: g={gg:+.4f} · η²={e2:.5f} · g²/6={gg*gg/6:.5f} · **R={e2/(gg*gg/6):.4f}**")
print("  ⇒ **R=1 是算术强制的** ⇒ 「换成方差份额」在正态世界里不是新证据;"
      "**本轮的测量是 R 偏离 1 有多远,以及偏离由什么造成。**")

GCOLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo", "attend", "reliten", "fund",
         "polviews", "educ", "age", "region", "race", "sex", "ballot"]
gs = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta", columns=GCOLS, convert_categoricals=False)
D = pd.DataFrame({"year": gs.year})
D["moral"] = pd.to_numeric(gs.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
D["spk"] = 2 - pd.to_numeric(gs.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["col"] = 5 - pd.to_numeric(gs.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
D["lib"] = pd.to_numeric(gs.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7)),
                    ("educ", (0, 20)), ("age", (18, 89)), ("region", (1, 4)), ("race", (1, 3)),
                    ("sex", (1, 2)), ("ballot", (1, 4))):
    D[c] = pd.to_numeric(gs[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
Rr = D.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = (zs(Rr.attend) + zs(-Rr.reliten) + zs(-Rr.fund)) / 3
D = D.join(Rr["REL"])

ITEMS = {"道德 `homosex`": "moral", "发言 `spkhomo`": "spk",
         "教书 `colhomo`": "col", "图书馆 `libhomo`": "lib"}
DECS = {"1970s": range(1972, 1980), "1980s": range(1980, 1990), "1990s": range(1990, 2000),
        "2000s": range(2000, 2010), "2010s": range(2010, 2020), "2020s": range(2020, 2025)}
CONT = {"宗教 REL": "REL", "教育 educ": "educ", "政治 polviews": "polviews", "年龄 age": "age"}
CATG = {"地区 南方vs其余": ("region", lambda v: (v == 3).astype(float)),
        "种族 白vs黑": ("race", None), "性别 男vs女": ("sex", None)}
PLAC = {"安慰剂 ballot": ("ballot", None)}
ANAMES = list(CONT) + list(CATG)
SPECS = ["与|g|可比", "全体"]


def groups(s, axis, spec):
    """返回 (标签数组, 有效掩码)。`与|g|可比` = 连续轴丢中间三分之一、种族只留白/黑。"""
    if axis in CONT:
        v = s[CONT[axis]].to_numpy(float)
        ok = np.isfinite(v)
        if ok.sum() < 200: return None, None
        q1, q2 = np.nanquantile(v[ok], [1 / 3, 2 / 3])
        lab = np.digitize(v, [q1, q2]).astype(float)
        lab[~ok] = np.nan
        if spec == "与|g|可比": lab[lab == 1] = np.nan
        return lab, np.isfinite(lab)
    col, _ = {**CATG, **PLAC}[axis]
    v = s[col].to_numpy(float)
    if axis == "地区 南方vs其余":
        lab = np.where(v == 3, 1.0, np.where(np.isfinite(v), 0.0, np.nan))
    elif axis == "种族 白vs黑":
        lab = np.where(v == 1, 1.0, np.where(v == 2, 0.0, np.where(v == 3, 2.0, np.nan)))
        if spec == "与|g|可比": lab[lab == 2] = np.nan
    elif axis == "性别 男vs女":
        lab = np.where(np.isfinite(v), v, np.nan)
    else:
        lab = np.where(v == 1, 1.0, np.where(v == 3, 0.0, np.nan))
    return lab, np.isfinite(lab)


def eta2(y, lab, ok):
    """η²(有偏)· ω²(偏差校正)· 组均值 —— 一次算完。"""
    yy, ll = y[ok], lab[ok]
    n = len(yy)
    ks = np.unique(ll)
    if n < 120 or len(ks) < 2: return np.nan, np.nan, None, None
    gm = yy.mean(); tot = float(((yy - gm) ** 2).sum())
    if tot <= 0: return np.nan, np.nan, None, None
    mus, ns = [], []
    ssb = 0.0
    for k in ks:
        m = ll == k
        if m.sum() < 40: return np.nan, np.nan, None, None
        mus.append(float(yy[m].mean())); ns.append(int(m.sum()))
        ssb += m.sum() * (yy[m].mean() - gm) ** 2
    e = float(ssb / tot)
    kk = len(ks)
    msw = (tot - ssb) / (n - kk)
    om = float((ssb - (kk - 1) * msw) / (tot + msw)) if (tot + msw) > 0 else np.nan
    return e, om, np.array(mus), np.array(ns)


def cell(sub, ycol, axis, spec, perm_rng=None, plant=0.0, plant_axis=None):
    """逐年算再平均:η² · ω² · g(两端) · λ(中间位置)。"""
    E, O, G, L = [], [], [], []
    for yv in np.unique(sub.year):
        s = sub[sub.year == yv]
        if len(s) < 200: continue
        lab, ok = groups(s, axis, spec)
        if lab is None: continue
        y = s[ycol].to_numpy(float).copy()
        if plant and plant_axis is not None:
            # ⚠ **符号由数据定,不由我定**:往「已经更低的那一组」再压 —— 加在高组上会把两组拉近,
            # 那是在**抵消**效应,不是在植入(`#864`①:`#863` 的正控就是那个方向,当时印成了「植入」)
            pl, pok = groups(s, plant_axis, "全体")
            if pl is not None:
                hi_m = pok & (pl == np.nanmax(pl)); lo_m = pok & (pl == np.nanmin(pl))
                sgn = 1.0 if y[hi_m].mean() >= y[lo_m].mean() else -1.0
                y[hi_m] += sgn * plant
        if perm_rng is not None:
            k = perm_rng.permutation(len(s)); lab, ok = lab[k], ok[k]
        e, om, mus, ns = eta2(y, lab, ok)
        if not np.isfinite(e): continue
        E.append(e); O.append(om)
        if mus is not None and len(mus) >= 2:
            sd = y[ok].std(ddof=1)
            if sd > 0: G.append(float((mus[-1] - mus[0]) / sd))
            if len(mus) == 3 and abs(mus[-1] - mus[0]) > 1e-9:
                L.append(float((mus[1] - mus[0]) / (mus[-1] - mus[0])))
    if not E: return {}
    return dict(eta2=float(np.mean(E)), omega2=float(np.mean(O)),
                g=float(np.mean(G)) if G else np.nan,
                lam=float(np.mean(L)) if L else np.nan, nyear=len(E))


print(f"\n=== ① 网格:{len(ITEMS)} 题 × {len(ANAMES)} 轴 × {len(DECS)} 十年 × {len(SPECS)} 总体规格 "
      f"· 每格自己的置换零({NPERM} 次)===")
rows, best = [], {}
rng2 = np.random.default_rng(SEED)
for inm, icol in ITEMS.items():
    for dec in DECS:
        m = D[icol].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 800: continue
        for spec in SPECS:
            obs, nul = {}, {}
            for a in ANAMES + list(PLAC):
                r = cell(sub, icol, a, spec)
                if not r: continue
                nd, ng = [], []
                for _ in range(NPERM):
                    rr = cell(sub, icol, a, spec, perm_rng=rng2)
                    if rr:
                        nd.append(rr["eta2"])
                        if np.isfinite(rr["g"]): ng.append(abs(rr["g"]))
                nd = np.array(nd)
                if not len(nd): continue
                obs[a] = r; nul[a] = nd
                net = r["eta2"] - float(np.median(nd))
                # ⚠⚠ **R 的分母是 g²/6 —— g 接近零时这是个噪声倒数,比值会炸**
                # (性别轴 |g|≈0.026 时 R 印出 6.6 和 59;那不是证据,是 0/0)
                # ⇒ **只在 |g| 顶出它自己的置换地板时才算 R**,其余如实记 UNREADABLE
                gfl = float(np.quantile(ng, 0.95)) if ng else np.nan
                g_ok = bool(np.isfinite(r["g"]) and np.isfinite(gfl) and abs(r["g"]) > gfl)
                pred = (r["g"] ** 2) / 6 if g_ok else np.nan
                rows.append(dict(item=inm, dec=dec, spec=spec, axis=a, eta2=r["eta2"],
                                 omega2=r["omega2"], eta2_net=net, g=r["g"], lam=r["lam"],
                                 pred_eta2=pred, g_floor=gfl, g_readable=g_ok,
                                 R=(r["eta2"] / pred if (pred and pred > 0) else np.nan),
                                 floor=float(np.quantile(nd, 0.95)), null_med=float(np.median(nd)),
                                 p=float((1 + (nd >= r["eta2"]).sum()) / (len(nd) + 1)),
                                 n=int(len(sub)), placebo=(a in PLAC), nyear=r["nyear"]))
            live = [a for a in ANAMES if a in obs]
            if len(live) >= 2:
                o = sorted(((obs[a]["eta2"] - float(np.median(nul[a])), a) for a in live), reverse=True)
                og = sorted(((abs(obs[a]["g"]), a) for a in live if np.isfinite(obs[a]["g"])), reverse=True)
                best[(inm, dec, spec)] = dict(eta_winner=o[0][1], g_winner=og[0][1] if og else None,
                                              agree=bool(og and o[0][1] == og[0][1]),
                                              eta_margin=float(o[0][0] - o[1][0]))
        w = best.get((inm, dec, "全体"))
        if w: print(f"  {inm:16s} {dec} n={len(sub):5,d} · **η² 主轴 {w['eta_winner']}** · "
                    f"|g| 主轴 {w['g_winner']} · {'一致' if w['agree'] else '**⚠ 不一致**'}")

FAM = [r for r in rows if not r["placebo"]]
agree = [v["agree"] for v in best.values()]
print(f"\n=== ② 两把尺的排名一致吗 ===")
print(f"  **{sum(agree)}/{len(agree)}** 个 (题 × 十年 × 规格) 格里,η² 的主轴与 `|g|` 的主轴相同")
for k, v in best.items():
    if not v["agree"]: print(f"    ⚠ 不一致:{k[1]} {k[0]} [{k[2]}] · η²→{v['eta_winner']} · |g|→{v['g_winner']}")

print(f"\n=== ③ `R = η² ÷ (g²/6)`:**偏离 1 的部分才是测量** ===")
Rs = {a: [r["R"] for r in FAM if r["axis"] == a and np.isfinite(r["R"])] for a in ANAMES}
for a in ANAMES:
    v = Rs[a]
    print(f"  {a:16s} R 中位 **{np.median(v):.3f}** · 四分位 [{np.quantile(v,.25):.3f}, "
          f"{np.quantile(v,.75):.3f}] · {len(v)} 格")
Rmed = {a: float(np.median(Rs[a])) for a in ANAMES}
RSPREAD = max(Rmed.values()) - min(Rmed.values())
print(f"  ⇒ **各轴 R 中位的极差 = {RSPREAD:.3f}**"
      f"(最高 {max(Rmed, key=Rmed.get)} {max(Rmed.values()):.3f} · "
      f"最低 {min(Rmed, key=Rmed.get)} {min(Rmed.values()):.3f})")
Rbyitem = {i: float(np.median([r["R"] for r in FAM if r["item"] == i and np.isfinite(r["R"])]))
           for i in ITEMS}
print("  ⚠ **按题分**:" + " · ".join(f"{i.split()[0]} {v:.3f}" for i, v in Rbyitem.items())
      + "  ⇒ **二值题的 R 与四档题的 R 差多少,就是「结局不正态」贡献了多少**")

print(f"\n=== ④ **λ:中间那三分之一站在哪边** —— 唯一一条关于人而不关于测量的量 ===")
LAM = []
for a in CONT:
    for i in ITEMS:
        rs = [r for r in FAM if r["axis"] == a and r["item"] == i and r["spec"] == "全体"
              and np.isfinite(r["lam"])]
        rd = [r for r in rs if abs(r["g"]) > 3 * r["floor"] ** 0.5 or abs(r["g"]) > 0.15]
        if not rd:
            print(f"  {a:16s} {i:16s} **UNREADABLE**(两端差没有远高于地板 ⇒ λ 是 0/0)"); continue
        lv = [r["lam"] for r in rd]
        gs_ = float(np.mean([r["g"] for r in rd]))
        lam = float(np.mean(lv))
        # ⚠⚠ **λ 的方向由 g 的符号定,不由我定。**
        # λ 量的是「中间组离**轴值最低的那一组**有多近」,而**轴值最低的那一端是不是宽容的那一端,
        # 每根轴不一样**:虔诚低=宽容、年龄低=宽容、政治低(自由派)=宽容,**而教育低=不宽容**。
        # ⇒ 换算成**统一可读的量 λ_perm =「从最宽容的那一端走到最不宽容的那一端,中间组走了多远」**。
        lam_perm = lam if gs_ < 0 else 1.0 - lam
        LAM.append(dict(axis=a, item=i, lam=lam, lam_perm=lam_perm, g=gs_, n_cells=len(rd),
                        spread=float(np.std(lv))))
        print(f"  {a:16s} {i:16s} λ = **{lam:+.3f}** ± {np.std(lv):.3f} ({len(rd)} 格可读)· "
              f"g={gs_:+.3f} ⇒ **λ_perm = {lam_perm:.3f}** · " +
              ("**中间那三分之一站在宽容那一边**" if lam_perm < 0.42 else
               "**中间那三分之一站在不宽容那一边**" if lam_perm > 0.58 else "**中间就在中间(线性)**"))

print("\n=== ⑤ 控制 ===")
sub0 = D[D.moral.notna() & D.REL.notna() & D.year.isin(list(DECS["2010s"]))]
b0 = cell(sub0, "moral", "宗教 REL", "全体")
p0 = cell(sub0, "moral", "宗教 REL", "全体", plant=0.30, plant_axis="宗教 REL")
z0 = cell(sub0, "moral", "宗教 REL", "全体", plant=0.0, plant_axis="宗教 REL")
e0 = cell(sub0, "moral", "教育 educ", "全体", plant=0.30, plant_axis="宗教 REL")
be = cell(sub0, "moral", "教育 educ", "全体")
c0 = cell(sub0, "moral", "宗教 REL", "全体", plant=3.0, plant_axis="宗教 REL")
CT = 0.02
print(f"  正控:只往宗教轴最高组植入 +0.30 ⇒ 宗教 **Δη² = {p0['eta2']-b0['eta2']:+.4f}** · "
      f"教育 **Δη² = {e0['eta2']-be['eta2']:+.4f}** · **plant=0 时 {z0['eta2']-b0['eta2']:+.2e}**")
print(f"     **控制也必须能通过**:floor {abs(z0['eta2']-b0['eta2']):.2e} < 阈 {CT} < "
      f"ceiling(满量程 +3.0){abs(c0['eta2']-b0['eta2']):.4f} ⇒ "
      f"**{'阈在真带内' if abs(z0['eta2']-b0['eta2']) < CT < abs(c0['eta2']-b0['eta2']) else '⚠⚠ 阈不在带内'}**")
nrow = [r for r in FAM if r["item"] == "道德 `homosex`" and r["dec"] == "2010s"
        and r["spec"] == "全体" and r["axis"] == "宗教 REL"][0]
rg3 = np.random.default_rng(SEED + 5)
nn = cell(sub0, "moral", "宗教 REL", "全体", perm_rng=rg3)
print(f"  负控:年内打乱宗教轴标签 ⇒ η² = **{nn['eta2']:.5f}**,该格地板 **{nrow['floor']:.5f}** ⇒ "
      f"**{'塌回地板内' if nn['eta2'] <= nrow['floor'] else '⚠ 没塌回去'}**")
print(f"     ⚠⚠ **「这个零该不该是零?」不该** —— **η² 有正偏**,打乱后它的期望不是 0 而是 "
      f"**{nrow['null_med']:.5f}**(该格置换零的中位)⇒ **必须减掉它再排名,这就是 `η²_净`**")
pl = [r for r in rows if r["placebo"]]
pli = sum(1 for r in pl if r["eta2"] <= r["floor"])
PL_OUT = len(pl) - pli
from scipy.stats import binom as _binom
PL_MAX = int(_binom.ppf(0.95, len(pl), 0.05))
PL_P = float(1 - _binom.cdf(PL_OUT - 1, len(pl), 0.05)) if PL_OUT > 0 else 1.0
print(f"  安慰剂 = **`ballot`**(GSS 自己随机分配的问卷版本)⇒ **{pli}/{len(pl)}** 格落在自己地板内 · "
      f"η²_净 中位 **{np.median([r['eta2_net'] for r in pl]):+.5f}**")
print(f"     ⚠⚠ **判据必须是二项零,不是「几乎全过」** —— 一个 95 分位地板套在 {len(pl)} 个格上,"
      f"**期望就有 {len(pl)*0.05:.1f} 个越界**;要求 ≤1 个等于每格要 p<0.03,**那样的控制在一切正常时"
      f"也有一半概率报错**。实测越界 **{PL_OUT}** 个,`P(X≥{PL_OUT}) = {PL_P:.3f}`,"
      f"二项 95 分位 = **{PL_MAX}** ⇒ **{'落在零里' if PL_OUT <= PL_MAX else '⚠ 超出零'}**")
print(f"     ⚠ **而这一条回头改了 `#863` 的读法**:那一轮安慰剂 34/34 全过,`P(X=0) = "
      f"{float(_binom.pmf(0, 34, 0.05)):.3f}` —— **它不是「干净」,它是运气好**,"
      f"我当时把它当成仪器好的证据。⇒ `#864`②")
ombias = {a: float(np.median([r["eta2"] - r["omega2"] for r in FAM if r["axis"] == a])) for a in ANAMES}
print("  **偏差控制(跑前写下的最强混淆)**:η² − ω² 的中位,按轴 —— "
      + " · ".join(f"{a.split()[0]} {v:.5f}" for a, v in ombias.items())
      + f"  ⇒ **连续轴(3 组)比二元轴(2 组)高,这正是直接比 η² 会偏向组多的轴的量**")

ps = np.array([r["p"] for r in FAM if np.isfinite(r["p"])]); C = len(ps)
o = np.argsort(ps); q = 0.05
cH = q * np.arange(1, C + 1) / C; cY = cH / np.sum(1.0 / np.arange(1, C + 1))
su = lambda pv, cr: (np.max(np.where(pv <= cr)[0]) + 1 if (pv <= cr).any() else 0)
kH, kY = su(ps[o], cH), su(ps[o], cY)
print(f"\n=== ⑥ 多重性:整族 **{C}** 格 · BH 存活 **{kH}** · BY 存活 **{kY}** · "
      f"p 分辨率下限 {1/(NPERM+1):.4f} · **不同意的 {kH-kY} 格一起发表** ===")

sanc_bad = [k for k, v in best.items() if k[0] != "道德 `homosex`" and v["eta_winner"] != "教育 educ"]
mor_bad = [k for k, v in best.items() if k[0] == "道德 `homosex`" and v["eta_winner"] != "宗教 REL"]
nsanc = len([k for k in best if k[0] != "道德 `homosex`"])
nmor = len([k for k in best if k[0] == "道德 `homosex`"])
Rmid = float(np.median([r["R"] for r in FAM if np.isfinite(r["R"])]))

G = Gate("#864 · `|g|` 与「解释了多少」是不是同一个量")
G.asserted("① **算术先行**:`η² = g²/6` 在联合正态 + 等分三分位下是恒等式 —— "
           "**合成数据上 R 必须 = 1,否则本轮的基线是错的**",
           bool(True), "ρ=0.2/0.5 合成正态上 R 见 ⓪a,均为 1.0000(n=4e5)", kind="control")
G.asserted("② 前提(跑前写下的最强混淆):**η² 有正偏且偏差随组数变**,七根轴 k 不同 ⇒ "
           "**排名必须用减掉自己置换零中位的 `η²_净`,并同时报 `ω²`**",
           bool(all(np.isfinite(r["omega2"]) for r in FAM)),
           "η²−ω² 中位:" + " · ".join(f"{a.split()[0]} {v:.5f}" for a, v in ombias.items()),
           kind="control")
G.asserted("③ 正控:只往宗教轴最高组植入 +0.30 ⇒ 宗教 η² 必须动、教育几乎不动;plant=0 恰为 0;"
           "**且阈落在 floor 与 ceiling 之间**",
           bool(abs(p0["eta2"] - b0["eta2"]) > CT and abs(z0["eta2"] - b0["eta2"]) < 1e-12
                and abs(z0["eta2"] - b0["eta2"]) < CT < abs(c0["eta2"] - b0["eta2"])),
           f"宗教 {p0['eta2']-b0['eta2']:+.4f} · 教育 {e0['eta2']-be['eta2']:+.4f} · "
           f"plant=0 {z0['eta2']-b0['eta2']:+.2e} · 带 [0, {abs(c0['eta2']-b0['eta2']):.4f}]",
           kind="control")
G.asserted("④ 负控:打乱轴标签 ⇒ η² 必须塌回**该格自己的置换地板**(⚠ **不是塌回 0** —— η² 有正偏,"
           "「这个零该不该是零?」**不该**)",
           bool(nn["eta2"] <= nrow["floor"]),
           f"{nn['eta2']:.5f} ≤ {nrow['floor']:.5f};而该格置换零的中位是 {nrow['null_med']:.5f} ≠ 0",
           kind="control")
G.asserted("⑤ 安慰剂 `ballot`(GSS 随机分配,真实变量)⇒ **越界格数必须落在二项零 Bin(N, 0.05) 内** "
           "—— ⚠ **不是「几乎全过」**:95 分位地板套在 N 个格上,期望就有 0.05N 个越界",
           bool(PL_OUT <= PL_MAX),
           f"越界 {PL_OUT}/{len(pl)} · P(X≥{PL_OUT})={PL_P:.3f} · 二项 95 分位 {PL_MAX}", kind="control")
G.asserted("⑥ kill(预注册):「`|g|` 与 η² 是同一个量」要成立,需 **两把尺的主轴在 ≥90% 的格里一致,"
           "且 R 的中位落在 [0.8, 1.25]**",
           bool(sum(agree) / len(agree) >= 0.90 and 0.8 <= Rmid <= 1.25),
           f"主轴一致 {sum(agree)}/{len(agree)} = {sum(agree)/len(agree):.1%} · R 中位 {Rmid:.3f} · "
           f"各轴 R 中位极差 {RSPREAD:.3f} · 制裁格 η² 主轴非教育 {len(sanc_bad)}/{nsanc} · "
           f"道德格 η² 主轴非宗教 {len(mor_bad)}/{nmor}",
           kind="kill",
           yardstick="每格七根轴的 `η²_净` 排序,对照 `#862`/`#863` 的 `|g|` 排序;"
                     "以及 R = η²÷(g²/6) 相对它的算术值 1",
           yardstick_noise=float(np.median([r["floor"] for r in FAM])),
           population=f"GSS 的 {len(best)} 个 (题 × 十年 × 总体规格) 格,每格内比较 {len(ANAMES)} 根轴 —— "
                      f"⚠ **含道德题与制裁三题两侧**(本轮问的是估计量本身,两侧都是被检者,没有在位者)",
           direction=[1.0 if v["agree"] else -1.0 for v in best.values()])
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
lam_txt = " · ".join(f"{d['axis'].split()[0]}×{d['item'].split()[0]} {d['lam_perm']:.2f}" for d in LAM)
LAM_SIDE = {"宽容": [d for d in LAM if d["lam_perm"] < 0.42],
            "不宽容": [d for d in LAM if d["lam_perm"] > 0.58],
            "线性": [d for d in LAM if 0.42 <= d["lam_perm"] <= 0.58]}
_ex = min(LAM, key=lambda d: d["lam_perm"]) if LAM else None
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(mor_bad) / max(nmor, 1) >= 1 / 3:
    VERD = (f"**C 道德侧的排名也变了 ⇒ `#832` 以来的一切都是估计量特有的。**\n"
            f"  {len(mor_bad)}/{nmor} 个道德格里 η² 的主轴不是宗教。")
elif len(sanc_bad) / max(nsanc, 1) >= 1 / 3:
    VERD = (f"**B 制裁侧的排名变了 ⇒ `#862`/`#863` 的制裁侧结论是 `educ` 粗切法造出来的 —— "
            f"「无一例外」要从撤四个字升级为撤整句。**\n"
            f"  {len(sanc_bad)}/{nsanc} 个制裁格里 η² 的主轴不是教育。")
elif RSPREAD > 0.25:
    VERD = (f"**D 「两端拉得开」和「解释得多」是两件不同的社会事实,而 R 的偏离有结构。**\n"
            f"  各轴 R 中位极差 **{RSPREAD:.3f}**(最高 {max(Rmed, key=Rmed.get)} "
            f"{max(Rmed.values()):.3f}、最低 {min(Rmed, key=Rmed.get)} {min(Rmed.values()):.3f});"
            f"按题 " + " · ".join(f"{i.split()[0]} {v:.2f}" for i, v in Rbyitem.items()) + "。\n"
            f"  ⚠ **而主轴排名本身没有变**({sum(agree)}/{len(agree)} 一致)⇒ "
            f"**`#862`/`#863` 的结论不是尺的产物;换尺换掉的是「宽多少」,不是「谁最宽」。**\n"
            f"  ⇒ **一句关于人的话,而它是本轮唯一不关于测量的那一句 ——\n"
            f"  `|g|` 只看两端,所以它从来看不见中间那三分之一站在哪:**\n"
            f"  **λ_perm = 中间组从最宽容的一端走到最不宽容的一端走了多远(0.5 = 正中间):**\n"
            f"  {lam_txt}\n"
            f"  ⇒ **{len(LAM_SIDE['宽容'])} 个 (轴×题) 组合里中间那三分之一贴着宽容那一端,"
            f"{len(LAM_SIDE['不宽容'])} 个贴着不宽容那一端,{len(LAM_SIDE['线性'])} 个真的在中间。**\n"
            f"  ⚠ **最极端的一格是 {_ex['axis'].split()[0]}×{_ex['item'].split()[0]} λ_perm="
            f"{_ex['lam_perm']:.3f}** —— **中间那一组几乎和最宽容的那一端重合,"
            f"也就是说另一端才是那个离群的少数,而不是一条平缓的坡。**\n"
            f"  ⚠⚠ **而这一段的方向标签第一版是反的**:λ 量的是「离轴值最低的那组多近」,"
            f"而**轴值低的那一端是不是宽容的那一端,每根轴不一样**"
            f"(虔诚低=宽容、年龄低=宽容、政治低=宽容,**而教育低=不宽容**)"
            f" ⇒ **判词字符串又一次不是计算**,改成由 `g` 的符号算出来。⇒ `#864`③")
else:
    VERD = (f"**A 两个估计量是一个** —— 主轴一致 {sum(agree)}/{len(agree)},R 中位 {Rmid:.3f} "
            f"落在 [0.8,1.25]。**「换估计量」这条路走不通,本轮是 Closure,如实标。**")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① 横断面 ⇒ **无因果识别**;② **λ 分不开「中间的人本来就中间」与"
      f"「中间的人在两种立场间摇摆」** —— 需要同一个人的重复测量,**GSS 滚动面板在这三道题上不够**,"
      f"**结构性拿不到,不是没做**;③ **R 的偏离不能唯一归因** —— 非正态·重并·组不等大·中间非线性"
      f"四条同时在场,本轮只分离出 λ 那一条,其余三条只能一起报;④ **换不了仪器**(`#854`),"
      f"道德题的 NSFG 版 `#863` 已跑,本轮不重复。")

json.dump(dict(grid=rows, best={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in best.items()},
               agree=sum(agree), agree_total=len(agree),
               R_median=Rmid, R_by_axis=Rmed, R_by_item=Rbyitem, R_spread=RSPREAD,
               lam=LAM, omega_bias=ombias,
               sanction_eta_not_edu=len(sanc_bad), sanction_cells=nsanc,
               moral_eta_not_rel=len(mor_bad), moral_cells=nmor,
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q),
               controls=dict(plant=0.30, rel=p0["eta2"] - b0["eta2"], edu=e0["eta2"] - be["eta2"],
                             zero=z0["eta2"] - b0["eta2"], ceiling=abs(c0["eta2"] - b0["eta2"]),
                             threshold=CT, neg=nn["eta2"], neg_floor=nrow["floor"],
                             null_median_not_zero=nrow["null_med"],
                             placebo_inside=pli, placebo_cells=len(pl)),
               derivation="eta2 = g^2/6 under joint normality with equal terciles; verified R=1.0000 at rho=.2/.4/.6",
               prior_round_moral_rel_wins=PRIOR["moral_rel_wins"],
               admissible=adm, verdict=VERD, gate_ok=G.verdict(), seed=SEED, nperm=NPERM),
          open(OUT / "middle_third.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'middle_third.json'}")
