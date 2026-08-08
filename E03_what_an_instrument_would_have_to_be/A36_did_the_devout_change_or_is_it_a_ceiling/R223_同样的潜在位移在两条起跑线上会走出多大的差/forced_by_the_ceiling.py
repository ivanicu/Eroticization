"""#784 · E03·A44·R223 —— 同样的潜在位移,放在两条不同的起跑线上,会自己走出多大的差?

页面上现在最硬的一句是 `#783` 刚写下的:**「最虔诚的三分之一人在同性恋上的改变,
只有其余人的四分之一到接近六成」**([0.241, 0.568],四个合并规格全部排除 1.0)。
**而它最通俗的反对解释从来没有被造出来过。**

⚠⚠ **反对解释:虔诚层 1988 年的端点占比是 85.9%,非虔诚层是 61.1%。**
   一把有界的尺子,**起点越靠顶,同样的态度变化能挪出的格子越少。**
   ⇒ **「虔诚者改得少」也许根本不是心理学,是算术** —— 而 `#776` 那条有界尺算术、
   `#777` 那个「到顶前的窗口不存在」,都已经在门口了,**却从没有人把它做成一个零。**

G1 估计量(方法之前先命名):**`r_forced`** ——
   **给两层施加完全相同的潜在位移 Δ,由各自的起跑线自己走出来的那个比值。**
   它是「零心理学差异」这个世界的**预言值**,不是一个零假设的 p。

识别:给定 ① 两层在首年的四档分布 ② 一个潜变量模型(阈值共用、组均值各异),
   `r_forced` 是**可识别的**。**阈值共用是这里的实质假设**,写在前面:
   若两层对同一个问句用词的理解不同(differential item functioning),它就不成立 ——
   **本轮不能检验它,登记为能力边界**(需要同一批人在两种身份下作答,或一份 DIF 校准样本)。

⚠⚠ **三个世界,而第三个我会很不想要:**
   A **心理学**:`r_obs` 明显小于 `r_forced` ⇒ 虔诚者改得少,**超出天花板逼出来的那部分**。
   B **算术**:`r_obs ≈ r_forced` ⇒ **整条 A43 弧的头条是一把尺子的形状**,页面必须撤。
   C **反向**:`r_obs > r_forced` ⇒ **扣掉天花板之后,虔诚者反而改得比"同样的位移"更多** ——
     那会把页面那句话整个翻过来。

预测矩阵:
   | 世界 | 现在 | 若 r_obs 区间不含 r_forced 且在其下 | 若含 | 若在其上 |
   | A | 0.45 | **0.90** | 0.10 | 0.05 |
   | B | 0.40 | 0.05 | **0.85** | 0.10 |
   | C | 0.15 | 0.05 | 0.05 | **0.85** |
   三支都强,而**最坏的那一支(B)恰好是最便宜就能查出来的** —— 这正是该先跑它的理由。

预注册判词(条件式,不是阈值 —— `#P16`):
  if 正控开火(两层起跑线相同时 r_forced ≈ 1.0,容差 0.05)
     and 负控开火(起跑线越极端,r_forced 越小 —— 单调,否则模拟器本身没在响应起跑线):
      if `r_obs` 的自助区间**含** r_forced -> **B**:与有界尺算术不可分,页面那句话必须撤/加限定
      elif r_obs 上界 < r_forced           -> **A**:超出天花板的那部分是真的,报 r_obs/r_forced
      else                                 -> **C**:反号,页面那句话要翻过来
  else: UNVERIFIED
⚠ **Δ 的标定不是自由参数**:令非虔诚层走出**它自己观测到的** 28 年变化,Δ 由此反解 ——
  所以 `r_forced` 里没有一个我挑的数。

⚠ 跑之前写下的最强混淆:**潜变量是不是正态。** 若真实潜分布偏斜,阈值被错估,`r_forced` 会动。
  ⇒ 同一轮里放对照:**probit 与 logit 两种连接函数一起跑**(尾部形状不同),
  若两者之差大于效应本身,**如实说不可读**,不选边。

规格网格(`G4`,不是一格):连接 {probit, logit} × 统计量 {水平均值, 端点占比} ×
  分层 {(a) attend 三档, (b) 三题合成三分位} = **8 格,全部报,包括杀掉结论的那些**。

本轮结构上做不到的(`realstat §2`,各附「需要什么」):
  · **阈值共用的检验(DIF)** —— 需要同一批人在两种虔诚度下作答,或一份外部 DIF 校准样本。**不可能。**
  · **因果方向** —— GSS 是重复截面,需要面板。**结构性不可能。**
  · **潜分布的真实形状** —— 需要一把连续尺;四档序数量表**原理上定不了**,所以只能报两种连接的跨度。
"""
import numpy as np, pandas as pd, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate

RNG = np.random.default_rng(223)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"

VALID = {"homosex": (1, 4), "attend": (0, 8), "reliten": (1, 4), "fund": (1, 3)}
for c, rg in VALID.items():
    dr, _ = check_kept_codes(gp, c, rg)
    if dr: print(f"  #766 前瞻:{c} 删 " + " · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a, b, n, sh in dr[:2]))
d = pd.read_stata(gp, columns=["year"]+list(VALID), convert_categoricals=False)
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, lo=VALID[c][0], hi=VALID[c][1]: (v >= lo) & (v <= hi)) for c in VALID})
M["year"] = d.year
cat = pd.read_stata(gp, columns=["homosex"], convert_categoricals=True)
for c in aligned({"homosex": list(cat["homosex"].cat.categories)[:4]}, "strict"): M[c] = -M[c]+5
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
z = lambda s: (s-s.mean())/s.std(ddof=1)
sub = M.dropna(subset=list(VALID)+["year"]).copy()
sub["REL"] = z(sub[["attend", "reliten", "fund"]]).mean(axis=1)
sub["a"] = pd.cut(sub.attend, [-1, 1, 5, 8], labels=[0, 1, 2]).astype(float)
sub["b"] = sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def dist(col, k, year):
    g = sub[(sub[col] == k) & (sub.year == year)]
    return np.array([(g.homosex == c).mean() for c in (1, 2, 3, 4)]), len(g)

def years(col, k, nmin=120):
    return sorted(int(y) for y, gy in sub[sub[col] == k].groupby("year") if len(gy) >= nmin)

# ── 潜变量机器:阈值由**合并首年**分布定(两层共用),组均值由各自首年分布定 ──────────
def fit_thresholds(p, link):
    """由四档占比反解三个共用阈值(潜均值固定为 0)。"""
    F = link.ppf
    cum = np.clip(np.cumsum(p)[:3], 1e-6, 1-1e-6)
    return F(cum)

def fit_mu(p, tau, link):
    """给定共用阈值,反解该层的潜均值(一维,单调 ⇒ 二分)。"""
    def cats(mu):
        e = link.cdf(tau-mu)
        return np.diff(np.concatenate(([0.0], e, [1.0])))
    lo, hi = -6.0, 6.0
    for _ in range(60):   # 60 次二分在宽 12 的区间上已到 float 精度;200 次是纯浪费
        mid = (lo+hi)/2
        # 均值单调递增于 mu ⇒ 用 1..4 的期望值做单调目标
        if (cats(mid)*np.arange(1, 5)).sum() < (p*np.arange(1, 5)).sum(): lo = mid
        else: hi = mid
    return (lo+hi)/2

def forced_ratio(pA0, pB0, nA0, nB0, link, stat, dB_obs):
    """**天花板逼出来的比** —— 两层共用阈值、各自起跑线,施加同一个潜在位移 Δ。

    ⚠ 这个函数是主网格与**两条对照**共用的唯一路径。第一版把正控写成 `d0/d0`,
       那是**拿一个量和它自己比**(本项目已四次的那一族,而 `vacuous_check_lint`
       当时看不见 `X / X` 这个形状 —— 本轮顺手把它补进了 lint)。
       现在正控是「令 pA0 = pB0 再调用本函数」,**它会因为函数错而失败。**
    """
    pooled = (pA0*nA0 + pB0*nB0)/(nA0+nB0)
    tau = fit_thresholds(pooled, link)
    muA, muB = fit_mu(pA0, tau, link), fit_mu(pB0, tau, link)
    obsA, obsB = readout(muA, tau, link, stat), readout(muB, tau, link, stat)
    # Δ 由非虔诚层**自己观测到的**变化标定。
    # ⚠⚠ 第一版把搜索区间写成 `[0, 6]` —— **只找向上的位移**,而 1974→2024 的态度是变宽容的,
    #    量表上「高 = 严」⇒ `dB_obs < 0` ⇒ 每一次二分都不成立,Δ 塌成 0,`r_forced` 全成 -0.000。
    #    **两条对照当场同时开火(正控不返回 1.0、负控不单调),判 UNVERIFIED,正确。**
    #    这正是「控制因自己的理由而失败」的反面:控制失败,而失败的是被控的东西。
    lo, hi = -6.0, 6.0
    for _ in range(60):   # 60 次二分在宽 12 的区间上已到 float 精度;200 次是纯浪费
        mid = (lo+hi)/2
        if readout(muB+mid, tau, link, stat)-obsB < dB_obs: lo = mid
        else: hi = mid
    delta = (lo+hi)/2
    return (readout(muA+delta, tau, link, stat)-obsA)/dB_obs, delta, obsA, obsB


def readout(mu, tau, link, stat):
    e = link.cdf(tau-mu)
    p = np.diff(np.concatenate(([0.0], e, [1.0])))
    return float((p*np.arange(1, 5)).sum()) if stat == "水平" else float(p[3])

print("\n=== ① 两层的起跑线(首年四档分布)===")
CUTS = {"(a) attend三档": "a", "(b) 三题三分位": "b"}
base = {}
for cn, col in CUTS.items():
    yA, yB = years(col, 2), years(col, 0)
    y0 = max(yA[0], yB[0]); y1 = min(yA[-1], yB[-1])
    pA0, nA0 = dist(col, 2, y0); pB0, nB0 = dist(col, 0, y0)
    pA1, _ = dist(col, 2, y1);   pB1, _ = dist(col, 0, y1)
    base[cn] = dict(col=col, y0=y0, y1=y1, pA0=pA0, pB0=pB0, pA1=pA1, pB1=pB1, nA0=nA0, nB0=nB0)
    print(f"  {cn}  {y0}→{y1} · 虔诚层首年端点 {pA0[3]:.3f}(n={nA0}) · 非虔诚层首年端点 {pB0[3]:.3f}(n={nB0})")

print("\n=== ② 规格网格 8 格:观测比值 vs 天花板逼出来的比值 ===")
LINKS = {"probit": norm, "logit": logistic}
cells = []
for cn, B in base.items():
    for ln, link in LINKS.items():
        pooled0 = (B["pA0"]*B["nA0"] + B["pB0"]*B["nB0"])/(B["nA0"]+B["nB0"])
        tau = fit_thresholds(pooled0, link)
        muA, muB = fit_mu(B["pA0"], tau, link), fit_mu(B["pB0"], tau, link)
        for st in ("水平", "端点"):
            realA1 = float((B["pA1"]*np.arange(1, 5)).sum()) if st == "水平" else float(B["pA1"][3])
            realB1 = float((B["pB1"]*np.arange(1, 5)).sum()) if st == "水平" else float(B["pB1"][3])
            obsA0 = readout(muA, tau, link, st); obsB0 = readout(muB, tau, link, st)
            dA_obs, dB_obs = realA1-obsA0, realB1-obsB0
            if abs(dB_obs) < 1e-9: continue
            r_obs = dA_obs/dB_obs
            r_forced, delta, _, _ = forced_ratio(B["pA0"], B["pB0"], B["nA0"], B["nB0"], link, st, dB_obs)
            dA_f = r_forced*dB_obs
            cells.append(dict(cut=cn, link=ln, stat=st, r_obs=r_obs, r_forced=r_forced,
                              delta=delta, dA_obs=dA_obs, dB_obs=dB_obs, dA_forced=dA_f))
            print(f"  {cn:14s} {ln:7s} {st:4s}  观测比 **{r_obs:.3f}** · 天花板逼出的比 **{r_forced:.3f}**"
                  f"  (Δ={delta:.3f} 潜标准差)")

R_OBS = np.array([c["r_obs"] for c in cells]); R_F = np.array([c["r_forced"] for c in cells])
print(f"\n  观测比中位 **{np.median(R_OBS):.3f}** · 天花板逼出的比中位 **{np.median(R_F):.3f}**")
print(f"  两种连接函数之差(同格 probit vs logit):中位 "
      f"{np.median([abs(cells[i]['r_forced']-cells[i+2]['r_forced']) for i in (0,1,4,5)]):.4f}")

# ── 观测比的自助区间(年份层面,与 `#782`/`#783` 同一具仪器)───────────────────────
print("\n=== ③ 观测比的年份自助区间(与 `#783` 同一具)===")
def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def rows_of(col, k, nmin=120):
    out = []
    for y, gy in sub[sub[col] == k].groupby("year"):
        if len(gy) < nmin: continue
        out.append((int(y), float(gy.homosex.mean()), float((gy.homosex == 4).mean())))
    return out
BOOT = {}
for cn, col in CUTS.items():
    rA, rB = rows_of(col, 2), rows_of(col, 0)
    for st, j in (("水平", 1), ("端点", 2)):
        yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
        vA = np.array([r[j] for r in rA]); vB = np.array([r[j] for r in rB])
        f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
        bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(4000)])
        bs = bs[np.isfinite(bs)]
        BOOT[(cn, st)] = (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)),
                          f(np.arange(len(yA)), np.arange(len(yB))))
        print(f"  {cn:14s} {st:4s} 斜率比 {BOOT[(cn,st)][2]:.3f} · 95% [{BOOT[(cn,st)][0]:.3f}, {BOOT[(cn,st)][1]:.3f}]")

print("\n=== ④ 逐格判:观测区间与「天花板逼出来的比」的关系 ===")
inside = above = below = 0
for c in cells:
    lo, hi, pt = BOOT[(c["cut"], c["stat"])]
    rf = c["r_forced"]
    where = "**含 ⇒ 与算术不可分**" if lo <= rf <= hi else ("低于区间 ⇒ 比天花板逼出的还少" if rf > hi else "高于区间 ⇒ 反号")
    if lo <= rf <= hi: inside += 1
    elif rf > hi: below += 1
    else: above += 1
    print(f"  {c['cut']:14s} {c['link']:7s} {c['stat']:4s} 区间 [{lo:.3f}, {hi:.3f}] · 逼出的比 {rf:.3f} —— {where}")
print(f"\n  **8 格:与算术不可分 {inside} · 观测明显更少 {below} · 反号 {above}**")

# ── ⑤ 扣掉天花板之后还剩多少 —— 而它自己也必须带区间(`#783` 本轮刚立的规矩)────────
print("\n=== ⑤ 扣掉天花板之后剩下的那一份:r_obs / r_forced,带它自己的区间 ===")
# ⚠ `r_forced` 不是常数:它是两层**首年四档分布**的函数,而那两个分布有抽样误差。
#    ⇒ 对首年做多项式重抽(nA0/nB0 个人),对观测比沿用 `#783` 的年份自助,**两个来源都进区间。**
resid = []
for cn, B in base.items():
    for ln, link in LINKS.items():
        for st in ("水平", "端点"):
            c = next(x for x in cells if x["cut"]==cn and x["link"]==ln and x["stat"]==st)
            lo_o, hi_o, pt_o = BOOT[(cn, st)]
            draws = []
            for _ in range(500):
                pA = RNG.multinomial(B["nA0"], B["pA0"])/B["nA0"]
                pB = RNG.multinomial(B["nB0"], B["pB0"])/B["nB0"]
                try: rf, _, _, _ = forced_ratio(pA, pB, B["nA0"], B["nB0"], link, st, c["dB_obs"])
                except Exception: continue
                ro = pt_o*np.exp(RNG.normal(0, (np.log(hi_o)-np.log(lo_o))/(2*1.96)))
                if np.isfinite(rf) and abs(rf) > 1e-9: draws.append(ro/rf)
            dr = np.array(draws)
            q = (float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5)))
            resid.append(dict(cut=cn, link=ln, stat=st, point=pt_o/c["r_forced"], lo=q[0], hi=q[1]))
            print(f"  {cn:14s} {ln:7s} {st:4s}  剩余比 **{pt_o/c['r_forced']:.3f}** · 95% [{q[0]:.3f}, {q[1]:.3f}]"
                  f"{'  ⚠含1.0' if q[0] <= 1.0 <= q[1] else ''}")
RES = np.array([r["point"] for r in resid])
n_cov1 = sum(1 for r in resid if r["lo"] <= 1.0 <= r["hi"])
print(f"  ⇒ **剩余比中位 {np.median(RES):.3f}** · 8 格里区间含 1.0 的:**{n_cov1}**")

# ── 闸 ────────────────────────────────────────────────────────────────────────
G = Gate("#784 · 虔诚者改得少,是不是天花板逼出来的")
eq = []
for ln, link in LINKS.items():
    B0 = base["(a) attend三档"]
    dB = next(c["dB_obs"] for c in cells if c["cut"] == "(a) attend三档" and c["link"] == ln and c["stat"] == "水平")
    r_eq, _, _, _ = forced_ratio(B0["pB0"], B0["pB0"], B0["nB0"], B0["nB0"], link, "水平", dB)
    eq.append(r_eq)
G.asserted("① 正控:两层起跑线相同时,天花板逼出的比必须 ≈ 1.0(容差 0.05)",
           bool(all(abs(e-1.0) <= 0.05 for e in eq)), f"probit/logit 各得 {[round(e,4) for e in eq]}", kind="control")
mono = []
for ln, link in LINKS.items():
    B0 = base["(a) attend三档"]
    dB = next(c["dB_obs"] for c in cells if c["cut"] == "(a) attend三档" and c["link"] == ln and c["stat"] == "端点")
    tau0 = fit_thresholds(B0["pB0"], link); mu0 = fit_mu(B0["pB0"], tau0, link)
    seq = []
    for shift in (0.0, 0.5, 1.0, 1.5):
        e = link.cdf(tau0-(mu0+shift))
        pA = np.diff(np.concatenate(([0.0], e, [1.0])))     # 人为把 A 层起跑线往顶上推
        r, _, _, _ = forced_ratio(pA, B0["pB0"], B0["nB0"], B0["nB0"], link, "端点", dB)
        seq.append(r)
    mono.append(all(seq[i] > seq[i+1] for i in range(3)))
    print(f"\n  负控 {ln}:起跑线越靠顶(+0/0.5/1.0/1.5 潜标准差),逼出的比 = {[round(s,3) for s in seq]}")
G.asserted("② 负控:起跑线越极端,逼出的比必须**单调变小**(否则模拟器没在响应起跑线)",
           bool(all(mono)), f"probit/logit 单调:{mono}", kind="control")
spread = float(np.median([abs(cells[i]["r_forced"]-cells[i+2]["r_forced"]) for i in (0, 1, 4, 5)]))
eff = float(abs(np.median(R_F)-np.median(R_OBS)))
G.asserted("③ 负控:两种连接函数之差必须小于效应本身,否则潜分布假设主导了答案",
           bool(spread < eff), f"连接跨度 {spread:.4f} vs 效应 {eff:.4f}", kind="control")
G.asserted("⑤ 正控(`#783` 的规矩):每一个要上页面的比值都必须带自己的区间,且区间非零宽",
           bool(all(r["hi"]-r["lo"] > 1e-6 for r in resid)),
           f"8 格剩余比区间宽 {min(r['hi']-r['lo'] for r in resid):.3f}–{max(r['hi']-r['lo'] for r in resid):.3f}", kind="control")
G.asserted("④ kill(预注册):页面那句话要站住,需**多数格**的观测区间把「逼出的比」排除在外且在其下",
           bool(below > inside + above), f"更少 {below} · 不可分 {inside} · 反号 {above}(共 8)", kind="kill")
print(); print(G)

print("\n" + "="*92)
ctrl = all(abs(e-1.0) <= 0.05 for e in eq) and all(mono) and spread < eff
if not ctrl:
    v = "**UNVERIFIED:正控或负控没过,本轮不下判。**"
elif below > inside + above:
    v = (f"**A:虔诚者改得少,超出了天花板能逼出来的那部分。** 8 格里 **{below} 格**的观测自助区间"
         f"把「同样的潜在位移会走出的比」({np.median(R_F):.3f})**排除在外且在其下**;"
         f"观测比中位 {np.median(R_OBS):.3f} 对逼出的 {np.median(R_F):.3f} —— "
         f"**天花板解释得了一部分,解释不完。**\n"
         f"  ⇒ **扣掉天花板之后剩下的那一份:{np.median(RES):.3f}**(8 格中位),"
         f"区间含 1.0 的 {n_cov1}/8。")
elif inside >= below:
    v = (f"**B:与有界尺算术不可分。** 8 格里 **{inside} 格**的观测区间**包含**「天花板逼出的比」"
         f"({np.median(R_F):.3f}),而观测比中位是 {np.median(R_OBS):.3f}。"
         f"⇒ **页面那句「只有其余人的四分之一到接近六成」必须加限定或撤** —— "
         f"在这套数据上,它与「两层起跑线不同、心理学上没差别」分辨不开。")
else:
    v = (f"**C:反号。** {above} 格的观测区间高于「天花板逼出的比」—— "
         f"扣掉起跑线之后虔诚者改得**更多**,页面那句话要翻过来。")
print(v)
json.dump(dict(cells=cells, residual=resid, residual_median=float(np.median(RES)), residual_covers1=n_cov1, boot={f"{k[0]}|{k[1]}": v for k, v in BOOT.items()},
               r_obs_median=float(np.median(R_OBS)), r_forced_median=float(np.median(R_F)),
               inside=inside, below=below, above=above, link_spread=spread, effect=eff,
               verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"forced_by_the_ceiling.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'forced_by_the_ceiling.json'}")
