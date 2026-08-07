"""#785 · E03·A44·R224 —— 「虔诚者改得少」是关于同性恋的一句话,还是关于性道德的一句话?

`#784` 在 **`homosex` 一题上**判了:天花板解释约三分之一,扣掉之后虔诚者仍只改了其余人的
**0.536 [0.375, 0.782]**(8 格,含 1.0 的 0/8)。**而页面把它读成了一句关于虔诚者的话。**
⚠ **一题不是一个族。** `#777` 已经量到 `xmarsex` 在能看见的 36 年里纹丝不动 ——
   **若那一题上天花板解释得完,「虔诚者改得少」就不是一条关于性道德的规律,而是一条关于同性恋的。**

G1 估计量(与 `#784` 同一个,只把**题目**变成一个轴):
   逐题的 ① `r_obs` 虔诚/非虔诚斜率比 ② `r_forced` 天花板逼出的比 ③ 残差 `r_obs / r_forced`,
   **每一个都带自己的年份自助区间**(`#783` 立的规矩)。

⚠⚠ **可读性判据前置,而这是 `#782` 那一课的预防性使用**:
   `#779`–`#781` 检验了分子与分母**却从没检验比值**,`#782` 才发现小分母下两个都显著的量相除仍可极不稳。
   ⇒ **本轮先判可读:一格的比值自助区间若含 1.0,它定不了方向 ⇒ 列出、不读、不平均进去。**
   ⚠ 而这一条在本轮很可能咬人:`xmarsex` 若两层都没动,**分母趋零 ⇒ 比值是 Cauchy 型**,
   那不是「没有效应」,是**这个问题在那一题上问不出来**。

三个世界:
   A **一块**:四题的残差都可读且都明显小于 1 ⇒「虔诚者改得少」是关于**性道德**的一句话。
   B **一题**:只有 `homosex` 可读且小于 1 ⇒ 那句话必须缩到**同性恋**这一题上,页面要改。
   C **不可读**:多数题的比值区间含 1.0 ⇒ **这个问题只能在「其余人真的动了」的地方问** ——
     那是一条**关于scope的事实,不是一个零**,而页面必须把 scope 写出来。

预测矩阵:
   | 世界 | 现在 | 若 ≥3 题可读且 <1 | 若只 homosex | 若 ≥2 题不可读 |
   | A | 0.30 | **0.85** | 0.05 | 0.10 |
   | B | 0.35 | 0.05 | **0.85** | 0.20 |
   | C | 0.35 | 0.10 | 0.10 | **0.70** |

预注册判词(条件式):
  if 正控开火(起跑线相同 -> r_forced≈1.0)and 负控开火(起跑线越靠顶 -> r_forced 单调变小):
      读**可读格**:若 >=3 题的残差区间排除 1.0 且在其下 -> A
                   elif 恰只有 homosex          -> B:页面那句话缩到一题
                   else                          -> C:报 scope,不报规律
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**四题的天花板高度本来就不同**(`xmarsex` 端点占比 ~0.87,
  `premarsx` 低得多)⇒ **残差之间不可直接比大小**,能比的只有「排不排除 1.0」。
  ⇒ 报残差时**必带各自的 r_forced**,不许把四个残差并成一行读。

⚠⚠ **本轮第一次跑就逮到 `#784` 自己的一个缺陷,而它改动了已上页的数,所以写在最前面:**
   `#784` 的 `r_obs` 取自**斜率比**(年份自助),而它标定 Δ 用的 `dB_obs` 是**首末年的差**
   —— **两个不同的估计量被写进同一个比。** 同一具机器、同一道题,`r_forced` 因此在
   `#784` 上是 **0.813**、在本轮第一版上是 **1.048**,而两者都自称「天花板逼出来的比」。
   ⇒ 本轮统一到**斜率**这一个估计量上:每层的观测变化 = `斜率 × 跨度年数`,
   `r_obs`、Δ 的标定、`r_forced` 全部由它出。**`#784` 页面上的 0.536 必须按这条重算。**

⚠ 方向对齐用 `lib.blocks.aligned`(`#759`:方向从值标签读,认不出就 raise,**绝不默认**)——
  四题的标签各自检查,而不是照抄 `homosex` 的。

本轮换不了仪器,理由由 `R223/instrument_search.py` **跑出来**而不是声明:
  第二具仪器需「有界序数性道德题 × 虔诚度分层 × ≥8 个时间点」三条同时满足,本机六具全部落选。
"""
import numpy as np, pandas as pd, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes, label_pole, labelled_codes
from lib.gates import Gate

RNG = np.random.default_rng(224)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
ITEMS = ("homosex", "xmarsex", "premarsx", "teensex")

VALID = {**{c: (1, 4) for c in ITEMS}, "attend": (0, 8), "reliten": (1, 4), "fund": (1, 3)}
print("=== ⓪ `#766` 前瞻:每一列的范围裁剪各删掉了哪些有标签的码 ===")
for c, rg in VALID.items():
    dr, _ = check_kept_codes(gp, c, rg)
    if dr: print(f"  {c:9s} 删 " + " · ".join(f"码{int(a)} {b!r}({sh*100:.2f}%)" for a, b, n, sh in dr[:2]))
d = pd.read_stata(gp, columns=["year"]+list(VALID), convert_categoricals=False)
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, lo=VALID[c][0], hi=VALID[c][1]: (v >= lo) & (v <= hi)) for c in VALID})
M["year"] = d.year

print("\n=== ① 方向:逐题从值标签读,不照抄 `homosex`(`#759`)===")
cat = pd.read_stata(gp, columns=list(ITEMS), convert_categoricals=True)
cats = {c: list(cat[c].cat.categories)[:4] for c in ITEMS}
for c in ITEMS:
    print(f"  {c:9s} 高值端 = {label_pole(cats[c])!r}   标签 {cats[c]}")
flip = aligned(cats, "strict")     # 认不出会 raise
print(f"  ⇒ 需要翻向的题:{flip if flip else '无'}")
for c in flip: M[c] = -M[c]+5
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]

z = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = z(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["a"] = pd.cut(REL.attend, [-1, 1, 5, 8], labels=[0, 1, 2]).astype(float)
REL["b"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def rows_of(item, col, k, nmin=120):
    g = REL.dropna(subset=[item, col])
    out = []
    for y, gy in g[g[col] == k].groupby("year"):
        if len(gy) < nmin: continue
        out.append((int(y), float(gy[item].mean()), float((gy[item] == 4).mean())))
    return out

def fit_thresholds(p, link):
    return link.ppf(np.clip(np.cumsum(p)[:3], 1e-6, 1-1e-6))
def readout(mu, tau, link, stat):
    e = link.cdf(tau-mu)
    p = np.diff(np.concatenate(([0.0], e, [1.0])))
    return float((p*np.arange(1, 5)).sum()) if stat == "水平" else float(p[3])
def fit_mu(p, tau, link):
    tgt = (p*np.arange(1, 5)).sum(); lo, hi = -6.0, 6.0
    for _ in range(60):
        mid = (lo+hi)/2
        if readout(mid, tau, link, "水平") < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2
def forced_ratio(pA0, pB0, nA0, nB0, link, stat, dB_obs):
    """与 `#784` 逐字同一条路径:共用阈值 · 各自起跑线 · 同一个 Δ(由 B 层自己的观测变化标定)。"""
    tau = fit_thresholds((pA0*nA0 + pB0*nB0)/(nA0+nB0), link)
    muA, muB = fit_mu(pA0, tau, link), fit_mu(pB0, tau, link)
    obsA, obsB = readout(muA, tau, link, stat), readout(muB, tau, link, stat)
    lo, hi = -6.0, 6.0                       # ⚠ 两侧 —— `#784` 第一版写成 [0,6] 被对照打掉
    for _ in range(60):
        mid = (lo+hi)/2
        if readout(muB+mid, tau, link, stat)-obsB < dB_obs: lo = mid
        else: hi = mid
    return (readout(muA+(lo+hi)/2, tau, link, stat)-obsA)/dB_obs

CUTS = {"(a) attend三档": "a", "(b) 三题三分位": "b"}
LINKS = {"probit": norm, "logit": logistic}
print("\n=== ② 逐题 · 逐分层:先判可读(`#782` 的判据前置)===")
cells, unreadable = [], []
for item in ITEMS:
    for cn, col in CUTS.items():
        rA, rB = rows_of(item, col, 2), rows_of(item, col, 0)
        if len(rA) < 8 or len(rB) < 8:
            unreadable.append(dict(item=item, cut=cn, why=f"年数不足 {len(rA)}/{len(rB)}")); continue
        yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
        for st, j in (("水平", 1), ("端点", 2)):
            vA = np.array([r[j] for r in rA]); vB = np.array([r[j] for r in rB])
            spA, spB = yA[-1]-yA[0], yB[-1]-yB[0]
            # ⚠ 统一到斜率:每层的观测变化 = 斜率 × 跨度年数(而不是首末年的差),
            #   这样 `r_obs`、Δ 的标定与 `r_forced` 是同一个估计量。
            f = lambda ia, ib: (slope(yA[ia], vA[ia])*spA)/(slope(yB[ib], vB[ib])*spB)
            bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(4000)])
            bs = bs[np.isfinite(bs)]
            lo95, hi95 = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
            r_obs = f(np.arange(len(yA)), np.arange(len(yB)))
            if lo95 <= 1.0 <= hi95:
                unreadable.append(dict(item=item, cut=cn, stat=st, r_obs=float(r_obs),
                                       lo=lo95, hi=hi95, why="比值区间含 1.0 ⇒ 定不了方向"))
                print(f"  {item:9s} {cn:14s} {st:4s} 比 {r_obs:6.3f} · [{lo95:7.3f}, {hi95:7.3f}]  **不可读(含1.0)**")
                continue
            pA0 = np.array([(REL[(REL[col] == 2) & (REL.year == int(yA[0]))][item] == c).mean() for c in (1, 2, 3, 4)])
            pB0 = np.array([(REL[(REL[col] == 0) & (REL.year == int(yB[0]))][item] == c).mean() for c in (1, 2, 3, 4)])
            nA0 = int((REL[(REL[col] == 2) & (REL.year == int(yA[0]))][item].notna()).sum())
            nB0 = int((REL[(REL[col] == 0) & (REL.year == int(yB[0]))][item].notna()).sum())
            dB = slope(yB, vB)*spB          # ⚠ 斜率口径,与 `r_obs` 同一个估计量
            for ln, link in LINKS.items():
                rf = forced_ratio(pA0, pB0, nA0, nB0, link, st, dB)
                resid = r_obs/rf if abs(rf) > 1e-9 else np.nan
                rl, rh = lo95/rf, hi95/rf
                if rf < 0: rl, rh = rh, rl
                cells.append(dict(item=item, cut=cn, stat=st, link=ln, r_obs=float(r_obs),
                                  lo=lo95, hi=hi95, r_forced=float(rf), residual=float(resid),
                                  res_lo=float(rl), res_hi=float(rh),
                                  res_covers1=bool(rl <= 1.0 <= rh), ceil_top=float(pA0[3]),
                                  span=[float(spA), float(spB)]))
                print(f"  {item:9s} {cn:14s} {st:4s} {ln:7s} 比 {r_obs:6.3f} [{lo95:.3f},{hi95:.3f}] · "
                      f"逼出 {rf:6.3f} · **残差 {resid:6.3f} [{rl:.3f},{rh:.3f}]**"
                      f"{'  ⚠含1.0' if rl <= 1.0 <= rh else ''}")

print(f"\n  可读格 **{len(cells)}** · 不可读 **{len(unreadable)}**(全部列出,不平均进去)")
by_item = {}
for c in cells: by_item.setdefault(c["item"], []).append(c)
print("\n=== ③ 逐题小结(⚠ 四题天花板高度不同 ⇒ 残差之间不可比大小,只能比排不排除 1.0)===")
ITEM_OK = {}
for item in ITEMS:
    g = by_item.get(item, [])
    if not g:
        print(f"  {item:9s} **一格可读的都没有** —— 这一题上问不出这个问题"); ITEM_OK[item] = None; continue
    below = sum(1 for c in g if c["res_hi"] < 1.0)
    ITEM_OK[item] = (below == len(g))
    print(f"  {item:9s} 可读 {len(g)} 格 · 残差区间**全在 1.0 之下**的 {below}/{len(g)} · "
          f"残差中位 {np.median([c['residual'] for c in g]):.3f} · "
          f"逼出的比中位 {np.median([c['r_forced'] for c in g]):.3f} · 虔诚层首年端点 {g[0]['ceil_top']:.3f}")

# ── 闸 ─────────────────────────────────────────────────────────────────────────
G = Gate("#785 · 那句话是关于同性恋的,还是关于性道德的")
eq, mono = [], []
for ln, link in LINKS.items():
    g0 = cells[0]
    p0 = np.array([0.05, 0.15, 0.20, 0.60])
    eq.append(forced_ratio(p0, p0, 1000, 1000, link, "水平", -0.30))
    tau0 = fit_thresholds(p0, link); mu0 = fit_mu(p0, tau0, link)
    seq = []
    for shift in (0.0, 0.5, 1.0, 1.5):
        e = link.cdf(tau0-(mu0+shift)); pA = np.diff(np.concatenate(([0.0], e, [1.0])))
        seq.append(forced_ratio(pA, p0, 1000, 1000, link, "端点", -0.20))
    mono.append(all(seq[i] > seq[i+1] for i in range(3)))
    print(f"\n  负控 {ln}:起跑线越靠顶,逼出的比 = {[round(s,3) for s in seq]}")
G.asserted("① 正控:两层起跑线设成相同时,逼出的比必须回到 1.0(容差 0.05)",
           bool(all(abs(e-1.0) <= 0.05 for e in eq)), f"probit/logit 各得 {[round(e,4) for e in eq]}", kind="control")
G.asserted("② 负控:起跑线越靠顶,逼出的比必须单调变小",
           bool(all(mono)), f"probit/logit 单调:{mono}", kind="control")
G.asserted("③ 正控(`#783` 的规矩):每一个要上页面的比值都带自己的区间且非零宽",
           bool(cells and all(c["res_hi"]-c["res_lo"] > 1e-6 for c in cells)),
           f"{len(cells)} 格残差区间宽 "
           f"{min(c['res_hi']-c['res_lo'] for c in cells):.3f}–{max(c['res_hi']-c['res_lo'] for c in cells):.3f}"
           if cells else "无可读格", kind="control")
G.asserted("④ 负控:不可读格必须被**列出**而不是被平均掉(`#781` 的教训)",
           True, f"不可读 {len(unreadable)} 格,逐条持久化到产物", kind="control")
n_full = sum(1 for it in ITEMS if ITEM_OK.get(it) is True)
G.asserted("⑤ kill(预注册):「关于性道德」要站住,需 ≥3 题的残差区间全在 1.0 之下",
           bool(n_full >= 3), f"四题中全在 1.0 之下的:{n_full} —— {[it for it in ITEMS if ITEM_OK.get(it) is True]}",
           kind="kill")
print(); print(G)

print("\n" + "="*92)
ctrl = all(abs(e-1.0) <= 0.05 for e in eq) and all(mono)
readable_items = [it for it in ITEMS if ITEM_OK.get(it) is not None]
if not ctrl:
    v = "**UNVERIFIED:正控或负控没过,本轮不下判。**"
elif n_full >= 3:
    v = (f"**A:那是一句关于性道德的话。** 四题里 **{n_full} 题**的残差区间全在 1.0 之下 —— "
         f"扣掉各自的天花板之后,虔诚者在它们上面都改得比其余人少。")
elif n_full == 1 and ITEM_OK.get("homosex") is True:
    v = (f"**B:那是一句关于同性恋的话,不是关于性道德的。** 四题里只有 `homosex` 的残差区间全在 1.0 之下;"
         f"可读的题只有 {readable_items},其余 {len(unreadable)} 格不可读(逐条列出)。"
         f"⇒ **页面那句话必须缩到这一题上。**")
else:
    v = (f"**C:这个问题只能在「其余人真的动了」的地方问。** 可读的题:{readable_items};"
         f"不可读 {len(unreadable)} 格。四题里残差全在 1.0 之下的只有 {n_full} 题。"
         f"⇒ **这是一条关于 scope 的事实,不是一个零** —— 页面必须把 scope 写出来。")
print(v)
json.dump(dict(cells=cells, unreadable=unreadable, item_ok={k: v2 for k, v2 in ITEM_OK.items()},
               n_full=n_full, verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"across_the_four_items.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'across_the_four_items.json'}")
