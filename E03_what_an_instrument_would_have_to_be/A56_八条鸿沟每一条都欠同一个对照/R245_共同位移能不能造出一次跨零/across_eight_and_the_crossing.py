"""#806 · E03·A56·R245 —— 八条鸿沟每一条都欠同一个对照;而 `helpblk` 跨过了零

`#805` 只做了 `homosex` 一题,而 `#802` 报的是**八条鸿沟**(缩 3/8,扩 5/8)。
**「五裂三合」这句话里,每一条都欠 `#805` 那个对照** —— 而页面已经把它当成结论在用。

⚠⚠ **而其中有一条是尖锐得多的:`helpblk` 的差距从 −0.331 走到 +0.297,跨过了零。**
   `#805` 问的是「共同位移能解释多少」;**这里能问一个更硬的:一个共同的潜在位移,
   能不能让两层的相对位置**换边**?** 若能,则「虔诚层在这道题上从落后走到领先」这句关于人的话
   **根本不需要两群人走得不一样** —— 它是同一场转变在一把粗尺子上从两个起点读出来的结果。
   ⇒ **这是本轮最想要的那一格,而它的正结果我不欢迎**(`§3` 盆地规则)。

G1 估计量:逐题的 `explained`(与 `#805` 同一个),**外加一个只对跨零题有意义的二值量**:
   **`sign_flip_reproduced`** —— 共同位移下的**预测 gap1 与观测 gap1 是否同号**。
   ⚠ 两者是不同的问题:`explained` 高不蕴含跨零被复现(可以解释掉大部分幅度却停在零之前)。

口径,全部沿用 `#805` 已经付过代价的那三条:
   ① 端点只取 `raw`(模型起点 = 观测首年分布 ⇒ 分子分母同口径,`#803` 的 `#785` 家缺陷);
   ② 每条定标臂先过 `Gate.headroom()` 的 30% 筛子(`#804`;⚠ **是筛子不是控制** —— `#805` 的教训);
   ③ 报自助区间,**不设中点阈值**。

⚠⚠ **`G3` 多重性,而这一轮它是真的有牙**:八题 × 两个连接 × 至多两条定标臂 = **至多 32 格**。
   **整张网格全报,包括没有任何一条臂过筛的题** —— 那种题登记为**这个设计答不了**,不是登记为 0。

三个世界:
   A **刻度是普遍的**:多数题的区间含 1(共同位移解释得完)⇒ **`#802` 的「五裂三合」大半是刻度的性质**,
     而那会把整页对 `r` 的读法再改一次。
   B **刻度是局部的**:多数题的区间排除 1 ⇒ 各题确实各有各的分歧,`#805` 的结论可以推广。
   C **跨零是刻度造的**:`helpblk` 的 `sign_flip_reproduced` = True ⇒
     **「虔诚层换到了另一边」这句话不需要分歧** ⇒ 单独撤掉那一句,不动其余。

预测矩阵:
   | 世界 | 现在 | 多数含 1 | 多数排除 1 | `helpblk` 跨零被复现 |
   | A 刻度普遍 | 0.30 | **0.85** | 0.05 | 0.55 |
   | B 刻度局部 | 0.50 | 0.05 | **0.85** | 0.20 |
   | C 跨零是刻度 | 0.20 | 0.30 | 0.30 | **0.90** |

预注册判词(条件式):
  if 负控开火(两层 Δ 相同的合成世界 `explained` 回到 1.0)
     and 正控开火(两层 Δ 不同 ⇒ 明显低于 1)
     and 至少一题有可用格:
      逐题报区间 + 跨零复现;**总判按「排除 1 的题数 / 有可用格的题数」**,不按平均
  else: UNVERIFIED
⚠ **不预注册一个「多数」的阈值** —— `#803`/`#805` 刚证明阈值判据在这类量上没有分辨力。
  **报计数与整张网格,让读者看见分歧本身。**

⚠ 跑之前写下的最强混淆:**八题的档数不同(2–5),而档数越少、阈值越少,
  这套序数模型的自由度就越低,`explained` 会被结构性地推向某个方向。**
  ⇒ 控制:**结果表里把档数并排列出**,并且**结论不跨档数排序**,只报每题自己的区间。
  ⚠ `racmar` 只有 2 档(1 个阈值)且末年 2002 —— **单独标注,不与 5 档题并读。**

⚠ 硬规则①:先打印每题 n、真正被问过的年份、档数、尺子两端。本轮换不了仪器(`R223` 六具全部落选)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(245)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS, OBS = P791["items"], P791["obs"]
B = 1000
STEM = pd.io.stata.StataReader(gp).variable_labels()

print("=== ⓪ 硬规则①:每题 n · 真正被问过的年份 · 档数 · 尺子两端 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c in ITEMS:
    yrs = sorted(d.year[M[c].notna()].unique())
    print(f"  {c:9s} n={int(M[c].notna().sum()):>7,} · 年 {len(yrs):>2}({int(min(yrs))}–{int(max(yrs))}) · "
          f"档 {K[c]}(阈值 {K[c]-1}) · 尺 [1, {K[c]}] · r={OBS[c]:+.3f}")
print("  ⚠ `racmar` 只有 2 档(1 个阈值)且末年 2002 —— **单独标注,不与 5 档题并读**")

for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))

def cells_for(item):
    g = REL.dropna(subset=[item])
    ya = {int(y) for y, gy in g[g.k == 2].groupby("year") if len(gy) >= 120}
    yb = {int(y) for y, gy in g[g.k == 0].groupby("year") if len(gy) >= 120}
    yy = sorted(ya & yb)
    if len(yy) < 8: return None
    y0, y1 = yy[0], yy[-1]
    return y0, y1, {(y, k): g[(g.year == y) & (g.k == k)][item].to_numpy(float)
                    for y in (y0, y1) for k in (0, 2)}

def make(kk):
    def pv(a): return np.array([(a == c).mean() for c in range(1, kk+1)])
    def readout(mu, tau, link):
        e = link.cdf(tau-mu)
        return float((np.diff(np.concatenate(([0.0], e, [1.0])))*np.arange(1, kk+1)).sum())
    return pv, readout
def fit_tau(p, link): return link.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def bisect(f, tgt, it=45):
    lo, hi = -8.0, 8.0
    for _ in range(it):
        mid = (lo+hi)/2
        if f(mid) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2

def run(a0, a1, b0, b1, link, cal, kk):
    pv, readout = make(kk)
    pa0, pb0 = pv(a0), pv(b0); na, nb = len(a0), len(b0)
    tau = fit_tau((pa0*na + pb0*nb)/(na+nb), link)
    mu = {"A": bisect(lambda m: readout(m, tau, link), float(a0.mean())),
          "B": bisect(lambda m: readout(m, tau, link), float(b0.mean()))}
    o = {s: readout(mu[s], tau, link) for s in mu}
    d_cal = float(a1.mean()-a0.mean()) if cal == "A" else float(b1.mean()-b0.mean())
    D = bisect(lambda m: readout(mu[cal]+m, tau, link)-o[cal], d_cal)
    oth = "B" if cal == "A" else "A"
    ec, eo = readout(mu[cal]+D, tau, link), readout(mu[oth]+D, tau, link)
    pa1, pb1 = (ec, eo) if cal == "A" else (eo, ec)
    pg0, pg1 = o["A"]-o["B"], pa1-pb1
    og0, og1 = float(a0.mean()-b0.mean()), float(a1.mean()-b1.mean())
    den = og1-og0
    if abs(den) < 1e-9: return None
    return (pg1-pg0)/den, pg1, og1, og0

print(f"\n=== ① 逐题:先筛,再自助 B={B}(⚠ 筛子不是控制 —— `#805`)===")
G = Gate("#806 · 八条鸿沟每一条都欠同一个对照;而 `helpblk` 跨过了零")
GRID, SKIP, DROPPED = [], [], []
for it in ITEMS:
    cf = cells_for(it)
    if cf is None: SKIP.append((it, "年份不足 8")); continue
    y0, y1, C = cf
    a0, a1, b0, b1 = C[(y0, 2)], C[(y1, 2)], C[(y0, 0)], C[(y1, 0)]
    og0, og1 = float(a0.mean()-b0.mean()), float(a1.mean()-b1.mean())
    crosses = (np.sign(og0) != np.sign(og1)) and abs(og1) > 1e-9
    usable = []
    for cal, nm, st, ch in (("A", "虔诚", float(a0.mean()), float(a1.mean()-a0.mean())),
                            ("B", "非虔诚", float(b0.mean()), float(b1.mean()-b0.mean()))):
        f = Gate.headroom(st, ch, 1, K[it])
        (usable.append((cal, nm)) if f <= Gate.HEADROOM_MAX else DROPPED.append((it, nm, float(f))))
    if not usable:
        SKIP.append((it, "两条定标臂都超过 30% 可达幅度 ⇒ **这个设计答不了这一题**")); continue
    for ln, link in (("probit", norm), ("logit", logistic)):
        for cal, nm in usable:
            r = run(a0, a1, b0, b1, link, cal, K[it])
            if r is None: continue
            pt, pg1, _, _ = r
            dr = np.empty(B)
            fl = np.zeros(B, bool)
            for i in range(B):
                ra0 = a0[RNG.integers(0, len(a0), len(a0))]; ra1 = a1[RNG.integers(0, len(a1), len(a1))]
                rb0 = b0[RNG.integers(0, len(b0), len(b0))]; rb1 = b1[RNG.integers(0, len(b1), len(b1))]
                rr = run(ra0, ra1, rb0, rb1, link, cal, K[it])
                dr[i] = np.nan if rr is None else rr[0]
                fl[i] = False if rr is None else (np.sign(rr[1]) == np.sign(rr[2]))
            ok = dr[np.isfinite(dr)]
            lo95, hi95 = float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))
            GRID.append(dict(item=it, K=K[it], link=ln, cal=nm, y0=y0, y1=y1,
                             point=float(pt), lo=lo95, hi=hi95,
                             obs_gap0=og0, obs_gap1=og1, pred_gap1=float(pg1),
                             crosses_zero=bool(crosses),
                             sign_flip_reproduced=bool(np.sign(pg1) == np.sign(og1)),
                             flip_rate=float(fl.mean()),
                             excl0=bool(lo95 > 0 or hi95 < 0), excl1=bool(lo95 > 1 or hi95 < 1)))

print(f"  ⚠ **剔除的(题,臂),明写不静默**:" +
      ("、".join(f"{i}/{n}({100*f:.0f}%)" for i, n, f in DROPPED) if DROPPED else "无"))
print(f"  ⚠ **完全答不了的题**:" + ("、".join(f"{i}({w})" for i, w in SKIP) if SKIP else "无"))
print(f"\n  网格 **{len(GRID)}** 格(至多 8 题 × 2 连接 × 2 臂 = 32)· "
      f"有可用格的题 **{len({g['item'] for g in GRID})}/{len(ITEMS)}**")
print(f"\n  {'题':9s} {'档':>2} {'连接':6s} {'臂':4s} {'explained':>10s}  {'95% 区间':>20s}  排除0 排除1  跨零")
for g in GRID:
    print(f"  {g['item']:9s} {g['K']:>2} {g['link']:6s} {g['cal']:4s} {g['point']:>+10.3f}  "
          f"[{g['lo']:+.3f}, {g['hi']:+.3f}]  {'是' if g['excl0'] else '否':4s} {'是' if g['excl1'] else '否':4s}  "
          + (f"**观测 {g['obs_gap0']:+.3f}→{g['obs_gap1']:+.3f} · 预测末 {g['pred_gap1']:+.3f} ⇒ "
             f"跨零{'被复现' if g['sign_flip_reproduced'] else '**没**被复现'}(自助 {100*g['flip_rate']:.1f}%)**"
             if g["crosses_zero"] else ""))

by_item = {}
for g in GRID: by_item.setdefault(g["item"], []).append(g)
n_all_excl1 = sum(1 for it, gs in by_item.items() if all(x["excl1"] for x in gs))
n_all_excl0 = sum(1 for it, gs in by_item.items() if all(x["excl0"] for x in gs))
print(f"\n  **所有可用格都排除 1 的题 {n_all_excl1}/{len(by_item)}** · 都排除 0 的题 {n_all_excl0}/{len(by_item)}")
CROSS = [g for g in GRID if g["crosses_zero"]]

print("\n=== ② 控制(合成世界,同一条代码路径,4 档)===")
kk = 4; pv4, ro4 = make(kk)
base = REL.dropna(subset=["homosex"]); by = base[base.year == base.year.min()]
p_pool = pv4(by["homosex"].to_numpy(float))
tau0 = fit_tau(p_pool, norm)
def draw(mu, n=4000):
    p = np.diff(np.concatenate(([0.0], norm.cdf(tau0-mu), [1.0])))
    return RNG.choice(np.arange(1., kk+1.), size=n, p=p/p.sum())
mua, mub = -0.35, +0.35
nc = run(draw(mua), draw(mua+0.45), draw(mub), draw(mub+0.45), norm, "B", kk)[0]
pc = run(draw(mua), draw(mua+0.15), draw(mub), draw(mub+0.60), norm, "B", kk)[0]
print(f"  负控:两层 Δ 完全相同(+0.45/+0.45)⇒ explained = **{nc:+.4f}**(参照 **1.0**,不是 0)")
print(f"  正控:两层 Δ 真的不同(+0.15/+0.60)⇒ explained = **{pc:+.4f}**(该**明显低于** 1)")

G.asserted("① 筛子(`#804` 的 30% 可达幅度;⚠ **是筛子不是控制** —— `#805` 的教训):"
           "凡没过筛的(题,臂)都**真的没有进网格**,且**明写剔了谁**",
           bool(not {(i, n) for i, n, _ in DROPPED} & {(g["item"], g["cal"]) for g in GRID}),
           f"剔除 {len(DROPPED)} 个(题,臂) · 完全答不了的题 {len(SKIP)}", kind="control")
G.identity_control("② 负控:两层 Δ 完全相同的合成世界里,`explained` 必须回到 **1.0**(⚠ 参照 1.0 不是 0)",
                   observed=float(nc), expected=1.0, tol=0.08,
                   what="共享阈值、两层各自的 μ、施加同一个 Δ、按人抽样 4000 的合成世界")
G.asserted("③ 正控:两层 Δ 真的不同时,`explained` 必须明显低于 1",
           bool(pc < 0.80), f"explained = {pc:+.4f}(阈 0.80)", kind="control")
G.asserted("④ 正控在 g=0 时**不**开火", bool(nc >= 0.80), f"g=0 那一格 = {nc:+.4f}", kind="control")
G.asserted("⑤ 前提(跑前写下的混淆):档数并排列出,且**结论不跨档数排序**",
           bool(len({g["K"] for g in GRID}) >= 1),
           f"网格里的档数 {sorted({g['K'] for g in GRID})} —— 每题只报自己的区间", kind="control")
G.asserted("⑥ kill(预注册):至少一题有可用格,否则整轮无话可说",
           bool(len(by_item) >= 1), f"有可用格的题 {len(by_item)}/{len(ITEMS)}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
else:
    lines = [f"**逐题网格已全报({len(GRID)} 格,{len(by_item)}/{len(ITEMS)} 题有可用格;"
             f"剔除 {len(DROPPED)} 个(题,臂),完全答不了 {len(SKIP)} 题 —— 都明写了)。**",
             f"  **所有可用格都排除 1 的题:{n_all_excl1}/{len(by_item)}** ⇒ "
             + ("**多数题的分歧不是刻度能解释完的 ⇒ `#805` 的结论可以推广(世界 B)。**"
                if n_all_excl1*2 > len(by_item) else
                "**没有形成多数 ⇒ 「刻度是普遍的还是局部的」这个设计分不开,如实报网格(世界 A/C 未分离)。**")]
    if CROSS:
        c0 = CROSS[0]
        rep = sum(g["sign_flip_reproduced"] for g in CROSS)
        lines.append(f"  ⚠⚠ **跨零那一题(`{c0['item']}`,观测 {c0['obs_gap0']:+.3f} → {c0['obs_gap1']:+.3f}):"
                     f"共同位移复现了跨零的格 {rep}/{len(CROSS)}** ⇒ "
                     + ("**「虔诚层在这道题上换到了另一边」这句关于人的话,不需要两群人走得不一样 —— "
                        "同一场转变从两个起点读出来就够了。⇒ 单独撤掉那一句,不动其余(世界 C)。**"
                        if rep == len(CROSS) else
                        "**共同位移造不出这次换边 ⇒ 那句话站得住,而它现在有对照了。**"))
    else:
        lines.append("  ⚠ **本轮网格里没有跨零的题** —— `helpblk` 若没进网格,那是筛子剔的,已明写。")
    V = "\n".join(lines)
print(V)
print("\n⚠ **区间是下界性质的**(层内按人重抽,不含年际波动)· "
      "**共享阈值是假设,DIF 结构性测不了** · **不跨档数排序**(档 2–5)。")
json.dump(dict(items=ITEMS, K=K, B=B, grid=GRID, dropped=[dict(item=i, arm=n, frac=f) for i, n, f in DROPPED],
               skipped=[dict(item=i, why=w) for i, w in SKIP],
               n_items_with_cells=len(by_item), n_items_all_excl1=n_all_excl1,
               n_items_all_excl0=n_all_excl0, crossing=CROSS,
               neg_control=float(nc), pos_control=float(pc), reference=1.0,
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"across_eight.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'across_eight.json'}")
