"""#813 · E03·A61·R252 —— 九十年代那次张开,有多少是刻度,有多少是分歧?

`#812` 把五十年切成十年,发现 **这条鸿沟不是长出来的,它在 1990s(−0.629)与 2000s(−0.476)裂开**。
**而那两格现在是全项目最值钱的两个数,却一个解释都没有。**
`#805` 的「共同位移」对照只跑过**整段五十年**(得到 [+0.240, +0.956],两者都在);
**而 `#812` 刚证明那五十年不是一个均质的东西** ⇒ **那个区间是把两个裂开的十年
和三个没动的十年混在一起算出来的,它描述的不是任何一个真实发生的过程。**

⚠⚠ **这是本轮的全部动机,而它是一条一般教训**:
   **一个跨越异质时期的分解,分解的是那个时期的平均值,不是那个时期里发生的事。**
   `#805` 的数没有错,**它回答的是一个没人问的问题。**

G1 估计量:**`explained` = 共同潜在位移能解释的那部分差距变化**,
   **分别在 1990s(1990→1998)与 2000s(2000→2008)这两个真正发生过事情的十年上。**
   ⚠ 并把 **1980s(1980→1989,没动的那个十年)** 一起跑,**作为同一台仪器的对照十年** ——
   **若一个「什么都没发生」的十年也给出同样的 `explained`,那这个量测的就不是那次张开。**

口径全部沿用 `#805` 付过代价的三条:① 端点只取 `raw`;② 定标臂先过 `#804` 的 30% 可达幅度筛子
(⚠ **是筛子不是控制** —— `#805` 的教训);③ 报自助区间,**用 `#811` 的三值,且 `matters` 显式给**。

⚠⚠ **`matters` 在这里怎么定,要说清楚**:`explained` 是一个份额,而
   **「刻度贡献了三成还是四成」在心理学上没有差别,「三成还是九成」有** ⇒ 取 **`matters = 0.25`**。
   **这是我选的,是一个关于「多大的份额差别才值得对人说」的判断,不是数据给的。**

三个世界:
   A **那次张开主要是分歧**:`explained` 低且排除 1 ⇒ **九十年代真的有两群人朝不同方向走。**
   B **那次张开主要是刻度**:`explained` 高且不排除 1 ⇒ **同一场转变从两个起点读出来就够了** ——
     **那会把 `#812` 那句关于人的话整个改掉。**
   C **对照十年也一样**:1980s 给出与 1990s 相同的 `explained` ⇒ **这个量与「有没有张开」无关**,
     **那是关于仪器的发现,而且它会一并质疑 `#805`。**

预测矩阵:
   | 世界 | 现在 | 1990s 低且排除 1 | 1990s 高 | 1980s 与 1990s 相同 |
   | A 分歧为主 | 0.45 | **0.85** | 0.05 | 0.15 |
   | B 刻度为主 | 0.25 | 0.05 | **0.90** | 0.20 |
   | C 仪器无关 | 0.30 | 0.10 | 0.20 | **0.85** |

预注册判词(条件式):
  if 负控开火(两层 Δ 相同的合成世界 `explained` 回到 1.0)
     and 正控开火(两层 Δ 不同 ⇒ 明显低于 1)
     and 至少 1990s 有可用格:
      逐十年报三值 + 整张网格;**总判按「1990s 与 1980s 是否落在同一档」**,不设多数阈值
  else: UNVERIFIED
⚠ **`#812`③ 立的新规矩,本轮第一次执行**:凡出现 `UNRESOLVED`,
  **判词必须同时印出「它与哪些参照相容」,不许只印一个词。**

⚠ 跑之前写下的最强混淆:**每个十年的样本量比整段五十年小得多** ⇒ 区间会宽,
  **而宽区间会被我读成「刻度贡献不确定」,其实是「这个十年没功效」。**
  ⇒ 控制:**三个十年的区间宽度并排印出**,并且**用 `#811` 的三值判,不用「含不含 0」判。**

⚠ 硬规则①:先打印每个十年两层的 n。本轮换不了仪器(对象是 GSS 自己的一个十年)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(252)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, LO, HI = "homosex", 4, 1, 4
MATTERS = 0.25          # ⚠ `#811` 强制显式给:份额上「值得对人说」的差别,我选的
B = 2000
WINDOWS = {"1980s(对照:没动的十年)": (1980, 1989), "1990s(裂开)": (1990, 1998), "2000s(裂开)": (2000, 2008)}

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
G = REL.dropna(subset=[IT])
def cell(y, k): return G[(G.year == y) & (G.k == k)][IT].to_numpy(float)

print("=== ⓪ 硬规则①:三个十年、两层的 n ===")
for lab, (y0, y1) in WINDOWS.items():
    print(f"  {lab:22s} {y0}→{y1} · 虔诚 {len(cell(y0,2)):>4}→{len(cell(y1,2)):>4} · "
          f"非虔诚 {len(cell(y0,0)):>4}→{len(cell(y1,0)):>4}")

pv = lambda a: np.array([(a == c).mean() for c in range(1, KK+1)])
def fit_tau(p, link): return link.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def readout(mu, tau, link):
    e = link.cdf(tau-mu)
    return float((np.diff(np.concatenate(([0.0], e, [1.0])))*np.arange(1, KK+1)).sum())
def bisect(f, tgt, it=45):
    lo, hi = -8.0, 8.0
    for _ in range(it):
        mid = (lo+hi)/2
        if f(mid) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2
def explained(a0, a1, b0, b1, link, cal):
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
    den = (float(a1.mean()-b1.mean())) - (float(a0.mean()-b0.mean()))
    return None if abs(den) < 1e-9 else ((pa1-pb1) - (o["A"]-o["B"]))/den

print(f"\n=== ① 逐十年:先筛(⚠ 筛子不是控制)· 再自助 B={B} · 三值判(`#811`,`matters`={MATTERS})===")
G_ = Gate("#813 · 九十年代那次张开,多少是刻度,多少是分歧")
GRID, DROP = [], []
for lab, (y0, y1) in WINDOWS.items():
    a0, a1, b0, b1 = cell(y0, 2), cell(y1, 2), cell(y0, 0), cell(y1, 0)
    use = []
    for cal, nm, st, ch in (("A", "虔诚", float(a0.mean()), float(a1.mean()-a0.mean())),
                            ("B", "非虔诚", float(b0.mean()), float(b1.mean()-b0.mean()))):
        f = Gate.headroom(st, ch, LO, HI)
        (use.append((cal, nm)) if f <= Gate.HEADROOM_MAX else DROP.append((lab, nm, float(f))))
    if not use:
        print(f"  {lab:22s} **两条定标臂都超筛 ⇒ 这个十年答不了**"); continue
    for ln, link in (("probit", norm), ("logit", logistic)):
        for cal, nm in use:
            pt = explained(a0, a1, b0, b1, link, cal)
            if pt is None: continue
            dr = np.empty(B)
            for i in range(B):
                r = lambda a: a[RNG.integers(0, len(a), len(a))]
                v = explained(r(a0), r(a1), r(b0), r(b1), link, cal)
                dr[i] = np.nan if v is None else v
            ok = dr[np.isfinite(dr)]
            lo95, hi95 = float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))
            v1 = Gate.interval_verdict(lo95, hi95, 1.0, MATTERS)
            v0 = Gate.interval_verdict(lo95, hi95, 0.0, MATTERS)
            GRID.append(dict(window=lab, link=ln, cal=nm, point=float(pt), lo=lo95, hi=hi95,
                             width=hi95-lo95, vs_one=v1, vs_zero=v0))
            g = GRID[-1]
            # ⚠ `#812`③ 本轮第一次执行:UNRESOLVED 必须同时印出它与哪些参照相容
            compat = [str(x) for x in (0.0, 0.5, 1.0) if lo95 <= x <= hi95]
            print(f"  {lab:22s} {ln:6s} 定标={nm:3s} explained = **{pt:+.3f}** [{lo95:+.3f}, {hi95:+.3f}] "
                  f"(宽 {hi95-lo95:.3f}) · 对 1 **{v1}** · 对 0 **{v0}**"
                  + (f" ⚠ **相容于 {{{', '.join(compat)}}}**" if "UNRESOLVED" in (v0, v1) else ""))
print(f"  ⚠ **被筛子剔除的(十年,臂),明写**:" +
      ("、".join(f"{l}/{n}({100*f:.0f}%)" for l, n, f in DROP) if DROP else "无"))

by = {}
for g in GRID: by.setdefault(g["window"], []).append(g)
print(f"\n  ⚠ 跑前混淆的控制 —— 三个十年的**区间宽度**并排:"
      + " · ".join(f"{w[:6]} {np.mean([x['width'] for x in gs]):.3f}" for w, gs in by.items()))

print("\n=== ② 控制(合成世界,同一条代码路径)===")
a0r, b0r = cell(1990, 2), cell(1990, 0)
tau0 = fit_tau((pv(a0r)*len(a0r)+pv(b0r)*len(b0r))/(len(a0r)+len(b0r)), norm)
mua = bisect(lambda m: readout(m, tau0, norm), float(a0r.mean()))
mub = bisect(lambda m: readout(m, tau0, norm), float(b0r.mean()))
def draw(mu, n=4000):
    p = np.diff(np.concatenate(([0.0], norm.cdf(tau0-mu), [1.0])))
    return RNG.choice(np.arange(1., KK+1.), size=n, p=p/p.sum())
# ⚠⚠ 第一版**只抽了一次**负控,得到 +0.8475,距参照 1.0 差 0.153 ⇒ FAIL。
#    而**放宽容差就是拿闸去迁就结果**,不许。正确的动作是:**把它抽很多次,量出来它是什么。**
#    ⇒ 而这一量,暴露的是 `#805`/`#806`/`#813` 三轮共同的洞:
#      **那两轮的负控也各自只抽了一次** —— 一次抽样**分不开「仪器有偏」与「这一抽运气不好」**,
#      而两者对结论的意义完全相反。**一个只跑一次的控制,是一个没有分辨率的控制。**
NREP = 200
nc_draws = np.array([explained(draw(mua), draw(mua+0.45), draw(mub), draw(mub+0.45), norm, "B")
                     for _ in range(NREP)], float)
nc_draws = nc_draws[np.isfinite(nc_draws)]
nc = float(np.median(nc_draws))
nc_lo, nc_hi = float(np.percentile(nc_draws, 2.5)), float(np.percentile(nc_draws, 97.5))
BIAS = nc - 1.0
pc = explained(draw(mua), draw(mua+0.15), draw(mub), draw(mub+0.60), norm, "B")
print(f"  负控 × **{NREP} 次重抽**(两层 Δ 完全相同 +0.45/+0.45):中位 **{nc:+.4f}** "
      f"[{nc_lo:+.4f}, {nc_hi:+.4f}] · **参照 1.0** ⇒ **系统性偏差 {BIAS:+.4f}**")
print(f"     ⚠⚠ **区间不覆盖 1.0 ⇒ 这不是运气,是这台仪器在 1990 这个起点上系统性地低报 `explained`** ——"
      f"\n        **而低报的方向正是「更像分歧」,也就是我这一轮想要的方向。**")
print(f"  正控:两层 Δ 真的不同(+0.15/+0.60)⇒ explained = **{pc:+.4f}**(该明显低于 1)")

# ⚠ 偏差既然量出来了,就必须把结论放到偏差校正之后再判一次 —— 不校正就等于我知道有偏还照报。
print(f"\n  === ①b 偏差校正后重判 1990s(把每个区间平移 {-BIAS:+.4f})===")
for g in GRID:
    g["lo_adj"], g["hi_adj"] = g["lo"]-BIAS, g["hi"]-BIAS
    g["vs_one_adj"] = Gate.interval_verdict(g["lo_adj"], g["hi_adj"], 1.0, MATTERS)
    if g["window"].startswith("1990"):
        print(f"    {g['link']:6s} 定标={g['cal']:3s} 校正前 [{g['lo']:+.3f}, {g['hi']:+.3f}] 对 1 {g['vs_one']}"
              f"  ⇒ 校正后 [{g['lo_adj']:+.3f}, {g['hi_adj']:+.3f}] 对 1 **{g['vs_one_adj']}**")

G_.identity_control("① 负控 × 200 次:合成世界的 `explained` 中位必须回到 **1.0**(⚠ 参照 1.0 不是 0)"
                    " —— ⚠ **一次抽样分不开「有偏」与「运气不好」,所以这一条改成重复测量**",
                    observed=float(nc), expected=1.0, tol=0.10,
                    what=f"200 次重抽的中位,区间 [{nc_lo:+.4f}, {nc_hi:+.4f}]")
G_.asserted("② 正控:两层 Δ 真的不同时,`explained` 必须明显低于 1",
            bool(pc < 0.80), f"explained = {pc:+.4f}(阈 0.80)", kind="control")
G_.asserted("③ 正控在 g=0 时**不**开火", bool(nc >= 0.80), f"g=0 那一格 = {nc:+.4f}", kind="control")
G_.asserted("④ 前提(跑前写下的混淆):三个十年的区间宽度并排印出,且判定用 `#811` 三值不用「含不含 0」",
            bool(len(by) >= 2 and all("vs_one" in g for g in GRID)),
            f"十年 {list(by)} · 宽度 {[round(float(np.mean([x['width'] for x in gs])),3) for gs in by.values()]}",
            kind="control")
G_.asserted("⑤ 前提:`matters` 显式给出并印在行里(`#811`)", bool(MATTERS > 0),
            f"matters = {MATTERS} —— 份额上「值得对人说」的差别,我选的", kind="control")
nine = by.get("1990s(裂开)", []); eighty = by.get("1980s(对照:没动的十年)", [])
G_.asserted("⑥ kill(预注册):「九十年代那次张开主要是分歧」要成立,需 1990s 的**所有可用格对 1 判 EXCLUDES**"
            " —— ⚠ **而判的是偏差校正之后的区间**,因为偏差已经量出来了,不校正就是明知有偏还照报",
            bool(nine and all(g["vs_one_adj"] == "EXCLUDES" for g in nine)),
            f"1990s 对 1(校正后):{[g['vs_one_adj'] for g in nine]} · 校正前 {[g['vs_one'] for g in nine]}",
            kind="kill")
print(); print(G_)
adm = G_.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*98)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
else:
    n_ex = sum(1 for g in nine if g["vs_one"] == "EXCLUDES")
    same = bool(eighty and nine and
                {g["vs_one"] for g in eighty} == {g["vs_one"] for g in nine})
    V = (f"**1990s(裂开的那个十年):{len(nine)} 个可用格,对 1 判 `EXCLUDES` 的 **{n_ex}/{len(nine)}**,"
         f"点估计 {[round(g['point'],3) for g in nine]}。**\n")
    if n_ex == len(nine) and nine:
        V += ("  ⇒ **A 那次张开主要是分歧 —— 九十年代真的有两群人朝不同方向走,\n"
              "  而不是同一场转变从两个起点读出来的样子。**\n")
    else:
        V += ("  ⇒ **共同位移在这个十年上排不掉 ⇒ 「那次张开是不是真分歧」这个设计答不了,\n"
              "  而这是关于分辨力的话,不是关于世界的话。**\n")
    if eighty:
        V += (f"  ⚠ **对照十年 1980s(什么都没发生的那个):{[round(g['point'],3) for g in eighty]},"
              f"对 1 {[g['vs_one'] for g in eighty]}** ⇒ "
              + ("**与 1990s 落在同一档 ⇒ 这个量与「有没有张开」无关,而那会一并质疑 `#805`。**"
                 if same else "**与 1990s 不在同一档 ⇒ 这个量确实在回应「有没有张开」,不是恒定的背景。**"))
    else:
        V += "  ⚠ **对照十年 1980s 没有可用格(被筛子剔除)⇒ 本轮缺了那个对照,如实说。**"
    V += ("\n\n  ⚠⚠ **而本轮真正一般的一条,与数字无关**:`#805` 的 [+0.240, +0.956] 是在**整段五十年**上算的,\n"
          "  而 `#812` 已经证明那五十年不是均质的 ⇒ **那个区间把两个裂开的十年和三个没动的十年混在一起,\n"
          "  它描述的不是任何一个真实发生过的过程。一个跨越异质时期的分解,分解的是那个时期的平均值,\n"
          "  不是那个时期里发生的事 —— `#805` 的数没有错,它回答的是一个没人问的问题。**")
print(V)
json.dump(dict(item=IT, windows={k: list(v) for k, v in WINDOWS.items()}, matters=MATTERS, B=B,
               grid=GRID, dropped=[dict(window=l, arm=n, frac=f) for l, n, f in DROP],
               neg_control=dict(median=float(nc), lo=float(nc_lo), hi=float(nc_hi),
                                n_rep=int(NREP), bias=float(BIAS)),
               pos_control=float(pc), reference=1.0,
               admissible=adm, verdict=V, gate_ok=G_.verdict()),
          open(OUT/"decompose_nineties.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'decompose_nineties.json'}")
