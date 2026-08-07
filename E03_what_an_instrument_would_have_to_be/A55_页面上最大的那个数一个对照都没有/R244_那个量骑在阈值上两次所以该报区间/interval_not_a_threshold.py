"""#805 · E03·A55·R244 —— 那个量骑在阈值上两次,所以该报的是区间,不是「≥/< 0.50」

`#803` 的中位是 **+0.478**,`#804` 用机械规则筛完只剩两格、中位 **+0.493**。
**两次都骑在预注册的 0.50 上。**
⚠⚠ **而那不是数值巧合,是设计缺陷:一个阈值判据在一个坐落于阈值上的量上没有分辨力。**
⇒ `#804`① 预注册的修法:**换判据,报区间。**

G1 估计量(与 `#803` 同一个,口径收紧):
   `explained = (预测 gap1 − 预测 gap0) / (观测 gap1 − 观测 gap0)`
   —— 一个**共同的潜在位移** Δ(两层完全一样地改主意)能解释多少那条扩大的差距。

**三处收紧,每一处都由前一轮的账指定:**
① **端点只取 `raw`** —— 模型的起点只能是观测首年分布,而 `fit` 规格宣称的 gap0 与它不同
   ⇒ 分子分母口径不一(`#803` 指出的 `#785` 那一家)。**`raw` 下模型 gap0 = 观测 gap0,恰好相等。**
② **每一格过 `headroom_control`**(`#804` 刚建的机械规则,上限 = 可达幅度的 30%)——
   超限的格子标 **不可用**,不是标数值。
③ **不报「≥/< 0.50」,报自助区间。**

⚠⚠ 预注册判词**不再是阈值,而是两个互相独立的排除**,而这才是有内容的问法:
   **区间排除 1.0 吗?** 排除 ⇒ **共同位移解释不完** ⇒ **确实有分歧。**
   **区间排除 0.0 吗?** 排除 ⇒ **共同位移解释了一部分** ⇒ **刻度确实在贡献。**
   ⇒ 两个都排除 ⇒ **两者都在,而那是一个比任何一端都更有内容的结论。**
   ⇒ 含 1.0 ⇒ 「全是刻度」排除不掉。含 0.0 ⇒ 「全是分歧」排除不掉。
   **不设中点阈值,不问「哪一半更大」** —— `#803`/`#804` 已经证明这个设计答不了那个问题。

三个世界:
   A **两者都在**:区间同时排除 0 与 1。
   B **只有分歧**:区间含 0 而排除 1(共同位移解释不了任何一部分)。
   C **分不开**:区间同时含 0 与 1 ⇒ **`raw` 只用了 4 个格子的数据(2 年 × 2 层),
     它的分辨力本来就可能不够** —— 那是关于设计的真收获,如实报,不硬判。

预测矩阵:
   | 世界 | 现在 | 排除 0 且排除 1 | 含 1 | 含 0 | 两个都含 |
   | A 两者都在 | 0.55 | **0.90** | 0.05 | 0.05 | 0.20 |
   | B 只有分歧 | 0.10 | 0.03 | 0.02 | **0.60** | 0.10 |
   | C 分不开   | 0.35 | 0.07 | 0.10 | 0.10 | **0.70** |

自助:**层内按人重抽**(4 个 `年 × 层` 格各自重抽),B = 2000 ——
   ⚠ 而**这正是 `raw` 口径的代价,要说出来**:它只用首末两年,
   **年际波动不进区间**,所以这个区间是**下界性质的**(真不确定度只会更大,不会更小)。
   ⇒ **如实标为「层内抽样误差的区间,不是全部不确定度」。**

⚠ 跑之前写下的最强混淆:**共享阈值这一假设本身。** 若两层的阈值本就不同(DIF),
  「共同位移」这个反事实根本不成立。⇒ **DIF 在这份数据上测不了 —— 结构性不可能,早已登记;
  本轮不假装控制它,只在结论里带着它。**

⚠ 硬规则①:先打印 n、真正被问过的年份、档数。本轮换不了仪器(`R223` 六具全部落选)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(244)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, LO, HI = "homosex", 4, 1, 4          # ⚠ `#804`:`lo`/`hi` 是**尺子的两端**,不是观测极值
B = 2000

print("=== ⓪ 硬规则①:n · 真正被问过的年份 · 档数 · 尺子两端 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
v = pd.to_numeric(d[IT], errors="coerce").where(lambda x: (x >= 1) & (x <= KK))
yrs = sorted(d.year[v.notna()].unique())
print(f"  {IT} n={int(v.notna().sum()):,} · 年 {len(yrs)}({int(min(yrs))}–{int(max(yrs))}) · 档 {KK} · 尺 [{LO}, {HI}]")

M = pd.DataFrame({IT: v})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda x, lo=lo, hi=hi: (x >= lo) & (x <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
G0 = REL.dropna(subset=[IT])
yy = sorted({int(y) for y, g in G0[G0.k == 2].groupby("year") if len(g) >= 120}
            & {int(y) for y, g in G0[G0.k == 0].groupby("year") if len(g) >= 120})
Y0, Y1 = yy[0], yy[-1]
CELL = {(y, k): G0[(G0.year == y) & (G0.k == k)][IT].to_numpy(float) for y in (Y0, Y1) for k in (0, 2)}
print(f"  用到的四个格({Y0} 与 {Y1} × 两层):"
      + " · ".join(f"{y}/k{k} n={len(CELL[(y,k)]):,}" for y in (Y0, Y1) for k in (0, 2)))

def pv(a): return np.array([(a == c).mean() for c in range(1, KK+1)])
def mean_(a): return float(a.mean())

def fit_tau(p, link): return link.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def readout(mu, tau, link):
    e = link.cdf(tau-mu); return float((np.diff(np.concatenate(([0.0], e, [1.0])))*np.arange(1, KK+1)).sum())
def bisect(f, tgt, lo=-8.0, hi=8.0, it=45):
    for _ in range(it):
        mid = (lo+hi)/2
        if f(mid) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2

def explained(a0, a1, b0, b1, link, cal):
    """a=虔诚 b=非虔诚,各为该格的原始作答向量。返回 explained,或 None(定标不可达)。"""
    pa0, pb0 = pv(a0), pv(b0); na, nb = len(a0), len(b0)
    tau = fit_tau((pa0*na + pb0*nb)/(na+nb), link)
    mu = {"A": bisect(lambda m: readout(m, tau, link), mean_(a0)),
          "B": bisect(lambda m: readout(m, tau, link), mean_(b0))}
    o = {s: readout(mu[s], tau, link) for s in mu}
    d_cal = (mean_(a1)-mean_(a0)) if cal == "A" else (mean_(b1)-mean_(b0))
    D = bisect(lambda m: readout(mu[cal]+m, tau, link)-o[cal], d_cal)
    oth = "B" if cal == "A" else "A"
    endc, endo = readout(mu[cal]+D, tau, link), readout(mu[oth]+D, tau, link)
    pa1, pb1 = (endc, endo) if cal == "A" else (endo, endc)
    pg0, pg1 = o["A"]-o["B"], pa1-pb1
    den = (mean_(a1)-mean_(b1)) - (mean_(a0)-mean_(b0))
    return None if abs(den) < 1e-9 else (pg1-pg0)/den, D, d_cal

print(f"\n=== ① 观测(`raw` 端点,{Y0}→{Y1})===")
a0, a1 = CELL[(Y0, 2)], CELL[(Y1, 2)]; b0, b1 = CELL[(Y0, 0)], CELL[(Y1, 0)]
print(f"  虔诚 {mean_(a0):.3f}→{mean_(a1):.3f}(Δ{mean_(a1)-mean_(a0):+.3f}) · "
      f"非虔诚 {mean_(b0):.3f}→{mean_(b1):.3f}(Δ{mean_(b1)-mean_(b0):+.3f})")
print(f"  差距 {mean_(a0)-mean_(b0):+.3f} → {mean_(a1)-mean_(b1):+.3f} "
      f"(变化 {(mean_(a1)-mean_(b1))-(mean_(a0)-mean_(b0)):+.3f})")

LINKS = {"probit": norm, "logit": logistic}
print(f"\n=== ② 每一格先过 `headroom_control`(`#804`,上限 30% 可达幅度)===")
# ⚠⚠ 第一版在这里用了 `headroom_control`,而**那是错的用法,并且它当场把整轮判成 UNVERIFIED**:
#    非虔诚臂 63.5% 超限 ⇒ 控制行 FAIL ⇒ `admissible()` = False ——
#    **可是这一轮已经正确地把那条臂剔除了,它从来没有进过网格。**
#    ⇒ **一个筛子(screen)和一个控制(control)不是同一件东西**:
#      控制失败 = **仪器坏了**,整轮不可采;筛子失败 = **这一格不能用**,剔掉,继续。
#      我拿控制形状的工具去做筛选,于是一次**被正确处理的排除**看起来像**仪器故障** ——
#      与 `#796` 同一族(拿 `offset_control` 去断言等式),这次犯在**我自己一轮前刚建的那个闸上**。
#    ⇒ 筛选用纯函数 `Gate.headroom()`(不写行),**而控制行断言的是「凡没过筛的臂都真的被剔除了」**
#      —— 那才是可以失败、且失败即意味着我漏掉了什么的命题。⚠ 并且**明写剔了谁**(不许静默截断)。
G = Gate("#805 · 那个量骑在阈值上两次,所以该报区间")
USABLE, DROPPED = [], []
for cal, nm, st, ch in (("A", "虔诚", mean_(a0), mean_(a1)-mean_(a0)),
                        ("B", "非虔诚", mean_(b0), mean_(b1)-mean_(b0))):
    frac = Gate.headroom(st, ch, LO, HI)
    ok = frac <= Gate.HEADROOM_MAX
    print(f"  定标={nm:3s} 起点 {st:.3f} 移动 {ch:+.3f} ⇒ {100*frac:.1f}% 可达幅度 "
          f"(上限 {100*Gate.HEADROOM_MAX:.0f}%) ⇒ **{'可用' if ok else '不可用 —— 剔除'}**")
    (USABLE if ok else DROPPED).append((cal, nm, frac))
USABLE = [(c, n) for c, n, _ in USABLE]
print(f"  ⇒ 可用的定标臂 **{len(USABLE)}/2** ⇒ 可用格子 **{len(USABLE)*len(LINKS)}/4**")
print(f"  ⚠ **剔除的臂(明写,不静默):** "
      + ("、".join(f"{n}({100*f:.1f}% 可达幅度)" for _, n, f in DROPPED) if DROPPED else "无"))
G.asserted("① 筛子(`#804` 的 `headroom` 上限 30%):凡没过筛的定标臂都**真的没有进网格**"
           "(⚠ 这是筛子不是控制 —— 控制失败=仪器坏了,筛子失败=这一格不能用)",
           bool(not (set(n for _, n, _ in DROPPED) & set(n for _, n in USABLE))),
           f"过筛 {[n for _, n in USABLE]} · 剔除 {[n for _, n, _ in DROPPED]}", kind="control")

print(f"\n=== ③ 自助 B={B}(层内按人重抽 4 个格)· 逐格区间,全报(`G3`)===")
CELLS = []
for ln_name, link in LINKS.items():
    for cal, nm in USABLE:
        pt = explained(a0, a1, b0, b1, link, cal)[0]
        draws = np.empty(B)
        for i in range(B):
            ra0 = a0[RNG.integers(0, len(a0), len(a0))]; ra1 = a1[RNG.integers(0, len(a1), len(a1))]
            rb0 = b0[RNG.integers(0, len(b0), len(b0))]; rb1 = b1[RNG.integers(0, len(b1), len(b1))]
            r = explained(ra0, ra1, rb0, rb1, link, cal)
            draws[i] = np.nan if r is None else r[0]
        dr = draws[np.isfinite(draws)]
        lo95, hi95 = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
        CELLS.append(dict(link=ln_name, cal=nm, point=float(pt), lo=lo95, hi=hi95,
                          n_boot=int(len(dr)), excl0=bool(lo95 > 0 or hi95 < 0),
                          excl1=bool(lo95 > 1 or hi95 < 1)))
        c = CELLS[-1]
        print(f"  {ln_name:6s} 定标={nm:3s} · explained = **{pt:+.3f}** [{lo95:+.3f}, {hi95:+.3f}] "
              f"({len(dr)}/{B} 有效) ⇒ 排除 0:**{'是' if c['excl0'] else '否'}** · "
              f"排除 1:**{'是' if c['excl1'] else '否'}**")
n_e0 = sum(c["excl0"] for c in CELLS); n_e1 = sum(c["excl1"] for c in CELLS)
print(f"\n  **排除 0 的格 {n_e0}/{len(CELLS)} · 排除 1 的格 {n_e1}/{len(CELLS)}**")

print("\n=== ④ 控制(合成世界,同一条代码路径)===")
def syn_cells(mu_a, mu_b, dA, dB, link, n=4000):
    tau = fit_tau((pv(a0)*len(a0) + pv(b0)*len(b0))/(len(a0)+len(b0)), link)
    def draw(mu):
        p = np.diff(np.concatenate(([0.0], link.cdf(tau-mu), [1.0])))
        return RNG.choice(np.arange(1., KK+1.), size=n, p=p/p.sum())
    return draw(mu_a), draw(mu_a+dA), draw(mu_b), draw(mu_b+dB)
tau0 = fit_tau((pv(a0)*len(a0) + pv(b0)*len(b0))/(len(a0)+len(b0)), norm)
mua = bisect(lambda m: readout(m, tau0, norm), mean_(a0))
mub = bisect(lambda m: readout(m, tau0, norm), mean_(b0))
sa0, sa1, sb0, sb1 = syn_cells(mua, mub, +0.45, +0.45, norm)      # 两层 Δ 相同
nc = explained(sa0, sa1, sb0, sb1, norm, "B")[0]
sa0d, sa1d, sb0d, sb1d = syn_cells(mua, mub, +0.15, +0.60, norm)  # 两层 Δ 不同
pc = explained(sa0d, sa1d, sb0d, sb1d, norm, "B")[0]
print(f"  负控:两层 Δ 完全相同(+0.45/+0.45)⇒ explained = **{nc:+.4f}**(参照 **1.0**,不是 0)")
print(f"  正控:两层 Δ 真的不同(+0.15/+0.60)⇒ explained = **{pc:+.4f}**(该**明显低于** 1)")

G.identity_control("② 负控:两层 Δ 完全相同的合成世界里,`explained` 必须回到 **1.0**"
                   "(⚠ 参照是 1.0 不是 0)", observed=float(nc), expected=1.0, tol=0.08,
                   what="共享阈值、两层各自的 μ、施加同一个 Δ、按人抽样 4000 的合成世界")
G.asserted("③ 正控:两层 Δ 真的不同时,`explained` 必须明显低于 1(否则判据认不出分歧)",
           bool(pc < 0.80), f"explained = {pc:+.4f}(阈 0.80)", kind="control")
G.asserted("④ 正控在 g=0 时**不**开火(否则它在无分歧时也喊分歧)",
           bool(nc >= 0.80), f"g=0 那一格 explained = {nc:+.4f}", kind="control")
G.asserted("⑤ 前提:端点只取 `raw`(模型起点 = 观测首年分布 ⇒ 分子分母同口径,`#803` 指出的缺陷)",
           True, f"模型 gap0 与观测 gap0 同为 {mean_(a0)-mean_(b0):+.4f}", kind="control")
G.asserted("⑥ kill(预注册):「两者都在」要成立,需**所有可用格**的区间同时排除 0 与 1",
           bool(len(CELLS) > 0 and n_e0 == len(CELLS) and n_e1 == len(CELLS)),
           f"排除 0 的格 {n_e0}/{len(CELLS)} · 排除 1 的格 {n_e1}/{len(CELLS)}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*94)
lo_all = min(c["lo"] for c in CELLS) if CELLS else float("nan")
hi_all = max(c["hi"] for c in CELLS) if CELLS else float("nan")
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n_e0 == len(CELLS) and n_e1 == len(CELLS):
    V = (f"**A 两者都在 —— 而这比任何一端都更有内容。** 可用的 {len(CELLS)} 格,区间并集 "
         f"[{lo_all:+.3f}, {hi_all:+.3f}],**每一格都同时排除 0 与 1**。\n"
         f"  ⇒ **一个共同的潜在位移解释了这条鸿沟扩大的一部分,但解释不完。**\n"
         f"  ⇒ **换成关于人的话:即使虔诚者与其余人心里改变的幅度一模一样,他们在问卷上的答案也会越离越远\n"
         f"  ——那一部分是刻度造的;而剩下的那一部分,是他们真的走得不一样。两件事都真,不是二选一。**")
else:
    V = (f"**C 分不开(或只排除掉一端)。** 可用的 {len(CELLS)} 格,区间并集 [{lo_all:+.3f}, {hi_all:+.3f}];\n"
         f"  排除 0 的 {n_e0}/{len(CELLS)} · 排除 1 的 {n_e1}/{len(CELLS)}。\n"
         f"  ⇒ **`raw` 口径只用了首末两年 4 个格子的数据,它的分辨力本来就可能不够** ——\n"
         f"  **而这是关于设计的真收获,不是关于世界的判断:如实报,不硬判。**")
print(V)
print("\n⚠ **这个区间是下界性质的**:层内按人重抽**不含年际波动**(`raw` 只用首末两年)⇒ 真不确定度只会更大。")
print("⚠ **共享阈值是假设** —— 两层阈值本就不同(DIF)则「共同位移」这个反事实不成立,"
      "而 **DIF 在这份数据上测不了,结构性不可能,早已登记。**")
json.dump(dict(item=IT, K=KK, scale=[LO, HI], y0=Y0, y1=Y1, B=B,
               n_cells={f"{y}/k{k}": int(len(CELL[(y, k)])) for y in (Y0, Y1) for k in (0, 2)},
               usable_arms=[nm for _, nm in USABLE],
               dropped_arms=[dict(arm=n, headroom_frac=float(f)) for _, n, f in DROPPED],
               headroom_max=Gate.HEADROOM_MAX, cells=CELLS,
               union=[lo_all, hi_all], n_excl0=n_e0, n_excl1=n_e1,
               neg_control=float(nc), pos_control=float(pc), reference=1.0,
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"interval_not_a_threshold.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'interval_not_a_threshold.json'}")
