"""#794 · E03·A47·R233 —— 「最虔诚的三分之一」是**逐年**定义的,而美国在世俗化

**⚠⚠ 一条从 `#776` 起就在每一轮里、而账本里一次都没被命名过的设计选择:**
`REL["b"] = REL.groupby("year")["REL"].transform(qcut(3))` —— **三分位是逐年算的。**
⇒ 1974 年的「最虔诚三分之一」与 2024 年的「最虔诚三分之一」**不是同一个绝对虔诚度**。
   五十年里美国在世俗化 ⇒ **那条线自己在往下漂。**
⇒ **「虔诚者改了主意」也许有一部分根本不是改主意,是「谁算虔诚者」变了。**
⚠ 而这一条**不在页面的任何一处限定里**,也不在 `#779` 那条世代限定里 ——
   世代钉住的是出生年,**钉不住分层线。**

**⚠ 先说清楚已经有的那一半,免得把已知当新发现**(`#778` 报过两种分层):
   **`(a) attend 三档` 是一个绝对门槛**(`attend ≤1 / 2–5 / 6–8`),**`(b) 三题三分位` 是逐年相对的**,
   两者给 0.431 与 0.409,只差 5%。**⇒ 这看起来已经控制住了。**
   ⚠ **但那是一个被混淆的对照**:(a) 与 (b) **同时差了两件事** ——
   「绝对 vs 相对」**以及**「一道题 vs 三题合成」。**⇒ 分不开,所以它不构成对这条混淆的控制。**

G1 估计量(两个,方法之前先命名):
   (a) **那条线漂了多远** —— 逐年第 67 百分位的 `REL` 值,以**合并样本 SD** 为单位的漂移量。
       ⚠ `REL` 是在**合并样本**上 z 化的(不是逐年),所以它是一把绝对的尺,这个量有定义。
   (b) **把线钉死之后,`#791` 的三堆还在不在** —— 用**首年那条线的绝对值**重切,重算八题的比值。

⚠⚠ 两个世界,而第二个我会非常不想要(它会削掉 `#776` 之后的一大片):
   A **改主意**:线几乎没漂,或钉死之后比值与三堆原样存活 ⇒ 那些人真的改了想法。
   B **成分**:线漂得多,且钉死之后比值移动超过它们自己的自助宽度
     ⇒ **「虔诚者变了」里有一块是「谁算虔诚者」变了**,页面上从 `#776` 到 `#791` 的每一条都要挂限定。

预测矩阵:
   | 世界 | 现在 | 若漂 <0.25 SD 且三堆存活 | 若漂 ≥0.25 SD 且比值移动超宽度 | 若漂大但比值不动 |
   | A 改主意 | 0.5 | **0.90** | 0.05 | 0.60 |
   | B 成分   | 0.5 | 0.10 | **0.90** | 0.40 |
   ⚠ 第四列是**真实的第三种结果**:线漂了而比值不动 —— 那说明**这个估计量对分层线不敏感**,
     是一个比 A 更强的结论(不是「线没漂所以没事」,是「线漂了也没事」)。

预注册判词(条件式):
  if 正控开火(把「绝对门槛」逐年重算 ⇒ 必须**逐字复现**逐年三分位的那一版)
     and 每年绝对门槛下两层的 n 都够(≥120,否则不是功效问题是没有样本):
      漂移 < 0.25 SD  -> A(线本来就没怎么动)
      漂移 >= 0.25 SD 且 ≥半数题的比值移动 > 它自己的自助宽度 -> B
      漂移 >= 0.25 SD 且 比值移动 <= 宽度 -> **A+**(线漂了,而估计量不敏感 —— 更强)
  else: UNVERIFIED
⚠ 0.25 SD 的理由写在跑之前:`REL` 是三题 z 均值,**四分之一个合并 SD 是一个人从「每月去」
  掉到「每年去」量级的移动** —— 小于它,「同一批人」这个说法还站得住。

⚠ **「这个零该不该是零?」** —— 正控那个零**该是零**:把绝对门槛**逐年重算**,
  它在构造上**就是**逐年三分位 ⇒ 两者之差的期望**恰好是 0** ⇒ `negative_control` 不是 `offset_control`。
  ⚠ 而这条正控**能失败**:若我的绝对切法实现有 off-by-one 或边界处理不同,差就不是 0。

⚠ 跑之前写下的最强混淆:**绝对门槛下,虔诚层会随时间萎缩** ⇒ 后期年份 n 变小 ⇒
  斜率变噪 ⇒ **比值移动可能只是噪声变大,不是位置变了。**
  ⇒ 控制:**逐年报两层的 n**,并且**比值的移动必须与它自己的自助宽度比**,不与 0 比。

本轮换不了仪器(对象是世界,而第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(233)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
PREV = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS = PREV["items"]; CLUMP = PREV["clumps"]; OBS0 = PREV["obs"]
print(f"=== ⓪ 对象:`#791` 的 8 题与它们的三堆 —— {ITEMS} ===")

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
KMAX = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= KMAX[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
# ⚠ 硬规则①:三道分层题各自的 n 与年份,先打印再用
for c in ("attend", "reliten", "fund"):
    v = REL[c]; yrs = sorted(REL.year[v.notna()].unique())
    print(f"  {c:9s} n={int(v.notna().sum()):>7,} · 年 {len(yrs)} ({int(min(yrs))}–{int(max(yrs))})")
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)   # ⚠ 合并样本 z 化 ⇒ 绝对尺
print(f"  合并 REL:均值 {REL.REL.mean():+.4f} · SD {REL.REL.std(ddof=1):.4f}(**这把尺是绝对的**)")

# ── ① 那条线漂了多远 ─────────────────────────────────────────────────────────
print("\n=== ① 逐年第 67 百分位的 `REL`(那条分层线自己)===")
SD = float(REL.REL.std(ddof=1))
byyear = REL.groupby("year")["REL"]
cut67 = byyear.quantile(2/3).dropna()
n_year = byyear.size()
keep = [int(y) for y in cut67.index if n_year[y] >= 300]
c67 = cut67.loc[keep]
drift = float(c67.iloc[-1] - c67.iloc[0]); drift_sd = drift/SD
rng_sd = float(c67.max()-c67.min())/SD
print(f"  {keep[0]} 年:{c67.iloc[0]:+.4f}  →  {keep[-1]} 年:{c67.iloc[-1]:+.4f}")
print(f"  **首末漂移 {drift:+.4f} = {drift_sd:+.3f} 个合并 SD** · 全期极差 {rng_sd:.3f} SD · 年数 {len(keep)}")
slope_c = float(np.cov(np.array(keep, float), c67.to_numpy(), ddof=1)[0, 1]/np.var(np.array(keep, float), ddof=1))
print(f"  线性趋势 {slope_c:+.5f}/年 = {slope_c*(keep[-1]-keep[0])/SD:+.3f} SD / 50 年")

# ── ② 把线钉死,重算八题的比值 ────────────────────────────────────────────────
def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def strata(mode):
    if mode == "相对(逐年三分位)":
        return REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))
    if mode == "绝对(钉在首年的两条线)":
        lo, hi = REL[REL.year == keep[0]]["REL"].quantile([1/3, 2/3])
        return pd.Series(np.where(REL.REL > hi, 2, np.where(REL.REL <= lo, 0, 1)), index=REL.index)
    if mode == "正控:绝对但逐年重算":     # 构造上等于逐年三分位 ⇒ 差必须恰为 0
        out = pd.Series(index=REL.index, dtype=float)
        for y, g in REL.groupby("year"):
            lo, hi = g["REL"].quantile([1/3, 2/3])
            # ⚠⚠ 第一版写成 `>= hi` —— 而 `pd.qcut` 的分箱是**右闭** `(a, b]`,顶箱是 `> q67`。
            #    `REL` 是三道序数题的 z 均值 ⇒ **边界上有大量并列**,于是并列的人在两种实现里进了不同的箱,
            #    正控给出 |差| = 4.94e-02 而不是 0。**控制开火了,而开火的是被控的那个东西。**
            out.loc[g.index] = np.where(g.REL > hi, 2, np.where(g.REL <= lo, 0, 1))
        return out
    raise ValueError(mode)

def ratios(col, nmin=120):
    out, ns = {}, {}
    for c in ITEMS:
        g = REL.dropna(subset=[c])
        s = {}
        for k in (2, 0):
            rows = [(int(y), float(gy[c].mean()), len(gy)) for y, gy in g[g[col] == k].groupby("year") if len(gy) >= nmin]
            s[k] = rows
        ns[c] = (len(s[2]), len(s[0]))
        if len(s[2]) < 8 or len(s[0]) < 8: out[c] = np.nan; continue
        sA = slope([r[0] for r in s[2]], [r[1] for r in s[2]])
        sB = slope([r[0] for r in s[0]], [r[1] for r in s[0]])
        out[c] = sA/sB if abs(sB) > 1e-12 else np.nan
    return out, ns

print("\n=== ② 三种分层,同一批题(`G4`:分层定义是一个规格轴)===")
RES = {}
for mode in ("相对(逐年三分位)", "正控:绝对但逐年重算", "绝对(钉在首年的两条线)"):
    REL["k"] = strata(mode)
    r, ns = ratios("k")
    RES[mode] = dict(r=r, ns=ns)
    print(f"  {mode}")
    print("    " + " · ".join(f"{c} {r[c]:+.3f}" if np.isfinite(r[c]) else f"{c} n/a" for c in ITEMS))
    print(f"    每题(虔诚年数/非虔诚年数):" + " · ".join(f"{c} {ns[c][0]}/{ns[c][1]}" for c in ITEMS))

# 正控:逐年重算的绝对切法必须逐字复现逐年三分位
pc_pairs = [(RES["相对(逐年三分位)"]["r"][c], RES["正控:绝对但逐年重算"]["r"][c])
            for c in ITEMS if np.isfinite(RES["相对(逐年三分位)"]["r"][c])]
pc_max = max(abs(a-b) for a, b in pc_pairs) if pc_pairs else np.inf
print(f"\n  正控:逐年重算的绝对切法 vs 逐年三分位 —— 最大 |差| = **{pc_max:.3e}** "
      f"(**该恰为 0**:两者在构造上是同一个切法)")

# ── ③ 钉死之后,移动了多少 —— 与它自己的自助宽度比,不与 0 比 ────────────────────
print("\n=== ③ 钉死之后每题移动了多少(⚠ 与它**自己的自助宽度**比,不与 0 比)===")
REL["k"] = strata("绝对(钉在首年的两条线)")
def boot_width(c, col="k", B=1500):
    g = REL.dropna(subset=[c])
    s = {k: [(int(y), float(gy[c].mean())) for y, gy in g[g[col] == k].groupby("year") if len(gy) >= 120]
         for k in (2, 0)}
    if len(s[2]) < 8 or len(s[0]) < 8: return np.nan
    yA = np.array([r[0] for r in s[2]], float); vA = np.array([r[1] for r in s[2]])
    yB = np.array([r[0] for r in s[0]], float); vB = np.array([r[1] for r in s[0]])
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(B)])
    bs = bs[np.isfinite(bs)]
    return float(np.percentile(bs, 97.5)-np.percentile(bs, 2.5))
rows, moved = [], 0
for c in ITEMS:
    r0, r1 = OBS0[c], RES["绝对(钉在首年的两条线)"]["r"][c]
    w = boot_width(c)
    if not np.isfinite(r1) or not np.isfinite(w):
        rows.append(dict(item=c, clump=CLUMP[c], r_rel=r0, r_abs=None, width=None, moved=None));
        print(f"  {c:9s} {CLUMP[c]}  {r0:+.3f} → **算不出**(绝对门槛下年数不足)"); continue
    mv = abs(r1-r0) > w
    moved += int(mv)
    rows.append(dict(item=c, clump=CLUMP[c], r_rel=float(r0), r_abs=float(r1), width=float(w), moved=bool(mv)))
    print(f"  {c:9s} {CLUMP[c]}  {r0:+.3f} → {r1:+.3f}  |移动| {abs(r1-r0):.3f} vs 自助宽 {w:.3f}"
          f"  {'**移动超过宽度**' if mv else '在宽度内'}")
ok_rows = [x for x in rows if x["r_abs"] is not None]
print(f"\n  可算 {len(ok_rows)}/{len(ITEMS)} 题 · **移动超过自己宽度的 {moved}/{len(ok_rows)}**")

G = Gate("#794 · 那条分层线自己漂了多远,钉死之后三堆还在不在")
# ⚠⚠ 库缺陷(本轮当场发现,见账本):`_control_rows` **只由 `asserted(kind="control")` 填充**,
#    而 `negative_control`/`positive_control` 等具名控制方法**不进那个集合** ⇒
#    `three_valued()` 把一条**控制**的失败读成 **kill 开火**,打印 `OVERTURNED`。
#    **一个假的 OVERTURNED 就是一次假撤回,而 `#782` 已经立过:假撤回与假无罪一样永久。**
#    ⇒ 本轮的总判**由脚本自己的 `ctrl` 算**,不采信库的那一行;库的修法与爆炸半径记在账本里。
G.negative_control("① 正控/负控合一:绝对切法**逐年重算**必须逐字复现逐年三分位(构造上同一个切法 ⇒ 差恰为 0)",
                   null=float(pc_max), effect=float(np.nanmax(np.abs(list(OBS0.values())))), ratio=1e-6,
                   null_kind="同一个切法的两种实现之差(期望恰为 0,非偏移)")
G.asserted("② 前提:绝对门槛下必须至少 6 题算得出来(否则不是功效问题,是没有样本)",
           bool(len(ok_rows) >= 6), f"可算 {len(ok_rows)}/{len(ITEMS)}", kind="control")
big_drift = bool(abs(drift_sd) >= 0.25)
G.asserted("③ 观测(非判据):那条线首末漂移是否 ≥ 0.25 个合并 SD",
           True, f"漂移 {drift_sd:+.3f} SD · 极差 {rng_sd:.3f} SD ⇒ {'大' if big_drift else '小'}", kind="control")
half = bool(ok_rows and moved > len(ok_rows)/2)
G.asserted("④ kill(预注册):「成分说」要成立,需**过半数**题的比值移动超过它自己的自助宽度",
           half, f"移动超宽度 {moved}/{len(ok_rows)}", kind="kill")
print(); print(G)

print("\n"+"="*92)
ctrl = bool(pc_max < 1e-9 and len(ok_rows) >= 6)
if not ctrl:
    v = f"**UNVERIFIED:正控 |差|={pc_max:.2e}(需 <1e-9)或可算题 {len(ok_rows)}<6 ⇒ 本轮不下判。**"
elif half:
    v = (f"**B 成分:「虔诚者变了」里有一块是「谁算虔诚者」变了。** 那条线五十年漂了 **{drift_sd:+.3f} 个 SD**,"
         f"钉死之后 **{moved}/{len(ok_rows)}** 题的比值移动超过它自己的自助宽度 ⇒ "
         f"**从 `#776` 到 `#791` 的每一条都要挂上这个限定。**")
elif big_drift:
    v = (f"**A+ 线漂了,而这个估计量对它不敏感 —— 这比「线没漂」更强。** 那条分层线五十年漂了 "
         f"**{drift_sd:+.3f} 个合并 SD**(极差 {rng_sd:.3f} SD),**确实在动**;"
         f"而把它钉死在首年之后,**{len(ok_rows)-moved}/{len(ok_rows)}** 题的比值移动**没有**超过它自己的自助宽度。\n"
         f"  ⇒ **「最虔诚的三分之一」这个说法里的成分漂移是真的,但它没有制造那些比值。**\n"
         f"  ⇒ 三堆:" + " · ".join(f"{x['item']}({x['clump']}) {x['r_rel']:+.2f}→{x['r_abs']:+.2f}"
                                    for x in ok_rows))
else:
    v = (f"**A 线本来就没怎么动。** 首末漂移只有 **{drift_sd:+.3f} 个 SD** < 0.25 ⇒ "
         f"「最虔诚的三分之一」在五十年里大体是同一个绝对虔诚度,成分说不成立。")
print(v)
json.dump(dict(drift=drift, drift_sd=drift_sd, range_sd=rng_sd, slope_per_year=slope_c,
               years=keep, cut67={int(y): float(v2) for y, v2 in c67.items()},
               positive_control_maxdiff=float(pc_max), rows=rows, n_moved=moved,
               n_computable=len(ok_rows), verdict=v, gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"boundary_drift.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'boundary_drift.json'}")
