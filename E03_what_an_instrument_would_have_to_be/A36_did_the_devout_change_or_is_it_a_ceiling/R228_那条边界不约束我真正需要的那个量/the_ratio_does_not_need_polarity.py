"""#789 · E03·A45·R228 —— `#788` 记下的那条边界,不约束我真正需要的那个量

`#788` 判 UNVERIFIED,并把原因记成一条**能力边界**:
   「这具自动定向仪器覆盖 Likert 式道德题,**结构性地覆盖不了 yes/no 政策题**;
    而『一个总闸还是一个个开关』恰恰需要非性的政策题,所以这个问题是被我的仪器挡住的。」

**⚠⚠ 而那句话是错的,并且错得可以用三行代数看出来 —— 这一轮先杀它。**

我要的量是 **`r = slope_虔诚 / slope_非虔诚`**。把任一题的编码翻向(`v → K+1−v`):
   `slope` 是 `cov(year, v)/var(year)`,而 `cov(year, K+1−v) = −cov(year, v)`
   ⇒ **两层的斜率同时变号 ⇒ 比值不变。**
⇒ **`r` 对题目的极性是不变的。我从来不需要知道哪一端是「宽容」。**

⚠ 这正是 `frontier §3` 里那条最便宜的攻击(gauge test):
   **说出让行为完全相同的那些变换,再问「测量」在它们之下变不变。**
   这里测量不变,**而我上一轮却因为一个只影响「怎么描述」的量而停了整轮。**
⇒ **`#788` 的能力边界要缩**:它对**依赖极性的量**(均值方向、端点占比是哪一端)成立,
   **对本轮的比值不成立。** 一条能力边界被写得比它实际的范围大,**和一个被写得太大的结论一样危险** ——
   而它更隐蔽,因为**一条「做不到」从来不会有人回去验它**(`realstat`:「一堵从没被查过的墙」)。

G1 估计量:与 `#784`–`#787` 逐字同一个 —— 逐题的「虔诚/非虔诚」年代斜率比,带年份自助区间。
   **本轮唯一的变化:题目集合扩到非性政策题,而极性一栏被删掉了,因为它不该在那里。**

⚠ 极性仍然出现在**一个**地方,而只出现在那里:**把结果翻译成一句话时**。
   `r > 0` = 两层同向动;`r < 0` = 反向。**方向的名字只用于叙述,不进任何计算。**
   ⇒ 产物里把 `题干原文`(从 `.dta` 的变量标签读出,不是我记的)与 `r` 并排存,
   **让读的人自己判方向,而不是让我替他判。**

预注册判词(条件式):
  if gauge 正控开火(翻向后 `r` 逐格不变,容差 1e-9)and 可读题 >= 6:
      设 gap = 非性题 `|r|` 中位 − 性题 `|r|` 中位
      if |gap| <= 0.165 -> **A 一个总闸**(0.165 = `#786` 同题跨三窗口的跨度 = 本设计的分辨极限)
      elif gap > 0.165  -> **B 一个个开关**:性题的比更小,刹车是教义性的
      else              -> **反号**:非性题上的比反而更小 —— 那会把「性道德特殊」整个翻过来
  else: UNVERIFIED
⚠ 用 `|r|`:因为跨题比较的是**幅度**,而 `r<0` 的题(两层反向动)属于另一种现象,单独列出不并入。

⚠ 跑之前写下的最强混淆(与 `#788` 同):非性题可能整体没怎么动 ⇒ 分母小 ⇒ 比值不稳。
  ⇒ 控制不变:同一轮报每题**非虔诚层自己的斜率与它自己的置换零**,把「比值不稳」与「社会没动」分开。
⚠ `#787` 的规矩:每一个比例都同时报「可读 / 尝试」两个数。

本轮换不了仪器,理由同 `R223/instrument_search.py`(对象是世界,第二具仪器本机六具全部落选)。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(228)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"

SEX = ["homosex", "xmarsex", "premarsx", "teensex", "pornlaw", "sexeduc", "divlaw"]
NONSEX = ["abany", "abdefect", "abnomore", "abpoor", "suicide1", "suicide2", "grass",
          "cappun", "gunlaw", "letdie1", "spanking", "prayer", "racmar", "natfare",
          "helpblk", "fechld", "fepresch"]
CAND = SEX + NONSEX
GRP = {c: ("性" if c in SEX else "非性") for c in CAND}

# 题干原文从 `.dta` 的变量标签读出 —— **不是我记的**(硬规则①)
STEM = pd.io.stata.StataReader(gp).variable_labels()

print("=== ① 普查:每列的 n · 真正被问过的年份 · 档数 · **题干原文**(硬规则①)===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+CAND, convert_categoricals=False)
cat = pd.read_stata(gp, columns=CAND, convert_categoricals=True)
usable, rejected = [], []
KMAX = {}
for c in CAND:
    v = pd.to_numeric(d[c], errors="coerce")
    cs = list(cat[c].cat.categories) if hasattr(cat[c], "cat") else []
    if c == "homosex": cs = cs[:4]          # 末档 `other` 不是作答档(`#788`①)
    KMAX[c] = len(cs)
    yrs = sorted(d.year[v.notna()].unique())
    ok = (2 <= len(cs) <= 5) and len(yrs) >= 12
    (usable if ok else rejected).append(c)
    print(f"  {c:9s} {GRP[c]:2s} n={int(v.notna().sum()):>7,} 年{len(yrs):>3} 档{len(cs):>2} "
          f"{'用' if ok else '弃'}  「{STEM.get(c,'?')}」")
print(f"\n  候选 {len(CAND)} · **可用 {len(usable)} / 尝试 {len(CAND)}** · 弃 {len(rejected)}:{rejected}")

M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= KMAX[c])) for c in usable})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
z = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = z(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["b"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def series(item, k, flip=False, nmin=120):
    g = REL.dropna(subset=[item]); out = []
    for y, gy in g[g["b"] == k].groupby("year"):
        if len(gy) < nmin: continue
        v = gy[item]
        out.append((int(y), float((KMAX[item]+1-v).mean() if flip else v.mean())))
    return out
def ratio(item, flip=False):
    rA, rB = series(item, 2, flip), series(item, 0, flip)
    if len(rA) < 10 or len(rB) < 10: return None
    yA = np.array([r[0] for r in rA], float); yB = np.array([r[0] for r in rB], float)
    vA = np.array([r[1] for r in rA]); vB = np.array([r[1] for r in rB])
    return slope(yA, vA)/slope(yB, vB), (yA, vA, yB, vB)

# ── ② gauge 正控:翻向后比值必须逐格不变 ────────────────────────────────────────
print("\n=== ② gauge 正控:把每一题的编码翻向,比值必须**逐格不变**(三行代数,现在把它变成测量)===")
gauge = []
for c in usable:
    a, b = ratio(c, False), ratio(c, True)
    if a is None or b is None: continue
    gauge.append(dict(item=c, r=a[0], r_flipped=b[0], absdiff=abs(a[0]-b[0])))
    print(f"  {c:9s} r={a[0]:+8.4f}  翻向后 r={b[0]:+8.4f}  |差|={abs(a[0]-b[0]):.3e}")
gauge_ok = bool(gauge and max(g["absdiff"] for g in gauge) < 1e-9)
print(f"  ⇒ 最大差 **{max(g['absdiff'] for g in gauge):.3e}** —— "
      f"{'**比值确实与极性无关 ⇒ `#788` 的边界不约束这个量**' if gauge_ok else '⛔ 不变性不成立,停'}")

# ── ③ 逐题 ────────────────────────────────────────────────────────────────────
print("\n=== ③ 逐题:比值 · 自助区间 · 非虔诚层自己动了没有 ===")
rows, attempted = [], 0
for c in usable:
    attempted += 1
    out = ratio(c)
    if out is None:
        rows.append(dict(item=c, grp=GRP[c], stem=STEM.get(c, ""), readable=False, why="年数不足")); continue
    r0, (yA, vA, yB, vB) = out
    sB = slope(yB, vB)
    q95 = float(np.quantile(np.abs([slope(RNG.permutation(yB), vB) for _ in range(2000)]), .95))
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(4000)])
    bs = bs[np.isfinite(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    ok = not (lo <= 1.0 <= hi)
    rows.append(dict(item=c, grp=GRP[c], stem=STEM.get(c, ""), readable=bool(ok), r=float(r0),
                     lo=lo, hi=hi, slopeB=sB, nullB95=q95, moved=bool(abs(sB) > q95),
                     nyr=[len(yA), len(yB)], why=None if ok else "比值区间含 1.0"))
    print(f"  {c:9s} {GRP[c]:2s} r={r0:+7.3f} [{lo:8.3f},{hi:8.3f}] 非虔诚斜率 {sB:+.5f} vs 零 {q95:.5f}"
          f" {'动了' if abs(sB) > q95 else '**没动**'}  {'可读' if ok else '不可读'}  「{STEM.get(c,'')[:34]}」")

R = [x for x in rows if x.get("readable")]
neg = [x for x in R if x["r"] < 0]
pos = [x for x in R if x["r"] >= 0]
sx = [abs(x["r"]) for x in pos if x["grp"] == "性"]; ns = [abs(x["r"]) for x in pos if x["grp"] == "非性"]
print(f"\n  **可读 {len(R)} / 尝试 {attempted}**(`#787`) · 同向 {len(pos)} · **反向(r<0){len(neg)}:"
      f"{[x['item'] for x in neg]}** —— 反向是另一种现象,单独列出,不并入幅度比较")
print(f"  性 {len(sx)} 题 · 非性 {len(ns)} 题")

# ── ④ ⚠⚠ 跑完之后补的一条,而它推翻了我自己的预注册判据 ────────────────────────
# 预注册拿「两组的**中位数之差**」去比 0.165,而 0.165 是**同一题跨三个窗口**的跨度 ——
# 那是**题内**的抖动。可这里比的是**题间**的中位数,而题间的散布是另一回事:
#   性 |r| = [0.364, 0.409, 2.192] · 非性 |r| = [0.105, 0.126, 0.410, 0.462, 1.841]
# **跨度 ~2.0,是那个阈值的十二倍。** ⇒ 判据用错了噪声来源,与 `#782` 同一族
# (那次是「检验了分子与分母,没检验比值」;这次是「用题内噪声去判题间差」)。
# ⇒ 补一个**对题目重抽**的自助,给中位数之差它**自己的**区间,再看它是不是分辨得了。
if len(pos) >= 6 and sx and ns:
    _sx, _ns = np.array(sx), np.array(ns)
    _bs = np.array([np.median(RNG.choice(_ns, len(_ns), replace=True))
                    - np.median(RNG.choice(_sx, len(_sx), replace=True)) for _ in range(4000)])
    gap_lo, gap_hi = float(np.percentile(_bs, 2.5)), float(np.percentile(_bs, 97.5))
    resolvable = bool(gap_hi - gap_lo <= 2*0.165)
    print(f"\n=== ④ 中位数之差**自己的**区间(对题目重抽,而不是对年份)===")
    print(f"  差 {float(np.median(_ns)-np.median(_sx)):+.3f} · **95% 区间 [{gap_lo:+.3f}, {gap_hi:+.3f}]**"
          f" 宽 {gap_hi-gap_lo:.3f}")
    print(f"  ⇒ 要分辨「差 ≤ 0.165」与「差 > 0.165」,区间宽度必须 ≤ 0.330 —— "
          f"**实测 {gap_hi-gap_lo:.3f} ⇒ {'够' if resolvable else '**不够,这个比较本设计做不了**'}**")
else:
    gap_lo = gap_hi = float("nan"); resolvable = False

G = Gate("#789 · 那条边界不约束我真正需要的那个量")
G.asserted("① gauge 正控:翻向后比值逐格不变(容差 1e-9)——**这一条同时是对 `#788` 那条边界的攻击**",
           gauge_ok, f"{len(gauge)} 题 · 最大 |差| {max(g['absdiff'] for g in gauge):.3e}", kind="control")
G.asserted("② 前提:可读且同向的题 >= 6,且性/非性两边都非空",
           bool(len(pos) >= 6 and sx and ns), f"可读 {len(R)}/{attempted} · 同向 {len(pos)} · 性 {len(sx)} · 非性 {len(ns)}",
           kind="control")
if len(pos) >= 6 and sx and ns:
    gap = float(np.median(ns)-np.median(sx))
    one_brake = bool(abs(gap) <= 0.165)
    G.asserted("③ 前提(跑完才想到,而它先于 kill):中位数之差必须先**可分辨** —— "
               "它自己的区间宽度 ≤ 2×0.165,否则 kill 无论落哪一边都不是测量",
               resolvable, f"差的 95% 区间 [{gap_lo:+.3f}, {gap_hi:+.3f}] 宽 {gap_hi-gap_lo:.3f} vs 需 ≤0.330",
               kind="control")
    G.asserted("④ kill(预注册):「一个总闸」要站住,需 |非性中位 − 性中位| ≤ 0.165"
               "(= `#786` 同题跨三窗口的跨度)",
               one_brake, f"性 {np.median(sx):.3f} vs 非性 {np.median(ns):.3f}(差 {gap:+.3f})", kind="kill")
else:
    gap = float("nan"); one_brake = False
    G.asserted("③ kill:样本不足以开火", False, f"同向 {len(pos)} · 性 {len(sx)} · 非性 {len(ns)}", kind="kill")
print(); print(G)

print("\n"+"="*92)
if gauge_ok and len(pos) >= 6 and sx and ns and not resolvable:
    v = (f"**不可分,而这推翻的是我自己的预注册判据,不是数据。** 两组的 |r| 中位确实几乎相等"
         f"(性 {np.median(sx):.3f} vs 非性 {np.median(ns):.3f},差 {gap:+.3f}),"
         f"**但那个差自己的 95% 区间是 [{gap_lo:+.3f}, {gap_hi:+.3f}],宽 {gap_hi-gap_lo:.3f}** —— "
         f"要分辨「≤0.165」与「>0.165」需要 ≤0.330。\n"
         f"  ⚠ 预注册拿 0.165 当阈值,而 0.165 是**同一题跨窗口**的抖动;这里比的是**题间**中位数,"
         f"题间散布 {min(sx+ns):.3f}–{max(sx+ns):.3f},约十二倍于它。**判据用错了噪声来源** —— "
         f"与 `#782` 同一族。\n"
         f"  ⇒ **而数据本身给出了一件比 A/B 都重要的事:变异在组内,不在组间。** "
         f"性题里 `sexeduc` 是 {max(sx):.2f}(虔诚者改得**更多**),`teensex` 是 {min(sx):.2f};"
         f"非性题里 `racmar` 是 {max(ns):.2f},`helpblk` 是 {min(ns):.2f}。"
         f"**「性 / 非性」这个切法解释不了任何东西 —— 我切错了维度。**")
elif not (gauge_ok and len(pos) >= 6 and sx and ns):
    v = f"**UNVERIFIED:gauge 正控或前提没过(可读 {len(R)}/{attempted},同向 {len(pos)},性 {len(sx)}、非性 {len(ns)})。**"
elif one_brake:
    v = (f"**A 一个总闸:那个刹车不挑题。** 可读 {len(R)}/{attempted},同向 {len(pos)} 题;"
         f"性题 |r| 中位 **{np.median(sx):.3f}**、非性题 **{np.median(ns):.3f}**,差 {gap:+.3f} "
         f"在本设计的分辨极限 0.165 之内 ⇒ **虔诚是一种对任何道德变迁的普遍阻力,性道德不特殊。**")
elif gap > 0.165:
    v = (f"**B 一个个开关:刹车是教义性的。** 性题 |r| 中位 **{np.median(sx):.3f}**、"
         f"非性题 **{np.median(ns):.3f}**,**差 {gap:+.3f} 超过分辨极限 0.165** ⇒ "
         f"**虔诚者不是对一切变迁都踩刹车,他们只在传统点过名的地方踩。**")
else:
    v = (f"**反号 —— 而这会把「性道德特殊」整个翻过来。** 非性题 |r| 中位 {np.median(ns):.3f} "
         f"**小于**性题的 {np.median(sx):.3f}(差 {gap:+.3f})⇒ 虔诚者在非性的道德题上刹得**更死**。")
print(v)
json.dump(dict(rows=rows, gauge=gauge, gauge_ok=gauge_ok, usable=usable, rejected=rejected,
               n_readable=len(R), n_attempted=attempted, negatives=[x["item"] for x in neg],
               sex_absr=sx, nonsex_absr=ns, gap=gap, verdict=v,
               gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"ratio_needs_no_polarity.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'ratio_needs_no_polarity.json'}")
