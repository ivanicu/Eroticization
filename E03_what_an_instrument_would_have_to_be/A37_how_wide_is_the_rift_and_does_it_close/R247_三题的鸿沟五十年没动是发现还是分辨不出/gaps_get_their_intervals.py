"""#808 · E03·A58·R247 —— 「五裂三合」一直只是方向;而三题的鸿沟五十年没动,是发现还是分辨不出?

**两笔账在这里合流,而它们是同一件事。**
`#802`①(**页面上最大的一笔欠账**):「给这八条差距挂上区间 —— 在挂上之前,
『五裂三合』只是方向,不许写成『显著』。」**它欠了六轮,而页面一直在用那句话。**
`#806`③ / `#807`②:`teensex`·`suicide2`·`prayer` 三题的 `explained` 区间跨到 ±7 ~ ±31,
**根因是它们的分母 `gap1 − gap0` 本来就小** ——
**而「这三题上宗教鸿沟五十年几乎没动」本身是一句关于人的话,页面把它混在「五裂三合」里报了。**

⚠⚠ **而这两笔合起来指向一个我一直没问的问题:**
   **「差距没变」到底是一个发现,还是这个设计分辨不出?**
   `#800`/`零的分辨率`那条规矩:**永远不要解释一个落在分辨率地板以下的差。**
   ⇒ **在给出区间之前,「缩 3/8、扩 5/8」这句话里可能有几条根本不该被数进任何一边。**

G1 估计量:**逐题的 `Δgap = gap1 − gap0`,带它自己的区间**,以及三值判定:
   **动了**(区间排除 0)· **没动**(区间含 0,**而且窄** ⇒ 真的没动)· **分辨不出**(区间含 0 **且宽**)。
⚠⚠ **后两者是完全不同的两句话,而「缩/扩」的二分把它们全压成了一个方向** ——
   **这正是 `#802` 那张表的结构性缺陷:它没有「没动」这一档,所以每一题都被迫站队。**

`G4` 两个规格,而**它们的重抽方案必须不同,这一条要写明**(`#797` 指出页面上并存两套不可比自助):
   ① **`fit`**(用上全部年份的线性趋势)⇒ **按年份聚类重抽**(年是抽样单位);
   ② **`raw`**(首末年观测均值)⇒ **在四个 `年 × 层` 格内按人重抽**。
   **两套不是同一个东西,不合并、不平均,各报各的,并在结果里标明是哪一套。**

「宽」的判据必须在跑之前定死,且**不能是我看了结果之后挑的**:
   **宽 = 区间宽度 > |最大那题的 Δgap|** —— 即「这条区间宽到能装下整张表里最大的一次移动」。
   ⚠ 这是一个**相对于本设计自身分辨力**的定义,不是相对于零。

三个世界:
   A **八题都动了**:全部区间排除 0 ⇒ 「五裂三合」站得住,只是现在带上了区间。
   B **有几题根本没动**:窄且含 0 ⇒ **「宗教鸿沟在这几道题上五十年没动」是一句新的、关于人的话。**
   C **有几题分辨不出**:宽且含 0 ⇒ **它们不该被数进「五裂三合」的任何一边**,
     而那会**直接改掉 `#802` 那句已发表的结论的分母。**

预测矩阵:
   | 世界 | 现在 | 全排除 0 | 有窄且含 0 | 有宽且含 0 |
   | A 都动了     | 0.35 | **0.90** | 0.05 | 0.05 |
   | B 有的没动   | 0.35 | 0.05 | **0.80** | 0.15 |
   | C 有的分辨不出 | 0.30 | 0.05 | 0.15 | **0.85** |

预注册判词(条件式):
  if 正控开火(**合成一个已知大小的 Δgap,两个规格都要取回来,且区间必须排除 0**)
     and 负控开火(**两层平行、差距按构造恒定的世界,Δgap 的区间必须含 0**):
      逐题三值判定 + 整张网格全报;**总判按计数,不设「多数」阈值**(`#805` 的教训)
  else: UNVERIFIED
⚠ **「这个零该不该是零?」** —— 负控问的是「差距恒定时 Δgap 会不会被判成动了」,
  **参照值是 0,而这一次它真的该是 0**(与 `#801`/`#805` 参照 1.0 的情形相反,所以要明确说)。

⚠ 跑之前写下的最强混淆:**八题的年份跨度不同(28–50 年)且档数不同(2–5)** ——
  一条 28 年的题即使每年动得一样多,总位移也更小。
  ⇒ 控制:**同时报「Δgap ÷ 年数 ÷ 量表跨度」**,而**结论只用三值判定,不用大小跨题排序。**

⚠ 硬规则①:先打印每题 n、真正被问过的年份、档数。本轮换不了仪器(`R223` 六具全部落选)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(247)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS, OBS = P791["items"], P791["obs"]
B = 3000

print("=== ⓪ 硬规则①:每题 n · 真正被问过的年份 · 档数 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c in ITEMS:
    yrs = sorted(d.year[M[c].notna()].unique())
    print(f"  {c:9s} n={int(M[c].notna().sum()):>7,} · 年 {len(yrs):>2}({int(min(yrs))}–{int(max(yrs))}) · 档 {K[c]}")
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))

def panel(item):
    g = REL.dropna(subset=[item])
    A = {int(y): gy[item].to_numpy(float) for y, gy in g[g.k == 2].groupby("year") if len(gy) >= 120}
    Bp = {int(y): gy[item].to_numpy(float) for y, gy in g[g.k == 0].groupby("year") if len(gy) >= 120}
    yy = sorted(set(A) & set(Bp))
    return (yy, A, Bp) if len(yy) >= 8 else None

def dgap_fit(yy, A, Bp):
    x = np.array(yy, float)
    out = []
    for D in (A, Bp):
        y = np.array([D[t].mean() for t in yy])
        b = np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1); a = y.mean()-b*x.mean()
        out.append((a+b*x[0], a+b*x[-1]))
    return (out[0][1]-out[1][1]) - (out[0][0]-out[1][0])
def dgap_raw(yy, A, Bp):
    y0, y1 = yy[0], yy[-1]
    return (A[y1].mean()-Bp[y1].mean()) - (A[y0].mean()-Bp[y0].mean())

def boot_fit(yy, A, Bp, rng):
    idx = rng.integers(0, len(yy), len(yy))            # ⚠ 年份聚类重抽(年是抽样单位)
    ys = sorted({yy[i] for i in idx})
    return dgap_fit(ys, A, Bp) if len(ys) >= 5 else np.nan
def boot_raw(yy, A, Bp, rng):
    y0, y1 = yy[0], yy[-1]                              # ⚠ 四个 `年 × 层` 格内按人重抽
    r = lambda a: a[rng.integers(0, len(a), len(a))]
    return (r(A[y1]).mean()-r(Bp[y1]).mean()) - (r(A[y0]).mean()-r(Bp[y0]).mean())

print(f"\n=== ① 逐题 Δgap 与它自己的区间(B={B};⚠ **两个规格用不同的重抽方案,不合并**)===")
ROWS = []
for it in ITEMS:
    p = panel(it)
    if p is None: continue
    yy, A, Bp = p
    for spec, pt_f, bo_f, scheme in (("fit", dgap_fit, boot_fit, "年份聚类重抽"),
                                     ("raw", dgap_raw, boot_raw, "格内按人重抽")):
        pt = pt_f(yy, A, Bp)
        dr = np.array([bo_f(yy, A, Bp, RNG) for _ in range(B)])
        dr = dr[np.isfinite(dr)]
        lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
        ROWS.append(dict(item=it, K=K[it], spec=spec, scheme=scheme, y0=yy[0], y1=yy[-1],
                         n_years=len(yy), dgap=float(pt), lo=lo, hi=hi, width=hi-lo,
                         per_year_per_span=float(pt)/(yy[-1]-yy[0])/(K[it]-1)))
MAXD = max(abs(r["dgap"]) for r in ROWS)
for r in ROWS:
    r["excl0"] = bool(r["lo"] > 0 or r["hi"] < 0)
    r["wide"] = bool(r["width"] > MAXD)          # ⚠ 判据在跑前定死:宽到能装下整张表最大的一次移动
    r["verdict"] = "动了" if r["excl0"] else ("**分辨不出**" if r["wide"] else "**没动**")
print(f"  ⚠ 「宽」的门槛 = 整张表里最大的 |Δgap| = **{MAXD:.3f}**(跑前定死,不是看了结果挑的)\n")
print(f"  {'题':9s} {'档':>2} {'规格':4s} {'重抽':6s} {'Δgap':>8s} {'95% 区间':>20s} {'宽':>7s}  判定")
for r in ROWS:
    print(f"  {r['item']:9s} {r['K']:>2} {r['spec']:4s} {r['scheme']:6s} {r['dgap']:>+8.3f} "
          f"[{r['lo']:+.3f}, {r['hi']:+.3f}] {r['width']:>7.3f}  **{r['verdict']}**")

by = {}
for r in ROWS: by.setdefault(r["item"], []).append(r)
agree = {it: (rs[0]["verdict"] == rs[1]["verdict"]) for it, rs in by.items() if len(rs) == 2}
n_move = sum(1 for it, rs in by.items() if all(x["excl0"] for x in rs))
n_still = sum(1 for it, rs in by.items() if all(not x["excl0"] and not x["wide"] for x in rs))
n_unres = sum(1 for it, rs in by.items() if any(not x["excl0"] and x["wide"] for x in rs))
print(f"\n  **两个规格判定一致的题 {sum(agree.values())}/{len(agree)}** · "
      f"两规格都判「动了」**{n_move}** · 都判「没动」**{n_still}** · 有一格判「分辨不出」**{n_unres}**")

print("\n=== ② 控制(合成世界,同一条代码路径)===")
yy_s = list(range(1974, 2025, 2))
def syn(g0, g1, n=400):
    A = {y: RNG.normal(2.0 + (g0 + (g1-g0)*(y-1974)/50), 0.9, n) for y in yy_s}
    Bp = {y: RNG.normal(2.0, 0.9, n) for y in yy_s}
    return yy_s, A, Bp
ny, nA, nB = syn(0.30, 0.30)                       # 负控:差距按构造恒定 ⇒ Δgap 该含 0
nc_pt = dgap_fit(ny, nA, nB)
nc_dr = np.array([boot_fit(ny, nA, nB, RNG) for _ in range(1500)]); nc_dr = nc_dr[np.isfinite(nc_dr)]
nc_lo, nc_hi = np.percentile(nc_dr, [2.5, 97.5])
py, pA, pB = syn(0.00, 0.50)                       # 正控:已知 Δgap = +0.50
pc_pt = dgap_fit(py, pA, pB)
pc_dr = np.array([boot_fit(py, pA, pB, RNG) for _ in range(1500)]); pc_dr = pc_dr[np.isfinite(pc_dr)]
pc_lo, pc_hi = np.percentile(pc_dr, [2.5, 97.5])
print(f"  负控:差距按构造恒定(+0.30 全程)⇒ Δgap = {nc_pt:+.4f} [{nc_lo:+.4f}, {nc_hi:+.4f}] "
      f"⇒ 含 0:**{'是' if nc_lo <= 0 <= nc_hi else '否'}**(⚠ **这一次参照真的是 0**)")
print(f"  正控:已知 Δgap = +0.500 ⇒ 取回 {pc_pt:+.4f} [{pc_lo:+.4f}, {pc_hi:+.4f}] "
      f"⇒ 排除 0:**{'是' if pc_lo > 0 else '否'}**")

G = Gate("#808 · 「五裂三合」一直只是方向;三题的鸿沟五十年没动?")
G.asserted("① 正控:合成一个已知 Δgap = +0.500 的世界,必须取回来且区间排除 0",
           bool(pc_lo > 0 and abs(pc_pt-0.50) < 0.15), f"取回 {pc_pt:+.4f} [{pc_lo:+.4f}, {pc_hi:+.4f}]",
           kind="control")
G.asserted("② 负控:差距按构造恒定的世界里,Δgap 的区间必须**含 0**"
           "(⚠ 参照真的是 0 —— 与 `#801`/`#805` 参照 1.0 的情形相反)",
           bool(nc_lo <= 0 <= nc_hi), f"Δgap = {nc_pt:+.4f} [{nc_lo:+.4f}, {nc_hi:+.4f}]", kind="control")
G.asserted("③ 正控在 g=0 时**不**开火(负控那一格不许被判成「动了」)",
           bool(not (nc_lo > 0 or nc_hi < 0)), f"负控区间 [{nc_lo:+.4f}, {nc_hi:+.4f}] 含 0", kind="control")
G.asserted("④ 前提(跑前写下的混淆):年份跨度 28–50 年、档数 2–5 不同 ⇒ 同时报「÷年数÷跨度」,"
           "且判定只用三值,不用大小跨题排序",
           bool(all("per_year_per_span" in r for r in ROWS)),
           f"年数 {sorted({r['n_years'] for r in ROWS})} · 档数 {sorted({r['K'] for r in ROWS})}", kind="control")
G.asserted("⑤ 前提:两个规格用**不同的重抽方案**,分别标明,不合并不平均(`#797` 指出的两套不可比自助)",
           bool(len({r["scheme"] for r in ROWS}) == 2),
           f"方案 {sorted({r['scheme'] for r in ROWS})}", kind="control")
G.asserted("⑥ kill(预注册):「八题都动了」要成立,需**每一题的两个规格都排除 0**",
           bool(n_move == len(by)), f"两规格都排除 0 的题 {n_move}/{len(by)}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*98)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
else:
    still = [it for it, rs in by.items() if all(not x["excl0"] and not x["wide"] for x in rs)]
    unres = [it for it, rs in by.items() if any(not x["excl0"] and x["wide"] for x in rs)]
    moved = [it for it, rs in by.items() if all(x["excl0"] for x in rs)]
    V = (f"**`#802` 那张表缺一档,而缺的那一档正是这几题。**\n"
         f"  **动了 {len(moved)}**:{moved}\n"
         f"  **没动(区间含 0 且窄){len(still)}**:{still}\n"
         f"  **分辨不出(区间含 0 且宽){len(unres)}**:{unres}\n"
         f"  ⇒ **「缩 3/8、扩 5/8」的分母里,有 {len(still)+len(unres)} 题本来就不该被数进任何一边** ——\n"
         f"  **`#802` 的二分把「没动」与「分辨不出」全压成了一个方向,而它们是三句不同的话。**")
    if still:
        V += (f"\n  ⚠⚠ **而「没动」是一句关于人的话,页面从没说过:在 {still} 这几道题上,\n"
              f"  虔诚者与其余人之间的距离,五十年里在这个设计能分辨的范围内没有移动 ——\n"
              f"  美国在这些问题上整体改了主意,而两群人是并排走的,没有靠近也没有拉开。**")
print(V)
json.dump(dict(items=ITEMS, K=K, B=B, rows=ROWS, wide_threshold=float(MAXD),
               spec_agree=int(sum(agree.values())), n_items=len(by),
               n_moved=n_move, n_still=n_still, n_unresolved=n_unres,
               neg_control=dict(point=float(nc_pt), lo=float(nc_lo), hi=float(nc_hi), reference=0.0),
               pos_control=dict(point=float(pc_pt), lo=float(pc_lo), hi=float(pc_hi), planted=0.50),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"gaps_get_intervals.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'gaps_get_intervals.json'}")
