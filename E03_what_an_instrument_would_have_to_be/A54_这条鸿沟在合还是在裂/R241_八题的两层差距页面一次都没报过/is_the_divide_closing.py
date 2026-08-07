"""#802 · E03·A54·R241 —— 八题的两层差距,页面一次都没报过:这条鸿沟在合还是在裂?

`#801` 在顶对上报了四个层的首末水平,发现两题都是**追上**。
**⚠⚠ 而那立刻暴露一件更大的事:整页发表过八题的比值 `r`,却从来没有报过任何一题的「两层差距」。**
`r` 说的是「谁动得快」;**差距说的是「他们离得近了还是远了」——那是两个不同的量,而后者才是
`E02` 那个框架(「社会拿它怎么办」)在年代这个单位上的问法。**

⚠ **而 `r` 与差距的关系不是显然的**:一个 `r < 1` 的题(虔诚者动得少)既可能差距在缩
(他们本来就落后、现在慢慢跟上),也可能差距在扩(别人跑了、他们没跟上)。
**同一个 `r` 兼容两种关于人的相反故事,而页面只发表了 `r`。**

G1 估计量:**逐题的两层差距(虔诚 − 非虔诚)在首年与末年**,以及它是缩是扩。
   ⚠ **两个规格都报**(`G4`):① **拟合端点**(各层线性趋势在首/末年的取值,用上全部年份)
   ② **原始端点**(首/末年的观测均值)。**不合并成一个数。**

⚠⚠ **而有一件事我明确拒绝做,并且要说明为什么:不报 `Δgap` 对 `gap0` 的回归。**
`Δgap = gap1 − gap0` 里**含着 `gap0` 本身** ⇒ 拿它对 `gap0` 回归是
**`realstat` 点名的那一行:「conditioning on the outcome —— binning a change score on one of its own
arms(Oldham 1962)」**,它**必然**给出负斜率,而那个负斜率**不是回归到中间的证据,是代数。**
⇒ **只报 `gap0` 与 `gap1` 两个数,让读者自己看。**

⚠⚠ 三个世界,而第三个是元分离:
   A **在合**:多数题差距在缩 ⇒ **这条宗教鸿沟在关上**,而 `r<1` 的题是「落后但在跟」。
   B **在裂**:多数题差距在扩 ⇒ **`r<1` 与「鸿沟在长」是同一件事**,
     而页面把 `r` 说成「改主意」其实说的是**极化**。
   C **两种都有,而分界由起点决定** ⇒ **差距的变化根本不是关于宗教的,是关于「谁离终点远」的** ——
     那会把整条线索的对象换掉。⚠ 但 **C 不能用 `Δgap ~ gap0` 去测**(上面那条 Oldham),
     **只能靠把 `gap0` 与 `gap1` 并排列出来让读者判**,而**我如实说这是本设计测不了的。**

预测矩阵:
   | 世界 | 现在 | 若 ≥6/8 在缩 | 若 ≥6/8 在扩 | 若混合 |
   | A 在合 | 0.45 | **0.85** | 0.03 | 0.15 |
   | B 在裂 | 0.25 | 0.05 | **0.90** | 0.20 |
   | C 起点决定 | 0.30 | 0.10 | 0.07 | **0.65**(而它测不了,只能报) |

预注册判词(条件式):
  if 正控开火(合成一个**已知**的差距变化,两个规格都要取回来)
     and 负控开火(两层轨迹**完全相同**时,`gap1` 必须等于 `gap0`):
      拟合端点下 >=6/8 在缩 -> A · >=6/8 在扩 -> B · 否则混合,报整张网格
  else: UNVERIFIED
⚠ **「这个零该不该是零?」** —— 负控问的是「轨迹相同时差距变不变」,
  而**正确的写法是比较 `gap1` 与 `gap0` 这两个值本身,不是拿它们的差与 0 比**
  (`#770` 立的规矩:等式检查要比两个值,不要比它们的差与零)⇒ **`identity_control`。**

⚠ 跑之前写下的最强混淆:**八题的量表长度不同**(2–4 档)⇒ **差距的绝对值跨题不可比。**
  ⇒ 控制:**同时报「差距 ÷ 量表跨度」**,并且**结论只用「缩/扩」的方向,不用大小跨题排序。**

⚠ 硬规则①:先打印每题的 n、真正被问过的年份、档数。
本轮换不了仪器(对象是世界;第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
⚠ 总判由 `Gate.admissible()` 决定(第五次用)。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(241)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS, OBS = P791["items"], P791["obs"]
STEM = pd.io.stata.StataReader(gp).variable_labels()

print("=== ⓪ 硬规则①:n · 真正被问过的年份 · 档数 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {}
for c in ITEMS:
    cs = list(cat[c].cat.categories); K[c] = 4 if c == "homosex" else len(cs)
    v = pd.to_numeric(d[c], errors="coerce").where(lambda x, k=K[c]: (x >= 1) & (x <= k))
    yrs = sorted(d.year[v.notna()].unique())
    print(f"  {c:9s} n={int(v.notna().sum()):>7,} · 年 {len(yrs):>2}({int(min(yrs))}–{int(max(yrs))}) · "
          f"档 {K[c]} · 「{STEM.get(c,'?')[:42]}」")

M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def series(item, k, nmin=120):
    g = REL.dropna(subset=[item])
    return [(int(y), float(gy[item].mean())) for y, gy in g[g.k == k].groupby("year") if len(gy) >= nmin]

def gaps(item, spec, rowsA=None, rowsB=None):
    """spec='fit' 用各层线性趋势在共同首/末年的取值;spec='raw' 用观测均值。"""
    A = rowsA if rowsA is not None else series(item, 2)
    B = rowsB if rowsB is not None else series(item, 0)
    if len(A) < 8 or len(B) < 8: return None
    y0 = max(A[0][0], B[0][0]); y1 = min(A[-1][0], B[-1][0])
    out = {}
    for lab, R in (("A", A), ("B", B)):
        x = np.array([r[0] for r in R], float); y = np.array([r[1] for r in R])
        if spec == "fit":
            b = np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1); a = y.mean()-b*x.mean()
            out[lab] = (a+b*y0, a+b*y1)
        else:
            out[lab] = (float(next(r[1] for r in R if r[0] == y0)),
                        float(next(r[1] for r in R if r[0] == y1)))
    return dict(y0=y0, y1=y1, gap0=out["A"][0]-out["B"][0], gap1=out["A"][1]-out["B"][1])

print("\n=== ① 逐题两层差距(虔诚 − 非虔诚)· 两个规格都报(`G4`)· ⚠ 量表长度不同 ⇒ 同时报「÷ 跨度」 ===")
TBL = {}
for spec in ("fit", "raw"):
    TBL[spec] = {}
    print(f"  —— 规格「{'拟合端点' if spec=='fit' else '原始端点'}」")
    for c in ITEMS:
        g = gaps(c, spec)
        if g is None: continue
        span = K[c]-1
        narrowing = abs(g["gap1"]) < abs(g["gap0"])
        crossed = np.sign(g["gap0"]) != np.sign(g["gap1"]) and abs(g["gap1"]) > 1e-9
        TBL[spec][c] = dict(**g, span=span, gap0_n=g["gap0"]/span, gap1_n=g["gap1"]/span,
                            narrowing=bool(narrowing), crossed=bool(crossed), r=OBS[c])
        t = TBL[spec][c]
        print(f"    {c:9s} r={OBS[c]:+6.3f} · {g['y0']}→{g['y1']} · 差距 {g['gap0']:+.3f} → {g['gap1']:+.3f}"
              f"  (÷跨度 {t['gap0_n']:+.3f} → {t['gap1_n']:+.3f})  "
              f"**{'缩' if narrowing else '扩'}**{' · 跨过 0' if crossed else ''}")
FIT = TBL["fit"]
n_nar = sum(v["narrowing"] for v in FIT.values()); n_tot = len(FIT)
agree = sum(1 for c in FIT if c in TBL["raw"] and FIT[c]["narrowing"] == TBL["raw"][c]["narrowing"])
print(f"\n  拟合端点:**缩 {n_nar}/{n_tot}** · 两个规格方向一致的题 **{agree}/{n_tot}**")

# ── ② 控制 ────────────────────────────────────────────────────────────────────
print("\n=== ② 控制 ===")
# ⚠⚠ 第一版把负控写成「两层轨迹**完全相同**」⇒ `gap0 = gap1 = 0.0` ⇒
#    `identity_control(0.0, 0.0)` **两个输入都恰好是零,库当场判 DEGENERATE**。
#    ⇒ 我以为在守 `#770`(比两个值而不是比差与零),**而两个值恰好都是零时,那条规矩没被守住**。
#    ⇒ 改成:两层**斜率相同但截距差 0.3**(差距按构造恒为 0.3,非零),
#      **并且给两层不同的年份集合**(A 缺若干年)—— 这样年份对齐那段逻辑**真的被走到**,
#      控制能因为实现错误而失败,而不是恒真。
yy = list(range(1974, 2025, 2))
base = [(y, 2.0 - 0.01*(y-1974)) for y in yy]
offs = [(y, 2.0 - 0.01*(y-1974) + 0.30) for y in yy if y % 6 != 0]   # ⚠ 年份集合不同
nc = gaps("_syn", "fit", rowsA=offs, rowsB=list(base))
print(f"  负控:两层斜率相同、截距差 0.30、**年份集合不同**(A {len(offs)} 年 · B {len(base)} 年)")
print(f"        → gap0 = {nc['gap0']:+.6f} · gap1 = {nc['gap1']:+.6f}  "
      f"(⚠ 两个值都非零 ⇒ 不退化;**它能因年份对齐写错而失败**)")
wide = [(y, 2.0 - 0.004*(y-1974)) for y in yy]     # 走得慢 ⇒ 差距必然扩大
pc = gaps("_syn", "fit", rowsA=wide, rowsB=list(base))
print(f"  正控:B 层走得快、A 层走得慢 → gap0 = {pc['gap0']:+.4f} → gap1 = {pc['gap1']:+.4f} "
      f"⇒ 判「{'扩' if abs(pc['gap1'])>abs(pc['gap0']) else '缩'}」(该是**扩**)")
pc_ok = abs(pc["gap1"]) > abs(pc["gap0"]) + 1e-9

G = Gate("#802 · 这条鸿沟在合还是在裂")
G.identity_control("① 负控:两层同斜率、截距差 0.30、年份集合不同 ⇒ `gap1` 必须等于 `gap0`(且两者非零)",
                   observed=float(nc["gap1"]), expected=float(nc["gap0"]), tol=1e-9,
                   what="同斜率不同截距、且年份集合不同的合成世界 —— 差距按构造恒为 0.30,"
                        "而年份对齐若写错就会破坏它")
G.asserted("② 正控:一层走得慢时,判据必须判「扩」(否则它连已知在裂的都认不出)",
           pc_ok, f"gap {pc['gap0']:+.4f} → {pc['gap1']:+.4f}", kind="control")
G.asserted("③ 前提(跑前写下的混淆):量表长度 2–4 档不同 ⇒ 同时报「÷ 跨度」,"
           "且结论只用方向不用跨题排序",
           bool(len({K[c] for c in FIT}) > 1),
           f"档数 {sorted({K[c] for c in FIT})} —— 已同时报归一化差距", kind="control")
G.asserted("④ 前提:两个规格(拟合/原始)的方向必须多数一致,否则「缩/扩」是规格的产物",
           bool(agree >= n_tot-1), f"方向一致 {agree}/{n_tot}", kind="control")
G.asserted("⑤ kill(预注册):「鸿沟在合」要站住,需拟合端点下 ≥6/8 在缩",
           bool(n_nar >= 6), f"缩 {n_nar}/{n_tot}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 仪器没资格下判。**"
elif n_nar >= 6:
    cross = [c for c in FIT if FIT[c]["crossed"]]
    v = (f"**A 在合。** 拟合端点下 **{n_nar}/{n_tot}** 题的两层差距在缩(两个规格方向一致 {agree}/{n_tot})"
         f"{f';其中 {cross} 跨过了 0 —— 虔诚层从一侧走到了另一侧' if cross else ''}。\n"
         f"  ⇒ **而这改写了 `r` 的读法:`r<1` 的那些题不是「鸿沟在长」,是「他们落后但在跟」。**")
elif (n_tot-n_nar) >= 6:
    v = (f"**B 在裂。** 拟合端点下只有 {n_nar}/{n_tot} 在缩 ⇒ **多数题的两层差距在扩** ——\n"
         f"  **`r<1` 与「鸿沟在长」是同一件事,而页面把 `r` 说成「改主意」其实说的是极化。**")
else:
    v = (f"**混合:缩 {n_nar}/{n_tot}。** 报整张网格,不选边。\n"
         f"  ⚠ 而**是不是「起点决定了缩扩」(世界 C)本设计测不了** —— "
         f"`Δgap` 里含着 `gap0`,拿它对 `gap0` 回归是 Oldham 的条件于结果,**必然给负斜率**。\n"
         f"  **只报两个数,不做那个回归。**")
print(v)
print(f"\n⚠ 无论哪一支:**`Δgap ~ gap0` 的回归本轮不做,也不许后来者拿本轮的数去做** —— "
      f"那是 `realstat` 点名的 Oldham 1962,**它给出的负斜率是代数不是发现。**")
json.dump(dict(items=ITEMS, K=K, stems={c: STEM.get(c, "") for c in ITEMS},
               fit=TBL["fit"], raw=TBL["raw"], n_narrowing=n_nar, n_total=n_tot, spec_agree=agree,
               neg_control=nc, pos_control=pc, admissible=adm, verdict=v,
               refused="Δgap ~ gap0 回归(Oldham 1962:条件于结果)", gate_ok=G.verdict()),
          open(OUT/"is_the_divide_closing.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'is_the_divide_closing.json'}")
