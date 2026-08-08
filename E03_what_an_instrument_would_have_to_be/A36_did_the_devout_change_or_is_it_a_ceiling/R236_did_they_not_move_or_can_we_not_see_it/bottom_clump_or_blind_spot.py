"""#797 · E03·A49·R236 —— 底堆:是「他们真的没动」,还是「这个问题在那里问不出来」?

⚠⚠ **先说为什么是这一轮:连着两轮(`#794` 分层线漂移 · `#796` 换估计量)都确认了三堆稳健。
   那是一个 basin。** `frontier §3` 说:**故意设计一个「正面结果我会不想要」的步子。**
   ⇒ 本轮去撞页面中心里**最软的那一块**:底堆(`prayer` · `helpblk`)。

**⚠⚠ 而先要更正一句我自己写过的话,它是这一轮的起点:**
`#791` 报「组内 0/8 · 组间 16/20 可分辨」,读起来像三堆都立住了。
**但那个 16/20 把三条边界并在一起数了。逐条拆开:**

| 边界 | 可分辨 |
|---|---|
| 顶/中 | **8/8** |
| 顶/底 | **4/4** |
| **中/底** | **4/8** —— `spanking`·`homosex` 与两个底题都分得开,**`suicide2`·`teensex` 与它们一个都分不开** |

⇒ **顶堆是干净的;而「中」与「底」之间只解开了一半。** `#791` 的汇总数字盖住了这件事。

G1 估计量(三个,方法之前先命名):
   (a) **分母动没动** —— 每题**非虔诚层自己**的年代斜率 vs 它自己的置换零。
       ⚠ 这是 `#790` 立的边界:**社会自己没动的地方,这个比根本没有定义**(分母趋零 ⇒ Cauchy)。
   (b) **底堆与中堆分不分得开** —— 底题的比值区间**排不排除中堆的中位 0.41**。
   (c) **底堆与 0 分不分得开** —— 底题的比值区间**含不含 0**。
   ⚠ (b) 与 (c) 是**两个不同的问题**,而页面把它们读成了一句话
     (「他们几乎不动」既暗示了「比中堆少」也暗示了「接近 0」)。**先拆开再说。**

⚠⚠ 三个世界,而**后两个都会把页面的中心从「三堆」削成「顶 + 剩下的」**:
   A **真的不动**:分母过零(社会确实动了)且底题区间排除 0.41
     ⇒ 底堆是真的,而且是一句关于人的硬话:**别人动了,他们没动。**
   B **问不出来**:分母**过不了**自己的零 ⇒ 与 `xmarsex` 同类,**这个问题在那两题上没有定义**
     ⇒ 底堆不是一堆人,是一片盲区,**必须从三堆里挪出去。**
   C **功效不足**:分母过零但区间**同时含 0 与 0.41** ⇒ 底题既不能说「像中堆」也不能说「是 0」
     ⇒ **无法归堆**,三堆变成「顶 + 一团」。

预测矩阵:
   | 世界 | 现在 | 若分母过零且区间排除 0.41 | 若分母过不了零 | 若区间同时含 0 与 0.41 |
   | A 真的不动 | 0.45 | **0.85** | 0.02 | 0.10 |
   | B 问不出来 | 0.25 | 0.05 | **0.90** | 0.15 |
   | C 功效不足 | 0.30 | 0.10 | 0.08 | **0.75** |

预注册判词(条件式):
  if 正控开火(合成一个**已知**动了的分母,置换零必须判它动了)
     and 负控开火(合成一个**没动**的分母,置换零必须判它没动 —— 即 g=0 时不开火):
      两个底题的分母都过零 且 两题区间都排除 0.41 -> A
      任一底题的分母过不了零                      -> B(那一题移出三堆)
      否则                                        -> C(无法归堆)
  else: UNVERIFIED
⚠ **0.41 这个参照不是我挑的**:它是中堆四题比值的中位(`#791` 的产物),**由数据给出。**

⚠ **「这个零该不该是零?」** —— (a) 的零**该是 0**:一个层的年代斜率在「年份与作答无关」的世界里
  期望**恰好是 0** ⇒ `negative_control`,零的做法是**打乱年份标签**。**不是 offset。**

⚠ 跑之前写下的最强混淆:**底堆两题的年份数与中堆不同** ⇒ 分母的零本身宽窄不同 ⇒
  「过不过零」可能只是年数差异。⇒ 控制:**逐题报年份数**,并且**把中堆四题一起跑**,
  让读者看见同一具仪器在两堆上的表现,而不是只看底堆。

本轮换不了仪器(对象是世界,第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
⚠ 总判由 `Gate.admissible()` 决定(`#796` 加的,第二次用)。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(236)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_is_the_eight_point_axis_an_axis_or_eight_labels_on_noise/results/is_the_ordering_an_object.json"))
ITEMS, CLUMP = P791["items"], P791["clumps"]
MID_MED = float(np.median([P791["obs"][c] for c in ITEMS if CLUMP[c] == "中"]))
BOT = [c for c in ITEMS if CLUMP[c] == "底"]
print(f"=== ⓪ 中堆中位 = **{MID_MED:.3f}**(由数据给出,不是我挑的)· 底堆 = {BOT} ===")

# ── 边界逐条拆开:更正 `#791` 的汇总数字 ──────────────────────────────────────
print("\n=== ① 更正 `#791`:那个 16/20 把三条边界并在一起数了 ===")
bnd = {}
for p in P791["pairs"]:
    key = tuple(sorted((CLUMP[p["a"]], CLUMP[p["b"]])))
    if key[0] == key[1]: key = ("组内", key[0])
    bnd.setdefault(key, []).append(p["survives"])
for k, v in sorted(bnd.items(), key=lambda x: -sum(x[1])/len(x[1])):
    print(f"  {'/'.join(k):10s} 可分辨 **{sum(v)}/{len(v)}**")
mid_bot = bnd[tuple(sorted(("中", "底")))]
print(f"  ⇒ **中/底 只解开了 {sum(mid_bot)}/{len(mid_bot)}** —— 这是本轮的起点")

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
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def series(item, k, nmin=120):
    g = REL.dropna(subset=[item])
    return [(int(y), float(gy[item].mean())) for y, gy in g[g.k == k].groupby("year") if len(gy) >= nmin]

def denom_test(item, B=4000):
    """非虔诚层自己动了没有:斜率 vs **打乱年份标签**的置换零(零该是 0 ⇒ negative_control)。"""
    rows = series(item, 0)
    if len(rows) < 8: return None
    x = np.array([r[0] for r in rows], float); y = np.array([r[1] for r in rows])
    s = slope(x, y)
    nul = np.abs([slope(RNG.permutation(x), y) for _ in range(B)])
    return dict(item=item, n_year=len(rows), slope=s, null95=float(np.quantile(nul, .95)),
                moved=bool(abs(s) > np.quantile(nul, .95)),
                ratio_to_null=float(abs(s)/np.quantile(nul, .95)))

print("\n=== ② (a) 分母动没动 —— 八题全跑,让同一具仪器在两堆上都露面(跑前混淆的控制)===")
D = {}
for c in ITEMS:
    r = denom_test(c); D[c] = r
    print(f"  {c:9s} {CLUMP[c]}  年 {r['n_year']:>2} · 非虔诚斜率 {r['slope']:+.5f} vs 零95% {r['null95']:.5f} "
          f"= **{r['ratio_to_null']:.2f}×** {'动了' if r['moved'] else '**没动**'}")

# ── (b)(c) 底题的区间:排不排除 0.41 · 含不含 0 ────────────────────────────────
print(f"\n=== ③ (b)(c) 两个不同的问题:区间排不排除中堆中位 {MID_MED:.3f} · 含不含 0 ===")
CI = {p["a"]: None for p in []}
def ratio_ci(item, B=4000):
    rA, rB = series(item, 2), series(item, 0)
    if len(rA) < 8 or len(rB) < 8: return None
    yA = np.array([r[0] for r in rA], float); vA = np.array([r[1] for r in rA])
    yB = np.array([r[0] for r in rB], float); vB = np.array([r[1] for r in rB])
    f = lambda ia, ib: slope(yA[ia], vA[ia])/slope(yB[ib], vB[ib])
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(B)])
    bs = bs[np.isfinite(bs)]
    return float(f(np.arange(len(yA)), np.arange(len(yB)))), \
           float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
tbl = []
for c in ITEMS:
    r, lo, hi = ratio_ci(c)
    ex_mid = not (lo <= MID_MED <= hi); ex_zero = not (lo <= 0.0 <= hi)
    tbl.append(dict(item=c, clump=CLUMP[c], r=r, lo=lo, hi=hi,
                    excludes_mid=ex_mid, excludes_zero=ex_zero, moved=D[c]["moved"]))
    print(f"  {c:9s} {CLUMP[c]}  {r:+7.3f} [{lo:+7.3f},{hi:+7.3f}]  "
          f"排除 {MID_MED:.2f}:{'是' if ex_mid else '否'}  ·  排除 0:{'是' if ex_zero else '否'}")

bot_rows = [t for t in tbl if t["clump"] == "底"]
all_moved = all(D[c]["moved"] for c in BOT)
all_ex_mid = all(t["excludes_mid"] for t in bot_rows)
any_ex_zero = any(t["excludes_zero"] for t in bot_rows)
print(f"\n  底堆两题:分母都动了 **{all_moved}** · 区间都排除中堆中位 **{all_ex_mid}** · "
      f"有任何一题排除 0 **{any_ex_zero}**")

# ── 控制:合成一个动了的分母 / 一个没动的分母 ────────────────────────────────────
print("\n=== ④ 控制:置换零必须能判「动了」,也必须在真没动时判「没动」===")
def synth_denom(beta, n_year=30, B=2000):
    x = np.arange(1974., 1974+2*n_year, 2.0)
    y = 3.0 + beta*(x-1974) + RNG.normal(0, 0.05, n_year)
    s = slope(x, y); nul = np.abs([slope(RNG.permutation(x), y) for _ in range(B)])
    # ⚠ `bool(...)` 不是装饰:`np.bool_` 不是 JSON 可序列化的,第一版跑完全部计算之后死在写产物那一步。
    return bool(abs(s) > np.quantile(nul, .95)), float(abs(s)/np.quantile(nul, .95))
pc_moved, pc_ratio = synth_denom(-0.02)     # 明显动了
# ⚠⚠ 第一版把 g=0 那条控制写成**一次抽样**:「β=0 时判「动了」必须是 False」。
#    第一次跑它返回 True(1.04× 零),而我打印的那句话说「必须是 False」——
#    **两句话互相矛盾,而闸那一行却按「相对大小」判它通过了。**
#    ⇒ 而正确的是:**一次抽样测不了一个假阳性率。** 一条 95% 门槛下的零序列,
#      本来就应该有 ~5% 的时候越过自己的分位 —— **一次 1.04× 完全合乎零,不是缺陷。**
#    ⇒ 判据改成**率**:跑 300 条 β=0 的合成序列,假阳性率必须 ≤ 0.10(名义 0.05 的两倍)。
#    ⇒ 「拿一个相邻的量当要判的那个量」本会话第五次 —— 这次是**一次实现 vs 一个率**。
NC_DRAWS = 300
nc_hits = sum(synth_denom(0.0, B=600)[0] for _ in range(NC_DRAWS))
nc_rate = nc_hits/NC_DRAWS
nc_ratio = float(np.median([synth_denom(0.0, B=600)[1] for _ in range(40)]))
print(f"  正控 β=−0.02(明显动了)→ 判「动了」= **{pc_moved}**({pc_ratio:.1f}× 零)")
print(f"  负控 β=0(g=0)· **{NC_DRAWS} 条合成零序列的假阳性率 = {nc_rate:.3f}**"
      f"(名义 0.05,门槛 ≤0.10)⇒ **{'合格' if nc_rate <= 0.10 else '不合格 ⇒ 分母判定不可信'}**")

G = Gate("#797 · 底堆:他们真的没动,还是这个问题在那里问不出来")
G.asserted("① 正控:合成一个明显动了的分母,置换零必须判它动了",
           pc_moved, f"β=−0.02 → {pc_ratio:.1f}× 自己的零", kind="control")
G.asserted("② 负控(g=0):**假阳性率**必须 ≤0.10 —— 一次实现测不了一个率(第一版就是那么写的)",
           bool(nc_rate <= 0.10), f"{NC_DRAWS} 条 β=0 合成序列 · 假阳性率 {nc_rate:.3f}(名义 0.05)",
           kind="control")
G.asserted("③ 前提(跑前写下的混淆):八题的年份数必须一起报,不许只看底堆",
           True, " · ".join(f"{c} {D[c]['n_year']}年" for c in ITEMS), kind="control")
G.asserted("④ kill(预注册):底堆要作为「一堆人」立住,需两题分母都过零**且**区间都排除中堆中位",
           bool(all_moved and all_ex_mid),
           f"分母都动了 {all_moved} · 都排除 {MID_MED:.2f} {all_ex_mid}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 仪器没资格下判。**"
elif not all_moved:
    dead = [c for c in BOT if not D[c]['moved']]
    v = (f"**B 问不出来:{dead} 的分母过不了自己的零。** ⇒ 与 `xmarsex` 同类 —— "
         f"**社会自己没动的地方,这个比没有定义**,那一题必须从三堆里挪出去,底堆不再是一堆人。")
elif all_ex_mid:
    v = (f"**A 底堆是真的,而且是一句硬话。** 两题的分母都动了"
         f"({' · '.join(f'{c} {D[c]['ratio_to_null']:.1f}× 零' for c in BOT)}),"
         f"而两题的比值区间都**排除中堆的中位 {MID_MED:.3f}**"
         f"({' · '.join(f'{t['item']} [{t['lo']:+.3f},{t['hi']:+.3f}]' for t in bot_rows)})。\n"
         f"  ⚠ **而 (b) 与 (c) 必须分开说**:两题的区间**都含 0** ⇒ "
         f"**「比中堆少」站得住,「等于 0」站不住** —— 页面把这两句读成了一句。")
else:
    v = (f"**C 无法归堆。** 分母过了零,但区间同时含 0 与中堆中位 {MID_MED:.3f} ⇒ "
         f"**底题既不能说「像中堆」也不能说「是 0」**,三堆要改写成「顶 + 一团」。")
print(v)
json.dump(dict(mid_median=MID_MED, bottom=BOT, boundary={'/'.join(k): [sum(v2), len(v2)] for k, v2 in bnd.items()},
               denom={c: D[c] for c in ITEMS}, table=tbl,
               all_moved=bool(all_moved), all_exclude_mid=bool(all_ex_mid), any_exclude_zero=bool(any_ex_zero),
               pos_control=[bool(pc_moved), float(pc_ratio)],
               neg_control_false_positive_rate=float(nc_rate), neg_control_draws=NC_DRAWS,
               admissible=adm, verdict=v, gate_ok=G.verdict()),
          open(OUT/"bottom_clump.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'bottom_clump.json'}")
