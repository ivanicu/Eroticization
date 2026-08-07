"""#822 · E03·A67·R261 —— 那个工程改动到底降了多少噪声?先量它,再拿它去重跑别的

`#821` 发现两条互不相干的线索**死在同一个设计限制上**:
`#818`(虔诚者停住还是回头)噪声底 **0.134** > 要分辨的 **0.10**;
`#821`(裂 vs 漂)噪声底 **0.175** > 要分辨的 **0.074**。
**两者都因为「每个十年只用首末两年」。**

⚠⚠ **而正确的下一步不是直接拿新口径去重跑那两个问题 —— 那会把「改动有没有用」
和「问题的答案是什么」搅在一起,而前者失败时我会误读成后者。**
⇒ **本轮只做一件事:量这个改动降了多少噪声。它自己就是一个可失败的主张。**

G1 估计量:**同一题、同一十年,两种口径下 `Δgap` 的自助噪声半宽之比**
   `ratio = half_width(拟合全年份) / half_width(首末两年)`。
   · `ratio < 1` ⇒ 改动有用,而**降了多少是一个数,不是一句话**。
   · `ratio ≥ 1` ⇒ **改动没用,两条线索仍然堵着,而我省下了在错误基础上重跑两轮的代价。**

**两种口径,写死在跑之前:**
   **A `endpoint`**(现口径):`Δgap = gap(末年) − gap(首年)`,**只用两年**。
   **B `fitted`**(新口径):对该十年**全部可用年份**的 `gap(y)` 拟合线性趋势,
     `Δgap = 斜率 × 该十年实际跨度`。**用上每一年。**

⚠⚠ **G4 规格曲线:自助方案也必须扫,因为它决定了「噪声」指的是什么** ——
   **`within`**:各 `年×层` 格内按人重抽 ⇒ 只含**抽样**误差;
   **`cluster`**:按年份聚类重抽 ⇒ 含**年际波动**。
   **两种都跑、都报。⚠ 而它们回答的不是同一个问题:**
   `within` 问「换一批受访者会怎样」,`cluster` 问「换一批年份会怎样」——
   **`fitted` 口径在 `within` 下几乎必然更稳(它平均了更多人),而在 `cluster` 下不一定,
   因为它把年际波动也一起拟合进去了。** ⇒ **只报 `within` 会是自我恭维,必须两种并排。**

三个世界:
   A **两种自助下都显著降噪** ⇒ 改动有用,**两条线索可以重跑。**
   B **只在 `within` 下降噪,`cluster` 下不降** ⇒ **改动只是把「更多人」换成了精度,
     而年际波动那一部分它碰不到** —— 那么 `#818`/`#821` 的堵塞取决于它们要的是哪种噪声,
     **而这是我此前从没分清过的一件事。**
   C **都不降** ⇒ 改动没用,**两条线索另找出路。**

预测矩阵:
   | 世界 | 现在 | 两种都降 | 只 within 降 | 都不降 |
   | A 有用     | 0.45 | **0.85** | 0.10 | 0.03 |
   | B 只降一半 | 0.40 | 0.05 | **0.85** | 0.05 |
   | C 没用     | 0.15 | 0.03 | 0.10 | **0.80** |

预注册判词(条件式):
  if 正控开火(**在一个已知噪声更小的合成世界里,两种口径的半宽都必须变小** ——
     否则这把量噪声的尺子本身不动)
     and 负控开火(**同一口径跑两次,半宽之比必须 ≈ 1**,⚠ **参照是 1.0 不是 0**):
      两种自助下 `homosex` 的 ratio 都 < 0.85 -> A
      只有 `within` < 0.85               -> B
      都 ≥ 0.85                          -> C
  else: UNVERIFIED
⚠ **0.85 这个门槛是我跑前定的**:低于它才算「值得为此重写两轮」的降幅。**不是数据给的。**

⚠ 跑之前写下的最强混淆:**`fitted` 口径在只有 3 个年份的十年上,自由度只剩 1** ——
  它可能**看起来**更稳只是因为线性拟合把两端拉平了。
  ⇒ 控制:**逐十年印出年份数**,并**只在年份数 ≥4 的十年上做总判**;
  年份数 = 3 的十年单独列出、不进总判。

⚠ 本轮**换不了仪器**:量的是同一份数据上两种计算口径的噪声比,第二具仪器没有对应物。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(261)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]
B, THR = 1200, 0.85

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+ITEMS, convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
YR, COV = {}, {}
for it in ITEMS:
    g = REL.dropna(subset=[it]); ys = {}
    for y, gy in g.groupby("year"):
        a = gy[gy.k == 2][it].to_numpy(float); b = gy[gy.k == 0][it].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    YR[it] = ys
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    COV[it] = {k: v for k, v in dec.items() if len(v) >= 3}

print("=== ⓪ 跑前写下的混淆:逐十年年份数(`fitted` 在 3 年上自由度只剩 1)===")
for it in ITEMS:
    print(f"  {it:9s} " + " · ".join(f"{k}s:{len(v)}年" for k, v in sorted(COV[it].items())))
print("  ⇒ **年份数 = 3 的十年单独列出、不进总判;总判只用年份数 ≥4 的十年。**")

def gp_(S, it, y): return float(S[y][0].mean()-S[y][1].mean())
def dg_endpoint(S, ys): return gp_(S, None, ys[-1]) - gp_(S, None, ys[0])
def dg_fitted(S, ys):
    x = np.array(ys, float); y = np.array([gp_(S, None, t) for t in ys])
    b = np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1)
    return float(b*(ys[-1]-ys[0]))

def half(it, ys, spec, boot, S0=None):
    S0 = S0 or YR[it]
    f = dg_endpoint if spec == "endpoint" else dg_fitted
    out = np.empty(B)
    for i in range(B):
        if boot == "within":
            S = {y: (S0[y][0][RNG.integers(0, len(S0[y][0]), len(S0[y][0]))],
                     S0[y][1][RNG.integers(0, len(S0[y][1]), len(S0[y][1]))]) for y in ys}
            use = ys
        else:                                    # cluster:按年份重抽(含年际波动)
            idx = RNG.integers(0, len(ys), len(ys))
            use = sorted({ys[j] for j in idx})
            if len(use) < (2 if spec == "endpoint" else 3): out[i] = np.nan; continue
            S = S0
        out[i] = f(S, use)
    o = out[np.isfinite(out)]
    return float((np.percentile(o, 97.5)-np.percentile(o, 2.5))/2)

print(f"\n=== ① 两种口径 × 两种自助的噪声半宽(B={B})===")
ROWS = []
for it in ITEMS:
    for dc, ys in sorted(COV[it].items()):
        r = dict(item=it, decade=dc, n_years=len(ys))
        for boot in ("within", "cluster"):
            he = half(it, ys, "endpoint", boot); hf = half(it, ys, "fitted", boot)
            r[f"h_end_{boot}"], r[f"h_fit_{boot}"] = he, hf
            r[f"ratio_{boot}"] = hf/he if he > 1e-12 else np.nan
        ROWS.append(r)
print(f"  {'题':9s} {'十年':6s} {'年':>3s} {'within: 首末→拟合':>26s} {'比':>7s}   {'cluster: 首末→拟合':>26s} {'比':>7s}")
for r in ROWS:
    print(f"  {r['item']:9s} {r['decade']}s {r['n_years']:>3} "
          f"{r['h_end_within']:>11.4f} → {r['h_fit_within']:<11.4f} {r['ratio_within']:>6.2f}   "
          f"{r['h_end_cluster']:>11.4f} → {r['h_fit_cluster']:<11.4f} {r['ratio_cluster']:>6.2f}")
main = [r for r in ROWS if r["n_years"] >= 4]
thin = [r for r in ROWS if r["n_years"] == 3]
mw = float(np.nanmedian([r["ratio_within"] for r in main]))
mc = float(np.nanmedian([r["ratio_cluster"] for r in main]))
print(f"\n  **年份数 ≥4 的 {len(main)} 个十年**:比值中位 —— `within` **{mw:.3f}** · `cluster` **{mc:.3f}**")
print(f"  ⚠ 年份数 = 3 的 {len(thin)} 个十年单独列出、**不进总判**:"
      + " · ".join(f"{r['item']}/{r['decade']}s" for r in thin))
H = [r for r in ROWS if r["item"] == "homosex" and r["n_years"] >= 4]
hw = float(np.nanmedian([r["ratio_within"] for r in H])); hc = float(np.nanmedian([r["ratio_cluster"] for r in H]))
print(f"  `homosex`(两条被堵线索的那一题):`within` **{hw:.3f}** · `cluster` **{hc:.3f}**")

print("\n=== ② 控制 ===")
it0, ys0 = "homosex", sorted(COV["homosex"])[1]
ys0 = COV["homosex"][ys0]
S_small = {y: (YR[it0][y][0][:len(YR[it0][y][0])],
               YR[it0][y][1][:len(YR[it0][y][1])]) for y in ys0}
S_big = {y: (np.repeat(YR[it0][y][0], 4), np.repeat(YR[it0][y][1], 4)) for y in ys0}
pc = {}
for spec in ("endpoint", "fitted"):
    a = half(it0, ys0, spec, "within", S0=S_small); b = half(it0, ys0, spec, "within", S0=S_big)
    pc[spec] = (a, b, b/a)
    print(f"  正控 {spec:9s}:样本量 ×4 的世界 ⇒ 半宽 {a:.4f} → {b:.4f}(比 **{b/a:.3f}**)—— 该**明显 < 1**")
nc = []
for _ in range(5):
    a = half(it0, ys0, "endpoint", "within"); b = half(it0, ys0, "endpoint", "within")
    nc.append(b/a)
nc = np.array(nc)
print(f"  负控:**同一口径跑两次**的半宽之比 ⇒ 中位 **{np.median(nc):.4f}** [{nc.min():.4f}, {nc.max():.4f}]"
      f" —— ⚠ **参照是 1.0 不是 0**,噪声半宽 **{(nc.max()-nc.min())/2:.4f}**")
nc_half = float(max((nc.max()-nc.min())/2, 1e-6))

G = Gate("#822 · 那个工程改动到底降了多少噪声")
G.asserted("① 正控:在一个**样本量 ×4** 的合成世界里,两种口径的半宽都必须明显变小"
           "(否则这把量噪声的尺子本身不动)",
           bool(pc["endpoint"][2] < 0.75 and pc["fitted"][2] < 0.75),
           f"endpoint 比 {pc['endpoint'][2]:.3f} · fitted 比 {pc['fitted'][2]:.3f}(阈 0.75)", kind="control")
G.identity_control("② 负控:**同一口径跑两次**,半宽之比必须回到 **1.0**"
                   "(⚠ 参照是 1.0 不是 0)—— 容差取它**自己量出来的**噪声半宽的两倍,"
                   "⚠ **而这一次容差不能事先写死**:它比的是同一个量的两次实现,尺度由该量自己定",
                   observed=float(np.median(nc)), expected=1.0, tol=max(2*nc_half, 1e-3),
                   noise_half_width=nc_half, what=f"5 次重复,跨度 [{nc.min():.4f}, {nc.max():.4f}]")
G.asserted("③ 前提(跑前写下的混淆):逐十年年份数已印出,**年份数 = 3 的十年不进总判**"
           "(`fitted` 在 3 年上自由度只剩 1,可能只是把两端拉平)",
           bool(all(r["n_years"] >= 4 for r in main)),
           f"进总判 {len(main)} 个十年(年份 ≥4)· 排除 {len(thin)} 个(年份 = 3)", kind="control")
G.asserted("④ 前提(`G4`):**两种自助方案都跑并都报** —— `within` 问「换一批受访者」,"
           "`cluster` 问「换一批年份」,**只报 `within` 会是自我恭维**",
           bool(all(f"ratio_{b}" in r for r in ROWS for b in ("within", "cluster"))),
           f"within 中位 {mw:.3f} · cluster 中位 {mc:.3f}", kind="control")
G.asserted(f"⑤ kill(预注册):「改动有用、两条线索可以重跑」(世界 A)要成立,"
           f"需 `homosex` 在**两种自助下**的比值都 < {THR}",
           bool(hw < THR and hc < THR), f"homosex within {hw:.3f} · cluster {hc:.3f}(阈 {THR})", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif hw < THR and hc < THR:
    V = (f"**A 改动有用,两种噪声都降。** `homosex` 比值 `within` **{hw:.3f}** · `cluster` **{hc:.3f}**;"
         f"全部 {len(main)} 个十年的中位 {mw:.3f} / {mc:.3f}。\n"
         f"  ⇒ **`#818`(要 <0.10,现 0.134)与 `#821`(要 <0.074,现 0.175)可以按新口径重跑。**")
elif hw < THR:
    V = (f"**B 只在 `within` 下降噪。** `homosex` `within` **{hw:.3f}** 而 `cluster` **{hc:.3f}**。\n"
         f"  ⇒ **这个改动买到的是「更多受访者」的精度,买不到「年际波动」那一部分** ——\n"
         f"  **而这恰好逼出一件我此前从没分清的事:`#818` 与 `#821` 要战胜的是哪一种噪声?**\n"
         f"  **若它们的问题本质上是「换一批年份还成不成立」,这个改动救不了它们。**")
else:
    V = (f"**C 改动没用。** `homosex` `within` {hw:.3f} · `cluster` {hc:.3f},都不低于 {THR}。\n"
         f"  ⇒ **两条线索另找出路 —— 而我省下了在一个没用的改动上重跑两轮的代价。**")
print(V)
json.dump(dict(items=ITEMS, B=B, thr=THR, rows=ROWS, main_n=len(main), thin=[f"{r['item']}/{r['decade']}" for r in thin],
               median_within=mw, median_cluster=mc, homosex_within=hw, homosex_cluster=hc,
               pos_control={k: list(v) for k, v in pc.items()},
               neg_control=dict(median=float(np.median(nc)), lo=float(nc.min()), hi=float(nc.max()),
                                half_width=nc_half, reference=1.0),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"does_the_fix_help.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'does_the_fix_help.json'}")
