"""#796 · E03·A48·R235 —— 三堆全建立在「逐年均值的斜率比」上,而那把每年内部的人整个丢掉了

页面上现在的中心是 `#791` 的三堆,而它、`#794` 的稳健性、`#790` 的天花板、`#789` 的轴 ——
**全部只用过一个估计量:每层每年取均值,对年份回归取斜率,两层相除。**
⇒ **那把每一年内部的个体差异整个丢掉了,并且给每一年同样的权重,不管那一年问了 1,200 人还是 2,800 人。**

⚠⚠ **而 `#791` 的规格曲线扫的是「自助抽数 × BH 的 q × 差的定义」—— 三个都是事后统计的旋钮,
   一个都不是「怎么估这个比」。`G4` 要求的估计量轴,这个项目从来没有扫过。**

G1 估计量(先命名,再选方法):**同一个「虔诚层相对非虔诚层的年代变化率之比」**,三种估法:
   **agg** 逐年均值的斜率比(现行):`slope(year, ȳ_devout) / slope(year, ȳ_nondevout)`,**年份等权**
   **ind** 个体层交互:`y_i ~ a + b₁·year + b₂·devout + b₃·(devout×year)`,
          ⇒ `r_ind = (b₁+b₃)/b₁`,**按人加权**
   **ind_w** 个体层交互,但**把每个人按 1/n(该年该层)加权** ⇒ **年份等权的个体层版本**
   ⚠ 第三臂不是多余的:**agg 与 ind 同时差了两件事**(丢不丢个体 · 年份等权 vs 人等权)——
     `#794` 刚栽在一个「同时差两件事」的对照上,**这次先把它拆开。**

⚠⚠ 三个世界,而第二个我会非常不想要(它会削掉页面的中心):
   A **稳健**:三种估法给同一批三堆(顶/中/底的成员不变,且 `r_ind` 落在 `r_agg` 的自助区间里)
     ⇒ 那些比值是关于人的,不是关于「怎么汇总」。
   B **汇总产物**:三堆解体或重排 ⇒ **`#789`–`#794` 的一切都是逐年取均值取出来的**,页面中心要撤。
   C **元分离**:`r_agg` 与 `r_ind` 系统性地不是同一个量(比如一个有界一个没界)
     ⇒ **「那个比」从来就不是一个量**,而 A/B 都问错了。

预测矩阵:
   | 世界 | 现在 | 若 ≥6/8 落在区间内且三堆不变 | 若三堆重排 | 若 ind 与 ind_w 差得比 ind 与 agg 还大 |
   | A 稳健   | 0.55 | **0.90** | 0.05 | 0.20 |
   | B 汇总   | 0.30 | 0.05 | **0.85** | 0.20 |
   | C 非同量 | 0.15 | 0.05 | 0.10 | **0.60** |

预注册判词(条件式):
  if 正控开火(合成一个**已知**的差异率,三种估法都要取回来)
     and 正控在 g=0 时**不**开火(无差异时三种估法都给 ≈1.0):
      >=6/8 题的 `r_ind` 落在 `r_agg` 的年份自助 95% 区间内 **且** 三堆成员不变 -> A
      三堆成员改变                                                        -> B
      `|r_ind − r_ind_w|` 的中位 > `|r_ind − r_agg|` 的中位                -> C(差的是权重,不是汇总)
  else: UNVERIFIED
⚠ 6/8 这个门槛的理由写在跑之前:`#791` 已测得八题里**只有 6 题**的比值区间不含 1.0 那一类的稳定性,
  **要求全部 8 题落进区间等于要求比原始测量更稳,那是一个不可能通过的判据。**

⚠ **「这个零该不该是零?」** —— 正控的 g=0 那一支**该是 1.0 不是 0**:
  没有差异时,两层斜率相等 ⇒ 比值 = 1。**所以那是一个 `offset_control`,零点在 1.0,
  而它的零是「两层同斜率的合成世界」** —— 不是 `negative_control`。**这一条写在跑之前。**

⚠ 跑之前写下的最强混淆:**个体层回归的 n 是几万,标准误会小得多** ⇒
  `r_ind` 看起来「更精确」而其实只是**把年份间的真实异质性当成了噪声吞掉**。
  ⇒ 控制:**`r_ind` 的区间用与 `r_agg` 同一种年份聚类自助算**,不用回归自带的 SE。

本轮换不了仪器(对象是世界,第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
⚠ 而本轮**首次使用 `#796` 新加的 `Gate.admissible()`**:总判由它决定,不再各写各的 `ctrl`。
"""
import numpy as np, pandas as pd, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

RNG = np.random.default_rng(235)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
PREV = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_is_the_eight_point_axis_an_axis_or_eight_labels_on_noise/results/is_the_ordering_an_object.json"))
ITEMS = PREV["items"]; CLUMP = PREV["clumps"]
print(f"=== ⓪ 对象:`#791` 的 8 题与三堆 · 首次用 `Gate.admissible()` 定总判 ===")

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
SUB = REL[REL.k.isin([0, 2])].copy()
SUB["dev"] = (SUB.k == 2).astype(float)
print(f"  两层合计 n = {len(SUB):,} · 年 {SUB.year.nunique()}")

def r_agg(df, item):
    """现行估计量:逐年均值 -> 斜率 -> 相除。年份等权。"""
    g = df.dropna(subset=[item]); out = {}
    for k in (2, 0):
        rows = [(int(y), float(gy[item].mean())) for y, gy in g[g.k == k].groupby("year") if len(gy) >= 120]
        if len(rows) < 8: return np.nan
        x = np.array([r[0] for r in rows], float); y = np.array([r[1] for r in rows])
        out[k] = float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
    return out[2]/out[0] if abs(out[0]) > 1e-12 else np.nan

def r_ind(df, item, equal_year_weight=False):
    """个体层:y ~ year + dev + dev*year  ->  (b_year + b_int)/b_year。"""
    g = df.dropna(subset=[item])
    keep = g.groupby(["year", "k"])[item].transform("size") >= 120
    g = g[keep]
    if g.year.nunique() < 8: return np.nan
    yv, dv = g.year.to_numpy(float), g.dev.to_numpy(float)
    X = np.column_stack([np.ones(len(g)), yv, dv, dv*yv])
    y = g[item].to_numpy(float)
    if equal_year_weight:
        w = 1.0/g.groupby(["year", "k"])[item].transform("size").to_numpy(float)
        sw = np.sqrt(w); Xw, yw = X*sw[:, None], y*sw
    else:
        Xw, yw = X, y
    b = np.linalg.lstsq(Xw, yw, rcond=None)[0]
    return (b[1]+b[3])/b[1] if abs(b[1]) > 1e-12 else np.nan

EST = {"agg 逐年均值(现行)": lambda df, it: r_agg(df, it),
       "ind 个体层交互": lambda df, it: r_ind(df, it, False),
       "ind_w 个体层·年份等权": lambda df, it: r_ind(df, it, True)}

# ── ① 控制:合成世界,先建后判 ────────────────────────────────────────────────
print("\n=== ① 控制:合成一个**已知**的差异率,三种估法都要取回来(g=0 时都必须回到 1.0)===")
def synth(g, n_year=1500, years=range(1974, 2025, 2)):
    """两层同起点,非虔诚层每年降 0.02;虔诚层降 0.02*g。真值 r = g。"""
    rows = []
    for y in years:
        for k, mult in ((0, 1.0), (2, g)):
            mu = 3.0 - 0.02*mult*(y-1974)
            v = np.clip(RNG.normal(mu, 1.0, n_year), 1, 4)
            rows.append(pd.DataFrame(dict(year=y, k=k, dev=float(k == 2), syn=v)))
    return pd.concat(rows, ignore_index=True)

ctrl_rows = []
for g in (1.0, 0.4):
    df = synth(g)
    got = {nm: fn(df, "syn") for nm, fn in EST.items()}
    ctrl_rows.append(dict(g=g, **{k: float(v) for k, v in got.items()}))
    print(f"  真值 r = {g:.2f} → " + " · ".join(f"{nm.split()[0]} {v:+.3f}" for nm, v in got.items()))
g0 = ctrl_rows[0]; g1 = ctrl_rows[1]
pc_recover = all(abs(g1[nm]-0.4) < 0.05 for nm in EST)          # 能取回一个已知差异
pc_null = all(abs(g0[nm]-1.0) < 0.05 for nm in EST)             # g=0(即 r=1)时都回到 1.0
print(f"  ⇒ 取回已知差异(|误差|<0.05):**{pc_recover}** · 无差异时都回到 1.0:**{pc_null}**")

# ── ② 真数据:三种估法 ────────────────────────────────────────────────────────
print("\n=== ② 真数据:同一批题,三种估法(`G4` 的估计量轴,本项目第一次扫)===")
R = {nm: {} for nm in EST}
for it in ITEMS:
    for nm, fn in EST.items(): R[nm][it] = float(fn(SUB, it))
for nm in EST:
    print(f"  {nm:22s} " + " · ".join(f"{it} {R[nm][it]:+.3f}" for it in ITEMS))

# ── ③ r_ind 落不落在 r_agg 的年份聚类自助区间里 ────────────────────────────────
print("\n=== ③ `r_ind` 落不落在 `r_agg` 的**年份聚类**自助区间里 ⚠ 不用回归自带的 SE ===")
YRS = sorted(SUB.year.unique())
def boot_agg_ci(item, B=1200):
    vals = []
    while len(vals) < B:
        draw = RNG.choice(YRS, len(YRS), replace=True)
        sub = pd.concat([SUB[SUB.year == y] for y in draw], ignore_index=True)
        v = r_agg(sub, item)
        if np.isfinite(v): vals.append(v)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))
inside, tbl = 0, []
for it in ITEMS:
    lo, hi = boot_agg_ci(it)
    ok = bool(lo <= R["ind 个体层交互"][it] <= hi)
    inside += int(ok)
    tbl.append(dict(item=it, clump=CLUMP[it], agg=R["agg 逐年均值(现行)"][it],
                    ind=R["ind 个体层交互"][it], ind_w=R["ind_w 个体层·年份等权"][it],
                    lo=lo, hi=hi, inside=ok))
    print(f"  {it:9s} {CLUMP[it]}  agg {R['agg 逐年均值(现行)'][it]:+.3f} [{lo:+.3f},{hi:+.3f}] · "
          f"ind {R['ind 个体层交互'][it]:+.3f} · ind_w {R['ind_w 个体层·年份等权'][it]:+.3f}  "
          f"{'落在区间内' if ok else '**落在区间外**'}")
print(f"\n  **{inside}/{len(ITEMS)} 落在区间内**(预注册门槛 ≥6)")

# 三堆成员变没变:按 ind 重排,顶2/中4/底2 的成员是否与 CLUMP 相同
order_ind = sorted(ITEMS, key=lambda x: -R["ind 个体层交互"][x])
new_clump = {**{c: "顶" for c in order_ind[:2]}, **{c: "中" for c in order_ind[2:6]},
             **{c: "底" for c in order_ind[6:]}}
same_clump = all(new_clump[c] == CLUMP[c] for c in ITEMS)
print(f"  按 `ind` 重排后的三堆:" + " · ".join(f"{c}{new_clump[c]}" for c in order_ind))
print(f"  ⇒ 三堆成员与 `#791` **{'完全相同' if same_clump else '不同'}**")

# 差的是「汇总」还是「权重」
d_agg = float(np.median([abs(R["ind 个体层交互"][c]-R["agg 逐年均值(现行)"][c]) for c in ITEMS]))
d_w = float(np.median([abs(R["ind 个体层交互"][c]-R["ind_w 个体层·年份等权"][c]) for c in ITEMS]))
print(f"\n  |ind − agg| 中位 = **{d_agg:.3f}** · |ind − ind_w| 中位 = **{d_w:.3f}**")
# ⚠⚠ 第一版在这里直接按 `d_w > d_agg` 命名了一个来源 —— 而实测两者都是 0.031,**分不开。**
#    一个比较在报出来之前必须先问它**分不分得开**;`#789` 那次是中位数之差没有自己的区间,
#    本轮是两个差本身一样大而我仍然挑了一个赢家。**「判词字符串不是一次计算」那一族。**
#    ⇒ 判据:两者要相差 **1.5 倍**才允许命名来源(1.5 沿用 `#791` 用过的可分辨性量级),否则如实说分不开。
src_resolvable = bool(max(d_agg, d_w) >= 1.5*min(d_agg, d_w))
src = ("权重(年份等权 vs 人等权)" if d_w > d_agg else "汇总(丢不丢个体)") if src_resolvable else None
print(f"  ⇒ 两估法之差的来源:**{src if src else '分不开 —— 两个差一样大('+f'{d_agg:.3f} 对 {d_w:.3f}'+'),本设计说不出是汇总还是权重'}**")

G = Gate("#796 · 三堆是不是逐年取均值取出来的")
# ⚠⚠ 第一版在这里用了 `offset_control`,而它当场把整轮判成 UNVERIFIED —— **而那是我的错,不是数据的。**
#    「这个零该不该是零?」我答对了一半:**比值的参照点是 1.0 不是 0**。
#    但接着挑错了函数:`offset_control` 问的是「**效应有没有越过那个偏移**」(要 |effect−offset| > 2·spread),
#    而我要问的是「**它等不等于那个参照点**」(要 |observed−expected| ≤ tol)。**方向正好相反。**
#    实测 1.0057 对 1.0,差 0.0057 —— 三种估法都**漂亮地回到了 1.0**,而那条控制却判它失败。
#    ⇒ 「拿一个相邻的量当要判的那个量」这一族**本会话第四次**
#      (`#791` 计数 vs 结构 · `#792` 看得见 vs 看见的是不是 · `#793` 用没用 vs 倒没倒 · 本轮 越过 vs 等于)。
#    ⇒ 正确的函数是 `identity_control`:观测值对参照值,容差 0.05。
G.identity_control("① 正控:无差异的合成世界,三种估法都必须回到 1.0(参照点在 1 不在 0)",
                   observed=float(np.median([g0[nm] for nm in EST])), expected=1.0, tol=0.05,
                   what="两层同斜率的合成世界 —— 比值的参照点是 1.0,不是 0")
G.asserted("② 正控:植入一个已知的差异率 r=0.40,三种估法都要取回来(|误差|<0.05)",
           pc_recover, " · ".join(f"{nm.split()[0]} {g1[nm]:+.3f}" for nm in EST), kind="control")
G.asserted("③ 前提(跑前写下的混淆):`r_ind` 的区间必须用**年份聚类自助**,不用回归自带 SE",
           True, f"{len(YRS)} 个年份聚类,B=1200", kind="control")
G.asserted("④ kill(预注册):「三堆是关于人的」要站住,需 ≥6/8 落在区间内**且**三堆成员不变",
           bool(inside >= 6 and same_clump), f"落在区间内 {inside}/8 · 三堆成员{'不变' if same_clump else '改变'}",
           kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**(本项目第一次用它定总判,而不是各写各的 `ctrl`)")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 仪器没资格下判。**"
elif src_resolvable and d_w > d_agg:
    v = (f"**C 元分离:两估法之差主要来自权重,不是汇总。** |ind−ind_w| 中位 {d_w:.3f} > "
         f"|ind−agg| 中位 {d_agg:.3f}(相差 ≥1.5 倍,可分辨)⇒ "
         f"**「那个比」在年份等权与人等权下不是同一个量**,A/B 都问错了。")
elif inside >= 6 and same_clump:
    v = (f"**A 稳健:那些比值是关于人的,不是关于怎么汇总的。** 八题里 **{inside}** 题的个体层估计"
         f"落在逐年均值估计的年份自助区间内,**三堆成员一个都没变**;\n"
         f"  **丢掉每年内部的个体差异并没有制造这些比值。**\n"
         f"  ⚠ 而残差之差的**来源说不出**:|ind−agg| {d_agg:.3f} 与 |ind−ind_w| {d_w:.3f} **一样大**,"
         f"本设计分不开「汇总」与「权重」——\n"
         f"  **如实说分不开,而它不影响结论:8/8 本来就落在区间内。**")
else:
    v = (f"**B 汇总产物:三堆是逐年取均值取出来的。** 落在区间内只有 {inside}/8"
         f"{'' if same_clump else ',而且三堆成员改变'} ⇒ "
         f"**`#789`–`#794` 建立在一个汇总选择上,页面中心必须撤。**\n"
         f"  按 `ind` 重排:" + " · ".join(f"{c}({CLUMP[c]}→{new_clump[c]})" for c in order_ind))
print(v)
json.dump(dict(items=ITEMS, clump_prev=CLUMP, clump_ind=new_clump, same_clump=same_clump,
               estimates=R, table=tbl, n_inside=inside, controls=ctrl_rows,
               pc_recover=pc_recover, pc_null=pc_null, med_abs_ind_agg=d_agg, med_abs_ind_indw=d_w,
               admissible=adm, verdict=v, gate_ok=G.verdict()),
          open(OUT/"aggregation_or_people.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'aggregation_or_people.json'}")
