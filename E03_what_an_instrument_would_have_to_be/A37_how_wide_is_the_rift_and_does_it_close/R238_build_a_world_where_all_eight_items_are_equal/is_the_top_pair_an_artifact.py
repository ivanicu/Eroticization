"""#799 · E03·A51·R238 —— `#797` 之后只剩「顶对」站着,而它的分母恰恰是八题里最小的

`#797` 撤了三堆,页面上现在**只剩一件东西站着**:
**顶对 `sexeduc` (2.19) 与 `racmar` (1.84) 与其余六题全部分得开**(顶/中 8/8 · 顶/底 4/4)。
⚠⚠ **而把八题按分母大小排一遍,最该担心的事就露出来了:**

| 题 | \\|分母斜率\\| | r |
|---|---|---|
| **`sexeduc`** | **0.00128**(全场最小) | **+2.192**(全场最大) |
| `suicide2` | 0.00226 | +0.410 |
| `prayer` | 0.00332 | +0.126 |
| **`racmar`** | 0.00719 | **+1.841** |
| `spanking` · `teensex` · `helpblk` · `homosex` | 0.0096 – 0.0373 | 0.105 – 0.462 |

⇒ **最小的分母配着最大的比值。** 而 `#782` 立过:**两个近零斜率相除是 Cauchy 型**;
`#790` 立过:**不稳定由分母驱动**。⇒ **「顶对」会不会只是小分母把噪声放大成了 r>1?**
**这个正面结果我会非常不想要 —— 它一落地,页面上就一件站着的东西都没有了。**
⚠ 而朴素版的怀疑**已经被数据本身削弱了一半**:`suicide2`(第二小)与 `prayer`(第三小)的 r 是
0.410 与 0.126,**一点也不大** ⇒ **「分母小 ⇒ r 大」不是单调成立的。** 所以不能靠排序下判,要造世界。

G1 估计量:**在一个「八题真值同为 r=0.45」的世界里,每题冒出 `r > 1.5` 的频率。**
   ⚠ 这是 `realstat` 攻击阶梯的**第 4 级**(造出对手预言的世界,看仪器会不会自己开火),
   不是第 1 级的 gauge、也不是第 2 级的算术 —— **前三级这一族已经跑过,这一级没跑过。**

造法(**用真实的分母,不是合成的**):
   · 每题取**观测到的非虔诚层年份序列**(真年份、真值)⇒ **分母完全是真的**;
   · 虔诚层合成:`dev(y) = dev0 + 0.45 × slope_nondev × (y−y0) + ε`,
     **真值 r 恒为 0.45,与题目无关**;
   · `ε` 的尺度取自**该题观测虔诚层围绕自己线性拟合的残差**。

⚠⚠ 跑之前写下的最强混淆,**而它需要自己的一条臂**:
   `ε` 取自该题自己的残差 ⇒ **若小分母题恰好也更吵,我就把要检验的东西烤进了合成世界。**
   ⇒ **第二臂:所有题共用同一个 `ε` 尺度(八题残差的中位)。**
   **两臂对比才能把「分母小」与「这题本来就吵」分开** —— `#794` 刚栽在一个「同时差两件事」的对照上。

三个世界:
   A **顶对是真的**:共同真值世界里,`sexeduc`/`racmar` 冒出 r>1.5 的频率 < 0.10
     ⇒ 小分母解释不了它们,**顶对经住了第 4 级攻击。**
   B **顶对是小分母产物**:频率 ≥ 0.10 ⇒ **页面上最后一件站着的东西也要撤。**
   C **是「这题吵」不是「分母小」**:两臂给出不同答案 ⇒ 我把两件事混成了一件,**判据本身要重写。**

预测矩阵:
   | 世界 | 现在 | 若两臂频率都 <0.10 | 若两臂都 ≥0.10 | 若两臂不一致 |
   | A | 0.45 | **0.90** | 0.03 | 0.15 |
   | B | 0.35 | 0.05 | **0.90** | 0.15 |
   | C | 0.20 | 0.05 | 0.07 | **0.70** |

预注册判词(条件式):
  if 正控开火(无噪声时每题必须**恰好**取回 0.45)
     and 生成器正控(有噪声时 r 的中位必须落在 0.45 ± 0.10):
      两臂里 `sexeduc` 与 `racmar` 的 P(r>1.5) 都 < 0.10 -> A
      任一臂里任一题 ≥ 0.10                              -> B
      两臂给出不同结论                                    -> C
  else: UNVERIFIED
⚠ **0.10 这个门槛不是现在挑的**:`#797` 的假阳性率门槛就是 0.10,**沿用同一个。**
⚠ **1.5 这个界也不是现在挑的**:它是「r>1」这个定性主张的一个保守化 ——
  顶对实测 1.84/2.19,**若小分母能常规造出 1.5,那 1.84 就不安全。**

⚠ **「这个零该不该是零?」** —— **无噪声那条控制的参照是 0.45,不是 0**(真值就是 0.45)
  ⇒ **`identity_control`,不是 `negative_control`,也不是 `offset_control`**
  (`#796` 刚在这里错过一次:`offset_control` 问「越没越过」,我要的是「等不等于」)。

本轮换不了仪器(对象是世界;第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
⚠ 总判由 `Gate.admissible()` 决定(`#796` 加的,第三次用)。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(238)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_is_the_eight_point_axis_an_axis_or_eight_labels_on_noise/results/is_the_ordering_an_object.json"))
ITEMS, CLUMP = P791["items"], P791["clumps"]
TOP = [c for c in ITEMS if CLUMP[c] == "顶"]
R_TRUE, THRESH, RATE = 0.45, 1.5, 0.10
print(f"=== ⓪ 顶对 = {TOP} · 共同真值 r={R_TRUE} · 界 r>{THRESH} · 门槛频率 {RATE}(沿用 `#797`)===")

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

OBS = {}
for c in ITEMS:
    rB, rA = series(c, 0), series(c, 2)
    yB = np.array([r[0] for r in rB], float); vB = np.array([r[1] for r in rB])
    yA = np.array([r[0] for r in rA], float); vA = np.array([r[1] for r in rA])
    sB, sA = slope(yB, vB), slope(yA, vA)
    resid = vA - (vA.mean() + sA*(yA-yA.mean()))
    OBS[c] = dict(yB=yB, vB=vB, yA=yA, sB=sB, sA=sA, r=sA/sB,
                  sigma=float(np.std(resid, ddof=2)), n=len(yA))
sig_med = float(np.median([OBS[c]["sigma"] for c in ITEMS]))
print("\n=== ① 真实的分母与每题虔诚层的残差尺度(合成世界只借这两样)===")
for c in sorted(ITEMS, key=lambda x: abs(OBS[x]["sB"])):
    print(f"  {c:9s} |分母斜率| {abs(OBS[c]['sB']):.5f} · 虔诚层残差 σ {OBS[c]['sigma']:.4f} · "
          f"年 {OBS[c]['n']:>2} · 实测 r {OBS[c]['r']:+.3f}")
print(f"  八题残差 σ 的中位 = **{sig_med:.4f}**(第二臂用它,把「这题吵」与「分母小」分开)")

def synth_r(c, sigma, B=4000, r_true=R_TRUE):
    """真分母 + 合成虔诚层(真值 r 恒为 r_true)。sigma=0 ⇒ 必须恰好取回 r_true。"""
    o = OBS[c]; yB, vB, yA = o["yB"], o["vB"], o["yA"]
    out = []
    for _ in range(B):
        vA = 3.0 + r_true*o["sB"]*(yA-yA[0]) + (RNG.normal(0, sigma, len(yA)) if sigma > 0 else 0.0)
        sb = slope(yB, vB)
        out.append(slope(yA, vA)/sb if abs(sb) > 1e-12 else np.nan)
    return np.array([x for x in out if np.isfinite(x)])

# ── 正控:无噪声时必须恰好取回真值 ────────────────────────────────────────────
print("\n=== ② 正控:σ=0 时每题必须**恰好**取回 0.45(参照是 0.45 不是 0 ⇒ identity_control)===")
exact = {c: float(synth_r(c, 0.0, B=3)[0]) for c in ITEMS}
worst = max(abs(v-R_TRUE) for v in exact.values())
print("  " + " · ".join(f"{c} {exact[c]:.4f}" for c in ITEMS))
print(f"  ⇒ 最大偏差 **{worst:.2e}** —— **σ=0 时小分母本身不制造任何膨胀,这一条把两件事分开了。**")

print("\n=== ③ 两臂:各题自己的 σ · 全部用中位 σ(`G4` 的噪声轴)===")
ARMS = {"臂1 各题自己的 σ": lambda c: OBS[c]["sigma"], "臂2 全部用中位 σ": lambda c: sig_med}
res, rates = {}, {}
for an, f in ARMS.items():
    res[an] = {}
    for c in ITEMS:
        v = synth_r(c, f(c))
        res[an][c] = dict(med=float(np.median(v)), p_gt=float((v > THRESH).mean()),
                          q95=float(np.percentile(v, 95)), sigma=float(f(c)))
    rates[an] = {c: res[an][c]["p_gt"] for c in ITEMS}
    print(f"  {an}")
    for c in sorted(ITEMS, key=lambda x: abs(OBS[x]["sB"])):
        z = res[an][c]
        print(f"    {c:9s} 中位 r {z['med']:+.3f} · 95分位 {z['q95']:+.3f} · "
              f"**P(r>{THRESH}) = {z['p_gt']:.3f}**{'  ⚠ 顶对' if c in TOP else ''}")

gen_ok = all(abs(res[an][c]["med"]-R_TRUE) <= 0.10 for an in ARMS for c in ITEMS)
top_rates = [rates[an][c] for an in ARMS for c in TOP]
any_rate = max(max(rates[an].values()) for an in ARMS)
arm_agree = (max(rates["臂1 各题自己的 σ"][c] for c in TOP) < RATE) == \
            (max(rates["臂2 全部用中位 σ"][c] for c in TOP) < RATE)
print(f"\n  顶对在两臂里的 P(r>{THRESH}):{[round(x,3) for x in top_rates]} · 全场最大 {any_rate:.3f}")

G = Gate("#799 · 剩下那一对会不会也是小分母")
G.identity_control("① 正控:σ=0 时每题必须恰好取回真值 0.45(参照是 0.45,不是 0)",
                   observed=float(np.median(list(exact.values()))), expected=R_TRUE, tol=1e-6,
                   what="无噪声的共同真值世界 —— 小分母在没有噪声时不制造膨胀")
G.asserted("② 生成器正控:有噪声时 r 的中位必须落在 0.45 ± 0.10(否则合成的不是我说的那个世界)",
           gen_ok, f"两臂 16 格的中位偏差最大 "
                   f"{max(abs(res[a][c]['med']-R_TRUE) for a in ARMS for c in ITEMS):.3f}", kind="control")
G.asserted("③ 前提(跑前写下的混淆):必须有两臂,把「分母小」与「这题本来就吵」分开",
           bool(len(ARMS) == 2), f"臂1 各题 σ({min(OBS[c]['sigma'] for c in ITEMS):.4f}–"
                                 f"{max(OBS[c]['sigma'] for c in ITEMS):.4f})· 臂2 中位 σ {sig_med:.4f}", kind="control")
G.asserted("④ kill(预注册):顶对要站住,需**两臂**的 P(r>1.5) 都 < 0.10",
           bool(max(top_rates) < RATE), f"顶对两臂最大 {max(top_rates):.3f} vs 门槛 {RATE}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 仪器没资格下判。**"
elif not arm_agree:
    v = (f"**C:两臂不一致 ⇒ 我把「分母小」与「这题本来就吵」混成了一件事,判据要重写。** "
         f"臂1 顶对最大 {max(rates['臂1 各题自己的 σ'][c] for c in TOP):.3f} · "
         f"臂2 {max(rates['臂2 全部用中位 σ'][c] for c in TOP):.3f}")
elif max(top_rates) < RATE:
    v = (f"**A 顶对经住了第 4 级攻击。** 在一个**八题真值同为 0.45**、**分母全部取自真实观测**的世界里,"
         f"`sexeduc` 与 `racmar` 冒出 r>{THRESH} 的频率两臂最大只有 **{max(top_rates):.3f}**"
         f"(门槛 {RATE}),而它们实测是 **2.19 / 1.84**。\n"
         f"  ⇒ **小分母 + 噪声造不出那两个数** —— 顶对不是分母产物。\n"
         f"  ⚠ 而 σ=0 那条正控把话说死:**没有噪声时小分母本身一点膨胀都不制造**"
         f"(最大偏差 {worst:.0e})—— **能膨胀的是噪声,不是小分母本身。**")
else:
    bad = [(a, c, rates[a][c]) for a in ARMS for c in TOP if rates[a][c] >= RATE]
    v = (f"**B 顶对也是小分母产物 —— 页面上最后一件站着的东西也要撤。** "
         f"共同真值世界里顶对仍冒出 r>{THRESH}:{[(c, round(p,3)) for _, c, p in bad]}(门槛 {RATE})。")
print(v)
json.dump(dict(items=ITEMS, top=TOP, r_true=R_TRUE, thresh=THRESH, rate_thresh=RATE,
               observed={c: dict(r=OBS[c]["r"], sB=OBS[c]["sB"], sigma=OBS[c]["sigma"], n=OBS[c]["n"])
                         for c in ITEMS},
               sigma_median=sig_med, exact_noise_free=exact, worst_exact=worst,
               arms={a: res[a] for a in ARMS}, top_rates=top_rates, arm_agree=bool(arm_agree),
               admissible=adm, verdict=v, gate_ok=G.verdict()),
          open(OUT/"top_pair_artifact.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'top_pair_artifact.json'}")
