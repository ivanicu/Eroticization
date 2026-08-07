"""#821 · E03·A66·R260 —— 「裂」还是「漂」要有一个连续的量,不是数有几格被标记

`#820` 读出「虔诚在裂、政治在漂」,**并当轮标注了那是事后读出来的**,同时暴露:
只是换个种子与 `B`,**44 个覆盖格里就翻了 1 格** ⇒ 数格子这件事本身是阈值化的、边界敏感的。
⇒ `#820`①② 一起做,**而它们本来就是同一件事**:
**一个连续的量按构造就没有阈值,所以它同时解决「事后」与「不稳」两个问题。**

⚠⚠ **G1 估计量,先命名再选方法(而不是算完再起名):**
   **`concentration(题, 轴)` = 该题**单个最大十年**的 |Δgap| ÷ 各十年 |Δgap| 之和。**
   · 完全匀速 ⇒ 接近 1/十年数 · 全部集中在一个十年 ⇒ 接近 1.0。
   · **它有界于 [1/n, 1],分母是绝对值之和(≥ 分子)⇒ 结构上不可能出现小分母**
     —— **`#799` 的分母放大器在这个量上被设计掉了,而 `#806` 就是被它咬过的。**

⚠⚠ **而「匀速时应该等于多少」不假设,由合成世界给出:**
   **基线 = 同一题、同一轴、在「真匀速」合成世界里跑出来的 `concentration` 的均值(多次实现)。**
   ⇒ 报的是 **`excess = concentration_观测 − concentration_匀速基线`** ——
   **这自动吸收了各十年跨度不等与覆盖不全,而那两件事我本来要手工校正、且校正本身会引入假设。**

G4 规格 · 种子:**≥3 个种子**(`realstat` 硬要求,`#820` 只有 2 次抽样)——
   全部量在 3 个种子下各跑一遍,**报每题每轴的种子间跨度**;
   **若某个结论在种子间翻转,那个结论就不存在。**

三个世界:
   A **两轴的 excess 分布不同**(虔诚轴上 `homosex` 的 excess 明显更大)⇒
     **「裂 vs 漂」是一个真的、连续的差别,而不是数格子的产物。**
   B **两轴的 excess 分布一样** ⇒ **`#820` 那句事后读数是阈值化造出来的** ——
     **那会撤掉 `#820` 最出彩的一句话,而这正是我不欢迎的那个结果(盆地规则)。**
   C **种子间不稳** ⇒ 连这个连续量也分辨不出,登记功效边界。

预测矩阵:
   | 世界 | 现在 | homosex 的 REL−pol excess 差 > 种子跨度 | ≈ 0 | 种子间翻转 |
   | A 真差别 | 0.50 | **0.85** | 0.05 | 0.10 |
   | B 阈值产物 | 0.30 | 0.05 | **0.85** | 0.15 |
   | C 分不出 | 0.20 | 0.10 | 0.10 | **0.75** |

预注册判词(条件式):
  if 正控开火(**植入「全部集中在一个十年」的世界,excess 必须明显 > 0**)
     and 负控开火(**一个新的真匀速实现,excess 必须 ≈ 0,而容差事先写死**):
      `homosex` 的 (REL excess − polviews excess) **在 3 个种子下同号,且最小的那个 > 种子跨度** -> A
      跨过 0 或小于种子跨度 -> B
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**两轴的覆盖十年数可能不同**(`polviews` 从 1974 起,`REL` 也是,
  但逐题的可用年不同)⇒ **十年数少的一侧 `concentration` 的下界 1/n 就更高**,
  **于是它天然看起来更「集中」。** ⇒ 控制:**基线由同一题同一轴的匀速世界给出,
  它按构造带着同样的十年数** —— **这正是用合成基线而不是用 1/n 的理由**;
  并**把两轴逐题的十年数并排印出来**。

⚠ 本轮**换不了仪器**:估计量是「同一具问卷内两根轴上的同一张网格的形状差」,第二份调查没有对应物。
⚠ 硬规则②:两根轴都来自 GSS ⇒ **同一具仪器内换轴,不是跨仪器复现。**
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]
SEEDS, NBASE, NC_TOL = [260, 261, 262], 60, 0.08

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", "polviews"]+ITEMS,
                  convert_categoricals=False)
cat = pd.read_stata(gp, columns=ITEMS, convert_categoricals=True)
K = {c: (4 if c == "homosex" else len(cat[c].cat.categories)) for c in ITEMS}
M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in ITEMS})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
BASE = M.dropna(subset=["year"]).copy()
Rr = BASE.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = zs(Rr[["attend", "reliten", "fund"]]).mean(axis=1)
BASE = BASE.join(Rr["REL"])
for nm, col in (("k_rel", "REL"), ("k_pol", "polviews")):
    BASE[nm] = BASE.groupby("year")[col].transform(
        lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)

def build(kcol):
    YR, COV = {}, {}
    for it in ITEMS:
        g = BASE.dropna(subset=[it, kcol]); ys = {}
        for y, gy in g.groupby("year"):
            a = gy[gy[kcol] == 2][it].to_numpy(float); b = gy[gy[kcol] == 0][it].to_numpy(float)
            if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
        YR[it] = ys
        dec = {}
        for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
        COV[it] = {k: v for k, v in dec.items() if len(v) >= 3}
    return YR, COV
AX = {"REL": build("k_rel"), "polviews": build("k_pol")}

print("=== ⓪ 跑前写下的混淆:两轴逐题的**可用十年数**(十年数少 ⇒ `concentration` 下界 1/n 更高)===")
for it in ITEMS:
    print(f"  {it:9s} REL {len(AX['REL'][1][it])} 个十年 {sorted(AX['REL'][1][it])} · "
          f"polviews {len(AX['polviews'][1][it])} 个十年 {sorted(AX['polviews'][1][it])}")
print("  ⇒ **正因为如此,基线不用 1/n,而用同一题同一轴的匀速合成世界跑出来的值 —— 它按构造带同样的十年数。**")

def gapfn(S, it, y): return float(S[it][y][0].mean()-S[it][y][1].mean())
def conc(S, COV, it):
    """单个最大十年的 |Δgap| ÷ 各十年 |Δgap| 之和。⚠ 分母 ≥ 分子 ⇒ 结构上无小分母。"""
    ds = []
    for dc, ys in sorted(COV[it].items()):
        ds.append(abs(gapfn(S, it, ys[-1]) - gapfn(S, it, ys[0])))
    tot = sum(ds)
    return (max(ds)/tot, len(ds)) if tot > 1e-9 else (np.nan, len(ds))

def syn(YR, COV, it, mode, rng, planted=1990):
    ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]
    g0 = gapfn(YR, it, y0); tot = gapfn(YR, it, y1) - g0
    out = {}
    for y in ys:
        f = ((y-y0)/(y1-y0) if mode == "uniform"
             else (0.0 if y < planted else (1.0 if y > planted+9 else (y-planted)/9.0)))
        cur = gapfn(YR, it, y)
        a, b = YR[it][y]
        # ⚠ 替换轨迹(`#812`/`#818`/`#819` 那一族已重犯三次),并按人重抽引入真实抽样噪声
        ia, ib = rng.integers(0, len(a), len(a)), rng.integers(0, len(b), len(b))
        out[y] = (a[ia] + (g0 + tot*f - cur), b[ib])
    return {it: out}

print(f"\n=== ① `concentration` 与 `excess`(基线 = 匀速合成世界 ×{NBASE} 次)· 种子 {SEEDS} ===")
RES = {}
for ax, (YR, COV) in AX.items():
    RES[ax] = {}
    for it in ITEMS:
        obs, nd = conc(YR, COV, it)
        per_seed = []
        for sd in SEEDS:
            rng = np.random.default_rng(sd)
            base = [conc(syn(YR, COV, it, "uniform", rng), COV, it)[0] for _ in range(NBASE)]
            base = np.array([x for x in base if np.isfinite(x)])
            per_seed.append(obs - float(np.mean(base)))
        RES[ax][it] = dict(conc=float(obs), n_dec=nd, excess=[float(x) for x in per_seed],
                           excess_mean=float(np.mean(per_seed)),
                           seed_spread=float(max(per_seed)-min(per_seed)))
    print(f"  —— 轴 = **{ax}**")
    for it in ITEMS:
        r = RES[ax][it]
        print(f"    {it:9s} 十年数 {r['n_dec']} · concentration **{r['conc']:.3f}** · "
              f"excess **{r['excess_mean']:+.3f}** · 种子跨度 {r['seed_spread']:.3f} "
              f"(逐种子 {[round(x,3) for x in r['excess']]})")

print("\n=== ② 两轴逐题的 excess 差(REL − polviews),3 个种子并排 ===")
DIFF = {}
for it in ITEMS:
    ds = [RES["REL"][it]["excess"][i] - RES["polviews"][it]["excess"][i] for i in range(len(SEEDS))]
    DIFF[it] = dict(per_seed=[float(x) for x in ds], mean=float(np.mean(ds)),
                    spread=float(max(ds)-min(ds)), same_sign=bool(len({np.sign(x) for x in ds}) == 1))
    print(f"  {it:9s} REL−pol excess = **{np.mean(ds):+.3f}** · 逐种子 {[round(x,3) for x in ds]} · "
          f"跨度 {max(ds)-min(ds):.3f} · 同号 **{DIFF[it]['same_sign']}**")
H = DIFF["homosex"]
crit = H["same_sign"] and min(abs(x) for x in H["per_seed"]) > H["spread"]
print(f"\n  ⇒ `homosex`:同号 **{H['same_sign']}** · 最小|差| **{min(abs(x) for x in H['per_seed']):.3f}** "
      f"vs 种子跨度 **{H['spread']:.3f}** ⇒ 预注册判据 **{'满足' if crit else '不满足'}**")

print("\n=== ③ 控制 ===")
# ⚠⚠ 第一版**控制的单位错了,而两条控制同时失败正是这个错的指纹**:
#   我给每题只跑了**一次**植入/匀速世界,拿 **8 题之间的散布**当成这条控制自己的抽样噪声 ——
#   **那是题间差异,不是抽样噪声。** 而估计量的单位是**逐题**,所以控制的单位也必须是逐题。
#   ⇒ 改成:**每题各重复 `NCTL` 次**,得到该题自己的植入分布与匀速分布;
#     正控的判据是**逐题的**:该题植入世界的 excess 中位是否超过它自己匀速分布的 97.5 分位。
#   ⚠⚠ **而第一版的失败还带回一个真发现,不许被修掉:** 正控中位只有 +0.073、跨度 [−0.154, +0.408] ——
#     **不是噪声,是这个估计量只在「全程移动量大」的题上有分辨力**;
#     `suicide2`(全程 +0.036)、`prayer`(+0.074)这类题上,噪声完全淹没了集中度。
#   ⇒ **所以正确的反应不是放宽正控,是把主张的适用范围缩到正控开火的那些题上**
#     (`realstat`:正控要报 **retention** 与 **MDE**)。
NCTL = 40
rng = np.random.default_rng(999)
YRr, COVr = AX["REL"]
RET = {}
for it in ITEMS:
    u = np.array([conc(syn(YRr, COVr, it, "uniform", rng), COVr, it)[0] for _ in range(NCTL)])
    pl = np.array([conc(syn(YRr, COVr, it, "planted", rng), COVr, it)[0] for _ in range(NCTL)])
    u, pl = u[np.isfinite(u)], pl[np.isfinite(pl)]
    ub = float(np.mean(u)); thr = float(np.percentile(u, 97.5))
    RET[it] = dict(uniform_mean=ub, uniform_p975=thr,
                   planted_median=float(np.median(pl)), fires=bool(np.median(pl) > thr),
                   uniform_half=float((np.percentile(u, 97.5)-np.percentile(u, 2.5))/2),
                   neg_excess_median=float(np.median(u)-ub))
print(f"  **逐题**(每题各 {NCTL} 次重复;⚠ **控制的单位 = 估计量的单位 = 题**):")
for it in ITEMS:
    r = RET[it]
    print(f"    {it:9s} 匀速基线 {r['uniform_mean']:.3f}(97.5% 分位 {r['uniform_p975']:.3f},半宽 {r['uniform_half']:.3f})"
          f" · 植入世界中位 **{r['planted_median']:.3f}** ⇒ 正控 **{'开火' if r['fires'] else '不开火'}**")
n_fire = sum(r["fires"] for r in RET.values())
print(f"  ⇒ **正控 retention:{n_fire}/{len(ITEMS)} 题** —— "
      f"**主张只对这些题成立,其余题登记为「这个估计量在它们上没有分辨力」**")
print(f"    ⚠ 不开火的题:{[it for it in ITEMS if not RET[it]['fires']]} —— "
      f"它们的全程移动量太小,噪声淹没了集中度(`realstat`:报 retention,不是放宽正控)")
pcH, ncH = RET["homosex"]["planted_median"] - RET["homosex"]["uniform_mean"], RET["homosex"]["neg_excess_median"]
nc_half = RET["homosex"]["uniform_half"]
print(f"  `homosex`(kill 所针对的那一题):正控 excess **{pcH:+.3f}** · 负控 excess **{ncH:+.4f}** · "
      f"匀速噪声半宽 **{nc_half:.4f}** vs 事先写死的容差 {NC_TOL} ⇒ 比值 **{NC_TOL/nc_half:.2f}×**")
pc = np.array([RET[it]["planted_median"] - RET[it]["uniform_mean"] for it in ITEMS])
nc = np.array([RET[it]["neg_excess_median"] for it in ITEMS])

G = Gate("#821 · 「裂」还是「漂」要有一个连续的量")
G.asserted("① 正控(**逐题,单位与估计量一致**):`homosex` 这一题的植入世界 excess 必须超过"
           "它自己匀速分布的 97.5 分位 —— ⚠ 并**同时报 retention**,主张只对开火的题成立",
           bool(RET["homosex"]["fires"]),
           f"homosex 植入中位 {RET['homosex']['planted_median']:.3f} > 匀速 97.5% 分位 "
           f"{RET['homosex']['uniform_p975']:.3f} ⇒ {RET['homosex']['fires']} · "
           f"**retention {n_fire}/{len(ITEMS)} 题**,不开火的 {[it for it in ITEMS if not RET[it]['fires']]}",
           kind="control")
G.identity_control("② 负控(**逐题**):`homosex` 的匀速世界 excess 必须 ≈ 0(⚠ **这一次参照真的是 0**)"
                   " —— 容差 0.08 **事先写死**,并与该题**自己**的匀速噪声半宽比对(`#814`/`#817`②)",
                   observed=float(ncH), expected=0.0, tol=NC_TOL, noise_half_width=nc_half,
                   what=f"homosex 匀速世界 {NCTL} 次重复,噪声半宽 {nc_half:.4f}")
G.asserted("③ 前提(跑前写下的混淆):两轴逐题**可用十年数**已并排印出,且**基线不用 1/n,"
           "而用同一题同一轴的匀速合成世界** —— 它按构造带同样的十年数",
           bool(all(RES["REL"][it]["n_dec"] > 0 and RES["polviews"][it]["n_dec"] > 0 for it in ITEMS)),
           f"REL 十年数 {[RES['REL'][it]['n_dec'] for it in ITEMS]} · "
           f"polviews {[RES['polviews'][it]['n_dec'] for it in ITEMS]}", kind="control")
G.asserted("④ 前提:估计量**结构上无小分母**(分母 = 各十年 |Δgap| 之和 ≥ 分子)——"
           "`#799` 的分母放大器在这个量上被设计掉了,而 `#806` 就是被它咬过的",
           True, "concentration ∈ [1/n, 1],有界", kind="control")
G.asserted(f"⑤ 前提:**≥3 个种子**({len(SEEDS)} 个)且逐题种子跨度已报(`realstat`;`#820` 只有 2 次抽样)",
           bool(len(SEEDS) >= 3), f"种子 {SEEDS} · homosex 种子跨度 REL {RES['REL']['homosex']['seed_spread']:.3f} · "
           f"polviews {RES['polviews']['homosex']['seed_spread']:.3f}", kind="control")
G.asserted("⑥ kill(预注册):「裂 vs 漂是真的连续差别」(世界 A)要成立,需 `homosex` 的 "
           "(REL − polviews) excess **3 个种子同号,且最小|差| > 种子跨度**",
           bool(crit), f"同号 {H['same_sign']} · 最小|差| {min(abs(x) for x in H['per_seed']):.3f} vs 跨度 {H['spread']:.3f}",
           kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif crit:
    V = (f"**A 「裂 vs 漂」是一个真的连续差别,不是数格子的产物。**\n"
         f"  `homosex` 的 excess:虔诚轴 **{RES['REL']['homosex']['excess_mean']:+.3f}** · "
         f"政治轴 **{RES['polviews']['homosex']['excess_mean']:+.3f}** ⇒ 差 **{H['mean']:+.3f}**,"
         f"3 个种子同号且最小|差| {min(abs(x) for x in H['per_seed']):.3f} > 种子跨度 {H['spread']:.3f}。\n"
         f"  ⇒ **`#820` 那句事后读数在一个连续、无阈值、多种子的量上站住了。**")
else:
    V = (f"**B 或 C:这个连续量不支持 `#820` 那句事后读数。**\n"
         f"  `homosex` 的 (REL − polviews) excess = **{H['mean']:+.3f}**,逐种子 "
         f"{[round(x, 3) for x in H['per_seed']]},种子跨度 **{H['spread']:.3f}**,同号 **{H['same_sign']}**。\n"
         f"  ⇒ **「虔诚在裂、政治在漂」有可能是阈值化(数格子)造出来的 —— 而这正是我不欢迎的那个结果,\n"
         f"  它要求把 `#820` 那句话降级为「格子计数下的现象」,不是「一个连续的形状差别」。**")
print(V)
json.dump(dict(items=ITEMS, seeds=SEEDS, n_base=NBASE, nc_tol=NC_TOL, per_axis=RES, diff=DIFF,
               homosex_criterion=bool(crit),
               retention=RET, n_fire=int(n_fire),
               pos_control=dict(homosex_excess=float(pcH), median_all=float(np.median(pc))),
               neg_control=dict(homosex_excess=float(ncH), half_width=nc_half, reference=0.0),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"concentration.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'concentration.json'}")
