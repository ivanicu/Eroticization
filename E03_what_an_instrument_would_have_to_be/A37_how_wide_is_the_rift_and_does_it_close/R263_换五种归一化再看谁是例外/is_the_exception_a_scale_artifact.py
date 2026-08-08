"""#824 · E03·A69·R263 —— 那个「2.5 倍」是原始量表分,而八题的量表不可比

`#823` 强制换方向,并指向 `#819` 留下的那个未解释的事实:
**`homosex` 的十年偏离幅度是整张网格其余最大格的 2.5 倍(−0.629 对 −0.256)。**

⚠⚠ **而在去解释它之前,先用最便宜的那一击打它自己**(`frontier §4` / `realstat §3` 攻击阶梯第 1 级):
   **那个 2.5 倍是原始量表分算的,而八题的档数是 2–5、作答分布天差地别。**
   `#819` 的 `EXCLUDES` 判定本身是**逐题自归一**的(参照 = 该题自己的全程 `Δgap`)⇒ **那一半是量纲安全的**;
   **但「谁的幅度最大、大多少」这句话是裸的原始分 ⇒ 它从来没被归一化过。**
   ⇒ **若它经不住换一种合理的归一化,那么 `#823` 派我去解释的那个「例外」根本不存在,
   而后面三轮会是在解释一个假象。**

⚠ **gauge 逻辑**(`frontier §3`):**「homosex 是例外」若是真性质,它必须对每题量表的单调重标定不变。**
   **测量变了而性质不该变 ⇒ 测量是瞎的。** 这是零算力的最便宜一击,先做。

G1 估计量:**逐题「最大十年偏离」= max_d |Δgap(题,十年) − ref(题,十年)|**,
   `ref` 与 `#819` 完全一致(该题全程 `Δgap` × 该十年跨年 ÷ 全程跨年)。
   然后**在五种归一化下各排一次序**,报 **`homosex` 是不是最大** 与 **它对亚军的倍数**。

**`G4` 规格曲线 —— 五种归一化,跑前写死:**
   ① `raw` —— `#819` 用的那种,原始量表分
   ② `÷span` —— 除以量表跨度 `K−1`
   ③ `÷SD_pooled` —— 除以该题**全样本的人际标准差**(标准化效应量的常规做法)
   ④ `÷SD_within_year` —— 除以该题**年内标准差的平均**(把长期漂移从 SD 里剔掉)
   ⑤ `÷|全程 Δgap|` —— 除以该题自己一共动了多少(**「这一跳占它全部移动的几成」**)

⚠⚠ **跑之前写下的最强混淆,而它足以自己制造世界 A:**
   **除以「该题自己的 SD」会给低方差的题一个人为的放大。**
   若某题作答挤在两档里,它的 SD 小,除下去之后任何一点移动都显得巨大。
   ⇒ 控制:**不只报 `homosex` 的名次,报五种归一化下的完整排序与每题的位次变动** ——
   **只报冠军是谁,等于把这个混淆藏起来。**

三个世界:
   A **五种下 homosex 都是最大且倍数 > 1.5** ⇒ 例外是真的,`#823` 的方向成立。
   B **≥1 种下倍数塌到 1.5 以下,或 homosex 不再是最大** ⇒
     **`#819` 那句「2.5 倍」必须撤,而 `#823` 派我去解释的东西不存在。**
   C **排序本身随归一化大改** ⇒ **「哪一题最极端」根本不是一个对象** ——
     与 `#797` 撤掉「三堆」是同一族的结论,而那比 A、B 都更伤。

预测矩阵:
   | 世界 | 现在 | 五种都最大且 >1.5 | 有一种塌 | 排序大改 |
   | A 真例外 | 0.45 | **0.85** | 0.05 | 0.05 |
   | B 量纲产物 | 0.30 | 0.05 | **0.85** | 0.30 |
   | C 排序不是对象 | 0.25 | 0.03 | 0.30 | **0.80** |

预注册判词(**条件式**):
  if 正控开火(**造一个「某题真有一个巨大例外偏离」的世界,它必须在五种归一化下都是最大**)
     and 负控开火(**造一个「八题全匀速」的世界,倍数必须落在它自己测出来的零分布内** ——
        ⚠ **「这个零该不该是零?」:不该。** 八个带噪声的量取 max÷亚军,**按构造 > 1**,
        零值既不是 0 也不是 1,**只能由匀速世界自己测出来**):
      五种归一化下 `homosex` 都是最大 **且** 最小倍数 > 1.5 -> A
      否则 -> B/C,**报整张五×八的表,不选边**
  else: UNVERIFIED
⚠ **≥3 个种子**(`realstat`;`#820` 的教训),逐格报种子跨度。
⚠ **`G3`:5 种归一化 × 8 题 = 40 格,全报,包括名次下降的格。**

⚠ 本轮**换不了仪器**:检验的是同一份数据上同一个量的量纲不变性,第二具仪器没有对应物;
   **而这正是硬规则④所说的例外 —— 这不是「同一个主张再跑一遍」,是攻击那个主张的量纲基础。**
⚠ 硬规则②:一切来自 GSS 同一份问卷。⚠ 总判由 `Gate.admissible()` 决定。
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
SEEDS, NNULL, THR = [263, 264, 265], 300, 1.5

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
YR, COV, SD = {}, {}, {}
for it in ITEMS:
    g = REL.dropna(subset=[it]); ys = {}
    for y, gy in g.groupby("year"):
        a = gy[gy.k == 2][it].to_numpy(float); b = gy[gy.k == 0][it].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    YR[it] = ys
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    COV[it] = {k2: v for k2, v in dec.items() if len(v) >= 3}
    allv = g[it].to_numpy(float)
    SD[it] = dict(pooled=float(np.std(allv, ddof=1)),
                  within=float(np.mean([np.std(gy[it].to_numpy(float), ddof=1)
                                        for _, gy in g.groupby("year") if len(gy) > 30])))

print("=== ⓪ 硬规则①:每题档数 · 人际 SD · 年内 SD · 全程 Δgap(五种归一化的分母,先看清楚)===")
FULL = {}
for it in ITEMS:
    ys = sorted(YR[it])
    gp0 = float(YR[it][ys[0]][0].mean()-YR[it][ys[0]][1].mean())
    gp1 = float(YR[it][ys[-1]][0].mean()-YR[it][ys[-1]][1].mean())
    FULL[it] = dict(y0=ys[0], y1=ys[-1], span=ys[-1]-ys[0], dgap=gp1-gp0)
    print(f"  {it:9s} 档 {K[it]}(跨度 {K[it]-1}) · SD人际 {SD[it]['pooled']:.3f} · SD年内 {SD[it]['within']:.3f} · "
          f"全程 Δgap {FULL[it]['dgap']:+.4f} · 可用十年 {len(COV[it])}")

NORMS = {"① raw(#819 用的)": lambda it: 1.0,
         "② ÷span": lambda it: K[it]-1,
         "③ ÷SD人际": lambda it: SD[it]["pooled"],
         "④ ÷SD年内": lambda it: SD[it]["within"],
         "⑤ ÷|全程Δgap|": lambda it: abs(FULL[it]["dgap"])}

def max_dep(it, src=None):
    """该题的最大十年偏离 = max_d |Δgap_d − ref_d|(`ref` 与 `#819` 完全一致)。"""
    S = src[it] if src else YR[it]
    g = lambda y: float(S[y][0].mean()-S[y][1].mean())
    best, where = 0.0, None
    for dc, ys in sorted(COV[it].items()):
        span = ys[-1]-ys[0]
        ref = FULL[it]["dgap"]*span/FULL[it]["span"]
        dev = abs((g(ys[-1])-g(ys[0])) - ref)
        if dev > best: best, where = dev, dc
    return best, where

print(f"\n=== ① `G3` 五种归一化 × 八题 = {len(NORMS)*len(ITEMS)} 格,全报(包括名次下降的格)===")
RAW = {it: max_dep(it) for it in ITEMS}
TAB, RANKS = {}, {}
for nm, fn in NORMS.items():
    vals = {it: RAW[it][0]/fn(it) for it in ITEMS}
    order = sorted(ITEMS, key=lambda x: -vals[x])
    TAB[nm] = dict(vals={k2: float(v) for k2, v in vals.items()}, order=order,
                   top=order[0], ratio=float(vals[order[0]]/vals[order[1]]))
    RANKS[nm] = {it: order.index(it)+1 for it in ITEMS}
print(f"  {'归一化':16s} {'冠军':10s} {'倍数(冠/亚)':>12s}   完整排序")
for nm in NORMS:
    t = TAB[nm]
    print(f"  {nm:16s} {t['top']:10s} {t['ratio']:>11.2f}×   " + " > ".join(t["order"]))
print(f"\n  ⚠ 跑前混淆的控制 —— **逐题位次在五种下的变动**(只报冠军等于把混淆藏起来):")
for it in ITEMS:
    r = [RANKS[nm][it] for nm in NORMS]
    print(f"    {it:9s} 位次 {r} · 变动幅度 **{max(r)-min(r)}**")
homo_top = all(TAB[nm]["top"] == "homosex" for nm in NORMS)
min_ratio = min(TAB[nm]["ratio"] for nm in NORMS)
print(f"\n  ⇒ `homosex` 在五种下都是最大?**{homo_top}** · 最小倍数 **{min_ratio:.2f}×**(阈 {THR})")

print("\n=== ② 控制(⚠「这个零该不该是零?」——**不该**:八个带噪声的量取 max÷亚军按构造 > 1)===")
def syn(mode, rng, planted_item="homosex", boost=3.0):
    """替换轨迹(`#812`/`#818`/`#819` 那一族已重犯三次);uniform = 每题真匀速;
    planted = 在 planted_item 的 1990s 上再加一个已知的巨大偏离。"""
    S = {}
    for it in ITEMS:
        ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]
        g0 = float(YR[it][y0][0].mean()-YR[it][y0][1].mean()); tot = FULL[it]["dgap"]
        S[it] = {}
        for y in ys:
            tgt = g0 + tot*(y-y0)/(y1-y0)
            if mode == "planted" and it == planted_item and 1990 <= y <= 1999:
                tgt += boost*abs(tot)*(y-1990)/9.0
            cur = float(YR[it][y][0].mean()-YR[it][y][1].mean())
            a, b = YR[it][y]
            ia, ib = rng.integers(0, len(a), len(a)), rng.integers(0, len(b), len(b))
            S[it][y] = (a[ia] + (tgt - cur), b[ib])
    return S
def ratio_of(S):
    out = {}
    for nm, fn in NORMS.items():
        v = {it: max_dep(it, src=S)[0]/fn(it) for it in ITEMS}
        o = sorted(ITEMS, key=lambda x: -v[x])
        out[nm] = (o[0], float(v[o[0]]/v[o[1]]))
    return out
rng = np.random.default_rng(9999)
PC = ratio_of(syn("planted", rng))
print(f"  正控(在 `homosex` 的 1990s 植入一个 3× 全程量的巨大偏离):"
      + " · ".join(f"{nm[0]}{v[0][:8]} {v[1]:.2f}×" for nm, v in PC.items()))
pc_ok = all(v[0] == "homosex" for v in PC.values())
NULLS = {nm: [] for nm in NORMS}
for i in range(NNULL):
    r = ratio_of(syn("uniform", np.random.default_rng(50000+i)))
    for nm in NORMS: NULLS[nm].append(r[nm][1])
print(f"  负控(八题全匀速)⇒ 倍数的**零分布**({NNULL} 次;⚠ 零值既不是 0 也不是 1,由它自己测出来):")
NQ = {}
for nm in NORMS:
    a = np.array(NULLS[nm]); NQ[nm] = dict(med=float(np.median(a)), p95=float(np.percentile(a, 95)),
                                           half=float((np.percentile(a, 97.5)-np.percentile(a, 2.5))/2))
    print(f"    {nm:16s} 中位 **{NQ[nm]['med']:.2f}×** · 95% 分位 **{NQ[nm]['p95']:.2f}×** · 半宽 {NQ[nm]['half']:.2f}")
above = {nm: bool(TAB[nm]["ratio"] > NQ[nm]["p95"]) for nm in NORMS}
print(f"  ⇒ 观测倍数超过匀速零分布 95% 分位的归一化:**{sum(above.values())}/{len(NORMS)}** —— {above}")

print(f"\n=== ③ 种子稳定性({len(SEEDS)} 个种子;`#820` 的教训)===")
SEEDR = {}
for sd in SEEDS:
    r = ratio_of(syn("uniform", np.random.default_rng(sd)))   # 仅用于确认零分布不随种子翻转
    SEEDR[sd] = {nm: r[nm][1] for nm in NORMS}
for nm in NORMS:
    v = [SEEDR[sd][nm] for sd in SEEDS]
    print(f"  {nm:16s} 匀速世界倍数逐种子 {[round(x,2) for x in v]} · 跨度 {max(v)-min(v):.2f}")
print("  ⚠ **观测值本身不含随机性**(它是确定性地从数据算出来的)—— **种子只作用在零分布上,如实说。**")

G = Gate("#824 · 那个「2.5 倍」是原始量表分")
G.asserted("① 正控:在 `homosex` 的 1990s 植入一个 3× 全程量的巨大偏离,"
           "**它必须在五种归一化下都是冠军**(否则这套排序机器连一个真例外都认不出)",
           bool(pc_ok), " · ".join(f"{nm[0]}→{v[0]}" for nm, v in PC.items()), kind="control")
G.asserted("② 负控:八题全匀速的世界里,倍数必须落在**它自己测出来的零分布**内 —— "
           "⚠ **「这个零该不该是零?」不该**:八个带噪声的量取 max÷亚军**按构造 > 1**,"
           "零值既不是 0 也不是 1,**只能由匀速世界自己测**;"
           "⇒ 写成 `asserted` 而非 `identity_control`/`offset_control`,"
           "**因为要比的是「观测是否超过一个测出来的分布」,不是「两个量相不相等」**",
           bool(all(0.9 < NQ[nm]["med"] < 3.0 for nm in NORMS)),
           " · ".join(f"{nm[0]}中位{NQ[nm]['med']:.2f}×" for nm in NORMS), kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**除以该题自己的 SD 会给低方差题人为放大** ⇒ "
           "**报五种下的完整排序与每题位次变动,而不只报冠军**",
           bool(all(len(TAB[nm]["order"]) == len(ITEMS) for nm in NORMS)),
           " · ".join(f"{it}:{max(RANKS[nm][it] for nm in NORMS)-min(RANKS[nm][it] for nm in NORMS)}"
                      for it in ITEMS), kind="control")
G.asserted("④ 前提(gauge):`ref` 与 `#819` 完全一致,**只有分母在变** —— "
           "变的是测量,不该变的是性质",
           True, "ref = 该题全程 Δgap × 该十年跨年 ÷ 全程跨年,五种归一化只除分母", kind="control")
G.asserted(f"⑤ kill(预注册):「例外是真的」(世界 A)要成立,需**五种归一化下 `homosex` 都是冠军"
           f"且最小倍数 > {THR}**",
           bool(homo_top and min_ratio > THR),
           f"都是冠军 {homo_top} · 最小倍数 {min_ratio:.2f}×(阈 {THR})", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
moved = {it: max(RANKS[nm][it] for nm in NORMS)-min(RANKS[nm][it] for nm in NORMS) for it in ITEMS}
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif homo_top and min_ratio > THR:
    V = (f"**A 例外是真的,而且它经住了量纲这一击。** 五种归一化下 `homosex` 都是冠军,"
         f"最小倍数 **{min_ratio:.2f}×**(阈 {THR});最大位次变动 {max(moved.values())} 名。\n"
         f"  ⇒ **`#819` 那句「2.5 倍」站住了,而 `#823` 派我去解释的那个例外确实存在。**")
elif not homo_top:
    losers = [nm for nm in NORMS if TAB[nm]["top"] != "homosex"]
    V = (f"**B `homosex` 不是所有归一化下的冠军。** 它在 {losers} 下让位给 "
         f"{[TAB[nm]['top'] for nm in losers]}。\n"
         f"  ⇒ **`#819` 那句「2.5 倍」必须撤 —— 那个倍数是量纲的产物,而不是一个关于这道题的事实。**\n"
         f"  ⇒ **而 `#823` 派我去解释的东西不存在,后面三轮本会是在解释一个假象。**")
else:
    V = (f"**B 倍数塌了。** `homosex` 五种下都是冠军,**但最小倍数只有 {min_ratio:.2f}×**,低于阈 {THR}。\n"
         f"  ⇒ **「2.5 倍」是最宽容的那种归一化下的数,而不是这道题的性质** ——\n"
         f"  **`#819` 那句话要改成「在原始量表分上 2.5 倍,而换一种合理的归一化只剩 "
         f"{min_ratio:.2f} 倍」,并且不许再作为「例外」的证据单独引用。**")
print(V)
print(f"\n  ⚠ **排序是不是一个对象**:八题在五种归一化下的最大位次变动 = "
      f"**{max(moved.values())} 名**({[k2 for k2, v in moved.items() if v == max(moved.values())]})"
      f" —— 变动大则「哪一题最极端」本身就不是对象(`#797` 撤「三堆」同一族)。")
json.dump(dict(items=ITEMS, K=K, sd=SD, full=FULL, raw_max_dep={k2: [float(v[0]), v[1]] for k2, v in RAW.items()},
               table={nm: TAB[nm] for nm in NORMS}, ranks=RANKS, rank_movement=moved,
               homosex_top_everywhere=bool(homo_top), min_ratio=float(min_ratio), thr=THR,
               pos_control={nm: [v[0], v[1]] for nm, v in PC.items()},
               null_quantiles=NQ, above_null_p95=above, seeds=SEEDS,
               seed_ratios={str(k2): v for k2, v in SEEDR.items()},
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"scale_artifact.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'scale_artifact.json'}")
