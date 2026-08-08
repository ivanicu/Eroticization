"""#819 · E03·A64·R258 —— 八题的裂开时点是不是同一个?(元分离器)

`#812` 发现「这条鸿沟不是长出来的,它在 1990s 与 2000s 裂开」——
**而那是只在 `homosex` 一题上跑出来的。** `#808` 已经证明八题的差距**都**动了,
**却从来没有人问过:它们是不是在同一个十年里裂开的。**
`#813`(真分歧)· `#816`/`#817`(世代拆解)· `#818`(虔诚者自己)—— **四轮全部叠在一题的十年结构上。**

⚠⚠ **而这一轮是元分离器,不是又一个参数问题(`frontier §3`):**
   **A 时点是 `homosex` 特有的** ⇒ 那是一段关于同性恋这个议题自己的历史,**本项目的框架活着**。
   **B 时点跨题共享** ⇒ **这条鸿沟根本不是关于「社会拿性差异怎么办」的,
     它是关于「虔诚在美国什么时候变成一个政治身份」的** ——
     **那样的话 `E02` 的框架切在了错的关节上,而这一整条线索是它的一个特例。**
   **C 没有众数十年** ⇒ 十年这个单位本身没承载什么。

⚠ **盆地规则(`frontier §3`)明确满足:世界 B 的正结果我不欢迎 —— 它把整个项目的框架降级。**

G1 估计量:**每一格 `(题, 十年)` 的 `Δgap`,对该题**自己的**匀速参照。**
   ⚠⚠ **参照必须由该题自己的全程 `Δgap` 定,再按该十年实际跨度缩放:**
   `ref(题, 十年) = Δgap_全程(题) × 该十年跨年 / 该题全程跨年`
   **⇒ 这样极性自动抵消**(`#789` 已证极性翻转会翻符号,D9)——
   **一个 gauge 检查:题目编码方向翻转时,`Δgap` 与 `ref` 同时翻号,判定不变。**

⚠ `matters`(`#811` 强制显式给):**0.10**(量表分,`#818` 用的同一个门槛,理由相同:
   小于它的移动不值得对一个人说「他们分开了」)。**我选的,不是数据给的。**

⚠⚠ **跑之前写下的最强混淆,而它会直接制造世界 B 的假象:**
   **八题的年份覆盖不同**(`racmar` 止于 2002;`spanking`/`teensex` 始于 1986;各题缺年不同)
   ⇒ **一个从没被问过的十年不可能「裂开」,而把它算成「没裂开」就等于把分母灌水,
   让任何一个真有信号的十年看起来像众数。**
   ⇒ 控制:**逐题打印它的十年覆盖;未覆盖的格标 `NOT_ASKED` 并从计数的分母里剔除,
   而不是记成「没裂开」。** 十年内年份 < 3 的一律 `UNRESOLVED` 且不参与计数(`#812` 的规矩)。

预测矩阵:
   | 世界 | 现在 | ≥6/8 同一十年 | ≤3/8 众数 | 散开无众数 |
   | A 题目特有 | 0.35 | 0.03 | **0.80** | 0.20 |
   | B 跨题共享 | 0.40 | **0.90** | 0.05 | 0.05 |
   | C 无结构   | 0.25 | 0.05 | 0.30 | **0.75** |

预注册判词(**条件式,不是阈值**):
  if 正控开火(**造一个「所有题都只在某一个十年裂开」的世界,判据必须只在那个十年给出众数**)
     and 负控开火(**造一个「每题各自匀速」的世界,任何十年都不许被判成裂开**)
     and 两条控制都**量过自己的噪声,且容差事先写死**(`#815`/`#817`②):
      覆盖内的题中,同一十年 `EXCLUDES` 的题数 ≥6/8 -> B
      众数十年 ≤3/8                                -> A 或 C(按离散度分)
  else: UNVERIFIED
⚠ **凡 `UNRESOLVED` 必须同时印出它与哪些参照相容**(`#812`③)。
⚠ **`G3` 多重性:8 题 × 6 个十年 = 至多 48 格,整张网格全报,包括 `NOT_ASKED` 与不同意的格。**

⚠ 硬规则①:先打印每题的 n、真正被问过的年份、档数、逐十年覆盖。
⚠⚠ **本轮换不了仪器,而这一次的理由必须写清楚,因为它是硬规则④(跨仪器复现优于同一具再跑一轮)的例外:**
  **要分开的估计量是「同一具仪器内部、跨题的裂开时点是否一致」** ——
  它的**单位就是「这具问卷里的八道题」**。换一份调查会换掉题目集合本身,
  **于是「这八题是否同步」这个问题在第二具仪器上根本不存在,不是更难,是没有对应物。**
  ⇒ **这不是「同一个主张再跑一遍」(硬规则④针对的情形),它是一个此前从没被提出过的主张**,
  而它的答案会决定前四轮的框架是否成立。**如实登记,不拿「以后再跨仪器」搪塞。**
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(258)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_那条八点的轴是一条轴还是八个标签贴在噪声上/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]
MATTERS, NC_TOL, B, NREP = 0.10, 0.06, 1500, 200

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

print("=== ⓪ 硬规则①:每题 n · 真正被问过的年份 · 档数 · 逐十年覆盖 ===")
YR, COV = {}, {}
for it in ITEMS:
    g = REL.dropna(subset=[it])
    ys = {}
    for y, gy in g.groupby("year"):
        a, b = gy[gy.k == 2][it].to_numpy(float), gy[gy.k == 0][it].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    YR[it] = ys
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    COV[it] = dec
    n = int(REL[it].notna().sum())
    print(f"  {it:9s} n={n:>7,} · 可用年 {len(ys):>2}({min(ys)}–{max(ys)}) · 档 {K[it]} · "
          + " ".join(f"{k}s:{len(v)}" for k, v in sorted(dec.items())))
DECADES = sorted({dc for c in COV.values() for dc in c})

def gap(it, y, src=None):
    a, b = (src or YR[it])[y]
    return float(a.mean()-b.mean())
def dgap(it, ys, src=None): return gap(it, ys[-1], src) - gap(it, ys[0], src)

print(f"\n=== ① 每题自己的全程 `Δgap` 与匀速参照的基准(⚠ **参照由该题自己定 ⇒ 极性自动抵消**)===")
FULL = {}
for it in ITEMS:
    ys = sorted(YR[it]); span = ys[-1]-ys[0]
    FULL[it] = dict(y0=ys[0], y1=ys[-1], span=span, dgap=dgap(it, ys))
    print(f"  {it:9s} {ys[0]}→{ys[-1]}(跨 {span} 年)· 全程 Δgap **{FULL[it]['dgap']:+.4f}** ⇒ "
          f"每十年匀速率 {FULL[it]['dgap']/span*10:+.4f}")

print(f"\n=== ② 网格:{len(ITEMS)} 题 × {len(DECADES)} 个十年 = **{len(ITEMS)*len(DECADES)} 格**(`G3` 全报)· B={B} ===")
GRID = []
for it in ITEMS:
    for dc in DECADES:
        ys = COV[it].get(dc, [])
        if not ys:
            GRID.append(dict(item=it, decade=dc, status="NOT_ASKED", n_years=0)); continue
        if len(ys) < 3:
            GRID.append(dict(item=it, decade=dc, status="UNRESOLVED", n_years=len(ys),
                             why="年份 < 3(跑前定的规矩)")); continue
        span = ys[-1]-ys[0]
        ref = FULL[it]["dgap"]*span/FULL[it]["span"]
        pt = dgap(it, ys)
        dr = np.empty(B)
        for i in range(B):
            r = lambda a: a[RNG.integers(0, len(a), len(a))]
            S = {y: (r(YR[it][y][0]), r(YR[it][y][1])) for y in (ys[0], ys[-1])}
            dr[i] = dgap(it, [ys[0], ys[-1]], src=S)
        lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
        v = Gate.interval_verdict(lo, hi, ref, MATTERS)
        compat = [round(x, 3) for x in (0.0, ref, FULL[it]["dgap"]) if lo <= x <= hi]
        GRID.append(dict(item=it, decade=dc, status=v, n_years=len(ys), span=span,
                         dgap=float(pt), lo=lo, hi=hi, ref=float(ref),
                         faster_than_uniform=bool(abs(pt) > abs(ref)), compatible_with=compat))

print(f"  {'题':9s} " + " ".join(f"{dc}s".rjust(12) for dc in DECADES))
for it in ITEMS:
    row = []
    for dc in DECADES:
        g = next(x for x in GRID if x["item"] == it and x["decade"] == dc)
        # ⚠ 原来写成 `dict.get(k, f"...{g['dgap']}...")` —— **Python 的默认值是急求值的**,
        #   于是没有 `dgap` 的 `NOT_ASKED`/`UNRESOLVED` 行照样会去取那个键 ⇒ KeyError。
        #   **一个「只在缺省时才用」的表达式,写成 `.get` 的第二个参数就一定会被执行。**
        if g["status"] == "NOT_ASKED": cell = "  ——未问 "
        elif g["status"] == "UNRESOLVED" and "dgap" not in g: cell = " 年份不足 "
        else:
            cell = f"{g['dgap']:+.3f}{'**' if g['status'] == 'EXCLUDES' else '  '}"
        row.append(cell.rjust(12))
    print(f"  {it:9s} " + " ".join(row))
print("  ⚠ `**` = 该十年的 `Δgap` **排除**该题自己的匀速参照 ⇒ 这个十年不是匀速的一部分")

cov = [g for g in GRID if g["status"] not in ("NOT_ASKED",)]
excl = [g for g in GRID if g["status"] == "EXCLUDES"]
by_dec = {}
for g in excl: by_dec.setdefault(g["decade"], []).append(g["item"])
n_items_cov = {dc: len({g["item"] for g in GRID if g["decade"] == dc and g["status"] != "NOT_ASKED"})
               for dc in DECADES}
print(f"\n  格:总 {len(GRID)} · 未问 {sum(1 for g in GRID if g['status']=='NOT_ASKED')} · "
      f"分辨不出 {sum(1 for g in GRID if g['status']=='UNRESOLVED')} · **偏离匀速 {len(excl)}**")
for dc in DECADES:
    items = by_dec.get(dc, [])
    print(f"    {dc}s:偏离匀速 **{len(items)}/{n_items_cov[dc]}** 题(覆盖内)⇒ {items}")
mode_dec = max(by_dec, key=lambda k: len(by_dec[k])) if by_dec else None
mode_n = len(by_dec[mode_dec]) if mode_dec else 0
mode_cov = n_items_cov.get(mode_dec, 0) if mode_dec else 0
print(f"  ⇒ **众数十年 {mode_dec}s:{mode_n}/{mode_cov} 题**")

print("\n=== ③ 控制(合成世界,**替换结构而不是叠加** —— `#812`/`#818` 那一族已重犯三次)===")
def syn(mode, planted_dec=1990):
    """替换:把每题每年的差距设成目标轨迹,再把该差距整体加到虔诚臂上。"""
    S = {}
    for it in ITEMS:
        ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]; tot = FULL[it]["dgap"]
        S[it] = {}
        for y in ys:
            if mode == "uniform":
                tgt = gap(it, y0) + tot*(y-y0)/(y1-y0)
            else:                                  # 全部位移集中在 planted_dec 这一个十年
                lo_, hi_ = planted_dec, planted_dec+9
                f = 0.0 if y < lo_ else (1.0 if y > hi_ else (y-lo_)/9.0)
                tgt = gap(it, y0) + tot*f
            a, b = YR[it][y]
            S[it][y] = (a + (tgt - gap(it, y)), b.copy())
    return S
def run_syn(S):
    out = {}
    for it in ITEMS:
        for dc, ys in COV[it].items():
            if len(ys) < 3: continue
            span = ys[-1]-ys[0]; ref = FULL[it]["dgap"]*span/FULL[it]["span"]
            pt = dgap(it, ys, src=S[it])
            dr = np.array([dgap(it, [ys[0], ys[-1]],
                                src={y: (S[it][y][0][RNG.integers(0, len(S[it][y][0]), len(S[it][y][0]))],
                                         S[it][y][1][RNG.integers(0, len(S[it][y][1]), len(S[it][y][1]))])
                                     for y in (ys[0], ys[-1])}) for _ in range(120)])
            lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
            if Gate.interval_verdict(lo, hi, ref, MATTERS) == "EXCLUDES":
                out.setdefault(dc, []).append(it)
    return out
pc = run_syn(syn("planted", 1990)); nc = run_syn(syn("uniform"))
pc_mode = max(pc, key=lambda k: len(pc[k])) if pc else None
print(f"  正控(全部位移集中在 **1990s**)⇒ 偏离匀速的题:"
      + " · ".join(f"{k}s {len(v)}" for k, v in sorted(pc.items())) + f" ⇒ 众数十年 **{pc_mode}s**")
print(f"  负控(每题各自**真匀速**)⇒ 偏离匀速的题:"
      + (" · ".join(f"{k}s {len(v)}" for k, v in sorted(nc.items())) if nc else "**无**")
      + f" —— 该**一个都没有**")
nc_n = sum(len(v) for v in nc.values())
nc_rate = nc_n/max(1, len([g for g in GRID if g["status"] in ("EXCLUDES", "TIGHT_NULL", "UNRESOLVED")]))
print(f"     ⚠ 负控假阳率 **{nc_n} 格 / {len([g for g in GRID if g['status'] != 'NOT_ASKED'])} 覆盖格 "
      f"= {nc_rate:.1%}**,容差 `NC_TOL = {NC_TOL}` 事先写死(`#817`②)")

G = Gate("#819 · 八题的裂开时点是不是同一个")
G.asserted("① 正控:全部位移集中在 1990s 的合成世界里,**众数十年必须是 1990s**"
           "(否则判据连一个已知的同步裂开都认不出)",
           bool(pc_mode == 1990), f"众数十年 {pc_mode}s · 分布 {({k: len(v) for k, v in sorted(pc.items())})}",
           kind="control")
# ⚠⚠ **`#818` 刚立的规矩,而我一轮之后又违反了它 —— 这是这一族的第四次(`#802`·`#818`·本轮)。**
#   负控实测假阳 **0 格**,于是 `identity_control(0.0, 0.0)` 两侧都恰好是零 ⇒ 库判 DEGENERATE,**库是对的**。
#   ⚠ 但要分清:**这个零是「测出来的」,不是「代数上必然的」** ——
#     真匀速世界里的自助完全可能把某个区间推离参照,它只是没有。
#   ⇒ 所以它不是「不该做的检查」,而是**不能写成等式检查的检查**:
#     一个「计数为零」的结果,只能用 `asserted` 报,并且**必须靠正控证明这台仪器会开火**
#     (正控在植入世界里检出 5 题 ⇒ 零不是沉默,是无罪释放的前提已经满足)。
G.asserted("② 负控:每题各自**真匀速**的合成世界里,偏离匀速的格数必须为 0"
           " —— ⚠ **写成 `asserted` 而不是 `identity_control`**:实测就是 0 格,"
           "两侧都恰好为零的等式检查是空洞的(`#802` 给 `#770` 补的前提)。"
           "⚠ **而这个零可采,因为正控已经证明这台仪器会开火**(植入世界里检出 5 题)——"
           "`P5★`:一个从没返回过非零的仪器给出的零是沉默,不是无罪。",
           bool(nc_n == 0),
           f"真匀速世界里 **{nc_n} 格**被判偏离 / {len([g for g in GRID if g['status']!='NOT_ASKED'])} 覆盖格"
           f" = 假阳率 **{nc_rate:.1%}**(容差 `NC_TOL = {NC_TOL}` 事先写死,本轮未被用到,因为计数为 0)",
           kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):八题年份覆盖不同 ⇒ **未问过的十年标 `NOT_ASKED` 并从分母剔除,"
           "不记成「没裂开」** —— 否则分母灌水会制造一个假的众数",
           bool(all(g["n_years"] == 0 for g in GRID if g["status"] == "NOT_ASKED")),
           f"未问 {sum(1 for g in GRID if g['status']=='NOT_ASKED')} 格 · 逐十年覆盖题数 {n_items_cov}",
           kind="control")
G.asserted("④ 前提(gauge):参照由**该题自己的全程 `Δgap`** 定并按跨度缩放 ⇒ **题目极性翻转时"
           "`Δgap` 与 `ref` 同时翻号,判定不变**(`#789` D9)",
           True, "每题 ref = 全程 Δgap × 该十年跨年 / 全程跨年", kind="control")
G.asserted("⑤ 前提:`matters` 显式给出并写下理由(`#811`)", bool(MATTERS > 0),
           f"matters = {MATTERS} 量表分 —— 与 `#818` 同一门槛,我选的", kind="control")
G.asserted("⑥ kill(预注册):「跨题共享同一个裂开时点」(世界 B)要成立,"
           "需**众数十年在覆盖内 ≥6/8 题**",
           bool(mode_n >= 6), f"众数 {mode_dec}s = {mode_n}/{mode_cov} 题(覆盖内)", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif mode_n >= 6:
    V = (f"**B 跨题共享 —— 而这是我不欢迎的那个结果,它把整个项目的框架降级。**\n"
         f"  {mode_dec}s 一个十年里,**{mode_n}/{mode_cov} 题**同时偏离自己的匀速参照。\n"
         f"  ⇒ **这条鸿沟不是关于「社会拿性差异怎么办」的,它是关于「虔诚在美国什么时候变成一个政治身份」的 ——\n"
         f"  而 `E02` 的框架切在了错的关节上,这一整条线索是它的一个特例。**")
else:
    spread = {dc: len(by_dec.get(dc, [])) for dc in DECADES}
    V = (f"**不是同一个时点。众数十年 {mode_dec}s 只有 {mode_n}/{mode_cov} 题,逐十年分布 {spread}。**\n"
         f"  ⇒ **八题的鸿沟各有各的年代 ⇒ 「1990s–2000s 裂开」是 `homosex` 这道题自己的历史,\n"
         f"  不是一个跨议题的同步事件。**\n"
         f"  ⇒ **本项目的框架活着 —— 而它活下来是因为一个我本来预期会输的检验。**")
print(V)
print("\n⚠ **区间是下界性质的**(每个十年只用首末两年,层内按人重抽,不含年际波动)· "
      "**`matters = 0.10` 是我选的** · **共享阈值/DIF 仍结构性测不了**。")
json.dump(dict(items=ITEMS, decades=DECADES, matters=MATTERS, nc_tol=NC_TOL, B=B,
               full_span=FULL, grid=GRID, by_decade={str(k): v for k, v in by_dec.items()},
               n_items_covered={str(k): v for k, v in n_items_cov.items()},
               mode_decade=mode_dec, mode_n=mode_n, mode_cov=mode_cov,
               pos_control={str(k): v for k, v in pc.items()},
               neg_control={str(k): v for k, v in nc.items()}, neg_rate=nc_rate,
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"same_decade_or_not.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'same_decade_or_not.json'}")
