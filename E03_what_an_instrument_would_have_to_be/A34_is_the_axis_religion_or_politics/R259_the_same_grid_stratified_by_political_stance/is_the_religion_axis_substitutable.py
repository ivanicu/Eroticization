"""#820 · E03·A65·R259 —— 「宗教」这个轴是不是可替换的?换成政治立场,再跑同一张网格

`#819` 的结论:**八题里只有 `homosex` 在猛地裂,其余七题各按自己的匀速走。**
⚠⚠ **而那句话有一个我没检验过的前提:它是按「虔诚度」分的层。**
   **若把分层换成「政治立场」,同一张网格还长这样吗?**
   **A 长得一样** ⇒ **「宗教」这个轴是可替换的** —— 那么 `#819` 说的根本不是关于宗教的事,
     是关于**同性恋这个议题本身**的事,而这条线索一直挂在「虔诚」上是个误挂。
   **B 长得不一样(多题同步裂)** ⇒ **`#819` 只对「宗教」这个切法成立** ——
     **宗教与政治在美国是按不同的节奏分开的,而那个差别本身就是发现。**
   **C 政治那张网格分辨不出** ⇒ 登记功效边界,不硬判。

⚠ **两个结果都会改写对象,所以这是又一个元分离器,不是参数问题。**

G1 估计量:**同一张 `(题 × 十年)` 网格,在两个分层变量下各跑一遍**,比较
   ① 各自的**众数十年题数** ② **逐格判定的一致率**。
   ⚠ **比较本身就是估计量** —— 不是跑两遍然后用眼睛看(`#819` 与本轮若分两轮,就只能用眼睛比)。

⚠ 硬规则①(已先跑过一次,结论抄在这里并在下面复打):
   `polviews` n=**65,878** · 年 **33**(1974–2024)· 取值 1–7 · **它的年份范围内没有缺年**。
⚠ 硬规则②(**这条主张路由过哪具仪器**):**两个分层都来自同一份问卷 GSS** ——
   **这不是跨仪器复现,是同一具仪器内换一根轴。** 若两张网格一致,
   **它排除的是「轴选错了」,不是「仪器选错了」。如实写死,不许升格。**

⚠⚠ **跑之前写下的最强混淆,而它足以自己制造世界 A:**
   **虔诚与政治保守在美国高度相关** ⇒ 两个「顶层」可能就是同一批人 ⇒
   **两张网格一致什么也不证明,它只证明我分了两次同一群人。**
   ⇒ 控制:**逐年测两个顶层的重叠率(Jaccard)与两个连续变量的相关**,并**印出来**;
   **重叠率高就明说这一轮的分辨力被它吃掉了多少 —— 不假装控制了它。**

`matters` = **0.10**(与 `#818`/`#819` 同一门槛,我选的)。

预注册判词(条件式):
  if 正控开火(**两个分层下,「全部位移集中在 1990s」的合成世界都必须只在 1990s 检出**)
     and 负控开火(**两个分层下,真匀速世界的假阳都必须为 0**):
      `polviews` 网格的众数十年 ≥6/8 题 -> B(政治是同步的,宗教不是)
      两张网格逐格一致率 ≥ 0.85 且众数都 ≤3/8 -> A(轴可替换)
      否则 -> 报整张双网格,不选边
  else: UNVERIFIED
⚠ **`G3`:两张网格共 96 格,全报,包括 `NOT_ASKED` 与不一致的格。**

⚠ 本轮**换不了仪器**:估计量是「同一具问卷内部换一根分层轴」,
  第二份调查会同时换掉题目集合与分层变量,**于是「同一张网格换一根轴」在那里没有对应物。**
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(259)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
P791 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A45_虔诚者踩的是一个总闸还是一个个开关/"
                      "R230_is_the_eight_point_axis_an_axis_or_eight_labels_on_noise/results/is_the_ordering_an_object.json"))
ITEMS = P791["items"]
MATTERS, B = 0.10, 1200

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
R = BASE.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1)
BASE = BASE.join(R["REL"])
print("=== ⓪ 硬规则①:两个分层变量的 n · 真正被问过的年份 ===")
for c in ("REL", "polviews"):
    v = BASE[c]; ys = sorted(BASE.year[v.notna()].unique())
    print(f"  {c:9s} n={int(v.notna().sum()):>7,} · 年 {len(ys):>2}({int(min(ys))}–{int(max(ys))})")
print("  ⚠ 硬规则②:**两个分层都来自同一份问卷 GSS ⇒ 这不是跨仪器复现,是同一具仪器内换一根轴。**")

for name, col in (("k_rel", "REL"), ("k_pol", "polviews")):
    BASE[name] = BASE.groupby("year")[col].transform(
        lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)

print("\n=== ⓪b 跑前写下的最强混淆:两个顶层是不是同一批人 ===")
ov = BASE.dropna(subset=["k_rel", "k_pol"])
jac, cor = [], []
for y, gy in ov.groupby("year"):
    a, b = gy.k_rel == 2, gy.k_pol == 2
    inter, union = int((a & b).sum()), int((a | b).sum())
    if union: jac.append(inter/union)
    if gy.REL.notna().sum() > 30: cor.append(float(gy.REL.corr(gy.polviews)))
print(f"  两个顶层的 Jaccard 重叠(逐年):中位 **{np.median(jac):.3f}** · 范围 [{min(jac):.3f}, {max(jac):.3f}]")
print(f"  `REL` 与 `polviews` 的逐年相关:中位 **{np.median(cor):+.3f}** · 范围 [{min(cor):+.3f}, {max(cor):+.3f}]")
print(f"  ⇒ ⚠ **重叠约 {np.median(jac):.0%} ⇒ 两张网格并非独立;一致时它排除的东西比看起来少,如实说。**")

def build(kcol):
    YR, COV = {}, {}
    for it in ITEMS:
        g = BASE.dropna(subset=[it, kcol])
        ys = {}
        for y, gy in g.groupby("year"):
            a = gy[gy[kcol] == 2][it].to_numpy(float); b = gy[gy[kcol] == 0][it].to_numpy(float)
            if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
        YR[it] = ys
        dec = {}
        for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
        COV[it] = dec
    return YR, COV

def grid(YR, COV, tag, src_override=None):
    S = src_override or YR
    gp_ = lambda it, y: float(S[it][y][0].mean()-S[it][y][1].mean())
    FULL = {}
    for it in ITEMS:
        ys = sorted(YR[it]); FULL[it] = dict(y0=ys[0], y1=ys[-1], span=ys[-1]-ys[0],
                                             dgap=gp_(it, ys[-1])-gp_(it, ys[0]))
    G_, by = [], {}
    for it in ITEMS:
        for dc in DECADES:
            ys = COV[it].get(dc, [])
            if not ys: G_.append(dict(item=it, decade=dc, status="NOT_ASKED")); continue
            if len(ys) < 3: G_.append(dict(item=it, decade=dc, status="UNRESOLVED")); continue
            span = ys[-1]-ys[0]; ref = FULL[it]["dgap"]*span/FULL[it]["span"]
            pt = gp_(it, ys[-1])-gp_(it, ys[0])
            dr = np.empty(B)
            for i in range(B):
                r = lambda a: a[RNG.integers(0, len(a), len(a))]
                s0, s1 = (r(S[it][ys[0]][0]), r(S[it][ys[0]][1])), (r(S[it][ys[-1]][0]), r(S[it][ys[-1]][1]))
                dr[i] = (s1[0].mean()-s1[1].mean()) - (s0[0].mean()-s0[1].mean())
            lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
            v = Gate.interval_verdict(lo, hi, ref, MATTERS)
            G_.append(dict(item=it, decade=dc, status=v, dgap=float(pt), lo=lo, hi=hi, ref=float(ref)))
            if v == "EXCLUDES": by.setdefault(dc, []).append(it)
    return G_, by, FULL

YRr, COVr = build("k_rel"); YRp, COVp = build("k_pol")
DECADES = sorted({(y//10)*10 for C in (COVr, COVp) for c in C.values() for y in sum(c.values(), [])})
Gr, BYr, FULLr = grid(YRr, COVr, "REL")
Gp, BYp, FULLp = grid(YRp, COVp, "polviews")

def show(G_, BY, tag):
    covn = {dc: len({g["item"] for g in G_ if g["decade"] == dc and g["status"] != "NOT_ASKED"})
            for dc in DECADES}
    print(f"\n  —— 分层 = **{tag}**")
    print("  " + "题".ljust(10) + "".join(f"{dc}s".rjust(12) for dc in DECADES))
    for it in ITEMS:
        row = ""
        for dc in DECADES:
            g = next(x for x in G_ if x["item"] == it and x["decade"] == dc)
            if g["status"] == "NOT_ASKED": c = "——未问"
            elif "dgap" not in g: c = "年份不足"
            else: c = f"{g['dgap']:+.3f}{'**' if g['status'] == 'EXCLUDES' else ''}"
            row += c.rjust(12)
        print("  " + it.ljust(10) + row)
    md = max(BY, key=lambda k: len(BY[k])) if BY else None
    print(f"    逐十年偏离匀速:{ {k: v for k, v in sorted(BY.items())} }")
    print(f"    ⇒ **众数十年 {md}s:{len(BY.get(md, []))}/{covn.get(md, 0)} 题**")
    return md, len(BY.get(md, [])), covn.get(md, 0), covn

print(f"\n=== ① 两张网格({len(ITEMS)}题 × {len(DECADES)}十年 × 2 分层 = **{len(ITEMS)*len(DECADES)*2} 格**,`G3` 全报)===")
mdr, mnr, mcr, covnr = show(Gr, BYr, "虔诚度 `REL`(`#819` 的那张)")
mdp, mnp, mcp, covnp = show(Gp, BYp, "政治立场 `polviews`")

pairs = [(a, b) for a in Gr for b in Gp if a["item"] == b["item"] and a["decade"] == b["decade"]]
both = [(a, b) for a, b in pairs if a["status"] not in ("NOT_ASKED",) and b["status"] not in ("NOT_ASKED",)]
agree = sum(1 for a, b in both if a["status"] == b["status"])
rate = agree/len(both) if both else float("nan")
print(f"\n  两张网格逐格一致:**{agree}/{len(both)} = {rate:.1%}**(两侧都有覆盖的格)")
diff = [(a["item"], a["decade"], a["status"], b["status"]) for a, b in both if a["status"] != b["status"]]
print(f"  ⚠ **不一致的格全列(`G3`)**:{diff if diff else '无'}")

print("\n=== ② 控制(两个分层各跑一遍;合成世界**替换轨迹,不叠加**)===")
def syn(YR, FULL, mode, planted=1990):
    S = {}
    for it in ITEMS:
        ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]; tot = FULL[it]["dgap"]
        g0 = float(YR[it][y0][0].mean()-YR[it][y0][1].mean())
        S[it] = {}
        for y in ys:
            f = ((y-y0)/(y1-y0) if mode == "uniform"
                 else (0.0 if y < planted else (1.0 if y > planted+9 else (y-planted)/9.0)))
            cur = float(YR[it][y][0].mean()-YR[it][y][1].mean())
            S[it][y] = (YR[it][y][0] + (g0 + tot*f - cur), YR[it][y][1].copy())
    return S
CTL = {}
for tag, YR, COV, FULL in (("REL", YRr, COVr, FULLr), ("polviews", YRp, COVp, FULLp)):
    gp1, by1, _ = grid(YR, COV, tag, src_override=syn(YR, FULL, "planted"))
    gp0, by0, _ = grid(YR, COV, tag, src_override=syn(YR, FULL, "uniform"))
    md1 = max(by1, key=lambda k: len(by1[k])) if by1 else None
    CTL[tag] = dict(pos_mode=md1, pos={str(k): v for k, v in by1.items()},
                    neg_n=sum(len(v) for v in by0.values()), neg={str(k): v for k, v in by0.items()})
    print(f"  分层 {tag:9s} 正控(全部集中在 1990s)⇒ 众数 **{md1}s**,分布 { {k: len(v) for k, v in sorted(by1.items())} }"
          f" · 负控(真匀速)⇒ 假阳 **{CTL[tag]['neg_n']} 格**")

G = Gate("#820 · 「宗教」这个轴是不是可替换的")
G.asserted("① 正控:**两个分层下**,「全部位移集中在 1990s」的合成世界都必须以 1990s 为众数",
           bool(CTL["REL"]["pos_mode"] == 1990 and CTL["polviews"]["pos_mode"] == 1990),
           f"REL ⇒ {CTL['REL']['pos_mode']}s · polviews ⇒ {CTL['polviews']['pos_mode']}s", kind="control")
G.asserted("② 负控:**两个分层下**,真匀速世界的假阳都必须为 0 —— ⚠ **写成 `asserted` 不是 `identity_control`**"
           "(`#819`:两侧恰好为零的等式检查是空洞的;而这个零可采,因为①已证明仪器会开火)",
           bool(CTL["REL"]["neg_n"] == 0 and CTL["polviews"]["neg_n"] == 0),
           f"REL 假阳 {CTL['REL']['neg_n']} 格 · polviews 假阳 {CTL['polviews']['neg_n']} 格", kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):两个顶层的重叠率与两变量相关**已测并印出**,"
           "**且不假装控制了它** —— 重叠高则一致所排除的东西比看起来少",
           bool(len(jac) > 0), f"Jaccard 中位 {np.median(jac):.3f} · corr 中位 {np.median(cor):+.3f}",
           kind="control")
G.asserted("④ 前提(gauge):参照由该题在**该分层下**自己的全程 `Δgap` 定 ⇒ 极性与分层方向都自动抵消",
           True, "每格 ref = 该分层下该题全程 Δgap × 十年跨年 / 全程跨年", kind="control")
G.asserted("⑤ 前提:`matters` 显式给出(`#811`)", bool(MATTERS > 0),
           f"matters = {MATTERS} —— 与 `#818`/`#819` 同一门槛,我选的", kind="control")
G.asserted("⑥ kill(预注册):「政治立场那一根轴是同步的」(世界 B)要成立,需 `polviews` 众数十年 ≥6/8 题",
           bool(mnp >= 6), f"polviews 众数 {mdp}s = {mnp}/{mcp} 题", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif mnp >= 6:
    V = (f"**B 政治那根轴是同步的,而宗教那根不是。** `polviews` 众数 {mdp}s = **{mnp}/{mcp}** 题,"
         f"而 `REL` 众数 {mdr}s = {mnr}/{mcr} 题。\n"
         f"  ⇒ **`#819` 的「各有各的年代」只对「宗教」这个切法成立 ——\n"
         f"  美国的宗教分歧与政治分歧是按不同节奏分开的,而那个差别本身就是发现。**")
elif rate >= 0.85 and mnr <= 3 and mnp <= 3:
    V = (f"**A 两根轴给出同一张网格(逐格一致 {rate:.1%},两边众数都 ≤3/8)。**\n"
         f"  ⇒ **「宗教」这个轴是可替换的 ⇒ `#819` 说的不是关于宗教的事,是关于同性恋这个议题本身的事。**\n"
         f"  ⚠⚠ **而这个结论要打一个大折扣,折扣是我跑前就写下的**:两个顶层的 Jaccard 重叠中位 "
         f"**{np.median(jac):.0%}** ⇒ **两张网格并非独立,一致所排除的东西比它看起来少。**")
else:
    V = (f"**两张网格既不同步也不完全一致 —— 报双网格,不选边。**\n"
         f"  `REL` 众数 {mdr}s {mnr}/{mcr} · `polviews` 众数 {mdp}s {mnp}/{mcp} · 逐格一致 {rate:.1%}\n"
         f"  ⚠ 不一致的格已全列(`G3`),而它们是下一轮该看的地方。")
print(V)
json.dump(dict(items=ITEMS, decades=DECADES, matters=MATTERS, B=B,
               grid_rel=Gr, grid_pol=Gp, by_rel={str(k): v for k, v in BYr.items()},
               by_pol={str(k): v for k, v in BYp.items()},
               mode_rel=[mdr, mnr, mcr], mode_pol=[mdp, mnp, mcp],
               agree=agree, n_both=len(both), agree_rate=rate, disagreements=diff,
               jaccard_median=float(np.median(jac)), corr_median=float(np.median(cor)),
               controls=CTL, admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"religion_axis_substitutable.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'religion_axis_substitutable.json'}")
