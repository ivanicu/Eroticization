"""#825 · E03·A70·R264 —— 那条幸存下来的事实,自己稳不稳?

`#824` 撤掉了「2.5 倍」,并写下**唯一还站着、值得追的一条**:
**`homosex` 是八题里唯一在**两个**十年偏离自己匀速率的题 —— 一个关于「次数」的事实。**

⚠⚠ **而 `#820` 早就量到过一件直接打在它身上的事,我在写 `#824`② 时没有回头看:**
   **只是换了种子与 `B`,`teensex` 2000s 那一格就从「未标记」翻成 `EXCLUDES`。**
   **而那一格恰好决定 `teensex` 是「一个十年」还是「两个十年」** ——
   **也就是说,`#824`② 那条「唯一」很可能是一次抽样的产物。**
   `#820` 当轮已登记:**「在补上多种子之前,单个格子的标记不许当硬事实」** ——
   **而我随即就用单个格子拼出了一个计数,并把它写进了 NEXT。这笔债现在到期。**

G1 估计量(**两个,而第二个才是结构性的**):
   ① **逐格标记稳定率** = 该格在 `S` 个种子里被判 `EXCLUDES` 的比例。
   ② **边界距离** = `|区间最近端 − 参照| ÷ 区间宽度` ——
      **一个区间端点离参照 0.001 的格子,不管 `B` 多大都会永远翻。**
   ⚠⚠ **这两个量回答的不是同一个问题,而分清它们正是 `#821` 的教训:**
      **①(翻转率)是蒙特卡洛误差,加大 `B` 就能压下去;
      ②(边界距离)是抽样不确定性,`B` 再大也压不动。**
      **只报 ① 会让人以为「多跑几次就稳了」,而那是可复现性不是精度。**

三个世界:
   A **计数稳**:`homosex` 唯一有两个十年,在 ≥4/5 种子下成立 ⇒ `#824`② 那条可以追。
   B **计数不稳**:`teensex`(或别的题)在部分种子下也有两个 ⇒
     **`#824`② 那条「唯一」撤回,而 `#819` 剩下的只有「众数十年远小于 6/8」这个粗结论。**
   C **边界距离说它永远不会稳** ⇒ **即使加大 `B`,这个计数也不是一个可测的量** ——
     那是关于设计的最终判词,比 B 更彻底。

预测矩阵:
   | 世界 | 现在 | ≥4/5 种子成立 | 有题也有两个 | 边界距离 < 0.1 的格 ≥1 |
   | A 稳     | 0.35 | **0.85** | 0.05 | 0.15 |
   | B 不稳   | 0.40 | 0.05 | **0.90** | 0.55 |
   | C 永不稳 | 0.25 | 0.10 | 0.50 | **0.85** |

预注册判词(条件式):
  if 正控开火(**造一个「某题只在一个十年裂开」的世界,该格必须在 5/5 种子下都被标记**)
     and 负控开火(**造一个全匀速世界,标记率必须 ≈ 0**,
        ⚠ **「这个零该不该是零?」——该。** 匀速世界里就是不该有任何格被标记,
        参照真的是 0;而**实测若恰好为 0,写成 `asserted` 不写 `identity_control`**,`#819` 的教训):
      `homosex` 唯一有 ≥2 个十年,在 ≥4/5 种子下成立 -> A
      否则 -> B/C,**报整张稳定率网格,不选边**
  else: UNVERIFIED
⚠ **`G3`:8 题 × 6 十年 × 5 种子,整张稳定率网格全报,包括从不翻的格。**

⚠ 跑之前写下的最强混淆:**加大 `B` 会让翻转率下降,而那容易被读成「问题解决了」。**
  ⇒ 控制:**在两个 `B`(1500 与 6000)下各跑一遍**,并**同时报边界距离** ——
  **若翻转率随 `B` 降但边界距离不变,那就证明稳定性是买来的,不是测出来的。**

⚠ 本轮**换不了仪器**:检验的是同一份数据上同一张网格的抽样稳定性。
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
SEEDS, BS, MATTERS = [301, 302, 303, 304, 305], [1500, 6000], 0.10

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
YR, COV, FULL = {}, {}, {}
for it in ITEMS:
    g = REL.dropna(subset=[it]); ys = {}
    for y, gy in g.groupby("year"):
        a = gy[gy.k == 2][it].to_numpy(float); b = gy[gy.k == 0][it].to_numpy(float)
        if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
    YR[it] = ys
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    COV[it] = {k2: v for k2, v in dec.items() if len(v) >= 3}
    s = sorted(ys)
    g0 = float(ys[s[0]][0].mean()-ys[s[0]][1].mean()); g1 = float(ys[s[-1]][0].mean()-ys[s[-1]][1].mean())
    FULL[it] = dict(span=s[-1]-s[0], dgap=g1-g0)
DECADES = sorted({dc for c in COV.values() for dc in c})

def cell(it, ys, rng, B, src=None):
    S = src[it] if src else YR[it]
    span = ys[-1]-ys[0]; ref = FULL[it]["dgap"]*span/FULL[it]["span"]
    dr = np.empty(B)
    for i in range(B):
        r = lambda a: a[rng.integers(0, len(a), len(a))]
        a0, b0 = r(S[ys[0]][0]), r(S[ys[0]][1]); a1, b1 = r(S[ys[-1]][0]), r(S[ys[-1]][1])
        dr[i] = (a1.mean()-b1.mean()) - (a0.mean()-b0.mean())
    lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
    v = Gate.interval_verdict(lo, hi, ref, MATTERS)
    # ⚠ 边界距离:区间最近端离参照多远,以区间宽度为单位。B 再大也压不动它。
    near = min(abs(lo-ref), abs(hi-ref)); width = hi-lo
    inside = lo <= ref <= hi
    return v, float(near/width if width > 1e-12 else np.nan), bool(inside), lo, hi, float(ref)

print(f"=== ① 逐格标记稳定率 · {len(SEEDS)} 个种子 × {len(BS)} 个 B(`G3` 全网格)===")
STAB = {}
for Bv in BS:
    for it in ITEMS:
        for dc, ys in sorted(COV[it].items()):
            marks, dists = [], []
            for sd in SEEDS:
                v, dist, inside, lo, hi, ref = cell(it, ys, np.random.default_rng(sd), Bv)
                marks.append(v == "EXCLUDES"); dists.append(dist)
            STAB[(Bv, it, dc)] = dict(rate=float(np.mean(marks)), dist=float(np.mean(dists)),
                                      n_seed=len(SEEDS))
for Bv in BS:
    print(f"\n  —— B = **{Bv}**   (格式:标记率 / 边界距离)")
    print("  " + "题".ljust(10) + "".join(f"{dc}s".rjust(15) for dc in DECADES))
    for it in ITEMS:
        row = ""
        for dc in DECADES:
            k2 = (Bv, it, dc)
            row += ("——".rjust(15) if k2 not in STAB
                    else f"{STAB[k2]['rate']:.1f} / {STAB[k2]['dist']:.3f}".rjust(15))
        print("  " + it.ljust(10) + row)
UNSTABLE = {k2: v for k2, v in STAB.items() if 0 < v["rate"] < 1}
print(f"\n  ⚠ **在种子间翻转的格(0 < 标记率 < 1)**:{len(UNSTABLE)} 个 —— "
      + (" · ".join(f"B{k2[0]}/{k2[1]}/{k2[2]}s(率{v['rate']:.1f},边界距{v['dist']:.3f})"
                    for k2, v in UNSTABLE.items()) if UNSTABLE else "**无**"))
CLOSE = {k2: v for k2, v in STAB.items() if v["dist"] < 0.10}
print(f"  ⚠ **边界距离 < 0.10 的格**(`B` 再大也压不动):{len(CLOSE)} 个 —— "
      + (" · ".join(f"{k2[1]}/{k2[2]}s({v['dist']:.3f})" for k2, v in CLOSE.items()) if CLOSE else "**无**"))

print(f"\n=== ② 逐种子的「departing 十年计数」 —— `#824`② 那条「唯一」的直接检验 ===")
COUNTS = {}
for Bv in BS:
    for sd in SEEDS:
        rng = np.random.default_rng(sd)
        c = {}
        for it in ITEMS:
            n = 0
            for dc, ys in sorted(COV[it].items()):
                v, *_ = cell(it, ys, np.random.default_rng(sd*1000+hash(it) % 997), Bv)
                n += (v == "EXCLUDES")
            c[it] = n
        COUNTS[(Bv, sd)] = c
for Bv in BS:
    print(f"  —— B = {Bv}")
    for sd in SEEDS:
        c = COUNTS[(Bv, sd)]
        two = [it for it in ITEMS if c[it] >= 2]
        print(f"    种子 {sd}:{ {k2: v for k2, v in c.items() if v > 0} } ⇒ ≥2 个十年的题:**{two}**")
uniq = [all(([it for it in ITEMS if COUNTS[(Bv, sd)][it] >= 2] == ["homosex"]) for sd in SEEDS) for Bv in BS]
n_uniq = sum(sum(1 for sd in SEEDS if [it for it in ITEMS if COUNTS[(Bv, sd)][it] >= 2] == ["homosex"])
             for Bv in BS)
tot = len(BS)*len(SEEDS)
print(f"  ⇒ **`homosex` 是唯一有 ≥2 个十年的题:{n_uniq}/{tot} 个(B × 种子)组合**")

print("\n=== ③ 控制 ===")
def syn(mode, rng, it_p="homosex", dec_p=1990):
    S = {}
    for it in ITEMS:
        ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]
        g0 = float(YR[it][y0][0].mean()-YR[it][y0][1].mean()); tot_ = FULL[it]["dgap"]
        S[it] = {}
        for y in ys:
            f = ((y-y0)/(y1-y0) if mode == "uniform"
                 else (0.0 if y < dec_p else (1.0 if y > dec_p+9 else (y-dec_p)/9.0)))
            cur = float(YR[it][y][0].mean()-YR[it][y][1].mean())
            a, b = YR[it][y]
            S[it][y] = (a[rng.integers(0, len(a), len(a))] + (g0 + tot_*f - cur),
                        b[rng.integers(0, len(b), len(b))])
    return S
Sp = syn("planted", np.random.default_rng(777)); Su = syn("uniform", np.random.default_rng(778))
pc_marks = [cell("homosex", COV["homosex"][1990], np.random.default_rng(sd), 1500, src=Sp)[0] == "EXCLUDES"
            for sd in SEEDS]
nc_n = 0
for it in ITEMS:
    for dc, ys in sorted(COV[it].items()):
        if cell(it, ys, np.random.default_rng(779), 1500, src=Su)[0] == "EXCLUDES": nc_n += 1
n_cells = sum(len(COV[it]) for it in ITEMS)
print(f"  正控(位移全集中在 `homosex` 的 1990s)⇒ 该格在 **{sum(pc_marks)}/{len(SEEDS)}** 个种子下被标记 —— 该 **5/5**")
print(f"  负控(八题全匀速)⇒ 被标记的格 **{nc_n}/{n_cells}** —— 该 **0**"
      f"(⚠ **「这个零该不该是零?」——该**:匀速世界里就不该有任何格被标记)")

G = Gate("#825 · 那条幸存下来的事实,自己稳不稳")
G.asserted("① 正控:位移全集中在 `homosex` 1990s 的合成世界里,该格必须在 **5/5** 种子下都被标记"
           "(否则这套标记机器连一个确定的裂开都认不稳)",
           bool(sum(pc_marks) == len(SEEDS)), f"{sum(pc_marks)}/{len(SEEDS)} 个种子", kind="control")
G.asserted("② 负控:八题全匀速的世界里,被标记的格必须为 0 —— ⚠ **「这个零该不该是零?」**"
           "**该**,匀速世界里就不该有任何格被标记,参照真的是 0;"
           "⚠ 而**实测若恰好为 0,写成 `asserted` 不写 `identity_control`**(`#819`/`#802` 的教训:"
           "两侧恰好为零的等式检查是空洞的),**而这个零可采,因为①已证明仪器会开火**",
           bool(nc_n == 0), f"被标记 {nc_n}/{n_cells} 格", kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**加大 `B` 会让翻转率下降而边界距离不变** ⇒ "
           "**两个 `B` 都跑并同时报边界距离** —— 否则「稳定」是买来的,不是测出来的",
           bool(len(BS) == 2 and all("dist" in v for v in STAB.values())),
           f"B = {BS} · 翻转格 {len(UNSTABLE)} · 边界距离 <0.10 的格 {len(CLOSE)}", kind="control")
G.asserted("④ 前提:两个估计量分开报 —— **①翻转率是蒙特卡洛误差(`B` 可压),②边界距离是抽样不确定性(`B` 压不动)**"
           "(`#821` 的教训:可复现性 ≠ 精度)", True,
           f"翻转率与边界距离逐格并列,共 {len(STAB)} 格", kind="control")
G.asserted("⑤ kill(预注册):「`homosex` 是唯一有 ≥2 个十年的题」(`#824`②)要成立,"
           "需在 **≥4/5 种子**下成立(两个 `B` 都要)",
           bool(n_uniq >= 4*len(BS)), f"成立 {n_uniq}/{tot} 个(B × 种子)组合", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n_uniq >= 4*len(BS):
    V = (f"**A 那条计数是稳的。** `homosex` 是唯一有 ≥2 个十年的题,在 **{n_uniq}/{tot}** 个"
         f"(B × 种子)组合下成立;翻转格 **{len(UNSTABLE)}**,边界距离 <0.10 的格 **{len(CLOSE)}**。\n"
         f"  ⇒ **`#824`② 那条可以追,而它现在带着自己的稳定性证据。**")
else:
    V = (f"**B 那条「唯一」不稳。** 只在 **{n_uniq}/{tot}** 个(B × 种子)组合下成立。\n"
         f"  ⇒ **`#824`② 撤回 —— 而 `#819` 剩下的只有「众数十年远小于 6/8」这个粗结论,\n"
         f"  它不依赖任何单个格子,所以它是这一串里唯一没被自己的稳定性打掉的东西。**\n"
         f"  ⚠ **而 `#820` 早就登记过「在补上多种子之前,单个格子的标记不许当硬事实」——\n"
         f"  我随即用单个格子拼出一个计数并写进 NEXT。这一轮是那笔债到期,而它到期时是负的。**")
print(V)
json.dump(dict(items=ITEMS, decades=DECADES, seeds=SEEDS, Bs=BS, matters=MATTERS,
               stability={f"{k2[0]}|{k2[1]}|{k2[2]}": v for k2, v in STAB.items()},
               unstable=[f"{k2[0]}|{k2[1]}|{k2[2]}" for k2 in UNSTABLE],
               close_to_boundary=[f"{k2[0]}|{k2[1]}|{k2[2]}" for k2 in CLOSE],
               counts={f"{k2[0]}|{k2[1]}": v for k2, v in COUNTS.items()},
               n_unique=n_uniq, n_total=tot,
               pos_control=int(sum(pc_marks)), neg_control=int(nc_n), n_cells=n_cells,
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"count_stability.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'count_stability.json'}")
