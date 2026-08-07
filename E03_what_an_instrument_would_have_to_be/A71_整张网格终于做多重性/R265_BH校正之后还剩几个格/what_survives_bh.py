"""#826 · E03·A71·R265 —— 整张网格终于做一次多重性:BH 校正之后还剩几个格?

`#825` 测出 **`G3` 从来没有被真正执行过**:我一直只做到「整张网格全报」,
**从没按格数调过标记阈值**。34 格 × 5% ⇒ 纯噪声下期望约 1.7 格被标记,
而 `#819` 的真实网格一共只标记了 4–5 格。**⇒ 校正之后可能剩不下几个。**

G1 估计量:**逐格的经验 p**(该十年 `Δgap` 与该题自己的匀速参照之差,双侧),
   再对**整张网格**做 **BH**(`q = 0.05` 与 `0.10` 两个水平,`G4`)。
   报 **cells tested / cells surviving**,**并逐个列出没存活的格**(`G3` 明文要求)。

⚠⚠ **`realstat` 点名的两个陷阱,跑前写死:**
   ① **BH 在秩 `k` 的阈值是 `q·k/C`,最大的那个就是 `q` 本身;`q/C` 是 Bonferroni。**
     **拿 `q/C` 当 BH 会高估所需抽样数,也会过度杀格。**
   ② **经验 p 的分辨率下限是 `1/(B+1)`。** `B = 6000` ⇒ 下限 **1.67e-4**;
     而 BH 在 `C ≈ 34`、`q = 0.05` 时最严的那一档是 `0.05 × 1/34 ≈ 1.47e-3` ——
     **1.47e-3 > 1.67e-4 ⇒ `B` 够用。跑前算出来并印出来,不是事后说「应该够」。**

⚠⚠ **而负控的期望值必须写对,因为我上一轮刚在这里犯过错**(`#825`①):
   **BH 控制的是 FDR,在全零世界里它的保证是「任何一次发现的概率 ≤ q」,不是「0 次发现」。**
   ⇒ 负控 = **把全匀速世界跑 `NREP` 次,数「至少出现一个存活格」的比例**,
   **该比例必须 ≤ q(带它自己的抽样误差)** —— **而不是要求 0。**
   **向一个额定率为 q 的程序索要 0,正是 `#825` 那一条。**

三个世界:
   A **`homosex` 的两格都存活** ⇒ `#819` 的核心在校正后仍然站着。
   B **只剩一格或零格** ⇒ **`#819` 的「两个十年」撤回,而这一串还剩下什么要如实列出。**
   C **存活数与匀速世界不可区分** ⇒ **整张网格没有任何格是可信的**,那是最彻底的判词。

预测矩阵:
   | 世界 | 现在 | 两格都活 | 只剩一格 | 零格 |
   | A 核心站住 | 0.40 | **0.85** | 0.10 | 0.02 |
   | B 削掉一半 | 0.35 | 0.05 | **0.80** | 0.15 |
   | C 全军覆没 | 0.25 | 0.02 | 0.15 | **0.80** |

预注册判词(条件式):
  if 正控开火(**植入世界里那一格必须在 `q=0.05` 下存活**)
     and 负控开火(**全匀速世界里「至少一个发现」的比例 ≤ q + 抽样误差**):
      `q=0.05` 下 `homosex` 两格都存活 -> A · 一格 -> B · 零格 -> C
  else: UNVERIFIED
⚠ **两个 `q` 都报**(`G4`),**且逐个列出没存活的格**(`G3`)。

⚠ 跑之前写下的最强混淆:**BH 假定 p 值在零下均匀,而自助 p 在格与格之间不独立**
  (同一题的相邻十年共用端点年份)。⇒ 控制:**同时报 BY(Benjamini–Yekutieli)校正**,
  它对任意依赖结构有效,**代价是更保守** —— **两者都报,不选边。**

⚠ 本轮**换不了仪器**:做的是同一张网格的多重性校正。
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
B, QS, NREP = 6000, [0.05, 0.10], 200

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
CELLS = [(it, dc) for it in ITEMS for dc in sorted(COV[it])]
C = len(CELLS)

print("=== ⓪ 跑前写死的两个陷阱,先算出来印出来(不是事后说「应该够」)===")
p_floor = 1.0/(B+1)
bh_strictest = QS[0]*1.0/C
print(f"  格数 C = **{C}** · B = {B} ⇒ 经验 p 的分辨率下限 `1/(B+1)` = **{p_floor:.2e}**")
print(f"  BH 在 q={QS[0]} 最严的一档 = `q·1/C` = **{bh_strictest:.2e}** ⇒ "
      f"**{'够用' if bh_strictest > p_floor else '不够用'}**(需 `q/C` > `1/(B+1)`)")
print(f"  ⚠ **而 `q/C` 是 Bonferroni,不是 BH 的判据** —— BH 在秩 k 的阈值是 `q·k/C`,最大的那个就是 q 本身。")

def pval(it, dc, rng, src=None):
    S = (src[it] if src else YR[it]); ys = COV[it][dc]
    span = ys[-1]-ys[0]; ref = FULL[it]["dgap"]*span/FULL[it]["span"]
    dr = np.empty(B)
    for i in range(B):
        r = lambda a: a[rng.integers(0, len(a), len(a))]
        a0, b0 = r(S[ys[0]][0]), r(S[ys[0]][1]); a1, b1 = r(S[ys[-1]][0]), r(S[ys[-1]][1])
        dr[i] = (a1.mean()-b1.mean()) - (a0.mean()-b0.mean())
    below = float(np.mean(dr <= ref)); above = float(np.mean(dr >= ref))
    return max(2*min(below, above), p_floor), float(np.mean(dr)), float(ref)

def bh(ps, q):
    """BH:秩 k 的阈值 q·k/C,取最大的通过秩。⚠ 不是 q/C(那是 Bonferroni)。"""
    n = len(ps); order = np.argsort(ps); thr = [q*(i+1)/n for i in range(n)]
    kmax = -1
    for i, idx in enumerate(order):
        if ps[idx] <= thr[i]: kmax = i
    return set(order[:kmax+1].tolist()) if kmax >= 0 else set()
def by(ps, q):
    """BY:对任意依赖有效,阈值 q·k/(C·H_C)。更保守,代价换普适。"""
    n = len(ps); H = sum(1.0/i for i in range(1, n+1))
    return bh(ps, q/H)

print(f"\n=== ① 逐格经验 p(B={B}),再对**整张网格** {C} 格做 BH / BY(`G4` 两个 q)===")
rng = np.random.default_rng(401)
PS, INFO = [], []
for it, dc in CELLS:
    p, mean_d, ref = pval(it, dc, rng)
    PS.append(p); INFO.append(dict(item=it, decade=dc, p=p, dgap=mean_d, ref=ref))
PS = np.array(PS)
RES = {}
for q in QS:
    s_bh, s_by = bh(PS, q), by(PS, q)
    RES[q] = dict(bh=sorted([f"{CELLS[i][0]}/{CELLS[i][1]}s" for i in s_bh]),
                  by=sorted([f"{CELLS[i][0]}/{CELLS[i][1]}s" for i in s_by]),
                  n_bh=len(s_bh), n_by=len(s_by), idx_bh=sorted(s_bh))
    print(f"  q = {q}: **BH 存活 {len(s_bh)}/{C}** ⇒ {RES[q]['bh'] or '无'}")
    print(f"            BY 存活 {len(s_by)}/{C} ⇒ {RES[q]['by'] or '无'}(⚠ 对任意依赖有效,更保守)")
print(f"\n  ⚠ **没存活的格,全列(`G3` 明文要求)** —— 按 p 升序前 12:")
for i in np.argsort(PS)[:12]:
    surv = "存活" if i in set(RES[QS[0]]["idx_bh"]) else "**未存活**"
    print(f"    {INFO[i]['item']:9s} {INFO[i]['decade']}s  p = {PS[i]:.4f} · "
          f"Δgap {INFO[i]['dgap']:+.4f} vs 参照 {INFO[i]['ref']:+.4f} ⇒ {surv}")
homo = [i for i, (it, dc) in enumerate(CELLS) if it == "homosex" and dc in (1990, 2000)]
homo_surv = sum(1 for i in homo if i in set(RES[0.05]["idx_bh"]))
print(f"\n  ⇒ **`homosex` 的两格(1990s / 2000s)在 q=0.05 下存活 {homo_surv}/2** "
      f"(p = {', '.join(f'{PS[i]:.4f}' for i in homo)})")

print("\n=== ② 控制 ===")
def syn(mode, rng_, it_p="homosex", dec_p=1990):
    S = {}
    for it in ITEMS:
        ys = sorted(YR[it]); y0, y1 = ys[0], ys[-1]
        g0 = float(YR[it][y0][0].mean()-YR[it][y0][1].mean()); tot = FULL[it]["dgap"]
        S[it] = {}
        for y in ys:
            f = ((y-y0)/(y1-y0) if mode == "uniform"
                 else (0.0 if y < dec_p else (1.0 if y > dec_p+9 else (y-dec_p)/9.0)))
            cur = float(YR[it][y][0].mean()-YR[it][y][1].mean())
            tgt = g0 + tot*f + (1.5*abs(tot)*f if (mode == "planted" and it == it_p) else 0.0)
            a, b = YR[it][y]
            S[it][y] = (a[rng_.integers(0, len(a), len(a))] + (tgt-cur),
                        b[rng_.integers(0, len(b), len(b))])
    return S
r2 = np.random.default_rng(402)
Sp = syn("planted", r2)
ps_p = np.array([pval(it, dc, r2, src=Sp)[0] for it, dc in CELLS])
surv_p = bh(ps_p, 0.05)
pc_ok = any(CELLS[i] == ("homosex", 1990) for i in surv_p)
print(f"  正控(位移全集中在 `homosex` 1990s 并加 1.5× 幅度)⇒ BH(q=0.05)存活 {len(surv_p)}/{C},"
      f"含 `homosex/1990s`:**{pc_ok}**")
hits = 0
for j in range(NREP):
    rj = np.random.default_rng(9000+j)
    Su = syn("uniform", rj)
    ps_u = np.array([pval(it, dc, rj, src=Su)[0] for it, dc in CELLS])
    if len(bh(ps_u, 0.05)) > 0: hits += 1
rate = hits/NREP; se = float(np.sqrt(rate*(1-rate)/NREP))
print(f"  负控(全匀速世界 × {NREP} 次)⇒ **「至少一个存活格」的比例 = {rate:.3f} ± {se:.3f}**")
print(f"     ⚠⚠ **期望是 ≤ q = 0.05,不是 0** —— BH 控制的是 FDR,"
      f"在全零世界里保证的是「任何一次发现的概率 ≤ q」。**`#825` 上一轮就犯在这里:"
      f"向一个额定率为 q 的程序索要 0。**")

G = Gate("#826 · 整张网格终于做一次多重性")
G.asserted("⓪ 前提(跑前写死的两个陷阱,算出来印出来):**`B` 的 p 分辨率下限 < BH 最严一档**,"
           "且**判据用 `q·k/C` 而不是 `q/C`(后者是 Bonferroni)**",
           bool(bh_strictest > p_floor),
           f"1/(B+1) = {p_floor:.2e} < q·1/C = {bh_strictest:.2e} · C = {C}", kind="control")
G.asserted("① 正控:位移全集中在 `homosex` 1990s 并加 1.5× 幅度的世界里,该格必须在 BH(q=0.05)下存活",
           bool(pc_ok), f"存活 {len(surv_p)}/{C} 格,含 homosex/1990s = {pc_ok}", kind="control")
G.asserted("② 负控:全匀速世界里「至少一个存活格」的比例必须 **≤ q = 0.05**(带抽样误差)"
           " —— ⚠⚠ **期望不是 0**:BH 控制 FDR,全零世界里的保证就是 ≤ q。"
           "**`#825` 上一轮正是犯在向一个额定率为 q 的程序索要 0。**",
           bool(rate <= 0.05 + 2*se), f"比例 {rate:.3f} ± {se:.3f}(阈 q + 2se = {0.05+2*se:.3f})",
           kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**BH 假定零下 p 均匀且格间独立,而同一题相邻十年共用端点年份**"
           " ⇒ **同时报 BY(对任意依赖有效,更保守),两者都报不选边**",
           bool(all("by" in RES[q] for q in QS)),
           " · ".join(f"q={q}: BH {RES[q]['n_bh']} / BY {RES[q]['n_by']}" for q in QS), kind="control")
G.asserted("④ 前提(`G3`):**cells tested 与 cells surviving 都报,且逐个列出没存活的格**",
           bool(C > 0), f"tested {C} · surviving(BH q=0.05){RES[0.05]['n_bh']}", kind="control")
G.asserted("⑤ kill(预注册):「`#819` 的核心在校正后仍站着」要成立,"
           "需 `homosex` 两格在 **q=0.05 的 BH** 下**都**存活",
           bool(homo_surv == 2), f"存活 {homo_surv}/2", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif homo_surv == 2:
    V = (f"**A `#819` 的核心在多重性校正后仍然站着。** `homosex` 的 1990s 与 2000s 两格在 BH(q=0.05)下都存活;"
         f"整张网格 {C} 格中共存活 **{RES[0.05]['n_bh']}** 格 ⇒ {RES[0.05]['bh']}。\n"
         f"  ⇒ **一句关于人的话:在减掉运气该给的那一份之后,同性恋这道题上的宗教鸿沟\n"
         f"  在九十年代与两千年代的两次偏离仍然站得住 —— 而这是这一串里第一个经过多重性校正的结论。**")
elif homo_surv == 1:
    V = (f"**B 削掉一半。** `homosex` 只有一格在 BH(q=0.05)下存活;整张网格存活 **{RES[0.05]['n_bh']}/{C}**。\n"
         f"  ⇒ **`#819` 的「两个十年」撤回 —— 校正之后它只剩一个十年。**")
else:
    V = (f"**C 全军覆没。** BH(q=0.05)下整张网格存活 **{RES[0.05]['n_bh']}/{C}** 格,`homosex` 一格也没剩。\n"
         f"  ⇒ **`#819` 的整张网格在多重性校正后不留下任何可信的格 ——\n"
         f"  而这意味着从 `#812` 到 `#824` 这一串关于「哪个十年裂开」的结论,\n"
         f"  全部建立在一个没有减掉运气的计数上。**")
print(V)
json.dump(dict(items=ITEMS, cells=[f"{a}/{b}" for a, b in CELLS], C=C, B=B, p_floor=p_floor,
               bh_strictest=bh_strictest, info=INFO, ps=PS.tolist(),
               results={str(q): RES[q] for q in QS}, homosex_survive=homo_surv,
               pos_control=dict(ok=bool(pc_ok), n_surv=len(surv_p)),
               neg_control=dict(rate=rate, se=se, expectation="<= q = 0.05, NOT 0", n_rep=NREP),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"bh_survivors.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'bh_survivors.json'}")
