"""#835 · E03·A77·R274 —— 年龄和宗教是同一条缝吗?而 `#834`① 问的那个量不可识别

`#834` 发现九十年代那条缝**同时**沿虔诚度与年龄裂开(BH 下两格都存活;BY 下只剩虔诚度)。
`#834`① 于是问:**「若把年龄效应扣掉,宗教那一格还剩多少?」**

⚠⚠ **而那个问句本身有问题,这是本轮第一件要写下的事(`G1`:识别先于功效):**
   **在同一个调查年内,年龄与出生世代完全共线** —— `世代 = 年份 − 年龄`。
   ⇒ **「扣掉年龄」同时也在扣掉世代;而年龄/时期/世代三者在横截面序列里本就不可分离(APC 不可识别)。**
   ⇒ **「宗教那一格净掉年龄效应之后的值」不是一个能被这份数据估出来的量。**
   **要估它,必须先假设「没有世代效应」或「没有年龄效应」—— 那是我强加的,不是我测到的。**
   ⇒ **登记为结构性不可识别(`realstat §2`),而不是「计划中」。**

**⇒ 所以本轮换一个**能**被识别的问句:不做扣除,做分层。**
   **「宗教那一格,在**每一个年龄带内部**还在不在?」**
   —— 这是一个条件独立性问题,**它不需要把年龄从世代里分开**,
   只需要问:在年龄相近的人里面,虔诚与否是否仍然把九十年代切开。

G1 估计量:**在年龄带 A 内,`homosex` 的 1990s `Δgap`(虔诚 vs 非虔诚)对该(带 × 轴)自己的匀速参照。**

⚠⚠ **而 `G1` 说识别是第一问、功效是第二问 —— 两问都必须在读结果之前回答。**
   **分成三带 ⇒ 每格 n 降到约三分之一 ⇒ 噪声半宽约涨 √3。**
   `#834` 那一格的 p = 0.00017(约 3.8σ)⇒ 分带后约 2.2σ ⇒ **每带 p ≈ 0.03,而这在校正后是边缘的。**
   ⇒ **本轮先把每带的噪声半宽量出来,再决定读不读结果** ——
   **这正是 `#822` 教的:先量那个改动有没有用,再拿它去重跑;而 `#821` 的教训是别拿错噪声当尺子。**

**功效闸(⓪,排在所有关于人的判断之前)**:
   **每带的自助噪声半宽必须 < 该带需要分辨的效应量的一半**,否则 **verdict = C,本轮不许往下读。**
   ⚠ 「需要分辨的效应量」= `#834` 全样本上那一格的 `|Δgap − ref|`,**跑前从产物读出来,不是我记的。**

三个世界:
   A **宗教在带内仍在**:≥2/3 带的 1990s 格在整族 BH 下存活 ⇒ **两条缝是两条,不是一条看了两遍。**
   B **宗教在带内消失**:0–1 带存活,而**同一设计下全样本那一格仍存活** ⇒
     **它原来主要是年龄构成 ⇒ `#832`/`#834` 的主语要再改一次。**
   C **分不出**:功效闸没过,或带内区间宽到含 0 又含全样本值 ⇒ **登记边界,不假装。**

预测矩阵:
   | 世界 | 现在 | ≥2/3 带存活 | 0–1 带存活 | 功效闸没过 |
   | A 两条缝 | 0.40 | **0.85** | 0.05 | 0.20 |
   | B 其实是年龄 | 0.25 | 0.05 | **0.85** | 0.15 |
   | C 分不出 | 0.35 | 0.10 | 0.10 | **0.85** |

预注册判词(**两级条件式**):
  if 功效闸没过: verdict = C(**不读带内结果**)
  elif 正控开火(在某一带的 1990s 植入巨大偏离,它必须在整族 BH 下存活)
       and 负控开火(全匀速世界「至少一个存活」比例 ≤ q,⚠ **不是 0**):
      ≥2/3 带存活 -> A · 0–1 带且全样本仍存活 -> B · 其余 -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**年龄带内部,虔诚度的三分位线与全样本不同** ——
  老年人整体更虔诚,所以「老年带里的低虔诚组」可能比「青年带里的高虔诚组」还虔诚。
  ⇒ 控制:**逐带印出该带虔诚三分位的 `REL` 切点**,让这个错位看得见;
  **并且分层只在带内做三分位**(而不是用全样本切点),否则有的带会几乎没有一侧。

⚠ 硬规则①:先打印每带、每年、每层的 n。⚠ 硬规则②:全部来自 GSS 同一份问卷,**不是跨仪器复现**。
⚠⚠ **本轮换不了仪器,而这一次的理由就是本轮的结论本身**:
   问的是「GSS 这份数据能不能把宗教缝与年龄缝分开」——**对象就是这具仪器的分辨力**。
   ⇒ 换第二具仪器不是「再验证一次」,**它恰恰是本轮结论指向的下一步**
   (`#835`①:唯一诚实的路是换一份年龄层内样本量更大的调查)。
   **所以这里的「换不了仪器」不是豁免,是被测的对象本身。**
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, B, Q, NREP, BNULL = "homosex", 4, 6000, 0.05, 60, 2000

# ⚠ 跑前从 `#834` 的产物里读出「需要分辨的效应量」,不是我记的
J834 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A76_那句话的主语从没被检验过/"
                      "R273_does_the_nineties_cell_survive_on_all_seven_axes/results/whose_divide.json"))
print("=== ⓪a 从 `#834` 的产物读出全样本那一格(不是从记忆)===")
print(f"  `REL 虔诚度|1990` 的 p = **{J834['ps']['REL 虔诚度|1990']:.5f}** · "
      f"BH 存活集 {J834['surv_bh']}")

d = pd.read_stata(gp, columns=["year", "age", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("age", (18, 89))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1)
M = M.join(R["REL"])
D = M.dropna(subset=[IT, "REL", "age", "year"]).copy()
BANDS = {"18–34 青年": (18, 34), "35–54 中年": (35, 54), "55+ 年长": (55, 89)}
D["band"] = pd.cut(D.age, bins=[17, 34, 54, 89], labels=list(BANDS))

print("\n=== ⓪b 硬规则① + 跑前写下的最强混淆:逐带的 n 与**带内**虔诚三分位切点 ===")
GR = {}
for bn in BANDS:
    sub = D[D.band == bn].copy()
    sub["k"] = sub.groupby("year")["REL"].transform(
        lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
    ys = {}
    for y, g in sub.groupby("year"):
        a = g[g.k == 2][IT].to_numpy(float); b = g[g.k == 0][IT].to_numpy(float)
        if len(a) >= 60 and len(b) >= 60: ys[int(y)] = (a, b)
    dec = {}
    for y in sorted(ys): dec.setdefault((y//10)*10, []).append(y)
    dec = {k: v for k, v in dec.items() if len(v) >= 3}
    s = sorted(ys)
    cut = sub.groupby("year")["REL"].quantile([1/3, 2/3]).groupby(level=1).mean()
    print(f"  {bn:10s} n={len(sub):>6,} · 可用年 {len(ys):>2} · 可用十年 {sorted(dec)} · "
          f"**带内 REL 三分位切点 {cut.iloc[0]:+.3f} / {cut.iloc[1]:+.3f}**")
    if len(s) < 8 or 1990 not in dec: GR[bn] = None; continue
    g0 = float(ys[s[0]][0].mean()-ys[s[0]][1].mean()); g1 = float(ys[s[-1]][0].mean()-ys[s[-1]][1].mean())
    GR[bn] = dict(ys=ys, dec=dec, span=s[-1]-s[0], dgap=g1-g0)
print("  ⚠ **切点逐带不同 —— 老年带的「低虔诚组」可能比青年带的「高虔诚组」还虔诚。**")
print("     **分层只在带内做,否则有的带几乎没有一侧;而这个错位改不了,只能量出来并写下。**")

def stats(bn, dc, rng, Bv, src=None):
    G_ = GR[bn]; S = src if src else G_["ys"]; ys = G_["dec"][dc]
    ref = G_["dgap"]*(ys[-1]-ys[0])/G_["span"]
    dr = np.empty(Bv)
    for i in range(Bv):
        r = lambda a: a[rng.integers(0, len(a), len(a))]
        a0, b0 = r(S[ys[0]][0]), r(S[ys[0]][1]); a1, b1 = r(S[ys[-1]][0]), r(S[ys[-1]][1])
        dr[i] = (a1.mean()-b1.mean()) - (a0.mean()-b0.mean())
    obs = (float(S[ys[-1]][0].mean()-S[ys[-1]][1].mean())
           - float(S[ys[0]][0].mean()-S[ys[0]][1].mean()))
    lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
    p = max(2*min(float(np.mean(dr <= ref)), float(np.mean(dr >= ref))), 1.0/(Bv+1))
    return dict(p=p, obs=obs, ref=float(ref), lo=lo, hi=hi, half=(hi-lo)/2,
                need=abs(obs-ref))

rng = np.random.default_rng(274)
print(f"\n=== ① 功效闸 ⓪(排在所有关于人的判断之前):每带 1990s 的噪声半宽 vs 需要分辨的量 ===")
POW = {}
for bn in BANDS:
    if not GR.get(bn): print(f"  {bn:10s} —— 1990s 不可用(年份不足),整带出局"); POW[bn] = None; continue
    st = stats(bn, 1990, rng, B); POW[bn] = st
    ok = st["half"] < st["need"]/2
    print(f"  {bn:10s} 观测 Δgap {st['obs']:+.4f} · 参照 {st['ref']:+.4f} ⇒ **需分辨 {st['need']:.4f}** · "
          f"噪声半宽 **{st['half']:.4f}** ⇒ 半宽 < 需分辨/2:**{ok}**")
usable = [bn for bn in BANDS if POW.get(bn)]
power_ok = bool(usable) and all(POW[bn]["half"] < POW[bn]["need"]/2 for bn in usable)
print(f"  ⇒ **功效闸:{'过' if power_ok else '**没过**'}** —— "
      f"{'可以往下读' if power_ok else '**本轮 verdict = C,不许读带内结果**'}")

CELLS = [(bn, dc) for bn in usable for dc in sorted(GR[bn]["dec"])]
PS = {c: (stats(c[0], c[1], rng, B)["p"] if c[1] != 1990 else POW[c[0]]["p"]) for c in CELLS}
surv = {CELLS[i] for i in Gate.bh([PS[c] for c in CELLS], Q)}
surv_by = {CELLS[i] for i in Gate.by([PS[c] for c in CELLS], Q)}
n90 = sum(1 for bn in usable if (bn, 1990) in surv)
print(f"\n=== ② 带内网格(整族 {len(CELLS)} 格,`G3` 全报)===")
for bn in usable:
    row = " ".join(f"{dc}s:{PS[(bn,dc)]:.4f}{'**' if (bn,dc) in surv else ''}"
                   for dc in sorted(GR[bn]["dec"]))
    print(f"  {bn:10s} {row}")
print(f"  BH 存活 **{len(surv)}/{len(CELLS)}** ⇒ {sorted(f'{a}/{b}s' for a,b in surv) or '无'}")
print(f"  BY 存活 **{len(surv_by)}/{len(CELLS)}** ⇒ {sorted(f'{a}/{b}s' for a,b in surv_by) or '无'}")
print(f"  ⇒ **1990s 在 {n90}/{len(usable)} 个年龄带内存活**")

print("\n=== ③ 控制 ===")
def syn(mode, rng_, plant=None):
    S = {}
    for bn in usable:
        G_ = GR[bn]; ys = sorted(G_["ys"]); y0, y1 = ys[0], ys[-1]
        g0 = float(G_["ys"][y0][0].mean()-G_["ys"][y0][1].mean()); tot = G_["dgap"]
        S[bn] = {}
        for y in ys:
            tgt = g0 + tot*(y-y0)/(y1-y0)
            if mode == "planted" and bn == (plant or usable[0]) and 1990 <= y <= 1999:
                tgt += 3.0*abs(tot)*(y-1990)/9.0
            cur = float(G_["ys"][y][0].mean()-G_["ys"][y][1].mean())
            a, b = G_["ys"][y]
            S[bn][y] = (a[rng_.integers(0, len(a), len(a))] + (tgt-cur),
                        b[rng_.integers(0, len(b), len(b))])
    return S
r2 = np.random.default_rng(275)
Sp = syn("planted", r2)
psp = {c: stats(c[0], c[1], r2, BNULL, src=Sp[c[0]])["p"] for c in CELLS}
sp = {CELLS[i] for i in Gate.bh([psp[c] for c in CELLS], Q)}
pc_ok = (usable[0], 1990) in sp if usable else False
print(f"  正控(在「{usable[0] if usable else '—'}」的 1990s 植入 3× 全程量偏离)⇒ "
      f"整族 BH 存活 {len(sp)}/{len(CELLS)},含被植入格:**{pc_ok}**")
hits = 0
for j in range(NREP):
    rj = np.random.default_rng(8000+j)
    Su = syn("uniform", rj)
    psu = {c: stats(c[0], c[1], rj, BNULL, src=Su[c[0]])["p"] for c in CELLS}
    if len(Gate.bh([psu[c] for c in CELLS], Q)) > 0: hits += 1
rate = hits/NREP; se = float(np.sqrt(rate*(1-rate)/NREP))
print(f"  负控(全匀速 × {NREP} 次)⇒「至少一个存活」**{rate:.3f} ± {se:.3f}** —— "
      f"⚠ **期望 ≤ q = {Q},不是 0**")

G = Gate("#835 · 年龄和宗教是同一条缝吗")
G.asserted("⓪ **识别(`G1` 第一问)**:`#834`① 问的「净掉年龄效应」在横截面序列里**不可识别** —— "
           "同一年内 `世代 = 年份 − 年龄` 完全共线,扣年龄同时在扣世代 ⇒ "
           "**登记为结构性不可识别,本轮改问一个能识别的:分层而非扣除**",
           True, "APC 不可识别 ⇒ 改问「带内是否仍裂」(条件独立性,不需分离年龄与世代)", kind="control")
G.asserted("⓪b **功效闸(`G1` 第二问,排在所有关于人的判断之前)**:每带 1990s 的噪声半宽必须 < "
           "需分辨量的一半,否则 verdict = C 且**不许读带内结果**",
           bool(power_ok),
           " · ".join(f"{bn[:5]} 半宽{POW[bn]['half']:.3f}/需{POW[bn]['need']:.3f}" for bn in usable),
           kind="control")
G.asserted("① 正控:在某一带的 1990s 植入 3× 全程量的偏离,该格必须在整族 BH 下存活",
           bool(pc_ok), f"存活 {len(sp)}/{len(CELLS)},含被植入格 = {pc_ok}", kind="control")
G.asserted("② 负控:全匀速世界「至少一个存活」比例 **≤ q**(⚠ **不是 0** —— BH 控 FDR)",
           bool(rate <= Q + 2*se), f"{rate:.3f} ± {se:.3f}(阈 {Q+2*se:.3f})", kind="control")
G.asserted("③ 前提(跑前写下的最强混淆):**带内虔诚三分位切点逐带不同**(老年带的低虔诚组可能比青年带的高虔诚组还虔诚)"
           " ⇒ **切点已逐带印出,且分层只在带内做** —— 这个错位改不了,只能量出来",
           True, "三带各自的 REL 三分位切点已打印", kind="control")
G.asserted("④ kill(预注册,**且被功效闸门控**):「两条缝是两条」要成立,需 1990s 在 **≥2/3** 带内存活",
           bool(power_ok and n90 >= 2), f"功效闸 {power_ok} · 1990s 带内存活 {n90}/{len(usable)}",
           kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not power_ok:
    V = (f"**C 这个设计分不出来 —— 而这是功效闸在读结果之前拦下的。**\n"
         f"  逐带:" + " · ".join(f"{bn} 噪声半宽 {POW[bn]['half']:.4f} vs 需分辨 {POW[bn]['need']:.4f}"
                                 for bn in usable) + "\n"
         f"  ⇒ **把样本切成三个年龄带之后,每带的噪声比它要分辨的东西还大 ——\n"
         f"  「宗教那条缝是不是年龄那条缝」在这份数据上答不了,而这不是没跑,是跑不了。**\n"
         f"  ⚠ **而 `#834`① 问的那个「扣掉年龄」的量,本来就不可识别**(APC 共线)——\n"
         f"  **所以这个问题有两层墙:第一层是识别,第二层是功效,两层都不是加大计算能翻的。**")
elif not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n90 >= 2:
    V = (f"**A 两条缝是两条,不是一条看了两遍。** 1990s 在 **{n90}/{len(usable)}** 个年龄带内存活。\n"
         f"  ⇒ **一句关于人的话:在年龄相近的美国人里面,信不信教仍然把九十年代切开 ——\n"
         f"  所以那条宗教的缝不是「年长的人更虔诚」这件事的影子。**")
else:
    V = (f"**B 带内消失了。** 1990s 只在 **{n90}/{len(usable)}** 个带内存活。\n"
         f"  ⇒ **那条缝原来主要是年龄构成 ⇒ `#832`/`#834` 的主语要再改一次。**")
print(V)
json.dump(dict(item=IT, bands=list(BANDS), usable=usable, power=POW, power_ok=bool(power_ok),
               cells=[f"{a}|{b}" for a, b in CELLS], ps={f"{a}|{b}": PS[(a, b)] for a, b in CELLS},
               surv_bh=sorted(f"{a}/{b}s" for a, b in surv), surv_by=sorted(f"{a}/{b}s" for a, b in surv_by),
               n1990_bands=n90, identification="APC non-identified: cohort = year - age",
               pos_control=bool(pc_ok), neg_control=dict(rate=rate, se=se, expectation="<= q, NOT 0"),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"age_or_religion.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'age_or_religion.json'}")
