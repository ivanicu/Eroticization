"""#818 · E03·A63·R257 —— 虔诚者是停住了,还是往回走了?

`#817` 并排看两个十年时冒出一个数:**两千年代虔诚层的世代内态度项是 −0.063** ——
而九十年代是 **+0.166**。**我当轮就写了「它没有区间,登记为下一轮要测的,本轮不当结论」。** 本轮测它。

⚠⚠ **而这一轮问的东西,和这条线索之前所有轮次都不同,这一点要先说清楚:**
   **`#802` 以来每一个数都是关于「两层之差」的** —— `r`、`Δgap`、`explained`、构成项份额,全是差。
   **一个差可以在两群人谁都没动的情况下变化,也可以在两群人都在动的情况下不变。**
   **而「虔诚的人自己做了什么」是一句独立的话,这条线索从来没有直接说过。**
   ⇒ **本轮的估计量只有一层:虔诚层自己的、世代内的态度变化,逐十年,带区间。**

G1 估计量:**`Σ_c w̄_c·(m_c1 − m_c0)`,只在虔诚层内,逐十年。**
   ⚠ **用世代内的项而不是原始均值差**,因为 `#816` 已经证明原始均值差里混着世代替换 ——
   **而「他们改了主意没有」问的正是世代内那一项。**

⚠ `matters`(`#811` 强制显式给):**0.10**。理由写下来:
   量表 1–4(跨度 3),而同期非虔诚层走了 +0.667 与 +0.360。
   **0.10 大约是「其余人一个十年走的六分之一」** —— **小于它的移动,不值得对一个人说「他们变了」。**
   **这是我选的,是一个关于「多大才算数」的判断,不是数据给的。**

⚠⚠ **`#817`② 的修正,本轮第一次执行:负控的容差必须事先定死,不许等于量出来的噪声。**
   `#817` 里我写了 `tol = max(噪声半宽, 1e-3)` ⇒ **比值恒等于 1.00×,那条检查永远不会开火。**
   ⇒ 本轮**先写死 `NC_TOL = 0.05`**,再去量噪声;
   **若量出来的噪声半宽 > 0.05,这条负控就判失败 —— 那是它应该做的事,不是我该去调的数。**

三个世界:
   A **停住了**:2000s 的区间含 0 且窄(`TIGHT_NULL`)⇒ **虔诚者在两千年代没有再改主意。**
   B **往回走了**:2000s 的区间排除 0 且为负 ⇒ **他们变得更保守了** ——
     **那是这条线索里第一句「有人朝反方向走」的话,而它比「鸿沟张开」强得多。**
   C **分辨不出**:区间含 0 且宽 ⇒ 登记功效不足,**不硬判**(`#811` 的三值)。

预测矩阵:
   | 世界 | 现在 | 2000s 排除 0 且负 | 2000s 窄零 | 2000s 宽零 |
   | A 停住 | 0.40 | 0.05 | **0.85** | 0.20 |
   | B 回走 | 0.25 | **0.90** | 0.03 | 0.10 |
   | C 分不出 | 0.35 | 0.05 | 0.12 | **0.70** |

预注册判词(条件式):
  if 正控开火(**在虔诚层里植入一个已知的世代内位移,必须取回**)
     and 负控开火(**世代内态度冻结的世界,该项必须 ≈ 0,容差 `NC_TOL = 0.05` 事先定死**):
      两个十年各自用 `#811` 三值判(参照 0,`matters` = 0.10);**并排报,不合并**
  else: UNVERIFIED
⚠ **凡 `UNRESOLVED` 必须同时印出它与哪些参照相容**(`#812`③)。

⚠ 跑之前写下的最强混淆:**虔诚层是按年份内三分位定义的** ⇒ 两个年份的「虔诚层」不是同一群人,
  **而「他们改了主意没有」这句话预设了它是。**
  ⇒ 控制:**同时报虔诚层自己的 `REL` 分布**(均值与两个四分位)在两个年份的差异 ——
  **若那个分布本身在漂,这句话的主语就在漂,而本设计修不了它,只能量它。**

⚠ 硬规则①:先打印两个十年、虔诚层、每个世代箱的 n。⚠ 换不了仪器。⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(257)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK = "homosex", 4
WINDOWS = {"1990s": (1990, 1998), "2000s": (2000, 2008)}
MATTERS, NC_TOL, B, NREP = 0.10, 0.05, 3000, 200

d = pd.read_stata(gp, columns=["year", "cohort", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("cohort", (1880, 2010))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
BASE = REL.dropna(subset=[IT, "cohort"]).copy()
BASE["gen"] = (BASE.cohort//10*10).astype(int)

def win(y0, y1, k=2):
    W = BASE[BASE.year.isin([y0, y1]) & (BASE.k == k)].copy()
    keep = [g for g, n in W.groupby("gen").size().items() if n >= 40]
    return W[W.gen.isin(keep)], sorted(keep)

def within(df, y0, y1):
    """只在这一层内:Σ w̄_c·(m_c1 − m_c0)。"""
    a, b = df[df.year == y0], df[df.year == y1]
    gens = sorted(set(a.gen) & set(b.gen))
    if not gens or not len(a) or not len(b): return None
    w0 = np.array([len(a[a.gen == g])/len(a) for g in gens])
    w1 = np.array([len(b[b.gen == g])/len(b) for g in gens])
    m0 = np.array([a[a.gen == g][IT].mean() for g in gens])
    m1 = np.array([b[b.gen == g][IT].mean() for g in gens])
    if np.isnan(m0).any() or np.isnan(m1).any(): return None
    return float((((w0+w1)/2)*(m1-m0)).sum())

print("=== ⓪ 硬规则①:两个十年 · 虔诚层 · 世代箱 n ===")
DATA = {}
for lab, (y0, y1) in WINDOWS.items():
    W, keep = win(y0, y1)
    n = W.groupby(["year", "gen"]).size()
    DATA[lab] = (W, keep, y0, y1)
    print(f"  {lab} {y0}→{y1} · 世代箱 {len(keep)} 个 · {y0} n={len(W[W.year==y0]):,} · {y1} n={len(W[W.year==y1]):,}")
    print(f"        每箱 {y0}/{y1}:" + " · ".join(f"{g}s {n.get((y0,g),0)}/{n.get((y1,g),0)}" for g in keep))

print(f"\n=== ① 虔诚层自己的世代内态度变化(⚠ **只有一层,不是两层之差** · B={B} · `matters`={MATTERS})===")
ROWS = []
for lab, (W, keep, y0, y1) in DATA.items():
    pt = within(W, y0, y1)
    dr = np.array([within(W.iloc[RNG.integers(0, len(W), len(W))], y0, y1) for _ in range(B)], float)
    dr = dr[np.isfinite(dr)]
    lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
    v = Gate.interval_verdict(lo, hi, 0.0, MATTERS)
    compat = [x for x in (-0.10, 0.0, 0.10, 0.166) if lo <= x <= hi]
    ROWS.append(dict(window=lab, y0=y0, y1=y1, point=float(pt), lo=lo, hi=hi, verdict=v,
                     compatible_with=compat, n_gens=len(keep)))
    print(f"  {lab}: **{pt:+.4f}** [{lo:+.4f}, {hi:+.4f}] ⇒ **{v}**"
          + (f"  ⚠ **相容于 {compat}**(`#812`③)" if v == "UNRESOLVED" else ""))

print("\n  ⚠ 跑前混淆的控制 —— 虔诚层自己的 `REL` 分布在两个年份漂了多少:")
for lab, (W, keep, y0, y1) in DATA.items():
    a, b = W[W.year == y0].REL, W[W.year == y1].REL
    print(f"    {lab}: {y0} 均值 {a.mean():+.3f}(Q1 {a.quantile(.25):+.3f} / Q3 {a.quantile(.75):+.3f}) → "
          f"{y1} 均值 {b.mean():+.3f}(Q1 {b.quantile(.25):+.3f} / Q3 {b.quantile(.75):+.3f})")
print("    ⚠ **若这个分布本身在漂,「他们改了主意没有」这句话的主语就在漂 —— 本设计修不了它,只能量它。**")

print("\n=== ② 控制(合成世界,同一条代码路径)===")
W9, _, y0_, y1_ = DATA["1990s"]
# ⚠⚠ 第一版**两条控制都错了,而两条错都是我自己的,不是仪器的**,并且两条都被闸当场抓住:
#  ① 正控**把 −0.300 叠加在真实数据上**,而真实数据本身已经有 +0.166 ⇒
#    期望取回的是 **−0.300 + 0.166 = −0.134**,实测 −0.1368 —— **仪器是对的,是我的期望值写错了。**
#    ⇒ 这是 `#812` 那条「destroy the structure under test, preserve everything else」的**第三次重犯**,
#    而 `#812` 已经把它写进账本了。**教训写下来并没有转移过来。**
#    ⇒ 改成**替换**:把 y1 每个人的值设成「自己世代在 y0 的均值 + 植入量」⇒ 世代内项**按构造恰好等于植入量**。
#  ② 负控「世代内冻结」使该项**代数上恰好为 0**(m_c1 ≡ m_c0 ⇒ Σw̄·0 = 0)⇒ 与 0 比 0,库判 DEGENERATE。
#    ⇒ **它不是一条随机控制,是一条代数恒等式的实现检查**(`realstat` 的算术陷阱)。
#    ⇒ 声明 `deterministic=True`,并**另加一条真正随机的负控**:y1 的值从**自己世代在 y0 的分布里重抽**
#      ⇒ 真值为 0 **而带抽样噪声** —— 那一条才量得出噪声半宽。
def syn(mode, shift=-0.30):
    H = W9.copy()
    gm = H[H.year == y0_].groupby("gen")[IT].mean()
    if mode == "frozen":                      # 代数恒等式:m_c1 ≡ m_c0
        H[IT] = [gm.get(g, np.nan) for g in H.gen]
        H = H.dropna(subset=[IT])
    elif mode == "planted":                   # ⚠ 替换,不是叠加
        m = H.year == y1_
        H.loc[m, IT] = [gm.get(g, np.nan) + shift for g in H.loc[m, "gen"]]
        H = H.dropna(subset=[IT])
    else:                                     # 'resampled':真值 0,带抽样噪声
        m = H.year == y1_
        pool = {g: H[(H.year == y0_) & (H.gen == g)][IT].to_numpy(float) for g in H.gen.unique()}
        H.loc[m, IT] = [RNG.choice(pool[g]) if len(pool.get(g, [])) else np.nan for g in H.loc[m, "gen"]]
        H = H.dropna(subset=[IT])
    return H
def ctl(mode, rep=NREP, shift=-0.30):
    v = []
    for _ in range(rep):
        H = syn(mode, shift)
        x = within(H.iloc[RNG.integers(0, len(H), len(H))], y0_, y1_)
        if x is not None: v.append(x)
    v = np.array(v, float)
    return float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
pc_m, pc_lo, pc_hi = ctl("planted")
id_m, id_lo, id_hi = ctl("frozen")            # 代数恒等式检查(确定性)
nc_m, nc_lo, nc_hi = ctl("resampled")         # 真正的随机负控
nc_half = (nc_hi-nc_lo)/2
print(f"  正控:在 {y1_} 年植入 −0.300 ⇒ 取回中位 **{pc_m:+.4f}** [{pc_lo:+.4f}, {pc_hi:+.4f}]")
print(f"  恒等式检查(世代内冻结,**代数上恰好为 0**)⇒ {id_m:+.6f} [{id_lo:+.6f}, {id_hi:+.6f}] —— **确定性,不是控制噪声的度量**")
print(f"  负控(y1 从自己世代在 y0 的分布重抽 ⇒ **真值 0 且带抽样噪声**)⇒ 中位 **{nc_m:+.4f}** "
      f"[{nc_lo:+.4f}, {nc_hi:+.4f}],噪声半宽 **{nc_half:.4f}**")
print(f"     ⚠⚠ **容差 `NC_TOL = {NC_TOL}` 是跑之前写死的(`#817`② 的修正)** —— "
      f"比值 **{NC_TOL/nc_half:.2f}×** ⇒ "
      f"{'过' if NC_TOL >= nc_half else '**这条负控判失败,而那正是它该做的事,不是我该去调的数**'}")

Gg = Gate("#818 · 虔诚者是停住了,还是往回走了")
Gg.asserted("① 正控:**替换式**植入一个已知的世代内位移 −0.300(而不是叠加在真实结构上)必须取回"
            " —— ⚠ **第一版是叠加,期望值写错了:−0.300 + 真实的 +0.166 = −0.134,实测 −0.1368,"
            "仪器是对的**(`#812` 那条教训的第三次重犯)",
            bool(abs(pc_m+0.30) < 0.10), f"取回 {pc_m:+.4f} [{pc_lo:+.4f}, {pc_hi:+.4f}]", kind="control")
# ⚠⚠ 这一条**不能写成 `identity_control`**,而这本身是本轮的一个发现:
#   世代内冻结时该项**代数上恰好为 0**,所以「观测」与「期望」都是 0 ⇒ 库判 DEGENERATE,而**库是对的**:
#   `#802` 已经给 `#770` 补过前提——**「比两个值而不是比差与零」要求那两个值非零**。
#   ⇒ 这里我想要的检查根本不是「两个量是同一个」,是「**这个量恰好是零**」——
#     **而一个恒等于零的量,它的检查只能是代码实现检查,不能是控制。**
#   ⇒ 写成 `asserted`,并**在行里说清楚它是什么**:不是证据,是代数被正确实现了。
Gg.asserted("①b 代码实现检查(**不是控制**):世代内冻结时该项**代数上恰好为 0**(m_c1 ≡ m_c0 ⇒ Σw̄·0 = 0)"
            " —— ⚠ **它不能写成 `identity_control`:两侧都恒为 0,库会正确地判它空洞("
            "`#802` 给 `#770` 补的前提是那两个值必须非零)。一个恒等于零的量,只能做实现检查,不能当控制。**",
            bool(abs(id_m) < 1e-12 and abs(id_lo) < 1e-12 and abs(id_hi) < 1e-12),
            f"{NREP} 次重复全部恰好为 0(跨度 [{id_lo:+.2e}, {id_hi:+.2e}])—— **代数被正确实现,这不是证据**",
            kind="control")
Gg.identity_control("② 负控(**真正随机的那一条**):y1 从自己世代在 y0 的分布重抽 ⇒ 真值 0 且带抽样噪声"
                    " —— ⚠⚠ **容差 0.05 事先写死,不等于量出来的噪声(`#817`② 的修正:"
                    "令 `容差 := 噪声` 会让比值恒等于 1.00×,那条检查永远不会开火)**",
                    observed=nc_m, expected=0.0, tol=NC_TOL, noise_half_width=nc_half,
                    what=f"{NREP} 次重复,95% 跨度 [{nc_lo:+.4f}, {nc_hi:+.4f}]")
Gg.asserted("③ 前提(跑前写下的混淆):虔诚层自己的 `REL` 分布在两个年份的漂移量并排印出,"
            "**如实登记本设计修不了它,只能量它**", True,
            " · ".join(f"{l} {DATA[l][0][DATA[l][0].year==DATA[l][2]].REL.mean():+.3f}→"
                       f"{DATA[l][0][DATA[l][0].year==DATA[l][3]].REL.mean():+.3f}" for l in DATA), kind="control")
Gg.asserted("④ 前提:`matters` 显式给出并写下理由(`#811`)", bool(MATTERS > 0),
            f"matters = {MATTERS} —— 约为其余人一个十年移动的六分之一,我选的", kind="control")
two = next(r for r in ROWS if r["window"] == "2000s")
Gg.asserted("⑤ kill(预注册):「虔诚者在两千年代往回走了」要成立,需 2000s 的区间**排除 0 且为负**",
            bool(two["verdict"] == "EXCLUDES" and two["hi"] < 0),
            f"2000s {two['point']:+.4f} [{two['lo']:+.4f}, {two['hi']:+.4f}] ⇒ {two['verdict']}", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
one = next(r for r in ROWS if r["window"] == "1990s")
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif two["verdict"] == "EXCLUDES" and two["hi"] < 0:
    V = (f"**B 往回走了。** 2000s **{two['point']:+.4f}** [{two['lo']:+.4f}, {two['hi']:+.4f}] ——\n"
         f"  ⇒ **这是这条线索里第一句「有人朝反方向走」的话,而它比「鸿沟张开」强得多:\n"
         f"  两千年代,虔诚的美国人在同性恋这道题上不是跟得慢,是自己往回退了。**")
elif two["verdict"] == "TIGHT_NULL":
    V = (f"**A 停住了。** 2000s **{two['point']:+.4f}** [{two['lo']:+.4f}, {two['hi']:+.4f}] —— **窄零。**\n"
         f"  ⇒ **一句关于人的话:九十年代虔诚者也在变宽容({one['point']:+.4f});\n"
         f"  到两千年代他们停住了 —— 不是往回走,是不再动。而其余人没有停。**\n"
         f"  ⚠ **「停住」与「往回走」是两句话,而本轮的区间只支持前一句。**")
else:
    V = (f"**C 分辨不出。** 2000s **{two['point']:+.4f}** [{two['lo']:+.4f}, {two['hi']:+.4f}] —— "
         f"**相容于 {two['compatible_with']}**(`#812`③)。\n"
         f"  ⇒ **「停住了」与「往回走了」这个设计分不开,如实说 —— 而 `#817` 那个 −0.063 因此\n"
         f"  仍然只是一次拆解里的读数,不是一个站得住的量。**")
print(V)
print(f"\n  ⚠ 对照:1990s **{one['point']:+.4f}** [{one['lo']:+.4f}, {one['hi']:+.4f}] ⇒ **{one['verdict']}**")
json.dump(dict(item=IT, matters=MATTERS, nc_tol=NC_TOL, B=B, n_rep=NREP, rows=ROWS,
               rel_drift={l: dict(y0=float(DATA[l][0][DATA[l][0].year==DATA[l][2]].REL.mean()),
                                  y1=float(DATA[l][0][DATA[l][0].year==DATA[l][3]].REL.mean())) for l in DATA},
               pos_control=dict(median=pc_m, lo=pc_lo, hi=pc_hi, planted=-0.30),
               identity_check=dict(median=id_m, lo=id_lo, hi=id_hi, deterministic=True),
               neg_control=dict(median=nc_m, lo=nc_lo, hi=nc_hi, half_width=nc_half,
                                tol_fixed_in_advance=NC_TOL, ratio=NC_TOL/nc_half, reference=0.0),
               admissible=adm, verdict=V, gate_ok=Gg.verdict()),
          open(OUT/"did_the_devout_turn_back.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'did_the_devout_turn_back.json'}")
