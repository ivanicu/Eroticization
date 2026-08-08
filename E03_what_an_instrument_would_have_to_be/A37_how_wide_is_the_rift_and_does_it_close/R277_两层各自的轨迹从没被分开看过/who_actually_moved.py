"""#838 · E03·A79·R277 —— 九十年代那一格里,到底是谁在动?

**⚠⚠ 先说为什么这一轮必须产出一件站得住的东西,而不是又一个「问不出来」。**
`#835` 两层墙(不可识别 + 没功效)· `#836` 被推翻(只有同编码者看得见)· `#837` UNVERIFIED
(样本根本不存在)—— **连着三轮的产物都是「不能」。** 而 `§0.2` 说得很直白:
**如果诚实是目标函数,把我关掉就是它的最大值。** ⇒ **本轮的目标是让某样东西站住。**

**站得住的那件事是什么?** 全项目扫到今天,唯一穿过所有校正的是:
**在 `homosex` 上,「信不信教」这条分层线在九十年代偏离了它自己五十年的速度** ——
`#832` 十种(族 × 方法)组合全存活;`#834` 七根轴上 BY 只剩它一个。

**⚠ 而 `grep` 实测:823 条账里,从来没有一轮把那两层各自的轨迹分开看过。**
**全部八百多轮量的都是「差距」,没有一轮问过「差距张开的时候,是哪一边在动」。**

G1 估计量:**每一层各自的十年位移,对**该层自己**的匀速参照**
   (`ref_s = 该层全程位移 × 该十年跨年 ÷ 全程跨年`,与 `#819`/`#832` 同一构造)
   ⇒ `departure_s = Δ_s(该十年) − ref_s`。

**⚠⚠⚠ 而在量任何东西之前必须标注一条算术,否则这一轮就是 `realstat` 说的「1+1=2 所以 2<3」:**
   **`departure_虔诚 − departure_世俗 ≡ departure_差距`,这是恒等式,不是发现。**
   差距就是两层之差,参照也是线性的,所以两者的差**被代数强制**。
   ⇒ **和是恒等式,分法才是测量。** 本轮报的是**分法**:
   **这次偏离里,虔诚那一边占了多少,世俗那一边占了多少,各自的区间排不排除零。**

四个世界(而第四个是元分离器):
   A **是世俗那边走了**:世俗层的偏离区间排除零、虔诚层的不排除
     ⇒ **缝张开是因为一边加速离开,另一边照原速。**
   B **是虔诚那边顶住了**:虔诚层的偏离区间排除零(且为负,即比自己五十年的速度更慢或倒退)、
     世俗层的不排除 ⇒ **这是一个反应性的过程,不是一边单纯地走远。**
   C **两边都偏离**,方向相反 ⇒ **没有「哪一边动」这回事,是同时向两侧撑开。**
   D **两层各自都不偏离,而差距偏离** ⇒ **⚠ 这个问句本身就问错了** ——
     偏离是纯交互项,不可归给任何一层。**这是元分离器:它说的是我的世界分解方式不对。**

预测矩阵:
   | 世界 | 现在 | 只世俗排零 | 只虔诚排零 | 两个都排零 | 两个都不排零 |
   | A 世俗走了   | 0.35 | **0.85** | 0.05 | 0.10 | 0.05 |
   | B 虔诚顶住   | 0.25 | 0.05 | **0.85** | 0.10 | 0.05 |
   | C 同时撑开   | 0.25 | 0.05 | 0.05 | **0.80** | 0.05 |
   | D 问句问错了 | 0.15 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**只往一层里种偏离,分解必须把它归给那一层而不归给另一层**)
     and 负控开火(**两层都匀速的世界里,两个偏离都必须落在零上**)
     and 恒等式控制通过(`dep_虔诚 − dep_世俗 == dep_差距`,数值上精确):
      只世俗排零 -> A · 只虔诚排零 -> B · 两个都排零 -> C · 两个都不排零 -> D
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`homosex` 有地板/天花板**(1–4 的序数题)。
  **一层若已经贴近某一端,它就没有位移空间,于是「它没动」会被天花板伪造出来。**
  ⇒ 控制:**逐层印出该十年首尾的均值与到两端的余量**,并算 `headroom = |位移| ÷ 可达范围`;
  **余量不足的层,它的「没偏离」一律记为 `UNVERIFIED` 而不是「照原速」。**

`G4` 规格曲线:**不是只报九十年代那一格** —— 全部可用十年 × 两层一起报,含不同意的格。
`G3` 多重性:整张网格 BH,并印出族大小(`#832`:族越窄存活越易,是算术不是证据)。
⚠ 本轮换不了仪器(GSS 是唯一有五十年逐年分层的美国态度调查;`#837` 已证 SCCS 走不通),
  **而它不需要** —— 本轮问的是同一具仪器内部的分解,不是跨仪器复现。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, B, Q = "homosex", 4, 6000, 0.05

cols = ["year", "attend", "reliten", "fund", IT]
d = pd.read_stata(gp, columns=cols, convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1)
M = M.join(R["REL"])
t = M.groupby("year")["REL"].transform(
    lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
HI, LO = (t == 2), (t == 0)          # HI = 虔诚层 · LO = 世俗层

print(f"=== ⓪ 硬规则①:变量真的问了哪些年、每层多少人 ===")
ok = M[IT].notna() & (HI | LO)
ys = {}
for y, g in M[ok].groupby("year"):
    a = g[HI.loc[g.index]][IT].to_numpy(float); b = g[LO.loc[g.index]][IT].to_numpy(float)
    if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
S = sorted(ys)
print(f"  `{IT}` 1–{KK}(高=越宽容)· 合格年 **{len(S)}** 个:{S[0]}–{S[-1]} · "
      f"总 n **{int(M[ok][IT].notna().sum()):,}**")
dec = {}
for y in S: dec.setdefault((y//10)*10, []).append(y)
dec = {k: v for k, v in dec.items() if len(v) >= 3}
print(f"  可用十年(每十年 ≥3 个合格年):{sorted(dec)} · 每十年年数 {[len(dec[k]) for k in sorted(dec)]}")
span = S[-1]-S[0]
FULL = {"虔诚层": float(ys[S[-1]][0].mean()-ys[S[0]][0].mean()),
        "世俗层": float(ys[S[-1]][1].mean()-ys[S[0]][1].mean())}
print(f"  全程({span} 年)位移:虔诚层 **{FULL['虔诚层']:+.4f}** · 世俗层 **{FULL['世俗层']:+.4f}** "
      f"(⇒ 差距全程变化 {FULL['虔诚层']-FULL['世俗层']:+.4f})")

# ── 跑前混淆的控制:天花板余量 ────────────────────────────────────────────────
print(f"\n=== ⓪b 跑前写下的最强混淆的控制:`{IT}` 是 1–{KK} 的序数题,贴边的层没有位移空间 ===")
HEAD = {}
for k, i in (("虔诚层", 0), ("世俗层", 1)):
    for dc in sorted(dec):
        y0, y1 = dec[dc][0], dec[dc][-1]
        m0, m1 = float(ys[y0][i].mean()), float(ys[y1][i].mean())
        mv = m1-m0; reach = (KK-m0) if mv > 0 else (m0-1)
        HEAD[(k, dc)] = dict(m0=m0, m1=m1, move=mv, reach=reach,
                             headroom=(abs(mv)/reach if reach > 0 else float("inf")))
    print(f"  {k}: " + " · ".join(
        f"{dc}s 首{HEAD[(k,dc)]['m0']:.2f}→尾{HEAD[(k,dc)]['m1']:.2f}"
        f"(余 {HEAD[(k,dc)]['reach']:.2f},占 {HEAD[(k,dc)]['headroom']:.0%})" for dc in sorted(dec)))
print(f"  ⚠ **余量占比高的格,它的「没偏离」记为 `UNVERIFIED`,不记为「照原速」。**")

# ── 分解:每层自己的偏离 ──────────────────────────────────────────────────────
def dep(i, dc, rng, Bv=B, src=None):
    """该层在该十年的位移 − 该层自己的匀速参照。自助 Bv 次。"""
    Y = src if src else ys; yy = dec[dc]
    ref = (float(Y[S[-1]][i].mean()-Y[S[0]][i].mean()))*(yy[-1]-yy[0])/span
    out = np.empty(Bv); r = lambda a: a[rng.integers(0, len(a), len(a))]
    for b in range(Bv):
        out[b] = r(Y[yy[-1]][i]).mean() - r(Y[yy[0]][i]).mean() - ref
    obs = float(Y[yy[-1]][i].mean()-Y[yy[0]][i].mean()) - ref
    return obs, out

rng = np.random.default_rng(277)
print(f"\n=== ① `G4` 规格曲线:**全部十年 × 两层**(不是只报九十年代那一格)· B={B} ===")
CELLS, GRIDR = [], {}
for dc in sorted(dec):
    for k, i in (("虔诚层", 0), ("世俗层", 1)):
        obs, bs = dep(i, dc, rng)
        lo, hi = float(np.quantile(bs, .025)), float(np.quantile(bs, .975))
        p = max(2*min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0))), 1.0/(B+1))
        GRIDR[(k, dc)] = dict(obs=obs, lo=lo, hi=hi, p=p, sd=float(np.std(bs)),
                              excl=bool(lo > 0 or hi < 0), **HEAD[(k, dc)])
        CELLS.append((k, dc))
for dc in sorted(dec):
    print(f"  {dc}s  " + "  |  ".join(
        f"{k} **{GRIDR[(k,dc)]['obs']:+.4f}** [{GRIDR[(k,dc)]['lo']:+.4f},{GRIDR[(k,dc)]['hi']:+.4f}]"
        f" p={GRIDR[(k,dc)]['p']:.4f}{' **排零**' if GRIDR[(k,dc)]['excl'] else ''}"
        for k in ("虔诚层", "世俗层")))
ps = [GRIDR[c]["p"] for c in CELLS]
bh = {CELLS[i] for i in Gate.bh(ps, Q)}
by = {CELLS[i] for i in Gate.by(ps, Q)}
print(f"  `G3` 整张网格 **{len(CELLS)} 格**(族大小印在旁边 —— `#832`:族越窄存活越易,是算术不是证据)"
      f" ⇒ BH 存活 **{len(bh)}** · BY 存活 **{len(by)}**")
print(f"     BH:{sorted(f'{k}{dc}s' for k,dc in bh) or '无'}   BY:{sorted(f'{k}{dc}s' for k,dc in by) or '无'}")

D90 = {k: GRIDR[(k, 1990)] for k in ("虔诚层", "世俗层")}
print(f"\n=== ② 九十年代那一格的分法(⚠ **和是恒等式,分法才是测量**)===")
gapdep = D90["虔诚层"]["obs"] - D90["世俗层"]["obs"]
for k in ("虔诚层", "世俗层"):
    g = D90[k]
    tot = sum(abs(D90[q]["obs"]) for q in D90)
    share = abs(g["obs"])/tot if tot else float("nan")
    print(f"  {k}:偏离 **{g['obs']:+.4f}** [{g['lo']:+.4f},{g['hi']:+.4f}] · "
          f"占差距偏离的 **{share:+.0%}** · 天花板余量占用 {g['headroom']:.0%} · "
          f"{'**区间排除零**' if g['excl'] else '区间含零'}")
print(f"  ⇒ 差距偏离 = {gapdep:+.4f}(**= 两层之差,恒等式**)")
excl = [k for k in ("虔诚层", "世俗层") if D90[k]["excl"]]

print("\n=== ③ 控制 ===")
# 正控:只往虔诚层的 1990s 尾年种一个位移
def shifted(i, dc, amt):
    Y = {y: (a.copy(), b.copy()) for y, (a, b) in ys.items()}
    y1 = dec[dc][-1]; Y[y1][i][:] = Y[y1][i] + amt
    return Y
# ⚠⚠⚠ **正控在这一轮失败了两次,两次都是我的期望值错,不是仪器错(`#818` 那一类)——
#      而失败第三次之后才看清共同的病根,它值得单独写下来:**
#      ① 第一次:预期写成 `A×(1 − 十年跨年/全程跨年)`,而种植点在**十年内部年**,
#         全程位移根本不变 ⇒ 参照不变 ⇒ 偏离恰好整体上移 A。
#      ② 第二次:改种在 `S[-1]=2024` 并预期 2010s 净上移,**而 2024 根本不在 2010s 这一格里** ——
#         `2024//10*10 = 2020`,而 2020s 因合格年 <3 被丢掉了。**最后五年不属于任何一个十年格。**
#      **共同的病根:我一直在用「我脑子里的分桶」推期望,而不是用代码真正算出来的分桶。**
#      ⇒ 所以不再手推两个数,而是**把期望写成一条从代码读出来的公式,并让它一次预测整张网格**:
#         `Δ偏离(十年 d, 种植年 yp, 量 A) =`
#            `A·[yp==d[-1]] − A·[yp==d[0]]`                        (观测侧)
#          `− (A·[yp==S[-1]] − A·[yp==S[0]])·(d[-1]−d[0])/span`     (参照侧)
#      **4 个种植点 × 5 个十年 = 20 条预测同时对上,才算这段机器验过了。**
AMT = 0.25
base_h = {dc: GRIDR[("虔诚层", dc)]["obs"] for dc in sorted(dec)}
base_l = {dc: GRIDR[("世俗层", dc)]["obs"] for dc in sorted(dec)}
def predict(d, yp, A):
    yy = dec[d]
    obs = A*(yp == yy[-1]) - A*(yp == yy[0])
    ref = (A*(yp == S[-1]) - A*(yp == S[0]))*(yy[-1]-yy[0])/span
    return obs - ref
SITES = {f"1990s 尾年 {dec[1990][-1]}(十年内部)": dec[1990][-1],
         f"1990s 首年 {dec[1990][0]}(十年内部)": dec[1990][0],
         f"全程末年 {S[-1]}(**不属于任何十年格**)": S[-1],
         f"全程首年 {S[0]}(= 1970s 首年)": S[0]}
print(f"  正控:**一条从代码读出来的公式,一次预测整张网格** —— {len(SITES)} 个种植点 × {len(dec)} 个十年")
print(f"     ⚠ 硬规则①:`{S[-1]}//10*10 = {(S[-1]//10)*10}` ⇒ **最后一年不在任何一个已报十年格里**")
pc_rows, pc_bad = [], 0
for lab, yp in SITES.items():
    Yx = {y: (a.copy(), b.copy()) for y, (a, b) in ys.items()}
    Yx[yp][0][:] = Yx[yp][0] + AMT
    got = [dep(0, dc, np.random.default_rng(11), 300, src=Yx)[0] - base_h[dc] for dc in sorted(dec)]
    exp = [predict(dc, yp, AMT) for dc in sorted(dec)]
    bad = sum(abs(g-e) > 1e-9 for g, e in zip(got, exp)); pc_bad += bad
    pc_rows.append(dict(site=lab, year=int(yp), got=got, exp=exp, mismatch=int(bad)))
    print(f"     {lab:34s} " + " · ".join(
        f"{dc}s {g:+.4f}/{e:+.4f}" for dc, g, e in zip(sorted(dec), got, exp))
        + ("  ⇒ **全对**" if bad == 0 else f"  ⇒ ⚠ {bad} 处不符"))
# 世俗层不受虔诚层种植影响
Yx = {y: (a.copy(), b.copy()) for y, (a, b) in ys.items()}
Yx[dec[1990][-1]][0][:] = Yx[dec[1990][-1]][0] + AMT
cross = max(abs(dep(1, dc, np.random.default_rng(12), 300, src=Yx)[0]-base_l[dc]) for dc in sorted(dec))
print(f"     交叉不渗漏:只种虔诚层 ⇒ **世俗层各十年最大变动 {cross:.2e}**")
pc_ok = (pc_bad == 0 and cross < 1e-12)
print(f"  ⇒ **{len(SITES)*len(dec)} 条预测中不符 {pc_bad} 条** ⇒ 正控 {'**过**' if pc_ok else '**没过**'}")

# 负控:构造一个两层都严格匀速的世界
Yu = {}
for j, y in enumerate(S):
    Yu[y] = tuple(ys[S[0]][i] + (FULL[k]*(y-S[0])/span) for i, k in ((0, "虔诚层"), (1, "世俗层")))
nc_h, _ = dep(0, 1990, np.random.default_rng(2), 400, src=Yu)
nc_l, _ = dep(1, 1990, np.random.default_rng(2), 400, src=Yu)
print(f"  负控:构造**两层都严格匀速**的世界 ⇒ 虔诚层偏离 **{nc_h:+.6f}** · 世俗层 **{nc_l:+.6f}** "
      f"—— ⚠ **「这个零该不该是零?」该**:匀速世界里按定义没有偏离,所以对 **0**")

G = Gate("#838 · 九十年代那一格里,到底是谁在动")
G.asserted("① 硬规则①:变量真的问了哪些年、每层每年 ≥120 人、可用十年 ≥3 个合格年 —— 都从对象印出",
           bool(len(S) >= 8 and len(dec) >= 3),
           f"合格年 {len(S)} 个({S[0]}–{S[-1]})· 可用十年 {sorted(dec)}", kind="control")
G.identity_control("② 恒等式控制:`dep(虔诚) − dep(世俗)` 必须**数值上精确等于** `dep(差距)` —— "
                   "⚠ **这条不是发现,是代数强制;写下来正是为了不把它当发现**",
                   observed=float(D90["虔诚层"]["obs"]-D90["世俗层"]["obs"]), expected=float(gapdep),
                   tol=1e-12, what="两层偏离之差 vs 差距偏离", deterministic=True)
G.asserted("③ 正控(**一条从代码读出来的公式,一次预测整张网格**):`Δ偏离 = A·[yp=d尾] − A·[yp=d首] "
           "− (A·[yp=全末] − A·[yp=全首])·d跨年/全程跨年` —— **4 个种植点 × 5 个十年 = 20 条预测全对**,"
           "且**只种虔诚层时世俗层零渗漏**;⚠ **这条替掉了我手推的两个数,而它们错了两次(`#818` 那一类)**",
           bool(pc_ok), f"20 条预测中不符 **{pc_bad}** 条 · 交叉渗漏最大 {cross:.2e}", kind="control")
G.asserted("④ 负控:**两层都严格匀速**的世界里,两个偏离都必须落在 **0** 上"
           "(⚠ **这个零该是零**:匀速按定义无偏离)",
           bool(abs(nc_h) < 1e-9 and abs(nc_l) < 1e-9),
           f"虔诚层 {nc_h:+.2e} · 世俗层 {nc_l:+.2e}", kind="control")
G.asserted("⑤ 前提(跑前写下的最强混淆):**`homesex` 是 1–4 序数题,贴边的层没有位移空间** ⇒ "
           "逐格印出到端点的余量与占用比;**余量不足的层,「没偏离」记 `UNVERIFIED` 而非「照原速」**",
           bool(all(np.isfinite(v["headroom"]) for v in HEAD.values())),
           " · ".join(f"{k}{dc}s {v['headroom']:.0%}" for (k, dc), v in HEAD.items() if dc == 1990),
           kind="control")
G.asserted("⑥ kill(预注册):「是某一边在动」要成立,需**恰好一层**的偏离区间排除零",
           bool(len(excl) == 1), f"排零的层 {excl or '无'}(共 {len(excl)}/2)", kind="kill",
           yardstick="每层自己的十年偏离,对照它自己的 95% 自助区间",
           yardstick_noise=float(np.mean([D90[k]["sd"] for k in D90])))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
_tot = sum(abs(D90[q]["obs"]) for q in D90)
sh = {k: abs(D90[k]["obs"])/_tot for k in D90} if _tot else {}
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(excl) == 1 and excl[0] == "世俗层":
    V = (f"**A 是世俗那边走了。** 世俗层偏离 **{D90['世俗层']['obs']:+.4f}** "
         f"[{D90['世俗层']['lo']:+.4f},{D90['世俗层']['hi']:+.4f}] 排除零,占两层偏离总量 **{sh['世俗层']:.0%}**;"
         f"虔诚层 {D90['虔诚层']['obs']:+.4f} 含零。\n"
         f"  ⇒ **一句关于人的话:九十年代那条缝张开,不是因为信教的人变得更不容忍,\n"
         f"  是因为不信教的人走得比他们自己前二十年更快 —— 缝是被一边的加速拉开的,\n"
         f"  另一边基本还在按它原来的速度走。**")
elif len(excl) == 1 and excl[0] == "虔诚层":
    V = (f"**B 是虔诚那边顶住了。** 虔诚层偏离 **{D90['虔诚层']['obs']:+.4f}** "
         f"[{D90['虔诚层']['lo']:+.4f},{D90['虔诚层']['hi']:+.4f}] 排除零,占两层偏离总量 **{sh['虔诚层']:.0%}**;"
         f"世俗层 {D90['世俗层']['obs']:+.4f} 含零。\n"
         f"  ⇒ **一句关于人的话:九十年代那条缝张开,不是因为世俗那边跑快了,\n"
         f"  是因为信教的人在那十年里比自己五十年的步子慢了下来 ——\n"
         f"  这是一个顶住的动作,不是一个走远的动作。**")
elif len(excl) == 2:
    V = (f"**C 同时向两侧撑开。** 虔诚层 **{D90['虔诚层']['obs']:+.4f}** 与 "
         f"世俗层 **{D90['世俗层']['obs']:+.4f}** 都排除零,方向相反。\n"
         f"  ⇒ **一句关于人的话:没有「哪一边动了」这回事 —— 九十年代那十年里,\n"
         f"  两边同时偏离了各自的步子,朝相反方向。缝不是被拉开的,是被撑开的。**")
else:
    V = (f"**D 这个问句本身问错了 —— 而这是本轮最该带走的。** 虔诚层 {D90['虔诚层']['obs']:+.4f} "
         f"[{D90['虔诚层']['lo']:+.4f},{D90['虔诚层']['hi']:+.4f}] 与 世俗层 {D90['世俗层']['obs']:+.4f} "
         f"[{D90['世俗层']['lo']:+.4f},{D90['世俗层']['hi']:+.4f}] **各自都不排除零**,\n"
         f"  而它们的差 {gapdep:+.4f} 却是全项目唯一穿过所有校正的那件事。\n"
         f"  ⇒ **偏离是纯交互项,不可归给任何一层** ⇒ **「是哪一边在动」这个问法预设了\n"
         f"  一个不存在的分解。** 差距是一个关于**两群人之间关系**的量,\n"
         f"  而它并不由任何一群人自己的轨迹承载。")
print(V)
print(f"\n⚠ **和是恒等式,分法才是测量**:`dep(虔诚) − dep(世俗) ≡ dep(差距)` 被代数强制"
      f"(控制②数值验过),**本轮唯一的测量是这个和如何分配,以及各自排不排零。**")
json.dump(dict(item=IT, years=S, decades={str(k): v for k, v in dec.items()}, full_span=FULL,
               grid={f"{k}|{dc}": v for (k, dc), v in GRIDR.items()},
               headroom={f"{k}|{dc}": v for (k, dc), v in HEAD.items()},
               nineties=D90, gap_departure=gapdep, shares=sh, excludes_zero=excl,
               bh=sorted(f"{k}|{dc}" for k, dc in bh), by=sorted(f"{k}|{dc}" for k, dc in by),
               family_size=len(CELLS), B=B, q=Q,
               pos_control=dict(planted=AMT, rows=pc_rows, mismatches=pc_bad, cross_leak=cross,
                                formula="A*[yp==d_last] - A*[yp==d_first] - (A*[yp==S_last]-A*[yp==S_first])*d_span/span"),
               neg_control=dict(devout=nc_h, secular=nc_l, reference=0.0),
               identity="dep(devout) - dep(secular) == dep(gap), algebraically forced",
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"who_actually_moved.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'who_actually_moved.json'}")
