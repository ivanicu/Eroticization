"""#839 · E03·A79·R278 —— 二〇一〇年代虔诚层那次「追上」,是真的,还是量表到头了?

`#838`① 预注册了这一轮,而**在写它之前必须先说清它同时是在攻击 `#838` 自己**(`§3` basin:
上一轮的结论我很满意 ⇒ 这一轮专门去找它站不住的地方)。

**`#838` 量到的两件事,合起来是一个很好听的故事:**
九十年代**世俗层加速**(+0.6266,排零)拉开了缝;二〇一〇年代**虔诚层加速**(+0.3331,排零)
—— 缝开始合上。**而「好听」正是要去查的理由。**

**⚠⚠ 威胁在 `#838` 自己印出来的余量表里,而我上一轮没有追下去:**
**二〇一〇年代,世俗层已经走到 3.08 → 3.49,而量表上端只有 4 —— 只剩 0.92 的余量,用掉 44%。**
**一个贴着天花板的层,它的「慢下来」可能根本不是它慢了,是尺子到头了。**
若如此,「虔诚层追上来」就有一部分是**世俗层被截断**造出来的假象。

G1 估计量:**同一套分解,换成一个没有天花板的估计量上重跑** ——
   把 1–4 的有序作答看成一个**潜在连续量**被三个阈值切开(有序 probit,阈值全样本共用、方差固定),
   **逐年逐层估一个潜在均值 μ。μ 不受 1–4 的边界约束,所以天花板无法伪造出「慢下来」。**
   然后在 μ 上重算 `#838` 的每一格偏离:`dep_s(d) = Δμ_s(d) − Δμ_s(全程)×d跨年÷全程跨年`。

三个世界:
   A **真的收敛**:潜在尺上,二〇一〇年代虔诚层仍偏离、世俗层的放缓仍在 ⇒ **天花板不是解释。**
   B **天花板伪造**:换到潜在尺后,二〇一〇年代那格**消失或翻号** ⇒ **`#838` 的第二半要撤。**
   C **换了估计量,九十年代那格也变了** ⇒ **⚠ 那就不只是二〇一〇年代的问题,
     `#838` 的主结论本身是估计量依赖的**,必须整条降级。**这是我最不欢迎的那个,所以要写在前面。**

预测矩阵:
   | 世界 | 现在 | 2010s 仍在 & 1990s 仍在 | 2010s 没了 | 1990s 也变了 |
   | A 真收敛       | 0.45 | **0.85** | 0.05 | 0.05 |
   | B 天花板伪造   | 0.35 | 0.10 | **0.85** | 0.10 |
   | C 估计量依赖   | 0.20 | 0.05 | 0.10 | **0.85** |

预注册判词(条件式):
  if 正控开火(**在潜在尺上种一个已知位移,估计量必须按预期取回**)
     and 负控开火(**潜在尺上严格匀速的世界,偏离必须为 0**)
     and 阈值控制通过(**三个阈值严格递增,且拟合出的类别份额与观测份额吻合**):
      1990s 世俗层仍排零 且 2010s 虔诚层仍排零 -> A
      1990s 仍排零 而 2010s 虔诚层不再排零     -> B(**撤 `#838` 的第二半**)
      1990s 也不再排零                         -> C(**整条降级为估计量依赖**)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**有序 probit 的阈值若逐年重估,μ 与阈值就不可分离**
  (同一份作答分布可以由「人更宽容」或「阈值移动」解释,二者共线)。
  ⇒ 控制:**阈值由全样本汇总一次估定并固定,方差固定为 1;逐年逐层只估 μ 一个参数。**
  **这是一个必须写下来的假设,不是一个被检验的结论:它假设「非常不对 / 有点不对 / …」
  这几个词的含义五十年没变。⚠ 而那几乎肯定不完全对 —— 所以本轮判的是
  「在词义不变这个假设下,天花板还解不解释得掉」,不是「天花板解释不掉」。**

⚠ 本轮换不了仪器(同一份 GSS),而它**不需要** —— 本轮换的是**估计量**,
  这正是 `G4` 要求的规格曲线里的一根轴(estimator-robust),不是跨仪器复现。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, B, Q = "homosex", 4, 1500, 0.05

d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
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
HI, LO = (t == 2), (t == 0)
ok = M[IT].notna() & (HI | LO)
ys = {}
for y, g in M[ok].groupby("year"):
    a = g[HI.loc[g.index]][IT].to_numpy(float); b = g[LO.loc[g.index]][IT].to_numpy(float)
    if len(a) >= 120 and len(b) >= 120: ys[int(y)] = (a, b)
S = sorted(ys); span = S[-1]-S[0]
dec = {}
for y in S: dec.setdefault((y//10)*10, []).append(y)
dec = {k: v for k, v in dec.items() if len(v) >= 3}

# ── 有序 probit:阈值全样本估一次并固定,逐格只估 μ ─────────────────────────────
Phi = lambda z: 0.5*(1.0+np.vectorize(math.erf)(z/math.sqrt(2.0)))
def ppf(p):
    lo_, hi_ = -8.0, 8.0
    for _ in range(200):
        mid = (lo_+hi_)/2
        if Phi(np.array([mid]))[0] < p: lo_ = mid
        else: hi_ = mid
    return (lo_+hi_)/2
allv = np.concatenate([np.concatenate(ys[y]) for y in S])
cum = np.cumsum([np.mean(allv == k) for k in range(1, KK+1)])[:-1]
TAU = np.array([ppf(c) for c in cum])
GRIDMU = np.linspace(-4.0, 4.0, 1601)
edges = np.concatenate([[-np.inf], TAU, [np.inf]])
LOGP = np.log(np.clip(np.array([[Phi(np.array([edges[k+1]-m]))[0] - Phi(np.array([edges[k]-m]))[0]
                                 for k in range(KK)] for m in GRIDMU]), 1e-12, 1))
def mu(v):
    c = np.array([(v == k).sum() for k in range(1, KK+1)], float)
    return float(GRIDMU[np.argmax(LOGP @ c)])

print("=== ⓪ 硬规则① + 阈值控制:阈值从全样本估一次并固定,逐格只估 μ 一个参数 ===")
print(f"  合格年 {len(S)} 个 {S[0]}–{S[-1]} · 可用十年 {sorted(dec)} · 总作答 {len(allv):,}")
print(f"  全样本各类份额 {[round(float(np.mean(allv==k)),4) for k in range(1,KK+1)]} ⇒ "
      f"阈值 τ = {np.round(TAU,4).tolist()}(严格递增 {bool(np.all(np.diff(TAU)>0))})")
fit = np.exp(LOGP[np.argmin(np.abs(GRIDMU-0.0))])
print(f"  拟合复原(μ=0 处)份额 {np.round(fit,4).tolist()} ⇒ 最大偏差 "
      f"**{float(np.max(np.abs(fit-np.array([np.mean(allv==k) for k in range(1,KK+1)])))):.2e}**")
print(f"  ⚠ **这是一个假设不是结论**:它假设那几个选项词的含义五十年没变。"
      f"⇒ 本轮判的是「在词义不变的假设下,天花板还解不解释得掉」。")

def dep_lat(i, dc, rng, Bv=B, src=None):
    Y = src if src else ys; yy = dec[dc]
    ref = (mu(Y[S[-1]][i])-mu(Y[S[0]][i]))*(yy[-1]-yy[0])/span
    obs = mu(Y[yy[-1]][i]) - mu(Y[yy[0]][i]) - ref
    out = np.empty(Bv); r = lambda a: a[rng.integers(0, len(a), len(a))]
    for b in range(Bv):
        rf = (mu(r(Y[S[-1]][i]))-mu(r(Y[S[0]][i])))*(yy[-1]-yy[0])/span
        out[b] = mu(r(Y[yy[-1]][i])) - mu(r(Y[yy[0]][i])) - rf
    return obs, out

rng = np.random.default_rng(278)
print(f"\n=== ① 潜在尺上重跑 `#838` 的整张网格(B={B})· ⚠ 对照列是 `#838` 的 1–4 均值尺 ===")
OLD = {("虔诚层",1970):-0.0214,("世俗层",1970):+0.0102,("虔诚层",1980):-0.1039,("世俗层",1980):-0.2859,
       ("虔诚层",1990):+0.1137,("世俗层",1990):+0.6266,("虔诚层",2000):-0.0885,("世俗层",2000):+0.2716,
       ("虔诚层",2010):+0.3331,("世俗层",2010):+0.1947}
NEW, CELLS = {}, []
for dc in sorted(dec):
    for k, i in (("虔诚层", 0), ("世俗层", 1)):
        obs, bs = dep_lat(i, dc, rng)
        lo_, hi_ = float(np.quantile(bs, .025)), float(np.quantile(bs, .975))
        p = max(2*min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0))), 1.0/(B+1))
        NEW[(k, dc)] = dict(obs=obs, lo=lo_, hi=hi_, p=p, sd=float(np.std(bs)),
                            excl=bool(lo_ > 0 or hi_ < 0), old=OLD[(k, dc)])
        CELLS.append((k, dc))
for dc in sorted(dec):
    print(f"  {dc}s  " + "  |  ".join(
        f"{k} 潜在 **{NEW[(k,dc)]['obs']:+.4f}** [{NEW[(k,dc)]['lo']:+.4f},{NEW[(k,dc)]['hi']:+.4f}]"
        f"{' **排零**' if NEW[(k,dc)]['excl'] else '     '} (均值尺 {OLD[(k,dc)]:+.4f})"
        for k in ("虔诚层", "世俗层")))
ps = [NEW[c]["p"] for c in CELLS]
bh = {CELLS[i] for i in Gate.bh(ps, Q)}; by = {CELLS[i] for i in Gate.by(ps, Q)}
print(f"  `G3` 整张网格 **{len(CELLS)} 格** ⇒ BH 存活 **{len(bh)}** · BY 存活 **{len(by)}**")
print(f"     BH:{sorted(f'{k}{dc}s' for k,dc in bh) or '无'}")
agree = sum(1 for c in CELLS if NEW[c]["excl"] == (abs(OLD[c]) > 0 and c in
            {("世俗层",1980),("世俗层",1990),("世俗层",2000),("世俗层",2010),("虔诚层",2010)}))
print(f"  ⇒ **与 `#838` 均值尺结论一致的格:{agree}/{len(CELLS)}**")

k90 = NEW[("世俗层", 1990)]["excl"]; k10 = NEW[("虔诚层", 2010)]["excl"]
sec10 = NEW[("世俗层", 2010)]
print(f"\n=== ② 两个关键格 ===")
print(f"  1990s 世俗层:潜在 **{NEW[('世俗层',1990)]['obs']:+.4f}** "
      f"[{NEW[('世俗层',1990)]['lo']:+.4f},{NEW[('世俗层',1990)]['hi']:+.4f}] ⇒ "
      f"{'**仍排零**' if k90 else '**不再排零**'}")
print(f"  2010s 虔诚层:潜在 **{NEW[('虔诚层',2010)]['obs']:+.4f}** "
      f"[{NEW[('虔诚层',2010)]['lo']:+.4f},{NEW[('虔诚层',2010)]['hi']:+.4f}] ⇒ "
      f"{'**仍排零**' if k10 else '**不再排零**'}")
print(f"  2010s 世俗层(受天花板威胁的那一格):潜在 **{sec10['obs']:+.4f}** "
      f"[{sec10['lo']:+.4f},{sec10['hi']:+.4f}](均值尺 {sec10['old']:+.4f})")

print("\n=== ③ 控制 ===")
# 正控:在潜在尺上把某年整体推高(用阈值反解构造作答)
def latent_shift(v, amt):
    m = mu(v); pr = np.exp(LOGP[np.argmin(np.abs(GRIDMU-(m+amt)))])
    rg2 = np.random.default_rng(9)
    return rg2.choice(np.arange(1, KK+1), size=len(v), p=pr/pr.sum()).astype(float)
AMT = 0.30
Yp = {y: (a.copy(), b.copy()) for y, (a, b) in ys.items()}
Yp[dec[2010][-1]] = (latent_shift(ys[dec[2010][-1]][0], AMT), Yp[dec[2010][-1]][1])
pc, _ = dep_lat(0, 2010, np.random.default_rng(5), 200, src=Yp)
pc_l, _ = dep_lat(1, 2010, np.random.default_rng(5), 200, src=Yp)
base10 = NEW[("虔诚层", 2010)]["obs"]; base10l = NEW[("世俗层", 2010)]["obs"]
print(f"  正控:在**潜在尺**上把虔诚层 {dec[2010][-1]} 年整体推高 +{AMT} ⇒ "
      f"偏离 {base10:+.4f} → **{pc:+.4f}**(动 **{pc-base10:+.4f}**,预期 ≈ +{AMT});"
      f"世俗层动 {pc_l-base10l:+.4f}")
# 负控:潜在尺上严格匀速的世界
Yu = {}
for y in S:
    Yu[y] = tuple(latent_shift(ys[S[0]][i], (mu(ys[S[-1]][i])-mu(ys[S[0]][i]))*(y-S[0])/span)
                  for i in (0, 1))
nc, _ = dep_lat(0, 2010, np.random.default_rng(6), 200, src=Yu)
nc_l, _ = dep_lat(1, 1990, np.random.default_rng(6), 200, src=Yu)
print(f"  负控:构造**潜在尺上严格匀速**的世界 ⇒ 虔诚层 2010s **{nc:+.4f}** · 世俗层 1990s **{nc_l:+.4f}** "
      f"—— ⚠ **「这个零该不该是零?」该**(匀速按定义无偏离);"
      f"⚠ 但它经过一次**重抽样离散化**,所以有抽样噪声,不是解析零")

G = Gate("#839 · 二〇一〇年代那次追上,是真的还是量表到头了")
G.asserted("① 阈值控制:阈值由**全样本估一次并固定**,方差固定为 1,逐格只估 μ 一个参数 ⇒ "
           "μ 与阈值不共线;且阈值严格递增、拟合份额复原观测份额 "
           "（⚠ **这是写下来的假设不是被检验的结论**:它假设选项词义五十年不变）",
           bool(np.all(np.diff(TAU) > 0) and
                np.max(np.abs(fit-np.array([np.mean(allv == k) for k in range(1, KK+1)]))) < 1e-3),
           f"τ={np.round(TAU,4).tolist()} · 份额最大偏差 "
           f"{float(np.max(np.abs(fit-np.array([np.mean(allv==k) for k in range(1,KK+1)])))):.2e}",
           kind="control")
G.asserted("② 正控:在**潜在尺**上推高 +0.30,估计量必须按预期取回,且**不渗漏到另一层**",
           bool(abs((pc-base10)-AMT) < 0.06 and abs(pc_l-base10l) < 0.02),
           f"虔诚层动 {pc-base10:+.4f}(预期 +{AMT})· 世俗层动 {pc_l-base10l:+.4f}", kind="control")
G.asserted("③ 负控:**潜在尺上严格匀速**的世界里偏离必须落在 0 上"
           "(⚠ **这个零该是零**;⚠ 但经过一次重抽样离散化 ⇒ 允许抽样噪声,不是解析零)",
           bool(abs(nc) < 0.10 and abs(nc_l) < 0.10),
           f"虔诚层 2010s {nc:+.4f} · 世俗层 1990s {nc_l:+.4f}", kind="control")
G.asserted("④ 前提(跑前写下的最强混淆):**世俗层 2010s 贴天花板(3.08→3.49,余 0.92,用掉 44%)** ⇒ "
           "换到**不受 1–4 边界约束**的潜在尺上重跑整张网格,而不是只重跑那一格",
           bool(len(CELLS) == 10), f"整张 {len(CELLS)} 格均在潜在尺上重跑", kind="control")
G.asserted("⑤ kill(预注册):「`#838` 的结论不是天花板伪造的」要成立,需 **1990s 世俗层与 2010s 虔诚层"
           "在潜在尺上都仍排除零**",
           bool(k90 and k10), f"1990s 世俗 {'排零' if k90 else '不排零'} · "
           f"2010s 虔诚 {'排零' if k10 else '不排零'}", kind="kill",
           yardstick="潜在尺上每格自己的偏离,对照它自己的 95% 自助区间",
           yardstick_noise=float(np.mean([NEW[c]["sd"] for c in CELLS])))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif k90 and k10:
    V = (f"**A 真的收敛 —— 天花板解释不掉它。** 换到不受 1–4 边界约束的潜在尺上,"
         f"1990s 世俗层 **{NEW[('世俗层',1990)]['obs']:+.4f}** 与 2010s 虔诚层 "
         f"**{NEW[('虔诚层',2010)]['obs']:+.4f}** 都仍排除零。\n"
         f"  ⇒ **一句关于人的话:那条缝先被不信教的人拉开,又被信教的人追上来 ——\n"
         f"  而追上来那一段不是尺子到头造出来的假象:换一把没有上端的尺子,它还在。**")
elif k90 and not k10:
    V = (f"**B 天花板伪造 ⇒ 撤 `#838` 的第二半。** 潜在尺上 1990s 世俗层仍排零 "
         f"({NEW[('世俗层',1990)]['obs']:+.4f}),而 **2010s 虔诚层不再排零** "
         f"({NEW[('虔诚层',2010)]['obs']:+.4f} [{NEW[('虔诚层',2010)]['lo']:+.4f},"
         f"{NEW[('虔诚层',2010)]['hi']:+.4f}])。\n"
         f"  ⇒ **「信教的那一边在二〇一〇年代追上来」这句话要撤:换一把没有上端的尺子,\n"
         f"  它就不在了 —— 它是世俗那边贴着天花板、显得慢下来造出来的。**")
else:
    V = (f"**C 估计量依赖 ⇒ 整条降级,而这是我最不欢迎的那个。** 潜在尺上 **1990s 世俗层也不再排零** "
         f"({NEW[('世俗层',1990)]['obs']:+.4f} [{NEW[('世俗层',1990)]['lo']:+.4f},"
         f"{NEW[('世俗层',1990)]['hi']:+.4f}])。\n"
         f"  ⇒ **`#838` 的主结论是估计量依赖的:它在 1–4 均值尺上成立,在潜在尺上不成立。\n"
         f"  ⇒ 那不是一条关于人的结论,是一条关于「用哪把尺子量」的结论。**")
print(V)
json.dump(dict(item=IT, tau=TAU.tolist(), grid={f"{k}|{dc}": v for (k, dc), v in NEW.items()},
               mean_scale_reference=({f"{k}|{dc}": v for (k, dc), v in OLD.items()}),
               bh=sorted(f"{k}|{dc}" for k, dc in bh), by=sorted(f"{k}|{dc}" for k, dc in by),
               key=dict(sec_1990_excl=k90, dev_2010_excl=k10),
               pos_control=dict(planted=AMT, devout=pc-base10, secular=pc_l-base10l),
               neg_control=dict(devout_2010=nc, secular_1990=nc_l),
               assumption="thresholds fixed across 50 years = option wording assumed stable",
               B=B, q=Q, admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"latent_or_ceiling.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'latent_or_ceiling.json'}")
