"""#840 · E03·A80·R279 —— 「缝合上了」是关于差距的主张,而差距的区间从没算过

`#839`① 预注册了这一轮,理由是 `#770` 早就立过的一条:
**要比的是那个量本身,不是两个量各自与零的关系。**

`#838`/`#839` 报的是**两层各自**的偏离与各自的区间,然后我读出了「先拉开、后追上」。
**而「缝拉开了」「缝合上了」是关于差距的主张** —— `dep(差距) = dep(虔诚) − dep(世俗)`
是恒等式(`#838`② 已验到 1e−12),**但它的区间不是两个区间相减能得到的**:
两层在同一年里是同一次抽样切出来的,**它们的抽样误差相关,而相关的方向决定差的区间是宽还是窄。**
⇒ **本轮唯一的事:给 `dep(差距)` 算它自己的自助区间,逐十年,两把尺子都算。**

G1 估计量:**`dep_gap(d) = Δgap(d) − Δgap(全程)×d跨年÷全程跨年`**,
   自助时**同一年的两层一起重抽**(保留它们的相关),而不是各自独立重抽。

三个世界:
   A **两头都实**:1990s 的 `dep_gap` 排零(负,拉开)**且** 2010s 排零(正,合上)
     ⇒ **「先拉开、后合上」整句站得住。**
   B **只有拉开是实的**:1990s 排零而 2010s 含零
     ⇒ **「合上了」不能说** —— 而这正是 `#839` 已经把它降级到的位置,本轮把它钉死。
   C **连拉开都不排零** ⇒ ⚠ **那会动到全项目唯一穿过所有校正的那件事**,
     必须整条降级。**这是我最不欢迎的那个,所以写在前面。**

预测矩阵:
   | 世界 | 现在 | 两头都排零 | 只 1990s | 都不排零 |
   | A 先拉开后合上 | 0.25 | **0.85** | 0.10 | 0.05 |
   | B 只有拉开是实的 | 0.60 | 0.10 | **0.85** | 0.05 |
   | C 连拉开都不实 | 0.15 | 0.05 | 0.05 | **0.90** |

预注册判词(条件式):
  if 恒等式控制通过(`dep_gap` 的点值 == `dep(虔诚) − dep(世俗)`,精确)
     and 正控开火(**只往一层种位移,`dep_gap` 必须按该层的符号整体移动**)
     and 负控开火(**两层都匀速的世界里 `dep_gap` == 0**)
     and 相关控制开火(**联合重抽的区间必须与「两层独立重抽」的区间不同** ——
        若相同,说明我保留相关这件事根本没生效,那这一轮就白做):
      1990s 排零 且 2010s 排零 -> A
      只 1990s 排零            -> B
      1990s 不排零             -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`#839` 已证两把尺子给的收敛强度差 38%**
  (均值尺 +0.1384 vs 潜在尺 +0.0854)—— **所以只报一把尺子的区间是选择性报告。**
  ⇒ 控制:**两把尺子都算,两套结果并排报,不管它们同不同意。**
`G3` 多重性:5 个十年 × 2 把尺子 = 10 格,整族 BH/BY,族大小印在旁边。
⚠ 本轮换不了仪器(同一份 GSS),而它**不需要** —— 本轮是给一个已有的量补它自己的区间。
"""
import numpy as np, json, pathlib, sys, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, B, Q = "homosex", 4, 4000, 0.05

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

# 潜在尺(与 `#839` 同一套阈值构造)
Phi = lambda z: 0.5*(1.0+np.vectorize(math.erf)(z/math.sqrt(2.0)))
def ppf(p):
    lo_, hi_ = -8.0, 8.0
    for _ in range(200):
        m = (lo_+hi_)/2
        if Phi(np.array([m]))[0] < p: lo_ = m
        else: hi_ = m
    return (lo_+hi_)/2
allv = np.concatenate([np.concatenate(ys[y]) for y in S])
TAU = np.array([ppf(c) for c in np.cumsum([np.mean(allv == k) for k in range(1, KK+1)])[:-1]])
GM = np.linspace(-4.0, 4.0, 1601); ed = np.concatenate([[-np.inf], TAU, [np.inf]])
LOGP = np.log(np.clip(np.array([[Phi(np.array([ed[k+1]-m]))[0]-Phi(np.array([ed[k]-m]))[0]
                                 for k in range(KK)] for m in GM]), 1e-12, 1))
def mu(v):
    c = np.array([(v == k).sum() for k in range(1, KK+1)], float)
    return float(GM[np.argmax(LOGP @ c)])
STAT = {"均值尺 (1–4)": lambda v: float(v.mean()), "潜在尺 (有序 probit)": mu}

print(f"=== ⓪ 硬规则①:合格年 {len(S)} 个 {S[0]}–{S[-1]} · 可用十年 {sorted(dec)} · n={len(allv):,}")
print(f"  阈值 τ = {np.round(TAU,4).tolist()}(与 `#839` 同一套构造)")

def dep_gap(f, dc, rng, Bv=B, src=None, joint=True):
    """差距的十年偏离 + 它**自己的**自助区间。joint=True ⇒ 同一年两层一起重抽,保留相关。"""
    Y = src if src else ys; yy = dec[dc]
    g = lambda Z, y: f(Z[y][0]) - f(Z[y][1])
    ref = (g(Y, S[-1])-g(Y, S[0]))*(yy[-1]-yy[0])/span
    obs = (g(Y, yy[-1])-g(Y, yy[0])) - ref
    out = np.empty(Bv)
    for b in range(Bv):
        Z = {}
        for y in (S[0], S[-1], yy[0], yy[-1]):
            a, bb = Y[y]
            if joint:                                     # 同一次抽样切出的两层 ⇒ 一起重抽
                ia = rng.integers(0, len(a), len(a)); ib = rng.integers(0, len(bb), len(bb))
                Z[y] = (a[ia], bb[ib])
            else:                                         # 对照:两层独立重抽(不同的 rng 流)
                Z[y] = (a[rng.integers(0, len(a), len(a))], bb[rng.integers(0, len(bb), len(bb))])
        rf = (g(Z, S[-1])-g(Z, S[0]))*(yy[-1]-yy[0])/span
        out[b] = (g(Z, yy[-1])-g(Z, yy[0])) - rf
    return obs, out

print(f"\n=== ① 差距自己的偏离与**自己的** 95% 自助区间(B={B})· 两把尺子并排,不论同不同意 ===")
rng = np.random.default_rng(279); RES, CELLS = {}, []
for sn, f in STAT.items():
    for dc in sorted(dec):
        obs, bs = dep_gap(f, dc, rng)
        lo_, hi_ = float(np.quantile(bs, .025)), float(np.quantile(bs, .975))
        p = max(2*min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0))), 1.0/(B+1))
        RES[(sn, dc)] = dict(obs=obs, lo=lo_, hi=hi_, p=p, sd=float(np.std(bs)),
                             excl=bool(lo_ > 0 or hi_ < 0)); CELLS.append((sn, dc))
for sn in STAT:
    print(f"  {sn}")
    for dc in sorted(dec):
        r = RES[(sn, dc)]
        print(f"     {dc}s  **{r['obs']:+.4f}** [{r['lo']:+.4f},{r['hi']:+.4f}] p={r['p']:.4f}"
              f"{'  **排零**' if r['excl'] else ''}")
ps = [RES[c]["p"] for c in CELLS]
bh = {CELLS[i] for i in Gate.bh(ps, Q)}; by = {CELLS[i] for i in Gate.by(ps, Q)}
print(f"  `G3` 整族 **{len(CELLS)} 格**(5 十年 × 2 尺;⚠ 族大小印在旁边 —— `#832`:族越窄存活越易)"
      f" ⇒ BH **{len(bh)}** · BY **{len(by)}**")
print(f"     BH:{sorted(f'{s[:3]}{dc}s' for s,dc in bh) or '无'}   BY:{sorted(f'{s[:3]}{dc}s' for s,dc in by) or '无'}")
k90 = all(RES[(sn, 1990)]["excl"] for sn in STAT)
k10 = all(RES[(sn, 2010)]["excl"] for sn in STAT)
a90 = any(RES[(sn, 1990)]["excl"] for sn in STAT); a10 = any(RES[(sn, 2010)]["excl"] for sn in STAT)
print(f"  ⇒ 1990s(拉开):两把尺**都**排零 {k90}(至少一把 {a90})· "
      f"2010s(合上):两把尺**都**排零 {k10}(至少一把 {a10})")

print("\n=== ② 控制 ===")
f0 = STAT["均值尺 (1–4)"]
# ⚠⚠ **第一版这里手抄了 `#838` 的四位小数,于是恒等式控制在 2.3e−5 上失败** ——
#    失败的是**转抄的精度**,不是那个恒等式。**教训是一条可以机械执行的规矩:
#    绝不把上一轮的数字打进下一轮的脚本,从它的产物里读。**
_p838 = (ROOT/"E03_what_an_instrument_would_have_to_be/A79_那一格里到底是谁在动/"
         "R277_两层各自的轨迹从没被分开看过/results/who_actually_moved.json")
_g838 = json.load(open(_p838, encoding="utf-8"))["grid"]
D838 = {dc: (_g838[f"虔诚层|{dc}"]["obs"], _g838[f"世俗层|{dc}"]["obs"]) for dc in (1990, 2010)}
print(f"  ⚠ `#838` 的数**从它的产物读**,不手抄:{ {k: (round(v[0],6), round(v[1],6)) for k,v in D838.items()} }")
for dc in (1990, 2010):
    idn = RES[("均值尺 (1–4)", dc)]["obs"] - (D838[dc][0]-D838[dc][1])
    print(f"  恒等式控制 {dc}s:`dep_gap` = {RES[('均值尺 (1–4)',dc)]['obs']:+.6f} vs "
          f"`#838` 两层之差 {D838[dc][0]-D838[dc][1]:+.6f} ⇒ 差 **{idn:+.2e}**")
idmax = max(abs(RES[("均值尺 (1–4)", dc)]["obs"]-(D838[dc][0]-D838[dc][1])) for dc in (1990, 2010))
Yp = {y: (a.copy(), b.copy()) for y, (a, b) in ys.items()}
Yp[dec[2010][-1]] = (Yp[dec[2010][-1]][0]+0.25, Yp[dec[2010][-1]][1])
pc, _ = dep_gap(f0, 2010, np.random.default_rng(7), 200, src=Yp)
base = RES[("均值尺 (1–4)", 2010)]["obs"]
print(f"  正控:只往**虔诚层** 2010s 尾年种 +0.25 ⇒ `dep_gap` {base:+.4f} → **{pc:+.4f}** "
      f"(动 **{pc-base:+.4f}**,预期 **+0.2500** —— 该年不是全程端点 ⇒ 参照不变,"
      f"`#838`③ 那条公式已验过整张网格)")
FULL = {i: float(np.mean(ys[S[-1]][i])-np.mean(ys[S[0]][i])) for i in (0, 1)}
Yu = {y: tuple(ys[S[0]][i] + FULL[i]*(y-S[0])/span for i in (0, 1)) for y in S}
nc, _ = dep_gap(f0, 2010, np.random.default_rng(8), 200, src=Yu)
print(f"  负控:两层都严格匀速的世界 ⇒ `dep_gap` = **{nc:+.2e}** —— "
      f"⚠ **「这个零该不该是零?」该**(匀速按定义无偏离,且这一步没有离散化 ⇒ 是解析零)")
# ⚠⚠⚠ **第四条控制的前提本身是错的,而这才是本轮最该带走的一条。**
#    我写这一轮的理由是:「两层在同一年是同一次抽样切出来的 ⇒ 误差相关 ⇒
#    差的区间不能由两个区间相减得到」。**而虔诚三分位与世俗三分位是不相交的两群人** ——
#    没有任何一个受访者同时属于两层 ⇒ **它们的抽样误差按构造就是独立的。**
#    第一版的 `joint=True/False` 两个分支写出来是**同一件事**,所以 SD 之比精确是 1.0000× ——
#    **控制没有失败,是它在检验一个先验为假的说法。**
#    ⇒ 改成一条**能独立推导、可能不成立**的预测:不相交 ⇒
#       `SD(dep_gap)² == SD(dep_虔诚)² + SD(dep_世俗)²`。
#    **若成立 ⇒ 不相交被证实,而「必须联合重抽」这个前提正式撤回;
#      若不成立 ⇒ 两层之间确实有共享结构(唯一的候选是逐年三分位切点由同一年数据估出),
#      那才是一个真发现。** ⚠ 这不是看了结果再改判据:判据换成了一条我在跑之前就能推出来的算式。
def dep_stratum(f, i, dc, rng, Bv):
    yy = dec[dc]
    ref = (f(ys[S[-1]][i])-f(ys[S[0]][i]))*(yy[-1]-yy[0])/span
    out = np.empty(Bv); r = lambda a: a[rng.integers(0, len(a), len(a))]
    for b in range(Bv):
        rf = (f(r(ys[S[-1]][i]))-f(r(ys[S[0]][i])))*(yy[-1]-yy[0])/span
        out[b] = f(r(ys[yy[-1]][i])) - f(r(ys[yy[0]][i])) - rf
    return float(np.std(out))
BV = 3000
sd_gap = float(np.std(dep_gap(f0, 2010, np.random.default_rng(21), BV)[1]))
sd_h = dep_stratum(f0, 0, 2010, np.random.default_rng(22), BV)
sd_l = dep_stratum(f0, 1, 2010, np.random.default_rng(23), BV)
pred = math.sqrt(sd_h**2 + sd_l**2)
ratio = sd_gap/pred
print(f"  不相交控制(**替掉第一版那条前提为假的「相关控制」**):")
print(f"     两层是**不相交的两群人** ⇒ 误差按构造独立 ⇒ 预测 "
      f"`SD(dep_gap) = √(SD虔诚² + SD世俗²)` = √({sd_h:.5f}² + {sd_l:.5f}²) = **{pred:.5f}**")
print(f"     实测 `SD(dep_gap)` = **{sd_gap:.5f}** ⇒ 比 **{ratio:.4f}×**")
print(f"     ⇒ {'**不相交成立** ⇒ 「必须联合重抽」这个前提正式撤回' if abs(ratio-1)<0.05 else '**不相交不成立** ⇒ 两层间有共享结构(候选:逐年三分位切点由同一年数据估出)'}")

G = Gate("#840 · 「缝合上了」是关于差距的主张,而差距的区间从没算过")
# ⚠⚠ **第二版这条又失败了,而这次是库替我抓的,理由完全正确(`#773`):**
#    我把 `max|差|` 拿去和字面 0 比 —— **两边都恰好是 0,库判 DEGENERATE:
#    一个算出来的 0 和一个写死的 0 相等,不能区分「算对了」与「根本没算」。**
#    ⇒ 正确写法是 `#770` 那条:**比两个非零的值本身**,不是比它们的差与零。
_obs90 = float(RES[("均值尺 (1–4)", 1990)]["obs"])
_exp90 = float(D838[1990][0]-D838[1990][1])
G.identity_control("① 恒等式控制:`dep_gap(1990s)` 的**点值**必须等于 `#838` 两层偏离之差 —— "
                   "⚠ **比的是两个非零值本身,不是它们的差与零**(`#770`/`#773`:"
                   "一个算出来的 0 和一个写死的 0 相等,证明不了任何事);"
                   "⚠ 而这是恒等式不是发现,写下来是为了确认我算的是同一个量",
                   observed=_obs90, expected=_exp90, tol=1e-9,
                   what="dep_gap(1990s) vs `#838` 的 dep(虔诚)−dep(世俗)", deterministic=True)
G.asserted("② 正控:只往虔诚层的十年内部年种 +0.25 ⇒ `dep_gap` 必须整体上移 +0.25"
           "(该年非全程端点 ⇒ 参照不变;这条公式 `#838`③ 已在整张网格上验过)",
           bool(abs((pc-base)-0.25) < 0.02), f"动 {pc-base:+.4f}(预期 +0.2500)", kind="control")
G.asserted("③ 负控:两层都严格匀速的世界里 `dep_gap` 必须 == 0"
           "(⚠ **这个零该是零**;且这一步无离散化 ⇒ **是解析零**,可以用极严容差)",
           bool(abs(nc) < 1e-9), f"{nc:+.2e}", kind="control")
G.asserted("④ 不相交控制(**替掉第一版那条前提为假的「相关控制」**):两层是**不相交的两群人** ⇒ "
           "误差按构造独立 ⇒ **可独立推导的预测** `SD(dep_gap) = √(SD虔诚²+SD世俗²)`;"
           "⚠ **这条可能不成立** —— 若不成立,说明两层间确有共享结构(候选:逐年三分位切点由同一年数据估出)",
           bool(abs(ratio-1.0) < 0.05),
           f"实测 {sd_gap:.5f} vs 预测 {pred:.5f} = {ratio:.4f}×(虔诚 {sd_h:.5f} · 世俗 {sd_l:.5f})",
           kind="control")
G.asserted("⑤ 前提(跑前写下的最强混淆):`#839` 已证两把尺子给的收敛强度差 38% ⇒ "
           "**只报一把尺子的区间是选择性报告** ⇒ 两把都算、并排报,不论同不同意",
           bool(len(STAT) == 2 and len(CELLS) == 2*len(dec)),
           f"{len(STAT)} 把尺 × {len(dec)} 个十年 = {len(CELLS)} 格全报", kind="control")
G.asserted("⑥ kill(预注册):「先拉开、后合上」整句要成立,需 **1990s 与 2010s 的 `dep_gap` "
           "在两把尺子上都排除零**",
           bool(k90 and k10), f"1990s 都排零 {k90} · 2010s 都排零 {k10}", kind="kill",
           yardstick="`dep_gap` 每格自己的 95% 联合自助区间",
           yardstick_noise=float(np.mean([RES[c]["sd"] for c in CELLS])))
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
r90 = {sn: RES[(sn, 1990)] for sn in STAT}; r10 = {sn: RES[(sn, 2010)] for sn in STAT}
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif k90 and k10:
    V = (f"**A 先拉开、后合上 —— 整句站得住。** 1990s `dep_gap` "
         f"{' / '.join(f'{r90[sn]['obs']:+.4f}' for sn in STAT)}(两把尺都排零);"
         f"2010s {' / '.join(f'{r10[sn]['obs']:+.4f}' for sn in STAT)}(两把尺都排零)。\n"
         f"  ⇒ **一句关于人的话:那条缝是先被拉开、后来又真的合回去了一截 ——\n"
         f"  而「合回去」这一半现在有它自己的区间,不再是从两边各自的动作里读出来的。**")
elif a90 and not k10:
    V = (f"**B 只有「拉开」是实的 —— 「合上了」不能说,而这一轮把它钉死了。**\n"
         f"  1990s `dep_gap` {' / '.join(f'{r90[sn]['obs']:+.4f}' for sn in STAT)} —— "
         f"{'两把尺都排零' if k90 else '至少一把排零'};\n"
         f"  2010s {' / '.join(f'{r10[sn]['obs']:+.4f} [{r10[sn]['lo']:+.4f},{r10[sn]['hi']:+.4f}]' for sn in STAT)}"
         f" —— **含零。**\n"
         f"  ⇒ **一句关于人的话:那条缝确实是在九十年代被拉开的,这一点经得住;\n"
         f"  但「后来又合上了」我说不出口 —— 两边在二〇一〇年代都往同一个方向走了,\n"
         f"  而他们之间的距离有没有真的缩短,这份数据给不出一个排除零的答案。**")
else:
    V = (f"**C 连「拉开」都不排零 ⇒ 整条降级,而这是我最不欢迎的那个。**\n"
         f"  1990s `dep_gap` "
         f"{' / '.join(f'{r90[sn]['obs']:+.4f} [{r90[sn]['lo']:+.4f},{r90[sn]['hi']:+.4f}]' for sn in STAT)}。\n"
         f"  ⇒ **给差距算它自己的区间之后,连全项目唯一穿过所有校正的那件事都不再排零 ——\n"
         f"  之前的「排零」是两层各自与零比出来的,而那从来不是这条主张说的量。**")
print(V)
json.dump(dict(item=IT, tau=TAU.tolist(), family_size=len(CELLS), B=B, q=Q,
               grid={f"{sn}|{dc}": v for (sn, dc), v in RES.items()},
               bh=sorted(f"{s}|{dc}" for s, dc in bh), by=sorted(f"{s}|{dc}" for s, dc in by),
               key=dict(widen_1990_both=k90, close_2010_both=k10, widen_any=a90, close_any=a10),
               controls=dict(identity_max=idmax, pos=pc-base, neg=nc,
                             sd_gap=sd_gap, sd_pred=pred, sd_ratio=ratio, sd_h=sd_h, sd_l=sd_l),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"gap_own_interval.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'gap_own_interval.json'}")
