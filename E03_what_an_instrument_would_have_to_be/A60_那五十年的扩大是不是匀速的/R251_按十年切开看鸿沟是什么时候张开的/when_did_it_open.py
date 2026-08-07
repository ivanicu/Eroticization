"""#812 · E03·A60·R251 —— 那五十年的扩大是不是匀速的?按十年切开,看这条鸿沟是什么时候张开的

`#810` 留下一条**明确标为事后**的观察:`#808` 测得 `homosex` 五十年的 `Δgap` = **−1.103**
⇒ 摊到每年约 −0.022 ⇒ NSFG 那个 6 年窗口该看到约 −0.17,**而 NSFG 的区间把它排除在外**。
⇒ 三个世界并存,而 `#810` 分不开它们:
   **① 扩大不是匀速的(发生在更早)· ② 它对 15–49 岁女性不成立 · ③ 两具仪器真的不一致。**

⚠⚠ **而世界①可以只用 GSS 自己分开,不需要第二具仪器 —— 这才是本轮值钱的地方:**
   **如果 GSS 自己的 2010 年代切片也几乎没动,那么两具仪器是一致的,而「匀速」是错的。**
   **如果 GSS 的 2010 年代仍在扩大而 NSFG 没有,那才是仪器不一致。**
   ⇒ **一个只用一具仪器的动作,把「仪器不一致」这个世界杀掉或留下 —— 这是本轮的分离器。**

G1 估计量:**逐十年的 `Δgap`(该十年末的两层差距 − 该十年初的两层差距)**,
   与**匀速模型的预测值**比较。

⚠⚠ **参照值不是 0,是 −0.2206** —— 即五十年总量 −1.103 摊到每个十年。
   **问的不是「这个十年动了没有」,而是「这个十年的移动是不是等于匀速模型说的那么多」。**
   ⇒ 用 `#811` 刚建的 `Gate.interval_verdict(lo, hi, reference=-0.2206, matters=…)`,
   **`EXCLUDES` = 这个十年偏离匀速。**
   ⚠ **而 `matters` 必须显式给**(`#811` 立的规矩):这里取 **0.11 = 匀速率的一半** ——
   **「一个十年走的量与匀速模型差了半个十年的量」,这是一个能被读者理解的偏离尺度。**
   **它是我选的,写下来,不藏。**

三个世界:
   A **匀速**:所有十年的区间都覆盖 −0.2206 ⇒ **这条鸿沟是稳定地、一代人接一代人地张开的** ——
     而那样 NSFG 的窄零就**真的**与 GSS 冲突,世界③活。
   B **集中在某几个十年**:至少一个十年 `EXCLUDES` ⇒ **扩大有它自己的年代** ——
     **那是一句关于人的话:这条鸿沟不是一直在长,它是在某个时候裂开的。**
   C **2010 年代本身就是一个窄零**:该十年 `TIGHT_NULL` ⇒ **GSS 与 NSFG 一致** ⇒
     **世界③(仪器不一致)被杀掉,而 `#810` 那条事后观察得到解释。**

预测矩阵:
   | 世界 | 现在 | 全覆盖 −0.22 | 有十年 EXCLUDES | 2010s 是窄零 |
   | A 匀速         | 0.30 | **0.85** | 0.05 | 0.05 |
   | B 有年代       | 0.45 | 0.05 | **0.90** | 0.45 |
   | C 两仪器一致   | 0.25 | 0.03 | 0.40 | **0.90** |

预注册判词(条件式):
  if 正控开火(合成一个**只在某一个十年**发生的位移,判据必须只在那个十年 `EXCLUDES`)
     and 负控开火(合成一个**真匀速**的世界,所有十年必须都覆盖参照):
      逐十年报三值 + 整张网格;**总判按计数,不设「多数」阈值**(`#805` 的教训)
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**每个十年的年份数不同**(GSS 有些年没问 `homosex`,2020 年后更稀)
  ⇒ 一个只有 2 个年份的十年,它的区间宽是**年份少**造成的,不是**没变化**。
  ⇒ 控制:**每个十年并排印出它实际用到的年份数与年份表**,而**年份数 < 3 的十年一律标
  `UNRESOLVED` 且不参与计数** —— 这一条写在跑之前。

⚠ 硬规则①:先打印每个十年的 n、真正被问过的年份、两层的样本量。
⚠ 本轮**只用 GSS 一具仪器,而这正是设计的一部分**(要杀掉的正是「仪器不一致」那个世界)。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(251)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK = "homosex", 4
UNIFORM = -1.103/5.0          # 五十年 −1.103 摊到每个十年(`#808` 的 `fit` 口径)
MATTERS = 0.11                # ⚠ `#811` 强制显式给:匀速率的一半,我选的,写下来
B = 3000

print("=== ⓪ 硬规则①:每个十年的年份 · 两层样本量 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
G = REL.dropna(subset=[IT])
YR = {}
for y, gy in G.groupby("year"):
    a, b = gy[gy.k == 2], gy[gy.k == 0]
    if len(a) >= 120 and len(b) >= 120:
        YR[int(y)] = (a[IT].to_numpy(float), b[IT].to_numpy(float))
DEC = {}
for y in sorted(YR):
    DEC.setdefault((y//10)*10, []).append(y)
for dc in sorted(DEC):
    ys = DEC[dc]
    print(f"  {dc}s: 年份 {len(ys)} 个 {ys} · 虔诚 n={sum(len(YR[y][0]) for y in ys):,} · "
          f"非虔诚 n={sum(len(YR[y][1]) for y in ys):,}")

def gap(y): a, b = YR[y]; return float(a.mean()-b.mean())
def dec_dgap(ys, src=None):
    """该十年内:末年差距 − 首年差距。⚠ 年份 < 3 的十年由调用方标 UNRESOLVED。"""
    S = src if src is not None else YR
    f = lambda y: float(S[y][0].mean()-S[y][1].mean())
    return f(ys[-1]) - f(ys[0])
def boot_dec(ys):
    r = lambda a: a[RNG.integers(0, len(a), len(a))]
    S = {y: (r(YR[y][0]), r(YR[y][1])) for y in (ys[0], ys[-1])}
    return dec_dgap([ys[0], ys[-1]], src=S)

print(f"\n=== ① 逐十年 `Δgap` 对**匀速模型**的参照 {UNIFORM:+.4f}(不是对 0!)· "
      f"「多大才算数」= {MATTERS}(`#811` 强制显式给,我选的)===")
ROWS = []
for dc in sorted(DEC):
    ys = DEC[dc]
    if len(ys) < 3:
        ROWS.append(dict(decade=dc, n_years=len(ys), years=ys, verdict="UNRESOLVED",
                         why="年份 < 3 —— 跑前写下的规矩,不参与计数", dgap=None, lo=None, hi=None))
        print(f"  {dc}s  年份 {len(ys)} 个 ⇒ **UNRESOLVED(年份 < 3,跑前定的规矩,不参与计数)**")
        continue
    pt = dec_dgap(ys)
    dr = np.array([boot_dec(ys) for _ in range(B)]); dr = dr[np.isfinite(dr)]
    lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
    # ⚠⚠ 第二处错,而且它同时污染了真实分析,不只是控制:
    #    每个十年**实际跨的年数不同**(有的十年只有 1974→1978),
    #    而我拿全部十年去比同一个「一个十年的匀速量」⇒ **跨得短的十年必然判成偏离。**
    #    ⇒ 参照必须按实际跨度缩放:`UNIFORM × 实际跨年 / 10`。
    ref = UNIFORM*(ys[-1]-ys[0])/10.0
    v = Gate.interval_verdict(lo, hi, reference=ref, matters=MATTERS)
    v0 = Gate.interval_verdict(lo, hi, reference=0.0, matters=MATTERS)
    ROWS.append(dict(decade=dc, n_years=len(ys), years=ys, dgap=float(pt), lo=lo, hi=hi,
                     span=ys[-1]-ys[0], reference=float(ref), verdict=v, verdict_vs_zero=v0, why=""))
    print(f"  {dc}s  {ys[0]}→{ys[-1]}(跨 {ys[-1]-ys[0]} 年 · {len(ys)} 个年份) · "
          f"Δgap = **{pt:+.4f}** [{lo:+.4f}, {hi:+.4f}] · 匀速参照 **{ref:+.4f}**(按实际跨度缩放)"
          f" ⇒ 对匀速 **{v}** · 对零 **{v0}**")

USE = [r for r in ROWS if r["dgap"] is not None]
n_excl = sum(1 for r in USE if r["verdict"] == "EXCLUDES")
d2010 = next((r for r in ROWS if r["decade"] == 2010), None)
print(f"\n  参与计数的十年 **{len(USE)}/{len(ROWS)}** · 偏离匀速(`EXCLUDES`)的十年 **{n_excl}**")

print("\n=== ② 控制(合成世界,同一条代码路径)===")
# ⚠⚠ 第一版的合成器**把结构叠加在观测数据上,而不是替换掉它** ——
#    「匀速世界」于是等于「真实世界 + 一条匀速趋势」,**它按构造就不是匀速的**,
#    负控因此不可能通过(实测 1990s/2000s 判 EXCLUDES)。
#    `realstat` 的原话是 **destroy the structure under test, preserve everything else** ——
#    而我把被测的结构**保留**了。⇒ 改成:直接把每年的差距**设成**目标值。
def syn(kind):
    """kind='uniform' 差距逐年按匀速率走;kind='onedecade' 全部位移集中在 1990s。"""
    S = {}; y0 = min(YR); g0 = gap(y0)
    for y in sorted(YR):
        a, b = YR[y]
        if kind == "uniform": target = g0 + UNIFORM*(y-y0)/10.0
        else: target = g0 + (-1.103)*(0.0 if y < 1990 else min((y-1990)/9.0, 1.0))
        S[y] = (a + (target - gap(y)), b.copy())      # ⚠ 替换,不是叠加
    return S
for kind, want in (("uniform", "全覆盖参照"), ("onedecade", "只有 1990s 偏离")):
    S = syn(kind); vs = {}
    for dc in sorted(DEC):
        ys = DEC[dc]
        if len(ys) < 3: continue
        pt = dec_dgap(ys, src=S)
        r = lambda a: a[RNG.integers(0, len(a), len(a))]
        dr = np.array([dec_dgap([ys[0], ys[-1]], src={y: (r(S[y][0]), r(S[y][1])) for y in (ys[0], ys[-1])})
                       for _ in range(800)])
        lo, hi = float(np.percentile(dr, 2.5)), float(np.percentile(dr, 97.5))
        vs[dc] = (float(pt), Gate.interval_verdict(lo, hi, UNIFORM*(ys[-1]-ys[0])/10.0, MATTERS))
    print(f"  {kind:10s}(该 {want}):" + " · ".join(f"{k}s {v[0]:+.3f}/{v[1][:4]}" for k, v in vs.items()))
    if kind == "uniform": nc_ok = all(v[1] != "EXCLUDES" for v in vs.values()); NCV = vs
    else: pc_ok = vs.get(1990, (0, ""))[1] == "EXCLUDES"; PCV = vs

Gg = Gate("#812 · 那五十年的扩大是不是匀速的")
Gg.asserted("① 正控:合成一个**只在 1990s** 发生全部位移的世界,1990s 必须判 `EXCLUDES`",
            bool(pc_ok), f"1990s ⇒ {PCV.get(1990)}", kind="control")
Gg.asserted("② 负控:合成一个**真匀速**的世界,所有十年都不许判 `EXCLUDES`"
            "(⚠ 参照是匀速率 −0.2206,**不是 0** —— `#811` 那条「这个零该不该是零」在这里是「这个参照该不该是零」)",
            bool(nc_ok), " · ".join(f"{k}s {v[1]}" for k, v in NCV.items()), kind="control")
Gg.asserted("③ 前提(跑前写下的混淆):每个十年并排印出年份数,**年份 < 3 一律 UNRESOLVED 且不参与计数**",
            bool(all(r["verdict"] == "UNRESOLVED" for r in ROWS if r["n_years"] < 3)),
            f"参与计数 {len(USE)}/{len(ROWS)} · 被排除的十年 {[r['decade'] for r in ROWS if r['dgap'] is None]}",
            kind="control")
Gg.asserted("④ 前提:「多大才算数」显式给出并印在行里(`#811`)",
            bool(MATTERS > 0), f"matters = {MATTERS}(匀速率 {UNIFORM:+.4f} 的一半,我选的)", kind="control")
Gg.asserted("⑤ kill(预注册):「匀速」要成立,需**所有参与计数的十年都覆盖它自己那个按跨度缩放的匀速参照**",
            bool(n_excl == 0), f"偏离匀速的十年 {n_excl}/{len(USE)}", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*98)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
else:
    ex = [f"{r['decade']}s({r['dgap']:+.3f})" for r in USE if r["verdict"] == "EXCLUDES"]
    if n_excl == 0:
        V = ("**A 匀速。所有参与计数的十年都覆盖 −0.2206。**\n"
             "  ⇒ **这条鸿沟是稳定地、一代人接一代人地张开的。**\n"
             "  ⚠⚠ **而那意味着 NSFG 那个窄零真的与 GSS 冲突 ——「两具仪器不一致」这个世界活下来了,\n"
             "  它现在是最大的未解项。**")
    else:
        V = (f"**B 这条鸿沟有它自己的年代。** 偏离匀速的十年:**{ex}**。\n"
             f"  ⇒ **一句关于人的话:虔诚者与其余人在同性恋这道题上的距离,不是五十年里一点点长出来的 ——\n"
             f"  它在某个时候裂开,而其余的时间它走的和匀速模型说的不一样。**")
    if d2010 and d2010["dgap"] is not None:
        V += (f"\n  ⚠⚠ **而 2010 年代那一格是本轮的分离器**:Δgap = **{d2010['dgap']:+.4f}** "
              f"[{d2010['lo']:+.4f}, {d2010['hi']:+.4f}] ⇒ 对匀速 **{d2010['verdict']}** · 对零 **{d2010['verdict_vs_zero']}**\n  ⇒ "
              + ("**GSS 自己的 2010 年代对零也是一个窄零 ⇒ 与 NSFG 一致 ⇒ 「两具仪器不一致」那个世界被杀掉,\n"
                 "  而 `#810` 那条事后观察得到解释:扩大发生在更早。**"
                 if d2010["verdict_vs_zero"] == "TIGHT_NULL" else
                 "**GSS 的 2010 年代对零不是窄零 ⇒ 它与 NSFG 的窄零并不一致 ⇒\n"
                 "  「两具仪器不一致」仍然活着,而下一轮该去查人群(15–49 岁女性)那一条。**"))
    elif d2010:
        V += f"\n  ⚠ **2010 年代那一格不参与计数**({d2010['why']})⇒ **本轮分不开 `#810` 那三个世界,如实说。**"
print(V)
json.dump(dict(item=IT, uniform_reference=UNIFORM, matters=MATTERS, B=B, rows=ROWS,
               n_used=len(USE), n_excludes=n_excl, decade_2010=d2010,
               pos_control=PCV, neg_control=NCV, admissible=adm, verdict=V, gate_ok=Gg.verdict()),
          open(OUT/"when_did_it_open.json", "w"), ensure_ascii=False, indent=1, default=str)
print(f"\n  产物 → {OUT/'when_did_it_open.json'}")
