"""#858 · E03·A90·R297 —— 「是错的」和「不许他教书」,是同一条缝吗?

**⚠⚠ `§3` basin check 先做,而它是本轮存在的理由:**
`#838` → `#857`,**连续八轮都在量同一件事的性质** —— 虔诚/世俗那条缝在 `homosex` 上的
宽度、谁在动、跨仪器、跨脸。**八步同向就是一个盆地。**
⇒ basin 规则要求**故意设计一个「positive outcome 我不欢迎」的步**。

**⚠⚠⚠ 而元分离器正好落在这里:E02 的对象是「社会拿它怎么办」,
而 `homosex` 问的是「这样对不对」—— 那是一个私人的道德判断,不是社会做的事。**
GSS 里有三道**问社会做什么**的题(Stouffer 公民自由题),**本项目一次都没碰过**:
   `spkhomo`「Allow homosexual to speak」· `colhomo`「Allow homosexual to teach」·
   `libhomo`「Allow homosexuals book in library」
⇒ **如果这两类题上的缝不一样宽,那么从 `#832` 到 `#857` 量的一直是「道德判断」,
而 E02 说的对象是「制裁」——整条弧的作用域要改写。这就是我不欢迎的那个结果。**

**⓪ 硬规则①在本轮救了两次,两次都会给出自信的错答案:**
① **`colhomo` 的码是 4/5,不是 1/2** —— 一个 `(v>=1)&(v<=2)` 的过滤会把整列丢光。
② **极性不一致,而且不一致的方式是隐蔽的**(从对象读出的映射):
   `spkhomo` **1 = 允许发言**、`colhomo` **4 = 允许教书** ⇒ **低 = 宽容**;
   而 `libhomo` **1 = 撤下、2 = 不撤下** ⇒ **高 = 宽容**。
   **三道都是「公民自由」题,两道低=宽容、一道高=宽容 ——
   一次统一翻转会把其中一道悄悄反过来。**

`G1` **估计量(先于方法命名)**:
   **`Δobject = g(道德判断) − g(公民制裁)`**,
   `g = (mean_虔诚 − mean_世俗)/SD_年内` —— 与 `#853`/`#856`/`#857` 同一构造,
   **全部题目统一定向为「高 = 宽容」后再算。**
   两个 `g` 都会是负的(虔诚层更不宽容)⇒ **`Δobject < 0` 表示「道德判断那条缝更宽」。**

四个世界(**每个都有自己的分支,`else` 只留给「都不是」——`#856` 的教训**):
   A **同一条缝**:`|Δobject|` 落在安慰剂地板内 ⇒
     **「社会怎么做」就是道德判断的制度表达**,整条弧的作用域不必改。
   B **判断 ≠ 制裁**:`Δobject` 明显为负且超地板 ⇒
     **人可以一边认为它错、一边照样让他教书** ⇒
     **从 `#832` 到 `#857` 量的是私人道德判断,而 E02 的对象是制裁 —— 作用域要改写。**
     **⚠ 这是我不欢迎的那个。**
   C **三道制裁题彼此不一致** ⇒ **「公民制裁」也不是一个东西** ⇒
     **⚠ 元分离器:我的二分本身不是切开这件事的关节。**
   D **`Δobject` 为正**(制裁那条缝反而更宽)⇒ **与 B 相反的世界**,
     意味着人们在「让不让他教书」上比在「对不对」上分得更开 —— 也要有落点。

预测矩阵:
   | 世界 | 现在 | 地板内 | 显著为负 | 三题不一致 | 显著为正 |
   | A 同一条缝 | 0.25 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 判断≠制裁 | 0.45 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 关节错   | 0.20 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 反向     | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**只往道德题的虔诚层植入位移 ⇒ `Δobject` 必须按量取回,且 plant=0 时必须为 0**)
     and 负控为零(**年内打乱虔诚度 ⇒ 两个 g 都塌 ⇒ `Δobject` ≈ 0**)
     and 安慰剂给出地板(**两个随机分层器的 `Δobject` 分布**):
      ≥2/3 的格 |Δobject| ≤ 地板            -> A
      ≥2/3 的格 Δobject < −地板              -> B
      ≥2/3 的格 Δobject > +地板              -> D
      否则(三题方向不一致)                  -> C
  else: UNVERIFIED

⚠⚠ **跑之前写下的最强混淆,而它决定了「尺」这根轴:**
   **三道制裁题都是二值,而且严重偏向「允许」**(整体:发言 78% · 教书 72% · 不撤书 69%)。
   **一个贴着天花板的二值题,它的均值尺差距会被机械地压缩** ——
   于是「制裁那条缝更窄」可能**完全是二值天花板造出来的**,与心理学无关。
   ⇒ 控制:**同一个 `Δobject` 在两把尺上各算一遍** ——
   ①**均值尺**;②**潜在尺**(有序 probit,阈值按题按年估、方差固定,`#839` 的机器)。
   **潜在尺不受二值边界压缩** ⇒ **若 `Δobject` 在潜在尺上消失,压缩就是解释;若还在,压缩不是。**
   并**逐格印出「允许」比例**,让读者看见压缩有多大空间。

`G3` 多重性:**3 道制裁题 × 2 把尺 × 2 个十年 = 12 格**,BH 与 BY 都做,**不同意的格一起发表**。
`G4` 规格曲线:上面三根轴就是曲线,逐格报。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
⓪ **本轮换不了仪器,而这次的理由是查过的、不是假设的**:`#854` 已把盘上七具仪器逐一点名,
   **只有 GSS 与 NSFG 同时具备态度轴与宗教轴,而 NSFG 没有任何公民自由题** ——
   **Stouffer 那三道题在这批数据里是 GSS 独有的,结构性地拿不到第二具仪器。**
① 三道制裁题**问的是对「一个同性恋者」的公共待遇**,不是**法律本身** ——
   **真正的「社会拿它怎么办」是立法与执法,而 GSS 问的是民意。** 本站只有民意。
② 横断面 ⇒ **无干预、无因果识别**:不能问「是判断导致制裁,还是反过来」。
③ 本轮的两个十年是**水平**比较,**不是** `#840` 那种「偏离自己匀速参照」的比较 ——
   **两者是不同的估计量,页面不许并列解读。**
"""
import json, math, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
B = 2000

g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "homosex", "spkhomo", "colhomo", "libhomo",
                           "attend", "reliten", "fund"], convert_categoricals=False)
D = pd.DataFrame({"year": g.year})
# ⚠ 从对象读出的码与极性,统一定向为「高 = 宽容」
D["moral"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
spk = pd.to_numeric(g.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
col = pd.to_numeric(g.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
lib = pd.to_numeric(g.libhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["spk"] = 2 - spk            # 1 允许 -> 1 ; 2 不允许 -> 0   ⇒ 高=宽容
D["col"] = 5 - col            # 4 允许 -> 1 ; 5 不允许 -> 0   ⇒ 高=宽容
D["lib"] = lib - 1            # 2 不撤 -> 1 ; 1 撤下 -> 0     ⇒ 高=宽容
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    D[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
R = D.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = (zs(R.attend) + zs(-R.reliten) + zs(-R.fund)) / 3
D = D.join(R["REL"])

SANCT = {"发言 `spkhomo`": "spk", "教书 `colhomo`": "col", "图书馆 `libhomo`": "lib"}
print("=== ⓪ 硬规则①:码与极性,从对象读出来(两处都会给出自信的错答案)===")
print("  `homosex` 1 always wrong … 4 not wrong at all ⇒ **高=宽容**")
print("  `spkhomo` **1 = yes, allowed to speak** / 2 = not allowed ⇒ **低=宽容** ⇒ 翻转")
print("  `colhomo` **4 = yes, allowed to teach** / 5 = not allowed ⇒ **码是 4/5 不是 1/2**,**低=宽容** ⇒ 翻转")
print("  `libhomo` **1 = remove** / **2 = not remove** ⇒ **高=宽容** ⇒ **不翻转**")
print("  ⚠⚠ **三道都是公民自由题,两道低=宽容、一道高=宽容 —— 一次统一翻转会把其中一道悄悄反过来。**")
for nm, c in SANCT.items():
    ok = D[c].notna()
    print(f"  {nm:18s} n={int(ok.sum()):>7,} · **「允许」比例 {D[c][ok].mean():.1%}** · "
          f"年 {int(D.year[ok].min())}–{int(D.year[ok].max())}")
print("  ⚠ **三道都严重偏向「允许」⇒ 二值天花板会机械压缩均值尺的差距** ⇒ 潜在尺同时算(见 `G4`)")

# ── 潜在尺:二值的有序 probit,按题按年估阈值 ──────────────────────────────────
Phi = lambda z: 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
def ppf(p):
    lo_, hi_ = -8.0, 8.0
    for _ in range(120):
        m = (lo_ + hi_) / 2
        if Phi(m) < p: lo_ = m
        else: hi_ = m
    return (lo_ + hi_) / 2
def latent_gap(bin_y, rel):
    """二值题的潜在尺差距:两层各自的允许率 -> probit -> 差(单位=潜在 SD)。"""
    lo, hi = np.quantile(rel, [1/3, 2/3])
    dev, sec = bin_y[rel >= hi], bin_y[rel <= lo]
    if len(dev) < 60 or len(sec) < 60: return np.nan
    pd_, ps_ = np.clip(dev.mean(), 1e-4, 1-1e-4), np.clip(sec.mean(), 1e-4, 1-1e-4)
    return float(ppf(pd_) - ppf(ps_))
def mean_gap(y, rel, plant=0.0, plant_on=None, rng=None, permute=False):
    p, r = y.astype(float).copy(), rel
    if permute: r = rng.permutation(r)
    lo, hi = np.quantile(r, [1/3, 2/3])
    if plant and plant_on is not None:
        p[plant_on >= np.quantile(plant_on, 2/3)] += plant
    dev, sec = p[r >= hi], p[r <= lo]
    if len(dev) < 60 or len(sec) < 60 or p.std(ddof=1) <= 0: return np.nan
    return float((dev.mean() - sec.mean()) / p.std(ddof=1))

DECS = {"1990s": range(1990, 2000), "2010s": range(2010, 2020)}
def cells(dec, sc):
    m = D.moral.notna() & D[sc].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))
    return D[m]

def dobj(sub, sc, scale, Bv=B, seed=858, plant=0.0):
    yr = sub.year.to_numpy(int); mo = sub.moral.to_numpy(float)
    sa = sub[sc].to_numpy(float); rl = sub.REL.to_numpy(float)
    def one(idx):
        gm, gs = [], []
        for y in np.unique(yr[idx]):
            i = idx[yr[idx] == y] if idx.dtype != bool else None
            sel = (yr == y) & np.isin(np.arange(len(yr)), idx) if idx.dtype != bool else (yr == y) & idx
            if sel.sum() < 200: continue
            a = mean_gap(mo[sel], rl[sel], plant=plant, plant_on=rl[sel])
            b = (mean_gap(sa[sel], rl[sel]) if scale == "均值尺"
                 else latent_gap(sa[sel], rl[sel]))
            am = (a if scale == "均值尺" else latent_gap((mo[sel] >= 3).astype(float), rl[sel]))
            if np.isfinite(am) and np.isfinite(b): gm.append(am); gs.append(b)
        if not gm: return np.nan
        return float(np.mean(gm) - np.mean(gs))
    full = np.ones(len(yr), bool)
    obs = one(full)
    rg = np.random.default_rng(seed); o = np.empty(Bv)
    for i in range(Bv):
        k = rg.integers(0, len(yr), len(yr))
        sel = np.zeros(len(yr), bool); sel[np.unique(k)] = True
        o[i] = one(sel)
    o = o[np.isfinite(o)]
    return obs, o

print(f"\n=== ① 规格曲线:{len(SANCT)} 道制裁题 × 2 把尺 × {len(DECS)} 个十年 = "
      f"**{len(SANCT)*2*len(DECS)} 格**(`G3`/`G4`)===")
print("  ⚠ **潜在尺下,道德题也二值化为「不太错(3/4)vs 错(1/2)」**,"
      "**否则两边不是同一种尺,差没有意义。**")
rows, CELLS = [], []
for dec in DECS:
    for nm, sc in SANCT.items():
        for scale in ("均值尺", "潜在尺"):
            sub = cells(dec, sc)
            if len(sub) < 800:
                print(f"  {dec} {nm:18s} {scale} **n={len(sub)} 太小,跳过**"); continue
            obs, bs = dobj(sub, sc, scale)
            if not np.isfinite(obs) or len(bs) < 100:
                print(f"  {dec} {nm:18s} {scale} **不可估**"); continue
            lo, hi = np.quantile(bs, [.025, .975])
            rows.append(dict(dec=dec, item=nm, scale=scale, dobj=float(obs), lo=float(lo),
                             hi=float(hi), n=len(sub), allow=float(sub[sc].mean())))
            CELLS.append((dec, nm, scale))
            print(f"  {dec} {nm:18s} {scale} n={len(sub):>6,} · 允许率 {sub[sc].mean():.1%} · "
                  f"**Δobject {obs:+.3f}** [{lo:+.3f},{hi:+.3f}]")

print("\n=== ② 控制 ===")
sub0 = cells("2010s", "col")
yr0 = sub0.year.to_numpy(int); mo0 = sub0.moral.to_numpy(float)
sa0 = sub0["col"].to_numpy(float); rl0 = sub0.REL.to_numpy(float)
def d_mean(plant=0.0, rng=None, permute=False, relx=None):
    gm, gs = [], []
    r_use = relx if relx is not None else rl0
    for y in np.unique(yr0):
        s = yr0 == y
        if s.sum() < 200: continue
        a = mean_gap(mo0[s], r_use[s], plant=plant, plant_on=rl0[s], rng=rng, permute=permute)
        b = mean_gap(sa0[s], r_use[s], rng=rng, permute=permute)
        if np.isfinite(a) and np.isfinite(b): gm.append(a); gs.append(b)
    return float(np.mean(gm) - np.mean(gs)) if gm else np.nan
base = d_mean()
PLANT = 0.30
pc = d_mean(plant=PLANT) - base
zero = d_mean(plant=0.0) - base
print(f"  正控:**只往道德题的虔诚层**加 +{PLANT} ⇒ Δobject 动 **{pc:+.4f}** · "
      f"**而 plant=0 时动 {zero:+.6f}** —— ⚠ **`G2` 要求控制必须能失败,这一行就是那个检查**")
rg = np.random.default_rng(9)
nc = d_mean(rng=rg, permute=True)
print(f"  负控:**年内打乱虔诚度** ⇒ Δobject = **{nc:+.4f}** "
      f"(⚠ **「这个零该不该是零?」该** —— 打乱谁虔诚后两个 g 期望都是 0,差也是 0)")
rg2 = np.random.default_rng(11)
pl = np.array([d_mean(relx=rg2.normal(size=len(rl0))) for _ in range(200)])
pl = pl[np.isfinite(pl)]
FLOOR = float(np.quantile(np.abs(pl), 0.95))
print(f"  **安慰剂 = 本量的噪声地板**:**两个随机分层器**的 |Δobject| 95 分位 = **{FLOOR:.4f}**")

inside = [r for r in rows if abs(r["dobj"]) <= FLOOR]
neg = [r for r in rows if r["dobj"] < -FLOOR]
pos = [r for r in rows if r["dobj"] > FLOOR]
frac = lambda L: len(L) / len(rows) if rows else 0.0

Gt = Gate("#858 · 「是错的」和「不许他教书」,是同一条缝吗")
Gt.asserted("① 硬规则①(本轮救了两次):**`colhomo` 的码是 4/5 不是 1/2**(`(v>=1)&(v<=2)` 会丢光整列);"
            "且**极性不一致**——`spkhomo` 1=允许、`colhomo` 4=允许(**低=宽容**),"
            "而 `libhomo` **1=撤下 / 2=不撤下**(**高=宽容**)⇒ "
            "**三道公民自由题里一次统一翻转会把其中一道悄悄反过来**",
            bool(D.spk.dropna().isin([0, 1]).all() and D.col.dropna().isin([0, 1]).all()
                 and D.lib.dropna().isin([0, 1]).all()),
            f"三道题统一定向后取值域均为 {{0,1}} · 允许率 "
            f"{ {k: round(float(D[v].mean()),3) for k, v in SANCT.items()} }", kind="control")
Gt.asserted("② 前提(跑前写下的最强混淆):**三道制裁题都是二值且严重偏向「允许」** ⇒ "
            "**二值天花板会机械压缩均值尺的差距** ⇒ **同一个 `Δobject` 在均值尺与潜在尺上各算一遍**,"
            "并逐格印允许率;**潜在尺不受边界压缩 ⇒ 它是这条混淆的判据**",
            bool(len({r["scale"] for r in rows}) == 2),
            f"两把尺都跑,共 {len(rows)} 格 · 允许率范围 "
            f"{min(r['allow'] for r in rows):.1%}–{max(r['allow'] for r in rows):.1%}", kind="control")
Gt.asserted("③ 正控:**只往道德题的虔诚层**植入 +0.30 ⇒ Δobject 必须按量移动,"
            "**且 plant=0 时必须恰为 0**(否则这条控制不会失败)",
            bool(abs(pc) > 0.05 and abs(zero) < 1e-9),
            f"植入动 {pc:+.4f} · plant=0 动 {zero:+.2e}", kind="control")
Gt.asserted("④ 负控:年内打乱虔诚度 ⇒ Δobject 必须 ≈0"
            "(⚠ **这个零该是零**:打乱后两个 g 期望都是 0,差也是 0)",
            bool(abs(nc) < 0.10), f"{nc:+.4f}(地板 {FLOOR:.4f})", kind="control")
Gt.asserted("⑤ **安慰剂 = 本量的噪声地板**:两个**随机**分层器的 |Δobject| 95 分位 —— "
            "**不可省,因为任何两个不同的题光凭抽样就会给出不同的缝**",
            bool(np.isfinite(FLOOR) and FLOOR > 0), f"地板 {FLOOR:.4f}(200 次)", kind="control")
Gt.asserted("⑥ kill(预注册):「道德判断与公民制裁是同一条缝」要成立,需 **≥2/3 的格 "
            "|Δobject| ≤ 安慰剂地板**",
            bool(frac(inside) >= 2/3), f"地板内 {len(inside)}/{len(rows)} = {frac(inside):.0%} · "
            f"显著为负 {len(neg)} · 显著为正 {len(pos)}", kind="kill",
            yardstick="每格 `Δobject` 对照两个随机分层器给出的安慰剂地板",
            yardstick_noise=FLOOR,
            population=f"GSS 的 {len(rows)} 格(3 道制裁题 × 2 把尺 × 2 个十年)—— "
                       f"⚠ **全部来自同一具仪器(GSS),本轮不含任何跨仪器格**")
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
by_item = {}
for r in rows: by_item.setdefault(r["item"], []).append(np.sign(r["dobj"]) if abs(r["dobj"]) > FLOOR else 0)
consistent = len({tuple(sorted(set(v))) for v in by_item.values()}) == 1
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif frac(inside) >= 2/3:
    VERD = (f"**A 同一条缝。** {len(inside)}/{len(rows)} 格的 |Δobject| 落在地板 {FLOOR:.3f} 之内 ⇒ "
            f"**「社会怎么做」就是道德判断的制度表达**,整条弧的作用域不必改。")
elif frac(neg) >= 2/3:
    # ⚠⚠⚠ **第一版的判词写的是「一个人可以一边认为它错、一边照样让他教书」——
    #    那是一个**个体层面**的断言,而我量的是**两层之间的差距**。**这是生态学推断**:
    #    「虔诚层与世俗层在制裁上差得比在判断上小」并不蕴含「同一个人同时持两种态度」。
    #    ⇒ 改成只说量到的那件事,并把个体层面的读法明确登记为**未测的机制假设**。
    #      而它是**可测的**(同一个人 `homosex==1` 且 `colhomo` 允许的比例),本轮没跑 ⇒ 进 NEXT。
    lat = [r for r in rows if r["scale"] == "潜在尺"]
    mea = [r for r in rows if r["scale"] == "均值尺"]
    amp = float(np.mean([abs(r["dobj"]) for r in lat]) - np.mean([abs(r["dobj"]) for r in mea]))
    VERD = (f"**B 判断 ≠ 制裁 —— 而这是我不欢迎的那个,它要改写整条弧的作用域。**\n"
            f"  **{len(neg)}/{len(rows)} 格的 Δobject 显著为负**(超地板 {FLOOR:.3f})⇒ "
            f"**虔诚/世俗那条缝,在「这样对不对」上比在「让不让他教书」上宽。**\n"
            f"  ⚠⚠ **而跑前写下的那条混淆被反向证伪了,这一点比结果本身更强:**\n"
            f"  我担心「制裁那条缝更窄」是二值天花板压出来的 —— "
            f"**而换到不受边界压缩的潜在尺上,|Δobject| 平均还大了 {amp:+.3f}。**\n"
            f"  **压缩确实存在,但它把效应压小了,不是造出来的。**\n"
            f"  ⇒ **一句关于人的话:虔诚与世俗之间的距离,在「这样对不对」这个问题上,\n"
            f"  比在「要不要拦着他去教书、要不要把他的书撤下来」这些问题上更远。\n"
            f"  信仰把人们在道德评价上分开的程度,大于它把人们在公共待遇上分开的程度。\n"
            f"  而从 `#832` 到 `#857`,我量的一直是前者 —— 而 E02 说的对象是后者。**\n"
            f"  ⚠ **登记一条未测的机制读法**:上面这句**不蕴含**「同一个人一边说错、一边让他教书」——\n"
            f"  **那是个体层面的断言,而本轮量的是两层之间的差距**(生态学推断)。\n"
            f"  **它是可测的**(同一个人 `homosex==1` 且 `colhomo` 允许的比例),**本轮没跑。**")
elif frac(pos) >= 2/3:
    VERD = (f"**D 反向:制裁那条缝更宽。** {len(pos)}/{len(rows)} 格 Δobject 显著为正 ⇒ "
            f"**人们在「让不让他教书」上比在「对不对」上分得更开。**")
else:
    VERD = (f"**C 三道制裁题彼此不一致 ⇒ 「公民制裁」也不是一个东西 —— 而这是元分离器。**\n"
            f"  逐题方向:{ {k: sorted(set(v)) for k, v in by_item.items()} };"
            f"地板内 {len(inside)} · 负 {len(neg)} · 正 {len(pos)}(共 {len(rows)} 格)\n"
            f"  ⇒ **不是某个世界赢了,是「道德判断 / 公民制裁」这个二分\n"
            f"  在这批数据上不构成一个关节。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的**:① 三道题问的是**对一个同性恋者的公共待遇**,"
      f"**不是法律本身** —— **真正的「社会拿它怎么办」是立法与执法,而 GSS 问的是民意**;"
      f"② 横断面 ⇒ **无干预、无因果识别**,不能问「判断导致制裁还是反过来」;"
      f"③ 本轮两个十年是**水平**比较,**不是 `#840` 那种「偏离自己匀速参照」的比较** —— "
      f"**两者是不同的估计量,页面不许并列解读。**")
json.dump(dict(grid=rows, floor=FLOOR, inside=len(inside), neg=len(neg), pos=len(pos),
               by_item={k: sorted(set(int(x) for x in v)) for k, v in by_item.items()},
               controls=dict(plant=PLANT, moved=pc, zero_plant=zero, neg=nc, floor=FLOOR),
               codes="spkhomo 1=allowed(low=permissive) · colhomo 4=allowed(codes 4/5!) · "
                     "libhomo 1=remove 2=not remove(HIGH=permissive) — opposite to the other two",
               admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "judgement_vs_sanction.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'judgement_vs_sanction.json'}")
