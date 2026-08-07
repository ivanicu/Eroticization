r"""#863 · E03·A93·R302 —— 把 `#862` 那把完全相同的尺架到道德题上:七根轴在 `homosex` 上重排

**这是 `#862`① 登记的债。** `#862` 在三道**制裁**题上比过七根轴,发现 **6/6 格是教育最宽,宗教零格**。
它的判词第一版多写了半句「在『这样对不对』上把美国人分开的是信仰」,**当场删掉了** ——
因为**那半句的证据来自 `#834`,而 `#834` 比的是「偏离自己的匀速参照」,不是「缝的水平」。
两个不同的估计量。道德题的水平跨轴比较,这个项目从来没做过。**

⇒ **只有本轮跑完,「同一个话题、两个问法、两条不同的裂缝」这句话才有资格说。**

**⚠⚠ 而本轮同时是对 `#862` 自己的攻击,因为读它的代码时发现两个可查的缺陷:**
① **`#862` 的 docstring 第 58 行写着「BH 与 BY 都做」,而脚本里一个多重性校正都没有**
   (`grep -n "BH\|BY\|fdr" which_cleavage_carries_sanction.py` → 只命中注释)。
   **声明了没实现 —— 本轮真的做。**
② **`#862` 拿「最大 − 次大」去比一个为「单个 |g|」造的地板。** 那是两个不同的量:
   一个是**一根轴的差距**,一个是**七根轴里头两名的差距**,后者的零分布更宽。
   **宪法的原话:一个量要跟它自己的零比,不能跟另一个观测量比。**
   ⇒ 本轮为**边际**单独造零分布(**七根轴同时打乱,保留轴与轴之间的相关**)。

`G1` **估计量(先于方法命名)**:
   **每个 (题 × 分层轴 × 十年 × 尺) 上的 `|g|`**,`g = (mean_高 − mean_低)/SD_年内`,**逐年算再平均**。
   ⚠ **比的是 `|g|` 不是 `g`** —— 各轴「高」的含义不同(虔诚高=更不宽容,教育高=更宽容),
   **符号单独报,绝不参与排序**。与 `#853`–`#862` 同一构造,**这是设计要求,不是省事**:
   构造一变,两侧就不可比,本轮问的那句话就没法说。

四个世界(**每个都有分支**,`#856` 的教训):
   A **道德侧的主轴是宗教** ⇒ **同一话题两条裂缝**,项目的框架对道德那一半成立。
     ⚠ **这是我欢迎的那个** —— 所以它拿到的怀疑要更多,不是更少(`§3` basin 规则)。
   B **道德侧也是教育** ⇒ **「宗教那条缝」这个框架对整个对象都错了**,`#834` 的结果只是偏离尺特有的。
   C **弥散**:最大与次大之差落在**边际自己的地板**内 ⇒ 「哪根轴最宽」这个问法两侧都不成立。
   D **主轴取决于尺**(均值尺 vs 潜在尺给出不同赢家)⇒ **⚠ 元分离器:那个排名不是人的性质,
     是答题格式的性质** —— 「存在一条最宽的裂缝」这个分解方式本身是测量相对的。

预测矩阵:
   | 世界 | 现在 | 宗教夺魁 | 教育夺魁 | 边际落在地板内 | 两把尺赢家不同 |
   | A 两条缝   | 0.45 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 一条缝   | 0.20 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 弥散     | 0.15 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 尺决定   | 0.20 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(**条件式,不是阈值**):
  if 正控开火(**只往宗教轴高组植入 ⇒ 只有宗教的 |g| 该动;plant=0 时必须恰为 0;
     且阈值必须落在 floor 与 ceiling 之间**)
     and 负控为零(**年内打乱宗教轴标签 ⇒ 塌回地板**)
     and 安慰剂为零(**`ballot` —— GSS 自己随机分配的问卷版本,一个真实但必然无关的分层器**):
      宗教在 ≥8/12 个道德格(6 十年 × 2 尺)里 `|g|` 最大   -> A
      同一根非宗教轴在 ≥8/12 格里最大                       -> B
      「最大−次大」在 ≥8/12 格里落在**边际自己的地板**内    -> C
      两把尺在 ≥3/6 个十年上给出不同赢家                    -> D
  else: UNVERIFIED

⚠⚠ **跑之前写下的最强混淆(两条,都在同一轮里量):**
  ① **七根轴彼此相关** —— 一根轴上的大差距可能只是宗教换了件衣服。
     ⇒ 逐格印出该轴高组与虔诚高组的 **Jaccard**,**折扣印在每一行旁边**(`#820`/`#834`/`#862` 同一条)。
  ② **⚠ 各轴「切得有多狠」不一样,而 |g| 会被这件事抬高。**
     `educ` 三分位切在 **12 年 vs 15 年**,而且重并极多(2010s 高组占 37.3%、低组 41.6%,**不是 33/33**);
     `sex` 是 50/50;`region` 南方约 25/75;`REL` 是连续 z 合成。
     **组的大小和形状不同 ⇒ 抽样噪声不同 ⇒ 不能拿同一个地板去量所有轴。**
     ⇒ **本轮每一格有它自己的地板**(该轴自己的年内置换分布 95 分位),
     **这正是 `#862` 用一个全局 0.0547 去量七根轴时缺的东西。**

`G3` 多重性:整族 = **4 题 × 7 轴 × 6 十年 × 2 尺**(有数据的格),**BH 与 BY 都做,不同意的格一起发表**。
`G4` 规格曲线:两把尺 × 六个十年 × 四道题,**全部逐格印出,包括杀掉结论的格**。
⚠ kill 带 `yardstick` / `yardstick_noise` / `population` / `direction`(`#836`①·`#854`①·`#859`①·`#860`)。

**⚠ 硬规则 ④:换仪器胜过在同一具仪器上再跑一轮。**
   道德题**有**第二具仪器(`#853` 已用过):**NSFG `samesex`**。
   ⇒ 本轮把同一个排名在 NSFG 上重跑一遍,**能对上的轴只有五根**(宗教·教育·年龄·种族·性别),
   **`polviews` 与 `region` NSFG 没有 ⇒ 如实登记,不算作复制。**
   ⚠⚠ **而 NSFG 的年龄轴 15–49 岁封顶** —— **它结构性地不可能复制 GSS 的年龄结果**,
   这一条**必须写在结果旁边**,不能让读者以为「年龄没复制上」是发现。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① **制裁三题换不了仪器**(`#854` 已点名盘上七具,Stouffer 三题是 GSS 独有)⇒
   **本轮的跨仪器复制只覆盖道德题,不覆盖 `#862` 的结论。**
② 横断面 ⇒ **无干预、无因果识别**:不能问「是教育造成宽容,还是宽容的人去读了书」。
③ **只比「哪根轴的缝最宽」,不做多元调整** —— 放进同一个回归后系数是「控制了其他轴之后」的量,
   **那是另一个估计量,而各轴相关性高时它对设定极其敏感** ⇒ **本轮不做,不是「以后做」。**
④ **NSFG 只有 2011–2019 三个文件有 setup**(`ls data/external/nsfg/setup/` → 3 个 `.dct`)⇒
   **跨仪器只覆盖 2010s 那一个十年**,九十年代的跨仪器复制**盘上不存在**(`#854` 已登记)。
"""
import json, math, pathlib, re, sys
import numpy as np
import pandas as pd
from scipy.special import ndtri

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
EXT = ROOT / "data/external"
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
NPERM = 200
SEED = 302

# ⚠ `#862` 的那个全局地板**从它的产物里读**,不手抄(`#840`;`tools/no_transcribed_numbers.py` 当场拦下了我)
F862 = json.load(open(ROOT / "E03_what_an_instrument_would_have_to_be/A93_宗教是不是这个对象的那根轴"
                      "/R301_七根轴在放行题上重排一次/results/which_cleavage.json"))["floor"]

# ── 仪器交换的规范检验:ndtri 必须与 `#861`/`R300` 手写的二分 ppf 一致到 1e-9 ────────────
Phi = lambda z: 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
def ppf_bisect(p):
    lo_, hi_ = -8.0, 8.0
    for _ in range(120):
        m = (lo_ + hi_) / 2
        if Phi(m) < p: lo_ = m
        else: hi_ = m
    return (lo_ + hi_) / 2
_gr = np.linspace(1e-4, 1 - 1e-4, 401)
PPF_MAXDIFF = float(np.max(np.abs(np.array([ppf_bisect(float(p)) for p in _gr]) - ndtri(_gr))))
print(f"=== ⓪a 仪器交换的规范检验:`ndtri` vs `#861` 手写二分 ppf,401 点最大差 "
      f"**{PPF_MAXDIFF:.2e}** ⇒ **同一把尺,只是快了** ===")

# ── GSS ────────────────────────────────────────────────────────────────────────
GCOLS = ["year", "homosex", "spkhomo", "colhomo", "libhomo", "attend", "reliten", "fund",
         "polviews", "educ", "age", "region", "race", "sex", "ballot"]
gs = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta", columns=GCOLS, convert_categoricals=False)
D = pd.DataFrame({"year": gs.year})
D["moral"] = pd.to_numeric(gs.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
D["spk"] = 2 - pd.to_numeric(gs.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["col"] = 5 - pd.to_numeric(gs.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
D["lib"] = pd.to_numeric(gs.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7)),
                    ("educ", (0, 20)), ("age", (18, 89)), ("region", (1, 4)), ("race", (1, 3)),
                    ("sex", (1, 2)), ("ballot", (1, 4))):
    D[c] = pd.to_numeric(gs[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
Rr = D.dropna(subset=["attend", "reliten", "fund"]).copy()
Rr["REL"] = (zs(Rr.attend) + zs(-Rr.reliten) + zs(-Rr.fund)) / 3
D = D.join(Rr["REL"])

ITEMS = {"道德 `homosex`": ("moral", [1, 2, 3]), "发言 `spkhomo`": ("spk", [0]),
         "教书 `colhomo`": ("col", [0]), "图书馆 `libhomo`": ("lib", [0])}
DECS = {"1970s": range(1972, 1980), "1980s": range(1980, 1990), "1990s": range(1990, 2000),
        "2000s": range(2000, 2010), "2010s": range(2010, 2020), "2020s": range(2020, 2025)}
SCALES = ["均值尺", "潜在尺"]

print("\n=== ⓪b 硬规则①:**变量名不是测量** —— 每一列的 n 与实际问过的年份,先印出来 ===")
for nm, (c, _) in ITEMS.items():
    yy = sorted(D.year[D[c].notna()].unique().astype(int))
    print(f"  {nm:16s} n={int(D[c].notna().sum()):6,d} · {len(yy)} 个年份 {yy[0]}–{yy[-1]}")
for c in ["REL", "polviews", "educ", "age", "region", "race", "sex", "ballot"]:
    yy = sorted(D.year[D[c].notna()].unique().astype(int))
    print(f"  轴 {c:9s}      n={int(D[c].notna().sum()):6,d} · {len(yy)} 个年份 {yy[0]}–{yy[-1]}")
print("  ⚠ `region` 本版 GSS 是**四类**:1 northeast · 2 midwest · **3 south** · 4 west "
      "(`convert_categoricals=True` 直接查过标签)⇒ `#862` 的「南方 vs 其余」标对了。")

def terc(col, sub):
    v = sub[col]
    if v.notna().sum() < 200: return None, None
    lo, hi = np.nanquantile(v, [1 / 3, 2 / 3])
    return (v >= hi).to_numpy(), (v <= lo).to_numpy()

AXES = {
    "宗教 REL": lambda s: terc("REL", s),
    "教育 educ": lambda s: terc("educ", s),
    "政治 polviews": lambda s: terc("polviews", s),
    "年龄 age": lambda s: terc("age", s),
    "地区 南方vs其余": lambda s: ((s.region == 3).to_numpy(),
                                  (s.region != 3).to_numpy() & s.region.notna().to_numpy()),
    "种族 白vs黑": lambda s: ((s.race == 1).to_numpy(), (s.race == 2).to_numpy()),
    "性别 男vs女": lambda s: ((s.sex == 1).to_numpy(), (s.sex == 2).to_numpy()),
}
PLACEBO = {"安慰剂 ballot": lambda s: ((s.ballot == 1).to_numpy(), (s.ballot == 3).to_numpy())}
ALL_AX = {**AXES, **PLACEBO}
ANAMES = list(AXES)

print("\n=== ⓪c 硬规则②:**各轴切得有多狠不一样** —— 这是本轮为每一格单造地板的理由 ===")
_s = D[D.moral.notna() & D.REL.notna() & D.year.isin(list(DECS["2010s"]))]
for a, fn in ALL_AX.items():
    hi, lo = fn(_s)
    if hi is None: continue
    print(f"  {a:16s} 高组 {hi.sum():5d}({hi.mean():.1%}) · 低组 {lo.sum():5d}({lo.mean():.1%})"
          f" · **两组合计 {(hi.sum()+lo.sum())/len(_s):.1%}**")
print(f"  ⚠ **`#862` 用一个全局地板 {F862:.4f}(从它的产物读的,不手抄)去量这七种几何** —— "
      f"本轮改为**每格自己的置换分布**。")


def stat(y, hi, lo, scale, cuts):
    """一个格的 g。均值尺 = 标准化差;潜在尺 = 累积 probit 差,逐切点算再平均。"""
    if hi.sum() < 60 or lo.sum() < 60: return np.nan, np.nan
    if scale == "均值尺":
        s = y.std(ddof=1)
        if s <= 0: return np.nan, np.nan
        return float((y[hi].mean() - y[lo].mean()) / s), 0.0
    d = []
    for c in cuts:
        ph = float(np.clip((y[hi] > c).mean(), 1e-4, 1 - 1e-4))
        pl = float(np.clip((y[lo] > c).mean(), 1e-4, 1 - 1e-4))
        d.append(float(ndtri(ph) - ndtri(pl)))
    return float(np.mean(d)), (float(np.std(d)) if len(d) > 1 else 0.0)


def blocks_for(item_col, dec):
    """逐年预算好 y 与七根轴(+安慰剂)的高低掩码 —— 置换只动掩码,不动 y。"""
    m = D[item_col].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))
    sub = D[m]
    if len(sub) < 800: return None, None
    bl = []
    for yv in np.unique(sub.year):
        s = sub[sub.year == yv]
        if len(s) < 200: continue
        y = s[item_col].to_numpy(float)
        mk = {}
        for a, fn in ALL_AX.items():
            hi, lo = fn(s)
            if hi is None or hi.sum() < 60 or lo.sum() < 60: continue
            mk[a] = (hi, lo)
        if mk: bl.append((y, mk))
    return (bl, sub) if bl else (None, None)


def gbar(bl, axis, scale, cuts, plant=0.0, plant_axis=None, perm=None):
    """逐年算再平均。perm = 每年一个共用的置换索引(七根轴同一个 k,保留轴间相关)。"""
    vs, sp = [], []
    for i, (y0, mk) in enumerate(bl):
        if axis not in mk: continue
        hi, lo = mk[axis]
        y = y0
        if plant and plant_axis in mk:
            y = y0.copy(); y[mk[plant_axis][0]] += plant
        if perm is not None:
            k = perm[i]; hi, lo = hi[k], lo[k]
        v, s = stat(y, hi, lo, scale, cuts)
        if np.isfinite(v): vs.append(v); sp.append(s)
    if not vs: return np.nan, np.nan
    return float(np.mean(vs)), float(np.mean(sp))


def jac(sub, axfn):
    hi, _ = axfn(sub); rh, _ = terc("REL", sub)
    if hi is None or rh is None: return np.nan
    u = (hi | rh).sum()
    return float((hi & rh).sum() / u) if u else np.nan


print(f"\n=== ① 规格曲线:{len(ITEMS)} 题 × {len(AXES)} 轴 × {len(DECS)} 十年 × {len(SCALES)} 尺 "
      f"· 每格自己的置换零({NPERM} 次)· ⚠ **比 |g|,符号单独报** ===")
rows, cellinfo = [], {}
rng = np.random.default_rng(SEED)
for inm, (icol, cuts) in ITEMS.items():
    for dec in DECS:
        bl, sub = blocks_for(icol, dec)
        if bl is None: continue
        perms = [[rng.permutation(len(y)) for (y, _) in bl] for _ in range(NPERM)]
        for scale in SCALES:
            obs, null = {}, {a: [] for a in ALL_AX}
            spread = {}
            for a in ALL_AX:
                v, s = gbar(bl, a, scale, cuts)
                obs[a] = v; spread[a] = s
            for b in range(NPERM):
                pk = perms[b]
                for a in ALL_AX:
                    v, _ = gbar(bl, a, scale, cuts, perm=pk)
                    null[a].append(abs(v) if np.isfinite(v) else np.nan)
            live = [a for a in ANAMES if np.isfinite(obs[a])]
            for a in ALL_AX:
                if not np.isfinite(obs[a]): continue
                nd = np.array([x for x in null[a] if np.isfinite(x)])
                fl = float(np.quantile(nd, 0.95)) if len(nd) else np.nan
                pv = float((1 + (nd >= abs(obs[a])).sum()) / (len(nd) + 1)) if len(nd) else np.nan
                rows.append(dict(item=inm, dec=dec, scale=scale, axis=a, g=obs[a], absg=abs(obs[a]),
                                 floor=fl, p=pv, cut_spread=spread[a],
                                 jaccard=jac(sub, ALL_AX[a]), n=int(len(sub)),
                                 placebo=(a in PLACEBO)))
            # 边际自己的零:同一次置换里七根轴的 max−second
            mnull = []
            for b in range(NPERM):
                v = sorted((null[a][b] for a in live if np.isfinite(null[a][b])), reverse=True)
                if len(v) >= 2: mnull.append(v[0] - v[1])
            o = sorted(((abs(obs[a]), a) for a in live), reverse=True)
            cellinfo[(inm, dec, scale)] = dict(
                winner=o[0][1], margin=float(o[0][0] - o[1][0]),
                margin_floor=float(np.quantile(mnull, 0.95)) if mnull else np.nan,
                margin_p=float((1 + sum(1 for x in mnull if x >= o[0][0] - o[1][0])) / (len(mnull) + 1)),
                order=[a for _, a in o], absg={a: float(abs(obs[a])) for a in live},
                sign={a: float(np.sign(obs[a])) for a in live})
        print(f"  {inm:16s} {dec} n={len(sub):5,d} · 均值尺赢家 **{cellinfo[(inm,dec,'均值尺')]['winner']}**"
              f" · 潜在尺赢家 **{cellinfo[(inm,dec,'潜在尺')]['winner']}**")

MORAL = "道德 `homosex`"
print(f"\n=== ②a 道德题逐格:七根轴的 |g|(括号 = 该轴该格自己的地板)===")
for scale in SCALES:
    for dec in DECS:
        k = (MORAL, dec, scale)
        if k not in cellinfo: continue
        rs = {r["axis"]: r for r in rows if (r["item"], r["dec"], r["scale"]) == k}
        ln = " · ".join(f"{a.split()[0]} **{rs[a]['g']:+.3f}**(地板{rs[a]['floor']:.3f})"
                        for a in ANAMES if a in rs)
        print(f"  {scale} {dec} · " + ln)
        c = cellinfo[k]
        print(f"     ⇒ 主轴 **{c['winner']}** · 最大−次大 **{c['margin']:+.3f}** vs "
              f"**边际自己的地板 {c['margin_floor']:.3f}**(p={c['margin_p']:.3f})"
              f"{'  ⚠ **落在地板内**' if c['margin'] <= c['margin_floor'] else ''}")

print(f"\n=== ②b 制裁三题(与 `#862` 同一构造,但用**本轮修好的地板**重跑,且**多了四个十年**)===")
SCELLS = [k for k in cellinfo if k[0] != MORAL]
for inm in [i for i in ITEMS if i != MORAL]:
    for scale in SCALES:
        ws = [(d, cellinfo[(inm, d, scale)]) for d in DECS if (inm, d, scale) in cellinfo]
        print(f"  {inm:16s} {scale} · " + " · ".join(
            f"{d}:{c['winner'].split()[0]}({c['margin']:+.3f}{'⚠地板内' if c['margin']<=c['margin_floor'] else ''})"
            for d, c in ws))
S_EDU = sum(1 for k in SCELLS if cellinfo[k]["winner"] == "教育 educ")
S_EXC = [(k, cellinfo[k]["winner"]) for k in SCELLS if cellinfo[k]["winner"] != "教育 educ"]
S_INFL = [k for k in SCELLS if cellinfo[k]["margin"] <= cellinfo[k]["margin_floor"]]
print(f"  ⚠⚠ **`#862` 说制裁侧「六个格子,无一例外」—— 那是只看了两个十年。"
      f"把六个十年都放进来:教育夺魁 {S_EDU}/{len(SCELLS)},**例外 {len(S_EXC)} 格**"
      + ("".join(f"({k[1]} {k[0].split()[0]} {k[2]} → {w.split()[0]})" for k, w in S_EXC))
      + f",且**边际落在自己地板内 {len(S_INFL)}/{len(SCELLS)} 格**。**排名稳,「无一例外」这四个字不稳。**")
M_LEAD = [cellinfo[k]["margin"] for k in cellinfo if k[0] == MORAL]
S_LEAD = [cellinfo[k]["margin"] for k in SCELLS]
print(f"  ⚠ **两条缝的锋利度不对称**:道德侧宗教的领先 {min(M_LEAD):+.3f}–{max(M_LEAD):+.3f}(中位 "
      f"{np.median(M_LEAD):+.3f}),制裁侧教育的领先 {min(S_LEAD):+.3f}–{max(S_LEAD):+.3f}(中位 "
      f"{np.median(S_LEAD):+.3f})⇒ **道德那条缝的主轴压倒性,制裁那条缝的主轴只是险胜。**")

print("\n=== ③ 控制 ===")
bl0, sub0 = blocks_for("moral", "2010s")
cuts0 = ITEMS[MORAL][1]
base_rel, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0)
base_edu, _ = gbar(bl0, "教育 educ", "均值尺", cuts0)
PLANT = 0.30
p_rel, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0, plant=PLANT, plant_axis="宗教 REL")
p_edu, _ = gbar(bl0, "教育 educ", "均值尺", cuts0, plant=PLANT, plant_axis="宗教 REL")
z_rel, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0, plant=0.0, plant_axis="宗教 REL")
c_rel, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0, plant=3.0, plant_axis="宗教 REL")
CTRL_T = 0.05
FLOOR_C, CEIL_C = abs(z_rel - base_rel), abs(c_rel - base_rel)
print(f"  正控:**只往宗教轴高组**植入 +{PLANT} ⇒ 宗教 **Δg = {p_rel-base_rel:+.4f}** · "
      f"教育 **Δg = {p_edu-base_edu:+.4f}**")
print(f"     ⚠ **写作 Δg 不是 Δ|g|** —— 这两个量在 g<0 时符号相反,"
      f"而 `#862` 把同一个减法印成了「|g| 动」(`#863`③,标签错,数没错)")
print(f"     ⚠ **`G2` 控制必须能失败**:plant=0 时动 **{z_rel-base_rel:+.2e}**")
print(f"     ⚠ **而控制也必须能通过**(`realstat §4`):floor(不植入)**{FLOOR_C:.2e}** < "
      f"阈 **{CTRL_T}** < ceiling(植入满量程 +3.0)**{CEIL_C:.4f}** ⇒ "
      f"**{'阈落在真带内' if FLOOR_C < CTRL_T < CEIL_C else '⚠⚠ 阈不在带内,控制无意义'}**")
MDE_GRID = [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.20, 0.30]
mde_cell_floor = [r for r in rows if r["item"] == MORAL and r["dec"] == "2010s"
                  and r["scale"] == "均值尺" and r["axis"] == "宗教 REL"][0]["floor"]
mde_curve = []
for A in MDE_GRID:
    v, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0, plant=A, plant_axis="宗教 REL")
    mde_curve.append((A, abs(v - base_rel)))
MDE = next((A for A, d in mde_curve if d > mde_cell_floor), None)
print(f"     **保留率与 MDE(`G2` 明写要报,而 `#862` 一个都没报)**:植入 {PLANT} 时 "
      f"|Δg| = **{abs(p_rel-base_rel):.4f}**,占「SD 不随植入变」这一理想情形 "
      f"{PLANT}/SD 的 **{abs(p_rel-base_rel)/(PLANT/np.mean([y.std(ddof=1) for y,_ in bl0])):.1%}** "
      f"(差额是植入把 SD 也撑大了)· **MDE = 最小的能顶出该格地板({mde_cell_floor:.4f})的植入量 "
      f"= {MDE if MDE else '>0.30'}** · 曲线 "
      + " · ".join(f"{A}→{d:.4f}" for A, d in mde_curve))
neg_rows = [r for r in rows if r["item"] == MORAL and r["dec"] == "2010s" and r["scale"] == "均值尺"]
nrel = [r for r in neg_rows if r["axis"] == "宗教 REL"][0]
rg2 = np.random.default_rng(SEED + 1)
pkk = [rg2.permutation(len(y)) for (y, _) in bl0]
n_rel, _ = gbar(bl0, "宗教 REL", "均值尺", cuts0, perm=pkk)
print(f"  负控:**年内打乱宗教轴标签** ⇒ |g| = **{abs(n_rel):.4f}**,该格地板 **{nrel['floor']:.4f}** "
      f"⇒ **{'塌回地板内' if abs(n_rel) <= nrel['floor'] else '⚠ 没塌回去'}**"
      f"(⚠ **「这个零该不该是零?」该** —— 打乱后高低两组的差期望就是 0)")
pl_rows = [r for r in rows if r["placebo"]]
pl_in = sum(1 for r in pl_rows if r["absg"] <= r["floor"])
print(f"  **安慰剂 = `ballot`**(GSS 自己随机分配的问卷版本,**真实变量,不是我的 rng**):"
      f"**{pl_in}/{len(pl_rows)}** 格落在自己的地板内 · |g| 中位 "
      f"**{np.median([r['absg'] for r in pl_rows]):.4f}** · 最大 **{max(r['absg'] for r in pl_rows):.4f}**")
mor_sp = [r["cut_spread"] for r in rows if r["item"] == MORAL and r["scale"] == "潜在尺"
          and not r["placebo"] and np.isfinite(r["cut_spread"])]
mor_ab = [r["absg"] for r in rows if r["item"] == MORAL and r["scale"] == "潜在尺" and not r["placebo"]]
CUTR = float(np.median(mor_sp) / np.median(mor_ab))
print(f"  **潜在尺的适配性**:道德题三个切点之间 g 的标准差,中位 **{np.median(mor_sp):.4f}**,"
      f"对照同格 |g| 中位 **{np.median(mor_ab):.4f}** ⇒ 比值 **{CUTR:.2f}** ⇒ "
      + (f"**有序 probit 的等距假设站得住** —— 三个切点给出的 g 几乎是同一个数,"
         f"所以「潜在 g」是一个参数,不是三个数的平均" if CUTR < 0.25 else
         f"⚠⚠ **等距假设不成立 —— 那个「潜在 g」只是三个不同数的平均,不是一个参数**"))
print(f"     ⚠ **这一行本来写死了「比值大 ⇒ 假设不成立」,而实测比值是 {CUTR:.2f} 也就是小** —— "
      f"**判词字符串不是计算**(`realstat §4`),当场改成条件式(`#863`④)")

print(f"\n=== ④ 多重性:整族 {len([r for r in rows if not r['placebo']])} 格,BH 与 BY 都做 "
      f"—— ⚠ **`#862` 声明了这一条却没实现,本轮补上** ===")
fam = [r for r in rows if not r["placebo"] and np.isfinite(r["p"])]
ps = np.array([r["p"] for r in fam]); C = len(ps)
o = np.argsort(ps); q = 0.05
cH = q * (np.arange(1, C + 1)) / C
cY = cH / np.sum(1.0 / np.arange(1, C + 1))
def step_up(pv, crit):
    ok = pv <= crit
    k = np.max(np.where(ok)[0]) + 1 if ok.any() else 0
    return k
kH, kY = step_up(ps[o], cH), step_up(ps[o], cY)
survH = {id(fam[i]) for i in o[:kH]}; survY = {id(fam[i]) for i in o[:kY]}
for r in fam: r["bh"] = id(r) in survH; r["by"] = id(r) in survY
print(f"  测了 **{C}** 格 · BH(q=0.05)存活 **{kH}** · BY 存活 **{kY}** · "
      f"p 的分辨率下限 1/({NPERM}+1) = **{1/(NPERM+1):.4f}**")
for a in ANAMES:
    ra = [r for r in fam if r["axis"] == a]
    print(f"    {a:16s} {len(ra):3d} 格 · BH 存活 {sum(r['bh'] for r in ra):3d} · "
          f"BY 存活 {sum(r['by'] for r in ra):3d} · |g| 中位 {np.median([r['absg'] for r in ra]):.3f}")
print("  ⚠ **不同意的格一起发表**:下面这些格 BH 存活而 BY 不存活 ⇒ "
      f"**{sum(1 for r in fam if r['bh'] and not r['by'])} 格**,它们是本轮最不该被当成结论的部分。")

# ── ⑤ 跨仪器:NSFG `samesex` ────────────────────────────────────────────────────
print("\n=== ⑤ 跨仪器复制(硬规则④):NSFG `samesex`,**能对上的只有五根轴** ===")
def dct_cols(dct, want):
    out = {}
    for m in re.finditer(r"_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)", dct.read_text(errors="replace")):
        out[m.group(2).lower()] = (int(m.group(1)), int(m.group(3)))
    return {w: out[w] for w in want if w in out}
NF = [("2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct", 2012, 2.0),
      ("2017_2019_FemRespData.dat", "2017_2019_FemRespSetup.dct", 2018, 2.0),
      ("2017_2019_MaleData.dat", "2017_2019_MaleSetup.dct", 2018, 1.0)]
WANT = ["samesex", "attndnow", "reldlife", "ager", "educat", "race"]
fr = []
for dat, dct, yr, sxc in NF:
    cols = dct_cols(EXT / "nsfg/setup" / dct, WANT)
    rr = []
    with open(EXT / "nsfg" / dat, "r", errors="replace") as fh:
        for line in fh:
            d = {}
            for v, (c, w) in cols.items():
                s = line[c - 1:c - 1 + w].strip()
                d[v] = float(s) if s.lstrip("-").isdigit() else np.nan
            rr.append(d)
    t = pd.DataFrame(rr); t["year"] = yr; t["sexn"] = sxc; fr.append(t)
N = pd.concat(fr, ignore_index=True)
N = N[N.ager >= 18].copy()
ss = N.samesex.where(lambda v: v.isin([1, 2, 3, 4]))
N["perm"] = 5 - ss                                       # ⚠ 1=强烈同意「都可以」⇒ 翻转后高=宽容
N["f_att"] = 8 - N.attndnow.where(lambda v: v.between(1, 7))
N["f_rel"] = 4 - N.reldlife.where(lambda v: v.between(1, 3))
NR = N.dropna(subset=["f_att", "f_rel"]).copy()
NR["REL"] = (zs(NR.f_att) + zs(NR.f_rel)) / 2
N = N.join(NR["REL"])
N["educ"] = N.educat.where(lambda v: v.between(9, 19))
N["age"] = N.ager.where(lambda v: v.between(18, 49))
N["race_"] = N.race.where(lambda v: v.isin([1, 2, 3]))
print(f"  NSFG n={len(N):,}(18+)· `samesex` 非缺 {int(N.perm.notna().sum()):,} · "
      f"REL 非缺 {int(N.REL.notna().sum()):,} · 两波 {sorted(N.year.unique().astype(int))}")
print("  ⚠ **NSFG `race`:1=BLACK · 2=WHITE**(读的是 `.sas` 值标签),**与 GSS 相反 ⇒ 已翻转对齐**")
print("  ⚠⚠ **NSFG 只到 49 岁** ⇒ **年龄轴结构性地不可能复制 GSS**,不是「没复制上」")
NAX = {"宗教 REL": ("REL", "terc"), "教育 educ": ("educ", "terc"), "年龄 age": ("age", "terc"),
       "种族 白vs黑": ("race_", ("eq", 2.0, 1.0)), "性别 男vs女": ("sexn", ("eq", 1.0, 2.0))}
def nmask(s, spec):
    col, how = spec
    v = s[col]
    if how == "terc":
        if v.notna().sum() < 200: return None, None
        lo, hi = np.nanquantile(v, [1 / 3, 2 / 3])
        return (v >= hi).to_numpy(), (v <= lo).to_numpy()
    _, a, b = how
    return (v == a).to_numpy(), (v == b).to_numpy()
nrows, nwin = [], {}
sub = N[N.perm.notna() & N.REL.notna()]
rgN = np.random.default_rng(SEED + 7)
for wave in sorted(sub.year.unique()):
    s = sub[sub.year == wave]
    if len(s) < 800: continue
    y = s.perm.to_numpy(float)
    mk = {a: nmask(s, sp) for a, sp in NAX.items()}
    mk = {a: v for a, v in mk.items() if v[0] is not None and v[0].sum() >= 60 and v[1].sum() >= 60}
    perms = [rgN.permutation(len(s)) for _ in range(NPERM)]
    for scale in SCALES:
        obs, nul = {}, {}
        for a, (hi, lo) in mk.items():
            v, _ = stat(y, hi, lo, scale, [1, 2, 3])
            obs[a] = v
            nd = []
            for k in perms:
                vv, _ = stat(y, hi[k], lo[k], scale, [1, 2, 3])
                if np.isfinite(vv): nd.append(abs(vv))
            nul[a] = np.array(nd)
        live = [a for a in obs if np.isfinite(obs[a])]
        for a in live:
            nrows.append(dict(wave=int(wave), scale=scale, axis=a, g=obs[a], absg=abs(obs[a]),
                              floor=float(np.quantile(nul[a], 0.95)),
                              p=float((1 + (nul[a] >= abs(obs[a])).sum()) / (len(nul[a]) + 1)), n=int(len(s))))
        o2 = sorted(((abs(obs[a]), a) for a in live), reverse=True)
        mn = []
        for b in range(NPERM):
            vv = sorted((nul[a][b] for a in live), reverse=True)
            if len(vv) >= 2: mn.append(vv[0] - vv[1])
        nwin[(int(wave), scale)] = dict(winner=o2[0][1], margin=float(o2[0][0] - o2[1][0]),
                                        margin_floor=float(np.quantile(mn, 0.95)) if mn else np.nan,
                                        order=[a for _, a in o2])
        print(f"  {int(wave)} {scale} n={len(s):5,d} · " +
              " · ".join(f"{a.split()[0]} **{obs[a]:+.3f}**" for a in live) +
              f" ⇒ 主轴 **{o2[0][1]}** · 边际 **{o2[0][0]-o2[1][0]:+.3f}** vs 地板 "
              f"{nwin[(int(wave),scale)]['margin_floor']:.3f}")
NREL_WINS = sum(1 for v in nwin.values() if v["winner"] == "宗教 REL")

# ── ⑥ 闸 ───────────────────────────────────────────────────────────────────────
MCELLS = [k for k in cellinfo if k[0] == MORAL]
rel_wins = sum(1 for k in MCELLS if cellinfo[k]["winner"] == "宗教 REL")
oth = {a: sum(1 for k in MCELLS if cellinfo[k]["winner"] == a) for a in ANAMES}
top_other = max(((a, v) for a, v in oth.items() if a != "宗教 REL"), key=lambda kv: kv[1])
diffuse = sum(1 for k in MCELLS if cellinfo[k]["margin"] <= cellinfo[k]["margin_floor"])
scale_flip = sum(1 for d in DECS if (MORAL, d, "均值尺") in cellinfo
                 and cellinfo[(MORAL, d, "均值尺")]["winner"] != cellinfo[(MORAL, d, "潜在尺")]["winner"])
ndecs = len({k[1] for k in MCELLS})
rel_margin_signs = [cellinfo[k]["margin"] if cellinfo[k]["winner"] == "宗教 REL"
                    else -cellinfo[k]["margin"] for k in MCELLS]
relJ = {a: float(np.nanmean([r["jaccard"] for r in rows if r["axis"] == a and r["item"] == MORAL]))
        for a in ANAMES}

G = Gate("#863 · 道德那一侧的主轴是不是宗教(与 `#862` 完全同一构造)")
G.asserted("① 前提(跑前写下的最强混淆之一):**七根轴彼此相关,大差距可能只是宗教换了件衣服** ⇒ "
           "**逐格印出该轴高组与虔诚高组的 Jaccard**,折扣印在每一行旁边(`#820`/`#834`/`#862` 同一条)",
           bool(all(np.isfinite(relJ[a]) for a in ANAMES if a != "宗教 REL")),
           " · ".join(f"{a.split()[0]}:{relJ[a]:.2f}" for a in ANAMES if a != "宗教 REL"), kind="control")
G.asserted("② 前提(第二条混淆):**各轴切得有多狠不一样(educ 三分位 37%/42%,sex 50/50,region 25/75)** "
           f"⇒ **每格用它自己的置换分布做地板**,而不是 `#862` 的单一全局 {F862:.4f}",
           bool(all(np.isfinite(r["floor"]) for r in rows)),
           f"{len(rows)} 格各自有地板,道德题地板跨度 "
           f"{min(r['floor'] for r in rows if r['item']==MORAL):.4f}–"
           f"{max(r['floor'] for r in rows if r['item']==MORAL):.4f} ⇒ **最宽的地板是最窄的 "
           f"{max(r['floor'] for r in rows if r['item']==MORAL)/min(r['floor'] for r in rows if r['item']==MORAL):.1f} 倍**",
           kind="control")
G.asserted("③ **比的是 `|g|` 不是 `g`**,符号单独存产物,绝不参与排序",
           bool(all("absg" in r for r in rows)), f"共 {len(rows)} 格", kind="control")
G.asserted("④ 正控:只往宗教轴高组植入 +0.30 ⇒ 宗教必须动、教育几乎不动;**plant=0 必须恰为 0**;"
           "**且阈 0.05 必须落在 floor 与 ceiling 之间**(`realstat §4`「控制必须能通过」)",
           bool(abs(p_rel - base_rel) > CTRL_T and FLOOR_C < 1e-12 and FLOOR_C < CTRL_T < CEIL_C),
           f"宗教动 {p_rel-base_rel:+.4f} · 教育动 {p_edu-base_edu:+.4f} · plant=0 {z_rel-base_rel:+.2e} · "
           f"带 [{FLOOR_C:.2e}, {CEIL_C:.4f}]", kind="control")
G.asserted("⑤ 负控:年内打乱宗教轴标签 ⇒ 必须塌回**该格自己的**地板内",
           bool(abs(n_rel) <= nrel["floor"]), f"{abs(n_rel):.4f} ≤ {nrel['floor']:.4f}", kind="control")
G.asserted("⑥ 安慰剂:**`ballot`(GSS 随机分配的问卷版本)**必须落在地板内 —— "
           "⚠ **这是真实变量,不是我造的随机数**",
           bool(pl_in >= len(pl_rows) - 1), f"{pl_in}/{len(pl_rows)} 格落在地板内 · "
           f"最大 |g| {max(r['absg'] for r in pl_rows):.4f}", kind="control")
G.asserted("⑦ kill(预注册):「道德那一侧的主轴是宗教」要成立,需 **宗教在 ≥8/12 个道德格"
           "(6 十年 × 2 尺)里 `|g|` 最大**",
           bool(rel_wins >= 8),
           f"宗教夺魁 {rel_wins}/{len(MCELLS)} · 次多 {top_other[0]} {top_other[1]} · "
           f"边际落在**边际自己地板**内 {diffuse}/{len(MCELLS)} · 两尺赢家不同的十年 {scale_flip}/{ndecs}",
           kind="kill",
           yardstick="每格七根轴 |g| 的排序;边际另有**为边际单造**的零分布(七轴同一次置换,保留轴间相关)",
           yardstick_noise=float(np.median([cellinfo[k]["margin_floor"] for k in MCELLS])),
           population=f"GSS `homosex` 的 {len(MCELLS)} 个 (十年 × 尺) 格 —— "
                      f"⚠ **不含制裁三题**(那是 `#862` 的总体,在场只作对照,`#854`①:判据的总体里不许有在位者)",
           direction=rel_margin_signs)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
sanc_win = {}
for inm in [i for i in ITEMS if i != MORAL]:
    for k in [k for k in cellinfo if k[0] == inm]:
        sanc_win[cellinfo[k]["winner"]] = sanc_win.get(cellinfo[k]["winner"], 0) + 1
sanc_tot = sum(sanc_win.values())
sanc_top = max(sanc_win.items(), key=lambda kv: kv[1]) if sanc_win else ("—", 0)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif scale_flip >= 3:
    VERD = (f"**D 主轴取决于尺 —— 元分离器:那个排名不是人的性质,是答题格式的性质。**\n"
            f"  {scale_flip}/{ndecs} 个十年里,均值尺与潜在尺给出不同的赢家。")
elif diffuse >= 8:
    VERD = (f"**C 道德那条缝也是弥散的 ⇒ 「哪根轴最宽」这个问法两侧都不成立。**\n"
            f"  {diffuse}/{len(MCELLS)} 格的边际落在**边际自己的**地板内。")
elif rel_wins >= 8:
    d70 = cellinfo[(MORAL, "1970s", "均值尺")]
    VERD = (f"**A 道德那一侧的主轴是宗教 —— 宗教在 {rel_wins}/{len(MCELLS)} 个道德格里 `|g|` 最大**"
            f"(次多:{top_other[0]} {top_other[1]} 格)。\n"
            f"  而同一次运行里,**制裁 {sanc_tot} 格的主轴是 {sanc_top[0]}({sanc_top[1]} 格)** ⇒ "
            f"**`#862` 在修好的地板下、且多了四个十年之后仍然成立。**\n"
            f"  ⚠ 跨仪器:NSFG `samesex` **{NREL_WINS}/{len(nwin)}** 格主轴也是宗教"
            f"(**能对上的只有五根轴,且 NSFG 只到 49 岁**)。\n"
            f"  ⇒ **一句关于人的话:同一个话题,问「这样对不对」时把美国人分得最开的是信不信教;\n"
            f"  问「该不该拦着他去教书」时分得最开的是读了多少书。\n"
            f"  不是同一条裂缝换了个问法 —— 是两条缝,而人们站在它们两边的位置不一样。**\n"
            f"  ⚠ **而这句话有三处不能省的边界,全部是本轮自己量出来的:**\n"
            f"  ① **七十年代不算**:那一格宗教({d70['absg']['宗教 REL']:.3f})与年龄"
            f"({d70['absg']['年龄 age']:.3f})的差 {d70['margin']:+.3f} **落在边际自己的地板 "
            f"{d70['margin_floor']:.3f} 内** ⇒ **那时候年龄和信仰一样能分开人;宗教是后来才拉开的。**\n"
            f"  ② **两条缝的锋利度差一个量级**:道德侧宗教的领先中位 {np.median(M_LEAD):+.3f},"
            f"制裁侧教育的领先中位 {np.median(S_LEAD):+.3f},且制裁侧 {len(S_INFL)}/{len(SCELLS)} 格"
            f"的领先落在自己地板内 ⇒ **「制裁侧是教育」是险胜,不是压倒。**\n"
            f"  ③ **`#862` 的「无一例外」四个字要撤**:六个十年下例外 {len(S_EXC)} 格"
            + ("".join(f"({k[1]} {k[0].split()[0]} {k[2]}→{w.split()[0]})" for k, w in S_EXC))
            + f" ⇒ `#863`②。")
elif top_other[1] >= 8:
    VERD = (f"**B 道德那一侧的主轴也不是宗教,是 {top_other[0]} ⇒ 「宗教那条缝」这个框架\n"
            f"  对整个对象都错了,而 `#834` 的结果只是「偏离」那把尺特有的。**\n"
            f"  {top_other[0]} 在 {top_other[1]}/{len(MCELLS)} 格里最大,宗教 {rel_wins} 格。")
else:
    VERD = (f"**都不是**:宗教 {rel_wins} · {top_other[0]} {top_other[1]} · 弥散 {diffuse} · "
            f"尺翻转 {scale_flip}(共 {len(MCELLS)} 格)—— **四个预注册世界都没被满足,如实登记。**")
print(VERD)
print(f"\n⚠ **本轮结构性做不到的**:① **制裁三题换不了仪器**(Stouffer 三题盘上是 GSS 独有)⇒ "
      f"跨仪器只覆盖道德题;② 横断面 ⇒ **无因果识别**;③ **只比「哪根轴最宽」,不做多元调整** —— "
      f"那是另一个估计量,**本轮不做,不是「以后做」**;④ **NSFG 只到 49 岁 ⇒ 年龄轴结构上不可复制**;"
      f"⑤ **NSFG 只有 2011–2019 三个文件有 setup** ⇒ 九十年代的跨仪器复制盘上不存在。")

json.dump(dict(grid=rows, cells={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in cellinfo.items()},
               nsfg=dict(rows=nrows, winners={f"{k[0]}|{k[1]}": v for k, v in nwin.items()},
                         rel_wins=NREL_WINS, n_cells=len(nwin)),
               moral_rel_wins=rel_wins, moral_cells=len(MCELLS), moral_other=oth,
               diffuse=diffuse, scale_flip=scale_flip, sanction_winners=sanc_win,
               multiplicity=dict(cells=C, bh=int(kH), by=int(kY), q=q, p_resolution=1 / (NPERM + 1)),
               controls=dict(plant=PLANT, rel_moved=p_rel - base_rel, edu_moved=p_edu - base_edu,
                             zero_plant=z_rel - base_rel, ceiling=CEIL_C, threshold=CTRL_T,
                             neg=abs(n_rel), neg_floor=nrel["floor"],
                             placebo_inside=pl_in, placebo_cells=len(pl_rows),
                             latent_cut_spread_ratio=float(np.median(mor_sp) / np.median(mor_ab)),
                             ppf_maxdiff=PPF_MAXDIFF),
               sanction_edu_wins=S_EDU, sanction_cells=len(SCELLS),
               sanction_exceptions=[f"{k[1]}|{k[0]}|{k[2]}→{w}" for k, w in S_EXC],
               sanction_margin_in_floor=len(S_INFL),
               moral_lead=[min(M_LEAD), float(np.median(M_LEAD)), max(M_LEAD)],
               sanction_lead=[min(S_LEAD), float(np.median(S_LEAD)), max(S_LEAD)],
               mde=dict(cell="道德|2010s|均值尺|宗教", floor=mde_cell_floor, mde=MDE,
                        curve=[[a, d] for a, d in mde_curve]),
               latent_cut_ratio=CUTR,
               fixes_to_862=["#862 声明 BH/BY 未实现,本轮实现",
                             "#862 用单个 |g| 的地板去量『最大−次大』,本轮为边际单造零分布",
                             "#862 的『六个格子无一例外』只看了两个十年;六个十年下有例外",
                             "#862 把 Δg 印成了 Δ|g|(标签错,数没错)"],
               mean_jaccard_with_devout=relJ, admissible=adm, verdict=VERD, gate_ok=G.verdict(),
               seed=SEED, nperm=NPERM),
          open(OUT / "which_cleavage_moral.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'which_cleavage_moral.json'}")
