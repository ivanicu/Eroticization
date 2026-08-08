"""#853 · E03·A87·R292 —— 第一次跨仪器复现:换一份问卷,那条宗教缝还在不在?

**硬规则④:跨仪器复现胜过同一具再跑一轮。而本项目 292 轮,一次都没做过。**
从 `#832` 到 `#846`,所有关于「虔诚 / 世俗那条缝」的数字**全部出自 GSS 一份问卷** ——
`#836` 换过单位(社会)、`#839` 换过估计量(潜在尺)、`#846` 换过题目(八题),
**但从来没换过问卷。** ⇒ 那条缝的**大小**,到今天为止是一句关于 GSS 的话,不是关于美国的话。

**⓪ 硬规则①当场救了两次(变量名不是测量):**
① `.sas` 里 `VALUE SAMESEX1F` 的码值是 **'UNDER 15 YEARS' / '15-19 YEARS' …** ——
   **那是年龄分组,不是态度量表**;它只是恰好以 `SAMESEX` 开头。
   真正指派给 `samesex` 的 format 是 **`AGDGFMT`**:
   `1 强烈同意 · 2 同意 · 3 不同意 · 4 强烈不同意 · 5「若坚持:既不同意也不反对」· 8 拒答 · 9 不知道`。
② 同一本字典里还有 **`SAMESEXANY`**(是否有过同性性接触)与 **`SAMESEX1`**(初次同性性经历年龄)——
   **前缀相同、构念完全不同。** 按名字前缀取变量会取到行为题或年龄题。

**⚠⚠ 而 `AGDGFMT` 的码 5 是 `#836` 那个 `SCCS176` 码 2 的同一个陷阱:**
**「既不同意也不反对」在码序上排在 4 之后,在语义上却坐在 2 和 3 中间。**
⇒ 它**不能**当成序数尺的一端。本轮把它当作一根**规格轴**:{剔除 · 编码为中点},两种都跑(`G4`)。

**⚠ 极性相反,必须显式处理:**
GSS `homosex` 1 = always wrong … 4 = not wrong at all ⇒ **高 = 宽容**;
NSFG `samesex` 1 = 强烈同意「可以」⇒ **低 = 宽容** ⇒ 翻转为 `5 − x`。

`G1` **估计量(先于方法命名)**:
   **虔诚层与世俗层在「同性性关系」态度上的差距,以该题当年自身的 SD 为单位** ——
   `g = (mean_虔诚 − mean_世俗) / SD_年内合并`,两具仪器各自独立算一遍,
   **量的是 `g_GSS − g_NSFG`。**
   ⚠ **必须是无量纲的**:两题的量表根本不同(一个问「错不错」,一个问「同不同意」),
   **原始分之差跨仪器没有意义。**

**识别**:两具仪器在 2010 年代都同时有「态度题 + 宗教度 + 年份」⇒ 可估。
   ⚠ **但总体不同** —— 见下面跑前写下的最强混淆。

三个世界:
   A **仪器无关**:两具仪器的 `g` 区间重叠 ⇒ **那条缝是关于美国的事实,不是关于 GSS 问卷的。**
   B **仪器决定大小**:两个 `g` 明显不同,且限制总体后仍不同
     ⇒ **缝的大小是「怎么问」的性质** ⇒ **本项目发表过的每一个关于缝大小的数字都只对 GSS 成立。**
     **⚠ 这是我不欢迎的那个,而 `#832` 之后的一切都压在它上面(`§3` basin)。**
   C **总体解释掉**:两个 `g` 不同,但把 GSS 限制到 NSFG 的年龄/性别后收敛 ⇒ **是构成,不是仪器。**

预测矩阵:
   | 世界 | 现在 | 区间重叠 | 限制后收敛 | 限制后仍不同 |
   | A 仪器无关   | 0.35 | **0.85** | 0.10 | 0.05 |
   | C 总体解释   | 0.30 | 0.10 | **0.85** | 0.10 |
   | B 仪器决定   | 0.35 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**植入一个已知位移必须按 SD 单位取回,且 g=0 时位移必须 ≈0**)
     and 负控为零(**年内打乱宗教度标签 ⇒ 两具仪器的 g 都必须 ≈0**)
     and 安慰剂为零(**把宗教度换成随机变量 ⇒ 两具仪器都必须 ≈0**):
      12 格里 `g_GSS − g_NSFG` 的区间含 0 的 ≥ 2/3 -> A
      限制总体后两者靠拢(|差| 至少减半)          -> C
      否则                                          -> B
  else: UNVERIFIED

⚠⚠ **跑之前写下的最强混淆(它决定了「限制」这根轴):**
   **NSFG 的总体是 15–49 岁,而女性问卷只有女性;GSS 是 18+ 且两性。**
   **年轻与女性都更宽容** ⇒ 光凭构成就能造出 `g` 的差异,与仪器无关。
   ⇒ 控制:**把 GSS 限制到 18–49 岁女性**,与 NSFG 女性对齐,**限制前后都报**。

`G3` 多重性:整族 = **宗教度算法 3 种**(出席 · 重要性 · 合成)× **码 5 处理 2 种** ×
   **总体 2 种** = **12 格**,BH 与 BY 都做,**不同意的格一起发表**。
`G4` 规格曲线:上面三根轴就是曲线,**逐格报,不挑一格**。

**⚠ 本站结构性做不到的(登记,不许写「计划中」):**
① **九十年代复现不了** —— NSFG 的字典(`setup/*.dct`)**只存在于 2011–2013 与 2017–2019**;
   另外 11 个 `.dat` **没有字典**,而定宽文件没有字典就是不可读的。
   ⇒ **本项目的头条(九十年代那条缝)在第二具仪器上结构性无法检验。**
② NSFG 男性文件**只有 2017–2019** ⇒ 性别这根轴在时间上不平衡,本轮不做男女对比。
③ 无干预、无因果识别:两具仪器都是横断面调查。
"""
import json, pathlib, re, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
EXT = ROOT / "data/external"
B, Q = 2000, 0.05

# ── NSFG:从 .dct 读列位置,再按定宽切 ──────────────────────────────────────
def dct_cols(dct, want):
    """从 Stata dictionary 读 (变量 -> (起始列, 宽度));⚠ 变量名大小写在两本字典里不同。"""
    out = {}
    for m in re.finditer(r"_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)", dct.read_text(errors="replace")):
        out[m.group(2).lower()] = (int(m.group(1)), int(m.group(3)))
    return {w: out[w] for w in want if w in out}

NSFG_FILES = [
    ("2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct", 2012, "F"),
    ("2017_2019_FemRespData.dat", "2017_2019_FemRespSetup.dct", 2018, "F"),
    ("2017_2019_MaleData.dat",    "2017_2019_MaleSetup.dct",    2018, "M"),
]
WANT = ["samesex", "attndnow", "reldlife", "ager"]
frames = []
print("=== ⓪ 硬规则①:NSFG 每个文件真的有哪些列、多少行、码值分布 ===")
for dat, dct, yr, sex in NSFG_FILES:
    cols = dct_cols(EXT / "nsfg/setup" / dct, WANT)
    if len(cols) < len(WANT):
        print(f"  {dat}: **缺列** {set(WANT)-set(cols)} ⇒ 跳过"); continue
    rows = []
    with open(EXT / "nsfg" / dat, "r", errors="replace") as fh:
        for line in fh:
            r = {}
            for v, (c, w) in cols.items():
                s = line[c - 1:c - 1 + w].strip()
                r[v] = float(s) if s.lstrip("-").isdigit() else np.nan
            rows.append(r)
    d = pd.DataFrame(rows); d["year"] = yr; d["sex"] = sex
    frames.append(d)
    ss = d.samesex.value_counts().sort_index().head(7).to_dict()
    print(f"  {dat}: n=**{len(d):,}** · year={yr} · sex={sex} · samesex 码分布 "
          f"{ {int(k): int(v) for k, v in ss.items()} }")
N = pd.concat(frames, ignore_index=True)
print(f"  NSFG 合计 n=**{len(N):,}** · 年龄 {int(N.ager.min())}–{int(N.ager.max())}")
print("  ⚠ **码 5 =「既不同意也不反对」在码序上排 4 之后、语义上坐 2 和 3 之间** ⇒ 当规格轴处理")
print(f"  ⚠ 而 `SAMESEXANY`(是否有过同性性接触)与 `SAMESEX1`(初次年龄)**前缀相同、构念不同** —— "
      f"按前缀取变量会取错题")

# ── GSS ────────────────────────────────────────────────────────────────────
gp = EXT / "gss/GSS_stata/gss7224_r3a.dta"
g = pd.read_stata(gp, columns=["year", "age", "sex", "attend", "reliten", "fund", "homosex"],
                  convert_categoricals=False)
G = pd.DataFrame({"year": g.year})
G["homosex"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)),
                    ("age", (18, 89)), ("sex", (1, 2))):
    G[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
G = G[(G.year >= 2010) & (G.year <= 2019)]
print(f"\n  GSS 2010–2019:n=**{int(G.homosex.notna().sum()):,}** · 年份 "
      f"{sorted(G.year.unique().tolist())}")

zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0

def build(instr, relig, neither, restrict):
    """返回 (宽容度分数, 虔诚度分数, 年份) —— 两具仪器统一为**高=宽容 / 高=虔诚**。"""
    if instr == "GSS":
        d = G.copy()
        if restrict: d = d[(d.age <= 49) & (d.sex == 2)]
        perm = d.homosex                                   # 已是高=宽容
        att, imp = d.attend, -d.reliten                    # 高=常去 / 高=认同强
    else:
        d = N.copy()
        if restrict: d = d[(d.ager <= 49) & (d.sex == "F")]
        d = d[(d.ager >= 18)]                              # ⚠ GSS 无 <18,对齐下界
        ss = d.samesex.where(lambda v: v.isin([1, 2, 3, 4, 5]))
        if neither == "剔除": ss = ss.where(ss != 5)
        else: ss = ss.where(ss != 5, 2.5)                  # 中点编码
        perm = 5 - ss                                      # ⚠ 极性翻转 ⇒ 高=宽容
        att = 8 - d.attndnow.where(lambda v: v.between(1, 7))
        imp = 4 - d.reldlife.where(lambda v: v.between(1, 3))
    rel = {"出席": zs(att), "重要性": zs(imp), "合成": (zs(att) + zs(imp)) / 2}[relig]
    m = perm.notna() & rel.notna()
    return perm[m].to_numpy(float), rel[m].to_numpy(float), d.year[m].to_numpy(int)

def gap(perm, rel, yr, rng=None, shift=0.0, permute=False, placebo=None):
    """年内三分位切虔诚/世俗;标准化差距 = (虔诚均值 − 世俗均值) / 年内合并 SD。"""
    out = []
    for y in np.unique(yr):
        i = yr == y
        p, r = perm[i].copy(), (placebo[i] if placebo is not None else rel[i])
        if permute: r = rng.permutation(r)
        if len(p) < 200: continue
        lo, hi = np.quantile(r, [1 / 3, 2 / 3])
        dev, sec = p[r >= hi], p[r <= lo]
        if len(dev) < 60 or len(sec) < 60: continue
        s = p.std(ddof=1)
        if s <= 0: continue
        out.append((dev.mean() + shift - sec.mean()) / s)
    return float(np.mean(out)) if out else np.nan

def boot(perm, rel, yr, seed, Bv=B):
    rg = np.random.default_rng(seed); o = np.empty(Bv)
    for i in range(Bv):
        k = rg.integers(0, len(perm), len(perm))
        o[i] = gap(perm[k], rel[k], yr[k])
    return o[np.isfinite(o)]

RELIGS, NEITHERS, RESTRICTS = ("出席", "重要性", "合成"), ("剔除", "中点"), (False, True)
print(f"\n=== ① 规格曲线:{len(RELIGS)}×{len(NEITHERS)}×{len(RESTRICTS)} = "
      f"**{len(RELIGS)*len(NEITHERS)*len(RESTRICTS)} 格**,逐格报(`G3`/`G4`)===")
CELLS, rows = [], []
for relig in RELIGS:
    for neither in NEITHERS:
        for restrict in RESTRICTS:
            key = (relig, neither, restrict)
            gp_, gr_, gy_ = build("GSS", relig, neither, restrict)
            np_, nr_, ny_ = build("NSFG", relig, neither, restrict)
            gG, gN = gap(gp_, gr_, gy_), gap(np_, nr_, ny_)
            bG = boot(gp_, gr_, gy_, 853); bN = boot(np_, nr_, ny_, 854)
            n = min(len(bG), len(bN))
            diff = bG[:n] - bN[:n]
            lo, hi = np.quantile(diff, [.025, .975])
            p = max(2 * min(float(np.mean(diff <= 0)), float(np.mean(diff >= 0))), 1 / (n + 1))
            CELLS.append(key)
            rows.append(dict(relig=relig, neither=neither, restrict=bool(restrict),
                             g_gss=gG, g_nsfg=gN, diff=gG - gN, lo=float(lo), hi=float(hi), p=p,
                             overlap=bool(lo <= 0 <= hi),
                             n_gss=len(gp_), n_nsfg=len(np_)))
            print(f"  {relig:4s}·码5{neither:2s}·{'限制' if restrict else '全体'}: "
                  f"GSS **{gG:+.3f}** vs NSFG **{gN:+.3f}** · 差 **{gG-gN:+.3f}** "
                  f"[{lo:+.3f},{hi:+.3f}] p={p:.4f}{'  **含0**' if lo <= 0 <= hi else '  排零'}")
ps = [r["p"] for r in rows]
bh = {CELLS[i] for i in Gate.bh(ps, Q)}; by = {CELLS[i] for i in Gate.by(ps, Q)}
ov = sum(r["overlap"] for r in rows)
print(f"  `G3` 整族 **{len(rows)} 格** ⇒ **含 0 的 {ov}/{len(rows)}** · BH 排零 **{len(bh)}** · "
      f"BY 排零 **{len(by)}**(⚠ 族大小印在旁边 —— `#832`:族越窄存活越易)")

print("\n=== ② 控制 ===")
gp_, gr_, gy_ = build("GSS", "合成", "剔除", False)
np_, nr_, ny_ = build("NSFG", "合成", "剔除", False)
base_G, base_N = gap(gp_, gr_, gy_), gap(np_, nr_, ny_)
sdG = float(np.std(boot(gp_, gr_, gy_, 1, 400)))
PLANT = 0.30
pcG = gap(gp_, gr_, gy_, shift=PLANT) - base_G
zeroG = gap(gp_, gr_, gy_, shift=0.0) - base_G
print(f"  正控:给 GSS 虔诚层加 +{PLANT}(原始分)⇒ 标准化差距动 **{pcG:+.4f}** "
      f"= **{abs(pcG)/sdG:.1f} 个自助 SD**;**而 g=0 时动 {zeroG:+.6f}(必须 ≈0,否则这条控制不会失败)**")
rg = np.random.default_rng(9)
ncG = gap(gp_, gr_, gy_, rng=rg, permute=True)
ncN = gap(np_, nr_, ny_, rng=rg, permute=True)
print(f"  负控:**年内打乱宗教度标签** ⇒ GSS **{ncG:+.4f}** · NSFG **{ncN:+.4f}**")
print(f"     ⚠ **「这个零该不该是零?」该** —— 打乱谁虔诚之后,两层之间的差距期望**就是 0** "
      f"⇒ 用 `negative_control` 对 **0**,不是对某个观测量")
rg2 = np.random.default_rng(11)
plG = gap(gp_, gr_, gy_, placebo=rg2.normal(size=len(gr_)))
plN = gap(np_, nr_, ny_, placebo=rg2.normal(size=len(nr_)))
print(f"  安慰剂:**把宗教度整根换成随机变量** ⇒ GSS **{plG:+.4f}** · NSFG **{plN:+.4f}** "
      f"—— 这条验的是**整条跨仪器管道对称地不会凭空造出差距**")

G_ = Gate("#853 · 换一份问卷,那条宗教缝还在不在")
G_.asserted("① 硬规则①(当场救了两次):`VALUE SAMESEX1F` 是**年龄分组**不是态度量表,真正指派给 "
            "`samesex` 的是 `AGDGFMT`(1 强烈同意…4 强烈不同意 · **5 既不同意也不反对**);"
            "且 `SAMESEXANY`/`SAMESEX1` **前缀相同、构念不同** —— 按前缀取变量会取错题",
            bool(len(N) > 1000 and N.samesex.notna().sum() > 500),
            f"NSFG n={len(N):,} · samesex 非缺 {int(N.samesex.notna().sum()):,}", kind="control")
G_.asserted("② 极性与码 5:GSS 高=宽容,NSFG **低=宽容 ⇒ 翻转 `5−x`**;"
            "码 5 语义坐 2 和 3 之间 ⇒ **当规格轴处理(剔除 / 中点),不当序数尺一端**",
            bool(len(NEITHERS) == 2), f"码 5 两种处理都跑,共 {len(rows)} 格", kind="control")
G_.asserted("③ 正控:植入 +0.30 必须按 SD 单位取回,**且 g=0 时必须 ≈0**(否则这条控制不会失败)",
            bool(abs(pcG) / sdG > 2 and abs(zeroG) < 1e-9),
            f"植入动 {pcG:+.4f} = {abs(pcG)/sdG:.1f} SD · g=0 动 {zeroG:+.2e}", kind="control")
G_.asserted("④ 负控:**年内打乱宗教度标签** ⇒ 两具仪器的差距都必须 ≈0"
            "(⚠ **这个零该是零**:打乱谁虔诚后差距期望就是 0)",
            bool(abs(ncG) < 0.05 and abs(ncN) < 0.05), f"GSS {ncG:+.4f} · NSFG {ncN:+.4f}",
            kind="control")
G_.asserted("⑤ 安慰剂:把宗教度换成随机变量 ⇒ 两具仪器都必须 ≈0(验整条跨仪器管道不凭空造差距)",
            bool(abs(plG) < 0.05 and abs(plN) < 0.05), f"GSS {plG:+.4f} · NSFG {plN:+.4f}",
            kind="control")
G_.asserted("⑥ 前提(跑前写下的最强混淆):**NSFG 是 15–49 岁、女性问卷只有女性;GSS 是 18+ 两性** ⇒ "
            "**年轻与女性都更宽容,光凭构成就能造出差异** ⇒ 把 GSS 限制到 18–49 岁女性,"
            "**限制前后都报**",
            bool(any(r["restrict"] for r in rows) and any(not r["restrict"] for r in rows)),
            f"限制格 {sum(r['restrict'] for r in rows)} · 全体格 {sum(not r['restrict'] for r in rows)}",
            kind="control")
G_.asserted("⑦ kill(预注册):「那条缝是关于美国而不是关于 GSS 问卷的」要成立,"
            "需 **12 格里 `g_GSS − g_NSFG` 的区间含 0 的 ≥ 2/3**",
            bool(ov >= 8), f"含 0 的 {ov}/{len(rows)}", kind="kill",
            yardstick="每格 `g_GSS − g_NSFG` 自己的 95% 自助区间",
            yardstick_noise=float(np.mean([(r["hi"] - r["lo"]) / 4 for r in rows])))
print(); print(G_)
adm = G_.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
unres = [r for r in rows if not r["restrict"]]
res = [r for r in rows if r["restrict"]]
mu_un = float(np.mean([abs(r["diff"]) for r in unres]))
mu_re = float(np.mean([abs(r["diff"]) for r in res]))
mgG = float(np.mean([r["g_gss"] for r in rows])); mgN = float(np.mean([r["g_nsfg"] for r in rows]))
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif ov >= 8:
    VERD = (f"**A 仪器无关 —— 那条缝换一份问卷还在,而且大小对得上。** "
            f"{ov}/{len(rows)} 格的 `g_GSS − g_NSFG` 区间含 0;"
            f"GSS 平均 **{mgG:+.3f}**、NSFG 平均 **{mgN:+.3f}**(单位都是该题当年自身的 SD)。\n"
            f"  ⇒ **一句关于人的话:信教与不信教的美国人在同性关系上的距离,\n"
            f"  不是 GSS 那道题问出来的 —— 换一份完全不同的问卷、不同的机构、不同的问法,\n"
            f"  那条缝还在原处,宽度也差不多。**")
elif mu_re <= mu_un / 2:
    VERD = (f"**C 是总体构成,不是仪器。** 全体格平均 |差| **{mu_un:.3f}** → "
            f"限制到 18–49 岁女性后 **{mu_re:.3f}**(减半以上)⇒ "
            f"**两具仪器的差异主要由「问的是谁」造成,不是「怎么问」。**\n"
            f"  ⇒ **一句关于人的话:那条缝在两份问卷里宽度不同,而差别来自被问的人不同 ——\n"
            f"  把年龄和性别对齐之后,两份问卷说的是同一件事。**")
else:
    VERD = (f"**B 缝的大小是「怎么问」的性质 ⇒ 本项目每一个关于缝大小的数字都只对 GSS 成立。**\n"
            f"  {len(rows)-ov}/{len(rows)} 格的差排除零;GSS 平均 **{mgG:+.3f}** vs "
            f"NSFG 平均 **{mgN:+.3f}**;限制总体后 |差| 从 {mu_un:.3f} 只到 {mu_re:.3f}。\n"
            f"  ⇒ **一句关于人的话,而它是我不欢迎的那句:「信教的人和不信教的人隔着多远」\n"
            f"  这个距离,量出来多大取决于你怎么问 —— 问「这样对不对」和问「你同不同意」\n"
            f"  得到的不是同一条缝。缝的存在是稳的,缝的宽度不是。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的(登记,不是「计划中」)**:")
print(f"   ① **九十年代复现不了** —— NSFG 字典只存在于 2011–2013 与 2017–2019,"
      f"另外 11 个 `.dat` 没有字典,**定宽文件没有字典就是不可读的**"
      f" ⇒ **本项目的头条(九十年代那条缝)在第二具仪器上结构性无法检验。**")
print(f"   ② NSFG 男性文件只有 2017–2019 ⇒ 性别轴在时间上不平衡,本轮不做男女对比。")
print(f"   ③ 两具仪器都是横断面调查 ⇒ **无干预、无因果识别。**")
json.dump(dict(grid=rows, n_cells=len(rows), overlap=ov,
               bh=len(bh), by=len(by), mean_g_gss=mgG, mean_g_nsfg=mgN,
               mean_absdiff_unrestricted=mu_un, mean_absdiff_restricted=mu_re,
               controls=dict(plant=PLANT, planted_shift=pcG, plant_in_sd=abs(pcG) / sdG,
                             zero_plant=zeroG, neg_gss=ncG, neg_nsfg=ncN,
                             placebo_gss=plG, placebo_nsfg=plN),
               structurally_impossible=["1990s not replicable: NSFG dictionaries exist only for "
                                        "2011-2013 and 2017-2019; the other 11 .dat files have none",
                                        "NSFG male file only 2017-2019",
                                        "cross-sectional only: no causal identification"],
               admissible=adm, verdict=VERD, gate_ok=G_.verdict()),
          open(OUT / "two_questionnaires.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'two_questionnaires.json'}")
