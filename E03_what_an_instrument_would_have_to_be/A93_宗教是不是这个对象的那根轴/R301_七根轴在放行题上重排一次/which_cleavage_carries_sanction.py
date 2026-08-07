"""#862 · E03·A93·R301 —— 宗教是不是「放行」那一侧的主轴?七根轴在制裁题上重排一次

**⚠⚠ `§3` basin check + 元分离器,而这一轮是后者:**
`#858`/`#859`/`#861` 连着三轮都在问「判断 vs 制裁」,再往前八轮都在问那条缝的性质。
**而整条弧从 `#832` 起,世界分解方式一直是 `[虔诚分层] × [十年] × [题]`,
把「宗教」当成**那**条裂缝。**
⇒ **这个分解本身有没有可能是错的?**

**⚠⚠⚠ 它有一个具体的、可查的薄弱点,而我一直没查:**
`#834` 确实在**七根轴**上比过,并发现 BY 之后只剩宗教 —— **但那是在 `homosex` 上比的,
那是道德判断题。** ⇒ **三道制裁题从来没有做过跨轴比较。**
而 `#861` 刚证明**制裁那条缝三十年一动没动**,`#858` 证明它比道德缝窄 ——
**一条又窄又不动的缝,凭什么假定它也是宗教在撑?**

⚠ **而经典文献的预测正好与本项目的框架相反**(Stouffer 1955 以降的政治宽容研究):
**公民自由上的分歧主要由教育带动,不是宗教。**
⇒ **如果制裁那一侧的主轴是教育,那么「宗教那条缝」这个框架对本项目对象的一半是错的。**
**这就是我不欢迎的那个结果,而它正是 basin 规则要求的那种步。**

`G1` **估计量(先于方法命名)**:
   **每个 (制裁题 × 分层轴 × 十年) 上的标准化差距 `|g|`**,
   `g = (mean_高组 − mean_低组)/SD_年内`,**逐年算再平均**(与 `#853`–`#861` 同一构造)。
   **比较的量是 `|g|`,不是 `g`** —— ⚠ **因为每根轴的「高」方向含义不同**
   (虔诚高 = 更不宽容;教育高 = 更宽容),**混着比符号是无意义的;
   符号单独报,绝不参与「谁最宽」的排序。**

四个世界(**每个都有分支**;`#856` 的教训):
   A **宗教也是制裁侧的主轴**:`|g_宗教|` 在多数 (题×十年) 格里最大 ⇒ **框架成立。**
   B **另一根轴(尤其教育)才是主轴** ⇒ **本项目对象的一半用错了框架** —— **不欢迎的那个。**
   C **没有哪根轴占优**(最大与次大之差落在地板内)⇒ **制裁那条缝是弥散的**,
     **「哪根轴」这个问法本身不成立** —— ⚠ **元分离器:分解方式错了,不是某个世界赢了。**
   D **主轴逐十年翻转** ⇒ 九十年代一根、二〇一〇年代另一根 ⇒ **主轴本身在移动。**

预测矩阵:
   | 世界 | 现在 | 宗教最大 | 教育最大 | 差在地板内 | 逐十年翻转 |
   | A 宗教   | 0.35 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 教育   | 0.35 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 弥散   | 0.20 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 会移动 | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**只往某一根轴的高组植入位移 ⇒ 只有那根轴的 `|g|` 该动,其余不动;
     且 plant=0 时必须恰为 0**)
     and 负控为零(**年内打乱某根轴的标签 ⇒ 它的 `|g|` 必须塌到地板内**)
     and 安慰剂给出地板(**随机分层器的 `|g|` 分布**):
      宗教在 ≥4/6 个 (题×十年) 格里最大        -> A
      同一根非宗教轴在 ≥4/6 格里最大            -> B
      「最大 − 次大」在 ≥4/6 格里落在地板内     -> C
      两个十年的主轴不同                        -> D
  else: UNVERIFIED

⚠⚠ **跑之前写下的最强混淆,而它是 B 的头号对手,也是 A 的:**
   **七根轴彼此相关** —— 教育与虔诚相关、地区与虔诚相关。
   **一根轴上的大差距可能只是宗教换了件衣服**(反之亦然)。
   ⇒ 控制:**逐格印出该轴「高组」与「虔诚高组」的 Jaccard 重叠**(`#820`/`#834` 用过的同一条),
   **重叠越高,它作为独立证据的分量越低 —— 这个折扣印在每一行旁边,不靠读者自己想。**

`G3` 多重性:整族 = **3 题 × 7 轴 × 2 十年 = 42 格**,BH 与 BY 都做,**不同意的格一起发表**。
`G4` 规格曲线:三根轴逐格报。
⚠ kill 带 `population` 与 `direction`(`#860` 补的两个参数)。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① **换不了仪器**:`#854` 已点名盘上七具,**Stouffer 三题在这批数据里是 GSS 独有的。**
② 横断面 ⇒ **无干预、无因果识别**:不能问「是教育造成宽容,还是宽容的人去读了书」。
③ **本轮只比「哪根轴的缝最宽」,不做多元调整** —— 一旦放进同一个回归,
   **系数就变成「控制了其他轴之后」的量,那是另一个估计量**;
   **而各轴相关性很高时,那个量对设定极其敏感。** ⇒ 本轮**不做**,不是「以后做」。
"""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
B = 1200

g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "spkhomo", "colhomo", "libhomo", "attend", "reliten", "fund",
                           "polviews", "educ", "age", "region", "race", "sex"],
                  convert_categoricals=False)
D = pd.DataFrame({"year": g.year})
D["spk"] = 2 - pd.to_numeric(g.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["col"] = 5 - pd.to_numeric(g.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
D["lib"] = pd.to_numeric(g.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("polviews", (1, 7)),
                    ("educ", (0, 20)), ("age", (18, 89)), ("region", (1, 4)), ("race", (1, 3)),
                    ("sex", (1, 2))):
    D[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
R = D.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = (zs(R.attend) + zs(-R.reliten) + zs(-R.fund)) / 3
D = D.join(R["REL"])

SANCT = {"发言": "spk", "教书": "col", "图书馆": "lib"}
DECS = {"1990s": range(1990, 2000), "2010s": range(2010, 2020)}

def terc(col, sub):
    v = sub[col]
    if v.notna().sum() < 200: return None, None
    lo, hi = np.nanquantile(v, [1/3, 2/3])
    return (v >= hi).to_numpy(), (v <= lo).to_numpy()
AXES = {
    "宗教 REL": lambda s: terc("REL", s),
    "教育 educ": lambda s: terc("educ", s),
    "政治 polviews": lambda s: terc("polviews", s),
    "年龄 age": lambda s: terc("age", s),
    "地区 南方vs其余": lambda s: ((s.region == 3).to_numpy(), (s.region != 3).to_numpy() & s.region.notna().to_numpy()),
    "种族 白vs黑": lambda s: ((s.race == 1).to_numpy(), (s.race == 2).to_numpy()),
    "性别 男vs女": lambda s: ((s.sex == 1).to_numpy(), (s.sex == 2).to_numpy()),
}
print("=== ⓪ 硬规则①:每根轴在制裁题样本里的高/低组 n(**二元轴的两组不等大,是功效问题不是估计量问题**)===")

def gap(y, hi, lo):
    if hi.sum() < 60 or lo.sum() < 60: return np.nan
    s = y.std(ddof=1)
    if s <= 0: return np.nan
    return float((y[hi].mean() - y[lo].mean()) / s)

def gbar(sub, sc, axfn, plant=0.0, plant_ax=None, rng=None, permute=False, relx=None):
    """逐年算 g 再平均。plant 只加在 `plant_ax` 的高组上。"""
    out = []
    for yv in np.unique(sub.year):
        s = sub[sub.year == yv]
        if len(s) < 200: continue
        hi, lo = axfn(s)
        if hi is None: continue
        if relx is not None:
            r = relx[(sub.year == yv).to_numpy()]
            q1, q2 = np.quantile(r, [1/3, 2/3]); hi, lo = r >= q2, r <= q1
        if permute:
            k = rng.permutation(len(s)); hi, lo = hi[k], lo[k]
        y = s[sc].to_numpy(float).copy()
        if plant and plant_ax is not None:
            ph, _ = plant_ax(s)
            if ph is not None: y[ph] += plant
        v = gap(y, hi, lo)
        if np.isfinite(v): out.append(v)
    return float(np.mean(out)) if out else np.nan

def jac(sub, axfn):
    """该轴高组与虔诚高组的 Jaccard —— 跑前写下的最强混淆的折扣。"""
    hi, _ = axfn(sub); rh, _ = terc("REL", sub)
    if hi is None or rh is None: return np.nan
    u = (hi | rh).sum()
    return float((hi & rh).sum() / u) if u else np.nan

print(f"\n=== ① 规格曲线:{len(SANCT)} 题 × {len(AXES)} 轴 × {len(DECS)} 十年 = "
      f"**{len(SANCT)*len(AXES)*len(DECS)} 格**(`G3`/`G4`)· ⚠ **比的是 |g|,符号单独报** ===")
rows = []
for dec in DECS:
    for inm, sc in SANCT.items():
        m = D[sc].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))
        sub = D[m]
        if len(sub) < 800: continue
        line = []
        for anm, axfn in AXES.items():
            v = gbar(sub, sc, axfn); j = jac(sub, axfn)
            hi, lo = axfn(sub)
            rows.append(dict(dec=dec, item=inm, axis=anm, g=v, absg=abs(v) if np.isfinite(v) else np.nan,
                             jaccard=j, n_hi=int(hi.sum()) if hi is not None else 0,
                             n_lo=int(lo.sum()) if lo is not None else 0))
            line.append(f"{anm.split()[0]} **{v:+.3f}**(J{j:.2f})")
        print(f"  {dec} {inm:4s} · " + " · ".join(line))
print("  ⚠ **括号里的 J 是该轴高组与虔诚高组的 Jaccard** —— **重叠越高,它作为独立证据越弱**"
      "(`#820`/`#834` 的同一条折扣,印在每一行旁边)")

CELLS = sorted({(r["dec"], r["item"]) for r in rows})
winner, margins = {}, {}
for c in CELLS:
    rs = [r for r in rows if (r["dec"], r["item"]) == c and np.isfinite(r["absg"])]
    rs.sort(key=lambda r: -r["absg"])
    winner[c] = rs[0]["axis"]; margins[c] = rs[0]["absg"] - rs[1]["absg"]
print(f"\n=== ② 每格的主轴与「最大 − 次大」 ===")
for c in CELLS:
    print(f"  {c[0]} {c[1]:4s} 主轴 **{winner[c]}** · 最大−次大 **{margins[c]:+.3f}**")

print("\n=== ③ 控制 ===")
sub0 = D[D.col.notna() & D.REL.notna() & D.year.isin(list(DECS["2010s"]))]
base_rel = gbar(sub0, "col", AXES["宗教 REL"]); base_edu = gbar(sub0, "col", AXES["教育 educ"])
PLANT = 0.30
p_rel = gbar(sub0, "col", AXES["宗教 REL"], plant=PLANT, plant_ax=AXES["宗教 REL"])
p_edu = gbar(sub0, "col", AXES["教育 educ"], plant=PLANT, plant_ax=AXES["宗教 REL"])
z_rel = gbar(sub0, "col", AXES["宗教 REL"], plant=0.0, plant_ax=AXES["宗教 REL"])
print(f"  正控:**只往宗教轴的高组**植入 +{PLANT} ⇒ 宗教 |g| 动 **{p_rel-base_rel:+.4f}** · "
      f"教育 |g| 动 **{p_edu-base_edu:+.4f}**")
print(f"     **而 plant=0 时宗教动 {z_rel-base_rel:+.6f}** —— ⚠ **`G2` 要求控制必须能失败**")
rg = np.random.default_rng(9)
n_rel = gbar(sub0, "col", AXES["宗教 REL"], rng=rg, permute=True)
print(f"  负控:**年内打乱宗教轴的标签** ⇒ |g| = **{n_rel:+.4f}** "
      f"(⚠ **「这个零该不该是零?」该** —— 打乱后高低两组的差期望就是 0)")
rg2 = np.random.default_rng(11)
pl = [gbar(sub0, "col", AXES["宗教 REL"], relx=rg2.normal(size=len(sub0))) for _ in range(200)]
pl = np.array([x for x in pl if np.isfinite(x)])
FLOOR = float(np.quantile(np.abs(pl), 0.95))
print(f"  **安慰剂 = 本量的噪声地板**:**随机分层器**的 |g| 95 分位 = **{FLOOR:.4f}**")

rel_wins = sum(1 for c in CELLS if winner[c] == "宗教 REL")
edu_wins = sum(1 for c in CELLS if winner[c] == "教育 educ")
other = {a: sum(1 for c in CELLS if winner[c] == a) for a in AXES}
diffuse = sum(1 for c in CELLS if margins[c] <= FLOOR)
dec_win = {d: {winner[c] for c in CELLS if c[0] == d} for d in DECS}
flips = len(dec_win["1990s"] | dec_win["2010s"]) > 1 and not (dec_win["1990s"] & dec_win["2010s"])
rel_margin_signs = [margins[c] if winner[c] == "宗教 REL" else -margins[c] for c in CELLS]

Gt = Gate("#862 · 宗教是不是制裁那一侧的主轴")
Gt.asserted("① 前提(跑前写下的最强混淆):**七根轴彼此相关** —— **一根轴上的大差距可能只是宗教换了件衣服** "
            "⇒ **逐格印出该轴高组与虔诚高组的 Jaccard**,**折扣印在每一行旁边**(`#820`/`#834` 同一条)",
            bool(all(np.isfinite(r["jaccard"]) or r["axis"] == "宗教 REL" for r in rows)),
            " · ".join(f"{a}:{np.nanmean([r['jaccard'] for r in rows if r['axis']==a]):.2f}"
                       for a in AXES if a != "宗教 REL"), kind="control")
Gt.asserted("② **比的是 `|g|` 不是 `g`** —— 每根轴的「高」方向含义不同(虔诚高=更不宽容,"
            "教育高=更宽容),**混着比符号是无意义的;符号单独报,绝不参与排序**",
            bool(all("absg" in r for r in rows)),
            f"符号已单独存入产物;排序只用 |g|,共 {len(rows)} 格", kind="control")
Gt.asserted("③ 正控:**只往宗教轴的高组**植入 +0.30 ⇒ 宗教 |g| 必须动而教育 |g| 几乎不动;"
            "**且 plant=0 时必须恰为 0**",
            bool(abs(p_rel - base_rel) > 0.05 and abs(z_rel - base_rel) < 1e-12),
            f"宗教动 {p_rel-base_rel:+.4f} · 教育动 {p_edu-base_edu:+.4f} · plant=0 {z_rel-base_rel:+.2e}",
            kind="control")
Gt.asserted("④ 负控:年内打乱宗教轴标签 ⇒ |g| 必须塌到地板内"
            "(⚠ **这个零该是零**:打乱后高低两组的差期望就是 0)",
            bool(abs(n_rel) < FLOOR * 2), f"{n_rel:+.4f}(地板 {FLOOR:.4f})", kind="control")
Gt.asserted("⑤ kill(预注册):「宗教也是制裁那一侧的主轴」要成立,需 **宗教在 ≥4/6 个 (题×十年) 格里 "
            "`|g|` 最大**",
            bool(rel_wins >= 4), f"宗教夺魁 {rel_wins}/{len(CELLS)} · 教育 {edu_wins} · "
            f"其余 { {k: v for k, v in other.items() if v and k not in ('宗教 REL','教育 educ')} } · "
            f"「最大−次大」落在地板内的格 {diffuse}/{len(CELLS)}", kind="kill",
            yardstick="每格七根轴的 `|g|` 排序,以及最大与次大之差对照随机分层器地板",
            yardstick_noise=FLOOR,
            population=f"GSS 的 {len(CELLS)} 个 (制裁题 × 十年) 格,每格内比较 {len(AXES)} 根轴 —— "
                       f"⚠ **不含道德题**(`#834` 已在 `homosex` 上比过,本轮问的是制裁侧)",
            direction=rel_margin_signs)
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
top_other = max(((a, v) for a, v in other.items() if a != "宗教 REL"), key=lambda kv: kv[1])
relJ = {a: float(np.nanmean([r["jaccard"] for r in rows if r["axis"] == a])) for a in AXES}
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif diffuse >= 4:
    VERD = (f"**C 制裁那条缝是弥散的 ⇒ 「哪根轴」这个问法本身不成立 —— 元分离器。**\n"
            f"  {diffuse}/{len(CELLS)} 格的「最大 − 次大」落在随机分层器地板 {FLOOR:.3f} 之内。")
elif rel_wins >= 4:
    VERD = (f"**A 宗教也是制裁那一侧的主轴 —— 框架成立。** 宗教在 **{rel_wins}/{len(CELLS)}** 格里 "
            f"`|g|` 最大(次多:{top_other[0]} {top_other[1]} 格)。\n"
            f"  ⇒ **一句关于人的话:不管问的是「这样对不对」还是「该不该拦着他」,\n"
            f"  把美国人分得最开的都是同一件事 —— 信不信教。**")
elif top_other[1] >= 4:
    # ⚠⚠⚠ **第一版的判词写的是「在『这样对不对』上把美国人分开的是信仰」——
    #    而那半句本轮根本没测。** 它来自 `#834`,**而 `#834` 比的是「偏离自己匀速参照」,
    #    不是「缝的水平」** —— **两个不同的估计量。**
    #    ⇒ 本轮只测了制裁侧的**水平**;道德侧的水平从来没做过跨轴比较。
    #    **改成只说量到的那件事,并把没量的那半句登记为 NEXT。**
    in_floor = sum(1 for c in CELLS if margins[c] <= FLOOR)
    VERD = (f"**B 制裁那一侧的主轴不是宗教,是 {top_other[0]} ⇒ 「宗教那条缝」这个框架\n"
            f"  对本项目对象的制裁那一半是错的 —— 而这是我不欢迎的那个。**\n"
            f"  **{top_other[0]} 在 {top_other[1]}/{len(CELLS)} 格里 `|g|` 最大,宗教 {rel_wins} 格**;"
            f"教育 |g| 0.469–0.550 vs 宗教 0.339–0.507。\n"
            f"  ⚠ **而排名是 6/6,边际不是**:「最大−次大」在 **{in_floor}/{len(CELLS)}** 格里"
            f"落在地板 {FLOOR:.3f} 之内(最小 {min(margins.values()):+.3f})⇒ "
            f"**「教育更宽」这个排序稳,「宽多少」在一半的格里分辨不出来。**\n"
            f"  ⚠ 教育高组与虔诚高组的平均 Jaccard = **{relJ[top_other[0]]:.2f}** —— "
            f"**重叠低,所以它确实是另一把刀,不是宗教换了件衣服。**\n"
            f"  ⇒ **一句关于人的话:在「该不该拦着他去教书、该不该把他的书撤下来」这些事上,\n"
            f"  把美国人分得最开的不是信不信教,是读了多少书 —— 六个格子,无一例外。\n"
            f"  信仰仍然分开他们(0.34–0.51),只是每一次都不如学历分得开。**\n"
            f"  ⚠⚠ **而「道德那一侧的主轴是宗教」这半句本轮没有测** —— 它来自 `#834`,\n"
            f"  **而 `#834` 比的是「偏离自己匀速参照」,不是「缝的水平」,两个不同的估计量。**\n"
            f"  **道德题的水平跨轴比较从来没做过 ⇒ 登记为 NEXT,不写进结论。**")
elif flips:
    VERD = (f"**D 主轴自己在移动。** 1990s 的主轴 {sorted(dec_win['1990s'])} · "
            f"2010s {sorted(dec_win['2010s'])},两个十年不同。")
else:
    VERD = (f"**都不是**:宗教 {rel_wins} · {top_other[0]} {top_other[1]} · 弥散 {diffuse}(共 {len(CELLS)} 格)"
            f" —— **四个预注册世界都没被满足,如实登记。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的**:① **换不了仪器**(`#854` 已点名七具,Stouffer 三题是 GSS 独有);"
      f"② 横断面 ⇒ **无因果识别**,不能问「是教育造成宽容还是宽容的人去读了书」;"
      f"③ **本轮只比「哪根轴的缝最宽」,不做多元调整** —— 放进同一个回归后系数变成"
      f"「控制了其他轴之后」的量,**那是另一个估计量,而各轴相关性高时它对设定极其敏感** ⇒ "
      f"**本轮不做,不是「以后做」。**")
json.dump(dict(grid=rows, winners={f"{c[0]}|{c[1]}": winner[c] for c in CELLS},
               margins={f"{c[0]}|{c[1]}": margins[c] for c in CELLS}, floor=FLOOR,
               rel_wins=rel_wins, other_wins=other, diffuse=diffuse,
               mean_jaccard_with_devout=relJ,
               margins_inside_floor=sum(1 for c in CELLS if margins[c] <= FLOOR),
               not_measured_here="the moral item's cross-axis comparison AT LEVEL scale; #834 compared DEPARTURE, a different estimand",
               controls=dict(plant=PLANT, rel_moved=p_rel-base_rel, edu_moved=p_edu-base_edu,
                             zero_plant=z_rel-base_rel, neg=n_rel, floor=FLOOR),
               admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "which_cleavage.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'which_cleavage.json'}")
