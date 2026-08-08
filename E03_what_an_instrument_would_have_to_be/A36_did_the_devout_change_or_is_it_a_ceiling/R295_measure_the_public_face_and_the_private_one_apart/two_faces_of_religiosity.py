"""#856 · E03·A89·R295 —— 虔诚有两张脸:公共的那张和私人的那张,量出来的缝不一样宽吗?

`#853`① 留下一条线索:12 格里两具仪器分歧最大的是 **`重要性 × 限制`**(+0.301/+0.292,
其余十格都在 ±0.09 内),含零但很勉强(p=0.084/0.102)。

**⓪ 硬规则①当场把那条线索改写成了一个更好的问题 —— 因为分歧的来源是我自己造的:**
**GSS 的 `reliten` 是「Strength of affiliation」(隶属强度),
而 NSFG 的 `reldlife` 是「How important is religion in R's daily life」(日常生活中的重要性)。**
**这是两道不同的题,而我在 `#853` 里把它们都贴上了「重要性」的标签。**
⇒ 那 +0.30 最可能是**我的构念替换**,不是关于女性 18–49 的事实。

**⚠⚠ 而这个错指向一个真问题,并且它比原来那个好:**
**「虔诚」不是一个东西。** `attend`(多久去一次礼拜)量的是**公共/群体**的那一面;
`relpersn`(是否自认为是个有宗教信仰的人)、`pray`(多久祷告一次)量的是**私人**的那一面。
**如果用哪一张脸来分层会改变缝的宽度,那么「信教的人和不信教的人隔多远」这句话
就缺了一个必须说出来的限定语。**

⇒ **而这个问题可以在单具仪器内部问,不需要跨仪器比较** ——
   这既避开了 `#853` 那个混淆,又让同一个量能在两具仪器上**各自独立**地算一遍(硬规则④)。

`G1` **估计量(先于方法命名)**:
   **`Δfacet = g(公共脸) − g(私人脸)`**,其中 `g = (mean_虔诚 − mean_世俗)/SD_年内合并`
   (与 `#853` 同一构造)。**在每具仪器内部各算一遍。**

**识别**:两张脸的题在同一份问卷、同一批人身上都问过 ⇒ 可估。
⚠ **但见下面的算术陷阱** —— 它决定这个量是否**根本能非零**。

**⚠⚠⚠ 算术陷阱(`realstat` 开篇那条),必须在读任何结果之前检查:**
**两张脸是相关的 —— 都是虔诚度。** 若它们切出的三分位**几乎相同**,
**`Δfacet` 就被代数逼近 0,与心理学无关。**
⇒ 控制:**逐格印出两张脸切出的「虔诚层」的重合度(Jaccard)**;
**重合度过高的格,其 `Δfacet ≈ 0` 一律记为 `DEGENERATE`,不计入判词** ——
**那是「统计量退化」,不是「两张脸一样」**(`realstat`:`floor == ceiling` 时无阈可用)。

三个世界:
   A **脸重要**:`Δfacet` 在两具仪器上同号且排零 ⇒
     **「信教的人和不信教的人隔多远」缺了一个限定语:按哪张脸分的。**
   B **脸不重要**:`Δfacet` 在两具仪器上都含零 ⇒ 对这个用途而言虔诚度实际上是一维的,
     **而 `#853` 那 +0.30 就是我的构念替换造出来的。**
   C **两具仪器对 `Δfacet` 不一致** ⇒ **脸效应本身是仪器绑定的**,说不了。

预测矩阵:
   | 世界 | 现在 | 两具同号排零 | 两具都含零 | 两具不一致 |
   | A 脸重要   | 0.40 | **0.85** | 0.05 | 0.10 |
   | B 脸不重要 | 0.35 | 0.05 | **0.85** | 0.10 |
   | C 仪器绑定 | 0.25 | 0.10 | 0.10 | **0.80** |

预注册判词(条件式):
  if 正控开火(**只在公共脸上植入的差距必须被 `Δfacet` 取回,且 plant=0 时必须 ≈0**)
     and 负控为零(**年内打乱虔诚度 ⇒ 两个 g 都 ≈0 ⇒ `Δfacet` ≈0**)
     and 安慰剂给出分辨率下限(**两个随机分层器的 `Δfacet` 分布** —— 这是本量的噪声地板)
     and 非退化(**两张脸的虔诚层 Jaccard < 0.80**):
      GSS 的**真私人脸**格 与 NSFG 格 同号且都超安慰剂地板 -> A
      两具都在地板内                                        -> B
      否则                                                  -> C
  else: UNVERIFIED

⚠⚠ **kill 的总体(`#855` 刚补的那个参数,这一轮第一次用):**
   **`reliten`(隶属强度)那两格留在网格里作对照,但不计入 kill** ——
   因为它正是 `#853` 里被我错贴成「重要性」的那一个,**它是构念替换的证据,不是私人脸的测量。**
   ⇒ kill 只在**真私人脸**(`relpersn` / `pray` / `reldlife`)的格上评。

`G3` 多重性:**4 对脸 × 2 总体 = 8 格**,BH 与 BY 都做,**不同意的格一起发表**。
`G4` 规格曲线:上面两根轴就是曲线,逐格报。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① NSFG 只有**一张**私人脸(`reldlife`),没有 `pray` 的对应题
   ⇒ **NSFG 侧的私人脸不能做规格扫描**,只有一个算法。
② `relpersn` 从 **1998** 年才问、`pray` 从 **1983** ⇒ **这两张脸都问不到九十年代以前**,
   本轮固定在 2010–2019 以与 `#853`/NSFG 对齐。
③ 两具都是横断面 ⇒ **无干预、无因果识别**;「哪张脸是原因」本轮不问也答不了。
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
B = 2000

# ── GSS ────────────────────────────────────────────────────────────────────
g = pd.read_stata(EXT / "gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "age", "sex", "homosex", "attend", "reliten", "pray",
                           "relpersn"], convert_categoricals=False)
G = pd.DataFrame({"year": g.year})
G["perm"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("pray", (1, 6)),
                    ("relpersn", (1, 4)), ("age", (18, 89)), ("sex", (1, 2))):
    G[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
G = G[(G.year >= 2010) & (G.year <= 2019)]
# 高 = 虔诚:attend 高=常去;reliten/pray/relpersn 低=更虔诚 ⇒ 取负
G["f_attend"] = G.attend
G["f_reliten"] = -G.reliten
G["f_pray"] = -G.pray
G["f_relpersn"] = -G.relpersn

print("=== ⓪ 硬规则①:每张脸真的问了什么、多少人(标签原文)===")
LAB = {"attend": "How often R attends religious services(**公共/群体**)",
       "reliten": "Strength of affiliation(**隶属强度 —— `#853` 里被我错贴成「重要性」**)",
       "pray": "How often Does R pray(**私人实践**)",
       "relpersn": "R consider self a religious person(**私人身份**)"}
for c, l in LAB.items():
    n = int(G[c].notna().sum()); ys = sorted(G.year[G[c].notna()].unique().tolist())
    print(f"  GSS `{c}` n={n:>6,} 年份 {ys} 「{l}」")
print(f"  ⚠ **`reliten` = 隶属强度 ≠ NSFG 的 `reldlife` = 日常生活中的重要性** —— "
      f"**两道不同的题,`#853` 把它们都贴上了「重要性」**")

# ── NSFG(与 `#853` 同一读法)────────────────────────────────────────────────
def dct_cols(dct, want):
    out = {}
    for m in re.finditer(r"_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)", dct.read_text(errors="replace")):
        out[m.group(2).lower()] = (int(m.group(1)), int(m.group(3)))
    return {w: out[w] for w in want if w in out}
NF = [("2011_2013_FemRespData.dat", "2011_2013_FemRespSetup.dct", 2012, "F"),
      ("2017_2019_FemRespData.dat", "2017_2019_FemRespSetup.dct", 2018, "F"),
      ("2017_2019_MaleData.dat", "2017_2019_MaleSetup.dct", 2018, "M")]
WANT = ["samesex", "attndnow", "reldlife", "ager"]
fr = []
for dat, dct, yr, sx in NF:
    cols = dct_cols(EXT / "nsfg/setup" / dct, WANT)
    rows = []
    with open(EXT / "nsfg" / dat, "r", errors="replace") as fh:
        for line in fh:
            r = {}
            for v, (c, w) in cols.items():
                s = line[c - 1:c - 1 + w].strip()
                r[v] = float(s) if s.lstrip("-").isdigit() else np.nan
            rows.append(r)
    d = pd.DataFrame(rows); d["year"] = yr; d["sex"] = sx; fr.append(d)
N = pd.concat(fr, ignore_index=True)
N = N[N.ager >= 18]
ss = N.samesex.where(lambda v: v.isin([1, 2, 3, 4]))     # 码 5「既不同意也不反对」剔除(`#853` 的规格轴之一)
N["perm"] = 5 - ss                                        # ⚠ 极性翻转 ⇒ 高=宽容
N["f_attnd"] = 8 - N.attndnow.where(lambda v: v.between(1, 7))
N["f_reldlife"] = 4 - N.reldlife.where(lambda v: v.between(1, 3))
print(f"  NSFG n={len(N):,} · `attndnow`(**公共**)非缺 {int(N.f_attnd.notna().sum()):,} · "
      f"`reldlife`(**私人:日常重要性**)非缺 {int(N.f_reldlife.notna().sum()):,}")

zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0

def gap(perm, rel, yr, rng=None, permute=False, plant=0.0, plant_on=None):
    """年内三分位;g = (虔诚均值 − 世俗均值)/年内合并 SD。plant_on:只在这个分层器的虔诚层加 plant。"""
    out = []
    for y in np.unique(yr):
        i = yr == y
        p, r = perm[i].astype(float).copy(), rel[i]
        if permute: r = rng.permutation(r)
        if len(p) < 200: continue
        lo, hi = np.quantile(r, [1/3, 2/3])
        if plant and plant_on is not None:
            po = plant_on[i]
            plo, phi = np.quantile(po, [1/3, 2/3])
            p[po >= phi] += plant                    # ⚠ 只按**公共脸**的层加,不按当前分层器
        dev, sec = p[r >= hi], p[r <= lo]
        if len(dev) < 60 or len(sec) < 60: continue
        s = p.std(ddof=1)
        if s <= 0: continue
        out.append((dev.mean() - sec.mean()) / s)
    return float(np.mean(out)) if out else np.nan

def jaccard(a, b, yr):
    """两张脸切出的**虔诚层**重合度 —— 算术陷阱的检查。"""
    js = []
    for y in np.unique(yr):
        i = yr == y
        if i.sum() < 200: continue
        ha = a[i] >= np.quantile(a[i], 2/3); hb = b[i] >= np.quantile(b[i], 2/3)
        u = (ha | hb).sum()
        if u: js.append(float((ha & hb).sum() / u))
    return float(np.mean(js)) if js else np.nan

def dfacet(perm, pub, priv, yr, Bv=B, seed=856):
    obs = gap(perm, pub, yr) - gap(perm, priv, yr)
    rg = np.random.default_rng(seed); o = np.empty(Bv)
    for i in range(Bv):
        k = rg.integers(0, len(perm), len(perm))
        o[i] = gap(perm[k], pub[k], yr[k]) - gap(perm[k], priv[k], yr[k])
    o = o[np.isfinite(o)]
    return obs, o

PAIRS = [("GSS", "f_attend", "f_relpersn", True,  "公共 vs **私人身份**(自认为有信仰的人)"),
         ("GSS", "f_attend", "f_pray",     True,  "公共 vs **私人实践**(祷告频率)"),
         ("GSS", "f_attend", "f_reliten",  False, "公共 vs **隶属强度** ⚠ `#853` 的构念替换,**对照,不计入 kill**"),
         ("NSFG", "f_attnd", "f_reldlife", True,  "公共 vs **私人重要性**(日常生活中多重要)")]
POPS = [("全体", False), ("女性18–49", True)]

print(f"\n=== ① 规格曲线:{len(PAIRS)} 对脸 × {len(POPS)} 总体 = **{len(PAIRS)*len(POPS)} 格**(`G3`/`G4`)===")
rows, CELLS = [], []
for instr, pub, priv, counts, desc in PAIRS:
    for pname, restrict in POPS:
        d = G if instr == "GSS" else N
        m = d.perm.notna() & d[pub].notna() & d[priv].notna()
        if restrict:
            m &= (d.age <= 49) & (d.sex == 2) if instr == "GSS" else (d.ager <= 49) & (d.sex == "F")
        sub = d[m]
        pe = sub.perm.to_numpy(float); pu = zs(sub[pub]).to_numpy(float)
        pr = zs(sub[priv]).to_numpy(float); yr = sub.year.to_numpy(int)
        if len(pe) < 500:
            print(f"  {instr:5s} {desc[:34]:36s} {pname:9s} **n={len(pe)} 太小,跳过**"); continue
        jac = jaccard(pu, pr, yr)
        obs, bs = dfacet(pe, pu, pr, yr)
        lo, hi = np.quantile(bs, [.025, .975])
        rows.append(dict(instr=instr, pub=pub, priv=priv, pop=pname, counts=counts, desc=desc,
                         g_pub=gap(pe, pu, yr), g_priv=gap(pe, pr, yr), dfacet=obs,
                         lo=float(lo), hi=float(hi), jaccard=jac, n=len(pe),
                         excl=bool(lo > 0 or hi < 0)))
        CELLS.append((instr, priv, pname))
        print(f"  {instr:5s} {desc[:34]:36s} {pname:9s} n={len(pe):>6,} · "
              f"g_公共 **{rows[-1]['g_pub']:+.3f}** · g_私人 **{rows[-1]['g_priv']:+.3f}** · "
              f"**Δfacet {obs:+.3f}** [{lo:+.3f},{hi:+.3f}] · Jaccard **{jac:.3f}**"
              f"{'  **排零**' if rows[-1]['excl'] else ''}")

print("\n=== ② 控制 ===")
sub = G[G.perm.notna() & G.f_attend.notna() & G.f_relpersn.notna()]
pe = sub.perm.to_numpy(float); pu = zs(sub.f_attend).to_numpy(float)
pr = zs(sub.f_relpersn).to_numpy(float); yr = sub.year.to_numpy(int)
PLANT = 0.40
d_plant = (gap(pe, pu, yr, plant=PLANT, plant_on=pu) - gap(pe, pr, yr, plant=PLANT, plant_on=pu))
d_base = gap(pe, pu, yr) - gap(pe, pr, yr)
d_zero = (gap(pe, pu, yr, plant=0.0, plant_on=pu) - gap(pe, pr, yr, plant=0.0, plant_on=pu)) - d_base
print(f"  正控:**只按公共脸的虔诚层**加 +{PLANT} ⇒ Δfacet {d_base:+.4f} → **{d_plant:+.4f}** "
      f"(动 **{d_plant-d_base:+.4f}**)· **而 plant=0 时动 {d_zero:+.6f}** —— "
      f"⚠ **`G2` 要求控制必须能失败,这一行就是那个检查**")
rg = np.random.default_rng(7)
nc = gap(pe, pu, yr, rng=rg, permute=True) - gap(pe, pr, yr, rng=rg, permute=True)
print(f"  负控:**年内打乱虔诚度** ⇒ 两个 g 都该塌 ⇒ Δfacet = **{nc:+.4f}** "
      f"(⚠ **「这个零该不该是零?」该** —— 打乱谁虔诚后两张脸各自的缝期望都是 0,差也是 0)")
rg2 = np.random.default_rng(11)
pl = np.array([gap(pe, rg2.normal(size=len(pe)), yr) - gap(pe, rg2.normal(size=len(pe)), yr)
               for _ in range(200)])
pl = pl[np.isfinite(pl)]
FLOOR = float(np.quantile(np.abs(pl), 0.95))
print(f"  **安慰剂(本量的噪声地板)**:**两个随机分层器**的 Δfacet 分布 ⇒ "
      f"中心 {np.mean(pl):+.4f} · |Δ| 的 95 分位 = **{FLOOR:.4f}**")
print(f"     ⚠ **这一条不可省** —— **任何两个不同的分层器,光凭抽样就会给出不同的缝**;"
      f"没有它,`Δfacet ≠ 0` 只是「它们是两个不同的变量」的同义反复")

DEG = [r for r in rows if r["jaccard"] >= 0.80]
kill_rows = [r for r in rows if r["counts"] and r["jaccard"] < 0.80]
gss_hit = [r for r in kill_rows if r["instr"] == "GSS" and abs(r["dfacet"]) > FLOOR and r["excl"]]
nsfg_hit = [r for r in kill_rows if r["instr"] == "NSFG" and abs(r["dfacet"]) > FLOOR and r["excl"]]
same_sign = bool(gss_hit and nsfg_hit and
                 np.sign(gss_hit[0]["dfacet"]) == np.sign(nsfg_hit[0]["dfacet"]))

Gt = Gate("#856 · 虔诚有两张脸:公共的和私人的")
Gt.asserted("① 硬规则①:**GSS `reliten` = Strength of affiliation(隶属强度),"
            "而 NSFG `reldlife` = 日常生活中的重要性 —— 两道不同的题**,"
            "而 `#853` 把它们都贴上了「重要性」⇒ 那 +0.30 最可能是**我的构念替换**",
            bool("reliten" in LAB and "relpersn" in LAB),
            "标签原文已逐条印出:attend=公共 · relpersn=私人身份 · pray=私人实践 · reliten=隶属强度",
            kind="control")
Gt.asserted("② **算术陷阱(`realstat` 开篇)**:两张脸都是虔诚度 ⇒ **若它们切出的三分位几乎相同,"
            "`Δfacet` 被代数逼近 0,与心理学无关** ⇒ 逐格印 Jaccard,**≥0.80 记 `DEGENERATE` 不计入判词**",
            bool(all(np.isfinite(r["jaccard"]) for r in rows)),
            " · ".join(f"{r['instr']}/{r['priv']}/{r['pop']}:{r['jaccard']:.3f}" for r in rows),
            kind="control")
Gt.asserted("③ 正控:**只按公共脸的虔诚层**植入 +0.40 ⇒ Δfacet 必须上移,"
            "**且 plant=0 时必须 ≈0**(否则这条控制不会失败)",
            bool((d_plant - d_base) > 0.05 and abs(d_zero) < 1e-9),
            f"植入动 {d_plant-d_base:+.4f} · plant=0 动 {d_zero:+.2e}", kind="control")
Gt.asserted("④ 负控:年内打乱虔诚度 ⇒ Δfacet 必须 ≈0"
            "(⚠ **这个零该是零**:打乱后两张脸各自的缝期望都是 0,差也是 0)",
            bool(abs(nc) < 0.10), f"{nc:+.4f}", kind="control")
Gt.asserted("⑤ **安慰剂 = 本量的噪声地板**:两个**随机**分层器的 |Δfacet| 95 分位 —— "
            "**不可省,因为任何两个不同的分层器光凭抽样就会给出不同的缝**",
            bool(np.isfinite(FLOOR) and FLOOR > 0), f"地板 {FLOOR:.4f}(200 次)", kind="control")
Gt.asserted("⑥ kill(预注册):「哪张脸重要」要成立,需 **GSS 的真私人脸格 与 NSFG 格 同号、"
            "都排零、且都超安慰剂地板**",
            same_sign, f"GSS 达标 {len(gss_hit)} · NSFG 达标 {len(nsfg_hit)} · 同号 {same_sign}",
            kind="kill", yardstick="每格 `Δfacet` 自己的 95% 自助区间,并须超过随机分层器的噪声地板",
            yardstick_noise=FLOOR,
            population=f"{len(kill_rows)} 个**真私人脸且非退化**的格 —— "
                       f"⚠ **`reliten`(隶属强度)那 {sum(1 for r in rows if not r['counts'])} 格留在网格里"
                       f"作对照但不计入**,因为它正是 `#853` 里被错贴成「重要性」的那一个")
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
# ⚠⚠⚠ **第一版的分支结构里没有 C，于是 `else` 把 C 的结果印成了 B。**
#    预测矩阵里我给 C(两具仪器不一致)配了 0.25 的先验,**却没有给它写一个分支** ——
#    实测正是 C:**GSS 四个真私人脸格全部排零且全部超地板(−0.136/−0.234/−0.215/−0.204),
#    而 NSFG 两格都在地板内外徘徊(+0.043 未超地板 · −0.041 含零)。**
#    ⇒ 这是 `#854` 那一类在**一轮之后**的复发,而位置又变了:
#    `#854` 错在 **kill 的总体**,这一次错在 **分支的结构** ——
#    **一个我在预测矩阵里明写过的世界,在代码里没有落点。**
#    ⇒ 规矩:**预测矩阵里的每一个世界,分支里都必须有一条与之对应的路;
#      `else` 不许充当任何一个已命名世界的落点。**
inside = [r for r in kill_rows if abs(r["dfacet"]) <= FLOOR]
gss_all = [r for r in kill_rows if r["instr"] == "GSS"]
nsfg_all = [r for r in kill_rows if r["instr"] == "NSFG"]
gss_unan = bool(gss_all) and all(abs(r["dfacet"]) > FLOOR and r["excl"] for r in gss_all)
nsfg_none = bool(nsfg_all) and all(not (abs(r["dfacet"]) > FLOOR and r["excl"]) for r in nsfg_all)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif same_sign:
    VERD = (f"**A 脸重要 —— 「隔多远」缺了一个限定语。** GSS {gss_hit[0]['desc'][:20]} "
            f"Δfacet **{gss_hit[0]['dfacet']:+.3f}**,NSFG **{nsfg_hit[0]['dfacet']:+.3f}**,同号且都超地板 "
            f"{FLOOR:.3f}。\n"
            f"  ⇒ **一句关于人的话:「信教的人和不信教的人隔多远」这句话本身不完整 ——\n"
            f"  按「多久去一次礼拜」分,和按「你自己算不算个有信仰的人」分,量出来的距离不一样宽。**")
elif gss_unan and nsfg_none:
    VERD = (f"**C 两具仪器不一致 —— 而这正是我在预测矩阵里给了 0.25、却忘了给它写分支的那个世界。**\n"
            f"  **GSS:{len(gss_all)}/{len(gss_all)} 个真私人脸格全部排零、全部超地板 {FLOOR:.3f}** —— "
            f"Δfacet {[round(r['dfacet'],3) for r in gss_all]},**一律为负**,\n"
            f"  即**公共脸(去不去礼拜)量出的缝,比私人脸(自认信仰 / 祷告频率)量出的更宽**;\n"
            f"  **NSFG:{len(nsfg_all)}/{len(nsfg_all)} 格都没达标**"
            f"（{[round(r['dfacet'],3) for r in nsfg_all]}，一个未超地板、一个含零）。\n"
            f"  ⚠⚠ **而「仪器不一致」这个说法本身有一个我分不开的混淆,必须写在同一句里:**\n"
            f"  **GSS 的私人脸是「自认是不是有信仰的人」和「多久祷告」,\n"
            f"  NSFG 的私人脸是「宗教在日常生活里多重要」—— 这不是同一张私人脸。**\n"
            f"  ⇒ **所以本轮分不清是「仪器不同」还是「私人脸不同」,两者在这里完全共线。**\n"
            f"  ⇒ **一句关于人的话,而它带着这条分不开:在 GSS 里,\n"
            f"  「你多久去一次礼拜」比「你算不算个有信仰的人」把人分得更开 ——\n"
            f"  公共的那张脸量出的缝,比私人的那张宽 0.14 到 0.23 个标准差。\n"
            f"  换到另一份问卷、也换到另一张私人脸,这个差就没了。**")
else:
    VERD = (f"**B 脸不重要(在本设计的分辨率上)。** 真私人脸且非退化的 {len(kill_rows)} 格里,"
            f"**{len(inside)} 格的 |Δfacet| 落在随机分层器地板 {FLOOR:.3f} 之内**。\n"
            f"  ⇒ **而这同时说明 `#853` 那 +0.30 最可能是我的构念替换** —— "
            f"`reliten` 是隶属强度,不是日常重要性。")
print(VERD)
if DEG:
    print(f"\n⚠ **退化格(Jaccard ≥ 0.80,`Δfacet` 被代数逼近 0,不计入判词)**:"
          f"{[(r['instr'], r['priv'], r['pop'], round(r['jaccard'],3)) for r in DEG]}")
print(f"\n⚠ **本站结构性做不到的**:① NSFG 只有**一张**私人脸(`reldlife`)⇒ **NSFG 侧不能做规格扫描**;"
      f"② `relpersn` 从 1998、`pray` 从 1983 才问 ⇒ **两张私人脸都够不到九十年代以前**;"
      f"③ 两具都是横断面 ⇒ **无干预、无因果识别**,「哪张脸是原因」本轮不问也答不了。")
json.dump(dict(grid=rows, degenerate=[r["desc"] for r in DEG], floor=FLOOR, n_cells=len(rows),
               kill_population=len(kill_rows), gss_hit=len(gss_hit), nsfg_hit=len(nsfg_hit),
               gss_unanimous=gss_unan, nsfg_none=nsfg_none,
               confound="instrument differs AND the private face differs (relpersn/pray vs reldlife) — collinear here",
               controls=dict(plant=PLANT, plant_moved=d_plant - d_base, zero_plant=d_zero,
                             neg=nc, placebo_floor=FLOOR),
               labels=LAB, admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "two_faces.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'two_faces.json'}")
