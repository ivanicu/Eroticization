"""#857 · E03·A89·R296 —— 同一年、同一批人,把三张私人脸并排放

`#856` 登记了一条**分不开的共线**:GSS 说「公共脸比私人脸把人分得更开」(四格一致,−0.14~−0.23),
NSFG 说没有 —— **而 GSS 的私人脸是「自认是不是有信仰的人」与「多久祷告」,
NSFG 的私人脸是「宗教在日常生活里多重要」。这不是同一张私人脸。**
⇒ **「仪器不同」与「私人脸不同」在 `#856` 里完全共线。**

**⓪ 硬规则①把这条共线拆开了,而拆它的东西一直在盘上:**
GSS 有 399 个宗教变量。逐个印 n 与**真的问过哪些年**之后:
   `religimp`「How important is religion」—— **只问过 1 年:2021**(n=3,609,与 `homosex` 同时非缺 **2,363**)
   `relidimp`「How important is R's religion」—— **只问过 2021**
   `godguide`「Feel guided by god in daily activities」—— **只问过 2004**
   `relpersn`「R consider self a religious person」—— 11 年,**含 2021**
   `pray`「How often Does R pray」—— **含 2021**
⇒ **2021 年这一份问卷里,「重要性」「身份」「实践」三张私人脸和「出席」这张公共脸同时存在,
问的是同一批人。** ⇒ **仪器、年份、人,全都固定;只有「脸」在变。**

`G1` **估计量(先于方法命名)**:
   **`Δfacet(f) = g(attend) − g(f)`**,`f ∈ {重要性, 身份, 实践}`,**全部在 GSS 2021 内部**。
   `g = (mean_虔诚 − mean_世俗)/SD` —— 与 `#853`/`#856` 同一构造。

三个世界:
   A **是脸,不是仪器**:`Δfacet(重要性)` 明显小于 `Δfacet(身份)`/`Δfacet(实践)`
     ⇒ **NSFG 之所以没有那个差,是因为它用的是「重要性」这张脸** ⇒ `#856` 的共线拆开,归因给脸。
   B **是仪器,不是脸**:三张私人脸给出的 `Δfacet` **彼此接近**(都为负、都超地板)
     ⇒ **在 GSS 里换哪张私人脸都一样,那 NSFG 的零就是仪器的事** ⇒ 共线拆开,归因给仪器。
   C **三张脸互相之间就不一致**,且没有「重要性 vs 其余」的清晰结构
     ⇒ **「公共/私人」这个二分本身不是切开这件事的正确关节** ——
     **⚠ 这是元分离器:它说的是我的世界分解方式错了,不是某个世界赢了。**

预测矩阵:
   | 世界 | 现在 | 重要性明显小 | 三张接近 | 三张互不一致 |
   | A 是脸   | 0.40 | **0.85** | 0.05 | 0.10 |
   | B 是仪器 | 0.35 | 0.05 | **0.85** | 0.10 |
   | C 关节错 | 0.25 | 0.10 | 0.10 | **0.80** |
**⚠ 三个世界各有一条分支,`else` 只留给「都不是」——`#856` 那个教训(预测矩阵里的世界
必须在代码里有落点)本轮当场执行。**

预注册判词(条件式,⓪ 排最前):
  ⓪ **功效闸**:n 从 `#856` 的 7,341 掉到约 2,363,**地板按 √n 会涨到约 0.14**,
     而 `#856` 观测到的 Δfacet 是 0.136–0.234 ⇒ **可能整段落进地板**。
     ⇒ **先在本年数据上植入一个等于 `#856` 观测量的 Δfacet,若捞不回 ≥2 个自助 SD
       ⇒ UNVERIFIED,不看结果**(`#835`/`#845`/`#849` 的做法:拦在读之前)。
  if 功效闸过 and 正控开火 and 负控为零 and 非退化(Jaccard<0.80):
      |Δ(重要性)| < 地板 且 至少一张其余脸 > 地板 -> A
      三张都 > 地板且同号                         -> B
      其余                                        -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**2021 年是 COVID 之后的一次特殊调查**(GSS 2021 改为网络/混合模式)。
  **模式变化会同时改变作答分布与谁回答** ⇒ 本年的任何绝对水平都不可与 2010s 直接比。
  ⇒ 控制:**本轮只做年内比较**(同一年、同一批人、只换脸),
  **绝不把 2021 的 `g` 与 `#856` 的 2010s `g` 并列解读**;并把 `attend` 在 2021 的 `g` 印出来,
  **让读者看见它与 2010s 的 −0.85 差多少,而不是让我在文字里替它解释。**

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① **「重要性」这张脸在 GSS 里只有 2021 一年** ⇒ **不能做跨年稳健性**,本轮是单年结果。
② 2021 的模式变化无法被本设计控制 —— **它对年内比较是共同的,对跨年比较是致命的**,
   所以跨年比较本轮**不做**,不是「以后做」。
③ 横断面 ⇒ **无干预、无因果识别**;「哪张脸是原因」问不了。
"""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
B, YEAR = 3000, 2021

g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "homosex", "attend", "religimp", "relidimp", "relpersn", "pray"],
                  convert_categoricals=False)
D = pd.DataFrame({"year": g.year})
D["perm"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
RANGES = {"attend": (0, 8), "religimp": (1, 5), "relidimp": (1, 5), "relpersn": (1, 4), "pray": (1, 6)}
for c, (lo, hi) in RANGES.items():
    D[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
D = D[D.year == YEAR]

print(f"=== ⓪ 硬规则①:GSS {YEAR} 年,每张脸真的问了多少人 ===")
for c in RANGES:
    n = int((D[c].notna() & D.perm.notna()).sum())
    vals = sorted(D[c].dropna().unique().tolist())[:7]
    print(f"  `{c}` 与 `homosex` 同时非缺 **{n:>5,}** · 取值 {vals}")
print(f"  ⚠ **`religimp`/`relidimp` 在整个 GSS 里只问过 {YEAR} 这一年** ⇒ 本轮是**单年**结果,"
      f"**不能做跨年稳健性**;而 `relpersn`/`pray` 也在这一年 ⇒ **三张私人脸同年同人并排。**")

# 高 = 虔诚
D["f_pub"] = D.attend
FACES = {"重要性 `religimp`": -D.religimp, "重要性2 `relidimp`": -D.relidimp,
         "身份 `relpersn`": -D.relpersn, "实践 `pray`": -D.pray}
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0

def gap(perm, rel, plant=0.0, plant_on=None, rng=None, permute=False):
    p, r = perm.astype(float).copy(), rel
    if permute: r = rng.permutation(r)
    lo, hi = np.quantile(r, [1/3, 2/3])
    if plant and plant_on is not None:
        phi = np.quantile(plant_on, 2/3)
        p[plant_on >= phi] += plant
    dev, sec = p[r >= hi], p[r <= lo]
    if len(dev) < 60 or len(sec) < 60 or p.std(ddof=1) <= 0: return np.nan
    return float((dev.mean() - sec.mean()) / p.std(ddof=1))

def jac(a, b):
    ha, hb = a >= np.quantile(a, 2/3), b >= np.quantile(b, 2/3)
    u = (ha | hb).sum()
    return float((ha & hb).sum() / u) if u else np.nan

def dfacet(perm, pub, priv, Bv=B, seed=857, plant=0.0):
    obs = gap(perm, pub, plant=plant, plant_on=pub) - gap(perm, priv, plant=plant, plant_on=pub)
    rg = np.random.default_rng(seed); o = np.empty(Bv)
    for i in range(Bv):
        k = rg.integers(0, len(perm), len(perm))
        o[i] = (gap(perm[k], pub[k], plant=plant, plant_on=pub[k])
                - gap(perm[k], priv[k], plant=plant, plant_on=pub[k]))
    o = o[np.isfinite(o)]
    return obs, o

# ── 安慰剂地板(本年 n 下的) ────────────────────────────────────────────────
base = D[D.perm.notna() & D.f_pub.notna() & D.relpersn.notna()]
pe0 = base.perm.to_numpy(float)
rg0 = np.random.default_rng(11)
pl = np.array([gap(pe0, rg0.normal(size=len(pe0))) - gap(pe0, rg0.normal(size=len(pe0)))
               for _ in range(300)])
pl = pl[np.isfinite(pl)]
FLOOR = float(np.quantile(np.abs(pl), 0.95))
# ⚠ **不手抄 `#856` 的地板与 n,从它的产物读**(`#840`/`#841` 立的规矩,而我的闸刚拦下了第一版)。
P856 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A89_虔诚有两张脸/"
                      "R295_measure_the_public_face_and_the_private_one_apart/results/two_faces.json", encoding="utf-8"))
F856 = P856["floor"]; N856 = max(r["n"] for r in P856["grid"] if r["instr"] == "GSS")
print(f"\n=== ⓪b 功效闸(排在所有关于人的判断之前)===")
print(f"  本年安慰剂地板(两个**随机**分层器的 |Δfacet| 95 分位,n={len(pe0):,}):**{FLOOR:.4f}**")
print(f"  ⚠ `#856` 在 n={N856:,} 上的地板是 {F856:.4f}(**从它的产物读,不手抄**)⇒ "
      f"**n 掉到 {len(pe0):,},地板按 √n 该涨到约 "
      f"{F856*np.sqrt(N856/max(len(pe0),1)):.4f}** —— 实测 {FLOOR:.4f}")
PLANT = 0.20
pu0 = zs(base.f_pub).to_numpy(float); pr0 = zs(-base.relpersn).to_numpy(float)
d_base, bs0 = dfacet(pe0, pu0, pr0)
d_pl, _ = dfacet(pe0, pu0, pr0, Bv=300, plant=PLANT)
sd0 = float(np.std(bs0))
got = d_pl - d_base
power_ok = bool(abs(got) / sd0 > 2)
print(f"  植入一个 Δfacet(只按公共脸加 +{PLANT})⇒ 动 **{got:+.4f}** = **{abs(got)/sd0:.2f} 个自助 SD** "
      f"⇒ {'**有功效**' if power_ok else '**没功效 ⇒ 不看结果**'}")

print(f"\n=== ① 同年同人,三张私人脸并排(B={B})===")
rows = []
for name, series in FACES.items():
    m = D.perm.notna() & D.f_pub.notna() & series.notna()
    sub = D[m]; sv = series[m]
    if m.sum() < 400:
        print(f"  {name:22s} **n={int(m.sum())} 太小,跳过**"); continue
    pe = sub.perm.to_numpy(float); pu = zs(sub.f_pub).to_numpy(float); pr = zs(sv).to_numpy(float)
    obs, bs = dfacet(pe, pu, pr)
    lo, hi = np.quantile(bs, [.025, .975])
    j = jac(pu, pr)
    rows.append(dict(face=name, n=int(m.sum()), g_pub=gap(pe, pu), g_priv=gap(pe, pr),
                     dfacet=obs, lo=float(lo), hi=float(hi), jaccard=j,
                     excl=bool(lo > 0 or hi < 0), above=bool(abs(obs) > FLOOR)))
    print(f"  {name:22s} n={int(m.sum()):>5,} · g_公共 **{rows[-1]['g_pub']:+.3f}** · "
          f"g_私人 **{rows[-1]['g_priv']:+.3f}** · **Δfacet {obs:+.3f}** [{lo:+.3f},{hi:+.3f}] · "
          f"Jaccard {j:.3f} · {'**超地板**' if rows[-1]['above'] else '地板内'}")
print(f"  ⚠ **`attend` 在 {YEAR} 的 g_公共 印在上面** —— 请与 `#856` 的 2010s 的 −0.85 对照看,"
      f"**但不要并列解读**:2021 是 GSS 改网络/混合模式的一年,**模式变化对年内比较是共同的,"
      f"对跨年比较是致命的。**")

print("\n=== ② 控制 ===")
rg = np.random.default_rng(7)
nc = gap(pe0, pu0, rng=rg, permute=True) - gap(pe0, pr0, rng=rg, permute=True)
print(f"  正控:见 ⓪b(植入 +{PLANT} ⇒ 动 {got:+.4f});**plant=0 时动 "
      f"{dfacet(pe0, pu0, pr0, Bv=50, plant=0.0)[0]-d_base:+.6f}** —— ⚠ **`G2` 要求控制必须能失败**")
print(f"  负控:打乱虔诚度 ⇒ Δfacet = **{nc:+.4f}**"
      f"(⚠ **「这个零该不该是零?」该** —— 打乱后两张脸各自的缝期望都是 0,差也是 0)")
print(f"  安慰剂地板 = **{FLOOR:.4f}**(两个随机分层器,300 次)")

imp = [r for r in rows if r["face"].startswith("重要性")]
oth = [r for r in rows if not r["face"].startswith("重要性")]
nondeg = [r for r in rows if r["jaccard"] < 0.80]
imp_small = bool(imp) and all(not r["above"] for r in imp)
oth_big = bool(oth) and any(r["above"] for r in oth)
all_big_same = bool(rows) and all(r["above"] for r in rows) and len({np.sign(r["dfacet"]) for r in rows}) == 1

Gt = Gate("#857 · 同一年同一批人,把三张私人脸并排放")
Gt.asserted("⓪ **功效闸(排在所有关于人的判断之前)**:n 从 7,341 掉到约 2,300,地板按 √n 会涨,"
            "**而 `#856` 观测到的 Δfacet 是 0.136–0.234,可能整段落进地板** ⇒ "
            "先植入再看,捞不回 ≥2 个自助 SD 就不看结果",
            power_ok, f"植入动 {got:+.4f} = {abs(got)/sd0:.2f} 个 SD · 地板 {FLOOR:.4f}", kind="control")
Gt.asserted("① 硬规则①:`religimp`/`relidimp` **整个 GSS 里只问过 2021 一年**,`godguide` 只问过 2004;"
            "而 `relpersn`/`pray` **也在 2021** ⇒ **三张私人脸同年同人并排,仪器/年份/人全固定,只有脸在变**",
            bool(len(rows) >= 3), f"进入分析的脸 {[r['face'] for r in rows]}", kind="control")
Gt.asserted("② 算术陷阱:两张脸切出的三分位若几乎相同,`Δfacet` 被代数逼近 0 ⇒ 逐格印 Jaccard,"
            "**≥0.80 记退化不计入**",
            bool(len(nondeg) == len(rows)),
            " · ".join(f"{r['face'][:6]}:{r['jaccard']:.3f}" for r in rows), kind="control")
Gt.asserted("③ 负控:打乱虔诚度 ⇒ Δfacet 必须 ≈0(⚠ **这个零该是零**)",
            bool(abs(nc) < 0.15), f"{nc:+.4f}(地板 {FLOOR:.4f})", kind="control")
Gt.asserted("④ 前提(跑前写下的最强混淆):**2021 是 GSS 改网络/混合模式的一年** ⇒ "
            "**本轮只做年内比较,绝不把 2021 的 `g` 与 `#856` 的 2010s `g` 并列解读**;"
            "并把 2021 的 `g_公共` 印出来让读者自己看差多少",
            bool(all(np.isfinite(r["g_pub"]) for r in rows)),
            f"2021 g_公共 = {rows[0]['g_pub']:+.3f}(2010s 为 −0.85,**不并列解读**)", kind="control")
Gt.asserted("⑤ kill(预注册):「那个差是**脸**造成的、不是仪器」要成立,需 "
            "**「重要性」两张全部落在地板内,而身份/实践至少一张超地板**",
            bool(imp_small and oth_big),
            f"重要性超地板 {sum(r['above'] for r in imp)}/{len(imp)} · "
            f"其余超地板 {sum(r['above'] for r in oth)}/{len(oth)}", kind="kill",
            yardstick="每张脸的 `Δfacet` 对照本年安慰剂地板(两个随机分层器的 95 分位)",
            yardstick_noise=FLOOR,
            population=f"GSS {YEAR} 年内、与 `homosex` 同时非缺、且 Jaccard<0.80 的 {len(nondeg)} 张脸 —— "
                       f"⚠ **不含任何跨年或跨仪器的格**")
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = ("**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
            + ("\n  ⚠ **而没过的是 ⓪ 功效闸** —— 本年 n 太小,"
               "**这不是「三张脸一样」,是「这一年分辨不了」。**" if not power_ok else ""))
elif imp_small and oth_big:
    VERD = (f"**A 是脸,不是仪器 —— `#856` 那条共线拆开了。**\n"
            f"  同一年、同一批人:**「重要性」两张脸的 Δfacet 都落在地板 {FLOOR:.3f} 之内**"
            f"({[round(r['dfacet'],3) for r in imp]}),\n"
            f"  **而身份/实践至少一张超地板**({[round(r['dfacet'],3) for r in oth]})。\n"
            f"  ⇒ **一句关于人的话:「你多久去一次礼拜」和「你自己算不算个有信仰的人」"
            f"把人分得不一样开;\n"
            f"  而「宗教对你有多重要」分出来的距离,和去不去礼拜分出来的一样宽。\n"
            f"  所以 NSFG 里那个「没有差」不是换了问卷的缘故,是换了那张脸。**")
elif all_big_same:
    VERD = (f"**B 是仪器,不是脸。** 三张私人脸的 Δfacet **全部超地板 {FLOOR:.3f} 且同号**"
            f"({[round(r['dfacet'],3) for r in rows]})\n"
            f"  ⇒ **在 GSS 里换哪张私人脸都一样** ⇒ NSFG 的零就是仪器的事,`#856` 的共线归因给仪器。")
else:
    VERD = (f"**C 「公共 / 私人」这个二分不是切开这件事的正确关节 —— 而这是元分离器。**\n"
            f"  三张私人脸彼此并不一致,也没有「重要性 vs 其余」的清晰结构:"
            f"{[(r['face'][:8], round(r['dfacet'],3), '超' if r['above'] else '内') for r in rows]}\n"
            f"  ⇒ **不是某个世界赢了,是我的分类方式在这批数据上不成立。**\n"
            f"  ⇒ **一句关于人的话:我把虔诚分成「公共的」和「私人的」两张脸,\n"
            f"  而同一年同一批人里,这几张脸给出的距离并不按这条线排队 ——\n"
            f"  「去礼拜 / 自认信仰 / 祷告 / 重要性」不是两组,是四把各不相同的尺。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的**:① **「重要性」这张脸在 GSS 里只有 {YEAR} 一年** ⇒ "
      f"**不能做跨年稳健性,本轮是单年结果**;② 2021 的模式变化**对年内比较共同、对跨年比较致命** ⇒ "
      f"**跨年比较本轮不做,不是「以后做」**;③ 横断面 ⇒ **无干预、无因果识别**。")
json.dump(dict(year=YEAR, rows=rows, floor=FLOOR, power=dict(plant=PLANT, moved=got,
               in_sd=abs(got)/sd0, ok=power_ok), neg_control=nc,
               structurally_impossible=["religimp/relidimp asked only in 2021 -> no cross-year robustness",
                                        "GSS 2021 changed to web/mixed mode -> cross-year comparison not done",
                                        "cross-sectional -> no causal identification"],
               admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "all_faces_one_year.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'all_faces_one_year.json'}")
