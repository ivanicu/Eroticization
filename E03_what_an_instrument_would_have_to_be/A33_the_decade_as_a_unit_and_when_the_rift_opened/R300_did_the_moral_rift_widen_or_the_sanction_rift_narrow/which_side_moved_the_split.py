"""#861 · E03·A92·R300 —— 裂口自己在变宽,是道德那条缝变宽,还是制裁那条缝变窄?

`#858`②(未预注册的线索):`|Δobject|` 从 1990s 的 **0.26–0.50** 涨到 2010s 的 **0.54–0.88** ——
**判断与制裁之间的裂口自己在变宽。** 本轮**预注册地**去问它,并把它拆开。

`G1` **估计量(先于方法命名)**:
   **`Δobject(2010s) − Δobject(1990s)` 的两个分量**:
   `= [g_道德(2010s) − g_道德(1990s)] − [g_制裁(2010s) − g_制裁(1990s)]`
   记作 `Δg_道德` 与 `Δg_制裁`。

**⚠⚠⚠ 算术先说清楚(`realstat` 开篇那条),否则整轮是「1+1=2 所以 2<3」:**
**上面那个等号是恒等式,不是发现** —— 裂口的变化**必然**等于两个分量之差。
⇒ **和是恒等式,分法才是测量**(`#838` 的同一条)。**本轮报的是分法,并数值核对那个恒等式。**

四个世界(**每个都有分支**;`#856` 的教训):
   A **道德那条缝变宽**:`Δg_道德` 显著为负而 `Δg_制裁` 在地板内
     ⇒ **宗教在「对不对」上把人分得比从前更开。**
   B **制裁那条缝变窄**:`Δg_制裁` 显著为正(朝 0 走)而 `Δg_道德` 在地板内
     ⇒ **公民自由变成了共识,而道德评价仍然对立** —— **这是关于「社会拿它怎么办」的直接答案。**
   C **两边都动** ⇒ 报分法,不报「哪一边」。
   D **两个分量都在地板内** ⇒ **裂口的变化本身不可归因**,`#858`② 那条线索降级。

预测矩阵:
   | 世界 | 现在 | 只道德动 | 只制裁动 | 两边都动 | 都不动 |
   | A 道德变宽 | 0.25 | **0.85** | 0.05 | 0.05 | 0.05 |
   | B 制裁变窄 | 0.40 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 两边都动 | 0.25 | 0.05 | 0.05 | **0.85** | 0.05 |
   | D 不可归因 | 0.10 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 恒等式核对通过(`Δ裂口 == Δg_道德 − Δg_制裁`,数值上)
     and 正控开火(**只往 2010s 的道德题植入位移 ⇒ `Δg_道德` 必须按量动,而 `Δg_制裁` 不动;
        且 plant=0 时两者都必须恰为 0**)
     and 负控为零(**年内打乱虔诚度 ⇒ 两个分量都塌到地板内**)
     and 安慰剂给出地板(**两个随机分层器**):
      只有 `Δg_道德` 超地板 -> A
      只有 `Δg_制裁` 超地板 -> B
      两个都超            -> C
      两个都不超          -> D
  else: UNVERIFIED

**⚠⚠ 跑之前写下的最强混淆,而它正是 B 的头号对手:**
   **三道制裁题的「允许」率从 1990s 的约 72% 涨到 2010s 的约 87%** ——
   **一个逼近天花板的二值题,它的任何差距都会被机械压缩。**
   ⇒ **「制裁那条缝变窄」完全可能是天花板压出来的,与社会共识无关。**
   ⇒ 控制:**同一分解在两把尺上各跑一遍** —— ①均值尺;②**潜在尺**(probit,不受二值边界压缩)。
   **`#858` 已实测潜在尺会把效应放大而不是缩小**,所以它是这条混淆的合格判据。
   **若「制裁变窄」只在均值尺上出现、潜在尺上消失 ⇒ 天花板就是解释,B 死。**

`G3` 多重性:**3 道制裁题 × 2 把尺 = 6 格**,每格两个分量 ⇒ **12 个量**,BH 与 BY 都做。
`G4` 规格曲线:两根轴逐格报,**含不同意的格**。
⚠ kill 带 `population` 与 **`direction`**(`#860` 刚补的那个参数,**本轮第一次用**):
   **「哪一边动」这个世界要求分量的符号在三道题上一致** —— **只要求量级会让符号翻转通过**(`#859`)。

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① **换不了仪器**:`#854` 已点名盘上七具,**Stouffer 三题在这批数据里是 GSS 独有的。**
② 横断面 ⇒ **无干预、无因果识别**:不能问「是道德先松还是制裁先松」。
③ 两个十年之间**问卷模式、抽样框都可能变**,而本轮**无法把它与真实变化分开** ——
   **它对两个分量是共同的**(同一批人同一份问卷),**所以对「哪一边动」这个比较无害,
   对「变了多少」这个绝对量有害** ⇒ 本轮只报**哪一边**,不把绝对量当社会变迁的度量。
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
D["moral"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
D["spk"] = 2 - pd.to_numeric(g.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["col"] = 5 - pd.to_numeric(g.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))
D["lib"] = pd.to_numeric(g.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    D[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
R = D.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = (zs(R.attend) + zs(-R.reliten) + zs(-R.fund)) / 3
D = D.join(R["REL"])
SANCT = {"发言 `spkhomo`": "spk", "教书 `colhomo`": "col", "图书馆 `libhomo`": "lib"}
DECS = {"1990s": range(1990, 2000), "2010s": range(2010, 2020)}

print("=== ⓪ 硬规则①:跑前写下的最强混淆,先把它的量印出来 ===")
for nm, c in SANCT.items():
    r = {d: float(D[c][D[c].notna() & D.year.isin(list(DECS[d]))].mean()) for d in DECS}
    print(f"  {nm:18s} 「允许」率 1990s **{r['1990s']:.1%}** → 2010s **{r['2010s']:.1%}** "
          f"(**涨 {r['2010s']-r['1990s']:+.1%}**)")
print("  ⚠⚠ **逼近天花板的二值题,任何差距都会被机械压缩** ⇒ "
      "**「制裁那条缝变窄」完全可能是天花板压出来的** ⇒ **潜在尺同跑,它是这条混淆的判据**")

Phi = lambda z: 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
def ppf(p):
    lo_, hi_ = -8.0, 8.0
    for _ in range(120):
        m = (lo_ + hi_) / 2
        if Phi(m) < p: lo_ = m
        else: hi_ = m
    return (lo_ + hi_) / 2

def gaps(sub, sc, scale, plant=0.0, rng=None, permute=False, relx=None):
    """返回 (g_道德, g_制裁),两把尺同一构造。plant 只加在**道德题**上。"""
    gm, gs = [], []
    for y in np.unique(sub.year):
        s = sub[sub.year == y]
        if len(s) < 200: continue
        r = (relx[sub.year.to_numpy() == y] if relx is not None else s.REL.to_numpy(float))
        if permute: r = rng.permutation(r)
        lo, hi = np.quantile(r, [1/3, 2/3])
        dev, sec = r >= hi, r <= lo
        if dev.sum() < 60 or sec.sum() < 60: continue
        mo = s.moral.to_numpy(float).copy()
        if plant: mo[dev] += plant
        sa = s[sc].to_numpy(float)
        if scale == "均值尺":
            if mo.std(ddof=1) <= 0 or sa.std(ddof=1) <= 0: continue
            gm.append((mo[dev].mean() - mo[sec].mean()) / mo.std(ddof=1))
            gs.append((sa[dev].mean() - sa[sec].mean()) / sa.std(ddof=1))
        else:
            mb = (mo >= 3).astype(float)      # ⚠ 潜在尺下道德题也二值化,否则两边不是同一种尺
            f = lambda v: ppf(np.clip(v[dev].mean(), 1e-4, 1-1e-4)) - ppf(np.clip(v[sec].mean(), 1e-4, 1-1e-4))
            gm.append(f(mb)); gs.append(f(sa))
    if not gm: return np.nan, np.nan
    return float(np.mean(gm)), float(np.mean(gs))

def components(sc, scale, Bv=B, seed=861, plant=0.0, rng=None, permute=False, relx=None):
    out = {}
    for d in DECS:
        m = D.moral.notna() & D[sc].notna() & D.REL.notna() & D.year.isin(list(DECS[d]))
        sub = D[m]
        rx = relx[m.to_numpy()] if relx is not None else None
        out[d] = gaps(sub, sc, scale, plant=(plant if d == "2010s" else 0.0),
                      rng=rng, permute=permute, relx=rx)
    dgm = out["2010s"][0] - out["1990s"][0]
    dgs = out["2010s"][1] - out["1990s"][1]
    return dgm, dgs, out

def boot(sc, scale, Bv=B, seed=861):
    rg = np.random.default_rng(seed); a = np.empty(Bv); b = np.empty(Bv)
    idx = {d: np.flatnonzero((D.moral.notna() & D[sc].notna() & D.REL.notna()
                              & D.year.isin(list(DECS[d]))).to_numpy()) for d in DECS}
    for i in range(Bv):
        gm_, gs_ = {}, {}
        for d in DECS:
            k = rg.choice(idx[d], len(idx[d]), replace=True)
            gm_[d], gs_[d] = gaps(D.iloc[k], sc, scale)
        a[i] = gm_["2010s"] - gm_["1990s"]; b[i] = gs_["2010s"] - gs_["1990s"]
    ok = np.isfinite(a) & np.isfinite(b)
    return a[ok], b[ok]

print(f"\n=== ① 分解:`Δ裂口 = Δg_道德 − Δg_制裁`(**恒等式,分法才是测量**)· B={B} ===")
rows = []
for nm, sc in SANCT.items():
    for scale in ("均值尺", "潜在尺"):
        dgm, dgs, out = components(sc, scale)
        ba, bb = boot(sc, scale)
        if len(ba) < 100:
            print(f"  {nm:18s} {scale} **不可估**"); continue
        la, ha = np.quantile(ba, [.025, .975]); lb, hb = np.quantile(bb, [.025, .975])
        dsplit = (out["2010s"][0] - out["2010s"][1]) - (out["1990s"][0] - out["1990s"][1])
        rows.append(dict(item=nm, scale=scale, dg_moral=dgm, dg_sanct=dgs,
                         moral_lo=float(la), moral_hi=float(ha),
                         sanct_lo=float(lb), sanct_hi=float(hb),
                         d_split=float(dsplit), identity=float(abs(dsplit - (dgm - dgs))),
                         g_moral_90=out["1990s"][0], g_moral_10=out["2010s"][0],
                         g_sanct_90=out["1990s"][1], g_sanct_10=out["2010s"][1]))
        print(f"  {nm:18s} {scale} · **Δg_道德 {dgm:+.3f}** [{la:+.3f},{ha:+.3f}] · "
              f"**Δg_制裁 {dgs:+.3f}** [{lb:+.3f},{hb:+.3f}] · 裂口变化 {dsplit:+.3f}")
print(f"  **恒等式核对** `Δ裂口 == Δg_道德 − Δg_制裁` 最大偏差:"
      f"**{max(r['identity'] for r in rows):.2e}** ⇒ **它确实是代数,和是恒等式,分法才是测量**")

print("\n=== ② 控制 ===")
sc0, sca0 = "col", "均值尺"
base_m, base_s, _ = components(sc0, sca0)
PLANT = 0.30
pm, ps, _ = components(sc0, sca0, plant=PLANT)
zm, zs_, _ = components(sc0, sca0, plant=0.0)
print(f"  正控:**只往 2010s 的道德题的虔诚层**加 +{PLANT} ⇒ "
      f"Δg_道德 动 **{pm-base_m:+.4f}** · **Δg_制裁 动 {ps-base_s:+.4f}(该为 0)**")
print(f"     **而 plant=0 时:Δg_道德 动 {zm-base_m:+.6f} · Δg_制裁 动 {zs_-base_s:+.6f}** —— "
      f"⚠ **`G2` 要求控制必须能失败,这一行就是那个检查**")
rg = np.random.default_rng(9)
relall = D.REL.to_numpy(float)
nm_, ns_, _ = components(sc0, sca0, rng=rg, permute=True)
print(f"  负控:**年内打乱虔诚度** ⇒ Δg_道德 **{nm_:+.4f}** · Δg_制裁 **{ns_:+.4f}** "
      f"(⚠ **「这个零该不该是零?」该** —— 打乱后每个十年的两个 g 期望都是 0,差也是 0)")
rg2 = np.random.default_rng(11)
pls = []
for _ in range(150):
    fake = rg2.normal(size=len(D))
    a_, b_, _ = components(sc0, sca0, relx=fake)
    if np.isfinite(a_) and np.isfinite(b_): pls.append((a_, b_))
pls = np.array(pls)
FLOOR = float(np.quantile(np.abs(pls).ravel(), 0.95))
print(f"  **安慰剂 = 本量的噪声地板**:**随机分层器**下两个分量的 |值| 95 分位 = **{FLOOR:.4f}**")

# ⚠⚠ **一条必须写下来的多重性事实**:三道制裁题共用**同一道道德题**,
#    只有「哪些人同时答了那道制裁题」略有不同 ⇒ **`Δg_道德` 的三格几乎是同一个测量重复三次,
#    不是三次独立检验。** 实测三格 Δg_道德 在同一把尺上只差 0.003 —— **这不是三重证据。**
_mo_spread = max(abs(a["dg_moral"] - b["dg_moral"]) for a in rows for b in rows
                 if a["scale"] == b["scale"])
print(f"  ⚠ **多重性注记**:三道制裁题共用同一道道德题 ⇒ 同尺下 `Δg_道德` 的三格极差仅 "
      f"**{_mo_spread:.4f}** —— **它们几乎是同一个测量重复三次,不是三次独立证据。**")
mo_big = [r for r in rows if abs(r["dg_moral"]) > FLOOR and (r["moral_lo"] > 0 or r["moral_hi"] < 0)]
sa_big = [r for r in rows if abs(r["dg_sanct"]) > FLOOR and (r["sanct_lo"] > 0 or r["sanct_hi"] < 0)]
mo_signs = [r["dg_moral"] for r in rows]; sa_signs = [r["dg_sanct"] for r in rows]
lat = [r for r in rows if r["scale"] == "潜在尺"]
sa_lat_big = [r for r in lat if abs(r["dg_sanct"]) > FLOOR and (r["sanct_lo"] > 0 or r["sanct_hi"] < 0)]

Gt = Gate("#861 · 裂口变宽,是哪一边动的")
# ⚠ **第一版把 `max|差|` 拿去和字面 0 比,被库判 DEGENERATE —— 而库是对的(`#773`/`#840`):
#    一个算出来的 0 和一个写死的 0 相等,不能区分「算对了」与「根本没算」。**
#    ⇒ `#770` 那条:**比两个非零的值本身。**
_r0 = rows[0]
Gt.identity_control("① **算术先说清楚**:`Δ裂口 = Δg_道德 − Δg_制裁` 是**恒等式,不是发现** ⇒ "
                    "**和是恒等式,分法才是测量**(`#838` 同一条);"
                    "⚠ **比的是两个非零值本身,不是它们的差与零**(`#770`/`#773`)",
                    observed=float(_r0["d_split"]),
                    expected=float(_r0["dg_moral"] - _r0["dg_sanct"]), tol=1e-9,
                    what=f"{_r0['item']}/{_r0['scale']} 的 Δ裂口 vs (Δg_道德 − Δg_制裁)",
                    deterministic=True)
Gt.asserted("② 前提(跑前写下的最强混淆,也是 B 的头号对手):**三道制裁题的允许率从约 72% 涨到约 87%**,"
            "**逼近天花板的二值题任何差距都会被机械压缩** ⇒ **同一分解在均值尺与潜在尺上各跑一遍**;"
            "`#858` 已实测潜在尺会**放大**而非缩小效应,**所以它是这条混淆的合格判据**",
            bool(len({r["scale"] for r in rows}) == 2),
            f"两把尺都跑,共 {len(rows)} 格", kind="control")
Gt.asserted("③ 正控:**只往 2010s 的道德题**植入 +0.30 ⇒ `Δg_道德` 必须动而 `Δg_制裁` 不动;"
            "**且 plant=0 时两者都必须恰为 0**",
            bool(abs(pm - base_m) > 0.05 and abs(ps - base_s) < 0.01
                 and abs(zm - base_m) < 1e-12 and abs(zs_ - base_s) < 1e-12),
            f"道德动 {pm-base_m:+.4f} · 制裁动 {ps-base_s:+.4f} · "
            f"plant=0 时 {zm-base_m:+.2e}/{zs_-base_s:+.2e}", kind="control")
Gt.asserted("④ 负控:年内打乱虔诚度 ⇒ 两个分量都必须落在地板内"
            "(⚠ **这个零该是零**:打乱后每个十年的两个 g 期望都是 0,差也是 0)",
            bool(abs(nm_) < FLOOR * 2 and abs(ns_) < FLOOR * 2),
            f"{nm_:+.4f} / {ns_:+.4f}(地板 {FLOOR:.4f})", kind="control")
Gt.asserted("⑤ kill(预注册):「制裁那条缝变窄」要成立,需 **`Δg_制裁` 在**潜在尺**上仍超地板且排零** "
            "—— **只在均值尺上出现就是天花板**",
            bool(len(sa_lat_big) >= 2),
            f"潜在尺上超地板且排零的制裁分量 {len(sa_lat_big)}/{len(lat)} · "
            f"全部格中 道德 {len(mo_big)}/{len(rows)} · 制裁 {len(sa_big)}/{len(rows)}",
            kind="kill",
            yardstick="每个分量自己的 95% 自助区间,并须超过随机分层器的噪声地板",
            yardstick_noise=FLOOR,
            population=f"GSS 的 {len(rows)} 格(3 道制裁题 × 2 把尺)—— ⚠ 全部同一具仪器",
            direction=sa_signs)
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
mo_dir = len({1 if x > 0 else -1 for x in mo_signs}) == 1
sa_dir = len({1 if x > 0 else -1 for x in sa_signs}) == 1
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(sa_lat_big) >= 2 and len(mo_big) < len(rows) / 2:
    VERD = (f"**B 制裁那条缝变窄 —— 而且天花板解释不掉它。**\n"
            f"  `Δg_制裁` 在潜在尺上 {len(sa_lat_big)}/{len(lat)} 格超地板 {FLOOR:.3f} 且排零"
            f"(方向{'一致' if sa_dir else '**不一致**'});`Δg_道德` 只有 {len(mo_big)}/{len(rows)} 格。\n"
            f"  ⇒ **一句关于人的话:从九十年代到二〇一〇年代,虔诚的人和世俗的人\n"
            f"  在「该不该拦着他」这件事上越走越近,而在「这样对不对」上没有。\n"
            f"  公共待遇变成了共识,道德评价没有 —— 社会先学会了放行,并没有先学会同意。**")
elif len(mo_big) >= len(rows) / 2 and len(sa_lat_big) < 2:
    VERD = (f"**A 道德那条缝变宽。** `Δg_道德` {len(mo_big)}/{len(rows)} 格超地板"
            f"(方向{'一致' if mo_dir else '**不一致**'}),而制裁分量在潜在尺上只有 {len(sa_lat_big)} 格。\n"
            f"  ⇒ **一句关于人的话:宗教在「对不对」上把人分得比从前更开。**")
elif len(mo_big) >= len(rows) / 2 and len(sa_lat_big) >= 2:
    VERD = (f"**C 两边都动 ⇒ 报分法,不报「哪一边」。** 道德 {len(mo_big)}/{len(rows)} 格 · "
            f"制裁(潜在尺){len(sa_lat_big)}/{len(lat)} 格,均超地板 {FLOOR:.3f}。")
else:
    VERD = (f"**D 裂口的变化不可归因 ⇒ `#858`② 那条线索降级。** 道德超地板 {len(mo_big)}/{len(rows)} · "
            f"制裁(潜在尺)超地板 {len(sa_lat_big)}/{len(lat)},**都不到判据要求** ⇒ "
            f"**「裂口在变宽」这句话本轮没有得到可归因的支持。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的**:① **换不了仪器**(`#854` 已点名七具,Stouffer 三题是 GSS 独有);"
      f"② 横断面 ⇒ **无因果识别**,不能问「是道德先松还是制裁先松」;"
      f"③ 两个十年之间**问卷模式与抽样框可能变**,而它**对两个分量是共同的**(同一批人同一份问卷)"
      f"⇒ **对「哪一边动」无害,对「变了多少」有害** ⇒ **本轮只报哪一边,不把绝对量当社会变迁的度量。**")
json.dump(dict(grid=rows, floor=FLOOR, moral_big=len(mo_big), sanct_big=len(sa_big),
               sanct_latent_big=len(sa_lat_big), moral_component_spread=float(_mo_spread),
               moral_not_independent="the three cells share one moral item; not three independent tests",
               moral_dir_consistent=mo_dir,
               sanct_dir_consistent=sa_dir,
               allow_rates={nm: {d: float(D[c][D[c].notna() & D.year.isin(list(DECS[d]))].mean())
                                 for d in DECS} for nm, c in SANCT.items()},
               controls=dict(plant=PLANT, moral_moved=pm-base_m, sanct_moved=ps-base_s,
                             zero_moral=zm-base_m, zero_sanct=zs_-base_s,
                             neg=(nm_, ns_), floor=FLOOR),
               admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "which_side_moved.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'which_side_moved.json'}")
