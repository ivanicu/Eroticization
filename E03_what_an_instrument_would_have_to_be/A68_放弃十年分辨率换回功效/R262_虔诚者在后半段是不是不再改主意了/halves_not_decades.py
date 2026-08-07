"""#823 · E03·A68·R262 —— 放弃十年分辨率,换回功效:虔诚者在后半段是不是不再改主意了?

`#818` 问「虔诚者在两千年代是停住了还是往回走了」,判 UNVERIFIED,理由是
**这个设计的噪声底 0.134 比它要分辨的 0.10 还大**。
`#822` 证明**那不是统计口径的问题,是 GSS 每个十年只跑四到六次** ——
**年份就是样本,而在十年这个分辨率上样本量是固定的、不可增加的。**

⚠⚠ **所以 `#822`① 的处方不是「重新跑」,是「重新提问」:**
   **问一个不需要十年分辨率的版本。** 把五十年切成**两半**而不是五个十年:
   每半段有 **15–20 个年份** 而不是 4–6 个 ⇒ **年际噪声按 √n 降**,
   预期从 0.134 降到 **≈0.08**,**而 0.08 < `#818` 要的 0.10。**
   ⚠ **它放弃的正是我已经证明拿不到的那部分(哪一个十年不寻常),
   换回的是我确实需要的那部分(前后两半有没有不同)。**

G1 估计量:**虔诚层自己的、世代内的态度变化**(`Σ_c w̄_c·(m_c1 − m_c0)`,与 `#818` 同一个量),
   **分别在前半段与后半段上**,并**报两者之差**。
   ⚠ **用世代内项而不是原始均值差** —— `#816` 已证原始均值差里混着世代替换。

**分段(跑前写死,不看结果调):** 前半 = 首年 → 1989 · 后半 = 1990 → 末年。
   **1990 这个切点不是挑出来的**:`#812` 的十年网格已独立地把 1990s 定为 `homosex` 的第一个偏离十年,
   **而那是在本轮之前、用另一个估计量得到的。⚠ 但它仍然是从同一份数据来的 ——
   如实登记:这个切点不是外生的,它借自本项目自己的一个先前结果。**

⚠⚠ **本轮第一个可失败的主张是关于仪器,不是关于人**(`#822` 的教训):
   **先验证噪声底真的降到了 0.10 以下,再读结论。若没降下来,本轮到此为止,不许读那个差。**

三个世界:
   A **后半段虔诚者不再改主意**:后半段的世代内变化含 0 且窄(`TIGHT_NULL`)⇒
     **一句关于人的话:虔诚的美国人在前半段跟着社会一起变宽容,在后半段停下了。**
   B **后半段往回走**:区间排除 0 且为负 ⇒ **他们变得更保守了**,比 A 强得多。
   C **两半段没有差别**:两段的差含 0 且窄 ⇒ **`#817` 那个「+0.166 → −0.063」是十年切法的产物**,
     **而那会撤掉 `#817` 最出彩的一句话 —— 这是我不欢迎的那个结果(盆地规则)。**

预测矩阵:
   | 世界 | 现在 | 后半窄零 | 后半排除 0 且负 | 两半之差含 0 且窄 |
   | A 停住 | 0.40 | **0.80** | 0.05 | 0.25 |
   | B 回走 | 0.20 | 0.05 | **0.90** | 0.05 |
   | C 无差别 | 0.40 | 0.20 | 0.03 | **0.85** |

预注册判词(条件式,**两级**):
  if 噪声底(按年份聚类)**未**降到 0.10 以下:
      verdict = UNVERIFIED —— **本轮的前提就是它降下来了,没降就没有资格往下读**
  elif 正控开火(植入一个已知的半段内位移必须取回)
       and 负控开火(世代内态度按年份重抽 ⇒ 真值 0 且带噪声,容差 **0.05 事先写死**):
      两半段各自用 `#811` 三值判(参照 0,`matters` = 0.10),**并报两半之差**
  else: UNVERIFIED
⚠ 凡 `UNRESOLVED` 必须同时印出它与哪些参照相容(`#812`③)。

⚠ 跑之前写下的最强混淆(与 `#818` 同一条,**重申因为它没被解决**):
  虔诚层按**年份内三分位**定义 ⇒ 两个半段的「虔诚层」不是同一群人。
  ⇒ 控制:**两段各自的 `REL` 均值与四分位并排印出**;**本设计修不了它,只能量它。**

⚠ 本轮**换不了仪器**;硬规则②:分层与题目都来自 GSS 同一份问卷。
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(262)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, KK, CUT = "homosex", 4, 1990
MATTERS, NC_TOL, B, NREP = 0.10, 0.05, 2000, 200

d = pd.read_stata(gp, columns=["year", "cohort", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("cohort", (1880, 2010))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
D = REL.dropna(subset=[IT, "cohort"]).copy()
D["gen"] = (D.cohort//10*10).astype(int)
D = D[D.k == 2]
YRS = sorted({int(y) for y, g in D.groupby("year") if len(g) >= 120})
EARLY = [y for y in YRS if y < CUT]; LATE = [y for y in YRS if y >= CUT]
keep = [g for g, n in D[D.year.isin(YRS)].groupby("gen").size().items() if n >= 60]
D = D[D.gen.isin(keep) & D.year.isin(YRS)]

print(f"=== ⓪ 硬规则①:虔诚层 · 两个半段的年份数(⚠ 这正是本轮想买回来的东西)===")
print(f"  前半 {EARLY[0]}–{EARLY[-1]}:**{len(EARLY)} 个年份** {EARLY}")
print(f"  后半 {LATE[0]}–{LATE[-1]}:**{len(LATE)} 个年份** {LATE}")
print(f"  ⇒ 对比 `#818` 每个十年只有 4–6 个年份 ⇒ **年际噪声预期按 √n 降到约 {0.134/np.sqrt(len(LATE)/5):.3f}**")
print(f"  世代箱 {len(keep)} 个:{sorted(keep)}")

def within(df, ys):
    """该层内:Σ w̄_c·(m_c末 − m_c首),用该半段的首末年;⚠ 年份重抽时首末年会变。"""
    a, b = df[df.year == ys[0]], df[df.year == ys[-1]]
    gens = sorted(set(a.gen) & set(b.gen))
    if not gens or not len(a) or not len(b): return None
    w0 = np.array([len(a[a.gen == g])/len(a) for g in gens])
    w1 = np.array([len(b[b.gen == g])/len(b) for g in gens])
    m0 = np.array([a[a.gen == g][IT].mean() for g in gens])
    m1 = np.array([b[b.gen == g][IT].mean() for g in gens])
    if np.isnan(m0).any() or np.isnan(m1).any(): return None
    return float((((w0+w1)/2)*(m1-m0)).sum())
def within_fit(df, ys):
    """⚠ `#822`:多用年份在 `cluster` 下不降噪 —— 所以这里仍用首末年,
    **本轮买到的功效来自「半段的年份多 ⇒ 重抽时首末年的组合多」,不是来自拟合。**"""
    return within(df, ys)

def boot_cluster(df, ys, rep=B):
    out = np.empty(rep)
    for i in range(rep):
        idx = RNG.integers(0, len(ys), len(ys))
        use = sorted({ys[j] for j in idx})
        out[i] = np.nan if len(use) < 2 else (within(df, use) or np.nan)
    o = out[np.isfinite(out)]
    return float(np.percentile(o, 2.5)), float(np.percentile(o, 97.5))

print(f"\n=== ① 先验证仪器:噪声底真的降到 0.10 以下了吗?(⚠ `#822` 的教训:先验仪器,再读人)===")
RES = {}
for lab, ys in (("前半", EARLY), ("后半", LATE)):
    pt = within(D, ys); lo, hi = boot_cluster(D, ys)
    RES[lab] = dict(years=ys, n_years=len(ys), point=float(pt), lo=lo, hi=hi, half=(hi-lo)/2)
    print(f"  {lab} {ys[0]}–{ys[-1]}({len(ys)} 年):**{pt:+.4f}** [{lo:+.4f}, {hi:+.4f}] · "
          f"噪声半宽 **{(hi-lo)/2:.4f}**")
floor_ok = all(RES[l]["half"] < MATTERS for l in RES)
print(f"  ⇒ 两段的噪声半宽都 < {MATTERS}?**{floor_ok}** —— "
      f"{'**前提成立,可以往下读**' if floor_ok else '**前提不成立,本轮到此为止**'}")

print(f"\n  ⚠ 跑前混淆的控制 —— 两段各自的 `REL` 分布(层的主语有没有漂):")
for lab, ys in (("前半", EARLY), ("后半", LATE)):
    s = D[D.year.isin(ys)].REL
    print(f"    {lab}:均值 {s.mean():+.3f} · Q1 {s.quantile(.25):+.3f} · Q3 {s.quantile(.75):+.3f}")

print("\n=== ② 控制 ===")
def syn(mode, ys, shift=-0.30):
    H = D[D.year.isin(ys)].copy()
    gm = H[H.year == ys[0]].groupby("gen")[IT].mean()
    if mode == "planted":
        m = H.year == ys[-1]
        H.loc[m, IT] = [gm.get(g, np.nan)+shift for g in H.loc[m, "gen"]]
    else:                                   # resampled:真值 0 且带抽样噪声
        m = H.year != ys[0]
        pool = {g: H[(H.year == ys[0]) & (H.gen == g)][IT].to_numpy(float) for g in H.gen.unique()}
        H.loc[m, IT] = [RNG.choice(pool[g]) if len(pool.get(g, [])) else np.nan for g in H.loc[m, "gen"]]
    return H.dropna(subset=[IT])
def ctl(mode, ys, rep=NREP, shift=-0.30):
    v = [within(syn(mode, ys, shift), ys) for _ in range(rep)]
    v = np.array([x for x in v if x is not None], float)
    return float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
pc_m, pc_lo, pc_hi = ctl("planted", LATE)
nc_m, nc_lo, nc_hi = ctl("resampled", LATE)
nc_half = (nc_hi-nc_lo)/2
print(f"  正控(后半段植入 −0.300)⇒ 取回中位 **{pc_m:+.4f}** [{pc_lo:+.4f}, {pc_hi:+.4f}]")
print(f"  负控(后半段各年从首年同世代重抽 ⇒ 真值 0)⇒ 中位 **{nc_m:+.4f}** [{nc_lo:+.4f}, {nc_hi:+.4f}],"
      f"噪声半宽 **{nc_half:.4f}** · 容差 `NC_TOL = {NC_TOL}` **事先写死** ⇒ 比值 **{NC_TOL/nc_half:.2f}×**")

G = Gate("#823 · 放弃十年分辨率,换回功效")
G.asserted("⓪ **前提(本轮第一个可失败的主张,而它是关于仪器不是关于人 —— `#822` 的教训)**:"
           f"两个半段的年份聚类噪声半宽都必须 < {MATTERS},否则本轮没有资格往下读",
           bool(floor_ok),
           f"前半 {RES['前半']['half']:.4f} · 后半 {RES['后半']['half']:.4f}(阈 {MATTERS})", kind="control")
G.asserted("① 正控:后半段植入 −0.300 必须取回(误差 < 0.10)",
           bool(abs(pc_m+0.30) < 0.10), f"取回 {pc_m:+.4f} [{pc_lo:+.4f}, {pc_hi:+.4f}]", kind="control")
G.identity_control("② 负控:后半段各年从首年同世代重抽(真值 0 且带噪声),该项必须 ≈ 0"
                   "(⚠ **参照真的是 0**)—— 容差 0.05 **事先写死**(`#817`②)",
                   observed=nc_m, expected=0.0, tol=NC_TOL, noise_half_width=nc_half,
                   what=f"{NREP} 次重复,95% 跨度 [{nc_lo:+.4f}, {nc_hi:+.4f}]")
G.asserted("③ 前提(跑前写下的混淆,重申因为没被解决):两段各自的 `REL` 分布已印出,"
           "**本设计修不了「层的主语在漂」,只能量它**", True,
           f"前半均值 {D[D.year.isin(EARLY)].REL.mean():+.3f} · 后半 {D[D.year.isin(LATE)].REL.mean():+.3f}",
           kind="control")
G.asserted("④ 前提:切点 1990 **跑前写死**,且来自 `#812` 的独立结果 —— "
           "⚠ **但它仍来自同一份数据,不是外生的,如实登记**", True,
           f"前半 {len(EARLY)} 年 · 后半 {len(LATE)} 年", kind="control")
late = RES["后半"]
v_late = Gate.interval_verdict(late["lo"], late["hi"], 0.0, MATTERS)
G.asserted("⑤ kill(预注册):「虔诚者在后半段往回走了」(世界 B)要成立,需后半段区间**排除 0 且为负**",
           bool(v_late == "EXCLUDES" and late["hi"] < 0),
           f"后半 {late['point']:+.4f} [{late['lo']:+.4f}, {late['hi']:+.4f}] ⇒ {v_late}", kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
early = RES["前半"]
v_early = Gate.interval_verdict(early["lo"], early["hi"], 0.0, MATTERS)
compat = [x for x in (-0.10, 0.0, 0.10, 0.166) if late["lo"] <= x <= late["hi"]]
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif v_late == "EXCLUDES" and late["hi"] < 0:
    V = (f"**B 往回走了。** 后半段 **{late['point']:+.4f}** [{late['lo']:+.4f}, {late['hi']:+.4f}]\n"
         f"  ⇒ **这条线索里第一句「有人朝反方向走」的话:后半段虔诚的美国人不是跟得慢,是自己往回退了。**")
elif v_late == "TIGHT_NULL":
    V = (f"**A 停住了。** 前半 **{early['point']:+.4f}** [{early['lo']:+.4f}, {early['hi']:+.4f}] ⇒ {v_early} · "
         f"后半 **{late['point']:+.4f}** [{late['lo']:+.4f}, {late['hi']:+.4f}] ⇒ **窄零**\n"
         f"  ⇒ **一句关于人的话:虔诚的美国人在前半段跟着社会一起变宽容;到后半段他们停下了 ——\n"
         f"  不是往回走,是不再动。而其余人没有停。**\n"
         f"  ⚠ **「停住」与「往回走」是两句话,而本轮的区间只支持前一句。**")
else:
    V = (f"**C 或分辨不出。** 后半 **{late['point']:+.4f}** [{late['lo']:+.4f}, {late['hi']:+.4f}] ⇒ **{v_late}**"
         f",相容于 {compat}(`#812`③)· 前半 {early['point']:+.4f} [{early['lo']:+.4f}, {early['hi']:+.4f}] ⇒ {v_early}\n"
         f"  ⇒ **即使把噪声底买下来了,这个问题在半段分辨率上仍然答不了 —— 如实说。**")
print(V)
json.dump(dict(item=IT, cut=CUT, matters=MATTERS, nc_tol=NC_TOL, B=B, halves=RES,
               verdict_early=v_early, verdict_late=v_late, floor_ok=bool(floor_ok),
               pos_control=dict(median=pc_m, lo=pc_lo, hi=pc_hi, planted=-0.30),
               neg_control=dict(median=nc_m, lo=nc_lo, hi=nc_hi, half_width=nc_half, reference=0.0),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"halves.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'halves.json'}")
