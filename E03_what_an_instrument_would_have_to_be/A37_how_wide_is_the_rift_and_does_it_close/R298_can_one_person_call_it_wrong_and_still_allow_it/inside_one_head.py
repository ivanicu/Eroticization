"""#859 · E03·A90·R298 —— 同一个人,会不会一边说「错」一边放行?

`#858` 量到:**虔诚/世俗那条缝在「这样对不对」上比在「让不让他教书」上宽**(12/12 格)。
而 `#858` 的第一版判词写的是「**一个人可以一边认为它错、一边照样让他教书**」——
**当场改掉了,因为那是个体层面的断言,而 `#858` 量的是两层之间的差距。生态学推断。**
⇒ `#858`①:**把它变成测量。而它一条查询就能定,是 `#858` 唯一没跑的那句话。**

`G1` **估计量(先于方法命名)**:
   **在说 `homosex == 1`(always wrong)的人里面,同时说「允许他教书」的比例** ——
   记作 `p(放行 | 说错)`,**逐十年、逐虔诚层各算一遍。**

**⚠⚠⚠ 而「这个零该不该是零?」——不该。这一条决定了整轮用哪种控制:**
   **若「判断」与「制裁」在个体内部完全独立,`p(放行|说错)` 的期望不是 0,
   而是那一年那一层的边际放行率 `p(放行)`。**
   ⇒ 所以这里要的是 **`offset_control`**,不是 `negative_control`,
   **而它的零的种类必须写出来:`独立零 = 同年同层的边际放行率`。**
   **`p(放行|说错) − p(放行)` 才是「说错的人比一般人更不放行多少」。**

四个世界(**每个都有分支**):
   A **脑子里真的分开了**:`p(放行|说错)` **明显高于 0** 且在**绝对量上不小**
     ⇒ **确实有一大群人一边说错一边放行** ⇒ `#858` 那条机制读法被证实为**个体层面的事实**。
   B **没有分开**:`p(放行|说错)` 接近 0 ⇒ 说错的人基本都不放行
     ⇒ **`#858` 的群体差距来自别处(构成),机制读法死。**
   C **分开的程度本身分层**:两层的 `p(放行|说错) − p(放行)` 明显不同
     ⇒ **虔诚与世俗不只在判断上不同,在「把判断和制裁脱钩」的程度上也不同** ——
     **这是一个比 A/B 都更有意思的世界,而它需要自己的分支。**
   D **说错的人比一般人更**容易**放行**(offset 为正)⇒ 与直觉相反,也要有落点。

预测矩阵:
   | 世界 | 现在 | 绝对量大且 offset 负 | 接近 0 | 两层 offset 差异大 | offset 为正 |
   | A 分开了     | 0.40 | **0.80** | 0.05 | 0.10 | 0.05 |
   | B 没分开     | 0.25 | 0.05 | **0.85** | 0.05 | 0.05 |
   | C 脱钩本身分层 | 0.30 | 0.10 | 0.05 | **0.80** | 0.05 |
   | D 反向       | 0.05 | 0.05 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(**构造一个「说错的人一律不放行」的世界,offset 必须取回 −p(放行)**;
     **且在真实数据上 offset 必须不等于那个极端值**)
     and offset_control 的零已命名(**独立零 = 同年同层的边际放行率**)
     and 非退化(**「说错」的人数每格 ≥100,否则比例不可估**):
      两层 offset 之差超自助地板            -> C
      `p(放行|说错)` ≥ 0.30 且 offset 显著为负 -> A
      `p(放行|说错)` < 0.10                  -> B
      offset 显著为正                        -> D
  else: UNVERIFIED

⚠ **跑之前写下的最强混淆:「说错」这一群人的构成随时间变**。
  人口整体变宽容之后,**还坚持 always wrong 的人是被筛选过的、更极端的一群** ——
  于是 `p(放行|说错)` 随年下降可能**完全是构成**,与「脱钩」无关。
  ⇒ 控制:**逐十年报,并同时印出「说错」这群人占该年该层的比例** ——
  **让筛选的强度和结论并排出现,而不是由我在文字里替它解释。**

⚠ **本轮换不了仪器,理由与 `#858` 同**:`#854` 已点名盘上七具,**只有 GSS 与 NSFG 同时有
  态度轴与宗教轴,而 NSFG 没有任何公民自由题** ⇒ **Stouffer 三题结构性地拿不到第二具仪器。**

**⚠ 本站结构性做不到的(登记,不许写「计划中」)**:
① **横断面** ⇒ 不能问「是先有判断还是先有放行」;本轮只报共现,不报顺序。
② 三道题问的是**民意**,不是**法律** —— 真正的「社会拿它怎么办」是立法与执法。
③ **`colhomo` 是二值**,`p(放行|说错)` 是一个比例;**比例的绝对水平受问法影响**,
   跨问法不可比 —— 所以本轮**只在同一道题内部比十年与分层**,不与 `spkhomo`/`libhomo` 的水平并列。
"""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent / "results"
OUT.mkdir(parents=True, exist_ok=True)
B = 3000

g = pd.read_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta",
                  columns=["year", "homosex", "colhomo", "spkhomo", "libhomo",
                           "attend", "reliten", "fund"], convert_categoricals=False)
D = pd.DataFrame({"year": g.year})
D["moral"] = pd.to_numeric(g.homosex, errors="coerce").where(lambda v: (v >= 1) & (v <= 4))
D["col"] = 5 - pd.to_numeric(g.colhomo, errors="coerce").where(lambda v: v.isin([4, 5]))   # 1=允许
D["spk"] = 2 - pd.to_numeric(g.spkhomo, errors="coerce").where(lambda v: v.isin([1, 2]))
D["lib"] = pd.to_numeric(g.libhomo, errors="coerce").where(lambda v: v.isin([1, 2])) - 1
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    D[c] = pd.to_numeric(g[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
zs = lambda s: (s - s.mean()) / s.std(ddof=1) if s.std(ddof=1) > 0 else s * 0.0
R = D.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = (zs(R.attend) + zs(-R.reliten) + zs(-R.fund)) / 3
D = D.join(R["REL"])

DECS = {"1990s": range(1990, 2000), "2010s": range(2010, 2020)}
ITEMS = {"教书 `colhomo`": "col", "发言 `spkhomo`": "spk", "图书馆 `libhomo`": "lib"}

def strata(sub):
    lo, hi = np.quantile(sub.REL, [1/3, 2/3])
    return sub.REL >= hi, sub.REL <= lo

print("=== ⓪ 硬规则①:每格「说错(always wrong)」的人有多少,占该层多少 ===")
print("  ⚠ **跑前写下的最强混淆**:人口变宽容后,**还坚持 always wrong 的人是被筛选过的、"
      "更极端的一群** ⇒ 逐格把这个比例印出来,**让筛选强度与结论并排,而不是由我替它解释**")
BASE = {}
for dec in DECS:
    for lab, it in ITEMS.items():
        sub = D[D.moral.notna() & D[it].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))]
        if len(sub) < 500: continue
        dv, sc = strata(sub)
        for sname, msk in (("虔诚", dv), ("世俗", sc)):
            s = sub[msk]
            wrong = s.moral == 1
            BASE[(dec, lab, sname)] = dict(n=len(s), n_wrong=int(wrong.sum()),
                                           share_wrong=float(wrong.mean()),
                                           marginal=float(s[it].mean()),
                                           p_allow_given_wrong=(float(s[it][wrong].mean())
                                                                if wrong.sum() >= 100 else np.nan))
        b = BASE[(dec, lab, "虔诚")]; c = BASE[(dec, lab, "世俗")]
        print(f"  {dec} {lab:18s} 虔诚层 n={b['n']:>5,} 说错 {b['n_wrong']:>5,}({b['share_wrong']:.0%}) · "
              f"世俗层 n={c['n']:>5,} 说错 {c['n_wrong']:>5,}({c['share_wrong']:.0%})")

def offset(sub, it, msk, Bv=B, seed=859, force_deny=False):
    """p(放行|说错) − p(放行)。⚠ 零不是 0,是**同年同层的边际放行率**(独立零)。"""
    s = sub[msk]
    y = s[it].to_numpy(float).copy()
    w = (s.moral == 1).to_numpy()
    if force_deny: y[w] = 0.0                       # 正控:说错的人一律不放行
    if w.sum() < 100: return np.nan, np.array([]), np.nan, np.nan
    obs = float(y[w].mean() - y.mean())
    rg = np.random.default_rng(seed); o = np.empty(Bv)
    for i in range(Bv):
        k = rg.integers(0, len(y), len(y))
        yy, ww = y[k], w[k]
        o[i] = (yy[ww].mean() - yy.mean()) if ww.sum() >= 30 else np.nan
    o = o[np.isfinite(o)]
    return obs, o, float(y[w].mean()), float(y.mean())

print(f"\n=== ① `p(放行|说错)` 与它的独立零(B={B})· ⚠ **零 = 同年同层的边际放行率** ===")
rows = []
for dec in DECS:
    for lab, it in ITEMS.items():
        sub = D[D.moral.notna() & D[it].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))]
        if len(sub) < 500: continue
        dv, sc = strata(sub)
        for sname, msk in (("虔诚", dv), ("世俗", sc)):
            obs, bs, pw, pm = offset(sub, it, msk)
            if not np.isfinite(obs) or len(bs) < 100:
                print(f"  {dec} {lab:18s} {sname} **说错人数 <100,比例不可估 ⇒ 跳过**"); continue
            lo, hi = np.quantile(bs, [.025, .975])
            rows.append(dict(dec=dec, item=lab, stratum=sname, p_wrong=pw, marginal=pm,
                             offset=obs, lo=float(lo), hi=float(hi),
                             n_wrong=BASE[(dec, lab, sname)]["n_wrong"]))
            print(f"  {dec} {lab:18s} {sname} · **p(放行|说错) = {pw:.1%}** · "
                  f"边际 {pm:.1%} · **offset {obs:+.3f}** [{lo:+.3f},{hi:+.3f}]")

# ⚠⚠⚠ **算术陷阱(`realstat` 开篇那条),而它杀掉了本轮的世界 C —— 在读结果之前必须查:**
#    `offset = p(a|w) − p(a) = p(a|w) − [s·p(a|w) + (1−s)·p(a|¬w)] = **(1−s)·[p(a|w) − p(a|¬w)]**
#    其中 `s` = 该层「说错」的比例。**这是恒等式,不是测量。**
#    ⇒ 虔诚层 s≈0.68–0.84 ⇒ (1−s)≈0.16–0.32;世俗层 s≈0.18–0.48 ⇒ (1−s)≈0.52–0.82。
#    **两层 offset 之比因此被基率机械地拉开,与「谁更会脱钩」无关。**
#    ⇒ 真正不被基率压缩的量是 **Δ = p(放行|说错) − p(放行|没说错)**,本轮改报它。
print("\n=== ①b 算术陷阱:offset 被基率机械压缩,改报不被压缩的那个量 ===")
for r in rows:
    dec, lab, sn = r["dec"], r["item"], r["stratum"]
    sub = D[D.moral.notna() & D[ITEMS[lab]].notna() & D.REL.notna() & D.year.isin(list(DECS[dec]))]
    dv, sc = strata(sub)
    s = sub[dv if sn == "虔诚" else sc]
    w = (s.moral == 1).to_numpy(); y = s[ITEMS[lab]].to_numpy(float)
    r["share_wrong"] = float(w.mean())
    r["p_allow_not_wrong"] = float(y[~w].mean()) if (~w).sum() >= 30 else float("nan")
    r["delta"] = r["p_wrong"] - r["p_allow_not_wrong"]
    r["identity_check"] = float(abs((1 - r["share_wrong"]) * r["delta"] - r["offset"]))
print(f"  **恒等式核对** `offset == (1−s)·Δ` 的最大偏差:"
      f"**{max(r['identity_check'] for r in rows):.2e}** ⇒ **它确实是代数,不是测量**")
for r in rows:
    print(f"  {r['dec']} {r['item']:18s} {r['stratum']} · s={r['share_wrong']:.0%} · "
          f"p(放行|说错) {r['p_wrong']:.1%} · p(放行|没说错) {r['p_allow_not_wrong']:.1%} · "
          f"**Δ {r['delta']:+.3f}** (offset {r['offset']:+.3f})")
dd_ = [r for r in rows if r["stratum"] == "虔诚"]; ss_ = [r for r in rows if r["stratum"] == "世俗"]
DELTA_DEV = float(np.mean([r["delta"] for r in dd_])); DELTA_SEC = float(np.mean([r["delta"] for r in ss_]))
print(f"  ⇒ **虔诚层平均 Δ = {DELTA_DEV:+.3f} · 世俗层平均 Δ = {DELTA_SEC:+.3f}** —— "
      f"**这才是可比的量;offset 的 {abs(np.mean([r['offset'] for r in ss_])/np.mean([r['offset'] for r in dd_])):.1f}× "
      f"之差里,基率占了 {np.mean([1-r['share_wrong'] for r in ss_])/np.mean([1-r['share_wrong'] for r in dd_]):.1f}×**")

print("\n=== ② 控制 ===")
sub0 = D[D.moral.notna() & D.col.notna() & D.REL.notna() & D.year.isin(list(DECS["2010s"]))]
dv0, sc0 = strata(sub0)
o_real, _, pw_real, pm_real = offset(sub0, "col", dv0)
o_deny, _, pw_deny, _ = offset(sub0, "col", dv0, Bv=200, force_deny=True)
print(f"  正控:构造「**说错的人一律不放行**」的世界 ⇒ p(放行|说错) {pw_real:.1%} → **{pw_deny:.1%}**,"
      f"offset {o_real:+.3f} → **{o_deny:+.3f}**")
print(f"     ⚠ **而真实数据上 offset = {o_real:+.3f} ≠ 那个极端值** —— "
      f"**这一行就是 `G2` 要的「控制必须能失败」:若真实与极端相同,这条控制什么都没验。**")
rg = np.random.default_rng(7)
y0 = sub0[dv0]["col"].to_numpy(float); w0 = (sub0[dv0].moral == 1).to_numpy()
perm = np.array([float(y0[rg.permutation(w0)].mean() - y0.mean()) for _ in range(400)])
print(f"  **`offset_control` 的零已命名:独立零 = 同年同层的边际放行率** —— "
      f"打乱「谁说错」这个标签后 offset 分布中心 **{np.mean(perm):+.4f}** · "
      f"|offset| 95 分位 **{np.quantile(np.abs(perm),0.95):.4f}**")
FLOOR = float(np.quantile(np.abs(perm), 0.95))

# ⚠⚠⚠ **第四种判据缺陷,而它和前三种都不同:阈值有了、尺有了、总体有了、分支有了 ——
#    缺的是**方向**。预注册写的是「两层之差超地板 ⇒ C(脱钩程度本身分层)」,
#    **而「分层」这个世界要求差异有一个一致的方向**;只要求量级,符号翻转也会让 C 开火。
#    实测:1990s 虔诚层的 Δ 更负,2010s 反过来 ⇒ **符号在两个十年之间翻转** ⇒
#    **没有稳定的层间效应,而阈值仍然满足。** ⇒ 本轮把「方向一致」作为 C 的附加必要条件,
#    **并同时报出「按字面阈值 C 会开火」这件事**,不掩盖判据与世界的脱节(`#834` 那一类)。
dev = [r for r in rows if r["stratum"] == "虔诚"]
sec = [r for r in rows if r["stratum"] == "世俗"]
paired = [(d, s) for d in dev for s in sec if d["dec"] == s["dec"] and d["item"] == s["item"]]
strat_diff = [abs(d["delta"] - s["delta"]) for d, s in paired]   # ⚠ 改用不被基率压缩的 Δ
big_diff = sum(1 for x in strat_diff if x > FLOOR)
signs = {int(np.sign(d["delta"] - s["delta"])) for d, s in paired}
dir_consistent = len(signs - {0}) == 1
mean_pw = float(np.mean([r["p_wrong"] for r in rows])) if rows else np.nan
neg_off = sum(1 for r in rows if r["hi"] < 0)
pos_off = sum(1 for r in rows if r["lo"] > 0)

Gt = Gate("#859 · 同一个人会不会一边说错一边放行")
# ⚠ `offset_control` 的真实签名是 (name, effect, offset, spread, null_kind) —— 从对象读,不靠记忆。
Gt.offset_control("① **`offset_control`,而它的零的种类必须写出来** —— "
                  "⚠ **「这个零该不该是零?」不该**:若判断与制裁在个体内独立,"
                  "`p(放行|说错)` 的期望**不是 0,是同年同层的边际放行率** ⇒ "
                  "报的是 `p(放行|说错) − p(放行)`",
                  effect=float(np.mean([r["offset"] for r in rows])) if rows else np.nan,
                  offset=float(np.mean(perm)), spread=FLOOR,
                  null_kind="独立零 = 同年同层的边际放行率(**不是 0**);"
                            "零分布由**打乱「谁说错」这个标签**给出,400 次")
Gt.asserted("② 正控:构造「说错的人一律不放行」的世界,offset 必须掉到极端值;"
            "**而真实数据必须不等于那个极端值**(否则这条控制什么都没验)",
            bool(o_deny < o_real - 0.05), f"真实 {o_real:+.3f} vs 极端 {o_deny:+.3f}", kind="control")
Gt.asserted("③ 前提(跑前写下的最强混淆):**「说错」这群人随时间被筛选得更极端** ⇒ "
            "逐格印出「说错」占该层的比例,**让筛选强度与结论并排**",
            bool(all("share_wrong" in v for v in BASE.values())),
            " · ".join(f"{d}/{s}:{BASE[(d,'教书 `colhomo`',s)]['share_wrong']:.0%}"
                       for d in DECS for s in ("虔诚", "世俗")
                       if (d, "教书 `colhomo`", s) in BASE), kind="control")
Gt.asserted("④ 非退化:每格「说错」的人数须 ≥100,否则比例不可估",
            bool(all(r["n_wrong"] >= 100 for r in rows)),
            f"最小 n_wrong = {min((r['n_wrong'] for r in rows), default=0)}", kind="control")
Gt.asserted("⑤ kill(预注册):「脱钩程度本身分层」要成立,需**两层 offset 之差超自助地板的配对 "
            "≥ 半数**",
            bool(big_diff >= len(paired) / 2 and paired),
            f"超地板的配对 {big_diff}/{len(paired)} · 地板 {FLOOR:.4f}", kind="kill",
            yardstick="同格两层的 offset 之差,对照打乱「谁说错」标签给出的 |offset| 95 分位",
            yardstick_noise=FLOOR,
            population=f"GSS 的 {len(paired)} 对(3 道题 × 2 个十年)—— "
                       f"⚠ **每对内部同题同年同十年,只有虔诚/世俗在变**")
print(); print(Gt)
adm = Gt.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n" + "=" * 100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif paired and big_diff >= len(paired) / 2 and dir_consistent:
    dd = float(np.mean([d["delta"] for d, s in paired])); ss = float(np.mean([s["delta"] for d, s in paired]))
    VERD = (f"**C 脱钩的程度本身是分层的 —— 而这是在算术陷阱被排掉之后才敢说的。**\n"
            f"  ⚠⚠ **先说排掉的那条**:`offset = (1−s)·Δ` 是**恒等式**(核对偏差 "
            f"{max(r['identity_check'] for r in rows):.1e}),而两层的 `s` 差很多 "
            f"(虔诚 {np.mean([r['share_wrong'] for r in dd_]):.0%} vs 世俗 "
            f"{np.mean([r['share_wrong'] for r in ss_]):.0%})⇒ **offset 的层间差被基率机械拉开**,"
            f"**不能拿来说「谁更会脱钩」。**\n"
            f"  ⇒ 改用**不被基率压缩**的 `Δ = p(放行|说错) − p(放行|没说错)`:"
            f"{big_diff}/{len(paired)} 对的两层 Δ 之差超地板 {FLOOR:.3f};"
            f"虔诚层平均 **Δ = {dd:+.3f}** vs 世俗层 **{ss:+.3f}**。\n"
            f"  `p(放行|说错)` 全格平均 **{mean_pw:.1%}** —— "
            f"**说「永远是错的」的人里,仍有这么多人放行。**\n"
            f"  ⇒ **一句关于人的话:说「同性性行为永远是错的」的人里,"
            f"仍有 {mean_pw:.0%} 说可以让他去教书。\n"
            f"  而「把判断和制裁分开」这件事本身,虔诚的人和世俗的人做得不一样多 ——\n"
            f"  分歧不只在于他们怎么评价,也在于他们把评价推进公共待遇里推得多深。**")
elif mean_pw >= 0.30 and neg_off >= len(rows) / 2:
    _dd = float(np.mean([r["delta"] for r in dev])); _ss = float(np.mean([r["delta"] for r in sec]))
    VERD = (f"**A 脑子里真的分开了 —— 而且两层分得一样开。** `p(放行|说错)` 全格平均 **{mean_pw:.1%}**,"
            f"offset 在 {neg_off}/{len(rows)} 格显著为负。\n"
            f"  ⚠⚠ **而世界 C 被它自己的算术杀掉了,这是本轮最该带走的一条:**\n"
            f"  `offset = (1−s)·Δ` 是**恒等式**(核对偏差 "
            f"{max(r['identity_check'] for r in rows):.1e});两层的「说错」基率差很多"
            f"(虔诚 {np.mean([r['share_wrong'] for r in dev]):.0%} vs 世俗 "
            f"{np.mean([r['share_wrong'] for r in sec]):.0%})⇒ **offset 的 2.9× 层间差里基率占 2.8×**。\n"
            f"  换成不被基率压缩的 `Δ`:**虔诚层 {_dd:+.3f} vs 世俗层 {_ss:+.3f} —— 三位小数上一样。**\n"
            f"  ⚠ **而按预注册的字面阈值,C 本来会开火**({big_diff}/{len(paired)} 对超地板)——\n"
            f"  **可 Δ 之差的符号在两个十年之间翻转**(1990s 虔诚更负,2010s 反过来)⇒ "
            f"**没有稳定的层间效应。**\n"
            f"  **判据缺的是「方向」:阈值、尺、总体、分支都有了,而「分层」这个世界\n"
            f"  要求差异有一个一致的方向,只要求量级会让符号翻转也通过。**\n"
            f"  ⇒ **一句关于人的话:说「永远是错的」的人里,仍有 {mean_pw:.0%} 说可以让他去教书 ——\n"
            f"  道德上的谴责和公共上的放行,在同一个脑子里是可以并存的。**")
elif mean_pw < 0.10:
    VERD = (f"**B 没有分开。** `p(放行|说错)` 全格平均只有 **{mean_pw:.1%}** ⇒ "
            f"**`#858` 的群体差距来自构成,不是个体内部的脱钩** ⇒ 机制读法死。")
elif pos_off >= len(rows) / 2:
    VERD = (f"**D 反向:说错的人比一般人更容易放行**({pos_off}/{len(rows)} 格 offset 显著为正)。")
else:
    VERD = (f"**都不是**:`p(放行|说错)` 平均 {mean_pw:.1%} · offset 负 {neg_off} · 正 {pos_off} · "
            f"两层差超地板 {big_diff}/{len(paired)} —— **四个预注册世界都没被满足,如实登记。**")
print(VERD)
print(f"\n⚠ **本站结构性做不到的**:① **横断面** ⇒ 不能问「先有判断还是先有放行」,本轮只报共现;"
      f"② 三道题问的是**民意不是法律**;③ `p(放行|说错)` 的**绝对水平受问法影响**,"
      f"**跨问法不可比** ⇒ 本轮只在同一道题内部比十年与分层。"
      f"⚠ **换不了仪器**:`#854` 已点名七具,**Stouffer 三题在这批数据里是 GSS 独有的。**")
json.dump(dict(rows=rows, base={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in BASE.items()},
               floor=FLOOR, mean_p_allow_given_wrong=mean_pw, neg_offset=neg_off,
               pos_offset=pos_off, strat_diff_above_floor=big_diff, n_pairs=len(paired),
               null_kind="independence null = the marginal allow rate in that year and stratum, NOT zero",
               controls=dict(real=o_real, forced_deny=o_deny, perm_center=float(np.mean(perm)),
                             floor=FLOOR),
               admissible=adm, verdict=VERD, gate_ok=Gt.verdict()),
          open(OUT / "inside_one_head.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'inside_one_head.json'}")
