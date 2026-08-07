"""#849 · E03·A85·R288 —— `#846` 顺手捞到的那条线索,经不经得住换一半人

`#846` 扫了 **8 题 × 2 层 = 16 格**,顺手捞到一格:
**虔诚层在 `sexeduc`(该不该在学校教性教育)上,九十年代偏离 +0.0720 [+0.015, +0.126]。**
`#846` 已如实登记它「没被预注册,所以现在只是一个假设」(`#111`)。

**⚠⚠ 而先要把它现在的地位说准,不能含糊:它的 `p = 0.0155`,
而 BH 在 16 格族上的秩-2 阈是 `0.05×2/16 = 0.00625` ⇒ 它没过 BH,也没过 BY。**
⇒ **所以它现在不是「一个较弱的发现」,是「一个被多重性校正拒绝了的格」。**
**唯一能让它翻身的,不是再算一遍同一批人,是拿它没被发现的那部分数据来看。**

G1 估计量:**同一格(`sexeduc` × 虔诚层 × 1990s)的偏离,在两种把样本切成两半的方式上各是多少** ——
   ① 按**受访者性别**切;② 按**九十年代内的奇数/偶数调查年**切。
   **两种切法都与「这一格是怎么被发现的」无关**(发现它的是跨 16 格的扫描),
   所以每一半都是它没被发现在其上的数据。

三个世界:
   A **两种切法的两半都同号且都超各自地板** ⇒ **它挺过了留出复现**,
     从「被 BH 拒绝的格」升为「值得预注册一轮去打的假设」。
   B **同号但只有部分超地板** ⇒ **不足以翻身**,仍是假设,但记下它没被证伪。
   C **有一半反号** ⇒ **它是扫描出来的噪声** ⇒ 明确写死,不再回来。

预注册判词(条件式,⓪ 排最前):
  ⓪ **功效闸**:在半样本上植入一个等于观测量的偏离,若捞不回 ≥2 个自助 SD
     ⇒ **UNVERIFIED,不看结果**(`#835`/`#845` 的做法:拦在读之前)。
  if 功效闸过 and 正控开火 and 负控开火:
      4 个半样本全同号且全超地板 -> A
      全同号但未全超地板         -> B
      任一反号                   -> C
  else: UNVERIFIED

⚠ 跑之前写下的最强混淆:**`sexeduc` 是二值题(1/2),而虔诚层在九十年代均值 1.16–1.22,
  贴着下端** —— 余量占用 `#846` 已量为 6.6%,不触发天花板标记;
  **但半样本会把噪声放大,而地板效应在小样本上更容易造出假的「同号」。**
  ⇒ 控制:**每一半的地板由它自己的置换零给出**(打乱十年内的年份标签),
  **不用全样本的地板去判半样本。**
⚠ 本轮换不了仪器(GSS);而它不需要 —— 本轮换的是**看哪一半人**,这正是留出要的那根轴。
"""
import json, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT, B, DEC = "sexeduc", 4000, 1990

d = pd.read_stata(gp, columns=["year", "sex", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= 2))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3)), ("sex", (1, 2))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
R = M.dropna(subset=["attend", "reliten", "fund"]).copy()
R["REL"] = zs(R[["attend", "reliten", "fund"]]).mean(axis=1); M = M.join(R["REL"])
t = M.groupby("year")["REL"].transform(
    lambda s: pd.qcut(s, 3, labels=False, duplicates="drop") if s.notna().sum() > 30 else np.nan)
M["dev"] = (t == 2)
W = M[M[IT].notna() & M.dev & M.sex.notna()].copy()

years = sorted(y for y, g in W.groupby("year") if len(g) >= 120)
dec = {}
for y in years: dec.setdefault((y//10)*10, []).append(y)
dec = {k: v for k, v in dec.items() if len(v) >= 3}
S = years; span = S[-1]-S[0]
print(f"=== ⓪ 硬规则①:`{IT}` 1–2 · 虔诚层每年 ≥120 的合格年 **{len(S)}** 个 {S[0]}–{S[-1]} · "
      f"n={len(W):,} · 可用十年 {sorted(dec)}")
# ⚠ **不手抄 `#846` 的数,从它的产物读**(`#840`/`#841` 立的规矩,而我的闸刚拦下了第一版)。
P846 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/"
                      "A84_是世俗那边走了还是只在这一题上走了/R285_八题上各自是谁在动/"
                      "results/whose_move_on_every_item.json", encoding="utf-8"))
_p = P846["grid"]["sexeduc|虔诚"]["p"]; _fam = P846["family"]; _q = P846["q"]
_thr = _q*2/_fam
print(f"  ⚠ **先把它现在的地位说准**(数从 `#846` 的产物读,不手抄):这一格 `p = {_p:.4f}`,"
      f"而 BH 在 {_fam} 格族上的秩-2 阈是 `{_q}×2/{_fam} = {_thr:.5f}` ⇒ "
      f"**{'它没过 BH,也没过 BY' if _p > _thr else '它过了 BH'}** —— "
      f"**不是「较弱的发现」,是「被多重性校正拒绝了的格」。**")

def dep(sub, Bv=B, seed=849):
    ys = {}
    for y, g in sub.groupby("year"):
        v = g[IT].to_numpy(float)
        if len(v) >= 60: ys[int(y)] = v
    have = [y for y in S if y in ys]
    if len(have) < 8 or DEC not in {(y//10)*10 for y in have}: return None
    yy = [y for y in have if (y//10)*10 == DEC]
    if len(yy) < 3: return None
    sp = have[-1]-have[0]
    ref = (ys[have[-1]].mean()-ys[have[0]].mean())*(yy[-1]-yy[0])/sp
    obs = float(ys[yy[-1]].mean()-ys[yy[0]].mean()) - ref
    rg = np.random.default_rng(seed); out = np.empty(Bv)
    r = lambda a: a[rg.integers(0, len(a), len(a))]
    for i in range(Bv):
        rf = (r(ys[have[-1]]).mean()-r(ys[have[0]]).mean())*(yy[-1]-yy[0])/sp
        out[i] = r(ys[yy[-1]]).mean() - r(ys[yy[0]]).mean() - rf
    lo, hi = np.quantile(out, [.025, .975])
    # 该半样本**自己的**地板:打乱十年内年份标签的置换零
    nul = np.empty(1500)
    pool = np.concatenate([ys[y] for y in yy]); sizes = [len(ys[y]) for y in yy]
    for i in range(1500):
        pm = rg.permutation(pool); k = 0; z = {}
        for y, n in zip(yy, sizes): z[y] = pm[k:k+n]; k += n
        rf = (r(ys[have[-1]]).mean()-r(ys[have[0]]).mean())*(yy[-1]-yy[0])/sp
        nul[i] = z[yy[-1]].mean() - z[yy[0]].mean() - rf
    floor = float(np.quantile(np.abs(nul), 0.95))
    return dict(obs=obs, lo=float(lo), hi=float(hi), sd=float(np.std(out)), n=int(len(sub)),
                floor=floor, above=bool(abs(obs) > floor), excl=bool(lo > 0 or hi < 0))

full = dep(W)
print(f"\n=== ① 全样本复算(与 `#846` 对齐)===")
print(f"  全样本:**{full['obs']:+.4f}** [{full['lo']:+.4f},{full['hi']:+.4f}] · "
      f"自身置换地板 **{full['floor']:.4f}** ⇒ {'超地板' if full['above'] else '**未超地板**'}")

SPLITS = {
 "性别=男": W[W.sex == 1], "性别=女": W[W.sex == 2],
 "九十年代奇数年": W[~((W.year//10*10 == DEC) & (W.year % 2 == 0))],
 "九十年代偶数年": W[~((W.year//10*10 == DEC) & (W.year % 2 == 1))],
}
print(f"\n=== ② 留出复现:两种与「它怎么被发现」无关的切法(B={B},每半用**它自己的**置换地板)===")
H = {}
for nm, sub in SPLITS.items():
    r = dep(sub, seed=849+len(nm))
    H[nm] = r
    if r is None: print(f"  {nm:16s} **不可用**(合格年或十年内年数不足)"); continue
    print(f"  {nm:16s} n={r['n']:>6,} **{r['obs']:+.4f}** [{r['lo']:+.4f},{r['hi']:+.4f}] · "
          f"自身地板 **{r['floor']:.4f}** ⇒ "
          f"{'**超地板**' if r['above'] else '未超地板'} · "
          f"{'同号' if np.sign(r['obs']) == np.sign(full['obs']) else '**反号**'}")
ok = [nm for nm, r in H.items() if r]
same = [nm for nm in ok if np.sign(H[nm]["obs"]) == np.sign(full["obs"])]
above = [nm for nm in ok if H[nm]["above"]]
print(f"  ⇒ 可用半样本 **{len(ok)}/4** · 同号 **{len(same)}/{len(ok)}** · 超自身地板 **{len(above)}/{len(ok)}**")

print(f"\n=== ③ 控制(⓪ 功效闸排在所有关于人的判断之前)===")
sub0 = SPLITS["性别=男"]
yy = [y for y in S if (y//10)*10 == DEC]
p0 = sub0.copy()
p0.loc[p0.year == yy[-1], IT] = p0.loc[p0.year == yy[-1], IT] + abs(full["obs"])
pw = dep(p0, Bv=800, seed=11)
got = pw["obs"] - H["性别=男"]["obs"]
print(f"  ⓪ **功效闸**:在半样本(性别=男)的九十年代尾年植入 +{abs(full['obs']):.4f} ⇒ "
      f"偏离动 **{got:+.4f}**(预期 +{abs(full['obs']):.4f})· "
      f"相对该半样本 SD **{abs(got)/H['性别=男']['sd']:.2f} 个 SD** ⇒ "
      f"{'**有功效**' if abs(got)/H['性别=男']['sd'] > 2 else '**没功效 ⇒ 不看结果**'}")
n0 = dep(sub0.assign(**{IT: np.random.default_rng(3).permutation(sub0[IT].values)}), Bv=800, seed=12)
print(f"  负控:把 `{IT}` 的作答**在半样本内整体打乱**(时间结构被毁,人不变)⇒ "
      f"**{n0['obs']:+.4f}**,自身地板 {n0['floor']:.4f} ⇒ "
      f"{'**落在地板内**' if abs(n0['obs']) <= n0['floor'] else '⚠ **超出地板**'}")
print(f"     ⚠ **「这个零该不该是零?」该** —— 打乱作答后,十年偏离的期望就是 0。")

G = Gate("#849 · `sexeduc` 那条线索经不经得住换一半人")
G.asserted("⓪ **功效闸(排在所有关于人的判断之前)**:半样本上植入一个等于观测量的偏离,"
           "必须捞回 ≥2 个该半样本自己的自助 SD",
           bool(abs(got)/H["性别=男"]["sd"] > 2),
           f"动 {got:+.4f} = {abs(got)/H['性别=男']['sd']:.2f} 个 SD", kind="control")
G.asserted("① 负控:半样本内把作答整体打乱(时间结构毁掉、人不变)⇒ 偏离必须落在它自己的置换地板内"
           "(⚠ **这个零该是零**:打乱作答后十年偏离期望为 0)",
           bool(abs(n0["obs"]) <= n0["floor"]),
           f"{n0['obs']:+.4f} vs 地板 {n0['floor']:.4f}", kind="control")
G.asserted("② 前提(跑前写下的最强混淆):`sexeduc` 二值且虔诚层贴下端,**半样本会放大噪声,"
           "地板效应在小样本上更容易造出假的「同号」** ⇒ **每一半的地板由它自己的置换零给出**,"
           "不用全样本的地板判半样本",
           bool(all(H[nm] and np.isfinite(H[nm]["floor"]) for nm in ok)),
           " · ".join(f"{nm}:{H[nm]['floor']:.4f}" for nm in ok), kind="control")
G.asserted("③ kill(预注册):「它挺过了留出复现」要成立,需**四个半样本全部同号且全部超自身地板**",
           bool(len(ok) == 4 and len(same) == 4 and len(above) == 4),
           f"可用 {len(ok)}/4 · 同号 {len(same)} · 超地板 {len(above)}", kind="kill",
           yardstick="每个半样本自己的十年偏离,对照它自己的置换地板(95 分位)",
           yardstick_noise=float(np.mean([H[nm]["sd"] for nm in ok])) if ok else 0.0)
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*100)
if not adm:
    VERD = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif len(same) == len(ok) == 4 and len(above) == 4:
    VERD = (f"**A 它挺过了留出复现。** 四个半样本全部同号、全部超自身置换地板 ⇒ "
            f"从「被 BH 拒绝的格」升为**值得预注册一轮去打的假设**。\n"
            f"  ⇒ **一句关于人的话:虔诚的那一层在九十年代唯一动过的地方,\n"
            f"  不是关于别人怎么活,是关于学校该不该教孩子性这件事。**")
elif len(same) < len(ok):
    VERD = (f"**C 它是扫描出来的噪声 ⇒ 写死,不再回来。** {len(ok)-len(same)} 个半样本**反号**"
            f"({[nm for nm in ok if nm not in same]})。\n"
            f"  ⇒ **一句关于方法的话:一格从 16 格里被捞出来、又没过 BH,\n"
            f"  换一半人就换了方向 —— 这正是多重性校正当初拒绝它的理由,而现在它被独立证实了。**")
else:
    VERD = (f"**B 同号但未全超地板 ⇒ 不足以翻身,仍是假设。** 同号 {len(same)}/{len(ok)},"
            f"超自身地板 **{len(above)}/{len(ok)}**({above})。\n"
            f"  ⇒ **登记为「没被证伪,也没被支持」** —— 而 `#846` 给它的地位不变:"
            f"**一个被多重性校正拒绝的格,不是一个较弱的发现。**")
print(VERD)
json.dump(dict(item=IT, decade=DEC, full=full, halves=H, usable=ok, same_sign=same, above_floor=above,
               power=dict(plant=abs(full["obs"]), got=got,
                          in_sd=abs(got)/H["性别=男"]["sd"]),
               neg_control=n0, bh_status=dict(p=_p, family=_fam, q=_q, rank2_threshold=_thr,
                                              rejected=bool(_p > _thr), source="read from #846 artifact"),
               admissible=adm, verdict=VERD, gate_ok=G.verdict()),
          open(OUT/"held_out_sexeduc.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'held_out_sexeduc.json'}")
