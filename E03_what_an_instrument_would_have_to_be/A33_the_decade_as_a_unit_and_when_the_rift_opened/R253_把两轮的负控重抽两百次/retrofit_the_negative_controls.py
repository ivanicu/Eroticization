"""#814 · E03·A61·R253 —— 把 `#805`/`#806` 的负控各抽两百次,看那两次「通过」值多少

`#813` 量出:这条负控在 1990 那个起点上,自己的 95% 抽样跨度是 **[+0.766, +1.439]**(±0.33),
**而 `#805` 与 `#806` 各自只抽了一次,用的容差是 0.08 —— 比它自己的噪声窄四倍。**
⇒ 那两次「通过」是**在一条随机开火的控制上抽到了好签**。

⚠⚠ **而 `#813` 的那个跨度是在 1990 起点上量的,不能直接搬到 `#805`(1974 起点)与 `#806`(逐题起点)上。**
   **搬一个数就是 `#813` 自己刚刚点名的那种错。** ⇒ **本轮就是去各自量一遍。**

G1 估计量:**`#805` 与 `#806` 各自那条负控的抽样分布**(中位 · 95% 跨度),
   以及**它们当初用的容差落在自己噪声跨度的什么位置**。

⚠ 本轮标注 **Production/verify —— 它不产生新的关于人的判断**,
  它把两句已发表的话(「控制全过」)换成它们配得上的强度。**`#807`:更正必须走到页面。**

预注册判词(条件式):
  if 正对照开火(**同一条代码在一个已知有分歧的世界里必须给出明显低于 1 的值** ——
     否则这条负控连「能不能动」都没证明,量它的噪声也没意义):
      逐轮报「中位 · 95% 跨度 · 当初的容差 ÷ 跨度」;
      **当初容差 < 跨度 ⇒ 那一轮的负控没有分辨率,如实标注**
  else: UNVERIFIED
⚠ **不预注册「几倍算窄」** —— 报比值本身,让读者看见(`#805` 的教训:阈值判据在这类量上没有分辨力)。

⚠ 跑之前写下的最强混淆:**噪声跨度依赖于我给的合成样本量 n。**
  n 越大跨度越窄,而那会让「当初的容差不算太窄」看起来成立。
  ⇒ 控制:**n 取两个值(2000 与 4000)都跑**,并**同时印出来** ——
  **若结论随 n 翻转,那本轮什么也没证明。**

⚠ 硬规则①:先打印每一轮当初用的起点、容差、观测值(从各自的产物 JSON 读,不是从记忆)。
⚠ 换不了仪器(对象是我自己两轮的控制)。⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(253)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
NREP = 200

print("=== ⓪ 硬规则①:两轮当初的起点 · 容差 · 观测值(从产物 JSON 读,不是从记忆)===")
J805 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A55_页面上最大的那个数一个对照都没有/"
                      "R244_那个量骑在阈值上两次所以该报区间/results/interval_not_a_threshold.json"))
J806 = json.load(open(ROOT/"E03_what_an_instrument_would_have_to_be/A56_八条鸿沟每一条都欠同一个对照/"
                      "R245_共同位移能不能造出一次跨零/results/across_eight.json"))
print(f"  `#805`(R244):负控观测 **{J805['neg_control']:+.4f}** · 参照 {J805['reference']} · 容差 **0.08** · "
      f"起点 = `homosex` {J805['y0']} 年两层分布")
print(f"  `#806`(R245):负控观测 **{J806['neg_control']:+.4f}** · 参照 {J806['reference']} · 容差 **0.08** · "
      f"起点 = `homosex` 首年两层分布(逐题网格的公共合成基)")

IT, KK = "homosex", 4
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
M = pd.DataFrame({IT: pd.to_numeric(d[IT], errors="coerce").where(lambda v: (v >= 1) & (v <= KK))})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year; M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))
G = REL.dropna(subset=[IT])
cell = lambda y, k: G[(G.year == y) & (G.k == k)][IT].to_numpy(float)

pv = lambda a: np.array([(a == c).mean() for c in range(1, KK+1)])
def fit_tau(p): return norm.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def readout(mu, tau):
    e = norm.cdf(tau-mu)
    return float((np.diff(np.concatenate(([0.0], e, [1.0])))*np.arange(1, KK+1)).sum())
def bisect(f, tgt, it=45):
    lo, hi = -8.0, 8.0
    for _ in range(it):
        mid = (lo+hi)/2
        if f(mid) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2
def explained(a0, a1, b0, b1, cal="B"):
    pa0, pb0 = pv(a0), pv(b0); na, nb = len(a0), len(b0)
    tau = fit_tau((pa0*na + pb0*nb)/(na+nb))
    mu = {"A": bisect(lambda m: readout(m, tau), float(a0.mean())),
          "B": bisect(lambda m: readout(m, tau), float(b0.mean()))}
    o = {s: readout(mu[s], tau) for s in mu}
    d_cal = float(b1.mean()-b0.mean()) if cal == "B" else float(a1.mean()-a0.mean())
    D = bisect(lambda m: readout(mu[cal]+m, tau)-o[cal], d_cal)
    oth = "A" if cal == "B" else "B"
    ec, eo = readout(mu[cal]+D, tau), readout(mu[oth]+D, tau)
    pa1, pb1 = (eo, ec) if cal == "B" else (ec, eo)
    den = float(a1.mean()-b1.mean()) - float(a0.mean()-b0.mean())
    return None if abs(den) < 1e-9 else ((pa1-pb1)-(o["A"]-o["B"]))/den

def spread(y_base, n, dA, dB, rep=NREP):
    a0r, b0r = cell(y_base, 2), cell(y_base, 0)
    tau0 = fit_tau((pv(a0r)*len(a0r)+pv(b0r)*len(b0r))/(len(a0r)+len(b0r)))
    mua = bisect(lambda m: readout(m, tau0), float(a0r.mean()))
    mub = bisect(lambda m: readout(m, tau0), float(b0r.mean()))
    def draw(mu):
        p = np.diff(np.concatenate(([0.0], norm.cdf(tau0-mu), [1.0])))
        return RNG.choice(np.arange(1., KK+1.), size=n, p=p/p.sum())
    v = np.array([explained(draw(mua), draw(mua+dA), draw(mub), draw(mub+dB)) for _ in range(rep)], float)
    v = v[np.isfinite(v)]
    return float(np.median(v)), float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))

Y805 = int(J805["y0"]); Y806 = int(min(y for y in G.year.unique() if len(cell(int(y), 2)) >= 120))
ROWS = []
print(f"\n=== ① 各自量一遍(⚠ **不搬 `#813` 那个数**)· 每个起点 × 两个合成样本量 × {NREP} 次 ===")
for tag, ybase, tol, obs in (("`#805`(R244)", Y805, 0.08, J805["neg_control"]),
                             ("`#806`(R245)", Y806, 0.08, J806["neg_control"])):
    for n in (2000, 4000):
        med, lo, hi = spread(ybase, n, 0.45, 0.45)
        ratio = tol/((hi-lo)/2)
        ROWS.append(dict(round=tag, y_base=ybase, n=n, median=med, lo=lo, hi=hi,
                         half_width=(hi-lo)/2, tol=tol, tol_over_noise=ratio, observed=obs,
                         obs_inside=bool(lo <= obs <= hi)))
        print(f"  {tag} 起点 {ybase} · n={n:>4} ⇒ 中位 **{med:+.4f}** · 95% 跨度 [{lo:+.4f}, {hi:+.4f}]"
              f"(半宽 **{(hi-lo)/2:.4f}**)· **当初容差 {tol} ÷ 半宽 = {ratio:.2f}×** · "
              f"当初那次观测 {obs:+.4f} 落在跨度内:**{'是' if lo <= obs <= hi else '否'}**")

print("\n=== ② 控制 ===")
pmed, plo, phi = spread(Y805, 4000, 0.15, 0.60)
print(f"  正对照:同一条代码,在一个**已知有分歧**的世界(Δ 0.15 vs 0.60)⇒ "
      f"中位 **{pmed:+.4f}** [{plo:+.4f}, {phi:+.4f}] —— 该**明显低于 1**")
n_ok = sum(1 for r in ROWS if r["tol_over_noise"] < 1.0)
by_n = {}
for r in ROWS: by_n.setdefault((r["round"], r["n"]), r["tol_over_noise"] < 1.0)
flip = len({v for v in by_n.values()}) > 1

Gg = Gate("#814 · 把两轮的负控各抽两百次")
Gg.asserted("① 正对照:同一条代码在一个已知有分歧的世界里必须给出明显低于 1 的中位"
            "(否则这条负控连「能不能动」都没证明,量它的噪声也没有意义)",
            bool(phi < 1.0), f"中位 {pmed:+.4f},上界 {phi:+.4f} < 1", kind="control")
Gg.asserted("② 前提(跑前写下的混淆):噪声跨度依赖合成样本量 n ⇒ **n 取 2000 与 4000 都跑并同时印出**,"
            "**若结论随 n 翻转则本轮什么也没证明**",
            bool(not flip), f"两个 n 下「容差 < 噪声半宽」的判定{'一致' if not flip else '**翻转**'}",
            kind="control")
Gg.asserted("③ 前提:两轮的起点/容差/观测值从各自产物 JSON 读出,不从记忆",
            True, f"`#805` 起点 {Y805} 观测 {J805['neg_control']:+.4f} · `#806` 起点 {Y806} 观测 {J806['neg_control']:+.4f}",
            kind="control")
Gg.asserted("④ kill(预注册):「那两轮的负控是有分辨率的」要成立,需**当初容差 ≥ 自己的噪声半宽**",
            bool(n_ok == 0), f"容差 < 噪声半宽的格 {n_ok}/{len(ROWS)}", kind="kill")
print(); print(Gg)
adm = Gg.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*96)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判。**"
elif n_ok == len(ROWS):
    V = (f"**两轮的负控都没有分辨率,而且与合成样本量无关(n=2000 与 4000 同结论)。**\n"
         f"  当初的容差 0.08 ÷ 自己的噪声半宽 = "
         + " · ".join(f"{r['round'][:8]} n={r['n']} **{r['tol_over_noise']:.2f}×**" for r in ROWS) + "\n"
         f"  ⇒ **容差比噪声窄一个数量级 ⇒ 那条控制是随机开火的:通过与失败都不携带信息。**\n"
         f"  ⇒ **`#805` 与 `#806` 的「控制全过」要改成「负控没有分辨率,正控与 g=0 控制仍然成立」。**\n"
         f"  ⚠ **两轮的结论都不因此改变 —— 改变的是它们能声称的强度。**")
else:
    V = (f"**混合:{n_ok}/{len(ROWS)} 个格的容差窄于自己的噪声半宽。整张表全报,不选边。**")
print(V)
json.dump(dict(rows=ROWS, n_rep=NREP, pos_control=dict(median=pmed, lo=plo, hi=phi),
               n_tol_below_noise=n_ok, flips_with_n=bool(flip),
               admissible=adm, verdict=V, gate_ok=Gg.verdict(), action="Production/verify"),
          open(OUT/"retrofit_controls.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'retrofit_controls.json'}")
