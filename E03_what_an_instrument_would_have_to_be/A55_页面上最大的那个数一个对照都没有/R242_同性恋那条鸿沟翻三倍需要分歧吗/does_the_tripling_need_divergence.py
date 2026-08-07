"""#803 · E03·A55·R242 —— 同性恋那条鸿沟翻三倍,需要「分歧」吗?

`#802` 刚发表了本页目前最大的一个数:`homosex` 的两层差距 **−0.549 → −1.651,约三倍**,
而我在同一条账里写下 **「它一个控制都还没有」**。本轮就是那个控制,**而它是故意挑一个
「正结果我不欢迎」的方向**(`§3` 盆地规则):**若它开火,我上一轮的头条就死。**

⚠⚠ 核心问题,而它不是「差距有没有扩大」(那是观测,已经在页上):
   **一个共同的潜在位移 Δ —— 两层完全一样地改主意 —— 在一把 1–4 的序数尺上,
   从两条不同的起跑线读出来,本来就会产生一条扩大的差距吗?**
   若会,则 **「三倍」不需要任何分歧就能出现**,它是刻度的性质,不是人的性质。

⚠ 而 `#801` 刚验证过这套机器的一条反直觉性质:**地板与天花板一样压缩** ——
   离阈值远的质量不跨过阈值 ⇒ **「起跑线远所以有更多空间」是假的。**
   ⇒ **本轮的预期方向因此是不确定的,这正是它值得跑的原因**:
     我既没有「共同 Δ 会造出三倍」的直觉,也没有相反的直觉。

G1 估计量:**共同 Δ 能解释的那部分差距变化的份额**
   `explained = (预测 gap1 − gap0) / (观测 gap1 − gap0)`
   —— 定标臂按构造复现自己的观测变化,读的是**另一臂在同一个 Δ 下会走到哪**。
   ⚠ **分母是观测的差距变化(−1.102),不是 0** ⇒ 这不是「差与零比」。

⚠⚠ 「这个零该不该是零?」——**这里根本没有零。共同 Δ 的世界里 `explained` 的参照值是 1.0。**
   ⇒ 负控必须用 `identity_control` 对着 **1.0**,不是对着 0(与 `#801` 同一个坑)。

三个世界:
   A **真分歧**:共同 Δ 解释不了(`explained` 小)⇒ 两层确实走得不一样,**头条站得住**。
   B **刻度产物**:共同 Δ 就能造出这条扩大的差距(`explained` 接近 1)⇒
     **「三倍」不需要分歧 ⇒ 我上一轮的头条必须降级为「在这把尺上,同样的改变读出来就是这样」。**
   C **模型不适配**:预测差距与观测差距差得离谱(两个方向都不像)⇒
     **这台机器回答不了差距问题**,而那是关于仪器边界的真收获(登记,不硬判)。

预测矩阵:
   | 世界 | 现在 | explained ≥ 0.50 | explained < 0.50 | 预测离谱 |
   | A 真分歧   | 0.45 | 0.05 | **0.85** | 0.10 |
   | B 刻度产物 | 0.40 | **0.90** | 0.05 | 0.10 |
   | C 不适配   | 0.15 | 0.05 | 0.10 | **0.80** |

预注册判词(条件式,**先问控制,再看阈值**):
  if 正控开火(**造一个两层 Δ 真的不同的世界,判据必须给出低 explained**)
     and 负控开火(**造一个两层 Δ 完全相同的世界,`explained` 必须等于 1.0**):
      八格中位 explained >= 0.50 -> B,**头条降级**
      < 0.50                     -> A,**头条站住,且现在它有对照了**
  else: UNVERIFIED
⚠ 「预测离谱」(世界 C)不设阈值硬判 —— 它由**预测差距落在观测差距的哪一侧、差多远**据实登记。

⚠ 跑之前写下的最强混淆:**定标臂的选择本身可能决定结论。**
  把 Δ 定在非虔诚臂上(它动得多)与定在虔诚臂上(它动得少),是两个不同的问题。
  ⇒ 控制:**两个方向都跑,进同一张规格曲线**,而**不取平均**(`realstat` §2.5:分歧本身是发现)。
`G4` 规格曲线 = 连接函数(probit/logit)× 定标臂(非虔诚/虔诚)× 端点(拟合/原始)= **8 格,全报。**

⚠ 硬规则①:先打印 n、真正被问过的年份、档数。
⚠ 硬规则③:本轮换不了仪器(六具候选全部落选 —— `R223/instrument_search.py`),
   **而它不需要第二具仪器:这是一个「同一份数据在一个模型下会不会自己造出这个数」的问题。**
⚠ 总判由 `Gate.admissible()` 决定。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(242)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
IT = "homosex"; KK = 4          # ⚠ `#788`:原始文件有 5 类,第 5 类是 'other',不是一档
STEM = pd.io.stata.StataReader(gp).variable_labels()

print("=== ⓪ 硬规则①:n · 真正被问过的年份 · 档数 · 题干原文 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund", IT], convert_categoricals=False)
cats = list(pd.read_stata(gp, columns=[IT], convert_categoricals=True)[IT].cat.categories)
v = pd.to_numeric(d[IT], errors="coerce").where(lambda x: (x >= 1) & (x <= KK))
yrs = sorted(d.year[v.notna()].unique())
print(f"  {IT}  n={int(v.notna().sum()):,} · 年 {len(yrs)}({int(min(yrs))}–{int(max(yrs))}) · "
      f"用档 {KK}(阈值 {KK-1}) · 原始 {len(cats)} 类 {cats}")
print(f"  题干:「{STEM.get(IT,'?')}」")
print(f"  ⚠ 第 5 类 'other' 不入模型(`#788` 那一次就是把它一起传进去被 `label_pole` 挡下的)")

M = pd.DataFrame({IT: v})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda x, lo=lo, hi=hi: (x >= lo) & (x <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda s: pd.qcut(s, 3, labels=False, duplicates="drop"))

def series(k, nmin=120):
    g = REL.dropna(subset=[IT])
    return [(int(y), float(gy[IT].mean()), len(gy)) for y, gy in g[g.k == k].groupby("year") if len(gy) >= nmin]
def dist(k, year):
    g = REL[(REL.k == k) & (REL.year == year)].dropna(subset=[IT])
    return np.array([(g[IT] == c).mean() for c in range(1, KK+1)]), len(g)

A, B = series(2), series(0)                      # A = 虔诚(top tercile) · B = 非虔诚
y0, y1 = max(A[0][0], B[0][0]), min(A[-1][0], B[-1][0])
def endp(R, spec):
    x = np.array([r[0] for r in R], float); y = np.array([r[1] for r in R])
    if spec == "raw":
        return (float(next(r[1] for r in R if r[0] == y0)), float(next(r[1] for r in R if r[0] == y1)))
    b = np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1); a = y.mean()-b*x.mean()
    return (a+b*y0, a+b*y1)

print(f"\n=== ① 观测:两层水平与差距({y0}→{y1},量表 1–{KK})===")
OBSV = {}
for spec in ("fit", "raw"):
    a0, a1 = endp(A, spec); b0, b1 = endp(B, spec)
    OBSV[spec] = dict(a0=a0, a1=a1, b0=b0, b1=b1, gap0=a0-b0, gap1=a1-b1,
                      dA=a1-a0, dB=b1-b0, dgap=(a1-b1)-(a0-b0))
    o = OBSV[spec]
    print(f"  [{spec}] 虔诚 {a0:.3f}→{a1:.3f}(Δ{o['dA']:+.3f}) · 非虔诚 {b0:.3f}→{b1:.3f}(Δ{o['dB']:+.3f})")
    print(f"        差距 {o['gap0']:+.3f} → {o['gap1']:+.3f}   **差距变化 {o['dgap']:+.3f}**")

# ── 共享阈值 · 各层自己的潜在均值 · 一个共同的 Δ ─────────────────────────────
def fit_tau(p, link): return link.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def readout(mu, tau, link):
    e = link.cdf(tau-mu); p = np.diff(np.concatenate(([0.0], e, [1.0])))
    return float((p*np.arange(1, KK+1)).sum())
def fit_mu(p, tau, link):
    tgt = float((p*np.arange(1, KK+1)).sum()); lo, hi = -8.0, 8.0
    for _ in range(90):
        mid = (lo+hi)/2
        if readout(mid, tau, link) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2

def common_delta(pA, nA, pB, nB, link, d_cal, cal):
    """共享阈值(首年合并分布)· 各层自己的 μ · 定标臂复现自己的观测变化 ⇒ 同一个 Δ 读另一臂。

    返回:定标臂与另一臂在同一个 Δ 下的**末年读数**,以及 Δ 本身。
    """
    tau = fit_tau((pA*nA + pB*nB)/(nA+nB), link)
    mu = {"A": fit_mu(pA, tau, link), "B": fit_mu(pB, tau, link)}
    o = {s: readout(mu[s], tau, link) for s in mu}
    lo, hi = -8.0, 8.0                          # ⚠ 两侧(`#784` 第一版写成 [0,6] 被对照打掉)
    for _ in range(90):
        mid = (lo+hi)/2
        if readout(mu[cal]+mid, tau, link) - o[cal] < d_cal: lo = mid
        else: hi = mid
    D = (lo+hi)/2
    oth = "B" if cal == "A" else "A"
    return dict(delta=D, cal_end=readout(mu[cal]+D, tau, link),
                oth_end=readout(mu[oth]+D, tau, link), cal_start=o[cal], oth_start=o[oth])

pA0, nA0 = dist(2, y0); pB0, nB0 = dist(0, y0)
LINKS = {"probit": norm, "logit": logistic}
print(f"\n=== ② `G4` 规格曲线:连接 × 定标臂 × 端点 = 8 格,全报(`G3` 多重性:8 格全部公布)===")
CELLS = []
for link_name, link in LINKS.items():
    for cal, cal_name in (("B", "非虔诚"), ("A", "虔诚")):
        for spec in ("fit", "raw"):
            o = OBSV[spec]
            d_cal = o["dB"] if cal == "B" else o["dA"]
            r = common_delta(pA0, nA0, pB0, nB0, link, d_cal, cal)
            # 预测的末年差距:定标臂按构造复现观测,另一臂由共同 Δ 决定
            if cal == "B":
                pred_a1, pred_b1 = r["oth_end"], r["cal_end"]
            else:
                pred_a1, pred_b1 = r["cal_end"], r["oth_end"]
            pred_gap1 = pred_a1 - pred_b1
            pred_gap0 = (r["oth_start"]-r["cal_start"]) if cal == "B" else (r["cal_start"]-r["oth_start"])
            expl = (pred_gap1 - pred_gap0)/(o["gap1"] - o["gap0"])
            CELLS.append(dict(link=link_name, cal=cal_name, spec=spec, delta=r["delta"],
                              pred_gap0=pred_gap0, pred_gap1=pred_gap1,
                              obs_gap0=o["gap0"], obs_gap1=o["gap1"], explained=expl))
            print(f"  {link_name:6s} 定标={cal_name:3s} {spec:3s} · Δ={r['delta']:+.3f} · "
                  f"预测差距 {pred_gap0:+.3f} → {pred_gap1:+.3f}(观测 {o['gap0']:+.3f} → {o['gap1']:+.3f}) "
                  f"⇒ **explained = {expl:+.3f}**")
EXPL = np.array([c["explained"] for c in CELLS])
med = float(np.median(EXPL))
print(f"\n  八格 explained:中位 **{med:+.3f}** · 全距 [{EXPL.min():+.3f}, {EXPL.max():+.3f}] · "
      f"≥0.50 的格 **{int((EXPL >= 0.5).sum())}/8**")

# ── ③ 控制:合成世界 ───────────────────────────────────────────────────────
print("\n=== ③ 控制(合成世界,同一条代码路径)===")
def synth(mu_a, mu_b, dA, dB, link):
    """造一个世界:共享阈值来自真实首年合并分布,两层各自的 μ 与各自的 Δ 由我指定。"""
    tau = fit_tau((pA0*nA0 + pB0*nB0)/(nA0+nB0), link)
    pa = np.diff(np.concatenate(([0.0], link.cdf(tau-mu_a), [1.0])))
    pb = np.diff(np.concatenate(([0.0], link.cdf(tau-mu_b), [1.0])))
    a0, b0 = readout(mu_a, tau, link), readout(mu_b, tau, link)
    a1, b1 = readout(mu_a+dA, tau, link), readout(mu_b+dB, tau, link)
    return pa, pb, dict(a0=a0, a1=a1, b0=b0, b1=b1, gap0=a0-b0, gap1=a1-b1)

def run_syn(dA, dB, link, cal="B"):
    tau = fit_tau((pA0*nA0 + pB0*nB0)/(nA0+nB0), link)
    mu_a, mu_b = fit_mu(pA0, tau, link), fit_mu(pB0, tau, link)
    pa, pb, o = synth(mu_a, mu_b, dA, dB, link)
    r = common_delta(pa, nA0, pb, nB0, link, o["b1"]-o["b0"] if cal == "B" else o["a1"]-o["a0"], cal)
    pred_a1, pred_b1 = (r["oth_end"], r["cal_end"]) if cal == "B" else (r["cal_end"], r["oth_end"])
    pg0 = (r["oth_start"]-r["cal_start"]) if cal == "B" else (r["cal_start"]-r["oth_start"])
    den = o["gap1"]-o["gap0"]
    return float(((pred_a1-pred_b1) - pg0)/den) if abs(den) > 1e-9 else float("nan"), o

# 负控:两层 Δ **完全相同** ⇒ explained 必须 = 1.0(⚠ 参照是 1.0,不是 0)
nc_expl, nc_o = run_syn(+0.60, +0.60, norm)
print(f"  负控:两层 Δ 完全相同(+0.600/+0.600)· 合成差距 {nc_o['gap0']:+.3f} → {nc_o['gap1']:+.3f} "
      f"⇒ explained = **{nc_expl:+.6f}**(参照 **1.0**,不是 0)")
# 正控:两层 Δ **真的不同** ⇒ explained 必须**低**(判据要认得出分歧)
pc_expl, pc_o = run_syn(+0.20, +0.90, norm)
print(f"  正控:两层 Δ 真的不同(虔诚 +0.200 / 非虔诚 +0.900)· 合成差距 "
      f"{pc_o['gap0']:+.3f} → {pc_o['gap1']:+.3f} ⇒ explained = **{pc_expl:+.6f}**(该**远低于** 1)")
# ⚠ 正控必须**在 g=0 时不开火** —— 上面那条负控就是 g=0 的那一格,两者互为补
print(f"  ⚠ 正控在 g=0(即负控那一格)不开火:{nc_expl:.6f} 不低于 0.50 ⇒ "
      f"**{'是' if nc_expl >= 0.5 else '否 —— 判据在无分歧时也喊分歧,不可用'}**")

G = Gate("#803 · 同性恋那条鸿沟翻三倍,需要分歧吗")
G.identity_control("① 负控:两层 Δ 完全相同的合成世界里,`explained` 必须**等于 1.0**"
                   "(⚠ 参照是 1.0 不是 0 —— 与 `#801` 同一个坑)",
                   observed=float(nc_expl), expected=1.0, tol=1e-3,
                   what="共享阈值、两层各自的 μ、两层施加同一个 Δ 的合成世界")
G.asserted("② 正控:两层 Δ 真的不同的合成世界里,判据必须给出**远低于 1** 的 explained",
           bool(pc_expl < 0.50), f"explained = {pc_expl:+.4f}(阈 0.50)", kind="control")
G.asserted("③ 正控必须在 g=0 时**不**开火(否则它在无分歧时也喊分歧)",
           bool(nc_expl >= 0.50), f"g=0 那一格 explained = {nc_expl:+.4f}", kind="control")
G.asserted("④ 前提(跑前写下的混淆):定标臂两个方向都跑,进同一张规格曲线,**不取平均**",
           bool(len({c['cal'] for c in CELLS}) == 2 and len(CELLS) == 8),
           f"{len(CELLS)} 格 = 连接 2 × 定标臂 {len({c['cal'] for c in CELLS})} × 端点 2", kind="control")
G.asserted("⑤ kill(预注册):「三倍是刻度产物、头条要降级」需八格中位 explained ≥ 0.50",
           bool(med >= 0.50), f"中位 explained = {med:+.3f} · ≥0.50 的格 {int((EXPL>=0.5).sum())}/8",
           kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*94)
if not adm:
    V = "**UNVERIFIED:控制行没有全过 ⇒ 判据没资格下判,而这一轮不许报那个 explained。**"
elif med >= 0.50:
    V = (f"**B 刻度产物 —— 而这杀掉的是我上一轮刚发表的头条。** 八格中位 explained = **{med:+.3f}**,\n"
         f"  ⇒ **一个共同的潜在位移 —— 两层完全一样地改主意 —— 在这把 1–4 的序数尺上,\n"
         f"  从两条不同的起跑线读出来,本来就会造出这条扩大的差距的大部分。**\n"
         f"  ⇒ **「同性恋那条鸿沟翻了三倍」必须降级:它不需要分歧就能出现。**")
else:
    V = (f"**A 真分歧 —— 头条站住了,而现在它有对照了。** 八格中位 explained = **{med:+.3f}**,\n"
         f"  ⇒ **一个共同的潜在位移造不出观测到的这条差距** ⇒ 两层确实走得不一样。\n"
         f"  ⇒ **`#802` 的那句话保留,并且从此带着这个对照:三倍不是刻度造的。**\n"
         f"  ⚠ 而 explained 仍然不是 0 ⇒ **刻度贡献了其中一部分,报这个数,不报「全是分歧」。**")
print(V)
json.dump(dict(item=IT, K=KK, y0=y0, y1=y1, n=int(v.notna().sum()), years=len(yrs),
               observed=OBSV, cells=CELLS, explained_median=med,
               explained_range=[float(EXPL.min()), float(EXPL.max())],
               n_ge_half=int((EXPL >= 0.5).sum()),
               neg_control=dict(explained=float(nc_expl), reference=1.0, **{k: float(x) for k, x in nc_o.items()}),
               pos_control=dict(explained=float(pc_expl), **{k: float(x) for k, x in pc_o.items()}),
               admissible=adm, verdict=V, gate_ok=G.verdict()),
          open(OUT/"tripling_needs_divergence.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'tripling_needs_divergence.json'}")
