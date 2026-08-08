"""#801 · E03·A53·R240 —— 页面说顶对「改得更多」,却从没说过他们从哪儿到哪儿

`#799` 之后,页面上最硬的一句是:**顶对(`sexeduc` 2.19 · `racmar` 1.84)上,虔诚者改得比其余人还多。**
**⚠⚠ 而那句话是一个比值,而比值分不开三个完全不同的故事:**
   ① **追上**:虔诚者从极端出发,只是赶到了其余人早就在的地方;
   ② **并肩**:两边一起动,虔诚者动得更快;
   ③ **超过**:虔诚者越过了其余人。
**页面从来没报过四个层的首末水平** ⇒ **读者手上只有一个比值,而那三个故事都能生成同一个比值。**

**⚠⚠ 而更硬的一件事:`#784`/`#785` 那套「天花板逼出来的比」对照,从来没跑到这两题上。**
实测(`#785` 的产物):它只跑过 `homosex`·`premarsx`·`teensex`·`xmarsex` 四题。
⇒ **页面上唯一还站着的那一对,恰恰是那条对照唯一没覆盖的地方。**
⚠ 而 `#784` 那套是问「虔诚者靠近**天花板**、动不了那么多」;
**这一对要问的是它的镜像:虔诚者的起跑线离得**远**,所以本来就有更多格子可走。**
**同一具机器,反过来用 —— 而它一次都没被反过来用过。**

G1 估计量(两个,方法之前先命名):
   (a) **四个层的首末水平** —— 虔诚/非虔诚 × 首年/末年,**在原始量表上**,不是比值。
   (b) **`r_forced`** —— 两层共用阈值、各自起跑线、施加同一个潜在位移 Δ,由起跑线自己走出的比。
       **Δ 由非虔诚层自己观测到的变化标定**(`#785` 修正后的**斜率口径**,不是首末差)。

⚠⚠ 三个世界,而第二个会把页面最响的那句话变成算术:
   A **真的多改了**:`r_obs` 的区间**排除** `r_forced` 且在其上
     ⇒ 他们改得比起跑线所能解释的还多,**那句话作为心理学主张站得住**。
   B **起跑线产物**:区间**含** `r_forced` ⇒ **「改得更多」是有界尺子在那条起跑线上必然给出的**,
     ⇒ 页面最响的一句要缩成「他们从更远的地方出发」。
   C **反号**:区间在 `r_forced` **之下** ⇒ 他们改得比起跑线预言的还**少** —— 那会把话整个翻过来。

预测矩阵:
   | 世界 | 现在 | 若区间排除且在上 | 若含 | 若在下 |
   | A | 0.45 | **0.90** | 0.05 | 0.05 |
   | B | 0.40 | 0.05 | **0.90** | 0.10 |
   | C | 0.15 | 0.05 | 0.05 | **0.85** |

预注册判词(条件式):
  if 正控开火(两层起跑线设成相同时 `r_forced` 必须回到 **1.0**)
     and 负控开火(起跑线拉得越远,`r_forced` 必须**单调变大** —— 否则机器没在响应起跑线):
      `r_obs` 区间排除 `r_forced` 且在其上 -> A
      区间含 `r_forced`                    -> B
      区间在其下                            -> C
  else: UNVERIFIED
⚠ **「这个零该不该是零?」** —— 正控的参照是 **1.0 不是 0**(起跑线相同 ⇒ 同样的位移走出同样的量)
  ⇒ **`identity_control`**,不是 `negative_control`,也不是 `offset_control`
  (`#796` 在这里错过一次:`offset_control` 问「越没越过」,我要的是「等不等于」)。

⚠ 跑之前写下的最强混淆:**`racmar` 只有两档(yes/no)、`sexeduc` 三档**,
  而潜变量模型的阈值数 = 档数 − 1 ⇒ **两题的模型自由度不同**,`r_forced` 的可比性受限。
  ⇒ 控制:**两题各自报,不合并成一个数**;并**报每题的档数与阈值数**,让读者看见这件事。

⚠ 硬规则①:两题的 **n 与真正被问过的年份**先打印再用;题干原文从 `.dta` 的变量标签读。

本轮换不了仪器(对象是世界;第二具仪器本机六具全部落选 —— `R223/instrument_search.py`)。
⚠ 总判由 `Gate.admissible()` 决定(第四次用)。
"""
import numpy as np, json, pathlib, sys
from scipy.stats import norm, logistic
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import pandas as pd
from lib.gates import Gate

RNG = np.random.default_rng(240)
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(parents=True, exist_ok=True)
gp = ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
TOP = ["sexeduc", "racmar"]
STEM = pd.io.stata.StataReader(gp).variable_labels()

print("=== ⓪ 硬规则①:先打印 n、真正被问过的年份、档数、题干原文 ===")
d = pd.read_stata(gp, columns=["year", "attend", "reliten", "fund"]+TOP, convert_categoricals=False)
cat = pd.read_stata(gp, columns=TOP, convert_categoricals=True)
K = {}
for c in TOP:
    cs = list(cat[c].cat.categories); K[c] = len(cs)
    v = pd.to_numeric(d[c], errors="coerce").where(lambda x: (x >= 1) & (x <= len(cs)))
    yrs = sorted(d.year[v.notna()].unique())
    print(f"  {c:8s} n={int(v.notna().sum()):>7,} · 年 {len(yrs):>2} ({int(min(yrs))}–{int(max(yrs))}) · "
          f"档 {len(cs)}(阈值 {len(cs)-1}) · {cs}")
    print(f"           题干:「{STEM.get(c,'?')}」")

M = pd.DataFrame({c: pd.to_numeric(d[c], errors="coerce").where(
    lambda v, c=c: (v >= 1) & (v <= K[c])) for c in TOP})
for c, (lo, hi) in (("attend", (0, 8)), ("reliten", (1, 4)), ("fund", (1, 3))):
    M[c] = pd.to_numeric(d[c], errors="coerce").where(lambda v, lo=lo, hi=hi: (v >= lo) & (v <= hi))
M["year"] = d.year
M["reliten"] = -M["reliten"]; M["fund"] = -M["fund"]
zs = lambda s: (s-s.mean())/s.std(ddof=1)
REL = M.dropna(subset=["attend", "reliten", "fund", "year"]).copy()
REL["REL"] = zs(REL[["attend", "reliten", "fund"]]).mean(axis=1)
REL["k"] = REL.groupby("year")["REL"].transform(lambda v: pd.qcut(v, 3, labels=False, duplicates="drop"))

def slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.cov(x, y, ddof=1)[0, 1]/np.var(x, ddof=1))
def series(item, k, nmin=120):
    g = REL.dropna(subset=[item])
    return [(int(y), float(gy[item].mean()), len(gy)) for y, gy in g[g.k == k].groupby("year") if len(gy) >= nmin]
def dist(item, k, year):
    g = REL[(REL.k == k) & (REL.year == year)].dropna(subset=[item])
    return np.array([(g[item] == c).mean() for c in range(1, K[item]+1)]), len(g)

# ── ① 四个层的首末水平 —— 页面从没报过的那件事 ──────────────────────────────
print("\n=== ① 四个层的首末水平(原始量表,不是比值)—— 页面从没报过 ===")
LEV = {}
for c in TOP:
    A, B = series(c, 2), series(c, 0)
    y0, y1 = max(A[0][0], B[0][0]), min(A[-1][0], B[-1][0])
    a0 = next(r for r in A if r[0] == y0); a1 = next(r for r in A if r[0] == y1)
    b0 = next(r for r in B if r[0] == y0); b1 = next(r for r in B if r[0] == y1)
    LEV[c] = dict(y0=y0, y1=y1, dev0=a0[1], dev1=a1[1], non0=b0[1], non1=b1[1],
                  n_dev0=a0[2], n_non0=b0[2], gap0=a0[1]-b0[1], gap1=a1[1]-b1[1],
                  d_dev=a1[1]-a0[1], d_non=b1[1]-b0[1], K=K[c])
    L = LEV[c]
    print(f"  {c:8s} {y0}→{y1}(量表 1–{K[c]})")
    print(f"    虔诚层  {L['dev0']:.3f} → {L['dev1']:.3f}   (Δ {L['d_dev']:+.3f})")
    print(f"    非虔诚层 {L['non0']:.3f} → {L['non1']:.3f}   (Δ {L['d_non']:+.3f})")
    print(f"    两层差距 {L['gap0']:+.3f} → {L['gap1']:+.3f}  ⇒ "
          f"**{'差距缩小(追上)' if abs(L['gap1']) < abs(L['gap0']) else '差距扩大或反转'}**"
          f"{' · **越过了对方**' if np.sign(L['gap1']) != np.sign(L['gap0']) and abs(L['gap1'])>1e-9 else ''}")

# ── ② r_forced ── #785 修正后的斜率口径 ───────────────────────────────────────
def fit_tau(p, link): return link.ppf(np.clip(np.cumsum(p)[:-1], 1e-6, 1-1e-6))
def readout(mu, tau, link, kk):
    e = link.cdf(tau-mu); p = np.diff(np.concatenate(([0.0], e, [1.0])))
    return float((p*np.arange(1, kk+1)).sum())
def fit_mu(p, tau, link, kk):
    tgt = float((p*np.arange(1, kk+1)).sum()); lo, hi = -8.0, 8.0
    for _ in range(80):
        mid = (lo+hi)/2
        if readout(mid, tau, link, kk) < tgt: lo = mid
        else: hi = mid
    return (lo+hi)/2
def forced(pA, pB, nA, nB, link, kk, dB):
    tau = fit_tau((pA*nA + pB*nB)/(nA+nB), link)
    muA, muB = fit_mu(pA, tau, link, kk), fit_mu(pB, tau, link, kk)
    oA, oB = readout(muA, tau, link, kk), readout(muB, tau, link, kk)
    lo, hi = -8.0, 8.0                      # ⚠ 两侧(`#784` 第一版写成 [0,6] 被对照打掉)
    for _ in range(80):
        mid = (lo+hi)/2
        if readout(muB+mid, tau, link, kk)-oB < dB: lo = mid
        else: hi = mid
    return (readout(muA+(lo+hi)/2, tau, link, kk)-oA)/dB

LINKS = {"probit": norm, "logit": logistic}
print("\n=== ② `r_forced`:同一个潜在位移,由两条起跑线自己走出的比(`#785` 修正后的斜率口径)===")
RES = {}
for c in TOP:
    A, B = series(c, 2), series(c, 0)
    yA = np.array([r[0] for r in A], float); vA = np.array([r[1] for r in A])
    yB = np.array([r[0] for r in B], float); vB = np.array([r[1] for r in B])
    spA, spB = yA[-1]-yA[0], yB[-1]-yB[0]
    r_obs = (slope(yA, vA)*spA)/(slope(yB, vB)*spB)
    dB = slope(yB, vB)*spB
    pA0, nA0 = dist(c, 2, int(yA[0])); pB0, nB0 = dist(c, 0, int(yB[0]))
    f = lambda ia, ib: (slope(yA[ia], vA[ia])*spA)/(slope(yB[ib], vB[ib])*spB)
    bs = np.array([f(RNG.integers(0, len(yA), len(yA)), RNG.integers(0, len(yB), len(yB))) for _ in range(4000)])
    bs = bs[np.isfinite(bs)]
    lo95, hi95 = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    rf = {ln: float(forced(pA0, pB0, nA0, nB0, lk, K[c], dB)) for ln, lk in LINKS.items()}
    RES[c] = dict(r_obs=float(r_obs), lo=lo95, hi=hi95, forced=rf, K=K[c],
                  above=bool(lo95 > max(rf.values())), contains=bool(any(lo95 <= v <= hi95 for v in rf.values())),
                  below=bool(hi95 < min(rf.values())))
    print(f"  {c:8s}(档 {K[c]})  r_obs **{r_obs:+.3f}** [{lo95:+.3f}, {hi95:+.3f}] · "
          f"r_forced probit {rf['probit']:+.3f} · logit {rf['logit']:+.3f}  ⇒ "
          f"**{'区间在其上' if RES[c]['above'] else ('区间含它' if RES[c]['contains'] else '区间在其下')}**")

# ── ③ 控制 ────────────────────────────────────────────────────────────────────
print("\n=== ③ 控制:起跑线相同 ⇒ 必须回到 1.0(参照 1.0 不是 0);起跑线越远 ⇒ 必须单调变大 ===")
# ⚠⚠ 第一版把控制的 `dB` 硬写成 **−0.20**,而 `sexeduc` 的非虔诚层坐在 **1.110**(量表 1–3)
#    ⇒ **它向下最多只能走 0.110。我要求这个模型产生一个尺子给不出的变化。**
#    Δ 的二分因此顶在下界饱和,比值退化 ⇒ 正控给 0.5502 而不是 1.0,负控也跟着不单调。
#    ⇒ **两条控制都是「因自己的理由而失败」** —— 失败的是我写的控制,不是被控的机器。
#    ⇒ 改成:`dB` 取**该起跑线在行进方向上可达幅度的 30%**,保证二分落在内部。
eq, mono = [], []
c = TOP[0]
pB0_, nB0_ = dist(c, 0, LEV[c]["y0"])
reach = float((pB0_*np.arange(1, K[c]+1)).sum()) - 1.0     # 向下(朝 1)可达的幅度
dB_ctl = -0.30*reach
print(f"  ⚠ 控制用的 dB = **{dB_ctl:+.4f}**(= 该起跑线向下可达幅度 {reach:.4f} 的 30%,保证不饱和)")
for ln, lk in LINKS.items():
    pB, nB = pB0_, nB0_
    eq.append(forced(pB, pB, nB, nB, lk, K[c], dB_ctl))
    tau0 = fit_tau(pB, lk); mu0 = fit_mu(pB, tau0, lk, K[c])
    seq = []
    for shift in (0.0, 0.4, 0.8, 1.2):
        e = lk.cdf(tau0-(mu0-shift)); pA = np.diff(np.concatenate(([0.0], e, [1.0])))
        seq.append(forced(pA, pB, nB, nB, lk, K[c], dB_ctl))
    mono.append(all(seq[i] < seq[i+1] for i in range(3)))
    print(f"  {ln}: 起跑线相同 → {eq[-1]:.4f} · 起跑线越拉越远 → {[round(s,3) for s in seq]} "
          f"{'单调变大' if mono[-1] else '**不单调**'}")

G = Gate("#801 · 顶对从哪儿到哪儿,而那条起跑线解释得了多少")
G.identity_control("① 正控:两层起跑线设成相同时,`r_forced` 必须回到 1.0(参照 1.0,不是 0)",
                   observed=float(np.median(eq)), expected=1.0, tol=0.02,
                   what="起跑线相同的合成世界 —— 同样的潜在位移走出同样的量,参照是 1.0")
G.asserted("② 负控:起跑线拉得越远,`r_forced` 必须**单调变大**(否则机器没在响应起跑线)",
           bool(all(mono)), f"probit/logit 单调:{mono}", kind="control")
G.asserted("③ 前提(跑前写下的混淆):两题档数不同(3 与 2 ⇒ 阈值 2 与 1)⇒ **各自报,不合并**",
           bool(K["sexeduc"] != K["racmar"]),
           f"sexeduc 档 {K['sexeduc']} · racmar 档 {K['racmar']} —— 两题分开报,不合成一个数", kind="control")
both_above = all(RES[c]["above"] for c in TOP)
G.asserted("④ kill(预注册):「改得更多」要作为心理学主张站住,需**两题**的区间都排除 `r_forced` 且在其上",
           both_above, " · ".join(f"{c} {'在上' if RES[c]['above'] else ('含它' if RES[c]['contains'] else '在下')}"
                                  for c in TOP), kind="kill")
print(); print(G)
adm = G.admissible()
print(f"\n  `Gate.admissible()` = **{adm}**")

print("\n"+"="*92)
if not adm:
    v = "**UNVERIFIED:控制行没有全过 ⇒ 仪器没资格下判。**"
elif both_above:
    v = (f"**A 他们改得比起跑线所能解释的还多。** 两题的观测区间都**排除**了「同样的位移由起跑线自己走出的比」"
         f"并落在其上:`sexeduc` [{RES['sexeduc']['lo']:+.3f}, {RES['sexeduc']['hi']:+.3f}] vs "
         f"forced {RES['sexeduc']['forced']['probit']:+.3f}/{RES['sexeduc']['forced']['logit']:+.3f};"
         f"`racmar` [{RES['racmar']['lo']:+.3f}, {RES['racmar']['hi']:+.3f}] vs "
         f"{RES['racmar']['forced']['probit']:+.3f}/{RES['racmar']['forced']['logit']:+.3f}。")
elif any(RES[c]["contains"] for c in TOP):
    bad = [c for c in TOP if RES[c]["contains"]]
    v = (f"**B 起跑线产物 —— 而这是我不想要的那一支。** {bad} 的观测区间**包含**「起跑线自己走出的比」"
         f"⇒ **在这一题上,「虔诚者改得更多」是有界尺子在那条起跑线上必然给出的**,\n"
         f"  页面最响的一句要缩成「他们从更远的地方出发」。")
else:
    v = (f"**C 反号:他们改得比起跑线预言的还少。** 区间落在 `r_forced` 之下 ⇒ 那句话要整个翻过来。")
print(v)
json.dump(dict(top=TOP, K=K, stems={c: STEM.get(c, "") for c in TOP}, levels=LEV, res=RES,
               eq_control=[float(x) for x in eq], mono=mono, dB_ctl=float(dB_ctl),
               reach=float(reach), admissible=adm,
               verdict=v, gate_ok=G.verdict()),
          open(OUT/"where_from_where_to.json", "w"), ensure_ascii=False, indent=1)
print(f"\n  产物 → {OUT/'where_from_where_to.json'}")
