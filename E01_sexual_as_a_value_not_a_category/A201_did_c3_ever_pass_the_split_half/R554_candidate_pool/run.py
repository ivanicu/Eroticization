"""E01·A201·R554 — 把候选池列出来,再从池子里随机抽参照

`#509` 的 NEXT,**它指向我不愿意看到的结果**:若跨池摆动明显大于池内摆动,
则 `#485a` 以来所有参照门槛都要带上「参照是从哪个池子挑的」。

⛔ 池子的定义**写在跑之前**,且「与性无关」用**预先声明的正则**做成机械筛,不是我逐个挑:
  ① 数值型,**3–12 个不同取值**(与各轮参照的形状一致);
  ② 在分析样本(有 `pornlaw` 与 `xmovie` 的人)上 **n ≥ 2000**;
  ③ 标签**不**匹配 `SEXRE`(性/道德/宗教/政治态度词表,下方写死);
  ④ 排除 `condemn` 自身的来源(`pornlaw`)与结局(`xmovie`)。
⚠ §P7:候选列数可能上百,而 `.dta` 是 598 MB ⇒ **先量读一批列的耗时**,再定池子上限。

G1 ESTIMAND:从池子里**无放回随机抽 12 个**当参照,算 q95 与中位;重复 2000 次 ->
  **跨池摆动**的 95% 区间宽度。与 `#509a` 的**池内**宽度(k=12:q95 0.0888 · 中位 0.0987)并排比。
判据(预注册):**跨池宽度 / 池内宽度 ≥ 1.5** 记为「明显更大」。

WORLDS:
  W-POOL-OK  跨池 ≈ 池内 -> 选择不确定性确实小,RULE-v3 稳固
  W-POOL-BIG 跨池明显更大 -> **门槛的真实不确定性一直被低估**,所有参照门槛加注
  | World      | now | 跨池≥1.5× | 相当 |
  | W-POOL-OK  | 0.4 | 0.15      | 0.85 |
  | W-POOL-BIG | 0.6 | 0.85      | 0.15 |
CONTROLS:正对照 = **池子必须非空且 ≥30 个变量**(否则「随机抽 12」没有意义);
  阴性 = 用**池子全体**算的 q95/中位,跨池抽样的分布应当**围绕它**(无系统偏移)。
IMPOSSIBLE:关键词表是我写的 ⇒ 「与性无关」的边界仍带我的判断,**但它是机械的、可复算的** ·
  未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, time, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
DTA = ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta"
SEXRE = re.compile(r"sex|porn|homosex|gay|lesbian|abort|marri|divorc|birth|contracep|"
                   r"relig|god|bible|pray|church|attend|moral|wrong|premar|xmar|teensex|"
                   r"polit|party|vote|liberal|conserv|gun|death penalt|capital pun", re.I)
POOL_CAP = 300     # §P7:池子上限,先量耗时再定

it = pd.read_stata(DTA, iterator=True); vl = it.variable_labels()
t0 = time.time()
base = pd.read_stata(DTA, columns=["pornlaw", "xmovie"], convert_categoricals=False)
print(f"读 2 列耗时 {time.time()-t0:.2f}s -> {POOL_CAP} 列预计 {POOL_CAP*(time.time()-t0)/2/60:.1f} 分钟")
ana = base.dropna(subset=["pornlaw", "xmovie"]).index
condemn_full = (base.pornlaw == 1).astype(float)

cands = [c for c, lab in vl.items()
         if not SEXRE.search(str(lab)) and c not in ("pornlaw", "xmovie")]
print(f"标签不匹配性/道德词表的变量:{len(cands)}(总 {len(vl)})")
cands = cands[:POOL_CAP]
t0 = time.time()
df = pd.read_stata(DTA, columns=cands, convert_categoricals=False)
print(f"读 {len(cands)} 列耗时 {time.time()-t0:.1f}s")

pool, rs = [], []
for c in cands:
    s = df[c]
    if not np.issubdtype(s.dtype, np.number): continue
    sub = s.loc[ana]
    m = sub.notna().values
    if m.sum() < 2000: continue
    u = sub[m].nunique()
    if not (3 <= u <= 12): continue
    r = float(np.corrcoef(condemn_full.loc[ana].values[m], sub[m].values)[0, 1])
    if np.isfinite(r): pool.append(c); rs.append(abs(r))
rs = np.array(rs)
print(f"\n池子大小 = {len(pool)};|r| 中位 {np.median(rs):.4f} q95 {np.quantile(rs,.95):.4f} "
      f"max {rs.max():.4f}")
print("  前 6 个:", pool[:6])

rng = np.random.default_rng(20260805)
def width(v):
    lo, hi = np.quantile(v, [.025, .975]); return float(hi - lo)
q_draws = [np.quantile(rng.choice(rs, 12, replace=False), .95) for _ in range(2000)]
m_draws = [np.median(rng.choice(rs, 12, replace=False)) for _ in range(2000)]
W_IN_Q, W_IN_M = 0.0888, 0.0987      # `#509a` 的池内宽度(k=12)
wq, wm = width(q_draws), width(m_draws)
print(f"\n跨池抽 12:q95 宽度={wq:.4f}(池内 {W_IN_Q:.4f},比 {wq/W_IN_Q:.2f})  "
      f"中位 宽度={wm:.4f}(池内 {W_IN_M:.4f},比 {wm/W_IN_M:.2f})")
print(f"  跨池 q95 中位={np.median(q_draws):.4f} · 全池 q95={np.quantile(rs,.95):.4f}")
print(f"  跨池 中位 中位={np.median(m_draws):.4f} · 全池 中位={np.median(rs):.4f}")

G = Gate("把候选池列出来,再从池子里随机抽参照")
G.asserted("正对照:池子 ≥30 个变量", len(pool) >= 30, f"池子 {len(pool)} 个", kind="control")
bias_q = abs(float(np.median(q_draws)) - float(np.quantile(rs, .95)))
G.negative_control("阴性:跨池抽样应围绕全池值(无系统偏移)", null=bias_q,
                   effect=wq, null_spread=float(np.std(q_draws)),
                   null_kind="全池值与跨池抽样中位之差")
big = (wq / W_IN_Q >= 1.5) or (wm / W_IN_M >= 1.5)
print("\n" + "=" * 70)
if len(pool) >= 30:
    verdict = ("**W-POOL-BIG:跨池摆动明显更大 -> 所有参照门槛都要带上「参照从哪个池子挑」**"
               if big else "**W-POOL-OK:跨池摆动不比池内大 -> 选择不确定性确实小**")
    print(f"控制齐备 ⇒ 评判。q95 比 {wq/W_IN_Q:.2f} · 中位比 {wm/W_IN_M:.2f} -> {verdict}")
    print("⚠ 通过的 KILL 会怎样失败:**关键词表是我写的** —— 「与性无关」的边界仍带我的判断,"
          f"而且池子被 `POOL_CAP={POOL_CAP}` 截断,**不是全问卷**。两者都会低估跨池摆动。")
else:
    verdict = f"UNVERIFIED —— 池子只有 {len(pool)} 个"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(pool_size=len(pool), pool_head=pool[:20], pool_r=rs.tolist(),
               w_cross_q95=wq, w_cross_med=wm, w_in_q95=W_IN_Q, w_in_med=W_IN_M,
               ratio_q95=wq/W_IN_Q, ratio_med=wm/W_IN_M, pool_cap=POOL_CAP,
               verdict=verdict, unchallenged=True),
          open(OUT / "candidate_pool.json", "w"), indent=1)
print(f"\nwrote {OUT/'candidate_pool.json'}")
