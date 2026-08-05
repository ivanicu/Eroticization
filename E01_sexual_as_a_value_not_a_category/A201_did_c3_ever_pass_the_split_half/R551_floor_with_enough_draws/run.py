"""E01·A201·R551 — 把那个 3 抽的底加到够,并回头核对所有「小 k 的 q95」

`#506` 的 NEXT。**行动类型:CLOSURE**(方法债,不开新世界)。
`#506c` 自认:打乱底只由 **3 抽**估得 `[0.781, 0.500, 0.844]`,q95≈0.838,
而 `c3` 是 0.875 —— **margin 0.037,且 3 个点的 q95 实质上就是最大值**
(`#490c`/`#491a` 刚教过,我又用了一次)。

G1 ESTIMAND:同一个符号一致率统计量,**打乱底抽数 3 -> 60**,给出真正的分位;
  并同时把 **`c3` 与 `c1` 的劈半也各加到 20 次**,让两边都有区间而不是三个点。
⚠ §P7 成本:单次 `loadings()` 要重算 32×32 相关矩阵(每格一次 `corrcoef`)。
  已先量:载入后单次 profile ≈ 0s,矩阵本身是 32²=1024 次相关 -> 单抽秒级。
  **60 + 40 抽在预算内**,不需要 pueue。

预注册的判据(写在跑之前):
  `c3` 的符号一致率**中位** vs 打乱底的 **q95(60 抽)**:
    高于 -> `#506c` 的「高于底」成立,可从「未建立」升级;
    不高于 -> **「高于底」正式撤回**;⚠ 而 `#506b` 的「三次相同」**独立成立,不随之倒**
    (这一点已写进页面措辞,所以撤回时**不必改页面**)。

WORLDS:
  W-ABOVE 高于底 -> 多一条独立支持
  W-BELOW 不高于 -> **撤回「高于底」**,页面不动
  | World   | now | 高于 | 不高于 |
  | W-ABOVE | 0.5 | 0.85 | 0.10 |
  | W-BELOW | 0.5 | 0.10 | 0.85 |

⛔ STRONGEST CONFOUND:符号一致率有 **0.5 的地板**,且 32 个符号 -> 取值只能是 k/32,
  **分辨率本身是 1/32 = 0.031** —— 与 `#506c` 的 margin 0.037 同量级。
  ⇒ **即使加抽,这个比较也只能分辨约一格。** 先写下,不事后解释。
CONTROLS:正对照 `c1` 必须比 `c3` 稳;阴性 = 打乱底本身(它就是零臂)。
KILL:if 正对照通过 -> 按上表读;else UNVERIFIED。
IMPOSSIBLE:统计量的取值粒度 1/32 ⇒ 精度有硬下限 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, warnings, time
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
_SRC = (ROOT / 'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/'
        'R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""', 2)[2].split('def fit_apply')[0])
ALLR = np.flatnonzero(ok)
print(f"NB={NB} NN={NN} 可用={len(ALLR)};统计量粒度 = 1/{NB} = {1/NB:.4f}")

def loadings(rows, k):
    m = np.zeros(NN, bool); m[rows] = True
    def prof_(X):
        F = np.isfinite(X); Z = np.where(F, X, 0.0); tot = Z.sum(0); ct = F.sum(0)
        R = np.full_like(X, np.nan)
        for b in range(NB):
            lo = np.where(ct - F[b] >= 6, (tot - Z[b]) / np.maximum(ct - F[b], 1), np.nan)
            R[b] = np.where(F[b], X[b] - lo, np.nan); R[b] = R[b] - np.nanmean(np.where(m, R[b], np.nan))
        return R
    Ra, Rb = prof_(A), prof_(B)
    C = np.zeros((NB, NB))
    for i in range(NB):
        for j in range(NB):
            mm = np.isfinite(Ra[i]) & np.isfinite(Rb[j]) & m
            if mm.sum() > 200: C[i, j] = np.corrcoef(Ra[i][mm], Rb[j][mm])[0, 1]
    C = (C + C.T) / 2; w, V = np.linalg.eigh(C)
    return V[:, np.argsort(-w)[k]]

def one(seed, k, shuffle=False):
    rng = np.random.default_rng(seed)
    r = rng.permutation(NN)[:len(ALLR)] if shuffle else ALLR.copy()
    p = rng.permutation(r); h = len(p) // 2
    v1, v2 = loadings(p[:h], k), loadings(p[h:], k)
    s = np.sign(v1) * np.sign(v2)
    return float(max(np.mean(s > 0), np.mean(s < 0)))

t0 = time.time(); _ = one(1, 2); dt = time.time() - t0
print(f"单抽耗时 = {dt:.2f}s -> 60+40 抽预计 {100*dt/60:.1f} 分钟(§P7 先量再定)")

NB_NULL, NB_REAL = 60, 20
null = np.array([one(10000 + i, 2, shuffle=True) for i in range(NB_NULL)])
c3 = np.array([one(2000 + i, 2) for i in range(NB_REAL)])
c1 = np.array([one(3000 + i, 0) for i in range(NB_REAL)])
q95 = float(np.quantile(null, .95))
print(f"\n打乱底({NB_NULL} 抽):中位 {np.median(null):.4f}  q95 {q95:.4f}  "
      f"范围 [{null.min():.3f}, {null.max():.3f}]")
print(f"  ⚠ 对比 `#506c` 的 3 抽 q95 = 0.838 —— 差 {abs(q95-0.838):.4f}")
print(f"c3({NB_REAL} 抽):中位 {np.median(c3):.4f}  范围 [{c3.min():.3f}, {c3.max():.3f}]")
print(f"c1({NB_REAL} 抽):中位 {np.median(c1):.4f}  范围 [{c1.min():.3f}, {c1.max():.3f}]")

G = Gate("把 3 抽的底加到够")
instr = float(np.median(c1)) > float(np.median(c3))
G.asserted("正对照:c1 必须比 c3 稳", instr,
           f"c1 {np.median(c1):.4f} vs c3 {np.median(c3):.4f}", kind="control")
above = float(np.median(c3)) > q95
print("\n" + "=" * 70)
if instr:
    verdict = ("**W-ABOVE:`c3` 中位高于 60 抽的底 q95 -> 「高于底」成立**" if above else
               "**W-BELOW:「高于底」正式撤回**(而 `#506b` 的「三次相同」独立成立,页面不动)")
    print(f"控制齐备 ⇒ 评判。c3 中位 {np.median(c3):.4f} vs 底 q95 {q95:.4f} -> {verdict}")
    print(f"⚠ 通过的 KILL 会怎样失败:统计量的取值只能是 k/32,**分辨率 {1/NB:.4f}**,"
          f"而这里的 margin 是 {abs(np.median(c3)-q95):.4f} —— **同量级,只能分辨约一格**。")
else:
    verdict = "UNVERIFIED —— 正对照失败"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(null=null.tolist(), c3=c3.tolist(), c1=c1.tolist(),
               null_q95=q95, null_q95_3draw_old=0.838, c3_med=float(np.median(c3)),
               c1_med=float(np.median(c1)), above=bool(above), granularity=1/NB,
               verdict=verdict, n_null=NB_NULL, n_real=NB_REAL, unchallenged=True),
          open(OUT / "floor_with_enough_draws.json", "w"), indent=1)
print(f"\nwrote {OUT/'floor_with_enough_draws.json'}")
