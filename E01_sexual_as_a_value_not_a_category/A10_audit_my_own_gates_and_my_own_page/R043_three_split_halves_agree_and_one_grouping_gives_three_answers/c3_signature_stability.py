"""E01·A201·R550 — `c3` 有没有通过过「三次劈半同一个答案」这一关

`#505` 的 NEXT。**行动类型:CLOSURE**(对已发表声明的稳定性审计,不开新世界)。
⚠ 这**不是**规则 ③ 禁止的「第四个候选」—— 它不为那 42% 提出任何解释,
   只问一个**已经在页面上的方向**有没有通过一个**我从未对它用过**的判据。

`#505c` 造出的判据:**三次随机劈半是否给出同一个答案。**
`#303` 已报过 `c3` 的可复现 `|cos| = 0.796`(`c1` 0.957),但**从未报过符号签名的稳定性**,
而 `#505` 的教训正是:**一个 |cos| 不低的方向,其符号模式仍可能每次都不一样。**

G1 ESTIMAND(先于方法):
  `c3` = 32×32 块相关矩阵的**第三个特征向量**(与 `R372`/`#303` 同一构造,直接复用其代码)。
  ① 三次随机劈半,每次算 `H1` 与 `H2` 的载荷,记 **|cos|** 与 **32 位符号签名的一致率**;
  ② **正对照 `c1`**(第一个特征向量)同样跑 —— 它必须明显更稳,否则这个仪器本身不能分辨稳与不稳;
  ③ 阴性/底:把**人**打乱后重跑(破坏人层结构),给出符号一致率的地板。

WORLDS:
  W-STABLE `c3` 是一个稳定方向 -> 符号一致率显著高于打乱底,且三种子彼此接近
  W-SHAKY  它一直不稳          -> 一致率接近打乱底 -> **页面上关于 `c3` 的每一句都要带这一条**
  | World    | now | 高于底 | 接近底 |
  | W-STABLE | 0.5 | 0.85   | 0.10   |
  | W-SHAKY  | 0.5 | 0.10   | 0.85   |

⛔ STRONGEST CONFOUND,写在跑之前:符号一致率**有 0.5 的地板**(`#502c` 刚踩过)——
  所以**不可用「零 < 0.5×效应」判**,必须与**打乱底**直接比分位。已写死。
KILL(条件式,预注册):
  if `c1` 的一致率 > `c3` 的(仪器能分辨):
      `c3` 一致率 > 打乱底 q95 -> W-STABLE
      否则                     -> **W-SHAKY,页面加注**
  else: UNVERIFIED(仪器分辨不了稳与不稳)
IMPOSSIBLE:32 维上的第三个特征向量本就脆弱 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
_SRC = (ROOT / 'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/'
        'R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""', 2)[2].split('def fit_apply')[0])
print(f"块数 NB={NB}  人数 NN={NN}  可用 {int(ok.sum())}")

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

ALLR = np.flatnonzero(ok)

def one_split(seed, k, shuffle=False):
    rng = np.random.default_rng(seed)
    r = ALLR.copy()
    if shuffle:
        r = rng.permutation(NN)[:len(ALLR)]      # 破坏「谁是有效样本」的结构
    p = rng.permutation(r); h = len(p) // 2
    v1, v2 = loadings(p[:h], k), loadings(p[h:], k)
    c = float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    s = np.sign(v1) * np.sign(v2)
    return abs(c), float(max(np.mean(s > 0), np.mean(s < 0)))

print("\n=== 三次劈半:|cos| 与 32 位符号签名一致率 ===")
res = {}
for nm, k in [("c1(正对照)", 0), ("c3", 2)]:
    rows = [one_split(s, k) for s in SEEDS]
    res[nm] = dict(cos=[r[0] for r in rows], sign=[r[1] for r in rows],
                   cos_med=float(np.median([r[0] for r in rows])),
                   sign_med=float(np.median([r[1] for r in rows])))
    print(f"  {nm:12s} |cos| = {[round(r[0],3) for r in rows]}  中位 {res[nm]['cos_med']:.3f}")
    print(f"  {'':12s} 符号一致 = {[round(r[1],3) for r in rows]}  中位 {res[nm]['sign_med']:.3f}")

null = [one_split(s + 555, 2, shuffle=True) for s in SEEDS]
null_sign = [r[1] for r in null]
nq95 = float(np.quantile(null_sign, .95)) if len(null_sign) > 2 else float(max(null_sign))
print(f"\n打乱底(c3,破坏有效样本结构):符号一致 = {[round(x,3) for x in null_sign]}  q95≈{nq95:.3f}")

G = Gate("`c3` 有没有通过过三次劈半这一关?")
instr_ok = res["c1(正对照)"]["sign_med"] > res["c3"]["sign_med"]
G.asserted("正对照:c1 必须比 c3 稳(否则仪器分辨不了)", instr_ok,
           f"c1 {res['c1(正对照)']['sign_med']:.3f} vs c3 {res['c3']['sign_med']:.3f}", kind="control")
above = res["c3"]["sign_med"] > nq95
print("\n" + "=" * 70)
if instr_ok:
    verdict = ("**W-STABLE:`c3` 的符号模式高于打乱底**" if above else
               "**W-SHAKY:`c3` 的符号模式没有高过打乱底 -> 页面上每一句关于 `c3` 的话都要带这一条**")
    print(f"控制齐备 ⇒ 评判。c3 符号一致 {res['c3']['sign_med']:.3f} vs 打乱底 q95 {nq95:.3f} -> {verdict}")
    print("⚠ 通过的 KILL 会怎样失败:符号一致率有 0.5 的地板(`#502c`),"
          "所以只能与打乱底比分位,不能用「零 < 0.5×效应」;而三次劈半只有三个点,"
          "**中位数是三个数的中间那个**,不是一个有区间的估计。")
else:
    verdict = "UNVERIFIED —— 正对照失败:c1 不比 c3 稳,仪器分辨不了稳与不稳"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(res=res, null_sign=null_sign, null_q95=nq95, instrument_ok=bool(instr_ok),
               above=bool(above), verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "c3_signature_stability.json", "w"), indent=1)
print(f"\nwrote {OUT/'c3_signature_stability.json'}")
