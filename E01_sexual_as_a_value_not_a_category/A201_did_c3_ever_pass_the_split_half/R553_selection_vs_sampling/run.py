"""E01·A201·R553 — 门槛的不确定性,有多少来自抽样,有多少来自「选了哪些变量」

`#508` 的 NEXT,**它指向我不愿意看到的结果**:若「选择」的摆动明显大于「抽样」的,
则 `#485a` 以来所有用参照分布定门槛的轮次(**包括我上一轮刚给出的那三行**)都要带上这一条。

G1 ESTIMAND(先于方法):在 `R529` 的 **12 个参照 |r|** 上,比较两种摆动的**宽度**:
  **抽样**:对同一批 12 个值**有放回**重抽(k 不变)-> q95 与中位的 95% 区间宽度;
  **选择**:从 12 个里**无放回**取 k' 个(k' = 11 与 6)-> 同样两个量的 95% 区间宽度。
⛔ **STRONGEST CONFOUND,写在跑之前:留一同时改变 k**,而 q95 的行为本身依赖 k。
  控制:**固定 k' 后,把「无放回子集」与「有放回重抽」在同一个 k' 上比** ——
  这样 k 被扣住,剩下的差别才是「选择」。

WORLDS:
  W-SAMPLING 选择的摆动 ≤ 抽样的 -> 参照构成不是主要风险,RULE-v3 可放心用
  W-SELECTION 选择的摆动明显更大 -> **门槛的真正不确定性一直被低估**,`#508a` 的三行要加注
  | World       | now | 选择≤抽样 | 选择明显更大 |
  | W-SAMPLING  | 0.4 | 0.85      | 0.10 |
  | W-SELECTION | 0.6 | 0.10      | 0.85 |
判据(预注册):同一 k' 上,**选择区间宽度 / 抽样区间宽度 ≥ 1.5** 记为「明显更大」。

CONTROLS:正对照 = **两种摆动都必须随 k' 增大而变窄**(否则这个仪器看不见 k 的作用);
  阴性 = 在**同一批值**上,两种做法在 k'=12 时应当**几乎相同**(无放回取 12 = 原样)。
IMPOSSIBLE:只有一个参照集(12 个) ⇒ **「选择」的总体本身就是这 12 个**,
  真正的选择不确定性(本该从更大的候选池里抽)**结构上量不到**,本轮只能量它的下界。
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
r529 = json.load(open(ROOT / "E02_condemnation_is_not_rarity/"
                      "A192_does_condemnation_suppress_the_report_or_the_act/"
                      "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))
ref = np.array([abs(x["r"]) for x in r529["reference"]])
print(f"R529 参照 k={len(ref)}  值={np.round(np.sort(ref),4).tolist()}")
rng = np.random.default_rng(20260805)

def width(vals):
    lo, hi = np.quantile(vals, [.025, .975]); return float(hi - lo)

rows = []
for kp in (11, 9, 6, 12):
    samp_q = [np.quantile(rng.choice(ref, kp, replace=True), .95) for _ in range(4000)]
    samp_m = [np.median(rng.choice(ref, kp, replace=True)) for _ in range(4000)]
    sel_q = [np.quantile(rng.choice(ref, kp, replace=False), .95) for _ in range(4000)]
    sel_m = [np.median(rng.choice(ref, kp, replace=False)) for _ in range(4000)]
    rows.append(dict(k=kp, samp_q95_w=width(samp_q), sel_q95_w=width(sel_q),
                     samp_med_w=width(samp_m), sel_med_w=width(sel_m),
                     ratio_q95=width(sel_q) / max(width(samp_q), 1e-12),
                     ratio_med=width(sel_m) / max(width(samp_m), 1e-12)))
    print(f"  k'={kp:2d}  q95 宽度 抽样={width(samp_q):.4f} 选择={width(sel_q):.4f} "
          f"比={rows[-1]['ratio_q95']:.2f}   |   中位 宽度 抽样={width(samp_m):.4f} "
          f"选择={width(sel_m):.4f} 比={rows[-1]['ratio_med']:.2f}")

G = Gate("门槛的不确定性,抽样 vs 选择")
w11 = [r for r in rows if r["k"] == 11][0]; w6 = [r for r in rows if r["k"] == 6][0]
narrows = bool(w11["samp_q95_w"] < w6["samp_q95_w"] and w11["sel_q95_w"] < w6["sel_q95_w"])
G.asserted("正对照:两种摆动都必须随 k' 增大而变窄", narrows,
           f"k'=11 {w11['samp_q95_w']:.4f}/{w11['sel_q95_w']:.4f} vs "
           f"k'=6 {w6['samp_q95_w']:.4f}/{w6['sel_q95_w']:.4f}", kind="control")
w12 = [r for r in rows if r["k"] == 12][0]
G.asserted("阴性:k'=12 时无放回取满 -> 选择宽度应为 0", w12["sel_q95_w"] < 1e-9,
           f"k'=12 选择宽度 = {w12['sel_q95_w']:.6f}", kind="control")

main = [r for r in rows if r["k"] in (11, 9, 6)]
big_q = sum(1 for r in main if r["ratio_q95"] >= 1.5)
big_m = sum(1 for r in main if r["ratio_med"] >= 1.5)
print(f"\n比值 ≥1.5 的 k':q95 {big_q}/3 · 中位 {big_m}/3")
print("\n" + "=" * 70)
if narrows and w12["sel_q95_w"] < 1e-9:
    verdict = ("**W-SELECTION:选择的摆动明显更大 -> `#508a` 的三行与所有参照门槛都要加注**"
               if big_q >= 2 else
               "**W-SAMPLING:选择的摆动不比抽样大 -> 参照构成不是主要风险**")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:**「选择」的总体在这里就是这 12 个变量本身** ——"
          "真正该量的是「从更大的候选池里会选出哪 12 个」,而那个池子没有被列过。"
          "⇒ **本轮量到的是选择不确定性的下界,不是它本身。**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(变窄={narrows} k12零宽={w12['sel_q95_w']:.2e})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(ref=ref.tolist(), rows=rows, big_q95=big_q, big_med=big_m,
               verdict=verdict, unchallenged=True),
          open(OUT / "selection_vs_sampling.json", "w"), indent=1)
print(f"\nwrote {OUT/'selection_vs_sampling.json'}")
