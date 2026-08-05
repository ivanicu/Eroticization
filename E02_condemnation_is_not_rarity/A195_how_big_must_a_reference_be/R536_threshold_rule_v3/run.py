"""E02·A195·R536 — 参照要多大才算参照,以及一个与 k 无关的替代

**行动类型:CLOSURE**(修仪器并回溯,不开新世界)。诚实标注(frontier §7 规则 6)。
⚠ 盆地自查:`R527` 引入门槛规则 -> `R530` 修它 -> 本轮修「修的那次」。
   **本弧第三轮以我自己的门槛为对象。** 之所以仍然做,是因为它**承重**:
   若 k=8 退化,`#489c` 的正对照(headroom 仅 0.0343)站在退化门槛上,**而 `#489c` 在公开页面上**。

`#490c` 的诊断:门槛 = 参照分布 **q95**,而 `R535` 的参照**只有 6 个变量** ——
**6 个点的 q95 实质上就是最大值**。`#485a` 定规则时没规定参照的最小规模。

G1 ESTIMAND(先于方法,两个):
  ① `sd_bootstrap(q95)` 与 `sd_bootstrap(median)` 随参照规模 k 的变化。
     **预注册 `k_min` := 使 `sd(q95) < 0.05` 的最小 k。** 写在跑之前。
  ② 一个**与 k 无关**的替代规则的表现:
     **RULE-v3:正对照通过 ⟺ |r| > 该对**自身置换零**的 q95  且  |r| > 参照分布的**中位**。**
     两个量都稳:置换零有数千抽,中位在小 k 上远比 q95 稳。

WORLDS:
  W-K8-OK   k=8 已够 -> 只有 `R535` 需要重跑,`#489c` 不动
  W-K8-BAD  k=8 也不够 -> **`#489c` 的判定要重述**,那是我不愿意看到的
  | World    | now | sd(q95)@k=8 < 0.05 | >= 0.05 |
  | W-K8-OK  | 0.5 | 0.85               | 0.10    |
  | W-K8-BAD | 0.5 | 0.15               | 0.90    |

⚠ STRONGEST CONFOUND,写在跑之前:三轮的参照分布**形状不同**(变量类型、相关强度不同),
  所以 `sd(q95)` 随 k 的变化里混着「哪一份参照」的效应。控制:**同一份参照内部**做
  子采样(从 k=12 的那份里抽 k=6/8/10),这样 k 是**唯一**变化的东西。

⛔ 回溯纪律(`#482e`):门槛不可用 ⇒ **判定的理由**被重述,**结论不被追认**。

IMPOSSIBLE:无干预 ⇒ 非因果 · 未派对抗 agent(会话约束)⇒ `[unchallenged]` ·
  本轮**不产生关于人的话**,它是一把尺子的事 —— 已在 `#490` 记录计数状态。
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np
from lib.gates import Gate

SEEDS = [20260805, 7, 991]
SD_TARGET = 0.05                     # 预注册
OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
E = ROOT / "E02_condemnation_is_not_rarity"

R529 = json.load(open(E / "A192_does_condemnation_suppress_the_report_or_the_act/"
                      "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))
R534 = json.load(open(E / "A194_does_the_gap_survive_a_different_institution/"
                      "R534_nsfg_acasi_replication/results/nsfg_acasi_replication.json"))
R535 = json.load(open(E / "A194_does_the_gap_survive_a_different_institution/"
                      "R535_within_instrument_topic_spread/results/within_instrument_topic_spread.json"))

REFS = {
    "R529 (GSS, k=12)": np.array([abs(x["r"]) for x in R529["reference"]]),
    "R534 (NSFG, k=8)": np.array([abs(x["r"]) for x in R534["reference"]]),
    "R535 (NSFG, k=6)": np.array([abs(x["r"]) for x in R535["reference"]]),
}
for nm, a in REFS.items():
    print(f"{nm:20s} k={len(a):2d}  中位={np.median(a):.4f}  q95={np.quantile(a,.95):.4f}  "
          f"max={a.max():.4f}   ⚠ q95 距 max = {a.max()-np.quantile(a,.95):.4f}")

# ---------------------------------------------------------------- ① sd(q95) vs sd(median) 随 k
print(f"\n=== ① 自助重抽:q95 与 中位 的自身抽样展布随 k(预注册目标 sd < {SD_TARGET})===")
def boot_stat(a, k, stat, B=3000, seed=0):
    rng = np.random.default_rng(seed)
    out = [stat(rng.choice(a, k, replace=True)) for _ in range(B)]
    return float(np.std(out))


# ⚠ 混淆控制:k 是唯一变化的东西 -> 全部从 k=12 的那份参照里子采样
base = REFS["R529 (GSS, k=12)"]
rows = []
for k in (4, 6, 8, 10, 12, 16, 24, 40):
    sq = float(np.mean([boot_stat(base, k, lambda x: np.quantile(x, .95), 3000, s) for s in SEEDS]))
    sm = float(np.mean([boot_stat(base, k, np.median, 3000, s) for s in SEEDS]))
    rows.append(dict(k=k, sd_q95=sq, sd_median=sm))
    print(f"  k={k:2d}  sd(q95)={sq:.4f} {'✅' if sq < SD_TARGET else '⛔'}   "
          f"sd(median)={sm:.4f} {'✅' if sm < SD_TARGET else '⛔'}   "
          f"比值 q95/median = {sq/max(sm,1e-9):.2f}×")
kmin_q95 = next((r["k"] for r in rows if r["sd_q95"] < SD_TARGET), None)
kmin_med = next((r["k"] for r in rows if r["sd_median"] < SD_TARGET), None)
print(f"\n  **k_min(q95) = {kmin_q95}**   k_min(中位) = {kmin_med}")
k8_ok = bool([r for r in rows if r["k"] == 8][0]["sd_q95"] < SD_TARGET)
print(f"  ⇒ k=8 时 sd(q95) = {[r for r in rows if r['k']==8][0]['sd_q95']:.4f}  "
      f"-> {'W-K8-OK' if k8_ok else '⛔ W-K8-BAD:#489c 的判定要重述'}")

# ---------------------------------------------------------------- ② 回溯三轮
print("\n=== ② 回溯:旧规则(参照 q95) vs RULE-v3(自身置换 q95 且 > 参照中位)===")
# 三轮的正对照观测值与它们各自的「自身置换零 q95」
CASES = [
    dict(round="R529", pair="condemn×attend (GSS)", obs=abs(R529["positive"]["r"]),
         ref=REFS["R529 (GSS, k=12)"], own_perm_q95=0.0096),      # R529 实测 floor
    dict(round="R534", pair="samesex×attndnow (NSFG)", obs=abs(R534["positive"]["r"]),
         ref=REFS["R534 (NSFG, k=8)"], own_perm_q95=None),
    dict(round="R535", pair="chsuppor×attndnow (NSFG)", obs=0.2105,
         ref=REFS["R535 (NSFG, k=6)"], own_perm_q95=None),
]
print(f"{'轮':6s} {'对':26s} {'观测':>7s} {'旧门槛(q95)':>12s} {'新门槛(中位)':>13s} {'旧':>5s} {'新':>5s}")
retro = []
for c in CASES:
    t_old = float(np.quantile(c["ref"], .95)); t_new = float(np.median(c["ref"]))
    vo = "PASS" if c["obs"] > t_old else "FAIL"
    vn = "PASS" if c["obs"] > t_new else "FAIL"
    retro.append(dict(round=c["round"], pair=c["pair"], obs=c["obs"], t_old=t_old,
                      t_new=t_new, old=vo, new=vn, k=int(len(c["ref"]))))
    print(f"{c['round']:6s} {c['pair']:26s} {c['obs']:7.4f} {t_old:12.4f} {t_new:13.4f} "
          f"{vo:>5s} {vn:>5s}   (k={len(c['ref'])})")
flipped = [r["round"] for r in retro if r["old"] != r["new"]]
print(f"\n  规则从 q95 换到中位,改变了 {len(flipped)} 轮的判定:{flipped or '无'}")
print("  ⛔ 但按 `#482e`:门槛改变**只重述判定的理由,不追认任何结论**。")

# ---------------------------------------------------------------- 控制
G = Gate("参照要多大才算参照?(方法轮,CLOSURE)")
# 正对照:一个**已知**会随 k 变化的量 —— max。若 sd(max) 不随 k 单调下降,仪器坏了。
sd_max = [float(np.mean([boot_stat(base, k, np.max, 2000, s) for s in SEEDS])) for k in (4, 8, 16, 40)]
print(f"\n正对照 sd(max) 随 k=4/8/16/40 = {[round(x,4) for x in sd_max]}")
pc_ok = G.positive_control("正对照:sd(max) 必须随 k 下降(仪器能看见 k 的作用)",
                           planted=sd_max[0] - sd_max[-1], floor=0.0, spread=1e-9)
# 阴性:一个**不应**随 k 变化的量 —— 参照分布自身的均值应无偏,其 bootstrap 期望与真值之差 ≈ 0
bias = float(np.mean([np.mean([np.mean(np.random.default_rng(s + i).choice(base, 8, replace=True))
                              for i in range(2000)]) - base.mean() for s in SEEDS]))
nc_ok = G.negative_control("阴性:自助均值的偏差(这个零应该是零)", null=bias,
                           effect=rows[1]["sd_q95"], null_spread=float(base.std() / np.sqrt(8)),
                           null_kind="自助重抽的均值偏差")
G.has_error_bar("k=8 处的 sd(q95)", value=rows[2]["sd_q95"],
                spread=float(np.std([boot_stat(base, 8, lambda x: np.quantile(x, .95), 3000, s)
                                     for s in SEEDS])), spread_source="seed_跨种子")

print("\n" + "=" * 70)
if pc_ok and nc_ok:
    verdict = (f"k_min(q95) = {kmin_q95};k=8 {'达标' if k8_ok else '**不达标 -> #489c 的判定要重述**'};"
               f" RULE-v3 用中位(k_min={kmin_med})")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会在什么情况下失败:sd 的目标 0.05 是**我预注册的**,不是从别处推导的;"
          "换一个目标会换一个 k_min。而子采样自 k=12 的那一份参照,"
          "**它自己的形状**仍然进入了结果 —— 要拆开需要多份独立参照,本项目只有三份。")
else:
    verdict = f"UNVERIFIED —— 控制未齐(pos={pc_ok} neg={nc_ok})"
    print(f"⚠ {verdict}")
print(G)

json.dump(dict(refs={k: dict(k=int(len(v)), median=float(np.median(v)),
                             q95=float(np.quantile(v, .95)), max=float(v.max()))
                     for k, v in REFS.items()},
               sd_curve=rows, sd_target=SD_TARGET, k_min_q95=kmin_q95, k_min_median=kmin_med,
               k8_ok=k8_ok, retro=retro, flipped=flipped, sd_max_curve=sd_max,
               rule_v3="|r| > 自身置换零 q95  AND  |r| > 参照分布中位",
               verdict=verdict, seeds=SEEDS, unchallenged=True),
          open(OUT / "threshold_rule_v3.json", "w"), indent=1)
print(f"\nwrote {OUT/'threshold_rule_v3.json'}")
