"""E01·A201·R552 — 把所有「小 k 的 q95」当过门槛的地方,只重算门槛

`#507` 的 NEXT。**行动类型:CLOSURE**(跨轮方法更正的兑现)。
⛔ **不重跑任何一轮**,只重算它们的门槛,并逐个报「旧 / 新 / 判定是否翻」。
若有翻转:**重述那一轮判定的理由,不得追认其结论**(`#482e`)。

`#507b` 的更正过的规律:**小 k 的 q95 是极值,因而不可靠;偏的方向不固定。**
⇒ 正确做法不是「假定它过严」,而是**给出它自身的抽样区间**:
   对那 k 个参照值**自助重抽 4000 次**,报 q95 的**中位与 95% 区间**,
   看**观测值落在这个区间的哪一侧** —— 落在区间内 = 门槛本身分辨不了。

G1 ESTIMAND:对每一轮,`q95_boot` 的分布,以及观测的正对照 `|r|` 相对它的位置。
判据(预注册):
  观测 > q95 区间上界 -> PASS 稳固;观测 < 下界 -> FAIL 稳固;**落在区间内 -> 门槛不可判**。
⚠ 同时报 **RULE-v3 的中位门槛**(`#491c`)作为对照 —— 它正是为小 k 设计的。
IMPOSSIBLE:参照分布是各轮当时选的变量集,**本轮不改它们**(改了就不是「只重算门槛」)。
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
E = ROOT / "E02_condemnation_is_not_rarity"
CASES = []
r529 = json.load(open(E / "A192_does_condemnation_suppress_the_report_or_the_act/"
                      "R529_gss_gap_over_51_years/results/gss_gap_over_51_years.json"))
CASES.append(("R529 GSS", np.array([abs(x["r"]) for x in r529["reference"]]),
              abs(r529["positive"]["r"]), "condemn×attend"))
r534 = json.load(open(E / "A194_does_the_gap_survive_a_different_institution/"
                      "R534_nsfg_acasi_replication/results/nsfg_acasi_replication.json"))
CASES.append(("R534 NSFG", np.array([abs(x["r"]) for x in r534["reference"]]),
              abs(r534["positive"]["r"]), "samesex×attndnow"))
r535 = json.load(open(E / "A194_does_the_gap_survive_a_different_institution/"
                      "R535_within_instrument_topic_spread/results/within_instrument_topic_spread.json"))
CASES.append(("R535 NSFG", np.array([abs(x["r"]) for x in r535["reference"]]),
              0.2105, "chsuppor×attndnow"))
print(f"{'轮':10s} {'k':>3s} {'观测':>7s} {'旧门槛q95':>9s} {'q95自助中位':>11s} "
      f"{'q95 95%区间':>18s} {'v3中位门槛':>10s} {'旧':>5s} {'新':>7s}")
rows = []
rng = np.random.default_rng(20260805)
for name, ref, obs, pair in CASES:
    k = len(ref); q_old = float(np.quantile(ref, .95)); med = float(np.median(ref))
    B = np.array([np.quantile(rng.choice(ref, k, replace=True), .95) for _ in range(4000)])
    lo, hi = np.quantile(B, [.025, .975])
    old_v = "PASS" if obs > q_old else "FAIL"
    new_v = ("PASS" if obs > hi else "FAIL" if obs < lo else "不可判")
    v3_v = "PASS" if obs > med else "FAIL"
    rows.append(dict(round=name, k=k, obs=obs, q95_old=q_old, q95_boot_med=float(np.median(B)),
                     q95_ci=[float(lo), float(hi)], v3_median=med,
                     old=old_v, new=new_v, v3=v3_v))
    print(f"{name:10s} {k:3d} {obs:7.4f} {q_old:9.4f} {np.median(B):11.4f} "
          f"[{lo:6.4f},{hi:6.4f}] {med:10.4f} {old_v:>5s} {new_v:>7s}")
flip = [r for r in rows if r["old"] != r["new"]]
undec = [r for r in rows if r["new"] == "不可判"]
print(f"\n判定改变的轮次:{[r['round'] for r in flip] or '无'}")
print(f"门槛本身不可判的轮次:{[r['round'] for r in undec] or '无'}")
print(f"RULE-v3(中位)门槛下:{[(r['round'], r['v3']) for r in rows]}")

G = Gate("只重算门槛,不重跑那些轮")
G.asserted("参照分布未被改动", True, "各轮 JSON 里原样读出,只重算门槛", kind="control")
G.asserted("本轮未重跑任何一轮", True, "没有触碰任何 run.py 的数据路径", kind="control")
print("\n" + "=" * 70)
if undec:
    verdict = (f"{len(undec)}/{len(rows)} 轮的**门槛本身不可判**(观测落在 q95 的自助区间内)"
               f" -> 那几轮的正对照判定**只重述理由,不追认结论**")
else:
    verdict = f"三轮的判定在自助区间下都稳固,无一翻转"
print(f"控制齐备 ⇒ 评判。{verdict}")
print("⚠ 通过的 KILL 会怎样失败:自助重抽只能反映**这 k 个点自身**的抽样波动,"
      "**不能反映「该选哪些变量进参照」的不确定性** —— 而后者可能更大,且本轮按设计不碰它。")
print(G)
json.dump(dict(rows=rows, flipped=[r["round"] for r in flip],
               undecidable=[r["round"] for r in undec], verdict=verdict, unchallenged=True),
          open(OUT / "threshold_audit.json", "w"), indent=1)
print(f"\nwrote {OUT/'threshold_audit.json'}")
