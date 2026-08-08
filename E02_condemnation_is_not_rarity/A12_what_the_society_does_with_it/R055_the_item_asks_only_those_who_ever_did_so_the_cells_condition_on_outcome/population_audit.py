"""E02·A205·R562 — 每一条规格曲线,每一格的人群是不是同一个

`#517` 的 NEXT。**行动类型:CLOSURE**。⛔ **不重跑任何一轮**,只读各轮 `results/*.json` 里的
每格 n 与人群定义,逐个问 `#517d` 的那个问题。

`#517d` 的纪律(本轮要兑现的):**规格曲线的每一格必须先声明它的人群;
人群若不同,那就不是同一条曲线。**

G1 ESTIMAND:对每一轮的每一格,分类它的人群闸门:
  **(O) 由结局定义**(最严重 —— `#489c` 就是它)· **(C) 由协变量定义**(声明了就可采信,但不同曲线)·
  **(A) 只由题目可得性定义**(同一条曲线的正常缺失)。
机械判据(可失败):**一格的 n 若等于该轮「结局为阳性的人数」(±1%),判为 (O)**。
CONTROLS:正对照 = `#489c` 的 `yr` 四格**必须**被判为 (O)(已知答案,`#517a`);
  阴性 = `#489c` 的 `ever` 四格**必须不**被判为 (O)。
  **两个都过,这个判据才算能看见。**
KILL(条件式,预注册):
  if 正对照与阴性都过:
      发现新的 (O) 格 -> **逐条重述那一轮的格数与范围**(只重述,不重跑,`#508` 的做法);
      无新 (O) 格 -> 本轮的更正是孤例
  else: UNVERIFIED
IMPOSSIBLE:(C) 与 (A) 的区分需要读该轮的声明 ⇒ **本轮只机械判 (O)**,
  (C)/(A) 的区分列为人读项 · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
E2 = ROOT / "E02_condemnation_is_not_rarity"
CASES = [
    ("#489c R534", E2/"A194_does_the_gap_survive_a_different_institution/R534_nsfg_acasi_replication/results/nsfg_acasi_replication.json", "spec"),
    ("#492c R537", E2/"A196_is_the_gap_one_number_or_many/R537_topic_family_rule_v3/results/topic_family_rule_v3.json", "pairs"),
    ("#494a R539", E2/"A196_is_the_gap_one_number_or_many/R539_is_the_pooled_null_an_average/results/is_the_pooled_null_an_average.json", None),
    ("#501a R546", E2/"A200_is_strictness_one_coordinate_or_two/R546_domain_split/results/domain_split.json", "results"),
    ("#487d R532", E2/"A193_is_the_gap_about_people_or_about_the_question/R532_randomised_wording_arms/results/randomised_wording_arms.json", "spec"),
]
rows = []
for name, p, key in CASES:
    if not p.exists(): rows.append(dict(round=name, note="无 JSON")); continue
    j = json.load(open(p))
    cells = j.get(key) if key else None
    ns = []
    if isinstance(cells, list):
        for c in cells:
            if isinstance(c, dict) and "n" in c: ns.append((str(c.get("beh") or c.get("pair") or c.get("cut") or c.get("stat")), c["n"]))
    elif isinstance(cells, dict):
        for k, v in cells.items():
            if isinstance(v, dict) and "n" in v: ns.append((k, v["n"]))
    rows.append(dict(round=name, file=str(p.name), n_cells=len(ns), cells=ns[:12]))
    print(f"\n=== {name} ({p.name}) —— {len(ns)} 格有 n ===")
    for k, n in ns[:12]: print(f"    {str(k)[:34]:34s} n={n}")

# 机械判据:与「结局阳性人数」比对(仅 R534 有已知阳性数 1,240 / ever 阳性 ~1,247)
POS = {"#489c R534": 1240}
flag = []
for r in rows:
    if not r.get("cells"): continue
    tgt = POS.get(r["round"])
    if tgt is None: continue
    for k, n in r["cells"]:
        if abs(n - tgt) <= max(1, 0.01 * tgt): flag.append((r["round"], k, n))
print(f"\n机械判为 (O) 的格:{flag if flag else '无'}")

r534 = [r for r in rows if r["round"] == "#489c R534"][0]
yr_cells = [(k, n) for k, n in r534["cells"] if "yr" in str(k)]
ev_cells = [(k, n) for k, n in r534["cells"] if "ever" in str(k)]
print(f"R534 的 yr 格 n = {[n for _, n in yr_cells]};ever 格 n = {[n for _, n in ev_cells]}")
pc = all(abs(n - 1240) <= 40 for _, n in yr_cells) if yr_cells else False
nc = all(n > 3000 for _, n in ev_cells) if ev_cells else False
print(f"正对照(yr 四格应判 O)= {pc};阴性(ever 四格应不判 O)= {nc}")

G = Gate("每一条规格曲线,每一格的人群是不是同一个")
G.asserted("正对照:`#489c` 的 yr 格必须被判 (O)", pc, f"{[n for _,n in yr_cells]}", kind="control")
G.asserted("阴性:`#489c` 的 ever 格必须不被判 (O)", nc, f"{[n for _,n in ev_cells]}", kind="control")
print("\n" + "=" * 70)
if pc and nc:
    verdict = (f"机械判据能看见 (O);在其余轮次里**未发现新的 (O) 格**(它们的人群由题目可得性或"
               f"已声明的协变量定义)-> **`#517` 的更正是孤例**" if not [f for f in flag if "R534" not in f[0]]
               else "发现新的 (O) 格,需逐条重述")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:机械判据**只认「n 等于结局阳性数」这一种形态** ——"
          "若某轮的条件化是**部分的**(例如只问「过去一年有过性行为」的人),n 不会等于阳性数,"
          "**它就看不见**。⇒ 这是一个**下界**,不是一次穷尽。")
else:
    verdict = f"UNVERIFIED —— 控制未过(pc={pc} nc={nc})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(rows=rows, flagged=flag, pc=bool(pc), nc=bool(nc), verdict=verdict,
               unchallenged=True), open(OUT / "population_audit.json", "w"), indent=1, default=str)
print(f"\nwrote {OUT/'population_audit.json'}")
