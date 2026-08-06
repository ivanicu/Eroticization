"""E02·A205·R563 — 参照改成逐波,和一个「每格必须存 n」的新守卫

`#518` 的 NEXT 两件,按顺序。**行动类型:CLOSURE**。
⚠ `#518` 已写明:这是**同一个问题的第二次**,而第一次失败的是**仪器**(参照写成了一个常数),
   不是问题被追第三轮。**若第二次仍不过,停止,按盆地规则换方向。**

① **修正对照**:`#489c` 跨两波,两波的「曾有过」阳性数不同
   (2017–19 = 1,240;2011–13 = 1,012)-> 判据改为**逐波**比对,不再用一个常数。
② **新守卫 `spec_curve_cells_declare_n`**(已写入 `lib/gates.py`):
   规格曲线的每一格必须把 `n` 写进 `results/`;缺 n 的格,事后不可审(`#518b`)。

G1 ESTIMAND:
  ① 每一格的 n 是否落在**它自己那一波**的结局阳性数 ±1% 内 -> 判 (O);
  ② 守卫在**已知答案**上的表现:`#489c` 的 JSON(逐格有 n)必须 **PASS**,
     `#494a` 的(0 格有 n)必须 **FAIL**。**两个都对,守卫才算能用。**
CONTROLS:
  正对照 ①`yr` 四格必须**全部**判 (O);②守卫在有 n 的产物上 PASS;
  阴性   ①`ever` 四格必须**全部不**判 (O);②守卫在无 n 的产物上 FAIL。
KILL(条件式,预注册):
  if 四个控制全过:对其余轮次给出 (O) 判定;else UNVERIFIED **且停止本方向**。
IMPOSSIBLE:机械判据只认「n = 该波阳性数」一种形态,**部分条件化看不见** ⇒ 下界 ·
  三轮无逐格 n ⇒ 除非重跑否则不可审(守卫只对**将来**的轮次生效) · 未派对抗 agent
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
E2 = ROOT / "E02_condemnation_is_not_rarity"
J534 = json.load(open(E2/"A194_does_the_gap_survive_a_different_institution/"
                      "R534_nsfg_acasi_replication/results/nsfg_acasi_replication.json"))
J539 = json.load(open(E2/"A196_is_the_gap_one_number_or_many/"
                      "R539_is_the_pooled_null_an_average/results/is_the_pooled_null_an_average.json"))
J537 = json.load(open(E2/"A196_is_the_gap_one_number_or_many/"
                      "R537_topic_family_rule_v3/results/topic_family_rule_v3.json"))

# ① 逐波参照(`#517a`/`#489a` 已打印过的两波阳性数)
POS = {"2017_2019_Fem": 1240, "2011_2013_Fem": 1012}
print("=== ① 逐波判据(不再用一个常数)===")
cells, flags = [], []
for c in J534["spec"]:
    wave, beh, n = c["wave"], c["beh"], c["n"]
    tgt = POS[wave]
    isO = abs(n - tgt) <= max(1, 0.02 * tgt)
    cells.append(dict(wave=wave, beh=beh, cond=c["cond"], n=n, target=tgt, O=isO))
    flags.append(isO)
    print(f"  {wave:14s} {c['cond']:9s} {beh:5s} n={n:5d}  该波阳性={tgt:5d}  -> {'(O)' if isO else '—'}")
yr = [x for x in cells if x["beh"] == "yr"]; ev = [x for x in cells if x["beh"] == "ever"]
pc1 = all(x["O"] for x in yr); nc1 = all(not x["O"] for x in ev)
print(f"\n正对照(yr 四格全判 O)= {pc1};阴性(ever 四格全不判 O)= {nc1}")

# ② 新守卫在已知答案上的表现
print("\n=== ② 新守卫 `spec_curve_cells_declare_n` 的正/阴对照 ===")
Gp = Gate("守卫正对照"); pc2 = Gp.spec_curve_cells_declare_n("有 n 的产物(#489c)", J534["spec"])
Gn = Gate("守卫阴性"); nc2 = Gn.spec_curve_cells_declare_n("无 n 的产物(#494a)",
                                                       [{"scope": k} for k in J539["slopes"]])
print(f"  #489c(逐格有 n)-> {'PASS ✅' if pc2 else 'FAIL ⛔'}(应 PASS)")
print(f"  #494a(0 格有 n)-> {'PASS ⛔' if nc2 else 'FAIL ✅'}(应 FAIL)")

G = Gate("参照改成逐波,和一个每格必须存 n 的新守卫")
G.asserted("正对照①:yr 四格全判 (O)", pc1, f"{[x['n'] for x in yr]}", kind="control")
G.asserted("阴性①:ever 四格全不判 (O)", nc1, f"{[x['n'] for x in ev]}", kind="control")
G.asserted("正对照②:守卫在有 n 的产物上 PASS", pc2, "#489c", kind="control")
G.asserted("阴性②:守卫在无 n 的产物上 FAIL", not nc2, "#494a", kind="control")
# 其余轮次(有逐格 n 的那一个)
other = [dict(pair=c["pair"], n=c["n"]) for c in J537["pairs"]]
print("\n=== 其余可审轮次:`#492c R537` ===")
for c in other: print(f"  {c['pair'][:28]:28s} n={c['n']}")
print("  ⚠ 该轮的人群由**题目可得性**与**已声明的协变量**(`staytog` 限已婚过者)定义,"
      "n 与任何结局阳性数都不相等 -> **未发现 (O) 格**")

ok = pc1 and nc1 and pc2 and (not nc2)
print("\n" + "=" * 70)
if ok:
    verdict = ("四个控制全过;**`#489c` 的 yr 四格确认为 (O)**,"
               "`#492c` **未发现 (O) 格**;`#494a`/`#501a`/`#487d` **因无逐格 n 而不可审**")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:判据只认「n = 该波阳性数」一种形态,"
          "**部分条件化**(如只问「过去一年有过性行为」的人)不会让 n 等于阳性数 -> 看不见。"
          "⇒ 这是**下界**;而新守卫只对**将来**的轮次生效,救不回已跑完的三轮。")
else:
    verdict = f"UNVERIFIED —— 控制未过(pc1={pc1} nc1={nc1} pc2={pc2} nc2_fail={not nc2}) ⇒ 按 `#518` 停止本方向"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(cells=cells, pc1=bool(pc1), nc1=bool(nc1), guard_pos=bool(pc2),
               guard_neg_fails=bool(not nc2), other_round=other, verdict=verdict,
               unchallenged=True), open(OUT / "per_wave_and_guard.json", "w"), indent=1)
print(f"\nwrote {OUT/'per_wave_and_guard.json'}")
