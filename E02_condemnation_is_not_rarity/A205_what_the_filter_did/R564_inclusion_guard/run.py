"""E01·A206·R564 — 产物必须记的不只是 n,还有「这一格由哪些条件共同定义」

`#519` 的 NEXT。**行动类型:PRODUCTION/仪器**(造一条新要求),**不是新实验**。

`#519a` 的三次失败给出根因:真正的纳入形态是「阳性数 **∩** 其余纳入条件」,
而**交集里有哪些条件,产物里一个字也没有**。⇒ 逐格 `n` 不够。

G1 ESTIMAND:新守卫 `spec_curve_cells_declare_inclusion` 在三种产物上的表现:
  ① **手工构造的、带 `inclusion` 的假产物** -> 必须 **PASS**;
  ② **真产物 `#489c`(有 n、无 inclusion)** -> 必须 **FAIL**(证明这是一条**新**要求);
  ③ **真产物 `#494a`(无 n、无 inclusion)** -> 必须 **FAIL**。
⚠ ② 是关键:若它 PASS,说明这条要求**本来就满足**,那这个守卫**什么也没加**。
CONTROLS:①=正对照(守卫能通过)· ②③=阴性(守卫能失败)。**三个都对,守卫才算能用。**
KILL:if 三个都对 -> 守卫可用,写入并记录;else -> 守卫不可用,不写入。
IMPOSSIBLE:守卫只检查**字段在不在**,不检查**写得对不对** ->
  一条写错的 `inclusion` 是一个**新的谎言源**(与 `#516d` 的锚同类风险),必须写明 ·
  它只对**将来**的轮次生效 · 未派对抗 agent ⇒ [unchallenged]
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

FAKE = [dict(cell="ever", n=5895,
             inclusion=["samesex ∈ {1,2,3,4}", "samesexany ∈ {1,5}"]),
        dict(cell="yr", n=1229,
             inclusion=["samesex ∈ {1,2,3,4}", "samyearnum 有效(⚠ 仅问曾有过者)"])]
print("=== 三种产物上的表现 ===")
g1 = Gate("①假产物"); a = g1.spec_curve_cells_declare_inclusion("带 inclusion 的假产物", FAKE)
g2 = Gate("②真产物#489c"); b = g2.spec_curve_cells_declare_inclusion("#489c(有 n 无 inclusion)", J534["spec"])
g3 = Gate("③真产物#494a"); c = g3.spec_curve_cells_declare_inclusion(
    "#494a(无 n 无 inclusion)", [{"scope": k} for k in J539["slopes"]])
print(f"  ① 假产物(带 inclusion) -> {'PASS ✅' if a else 'FAIL ⛔'}(应 PASS)")
print(f"  ② `#489c`(有 n 无 inclusion) -> {'PASS ⛔' if b else 'FAIL ✅'}(应 FAIL —— 证明这是新要求)")
print(f"  ③ `#494a`(两者皆无) -> {'PASS ⛔' if c else 'FAIL ✅'}(应 FAIL)")

G = Gate("产物必须记的不只是 n,还有纳入条件")
G.asserted("正对照:带 inclusion 的产物必须 PASS", a, "手工构造", kind="control")
G.asserted("阴性②:`#489c` 必须 FAIL(证明这是新要求)", not b, "有 n、无 inclusion", kind="control")
G.asserted("阴性③:`#494a` 必须 FAIL", not c, "两者皆无", kind="control")
ok = a and (not b) and (not c)
print("\n" + "=" * 70)
verdict = ("**守卫可用**:三个对照全对 -> 已写入 `lib/gates.py`;"
           "本项目**所有既有产物都会 FAIL**,这正是它可失败的证明" if ok else
           "**守卫不可用**,不写入")
print(f"评判:{verdict}")
print("⚠ 通过的 KILL 会怎样失败:守卫只检查**字段在不在**,不检查**写得对不对** ——"
      "一条写错的 `inclusion` 是一个**新的谎言源**(与 `#516d` 的来源锚同类风险);"
      "而它只对**将来**的轮次生效,救不回已跑完的任何一轮。")
print(G)
json.dump(dict(fake_pass=bool(a), r534_fail=bool(not b), r539_fail=bool(not c),
               guard_usable=bool(ok), verdict=verdict, unchallenged=True),
          open(OUT / "inclusion_guard.json", "w"), indent=1)
print(f"\nwrote {OUT/'inclusion_guard.json'}")
