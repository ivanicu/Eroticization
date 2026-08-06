"""E02·A209·R568 — 给那六条加来源锚,并写一个能失败的检查

`#523` 的 NEXT。**行动类型:PRODUCTION/仪器**。
`#516d` 提出「给承重断言加显式来源锚」,并同时写下它的代价:
**锚一旦写错就是一个新的谎言源。** ⇒ 加锚**必须**与一个能失败的检查同时落地。

G1 ESTIMAND:对新段里每一条带锚的断言,检查
  **锚指向的账本条目,其正文是否含该条引用的至少一个数**。
判据(预注册):全部指得到 -> 锚可用;有指不到的 -> **那一条的来源是我记错的**,当场修。
CONTROLS:
  ⛔ **正对照(先跑):把一条锚故意改成一个不相干的条目号,检查必须报错。**
     不报错 -> 这个检查不会失败 -> 它不是检查(realstat「不会失败的检查」)。
  阴性:真实的六条锚在**未被篡改**时应当全部通过。
IMPOSSIBLE:检查只验「数在不在那条里」,**不验那条说的是不是这件事** ->
  与 `#512e` 同一种单向性:**只能证伪** · 未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
led = pathlib.Path("RETRACTIONS.md").read_text()
marks = [(int(m.group(1)), m.start()) for m in re.finditer(r"^## Entry (\d+)", led, re.M)]
bodies = {}
for i, (n, s) in enumerate(marks):
    bodies[n] = led[s:(marks[i+1][1] if i+1 < len(marks) else len(led))]
NUM = re.compile(r"[−+-]?\d+\.\d+")

def check(page_text, tamper=None):
    """返回 [(bullet_idx, anchors, nums, ok)]。tamper=(idx, fake_anchor) 用于正对照。"""
    i = page_text.find("# What the second half was for")
    if i < 0: i = page_text.find("# 下半场是为了什么")
    j = page_text.find("\n# ", i + 5)
    sec = page_text[i:j if j > 0 else len(page_text)]
    out = []
    for k, line in enumerate([l for l in sec.split("\n- ") if "`[#" in l]):
        anch = [int(x) for x in re.findall(r"#(\d+)", line)]
        if tamper and tamper[0] == k: anch = [tamper[1]]
        ns = [t.replace("−", "-") for t in NUM.findall(line)]
        ok = any(any(f in bodies.get(a, "") for f in {t, t.replace("-", "−")}) for a in anch for t in ns) if ns else None
        out.append((k, anch, ns[:4], ok))
    return out

page = pathlib.Path("README.md").read_text()
print("=== ⛔ 正对照:把第 0 条的锚篡改成 #1(不相干)===")
tam = check(page, tamper=(0, 1))
print(f"  篡改后第 0 条 -> anchors={tam[0][1]} nums={tam[0][2]} ok={tam[0][3]}")
pc = (tam[0][3] is False)
print(f"  检查{'报错 ✅(它能失败)' if pc else '没报错 ⛔(它不是检查)'}")

print("\n=== 真实六条锚 ===")
real = check(page)
for k, a, ns, ok in real:
    print(f"  [{k}] 锚 {a}  引用数 {ns}  -> {'PASS ✅' if ok else 'FAIL ⛔' if ok is False else '无数字(跳过)'}")
bad = [r for r in real if r[3] is False]

G = Gate("给那六条加来源锚,并写一个能失败的检查")
G.asserted("⛔ 正对照:篡改锚之后检查必须报错", pc, "第 0 条锚改成 #1", kind="control")
G.asserted("阴性:真实六条锚全部指得到", not bad, f"FAIL {len(bad)} 条", kind="control")
print("\n" + "=" * 70)
verdict = ("**锚可用**:篡改能被抓,真实六条全部指得到 -> 可逐段推广" if pc and not bad else
           (f"⛔ **{len(bad)} 条锚指不到 -> 那些来源是我记错的,当场修**" if pc else
            "UNVERIFIED —— 检查不会失败,它不是检查"))
print(f"评判:{verdict}")
print("⚠ 通过的 KILL 会怎样失败:检查只验「数在不在那条里」,**不验那条说的是不是这件事** ——"
      "与 `#512e` 同一种单向性,**只能证伪**;一条指得到但语义错的锚,它抓不到。")
print(G)
json.dump(dict(tamper_detected=bool(pc), rows=[dict(idx=k, anchors=a, nums=ns, ok=ok)
                                               for k, a, ns, ok in real],
               n_fail=len(bad), verdict=verdict,
               inclusion=["README.md 的『下半场是为了什么』小节", "仅带 `[#NNN]` 锚的条目"],
               unchallenged=True), open(OUT / "anchor_check.json", "w"), indent=1)
print(f"\nwrote {OUT/'anchor_check.json'}")
