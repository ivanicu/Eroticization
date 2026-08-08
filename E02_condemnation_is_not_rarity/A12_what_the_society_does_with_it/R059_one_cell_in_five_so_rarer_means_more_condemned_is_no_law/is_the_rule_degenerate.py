"""E02·A12·R657 —— 那条规则的四格,是真不一致,还是退化命中?

`#620` 的 NEXT。**行动类型:CLOSURE**(如实标注)。
⛔ **判据在 `#620` 已写死,本轮只执行**:**只有四格全是退化命中,才允许改规则**;
   而改规则**必须先用 git 历史回测**,不许为了让计数变绿而放松。

规则实况(读源码,不是读它的名字):`internal_consistency` 是**行级**的 ——
按「该行含不含 CJK 字符」把同一引用标记的出现分成两侧,比较两侧数字的**并集**。
⇒ 它真正的语义是「**中文表**里写的数 vs **英文正文**里写的数」(`#144` 的原意)。

G1 ESTIMAND(先于方法):对当前每一格,算两侧并集 `S_cjk` / `S_lat`,判三值:
  **退化命中** = 有一侧为**空集** · **真不一致** = 两侧都非空**且**存在一个数只在一侧 ·
  **判不了** = 两侧都非空且差异只来自**四舍五入位数不同**(如 `0.797` vs `0.7966`)。

CONTROLS:
  正对照:**人为制造一处真不一致** —— 把中文页某处的 `0.432` 改成 `0.532` ——
    规则必须抓到 ⇒ **计数必须上升**。
  **g=0**:恢复后必须**回到原值**。⇒ 两者一起才说明规则**既能开火也能熄火**。
  安慰剂:在一个**不带任何数字**的位置加一个新锚 -> 若计数上升,**说明规则对退化命中敏感**(这正是被告的行为)。
KILL(条件式,预注册于 `#620`):
  if 正对照使计数上升 and g=0 使它回落:
      四格**全是**退化命中 -> 允许改规则,且**必须回测**
      否则 -> **不改规则**,把分类结果如实登记
  else: UNVERIFIED
G3:四格逐格发布两侧并集与判定。G4:分侧口径 {CJK/非CJK} × {行级/段级} 两档对比。
IMPOSSIBLE(不写 planned):这是**页面自洽性**核对,**不核对页面与账本**(那是 `#620` 做的) ·
  「两侧」是 CJK 的代理,而**中文短语锚会把英文行推到中文侧** ⇒ 代理本身会漂 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r'(?<![\w.])\d+\.\d{2,4}(?![\w])')


def sides(page="README.md"):
    lines = pathlib.Path(page).read_text().splitlines()
    cj = lambda s: bool(re.search(r'[一-鿿]', s))
    cites = {}
    for i, l in enumerate(lines, 1):
        for c in set(re.findall(r'#\d{1,3}\b|\bA\d{2}\b|Entry \d{1,3}', l)):
            cites.setdefault(c, []).append((i, cj(l), set(NUM.findall(l))))
    out = {}
    for c, rows in cites.items():
        S = {True: set(), False: set()}
        for i, k, ns in rows: S[k] |= ns
        if len(set(k for _, k, _ in rows)) == 2 and S[True] != S[False]:
            out[c] = (S[True], S[False], rows)
    return out


def classify(a, b):
    if not a or not b: return "退化命中(一侧空集)"
    only_a = a - b; only_b = b - a
    def rounded_pair(x, y):
        for p in x:
            d = len(p.split(".")[1])
            for q in y:
                try:
                    if f"{round(float(q), d):.{d}f}" == p: return True
                except Exception: pass
        return False
    if rounded_pair(only_a, b) and rounded_pair(only_b, a): return "判不了(只是四舍五入位数不同)"
    return "真不一致"


S = sides()
print(f"=== G3:当前 {len(S)} 格,逐格两侧并集与判定 ===")
rows = []
for c, (cjk, lat, rr) in sorted(S.items()):
    v = classify(cjk, lat)
    rows.append(dict(cite=c, n_cjk=len(cjk), n_lat=len(lat), verdict=v,
                     only_cjk=sorted(cjk - lat)[:6], only_lat=sorted(lat - cjk)[:6],
                     lines=[r[0] for r in rr]))
    print(f"\n  **{c}** -> {v}")
    print(f"    中文侧 {len(cjk)} 个数 · 英文侧 {len(lat)} 个数 · 出现在行 {[r[0] for r in rr]}")
    print(f"    只在中文侧:{sorted(cjk-lat)[:8]}")
    print(f"    只在英文侧:{sorted(lat-cjk)[:8]}")
T = pd.DataFrame(rows)
n_deg = int((T.verdict.str.startswith("退化")).sum())
n_real = int((T.verdict == "真不一致").sum())
n_und = int((T.verdict.str.startswith("判不了")).sum())
print(f"\n**退化命中 {n_deg} · 真不一致 {n_real} · 判不了 {n_und}(共 {len(T)})**")

# ── 控制 ─────────────────────────────────────────────────────
G = Gate("那条规则的四格,是真不一致,还是退化命中?")
# ⛔ 第一版正对照放在**中文页**,而它永远开不了火 —— 见下面这行实测:
#   规则按「行含不含 CJK」分两侧,中文页几乎每行都含 CJK ⇒ 只有一侧 ⇒ `len(set)==2` 永假。
#   **这条规则在中文页上结构性地失明**,而这本身就是本轮最硬的一件事。
print(f"\n  ⛔ 结构性失明:`README_zh.md` 给出 **{len(sides('README_zh.md'))} 格** ——"
      f" 规则按 CJK 分侧,中文页只有一侧,**它从来没有检查过中文页**。")
# 正对照必须放在规则**能**开火的地方(英文页),且改的必须是**已在两侧出现**的标记所在行的数
p2 = pathlib.Path("README.md"); o2 = p2.read_text()
base = len(sides())
target = None
for c, (cjk, lat, rr) in sides().items():
    for n in sorted(lat):
        if o2.count(n) >= 1: target = (c, n); break
    if target: break
mutn = None
if target:
    c, n = target
    d = len(n.split(".")[1]); mutn = f"{float(n)+0.1:.{d}f}"
    p2.write_text(o2.replace(n, mutn, 1))
after_pos = len(sides()); p2.write_text(o2); back = len(sides())
print(f"  正对照(英文页):把一处 `{target[1] if target else '?'}` 改成 `{mutn}` -> 格数 {base} -> **{after_pos}**")
print(f"  g=0:恢复后 -> **{back}**({'✅ 回落' if back == base else '⛔ 没有回落'})")
pos_ok = G.positive_control("正对照:在规则能开火的页面上人为制造一处不一致,计数必须上升",
                            planted=float(after_pos - base), floor=0.0, spread=0.4)
p3 = pathlib.Path("README.md"); o3 = p3.read_text()
p3.write_text(o3.replace("\n## ", "\n`[#999]`\n\n## ", 1))
after_pla = len(sides()); p3.write_text(o3)
print(f"  安慰剂:在一个**不带数字**的位置加锚 `[#999]` -> 格数 {base} -> **{after_pla}**")
pla_ok = G.negative_control("安慰剂:加一个不带数字的锚", null=float(abs(after_pla - base)),
                            effect=float(max(after_pos - base, 1)), null_spread=0.4,
                            null_kind="一个不携带任何数字的新锚")
base_zh, after_pos_zh = len(sides('README_zh.md')), 0

if pos_ok:
    if n_deg == len(T):
        verdict = "四格全是退化命中 -> **允许改规则,且必须回测**"
    else:
        verdict = (f"**不改规则**:退化 {n_deg}/{len(T)},真不一致 {n_real},判不了 {n_und} —— "
                   f"判据要求「全是退化」才允许改,条件**不成立**")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 正对照未触发({after_pos - base_zh})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4:两个页面各自的格数 ===")
spec = [dict(page=f, cells=len(sides(f))) for f in ("README.md", "README_zh.md")]
for s in spec: print(f"  {s['page']}: {s['cells']} 格")
json.dump(dict(cells=T.to_dict("records"), n_degenerate=n_deg, n_real=n_real, n_undecided=n_und,
               positive=[base, after_pos, back], zh_cells=base_zh, placebo=[base, after_pla], spec=spec,
               verdict=verdict, unchallenged=True),
          open(OUT/"is_the_rule_degenerate.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'is_the_rule_degenerate.json'}")
