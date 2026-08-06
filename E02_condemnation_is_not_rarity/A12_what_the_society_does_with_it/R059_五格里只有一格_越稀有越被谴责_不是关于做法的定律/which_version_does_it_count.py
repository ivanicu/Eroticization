"""E02·A12·R662 —— 每条规则,它数的是哪一版?(量出来,不是声明出来)

`#625` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
本次会话里「计数取决于看哪一版」已同型**三次**(`#607b` 分票零重叠 · `#622` 规则只在半张纸上工作 ·
`#625d` 删限定语只在中文页)。**第四次不该是再记一笔,该是把它变成规则的一个属性。**

⚠ **关键设计选择,写在前面:「它读了哪一版」必须是量出来的,不是声明出来的。**
   声明会错(`#602` 就是靠读代码发现两条规则写死了 `README.md`,而那是**读**出来的,不是**测**出来的)。
   ⇒ 本轮**拦截文件读取本身**:在规则运行期间替换 `pathlib.Path.read_text`,记录它实际打开了哪些文件。

G1 ESTIMAND(先于方法):对 `readme_gate` 的每一条规则,
  `覆盖版本数` = 该规则运行期间**实际读取**的 `README*.md` 文件个数(0 / 1 / 2)。
CONTROLS:
  正对照:一个**只读 `README.md`** 的合成规则 -> 必须报 **1**。
  **g=0**:一个**两版都读**的合成规则 -> 必须报 **2**(否则这个探针分不出 1 和 2)。
  安慰剂:一个**一个页面都不读**的合成规则(只读账本)-> 必须报 **0**。
KILL(条件式,预注册):
  if 正对照=1 and g=0=2 and 安慰剂=0:
      逐条报覆盖版本数,并把它**打进闸门输出**
  else: UNVERIFIED —— 探针不可信,任何覆盖数都不采信
G3:七条规则逐条发布。G4:探针 {只拦 `read_text` / 同时拦 `open`} 两档,比较是否漏读。
IMPOSSIBLE(不写 planned):**探针只看得见经过 `pathlib.Path.read_text` 的读取** ——
  经 `open()`/`subprocess`/缓存的读取它看不见(G4 的第二档就是为量这个) ·
  「覆盖两版」**不等于**「对两版同样有效」(`#622`:规则读了中文页,但在它上面结构性开不了火)·
  `[unchallenged]`
"""
import os, sys, pathlib, json, warnings, builtins
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_gate as G_, readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
_ORIG_RT = pathlib.Path.read_text
_ORIG_OPEN = builtins.open


class Probe:
    """在一段代码运行期间,记录它实际读了哪些 README*。"""
    def __init__(self, catch_open=False): self.seen = set(); self.catch_open = catch_open
    def __enter__(self):
        seen = self.seen
        def rt(self_, *a, **k):
            n = pathlib.Path(self_).name
            if n.startswith("README") and n.endswith(".md"): seen.add(n)
            return _ORIG_RT(self_, *a, **k)
        pathlib.Path.read_text = rt
        if self.catch_open:
            def op(f, *a, **k):
                try: n = pathlib.Path(f).name
                except Exception: n = ""
                if n.startswith("README") and n.endswith(".md"): seen.add(n)
                return _ORIG_OPEN(f, *a, **k)
            builtins.open = op
        return self
    def __exit__(self, *e):
        pathlib.Path.read_text = _ORIG_RT; builtins.open = _ORIG_OPEN


def cover(fn, catch_open=False):
    with Probe(catch_open) as p:
        try: fn()
        except Exception as e: return sorted(p.seen), f"错误 {type(e).__name__}"
    return sorted(p.seen), "ok"


# ── 控制先跑 ─────────────────────────────────────────────────
G = Gate("每条规则,它数的是哪一版?")
pos = cover(lambda: pathlib.Path("README.md").read_text())
g0 = cover(lambda: (pathlib.Path("README.md").read_text(), pathlib.Path("README_zh.md").read_text()))
pla = cover(lambda: pathlib.Path("RETRACTIONS.md").read_text())
print("=== 控制 ===")
print(f"  正对照(只读 README.md)      -> {len(pos[0])} 版 {pos[0]}(须 1)")
print(f"  g=0  (两版都读)             -> {len(g0[0])} 版 {g0[0]}(须 2)")
print(f"  安慰剂(一个页面都不读)       -> {len(pla[0])} 版 {pla[0]}(须 0)")
pos_ok = G.positive_control("探针必须能分出 1 和 2", planted=float(len(g0[0]) - len(pos[0])),
                            floor=0.0, spread=0.4)
pla_ok = G.negative_control("安慰剂:不读页面的规则必须报 0", null=float(len(pla[0])),
                            effect=float(len(g0[0])), null_spread=0.4,
                            null_kind="只读账本、不碰页面的规则")

# ── 主测:逐条规则 ───────────────────────────────────────────
PAGES = ("README.md", "README_zh.md")
RULES = {
 "named_defects":        lambda: A.named_defects(),
 "numbers_that_left":    lambda: A.numbers_that_left(rev="HEAD~1"),
 "uncited_numbers":      lambda: [A.uncited_numbers(p) for p in PAGES],
 "internal_consistency": lambda: [A.internal_consistency(p) for p in PAGES],
 "dangling_anchors":     lambda: G_.dangling_anchors(),
 "claims_without_anchor":lambda: G_.claims_page_edit_without_anchor(),
 "qualifiers_stripped(未接入)": lambda: A.qualifiers_stripped("HEAD~1"),
}
rows = []
print(f"\n=== G3:{len(RULES)} 条规则,逐条量它实际读了哪几版 ===")
for name, fn in RULES.items():
    f1, s1 = cover(fn, False); f2, s2 = cover(fn, True)
    rows.append(dict(rule=name, n_versions=len(f1), files=f1, status=s1,
                     n_with_open=len(f2), files_with_open=f2))
    mark = "✅ 两版" if len(f1) == 2 else ("⚠ 只一版" if len(f1) == 1 else "⛔ 零版")
    print(f"  {name:26s} {mark}  读到 {f1}  ·(加拦 open 后 {f2})  {s1}")
T = pd.DataFrame(rows)
n2 = int((T.n_versions == 2).sum()); n1 = int((T.n_versions == 1).sum()); n0 = int((T.n_versions == 0).sum())
print(f"\n**两版 {n2} 条 · 只一版 {n1} 条 · 零版 {n0} 条(共 {len(T)})**")
diff_open = T[T.n_versions != T.n_with_open]
print(f"  G4:加拦 `open` 后覆盖数变化的规则:**{len(diff_open)}** 条"
      f"{'' if not len(diff_open) else ' -> ' + str(list(diff_open.rule))}")

if pos_ok and pla_ok:
    verdict = f"两版 {n2} · 只一版 {n1} · 零版 {n0} ⇒ **把覆盖版本数打进闸门输出**"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 探针不可信(正对照 {pos_ok} · 安慰剂 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(rules=T.to_dict("records"), n_two=n2, n_one=n1, n_zero=n0,
               controls=dict(pos=pos[0], g0=g0[0], placebo=pla[0]),
               open_probe_diff=list(diff_open.rule), verdict=verdict, unchallenged=True),
          open(OUT/"which_version.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'which_version.json'}")
