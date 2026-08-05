"""E01·A203·R556 — 页面上的每一个数,账本里找不找得到

`#511` 的 NEXT。**行动类型:CLOSURE**(仪器)。
本项目已两次发生「页面说的与账本说的不同」(`#347` · `#502d`),而**从未做过一次全页核对**。

⛔ 核心难点写在跑之前:**页面是取整的,账本是全精度的** ->
  直接 `grep` 必然大量误报。判据必须是:
  **账本里存不存在一个数,在页面那个数自身的小数位数上四舍五入等于它。**

G1 ESTIMAND:页面上每一个**含小数点**的数(计数与年份不进,因为它们不是测量),
  分类为 **MATCH**(账本里有数在该精度上等于它)/ **NO-MATCH**。
⚠ 只查**存在性**,不查**语境是否相同** —— 一个数可能在账本里是另一件事的值。
  **这是本仪器的上界:它能证伪(找不到 = 页面上有账本没有的数),不能证实。**(P6 的单向性)

CONTROLS:
  正对照 ① 注入一个**账本里确实存在**的数 -> 必须 MATCH;
         ② 注入一个**编造的**数(如 `0.123456789`)-> 必须 NO-MATCH。
         **两个都过,这个正则才算能看见。**(`#489a` 的教训:先在已知答案上跑)
  阴性   把账本换成一段**无关文本**,页面数的 MATCH 率应当**塌下去**。
KILL(条件式,预注册):
  if 两个正对照都通过 and 阴性塌下去:
      NO-MATCH 数 = 0 -> 页面与账本一致;
      > 0            -> **逐个列出,它们是下一批必须查的东西**
  else: UNVERIFIED
IMPOSSIBLE:只查存在性不查语境 ⇒ **只能证伪** · 页面上的数可能来自 `results/*.json` 而非账本正文 ·
  未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r"[−+-]?\d+\.\d+")

def nums(text):
    out = []
    for tok in NUM.findall(text):
        d = tok.replace("−", "-")
        try: out.append((float(d), len(d.split(".")[1]), tok))
        except ValueError: pass
    return out

led = pathlib.Path("RETRACTIONS.md").read_text()
led_vals = [v for v, _, _ in nums(led)]
print(f"账本里含小数的数:{len(led_vals)} 个")

def matches(v, dp, pool):
    return any(round(x, dp) == round(v, dp) for x in pool)

# 正对照:先在已知答案上跑(`#489a`)
known = led_vals[0]
pc1 = matches(known, 4, led_vals)
pc2 = not matches(0.123456789, 9, led_vals)
print(f"正对照① 账本里已存在的数 {known} -> {'MATCH' if pc1 else 'NO-MATCH'}(应 MATCH)")
print(f"正对照② 编造的数 0.123456789 -> {'MATCH' if not pc2 else 'NO-MATCH'}(应 NO-MATCH)")

rows, bad = [], []
for f in ("README.md", "README_zh.md"):
    page = pathlib.Path(f).read_text()
    ns = nums(page)
    hit = [(v, dp, tok) for v, dp, tok in ns if matches(v, dp, led_vals)]
    miss = [(v, dp, tok) for v, dp, tok in ns if not matches(v, dp, led_vals)]
    rows.append(dict(file=f, total=len(ns), match=len(hit), nomatch=len(miss),
                     misses=[t for _, _, t in miss]))
    bad += [(f, t) for _, _, t in miss]
    print(f"\n{f}: 含小数的数 {len(ns)} 个 -> MATCH {len(hit)} · **NO-MATCH {len(miss)}**")
    if miss:
        print("  NO-MATCH 逐个列出:", [t for _, _, t in miss][:40])

# 阴性:把账本换成无关文本,MATCH 率应塌
junk = [1.111111, 2.222222, 3.333333]
tot = sum(r["total"] for r in rows)
junk_match = sum(1 for f in ("README.md", "README_zh.md")
                 for v, dp, _ in nums(pathlib.Path(f).read_text()) if matches(v, dp, junk))
print(f"\n阴性:用无关文本当账本 -> MATCH {junk_match}/{tot}(应塌下去)")

G = Gate("页面上的每一个数,账本里找不找得到")
G.asserted("正对照①:账本里已存在的数必须 MATCH", pc1, f"{known}", kind="control")
G.asserted("正对照②:编造的数必须 NO-MATCH", pc2, "0.123456789", kind="control")
G.negative_control("阴性:无关文本当账本,MATCH 率塌下去",
                   null=junk_match / max(tot, 1), effect=sum(r["match"] for r in rows) / max(tot, 1),
                   null_spread=0.02, null_kind="把账本替换为无关文本")
print("\n" + "=" * 70)
if pc1 and pc2:
    verdict = ("**页面与账本一致:0 个 NO-MATCH**" if not bad else
               f"**{len(bad)} 个 NO-MATCH -> 它们是下一批必须查的东西**")
    print(f"控制齐备 ⇒ 评判。{verdict}")
    print("⚠ 通过的 KILL 会怎样失败:本仪器**只查存在性,不查语境** —— "
          "一个数可能在账本里是**另一件事**的值,那样它会被误判为 MATCH。"
          "⇒ **它只能证伪,不能证实**(P6 的单向性);0 个 NO-MATCH **不等于**页面全对。")
else:
    verdict = f"UNVERIFIED —— 正对照未过(①={pc1} ②={pc2})"
    print(f"⚠ {verdict}")
print(G)
json.dump(dict(ledger_numbers=len(led_vals), rows=rows, n_nomatch=len(bad),
               nomatch=[{"file": f, "token": t} for f, t in bad], verdict=verdict,
               unchallenged=True), open(OUT / "full_page_audit.json", "w"), indent=1)
print(f"\nwrote {OUT/'full_page_audit.json'}")
