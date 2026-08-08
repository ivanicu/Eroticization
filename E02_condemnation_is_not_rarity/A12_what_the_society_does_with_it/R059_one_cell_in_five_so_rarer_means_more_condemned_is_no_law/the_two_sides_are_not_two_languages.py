"""E02·A12·R658 —— 那条规则的「两侧」,根本不是两种语言

`#621` 的 NEXT。**行动类型:CLOSURE**(如实标注)。
判据在 `#621` 已写死:逐个判定 `#129` 那格的数字是「中英两版同一主张的两次陈述」,
还是「整行并集带进来的邻居数字」;**若全部是邻居 ⇒ 记「`#129` 也是行并集伪影」,真不一致数改为 0**。
⚠ 判据还要求:**必须报「同一主张」是怎么判的 —— 我判的,不是数据判的。**

G1 ESTIMAND(先于方法,两个):
  E1 `#129` 两侧的行,**是不是同一条主张**?判法先写死:
     **取每一行第一个表格单元格(主张列)的文字**,两行文字不同 ⇒ **不是同一条主张**。
     ⚠ 这是**我定的判法**,不是数据给的;它可失败之处:一条主张可能被拆成两行。
  E2 **「CJK 侧」到底装的是什么?** 在 `README.md` 上统计:
     含 CJK 的行里,有多少是**整行中文**,有多少只是**英文行里嵌了一个中文短语锚**。
     ⇒ 若后者占压倒多数,则「两侧 = 两种语言」这个前提**在这一页上是假的**。

WORLDS:**A** 两侧确实是两种语言(前提成立,`#129` 可能是真缺陷)·
  **B** 两侧只是「有没有中文短语锚」(前提为假 ⇒ 整个计数是伪影)
CONTROLS:
  正对照:人为造一行**整行中文**且带某个已有标记 -> `含CJK且整行中文` 的计数必须 +1。
  **g=0**:造一行英文但嵌一个中文短语锚 -> **`整行中文` 计数必须不变**(否则这个分类器本身没用)。
KILL(条件式,预注册):
  if 正对照 +1 and g=0 不变:
      `#129` 两行主张列文字不同 **且** CJK 侧压倒多数是「嵌锚的英文行」-> **W-B,真不一致数改为 0**
      否则 -> 逐条对齐 `#129`
  else: UNVERIFIED
G3:`README.md` 全部含 CJK 行的分类全表。G4:判「整行中文」的阈值 {30%, 50%, 70%} 三档。
IMPOSSIBLE(不写 planned):**这是对规则的核对,不是对页面内容的核对** ·
  「同一主张」的判法是**我定的** · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r'(?<![\w.])\d+\.\d{2,4}(?![\w])')
CJK = re.compile(r'[一-鿿]')
ANCHOR = re.compile(r'`\[#\d+[^`]*?\]`')


def cjk_share(line, strip_anchor=True):
    s = ANCHOR.sub("", line) if strip_anchor else line
    if not s.strip(): return 0.0
    return len(CJK.findall(s)) / max(len(s), 1)


lines = pathlib.Path("README.md").read_text().splitlines()

# ── E1:`#129` 两行是不是同一条主张 ────────────────────────────
print("=== E1:`#129` 两行的**主张列**(每行第一个表格单元格)===")
rows129 = []
for i, l in enumerate(lines, 1):
    if re.search(r'#129\b', l):
        cells = [c.strip() for c in l.split("|") if c.strip()]
        claim = cells[0] if cells else l[:60]
        rows129.append(dict(line=i, cjk=bool(CJK.search(l)), nums=len(NUM.findall(l)),
                            claim=claim[:96]))
        print(f"  行 {i} · {'CJK侧' if CJK.search(l) else '英文侧'} · {len(NUM.findall(l))} 个数")
        print(f"    主张列:「{claim[:92]}」")
same_claim = len({r["claim"] for r in rows129}) == 1
print(f"\n  **两行主张列文字相同?{same_claim}** ⇒ "
      f"{'同一条主张' if same_claim else '**不是同一条主张 —— 是同一张表里的两行不同主张**'}")
print("  ⚠ 判法是**我定的**:取每行第一个表格单元格。可失败之处:一条主张若被拆成两行,这个判法会误判。")

# ── E2:CJK 侧装的是什么 ─────────────────────────────────────
print("\n=== E2 · G3:`README.md` 里含 CJK 的行,到底是什么 ===")
cat = []
for i, l in enumerate(lines, 1):
    if not CJK.search(l): continue
    sh = cjk_share(l)
    kind = "整行中文" if sh >= 0.30 else ("嵌了中文短语锚的英文行" if ANCHOR.search(l) and CJK.search(ANCHOR.search(l).group(0)) else "其他(少量中文)")
    cat.append(dict(line=i, share=round(sh, 4), kind=kind, n=len(NUM.findall(l))))
C = pd.DataFrame(cat)
print(C.kind.value_counts().to_string())
print(f"\n  含 CJK 的行共 **{len(C)}**;其中**整行中文 {int((C.kind=='整行中文').sum())} 行**,"
      f"**嵌锚的英文行 {int((C.kind=='嵌了中文短语锚的英文行').sum())} 行**")
print(f"  这些行携带的数字总数:整行中文 {int(C[C.kind=='整行中文'].n.sum())} · "
      f"嵌锚英文行 **{int(C[C.kind=='嵌了中文短语锚的英文行'].n.sum())}**")

# ── 控制 ────────────────────────────────────────────────────
G = Gate("那条规则的「两侧」,是不是两种语言?")
base = int((C.kind == "整行中文").sum())
pos = cjk_share("| **这是一整行中文的主张** | 它带着数字 0.123 与 0.456,并且引用 `[#129]` |")
neg = cjk_share("| **This row is English** | it carries 0.123 and cites `[#129「一个中文短语锚」]` |")
print(f"\n  正对照:一整行中文 -> CJK 占比 **{pos:.3f}**(须 ≥0.30)")
print(f"  g=0:英文行只嵌一个中文短语锚(锚已剥离)-> CJK 占比 **{neg:.3f}**(须 <0.30)")
pos_ok = G.positive_control("正对照:整行中文必须被判为整行中文", planted=float(pos), floor=0.30, spread=0.02)
pla_ok = G.negative_control("g=0:嵌锚的英文行必须不被判为整行中文", null=float(neg), effect=float(pos),
                            null_spread=0.02, null_kind="英文行里嵌一个中文短语锚")

dominated = int((C.kind == "嵌了中文短语锚的英文行").sum()) > int((C.kind == "整行中文").sum())
if pos_ok and pla_ok:
    if (not same_claim) and dominated:
        verdict = ("W-B:**两侧不是两种语言 —— `#129` 也是行并集伪影,四格的真不一致数改为 0**")
    else:
        verdict = f"逐条对齐 `#129`(same_claim={same_claim} · CJK侧被嵌锚行主导={dominated})"
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:判「整行中文」的阈值三档 ===")
spec = []
for th in (0.30, 0.50, 0.70):
    n_full = sum(1 for _, r in C.iterrows() if r.share >= th)
    spec.append(dict(th=th, full_cjk=int(n_full), embedded=int(len(C) - n_full)))
    print(f"  阈值 {th:.2f}: 整行中文 {n_full:3d} · 其余 {len(C)-n_full:3d}")
json.dump(dict(rows129=rows129, same_claim=bool(same_claim), categories=C.kind.value_counts().to_dict(),
               n_cjk_lines=len(C), positive=float(pos), g0=float(neg), spec=spec,
               verdict=verdict, unchallenged=True),
          open(OUT/"two_sides_not_two_languages.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'two_sides_not_two_languages.json'}")
