"""E02·A12·R656 —— 那些锚,引对了吗?(逐位比对页面的数与条目里的数)

`#619` 的 NEXT。**行动类型:CLOSURE**(如实标注)——它保护一个已有结论,不分离任何世界。
`#619e` 的 ⛔:我核对了七条的**原文**,但**从没反向验证过账本里的数与页面上的数逐位相同**。

⚠ **`#618` 刚证明:我造的搜索会失败,而且每次失败的位置不同。**
   ⇒ **本轮先造正对照,再造搜索。任何「一致」的结论,在正对照通过之前都不可采信。**

G1 ESTIMAND(先于方法):对页面上每一个锚 `[#N]` / `[#N「…」]` / `[#N, #M]`,
  取它**所在段落**里的每一个小数(形如 `\\d+\\.\\d{2,4}`),
  判它是否是**它所引条目里某个数的四舍五入**(`#571`:页面 `0.797`,账本 `+0.7966`)。
  **三值**:`一致` / `不一致` / `判不了`(条目里没有任何小数,或页面的数是年份/计数)。
  ⚠ **只查小数**;整数(n、年份、条目号)一律记 `判不了`,**不当成不一致**。

CONTROLS:
  正对照:三个我已手工确认对上的 —— **`#251` 的 0.432 · `#487` 的 0.320 · `#507` 的 0.838** ——
    喂进同一套判定,**必须全部判为一致**。
  **g=0**:把同一段的数字拿去和一个**错误的条目**(`#1`)比,**必须判为不一致**。
    ⇒ 正对照与 g=0 一起,才说明这套判定**既能通过也能失败**。
  安慰剂:把页面的数**各加 0.01** 再判,**一致数必须显著下降**(否则判定太松)。
KILL(条件式,预注册):
  if 正对照 3/3 一致 and g=0 判为不一致 and 安慰剂使一致数下降:
      报 `一致 / 不一致 / 判不了` 三个数,**不一致的逐条列出**
  else: UNVERIFIED —— 判定不合格,任何「一致」都不可采信
G3:全部锚 × 全部小数的判定结果发布,**不一致的逐条列出,不许只报计数**。
G4:容差 = {恰好四舍五入 / ±1 末位 / ±2 末位} 三档。
IMPOSSIBLE(不写 planned):**这是一致性核对,不是正确性核对** ——
  它只能说「页面的数在被引条目里找得到」,**不能说那个条目本身是对的** ·
  段落是锚的作用域,而**一个锚可能只为段落里的一部分数字负责** ⇒ 会有假「不一致」· `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
# 页面侧:与 `uncited_numbers` 同一个口径(2–4 位小数)
NUM = re.compile(r'(?<![\w.])(\d+\.\d{2,4})(?![\w])')
# ⚠ 账本侧必须更宽:`#501` 的表写的是 `−0.05181`(**五位**),而上面那个模式
#   因为负向前瞻 `(?![\w])`,对 5 位以上的小数**完全看不见** —— 于是 `#501` 被判成 0/6。
#   **正对照通过了,因为那三个恰好 ≤4 位** —— realstat:正对照只问「仪器看得见吗」,
#   从不问「它看见的是不是我要断言的那个东西」。
LNUM = re.compile(r'(?<![\w.])(\d+\.\d{2,8})(?![\w])')
ANCH = re.compile(r'\[#(\d+)(?:[^\]]*?)\]|\[#(\d+)\s*,\s*#(\d+)\]')

led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
BODY = {}
for j, l in enumerate(led):
    if ent[j]: BODY.setdefault(ent[j], []).append(l)
BODY = {k: "\n".join(v) for k, v in BODY.items()}
ENTNUM = {k: set(LNUM.findall(v)) for k, v in BODY.items()}
print(f"账本 {len(BODY)} 条条目 · 带小数的 {sum(1 for v in ENTNUM.values() if v)} 条")


def rounds_to(page, cand):
    """页面的数是不是 cand 的四舍五入?"""
    d = len(page.split(".")[1])
    try: return f"{round(float(cand), d):.{d}f}" == page
    except Exception: return False


def judge(page_num, entries, tol=0):
    """三值。tol = 允许末位差几个单位。"""
    pool = set()
    for e in entries: pool |= ENTNUM.get(e, set())
    if not pool: return "判不了", None
    for c in pool:
        if rounds_to(page_num, c): return "一致", c
    if tol:
        d = len(page_num.split(".")[1]); step = 10**(-d)
        for c in pool:
            try:
                if abs(float(c) - float(page_num)) <= tol*step + 1e-12: return "一致", c
            except Exception: pass
    return "不一致", None


def scan(page, tol=0):
    paras = re.split(r'\n\s*\n', pathlib.Path(page).read_text())
    rows = []
    for para in paras:
        es = sorted({int(g) for m in ANCH.finditer(para) for g in m.groups() if g})
        if not es: continue
        for n in NUM.findall(para):
            v, src = judge(n, es, tol)
            rows.append(dict(page=page, anchors=es, num=n, verdict=v, matched=src,
                             snip=para.strip().replace("\n", " ")[:70]))
    return pd.DataFrame(rows)


# ── 控制:先造对照,再信结果 ─────────────────────────────────
G = Gate("那些锚,引对了吗?")
POS = [("0.432", [251]), ("0.320", [487]), ("0.838", [507])]
pos_hits = [judge(n, e)[0] for n, e in POS]
print("\n=== 正对照:三个已手工确认对上的 ===")
for (n, e), v in zip(POS, pos_hits): print(f"  {n} vs #{e[0]} -> **{v}**")
g0 = [judge(n, [1])[0] for n, _ in POS]
print(f"  **g=0:同样三个数拿去和 `#1` 比 -> {g0}**(须全部「不一致」或「判不了」)")
pos_ok = G.positive_control("正对照:三个已确认的必须全部判为一致",
                            planted=float(sum(v == "一致" for v in pos_hits)), floor=1.0, spread=0.5)
g0_bad = sum(v == "一致" for v in g0)
print(f"  g=0 判为「一致」的个数 = {g0_bad}(须 0)")

T = scan("README.md"); TZ = scan("README_zh.md")
ALL = pd.concat([T, TZ], ignore_index=True)
# 安慰剂:把页面的数各加 0.01 再判
def bump(n):
    d = len(n.split(".")[1]); return f"{float(n)+0.01:.{d}f}"
pl = [judge(bump(r.num), r.anchors)[0] for r in ALL.itertuples()]
base_ok = int((ALL.verdict == "一致").sum()); pl_ok = sum(v == "一致" for v in pl)
print(f"\n=== 安慰剂:页面每个数 +0.01 再判 -> 一致数 {base_ok} -> **{pl_ok}**(须显著下降)===")
pla_ok = G.negative_control("安慰剂:+0.01 之后一致数必须下降",
                            null=float(pl_ok), effect=float(base_ok), null_spread=1.0,
                            null_kind="人为把页面的数移开一个末位量级")

print(f"\n=== G3:全部判定({len(ALL)} 个「锚×小数」)===")
print(ALL.verdict.value_counts().to_string())
bad = ALL[ALL.verdict == "不一致"]
print(f"\n--- **不一致的逐条列出({len(bad)} 条),不许只报计数** ---")
for r in bad.itertuples():
    print(f"  ⛔ {r.page:12s} 数 **{r.num}** 锚 {r.anchors}\n      「{r.snip}」")
und = ALL[ALL.verdict == "判不了"]
print(f"\n--- 判不了 {len(und)} 条(所引条目里没有任何小数)· 涉及的锚:"
      f"{sorted({e for r in und.itertuples() for e in r.anchors})[:14]} ---")

if pos_ok and pla_ok and g0_bad == 0:
    verdict = (f"一致 **{int((ALL.verdict=='一致').sum())}** · 不一致 **{len(bad)}** · 判不了 **{len(und)}**")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 判定不合格(正对照 {pos_ok} · 安慰剂 {pla_ok} · g=0 假一致 {g0_bad})"
    print(f"\n⚠ {verdict}")
print(G)

print("\n=== G4 规格曲线:容差三档 ===")
spec = []
for tol in (0, 1, 2):
    A2 = pd.concat([scan("README.md", tol), scan("README_zh.md", tol)], ignore_index=True)
    c = A2.verdict.value_counts().to_dict()
    spec.append(dict(tol=tol, **{k: int(v) for k, v in c.items()}))
    print(f"  容差 ±{tol} 末位: 一致 {c.get('一致',0):3d} · 不一致 {c.get('不一致',0):3d} · 判不了 {c.get('判不了',0):3d}")
json.dump(dict(rows=ALL.to_dict("records"), positive=pos_hits, g0=g0, placebo=[base_ok, pl_ok],
               spec_curve=spec, verdict=verdict, unchallenged=True),
          open(OUT/"do_the_anchors_point_true.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'do_the_anchors_point_true.json'}")
