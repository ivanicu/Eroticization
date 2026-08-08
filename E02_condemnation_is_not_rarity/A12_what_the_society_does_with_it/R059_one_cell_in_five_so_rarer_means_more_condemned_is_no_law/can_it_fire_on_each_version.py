"""E02·A12·R663 —— 每条规则,在每一版页面上开得了火吗?(植入检验)

`#626` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
`#626c` 的 ⛔:我量的是「它往哪里看」,而出错的是「它能看见什么」。本轮把后者变成检查。

⚠ **§3 梯度检查先定操作化,而这决定了对照怎么搭:**
   好几条规则**基线本来就非零**(`uncited_numbers=1` · `internal_consistency=4`)⇒
   **「开火」不能读成「计数>0」,只能读成「计数上升」**;
   而 **g=0 不能是「不注入」**(那必然不上升,是个不可能失败的对照)——
   **g=0 必须是「注入一个该规则不该抓的东西」**。

G1 ESTIMAND(先于方法):对每条规则 R × 每一版页面 P,
  `fires(R,P)` = 在 P 里注入一个**最小的、R 本应抓到的合成缺陷**后,R 的计数是否**上升**。
  **不上升 = R 在 P 上失明。**
CONTROLS:
  正对照 = 注入本身(每格自带)· **g=0** = 在同一版注入一个**该规则不该抓的**东西,计数**必须不上升**。
KILL(条件式,预注册):
  if 某格的 g=0 不上升:
      注入后上升 -> **该版有效** · 不上升 -> **该版失明,写进页面「做不到什么」**
  else: 该格记 **UNVERIFIED**(g=0 都会上升 ⇒ 这个注入检验对该格没有分辨力)
⛔ **不许因为发现失明就悄悄修规则** —— `#620` 的回测规矩仍然有效。
G3:全表发布,**含失明的格与 `UNCOMPUTED` 的格**。
IMPOSSIBLE(不写 planned):
  · `named_defects` / `numbers_that_left` / `claims_without_anchor` 的缺陷**不住在页面里**
    (分别住在「账本点名 × 页面残留」「与上一版的差」「账本条目有没有锚」)⇒
    **「在这一版上开得了火吗」对它们不是良定义的**,记 `UNCOMPUTED` 并说明它需要什么;
  · `qualifiers_stripped` **仍未接入闸门**(`#624a` 第三次)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_gate as G_, readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
PAGES = ("README.md", "README_zh.md")

# 每条可植入规则:(计数函数, 该抓的注入, 不该抓的注入)
INJ = {
 "uncited_numbers": (
    lambda p: len(A.uncited_numbers(p)),
    "\n\n这一段带着一个没有出处的数 +0.4242,而它不引用任何条目。\n",
    "\n\n这一段不带任何小数,只是一句普通的话。\n"),
 "internal_consistency": (
    lambda p: len(A.internal_consistency(p)),
    "\n\n中文一行 `#251` 写着 0.9999。\n`#251` in an English line says 0.1111.\n",
    "\n\n一句不带引用标记也不带小数的普通话。\n"),
 "dangling_anchors": (
    lambda p: sum(1 for x in G_.dangling_anchors() if x[3] == "missing"),
    "\n\n这一段带着一个指向账本里不存在句子的锚 `[#251「这句话账本里绝对没有」]`。\n",
    "\n\n这一段带着一个已经存在的锚 `[#251]`,不该被判为悬空。\n"),
}
UNCOMPUTED = {
 "named_defects": "缺陷 = 「账本点名过的缺陷仍在页面上原样活着」⇒ 需要同时伪造账本条目与页面残留,不是单版注入",
 "numbers_that_left": "缺陷 = 「与上一版相比某个数离开了页面」⇒ 需要伪造一次提交,不是单版注入",
 "claims_without_anchor": "缺陷 = 「账本条目没有锚」⇒ 它住在账本里,不住在任何一版页面里",
 "qualifiers_stripped(未接入)": "仍未接入闸门(`#624a` 第三次);且它的缺陷跨版本比较,不是单版注入",
}

rows = []
print("=== G3:可植入的规则 × 两版 ===")
for rule, (cnt, bad, good) in INJ.items():
    for P in PAGES:
        p = pathlib.Path(P); orig = p.read_text()
        try:
            c0 = cnt(P)
            p.write_text(orig + bad); c_bad = cnt(P)
            p.write_text(orig + good); c_good = cnt(P)
        finally:
            p.write_text(orig)
        fired = c_bad > c0; g0_clean = c_good <= c0
        v = ("UNVERIFIED(g=0 也上升)" if not g0_clean else ("✅ 该版有效" if fired else "⛔ 该版失明"))
        rows.append(dict(rule=rule, page=P, c0=c0, c_bad=c_bad, c_good=c_good,
                         fired=bool(fired), g0_clean=bool(g0_clean), verdict=v))
        print(f"  {rule:22s} {P:14s} 基线 {c0:2d} · 注入缺陷 {c_bad:2d} · 注入无害 {c_good:2d}  -> **{v}**")
for rule, why in UNCOMPUTED.items():
    rows.append(dict(rule=rule, page="—", verdict="UNCOMPUTED", why=why))
    print(f"  {rule:22s} {'—':14s} **UNCOMPUTED** —— {why}")
T = pd.DataFrame(rows)
M = T[T.verdict.isin(["✅ 该版有效", "⛔ 该版失明"])]
print(f"\n**可测格 {len(M)} · 有效 {int((M.verdict=='✅ 该版有效').sum())} · "
      f"失明 **{int((M.verdict=='⛔ 该版失明').sum())}** · UNCOMPUTED {int((T.verdict=='UNCOMPUTED').sum())}**")

G = Gate("每条规则,在每一版页面上开得了火吗?")
g0_all = bool(T.g0_clean.dropna().all()) if "g0_clean" in T else False
print(f"\n  g=0:所有可测格注入「不该抓的东西」后计数都不上升?**{g0_all}**")
pos_ok = G.positive_control("正对照:至少一格在注入缺陷后上升",
                            planted=float(M.fired.sum()) if len(M) else 0.0, floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:注入无害内容不得使计数上升",
                            # ⚠ 第一版这里写  -> object dtype 上 ,
                            #   于是 6 个干净格被算成 -12,**对照因它自己的算术而失败**(realstat 点名的那一行)。
                            null=float((T.g0_clean.dropna() == False).sum()),
                            effect=float(M.fired.sum()) if len(M) else 1.0,
                            null_spread=0.4, null_kind="注入一个该规则不该抓的东西")
blind = M[M.verdict == "⛔ 该版失明"]
if pos_ok and pla_ok:
    verdict = (f"可测 {len(M)} 格:有效 {int((M.verdict=='✅ 该版有效').sum())} · 失明 {len(blind)}"
               + (f" ⇒ 失明的格写进页面「做不到什么」:{list(zip(blind.rule, blind.page))}" if len(blind)
                  else " ⇒ 无失明格"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(table=T.to_dict("records"), n_testable=len(M),
               n_effective=int((M.verdict=="✅ 该版有效").sum()), n_blind=len(blind),
               n_uncomputed=int((T.verdict=="UNCOMPUTED").sum()), verdict=verdict, unchallenged=True),
          open(OUT/"can_it_fire.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'can_it_fire.json'}")
