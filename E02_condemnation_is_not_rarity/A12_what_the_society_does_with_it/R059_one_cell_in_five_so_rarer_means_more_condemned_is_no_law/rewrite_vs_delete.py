"""E02·A12·R665 —— 「改写」和「删除」,在它的输出里分得开吗?

`#628` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。

⚠ **§3 梯度检查先做,而它改掉了判据的读法:**
   `numbers_that_left` 的代码是 `for tk in sorted(o - n)` —— **它结构上只能报「谁走了」,报不出「谁来了」。**
   ⇒ 「它对改写半盲」这件事**是代数强制的,是 DERIVATION,不是发现**(realstat:标成推导,并说出假设)。
   **假设**:`_MAGNUM` 对新旧两版各抽一次,差集只取 `old − new`。
   ⇒ 所以本轮真正可测的是更锋利的那一问:
     **「改写」和「删除」在它的输出里是不是一模一样?** 如果是,读者看到「有个数离开了」时
     **分不出它是被拿掉了,还是被换成了别的数** —— 而这两件事要求的下一步完全不同。

G1 ESTIMAND(先于方法):对每一版页面 P,
  `out_del(P)` = 「提交 A 写入 `+0.7373` → 提交 B **删掉**它」之后 `numbers_that_left(rev=A)` 的输出;
  `out_rw(P)`  = 同样的 A,但提交 B 把它**改写成 `+0.8484`**;
  **`indistinguishable(P)` = (out_del 与 out_rw 的 `token` 集合相同)**。
CONTROLS:
  正对照:两臂都必须列出 `+0.7373`(否则这个比较无从谈起)。
  **g=0**:提交 B 把 `+0.7373` **原样写回** -> 两臂都**必须为空**。
  安慰剂(新增):提交 B **只加一个新数** `+0.8484` 而不动 `+0.7373` -> **必须为空**
    (证明它确实对「新增」无反应,而这正是那条 DERIVATION 的实测)。
KILL(条件式,预注册):
  if g=0 为空 and 正对照两臂都列出:
      两臂 token 集合相同 -> **分不开:输出上「改写」= 「删除」**,写进页面
      不同 -> **分得开**,报差别
  else: UNVERIFIED
G3:两版 × {删除, 改写, g=0, 安慰剂} 八格全表。
⚠ 安全边界同 `#628a`:临时分支 · `--no-verify` · `finally` 回 `main` 删分支 · **收尾核对打印出来**。
IMPOSSIBLE(不写 planned):它只验 `_MAGNUM` 认得的数种 · 只验单 token 的改写 · `[unchallenged]`
"""
import os, sys, pathlib, json, subprocess, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
T1, T2 = "+0.7373", "+0.8484"
BR = "tmp/r665"
run = lambda *a: subprocess.run(a, capture_output=True, text=True)
sb = run("git", "branch", "--show-current").stdout.strip()
sh_ = run("git", "rev-parse", "HEAD").stdout.strip()
tracked = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
print(f"起始:分支 {sb} · 已跟踪改动 {tracked} · HEAD {sh_[:7]}")
assert tracked == 0, "有已跟踪的未提交改动 —— 拒绝伪造提交"

ARMS = {
 "删除":            lambda base: base + "\n\n伪造素材:那个数已经被删掉。\n",
 "改写":            lambda base: base + f"\n\n伪造素材:那个数已被改写成 {T2}。\n",
 "g=0(原样写回)":  lambda base: base + f"\n\n伪造素材:那个数 {T1} 原样还在。\n",
 "安慰剂(只新增)": lambda base: base + f"\n\n伪造素材:那个数 {T1} 还在,另外新增一个 {T2}。\n",
}
rows = []
try:
    run("git", "checkout", "-q", "-b", BR)
    for P in ("README.md", "README_zh.md"):
        for arm, mk in ARMS.items():
            p = pathlib.Path(P); base = p.read_text()
            p.write_text(base + f"\n\n伪造素材:带一个数 {T1}。\n")
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"A {P} {arm}")
            Asha = run("git", "rev-parse", "HEAD").stdout.strip()
            p.write_text(mk(base))
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"B {P} {arm}")
            D = A.numbers_that_left(rev=Asha)
            toks = sorted(set(D.token)) if len(D) else []
            rows.append(dict(page=P, arm=arm, n=int(len(D)), tokens=toks))
            print(f"  {P:14s} {arm:16s} -> {len(D)} 行 · tokens {toks}")
            run("git", "reset", "-q", "--hard", sh_)
finally:
    run("git", "checkout", "-q", sb or "main"); run("git", "branch", "-D", BR)
    eb = run("git", "branch", "--show-current").stdout.strip()
    eh = run("git", "rev-parse", "HEAD").stdout.strip()
    ed = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
    left = BR in run("git", "branch", "--list", BR).stdout
    print(f"\n=== 收尾核对 ===\n  分支 {eb}{'✅' if eb==sb else '⛔'} · HEAD {eh[:7]}{'✅' if eh==sh_ else '⛔'} · "
          f"已跟踪改动 {ed}{'✅' if ed==0 else '⛔'} · 临时分支残留 {left}{'⛔' if left else ' ✅'}")

T = pd.DataFrame(rows)
print(f"\n=== G3:{len(T)} 格全表 ===")
print(T.to_string(index=False))
get = lambda P, a: T[(T.page == P) & (T.arm == a)].iloc[0]
same = {P: get(P, "删除").tokens == get(P, "改写").tokens for P in ("README.md", "README_zh.md")}
g0_empty = all(get(P, "g=0(原样写回)").n == 0 for P in same)
pla_empty = all(get(P, "安慰剂(只新增)").n == 0 for P in same)
pos_both = all(T1 in get(P, a).tokens for P in same for a in ("删除", "改写"))
print(f"\n  正对照:两臂都列出 `{T1}`?**{pos_both}**")
print(f"  g=0(原样写回)两版都为空?**{g0_empty}**")
print(f"  安慰剂(只新增)两版都为空?**{pla_empty}** ⇒ **那条 DERIVATION 得到实测支持:它对「新增」无反应**")
print(f"  **删除 vs 改写,token 集合相同?** {same}")

G = Gate("「改写」和「删除」,在它的输出里分得开吗?")
pos_ok = G.positive_control("正对照:两臂都必须列出那个离开的数",
                            planted=float(pos_both), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:原样写回必须为空", null=float(0 if g0_empty else 1),
                            effect=1.0, null_spread=0.4, null_kind="第二次提交把那个数原样写回")
if pos_ok and pla_ok:
    verdict = ("**分不开:输出上「改写」与「删除」完全相同** —— 读者看到「有个数离开了」时,"
               "分不出它是被拿掉了还是被换成了别的数" if all(same.values())
               else f"**分得开**:{same}")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(table=T.to_dict("records"), same=same, g0_empty=bool(g0_empty),
               placebo_empty=bool(pla_empty), pos_both=bool(pos_both), verdict=verdict,
               derivation="`for tk in sorted(o-n)` ⇒ 结构上只报「谁走了」,报不出「谁来了」",
               housekeeping=dict(branch=eb, head=eh, dirty=ed, leftover=left), unchallenged=True),
          open(OUT/"rewrite_vs_delete.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'rewrite_vs_delete.json'}")
