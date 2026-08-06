"""E02·A12·R667 —— 一段里同时走两个、来两个:它会不会编一个一一对应出来?

`#630` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
⚠ **我在写这一轮之前就知道答案**:`#630` 接进去的代码是 `cand[0]` ——
**它会挑第一个候选而不声明还有别的。那就是编造。**
**而我仍然先量,再修** —— 因为「我以为我知道」和「我量过」是两件事,
而这个项目已经有 `#620`/`#624` 两次记录:我以为的判据,不如仓库里躺着的代码。

G1 ESTIMAND(先于方法):在**同一段**里同时移除两个量、引入两个量之后,
  新列 `replaced_by` 给出的是什么:
  **唯一值(编造)** / **全部候选(诚实)** / **`UNMATCHED`(退回)**。
KILL(条件式,预注册于 `#630`):
  if g=0(一走一来)给出唯一且正确的 `replaced_by`:
      两走两来时若给出**任意一种一一对应而不声明还有别的** -> **记为缺陷并修**
      若列出全部候选或退回 `UNMATCHED` -> **合格**
  else: UNVERIFIED
修法(先写死,免得看到结果再设计):新增第四值 **`多候选`**,并把**全部**候选写进 `replaced_by`。
**回测**(`#620` 硬要求):修完必须在 `#623` 的候选提交上重放,**漏报任何旧 token ⇒ 不改**。
⚠ 安全边界同 `#628a`;收尾核对打印。
IMPOSSIBLE(不写 planned):它仍只处理 `_MAGNUM` 认得的数种 ·
  **「哪一个换了哪一个」在没有语义的前提下不可判定** —— 本轮的修法是**如实列出全部**,不是解决它 ·
  `[unchallenged]`
"""
import os, sys, pathlib, json, re, subprocess, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
sys.path.insert(0, str(ROOT/"tools"))
import readme_ledger_audit as A
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
BR = "tmp/r667"
run = lambda *a: subprocess.run(a, capture_output=True, text=True)
sb = run("git", "branch", "--show-current").stdout.strip()
sh_ = run("git", "rev-parse", "HEAD").stdout.strip()
tracked = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
print(f"起始:分支 {sb} · 已跟踪改动 {tracked} · HEAD {sh_[:7]}")
assert tracked == 0, "有已跟踪的未提交改动 —— 拒绝伪造提交"

SEED = "\n\n伪造素材段落:这里带着 {a} 与 {b},其余文字保持不变以便配对。\n"
ARMS = {
 "两走两来":  (("+0.7373", "+0.5151"), ("+0.8484", "+0.6262")),
 "g=0(一走一来)": (("+0.7373", "+0.5151"), ("+0.8484", "+0.5151")),
}
rows = []
try:
    run("git", "checkout", "-q", "-b", BR)
    for P in ("README.md", "README_zh.md"):
        for arm, (old2, new2) in ARMS.items():
            p = pathlib.Path(P); base = p.read_text()
            p.write_text(base + SEED.format(a=old2[0], b=old2[1]))
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"A {P} {arm}")
            Asha = run("git", "rev-parse", "HEAD").stdout.strip()
            p.write_text(base + SEED.format(a=new2[0], b=new2[1]))
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"B {P} {arm}")
            D = A.numbers_that_left(rev=Asha)
            d = D[D.file == P] if len(D) else D
            rows.append(dict(page=P, arm=arm, n=int(len(d)),
                             detail=[(r.token, r.kind, r.replaced_by) for r in d.itertuples()]))
            print(f"  {P:14s} {arm:14s} -> {len(d)} 行 · {rows[-1]['detail']}")
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
print(f"\n=== G3:{len(T)} 格 ===\n{T.to_string(index=False)}")
g = lambda P, a: T[(T.page == P) & (T.arm == a)].iloc[0]
PGS = ("README.md", "README_zh.md")
# 两走两来:是否给出了「唯一值」而不声明还有别的
fabricates = []
for P in PGS:
    for tok, kind, rep in g(P, "两走两来").detail:
        if kind == "替换" and isinstance(rep, str) and "," not in rep and "|" not in rep:
            fabricates.append((P, tok, rep))
g0_ok = all(any(k == "替换" and r == "+0.8484" for _, k, r in g(P, "g=0(一走一来)").detail) for P in PGS)
print(f"\n  g=0(一走一来):给出唯一且正确的 `+0.8484`?**{g0_ok}**")
print(f"  两走两来:给出**唯一值而不声明还有别的**的格数 = **{len(fabricates)}** -> {fabricates}")

G = Gate("一段里同时走两个、来两个:它会不会编一个一一对应出来?")
pos_ok = G.positive_control("正对照:g=0 必须给出唯一且正确的替换", planted=float(g0_ok), floor=0.0, spread=0.4)
verdict = ("**缺陷确认:两走两来时它挑了第一个候选,没有声明还有别的 ⇒ 必须修**" if fabricates
           else "**合格:它列出了全部候选或退回 UNMATCHED**")
print(f"\n{'控制齐备 ⇒ ' if pos_ok else '⚠ '}判定。**{verdict}**")
print(G)
json.dump(dict(table=T.to_dict("records"), fabricates=fabricates, g0_ok=bool(g0_ok),
               verdict=verdict, housekeeping=dict(branch=eb, head=eh, dirty=ed, leftover=left),
               unchallenged=True),
          open(OUT/"two_leave_two_arrive.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'two_leave_two_arrive.json'}")
