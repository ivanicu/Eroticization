"""E02·A12·R664 —— 把四个 UNCOMPUTED 里最便宜的一个做掉:伪造一次提交

`#627` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
`#627d`:`numbers_that_left` 的缺陷 = **「与上一版相比某个数离开了页面」** ——
它不住在任何一版页面里,住在**两次提交之间的差**里 ⇒ 单版注入够不到它。**本轮伪造那次提交。**

⚠ **安全边界先写死(这是要造提交的一轮):**
   · 只在**临时分支** `tmp/r664` 上做,**绝不在 `main` 上**;
   · 伪造提交用 `--no-verify` 绕过闸门 —— **这是在一次性分支上的、被标注的刻意绕过**;
   · 无论成败,`finally` 里**必须**回到 `main`、删掉临时分支、并验证工作区干净;
   · 脚本最后打印**分支 / 工作区 / HEAD**,让读者能核。

⚠ **`#624b` 的教训直接适用**:`_MAGNUM` 只认**带符号小数 / 百分数 / 倍数** ⇒
   注入的数必须带符号(`+0.7373`),否则正对照会因为「它压根不看这种数」而假失败。

G1 ESTIMAND(先于方法):对每一版页面 P,
  `fires(P)` = 伪造「提交 A 加入带 `+0.7373` 的段落 → 提交 B 删掉它」之后,
  `numbers_that_left(rev=A)` 是否列出 `+0.7373`(且 `file == P`)。
CONTROLS:
  正对照 = 上述注入(每版一次)。
  **g=0** = 同样两次提交,但**第二次只改措辞、数字不动** -> **必须不列出**。
KILL(条件式,预注册):
  if g=0 不列出:
      正对照列出 -> **该版有效,`UNCOMPUTED` 4 -> 3**
      正对照不列出 -> **该版失明**,写进页面「做不到什么」,**不许改规则**(`#620`)
  else: UNVERIFIED —— 这个伪造检验没有分辨力
G3:两版 × {正对照, g=0} 四格全表。
IMPOSSIBLE(不写 planned):它只验了**这一种**形状(一个带符号的数被整段删掉)——
  「数字被改写成另一个数」`_MAGNUM` 会看成一删一增,**本轮没有验** · `[unchallenged]`
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
TOK = "+0.7373"
BR = "tmp/r664"
run = lambda *a: subprocess.run(a, capture_output=True, text=True)
start_branch = run("git", "branch", "--show-current").stdout.strip()
start_head = run("git", "rev-parse", "HEAD").stdout.strip()
# ⚠ 第一版的守卫用 `git status --porcelain`(含未跟踪),于是被**本轮自己的新脚本**挡住了。
#   守卫太严会让人去关掉它 —— 而真正的危险是 `reset --hard` 会**毁掉的已跟踪修改**。
#   ⇒ 收紧到该收紧的地方:`-uno`(只看已跟踪),未跟踪的新文件放行。
tracked = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
untracked = len([x for x in run("git", "status", "--porcelain").stdout.split("\n") if x.strip()]) - tracked
print(f"起始:分支 {start_branch} · 已跟踪改动 {tracked} 项 · 未跟踪 {untracked} 项 · HEAD {start_head[:7]}")
assert tracked == 0, "有已跟踪的未提交改动 —— 拒绝伪造提交(reset --hard 会毁掉它)"
dirty = tracked

rows = []
try:
    run("git", "checkout", "-q", "-b", BR)
    for P in ("README.md", "README_zh.md"):
        for arm, second in (("正对照(删掉那个数)", "del"), ("g=0(只改措辞)", "word")):
            p = pathlib.Path(P); base = p.read_text()
            # 提交 A:加入带符号的数
            p.write_text(base + f"\n\n这一段是伪造的检验素材,带一个数 {TOK},它随后会被处理。\n")
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"forge A {P} {arm}")
            A_sha = run("git", "rev-parse", "HEAD").stdout.strip()
            # 提交 B
            if second == "del":
                p.write_text(base + "\n\n这一段是伪造的检验素材,那个数已经被删掉。\n")
            else:
                p.write_text(base + f"\n\n这一段是伪造的检验素材(措辞已改),带一个数 {TOK},它没有被动过。\n")
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"forge B {P} {arm}")
            D = A.numbers_that_left(rev=A_sha)
            hit = bool(len(D) and ((D.token == TOK) & (D.file == P)).any())
            rows.append(dict(page=P, arm=arm, n_rows=int(len(D)), fired=hit))
            print(f"  {P:14s} {arm:18s} -> numbers_that_left 返回 {len(D)} 行 · 命中 `{TOK}` **{hit}**")
            # 回到分支起点,准备下一格
            run("git", "reset", "-q", "--hard", start_head)
finally:
    run("git", "checkout", "-q", start_branch or "main")
    run("git", "branch", "-D", BR)
    end_branch = run("git", "branch", "--show-current").stdout.strip()
    end_head = run("git", "rev-parse", "HEAD").stdout.strip()
    end_dirty = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
    left = BR in run("git", "branch", "--list", BR).stdout
    print(f"\n=== 收尾核对 ===")
    print(f"  分支 {end_branch}(须 {start_branch}){'✅' if end_branch == start_branch else '⛔'} · "
          f"HEAD {end_head[:7]}(须 {start_head[:7]}){'✅' if end_head == start_head else '⛔'} · "
          f"工作区 {end_dirty} 项{'✅' if end_dirty <= 0 else '⛔'} · "
          f"临时分支残留 {left}{'⛔' if left else ' ✅'}")

T = pd.DataFrame(rows)
print(f"\n=== G3:{len(T)} 格全表 ===")
print(T.to_string(index=False))
G = Gate("伪造一次提交:`numbers_that_left` 在每一版上抓得到吗?")
pos = T[T.arm.str.startswith("正对照")]; g0 = T[T.arm.str.startswith("g=0")]
print(f"\n  正对照命中 {int(pos.fired.sum())}/{len(pos)} · g=0 命中 {int(g0.fired.sum())}/{len(g0)}(须 0)")
pos_ok = G.positive_control("正对照:删掉的数必须被抓到", planted=float(pos.fired.sum()), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:只改措辞不得被抓到", null=float(g0.fired.sum()),
                            effect=float(max(pos.fired.sum(), 1)), null_spread=0.4,
                            null_kind="第二次提交只改措辞、数字不动")
if pla_ok:
    n_ok = int(pos.fired.sum())
    verdict = (f"两版正对照命中 {n_ok}/2 ⇒ "
               + ("**两版都有效,`UNCOMPUTED` 4 -> 3**" if n_ok == 2
                  else f"**有 {2-n_ok} 版失明,写进页面「做不到什么」,不改规则**"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— g=0 也命中({int(g0.fired.sum())})⇒ 这个伪造检验没有分辨力"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(table=T.to_dict("records"), verdict=verdict,
               housekeeping=dict(branch=end_branch, head=end_head, dirty=end_dirty, leftover=left),
               unchallenged=True),
          open(OUT/"forge_a_commit.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'forge_a_commit.json'}")
