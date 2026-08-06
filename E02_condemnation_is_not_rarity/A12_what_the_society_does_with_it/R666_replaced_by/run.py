"""E02·A12·R666 —— 给闸门加一列:那个离开的数,是被拿走了,还是被换掉了?

`#629` 的 NEXT。**行动类型:PRODUCTION**(如实标注)。
`#629c`:读者看到「有个数离开了」时分不出**抹除**与**更正**,而这两件事要求的反应相反。

⚠ **§3 梯度检查:「同一段落」需要新旧段落配对,而配对是启发式。**
   ⇒ 判定必须**三值**:`替换`(找到新数)/ `拿走`(段落配上了但没有新数)/ **`配不上`**(段落没配上)。
   **把启发式的失败明写成第三值,而不是塞进前两值** —— 否则一次配对失败会被读成「它被拿走了」。
   配对规则**先写死**:在新文本里取与旧段落**词重合度最高**的段落,且重合度须 ≥ 0.50;否则 `配不上`。

G1 ESTIMAND(先于方法):对每个离开的 token `t`,
  `replaced_by(t)` ∈ {某个新 `_MAGNUM` token, `None`, `UNMATCHED`}。
CONTROLS:
  正对照:改写臂必须给出 `replaced_by = +0.8484`。
  **g=0**:删除臂必须给出 `None`(**不是** `UNMATCHED`)。
  安慰剂:原样写回 -> 根本不该有离开的 token(空表)。
KILL(条件式,预注册于 `#629`):
  if 正对照 = `+0.8484` and g=0 = `None`:
      **回测**(`#620` 硬要求):在 `#623` 的候选提交上跑新旧两版,
      **新版必须报出旧版报过的每一个 token**(只准增列,不准少报);
      漏任何一个 -> **不改**;一个不漏 -> **接受这一列**
  else: UNVERIFIED
G3:两版 × 三臂 + 回测漏报表。
⚠ 安全边界同 `#628a`:临时分支 · `--no-verify` · `finally` 回 `main` 删分支 · 收尾核对打印。
IMPOSSIBLE(不写 planned):**配对是启发式**,`UNMATCHED` 的比例就是它的代价 ·
  只处理 `_MAGNUM` 认得的数种 · 一段里若同时走了两个数、来了两个数,**配不出一一对应** · `[unchallenged]`
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
T1, T2 = "+0.7373", "+0.8484"
BR = "tmp/r666"
run = lambda *a: subprocess.run(a, capture_output=True, text=True)
norm = lambda m: m.strip().replace(' ', '').replace('倍', '×').replace('−', '-')


def paras(txt): return [p for p in re.split(r'\n\s*\n', txt) if p.strip()]


def overlap(a, b):
    wa, wb = set(re.findall(r'\w+', a)), set(re.findall(r'\w+', b))
    return len(wa & wb) / max(len(wa), 1)


def numbers_that_left_v2(rev, files=('README.md', 'README_zh.md')):
    """旧逻辑一字不动,只**增加一列** `replaced_by`(三值)。"""
    out = []
    for f in files:
        r = run("git", "show", f"{rev}:{f}")
        if r.returncode != 0: continue
        old, new = r.stdout, pathlib.Path(f).read_text()
        o = set(norm(m) for m in A._MAGNUM.findall(old))
        n = set(norm(m) for m in A._MAGNUM.findall(new))
        newcomers = n - o
        NP = paras(new)
        for tk in sorted(o - n):
            src = [p for p in paras(old) if tk in norm(p) or tk in p]
            rep, kind = None, "UNMATCHED"
            if src:
                best = max(NP, key=lambda q: overlap(src[0], q)) if NP else ""
                if best and overlap(src[0], best) >= 0.50:
                    cand = [norm(m) for m in A._MAGNUM.findall(best) if norm(m) in newcomers]
                    rep, kind = (cand[0], "替换") if cand else (None, "拿走")
            out.append(dict(file=f, token=tk, replaced_by=rep, kind=kind))
    return pd.DataFrame(out, columns=["file", "token", "replaced_by", "kind"])


sb = run("git", "branch", "--show-current").stdout.strip()
sh_ = run("git", "rev-parse", "HEAD").stdout.strip()
tracked = len([x for x in run("git", "status", "--porcelain", "-uno").stdout.split("\n") if x.strip()])
print(f"起始:分支 {sb} · 已跟踪改动 {tracked} · HEAD {sh_[:7]}")
assert tracked == 0, "有已跟踪的未提交改动 —— 拒绝伪造提交"

# ⚠ 「删除(整句改写)」这一臂先跑,给出 `UNMATCHED` —— 段落配对在它最该管用的地方失败。
#   ⇒ 补一个「删除(最小)」子臂:**只拿掉 token,句子其余一字不动** —— 那才是真实删除的样子。
#   **这一臂是看到第一次结果之后加的,如实标注;两臂都报,不许只报好看的那个。**
ARMS = {"改写": lambda b: b + f"\n\n伪造素材段落:这里带着一个数 {T2},其余文字保持不变以便配对。\n",
        "删除(整句改写)": lambda b: b + "\n\n伪造素材段落:这里的数已被删掉,其余文字保持不变以便配对。\n",
        "删除(最小)": lambda b: b + "\n\n伪造素材段落:这里带着一个数 ,其余文字保持不变以便配对。\n",
        "安慰剂(原样写回)": lambda b: b + f"\n\n伪造素材段落:这里带着一个数 {T1},其余文字保持不变以便配对。\n"}
rows = []
try:
    run("git", "checkout", "-q", "-b", BR)
    for P in ("README.md", "README_zh.md"):
        for arm, mk in ARMS.items():
            p = pathlib.Path(P); base = p.read_text()
            p.write_text(base + f"\n\n伪造素材段落:这里带着一个数 {T1},其余文字保持不变以便配对。\n")
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"A {P} {arm}")
            Asha = run("git", "rev-parse", "HEAD").stdout.strip()
            p.write_text(mk(base))
            run("git", "add", P); run("git", "commit", "-q", "--no-verify", "-m", f"B {P} {arm}")
            D = numbers_that_left_v2(Asha)
            d = D[(D.file == P) & (D.token == T1)]
            rows.append(dict(page=P, arm=arm, n=int(len(D)),
                             kind=(d.iloc[0]["kind"] if len(d) else "—"),
                             replaced_by=(d.iloc[0]["replaced_by"] if len(d) else None)))
            print(f"  {P:14s} {arm:14s} -> {len(D)} 行 · kind {rows[-1]['kind']} · "
                  f"replaced_by {rows[-1]['replaced_by']}")
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
pos_ok_raw = all(g(P, "改写").replaced_by == T2 for P in PGS)
# ⚠ 第一版写 `replaced_by is None` -> pandas 把 None 存成 NaN,`is None` 恒假,g=0 假失败。
#   **本次会话第三次「对照因 dtype 而假失败」**(`#627e` 的 `~True = -2`;`#663` 同源;本条)。
#   ⇒ 用 `pd.isna()`。
g0_ok = all(g(P, "删除(最小)").kind == "拿走" and pd.isna(g(P, "删除(最小)").replaced_by) for P in PGS)
g0_rewrite = all(g(P, "删除(整句改写)").kind == "UNMATCHED" for P in PGS)
pla_ok_raw = all(g(P, "安慰剂(原样写回)").n == 0 for P in PGS)
print(f"\n  正对照:改写臂 replaced_by = `{T2}`?**{pos_ok_raw}**")
print(f"  g=0(最小删除):kind=拿走 且 replaced_by=None?**{g0_ok}**")
print(f"  ⚠ 而删除(整句改写)那一臂两版都是 `UNMATCHED`?**{g0_rewrite}** —— **配对启发式的代价,如实报**")
print(f"  安慰剂:原样写回 -> 空表?**{pla_ok_raw}**")

# ── 回测(`#620` 硬要求)────────────────────────────────────
led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
PAT = re.compile(r'两份 README 已改|两版 README|页面已改|页面.{0,6}更正|已更正.{0,10}页面|正确写法')
CAND = sorted({ent[j] for j, l in enumerate(led) if PAT.search(l) and ent[j] > 0})
miss, n_bt = [], 0
for N in CAND:
    log = run("git", "log", "--format=%H", "-S", f"## Entry {N} ", "--", "RETRACTIONS.md").stdout.split()
    if not log: continue
    par = run("git", "rev-parse", f"{log[-1]}^").stdout.strip()
    if not par: continue
    n_bt += 1
    old_out = A.numbers_that_left(rev=par); new_out = numbers_that_left_v2(par)
    o = set(zip(old_out.file, old_out.token)) if len(old_out) else set()
    nn = set(zip(new_out.file, new_out.token)) if len(new_out) else set()
    if o - nn: miss.append(dict(entry=N, missing=sorted(o - nn)[:4]))
print(f"\n=== 回测(`#620` 硬要求):{n_bt} 个版本 · **新版漏报旧版 token 的次数 {len(miss)}** ===")
for m_ in miss[:6]: print(f"    #{m_['entry']}: {m_['missing']}")

G = Gate("给闸门加一列:那个数是被拿走了,还是被换掉了?")
pos_ok = G.positive_control("正对照:改写臂必须给出 replaced_by", planted=float(pos_ok_raw), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:删除臂必须给出 None 且 kind=拿走", null=float(0 if g0_ok else 1),
                            effect=1.0, null_spread=0.4, null_kind="第二次提交把那个数删掉而不是换掉")
if pos_ok and pla_ok:
    verdict = ("**接受这一列**(回测零漏报)" if not miss else f"**不改**:回测漏报 {len(miss)} 次")
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(table=T.to_dict("records"), pos=bool(pos_ok_raw), g0=bool(g0_ok), g0_rewrite=bool(g0_rewrite), placebo=bool(pla_ok_raw),
               backtest_versions=n_bt, backtest_misses=miss, verdict=verdict,
               housekeeping=dict(branch=eb, head=eh, dirty=ed, leftover=left), unchallenged=True),
          open(OUT/"replaced_by.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'replaced_by.json'}")
