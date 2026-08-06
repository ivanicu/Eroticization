"""E01·A204·R560 — 只审最承重的那一小节,顺着账本的交叉引用链走

`#515` 的 NEXT,预注册的分支已在 `#515d` 执行:77 条降级为背景,
**本轮只审页面「这说明了关于人的什么 / What this says about people」小节**。

⛔ 两条从 `#515` 带过来的硬约束:
  ① **排除 `A203_` 一族的账本条目**(`#515c`:账本已开始引用页面,独立性在流失);
  ② **不能只看相邻文本**(`#515b`:页面「标注而非重写」,限定可以离那个数任意远)
     -> 必须顺着账本里 `#NNN` 的**交叉引用**走。

G1 ESTIMAND:该小节里每一个含数的断言,给出
  ① **来源条目号**(哪些 `## Entry N` 的正文含这个数);
  ② **后续引用**(有没有编号**更大**的条目在正文里引用了那个来源号);
  ③ 页面那一句**是否已带**该后续引用所说的修改(**这一条本轮不自动判,只把两段摊开**)。
判据(预注册):
  **有后续引用且页面未带** -> 逐条修;**全部已带或无后续引用** -> 该小节无过期引用。
CONTROLS:正对照 = 已知有引用链的一对(`#243` -> `#245c` 改写了它的表述)必须被链路发现;
  阴性 = 随机取同样多的**非该小节**的数,其「有后续引用」的比例作基线。
IMPOSSIBLE:「页面是否已带」需人读 ⇒ 本轮只产出**可复核的并排材料**,不是一个可失败的门 ·
  未派对抗 agent ⇒ [unchallenged]
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
led = pathlib.Path("RETRACTIONS.md").read_text()
# 切成条目
ENT = re.compile(r"^## Entry (\d+)", re.M)
marks = [(int(m.group(1)), m.start()) for m in ENT.finditer(led)]
bodies = {}
for i, (n, s) in enumerate(marks):
    e = marks[i + 1][1] if i + 1 < len(marks) else len(led)
    bodies[n] = led[s:e]
# ⛔ 排除 A203_ 一族(本轮的审计条目自身)
EXCL = {n for n, b in bodies.items() if "A203_" in b or "A204_" in b}
print(f"账本条目 {len(bodies)} 条;排除 A203/A204 自身 {len(EXCL)} 条 -> 可用 {len(bodies)-len(EXCL)}")

page = pathlib.Path("README.md").read_text()
i = page.find("# What this says about people")
j = page.find("\n# ", i + 5)
sec = page[i:j if j > 0 else len(page)]
NUM = re.compile(r"[−+-]?\d+\.\d+")
toks = sorted(set(NUM.findall(sec)))
print(f"该小节长度 {len(sec):,} 字;含小数的数 {len(toks)} 个")

def forms(t): return {t, t.replace("−", "-"), t.replace("-", "−")}
def src_entries(t):
    return sorted(n for n, b in bodies.items() if n not in EXCL and any(f in b for f in forms(t)))

rows = []
for t in toks:
    src = src_entries(t)
    if not src: rows.append(dict(tok=t, src=[], later=[])); continue
    later = []
    for s in src:
        pat = re.compile(rf"#{s}[a-z]?\b")
        later += [n for n, b in bodies.items()
                  if n > s and n not in EXCL and pat.search(b)]
    rows.append(dict(tok=t, src=src, later=sorted(set(later))))
have_later = [r for r in rows if r["later"]]
no_src = [r for r in rows if not r["src"]]
print(f"\n有来源条目的:{len(rows)-len(no_src)} · 无来源:{len(no_src)} · **有后续引用的:{len(have_later)}**")

# 正对照:#243 -> #245 的链路必须被发现
pc = any(243 in r["src"] and any(n >= 245 for n in r["later"]) for r in rows)
print(f"正对照(`#243` 的后续 `#245` 链路)-> {'发现 ✅' if pc else '未发现 ⛔'}")

# 阴性:非该小节的数的「有后续引用」比例
other = sorted(set(NUM.findall(page)) - set(toks))
rng = np.random.default_rng(20260805)
samp = [other[k] for k in rng.choice(len(other), min(60, len(other)), replace=False)]
base = np.mean([1.0 if any(n > s for s in src_entries(x) for n in bodies
                           if n > s and n not in EXCL and re.search(rf"#{s}[a-z]?\b", bodies[n]))
                else 0.0 for x in samp])
print(f"阴性基线(非该小节 {len(samp)} 个数的「有后续引用」比例)= {base:.3f}")
rate = len(have_later) / max(len(rows) - len(no_src), 1)
print(f"该小节的比例 = {rate:.3f}")

print("\n" + "=" * 74)
print("有后续引用的条目 —— 来源 -> 后续(供逐条人读)")
print("=" * 74)
for r in have_later[:14]:
    print(f"  {r['tok']:>10s}  来源 {r['src']}  -> 后续 {r['later'][:6]}")

G = Gate("只审最承重的那一小节,顺着交叉引用链走")
G.asserted("正对照:已知链路 `#243`→`#245` 必须被发现", pc, "见上", kind="control")
G.negative_control("阴性:非该小节的数的后续引用比例", null=float(base), effect=rate,
                   null_spread=0.05, null_kind="页面其余部分的同一统计量")
print(G)
json.dump(dict(n_entries=len(bodies), excluded=sorted(EXCL), sec_tokens=toks,
               rows=rows, n_have_later=len(have_later), n_no_src=len(no_src),
               rate=rate, base=float(base), pc=bool(pc), unchallenged=True),
          open(OUT / "claim_chain_audit.json", "w"), indent=1)
print(f"\nwrote {OUT/'claim_chain_audit.json'}")
