"""E02·A12·R659 —— 回测:这条规则在历史上抓到过什么?

`#622` 的 NEXT。**行动类型:CLOSURE**(如实标注)。
`#620` 立的规矩:**改规则必须先用 git 历史回测**。本轮把它跑掉。

⚠ **§3 梯度检查(先做,它决定判据怎么写)**:`#622` 已证明现规则**当前**的四格全是伪影 ——
   若把判据写成「候选必须抓到现规则抓到的每一个真缺陷」,而现规则一个也没抓到,
   **那就是一个空真的判据,一个不可能失败的检查。**
   ⇒ 所以「真缺陷」的定义**取自账本,不取自现规则**(这一条在 `#622` 的 NEXT 里已写死)。

G1 ESTIMAND(先于方法):对每一个「账本记过页面被更正」的条目 `#N`:
  · **缺陷版本** = 首次引入 `## Entry N` 的那次提交的**父提交**;
  · **错的数字** = 该次提交从 `README.md` 里**删掉**的小数(diff 的 `-` 行减去 `+` 行);
  · **「抓到」的定义(先写死)**:规则在**父提交**上标记出的某一格,其**两侧差集**里
    至少含有一个「错的数字」。
  三条规则各算一次:**现规则** · **候选①跨文件**(比较两版页面同一标记的数字集合)·
  **候选②一侧空集不参与比较**。

KILL(条件式,预注册于 `#622`):
  if 缺陷版本集合非空 and 至少一条规则抓到过 >0 个(否则整个回测没有分辨力):
      候选漏掉现规则抓到的任何一个 -> **不换规则**
      候选抓到的是现规则的超集 -> 允许换,并报两者的命中/漏报
  else: **UNVERIFIED —— 回测没有分辨力**,不据它换规则
G3:36 个缺陷版本 × 3 条规则的完整命中表。G4:候选②的空集判定 {严格空 / ≤1 个数} 两档。
IMPOSSIBLE(不写 planned):**账本记「页面已改」的条目里,很多是措辞更正,不是数字更正** ——
  这些条目**结构上不可能被任何数字规则抓到**,必须单独计数,不能算作漏报 ·
  条目 -> 提交的映射靠 `git log -S`,**它是一个搜索**,可能失配 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, subprocess, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
NUM = re.compile(r'(?<![\w.])\d+\.\d{2,4}(?![\w])')
CJK = re.compile(r'[一-鿿]')
CITE = re.compile(r'#\d{1,3}\b|\bA\d{2}\b|Entry \d{1,3}')


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True).stdout


def rule_current(en, zh):
    cites = {}
    for i, l in enumerate(en.splitlines(), 1):
        for c in set(CITE.findall(l)):
            cites.setdefault(c, {True: set(), False: set()})[bool(CJK.search(l))] |= set(NUM.findall(l))
    return {c: (S[True], S[False]) for c, S in cites.items()
            if S[True] and S[False] or (S[True] != S[False] and (S[True] or S[False]))}


def rule_cross(en, zh):
    out = {}
    for c in set(CITE.findall(en)) | set(CITE.findall(zh)):
        a = {n for l in en.splitlines() if re.search(re.escape(c), l) for n in NUM.findall(l)}
        b = {n for l in zh.splitlines() if re.search(re.escape(c), l) for n in NUM.findall(l)}
        if a and b and a != b: out[c] = (a, b)
    return out


def rule_nonempty(en, zh):
    return {c: v for c, v in rule_current(en, zh).items() if v[0] and v[1]}


RULES = {"现规则": rule_current, "候选①跨文件": rule_cross, "候选②一侧空集不参与": rule_nonempty}

led = pathlib.Path("RETRACTIONS.md").read_text().splitlines()
ent, cur = [], 0
for l in led:
    m = re.match(r'## Entry (\d+)', l); cur = int(m.group(1)) if m else cur
    ent.append(cur)
PAT = re.compile(r'两份 README 已改|两版 README|页面已改|页面.{0,6}更正|已更正.{0,10}页面|正确写法')
CAND = sorted({ent[j] for j, l in enumerate(led) if PAT.search(l) and ent[j] > 0})
print(f"账本里记过页面更正的条目:{len(CAND)} 个")

rows = []
for N in CAND:
    log = sh("git", "log", "--format=%H", "-S", f"## Entry {N} ", "--", "RETRACTIONS.md").split()
    if not log: rows.append(dict(entry=N, status="找不到提交")); continue
    c = log[-1]
    par = sh("git", "rev-parse", f"{c}^").strip()
    if not par: rows.append(dict(entry=N, status="无父提交")); continue
    d = sh("git", "show", c, "--", "README.md")
    removed = {n for l in d.splitlines() if l.startswith("-") and not l.startswith("---") for n in NUM.findall(l)}
    added = {n for l in d.splitlines() if l.startswith("+") and not l.startswith("+++") for n in NUM.findall(l)}
    wrong = removed - added
    if not wrong: rows.append(dict(entry=N, status="非数字更正(结构上抓不到)")); continue
    en = sh("git", "show", f"{par}:README.md"); zh = sh("git", "show", f"{par}:README_zh.md")
    if not en: rows.append(dict(entry=N, status="父提交无页面")); continue
    r = dict(entry=N, status="可测", n_wrong=len(wrong))
    for name, fn in RULES.items():
        try: cells = fn(en, zh)
        except Exception: r[name] = "错误"; continue
        caught = any((a - b) & wrong or (b - a) & wrong for a, b in cells.values())
        r[name] = bool(caught)
    rows.append(r)
T = pd.DataFrame(rows)
print("\n=== G3:回测总体 ===")
print(T.status.value_counts().to_string())
M = T[T.status == "可测"]
print(f"\n**可测的缺陷版本:{len(M)} 个**")
if len(M):
    for name in RULES:
        if name in M.columns:
            print(f"  {name:18s} 抓到 **{int(M[name].sum())}/{len(M)}**")
    print("\n  逐条(条目 · 错的数字个数 · 三条规则):")
    for r in M.itertuples():
        print(f"    #{int(r.entry):4d} 错数 {int(r.n_wrong):2d} · 现 {getattr(r,'现规则','?')} · "
              f"候选① {getattr(r,'候选①跨文件','?')} · 候选② {getattr(r,'候选②一侧空集不参与','?')}")

G = Gate("回测:这条规则在历史上抓到过什么?")
cur_hit = int(M["现规则"].sum()) if len(M) and "现规则" in M else 0
c1 = int(M["候选①跨文件"].sum()) if len(M) and "候选①跨文件" in M else 0
c2 = int(M["候选②一侧空集不参与"].sum()) if len(M) and "候选②一侧空集不参与" in M else 0
any_hit = max(cur_hit, c1, c2)
pos_ok = G.positive_control("回测必须有分辨力:至少一条规则在历史上抓到过 >0 个",
                            planted=float(any_hit), floor=0.0, spread=0.4)
if pos_ok:
    miss1 = int(((M["现规则"]) & (~M["候选①跨文件"])).sum()) if len(M) else 0
    miss2 = int(((M["现规则"]) & (~M["候选②一侧空集不参与"])).sum()) if len(M) else 0
    verdict = (f"现 {cur_hit} · 候选① {c1}(漏现规则 {miss1})· 候选② {c2}(漏现规则 {miss2}) ⇒ "
               + ("**不换规则**" if (miss1 or miss2) or max(c1, c2) <= cur_hit
                  else "**允许换**"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = (f"**UNVERIFIED —— 回测没有分辨力**:三条规则在 {len(M)} 个可测缺陷版本上"
               f"**一个都没抓到**({cur_hit}/{c1}/{c2})⇒ **不据它换规则,也不据它保留规则**")
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(candidates=CAND, table=T.to_dict("records"), testable=len(M),
               hits=dict(现规则=cur_hit, 候选1=c1, 候选2=c2), verdict=verdict, unchallenged=True),
          open(OUT/"backtest.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'backtest.json'}")
