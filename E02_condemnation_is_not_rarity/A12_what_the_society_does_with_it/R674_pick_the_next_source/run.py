"""E02·A12·R674 —— 用那张边界图,选出这个项目还没碰过、而又足够密的那一篇

`#637` 的 NEXT。**行动类型:PRODUCTION + 人工判定**,如实标注。

⚠ **「本项目从未用过」不能凭记忆** —— 机械确定:
   把仓库里(账本 + 全部轮次脚本 + 两版页面)出现过的每一个 `SCCS\\d+` 抓出来,映射回它的 `source`。
   ⚠ 这是一个**搜索**,而 `#618` 记过我的搜索会失败 ⇒ **正对照**:
   它必须找回 `broude1976cross` 与 `ross1983political`(我确知用过的两篇)。

G1 ESTIMAND(先于方法):候选 = 44 篇合格 source 里同时满足
  (a) 对内中位 n >= 60 · (b) `k >= 10` · (c) **从未被本项目引用过**。
  然后**逐篇人工判定**它有没有「谴责 / 规范 / 制裁」类变量 —— **判定是我做的,逐条列出。**
CONTROLS:
  正对照(对**用过**的检测):必须找回 `broude1976cross` 与 `ross1983political`。
  **g=0**:一个**根本不存在**的变量号(如 `SCCS99999`)不得把任何 source 标成「用过」。
KILL(条件式,预注册于 `#637`):
  if 正对照找回两篇 and g=0 不误标:
      >= 1 篇候选带有谴责/规范/制裁类变量 -> **它是 E02 社会层的下一具仪器**
      0 篇 -> **如实记「社会层的谴责问题只能靠那两篇密度不足的性编码」**,写进页面
  else: UNVERIFIED
G3:候选全表 + 每篇的变量样例(供人工判定,让下一个人能推翻我)。
IMPOSSIBLE(不写 planned):「有没有谴责类变量」是**我读标题判的**,不是数据判的 ·
  密度够 **不等于** 能回答(`#635`:n=17 也可能因变量无方差而死)· `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import pandas as pd
from lib.gates import Gate

OUT = pathlib.Path(__file__).parent / "results"; OUT.mkdir(exist_ok=True)
S = ROOT / "data/external/dplace/repo/datasets/SCCS"
V = pd.read_csv(S/"variables.csv")
MAP = json.load(open(ROOT/"E02_condemnation_is_not_rarity/A12_what_the_society_does_with_it/"
                     "R673_source_boundary_map/results/source_map.json"))
T = pd.DataFrame(MAP["table"])

# ── 机械确定「用过」───────────────────────────────────────────
pat = re.compile(r'SCCS\d+(?:\.\d+)?')
seen = set()
for p in [ROOT/"RETRACTIONS.md", ROOT/"README.md", ROOT/"README_zh.md"]:
    seen |= set(pat.findall(p.read_text()))
for p in ROOT.glob("E0*/A*/R*/run.py"):
    seen |= set(pat.findall(p.read_text()))
v2s = dict(zip(V.id, V.source))
used_src = {str(v2s[x]) for x in seen if x in v2s}
print(f"仓库里出现过的 SCCS 变量号 {len(seen)} 个 · 映射到 **{len(used_src)}** 个 source")
pos = {"broude1976cross", "ross1983political"} <= used_src
g0 = "SCCS99999" in v2s
print(f"  正对照:找回 broude1976cross 与 ross1983political?**{pos}**")
print(f"  g=0:不存在的变量号 `SCCS99999` 会被映射吗?**{g0}**(须 False)")

# ── 候选 ────────────────────────────────────────────────────
T["used"] = T.source.apply(lambda s: any(str(s).startswith(u[:20]) for u in used_src))
C = T[(T.pair_med >= 60) & (T.k >= 10) & (~T.used)].sort_values("pair_med", ascending=False)
print(f"\n=== G3:候选(对内中位 n>=60 · k>=10 · 从未用过)= **{len(C)}** 篇 ===")
KEY = re.compile(r'punish|sanction|penalt|taboo|norm|approv|disapprov|condemn|shame|'
                 r'moral|deviance|crime|law|justice|restrict|prohibit|attitude', re.I)
rows = []
for r in C.itertuples():
    g = V[V.source.astype(str).str.startswith(r.source[:24])]
    hits = g[g.title.str.contains(KEY, na=False)]
    ex = list(g.title.head(4))
    rows.append(dict(source=r.source, k=int(r.k), pair_med=int(r.pair_med),
                     n_norm_vars=int(len(hits)), examples=[str(x)[:52] for x in ex],
                     norm_examples=[str(x)[:52] for x in hits.title.head(3)]))
    print(f"\n  {r.source:32s} k={int(r.k):3d} 对内中位 n={int(r.pair_med):3d} · "
          f"标题含规范/制裁词的变量 **{len(hits)}**")
    print(f"     变量样例: {ex[0][:60] if ex else '—'}")
    if len(hits): print(f"     规范类样例: " + " | ".join(str(x)[:44] for x in hits.title.head(3)))
R = pd.DataFrame(rows)
withnorm = R[R.n_norm_vars > 0] if len(R) else R
print(f"\n**候选里带规范/制裁类变量的:{len(withnorm)} / {len(R)} 篇**")

G = Gate("用边界图选出下一具社会层仪器")
pos_ok = G.positive_control("正对照:必须找回两篇确知用过的 source", planted=float(pos), floor=0.0, spread=0.4)
pla_ok = G.negative_control("g=0:不存在的变量号不得被映射", null=float(g0), effect=1.0,
                            null_spread=0.4, null_kind="一个根本不存在的变量号")
if pos_ok and pla_ok:
    verdict = (f"**{len(withnorm)} 篇候选带规范/制裁类变量** ⇒ "
               + ("**它们是 E02 社会层的下一具仪器候选**" if len(withnorm)
                  else "**0 篇 ⇒ 社会层的谴责问题只能靠那两篇密度不足的性编码**"))
    print(f"\n控制齐备 ⇒ 评判。**{verdict}**")
else:
    verdict = f"UNVERIFIED —— 控制未齐(正对照 {pos_ok} · g=0 {pla_ok})"
    print(f"\n⚠ {verdict}")
print(G)
json.dump(dict(used_sources=sorted(used_src), n_seen_vars=len(seen),
               candidates=R.to_dict("records"), n_candidates=int(len(R)),
               n_with_norm=int(len(withnorm)), verdict=verdict, unchallenged=True),
          open(OUT/"next_source.json","w"), indent=1, ensure_ascii=False)
print(f"\nwrote {OUT/'next_source.json'}")
