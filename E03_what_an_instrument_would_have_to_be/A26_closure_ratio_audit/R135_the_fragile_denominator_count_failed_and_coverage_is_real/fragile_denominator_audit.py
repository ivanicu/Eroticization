"""E03·A26·R135 —— 脆分母清点失败了,而它失败的原因本身是可报的

**类型:CLOSURE(诚实标注)。不推进对象。**
**裁决:预注册的脆分母清点 UNVERIFIED(④ 正对照失败,不得报);
而它失败的原因 —— 大多数条目从没报过自己的零 —— 用一具校准过的仪器量出来了。**

**这一轮没有关于人的那一句话,我不编一个。** 它保护的是页面的可读性:
**读者无法判断页面上哪个数是脆的,因为大多数数后面没有零。**

## ④ 正对照失败,先说这个
`#692` 写死:「`#677` 的 0.0423 / 0.0313 = 1.35 必须被这套清点复现,否则清点本身错了」。
我的清点用 `min(effs, key=|v − 2q|)` 去**猜**哪个效应配哪个零 —— **那是猜配对,不是测量**;
它把 `#675` 配成 `+0.4732 / 1.0000 = 0.47`(两个数根本不属于同一个检验),
而 `#677` 的 1.35 **没被复现**。**⇒ 清点作废,脆分母名单不发布。**

## 而覆盖率这个数不依赖配对,所以它可以报 —— 但仪器先自检
`#600` 起共 **93** 条。两个正对照:`#677` 必须被「打乱…的零」命中 ✅;
`#691` 必须被「零 95% 分位」命中 ✅ ⇒ **模式不是瞎的。**

| 模式 | 命中 |
|---|---|
| 紧:`零 (的) 95%` | **13 / 93 = 14%** |
| 松:紧 ∪ 打乱…的零 ∪ 置换零 ∪ 自助区间 | **20 / 93 = 22%** |
| 仅出现「零」字(**上界,不是测量**) | 55 / 93 = 59% |

**⇒ 报区间不报点:`#600` 起只有 14–22% 的条目报了自己的零,
其余 78–86% 结构性无法判定分母脆不脆。**

## 仪器(硬规则②)
**本轮唯一的仪器是账本文本自身**(`RETRACTIONS.md` 的 93 条)。**没有第二具仪器** ——
本轮不触碰任何数据集(GSS / SCCS / NSFG / MFQ / RWAS 一个都没读),
所以「跨仪器」在这里**不适用**,而不是被跳过。

## IMPOSSIBLE(不写 planned)
**78–86% 的条目没有可回算的零** ⇒ **脆分母清点在本账本上不是「难」,是数据不存在**;
本轮不重新检验任何声明。`[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
s=open("RETRACTIONS.md").read(); ent=re.split(r'\n## Entry ',s)
PATS={"strict":r'零\s*(?:的)?\s*95%',"shuffle":r'打乱[^。\n]{0,20}零',"perm":r'置换零',
      "boot":r'自助[^。\n]{0,10}区间|95%\s*(?:自助)?区间'}
rows={}
for e in ent[1:]:
    n=int(re.match(r'(\d+)',e).group(1))
    if n>=600: rows[n]={k:bool(re.search(v,e)) for k,v in PATS.items()}
N=len(rows)
strict=sum(1 for d in rows.values() if d["strict"])
loose=sum(1 for d in rows.values() if any(d.values()))
pc = rows[677]["shuffle"] and rows[691]["strict"]
print(f"#600 起 {N} 条 · 紧 {strict} ({100*strict/N:.0f}%) · 松 {loose} ({100*loose/N:.0f}%)")
print(f"④ 仪器正对照:#677 命中「打乱…的零」{rows[677]['shuffle']} · #691 命中「零 95%」{rows[691]['strict']}")
G=Gate("脆分母清点")
p1=G.positive_control("清点必须复现 #677 的 0.0423/0.0313 = 1.35",planted=0.0,floor=0.5,spread=0.01)
p2=G.positive_control("覆盖率仪器自检:两个已知条目必须被各自的模式命中",
                      planted=1.0 if pc else 0.0,floor=0.5,spread=0.01)
print(f"\n**脆分母清点:UNVERIFIED —— ④ 正对照失败(猜配对不是测量),名单不发布**")
print(f"**覆盖率:14–22%(紧/松两读数都报,不挑一个)⇒ 78–86% 无法判定**")
print(G)
json.dump(dict(n_entries=N,strict=strict,loose=loose,
               strict_pct=round(100*strict/N,1),loose_pct=round(100*loose/N,1),
               instrument_ok=bool(pc),fragility_list="UNVERIFIED — positive control failed",
               unchallenged=True),open(OUT/"coverage.json","w"),indent=1,ensure_ascii=False)
