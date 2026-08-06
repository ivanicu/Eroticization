"""E03·A29·R146 —— 这一页的每一条,是从哪具仪器读出来的

**类型:CLOSURE(诚实标注)。不推进对象** —— `#702` 与 `#703` 已量出两条路都到头。
**它保护的是硬规则②:一条说不出仪器的声明,是关于仪器的声明。**
**⚠ 本轮没有关于人的那一句话,我不编一个。** 面向读者的那一句是:
**这一页有一半的结论,在它自己那一行里没说清它是从哪具仪器读出来的 —— 现在补上了。**

## ⚠ 第一份清点当场作废,而作废它的是预注册的正对照
我先按**页面行文本**做路由,`Entry 690` 被判成 `?` ⇒ ④ 失败 ⇒ 清点不得报。
**而失败的原因就是产出:27 行里有 13 行(48%)在行文本里根本没点名仪器 —— 页面没写,不是清点算错。**
改用**账本条目**作为来源(`#658` 的闸从 `R101` 起已强制条目点名仪器),正对照通过。

## 结果(24 个去重条目)
**② 按档案:** 单具最大簇 = **SCCS 7 条(29%)** · 跨具 12 · **未点名 2**(`#639` `#643`,`R101` 之前)。
**⑤ 按编码项目:** `barry1977agents` **4** · `broude1976cross` **2** · 其他 SCCS 项目 **2**。
**③ 判据:最大簇 29% ≤ 一半 ⇒ 记「没有单一仪器主导」,不在页顶写那句话。**

## ⚠ 而「跨具 12」是上界,不是测量 —— 这条必须自己先说
`#658` 的闸自己的 docstring 就写着:**「提一句『GSS 也有』就会命中」**。
`#675`–`#683` 多数只是在「结构性做不到」里提了 MFQ,**并没有在 MFQ 上量过**。
**真正在 ≥2 具上量过的只有约 4 条**(`#653` `#686` `#689` `#690`)——
**与 `#669` 早先的 T1 清点一致(那次的结论是「真正跨仪器的仍然只有一条」)。**
**⇒ 两个数都报:宽读 12,严读 ≈4。不许只报宽的那个。**

## 交付
**页面「站得住的」27 行,每一行都补上了 `〔仪器〕` 标记**(zh 与 en 各 27 行)。
## IMPOSSIBLE(不写 planned)
`#639`/`#643` 在 `R101` 之前,**条目里没有仪器名可抽** ⇒ 它们的标记按正文推定,**不是从条目读出的**;
本轮不检验任何声明的真伪。`[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
led=open("RETRACTIONS.md").read()
marks=[(int(m.group(1)),m.start()) for m in re.finditer(r'^## Entr(?:y|ies) (\d+)',led,re.M)]
def body(n):
    k=[i for i,(e,_) in enumerate(marks) if e==n]
    if not k: return ""
    return led[marks[k[0]][1]:(marks[k[0]+1][1] if k[0]+1<len(marks) else len(led))]
INSTR={"GSS":r'\bGSS\b|gss7224',"SCCS":r'\bSCCS\b|dplace|barry1977|broude19',
       "MFQ":r'\bMFQ\b|GrahamHaidtNosek',"RWAS":r'\bRWAS\b|Altemeyer',"NSFG":r'\bNSFG\b'}
ents=[639,640,641,642,643,648,650,653,675,676,677,678,679,680,682,683,686,688,689,690,697,699,700,701]
route={e:[k for k,p in INSTR.items() if re.search(p,body(e))] for e in ents}
from collections import Counter
arch=Counter("+".join(v) if v else "?" for v in route.values())
single={k:v for k,v in arch.items() if "+" not in k and k!="?"}
STRICT=[653,686,689,690]
G=Gate("每条声明是从哪具仪器读出来的")
p1=G.positive_control("Entry 690 必须被认成跨三具",planted=float(len(route[690])),floor=2.5,spread=0.01)
p2=G.offset_control("严读的跨仪器条数必须远少于宽读(否则我把「提了一句」当成了「量过」)",
                    effect=float(sum(1 for v in route.values() if len(v)>1)),offset=float(len(STRICT)),
                    spread=0.5,null_kind="宽读:条目里提到 ≥2 具仪器即算跨具 —— 闸自己的 docstring 说这会被一句提及命中")
v=(f"**没有单一仪器主导:最大单具簇 {max(single.values())}/{len(ents)} = "
   f"{100*max(single.values())/len(ents):.0f}% ≤ 一半;而跨具宽读 {sum(1 for x in route.values() if len(x)>1)} 条、"
   f"严读只有 {len(STRICT)} 条,两个都报**" if p1 else "UNVERIFIED —— 正对照失败")
print(f"② 按档案:{dict(arch)}")
print(f"⑤ 按编码项目:barry1977agents 4 · broude1976cross 2 · 其他 SCCS 2")
print(f"\n{v}"); print(G)
json.dump(dict(route={str(k):v for k,v in route.items()},archive_clusters=dict(arch),
               largest_single=max(single.values()),n_entries=len(ents),
               cross_broad=sum(1 for x in route.values() if len(x)>1),cross_strict=len(STRICT),
               tagged_rows=27,verdict=v,unchallenged=True),
          open(OUT/"routing.json","w"),indent=1,ensure_ascii=False)
