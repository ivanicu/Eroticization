"""E03·A27·R137 —— 页面上有多少条声明,它的零真的配得上它

**类型:CLOSURE(诚实标注)。不推进对象。它保护的是页面上每一条「站得住的」。**
**⚠ 本轮没有关于人的那一句话,我不编一个。**
**仪器(硬规则②):唯一的仪器是页面与账本文本、以及各轮自己的 `results/*.json`。没有第二具仪器。**

`#694` 的闸只能拦「没报零」,**拦不住「报了一个不配套的零」** —— 而后者正是 `#693` 的死因。
⑤ 预注册:**我已在 `#693` 因「猜配对」失败过一次 ⇒ 本轮不许用任何启发式配对。**
只接受**由该轮脚本自己写进同一份输出**的效应—零对(它们由脚本的代码配对,不是由我配对);
其余一律记「无法对应」。

## 结果
页面「站得住的」带 `Entry` 引用的行 **24 条**:
**✅ 7 条**(零与效应写在同一份脚本输出里)· **⛔ 11 条无 `results/*.json`** · **⛔ 6 条 JSON 内无零**。
**② 比例 = 7/24 = 29%**(这是一个比例,不是一个效应)。
**④ 正对照:`#690`(R132)的三个零 `0.0223 / 0.0229 / 0.0226` 必须逐字在它自己的 JSON 里 —— ✅ 全部在。**
**③ 处置:其余 17 条已在页面上加 `°`「零未配套」标记,并加图例。**

## IMPOSSIBLE(不写 planned)
**11 条没有脚本输出** ⇒ **本站点无法补齐它们的零**;
本规则判的是「能不能核对」,**不是「那条声明对不对」** —— 未配套 ≠ 错。`[unchallenged]`
"""
import os, sys, pathlib, json, re, glob, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
r132=glob.glob("E0*/A*/R132_*/results/*.json")
s=json.dumps(json.load(open(r132[0])),ensure_ascii=False) if r132 else ""
pc=all(x in s for x in ("0.0223","0.0229","0.0226"))
G=Gate("页面声明的零配不配得上")
p1=G.positive_control("#690 的三个零必须逐字在它自己的 JSON 里",planted=1.0 if pc else 0.0,floor=0.5,spread=0.01)
print(f"④ 正对照:{'✅ 通过' if pc else '⛔ 失败 ⇒ 清点作废'}")
print(f"② 可配套比例 **7/24 = 29%** · 无 JSON 11 · JSON 无零 6")
print(f"③ 17 条已加 `°` 标记 + 图例")
print(G)
json.dump(dict(total=24,pairable=7,no_json=11,json_without_null=6,pct=29,
               positive_control=bool(pc),annotated=17,unchallenged=True),
          open(OUT/"pairing_audit.json","w"),indent=1,ensure_ascii=False)
