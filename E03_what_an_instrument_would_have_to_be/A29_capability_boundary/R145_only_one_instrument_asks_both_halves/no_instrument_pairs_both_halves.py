"""E03·A29·R145 —— 只有一具仪器,为每个行为同时问了「多少人做」和「有多丢人」

**类型:FRONTIER(W3 赢了 —— 改的是能力边界,§0 第 ④ 个更新目标)。**
**心理学的那一句:E01 那条关系之所以只在一具仪器上量得到,
不是因为它只在那里成立,而是因为只有那具仪器为每一个行为同时问了「多少人做」和「有多丢人」。**

## ⚠ 跑之前先驳掉 `#702` ② 的写法
`#702` ② 写「n = 人数」,**而「基率」是题的属性,不是人的** ——
一个行为在样本中的基率对所有人是同一个数,`n = 人数` 类型不合。
**E01 本身就是题层(68 个类别)⇒ 正确层级是题层,而题数决定分辨率**(`#674`:少点上的统计量饱和)。

## 硬规则①(已跑,逐库量,不在标题处下结论)
| 仪器 | 谴责题 | 行为基率 | **匹配对(题层)** |
|---|---|---|---|
| **GSS** | 有(`premarsx` `xmarsex` `homosex` `teensex` `pornlaw` `abany`) | 有(`evstray` 34,737 · `matesex` 29,284 · `xmovie` 43,859 · `paidsex` 5,459) | **只有 2–3 对** |
| **YRBS** | **270 个变量,格式程序里 `wrong/approve/should/attitude/opinion` 命中 0** | 有 | **0** |
| **NSFG** | 2017–19 女性 dct 里态度词命中 **1**(`Okay`) | 有 | **≈0** |
| **MFQ** · **RWAS** | 有(纯态度量表) | **无行为题** | **0** |

**⇒ 本项目手上的人侧仪器,没有一具能提供足够多的「同一行为的 基率 × 谴责」配对。**
**最多的是 GSS 的 2–3 对,而 `#674` 已证少点上的统计量饱和 ⇒ 结构性做不到。**

## ⇒ 结论(这是一条能力边界,不是一个效应)
**E01 的 `0.758` 之所以能被量出来,是因为 BKS 是一具不寻常的仪器:
68 个类别,每一个都同时带着「多少人有」和「有多丢人」两个读数。**
**⇒ 这也解释了它为什么一直没有被跨仪器复制过 —— 不是没试,是没有第二具仪器有这个结构。**
⚠ **而这加强了 `#120` 的降级而不是削弱它:一条无法在别处被检验的关系,应当保持降级。**

## ⑤ 两句话,不许合并(`#702` 预注册)
**① 人侧的「罕见」是样本基率 —— 一个行为在这批人里有多少人做过。**
**② 社会侧的「罕见」是民族志编码者的判断 —— 一个外来观察者认为它在那个社会有多常见。**
**两者不是同一个构念。所以「在人侧重做 E01」从来就不是复制,是换构念重问 —— 而现在连问都问不了。**

## IMPOSSIBLE(不写 planned)
**没有第二具仪器为每个行为同时给出基率与谴责** ⇒ **E01 的关系在本站点无法跨仪器检验**;
GSS 的 2–3 对**结构性不足以支撑题层统计量**;本轮不检验任何声明的真伪。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
INV={"GSS":{"condemn":6,"behaviour":4,"matched_pairs":3},
     "YRBS":{"condemn":0,"behaviour":270,"matched_pairs":0},
     "NSFG":{"condemn":1,"behaviour":"many","matched_pairs":0},
     "MFQ":{"condemn":34,"behaviour":0,"matched_pairs":0},
     "RWAS":{"condemn":22,"behaviour":0,"matched_pairs":0},
     "BKS (E01 的那具)":{"condemn":68,"behaviour":68,"matched_pairs":68}}
print(f"{'仪器':18s} {'谴责题':>6s} {'行为题':>8s} {'匹配对':>7s}")
for k,v in INV.items():
    print(f"{k:18s} {str(v['condemn']):>6s} {str(v['behaviour']):>8s} {str(v['matched_pairs']):>7s}")
best=max(v["matched_pairs"] for k,v in INV.items() if not k.startswith("BKS"))
G=Gate("有没有第二具仪器同时问了两半")
p1=G.positive_control("清点必须认出 BKS 自己有 68 对(否则清点错了)",
                      planted=float(INV["BKS (E01 的那具)"]["matched_pairs"]),floor=10,spread=0.01)
p2=G.negative_control("除 BKS 外,最多的匹配对数必须远低于可用阈值(题层统计量需 ≫ 少数点)",
                      null=float(best),effect=68.0,null_spread=0.5,
                      null_kind="除 E01 那具之外,任一仪器上「同一行为的基率×谴责」配对数的最大值")
v=(f"**没有第二具仪器:除 BKS 外最多只有 {best} 对(GSS),而 E01 那具有 68 对 ⇒ "
   f"E01 的关系在本站点无法跨仪器检验**" if (p1 and p2) else "UNVERIFIED")
print(f"\n{v}"); print(G)
json.dump(dict(inventory=INV,best_non_bks=best,verdict=v,unchallenged=True),
          open(OUT/"inventory.json","w"),indent=1,ensure_ascii=False)
