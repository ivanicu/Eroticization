"""E03·A29·R147 —— 跨仪器让一条结论更可信,但没有让它更精确

**类型:CLOSURE(诚实标注)。裁决:`#704` 预注册的那个比较 UNVERIFIED —— 检查不合用。**
**⚠ 本轮没有关于人的那一句话,我不编一个。** 面向读者的那一句是:
**跨仪器复制让一条结论更可信,但没有让它更精确 —— 而这一页此前把这两件事当成了一件。**

## ⚠ 跑之前抓到两处,第二处让预注册的比较作废
**① 跨仪器组只有 n = 2,不是 `#704` 写的 4** —— `#653` 与 `#686` 的三具行**没有可配套的(效应,零)**
(`#695` 的纪律:只认脚本自己写进同一份输出的),能配套的只有 `#690` 与 `#689`。
**② G1 失败,而写坏它的是我自己在 `#704` 里的估计量。**
我取的 `#690` 的 `r_tot / q_tot` **是 RWAS 那一具上自己的相关与自己的零**,不是一个跨仪器统计量;
**「效应 ÷ 自身零」量的是「这一具仪器上这次测量有多精确」,而不是「跨仪器让结论多硬」。**
**跨仪器的硬度根本不是一个比值,它是一个类别事实:这条在不在 ≥2 具上成立。**
**⇒ 预注册的比较不测它自己声称要测的东西 ⇒ 作废,不跑。**

## ⑤ 停止检查(`#704` 预注册)照样跑完并登记
随机分组零的可能取值档数 = **21**(而总组合数也是 21,C(7,2))⇒ **档数 ≥20 通过,
但最细的经验 p 只能到 1/21 = 0.048** —— **这个设计的 p 天花板就在 0.05 边上,记下来。**

## 而那张表本身可读(全部来自各轮自己的 JSON,效应与零同源)
| 条目 | 效应 | 零 95% | 比值 | 仪器 |
|---|---|---|---|---|
| `#686` 人内压缩 | 0.1680 | 0.0091 | **18.46** | GSS(原始) |
| `#690` 三具之一 | 0.1069 | 0.0183 | **5.84** | RWAS(复制) |
| `#678` 结构不动 | 0.2189 | 0.0720 | 3.04 | GSS |
| `#699` 施加轴 | 0.3243 | 0.1231 | 2.63 | SCCS |
| `#698` 跨手段 | 0.3241 | 0.1684 | 1.93 | SCCS |
| `#689` 二具之一 | 0.0542 | 0.0328 | **1.66** | MFQ(复制) |
| `#700` 轴×谴责 | 0.3631 | 0.2381 | 1.52 | SCCS |

**⇒ 这一页最精确的一次测量是单仪器的(GSS,18.46);
而第二、第三具仪器上的复制是 1.66 与 5.84 —— 复制出来的数总是更松,因为那些样本更小、零更宽。**
**⇒ 「更可信」与「更精确」是两件事,而这一页此前把它们混过。**

## IMPOSSIBLE(不写 planned)
**跨仪器硬度不是一个可比的连续量** ⇒ 本站点无法把它与精确度放在同一把尺子上;
可配套(效应,零)的条目只有 7 条(`#695`:24 条里 7 条可配套);
本轮不检验任何声明的真伪。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, itertools
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
T={"#686 人内压缩(GSS 原始)":(0.1680,0.0091,False),
   "#690 三具之一(RWAS 复制)":(0.1069,0.0183,True),
   "#678 结构不动(GSS)":(0.2189,0.0720,False),
   "#699 施加轴(SCCS)":(0.3243,0.1231,False),
   "#698 跨手段(SCCS)":(0.3241,0.1684,False),
   "#689 二具之一(MFQ 复制)":(0.0542,0.0328,True),
   "#700 轴×谴责(SCCS)":(0.3631,0.2381,False)}
for k,(e,n,c) in T.items(): print(f"  {k:26s} 效应 {e:.4f} · 零 {n:.4f} · 比值 **{e/n:5.2f}** · {'复制' if c else '原始'}")
allr=[e/n for e,n,_ in T.values()]; k=sum(1 for *_,c in T.values() if c)
combs=list(itertools.combinations(range(len(allr)),k))
diffs=sorted({float(np.median([allr[i] for i in c])-np.median([allr[i] for i in range(len(allr)) if i not in c])) for c in combs})
print(f"\n⑤ 随机分组零:可能取值 **{len(diffs)}** 档 · 总组合 {len(combs)} ⇒ 最细经验 p = **1/{len(combs)} = {1/len(combs):.3f}**")
rep=[e/n for e,n,c in T.values() if c]; ori=[e/n for e,n,c in T.values() if not c]
print(f"   复制组比值 {['%.2f'%v for v in rep]} · 原始组 {['%.2f'%v for v in ori]}")
G=Gate("跨仪器让结论更硬吗")
p1=G.positive_control("#690 的复制必须落在复制组里",planted=1.0 if T["#690 三具之一(RWAS 复制)"][2] else 0.0,floor=0.5,spread=0.01)
v=("**UNVERIFIED —— `#704` 预注册的估计量不测它声称要测的东西:"
   "「效应÷自身零」量的是这一具上的精确度,而跨仪器硬度是一个类别事实,不是比值。比较作废,不跑。**"
   "\n**而可报的是描述:最精确的一次测量是单仪器的 18.46,两次复制只有 1.66 与 5.84 ⇒ "
   "更可信 ≠ 更精确。**")
print(f"\n{v}"); print(G)
json.dump(dict(table={k:{"effect":e,"null":n,"ratio":e/n,"replication":c} for k,(e,n,c) in T.items()},
               null_levels=len(diffs),n_combinations=len(combs),p_floor=1/len(combs),
               verdict="UNVERIFIED — 预注册估计量不合用(G1)",unchallenged=True),
          open(OUT/"credible_vs_precise.json","w"),indent=1,ensure_ascii=False)
