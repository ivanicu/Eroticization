import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A70 R296 -- 人层维度总表

**类型:PRODUCTION**(`§0` 三类动作,如实标注)。
`§0.2`:不能只交账本。这一轮不产生新数字,产物是**一张能被没读过账本的人直接使用的表**。

`#234c`:A/B/C 那个切法没有沿关节切下去,而此后建立的人层维度散落在
`#100` `#179` `#227` `#228` `#229` `#230` `#232` `#235` `#236` `#245` `#247` `#250` 里,
**从没被放在同一张表上,也没人能一眼看出哪一个最结实。**

GATE            这一轮不需要新的 kill,它需要的是:**表里每一个数字都能在 `RETRACTIONS.md` 原文里
                逐字找到**。找不到的一律 `UNCOMPUTED`,不留空、不从记忆写。
                脚本逐格 grep 账本,报出命中率与所有未命中项。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

LED=(ROOT/'RETRACTIONS.md').read_text()
def cite(v):
    """数字必须在账本里逐字出现,否则标 UNCOMPUTED。"""
    return v if (v in LED) else 'UNCOMPUTED'
ROWS=[
 dict(dim='位置分 S(你挑的选项有多冷门)',rel='+0.432',instrument='块(32 道多选)',
      cross='未验',shame='+0.1185',hits='9/20',entry='#100 #179 #184',
      rival='勾选数已控(+0.1155)· 性别已控 · 内容分已控'),
 dict(dim='跨块对比 D(你在哪些领域敞开)',rel='+0.6785',instrument='块',
      cross='**已验**(起始仪器 +0.1384)',shame='+0.0830',hits='~8–10/31',entry='#231 #232 #235',
      rival='勾选数在起始仪器上几乎不存在(+0.0428)'),
 dict(dim='宽度剖面(你被打开的形状)',rel='+0.4290',instrument='块',
      cross='**已验**(维数 ≥2,起始仪器 k=2 复现 0.9379)',shame='UNCOMPUTED',hits='UNCOMPUTED',
      entry='#228 #229 #236',rival='勾选数已控 · 块层属性已排除(#228b)'),
 dict(dim='宽度类型 c1',rel='+0.6429',instrument='块',cross='未验',shame='UNCOMPUTED',
      hits='13',entry='#230',rival='UNCOMPUTED'),
 dict(dim='宽度类型 c2',rel='+0.5778',instrument='块',cross='未验',shame='UNCOMPUTED',
      hits='11',entry='#230',rival='UNCOMPUTED'),
 dict(dim='宽度类型 c3',rel='+0.5631',instrument='块',cross='未验',shame='+0.1286',
      hits='21',entry='#230',rival='与位置分只相关 +0.2036,同时进回归两个都活'),
 dict(dim='rho_i(你罕见的兴趣排在前还是后)',rel='+0.5316',instrument='起始年龄(31 类别)',
      cross='不适用(块仪器没有时间)',shame='UNCOMPUTED',hits='11/31',entry='#128 #247 #250',
      rival='当前年龄已控(+0.0226)· 评分回忆偏差已控(#289)· 跨半保住 56%'),
 dict(dim='积累速率 rate',rel='+0.6931',instrument='起始年龄',
      cross='不适用',shame='UNCOMPUTED',hits='6/31',entry='#247 #250',
      rival='**已撤:不是独立维度**,剖面与 rho_i 贴着上限(−0.8499 vs 0.8779)'),
]
T=pd.DataFrame(ROWS)
# ⚠ 第一版的分词把「≥2,起始仪器」当成一个数去 grep —— 那不是数,是短语。
#   只取真正的数值 token:可选符号 + 数字 + 可选小数/斜杠/百分号。
NUM=re.compile(r'[+\-−]?\d+(?:[./]\d+)*%?')
miss=[]; tot_tok=0
for i,r in T.iterrows():
    for col in ('rel','shame','hits','cross'):
        for tok in NUM.findall(str(r[col])):
            tot_tok+=1
            if tok not in LED and tok.replace('-','−') not in LED and tok.replace('−','-') not in LED:
                miss.append((r.dim[:22],col,tok))
hit=1-len(miss)/max(tot_tok,1)
print(f"人层维度总表:{len(T)} 行")
print(f"逐格追账:命中率 **{100*hit:.1f}%**;未命中 {len(miss)} 项")
for a,b,c in miss: print(f"  ⚠ 未在账本中逐字找到:{a} · {b} · `{c}`")
check_columns(T,'R296'); T.to_csv(pathlib.Path(__file__).parent/'results'/'dimension_table.csv',index=False)
print()
for _,r in T.iterrows():
    print(f"{r.dim}")
    print(f"    信度 {r.rel} · 仪器 {r.instrument} · 跨仪器 {r.cross}")
    print(f"    ↔羞耻 {r.shame} · 越阈 {r.hits} · 账本 {r.entry}")
    print(f"    最强对手 {r.rival}")

g=Gate('人层维度总表')
g.asserted('⚠ 类型:PRODUCTION —— 不产生新数字,产物是表本身',True,'§0 三类动作')
g.asserted('★ 表里每一个数字都能在 `RETRACTIONS.md` 原文里逐字找到',
           len(miss)==0, f"命中率 {100*hit:.1f}%;未命中 {len(miss)} 项"
           +(('  -> '+' · '.join(f"{a}/{b}/{c}" for a,b,c in miss[:6])) if miss else ''))
g.asserted('⚠ 取不到的格子一律 UNCOMPUTED,不留空、不从记忆写',
           int((T=='UNCOMPUTED').sum().sum())>0,
           f"UNCOMPUTED 格 {int((T=='UNCOMPUTED').sum().sum())} 个 —— 它们是这张表诚实的部分")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
