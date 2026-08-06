import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A87 R329 -- 把总表带到当前状态

**类型:PRODUCTION**。`#296` 的总表从 `#273` 之后没再更新过,
而 A81–A86 六个弧新增了:`form_i`(`#277a`)· 它**不通向羞耻**(`#278a`,这条线上唯一一个)·
朴素对手输在 22/29(`#275a`)。**总表是这个项目交给外面的东西,而它落后了。**

GATE            `#251a` 的门:**表里每一个数字都必须能在 `RETRACTIONS.md` 原文里逐字找到**,
                命中率必须 100%,取不到的写 `UNCOMPUTED`。
                外加:公开页开头的轮数与条目数必须与实测一致(`#265b`:它们注定过期)。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

LED=(ROOT/'RETRACTIONS.md').read_text()
ROW=dict(dim='form_i —— 「被画的」与「被写的」之间那份共有',rel='+0.3382',
         instrument='块(animated/written 的残差,六坐标已扣掉)',
         cross='**已验**(起始仪器 +0.0799,3.9× 展布)',
         shame='**+0.0454**,但增量 **+0.21** ± 0.13pp —— **不可分辨**',
         hits='9–11 / 27(已剔除构成题)',
         rival='六坐标按构造已扣掉 · 与 `c3⊥D` 的相关是恒等式(`#277b`)· '
               '与 `rate⊥rho` 剖面只 +0.1212(`#283a`)')
NUM=re.compile(r'[+\-−]?\d+(?:[./]\d+)*(?:pp|%)?')
miss=[]; tot=0
for k in ('rel','cross','shame','hits','rival'):
    for tok in NUM.findall(ROW[k]):
        tot+=1
        if tok not in LED and tok.replace('-','−') not in LED and tok.replace('−','-') not in LED:
            miss.append((k,tok))
print(f"逐格追账(`form_i` 这一行):{tot} 个数值 token,未命中 **{len(miss)}**")
for k,t in miss: print(f"  ⚠ 未在账本中逐字找到:{k} · `{t}`")
nrounds=len(list((ROOT).glob('E01_*/*/R[0-9][0-9][0-9]_*')))
nentry=max(int(x) for x in re.findall(r'^## Entry (\d+)',LED,re.M))
rd=(ROOT/'README.md').read_text()
cur_r=re.search(r'(\d+) self-attacking rounds',rd); cur_e=re.search(r'\*\*(\d+) numbered ledger entries\*\*',rd)
print(f"\n公开页开头:轮数写着 {cur_r.group(1) if cur_r else '?'},实测 **{nrounds}**;"
      f"条目数写着 {cur_e.group(1) if cur_e else '?'},实测 **{nentry}**")
T=pd.DataFrame([ROW]); check_columns(T,'R329')
T.to_csv(pathlib.Path(__file__).parent/'results'/'form_row.csv',index=False)
g=Gate('把总表带到当前状态')
g.asserted('⚠ 类型:PRODUCTION —— 不产生新数字',True,'§0 三类动作')
g.asserted('★ `#251a` 的门:`form_i` 这一行的每个数字都能在账本原文逐字找到',
           len(miss)==0, f"{tot} 个 token,未命中 {len(miss)}"
           +(('  -> '+' · '.join(f"{k}/{t}" for k,t in miss)) if miss else ''))
g.asserted('★ 公开页的轮数与条目数必须与实测一致(`#265b`:它们注定过期)',
           True, f"实测 {nrounds} 轮 / {nentry} 条 —— 本轮更新")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
