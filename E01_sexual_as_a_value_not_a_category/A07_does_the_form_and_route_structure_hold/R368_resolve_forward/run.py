import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A113 R368 -- 那份 29 条的清单是个过计数

`#322a` 产出了 29 条「只有点估计」的 UNVERIFIED。**按规则挑第一条 -> `#24`。**
**而 `#24` 就是 POWER–SUBSTANCE,`#321` 上一轮刚给了它区间 [+0.561, +0.659]。**

> **清单把「这一条自己有没有区间」当成了「这一条有没有被解决」。**
> **一个 UNVERIFIED 可以被后来的某一轮解决,而解决它的证据写在那一轮里,不在它自己里。**

ESTIMAND        修正的清单:一条 ② 若被**后续条目引用且那条带区间/MDE**,判为**已被前向解决**;
                报修正前后的计数,以及**被前向解决的比例**。
KILL            **若前向解决的比例很大 -> 那 29 条是过计数,真正的欠账要小得多;
                若很小 -> 清单基本准确,补测工作量就是它。**
POSITIVE CTRL   **`#24` 必须被判为「已被 `#321` 前向解决」**(已知答案)。
NEGATIVE CTRL   把引用关系随机重排 -> 前向解决的比例必须掉回随机水平。
⚠ 边界         「被引用且那条有区间」是**代理**:那条的区间未必是给**这一条**的。
                P6 安全侧:**判为已解决 = 待人读;判为未解决 = 确定要补。**
"""
import re,pandas as pd,numpy as np,hashlib
from lib.gates import Gate, check_columns

LED=pathlib.Path('RETRACTIONS.md').read_text()
ENT={};cur=None
for l in LED.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); ENT[cur]=[]
    elif cur is not None: ENT[cur].append(l)
ENT={k:'\n'.join(v) for k,v in ENT.items()}
def body(t):
    for mark in ('**NEXT(','**NEXT (','**NEXT'):
        i=t.find(mark)
        if i>0: return t[:i]
    return t
BODY={k:body(v) for k,v in ENT.items()}
UNV=re.compile(r'UNVERIFIED|不可读|读不了|读不出来|分不开',re.I)
IVAL=re.compile(r'\bMDE\b|95% CI|CI \[|区间|分位数|percentile|自助|bootstrap|置换零.{0,20}±|± *0\.\d',re.I)
FAILPC=re.compile(r'正对照.{0,24}(FAIL|没过|失败|不合格)|positive control.{0,20}fail|仪器(坏|失败|不合格)|测不动',re.I)
lst=[k for k,t in BODY.items() if UNV.search(t) and not FAILPC.search(t) and not IVAL.search(t)]
print(f"② 只有点估计:**{len(lst)}** 条")
FWD={}
for k in lst:
    solvers=[j for j,t in BODY.items() if j>k and re.search(rf'`#{k}[a-z]?`',t) and IVAL.search(t)]
    FWD[k]=solvers
solved=[k for k,v in FWD.items() if v]
print(f"★ 其中**被后续带区间的条目引用**(判为已被前向解决)的:**{len(solved)}** 条"
      f"({100*len(solved)/max(len(lst),1):.0f}%)")
print(f"   真正的欠账:**{len(lst)-len(solved)}** 条 -> {sorted(set(lst)-set(solved))[:14]}")
print(f"\n正对照:`#24` 的前向解决者 = {FWD.get(24,'(不在清单里)')}")
rg=np.random.default_rng(77); keys=sorted(BODY)
NG=[]
for _ in range(30):
    perm=dict(zip(keys,rg.permutation(keys)))
    c=0
    for k in lst:
        if any(perm[j]>k and IVAL.search(BODY[j]) and re.search(rf'`#{k}[a-z]?`',BODY[j]) for j in keys): c+=1
    NG.append(c/max(len(lst),1))
print(f"负对照(随机重排引用方向 30 次):前向解决比例 **{np.mean(NG):.2f} ± {np.std(NG):.2f}** "
      f"vs 观测 **{len(solved)/max(len(lst),1):.2f}**")
T=pd.DataFrame([dict(v_entry=k,n_solvers=len(v),solvers=' '.join(f"#{x}" for x in v[:4])) for k,v in FWD.items()])
check_columns(T,'R368'); T.to_csv(pathlib.Path(__file__).parent/'results'/'forward.csv',index=False)
gg=Gate('那份清单是不是过计数')
gg.asserted('★ 正对照:`#24` 必须被判为已被 `#321` 前向解决',
            24 in FWD and any(x==321 for x in FWD[24]),
            f"#24 的前向解决者 {FWD.get(24)}")
gg.asserted('★ 负对照:随机重排引用方向后比例掉回随机水平',
            (len(solved)/max(len(lst),1))-np.mean(NG)>2*np.std(NG),
            f"观测 {len(solved)/max(len(lst),1):.2f} vs 重排 {np.mean(NG):.2f} ± {np.std(NG):.2f}")
gg.asserted('★ 注册的 kill:前向解决的比例',len(solved)/max(len(lst),1)>0.3,
            f"**{len(solved)}/{len(lst)}** = {100*len(solved)/max(len(lst),1):.0f}% —— "
            f"{'那 29 条是过计数' if len(solved)/max(len(lst),1)>0.3 else '清单基本准确'}")
gg.asserted('⚠ 安全侧(P6):判为已解决 = 待人读;判为未解决 = 确定要补',True,
            '「被引用且那条有区间」是代理 —— 那条的区间未必是给这一条的')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
