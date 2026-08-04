import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A105 R357 -- guard 21:为「零」设计的判据

`#308c` 与 `#311d` 是同一条教训的两次出现,而两次我都是**事后用散文解释**为什么那个 FAIL 不算数。
**散文不是执行**(P9)。这一轮把它变成一个能自己开火的守卫。

⚠ **这是 Closure**(工具),不是 Frontier。

ESTIMAND        `null_claim_uses_null_criteria`:一轮声明结论是 `EFFECT` 还是 `NULL`;
                `NULL` 时强制**置换分位数 · MDE · 正对照证明的灵敏度**三样同时在场,
                并要求 **MDE < 一个有意义的效应量**,否则这个零没有内容。
POSITIVE CTRL   拿 `#311`(合格的零:分位数 58.3% · MDE 10–30% · 扫描灵敏度)-> 必须 PASS。
NEGATIVE CTRL   拿 `#310a`(不合格的零:MDE 55% 而 30% 的缓冲有意义)-> 必须 FAIL;
                以及三种缺件形式 -> 各自 FAIL。
⚠ 第三类        拿一个**有效应**的轮次(`#308a` 的加性)-> `EFFECT` 分支必须**不触发**。
IMPOSSIBLE      守卫检的是**证据是否到齐**,不检那些数**算得对不对** ——
                一个编造的 MDE 一样能通过。它挡的是遗漏,不是造假。
"""
import pandas as pd,hashlib
from lib.gates import Gate

CASES=[
 ('★ 正对照:`#311` 的零(合格)','NULL',0.583,0.20,'扫描 10–30% 抓到',0.30,True),
 ('★ 负对照:`#310a` 的零(MDE 55%,30% 缓冲有意义)','NULL',0.42,0.55,'植入 30% 读出 77%',0.30,False),
 ('★ 第三类:`#308a` 的加性(有效应)','EFFECT',None,None,None,None,True),
 ('负对照:缺置换分位数','NULL',None,0.20,'有',0.30,False),
 ('负对照:缺 MDE','NULL',0.58,None,'有',0.30,False),
 ('负对照:缺正对照灵敏度','NULL',0.58,0.20,None,0.30,False),
 ('负对照:没声明 claim_kind','(未声明)',0.58,0.20,'有',0.30,False),
]
g=Gate('guard 21 的验收')
res=[]
for nm,kind,q,mde,sens,mean_,want in CASES:
    got=g.null_claim_uses_null_criteria(nm,kind,perm_quantile=q,mde=mde,
                                        sensitivity_shown=sens,meaningful=mean_)
    res.append(dict(v_case=nm[:44],want=want,got=got,ok=(got==want)))
print(g)
T=pd.DataFrame(res)
print(f"\n七个用例:**{int(T.ok.sum())}/{len(T)}** 与预期一致")
for _,r in T[~T.ok].iterrows(): print(f"   ❌ {r.v_case}: 期望 {r.want} 得到 {r.got}")
T.to_csv(pathlib.Path(__file__).parent/'results'/'guard21.csv',index=False)
g2=Gate('guard 21 本身')
g2.asserted('★ 七个用例全部与预期一致(三缺一 · MDE 无内容 · EFFECT 不触发 · 未声明)',
            bool(T.ok.all()),f"{int(T.ok.sum())}/{len(T)}")
g2.asserted('⚠ 边界:守卫检的是证据是否到齐,不检那些数算得对不对',True,
            '一个编造的 MDE 一样能通过 —— 它挡的是遗漏,不是造假')
print(g2)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
