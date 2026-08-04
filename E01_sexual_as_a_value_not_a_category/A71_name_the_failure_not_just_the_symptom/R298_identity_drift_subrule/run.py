import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A71 R298 -- 守卫 13 的身份漂移子规则:说出失败的原因,而不只是症状

**类型:CLOSURE**。

`#252c`:守卫 13 面对「被测对象在扫描中途被替换」时,只会说「非单调」——**方向对但不够具体**。
一个由**数据定义**的量(特征向量、聚类中心、主成分),种得够强时**种入本身会重新定义它**。

GATE            四端自检(`#244b` 的规矩),**走独立的 `Gate` 实例**(`#249c` 定下的做法),
                只把「四端是否全对」这一条汇报进主门:
                ① 身份稳定 + 单调 → 放行
                ② 身份稳定 + 非单调 → 报「非单调」
                ③ **身份漂移 + 非单调 → 报「对象被替换」**(新)
                ④ **身份漂移 + 单调 → 也报「对象被替换」** —— 那时的单调是假的
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

SELF=Gate('守卫 13 身份子规则 —— 四端自检(独立实例,#249c)')
# ⚠ 第一版把「计数尺度」的灵敏度门槛(2.0)用在了**相关尺度**的扫描上(跨度只有 0.55),
#   于是 ① 因「灵敏度未证明」而 FAIL —— 那是我的测试用例参数错,不是守卫错。
#   每条用例自带它自己尺度上的门槛。
CASES=[
 ('①身份稳定+单调 → 放行',[(0,-0.26),(0.25,-0.07),(0.5,0.07),(1.0,0.29)],[1.0,0.99,0.98,0.97],True,0.13),
 ('②身份稳定+非单调 → 报「非单调」',[(0,-0.26),(0.25,-0.30),(0.5,-0.10),(1.0,-0.28)],[1.0,0.99,0.98,0.97],False,0.13),
 ('③身份漂移+非单调 → 报「对象被替换」',[(0,10.0),(0.05,12.1),(0.15,9.4),(0.40,17.4)],[1.0,0.95,0.62,0.20],False,2.0),
 ('④身份漂移+单调 → 也报「对象被替换」(单调是假的)',[(0,10.0),(0.05,12.0),(0.15,14.0),(0.40,17.4)],[1.0,0.93,0.55,0.18],False,2.0),
]
got=[SELF.plant_direction_from_sweep(n,s,baseline=s[0][1],baseline_spread=abs(s[0][1])*0.1+0.02,
                                     half_of=h,identity=idv)
     for n,s,idv,_,h in CASES]
allok=all(x==e for x,(_,_,_,e,_h) in zip(got,CASES))
print(SELF)
# ③④ 的说明文字必须真的说出「对象被替换」,而不只是 FAIL
# ⚠ 第一版只看 rows[i][3](说明字段),而「对象被替换」写在 rows[i][1](明细字段)里。
msg34=[f"{r[1]} || {r[3]}" for r in SELF.rows[2:4]]
named=all('替换' in m for m in msg34)
print(f"\n③④ 的诊断文字是否指名「对象被替换」:**{named}**")
for m in msg34: print(f"   {m}")
T=pd.DataFrame([dict(v_case=n,expected=e,v_got=x,detail=SELF.rows[i][1][:80])
                for i,((n,_,_,e,_h),x) in enumerate(zip(CASES,got))])
check_columns(T,'R298'); T.to_csv(pathlib.Path(__file__).parent/'results'/'four_ends.csv',index=False)

g=Gate('守卫 13 的身份漂移子规则')
g.asserted('⚠ 类型:CLOSURE —— 保护的是「诊断的具体性」,不产生新数字',True,'§0 三类动作')
g.asserted('★ 四端全对(自检走独立 Gate 实例,#249c;只把这一条汇报进主门)',
           allok, ' · '.join(f"{n[:14]}{'✅' if x==e else '❌'}" for x,(n,_,_,e,_h) in zip(got,CASES)))
g.asserted('★ ③④ 的诊断必须指名「对象被替换」,不能只说 FAIL —— 这一轮的全部意义就在这里',
           named, ' | '.join(m[:52] for m in msg34))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
