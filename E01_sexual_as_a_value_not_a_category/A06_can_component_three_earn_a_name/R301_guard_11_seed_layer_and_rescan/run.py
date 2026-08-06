import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A73 R301 -- 守卫 11 的种子层硬规则,并回头扫有多少条结论建在错的来源上

**类型:CLOSURE**。

`#255c`:`#233` 只量了**阈值重抽样**那一层,而 `#300` 实测**数据劈分种子**那一层
可以比它大 40%(±2.5 vs ±1.8),并且正是它让 `#254a` 得出了一个错的结论。

⚠ **我注册的做法是「降级为 WARN」,这里改成 `FAIL`** —— **严于注册,方向安全**:
一个未被测量的**主导**不确定性来源,等于**精度未知**,而不是精度尚可。

GATE            四端自检(走**独立 `Gate` 实例**,`#249c`),只把「四端是否全对」汇报进主门:
                ① 只给阈值重抽样、无种子展布 → **FAIL**(新)
                ② 两层都给 → 放行,并按**较大**的那层给区间
                ③ 来源是零臂 → 仍然 FAIL(旧规则不能被新规则吃掉)
                ④ 展布为 0 且重抽 ≥10 → 仍然「钉死」放行(旧规则不能被新规则吃掉)
                外加:普查有多少条 `count_needs_interval` 调用只给了阈值层。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

SELF=Gate('守卫 11 种子层 —— 四端自检(独立实例,#249c)')
e1=SELF.count_needs_interval('①只给阈值层',8,29,1.8,'threshold_resample_阈值重抽样',n_resamples=12)
e2=SELF.count_needs_interval('②两层都给',8,29,1.8,'threshold_resample_阈值重抽样',n_resamples=12,seed_spread=2.5)
e3=SELF.count_needs_interval('③来源是零臂(旧规则不能被吃掉)',8,29,1.8,'null_零臂')
e4=SELF.count_needs_interval('④钉死(旧规则不能被吃掉)',0,29,0.0,'seed_跨种子',n_resamples=12)
got=[e1,e2,e3,e4]; exp=[False,True,False,True]
allok=(got==exp)
print(SELF)
print(f"四端:期望 {exp} · 实测 {got} -> **{'全对' if allok else '不全对'}**")

runs=sorted(ROOT.glob('E01_*/*/*/run.py'))
only_thr=[]; both=[]
for f in runs:
    t=f.read_text()
    for m in re.finditer(r"count_needs_interval\((.{0,260}?)\)\n",t,re.S):
        s=m.group(1)
        if 'threshold_resample' in s:
            (both if 'seed_spread' in s else only_thr).append(f.parent.name[:52])
print(f"\n普查:{len(runs)} 个 run.py 里 `count_needs_interval` 的调用 —— "
      f"**只给阈值层 {len(only_thr)} 处**,两层都给 {len(both)} 处")
for x in sorted(set(only_thr)): print(f"  ⚠ {x}")
T=pd.DataFrame([dict(v_round=x,source='只给阈值层') for x in sorted(set(only_thr))])
check_columns(T,'R301'); T.to_csv(pathlib.Path(__file__).parent/'results'/'only_threshold.csv',index=False)

g=Gate('守卫 11 的种子层硬规则')
g.asserted('⚠ 类型:CLOSURE;且我把注册的 WARN 改成了 FAIL —— 严于注册,方向安全',
           True, '未被测量的主导不确定性来源 = 精度未知,不是精度尚可')
g.asserted('★ 四端全对(自检走独立 Gate 实例;新规则没有吃掉旧规则)',
           allok, f"期望 {exp} · 实测 {got}")
g.asserted('★ 普查给出计数与名单',
           True, f"只给阈值层 {len(set(only_thr))} 个轮次:{sorted(set(only_thr))}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
