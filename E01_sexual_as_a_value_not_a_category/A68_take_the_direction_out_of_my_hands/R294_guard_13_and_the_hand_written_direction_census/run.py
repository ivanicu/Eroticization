import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A68 R294 -- 第 13 个守卫:把方向从我手里拿走;并普查还有多少条方向是我手写的

**类型:CLOSURE**(如实标注)。

`#248c`:预注册一个**方向**,和预注册一个**阈值**,不是同一件事。
阈值我算得出来;**方向我算不出来,我是在猜** —— 而我已经猜错四次
(`#132b` 审查 · `#134f` 剂量 · `#146e` 罕见度 · `#248c` 跨半种入),**四次都是扫描纠正了我。**

ESTIMAND        ① 落 `plant_direction_from_sweep`:**不接受调用者给的方向**,
                   只接受一条 `(g, 统计量)` 扫描,自己判单调性、`g=0` 落基线、并报灵敏度;
                ② 普查:现存 `run.py` 里有多少条正对照的方向是**我手写**的。
KILL            守卫必须在 **5 端**上全对(`#244b` 的新规矩:至少一端是「它可能误伤的正常情形」);
                普查必须给出一个**计数与名单**,不是一句「基本没问题」。
POSITIVE CTRL   五端:① 真的单调上升 → 放行;② 真的单调**下降** → 也放行(**方向无关**);
                ③ 非单调(自干扰的种入)→ 报警;④ 平(灵敏度未证明)→ 报警;
                ⑤ `g=0` 不在基线上 → 报警。
IMPOSSIBLE      普查是**词面**的(找「必须上升/加深/变大」这类词),
                所以它给的是**候选**;一条写着「必须上升」而实际做了扫描的会被误算进来。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

g=Gate('守卫 13 + 手写方向普查')
S=[('①真的单调上升',[(0,-0.26),(0.25,-0.07),(0.5,0.07),(1.0,0.29)],-0.26,0.13,True),
   ('②真的单调下降(方向无关)',[(0,0.09),(0.25,-0.02),(0.5,-0.18),(1.0,-0.36)],0.09,0.05,True),
   ('③非单调(自干扰的种入)',[(0,-0.26),(0.25,-0.30),(0.5,-0.10),(1.0,-0.28)],-0.26,0.13,False),
   ('④平(灵敏度未证明)',[(0,0.43),(0.02,0.44),(0.05,0.45),(0.10,0.48)],0.43,0.20,False),
   ('⑤g=0 不在基线上',[(0,0.10),(0.5,0.30),(1.0,0.50)],-0.26,0.13,False)]
got=[g.plant_direction_from_sweep(n,s,baseline=b,half_of=h) for n,s,b,h,_ in S]
allok=all(x==e for x,(_,_,_,_,e) in zip(got,S))

runs=sorted(ROOT.glob('E01_*/*/*/run.py'))
DIR=('上升','加深','变大','更大','更强','单调','超过')
tot=pos=hand=swept=0; rows=[]
for f in runs:
    t=f.read_text()
    for m in re.finditer(r"g\.asserted\(\s*'([^']{0,140})'",t):
        s=m.group(1); tot+=1
        if '正对照' in s:
            pos+=1
            hd=any(d in s for d in DIR); sw=('扫描' in s or 'sweep' in s.lower())
            two=('两端' in s or '两边' in s)
            if hd and not sw:
                hand+=1; rows.append(dict(v_round=f.parent.name[:46],two_ended=two,text=s[:90]))
            elif sw: swept+=1
T=pd.DataFrame(rows); check_columns(T,'R294')
T.to_csv(pathlib.Path(__file__).parent/'results'/'hand_written_directions.csv',index=False)
single=int((~T.two_ended).sum()) if len(T) else 0
print(f"普查:{len(runs)} 个 run.py · `g.asserted` {tot} 条 · 其中正对照 **{pos}** 条")
print(f"  方向由我手写且未提扫描:**{hand}** 条({100*hand/max(pos,1):.0f}%);"
      f"其中**单方向**(非两端式)**{single}** 条 <- 真正的待办")
print(f"  明确提到扫描的:{swept} 条")
for _,r in T.iterrows():
    print(f"    {'两端' if r.two_ended else '单向'}  {r['v_round']:<48} {r.text[:62]}")

g.asserted('★ 守卫 13 必须在 5 端上全对(其中 ② 是「它可能误伤的正常情形」:方向相反的真扫描)',
           allok, ' · '.join(f"{n[:12]}{'✅' if x==e else '❌'}" for x,(n,_,_,_,e) in zip(got,S)))
g.asserted('★ 普查必须给出计数与名单,不是一句「基本没问题」',
           len(T)==hand and hand>0, f"{hand} 条(单方向 {single} 条),名单已写入 results/")
g.asserted('⚠ 普查是词面的:一条写着「必须上升」而实际做了扫描的会被误算进来',
           True, '所以它给的是候选,不是判决')
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
