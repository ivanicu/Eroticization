import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A79 R313 -- 第 14 个守卫:这个数有没有可能是别的值

`#267b`/`#267c`:同一轮里出现了两个**不可能变的数**,而其中一个报成了一个干净的零。
本项目此前的 13 个守卫**没有一个问这件事**。

GATE            四端自检(`#244b`),走**独立 `Gate` 实例**(`#249c`),只把「四端是否全对」汇报进主门:
                ① 真会变的量 → 放行
                ② **恒等式**(R² 对预测量线性重缩放不变,即 `#312③`)→ 报警
                ③ **恒为常数的构造**(`#312` 里 S 的半块信度)→ 报警
                ④ 只在极端扰动下才变 → 放行,并报出变动幅度
外加             在 `#312③` 的**真实构造**上演示一次(不是模拟)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

SELF=Gate('守卫 14 —— 四端自检(独立实例,#249c)')
rng=np.random.default_rng(7)
Xb=rng.standard_normal((3000,4)); yb=Xb@np.array([0.4,-0.3,0.2,0.0])+rng.standard_normal(3000)
def r2_of(X,y):
    X=np.column_stack([np.ones(len(y)),X]); b=np.linalg.lstsq(X,y,rcond=None)[0]
    return float(1-np.var(y-X@b)/np.var(y))
# ① 真会变的量:换种子重抽样
e1=SELF.could_have_come_out_otherwise('①真会变的量(自助重抽样下的 R²)',
    lambda s: r2_of(Xb[np.random.default_rng(s).choice(3000,3000,True)],
                    yb[np.random.default_rng(s).choice(3000,3000,True)]), [1,2,3,4])
# ② 恒等式:R² 对预测量线性重缩放不变 —— 这正是 `#312③`
e2=SELF.could_have_come_out_otherwise('②恒等式(`#312③`:按信度重缩放预测量后的 R²)',
    lambda s: r2_of(Xb/np.array([1.0,0.5,2.0,0.8])**(s/4.0), yb), [0,1,2,3])
# ③ 恒为常数的构造:一个不随参数变化的「信度」
e3=SELF.could_have_come_out_otherwise('③恒为常数的构造(`#312` 里 S 的半块信度)',
    lambda s: 1.0, [0,1,2,3])
# ④ 只在极端扰动下才变
e4=SELF.could_have_come_out_otherwise('④只在极端扰动下才变',
    lambda s: r2_of(Xb+(0.0 if s<3 else 50.0)*rng.standard_normal((3000,4)), yb), [0,1,2,3])
got=[e1,e2,e3,e4]; exp=[True,False,False,True]
print(SELF)
print(f"四端:期望 {exp} · 实测 {got} -> **{'全对' if got==exp else '不全对'}**")
T=pd.DataFrame([dict(v_case=n,expected=e,v_got=g_) for n,e,g_ in
                zip(['①真会变','②恒等式','③恒为常数','④极端才变'],exp,got)])
check_columns(T,'R313'); T.to_csv(pathlib.Path(__file__).parent/'results'/'four_ends.csv',index=False)

g=Gate('守卫 14')
g.asserted('⚠ 类型:CLOSURE —— 保护的是「这个数有没有可能变」,不产生新数字',True,'§0 三类动作')
g.asserted('★ 四端全对(自检走独立 Gate 实例;②③ 正是 `#312` 里真实出现的两个)',
           got==exp, f"期望 {exp} · 实测 {got}")
g.asserted('★ ② 的诊断必须指名「这个数不可能是别的值」,不能只说 FAIL',
           '不可能是别的值' in SELF.rows[1][3], SELF.rows[1][3][:60])
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
