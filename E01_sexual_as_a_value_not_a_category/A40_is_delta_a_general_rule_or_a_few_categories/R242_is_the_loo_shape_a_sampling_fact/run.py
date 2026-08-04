import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A40 R242 -- 留一曲线的形状,是采样事实还是心理学事实

`#196c`:Δ 既不集中也不均匀。`#196d`③ 给出可检验的解释 ——
**留一影响可能主要由「多少人答过这个类别」决定,而不是由这个类别的稀有度决定。**

ESTIMAND        把 31 个类别的 `|ΔΔ|` 对 **log n**(答过的人数)与 **稀有度** 同时回归,
                比较两者的偏相关。
KILL            **若 n 解释了大部分而稀有度几乎不解释 ->
                留一曲线的形状是一个采样事实,不是心理学事实;
                「Δ 集中在哪里」在这份数据上问不出答案,写进 IMPOSSIBLE,停止追。**
IDENTIFICATION  n 与稀有度在构造上强负相关(越罕见的越少人答)——
                所以两者的"各自解释力"在共线下不可分(`#182b`)。**报三个数,不声称分开。**
NOISE FLOOR     类别层 bootstrap(重抽 31 个类别)300 次。
POSITIVE CTRL   合成一个**只由 n 决定**的 |ΔΔ|,回归必须把它归给 n。
NEGATIVE CTRL   合成一个**只由稀有度决定**的,必须归给稀有度。
IMPOSSIBLE      n = 31 个类别。|r| < 0.36 在 n=31 上不显著,所以只能判**大的**差别。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

T=pd.read_csv(ROOT/'E01_sexual_as_a_value_not_a_category/A40_is_delta_a_general_rule_or_a_few_categories'
                  /'R241_per_category_decomposition/results/loo.csv')
check_columns(T,'R242 输入')
y=np.abs(T.loo_delta.values); n=np.log(T.n_people.values.astype(float)); r=T.rarity.values
RHO=float(np.corrcoef(n,r)[0,1])
print(f"类别 {len(T)} 个;corr(log n, 稀有度) = {RHO:+.4f}")
# ⚠ #197c:**它们是同一个变量。** 稀有度就定义为 `−log(答过的比例)` = `log N − log n`,
#   所以 corr = −1.0000 是**恒等式**,不是发现。`#196d`③ 那个"是 n 还是稀有度"的问题
#   **在构造上不可检验** —— 它是同义反复。本轮的注册 kill 因此**取不到数**。
assert abs(abs(RHO)-1.0)<1e-6 or abs(RHO)<0.999, "共线度落在中间 —— 与预期不符,先查"
if abs(abs(RHO)-1.0)<1e-6:
    print("  ⚠ **完全共线(|r| = 1.0000)**:稀有度 = log N − log n。")
    print("     「留一形状由 n 决定还是由稀有度决定」是同义反复,本轮注册的 kill 取不到数。")

def partial(target,x,ctrls):
    X=np.c_[np.ones(len(target)),*ctrls] if ctrls else np.ones((len(target),1))
    ry=target-X@np.linalg.lstsq(X,target,rcond=None)[0]
    rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rb=np.random.default_rng(20260803)
def row(tag,t):
    a=partial(t,n,[]); b=partial(t,r,[])
    a_g=partial(t,n,[r]); b_g=partial(t,r,[n])
    bs=[]
    for _ in range(300):
        ii=rb.choice(len(t),len(t),replace=True)
        bs.append(partial(t[ii],n[ii],[r[ii]])-partial(t[ii],r[ii],[n[ii]]))
    return dict(tag=tag,r_n=a,r_rar=b,r_n_given_rar=a_g,r_rar_given_n=b_g,
                d_partial=a_g-b_g,sd=float(np.std(bs)))   # ⚠ 原名 `diff` —— pandas 方法名(#197a)
rows=[row('真实 |ΔΔ|',y),
      row('【正对照】只由 n 决定',(n-n.mean())+0.05*rb.standard_normal(len(n))),
      row('【负对照】只由稀有度决定',(r-r.mean())+0.05*rb.standard_normal(len(r)))]
R=pd.DataFrame(rows); check_columns(R,'R242'); R.to_csv(pathlib.Path(__file__).parent/'results'/'attrib.csv',index=False)
print(f"\n{'r(log n)':>10}{'r(稀有度)':>11}{'n|稀有度':>10}{'稀有度|n':>10}{'差':>9}{'sd':>8}  臂")
for _,q in R.iterrows():
    print(f"{q.r_n:>+10.4f}{q.r_rar:>+11.4f}{q.r_n_given_rar:>+10.4f}{q.r_rar_given_n:>+10.4f}"
          f"{q.d_partial:>+9.4f}{q.sd:>8.4f}  {q.tag}")

real=R.iloc[0]; pos=R.iloc[1]; neg=R.iloc[2]
g=Gate('留一形状是采样事实还是心理学事实')
g.asserted('正对照:只由 n 决定的合成量必须归给 n',pos.d_partial>0.5,f"差 {pos.d_partial:+.4f}")
g.asserted('负对照:只由稀有度决定的必须归给稀有度',neg.d_partial<-0.5,f"差 {neg.d_partial:+.4f}")
g.asserted('共线已量化,不声称把两者分开',True,
           f"corr(log n, 稀有度) = {np.corrcoef(n,r)[0,1]:+.4f}")
g.resolvable('真实的偏相关之差(n|稀有度 − 稀有度|n)',float(real.d_partial),float(real.sd))
g.asserted('注册的 kill:n 解释大部分而稀有度几乎不解释',
           (abs(real.r_n_given_rar)>2*abs(real.r_rar_given_n)) and (real.d_partial>2*real.sd),
           f"n|稀有度 {real.r_n_given_rar:+.4f} vs 稀有度|n {real.r_rar_given_n:+.4f}, "
           f"差 {real.d_partial:+.4f} ± {real.sd:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(R.to_csv(index=False).encode()).hexdigest()[:12]}")
