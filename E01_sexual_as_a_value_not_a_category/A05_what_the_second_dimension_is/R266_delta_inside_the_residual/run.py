import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A53 R266 -- Δ 在那个残差矩阵里,是多大的一个成分

`#220`:双约束零与 Δ 在定义上纠缠,**因为 Δ 是两侧去均值之后的量,而去均值移除的正是两侧边际**。
这指向一个更前置的问题:**Δ 这个估计量,到底还剩下什么不是边际 —— 它有多大。**

ESTIMAND        两侧去均值后的残差 `D`(只在观测格上)里:
                ① **Δ 方向**(固定 `x_j = 稀有度 − 均值`,人载荷自由)解释的平方和占比;
                ② **最佳秩一方向**(x 也自由,带缺失的幂迭代)解释的占比;
                ③ 随机 `x` 方向的占比(基线)。
KILL            **若 Δ 方向占的比例远小于最佳秩一方向(<1/5),那么这个项目一直在讨论的
                是残差里一个次要方向 —— 那必须写进它的每一次陈述。**
POSITIVE CTRL   种一个强 Δ -> Δ 方向的占比必须明显上升。
NEGATIVE CTRL   随机 `x`(打乱稀有度标签)的占比 = 基线。
IMPOSSIBLE      "占比小"不等于"不重要" —— 一个小而稳的方向可以比一个大而杂的方向更有解释力。
                本轮只量**大小**,不判**重要性**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
KEEP0=(np.isfinite(V0).sum(1)>=8)

def demean_np(Aa,iters=200,tol=1e-10):
    D=np.where(np.isfinite(Aa),Aa,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D

def shares(Aa, x_fixed=None, seed=0):
    D=demean_np(Aa); W=np.isfinite(D)&KEEP0[:,None]
    Z=np.where(W,D,0.0); tot=float((Z*Z).sum())
    # ① 固定 x 的秩一(= 人内斜率):每人最优载荷 u_i = <D_i, x_i>/<x_i,x_i>
    def fixed_x_share(x):
        Xc=W*(x[None,:]-np.where(W.sum(1,keepdims=True)>0,
              (W*x[None,:]).sum(1,keepdims=True)/np.maximum(W.sum(1,keepdims=True),1),0.0))
        num=(Z*Xc).sum(1); den=(Xc*Xc).sum(1)
        ss=float((num[den>1e-12]**2/den[den>1e-12]).sum())
        return ss/tot
    s_delta=fixed_x_share(rar0-rar0.mean())
    # ② 最佳秩一(x 自由):带缺失的交替最小二乘
    rng=np.random.default_rng(seed); x=rng.standard_normal(Mc)
    for _ in range(200):
        Xc=W*x[None,:]; den=(Xc*Xc).sum(1); u=np.where(den>1e-12,(Z*Xc).sum(1)/np.maximum(den,1e-12),0.0)
        Uc=W*u[:,None]; den2=(Uc*Uc).sum(0); x=np.where(den2>1e-12,(Z*Uc).sum(0)/np.maximum(den2,1e-12),0.0)
        if np.linalg.norm(x)>0: x=x/np.linalg.norm(x)
    Xc=W*x[None,:]; den=(Xc*Xc).sum(1); num=(Z*Xc).sum(1)
    s_best=float((num[den>1e-12]**2/den[den>1e-12]).sum())/tot
    # ③ 随机 x 基线
    s_rand=[fixed_x_share(rng.permutation(rar0)-rar0.mean()) for _ in range(20)]
    return s_delta,s_best,float(np.mean(s_rand)),float(np.std(s_rand)),tot,x

sd_,sb_,sr_,sr_sd,tot,xbest=shares(V0)
print(f"残差总平方和(观测格,KEEP 内)= {tot:,.0f}\n")
print(f"  ① **Δ 方向**(x = 稀有度)      占 **{100*sd_:.2f}%**")
print(f"  ② 最佳秩一方向(x 自由)        占 **{100*sb_:.2f}%**")
print(f"  ③ 随机 x 基线(20 次)          占 {100*sr_:.2f}% ± {100*sr_sd:.2f}%")
print(f"\n  Δ 方向 / 最佳秩一 = {sd_/sb_:.3f}")
print(f"  corr(最佳秩一的 x, 稀有度) = {np.corrcoef(xbest,rar0)[0,1]:+.4f}")

rng=np.random.default_rng(20260803)
u=np.abs(rng.standard_normal(N))+0.5
Vp=V0+8.0*np.outer(u,rar0-rar0.mean())*np.isfinite(V0)
sd_p,sb_p,_,_,_,_=shares(Vp,seed=3)
print(f"\n正对照(种强 Δ):Δ 方向占比 {100*sd_:.2f}% -> **{100*sd_p:.2f}%**")

T=pd.DataFrame([dict(direction='Δ(x=稀有度)',share=sd_),
                dict(direction='最佳秩一(x 自由)',share=sb_),
                dict(direction='随机 x 基线',share=sr_)])
check_columns(T,'R266'); T.to_csv(pathlib.Path(__file__).parent/'results'/'shares.csv',index=False)
g=Gate('Δ 在残差里有多大')
g.asserted('正对照:种强 Δ 后 Δ 方向占比必须明显上升',sd_p>sd_+0.02,f"{100*sd_:.2f}% -> {100*sd_p:.2f}%")
g.offset_control('Δ 方向占比 vs 随机 x 基线',sd_,sr_,sr_sd,
                 null_kind='把稀有度标签打乱后的同型方向 —— 不是零假设,'
                           '是"一个任意的题目权重向量本来能占多少"')
g.asserted('⚠ 占比小 ≠ 不重要 —— 本轮只量大小,不判重要性',True,
           '一个小而稳的方向可以比一个大而杂的方向更有解释力')
g.asserted('注册的 kill:Δ 方向占比 < 最佳秩一的 1/5 -> 它是残差里一个次要方向',
           sd_ < sb_/5, f"Δ {100*sd_:.2f}% vs 最佳秩一 {100*sb_:.2f}% -> 比值 {sd_/sb_:.3f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 原始比值被共同基线抬高了,必须扣掉基线再比 --------------------------------
# 随机方向本来就占 8.43% —— 那是"任意一个题目权重向量"的免费额度。
# 所以 `Δ/最佳 = 0.721` 是虚高的;要比的是**超出基线的那部分**。
print("\n---- 扣掉基线之后 ----")
ex_d=sd_-sr_; ex_b=sb_-sr_
print(f"  超出基线:Δ 方向 {100*ex_d:.2f} 个百分点 · 最佳秩一 {100*ex_b:.2f} 个百分点")
print(f"  **净比值 = {ex_d/ex_b:.3f}**(原始 {sd_/sb_:.3f})")
print(f"  -> 稀有度这个权重向量,拿到了'最好的那个方向'超额解释力的 **{100*ex_d/ex_b:.0f}%**,"
      f"而随机权重拿到 0%")
g2=Gate('扣掉基线之后,Δ 方向还占多少')
g2.asserted('基线本身可分辨地非零(所以必须扣)',sr_>3*sr_sd,f"随机基线 {100*sr_:.2f}% ± {100*sr_sd:.2f}%")
g2.no_sign_crossing('三个占比同号',[sd_,sb_,sr_])
g2.asserted('原始比值被共同基线抬高',sd_/sb_ > ex_d/ex_b,f"原始 {sd_/sb_:.3f} vs 净 {ex_d/ex_b:.3f}")
g2.resolvable('Δ 方向超出基线的部分',float(ex_d),float(sr_sd))
g2.asserted('净比值 > 1/5 -> Δ 不是残差里一个次要方向',ex_d/ex_b>0.2,f"{ex_d/ex_b:.3f}")
print(g2)
