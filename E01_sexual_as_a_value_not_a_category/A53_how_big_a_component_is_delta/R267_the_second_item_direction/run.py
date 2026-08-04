import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A53 R267 -- 起始年龄矩阵里,有没有第二个稳定的题目维度

`#221c`:最佳秩一方向与稀有度只有 **+0.2935** 的相关 —— **六成不是稀有度**。
**而那个方向的题目权重是可以直接看的**(`#203` 对多选块做过同样的事)。

ESTIMAND        ① 带缺失 ALS 解出的题目权重 `x`(31 个类别),排出来看;
                ② `x` 对稀有度回归后的**残差方向**,在**两半人**上复不复现(`#204b` 同款)。
KILL            **若残差方向在两半人上复现(|r| 明显高于置换零)->
                起始年龄矩阵里存在第二个稳定的题目维度,而本项目从 `#128` 起只用了第一个。**
⚠ 符号不定       ALS 的 `x` 符号任意 -> 只能判 **|r|**,并对着**同样取绝对值的零**比(`#204b` 的教训)。
NEGATIVE CTRL   题内跨人置换后走同一条管道 -> 零的 |r|(它也取绝对值)。
POSITIVE CTRL   种一个已知的第二题目方向 -> 复现 |r| 必须明显上升。
IMPOSSIBLE      31 个类别 · 每半约 5,000 人 —— 权重本身噪声大,|r| 是复现度的下界。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

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
def als_x(Aa, rows=None, seed=0, iters=300):
    D=demean_np(Aa); W=np.isfinite(D)&KEEP0[:,None]
    if rows is not None:
        keep=np.zeros(N,bool); keep[rows]=True; W=W&keep[:,None]
    Z=np.where(W,D,0.0)
    rng=np.random.default_rng(seed); x=rng.standard_normal(Mc)
    for _ in range(iters):
        Xc=W*x[None,:]; den=(Xc*Xc).sum(1); u=np.where(den>1e-12,(Z*Xc).sum(1)/np.maximum(den,1e-12),0.0)
        Uc=W*u[:,None]; den2=(Uc*Uc).sum(0); x=np.where(den2>1e-12,(Z*Uc).sum(0)/np.maximum(den2,1e-12),0.0)
        n=np.linalg.norm(x)
        if n>0: x=x/n
    return x
def resid_on_rar(x):
    X=np.c_[np.ones(Mc),rar0-rar0.mean()]
    return x-X@np.linalg.lstsq(X,x,rcond=None)[0]

x_all=als_x(V0)
if np.corrcoef(x_all,rar0)[0,1]<0: x_all=-x_all      # 符号锚到稀有度正向,只为展示
r_all=resid_on_rar(x_all)
check_residualized(r_all,rar0-rar0.mean(),'R267 x 对稀有度的残差')
print(f"corr(x, 稀有度) = {np.corrcoef(x_all,rar0)[0,1]:+.4f}\n")
o=np.argsort(x_all)
print("题目权重最**高**的 6 个类别:")
for j in o[::-1][:6]: print(f"  w {x_all[j]:>+7.3f}  稀有度 {rar0[j]:.2f}  {str(ons[j])[:56]}")
print("题目权重最**低**的 6 个类别:")
for j in o[:6]: print(f"  w {x_all[j]:>+7.3f}  稀有度 {rar0[j]:.2f}  {str(ons[j])[:56]}")
o2=np.argsort(r_all)
print("\n**去掉稀有度之后**,残差权重最高/最低的各 4 个:")
for j in o2[::-1][:4]: print(f"  高 r {r_all[j]:>+7.3f}  稀有度 {rar0[j]:.2f}  {str(ons[j])[:52]}")
for j in o2[:4]:       print(f"  低 r {r_all[j]:>+7.3f}  稀有度 {rar0[j]:.2f}  {str(ons[j])[:52]}")

def half_rep(Aa,seed):
    rng=np.random.default_rng(seed); p=rng.permutation(np.flatnonzero(KEEP0)); h=len(p)//2
    xa=als_x(Aa,rows=p[:h],seed=seed); xb=als_x(Aa,rows=p[h:2*h],seed=seed+1)
    return (abs(float(np.corrcoef(xa,xb)[0,1])),
            abs(float(np.corrcoef(resid_on_rar(xa),resid_on_rar(xb))[0,1])))
real=[half_rep(V0,3000+s) for s in range(8)]
nul =[half_rep(perm_null(V0,np.random.default_rng(4000+s)),5000+s) for s in range(8)]
rf=float(np.mean([a for a,_ in real])); rr_=float(np.mean([b for _,b in real]))
nf=float(np.mean([a for a,_ in nul]));  nr=float(np.mean([b for _,b in nul]))
sdr=float(np.std([b for _,b in real])); sdn=float(np.std([b for _,b in nul]))
print(f"\n两半复现 |r|:整条 x 真实 {rf:.3f} vs 零 {nf:.3f}")
print(f"            **去稀有度残差 真实 {rr_:.3f} ± {sdr:.3f} vs 零 {nr:.3f} ± {sdn:.3f}**")

rng=np.random.default_rng(20260803)
w2=rng.standard_normal(Mc); w2=resid_on_rar(w2); w2/=np.linalg.norm(w2)
u2=rng.standard_normal(N)
Vp=V0+6.0*np.outer(u2,w2)*np.isfinite(V0)
pl=[half_rep(Vp,7000+s) for s in range(4)]
rp=float(np.mean([b for _,b in pl]))
print(f"正对照(种一个已知的第二题目方向):残差复现 |r| {rr_:.3f} -> **{rp:.3f}**")

T=pd.DataFrame(dict(cat_q=[str(c)[:52] for c in ons],rarity=rar0,w=x_all,w_resid=r_all))
check_columns(T,'R267'); T.to_csv(pathlib.Path(__file__).parent/'results'/'second_dir.csv',index=False)
g=Gate('有没有第二个稳定的题目维度')
g.asserted('正对照:种一个已知第二方向,残差复现必须上升',rp>rr_+0.10,f"{rr_:.3f} -> {rp:.3f}")
g.asserted('⚠ ALS 的 x 符号任意 -> 只判 |r|,而零也取绝对值(`#204b`)',True,
           f"零的残差 |r| = {nr:.3f},不是 0")
g.negative_control('题内跨人置换的残差复现 |r|',nr,rr_,null_spread=sdn,
                   null_kind='题内跨人置换后走同一条 ALS 管道,同样取绝对值')
g.resolvable('残差复现 |r| 减去零',float(rr_-nr),float(np.hypot(sdr,sdn)))
g.asserted('注册的 kill:残差方向在两半人上复现 -> 存在第二个稳定的题目维度',
           rr_-nr>2*np.hypot(sdr,sdn),f"真实 {rr_:.3f} vs 零 {nr:.3f},差 {rr_-nr:+.3f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
